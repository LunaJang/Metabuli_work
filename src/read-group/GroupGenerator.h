
#ifndef GROUP_GENERATOR_H
#define GROUP_GENERATOR_H

#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <set>
#include <cassert>
#include <thread>
#include <atomic>
#include <sys/sysinfo.h>
#include <sys/resource.h>
#include <sys/statvfs.h>
#include <cstdio>
#include <cstdlib>
#include <cerrno>
#include <cstring>
#include <memory>
#include <mutex>
#include <algorithm>
#include <zstd.h>
#include "IndexCreator.h"
#include "SeqIterator.h"
#include "NcbiTaxonomy.h"
#include "LocalParameters.h"
#include "QueryIndexer.h"
#include "ReducedKmerMatcher.h"
#include "KmerExtractor.h"
#include "KSeqWrapper.h"
#include "DeltaIdxReader.h"

#define BufferSize 16'777'216 //16 * 1024 * 1024 // 16 M
using namespace std;

struct Relation {
    uint32_t id1;
    uint32_t id2;
    uint16_t weight;

    Relation(uint32_t a = 0, uint32_t b = 0): id1(a), id2(b), weight(0) {}
    Relation(uint32_t a, uint32_t b, uint16_t w): id1(a), id2(b), weight(w) {}

    static bool compare(const Relation& a, const Relation& b) {
        if (a.id1 != b.id1) return a.id1 < b.id1;
        return a.id2 < b.id2;
    }

    bool operator<(const Relation& other) const {
        if (id1 != other.id1) return id1 < other.id1;
        return id2 < other.id2;
    }

    bool operator==(const Relation& other) const {
        return id1 == other.id1 && id2 == other.id2;
    }
};

// Relation is 12 bytes but only 10 carry payload; the compiler never initialises the trailing
// two. They are written to disk verbatim, so leaving them undefined makes the intermediate files
// differ byte for byte between runs even when every record is identical -- which is exactly what
// happened to relations.bin once the surrounding code changed and the stack slot stopped being
// reused. Build every record that reaches a file through here.
// The void* cast is the documented way to memset a class type without -Wclass-memaccess.
inline Relation makeRelation(uint32_t id1, uint32_t id2, uint16_t weight) {
    Relation rel;
    std::memset(static_cast<void *>(&rel), 0, sizeof(rel));
    rel.id1 = id1;
    rel.id2 = id2;
    rel.weight = weight;
    return rel;
}

// Sequential Relation streams for the subGraph_* files ONLY. relations_* keeps using
// ReadBuffer/WriteBuffer<Relation>: giving the two a different type is the point, so that a
// change to the subgraph representation cannot silently reach the files makeGroups reads.
//
// On-disk layout: one zstd stream per file, nothing else. Records go in as raw Relation bytes,
// so the two padding bytes of every 12-byte struct are zeros the compressor removes for free.
// Streaming (not one-shot) because a flush can carry 200 M records -- one-shot would need a
// 2.4 GB input and a matching output buffer, while the streaming API bounds both.
//
// Every failure is fatal. saveSubGraphToFile used to print to cerr and continue, which let a
// 21-hour run finish with a zero exit code and an empty result after 70,399 writes had failed.
static const int RELATION_ZSTD_LEVEL = 1;             // ~500 MB/s/core; emit runs near 100 MB/s
static const size_t RELATION_WRITE_ELEMS = 65'536;    // 768 KB of records per compressor feed

static inline void relationStreamFatal(const std::string & path, const char * action) {
    std::cerr << "Error " << action << " file: " << path
              << " (" << strerror(errno) << ")" << std::endl;
    exit(EXIT_FAILURE);
}

static inline void relationStreamZstdFatal(const std::string & path, const char * action,
                                           size_t code) {
    std::cerr << "Error " << action << " file: " << path
              << " (zstd: " << ZSTD_getErrorName(code) << ")" << std::endl;
    exit(EXIT_FAILURE);
}

// Fixed per-stream cost of the codec, on top of the caller's record buffer. Charged to the
// merge memory budget in getMergeBufferElems so the disk representation and the memory
// representation are not silently conflated.
inline size_t relationStreamFixedBytes() {
    return ZSTD_DStreamInSize();
}

struct ZstdCStreamDeleter {
    void operator()(ZSTD_CStream * s) const { ZSTD_freeCStream(s); }
};

struct ZstdDStreamDeleter {
    void operator()(ZSTD_DStream * s) const { ZSTD_freeDStream(s); }
};

class RelationWriter {
public:
    explicit RelationWriter(const std::string & filePath)
        : path(filePath), fp(fopen(filePath.c_str(), "wb")), cstream(ZSTD_createCStream()),
          records(), outBuf(ZSTD_CStreamOutSize()), bytesOnDisk(0) {
        if (fp == nullptr) {
            relationStreamFatal(path, "opening");
        }
        if (!cstream) {
            relationStreamFatal(path, "allocating a compressor for");
        }
        const size_t init = ZSTD_initCStream(cstream.get(), RELATION_ZSTD_LEVEL);
        if (ZSTD_isError(init)) {
            relationStreamZstdFatal(path, "initializing the compressor for", init);
        }
        records.reserve(RELATION_WRITE_ELEMS);
    }

    RelationWriter(const RelationWriter &) = delete;
    RelationWriter & operator=(const RelationWriter &) = delete;

    ~RelationWriter() { finish(); }

    void write(const Relation & rel) {
        records.push_back(rel);
        if (records.size() >= RELATION_WRITE_ELEMS) {
            feedCompressor();
        }
    }

    // Idempotent so the destructor can act as a safety net after an explicit call.
    void finish() {
        if (fp == nullptr) {
            return;
        }
        feedCompressor();

        size_t remaining = 0;
        do {
            ZSTD_outBuffer out = { outBuf.data(), outBuf.size(), 0 };
            remaining = ZSTD_endStream(cstream.get(), &out);
            if (ZSTD_isError(remaining)) {
                relationStreamZstdFatal(path, "finishing", remaining);
            }
            writeOut(out.pos);
        } while (remaining != 0);

        FILE * const closing = fp;
        fp = nullptr;
        if (fclose(closing) != 0) {
            relationStreamFatal(path, "closing");
        }
    }

    size_t compressedBytes() const { return bytesOnDisk; }

private:
    void feedCompressor() {
        if (records.empty()) {
            return;
        }
        ZSTD_inBuffer in = { records.data(), records.size() * sizeof(Relation), 0 };
        while (in.pos < in.size) {
            ZSTD_outBuffer out = { outBuf.data(), outBuf.size(), 0 };
            const size_t code = ZSTD_compressStream(cstream.get(), &out, &in);
            if (ZSTD_isError(code)) {
                relationStreamZstdFatal(path, "compressing into", code);
            }
            writeOut(out.pos);
        }
        records.clear();
    }

    void writeOut(size_t bytes) {
        if (bytes == 0) {
            return;
        }
        if (fwrite(outBuf.data(), 1, bytes, fp) != bytes) {
            relationStreamFatal(path, "writing to");
        }
        bytesOnDisk += bytes;
    }

    std::string path;
    FILE * fp;
    std::unique_ptr<ZSTD_CStream, ZstdCStreamDeleter> cstream;
    std::vector<Relation> records; // reserve() fixes the capacity, so push_back never reallocates
    std::vector<char> outBuf;
    size_t bytesOnDisk;
};

class RelationReader {
public:
    RelationReader(const std::string & filePath, size_t blockElems)
        : path(filePath), fp(fopen(filePath.c_str(), "rb")), dstream(ZSTD_createDStream()),
          inBuf(ZSTD_DStreamInSize()), inFilled(0), inPos(0),
          outBytes((blockElems == 0 ? 1 : blockElems) * sizeof(Relation)),
          outFilled(0), outPos(0), sourceDrained(false) {
        if (fp == nullptr) {
            relationStreamFatal(path, "opening");
        }
        if (!dstream) {
            relationStreamFatal(path, "allocating a decompressor for");
        }
        const size_t init = ZSTD_initDStream(dstream.get());
        if (ZSTD_isError(init)) {
            relationStreamZstdFatal(path, "initializing the decompressor for", init);
        }
    }

    RelationReader(const RelationReader &) = delete;
    RelationReader & operator=(const RelationReader &) = delete;

    ~RelationReader() {
        if (fp != nullptr) {
            fclose(fp);
        }
    }

    // Returns a default Relation at end of stream, matching ReadBuffer<Relation>'s contract so
    // the callers' `if (!(next == Relation()))` termination test still holds. Relation::operator==
    // ignores weight, so the sentinel is really (id1,id2) == (0,0) -- an edge that cannot exist,
    // because writeKmers globalizes query IDs as 1-based (sequenceID += processedReadCnt).
    Relation getNext() {
        if (outFilled - outPos < sizeof(Relation) && !fill()) {
            return Relation();
        }
        Relation rel;
        memcpy(&rel, outBytes.data() + outPos, sizeof(Relation));
        outPos += sizeof(Relation);
        return rel;
    }

private:
    // Decompress until at least one whole record is buffered. A record can straddle two
    // decompressor outputs, so the leftover tail is carried to the front instead of dropped.
    bool fill() {
        const size_t leftover = outFilled - outPos;
        if (leftover > 0 && outPos > 0) {
            memmove(outBytes.data(), outBytes.data() + outPos, leftover);
        }
        outFilled = leftover;
        outPos = 0;

        while (outFilled < sizeof(Relation)) {
            if (inPos >= inFilled) {
                if (sourceDrained) {
                    return false;
                }
                inFilled = fread(inBuf.data(), 1, inBuf.size(), fp);
                inPos = 0;
                if (inFilled == 0) {
                    sourceDrained = true;
                    return false;
                }
            }

            ZSTD_inBuffer in = { inBuf.data(), inFilled, inPos };
            ZSTD_outBuffer out = { outBytes.data(), outBytes.size(), outFilled };
            const size_t code = ZSTD_decompressStream(dstream.get(), &out, &in);
            if (ZSTD_isError(code)) {
                relationStreamZstdFatal(path, "decompressing", code);
            }
            inPos = in.pos;
            outFilled = out.pos;
        }
        return true;
    }

    std::string path;
    FILE * fp;
    std::unique_ptr<ZSTD_DStream, ZstdDStreamDeleter> dstream;
    std::vector<char> inBuf;
    size_t inFilled;
    size_t inPos;
    std::vector<char> outBytes;
    size_t outFilled;
    size_t outPos;
    bool sourceDrained;
};

class DisjointSet {
public:
    std::vector<uint32_t> parent;
    std::vector<uint32_t> rank;
    std::vector<bool> grouped;

    DisjointSet(size_t numQuery) 
        : parent(numQuery + 1), rank(numQuery + 1), grouped(numQuery + 1) {
        for (size_t i = 1; i <= numQuery; ++i) {
            parent[i] = i;
            rank[i] = 0;
            grouped[i] = false;
        }
    }

    uint32_t find(uint32_t element) {
        if (parent[element] != element) {
            parent[element] = find(parent[element]);
        }
        return parent[element];
    }

    void unionSets(uint32_t elem1, uint32_t elem2) {
        uint32_t root1 = find(elem1);
        uint32_t root2 = find(elem2);

        grouped[elem1] = true;
        grouped[elem2] = true;

        if (root1 != root2) {
            if (rank[root1] < rank[root2]) {
                parent[root1] = root2;
            } else if (rank[root1] > rank[root2]) {
                parent[root2] = root1;
            } else {
                if (root1 < root2) {
                    parent[root2] = root1;
                    rank[root1]++;
                }
                else{
                    parent[root1] = root2;
                    rank[root2]++;
                }
            }
        }
    }

    void flatten() {
        for (size_t i = 1; i < parent.size(); ++i) {
            parent[i] = find(i);
        }
    }

    DisjointSet& operator+=(const DisjointSet& rhs) {
        for (size_t i = 1; i < parent.size(); ++i) {
            uint32_t p = rhs.parent[i];
            if (p != i) unionSets(i, p);
        }
        return *this;
    }

    friend DisjointSet operator+(DisjointSet lhs, const DisjointSet& rhs) {
        lhs += rhs;
        return lhs;
    }
};

// Saturating accumulation into a uint16 edge weight: clamps at UINT16_MAX
// instead of wrapping. Relation::weight is uint16_t, so both the per-subgraph
// count and the cross-subgraph sum can exceed 65535; a wrapped weight would
// look small and be dropped by the Phase 1/2 threshold comparisons.
static inline void addSat16(uint16_t & dst, uint32_t add) {
    const uint32_t sum = static_cast<uint32_t>(dst) + add;
    dst = (sum > UINT16_MAX) ? UINT16_MAX : static_cast<uint16_t>(sum);
}

// Memory usable for buffers.
//
// Linux uses /proc/meminfo's MemAvailable, not sysinfo's freeram: freeram excludes the page
// cache even though the cache is reclaimable, so during a multi-hundred-GB write it collapses.
// That collapse used to shrink the flush threshold and split the graph into tens of thousands
// of tiny subgraph files, which then exhausted the process's file descriptors.
// macOS uses hw.memsize, which is TOTAL memory and therefore an over-estimate.
inline size_t getAvailableMemoryBytes() {
#if defined(__linux__)
    // Line-by-line: /proc/meminfo rows do not all carry a unit token, so streaming with >>
    // would read the next row's key as this row's unit.
    std::ifstream meminfo("/proc/meminfo");
    std::string line;
    while (std::getline(meminfo, line)) {
        unsigned long long kb = 0;
        if (sscanf(line.c_str(), "MemAvailable: %llu kB", &kb) == 1) {
            return static_cast<size_t>(kb) * 1024;
        }
    }
    // Pre-3.14 kernels have no MemAvailable; free + reclaimable buffers is the closest stand-in.
    struct sysinfo info;
    sysinfo(&info);
    return (size_t)(info.freeram + info.bufferram) * info.mem_unit;
#elif defined(__APPLE__)
    int64_t freeMemory;
    size_t len = sizeof(freeMemory);
    sysctlbyname("hw.memsize", &freeMemory, &len, nullptr, 0);
    return (size_t)freeMemory; // 근사치
#else
    return 8ULL * 1024 * 1024 * 1024; // fallback 8GB
#endif
}

// Memory budget: what the OS can spare, capped by the user's declared --max-ram. The rest of
// the pipeline already budgets against --max-ram (Buffer::calculateBufferSize); these two
// helpers used to ignore it and trust an instantaneous OS reading alone.
inline size_t getMemoryBudgetBytes(int maxRamGiB) {
    const size_t availableBytes = getAvailableMemoryBytes();
    if (maxRamGiB <= 0) {
        return availableBytes;
    }
    const size_t declared = static_cast<size_t>(maxRamGiB) * 1024 * 1024 * 1024;
    return std::min(availableBytes, declared);
}

// Free space at `path`, or 0 when it cannot be measured. Same platform-guard shape as
// getAvailableMemoryBytes: a value we cannot trust is reported as "unknown", not guessed.
inline size_t getFreeDiskBytes(const std::string & path) {
#if defined(__linux__) || defined(__APPLE__)
    struct statvfs info;
    if (statvfs(path.c_str(), &info) != 0) {
        return 0;
    }
    return static_cast<size_t>(info.f_bavail) * static_cast<size_t>(info.f_frsize);
#else
    (void) path;
    return 0;
#endif
}

// Returned when no cap applies -- either the user asked for none or free space is unknown.
static const size_t TMP_DISK_UNLIMITED = SIZE_MAX;

// Total capacity of the filesystem holding `path`, or 0 when it cannot be measured.
inline size_t getTotalDiskBytes(const std::string & path) {
#if defined(__linux__) || defined(__APPLE__)
    struct statvfs info;
    if (statvfs(path.c_str(), &info) != 0) {
        return 0;
    }
    return static_cast<size_t>(info.f_blocks) * static_cast<size_t>(info.f_frsize);
#else
    (void) path;
    return 0;
#endif
}

// Space the run refuses to consume. Folding bounds the footprint but cannot bound it below the
// merged graph, so a hard floor is still needed -- this is what actually keeps the run from
// filling somebody else's filesystem. 2% of the volume, clamped so it is neither trivial on a
// laptop nor wasteful on a 100 TB array.
inline size_t getDiskReserveBytes(const std::string & path) {
    const size_t MIN_RESERVE = 64ULL * 1024 * 1024;
    const size_t MAX_RESERVE = 1024ULL * 1024 * 1024;
    const size_t total = getTotalDiskBytes(path);
    if (total == 0) {
        return MIN_RESERVE;
    }
    return std::max(MIN_RESERVE, std::min(MAX_RESERVE, total / 50));
}

// Disk counterpart of getMemoryBudgetBytes. Note what this value is NOT: it is not a
// performance target. Anchoring it at a small fraction of free space would make a laptop fold
// constantly while leaving a 100 TB filesystem unprotected, and would make the same data behave
// differently depending on what else happens to sit on the disk. 80% is a safety ceiling: a run
// that fits under it never folds early and takes the same path it takes today.
inline size_t getTmpDiskBudgetBytes(const std::string & outDir, int maxTmpDiskMiB) {
    if (maxTmpDiskMiB > 0) {
        return static_cast<size_t>(maxTmpDiskMiB) * 1024 * 1024;
    }
    const size_t freeBytes = getFreeDiskBytes(outDir);
    if (freeBytes == 0) {
        return TMP_DISK_UNLIMITED;
    }
    return freeBytes / 100 * 80;
}

inline size_t getRelationThreshold(int numThreads, int maxRamGiB) {
    const size_t availableBytes = getMemoryBudgetBytes(maxRamGiB);

    const double safetyFactor = 0.6;
    const size_t bytesPerEntry = 48; // unordered_map node overhead

    size_t threshold = (size_t)(availableBytes * safetyFactor)
                       / (numThreads * bytesPerEntry);

    const size_t MIN_THRESHOLD = 1'000'000;
    const size_t MAX_THRESHOLD = 200'000'000;
    return std::max(MIN_THRESHOLD, std::min(threshold, MAX_THRESHOLD));
}

// Width of one id range in the relations_* routing scheme.
inline size_t getRouteRangeSize(size_t processedReadCnt, int numThreads) {
    const size_t threads = static_cast<size_t>(numThreads);
    return (processedReadCnt > threads) ? (processedReadCnt / threads) : processedReadCnt;
}

// Which relations_{r}.bin an edge belongs to. Depends only on (id1, id2), which is what lets
// the merge be split by route: a pair can never appear under two different routes.
// Routes 0..numThreads-1 and numThreads+1..2*numThreads are each read by a single union-find
// thread; numThreads is the cross bucket that spans partitions.
inline size_t routeOf(uint32_t id1, uint32_t id2, int numThreads, size_t rangeSize) {
    const uint32_t threads = static_cast<uint32_t>(numThreads);
    if (id1 % threads == id2 % threads) {
        return static_cast<size_t>(id1 % threads);
    }
    // rangeSize is 0 only when there are no reads at all, in which case there are no edges
    // either; guard anyway so the helper can never divide by zero.
    if (rangeSize != 0 && id1 / rangeSize == id2 / rangeSize) {
        return static_cast<size_t>(id1 / rangeSize) + static_cast<size_t>(numThreads);
    }
    return static_cast<size_t>(numThreads);
}

// How many shards a route is split into for merging. Only the cross bucket needs it: the
// routing rule sends ~88% of all edges there (measured), so merging it as one unit caps the
// whole parallel merge at ~1.13x no matter how many threads are available.
inline size_t shardsForRoute(size_t route, int numThreads) {
    if (numThreads < 1) {
        return 1;
    }
    return (route == static_cast<size_t>(numThreads)) ? static_cast<size_t>(numThreads) : 1;
}

// Which shard of its route an edge belongs to. Keyed on id1 ALONE -- keying on anything
// derived from both ids would let one pair land in two shards and split its weight.
inline size_t shardOf(uint32_t id1, size_t shardCnt) {
    if (shardCnt <= 1) {
        return 0;
    }
    return static_cast<size_t>(id1 % static_cast<uint32_t>(shardCnt));
}

// Soft cap on how many files this process may hold open.
inline size_t getOpenFileLimit() {
    struct rlimit lim;
    if (getrlimit(RLIMIT_NOFILE, &lim) == 0 && lim.rlim_cur != RLIM_INFINITY) {
        return static_cast<size_t>(lim.rlim_cur);
    }
    return 1024; // POSIX-typical default when the limit cannot be read
}

// How many subgraph files one merge pass may open at once. Bounded by the fd limit -- with a
// reserve for the relations_* writers, stdio, and the k-mer readers -- and by a cap that keeps
// each stream's read buffer large enough for sequential I/O. Anything beyond this is handled
// by merging in rounds rather than by opening more files.
// Descriptors one merge unit holds at its peak: `fanIn` readers plus the single file it is
// writing (its output, or a fold intermediate).
inline size_t mergeFdsPerUnit(size_t fanIn) {
    return fanIn + 1;
}

// Largest fan-in whose peak descriptor use fits the process limit.
//
// `concurrentMergers` units run at once, each holding mergeFdsPerUnit(fanIn). An earlier
// version reserved a flat slack and forgot the per-unit writer, which put the peak one
// descriptor over an 8192 limit at 16 threads and killed the merge with a bare
// "Error opening file". Only a fraction of the limit is spent so that descriptors held
// elsewhere -- and the writer this arithmetic now counts -- cannot tip it over.
inline size_t getMergeFanIn(int numThreads, size_t concurrentMergers) {
    (void) numThreads;
    const size_t FAN_IN_CAP = 512;
    const double utilisation = 0.75;
    const size_t mergers = std::max<size_t>(1, concurrentMergers);
    const size_t budget = static_cast<size_t>(getOpenFileLimit() * utilisation);
    // Solve mergers * (fanIn + 1) <= budget.
    const size_t perUnit = (budget > mergers) ? (budget / mergers - 1) : 0;
    return std::max<size_t>(2, std::min(FAN_IN_CAP, perUnit));
}

// m-distribution histogram: reads-per-k-mer bucketed by floor(log2(m)). m can reach the
// read count (tens of millions), so a flat array is not affordable per thread; log2 buckets
// give cap candidates at factor-2 spacing, which is enough to pick a disk budget.
static const size_t M_HIST_BUCKETS = 64;

static inline size_t mHistBucket(size_t m) {
    size_t bucket = 0;
    while ((m >> bucket) > 1 && bucket + 1 < M_HIST_BUCKETS) {
        ++bucket;
    }
    return bucket; // floor(log2(m)) for m >= 1; 0 for m == 0
}

// Floor for mergeGraph's per-stream read buffer. Reaching it means numOfGraph is so large
// that even the memory budget cannot give each stream a useful buffer; mergeGraph warns,
// and the real fix is multi-round merging (not implemented).
static const size_t MERGE_BUFFER_MIN_ELEMS = 4'096;

// Per-stream ReadBuffer size (in Relation elements) for the k-way subgraph merge.
// The merge opens one buffer per subGraph_* file, so a fixed per-stream size makes total
// memory scale with numOfGraph: the historical 1M elements is 12 MB per stream, i.e. 12 GB
// at 1000 streams. Budget the total against free memory and split it across streams.
inline size_t getMergeBufferElems(size_t numOfGraph, int maxRamGiB, size_t concurrentMergers) {
    const size_t MAX_ELEMS = 1'048'576; // 12 MB per stream; the historical fixed size
    if (numOfGraph == 0) {
        return MAX_ELEMS;
    }

    const double safetyFactor = 0.5;
    const size_t mergers = std::max<size_t>(1, concurrentMergers);
    const size_t budget = (size_t)(getMemoryBudgetBytes(maxRamGiB) * safetyFactor) / mergers;

    // Disk bytes and memory bytes are not the same thing here. Each compressed stream also
    // holds a fixed decompressor input window; charge that first, then split what is left
    // across the record buffers. Conflating the two is how the merge blew past its fd budget
    // once before (see research.md, 2026-08-06 (2)).
    const size_t fixedCost = relationStreamFixedBytes() * numOfGraph;
    const size_t forRecords = (budget > fixedCost) ? (budget - fixedCost) : 0;
    const size_t perStream = forRecords / (numOfGraph * sizeof(Relation));

    return std::max(MERGE_BUFFER_MIN_ELEMS, std::min(perStream, MAX_ELEMS));
}

class GroupGenerator {
protected:
    const LocalParameters & par;
    string commonKmerDB;
    string outDir;
    size_t matchPerKmer;
    int kmerFormat;
    int printLog;
    
    // Agents    
    GeneticCode * geneticCode = nullptr;
    QueryIndexer *queryIndexer = nullptr;
    KmerExtractor *kmerExtractor = nullptr;

    size_t numOfSplits = 0;
    size_t numOfGraph = 0;
    // Bytes the subGraph_* files actually occupy. saveSubGraphToFile runs on every emit
    // thread, hence the atomic. Compared against emittedEdgeCnt * sizeof(Relation) it gives
    // the compression ratio the run achieved.
    std::atomic<uint64_t> subGraphBytesOnDisk{0};

    // --- Incremental fold state -------------------------------------------------------------
    // The merge used to start only after every flush had been written, so peak disk was the
    // whole emit. Folding a completed prefix while emission continues bounds it instead.
    //
    // A flush index is handed out by counter.fetch_add BEFORE its files exist, so "written" is
    // tracked separately: emitDone[i] is set once all of flush i's unit files are closed, and
    // emitWatermark walks the contiguous prefix of those. Only that prefix may be folded --
    // folding across a hole would list a file that has not been written yet.
    std::mutex foldMutex;                                // guards the five members below
    std::vector<char> emitDone;
    size_t emitWatermark = 0;                            // flushes [0, emitWatermark) are complete
    size_t foldedUpTo = 0;                               // flushes [0, foldedUpTo) folded away
    size_t foldRounds = 0;
    std::vector<std::vector<std::string>> foldedOutputs; // [unit] -> folded file names

    // Live footprint of subGraph_* and its high-water mark. Signed because a fold subtracts its
    // inputs before adding its output.
    std::atomic<int64_t> subGraphLiveBytes{0};
    std::atomic<uint64_t> subGraphPeakBytes{0};

    // Emit progress, for projecting the total write volume. The [mhist] table cannot serve here:
    // it is only assembled after the parallel region joins. These two are updated as k-mers are
    // consumed, so a projection is available while the emit is still running.
    std::atomic<uint64_t> kmerValuesTotal{0};
    std::atomic<uint64_t> kmerValuesDone{0};
    std::atomic<int> volumeWarned{0};   // the projection warning is printed at most once

    std::vector<uint64_t> kmerBoundaries;
    bool boundariesInitialized = false;
    bool useOnlyTrueRelations = false; // for debug

    // Filter-surviving k-mers, accumulated across splits. Divided by the read count it gives
    // k-mers per read, which is the scale the Phase 1 edge-weight threshold is expressed in:
    // an edge weight is the number of k-mers two reads share, so weight/kmersPerRead is the
    // fraction of a read the two overlap by.
    size_t totalFilteredKmers = 0;

public:
    GroupGenerator(LocalParameters & par);

    ~GroupGenerator();

    void startGroupGeneration(const LocalParameters & par);
    
    void filterCommonKmers(Buffer<Kmer>& queryKmerBuffer,
                           Buffer<std::pair<uint32_t, uint32_t>> & matchBuffer,
                           const string & db="");

    void writeKmers(Buffer<Kmer>& queryKmerBuffer,
                    size_t processedReadCnt);

    std::vector<std::pair<size_t, size_t>> getKmerRanges(const Buffer<Kmer>& kmerBuffer,
                                                         size_t offset);

    void makeSubGraph(size_t processedReadCnt);
    
    // Writes one file per relations_* route (subGraph_{counter}_{route}) instead of one file
    // per flush, so the merge can process each route independently and in parallel.
    void saveSubGraphToFile(const unordered_map<uint64_t, uint16_t>& pair2weight,
                            const size_t counter_now,
                            size_t processedReadCnt);

    // No processedReadCnt: routing now happens in saveSubGraphToFile, which is where the
    // read count is needed. This function only folds and merges what is already partitioned.
    void mergeGraph(std::vector<uint64_t>& edgeWeightHist);

    // Merge one batch of sorted subgraph files into a single sorted file, summing the weights
    // of duplicate (id1,id2) pairs. Returns the number of records written. Used to fold many
    // subgraphs down to a mergeable fan-in; it does no partitioning or histogramming, which
    // only the final pass performs.
    size_t mergeSubGraphBatch(const std::vector<std::string>& inputs,
                              const std::string& output,
                              size_t bufElems);

    // Reduce `files` to at most maxFanIn entries by merging in rounds, deleting each batch's
    // inputs as soon as it is folded so peak disk grows by one batch rather than a full copy.
    // `route`/`shard` only name the intermediates: units fold concurrently, so a shared name
    // would let two of them write the same file and silently mix their edges.
    void reduceSubGraphFanIn(std::vector<std::string>& files, size_t maxFanIn,
                             size_t route, size_t shard);

    // Merge every subgraph belonging to one (route, shard) unit into `outPath`. Purely
    // sequential -- mergeGraph runs these concurrently, one unit per thread.
    // `files` comes from unitInputFiles(): incremental folding replaces a prefix of the flushes
    // with folded files, so the unit's inputs can no longer be derived from numOfGraph alone.
    void mergeUnit(size_t route, size_t shard, const std::string& outPath,
                   const std::vector<std::string>& files,
                   size_t bufElems, size_t maxFanIn,
                   std::vector<uint64_t>& histOut, size_t& mergedOut, size_t& ceilingOut);

    // --- Incremental fold ---------------------------------------------------------------
    // Units are numbered route-major with a fixed stride so the index is a pure function of
    // (route, shard); routes that carry one shard simply leave their extra slots unused.
    size_t shardStride() const {
        return shardsForRoute(static_cast<size_t>(par.threads), par.threads);
    }
    size_t unitCount() const {
        return (static_cast<size_t>(par.threads) * 2 + 1) * shardStride();
    }
    size_t unitIndexOf(size_t route, size_t shard) const {
        return route * shardStride() + shard;
    }
    std::string subGraphName(size_t flushIdx, size_t route, size_t shard) const {
        return outDir + "/subGraph_" + std::to_string(flushIdx) + "_" + std::to_string(route)
             + "_" + std::to_string(shard);
    }

    // Live-footprint bookkeeping. `delta` is positive for a file written, negative for one
    // removed; the high-water mark is what the [disk] summary reports.
    void noteSubGraphBytes(int64_t delta);

    // Called once every unit file of `flushIdx` is closed. Advances the contiguous watermark.
    void markEmitComplete(size_t flushIdx);

    // Two guards, checked at every flush.
    //  - projection: extrapolate the write volume from emit progress and warn once if it will
    //    not fit. Only a warning: folding routinely keeps the footprint far below the total,
    //    so a projected overflow does not mean the run must fail.
    //  - headroom: if the filesystem's actual free space falls to the reserve, stop. This one
    //    depends on no estimate, which is what makes it the real guarantee.
    void checkDiskHeadroom();

    // Fold the completed-but-unfolded prefix when the live footprint reaches `tmpDiskBudget`,
    // or when that prefix alone would exceed the merge's fan-in. Returns without doing anything
    // if another thread is already folding -- the caller just keeps emitting.
    void maybeFoldEmitted(size_t tmpDiskBudget, size_t maxFanIn);

    // Every file the merge must read for one unit: folded outputs first, then the flushes that
    // were never folded.
    std::vector<std::string> unitInputFiles(size_t route, size_t shard) const;

    // How many streams one unit contributes to the merge. Every unit folds on the same rounds
    // and keeps the same unfolded tail, so one unit's count answers for all of them. This --
    // not numOfGraph -- is what the fan-in and buffer budgets have to be sized against once
    // folding has replaced part of the flush range.
    size_t unitStreamCount() const {
        const size_t folded = foldedOutputs.empty() ? 0 : foldedOutputs[0].size();
        return folded + (numOfGraph - foldedUpTo);
    }

    static int kneeThreshold(const std::vector<uint64_t>& hist, int minWeight);

    void mergeGraph_one(size_t processedReadCnt);
    
    void makeGroups(int groupKmerThr,
                    size_t processedReadCnt,
                    unordered_map<uint32_t, unordered_set<uint32_t>>& groupInfo,
                    vector<uint32_t> &queryGroupInfo);
    // Phase 1.5: merge Phase-1 units that are joined by several independent weak edges.
    // A single weak edge (few shared k-mers) is indistinguishable from chance sharing across
    // repeats or conserved regions; several between the same pair of units are not. Counting
    // support at the unit level -- rather than testing node-level triangles -- is what makes
    // this affordable: two singletons can only ever have one edge between them, so every pair
    // worth counting involves at least one multi-read core group.
    void mergeBySupport(int coreThr,
                        int weakThr,
                        int minSupport,
                        size_t processedReadCnt,
                        std::unordered_map<uint32_t, std::unordered_set<uint32_t>> &groupInfo,
                        std::vector<uint32_t> &queryGroupInfo);


    void makeGroupsPhase2(int groupKmerThr,
                          size_t processedReadCnt,
                          const std::vector<bool>& isSingleton,
                          unordered_map<uint32_t, unordered_set<uint32_t>>& groupInfo,
                          vector<uint32_t>& queryGroupInfo);

    void saveGroupsToFile(const unordered_map<uint32_t, unordered_set<uint32_t>>& groupInfo,
                          const vector<uint32_t>& queryGroupInfo);
    
};


#endif // GROUP_GENERATOR_H