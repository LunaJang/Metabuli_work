
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
#include <cstdio>
#include <algorithm>
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

static inline bool keepEdgeGeo(uint16_t w, uint16_t tu, uint16_t tv) {
    // w >= sqrt(tu*tv)  <=>  w*w >= tu*tv
    return (uint64_t)w * (uint64_t)w >= (uint64_t)tu * (uint64_t)tv;
}

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
inline size_t getMergeFanIn(int numThreads) {
    const size_t FAN_IN_CAP = 512;
    const size_t reserved = static_cast<size_t>(numThreads) * 2 + 1 + 64; // relations_* + slack
    const size_t limit = getOpenFileLimit();
    const size_t byFd = (limit > reserved + 8) ? (limit - reserved) : 8;
    return std::max<size_t>(2, std::min(FAN_IN_CAP, byFd));
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
inline size_t getMergeBufferElems(size_t numOfGraph, int maxRamGiB) {
    const size_t MAX_ELEMS = 1'048'576; // 12 MB per stream; the historical fixed size
    if (numOfGraph == 0) {
        return MAX_ELEMS;
    }

    const double safetyFactor = 0.5;
    const size_t budget = (size_t)(getMemoryBudgetBytes(maxRamGiB) * safetyFactor);
    const size_t perStream = budget / (numOfGraph * sizeof(Relation));

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
    std::vector<uint64_t> kmerBoundaries;
    bool boundariesInitialized = false;
    bool useOnlyTrueRelations = false; // for debug

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
    
    void saveSubGraphToFile(const unordered_map<uint64_t, uint16_t>& pair2weight,
                            const size_t counter_now);

    void mergeGraph(size_t processedReadCnt, std::vector<uint64_t>& edgeWeightHist);

    // Merge one batch of sorted subgraph files into a single sorted file, summing the weights
    // of duplicate (id1,id2) pairs. Returns the number of records written. Used to fold many
    // subgraphs down to a mergeable fan-in; it does no partitioning or histogramming, which
    // only the final pass performs.
    size_t mergeSubGraphBatch(const std::vector<std::string>& inputs,
                              const std::string& output,
                              size_t bufElems);

    // Reduce `files` to at most maxFanIn entries by merging in rounds, deleting each batch's
    // inputs as soon as it is folded so peak disk grows by one batch rather than a full copy.
    void reduceSubGraphFanIn(std::vector<std::string>& files, size_t maxFanIn);

    static int otsuThreshold(const std::vector<uint64_t>& hist);

    static int kneeThreshold(const std::vector<uint64_t>& hist, int minWeight);

    void mergeGraph_one(size_t processedReadCnt);
    
    void computeNodeDegree(int groupKmerThr, 
                           size_t processedReadCnt, 
                           std::vector<uint32_t>& degree);
                           
    void computeGroupQuarterDegree(const std::vector<uint32_t>& queryGroupInfo,
                                  const std::vector<uint32_t>& degree,
                                  std::unordered_map<uint32_t, uint32_t>& groupQuarterDeg);
    
    void makeGroupsAdaptive(const vector<uint16_t>& nodeThr,
                            size_t processedReadCnt,
                            vector<uint32_t>& queryGroupInfo);         

    void makeGroups(int groupKmerThr,
                    size_t processedReadCnt,
                    unordered_map<uint32_t, unordered_set<uint32_t>>& groupInfo,
                    vector<uint32_t> &queryGroupInfo);

    void makeGroupsPhase2(int groupKmerThr,
                          size_t processedReadCnt,
                          const std::vector<bool>& isSingleton,
                          unordered_map<uint32_t, unordered_set<uint32_t>>& groupInfo,
                          vector<uint32_t>& queryGroupInfo);

    void saveGroupsToFile(const unordered_map<uint32_t, unordered_set<uint32_t>>& groupInfo,
                          const vector<uint32_t>& queryGroupInfo);
    
    uint16_t degreeToThr(uint32_t quarterDegree) const {
        float predCoverage = quarterDegree * 0.5f;
        float thr = predCoverage * 3.5f;
        return static_cast<uint16_t>(std::max(1.0f, std::min(thr, 150.0f))); 
    }

};


#endif // GROUP_GENERATOR_H