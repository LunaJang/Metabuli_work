#include "GroupGenerator.h"
#include "FileUtil.h"
#include "QueryIndexer.h"
#include "common.h"
#include "Kmer.h"
#include <sys/stat.h>
#include <queue>
#include <iomanip>
#include <functional>

// One stream head in mergeGraph's k-way merge: the relation plus the subGraph_* index it
// came from, so the merge knows which stream to refill after popping.
struct MergeHeapEntry {
    Relation rel;
    size_t stream;
};

// Min-heap ordering for MergeHeapEntry. std::priority_queue pops the "largest" element, so
// this comparator answers "does a come after b". Relation::operator< compares (id1,id2)
// only; equal keys tie-break on stream index to keep the pop order deterministic.
struct MergeHeapGreater {
    bool operator()(const MergeHeapEntry & a, const MergeHeapEntry & b) const {
        if (a.rel < b.rel) { return false; }
        if (b.rel < a.rel) { return true; }
        return a.stream > b.stream;
    }
};

using MergeHeap = std::priority_queue<MergeHeapEntry, std::vector<MergeHeapEntry>, MergeHeapGreater>;

// Human-readable byte size for disk-usage logs.
static std::string humanBytes(uint64_t bytes) {
    const char* unit[] = {"B", "KB", "MB", "GB", "TB"};
    double v = static_cast<double>(bytes);
    int i = 0;
    while (v >= 1024.0 && i < 4) { v /= 1024.0; ++i; }
    char buf[64];
    snprintf(buf, sizeof(buf), "%.2f %s", v, unit[i]);
    return std::string(buf);
}

// Sum sizes of existing files in `paths`, log "<label>: N files, X freed", then remove them.
// Non-existent paths are skipped (not counted).
static void reportAndRemoveFiles(const std::vector<std::string>& paths, const std::string& label) {
    uint64_t total = 0;
    size_t count = 0;
    for (const std::string& p : paths) {
        struct stat st;
        if (stat(p.c_str(), &st) == 0) {
            total += static_cast<uint64_t>(st.st_size);
            ++count;
            std::remove(p.c_str());
        }
    }
    std::cout << "[disk] " << label << ": " << count << " files, "
              << humanBytes(total) << " freed" << std::endl;
}

GroupGenerator::GroupGenerator(LocalParameters & par) : par(par) {
    commonKmerDB = par.filenames[1 + (par.seqMode == 2)];
    outDir       = par.filenames[2 + (par.seqMode == 2)];
    matchPerKmer = par.matchPerKmer;
    kmerFormat = par.kmerFormat;
    printLog = par.printLog;
    
    geneticCode = new GeneticCode(par.reducedAA == 1);
    queryIndexer = new QueryIndexer(par);
    queryIndexer->setKmerLen(12);
    kmerExtractor = new KmerExtractor(par, *geneticCode, kmerFormat);
}

GroupGenerator::~GroupGenerator() {
    delete queryIndexer;
    delete kmerExtractor;
    delete geneticCode;
}

void GroupGenerator::startGroupGeneration(const LocalParameters &par) {  
    Buffer<Kmer> queryKmerBuffer;
    Buffer<std::pair<uint32_t, uint32_t>> matchBuffer; // seq id, pos
    vector<Query> queryList;

    bool complete = false;
    size_t processedReadCnt = 0;
    size_t tries = 0;
    size_t totalSeqCnt = 0;

    // Extract k-mers from query sequences and compare them to target k-mer DB
    while (!complete) {
        tries++;

        // new code
        if (tries == 1) {
                cout << "Indexing query file ...";
        }
        queryIndexer->setBytesPerKmer(matchPerKmer);
        queryIndexer->indexQueryFile(processedReadCnt);
        const vector<QuerySplit> & queryReadSplit = queryIndexer->getQuerySplits();

        if (tries == 1) {
            totalSeqCnt = queryIndexer->getReadNum_1();
            cout << "Done" << endl;
            cout << "Total number of sequences: " << queryIndexer->getReadNum_1() << endl;
            cout << "Total read length: " << queryIndexer->getTotalReadLength() <<  "nt" << endl;
        }

        // Set up kseq
        KSeqWrapper* kseq1 = KSeqFactory(par.filenames[0].c_str());
        KSeqWrapper* kseq2 = nullptr;
        if (par.seqMode == 2) { kseq2 = KSeqFactory(par.filenames[1].c_str()); }

        // Move kseq to unprocessed reads
        for (size_t i = 0; i < processedReadCnt; i++) {
            kseq1->ReadEntry();
            if (par.seqMode == 2) { kseq2->ReadEntry(); }
        }

        for (size_t splitIdx = 0; splitIdx < queryReadSplit.size(); splitIdx++) {
            // Allocate memory for query list
            queryList.clear();
            queryList.resize(queryReadSplit[splitIdx].end - queryReadSplit[splitIdx].start);

            // Allocate memory for query k-mer buffer
            queryKmerBuffer.reallocateMemory(queryReadSplit[splitIdx].kmerCnt);
            queryKmerBuffer.init();
            matchBuffer.reallocateMemory(queryReadSplit[splitIdx].kmerCnt);
            matchBuffer.init();

            // Extract query k-mers
            kmerExtractor->extractQueryKmers(queryKmerBuffer,
                                             queryList,
                                             queryReadSplit[splitIdx],
                                             par,
                                             kseq1,
                                             kseq2);

            filterCommonKmers(queryKmerBuffer, matchBuffer, commonKmerDB);
            time_t t = time(nullptr);
            writeKmers(queryKmerBuffer, processedReadCnt);
            cout << "Writing query k-mer file: " << double(time(nullptr) - t) << " s" << endl;
            processedReadCnt += queryReadSplit[splitIdx].readCnt;
            cout << "The number of processed sequences: " << processedReadCnt << " (" << (double) processedReadCnt / (double) totalSeqCnt << ")" << endl;
            cout << "-----------------------------------" << endl;
        }
        delete kseq1;
        if (par.seqMode == 2) {
            delete kseq2;
        }
        if (processedReadCnt == totalSeqCnt) {
            complete = true;
        } 
    }   

    makeSubGraph(processedReadCnt);   

    unordered_map<uint32_t, unordered_set<uint32_t>> groupInfo;
    vector<uint32_t> queryGroupInfo;
    queryGroupInfo.resize(processedReadCnt + 1, 0); 

    if (printLog) {
        mergeGraph_one(processedReadCnt);
    } else {
        mergeGraph();

        // An edge weight counts the k-mers two reads share, so weight / kmersPerRead is the
        // fraction of a read the two overlap by. That makes the core threshold a property of read
        // geometry: one ratio means the same overlap on every dataset, which an absolute weight
        // does not. There is no fallback -- --min-overlap-ratio is validated to be > 0 before the
        // run starts, so a threshold that cannot be formed is a hard error, not a switch to some
        // other rule.
        const double kmersPerRead = (processedReadCnt > 0)
            ? static_cast<double>(totalFilteredKmers) / static_cast<double>(processedReadCnt)
            : 0.0;
        const int coreThr = static_cast<int>(par.minOverlapRatio * kmersPerRead + 0.5);
        if (coreThr < 2) {
            cerr << "Error: core threshold resolves to " << coreThr << " (overlap ratio "
                 << par.minOverlapRatio << " x " << kmersPerRead << " k-mers/read)." << endl;
            cerr << "       It must be at least 2, or every shared k-mer would form a core edge"
                 << " and no weak band would exist." << endl;
            cerr << "       Raise --min-overlap-ratio, or check that the common k-mer filter did"
                 << " not remove nearly every k-mer." << endl;
            exit(EXIT_FAILURE);
        }

        // Weak band lower bound, as a fraction of the core threshold. Keeping it proportional is
        // the point: an absolute 5 was 5/15 = 0.333 of the core on the species-inclusion benchmark
        // but 5/34 = 0.147 on CAMI2 marine, so marine's band was three times wider in absolute
        // terms and Phase 1.5 there took in far more chance links. The clamps only guard rounding:
        // the band has to contain at least one weight and stay below the core threshold.
        int weakThr = static_cast<int>(par.weakBandRatio * coreThr + 0.5f);
        if (weakThr < 1) { weakThr = 1; }
        if (weakThr >= coreThr) { weakThr = coreThr - 1; }

        cout << "Core threshold: " << coreThr << " (overlap ratio " << par.minOverlapRatio
             << " x " << kmersPerRead << " k-mers/read)" << endl;
        cout << "Weak band: (" << weakThr << ", " << coreThr << "] (ratio " << par.weakBandRatio
             << " x core " << coreThr << "); Phase 2 floor " << weakThr << endl;

        // Phase 1: form core groups with strong edges
        makeGroups(coreThr, processedReadCnt, groupInfo, queryGroupInfo);

        // Phase 1.5: merge units joined by several independent weak links
        mergeBySupport(coreThr, weakThr, par.mergeSupportRatio,
                       static_cast<size_t>(par.mergeMaxUnitReads),
                       processedReadCnt, groupInfo, queryGroupInfo);

        // Phase 2: link Phase-1 singletons among themselves with weak edges
        {
            std::vector<bool> isSingleton(processedReadCnt + 1, false);
            for (uint32_t i = 1; i <= processedReadCnt; i++) {
                if (queryGroupInfo[i] == 0) isSingleton[i] = true;
            }
            makeGroupsPhase2(weakThr, processedReadCnt, isSingleton, groupInfo, queryGroupInfo);
        }

        saveGroupsToFile(groupInfo, queryGroupInfo);

        // relations_*.bin are reparsed every refinement iteration but are no longer
        // needed once grouping is saved. Indices: relations_0 .. relations_{2*threads}.
        std::vector<std::string> relationFiles;
        for (int i = 0; i <= par.threads * 2; i++) {
            relationFiles.push_back(outDir + "/relations_" + std::to_string(i) + ".bin");
        }
        reportAndRemoveFiles(relationFiles, "relations");
    }
}

void GroupGenerator::filterCommonKmers(Buffer<Kmer> & qKmers,
                                       Buffer<std::pair<uint32_t, uint32_t>> & matchBuffer,
                                       const string & commonKmerDB) {
    string gtdbListDB;
    std::string diffIdxFileName = commonKmerDB +"/diffIdx";
    std::string diffIdxSplitFileName = commonKmerDB + "/split";

    size_t blankCnt = std::find_if(qKmers.buffer,
                                   qKmers.buffer + qKmers.startIndexOfReserve, 
                                   [](const auto& kmer) { return kmer.qInfo.sequenceID != 0;}
                                  ) - qKmers.buffer;

    size_t queryKmerNum = qKmers.startIndexOfReserve - blankCnt;
    std::cout << "Query k-mer number     : " << queryKmerNum << endl;

    // Filter out meaningless target splits
    MmapedData<DiffIdxSplit> diffIdxSplits = mmapData<DiffIdxSplit>(diffIdxSplitFileName.c_str(), 3);
    size_t numOfDiffIdxSplits = diffIdxSplits.fileSize / sizeof(DiffIdxSplit);
    size_t numOfDiffIdxSplits_use = numOfDiffIdxSplits;
    for (size_t i = 1; i < numOfDiffIdxSplits; i++) {
        if (diffIdxSplits.data[i].ADkmer == 0 || diffIdxSplits.data[i].ADkmer == UINT64_MAX) {
            numOfDiffIdxSplits_use--;
        }
    }

    // Divide query k-mer list into blocks for multi threading.
    std::vector<QueryKmerSplit> querySplits;
    size_t quotient = queryKmerNum / par.threads;
    size_t remainder = queryKmerNum % par.threads;
    size_t startIdx = blankCnt;
    size_t endIdx = 0; // endIdx is inclusive
    for (size_t i = 0; i < par.threads; i++) {
        endIdx = startIdx + quotient - 1;
        if (remainder > 0) {
            endIdx++;
            remainder--;
        }
        bool needLastTargetBlock = true;
        uint64_t queryAA = qKmers.buffer[startIdx].value;
        for (size_t j = 0; j < numOfDiffIdxSplits_use; j ++) {
            if (queryAA <= diffIdxSplits.data[j].ADkmer) {
                querySplits.emplace_back(startIdx, endIdx, diffIdxSplits.data[j - (j != 0)]);
                needLastTargetBlock = false;
                break;
            }
        }
        if (needLastTargetBlock) {
            querySplits.emplace_back(startIdx, endIdx, diffIdxSplits.data[numOfDiffIdxSplits_use - 2]);
        }
        startIdx = endIdx + 1;
    }
    munmap(diffIdxSplits.data, diffIdxSplits.fileSize);

    time_t beforeFilter = time(nullptr);
    std::cout << "Common k-mer searching : " << std::flush;
    #pragma omp parallel default(none), shared(matchBuffer, commonKmerDB, querySplits, qKmers, cout)
    {
        Buffer<std::pair<uint32_t, uint32_t>> localMatches(1024 * 1024 * 2);  // 16 Mb <queryID, pos>
        KmerDbReader * kmerDbReader
            = new KmerDbReader(commonKmerDB + "/diffIdx", 1024 * 1024, 1024 * 1024);
        std::vector<std::pair<uint32_t, uint32_t>> tempMatches;  
        bool hasOverflow = false;
    
        #pragma omp for schedule(dynamic, 1)
        for (size_t i = 0; i < querySplits.size(); i++) {
            kmerDbReader->setReadPosition(querySplits[i].diffIdxSplit);
            uint64_t tKmer = kmerDbReader->next();
            Kmer qKmer(UINT64_MAX, 0);
            for (size_t j = querySplits[i].start; j < querySplits[i].end + 1; j++) {
                // Reuse the AA matches if queries are identical
                if (qKmer.value == qKmers.buffer[j].value) {
                    if (unlikely(!localMatches.afford(tempMatches.size()))) {
                        if (!Buffer<std::pair<uint32_t, uint32_t>>::moveSmallToLarge(&localMatches, &matchBuffer)) {
                            hasOverflow = true;
                            break;
                        }
                    }
                    size_t posToWrite = localMatches.reserveMemory(tempMatches.size());
                    memcpy(localMatches.buffer + posToWrite, tempMatches.data(),
                           sizeof(std::pair<uint32_t, uint32_t>) * tempMatches.size());
                    for (size_t k = 0; k < tempMatches.size(); k++) {
                        localMatches.buffer[posToWrite + k].first = qKmers.buffer[j].qInfo.sequenceID;
                        localMatches.buffer[posToWrite + k].second = qKmers.buffer[j].qInfo.pos;
                    }
                    continue;
                }
                tempMatches.clear();
                // Get next query, and start to find
                qKmer = qKmers.buffer[j];

                // Skip target k-mers lexiocographically smaller
                while (!kmerDbReader->isCompleted() && qKmer.value > tKmer) {
                    tKmer = kmerDbReader->next();
                }

                // No match found - skip to the next query
                if (qKmer.value != tKmer) { continue; } 

                // Match found - load target k-mers matching at amino acid level
                while (!kmerDbReader->isCompleted() && qKmer.value == tKmer) {
                    tempMatches.emplace_back((uint32_t) qKmer.qInfo.sequenceID, (uint32_t) qKmer.qInfo.pos);
                    tKmer = kmerDbReader->next();                                      
                }

                if (unlikely(!localMatches.afford(tempMatches.size()))) {
                    if (!Buffer<std::pair<uint32_t, uint32_t>>::moveSmallToLarge(&localMatches, &matchBuffer)) {
                        hasOverflow = true;
                        break;
                    }
                }

                size_t posToWrite = localMatches.reserveMemory(tempMatches.size());
                memcpy(localMatches.buffer + posToWrite, tempMatches.data(),
                       sizeof(std::pair<uint32_t, uint32_t>) * tempMatches.size());
            } // End of one split

            // Move matches in the local buffer to the shared buffer
            if (!Buffer<std::pair<uint32_t, uint32_t>>::moveSmallToLarge(&localMatches, &matchBuffer)) {
                hasOverflow = true;
            }
        } // End of omp for (Iterating for splits)
    } // End of omp parallel
    std::cout << time(nullptr) - beforeFilter << " s (" << matchBuffer.startIndexOfReserve << " matches found)" << endl;

    time_t here = time(nullptr);
    SORT_PARALLEL(matchBuffer.buffer, matchBuffer.buffer + matchBuffer.startIndexOfReserve);    
    std::cout << "Sorting matches        : " << double(time(nullptr) - here) << " s" << std::endl;
 
    // Sort query k-mers by <seqID, pos>
    time_t firstSort = time(nullptr);
    SORT_PARALLEL(qKmers.buffer + blankCnt, qKmers.buffer + qKmers.startIndexOfReserve, Kmer::compareQKmerByIdAndPos);
    std::cout << "Query k-mer sorting (1): " << double(time(nullptr) - firstSort) << " s" << std::endl;

    // Filter neighbor k-mers
    here = time(nullptr);
    std::cout << "Filtering common k-mers: " << std::flush;
    size_t storePos = blankCnt;
    size_t lookingPos = blankCnt;
    size_t matchIdx = 0;
    while (lookingPos < qKmers.startIndexOfReserve) {
        if (matchIdx < matchBuffer.startIndexOfReserve) {
            // copy
            if (qKmers.buffer[lookingPos].qInfo.sequenceID < matchBuffer.buffer[matchIdx].first){
                qKmers.buffer[storePos++] = qKmers.buffer[lookingPos++];
            }
            // next target check
            else if(qKmers.buffer[lookingPos].qInfo.sequenceID > matchBuffer.buffer[matchIdx].first){
                matchIdx++;
            }
            // same seq
            else{
                // copy
                if (int64_t(qKmers.buffer[lookingPos].qInfo.pos) < int(matchBuffer.buffer[matchIdx].second) - COMMON_KMER_NEIGHBOR_SPAN){
                    qKmers.buffer[storePos++] = qKmers.buffer[lookingPos++];
                }
                // next target check
                else if(int(matchBuffer.buffer[matchIdx].second) + COMMON_KMER_NEIGHBOR_SPAN < int64_t(qKmers.buffer[lookingPos].qInfo.pos)){
                    matchIdx++;
                }
                // pass
                else{
                    lookingPos++;
                }
            }            
        }
        else{
            qKmers.buffer[storePos++] = qKmers.buffer[lookingPos++];
        }
    }
    qKmers.startIndexOfReserve = size_t(storePos);
    cout << double(time(nullptr) - here) << " s" << endl;
    
    // sort buffer by kmer
    time_t secondSort = time(nullptr);
    SORT_PARALLEL(qKmers.buffer, qKmers.buffer + qKmers.startIndexOfReserve, Kmer::compareQueryKmer);
    secondSort = time(nullptr) - secondSort;    
    cout << "Query k-mer sorting (2): " << double(secondSort) << " s" << endl;
    cout << "Filtered k-mer number  : " << storePos - blankCnt << endl;
    // Accumulate across splits; this point is outside the OpenMP region above, so no race.
    totalFilteredKmers += storePos - blankCnt;
}

void GroupGenerator::writeKmers(Buffer<Kmer>& queryKmerBuffer,
                                size_t processedReadCnt) {
    size_t blankCnt = std::find_if(queryKmerBuffer.buffer,
                                   queryKmerBuffer.buffer + queryKmerBuffer.startIndexOfReserve, 
                                   [](const auto& kmer) { return kmer.qInfo.sequenceID != 0;}
                                  ) - queryKmerBuffer.buffer;                                        
    size_t queryKmerNum = queryKmerBuffer.startIndexOfReserve - blankCnt;

    // Make k-mer boundaries based on actual distribution
    if (!boundariesInitialized) {
        size_t quotient = queryKmerNum / par.threads;
        size_t remainder = queryKmerNum % par.threads;
        size_t idx = blankCnt;
        for (size_t i = 0; i < par.threads - 1; i++) {
            idx = idx + quotient - (i == 0);
            if (remainder > 0) {
                idx++;
                remainder--;
            }
            if (idx >= queryKmerNum + blankCnt) {
                std::cout << "Warning: endIdx exceeded queryKmerNum, adjusting to max." << std::endl;
                idx = queryKmerNum + blankCnt - 1;
            }
            kmerBoundaries.emplace_back(queryKmerBuffer.buffer[idx].value);
        }
        boundariesInitialized = true;
    }

    // Make query k-mer ranges for each split
    std::vector<std::pair<size_t, size_t>> queryKmerRanges = getKmerRanges(queryKmerBuffer, blankCnt);
    #pragma omp parallel default(none), shared(queryKmerBuffer, queryKmerRanges, processedReadCnt)
    {
        size_t threadId = omp_get_thread_num();
        size_t startIdx = queryKmerRanges[threadId].first;
        size_t endIdx = queryKmerRanges[threadId].second;
        WriteBuffer<uint16_t> diffBuffer(this->outDir + "/kmer_delta_" + to_string(this->numOfSplits) + "_" + to_string(threadId), 1024 * 1024);
        WriteBuffer<uint32_t> infoBuffer(this->outDir + "/kmer_info_"  + to_string(this->numOfSplits) + "_" + to_string(threadId), 1024 * 1024);
        uint64_t lastKmer = 0;
        for (size_t i = startIdx; i < endIdx; i++) {
            queryKmerBuffer.buffer[i].qInfo.sequenceID += processedReadCnt;
            uint32_t id = static_cast<uint32_t>(queryKmerBuffer.buffer[i].qInfo.sequenceID);
            infoBuffer.write(&id);
            IndexCreator::getDiffIdx(lastKmer, queryKmerBuffer.buffer[i].value, diffBuffer);
        }
        infoBuffer.flush();
        diffBuffer.flush();
    }
    this->numOfSplits++;
}

std::vector<std::pair<size_t, size_t>> GroupGenerator::getKmerRanges(const Buffer<Kmer> & kmerBuffer, 
                                                                     size_t offset) {
    // Custom comparator to find the first k-mer that IS GREATER THAN the value (>)
    auto value_less_than_kmer = [](uint64_t val, const Kmer& kmer) {
        return val < kmer.value;
    };

    std::vector<std::pair<size_t, size_t>> ranges;
    size_t startIdx = offset;

    // For every boundary, find the first element > it. This marks the end of the current range.
    for (const uint64_t& boundary : kmerBoundaries) {
        // Find the first element that is > boundary.
        auto it_cut = std::upper_bound(kmerBuffer.buffer + startIdx,
                                       kmerBuffer.buffer + kmerBuffer.startIndexOfReserve,
                                       boundary, 
                                       value_less_than_kmer);

        // The distance gives the index relative to the beginning of the vector
        size_t cut_idx = std::distance(kmerBuffer.buffer, it_cut);

        ranges.push_back({startIdx, cut_idx});
        startIdx = cut_idx;
    }

    // Add the final range for all elements > the last boundary
    ranges.push_back({startIdx, kmerBuffer.startIndexOfReserve});
    return ranges;
}

void GroupGenerator::scanMHistogram(size_t processedReadCnt,
                                    std::vector<uint64_t> & mHistKmers,
                                    std::vector<uint64_t> & mHistPairs) {
    mHistKmers.assign(M_HIST_BUCKETS, 0);
    mHistPairs.assign(M_HIST_BUCKETS, 0);

    // Same partitioning as the emit pass: thread t owns k-mer value range t, so the per-thread
    // histograms are over disjoint k-mers and their sum is order-independent. Two runs with
    // different --threads therefore produce the same histogram, and the same cap.
    #pragma omp parallel num_threads(par.threads)
    {
        const int threadIdx = omp_get_thread_num();
        std::vector<uint64_t> localKmers(M_HIST_BUCKETS, 0);
        std::vector<uint64_t> localPairs(M_HIST_BUCKETS, 0);
        uint64_t localValues = 0;

        scanKmerRuns(threadIdx, processedReadCnt, localValues,
                     [&](uint64_t, const std::vector<uint32_t> & queryIds) {
            const size_t m = queryIds.size();
            if (m < 2) {
                return;
            }
            const size_t bucket = mHistBucket(m);
            localKmers[bucket]++;
            localPairs[bucket] += static_cast<uint64_t>(m) * static_cast<uint64_t>(m - 1) / 2;
        });

        #pragma omp critical
        {
            for (size_t b = 0; b < M_HIST_BUCKETS; ++b) {
                mHistKmers[b] += localKmers[b];
                mHistPairs[b] += localPairs[b];
            }
        }
    }
}

void GroupGenerator::makeSubGraph(size_t processedReadCnt) {
    cout << "Connecting reads with shared k-mer..." << endl;
    time_t beforeSearch = time(nullptr);

    const size_t RELATION_THRESHOLD = getRelationThreshold(par.threads, par.ramUsage);
    cout << "Flush threshold: " << RELATION_THRESHOLD << " pairs/thread (memory budget "
         << humanBytes(getMemoryBudgetBytes(par.ramUsage)) << ", --max-ram "
         << par.ramUsage << " GiB)" << endl;

    // Disk counterpart of the memory budget above. A safety ceiling, not a target: a run that
    // fits under it behaves exactly as it does without one.
    const size_t tmpDiskBudget = getTmpDiskBudgetBytes(outDir, par.maxTmpDiskMiB);
    if (par.maxTmpDiskMiB > 0) {
        cout << "Tmp disk budget: " << humanBytes(tmpDiskBudget)
             << " (--max-tmp-disk " << par.maxTmpDiskMiB << " MiB)" << endl;
    } else if (tmpDiskBudget == TMP_DISK_UNLIMITED) {
        cout << "[WARN] could not measure free space at " << outDir
             << "; running without a tmp disk cap. Pass --max-tmp-disk to set one." << endl;
    } else {
        cout << "Tmp disk budget: " << humanBytes(tmpDiskBudget) << " (auto: 80% of "
             << humanBytes(getFreeDiskBytes(outDir)) << " free at " << outDir << ")" << endl;
    }

    // Second fold trigger, independent of the budget: the merge can only open so many streams
    // per unit, so a prefix longer than the fan-in has to be folded regardless of disk. Note
    // that getMergeFanIn looks at the fd limit and the thread count only -- it knows nothing
    // about how much data is on disk, which is why the budget above is the primary trigger.
    const size_t foldFanIn = getMergeFanIn(par.threads, static_cast<size_t>(par.threads));
    std::atomic<int> counter(0);

    // High-frequency k-mers are dropped by --max-kmer-freq-ratio. Record which ones so the
    // information loss stays auditable instead of silent.
    size_t skippedKmerCnt = 0;
    size_t skippedMaxM = 0;
    size_t skippedSumM = 0;
    double skippedPairEst = 0.0; // sum of C(m,2) over skipped k-mers; double because it overflows uint64

    // Edge-volume instrumentation: how much the C(m,2) clique actually produced, and the
    // largest m that survived the skip -- that m sets the per-k-mer worst case.
    size_t emittedEdgeCnt = 0; // Relation records written to subGraph_*
    size_t maxKeptM = 0;       // largest reads-per-k-mer that was NOT skipped
    // Sum of C(m,2) over kept k-mers = pair OCCURRENCES the clique produced. Records written
    // are fewer, because pair2weight collapses a pair seen across several k-mers into one
    // entry. emittedEdgeCnt / keptPairSum is that collapse factor, measured on this run.
    uint64_t keptPairSum = 0;

    // m distribution over ALL k-mers, skipped or not. The [edges] summary only reports what
    // the chosen thresholds produced; this reports what any other threshold would have.
    std::vector<uint64_t> mHistKmers(M_HIST_BUCKETS, 0);
    std::vector<uint64_t> mHistPairs(M_HIST_BUCKETS, 0);

    // Skip threshold, resolved once. An explicit --max-kmer-reads wins and skips the pre-pass
    // entirely; otherwise --max-kmer-quantile picks the cap from the reads-per-k-mer
    // distribution, which has to be measured before any edge is emitted. The cap stays an
    // absolute m, only its choice becomes a property of the data: m does not scale with the
    // size of the dataset but with per-genome coverage, so neither a hand-picked constant nor a
    // fraction of the read count transfers (--max-kmer-freq-ratio failed on exactly that).
    size_t absThr = (par.maxKmerReads > 0) ? static_cast<size_t>(par.maxKmerReads) : 0;
    bool capIsAuto = false;
    if (par.maxKmerReads <= 0 && par.maxKmerQuantile > 0.0f && par.maxKmerQuantile < 1.0f) {
        cout << "[skip] pre-pass: measuring the reads-per-k-mer distribution..." << endl;
        const time_t preStart = time(nullptr);
        std::vector<uint64_t> preHistKmers;
        std::vector<uint64_t> preHistPairs;
        scanMHistogram(processedReadCnt, preHistKmers, preHistPairs);
        // The pre-pass counted every k-mer value once already; the emit pass counts them again.
        // Leaving both in would double the denominator the disk projection divides by.
        kmerValuesTotal.store(0, std::memory_order_relaxed);

        absThr = capFromQuantile(preHistKmers, par.maxKmerQuantile);
        capIsAuto = (absThr > 0);

        uint64_t totalKmers = 0, totalPairs = 0, keptKmers = 0, keptPairs = 0;
        for (size_t b = 0; b < M_HIST_BUCKETS; ++b) {
            totalKmers += preHistKmers[b];
            totalPairs += preHistPairs[b];
            // Bucket b covers m in [2^b, 2^(b+1) - 1]; it survives when its top is under the cap.
            const bool kept = (absThr == 0)
                || (b + 1 < 63 && ((static_cast<size_t>(1) << (b + 1)) - 1) <= absThr);
            if (kept) {
                keptKmers += preHistKmers[b];
                keptPairs += preHistPairs[b];
            }
        }
        cout << "[skip] pre-pass done: " << double(time(nullptr) - preStart) << " s, "
             << totalKmers << " k-mers with m >= 2" << endl;
        cout << "[skip] auto cap " << (absThr > 0 ? to_string(absThr) : string("none"))
             << " (quantile " << par.maxKmerQuantile << "): keeps "
             << (totalKmers ? 100.0 * static_cast<double>(keptKmers)
                                    / static_cast<double>(totalKmers) : 0.0)
             << "% of k-mers, "
             << (totalPairs ? 100.0 * static_cast<double>(keptPairs)
                                    / static_cast<double>(totalPairs) : 0.0)
             << "% of pair occurrences" << endl;
    }
    cout << "[skip] threshold: --max-kmer-reads "
         << (absThr > 0 ? to_string(absThr) : string("off"))
         << (capIsAuto ? " (auto)" : "") << endl;

    #pragma omp parallel num_threads(par.threads)
    {
        int threadIdx = omp_get_thread_num();
        std::unordered_map<uint64_t, uint16_t> pair2weight;
        pair2weight.reserve(RELATION_THRESHOLD);

        const string skippedPartName = outDir + "/skipped_kmers_" + to_string(threadIdx);
        ofstream skippedPart(skippedPartName);
        if (!skippedPart.is_open()) {
            cerr << "Error opening file: " << skippedPartName << endl;
        }
        size_t localSkippedCnt = 0;
        size_t localMaxM = 0;
        size_t localSumM = 0;
        double localPairEst = 0.0;
        size_t localEmittedEdges = 0;
        size_t localMaxKeptM = 0;
        std::vector<uint64_t> localHistKmers(M_HIST_BUCKETS, 0);
        std::vector<uint64_t> localHistPairs(M_HIST_BUCKETS, 0);
        uint64_t localKeptPairs = 0;

        uint64_t localValuesDone = 0;

        scanKmerRuns(threadIdx, processedReadCnt, localValuesDone,
                     [&](uint64_t minKmer, const std::vector<uint32_t> & currentQueryIds) {
            const size_t m = currentQueryIds.size();

            // Histogram every k-mer, skipped or not, so the summary table can answer
            // "what would a different cap have cost".
            if (m >= 2) {
                const size_t bucket = mHistBucket(m);
                localHistKmers[bucket]++;
                localHistPairs[bucket] += static_cast<uint64_t>(m) * static_cast<uint64_t>(m - 1) / 2;
            }

            if (absThr > 0 && m > absThr) {
                skippedPart << minKmer << "\t" << m << "\n";
                localSkippedCnt++;
                localSumM += m;
                if (m > localMaxM) { localMaxM = m; }
                localPairEst += static_cast<double>(m) * static_cast<double>(m - 1) * 0.5;
                return; // `continue` before the scan became a callback
            }
            // Every read pair sharing this k-mer gets an edge: the full C(m,2) clique.
            // Edge weight therefore counts shared k-mers directly, which is what the
            // Phase 1/2 thresholds and the knee detection assume. currentQueryIds is
            // sorted and deduplicated above, so id[i] < id[j] holds and pairKey needs no
            // min/max normalization. The quadratic term is bounded by the high-frequency
            // k-mer skip above, not by reducing the clique.
            if (m > localMaxKeptM) { localMaxKeptM = m; }
            if (m >= 2) {
                localKeptPairs += static_cast<uint64_t>(m) * static_cast<uint64_t>(m - 1) / 2;
            }
            for (size_t i = 0; i + 1 < m; ++i) {
                const uint64_t idHi = static_cast<uint64_t>(currentQueryIds[i]) << 32;
                for (size_t j = i + 1; j < m; ++j) {
                    const uint64_t pairKey = idHi | static_cast<uint64_t>(currentQueryIds[j]);
                    addSat16(pair2weight[pairKey], 1);
                }
            }
            if (pair2weight.size() >= RELATION_THRESHOLD) {
                size_t counter_now = counter.fetch_add(1, memory_order_relaxed);
                localEmittedEdges += pair2weight.size();
                saveSubGraphToFile(pair2weight, counter_now, processedReadCnt);
                pair2weight.clear();
                // Publish progress before the guards read it, so the projection matches the
                // bytes that have just been written.
                kmerValuesDone.fetch_add(localValuesDone, std::memory_order_relaxed);
                localValuesDone = 0;
                checkDiskHeadroom();
                // Fold what is already complete instead of letting the whole emit pile up.
                // No barrier: this thread folds only the finished prefix while the others keep
                // writing indices above it, so the two never touch the same file.
                maybeFoldEmitted(tmpDiskBudget, foldFanIn);
            }
        });

        if (!pair2weight.empty()) {
            size_t counter_now = counter.fetch_add(1, std::memory_order_relaxed);
            localEmittedEdges += pair2weight.size();
            saveSubGraphToFile(pair2weight, counter_now, processedReadCnt);
        } else {
            cout << "Thread " << threadIdx << " has no relations to write." << endl;
        }
        skippedPart.close();

        #pragma omp critical
        {
            skippedKmerCnt  += localSkippedCnt;
            skippedSumM     += localSumM;
            skippedPairEst  += localPairEst;
            emittedEdgeCnt  += localEmittedEdges;
            if (localMaxM > skippedMaxM) { skippedMaxM = localMaxM; }
            if (localMaxKeptM > maxKeptM) { maxKeptM = localMaxKeptM; }
            keptPairSum += localKeptPairs;
            for (size_t b = 0; b < M_HIST_BUCKETS; ++b) {
                mHistKmers[b] += localHistKmers[b];
                mHistPairs[b] += localHistPairs[b];
            }
        }
    }
    this->numOfGraph = counter.load(std::memory_order_relaxed);

    // Concatenate the per-thread skip logs. Thread t owns k-mer value range t (writeKmers
    // fixes kmerBoundaries once, so the ranges are disjoint and ordered), which makes thread
    // order equal to ascending k-mer order -- no merge sort needed. Unlike the other
    // intermediates this file is diagnostic output and is kept, not deleted.
    {
        const string skippedFileName = outDir + "/skipped_kmers";
        ofstream skippedAll(skippedFileName);
        if (!skippedAll.is_open()) {
            cerr << "Error opening file: " << skippedFileName << endl;
        } else {
            std::vector<std::string> skippedParts;
            for (int t = 0; t < par.threads; t++) {
                skippedParts.push_back(outDir + "/skipped_kmers_" + to_string(t));
            }
            for (const std::string & part : skippedParts) {
                ifstream in(part);
                // Streaming an empty rdbuf sets failbit on the destination, so guard it.
                if (in.is_open() && in.peek() != std::ifstream::traits_type::eof()) {
                    skippedAll << in.rdbuf();
                }
            }
            skippedAll.close();
            reportAndRemoveFiles(skippedParts, "skipped_kmers parts");
        }
    }

    if (skippedKmerCnt > 0) {
        cout << "[skip] " << skippedKmerCnt << " k-mers skipped (max m: " << skippedMaxM
             << ", sum m: " << skippedSumM
             << ", est. dropped read-pairs: " << skippedPairEst << ")" << endl;
        cout << "[skip] see " << outDir << "/skipped_kmers" << endl;
    } else if (absThr == 0) {
        cout << "[skip] 0 k-mers skipped (--max-kmer-reads off)" << endl;
    } else {
        cout << "[skip] 0 k-mers skipped (no k-mer exceeded the threshold)" << endl;
    }

    // k-mer files are fully consumed by the merge above; report + remove to free disk.
    {
        std::vector<std::string> deltaFiles, infoFiles;
        for (size_t i = 0; i < this->numOfSplits; i++) {
            for (int t = 0; t < par.threads; t++) {
                deltaFiles.push_back(outDir + "/kmer_delta_" + to_string(i) + "_" + to_string(t));
                infoFiles.push_back(outDir + "/kmer_info_"  + to_string(i) + "_" + to_string(t));
            }
        }
        reportAndRemoveFiles(deltaFiles, "kmer_delta");
        reportAndRemoveFiles(infoFiles, "kmer_info");
    }

    // m-distribution table. Each row is a candidate --max-kmer-reads value (the bucket's
    // upper bound) with the edge volume that cap would produce: the cumulative sum of
    // C(m,2) over every k-mer with m at or below it. Pick the row that fits the disk budget.
    // Counts cover all k-mers with m >= 2, whether or not the current run skipped them, so
    // the table stays valid for thresholds other than the one just used.
    {
        size_t topBucket = 0;
        uint64_t totalKmers = 0;
        uint64_t totalPairs = 0;
        for (size_t b = 0; b < M_HIST_BUCKETS; ++b) {
            totalKmers += mHistKmers[b];
            totalPairs += mHistPairs[b];
            if (mHistKmers[b] > 0) { topBucket = b; }
        }

        // Pair occurrences over-count what actually lands on disk, because pair2weight
        // collapses a pair seen across several k-mers into one record. Measure that factor
        // on this run and apply it, so the table predicts records rather than an upper bound.
        const double dedup = (keptPairSum > 0)
            ? static_cast<double>(emittedEdgeCnt) / static_cast<double>(keptPairSum)
            : 1.0;

        cout << "[mhist] reads-per-k-mer distribution (m >= 2): " << totalKmers
             << " k-mers, " << totalPairs << " pair occurrences uncapped" << endl;
        cout << "[mhist] measured collapse factor: " << emittedEdgeCnt << " records / "
             << keptPairSum << " kept pair occurrences = " << dedup << endl;
        cout << "[mhist] estimates below extrapolate that factor to other caps; it was"
             << " measured only on k-mers this run kept, so expect roughly +/-20%" << endl;
        if (totalKmers > 0) {
            cout << "[mhist] " << setw(16) << "cap (max m)"
                 << setw(14) << "k-mers"
                 << setw(18) << "pairs in range"
                 << setw(18) << "cum pairs"
                 << setw(16) << "est. records"
                 << setw(12) << "est. disk" << endl;
            uint64_t cumPairs = 0;
            for (size_t b = 0; b <= topBucket; ++b) {
                cumPairs += mHistPairs[b];
                if (mHistKmers[b] == 0) { continue; }
                const size_t lo = (size_t(1) << b);
                const size_t hi = (size_t(1) << (b + 1)) - 1;
                const uint64_t estRecords = static_cast<uint64_t>(static_cast<double>(cumPairs) * dedup);
                cout << "[mhist] " << setw(16) << (to_string(hi) + " (" + to_string(lo) + "-" + to_string(hi) + ")")
                     << setw(14) << mHistKmers[b]
                     << setw(18) << mHistPairs[b]
                     << setw(18) << cumPairs
                     << setw(16) << estRecords
                     << setw(12) << humanBytes(estRecords * sizeof(Relation)) << endl;
            }
            cout << "[mhist] read off the row whose est. disk fits the budget,"
                 << " then pass its cap as --max-kmer-reads" << endl;
        }
    }

    // Edge-volume summary. maxKeptM is the worst per-k-mer case that survived the skip:
    // that single k-mer contributed C(maxKeptM, 2) edges, so it is the knob to watch when
    // tuning --max-kmer-freq-ratio against a disk budget.
    const uint64_t rawEdgeBytes = static_cast<uint64_t>(emittedEdgeCnt) * sizeof(Relation);
    const uint64_t diskEdgeBytes = this->subGraphBytesOnDisk.load(std::memory_order_relaxed);
    cout << "[edges] emitted " << emittedEdgeCnt << " edge records into "
         << this->numOfGraph << " subgraphs ("
         << humanBytes(rawEdgeBytes) << " before merge, "
         << humanBytes(diskEdgeBytes) << " on disk, zstd level " << RELATION_ZSTD_LEVEL;
    if (diskEdgeBytes > 0) {
        cout << ", " << fixed << setprecision(2)
             << (static_cast<double>(rawEdgeBytes) / static_cast<double>(diskEdgeBytes)) << "x";
    }
    cout << ")" << endl;
    // What the run actually held at once, which is the number a shared filesystem cares about.
    // "0 folds" means the data fit under the budget and nothing was folded early -- the path
    // this build takes when there is room, identical to the one before incremental folding.
    cout << "[disk] peak subGraph footprint: "
         << humanBytes(this->subGraphPeakBytes.load(std::memory_order_relaxed)) << " (budget "
         << (getTmpDiskBudgetBytes(outDir, par.maxTmpDiskMiB) == TMP_DISK_UNLIMITED
             ? string("unlimited")
             : humanBytes(getTmpDiskBudgetBytes(outDir, par.maxTmpDiskMiB)))
         << ", " << this->foldRounds << " folds)" << endl;
    cout << "[edges] largest kept k-mer: m = " << maxKeptM << " -> "
         << (static_cast<double>(maxKeptM) * static_cast<double>(maxKeptM > 0 ? maxKeptM - 1 : 0) * 0.5)
         << " edges from that k-mer alone" << endl;

    cout << "Relations generated from files successfully." << endl;
    cout << "Time spent: " << double(time(nullptr) - beforeSearch) << " seconds." << endl;
}

void GroupGenerator::saveSubGraphToFile(const unordered_map<uint64_t, uint16_t>& pair2weight,
                                        const size_t counter_now,
                                        size_t processedReadCnt) {
    // Split by relations_* route at write time. Because the route is a function of
    // (id1, id2) alone, a pair can never land under two routes -- which is what makes the
    // per-route merges independent and lets them run in parallel.
    const size_t routeCnt = static_cast<size_t>(par.threads) * 2 + 1;
    const size_t rangeSize = getRouteRangeSize(processedReadCnt, par.threads);

    // byUnit[route][shard]. Every route carries a shard index so the merge's filename rule is
    // uniform; only the cross bucket actually has more than one.
    std::vector<std::vector<std::vector<Relation>>> byUnit(routeCnt);
    for (size_t route = 0; route < routeCnt; ++route) {
        byUnit[route].resize(shardsForRoute(route, par.threads));
    }
    for (const auto& [pairKey, weight] : pair2weight) {
        const uint32_t id1 = static_cast<uint32_t>(pairKey >> 32);
        const uint32_t id2 = static_cast<uint32_t>(pairKey & 0xFFFFFFFF);
        const size_t route = routeOf(id1, id2, par.threads, rangeSize);
        byUnit[route][shardOf(id1, byUnit[route].size())].push_back(makeRelation(id1, id2, weight));
    }

    // Empty units still get an empty file: the merge's input list is then simply
    // subGraph_{0..numOfGraph-1}_{route}_{shard} with no existence checks.
    uint64_t compressedBytes = 0;
    for (size_t route = 0; route < routeCnt; ++route) {
        for (size_t shard = 0; shard < byUnit[route].size(); ++shard) {
            const string subGraphFileName = subGraphName(counter_now, route, shard);
            std::vector<Relation> & relations = byUnit[route][shard];
            sort(relations.begin(), relations.end(), Relation::compare);
            // A failed write used to be reported and skipped. It is fatal now: the merge would
            // otherwise read a file that is missing or short and produce a silently wrong graph.
            RelationWriter out(subGraphFileName);
            for (size_t i = 0; i < relations.size(); ++i) {
                out.write(relations[i]);
            }
            out.finish();
            compressedBytes += out.compressedBytes();
        }
    }
    subGraphBytesOnDisk.fetch_add(compressedBytes, std::memory_order_relaxed);
    noteSubGraphBytes(static_cast<int64_t>(compressedBytes));
    // Only now, with every unit file of this flush closed, may the folder see this index.
    markEmitComplete(counter_now);
}

void GroupGenerator::noteSubGraphBytes(int64_t delta) {
    const int64_t live = subGraphLiveBytes.fetch_add(delta, std::memory_order_relaxed) + delta;
    if (live <= 0) {
        return;
    }
    const uint64_t liveBytes = static_cast<uint64_t>(live);
    uint64_t peak = subGraphPeakBytes.load(std::memory_order_relaxed);
    while (liveBytes > peak
           && !subGraphPeakBytes.compare_exchange_weak(peak, liveBytes,
                                                       std::memory_order_relaxed)) {
        // compare_exchange_weak refreshed `peak`; retry until it stops growing.
    }
}

void GroupGenerator::checkDiskHeadroom() {
    const size_t freeBytes = getFreeDiskBytes(outDir);
    if (freeBytes == 0) {
        return; // cannot measure; getTmpDiskBudgetBytes already warned about that
    }

    // Hard guard first. It needs no estimate, so it holds even when the projection is wrong.
    const size_t reserve = getDiskReserveBytes(outDir);
    if (freeBytes <= reserve) {
        cerr << "[ERROR] out of disk space at " << outDir << ": "
             << humanBytes(freeBytes) << " free, keeping " << humanBytes(reserve)
             << " in reserve.\n"
             << "        Intermediates hold " << humanBytes(subGraphBytesOnDisk.load())
             << " so far. Point the output directory at a larger filesystem, or lower\n"
             << "        --max-kmer-reads (currently " << par.maxKmerReads
             << ") to cut the edge volume -- see the [mhist] table for what each cap costs."
             << endl;
        exit(EXIT_FAILURE);
    }

    // Projection. Warn once, do not stop: folding decouples the footprint from the total.
    if (volumeWarned.load(std::memory_order_relaxed) != 0) {
        return;
    }
    const uint64_t total = kmerValuesTotal.load(std::memory_order_relaxed);
    const uint64_t done = kmerValuesDone.load(std::memory_order_relaxed);
    if (total == 0 || done == 0) {
        return;
    }
    const double progress = static_cast<double>(done) / static_cast<double>(total);
    if (progress < 0.02 || progress >= 1.0) {
        return; // too early to extrapolate, or already finished
    }
    const double projected = static_cast<double>(subGraphBytesOnDisk.load()) / progress;
    if (projected <= static_cast<double>(freeBytes)) {
        return;
    }
    int expected = 0;
    if (!volumeWarned.compare_exchange_strong(expected, 1)) {
        return; // another thread is printing it
    }
    cout << "[WARN] at " << fixed << setprecision(1) << (progress * 100.0)
         << "% of the emit, the intermediates project to "
         << humanBytes(static_cast<uint64_t>(projected)) << " in total but only "
         << humanBytes(freeBytes) << " is free at " << outDir << ".\n"
         << "       Folding keeps what is held at once well below that total, so the run may\n"
         << "       still finish. If it does not, lower --max-kmer-reads (currently "
         << par.maxKmerReads << ") -- the [mhist] table prints what each cap would cost --\n"
         << "       or use a larger filesystem. This early in the run the projection is rough\n"
         << "       and tends to overshoot; it is a heads-up, not a verdict." << endl;
}

void GroupGenerator::markEmitComplete(size_t flushIdx) {
    std::lock_guard<std::mutex> guard(foldMutex);
    if (emitDone.size() <= flushIdx) {
        emitDone.resize(flushIdx + 1, 0);
    }
    emitDone[flushIdx] = 1;
    while (emitWatermark < emitDone.size() && emitDone[emitWatermark] != 0) {
        ++emitWatermark;
    }
}

void GroupGenerator::maybeFoldEmitted(size_t tmpDiskBudget, size_t maxFanIn) {
    const int64_t live = subGraphLiveBytes.load(std::memory_order_relaxed);
    const bool overBudget = tmpDiskBudget != TMP_DISK_UNLIMITED
                            && live > 0
                            && static_cast<size_t>(live) >= tmpDiskBudget;

    // Cheap pre-check outside the lock. The authoritative test is repeated under it.
    if (!overBudget) {
        std::lock_guard<std::mutex> guard(foldMutex);
        if (emitWatermark - foldedUpTo < maxFanIn) {
            return;
        }
    }

    std::unique_lock<std::mutex> guard(foldMutex, std::try_to_lock);
    if (!guard.owns_lock()) {
        return; // another thread is folding; keep emitting instead of waiting on it
    }

    const size_t from = foldedUpTo;
    const size_t to = emitWatermark;
    // Nothing to gain from folding a single flush: it would rewrite the same records.
    if (to - from < 2) {
        return;
    }

    const size_t round = foldRounds++;
    foldedUpTo = to;
    if (foldedOutputs.empty()) {
        foldedOutputs.resize(unitCount());
    }

    const size_t routeCnt = static_cast<size_t>(par.threads) * 2 + 1;
    // One extra stream per unit for the previous fold's output, which is re-folded below.
    const size_t bufElems = getMergeBufferElems(to - from + 1, par.ramUsage, 1);
    size_t foldedFiles = 0;
    for (size_t route = 0; route < routeCnt; ++route) {
        const size_t shardCnt = shardsForRoute(route, par.threads);
        for (size_t shard = 0; shard < shardCnt; ++shard) {
            // Earlier folds' outputs go back in. Folding only the new flushes would let those
            // outputs pile up untouched, and the live footprint would climb past the budget one
            // fold at a time. Re-folding keeps exactly one file per unit, so the footprint
            // tracks the merged-so-far size -- the floor no amount of folding can go below.
            std::vector<std::string> & folded = foldedOutputs[unitIndexOf(route, shard)];
            std::vector<std::string> inputs = folded;
            inputs.reserve(inputs.size() + (to - from));
            for (size_t i = from; i < to; ++i) {
                inputs.push_back(subGraphName(i, route, shard));
            }
            // A distinct prefix from reduceSubGraphFanIn's subGraph_p*: that one numbers its
            // rounds independently during the final merge, and a shared name would let the two
            // overwrite each other.
            const std::string output = outDir + "/subGraph_f" + to_string(round)
                                     + "_" + to_string(route) + "_" + to_string(shard);
            mergeSubGraphBatch(inputs, output, bufElems);

            for (size_t i = 0; i < inputs.size(); ++i) {
                noteSubGraphBytes(-static_cast<int64_t>(FileUtil::getFileSize(inputs[i])));
                std::remove(inputs[i].c_str());
            }
            noteSubGraphBytes(static_cast<int64_t>(FileUtil::getFileSize(output)));
            folded.assign(1, output);
            ++foldedFiles;
        }
    }

    cout << "[disk] fold " << round << ": flushes [" << from << ", " << to << ") -> "
         << foldedFiles << " files, live "
         << humanBytes(static_cast<uint64_t>(std::max<int64_t>(0, subGraphLiveBytes.load())))
         << endl;
}

// No lock: this runs after makeSubGraph's parallel region has joined, so folding is over and
// foldedUpTo/foldedOutputs are stable.
std::vector<std::string> GroupGenerator::unitInputFiles(size_t route, size_t shard) const {
    std::vector<std::string> files;
    const size_t unit = unitIndexOf(route, shard);
    if (unit < foldedOutputs.size()) {
        files = foldedOutputs[unit];
    }
    files.reserve(files.size() + (this->numOfGraph - foldedUpTo));
    for (size_t i = foldedUpTo; i < this->numOfGraph; ++i) {
        files.push_back(subGraphName(i, route, shard));
    }
    return files;
}

size_t GroupGenerator::mergeSubGraphBatch(const std::vector<std::string> & inputs,
                                          const std::string & output,
                                          size_t bufElems) {
    std::vector<std::unique_ptr<RelationReader>> readers(inputs.size());
    for (size_t i = 0; i < inputs.size(); ++i) {
        readers[i].reset(new RelationReader(inputs[i], bufElems));
    }

    RelationWriter out(output);

    MergeHeap heap;
    for (size_t i = 0; i < readers.size(); ++i) {
        const Relation first = readers[i]->getNext();
        if (!(first == Relation())) {
            MergeHeapEntry entry;
            entry.rel = first;
            entry.stream = i;
            heap.push(entry);
        }
    }

    size_t written = 0;
    while (!heap.empty()) {
        const Relation minRelation = heap.top().rel;
        uint16_t totalWeight = 0;
        while (!heap.empty() && heap.top().rel == minRelation) {
            const size_t stream = heap.top().stream;
            addSat16(totalWeight, heap.top().rel.weight);
            heap.pop();
            const Relation next = readers[stream]->getNext();
            if (!(next == Relation())) {
                MergeHeapEntry entry;
                entry.rel = next;
                entry.stream = stream;
                heap.push(entry);
            }
        }
        out.write(makeRelation(minRelation.id1, minRelation.id2, totalWeight));
        ++written;
    }
    out.finish();

    return written;
}

void GroupGenerator::reduceSubGraphFanIn(std::vector<std::string> & files, size_t maxFanIn,
                                         size_t route, size_t shard) {
    if (files.size() <= maxFanIn) {
        return;
    }

    const size_t bufElems = getMergeBufferElems(maxFanIn, par.ramUsage, par.threads);
    size_t round = 0;
    while (files.size() > maxFanIn) {
        std::vector<std::string> next;
        size_t batchIdx = 0;
        for (size_t start = 0; start < files.size(); start += maxFanIn) {
            const size_t end = std::min(start + maxFanIn, files.size());
            const std::vector<std::string> batch(files.begin() + start, files.begin() + end);
            // Route AND shard are part of the name: units fold at the same time, and a shared
            // name would have two of them write the same file. Route alone is not enough --
            // the cross bucket's shards all carry the same route.
            const std::string out = outDir + "/subGraph_p" + to_string(route)
                                  + "_s" + to_string(shard)
                                  + "_r" + to_string(round) + "_" + to_string(batchIdx);
            mergeSubGraphBatch(batch, out, bufElems);
            // Delete this batch now, not after the whole round: peak disk then grows by one
            // batch's output instead of by a second copy of every subgraph.
            for (size_t i = 0; i < batch.size(); ++i) {
                std::remove(batch[i].c_str());
            }
            next.push_back(out);
            ++batchIdx;
        }
        files.swap(next);
        ++round;
    }
}

void GroupGenerator::mergeUnit(size_t route, size_t shard, const std::string & outPath,
                               const std::vector<std::string> & inputFiles,
                               size_t bufElems, size_t maxFanIn,
                               size_t & mergedOut, size_t & ceilingOut) {
    mergedOut = 0;
    ceilingOut = 0;

    std::vector<std::string> files = inputFiles;
    reduceSubGraphFanIn(files, maxFanIn, route, shard);

    const size_t streamCnt = files.size();
    std::vector<std::unique_ptr<RelationReader>> readers(streamCnt);
    for (size_t i = 0; i < streamCnt; ++i) {
        readers[i].reset(new RelationReader(files[i], bufElems));
    }

    // Output is relations_*, which makeGroups/mergeBySupport/makeGroupsPhase2 read with
    // ReadBuffer<Relation>. That format stays raw -- only the subGraph_* inputs above change.
    WriteBuffer<Relation> out(outPath, 1024 * 1024);

    // No routing decision here: every edge in these files already belongs to this unit.
    MergeHeap heap;
    for (size_t i = 0; i < streamCnt; ++i) {
        const Relation first = readers[i]->getNext();
        if (!(first == Relation())) {
            MergeHeapEntry entry;
            entry.rel = first;
            entry.stream = i;
            heap.push(entry);
        }
    }

    while (!heap.empty()) {
        const Relation minRelation = heap.top().rel;
        uint16_t totalWeight = 0;
        while (!heap.empty() && heap.top().rel == minRelation) {
            const size_t stream = heap.top().stream;
            addSat16(totalWeight, heap.top().rel.weight);
            heap.pop();
            const Relation next = readers[stream]->getNext();
            if (!(next == Relation())) {
                MergeHeapEntry entry;
                entry.rel = next;
                entry.stream = stream;
                heap.push(entry);
            }
        }
        if (totalWeight == UINT16_MAX) ceilingOut++;
        mergedOut++;

        // Zeroed in place, not via a helper's return value -- see mergeGraph_one for why.
        Relation rel;
        std::memset(static_cast<void *>(&rel), 0, sizeof(rel));
        rel.id1 = minRelation.id1;
        rel.id2 = minRelation.id2;
        rel.weight = totalWeight;
        out.write(&rel);
    }
    out.flush();

    readers.clear(); // close every input before unlinking it
    for (size_t i = 0; i < streamCnt; ++i) {
        std::remove(files[i].c_str());
    }
}

void GroupGenerator::mergeGraph() {
    cout << "Merging subgraphs" << endl;
    time_t before = time(nullptr);

    // Units merge independently: saveSubGraphToFile already split every flush by
    // (relations_* route, shard), and both are functions of the ids alone, so no pair spans
    // two units. Sharding exists because the routing rule funnels ~88% of all edges into the
    // cross bucket, which would otherwise cap the whole parallel merge at ~1.13x.
    const size_t routeCnt = static_cast<size_t>(par.threads) * 2 + 1;
    std::vector<std::pair<size_t, size_t>> units; // (route, shard)
    for (size_t r = 0; r < routeCnt; ++r) {
        const size_t shardCnt = shardsForRoute(r, par.threads);
        for (size_t s = 0; s < shardCnt; ++s) {
            units.emplace_back(r, s);
        }
    }
    const size_t unitCnt = units.size();
    const size_t concurrentMergers = std::min(static_cast<size_t>(par.threads), unitCnt);
    const size_t maxFanIn = getMergeFanIn(par.threads, concurrentMergers);
    const size_t mergeBufElems = getMergeBufferElems(std::min(this->unitStreamCount(), maxFanIn),
                                                     par.ramUsage, concurrentMergers);
    const size_t streamsPerUnit = std::min(this->unitStreamCount(), maxFanIn);

    // A shard writes to a temp; the cross bucket's shards are concatenated afterwards.
    // Concatenation is legal because the relations_* consumers scan the file end to end and
    // never depend on ordering.
    auto unitOutPath = [&](size_t route, size_t shard) {
        return (shardsForRoute(route, par.threads) == 1)
             ? (outDir + "/relations_" + to_string(route) + ".bin")
             : (outDir + "/relations_" + to_string(route) + "_s" + to_string(shard) + ".tmp");
    };

    const size_t peakFds = concurrentMergers * mergeFdsPerUnit(std::min(this->unitStreamCount(), maxFanIn));
    cout << "Merge: " << unitCnt << " units (" << routeCnt << " routes, cross bucket sharded x"
         << shardsForRoute(static_cast<size_t>(par.threads), par.threads) << ") x "
         << this->numOfGraph << " subgraphs, " << concurrentMergers
         << " concurrent (fan-in " << maxFanIn << ", peak fds " << peakFds
         << " of soft limit " << getOpenFileLimit() << ")" << endl;
    if (peakFds >= getOpenFileLimit()) {
        cout << "Warning: peak descriptor use " << peakFds << " meets the soft limit "
             << getOpenFileLimit() << "; raise it or lower --threads." << endl;
    }
    cout << "Merge read buffers: " << streamsPerUnit << " streams x " << mergeBufElems
         << " Relations per unit, "
         << humanBytes(static_cast<uint64_t>(streamsPerUnit) * mergeBufElems
                       * sizeof(Relation) * concurrentMergers)
         << " total" << endl;
    if (this->numOfGraph > maxFanIn) {
        cout << "Folding each unit's " << this->numOfGraph << " subgraphs down to <= "
             << maxFanIn << endl;
    }
    if (mergeBufElems == MERGE_BUFFER_MIN_ELEMS) {
        cout << "Warning: read buffers are at the " << MERGE_BUFFER_MIN_ELEMS
             << "-element floor. Merging will be I/O bound." << endl;
    }

    // Each unit deletes its own inputs as it finishes, so total subgraph bytes are measured
    // up front to keep the [disk] report the single-pass merge used to print.
    uint64_t subGraphBytes = 0;
    for (size_t i = 0; i < this->numOfGraph; ++i) {
        for (size_t u = 0; u < unitCnt; ++u) {
            struct stat st;
            const string name = outDir + "/subGraph_" + to_string(i) + "_"
                              + to_string(units[u].first) + "_" + to_string(units[u].second);
            if (stat(name.c_str(), &st) == 0) {
                subGraphBytes += static_cast<uint64_t>(st.st_size);
            }
        }
    }

    std::vector<size_t> unitMerged(unitCnt, 0);
    std::vector<size_t> unitCeiling(unitCnt, 0);
    std::vector<double> unitSeconds(unitCnt, 0.0);

    // schedule(dynamic): units still carry different loads even after sharding.
    #pragma omp parallel for schedule(dynamic) num_threads(concurrentMergers)
    for (size_t u = 0; u < unitCnt; ++u) {
        const time_t unitStart = time(nullptr);
        mergeUnit(units[u].first, units[u].second, unitOutPath(units[u].first, units[u].second),
                  unitInputFiles(units[u].first, units[u].second),
                  mergeBufElems, maxFanIn, unitMerged[u], unitCeiling[u]);
        unitSeconds[u] = double(time(nullptr) - unitStart);
    }

    // Fold the cross bucket's shard temps into its relations_*.bin.
    {
        const size_t crossRoute = static_cast<size_t>(par.threads);
        const size_t crossShards = shardsForRoute(crossRoute, par.threads);
        if (crossShards > 1) {
            const string finalPath = outDir + "/relations_" + to_string(crossRoute) + ".bin";
            ofstream dst(finalPath, std::ios::binary);
            if (!dst.is_open()) {
                cerr << "Error opening file: " << finalPath << endl;
            } else {
                for (size_t s = 0; s < crossShards; ++s) {
                    const string part = unitOutPath(crossRoute, s);
                    ifstream in(part, std::ios::binary);
                    // Streaming an empty rdbuf sets failbit on the destination, so guard it.
                    if (in.is_open() && in.peek() != std::ifstream::traits_type::eof()) {
                        dst << in.rdbuf();
                    }
                }
                dst.close();
                for (size_t s = 0; s < crossShards; ++s) {
                    std::remove(unitOutPath(crossRoute, s).c_str());
                }
            }
        }
    }

    size_t ceilingEdgeCnt = 0; // edges whose merged weight sits at the uint16 ceiling
    size_t mergedEdgeCnt = 0;  // distinct (id1,id2) pairs after merging
    size_t heaviest = 0;
    for (size_t u = 0; u < unitCnt; ++u) {
        mergedEdgeCnt += unitMerged[u];
        ceilingEdgeCnt += unitCeiling[u];
        if (unitMerged[u] > unitMerged[heaviest]) { heaviest = u; }
    }

    // Report the load split by EDGE COUNT, not just by time. On small inputs every unit
    // reports 0 s, which is how a 88%-of-everything unit stayed invisible before.
    {
        std::vector<size_t> order(unitCnt);
        for (size_t u = 0; u < unitCnt; ++u) { order[u] = u; }
        std::sort(order.begin(), order.end(),
                  [&](size_t a, size_t b) { return unitMerged[a] > unitMerged[b]; });
        cout << "[edges] unit load (top 3 of " << unitCnt << "), max share "
             << (mergedEdgeCnt ? 100.0 * static_cast<double>(unitMerged[heaviest])
                                     / static_cast<double>(mergedEdgeCnt) : 0.0)
             << "%:" << endl;
        for (size_t i = 0; i < order.size() && i < 3; ++i) {
            const size_t u = order[i];
            cout << "[edges]   route " << units[u].first << " shard " << units[u].second
                 << ": " << unitMerged[u] << " edges, " << unitSeconds[u] << " s" << endl;
        }
    }
    cout << "[disk] subGraph: " << (this->numOfGraph * unitCnt) << " files, "
         << humanBytes(subGraphBytes) << " freed" << endl;

    cout << "[edges] merged into " << mergedEdgeCnt << " distinct edges ("
         << humanBytes(static_cast<uint64_t>(mergedEdgeCnt) * sizeof(Relation)) << ")" << endl;

    if (ceilingEdgeCnt > 0) {
        cout << "Warning: " << ceilingEdgeCnt << " edges reached the weight ceiling ("
             << UINT16_MAX << "); their weights are saturated." << endl;
    } else {
        cout << "Edges at weight ceiling (" << UINT16_MAX << "): 0" << endl;
    }

    cout << "Query relation graph merged successfully" << endl;
    cout << "Time spent: " << double(time(nullptr) - before) << " seconds." << endl;
    return;
}

void GroupGenerator::mergeGraph_one(size_t processedReadCnt) {
    cout << "Merging subgraphs" << endl;
    time_t before = time(nullptr);

    // Subgraphs are now split by route, so this path walks the routes in turn and appends
    // each one's merge to relations.bin. Records are sorted within a route but the file as a
    // whole is segmented by route -- this dump exists for weight-histogram analysis, which
    // does not depend on global ordering. The merge itself stays a linear scan.
    const size_t routeCnt = static_cast<size_t>(par.threads) * 2 + 1;
    std::vector<std::pair<size_t, size_t>> units; // (route, shard)
    for (size_t r = 0; r < routeCnt; ++r) {
        const size_t shardCnt = shardsForRoute(r, par.threads);
        for (size_t s = 0; s < shardCnt; ++s) {
            units.emplace_back(r, s);
        }
    }
    const size_t maxFanIn = getMergeFanIn(par.threads, 1);
    const size_t mergeBufElems = getMergeBufferElems(std::min(this->unitStreamCount(), maxFanIn),
                                                     par.ramUsage, 1);
    cout << "Merge: " << units.size() << " units x " << this->numOfGraph
         << " subgraphs, sequential (fd soft limit " << getOpenFileLimit()
         << ", fan-in " << maxFanIn << ")" << endl;

    WriteBuffer<Relation> relationLog(outDir + "/relations.bin", 1024 * 1024);

    size_t ceilingEdgeCnt = 0; // edges whose merged weight sits at the uint16 ceiling

    for (size_t u = 0; u < units.size(); ++u) {
        const size_t route = units[u].first;
        const size_t shard = units[u].second;
        std::vector<std::string> subGraphFiles = unitInputFiles(route, shard);
        reduceSubGraphFanIn(subGraphFiles, maxFanIn, route, shard);

        const size_t streamCnt = subGraphFiles.size();
        std::vector<std::unique_ptr<RelationReader>> relationBuffers(streamCnt);
        std::vector<Relation> currentRelations(streamCnt);
        for (size_t i = 0; i < streamCnt; ++i) {
            relationBuffers[i].reset(new RelationReader(subGraphFiles[i], mergeBufElems));
            currentRelations[i] = relationBuffers[i]->getNext();
            if (currentRelations[i] == Relation()) {
                currentRelations[i] = Relation(UINT32_MAX, UINT32_MAX, UINT16_MAX);
            }
        }

        while (true) {
            Relation minRelation(UINT32_MAX, UINT32_MAX, 0);
            for (size_t i = 0; i < streamCnt; ++i) {
                if (currentRelations[i] < minRelation) {
                    minRelation = currentRelations[i];
                }
            }
            if (minRelation.id1 == UINT32_MAX) break;
            uint16_t totalWeight = 0;
            for (size_t i = 0; i < streamCnt; ++i) {
                if (currentRelations[i] == minRelation) {
                    addSat16(totalWeight, currentRelations[i].weight);
                    currentRelations[i] = relationBuffers[i]->getNext();
                    if (currentRelations[i] == Relation()) {
                        currentRelations[i] = Relation(UINT32_MAX, UINT32_MAX, UINT16_MAX);
                    }
                }
            }
            if (totalWeight == UINT16_MAX) ceilingEdgeCnt++;

            // Zero in place rather than copying from a helper's return value: a return-by-value
            // copy is free to leave the padding indeterminate, and this file is written byte
            // for byte.
            Relation rel;
            std::memset(static_cast<void *>(&rel), 0, sizeof(rel));
            rel.id1 = minRelation.id1;
            rel.id2 = minRelation.id2;
            rel.weight = totalWeight;
            relationLog.write(&rel);
        }

        relationBuffers.clear(); // close every input before unlinking it
        for (size_t i = 0; i < streamCnt; ++i) {
            std::remove(subGraphFiles[i].c_str());
        }
    }
    relationLog.flush();

    // This path exists to dump relations.bin for weight-histogram analysis, so a
    // silent saturation would quietly distort that histogram. Report it.
    if (ceilingEdgeCnt > 0) {
        cout << "Warning: " << ceilingEdgeCnt << " edges reached the weight ceiling ("
             << UINT16_MAX << "); their weights are saturated." << endl;
    } else {
        cout << "Edges at weight ceiling (" << UINT16_MAX << "): 0" << endl;
    }

    cout << "Query relation graph merged successfully" << endl;
    cout << "Time spent: " << double(time(nullptr) - before) << " seconds." << endl;
    return;
}

void GroupGenerator::makeGroups(int groupKmerThr,
                                size_t processedReadCnt,
                                unordered_map<uint32_t, unordered_set<uint32_t>> &groupInfo, 
                                vector<uint32_t> &queryGroupInfo) {
    cout << "Creating groups from relation file..." << endl;
    time_t beforeSearch = time(nullptr);

    DisjointSet ds(processedReadCnt);

    auto processFile = [&](const std::string& fname, DisjointSet& subDs) {
        ReadBuffer<Relation> rb(fname, 1024 * 1024);
        for (Relation r = rb.getNext(); !(r == Relation()); r = rb.getNext()) {
            uint32_t id1 = r.id1, id2 = r.id2;
            uint16_t weight = r.weight;
            if (static_cast<int>(weight) > groupKmerThr) {
                if (id1 == 0 || id2 == 0) continue;
                if (id1 > processedReadCnt || id2 > processedReadCnt) continue;
                subDs.unionSets(id1, id2);
            }
        }
    };


    #pragma omp parallel num_threads(par.threads)
    {
        int threadIdx = omp_get_thread_num();
        DisjointSet subDs(processedReadCnt);

        processFile(outDir + "/relations_" + std::to_string(threadIdx) + ".bin", subDs);

        subDs.flatten();

        #pragma omp critical
        {
            ds += subDs;
        }
    }

    #pragma omp parallel num_threads(par.threads)
    {
        int threadIdx = omp_get_thread_num();
        DisjointSet subDs(processedReadCnt);
        #pragma omp critical
        { subDs = ds;}

        processFile(outDir + "/relations_" + std::to_string(par.threads + threadIdx) + ".bin", subDs);

        subDs.flatten();

        #pragma omp critical
        { ds += subDs; }
    }

    {
        ReadBuffer<Relation> rb(outDir + "/relations_" + std::to_string(par.threads * 2) + ".bin", 1024 * 1024);
        for (Relation r = rb.getNext(); !(r == Relation()); r = rb.getNext()) {
            if (static_cast<int>(r.weight) > groupKmerThr) {
                ds.unionSets(r.id1, r.id2);
            }
        }
    }

    // Pass 1: map each component root -> min node ID in that component (canonical label)
    std::unordered_map<uint32_t, uint32_t> rootToMin;
    for (uint32_t i = 1; i < ds.parent.size(); ++i) {
        if (!ds.grouped[i]) continue;
        uint32_t root = ds.find(i);
        auto it = rootToMin.find(root);
        if (it == rootToMin.end() || i < it->second) {
            rootToMin[root] = i;
        }
    }
    // Pass 2: assign canonical group IDs
    for (uint32_t queryId = 1; queryId < ds.parent.size(); queryId++) {
        if (ds.grouped[queryId]){
            uint32_t groupId = rootToMin[ds.find(queryId)];
            groupInfo[groupId].insert(queryId);
            queryGroupInfo[queryId] = groupId;
        }
    }

    cout << "Query groups created successfully: " << groupInfo.size() << " groups." << endl;
    cout << "Time spent: " << double(time(nullptr) - beforeSearch) << " seconds." << endl;
}

// Phase 1.5 -- see the declaration in GroupGenerator.h for the rationale.
//
// A "unit" is a Phase-1 core group (labelled by its minimum read id) or, for a read that Phase 1
// left alone, the read itself. Group ids are read ids of grouped reads, so a singleton's own id can
// never collide with a group id. Weak edges are those in (weakThr, coreThr] -- above Phase 2's
// floor, below the core threshold. Pairs where both sides are singletons are skipped: distinct edges
// are merged upstream, so two reads share exactly one edge and can never reach a support of 2. That
// also keeps the counting map bounded to pairs involving at least one multi-read group.
//
// Pass 1 counts weak EDGES per unit pair. With supportRatio == 0 that count, against the floor, is
// the whole rule. With supportRatio > 0 the requirement becomes a fraction of the smaller unit's
// read count and pass 2 recounts the qualifying pairs by DISTINCT READS on the smaller side. Both
// changes are needed together:
//   - the requirement has to scale, because chance links between units u and v grow with
//     |u| * |v|; a fixed count is met automatically once coverage is high.
//   - the count has to be distinct reads, because one repeat-bearing read linked to 50 reads of the
//     other unit yields 50 edges but explains only one read -- exactly the case being excluded.
// Pass 1's edge count is an upper bound on pass 2's distinct-read count, so it is a sound prefilter:
// a pair that fails on edges cannot pass on distinct reads.
//
// Unit sizes are the STATIC Phase-1 sizes, never the disjoint set's evolving component sizes. The
// support map is an unordered_map, so the order in which pairs are merged is unspecified; a
// size test against a changing component would make the output depend on that order.
void GroupGenerator::mergeBySupport(int coreThr,
                                   int weakThr,
                                   float supportRatio,
                                   size_t maxUnitReads,
                                   size_t processedReadCnt,
                                   unordered_map<uint32_t, unordered_set<uint32_t>> &groupInfo,
                                   vector<uint32_t> &queryGroupInfo) {
    const bool useRatio = (supportRatio > 0.0f);
    if (useRatio) {
        cout << "Phase 1.5: merging units whose smaller side has >= max(" << MERGE_SUPPORT_FLOOR
             << ", " << supportRatio << " x its read count) distinct reads carrying a weak link"
             << " (weight in (" << weakThr << ", " << coreThr << "])..." << endl;
    } else {
        cout << "Phase 1.5: merging units with >= " << MERGE_SUPPORT_FLOOR
             << " weak links (weight in (" << weakThr << ", " << coreThr << "])..." << endl;
    }
    time_t t0 = time(nullptr);

    const size_t groupsBefore = groupInfo.size();

    // unit[i]: Phase-1 group id, or the read's own id when Phase 1 left it alone.
    std::vector<uint32_t> unit(processedReadCnt + 1, 0);
    for (uint32_t i = 1; i <= processedReadCnt; ++i) {
        unit[i] = queryGroupInfo[i] ? queryGroupInfo[i] : i;
    }

    // The disjoint set starts from the Phase-1 components and the supported pairs are added to it
    // further down. It is built here, before the counting passes, because `compSize` is indexed by
    // its roots.
    DisjointSet ds(processedReadCnt);
    for (uint32_t i = 1; i <= processedReadCnt; ++i) {
        if (queryGroupInfo[i] != 0 && queryGroupInfo[i] != i) {
            ds.unionSets(i, queryGroupInfo[i]);
        }
    }

    // Reads per component, indexed by disjoint-set root. One array serves both roles: until the
    // merge loop runs, every root is still a Phase-1 root, so a lookup returns the static unit size
    // the ratio rule needs; from then on the loop keeps it current so the size gate can see what a
    // component has already grown to. Keeping two arrays instead would cost another 4 B per read,
    // which is ~666 MB on a 166 M-read run.
    // The order matters and is load-bearing: `requiredSupport` and `smallSideOf` must see static
    // sizes, or the support rule would depend on the merge order. Both are only called while
    // building `candidates`, which finishes before the merge loop starts.
    std::vector<uint32_t> compSize(processedReadCnt + 1, 0u);
    for (uint32_t i = 1; i <= processedReadCnt; ++i) { compSize[ds.find(i)]++; }
    const auto sizeOf = [&](uint32_t node) -> uint32_t { return compSize[ds.find(node)]; };

    // Counting map memory guard. Entries are ~24 B in unordered_map, so the cap below is a few GB.
    // On overflow, existing keys keep accumulating and new ones are dropped -- support can then only
    // be under-counted, never over-counted, so a merge is never invented. The count is reported.
    const size_t SUPPORT_PAIR_CAP = 200000000;
    const size_t perThreadCap = std::max<size_t>(1, SUPPORT_PAIR_CAP / std::max(1, par.threads));

    std::unordered_map<uint64_t, uint32_t> support;
    size_t weakEdgeCnt = 0, droppedPairs = 0;

    auto packPair = [](uint32_t a, uint32_t b) -> uint64_t {
        return (a < b) ? ((uint64_t)a << 32 | b) : ((uint64_t)b << 32 | a);
    };

    auto countFile = [&](const std::string& fname,
                         std::unordered_map<uint64_t, uint32_t>& local,
                         size_t cap, size_t& weakCnt, size_t& dropped) {
        ReadBuffer<Relation> rb(fname, 1024 * 1024);
        for (Relation r = rb.getNext(); !(r == Relation()); r = rb.getNext()) {
            const uint32_t id1 = r.id1, id2 = r.id2;
            if (id1 == 0 || id2 == 0) continue;
            if (id1 > processedReadCnt || id2 > processedReadCnt) continue;
            const int w = static_cast<int>(r.weight);
            if (w <= weakThr || w > coreThr) continue;
            if (queryGroupInfo[id1] == 0 && queryGroupInfo[id2] == 0) continue;
            const uint32_t u = unit[id1], v = unit[id2];
            if (u == v) continue;
            weakCnt++;
            const uint64_t key = packPair(u, v);
            auto it = local.find(key);
            if (it != local.end()) {
                it->second++;
            } else if (local.size() < cap) {
                local.emplace(key, 1u);
            } else {
                dropped++;
            }
        }
    };

    #pragma omp parallel num_threads(par.threads)
    {
        const int threadIdx = omp_get_thread_num();
        std::unordered_map<uint64_t, uint32_t> local;
        size_t localWeak = 0, localDropped = 0;

        countFile(outDir + "/relations_" + std::to_string(threadIdx) + ".bin",
                  local, perThreadCap, localWeak, localDropped);
        countFile(outDir + "/relations_" + std::to_string(par.threads + threadIdx) + ".bin",
                  local, perThreadCap, localWeak, localDropped);

        #pragma omp critical
        {
            for (const auto& kv : local) support[kv.first] += kv.second;
            weakEdgeCnt += localWeak;
            droppedPairs += localDropped;
        }
    }
    countFile(outDir + "/relations_" + std::to_string(par.threads * 2) + ".bin",
              support, SUPPORT_PAIR_CAP, weakEdgeCnt, droppedPairs);

    // Support required of a pair. Without the ratio this is the bare floor, which is the rule the
    // benchmark operating point was measured with.
    auto requiredSupport = [&](uint32_t u, uint32_t v) -> uint32_t {
        if (!useRatio) { return MERGE_SUPPORT_FLOOR; }
        const uint32_t smaller = std::min(sizeOf(u), sizeOf(v));
        const double need = std::ceil(static_cast<double>(supportRatio) * static_cast<double>(smaller));
        const uint32_t needed = (need >= static_cast<double>(UINT32_MAX))
            ? UINT32_MAX : static_cast<uint32_t>(need);
        return std::max(MERGE_SUPPORT_FLOOR, needed);
    };

    // Which side of a pair the distinct reads are counted on. Sizes decide it; equal sizes fall back
    // to the smaller unit id so the choice does not depend on map order.
    auto smallSideOf = [&](uint32_t u, uint32_t v) -> uint32_t {
        const uint32_t su = sizeOf(u), sv = sizeOf(v);
        if (su != sv) { return (su < sv) ? u : v; }
        return std::min(u, v);
    };

    // Pairs that cleared pass 1. With the ratio on, this is only a prefilter -- pass 2 recounts them
    // by distinct reads, which can only be smaller.
    // The size gate is applied here as well as in the merge loop. Here it is the cheap half: a pair
    // whose Phase-1 units are already too big can never qualify, whatever happens later, so dropping
    // it now also keeps it out of pass 2's distinct-read counting.
    std::unordered_map<uint64_t, uint32_t> candidates; // pair -> required support
    size_t gatedStatic = 0, gatedStaticReads = 0;
    for (const auto& kv : support) {
        const uint32_t u = static_cast<uint32_t>(kv.first >> 32);
        const uint32_t v = static_cast<uint32_t>(kv.first & 0xFFFFFFFFull);
        const uint32_t need = requiredSupport(u, v);
        if (kv.second < need) { continue; }
        if (maxUnitReads > 0
            && (sizeOf(u) > maxUnitReads || sizeOf(v) > maxUnitReads)) {
            gatedStatic++;
            gatedStaticReads += static_cast<size_t>(sizeOf(u)) + sizeOf(v);
            continue;
        }
        candidates.emplace(kv.first, need);
    }
    if (maxUnitReads > 0) {
        cout << "Phase 1.5: size gate " << maxUnitReads << " reads/unit dropped " << gatedStatic
             << " supported pairs spanning " << gatedStaticReads << " unit-reads" << endl;
    }

    // Pass 2: distinct reads on the smaller side, for candidate pairs only.
    std::unordered_set<uint64_t> satisfied;
    size_t droppedReads = 0;
    if (useRatio && !candidates.empty()) {
        // Read-set memory guard, same contract as SUPPORT_PAIR_CAP: on overflow new reads are
        // dropped, so a count can only fall short and a merge is never invented.
        const size_t PASS2_READ_CAP = 200000000;
        const size_t perThreadReadCap = std::max<size_t>(1, PASS2_READ_CAP / std::max(1, par.threads));

        // A pair is satisfied as soon as `required` distinct reads are seen, so a thread that stops
        // inserting at that point cannot make any pair's union fall short: no thread stops earlier.
        auto distinctFile = [&](const std::string& fname,
                                std::unordered_map<uint64_t, std::unordered_set<uint32_t>>& local,
                                std::unordered_set<uint64_t>& localSatisfied,
                                size_t cap, size_t& held, size_t& dropped) {
            ReadBuffer<Relation> rb(fname, 1024 * 1024);
            for (Relation r = rb.getNext(); !(r == Relation()); r = rb.getNext()) {
                const uint32_t id1 = r.id1, id2 = r.id2;
                if (id1 == 0 || id2 == 0) continue;
                if (id1 > processedReadCnt || id2 > processedReadCnt) continue;
                const int w = static_cast<int>(r.weight);
                if (w <= weakThr || w > coreThr) continue;
                if (queryGroupInfo[id1] == 0 && queryGroupInfo[id2] == 0) continue;
                const uint32_t u = unit[id1], v = unit[id2];
                if (u == v) continue;
                const uint64_t key = packPair(u, v);
                const auto cand = candidates.find(key);
                if (cand == candidates.end()) continue;
                if (localSatisfied.count(key)) continue;
                const uint32_t small = smallSideOf(u, v);
                const uint32_t readId = (unit[id1] == small) ? id1 : id2;
                if (held >= cap) { dropped++; continue; }
                auto& reads = local[key];
                if (reads.insert(readId).second) { held++; }
                if (reads.size() >= cand->second) {
                    localSatisfied.insert(key);
                    held -= reads.size();
                    local.erase(key);
                }
            }
        };

        size_t heldReads = 0;
        std::unordered_map<uint64_t, std::unordered_set<uint32_t>> partial;
        #pragma omp parallel num_threads(par.threads)
        {
            const int threadIdx = omp_get_thread_num();
            std::unordered_map<uint64_t, std::unordered_set<uint32_t>> local;
            std::unordered_set<uint64_t> localSatisfied;
            size_t localHeld = 0, localDropped = 0;

            distinctFile(outDir + "/relations_" + std::to_string(threadIdx) + ".bin",
                         local, localSatisfied, perThreadReadCap, localHeld, localDropped);
            distinctFile(outDir + "/relations_" + std::to_string(par.threads + threadIdx) + ".bin",
                         local, localSatisfied, perThreadReadCap, localHeld, localDropped);

            #pragma omp critical
            {
                satisfied.insert(localSatisfied.begin(), localSatisfied.end());
                for (const auto& kv : local) {
                    partial[kv.first].insert(kv.second.begin(), kv.second.end());
                }
                heldReads += localHeld;
                droppedReads += localDropped;
            }
        }
        // Unions can have crossed the requirement even when no single thread's share did.
        for (auto it = partial.begin(); it != partial.end(); ) {
            if (it->second.size() >= candidates[it->first]) {
                satisfied.insert(it->first);
                it = partial.erase(it);
            } else {
                ++it;
            }
        }
        distinctFile(outDir + "/relations_" + std::to_string(par.threads * 2) + ".bin",
                     partial, satisfied, PASS2_READ_CAP, heldReads, droppedReads);

        cout << "Phase 1.5: " << candidates.size() << " candidate pairs, " << satisfied.size()
             << " reached the distinct-read requirement (" << heldReads
             << " reads still held for the rest)";
        if (droppedReads) {
            cout << " [read cap hit, " << droppedReads
                 << " observations dropped -- support undercounted]";
        }
        cout << endl;
    }

    // Merge the supported pairs into the Phase-1 components built above.
    //
    // The order is fixed on purpose. The gate below asks how big a component already is, so with
    // the plain `unordered_map` walk the outcome would depend on the map's bucket layout. Sorting
    // by support puts the best-evidenced pairs first, which is also where the size budget should go.
    const size_t LOG_BUCKETS = 8;
    std::vector<std::pair<uint64_t, uint32_t>> ordered; // (pair key, support observed in pass 1)
    ordered.reserve(candidates.size());
    for (const auto& kv : candidates) {
        if (useRatio && satisfied.count(kv.first) == 0) continue;
        const auto found = support.find(kv.first);
        ordered.emplace_back(kv.first, (found == support.end()) ? 0u : found->second);
    }
    std::sort(ordered.begin(), ordered.end(),
              [](const std::pair<uint64_t, uint32_t>& a, const std::pair<uint64_t, uint32_t>& b) {
                  if (a.second != b.second) { return a.second > b.second; }
                  return a.first < b.first;
              });

    // Units touched by an accepted merge. A bitmap, not a set of ids: the id space is the read
    // space, so this is one bit per read where a hash set would be tens of bytes per merged pair.
    std::vector<bool> unitMerged(processedReadCnt + 1, false);
    std::vector<size_t> mergeSizeHist(LOG_BUCKETS, 0); // sizes the gate actually saw, both sides
    size_t mergedPairs = 0, gatedDynamic = 0;

    const auto bucketOf = [LOG_BUCKETS](uint32_t sz) -> size_t {
        // bucket 0 is the singleton component; from 2 on it is 1 + floor(log10(sz)), capped.
        if (sz < 2) { return 0; }
        size_t bucket = 1;
        for (uint32_t bound = 10; bucket + 1 < LOG_BUCKETS && sz >= bound; ++bucket) {
            bound *= 10;
        }
        return bucket;
    };

    for (const auto& entry : ordered) {
        const uint32_t u = static_cast<uint32_t>(entry.first >> 32);
        const uint32_t v = static_cast<uint32_t>(entry.first & 0xFFFFFFFFull);
        const uint32_t ru = ds.find(u), rv = ds.find(v);
        if (ru == rv) { continue; } // an earlier merge already joined them
        // Re-checked here, not just on the Phase-1 sizes: a component that grew past the bound has
        // to stop growing, and that is what makes the bound hold. Both sides are at most
        // maxUnitReads, so no component can end up with 2 * maxUnitReads reads or more, and one
        // that reaches the bound is out of the running for every remaining pair.
        if (maxUnitReads > 0
            && (compSize[ru] > maxUnitReads || compSize[rv] > maxUnitReads)) {
            gatedDynamic++;
            continue;
        }
        mergeSizeHist[bucketOf(compSize[ru])]++;
        mergeSizeHist[bucketOf(compSize[rv])]++;
        const uint64_t merged = static_cast<uint64_t>(compSize[ru]) + compSize[rv];
        ds.unionSets(u, v);
        compSize[ds.find(u)] = static_cast<uint32_t>(merged); // read ids are uint32, so is any sum
        unitMerged[u] = true;
        unitMerged[v] = true;
        mergedPairs++;
    }

    // A Phase-1 group of exactly one read (its own label) never entered the disjoint set above, so
    // re-mark those reads as grouped to keep them out of Phase 2's singleton pass.
    for (uint32_t i = 1; i <= processedReadCnt; ++i) {
        if (queryGroupInfo[i] != 0) ds.grouped[i] = true;
    }

    // Relabel every component by its minimum read id, as Phase 1 does.
    groupInfo.clear();
    std::fill(queryGroupInfo.begin(), queryGroupInfo.end(), 0u);
    std::unordered_map<uint32_t, uint32_t> rootToMin;
    for (uint32_t i = 1; i < ds.parent.size(); ++i) {
        if (!ds.grouped[i]) continue;
        const uint32_t root = ds.find(i);
        auto it = rootToMin.find(root);
        if (it == rootToMin.end() || i < it->second) rootToMin[root] = i;
    }
    for (uint32_t i = 1; i < ds.parent.size(); ++i) {
        if (!ds.grouped[i]) continue;
        const uint32_t groupId = rootToMin[ds.find(i)];
        groupInfo[groupId].insert(i);
        queryGroupInfo[i] = groupId;
    }

    cout << "Phase 1.5: " << weakEdgeCnt << " weak edges over " << support.size()
         << " unit pairs, " << mergedPairs << " pairs merged"
         << "; groups " << groupsBefore << " -> " << groupInfo.size();
    if (droppedPairs) {
        cout << " [pair cap hit, " << droppedPairs << " observations dropped -- support undercounted]";
    }
    cout << endl;

    // Read mass this phase moved, and the sizes it moved it from. Pair counts alone hide the damage:
    // a few merges between large units cost more purity than many merges between small ones, because
    // every read of a wrongly merged unit becomes impure. The share below is the worst case that
    // purity can lose here, so it is readable before any grading run.
    {
        // Reads sitting in a component that this phase enlarged. Counted per component, from the
        // final sizes, so a component merged several times is counted once.
        std::unordered_set<uint32_t> mergedRoots;
        size_t mergedUnits = 0;
        for (uint32_t i = 1; i <= processedReadCnt; ++i) {
            if (!unitMerged[i]) continue;
            mergedUnits++;
            mergedRoots.insert(ds.find(i));
        }
        size_t mergedReads = 0, maxMergedComp = 0;
        for (const uint32_t root : mergedRoots) {
            mergedReads += compSize[root];
            if (compSize[root] > maxMergedComp) { maxMergedComp = compSize[root]; }
        }

        size_t groupedReads = 0, maxGroupSize = 0;
        std::vector<size_t> groupSizes;
        groupSizes.reserve(groupInfo.size());
        for (const auto& g : groupInfo) {
            groupSizes.push_back(g.second.size());
            groupedReads += g.second.size();
            if (g.second.size() > maxGroupSize) { maxGroupSize = g.second.size(); }
        }
        const size_t topN = std::min<size_t>(5, groupSizes.size());
        std::partial_sort(groupSizes.begin(), groupSizes.begin() + topN, groupSizes.end(),
                          std::greater<size_t>());

        const double massShare = groupedReads
            ? (100.0 * static_cast<double>(mergedReads) / static_cast<double>(groupedReads)) : 0.0;
        cout << "Phase 1.5: " << mergedUnits << " units merged, holding " << mergedReads
             << " reads (" << massShare << "% of grouped reads -- upper bound on the purity this"
             << " phase can cost)" << endl;

        // Sizes as the gate saw them: one entry per side of every accepted merge, taken from the
        // component size at that moment. This is the distribution --merge-max-unit-reads cuts, so
        // it is the table to read when picking a value for it.
        cout << "Phase 1.5: merged component sizes";
        static const char* BUCKET_NAME[] = {
            "1", "2-9", "10-99", "100-999", "1e3-1e4", "1e4-1e5", "1e5-1e6", "1e6+"
        };
        for (size_t b = 0; b < LOG_BUCKETS; ++b) {
            cout << " | " << BUCKET_NAME[b] << ": " << mergeSizeHist[b];
        }
        cout << endl;
        if (maxUnitReads > 0) {
            cout << "Phase 1.5: size gate also blocked " << gatedDynamic
                 << " pairs at merge time (component had grown past the bound)" << endl;
        }

        // Two different numbers on purpose. The first covers every group, including Phase-1 groups
        // this phase never touched -- those can be larger than the bound and are left alone. The
        // second is the one the bound governs: no component this phase enlarged reaches
        // 2 * maxUnitReads reads.
        cout << "Phase 1.5: max group size " << maxGroupSize << ", top " << topN << ":";
        for (size_t i = 0; i < topN; ++i) { cout << " " << groupSizes[i]; }
        cout << endl;
        cout << "Phase 1.5: largest component this phase enlarged: " << maxMergedComp;
        if (maxUnitReads > 0) {
            cout << " (bound 2 x " << maxUnitReads << " = " << (2 * maxUnitReads) << ")";
        }
        cout << endl;
    }

    cout << "Phase 1.5 done: " << double(time(nullptr) - t0) << " s" << endl;
}

void GroupGenerator::makeGroupsPhase2(
    int groupKmerThr,
    size_t processedReadCnt,
    const std::vector<bool>& isSingleton,
    unordered_map<uint32_t, unordered_set<uint32_t>>& groupInfo,
    vector<uint32_t>& queryGroupInfo)
{
    cout << "Phase 2: linking singletons..." << endl;
    time_t t0 = time(nullptr);

    DisjointSet ds(processedReadCnt);

    auto processFile = [&](const std::string& fname, DisjointSet& subDs) {
        ReadBuffer<Relation> rb(fname, 1024 * 1024);
        for (Relation r = rb.getNext(); !(r == Relation()); r = rb.getNext()) {
            uint32_t id1 = r.id1, id2 = r.id2;
            if (id1 == 0 || id2 == 0) continue;
            if (id1 > processedReadCnt || id2 > processedReadCnt) continue;
            if (!isSingleton[id1] || !isSingleton[id2]) continue;
            if (static_cast<int>(r.weight) > groupKmerThr) {
                subDs.unionSets(id1, id2);
            }
        }
    };

    #pragma omp parallel num_threads(par.threads)
    {
        int threadIdx = omp_get_thread_num();
        DisjointSet subDs(processedReadCnt);
        processFile(outDir + "/relations_" + std::to_string(threadIdx) + ".bin", subDs);
        subDs.flatten();
        #pragma omp critical
        { ds += subDs; }
    }

    #pragma omp parallel num_threads(par.threads)
    {
        int threadIdx = omp_get_thread_num();
        DisjointSet subDs(processedReadCnt);
        #pragma omp critical
        { subDs = ds; }
        processFile(outDir + "/relations_" + std::to_string(par.threads + threadIdx) + ".bin", subDs);
        subDs.flatten();
        #pragma omp critical
        { ds += subDs; }
    }

    {
        ReadBuffer<Relation> rb(outDir + "/relations_" + std::to_string(par.threads * 2) + ".bin", 1024 * 1024);
        for (Relation r = rb.getNext(); !(r == Relation()); r = rb.getNext()) {
            if (r.id1 == 0 || r.id2 == 0) continue;
            if (r.id1 > processedReadCnt || r.id2 > processedReadCnt) continue;
            if (!isSingleton[r.id1] || !isSingleton[r.id2]) continue;
            if (static_cast<int>(r.weight) > groupKmerThr) {
                ds.unionSets(r.id1, r.id2);
            }
        }
    }

    std::unordered_map<uint32_t, uint32_t> rootToMin;
    for (uint32_t i = 1; i < ds.parent.size(); ++i) {
        if (!ds.grouped[i]) continue;
        uint32_t root = ds.find(i);
        auto it = rootToMin.find(root);
        if (it == rootToMin.end() || i < it->second) rootToMin[root] = i;
    }
    for (uint32_t queryId = 1; queryId < ds.parent.size(); queryId++) {
        if (ds.grouped[queryId]) {
            uint32_t groupId = rootToMin[ds.find(queryId)];
            groupInfo[groupId].insert(queryId);
            queryGroupInfo[queryId] = groupId;
        }
    }

    cout << "Phase 2 done: " << double(time(nullptr) - t0) << " s" << endl;
}

void GroupGenerator::saveGroupsToFile(const unordered_map<uint32_t, unordered_set<uint32_t>> &groupInfo,
                                      const vector<uint32_t> &queryGroupInfo) {
    // save group in txt file
    const string& groupInfoFileName = outDir + "/groups";
    ofstream outFile1(groupInfoFileName);
    if (!outFile1.is_open()) {
        cerr << "Error opening file: " << groupInfoFileName << endl;
        return;
    }

    for (const auto& [groupId, queryIds] : groupInfo) {
        outFile1 << groupId << "\t";
        for (const auto& queryId : queryIds) {
            outFile1 << queryId << "\t";
        }
        outFile1 << endl;
    }
    outFile1.close();
    cout << "Query group saved to " << groupInfoFileName << " successfully." << endl;
    

    const string& queryGroupInfoFileName = outDir + "/groupMap";
    ofstream outFile2(queryGroupInfoFileName);
    if (!outFile2.is_open()) {
        cerr << "Error opening file: " << queryGroupInfoFileName << endl;
        return;
    }

    for (size_t i = 1; i < queryGroupInfo.size(); ++i) {
        outFile2 << i << "\t" << queryGroupInfo[i] << "\n";
    }
    outFile2.close();
    cout << "Query group saved to " << queryGroupInfoFileName << " successfully." << endl;

    return;
}