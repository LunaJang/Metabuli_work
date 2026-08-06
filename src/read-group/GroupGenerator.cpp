#include "GroupGenerator.h"
#include "FileUtil.h"
#include "QueryIndexer.h"
#include "common.h"
#include "Kmer.h"
#include <sys/stat.h>
#include <queue>
#include <iomanip>

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
        std::vector<uint64_t> edgeWeightHist;
        mergeGraph(edgeWeightHist);

        // The knee is always computed, even when unused, so its value can be compared against
        // whatever is actually applied. edgeWeightHist is built by mergeGraph regardless, so
        // this costs nothing.
        int kneeEdge = kneeThreshold(edgeWeightHist, par.minEdgeWeight);
        if (kneeEdge > par.minEdgeWeight && par.kneeScale != 1.0f) {
            // --knee-scale < 1.0 lowers the core threshold; clamp so Phase 2 stays enabled.
            kneeEdge = static_cast<int>(kneeEdge * par.kneeScale + 0.5f);
            if (kneeEdge < par.minEdgeWeight + 1) kneeEdge = par.minEdgeWeight + 1;
        }

        // An edge weight counts the k-mers two reads share, so weight / kmersPerRead is the
        // fraction of a read the two overlap by. That makes the threshold a property of read
        // geometry rather than of the weight distribution's shape -- unlike the knee, which
        // tracks the tail's extent and therefore shifts with read length and coverage.
        const double kmersPerRead = (processedReadCnt > 0)
            ? static_cast<double>(totalFilteredKmers) / static_cast<double>(processedReadCnt)
            : 0.0;
        const int ratioEdge = (par.minOverlapRatio > 0.0f && kmersPerRead > 0.0)
            ? static_cast<int>(par.minOverlapRatio * kmersPerRead + 0.5)
            : 0;

        int effectiveCoreEdge = par.coreEdgeWeight;
        if (ratioEdge > par.minEdgeWeight) {
            effectiveCoreEdge = ratioEdge;
            cout << "Core threshold: " << effectiveCoreEdge << " (overlap ratio "
                 << par.minOverlapRatio << " x " << kmersPerRead << " k-mers/read)"
                 << " [knee would give " << kneeEdge << "]" << endl;
        } else if (par.minOverlapRatio > 0.0f) {
            cout << "Core threshold: overlap ratio " << par.minOverlapRatio << " x "
                 << kmersPerRead << " k-mers/read = " << ratioEdge
                 << " is not above --min-edge " << par.minEdgeWeight
                 << "; falling back" << endl;
        }
        if (effectiveCoreEdge == par.coreEdgeWeight) {
            if (kneeEdge > par.minEdgeWeight) {
                effectiveCoreEdge = kneeEdge;
                cout << "Auto coreEdgeWeight (knee x" << par.kneeScale << "): "
                     << effectiveCoreEdge << endl;
            } else {
                cout << "Knee: insufficient data, using --core-edge: "
                     << effectiveCoreEdge << endl;
            }
        }

        // Phase 1: form core groups with strong edges
        makeGroups(effectiveCoreEdge, processedReadCnt, groupInfo, queryGroupInfo);

        // Phase 2: link Phase-1 singletons among themselves with weak edges
        if (effectiveCoreEdge > par.minEdgeWeight) {
            std::vector<bool> isSingleton(processedReadCnt + 1, false);
            for (uint32_t i = 1; i <= processedReadCnt; i++) {
                if (queryGroupInfo[i] == 0) isSingleton[i] = true;
            }
            makeGroupsPhase2(par.minEdgeWeight, processedReadCnt, isSingleton, groupInfo, queryGroupInfo);
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
                if (int64_t(qKmers.buffer[lookingPos].qInfo.pos) < int(matchBuffer.buffer[matchIdx].second) - par.neighborKmers){
                    qKmers.buffer[storePos++] = qKmers.buffer[lookingPos++];
                }
                // next target check
                else if(int(matchBuffer.buffer[matchIdx].second) + par.neighborKmers < int64_t(qKmers.buffer[lookingPos].qInfo.pos)){
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

void GroupGenerator::makeSubGraph(size_t processedReadCnt) {
    cout << "Connecting reads with shared k-mer..." << endl;
    time_t beforeSearch = time(nullptr);

    const size_t RELATION_THRESHOLD = getRelationThreshold(par.threads, par.ramUsage);
    cout << "Flush threshold: " << RELATION_THRESHOLD << " pairs/thread (memory budget "
         << humanBytes(getMemoryBudgetBytes(par.ramUsage)) << ", --max-ram "
         << par.ramUsage << " GiB)" << endl;
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

    // Skip thresholds, resolved once. The ratio threshold is floored at 2: readCnt * ratio
    // truncates to 0 whenever it is below 1, and "m > 0" skips every k-mer including m = 1,
    // which yields zero edges and zero groups. That is never the intent.
    const size_t ratioThr = (par.maxKmerFreqRatio > 0.0f)
        ? std::max<size_t>(2, static_cast<size_t>(static_cast<double>(processedReadCnt) * par.maxKmerFreqRatio))
        : 0;
    const size_t absThr = (par.maxKmerReads > 0) ? static_cast<size_t>(par.maxKmerReads) : 0;
    cout << "[skip] thresholds: --max-kmer-reads "
         << (absThr > 0 ? to_string(absThr) : string("off"))
         << ", --max-kmer-freq-ratio " << par.maxKmerFreqRatio
         << " (m > " << (ratioThr > 0 ? to_string(ratioThr) : string("off")) << ")" << endl;

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

        std::vector<DeltaIdxReader*> deltaIdxReaders;
        std::vector<Kmer> currentKmers;
        for (size_t i = 0; i < this->numOfSplits; i++) {
            string diffFile = outDir + "/kmer_delta_" + to_string(i) + "_" + to_string(threadIdx);
            string infoFile = outDir + "/kmer_info_"  + to_string(i) + "_" + to_string(threadIdx);
            DeltaIdxReader* reader = new DeltaIdxReader(diffFile, infoFile, 1024 * 1024, 1024 * 1024);
            deltaIdxReaders.push_back(reader);
            currentKmers.push_back(reader->next());
        }

        vector<uint32_t> currentQueryIds;
        currentQueryIds.reserve(1024);

        while (true) {
            // Find the smallest k-mer
            uint64_t minKmer = UINT64_MAX;
            for (size_t file = 0; file < this->numOfSplits; ++file) {
                if (!deltaIdxReaders[file]->isCompleted()) {
                    minKmer = min(minKmer, currentKmers[file].value);
                }
            }
            if (minKmer == UINT64_MAX) break;

            currentQueryIds.clear();
            for (size_t file = 0; file < this->numOfSplits; ++file) {
                while (currentKmers[file].value == minKmer) {
                    uint32_t seqId = currentKmers[file].tInfo.taxId; // query ID is stored in taxId field
                    if (seqId != UINT32_MAX && seqId <= processedReadCnt) {
                        currentQueryIds.emplace_back(seqId);
                    }
                    currentKmers[file] = deltaIdxReaders[file]->next();
                    if (deltaIdxReaders[file]->isCompleted()) {
                        currentKmers[file].value = UINT64_MAX;
                        break;
                    }
                }
            }

            std::sort(currentQueryIds.begin(), currentQueryIds.end());
            auto last = std::unique(currentQueryIds.begin(), currentQueryIds.end());
            currentQueryIds.erase(last, currentQueryIds.end());
            const size_t m = currentQueryIds.size();

            // Histogram every k-mer, skipped or not, so the summary table can answer
            // "what would a different cap have cost".
            if (m >= 2) {
                const size_t bucket = mHistBucket(m);
                localHistKmers[bucket]++;
                localHistPairs[bucket] += static_cast<uint64_t>(m) * static_cast<uint64_t>(m - 1) / 2;
            }

            if ((absThr > 0 && m > absThr) || (ratioThr > 0 && m > ratioThr)) {
                skippedPart << minKmer << "\t" << m << "\n";
                localSkippedCnt++;
                localSumM += m;
                if (m > localMaxM) { localMaxM = m; }
                localPairEst += static_cast<double>(m) * static_cast<double>(m - 1) * 0.5;
                continue;
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
            }
        }
        if (!pair2weight.empty()) {
            size_t counter_now = counter.fetch_add(1, std::memory_order_relaxed);
            localEmittedEdges += pair2weight.size();
            saveSubGraphToFile(pair2weight, counter_now, processedReadCnt);
        } else {
            cout << "Thread " << threadIdx << " has no relations to write." << endl;
        }
        for (size_t file = 0; file < this->numOfSplits; file++) {
            delete deltaIdxReaders[file];
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
    } else if (absThr == 0 && ratioThr == 0) {
        cout << "[skip] 0 k-mers skipped (both thresholds off)" << endl;
    } else {
        cout << "[skip] 0 k-mers skipped (no k-mer exceeded the thresholds)" << endl;
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
    cout << "[edges] emitted " << emittedEdgeCnt << " edge records into "
         << this->numOfGraph << " subgraphs ("
         << humanBytes(static_cast<uint64_t>(emittedEdgeCnt) * sizeof(Relation)) << " before merge)" << endl;
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
        byUnit[route][shardOf(id1, byUnit[route].size())].emplace_back(id1, id2, weight);
    }

    // Empty units still get an empty file: the merge's input list is then simply
    // subGraph_{0..numOfGraph-1}_{route}_{shard} with no existence checks.
    for (size_t route = 0; route < routeCnt; ++route) {
        for (size_t shard = 0; shard < byUnit[route].size(); ++shard) {
            const string subGraphFileName = outDir + "/subGraph_" + to_string(counter_now)
                                          + "_" + to_string(route) + "_" + to_string(shard);
            FILE * outFile = fopen(subGraphFileName.c_str(), "wb");
            if (!outFile) {
                cerr << "Error opening file: " << subGraphFileName << endl;
                continue;
            }
            std::vector<Relation> & relations = byUnit[route][shard];
            sort(relations.begin(), relations.end(), Relation::compare);
            fwrite(relations.data(), sizeof(Relation), relations.size(), outFile);
            fclose(outFile);
        }
    }
}

size_t GroupGenerator::mergeSubGraphBatch(const std::vector<std::string> & inputs,
                                          const std::string & output,
                                          size_t bufElems) {
    std::vector<ReadBuffer<Relation> *> readers(inputs.size());
    for (size_t i = 0; i < inputs.size(); ++i) {
        readers[i] = new ReadBuffer<Relation>(inputs[i], bufElems);
    }

    WriteBuffer<Relation> out(output, 1024 * 1024);

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
        Relation rel(minRelation.id1, minRelation.id2, totalWeight);
        out.write(&rel);
        ++written;
    }
    out.flush();

    for (size_t i = 0; i < readers.size(); ++i) {
        delete readers[i];
    }
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
                               size_t bufElems, size_t maxFanIn,
                               std::vector<uint64_t> & histOut, size_t & mergedOut,
                               size_t & ceilingOut) {
    histOut.assign(65536, 0);
    mergedOut = 0;
    ceilingOut = 0;

    std::vector<std::string> files;
    files.reserve(this->numOfGraph);
    for (size_t i = 0; i < this->numOfGraph; ++i) {
        files.push_back(outDir + "/subGraph_" + to_string(i) + "_" + to_string(route)
                        + "_" + to_string(shard));
    }
    reduceSubGraphFanIn(files, maxFanIn, route, shard);

    const size_t streamCnt = files.size();
    std::vector<ReadBuffer<Relation> *> readers(streamCnt);
    for (size_t i = 0; i < streamCnt; ++i) {
        readers[i] = new ReadBuffer<Relation>(files[i], bufElems);
    }

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
        if (totalWeight > 0) histOut[totalWeight]++;

        Relation rel(minRelation.id1, minRelation.id2, totalWeight);
        out.write(&rel);
    }
    out.flush();

    for (size_t i = 0; i < streamCnt; ++i) {
        delete readers[i];
    }
    for (size_t i = 0; i < streamCnt; ++i) {
        std::remove(files[i].c_str());
    }
}

int GroupGenerator::otsuThreshold(const std::vector<uint64_t>& hist) {
    uint64_t N = 0;
    double S = 0.0;
    for (size_t i = 0; i < hist.size(); i++) {
        N += hist[i];
        S += static_cast<double>(i) * static_cast<double>(hist[i]);
    }
    if (N == 0) return 0;

    double mean = S / static_cast<double>(N);
    double varTotal = 0.0;
    for (size_t i = 0; i < hist.size(); i++) {
        double d = static_cast<double>(i) - mean;
        varTotal += d * d * static_cast<double>(hist[i]);
    }
    varTotal /= static_cast<double>(N);
    if (varTotal == 0.0) return 0;

    double maxVarB = 0.0;
    int threshold = 0;
    uint64_t wB = 0;
    double sumB = 0.0;
    for (size_t t = 0; t < hist.size(); t++) {
        wB += hist[t];
        sumB += static_cast<double>(t) * static_cast<double>(hist[t]);
        uint64_t wF = N - wB;
        if (wB == 0 || wF == 0) continue;
        double muB = sumB / static_cast<double>(wB);
        double muF = (S - sumB) / static_cast<double>(wF);
        double varB = (static_cast<double>(wB) * static_cast<double>(wF)
                       / (static_cast<double>(N) * static_cast<double>(N)))
                      * (muB - muF) * (muB - muF);
        if (varB > maxVarB) {
            maxVarB = varB;
            threshold = static_cast<int>(t);
        }
    }

    return (maxVarB / varTotal >= 0.1) ? threshold : 0;
}

// Single-knee (Kneedle) threshold on the edge-weight CCDF.
// hist[w] = number of edges with total weight w (index 0 unused, size 65536).
// Unlike Otsu (which needs two comparable-mass classes), the knee only needs a
// curvature transition, so it survives the unimodal/power-law distributions of
// natural metagenomes. Returns 0 (= "no usable knee, fall back to --core-edge")
// on degenerate input, matching otsuThreshold()'s contract.
int GroupGenerator::kneeThreshold(const std::vector<uint64_t>& hist, int minWeight) {
    const int wLo = minWeight + 1;                 // Phase 1 only cuts above minWeight
 
    // Scan the usable domain: total edges, highest non-empty weight, distinct bins.
    uint64_t N = 0;
    int wmax = 0;
    int nonzero = 0;
    for (int w = wLo; w < static_cast<int>(hist.size()); w++) {
        if (hist[w] > 0) {
            N += hist[w];
            wmax = w;
            nonzero++;
        }
    }
    if (N == 0) return 0;                           // no edges above minWeight
    if (nonzero < 3 || wmax <= minWeight) return 0; // too few points to form a curve
 
    // --- robust wHi: anchor the chord's top end where enough edges actually
    // support that weight, instead of at the single largest (possibly
    // near-empty outlier) weight bin. minSamples floors at 100 edges or
    // 0.05% of N, whichever is larger, so it scales with dataset size.
    const uint64_t minSamples = std::max<uint64_t>(100, static_cast<uint64_t>(N * 0.0005));
    uint64_t survFromTop = 0;
    int wHi = wLo;
    for (int w = wmax; w >= wLo; w--) {
        survFromTop += hist[w];
        if (survFromTop >= minSamples) {
            wHi = w;
            break;
        }
    }
    if (wHi <= wLo) return 0;                       // not enough support anywhere above wLo
 
    // survFromTop was accumulated top-down from wmax to wHi, so it already
    // equals surv(wHi) = #edges with weight >= wHi.
    const uint64_t survHi = survFromTop;
    // --- end robust wHi ---
 
    const double xden = static_cast<double>(wHi - wLo);
    const double yden = static_cast<double>(N - survHi);
    if (xden <= 0.0 || yden <= 0.0) return 0;       // degenerate curve
 
    // CCDF surv(w) = #edges with weight >= w, monotone decreasing. After normalizing
    // x,y to [0,1] (endpoints (0,1) and (1,0)), the convex curve lies below the chord
    // y = 1 - x; the knee is the point of maximum distance below that chord.
    double bestDist = -1.0;
    int knee = 0;
    uint64_t cumBelow = 0;                          // edges with weight in [wLo, w)
    for (int w = wLo; w <= wHi; w++) {
        const uint64_t surv = N - cumBelow;
        const double xn = static_cast<double>(w - wLo) / xden;
        const double yn = (static_cast<double>(surv) - static_cast<double>(survHi)) / yden;
        const double dist = (1.0 - xn) - yn;
        if (dist > bestDist) {
            bestDist = dist;
            knee = w;
        }
        cumBelow += hist[w];
    }
 
    return (knee > minWeight) ? knee : 0;
}

void GroupGenerator::mergeGraph(std::vector<uint64_t>& edgeWeightHist) {
    cout << "Merging subgraphs" << endl;
    time_t before = time(nullptr);
    edgeWeightHist.assign(65536, 0);

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
    const size_t mergeBufElems = getMergeBufferElems(std::min(this->numOfGraph, maxFanIn),
                                                     par.ramUsage, concurrentMergers);
    const size_t streamsPerUnit = std::min(this->numOfGraph, maxFanIn);

    // A shard writes to a temp; the cross bucket's shards are concatenated afterwards.
    // Concatenation is legal because the relations_* consumers scan the file end to end and
    // never depend on ordering.
    auto unitOutPath = [&](size_t route, size_t shard) {
        return (shardsForRoute(route, par.threads) == 1)
             ? (outDir + "/relations_" + to_string(route) + ".bin")
             : (outDir + "/relations_" + to_string(route) + "_s" + to_string(shard) + ".tmp");
    };

    cout << "Merge: " << unitCnt << " units (" << routeCnt << " routes, cross bucket sharded x"
         << shardsForRoute(static_cast<size_t>(par.threads), par.threads) << ") x "
         << this->numOfGraph << " subgraphs, " << concurrentMergers
         << " concurrent (fd soft limit " << getOpenFileLimit()
         << ", fan-in " << maxFanIn << ")" << endl;
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

    std::vector<std::vector<uint64_t>> unitHist(unitCnt);
    std::vector<size_t> unitMerged(unitCnt, 0);
    std::vector<size_t> unitCeiling(unitCnt, 0);
    std::vector<double> unitSeconds(unitCnt, 0.0);

    // schedule(dynamic): units still carry different loads even after sharding.
    #pragma omp parallel for schedule(dynamic) num_threads(concurrentMergers)
    for (size_t u = 0; u < unitCnt; ++u) {
        const time_t unitStart = time(nullptr);
        mergeUnit(units[u].first, units[u].second, unitOutPath(units[u].first, units[u].second),
                  mergeBufElems, maxFanIn, unitHist[u], unitMerged[u], unitCeiling[u]);
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
        for (size_t w = 0; w < edgeWeightHist.size(); ++w) {
            edgeWeightHist[w] += unitHist[u][w];
        }
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
    const size_t mergeBufElems = getMergeBufferElems(std::min(this->numOfGraph, maxFanIn),
                                                     par.ramUsage, 1);
    cout << "Merge: " << units.size() << " units x " << this->numOfGraph
         << " subgraphs, sequential (fd soft limit " << getOpenFileLimit()
         << ", fan-in " << maxFanIn << ")" << endl;

    WriteBuffer<Relation> relationLog(outDir + "/relations.bin", 1024 * 1024);

    size_t ceilingEdgeCnt = 0; // edges whose merged weight sits at the uint16 ceiling

    for (size_t u = 0; u < units.size(); ++u) {
        const size_t route = units[u].first;
        const size_t shard = units[u].second;
        std::vector<std::string> subGraphFiles;
        subGraphFiles.reserve(this->numOfGraph);
        for (size_t i = 0; i < this->numOfGraph; ++i) {
            subGraphFiles.push_back(outDir + "/subGraph_" + to_string(i) + "_"
                                    + to_string(route) + "_" + to_string(shard));
        }
        reduceSubGraphFanIn(subGraphFiles, maxFanIn, route, shard);

        const size_t streamCnt = subGraphFiles.size();
        std::vector<ReadBuffer<Relation> *> relationBuffers(streamCnt);
        std::vector<Relation> currentRelations(streamCnt);
        for (size_t i = 0; i < streamCnt; ++i) {
            relationBuffers[i] = new ReadBuffer<Relation>(subGraphFiles[i], mergeBufElems);
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

            Relation rel(minRelation.id1, minRelation.id2, totalWeight);
            relationLog.write(&rel);
        }

        for (size_t i = 0; i < streamCnt; ++i) {
            delete relationBuffers[i];
        }
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

void GroupGenerator::computeNodeDegree(
    int threshold,
    size_t processedReadCnt,
    std::vector<uint32_t>& degree)
{
    degree.assign(processedReadCnt + 1, 0);

    auto processFile = [&](const std::string& fname) {
        ReadBuffer<Relation> rb(fname, 1024 * 1024);
        for (Relation r = rb.getNext(); !(r == Relation()); r = rb.getNext()) {
            uint32_t id1 = r.id1, id2 = r.id2;
            uint16_t w = r.weight;
            if (id1 == 0 || id2 == 0) continue;
            if (id1 > processedReadCnt || id2 > processedReadCnt) continue;
            if (w > threshold) {
                degree[id1]++;
                degree[id2]++;
            }
        }
    };

    for (int i = 0; i < par.threads * 2 + 1; i++) {
        processFile(outDir + "/relations_" + std::to_string(i) + ".bin");
    }
}

void GroupGenerator::computeGroupQuarterDegree(
    const std::vector<uint32_t>& queryGroupInfo,
    const std::vector<uint32_t>& degree,
    std::unordered_map<uint32_t, uint32_t>& groupQuarterDeg)
{
    std::unordered_map<uint32_t, std::vector<uint32_t>> groupDegrees;

    for (uint32_t i = 1; i < queryGroupInfo.size(); i++) {
        uint32_t groupId = queryGroupInfo[i];
        if (groupId == 0) continue; // skip ungrouped node
        groupDegrees[groupId].push_back(degree[i]);
    }

    groupQuarterDeg.clear();
    for (auto& [groupId, degrees] : groupDegrees) {
        size_t n = degrees.size();
        std::nth_element(degrees.begin(), degrees.begin() + n / 4, degrees.end());
        uint32_t p25 = degrees[n / 4];
        groupQuarterDeg[groupId] = p25;
    }
}

void GroupGenerator::makeGroupsAdaptive(
    const std::vector<uint16_t>& nodeThr,
    size_t processedReadCnt,
    std::vector<uint32_t>& queryGroupInfo
) {
    cout << "Creating groups (adaptive thresholds)..." << endl;
    time_t beforeSearch = time(nullptr);

    DisjointSet ds(processedReadCnt);

    auto processFile = [&](const std::string& fname, DisjointSet& subDs) {
        ReadBuffer<Relation> rb(fname, 1024 * 1024);
        for (Relation r = rb.getNext(); !(r == Relation()); r = rb.getNext()) {
            uint32_t id1 = r.id1, id2 = r.id2;
            uint16_t w = r.weight;
            if (id1 == 0 || id2 == 0) continue;
            if (id1 > processedReadCnt || id2 > processedReadCnt) continue;

            if (keepEdgeGeo(w, nodeThr[id1], nodeThr[id2])) {
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
            if (keepEdgeGeo(r.weight, nodeThr[r.id1], nodeThr[r.id2])) {
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
        if (ds.grouped[queryId]) {
            queryGroupInfo[queryId] = rootToMin[ds.find(queryId)];
        }
    }

    cout << "Adaptive grouping done." << endl;
    cout << "Time spent: " << double(time(nullptr) - beforeSearch) << " seconds." << endl;
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