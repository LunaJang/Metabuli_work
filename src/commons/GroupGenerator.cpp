#include "GroupGenerator.h"
#include "FileUtil.h"
#include "QueryIndexer.h"
#include "common.h"
#include "Kmer.h"

GroupGenerator::GroupGenerator(LocalParameters & par) {
    // Load parameters
    dbDir = par.filenames[1 + (par.seqMode == 2)];
    matchPerKmer = par.matchPerKmer;
    loadDbParameters(par, par.filenames[1 + (par.seqMode == 2)]);
    
    cout << "DB name: " << par.dbName << endl;
    cout << "DB creation date: " << par.dbDate << endl;
    
    // Taxonomy
    taxonomy = loadTaxonomy(dbDir, par.taxonomyPath);
    
    // Agents
    queryIndexer = new QueryIndexer(par);
    kmerExtractor = new KmerExtractor(par);
    if (par.reducedAA) {
        kmerMatcher = new ReducedKmerMatcher(par, taxonomy);
    } else {
        kmerMatcher = new KmerMatcher(par, taxonomy);
    }
    taxonomer = new Taxonomer(par, taxonomy);
    reporter = new Reporter(par, taxonomy);
    seqIterator = new SeqIterator(par);
    kmerFileHandler = new KmerFileHandler();
}

GroupGenerator::~GroupGenerator() {
    delete taxonomy;
    delete queryIndexer;
    delete kmerExtractor;
    delete kmerMatcher;
    delete taxonomer;
    delete reporter;
    delete seqIterator;
    delete kmerFileHandler;
}

void GroupGenerator::startGroupGeneration(const LocalParameters &par) {  

    QueryKmerBuffer queryKmerBuffer;
    vector<Query> queryList;
    size_t numOfTatalQueryKmerCnt = 0;
    size_t numOfThreads = par.threads;

    bool complete = false;
    size_t processedReadCnt = 0;
    size_t tries = 0;
    size_t totalSeqCnt = 0;
    size_t numOfSplits = 0;
    size_t numOfGraph = 0;

    string outDir;
    string jobId;
    if (par.seqMode == 2) {
        outDir = par.filenames[3];
        jobId = par.filenames[4];
    } else {
        outDir = par.filenames[2];
        jobId = par.filenames[3];
    }

    size_t voteMode = par.voteMode;
    float majorityThr = par.majorityThr;
    float groupScoreThr = par.groupScoreThr;
    size_t groupKmerThr = par.groupKmerThr;
    cout << "voteMode: " << voteMode << endl;
    cout << "majorityThr: " << majorityThr << endl;
    cout << "groupScoreThr: " << groupScoreThr << endl;
    cout << "groupKmerThr: " << groupKmerThr << endl;
    
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

            // Initialize query k-mer buffer and match buffer
            queryKmerBuffer.startIndexOfReserve = 0;

            // Extract query k-mers
            kmerExtractor->extractQueryKmers(queryKmerBuffer,
                                             queryList,
                                             queryReadSplit[splitIdx],
                                             par,
                                             kseq1,
                                             kseq2); // sync kseq1 and kseq2            
            // saveQueryIdToFile
            kmerFileHandler->writeQueryKmerFile(&queryKmerBuffer, outDir, numOfSplits, numOfThreads, processedReadCnt, jobId);
            processedReadCnt += queryReadSplit[splitIdx].readCnt;
            cout << "The number of processed sequences: " << processedReadCnt << " (" << (double) processedReadCnt / (double) totalSeqCnt << ")" << endl;

            numOfTatalQueryKmerCnt += queryKmerBuffer.startIndexOfReserve;
        }

        delete kseq1;
        if (par.seqMode == 2) {
            delete kseq2;
        }
        if (processedReadCnt == totalSeqCnt) {
            complete = true;
        }
    }   

    makeGraph(outDir, numOfSplits, numOfThreads, numOfGraph, jobId);   
    unordered_map<uint32_t, unordered_set<uint32_t>> groupInfo;
    vector<int> queryGroupInfo;
    queryGroupInfo.resize(processedReadCnt, -1);
    makeGroups(groupInfo, outDir, queryGroupInfo, groupKmerThr, numOfGraph, jobId);
    
    saveGroupsToFile(groupInfo, queryGroupInfo, outDir, jobId);
    
    vector<pair<int, float>> metabuliResult;       
    metabuliResult.resize(processedReadCnt, make_pair(-1, 0.0f));
    loadMetabuliResult(outDir, metabuliResult);

    unordered_map<uint32_t, int> repLabel; 
    getRepLabel(outDir, metabuliResult, groupInfo, repLabel, jobId, voteMode, majorityThr, groupScoreThr);

    applyRepLabel(outDir, outDir, queryGroupInfo, repLabel, groupScoreThr, jobId);
    
    cout << "Number of query k-mers: " << numOfTatalQueryKmerCnt << endl;
}

void GroupGenerator::makeGroupsFromBinning(const string &binningFileDir, 
                                           unordered_map<uint32_t, unordered_map<uint32_t, uint32_t>> &relation,
                                           unordered_map<uint32_t, unordered_set<uint32_t>> &groupInfo, 
                                           vector<int> &queryGroupInfo, 
                                           int groupKmerThr) {
    // Map to store the count of shared Y nodes between X nodes
    DisjointSet ds;
    string line;

    cout << "Creating groups based on graph..." << endl;
    time_t beforeSearch = time(nullptr);
    
    for (const auto& currentQueryRelation : relation) {
        uint32_t currentQueryId = currentQueryRelation.first;
        for (const auto& [otherQueryId, count] : currentQueryRelation.second) {
            if (count >= groupKmerThr) {
                if (ds.parent.find(currentQueryId) == ds.parent.end()) {
                    ds.makeSet(currentQueryId);
                }
                if (ds.parent.find(otherQueryId) == ds.parent.end()) {
                    ds.makeSet(otherQueryId);
                }
                ds.unionSets(currentQueryId, otherQueryId);
            }
        }
    }

    // Collect nodes into groups
    for (const auto& p : ds.parent) {
        uint32_t currentQueryId = p.first;
        uint32_t groupId = ds.find(currentQueryId);
        groupInfo[groupId].insert(currentQueryId);
        queryGroupInfo[currentQueryId] = groupId;
    }

    cout << "Query group created successfully : " << groupInfo.size() << " groups" << endl;
    cout << "Time spent for query groups: " << double(time(nullptr) - beforeSearch) << endl;

    return;
}

void GroupGenerator::makeGraph(const string &queryKmerFileDir, 
                               size_t &numOfSplits, 
                               size_t &numOfThreads, 
                               size_t &numOfGraph,
                               const string &jobId) {
    
    cout << "Creating graphs based on kmer-query relation..." << endl;
    time_t beforeSearch = time(nullptr);

    const size_t RELATION_THRESHOLD = 100'000'000;  // relation 크기 제한
    std::atomic<int> counter(0);

    #pragma omp parallel num_threads(numOfThreads)
    {
        int threadIdx = omp_get_thread_num();
        size_t processedRelationCnt = 0;
        map<uint32_t, map<uint32_t, uint32_t>> threadRelation;

        // 스레드별 파일 정보 할당
        struct MmapedData<uint16_t> *diffFileList = new struct MmapedData<uint16_t>[numOfSplits];
        struct MmapedData<QueryKmerInfo> *infoFileList = new struct MmapedData<QueryKmerInfo>[numOfSplits];

        uint64_t *lookingKmers = new uint64_t[numOfSplits];
        auto *lookingInfos = new QueryKmerInfo[numOfSplits];
        auto *diffFileIdx = new size_t[numOfSplits];
        memset(diffFileIdx, 0, numOfSplits * sizeof(size_t));
        auto *infoFileIdx = new size_t[numOfSplits];
        memset(infoFileIdx, 0, numOfSplits * sizeof(size_t));
        auto *maxIdxOfEachFiles = new size_t[numOfSplits];

        // 파일 초기화 및 mmap 처리
        for (size_t file = 0; file < numOfSplits; ++file) {
            string diffIdxFilename = queryKmerFileDir + "/" + jobId + "_queryKmerRelation" + to_string(file) + "_diffIdx" + to_string(threadIdx);
            string infoFilename = queryKmerFileDir + "/" + jobId + "_queryKmerRelation" + to_string(file) + "_info" + to_string(threadIdx);

            diffFileList[file] = mmapData<uint16_t>(diffIdxFilename.c_str());
            infoFileList[file] = mmapData<QueryKmerInfo>(infoFilename.c_str());

            if (!diffFileList[file].data || !infoFileList[file].data || infoFileList[file].fileSize == 0) {
                lookingKmers[file] = UINT64_MAX; // 해당 파일을 무시하도록 설정
                continue;
            }

            lookingKmers[file] = kmerFileHandler->getNextKmer(0, diffFileList[file], diffFileIdx[file]);
            lookingInfos[file] = infoFileList[file].data[0];
            infoFileIdx[file]++;
            maxIdxOfEachFiles[file] = diffFileList[file].fileSize / sizeof(uint16_t);
        }

        vector<uint32_t> currentQueryIds;
        currentQueryIds.reserve(1024);
        while (true) {
            uint64_t currentKmer = UINT64_MAX;

            // 가장 작은 k-mer 찾기 (스레드별 파일에서 찾음)
            for (size_t file = 0; file < numOfSplits; ++file) {
                if (lookingKmers[file] < currentKmer) {
                    currentKmer = lookingKmers[file];
                }
            }

            if (currentKmer == UINT64_MAX) {
                break;
            }

            currentQueryIds.clear();

            // 현재 k-mer에 해당하는 query ID들을 수집
            for (size_t file = 0; file < numOfSplits; ++file) {
                while (lookingKmers[file] == currentKmer) {
                    currentQueryIds.push_back(lookingInfos[file].sequenceID);

                    lookingKmers[file] = kmerFileHandler->getNextKmer(lookingKmers[file], diffFileList[file], diffFileIdx[file]);

                    if (!infoFileList[file].data || infoFileIdx[file] >= infoFileList[file].fileSize / sizeof(QueryKmerInfo)) {
                        lookingKmers[file] = UINT64_MAX; // 파일 끝 도달
                        break;
                    }
                    lookingInfos[file] = infoFileList[file].data[infoFileIdx[file]];
                    infoFileIdx[file]++;
                }
            }
            // 관계 저장 (스레드별)
            for (size_t i = 0; i < currentQueryIds.size(); ++i) {
                for (size_t j = i; j < currentQueryIds.size(); ++j) {
                    if (currentQueryIds[i] != UINT32_MAX && currentQueryIds[j] != UINT32_MAX){
                        if (currentQueryIds[i] < currentQueryIds[j]) {
                            threadRelation[currentQueryIds[i]][currentQueryIds[j]]++;
                            processedRelationCnt++;
                        } else if (currentQueryIds[i] > currentQueryIds[j]) {
                            threadRelation[currentQueryIds[j]][currentQueryIds[i]]++;
                            processedRelationCnt++;
                        }
                    }
                }
            }

            // relation이 threshold를 넘으면 저장 후 초기화
            if (processedRelationCnt > RELATION_THRESHOLD) {                
                size_t counter_now = counter.fetch_add(1, std::memory_order_relaxed);
                saveSubGraphToFile(threadRelation, queryKmerFileDir, counter_now, jobId);
                threadRelation.clear();
                processedRelationCnt = 0;
            }
        }

        // 남아 있는 관계 저장
        if (!threadRelation.empty()) {
            size_t counter_now = counter.fetch_add(1, std::memory_order_relaxed);
            saveSubGraphToFile(threadRelation, queryKmerFileDir, counter_now, jobId);
        }


        // 파일 닫기
        for (size_t file = 0; file < numOfSplits; file++) {
            if (diffFileList[file].data) {
                munmap(diffFileList[file].data, diffFileList[file].fileSize);
            }
            if (infoFileList[file].data) {
                munmap(infoFileList[file].data, infoFileList[file].fileSize);
            }
        }

        // 동적 메모리 해제
        delete[] diffFileList;
        delete[] infoFileList;
        delete[] lookingInfos;
        delete[] lookingKmers;
        delete[] diffFileIdx;
        delete[] infoFileIdx;
        delete[] maxIdxOfEachFiles;
    }

    
    numOfGraph = counter.fetch_add(1, std::memory_order_relaxed);

    cout << "Relations generated from files successfully." << endl;
    cout << "Time spent: " << double(time(nullptr) - beforeSearch) << " seconds." << endl;
}



void GroupGenerator::saveSubGraphToFile(const map<uint32_t, map<uint32_t, uint32_t>> &subRelation, 
                                        const string &subGraphFileDir, 
                                        const size_t counter_now,
                                        const string &jobId){
    const string& subGraphFileName = subGraphFileDir + "/" + jobId + "_subGraph" + to_string(counter_now);

    ofstream outFile(subGraphFileName);
    if (!outFile.is_open()) {
        cerr << "Error opening file: " << subGraphFileName << endl;
        return;
    }

    for (const auto &[queryid_1, inner_map] : subRelation) {
        for (const auto &[queryid_2, weight] : inner_map) {
            outFile << queryid_1 << "," << queryid_2 << "," << weight << "\n";
        }
    }
    
    outFile.close();
    #pragma omp critical
    {
        cout << "Query sub-graph saved to " << subGraphFileName << " successfully." << endl;
    }
    
}


void GroupGenerator::makeGroups(unordered_map<uint32_t, unordered_set<uint32_t>> &groupInfo, 
                                const string &subGraphFileDir, 
                                vector<int> &queryGroupInfo, 
                                int groupKmerThr, 
                                size_t &numOfGraph,
                                const string &jobId) {
    DisjointSet ds;
    time_t beforeSearch = time(nullptr);

    cout << "Creating groups based on graph..." << endl;

    vector<ifstream> files(numOfGraph);
    vector<tuple<uint32_t, uint32_t, uint32_t>> currentLines(numOfGraph);  // (id1, id2, weight)
    vector<bool> fileHasData(numOfGraph, true);

    // 1. open all files and read the first line
    for (size_t i = 0; i < numOfGraph; ++i) {
        string fileName = subGraphFileDir + "/" + jobId + "_subGraph" + to_string(i);
        files[i].open(fileName);
        if (!files[i].is_open()) {
            cerr << "Error opening file: " << fileName << endl;
            fileHasData[i] = false;
            continue;
        }
        string line;
        if (getline(files[i], line)) {
            istringstream ss(line);
            string token;
            getline(ss, token, ','); uint32_t id1 = stoi(token);
            getline(ss, token, ','); uint32_t id2 = stoi(token);
            getline(ss, token, ','); uint32_t weight = stoi(token);
            currentLines[i] = {id1, id2, weight};
        } else {
            fileHasData[i] = false;
        }
    }

    // 2. process until all files are done
    while (true) {
        // 가장 작은 (id1, id2) 찾기
        pair<uint32_t, uint32_t> minKey = {UINT32_MAX, UINT32_MAX};
        for (size_t i = 0; i < numOfGraph; ++i) {
            if (fileHasData[i]) {
                auto [id1, id2, _] = currentLines[i];
                if (make_pair(id1, id2) < minKey) {
                    minKey = {id1, id2};
                }
            }
        }
        if (minKey == make_pair(UINT32_MAX, UINT32_MAX)) break;

        // 동일한 relation에 대해 weight 누적
        uint32_t totalWeight = 0;
        for (size_t i = 0; i < numOfGraph; ++i) {
            if (fileHasData[i]) {
                auto [id1, id2, weight] = currentLines[i];
                if (make_pair(id1, id2) == minKey) {
                    totalWeight += weight;
                    string line;
                    if (getline(files[i], line)) {
                        istringstream ss(line);
                        string token;
                        getline(ss, token, ','); id1 = stoi(token);
                        getline(ss, token, ','); id2 = stoi(token);
                        getline(ss, token, ','); weight = stoi(token);
                        currentLines[i] = {id1, id2, weight};
                    } else {
                        fileHasData[i] = false;
                    }
                }
            }
        }

        // 조건에 따라 union
        if (totalWeight >= groupKmerThr) {
            uint32_t id1 = get<0>(minKey), id2 = get<1>(minKey);
            if (ds.parent.find(id1) == ds.parent.end()) ds.makeSet(id1);
            if (ds.parent.find(id2) == ds.parent.end()) ds.makeSet(id2);
            ds.unionSets(id1, id2);
        }
    }

    // 그룹 수집
    for (const auto& [queryId, _] : ds.parent) {
        uint32_t groupId = ds.find(queryId);
        groupInfo[groupId].insert(queryId);
        if (queryId >= queryGroupInfo.size()) queryGroupInfo.resize(queryId + 1, -1);
        queryGroupInfo[queryId] = groupId;
    }

    cout << "Query group created successfully : " << groupInfo.size() << " groups" << endl;
    cout << "Time spent: " << double(time(nullptr) - beforeSearch) << " seconds." << endl;
}


void GroupGenerator::saveGroupsToFile(const unordered_map<uint32_t, unordered_set<uint32_t>> &groupInfo, 
                                      const vector<int> &queryGroupInfo, 
                                      const string &groupFileDir, 
                                      const string &jobId) {
    // save group in txt file
    const string& groupInfoFileName = groupFileDir + "/" + jobId + "_groups";
    ofstream outFile1(groupInfoFileName);
    if (!outFile1.is_open()) {
        cerr << "Error opening file: " << groupInfoFileName << endl;
        return;
    }

    for (const auto& [groupId, queryIds] : groupInfo) {
        outFile1 << groupId << " ";
        for (const auto& queryId : queryIds) {
            outFile1 << queryId << " ";
        }
        outFile1 << endl;
    }
    outFile1.close();
    cout << "Query group saved to " << groupInfoFileName << " successfully." << endl;
    

    const string& queryGroupInfoFileName = groupFileDir + "/" + jobId + "_queryGroupMap";
    ofstream outFile2(queryGroupInfoFileName);
    if (!outFile2.is_open()) {
        cerr << "Error opening file: " << queryGroupInfoFileName << endl;
        return;
    }

    for (size_t i = 0; i < queryGroupInfo.size(); ++i) {
        outFile2 << queryGroupInfo[i] << "\n";
    }
    outFile2.close();
    cout << "Query group saved to " << queryGroupInfoFileName << " successfully." << endl;

    return;
}

void GroupGenerator::loadMetabuliResult(const string &resultFileDir, 
                                        vector<pair<int, float>> &metabuliResult) {
    ifstream inFile(resultFileDir + "/1_classifications.tsv");
    if (!inFile.is_open()) {
        cerr << "Error opening file: " << resultFileDir + "/1_classifications.tsv" << endl;
        return;
    }

    string line;
    uint32_t id = 0;
    while (getline(inFile, line)) {
        stringstream ss(line);
        int label, query_label, read_length;
        float score;
        string query_name; // query_name is not used.
        ss >> label >> query_name >> query_label >> read_length >> score;
        metabuliResult[id] = {query_label, score}; 
        id ++;
    }

    inFile.close();
    cout << "Original Metabuli result loaded from " << resultFileDir + "/1_classifications.tsv" << " successfully." << endl;
}

void GroupGenerator::getRepLabel(const string &groupRepFileDir, 
                                 vector<pair<int, float>> &metabuliResult, 
                                 const unordered_map<uint32_t, unordered_set<uint32_t>> &groupInfo, 
                                 unordered_map<uint32_t, int> &repLabel, 
                                 const string &jobId, 
                                 int voteMode, 
                                 float majorityThr,
                                 const float groupScoreThr) {
    for (const auto& group : groupInfo) {
        uint32_t groupId = group.first;
        const unordered_set<uint32_t>& queryIds = group.second;

        vector<WeightedTaxHit> setTaxa;

        for (const auto& queryId : queryIds) {
            int query_label = metabuliResult[queryId].first; 
            float score = metabuliResult[queryId].second;
            if (query_label != 0 && score >= groupScoreThr) {
                setTaxa.emplace_back(query_label, score, voteMode);
            }
        }

        WeightedTaxResult result = taxonomy->weightedMajorityLCA(setTaxa, majorityThr);

        if (result.taxon != 0) {
            repLabel[groupId] = result.taxon;
        }
    }

    const string& groupRepFileName = groupRepFileDir + "/" + jobId + "_groupRep";
    ofstream outFile(groupRepFileName);
    if (!outFile.is_open()) {
        cerr << "Error opening file: " << groupRepFileName << endl;
        return;
    }

    for (const auto& [groupId, groupRep] : repLabel) {
        outFile << groupId << "\t" << groupRep << "\n";
    }

    outFile.close();

    cout << "Query group representative label saved to " << groupRepFileName << " successfully." << endl;
}

void GroupGenerator::applyRepLabel(const string &resultFileDir, 
                                   const string &newResultFileDir, 
                                   const vector<int> &queryGroupInfo, 
                                   const unordered_map<uint32_t, int> &repLabel, 
                                   const float groupScoreThr,
                                   const string &jobId) {
    ifstream inFile(resultFileDir + "/" + jobId + "_classifications.tsv");
    if (!inFile.is_open()) {
        cerr << "Error opening file: " << newResultFileDir + "/" + jobId + "_classifications.tsv" << endl;
        return;
    }

    ofstream outFile(newResultFileDir + "/" + jobId + "_updated_classifications.tsv");
    if (!outFile.is_open()) {
        cerr << "Error opening file: " << newResultFileDir + "/" + jobId + "_updated_classifications.tsv" << endl;
        return;
    }

    string line;
    uint32_t queryIdx = 0;
    while (getline(inFile, line)) {
        stringstream ss(line);
        vector<string> fields;
        string field;

        while (getline(ss, field, '\t')) {
            fields.push_back(field);
        }

        if ((fields.size() > 2) && ((fields[0] == "0") || (std::stof(fields[4]) < groupScoreThr))) {
            uint32_t groupId = queryGroupInfo[queryIdx];
            auto repLabelIt = repLabel.find(groupId);
            if (repLabelIt != repLabel.end() && repLabelIt->second != 0) {
                fields[2] = to_string(repLabelIt->second);
                fields[0] = "1";
                // fields[5] = taxonomy->getString(taxonomy->taxonNode(repLabelIt->second)->rankIdx);
                fields[5] = "-";
            }
        }

        for (size_t i = 0; i < fields.size(); ++i) {
            outFile << fields[i]; 
            if (i < fields.size() - 1) {  
                outFile << "\t";
            }
        }
        outFile << endl;
        queryIdx++;
    }

    inFile.close();
    outFile.close();
    
    cout << "Result saved to " << newResultFileDir << " successfully." << endl;
}

// KmerFileHandler function implementations
void KmerFileHandler::flushKmerBuf(uint16_t *buffer, FILE *handleKmerTable, size_t &localBufIdx) {
    fwrite(buffer, sizeof(uint16_t), localBufIdx, handleKmerTable);
    localBufIdx = 0; // Reset buffer index
}

void KmerFileHandler::getDiffIdx(const uint64_t &lastKmer, const uint64_t &entryToWrite, FILE *handleKmerTable, uint16_t *buffer, size_t &localBufIdx, size_t bufferSize) {
    uint64_t kmerdiff = entryToWrite - lastKmer;
    uint16_t tempBuffer[5];
    int idx = 3;
    tempBuffer[4] = SET_END_FLAG(GET_15_BITS(kmerdiff));
    kmerdiff >>= 15U;
    while (kmerdiff) {
        uint16_t toWrite = GET_15_BITS(kmerdiff);
        kmerdiff >>= 15U;
        tempBuffer[idx] = toWrite;
        idx--;
    }
    // Buffer the result
    writeDiffIdx(buffer, handleKmerTable, tempBuffer + idx + 1, 4 - idx, localBufIdx, bufferSize);
}

void KmerFileHandler::writeDiffIdx(uint16_t *buffer, FILE *handleKmerTable, uint16_t *toWrite, size_t size, size_t &localBufIdx, size_t bufferSize) {
    if (localBufIdx + size >= bufferSize) {
        flushKmerBuf(buffer, handleKmerTable, localBufIdx); // Flush if buffer is full
    }
    memcpy(buffer + localBufIdx, toWrite, sizeof(uint16_t) * size); // Add new data to buffer
    localBufIdx += size; // Update buffer index
}

uint64_t KmerFileHandler::getNextKmer(uint64_t lookingTarget, const struct MmapedData<uint16_t> & diffList, size_t & idx) {
    uint16_t fragment = 0;
    uint64_t diffIn64bit = 0;
    uint16_t check = 32768; // 2^15

    if (idx >= diffList.fileSize / sizeof(uint16_t)) {
        return UINT64_MAX;  
    }

    fragment = diffList.data[idx++];
    while(!(fragment & check)) {
        diffIn64bit |= fragment;
        diffIn64bit <<= 15u;
        fragment = diffList.data[idx++];
    }

    fragment &= ~check;
    diffIn64bit |= fragment;

    return diffIn64bit + lookingTarget;
}

void KmerFileHandler::writeQueryKmerFile(QueryKmerBuffer * queryKmerBuffer, 
                                         const string& queryKmerFileDir, 
                                         size_t& numOfSplits, 
                                         size_t numOfThreads,
                                         size_t processedReadCnt, 
                                         const string & jobId) {

    size_t queryKmerNum = queryKmerBuffer->startIndexOfReserve;
    QueryKmer *queryKmerList = queryKmerBuffer->buffer;
    for (size_t checkN = queryKmerNum - 1; checkN > 0; checkN--) {
        if (queryKmerList[checkN].ADkmer != UINT64_MAX) {
            queryKmerNum = checkN + 1;
            break;
        }
    }
    
    vector<FILE*> diffIdxFiles(numOfThreads, nullptr);
    vector<FILE*> infoFiles(numOfThreads, nullptr);
    vector<uint16_t*> diffIdxBuffers(numOfThreads, nullptr);
    vector<size_t> localBufIdxs(numOfThreads, 0);

    size_t bufferSize = max(10'000'000 / numOfThreads, size_t(100'000)); 

    for (size_t i = 0; i < numOfThreads; i++) {
        string diffIdxFileName = queryKmerFileDir + "/" + jobId + "_queryKmerRelation" + to_string(numOfSplits) + "_diffIdx" + to_string(i);
        string infoFileName = queryKmerFileDir + "/" + jobId + "_queryKmerRelation" + to_string(numOfSplits) + "_info" + to_string(i);

        diffIdxFiles[i] = fopen(diffIdxFileName.c_str(), "wb");
        infoFiles[i] = fopen(infoFileName.c_str(), "wb");
        if (diffIdxFiles[i] == nullptr || infoFiles[i] == nullptr) {
            cout << "Cannot open file for writing: " << diffIdxFileName << " or " << infoFileName << endl;
            return;
        }
        
        diffIdxBuffers[i] = (uint16_t *)malloc(sizeof(uint16_t) * bufferSize);
        if (diffIdxBuffers[i] == nullptr ) {
            cout << "Memory allocation failed for diffIdxBuffer[" << i << "]" << endl;
            flushKmerBuf(diffIdxBuffers[i], diffIdxFiles[i], localBufIdxs[i]);
            free(diffIdxBuffers[i]); 
            return;
        }
    }

    uint64_t lastKmer[numOfThreads] = {0, };
    size_t write = 0;
    const uint64_t MAX_KMER = UINT64_MAX;
    size_t rangeSize = MAX_KMER / numOfThreads;
    
    for(size_t i = 0; i < queryKmerNum ; i++) {
        queryKmerList[i].info.sequenceID += processedReadCnt;
        size_t splitIdx = min(queryKmerList[i].ADkmer / rangeSize, numOfThreads - 1);   

        fwrite(&queryKmerList[i].info, sizeof (QueryKmerInfo), 1, infoFiles[splitIdx]); 
        getDiffIdx(lastKmer[splitIdx], queryKmerList[i].ADkmer, diffIdxFiles[splitIdx], diffIdxBuffers[splitIdx], localBufIdxs[splitIdx], bufferSize); 
        lastKmer[splitIdx] = queryKmerList[i].ADkmer;   
        write++; 
    }

    cout<<"total k-mer count  : "<< queryKmerNum << endl;
    cout<<"written k-mer count: "<< write <<endl;

    for (size_t i = 0; i < numOfThreads; i++) {
        flushKmerBuf(diffIdxBuffers[i], diffIdxFiles[i], localBufIdxs[i]);
        free(diffIdxBuffers[i]);     
        diffIdxBuffers[i] = nullptr;
        fclose(diffIdxFiles[i]);
        fclose(infoFiles[i]);
    }

    numOfSplits++;

    return;
}
