#include "GroupApplier.h"
#include "FileUtil.h"
#include "common.h"
#include "Kmer.h"
GroupApplier::GroupApplier(LocalParameters & par) : par(par) {
    groupFileDir    = par.filenames[1 + (par.seqMode == 2)];
    groupmapFileDir = par.filenames[2 + (par.seqMode == 2)];
    taxDbDir        = par.filenames[3 + (par.seqMode == 2)];
    orgRes          = par.filenames[4 + (par.seqMode == 2)];
    outDir          = par.filenames[5 + (par.seqMode == 2)];
    
    taxonomy = new TaxonomyWrapper(
                    taxDbDir + "/names.dmp",
                    taxDbDir + "/nodes.dmp",
                    taxDbDir + "/merged.dmp",
                    true);
    
    updatedResultFileName = outDir + "/updated_classifications.tsv";
    updatedReportFileName = outDir + "/updated_report.tsv";

    reporter = new Reporter(par, taxonomy);
}

GroupApplier::~GroupApplier() {
    delete taxonomy;
    delete reporter;
}

void GroupApplier::startGroupApplication(const LocalParameters &par) {  
    Buffer<Kmer> queryKmerBuffer;
    Buffer<std::pair<uint32_t, uint32_t>> matchBuffer; // seq id, pos
    vector<Query> queryList;

    bool complete = false;
    size_t processedReadCnt = 0;
    size_t tries = 0;
    size_t totalSeqCnt = 0;

    vector<OrgResult> orgResult;       
    loadOrgResult(orgResult, processedReadCnt);

    unordered_map<uint32_t, unordered_set<uint32_t>> groupInfo;
    vector<uint32_t> groupMappingInfo;
    groupMappingInfo.resize(processedReadCnt, 0);
    loadGroupsFromFile(groupInfo, groupMappingInfo, outDir);

    std::unordered_map<int, int> external2internalTaxId;
    taxonomy->getExternal2internalTaxID(external2internalTaxId);
    unordered_map<uint32_t, uint32_t> repLabel; 
    getRepLabel(orgResult, groupInfo, repLabel, external2internalTaxId);
    applyRepLabel(orgResult, groupMappingInfo, repLabel, external2internalTaxId);
}


void GroupApplier::loadOrgResult(vector<OrgResult>& orgResults, size_t& processedReadCnt) {
    ifstream inFile(orgRes);
    if (!inFile.is_open()) {
        cerr << "Error opening file: " << orgRes << endl;
        return;
    }

    int classificationCol = par.taxidCol - 1; 
    int scoreCol = par.scoreCol - 1; 
    int readNameCol = par.readIdCol - 1; 
    if (par.weightMode == 0 || scoreCol < 0) {
        string line;
        while (getline(inFile, line)) {
            if (line.empty()) continue;
            if (line.front() == '#') continue;
            std::vector<std::string> columns = TaxonomyWrapper::splitByDelimiter(line, "\t", 20);
            TaxID taxId = stoi(columns[classificationCol]);
            orgResults.push_back({taxId, 1.0, columns[readNameCol]});
            processedReadCnt++;
        }
    } else {
        string line;
        while (getline(inFile, line)) {
            if (line.empty()) continue;
            if (line.front() == '#') continue;
            std::vector<std::string> columns = TaxonomyWrapper::splitByDelimiter(line, "\t", 20);
            TaxID taxId = stoi(columns[classificationCol]);
            float score = stof(columns[scoreCol]);
            orgResults.push_back({taxId, score, columns[readNameCol]});
            processedReadCnt++;
        }
    }
    inFile.close();
    cout << "Original Metabuli result loaded from " << orgRes << " successfully." << endl;
    cout << "Number of query result: " << orgResults.size() << endl;
}

void GroupApplier::loadGroupsFromFile(unordered_map<uint32_t, unordered_set<uint32_t>> &groupInfo,
                                        vector<uint32_t> &groupMappingInfo,
                                        const string &groupFileDir) {
    const string groupInfoFileName = groupFileDir + "/groups";
    const string groupMappingInfoFileName = groupFileDir + "/groupsMapping";

    // 1. Load groupInfo
    ifstream inFile1(groupInfoFileName);
    if (!inFile1.is_open()) {
        cerr << "Error opening file: " << groupInfoFileName << endl;
        return;
    }

    string line;
    while (getline(inFile1, line)) {
        istringstream ss(line);
        uint32_t groupId;
        ss >> groupId;

        uint32_t queryId;
        while (ss >> queryId) {
            groupInfo[groupId].insert(queryId);
        }
    }
    inFile1.close();
    cout << "Group info loaded from " << groupInfoFileName << " successfully." << endl;

    // 2. Load groupMappingInfo
    ifstream inFile2(groupMappingInfoFileName);
    if (!inFile2.is_open()) {
        cerr << "Error opening file: " << groupMappingInfoFileName << endl;
        return;
    }

    groupMappingInfo.clear();
    int groupId;
    while (inFile2 >> groupId) {
        groupMappingInfo.emplace_back(groupId);
    }
    inFile2.close();
    cout << "Query-to-group map loaded from " << groupMappingInfoFileName << " successfully." << endl;
}

void GroupApplier::getRepLabel(vector<OrgResult> &orgResults, 
                               const unordered_map<uint32_t, unordered_set<uint32_t>> &groupInfo, 
                                 unordered_map<uint32_t, uint32_t> &repLabel,
                                 std::unordered_map<int, int>& external2internalTaxId) {
    cout << "Find query group representative labels..." << endl;    
    time_t beforeSearch = time(nullptr);

    for (const auto& group : groupInfo) {
        uint32_t groupId = group.first;
        const unordered_set<uint32_t>& queryIds = group.second;

        vector<WeightedTaxHit> setTaxa;

        for (const auto& queryId : queryIds) {
            int query_label = external2internalTaxId[orgResults[queryId].label]; 
            if (par.weightMode == 0) {
                float score = 1; 
                if (query_label != 0) {
                    setTaxa.emplace_back(query_label, score, 2);
                }
            } else if (par.weightMode == 1) {
                float score = orgResults[queryId].score;
                if (query_label != 0 && score >= par.minVoteScr) {
                    setTaxa.emplace_back(query_label, score, 2);
                }
            } else if (par.weightMode == 2) {
                float score = orgResults[queryId].score;
                if (query_label != 0 && score >= par.minVoteScr) {
                    setTaxa.emplace_back(query_label, score * score, 2);
                }
            }
        }

        WeightedTaxResult result = taxonomy->weightedMajorityLCA(setTaxa, 0.5);

        if (result.taxon != 0 && result.taxon != 1) {
            repLabel[groupId] = result.taxon;
        }
        else{
            repLabel[groupId] = 0;
        }
    }

    const string& groupRepFileName = outDir + "/groupRep";
    ofstream outFile(groupRepFileName);
    if (!outFile.is_open()) {
        cerr << "Error opening file: " << groupRepFileName << endl;
        return;
    }

    for (const auto& [groupId, groupRep] : repLabel) {
        outFile << groupId << "\t" << taxonomy->getOriginalTaxID(groupRep) << "\n";
    }

    outFile.close();

    cout << "Query group representative label saved to " << groupRepFileName << " successfully." << endl;    
    cout << "Time spent: " << double(time(nullptr) - beforeSearch) << " seconds." << endl;
}

void GroupApplier::applyRepLabel(const vector<OrgResult>& orgResults, 
                                   const vector<uint32_t>& groupMappingInfo, 
                                   const unordered_map<uint32_t, uint32_t>& repLabel,
                                   std::unordered_map<int, int>& external2internalTaxId) {
    cout << "Apply query group representative labels..." << endl;    
    time_t beforeSearch = time(nullptr);

    vector<Query> queryList;

    for (int queryIdx = 0; queryIdx < orgResults.size() ; queryIdx++){
        uint32_t groupId = groupMappingInfo[queryIdx];
        auto repLabelIt = repLabel.find(groupId);
        if (repLabelIt != repLabel.end() && repLabelIt->second != 0){
            queryList.emplace_back(Query(repLabelIt->second, orgResults[queryIdx].score, repLabelIt->second, orgResults[queryIdx].name));
        } else{
            queryList.emplace_back(Query(external2internalTaxId[orgResults[queryIdx].label], orgResults[queryIdx].score, orgResults[queryIdx].label, orgResults[queryIdx].name));           

        }
    }    

    reporter->openReadClassificationFile(updatedResultFileName);
    reporter->writeReadClassification(queryList, groupMappingInfo);    
    reporter->closeReadClassificationFile();
    
    cout << "Result saved to " << updatedResultFileName << " successfully." << endl;    
    cout << "Time spent: " << double(time(nullptr) - beforeSearch) << " seconds." << endl;
}