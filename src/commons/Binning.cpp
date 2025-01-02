#include "Binning.h"
#include "FileUtil.h"
#include "common.h"

Binning::Binning(LocalParameters & par) {
    // Load parameters
    dbDir = par.filenames[0];
    loadDbParameters(par);
    
    cout << "DB name: " << par.dbName << endl;
    cout << "DB creation date: " << par.dbDate << endl;
    
    // Taxonomy
    taxonomy = loadTaxonomy(dbDir, par.taxonomyPath);
    
}

Binning::~Binning() {
    delete taxonomy;
}

void Binning::startBinning(const LocalParameters &par) {  
    size_t processedReadCnt = 0;

    string outDir;
    string binningOutDir;
    string jobId;
    outDir = par.filenames[1];
    binningOutDir = par.filenames[1] + "/binning_result";
    jobId = par.filenames[2];

    size_t voteMode = par.voteMode;
    float majorityThr = par.majorityThr;
    cout << "voteMode: " << voteMode << endl;
    cout << "majorityThr: " << majorityThr << endl;

    unordered_map<uint32_t, unordered_set<string>> binningReadGroupInfo;
    unordered_map<string, uint32_t> queryGroupInfo;
    //queryGroupInfo.resize(processedReadCnt, -1);
    makeGroupsFromBinning(binningOutDir, binningReadGroupInfo, queryGroupInfo);

    vector<pair<string, float>> binnningMetabuliResult;       
    //binnningMetabuliResult.resize(processedReadCnt, make_pair(-1, 0.0f));
    loadMetabuliResult(binningOutDir, binnningMetabuliResult);

    unordered_map<uint32_t, int> repLabel; 
    getRepLabel(binningOutDir, binnningMetabuliResult, binningReadGroupInfo, repLabel, jobId, voteMode, majorityThr);

    applyRepLabel(outDir, queryGroupInfo, repLabel, jobId);
    
    return;
}

void Binning::makeGroupsFromBinning(const string &groupFileDir, 
                                    unordered_map<uint32_t, unordered_set<string>> &binningReadGroupInfo, 
                                    unordered_map<string, uint32_t> &queryGroupInfo) {
    // bin에 포함된 new reads를 읽어서 group으로 묶는다. group: reads
    const string& binningReadGroupFileName = groupFileDir + "/binning_read_ids.tsv";
    std::ifstream inFile(binningReadGroupFileName);
    if (!inFile.is_open()) {
        std::cerr << "Error opening file: " << binningReadGroupFileName << std::endl;
        return;
    }

    string line;
    while (getline(inFile, line)) {
        stringstream ss(line);
        uint32_t groupId;
        string binningReadId;
        ss >> groupId >> binningReadId;
        binningReadGroupInfo[groupId].insert(binningReadId);
    }

    inFile.close();
    std::cout << "Group info loaded from " << binningReadGroupFileName << " successfully." << std::endl;

    // bin에 포함된 query reads를 읽어서 group으로 묶는다. query: groups
    const string& queryGroupFileName = groupFileDir + "/binning_query_ids.tsv";
    std::ifstream inFile(queryGroupFileName);
    if (!inFile.is_open()) {
        std::cerr << "Error opening file: " << queryGroupFileName << std::endl;
        return;
    }

    string line;
    while (getline(inFile, line)) {
        stringstream ss(line);
        uint32_t groupId;
        string queryId;
        ss >> groupId >> queryId;
        queryGroupInfo[queryId] = groupId;
    }

    inFile.close();
    std::cout << "Query group info loaded from " << queryGroupFileName << " successfully." << std::endl;
}

void Binning::loadMetabuliResult(const string &resultFileDir, 
                                 vector<pair<string, float>> &binnningMetabuliResult) {
    // binning reads의 metabuli result를 읽어온다.
    ifstream inFile(resultFileDir + "/1_classifications.tsv");
    if (!inFile.is_open()) {
        cerr << "Error opening file: " << resultFileDir + "/1_classifications.tsv" << endl;
        return;
    }

    string line;
    while (getline(inFile, line)) {
        stringstream ss(line);
        int label, query_label, read_length;
        float score;
        string query_name;
        ss >> label >> query_name >> query_label >> read_length >> score;
        binnningMetabuliResult[query_name] = {query_label, score}; 
    }

    inFile.close();
    cout << "Metabuli result for binning reads loaded from " << resultFileDir + "/1_classifications.tsv" << " successfully." << endl;
}

void Binning::getRepLabel(const string &groupRepFileDir, 
                          vector<pair<int, float>> &binnningMetabuliResult, 
                          unordered_map<uint32_t, unordered_set<string>> &binningReadGroupInfo, 
                          unordered_map<uint32_t, int> &repLabel, 
                          const string &jobId, 
                          int voteMode, 
                          float majorityThr) {
    for (const auto& group : binningReadGroupInfo) {
        uint32_t groupId = group.first;
        const unordered_set<string>& binningReadIds = group.second;

        vector<WeightedTaxHit> setTaxa;

        for (const auto& binningReadId : binningReadIds) {
            int query_label = binnningMetabuliResult[binningReadId].first; 
            float score = binnningMetabuliResult[binningReadId].second;
            if (query_label != 0) {
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

// Todo: 하나의 query가 여러 bin에 포함 될 경우, 각 bin의 representative label에 LCA를 적용하여 query의 representative label을 결정한다.
void Binning::applyRepLabel(const string &originalResultFileDir, 
                            unordered_map<string, uint32_t> queryGroupInfo, 
                            const unordered_map<uint32_t, int> &repLabel, 
                            const string &jobId) {
    ifstream inFile(originalResultFileDir + "/" + jobId + "_classifications.tsv");
    if (!inFile.is_open()) {
        cerr << "Error opening file: " << originalResultFileDir + "/" + jobId + "_classifications.tsv" << endl;
        return;
    }

    ofstream outFile(originalResultFileDir + "/" + jobId + "_updated_classifications.tsv");
    if (!outFile.is_open()) {
        cerr << "Error opening file: " << originalResultFileDir + "/" + jobId + "_updated_classifications.tsv" << endl;
        return;
    }

    string line;
    while (getline(inFile, line)) {
        stringstream ss(line);
        vector<string> fields;
        string field;

        while (getline(ss, field, '\t')) {
            fields.push_back(field);
        }

        if (fields.size() > 2 && fields[0] == "0") {
            auto queryGroupInfoIt = queryGroupInfo.find(fields[1]);
            if (queryGroupInfoIt != queryGroupInfo.end() && queryGroupInfoIt->second != 0) {
                auto repLabelIt = repLabel.find(queryGroupInfoIt->second);
                if (repLabelIt != repLabel.end() && repLabelIt->second != 0) {
                    fields[2] = to_string(repLabelIt->second);
                    fields[0] = "1";
                    fields[5] = taxonomy->getString(taxonomy->taxonNode(repLabelIt->second)->rankIdx);
                }
            }
        }

        for (size_t i = 0; i < fields.size(); ++i) {
            outFile << fields[i]; 
            if (i < fields.size() - 1) {  
                outFile << "\t";
            }
        }
        outFile << endl;
    }

    inFile.close();
    outFile.close();
    
    cout << "Result saved to " << newResultFileDir << " successfully." << endl;
}