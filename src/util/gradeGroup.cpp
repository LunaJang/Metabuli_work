#include "NcbiTaxonomy.h"
#include "Parameters.h"
#include "LocalParameters.h"
#include "IndexCreator.h"

#include <string>
#include <iostream>
#include <regex>
#include <cstdint>
#include <fstream>
#include <unordered_map>
#include <vector>

using namespace std;

struct CountAtRank {
    int total;
    float purity;
    float recall;
    float f1;
    float ARI;
    void calculate() {
        f1 = 2 * purity * recall / (purity + recall);
        ARI = 0.0;
    }
};

struct GradeResult{
    unordered_map<string, CountAtRank> countsAtRanks;
    string path;
};

void countAtRank_CAMI(const unordered_map<int, vector<int>>& tax2groupsAtRank, 
                      const unordered_map<int, vector<int>>& group2taxs, 
                      const TaxonomyWrapper & ncbiTaxonomy, CountAtRank & count, const string & rank);

void setGradeGroupDefault(LocalParameters & par){
    par.testRank = "";
    par.testType = "gtdb";
}

int gradeGroup(int argc, const char **argv, const Command &command) {

    LocalParameters &par = LocalParameters::getLocalInstance();
    setGradeGroupDefault(par);
    par.parseParameters(argc, argv, command, false, Parameters::PARSE_ALLOW_EMPTY, 0);

    const string groupFileList     = par.filenames[0];
    const string readGroupFileList = par.filenames[1];
    const string mappingFileList   = par.filenames[2];
    const string queryNameFileList = par.filenames[3]; // NEW
    const string taxonomy          = par.filenames[4];

    // Parse ranks
    vector<string> ranks;
    if (!par.testRank.empty()) {
        ranks = Util::split(par.testRank, ",");
    } else {
        ranks = {"class", "order", "family", "genus", "species"};
    }

    // Load Taxonomy
    string names  = taxonomy + "/names.dmp";
    string nodes  = taxonomy + "/nodes.dmp";
    string merged = taxonomy + "/merged.dmp";

    string eachLine;

    // Load mapping file names
    ifstream mappingFileListFile(mappingFileList);
    vector<string> mappingFileNames;
    if (mappingFileListFile.is_open()) {
        while (getline(mappingFileListFile, eachLine)) {
            if (!eachLine.empty())
                mappingFileNames.push_back(eachLine);
        }
    } else {
        cerr << "Cannot open file for mapping file list" << endl;
    }
    cout << "Answer sheet list loaded" << endl;

    // Load group file names
    ifstream groupFileListFile(groupFileList);
    vector<string> groupFileNames;
    if (groupFileListFile.is_open()) {
        while (getline(groupFileListFile, eachLine)) {
            if (!eachLine.empty())
                groupFileNames.push_back(eachLine);
        }
    } else {
        cerr << "Cannot open file for group file list" << endl;
    }
    cout << "Grouping results list loaded" << endl;

    // Load readgroup file names
    ifstream readGroupFileListFile(readGroupFileList);
    vector<string> readGroupFileNames;
    if (readGroupFileListFile.is_open()) {
        while (getline(readGroupFileListFile, eachLine)) {
            if (!eachLine.empty())
                readGroupFileNames.push_back(eachLine);
        }
    } else {
        cerr << "Cannot open file for read group file list" << endl;
    }
    cout << "Read-group mapping results list loaded" << endl;

    // NEW: Load query-name file names
    ifstream queryNameFileListFile(queryNameFileList);
    vector<string> queryNameFileNames;
    if (queryNameFileListFile.is_open()) {
        while (getline(queryNameFileListFile, eachLine)) {
            if (!eachLine.empty())
                queryNameFileNames.push_back(eachLine);
        }
    } else {
        cerr << "Cannot open file for query name file list" << endl;
    }
    cout << "Query-name list files loaded" << endl;

    size_t numberOfFiles = groupFileNames.size();
    vector<GradeResult> results(numberOfFiles);

#ifdef OPENMP
    omp_set_num_threads(par.threads);
#endif

#pragma omp parallel default(none), \
    shared(results, ranks, numberOfFiles, mappingFileNames, groupFileNames, \
           readGroupFileNames, queryNameFileNames, \
           par, cout, cerr, names, nodes, merged)
    {
        unordered_map<string, int> assacc2taxid;
        string mappingFile;
        string groupFileName;
        string readGroupFileName;
        string queryNameFileName;

        vector<string> ranks_local = ranks;

        TaxonomyWrapper ncbiTaxonomy(names, nodes, merged, false);
        cout << "Taxonomy loaded" << endl;

        regex regexName("(GC[AF]_[0-9]+\\.[0-9]+)");

#pragma omp for schedule(dynamic)
        for (size_t i = 0; i < numberOfFiles; ++i) {
            assacc2taxid.clear();
            mappingFile       = mappingFileNames[i];
            groupFileName     = groupFileNames[i];
            readGroupFileName = readGroupFileNames[i];
            queryNameFileName = queryNameFileNames[i];

            {
                string key, value;
                ifstream map(mappingFile);
                if (map.is_open()) {
                    while (getline(map, key, '\t')) {
                        if (!getline(map, value, '\n')) break;
                        // remove version number
                        size_t pos = key.find('.');
                        if (pos != string::npos) {
                            key = key.substr(0, pos);
                        }
                        assacc2taxid[key] = stoi(value);
                    }
                } else {
                    cerr << "Cannot open file for answer: " << mappingFile << endl;
                }
            }
            cout << "Mapping loaded: " << mappingFile << endl;

            vector<int> queryTaxid;          // index 0 -> query index 1
            {
                ifstream qf(queryNameFileName);
                string qline;
                while (getline(qf, qline)) {
                    if (qline.empty()) {
                        queryTaxid.push_back(-1);
                        continue;
                    }
                    if (qline[0] == '>')
                        qline = qline.substr(1);

                    string id = qline;

                    if (par.testType == "gtdb") {
                        smatch m;
                        if (regex_search(qline, m, regexName)) {
                            id = m[0];
                            size_t pos = id.find('.');
                            if (pos != string::npos)
                                id = id.substr(0, pos);
                        } else {
                            queryTaxid.push_back(-1);
                            continue;
                        }
                    } else if (par.testType == "hiv" || par.testType == "hiv-ex") {
                        size_t pos = id.find('_');
                        if (pos != string::npos) id = id.substr(0, pos);
                    } else if (par.testType == "cami" || par.testType == "cami-long" || par.testType == "cami-euk") {
                        size_t pos = id.find('/');
                        if (pos != string::npos) id = id.substr(0, pos);
                    } else if (par.testType == "over") {
                        smatch m;
                        if (regex_search(qline, m, regexName)) {
                            id = m[0];
                        }
                    }

                    auto it = assacc2taxid.find(id);
                    if (it == assacc2taxid.end()) {
                        queryTaxid.push_back(-1);
                    } else {
                        queryTaxid.push_back(it->second);
                    }
                }
            }
            cout << "Query names loaded: " << queryNameFileName 
                 << " (size=" << queryTaxid.size() << ")" << endl;

            string groupResultLine;
            ifstream groupFile(groupFileName);            

            size_t numGroup = 0; 
            size_t numReadsInGroup = 0; 
            unordered_map<int, vector<int>> group2taxs;

            while (getline(groupFile, groupResultLine)) {
                if (groupResultLine.empty())
                    continue;

                vector<string> fields = Util::split(groupResultLine, "\t");
                if (fields.size() < 2)
                    continue;

                int groupId = stoi(fields[0]);

                for (size_t j = 1; j < fields.size(); j++) {
                    int qIdx = stoi(fields[j]);  // 1-based
                    size_t vecIdx = static_cast<size_t>(qIdx - 1);

                    int taxid = queryTaxid[vecIdx];
                    if (taxid < 0) {
                        continue;
                    }

                    group2taxs[groupId].push_back(taxid);
                    numReadsInGroup++;
                }
                numGroup++;
            }
            groupFile.close();

            string readGroupResultLine;
            ifstream readGroupFile(readGroupFileName);

            size_t numReads = 0;
            unordered_map<string, unordered_map<int, vector<int>>> tax2groups;

            while (getline(readGroupFile, readGroupResultLine)) {
                if (readGroupResultLine.empty())
                    continue;

                vector<string> fields = Util::split(readGroupResultLine, "\t");
                if (fields.size() < 2)
                    continue;

                int qIdx = stoi(fields[0]);   // 1-based query index
                int gId  = stoi(fields[1]);   // group id

                numReads++;

                if (gId == 0)
                    continue;

                size_t vecIdx = static_cast<size_t>(qIdx - 1);
                if (vecIdx >= queryTaxid.size() || qIdx <= 0) {
                    cerr << "Warning: qIdx " << qIdx << " out of range in read-group file "
                         << readGroupFileName << endl;
                    continue;
                }

                int taxid = queryTaxid[vecIdx];
                if (taxid < 0)
                    continue;

                for (const string &rank : ranks_local) {
                    int taxAtRank = ncbiTaxonomy.getTaxIdAtRank(taxid, rank);
                    tax2groups[rank][taxAtRank].push_back(gId);
                }
            }
            readGroupFile.close();

            for (const string &rank: ranks_local) {
                countAtRank_CAMI(tax2groups[rank], group2taxs, ncbiTaxonomy,
                                 results[i].countsAtRanks[rank], rank);
                results[i].countsAtRanks[rank].calculate();                 
            }

            cout << readGroupFileName << endl;
            cout << "The number of reads: " << numReads << endl;
            cout << "The number of groups: " << numGroup << endl;
            cout << "The number of reads in groups: " << numReadsInGroup << endl;
        }
    } // End of parallel region

    cout << "Rank\t";
    for (size_t i = 0; i < results.size(); i++) {
        cout << "Purity\tRecall\tF1\t";
    }
    cout << endl;

    for (const string &rank: ranks) {
        cout << rank << "\t";
        for (auto & result : results) {
            cout << result.countsAtRanks[rank].purity << "\t"
                 << result.countsAtRanks[rank].recall << "\t"
                 << result.countsAtRanks[rank].f1 << "\t";
        }
        cout << endl;
    }
    return 0;
}

void countAtRank_CAMI(const unordered_map<int, vector<int>>& tax2groupsAtRank, 
                      const unordered_map<int, vector<int>>& group2taxs, 
                      const TaxonomyWrapper& ncbiTaxonomy, CountAtRank& count, const string& rank) {
    // purity (number of reads with major taxonomy in each group)
    float numMajorTaxs = 0;
    float numGroupedReads = 0;
    for (const auto& [groupId, memberTaxs] : group2taxs){
        std::unordered_map<TaxID, float> freq;
        for (const auto& memberTax : memberTaxs){
            freq[ncbiTaxonomy.getTaxIdAtRank(memberTax, rank)]++;
        }
        
        float numMajorTaxsinGroup = 0;
        for (auto &p : freq) {
            numMajorTaxsinGroup = std::max(numMajorTaxsinGroup, p.second);
        }
        numMajorTaxs += numMajorTaxsinGroup;
        numGroupedReads += memberTaxs.size();
    }
    count.purity = numMajorTaxs / numGroupedReads;

    // recall (number of reads with same taxonomy from one group containing most of the reads / number of reads with same taxonomy)
    float numMajorGroup = 0;
    float numTaxGroups = 0;
    for (const auto& [tax, groupIds] : tax2groupsAtRank){
        std::unordered_map<TaxID, float> freq;
        for (const auto& groupId : groupIds){
            freq[groupId]++;
        }
        
        float numMajorGroupinTax = 0;
        for (auto &p : freq) {
            numMajorGroupinTax = max(numMajorGroupinTax, p.second);
        }
        numMajorGroup += numMajorGroupinTax;
        numTaxGroups += groupIds.size();
    }
    count.recall = (numTaxGroups == 0.0f) ? 0.0f : numMajorGroup / numTaxGroups;
}
