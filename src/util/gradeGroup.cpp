#include "NcbiTaxonomy.h"
#include "Parameters.h"
#include "LocalParameters.h"
#include "IndexCreator.h"

#include <string>
#include <iostream>
#include <regex>
#include <cstdint>

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
    const string groupFileList = par.filenames[0];
    const string readGroupFileList = par.filenames[1];
    const string mappingFileList = par.filenames[2];
    const string taxonomy = par.filenames[3];

    // Parse ranks
    vector<string> ranks;
    if (!par.testRank.empty()) {
        ranks = Util::split(par.testRank, ",");
    } else {
        ranks = {"class", "order", "family", "genus", "species"};
    }

    // Parse print columns
    vector<string> printColumns;
    vector<size_t> printColumnsIdx;
    if (!par.printColumns.empty()) {
        printColumns = Util::split(par.printColumns, ",");
        // stoi
        for (const auto &printColumn : printColumns) {
            printColumnsIdx.push_back(stoi(printColumn));
        }
    }

    // Load Taxonomy
    string names = taxonomy + "/names.dmp";
    string nodes = taxonomy + "/nodes.dmp";
    string merged = taxonomy + "/merged.dmp";

    // Load mapping file names
    ifstream mappingFileListFile;
    mappingFileListFile.open(mappingFileList);
    string eachLine;
    vector<string> mappingFileNames;
    if (mappingFileListFile.is_open()) {
        while (getline(mappingFileListFile, eachLine)) {
            mappingFileNames.push_back(eachLine);
        }
    } else {
        cerr << "Cannot open file for mapping file list" << endl;
    }
    cout << "Answer sheet loaded" << endl;

    // Load group file names
    ifstream groupFileListFile;
    groupFileListFile.open(groupFileList);
    vector<string> groupFileNames;
    if (groupFileListFile.is_open()) {
        while (getline(groupFileListFile, eachLine)) {
            groupFileNames.push_back(eachLine);
        }
    } else {
        cerr << "Cannot open file for read group file list" << endl;
    }
    cout << "Grouping results loaded" << endl;

    // Load readgroup file names
    ifstream readGroupFileListFile;
    readGroupFileListFile.open(readGroupFileList);
    vector<string> readGroupFileNames;
    if (readGroupFileListFile.is_open()) {
        while (getline(readGroupFileListFile, eachLine)) {
            readGroupFileNames.push_back(eachLine);
        }
    } else {
        cerr << "Cannot open file for read group file list" << endl;
    }
    cout << "Read-group mapping results loaded" << endl;

    size_t numberOfFiles = groupFileNames.size();
    vector<GradeResult> results;
    results.resize(numberOfFiles);

#ifdef OPENMP
    omp_set_num_threads(par.threads);
#endif

#pragma omp parallel default(none), shared(results, ranks, numberOfFiles, mappingFileNames, groupFileNames,readGroupFileNames,\
par, cout, printColumnsIdx, cerr, names, nodes, merged)
    {
        // Grade each file
        unordered_map<string, int> assacc2taxid;
        vector<int> classList;
        vector<float> scores;
        string mappingFile;
        string groupFileName;
        string readGroupFileName;

        vector<string> ranks_local = ranks;

        TaxonomyWrapper ncbiTaxonomy(names, nodes, merged, false);
        cout << "Taxonomy loaded" << endl;

#pragma omp for schedule(dynamic)
        for (size_t i = 0; i < numberOfFiles; ++i) {
            // Initialize
            assacc2taxid.clear();
            mappingFile = mappingFileNames[i];
            groupFileName = groupFileNames[i];
            readGroupFileName = readGroupFileNames[i];

            // Load the mapping file (answer sheet) (accession to taxID)
            string key, value;
            ifstream map;
            map.open(mappingFile);
            if (map.is_open()) {
                while (getline(map, key, '\t')) {
                    getline(map, value, '\n');
                    // remove version number
                    size_t pos = key.find('.');
                    if (pos != string::npos) {
                        key = key.substr(0, pos);
                    }
                    assacc2taxid[key] = stoi(value);
                }
            } else {
                cout << "Cannot open file for answer" << endl;
            }
            map.close();
            cout << "Mapping loaded" << endl;

            regex regexName("(GC[AF]_[0-9]+\\.[0-9]+)");

            // Load group results
            string groupResultLine;
            ifstream groupFile;
            groupFile.open(groupFileName);            

            smatch assacc1;
            size_t numGroup = 0; 
            size_t numReadsInGroup = 0; 
            unordered_map<int, vector<int>> group2taxs;

            while (getline(groupFile, groupResultLine, '\n')) {
                if (groupResultLine.empty()) {
                    continue;
                }

                vector<string> fields = Util::split(groupResultLine, "/t");     
                for (int i = 1; i < fields.size(); i++){
                    string id = fields[i];
                    if (par.testType == "gtdb") {
                        regex_search(id, assacc1, regexName);
                        id = assacc1[0];
                        // remove version number
                        size_t pos = id.find('.');
                        if (pos != string::npos) {
                            id = id.substr(0, pos);
                        }
                    } else if (par.testType == "hiv" || par.testType == "hiv-ex") {
                        size_t pos = id.find('_');
                        id = id.substr(0, pos);
                    } else if (par.testType == "cami" || par.testType == "cami-long" || par.testType == "cami-euk") {
                        size_t pos = id.find('/');
                        id = id.substr(0, pos);
                    } else if (par.testType == "over") {
                        regex_search(id, assacc1, regexName);
                        id = assacc1[0];
                    }
                    group2taxs[stoi(fields[0])].push_back(assacc2taxid[id]);   
                    numReadsInGroup++;                    
                }
                numGroup++;
            }
            groupFile.close();
            

            // Load read-group results
            string readGroupResultLine;
            ifstream readGroupFile;
            readGroupFile.open(readGroupFileName);

            smatch assacc2;
            size_t numReads = 0;    
            unordered_map<string, unordered_map<int, vector<int>>> tax2groups;

            while (getline(readGroupFile, readGroupResultLine, '\n')) {
                if (readGroupResultLine.empty()) {
                    continue;
                }
                
                // Parse grouping result
                vector<string> fields = Util::split(readGroupResultLine, "\t");    
                if (stoi(fields[1]) != -1){
                string id = fields[0];
                    if (par.testType == "gtdb") {
                        regex_search(id, assacc2, regexName);
                        id = assacc2[0];
                        // remove version number
                        size_t pos = id.find('.');
                        if (pos != string::npos) {
                            id = id.substr(0, pos);
                        }
                    } else if (par.testType == "hiv" || par.testType == "hiv-ex") {
                        size_t pos = id.find('_');
                        id = id.substr(0, pos);
                    } else if (par.testType == "cami" || par.testType == "cami-long" || par.testType == "cami-euk") {
                        size_t pos = id.find('/');
                        id = id.substr(0, pos);
                    } else if (par.testType == "over") {
                        regex_search(id, assacc2, regexName);
                        id = assacc2[0];
                    }

                    for (const string &rank: ranks_local) {
                        tax2groups[rank][ncbiTaxonomy.getTaxIdAtRank(assacc2taxid[id], rank)].push_back(stoi(fields[1]));     
                    }
                }
                numReads++;

            }
            readGroupFile.close();
            
            // Score the groupinh
            for (const string &rank: ranks_local) {
                countAtRank_CAMI(tax2groups[rank], group2taxs, ncbiTaxonomy,
                                 results[i].countsAtRanks[rank], rank);
                results[i].countsAtRanks[rank].calculate();                 
            }

            // Print Grade Result of each file
            cout << readGroupFileName << endl;
            cout << "The number of reads: " << numReads << endl;
            cout << "The number of groups: " << numGroup << endl;
            cout << "The number of reads in groups: " << numReadsInGroup << endl;
        }
    } // End of parallel region

    cout << "Rank\t";
    for (size_t i = 0; i < results.size(); i++) {
        // cout << "Purity\tRecall\tF1\tARI\t";
        cout << "Purity\tRecall\tF1\t";
    }
    cout << endl;
    for (const string &rank: ranks) {
        cout << rank << "\t";
        for (auto & result : results) {
            // cout << result.countsAtRanks[rank].purity << "\t" << result.countsAtRanks[rank].recall
            //      << "\t" << result.countsAtRanks[rank].f1 << "\t" << result.countsAtRanks[rank].ARI << "\t";
            cout << result.countsAtRanks[rank].purity << "\t" << result.countsAtRanks[rank].recall
                 << "\t" << result.countsAtRanks[rank].f1 << "\t";
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
        
        float groupRecall = 0;
        float numMajorGroupinTax = 0;
        for (auto &p : freq) {
            numMajorGroupinTax = max(numMajorGroupinTax, p.second);
        }
        numMajorGroup += numMajorGroupinTax;
        numTaxGroups += groupIds.size();
    }
    count.recall = numMajorGroup / numTaxGroups;

    auto comb2 = [](long long n) -> long long { return (n < 2) ? 0LL : (n * (n - 1)) / 2; };

    // contingency: groupId -> (labelAtRank -> count)
    // std::unordered_map<int, std::unordered_map<int, long long>> contingency;
    // std::unordered_map<int, long long> groupSizes;  // a_i
    // std::unordered_map<int, long long> labelSizes;  // b_j
    // long long N = 0;

    // for (const auto& [gid, memberTaxs] : group2taxs) {
    //     auto& row = contingency[gid];
    //     long long gi = 0;
    //     for (int t : memberTaxs) {
    //         int lab = ncbiTaxonomy.getTaxIdAtRank(t, rank);
    //         if (lab <= 0) continue;                 // unknown 제외(규칙 일관)
    //         ++row[lab];
    //         ++labelSizes[lab];
    //         ++gi;
    //     }
    //     groupSizes[gid] += gi;
    //     N += gi;
    // }

    // long double sum_nij2 = 0.0L, A = 0.0L, B = 0.0L;
    // for (const auto& [gid, row] : contingency)
    //     for (const auto& [lab, cnt] : row)
    //         sum_nij2 += static_cast<long double>(comb2(cnt));

    // for (const auto& [gid, a] : groupSizes) A += static_cast<long double>(comb2(a));
    // for (const auto& [lab, b] : labelSizes) B += static_cast<long double>(comb2(b));

    // const long double CN2 = static_cast<long double>(comb2(N));
    // double ari = 0.0;
    // if (CN2 > 0.0L) {
    //     const long double expected = (A * B) / CN2;
    //     const long double maxIndex = 0.5L * (A + B);
    //     const long double denom = maxIndex - expected;
    //     const long double numer = sum_nij2 - expected;
    //     ari = (denom != 0.0L) ? static_cast<double>(numer / denom) : 0.0;
    // }
    // count.ARI = ari;
}
