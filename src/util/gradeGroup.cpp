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
    }
};

struct GradeResult{
    unordered_map<string, CountAtRank> countsAtRanks;
    string path;
};

struct Score2{
    Score2(int tf, std::string rank, float score) : tf(tf), rank(rank), score(score) { }
    int tf; // 1 = t, 2 = f
    std::string rank;
    float score;
};

char compareTaxonAtRank_CAMI(TaxID shot, TaxID target, const TaxonomyWrapper & ncbiTaxonomy, CountAtRank & count,
                             const string & rank);

char compareTaxonAtRank_CAMI_euk(TaxID shot, TaxID target, TaxonomyWrapper & ncbiTaxonomy, CountAtRank & count,
                                 const string & rank);

char compareTaxon_overgroup(TaxID shot, TaxID target, TaxonomyWrapper & ncbiTaxonomy, CountAtRank & count,
                                     const string & rank);

char compareTaxon_hivExclusion(TaxID shot, TaxID target, CountAtRank & count);

void setGradeDefault(LocalParameters & par){
    par.groupIdCol = 0;
    par.memberIdCol = 1;
    // par.verbosity = 2;
    // par.scoreCol = 0;
    par.testRank = "";
    par.testType = "gtdb";
    // par.skipSecondary = 0;
}

int gradeGroup(int argc, const char **argv, const Command &command) {

    LocalParameters &par = LocalParameters::getLocalInstance();
    setGradeDefault(par);
    par.parseParameters(argc, argv, command, false, Parameters::PARSE_ALLOW_EMPTY, 0);
    // const string readClassificationFileList = par.filenames[0];
    const string readGroupFileList = par.filenames[0];
    const string mappingFileList = par.filenames[1];
    const string taxonomy = par.filenames[2];

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
    cout << "Grouping results loaded" << endl;

    size_t numberOfFiles = readGroupFileNames.size();
    vector<GradeResult> results;
    results.resize(numberOfFiles);

#ifdef OPENMP
    omp_set_num_threads(par.threads);
#endif

#pragma omp parallel default(none), shared(results, ranks, numberOfFiles, mappingFileNames, readGroupFileNames,\
par, cout, printColumnsIdx, cerr, names, nodes, merged)
    {
        // Grade each file
        unordered_map<string, int> assacc2taxid;
        vector<int> classList;
        vector<float> scores;
        string mappingFile;
        string readGroupFileName;

        vector<string> ranks_local = ranks;

        TaxonomyWrapper ncbiTaxonomy(names, nodes, merged, false);
        cout << "Taxonomy loaded" << endl;

        // Print scores of TP and FP
        unordered_map<string, vector<size_t>> rank2TpIdx;
        unordered_map<string, vector<size_t>> rank2FpIdx;
        unordered_map<string, vector<size_t>> rank2FnIdx;
        vector<vector<string>> idx2values;
        if (!printColumnsIdx.empty()){
            for (const auto & rank : ranks_local) {
                rank2TpIdx[rank] = vector<size_t>();
                rank2FpIdx[rank] = vector<size_t>();
                rank2FnIdx[rank] = vector<size_t>();
            }
        }

#pragma omp for schedule(dynamic)
        for (size_t i = 0; i < numberOfFiles; ++i) {
            // Initialize
            assacc2taxid.clear();
            vector<int> readTaxs;
            vector<string> readIds;
            
            if (!printColumnsIdx.empty()){
                for (const auto & rank : ranks_local) {
                    rank2TpIdx[rank].clear();
                    rank2FpIdx[rank].clear();
                    rank2FnIdx[rank].clear();
                }
            }
            mappingFile = mappingFileNames[i];
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

            // Load group results
            string resultLine;
            ifstream readGroupFile;
            readGroupFile.open(readGroupFileName);
            vector<string> fields;
            string field;

            vector<Score> tpOrFp;
            regex regex1("(GC[AF]_[0-9]+\\.[0-9]+)");
            smatch assacc;
            size_t numberOfGroups = 0;
            size_t numberOfReadinGroups = 0;
            while (getline(readGroupFile, resultLine, '\n')) {
                // skip line starting with '#'
                if (resultLine.empty()) {
                    continue;
                }
                
                // Parse grouping result
                fields = Util::split(resultLine, "\t");                     
                span<string> groupMembers(fields.begin() + memberIdCol, fields.end());
                for (int i = 0; i< groupMembers.size(); i++){
                    string fullId = id;
                    if (par.testType == "gtdb") {
                        regex_search(id, assacc, regex1);
                        id = assacc[0];
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
                        regex_search(id, assacc, regex1);
                        id = assacc[0];
                    }
                    readIds.push_back(fullId);
                    readTaxs.push_back(assacc2taxid[id]);   
                    numberOfReadinGroups++;
                }
                numberOfGroups++;
            }
            readGroup.close();

            // Score the grouping
            char p;
            for (size_t j = 0; j < classList.size(); j++) {
                if (par.verbosity == 3) cout << readIds[j] << " " << classList[j] << " " << rightAnswers[j];
                for (const string &rank: ranks_local) {
                    if (par.testType == "over") {
                        p = compareTaxon_overGroup(classList[j], rightAnswers[j], ncbiTaxonomy,
                                                            results[i].countsAtRanks[rank], rank);
                    } else if(par.testType == "hiv-ex"){
                        p = compareTaxon_hivExclusion(classList[j], 11676, results[i].countsAtRanks[rank]);
                    } else if (par.testType == "cami-euk"){
                        p = compareTaxonAtRank_CAMI_euk(classList[j], rightAnswers[j], ncbiTaxonomy,
                                                        results[i].countsAtRanks[rank], rank);
                    } else {
                        p = compareTaxonAtRank_CAMI(classList[j], rightAnswers[j], ncbiTaxonomy,
                                                         results[i].countsAtRanks[rank], rank);
                    }
                    if (!printColumnsIdx.empty()) {
                        if (p == 'O') rank2TpIdx[rank].push_back(j);
                        else if (p == 'X') rank2FpIdx[rank].push_back(j);
                        else if (p == 'N') rank2FnIdx[rank].push_back(j);
                    }
                    if (par.verbosity == 3) cout << " " << p;
                }
                if (par.verbosity == 3) cout << endl;
            }

            // Calculate the scores
            for (const string &rank: ranks_local) {
                results[i].countsAtRanks[rank].calculate();
            }

            // Print Grade Result of each file
            cout << readGroupFileName << endl;
            cout << "The number of reads: " << rightAnswers.size() << endl;
            cout << "The number of groups: " << numberOfGroups << endl;
            cout << "The number of reads in groups: " << numberOfReadinGroups << endl;
            for (const string &rank: ranks_local) {
                cout << rank << " " << results[i].countsAtRanks[rank].total << " "
                     << results[i].countsAtRanks[rank].TP + results[i].countsAtRanks[rank].FP << " "
                     << results[i].countsAtRanks[rank].TP << " " << results[i].countsAtRanks[rank].FP << " "
                     << results[i].countsAtRanks[rank].precision << " "
                     << results[i].countsAtRanks[rank].sensitivity << " " << results[i].countsAtRanks[rank].f1 << endl;
            }
            cout << endl;
        }
    } // End of parallel region

    cout << "Rank\t";
    for (size_t i = 0; i < results.size(); i++) {
        cout << "Precision\tSensitivity\tF1\t";
    }
    cout << endl;
    for (const string &rank: ranks) {
        cout << rank << "\t";
        for (auto & result : results) {
            cout << result.countsAtRanks[rank].precision << "\t" << result.countsAtRanks[rank].sensitivity
                 << "\t" << result.countsAtRanks[rank].f1 << "\t";
        }
        cout << endl;
    }
    return 0;
}

char compareTaxonAtRank_CAMI(TaxID shot, TaxID target, const TaxonomyWrapper & ncbiTaxonomy, CountAtRank & count,
                             const string & rank) {
    if (rank == "subspecies") {
        // Do not count if the rank of target is higher than current rank
        // current rank is subspecies
        // the rank of target is subspecies


        // False negative; no groups or meaningless groups
        if (shot == 1 || shot == 0) {
            count.FN ++;
            count.total ++;
            return 'N';
        }

        // False negative if the rank of shot is higher than current rank
        const TaxonNode * shotNode = ncbiTaxonomy.taxonNode(shot);
        if (strcmp(ncbiTaxonomy.getString(shotNode->rankIdx), "no rank") != 0) { // no rank is subspecies
            // cout << ncbiTaxonomy.getString(shotNode->rankIdx) << endl;
            count.FN ++;
            count.total ++;
            return 'N';
        }

        count.total++;
        if(shot == target){
            count.TP++;
            return 'O';
        } else {
            count.FP++;
            return 'X';
        }
    } else {
        // Do not count if the rank of target is higher than current rank
        TaxID targetTaxIdAtRank = ncbiTaxonomy.getTaxIdAtRank(target, rank);
        // cout << targetTaxIdAtRank << endl;
        const TaxonNode * targetNode = ncbiTaxonomy.taxonNode(targetTaxIdAtRank);
        int rankIdx = ncbiTaxonomy.findRankIndex2(rank);
        // cout << shot << " " << targetTaxIdAtRank << " " << targetNode->rankIdx << " " << endl;
        if (ncbiTaxonomy.findRankIndex2(ncbiTaxonomy.getString(targetNode->rankIdx)) > rankIdx) {
            return '-';
        }
        
        // False negative; no group or meaningless group
        if(shot == 1 || shot == 0) {
            count.FN ++;
            count.total ++;
            return 'N';
        }

        // False negative if the rank of shot is higher than current rank
        TaxID shotTaxIdAtRank = ncbiTaxonomy.getTaxIdAtRank(shot, rank);
        const TaxonNode * shotNode = ncbiTaxonomy.taxonNode(shotTaxIdAtRank);
        if (ncbiTaxonomy.findRankIndex2(ncbiTaxonomy.getString(shotNode->rankIdx)) > rankIdx) {
            count.FN ++;
            count.total ++;
            return 'N';
        }   
        count.total++;
        if(shotTaxIdAtRank == targetTaxIdAtRank){
            count.TP++;
            return 'O';
        } else {
            count.FP++;
            return 'X';
        }
    }
}

char compareTaxonAtRank_CAMI_euk(TaxID shot, TaxID target, TaxonomyWrapper & ncbiTaxonomy, CountAtRank & count,
                             const string & rank) {
    // Do not count if the rank of target is higher than current rank
    TaxID targetTaxIdAtRank = ncbiTaxonomy.getTaxIdAtRank(target, rank);
    const TaxonNode * targetNode = ncbiTaxonomy.taxonNode(targetTaxIdAtRank);
    int rankIdx = NcbiTaxonomy::findRankIndex(rank);
    if (NcbiTaxonomy::findRankIndex(ncbiTaxonomy.getString(targetNode->rankIdx)) > rankIdx) {
        return '-';
    }

    // Do not count if target is not eukaryote
    if (ncbiTaxonomy.getTaxIdAtRank(target, "superkingdom") != 2759) {
        return '-';
    }

    // False negative; no group or meaningless group
    if(shot == 1 || shot == 0) {
        count.FN ++;
        count.total ++;
        return 'N';
    }

    // False negative if the rank of shot is higher than current rank
    TaxID shotTaxIdAtRank = ncbiTaxonomy.getTaxIdAtRank(shot, rank);
    const TaxonNode * shotNode = ncbiTaxonomy.taxonNode(shotTaxIdAtRank);
    if (NcbiTaxonomy::findRankIndex(ncbiTaxonomy.getString(shotNode->rankIdx)) > rankIdx) {
        count.FN ++;
        count.total ++;
        return 'N';
    }

    count.total++;
    if(shotTaxIdAtRank == targetTaxIdAtRank){
        count.TP++;
        return 'O';
    } else {
        count.FP++;
        return 'X';
    }
}

char compareTaxon_overGroup(TaxID shot, TaxID target, TaxonomyWrapper & ncbiTaxonomy, CountAtRank & count,
                                     const string & rank){
    // Do not count if the rank of target is higher than current rank
//    TaxID targetTaxIdAtRank = ncbiTaxonomy.getTaxIdAtRank(target, rank);
    const TaxonNode * targetNode = ncbiTaxonomy.taxonNode(target);
    int rankIdx = NcbiTaxonomy::findRankIndex(rank);
    if (NcbiTaxonomy::findRankIndex(ncbiTaxonomy.getString(targetNode->rankIdx)) > rankIdx) {
        return '-';
    }


    // False negative; no group or meaningless group
    if(shot == 1 || shot == 0) {
        count.FN ++;
        count.total ++;
        return 'N';
    }

    // False negative if the rank of shot is higher than current rank
    const TaxonNode * shotNode = ncbiTaxonomy.taxonNode(shot);
    if (NcbiTaxonomy::findRankIndex(ncbiTaxonomy.getString(shotNode->rankIdx)) > rankIdx) {
        count.FN ++;
        count.total ++;
        return 'N';
    }

    count.total++;
    if(shot == target){
        count.TP++;
        return 'O';
    } else {
        count.FP++;
        return 'X';
    }
}

// TP: HIV-1 at species rank
// FP: Groups to other taxa
// FN: Not-classified
char compareTaxon_hivExclusion(TaxID shot, TaxID target, CountAtRank & count){
    // False negative; no group or meaningless group
    if(shot == 1 || shot == 0) {
        count.FN ++;
        count.total ++;
        return 'N';
    }
    count.total++;
    if(shot == target){
        count.TP++;
        return 'O';
    } else {
        count.FP++;
        return 'X';
    }
}
