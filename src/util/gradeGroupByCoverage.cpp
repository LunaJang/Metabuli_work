#include "NcbiTaxonomy.h"
#include "Parameters.h"
#include "LocalParameters.h"
#include "IndexCreator.h"

#include <string>
#include <iostream>
#include <fstream>
#include <regex>
#include <map>
#include <set>
#include <unordered_map>
#include <vector>
#include <algorithm>
#include <cstdint>

using namespace std;

struct CovGroupCount {
    float purity = 0.0f, recall = 0.0f, f1 = 0.0f;
    void calculate() {
        f1 = (purity + recall > 0.0f)
             ? 2.0f * purity * recall / (purity + recall)
             : 0.0f;
    }
};

struct GradeGroupByCovResult {
    map<float, unordered_map<string, CovGroupCount>> cov2rank;
    string path;
};

static void computeGroupPurityRecallByCov(
    const unordered_map<int, vector<TaxID>> &group2taxs,
    const unordered_map<TaxID, vector<int>> &tax2groups,
    CovGroupCount &count)
{
    float numMajor = 0.0f, totalGrouped = 0.0f;
    for (const auto &kv : group2taxs) {
        unordered_map<TaxID, int> freq;
        for (TaxID t : kv.second) freq[t]++;
        int mx = 0;
        for (auto &p : freq) mx = max(mx, p.second);
        numMajor     += (float)mx;
        totalGrouped += (float)kv.second.size();
    }
    count.purity = (totalGrouped > 0.0f) ? numMajor / totalGrouped : 0.0f;

    float numMajorGroup = 0.0f, totalTaxReads = 0.0f;
    for (const auto &kv : tax2groups) {
        unordered_map<int, int> freq;
        for (int g : kv.second) freq[g]++;
        int mx = 0;
        for (auto &p : freq) mx = max(mx, p.second);
        numMajorGroup += (float)mx;
        totalTaxReads += (float)kv.second.size();
    }
    count.recall = (totalTaxReads > 0.0f) ? numMajorGroup / totalTaxReads : 0.0f;
}

void setGradeGroupByCoverageDefault(LocalParameters &par) {
    par.testRank = "";
    par.testType = "gtdb";
}

int gradeGroupByCoverage(int argc, const char **argv, const Command &command) {
    LocalParameters &par = LocalParameters::getLocalInstance();
    setGradeGroupByCoverageDefault(par);
    par.parseParameters(argc, argv, command, false, Parameters::PARSE_ALLOW_EMPTY, 0);

    const string groupFileList     = par.filenames[0];
    const string readGroupFileList = par.filenames[1];
    const string mappingFileList   = par.filenames[2];
    const string queryNameFileList = par.filenames[3];
    const string covMappingFile    = par.filenames[4]; // accession<whitespace>coverage
    const string taxonomy          = par.filenames[5];

    vector<string> ranks = par.testRank.empty()
        ? vector<string>{"class", "order", "family", "genus", "species"}
        : Util::split(par.testRank, ",");

    string names  = taxonomy + "/names.dmp";
    string nodes  = taxonomy + "/nodes.dmp";
    string merged = taxonomy + "/merged.dmp";

    // accession (no version) -> coverage
    unordered_map<string, float> assacc2cov;
    {
        ifstream f(covMappingFile);
        if (!f.is_open()) { cerr << "Cannot open: " << covMappingFile << endl; exit(1); }
        string acc; float cov;
        while (f >> acc >> cov) {
            size_t dot = acc.find('.');
            if (dot != string::npos) acc = acc.substr(0, dot);
            assacc2cov[acc] = cov;
        }
    }
    cerr << "Coverage map loaded: " << assacc2cov.size() << " species" << endl;

    auto loadList = [](const string &path) {
        vector<string> v;
        ifstream f(path); string line;
        while (getline(f, line)) if (!line.empty()) v.push_back(line);
        return v;
    };
    vector<string> groupFiles     = loadList(groupFileList);
    vector<string> readGroupFiles = loadList(readGroupFileList);
    vector<string> mappingFiles   = loadList(mappingFileList);
    vector<string> queryNameFiles = loadList(queryNameFileList);

    size_t N = groupFiles.size();
    vector<GradeGroupByCovResult> results(N);

#ifdef OPENMP
    omp_set_num_threads(par.threads);
#endif

#pragma omp parallel default(none), \
    shared(results, ranks, N, mappingFiles, readGroupFiles, queryNameFiles, \
           assacc2cov, par, cerr, names, nodes, merged)
    {
        unordered_map<string, int> assacc2taxid;
        TaxonomyWrapper ncbiTax(names, nodes, merged, false);
        regex accRe("(GC[AF]_[0-9]+\\.[0-9]+)");
        smatch m;
        vector<string> ranks_local = ranks;
        cerr << "Taxonomy loaded" << endl;

#pragma omp for schedule(dynamic)
        for (size_t i = 0; i < N; ++i) {
            assacc2taxid.clear();

            {
                string k, v;
                ifstream mf(mappingFiles[i]);
                if (!mf.is_open()) { cerr << "Cannot open: " << mappingFiles[i] << endl; continue; }
                while (getline(mf, k, '\t') && getline(mf, v, '\n')) {
                    if (!v.empty() && v.back() == '\r') v.pop_back();
                    if (v.empty() || !isdigit((unsigned char)v[0])) continue;
                    size_t d = k.find('.');
                    if (d != string::npos) k = k.substr(0, d);
                    assacc2taxid[k] = stoi(v);
                }
            }

            // Build queryTaxid and queryCov indexed by 0-based query position
            vector<TaxID> queryTaxid;
            vector<float> queryCov;
            {
                ifstream qf(queryNameFiles[i]);
                if (!qf.is_open()) { cerr << "Cannot open: " << queryNameFiles[i] << endl; continue; }
                string qline;
                while (getline(qf, qline)) {
                    if (qline.empty()) { queryTaxid.push_back(-1); queryCov.push_back(-1.0f); continue; }
                    if (qline[0] == '>') qline = qline.substr(1);
                    TaxID taxid = -1; float cov = -1.0f;
                    if (regex_search(qline, m, accRe)) {
                        string acc = m[1];
                        size_t dot = acc.find('.');
                        if (dot != string::npos) acc = acc.substr(0, dot);
                        auto taxIt = assacc2taxid.find(acc);
                        if (taxIt != assacc2taxid.end()) taxid = taxIt->second;
                        auto covIt = assacc2cov.find(acc);
                        if (covIt != assacc2cov.end()) cov = covIt->second;
                    }
                    queryTaxid.push_back(taxid);
                    queryCov.push_back(cov);
                }
            }
            cerr << "Query names loaded: " << queryNameFiles[i]
                 << " (n=" << queryTaxid.size() << ")" << endl;

            // Phase 1: read the read-group file once into a flat per-coverage list.
            // covReads[cov] = vector of (taxid, groupId) — compact, no per-rank expansion yet.
            map<float, vector<pair<TaxID, int>>> covReads;
            {
                ifstream rgf(readGroupFiles[i]);
                if (!rgf.is_open()) { cerr << "Cannot open: " << readGroupFiles[i] << endl; continue; }
                string line;
                size_t badLines = 0;
                while (getline(rgf, line)) {
                    if (line.empty()) continue;
                    auto fields = Util::split(line, "\t");
                    if (fields.size() < 2) continue;
                    if (!isdigit((unsigned char)fields[0][0]) ||
                        !isdigit((unsigned char)fields[1][0])) {
                        if (++badLines <= 3)
                            cerr << "Warning: skipping non-numeric line in "
                                 << readGroupFiles[i] << ": " << line << endl;
                        continue;
                    }
                    int qIdx = stoi(fields[0]); // 1-based
                    int gId  = stoi(fields[1]);
                    if (gId == 0) continue;
                    size_t vecIdx = (size_t)(qIdx - 1);
                    if (vecIdx >= queryTaxid.size()) {
                        cerr << "Warning: qIdx " << qIdx << " out of range" << endl;
                        continue;
                    }
                    TaxID taxid = queryTaxid[vecIdx];
                    float cov   = queryCov[vecIdx];
                    if (taxid < 0 || cov < 0.0f) continue;
                    covReads[cov].emplace_back(taxid, gId);
                }
            }

            // Phase 2: for each rank, iterate over coverages, build group/tax maps on the fly.
            // Memory at any point: O(reads_at_one_cov) instead of O(ranks × covs × reads).
            for (const string &rank : ranks_local) {
                for (auto &kv : covReads) {
                    float cov = kv.first;
                    unordered_map<int, vector<TaxID>> group2taxs;
                    unordered_map<TaxID, vector<int>> tax2groups;
                    for (auto &entry : kv.second) {
                        TaxID taxAtRank = ncbiTax.getTaxIdAtRank(entry.first, rank);
                        group2taxs[entry.second].push_back(taxAtRank);
                        tax2groups[taxAtRank].push_back(entry.second);
                    }
                    computeGroupPurityRecallByCov(group2taxs, tax2groups,
                                                  results[i].cov2rank[cov][rank]);
                    results[i].cov2rank[cov][rank].calculate();
                }
            }

            cerr << readGroupFiles[i] << " done" << endl;
        }
    }

    // Summary table to stdout only
    std::set<float> allCovs;
    for (auto &r : results)
        for (auto &p : r.cov2rank)
            allCovs.insert(p.first);

    cout << "Coverage\tRank";
    for (size_t i = 0; i < results.size(); i++)
        cout << "\tPurity\tRecall\tF1";
    cout << "\n";

    for (float cov : allCovs) {
        for (const string &rank : ranks) {
            cout << cov << "\t" << rank;
            for (auto &r : results) {
                auto cIt = r.cov2rank.find(cov);
                if (cIt != r.cov2rank.end()) {
                    auto rIt = cIt->second.find(rank);
                    if (rIt != cIt->second.end()) {
                        cout << "\t" << rIt->second.purity
                             << "\t" << rIt->second.recall
                             << "\t" << rIt->second.f1;
                        continue;
                    }
                }
                cout << "\tN/A\tN/A\tN/A";
            }
            cout << "\n";
        }
    }
    cout.flush();
    return 0;
}
