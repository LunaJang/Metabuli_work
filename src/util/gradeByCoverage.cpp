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
#include <cstdint>

using namespace std;

struct CovCount {
    int total = 0, TP = 0, FP = 0, FN = 0;
    float precision = 0.0f, sensitivity = 0.0f, f1 = 0.0f;
    void calculate() {
        precision   = (TP + FP > 0) ? (float)TP / (float)(TP + FP) : 0.0f;
        sensitivity = (total  > 0)  ? (float)TP / (float)total      : 0.0f;
        f1 = (precision + sensitivity > 0.0f)
             ? 2.0f * precision * sensitivity / (precision + sensitivity)
             : 0.0f;
    }
};

struct GradeByCovResult {
    map<float, unordered_map<string, CovCount>> cov2rank;
    string path;
};

static char compareTaxonAtRank_ByCov(TaxID shot, TaxID target,
                                      const TaxonomyWrapper &ncbiTax,
                                      CovCount &count, const string &rank) {
    TaxID targetAtRank = ncbiTax.getTaxIdAtRank(target, rank);
    const TaxonNode *targetNode = ncbiTax.taxonNode(targetAtRank);
    int rankIdx = ncbiTax.findRankIndex2(rank);
    if (ncbiTax.findRankIndex2(ncbiTax.getString(targetNode->rankIdx)) > rankIdx)
        return '-';

    if (shot == 1 || shot == 0) { count.FN++; count.total++; return 'N'; }

    TaxID shotAtRank = ncbiTax.getTaxIdAtRank(shot, rank);
    const TaxonNode *shotNode = ncbiTax.taxonNode(shotAtRank);
    if (ncbiTax.findRankIndex2(ncbiTax.getString(shotNode->rankIdx)) > rankIdx) {
        count.FN++; count.total++; return 'N';
    }

    count.total++;
    if (shotAtRank == targetAtRank) { count.TP++; return 'O'; }
    else                            { count.FP++; return 'X'; }
}

void setGradeByCoverageDefault(LocalParameters &par) {
    par.readIdCol = 1;
    par.taxidCol  = 2;
    par.verbosity = 2;
    par.testRank  = "";
    par.testType  = "gtdb";
}

int gradeByCoverage(int argc, const char **argv, const Command &command) {
    LocalParameters &par = LocalParameters::getLocalInstance();
    setGradeByCoverageDefault(par);
    par.parseParameters(argc, argv, command, false, Parameters::PARSE_ALLOW_EMPTY, 0);

    const string classListPath   = par.filenames[0];
    const string mappingListPath = par.filenames[1];
    const string covMappingFile  = par.filenames[2]; // accession<whitespace>coverage
    const string taxonomy        = par.filenames[3];

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
    vector<string> classFiles   = loadList(classListPath);
    vector<string> mappingFiles = loadList(mappingListPath);

    size_t N = classFiles.size();
    vector<GradeByCovResult> results(N);

#ifdef OPENMP
    omp_set_num_threads(par.threads);
#endif

#pragma omp parallel default(none), \
    shared(results, ranks, N, mappingFiles, classFiles, assacc2cov, par, cerr, names, nodes, merged)
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

            ifstream cf(classFiles[i]);
            if (!cf.is_open()) { cerr << "Cannot open: " << classFiles[i] << endl; continue; }
            string line;
            size_t classified = 0;
            while (getline(cf, line)) {
                if (line.empty() || line[0] == '#') continue;
                auto fields = Util::split(line, "\t");
                if ((int)fields.size() <= max(par.readIdCol, par.taxidCol)) continue;
                if (!isdigit(fields[par.taxidCol][0])) continue;

                string id = fields[par.readIdCol];
                if (!regex_search(id, m, accRe)) continue;
                string acc = m[1];
                size_t dot = acc.find('.');
                if (dot != string::npos) acc = acc.substr(0, dot);

                int classInt = stoi(fields[par.taxidCol]);
                if (classInt != 0) classified++;

                auto taxIt = assacc2taxid.find(acc);
                if (taxIt == assacc2taxid.end()) continue;
                TaxID rightAnswer = taxIt->second;

                auto covIt = assacc2cov.find(acc);
                if (covIt == assacc2cov.end()) continue;
                float cov = covIt->second;

                for (const string &rank : ranks_local)
                    compareTaxonAtRank_ByCov(classInt, rightAnswer, ncbiTax,
                                             results[i].cov2rank[cov][rank], rank);
            }
            cf.close();

            for (auto &p : results[i].cov2rank)
                for (const string &rank : ranks_local)
                    p.second[rank].calculate();

            cerr << classFiles[i] << " done (classified=" << classified << ")" << endl;
        }
    }

    // Summary table to stdout only
    std::set<float> allCovs;
    for (auto &r : results)
        for (auto &p : r.cov2rank)
            allCovs.insert(p.first);

    cout << "Coverage\tRank";
    for (size_t i = 0; i < results.size(); i++)
        cout << "\tPrecision\tSensitivity\tF1";
    cout << "\n";

    for (float cov : allCovs) {
        for (const string &rank : ranks) {
            cout << cov << "\t" << rank;
            for (auto &r : results) {
                auto cIt = r.cov2rank.find(cov);
                if (cIt != r.cov2rank.end()) {
                    auto rIt = cIt->second.find(rank);
                    if (rIt != cIt->second.end()) {
                        cout << "\t" << rIt->second.precision
                             << "\t" << rIt->second.sensitivity
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
