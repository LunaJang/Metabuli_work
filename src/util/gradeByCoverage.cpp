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

static const float COV_BREAKS[] = {1, 5, 10, 20, 50};
static const char* COV_LABELS[] = {"<1", "1-5", "5-10", "10-20", "20-50", "50+"};
static const int N_COV_BINS = 6;

static int covToBin(float c) {
    for (int i = 0; i < N_COV_BINS - 1; i++)
        if (c < COV_BREAKS[i]) return i;
    return N_COV_BINS - 1;
}

struct GradeByCovResult {
    map<int, unordered_map<string, CovCount>> cov2rank;
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

    // No negative-value check on --readid-col / --taxid-col: their parameter regex is ^[0-9]+$,
    // so the argument parser rejects a negative before it reaches here.

    vector<string> ranks = par.testRank.empty()
        ? vector<string>{"class", "order", "family", "genus", "species"}
        : Util::split(par.testRank, ",");

    string names  = taxonomy + "/names.dmp";
    string nodes  = taxonomy + "/nodes.dmp";
    string merged = taxonomy + "/merged.dmp";

    // read id / accession (no version) -> coverage
    unordered_map<string, float> id2cov;
    {
        ifstream f(covMappingFile);
        if (!f.is_open()) { cerr << "Cannot open: " << covMappingFile << endl; exit(1); }
        string key; float cov;
        while (f >> key >> cov) {
            if (par.testType != "cami") {
                size_t dot = key.find('.');
                if (dot != string::npos) key = key.substr(0, dot);
            }
            id2cov[key] = cov;
        }
    }
    cerr << "Coverage map loaded: " << id2cov.size() << " entries" << endl;

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
    shared(results, ranks, N, mappingFiles, classFiles, id2cov, par, cerr, names, nodes, merged)
    {
        unordered_map<string, int> id2taxid;
        TaxonomyWrapper ncbiTax(names, nodes, merged, false);
        regex accRe("(GC[AF]_[0-9]+\\.[0-9]+)");
        smatch m;
        vector<string> ranks_local = ranks;
        cerr << "Taxonomy loaded" << endl;

#pragma omp for schedule(dynamic)
        for (size_t i = 0; i < N; ++i) {
            id2taxid.clear();

            {
                string k, v;
                ifstream mf(mappingFiles[i]);
                if (!mf.is_open()) { cerr << "Cannot open: " << mappingFiles[i] << endl; continue; }
                while (getline(mf, k, '\t') && getline(mf, v, '\n')) {
                    if (!v.empty() && v.back() == '\r') v.pop_back();
                    if (v.empty() || !isdigit((unsigned char)v[0])) continue;
                    if (par.testType == "cami") {
                        size_t slash = k.rfind('/');
                        if (slash != string::npos) k = k.substr(0, slash);
                    } else {
                        size_t d = k.find('.');
                        if (d != string::npos) k = k.substr(0, d);
                    }
                    id2taxid[k] = stoi(v);
                }
            }

            ifstream cf(classFiles[i]);
            if (!cf.is_open()) { cerr << "Cannot open: " << classFiles[i] << endl; continue; }
            string line;
            size_t classified = 0;
            while (getline(cf, line)) {
                if (line.empty() || line[0] == '#') continue;
                auto fields = Util::split(line, "\t");
                // Fatal, not skipped. Blank lines are already gone above, so a short line here
                // means the column indices do not match the file: every record would be skipped
                // and the run would report a score computed over nothing. Silently continuing is
                // what let a centrifuge run against the metabuli defaults produce output files
                // while `grade` on the same input segfaulted. exit rather than return because
                // this loop is inside the #pragma omp parallel above.
                if ((int) fields.size() <= max(par.readIdCol, par.taxidCol)) {
                    cerr << "Error: " << classFiles[i] << " has a line with " << fields.size()
                         << " tab-separated fields, but --readid-col " << par.readIdCol
                         << " --taxid-col " << par.taxidCol << " needs at least "
                         << (max(par.readIdCol, par.taxidCol) + 1) << "." << endl;
                    cerr << "       The line was: " << line << endl;
                    cerr << "       These are ZERO-based indices. Centrifuge writes"
                         << " readID/seqID/taxID, so it needs --readid-col 0 --taxid-col 2." << endl;
                    exit(EXIT_FAILURE);
                }
                if (fields[par.taxidCol].empty() || !isdigit(fields[par.taxidCol][0])) continue;

                string id = fields[par.readIdCol];
                int classInt = stoi(fields[par.taxidCol]);

                string lookupKey;
                if (par.testType == "cami") {
                    lookupKey = id;
                    size_t slash = lookupKey.rfind('/');
                    if (slash != string::npos) lookupKey = lookupKey.substr(0, slash);
                } else {
                    if (!regex_search(id, m, accRe)) continue;
                    string acc = m[1];
                    size_t dot = acc.find('.');
                    if (dot != string::npos) acc = acc.substr(0, dot);
                    lookupKey = acc;
                }

                if (classInt != 0) classified++;

                auto taxIt = id2taxid.find(lookupKey);
                if (taxIt == id2taxid.end()) continue;
                TaxID rightAnswer = taxIt->second;

                auto covIt = id2cov.find(lookupKey);
                if (covIt == id2cov.end()) continue;
                int bin = covToBin(covIt->second);

                for (const string &rank : ranks_local)
                    compareTaxonAtRank_ByCov(classInt, rightAnswer, ncbiTax,
                                             results[i].cov2rank[bin][rank], rank);
            }
            cf.close();

            for (auto &p : results[i].cov2rank)
                for (const string &rank : ranks_local)
                    p.second[rank].calculate();

            cerr << classFiles[i] << " done (classified=" << classified << ")" << endl;
        }
    }

    cout << "Coverage\tRank";
    for (size_t i = 0; i < results.size(); i++)
        cout << "\tPrecision\tSensitivity\tF1";
    cout << "\n";

    for (int bin = 0; bin < N_COV_BINS; bin++) {
        for (const string &rank : ranks) {
            cout << COV_LABELS[bin] << "\t" << rank;
            for (auto &r : results) {
                auto cIt = r.cov2rank.find(bin);
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
