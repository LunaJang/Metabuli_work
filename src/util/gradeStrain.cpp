#include "Parameters.h"
#include "LocalParameters.h"
#include "Util.h"

#include <string>
#include <iostream>
#include <fstream>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <cstdint>

using namespace std;

// Strain-level grading. gradeGroup stops at species rank because it resolves labels through the
// taxonomy; CAMI2 carries genome_id per read, which is strain identity, and on strain-madness that
// is 408 labels against 26 tax ids. Scoring the same groups on both label sets with one formula is
// what separates "the groups mix species" from "the groups cannot split strains".
//
// Deliberately different from gradeGroup in three places:
//   - no taxonomy is loaded; a label is an opaque string, and no rank is derived from it
//   - answer keys are matched whole. gradeGroup truncates them at the first '.' to drop a GTDB
//     accession version, which is not a rule that applies to read names
//   - reads whose name is not in the answer are counted and reported. gradeGroup drops them
//     silently, so a join that fails completely still prints a full-looking table

struct Score {
    double purity = 0.0;
    double recall = 0.0;
    double f1() const {
        return (purity + recall > 0.0) ? 2.0 * purity * recall / (purity + recall) : 0.0;
    }
};

// Dense codes for the label strings, so the per-read vectors stay int-sized. A 166 M-read sample
// keeps two of them, which is 1.3 GB; holding the strings per read instead would be far worse.
class LabelTable {
public:
    int32_t code(const string &label) {
        auto found = codes.find(label);
        if (found != codes.end()) { return found->second; }
        const int32_t next = static_cast<int32_t>(codes.size());
        codes.emplace(label, next);
        return next;
    }
    size_t size() const { return codes.size(); }

private:
    unordered_map<string, int32_t> codes;
};

static vector<string> readList(const string &listFile, const string &what) {
    vector<string> names;
    ifstream list(listFile);
    if (!list.is_open()) {
        cerr << "Cannot open file for " << what << ": " << listFile << endl;
        return names;
    }
    string line;
    while (getline(list, line)) {
        if (!line.empty()) { names.push_back(line); }
    }
    return names;
}

// answer file: <read name>\t<label>. The label stays a string, so no stoi here.
static void loadAnswer(const string &path,
                       unordered_map<string, int32_t> &name2code,
                       LabelTable &labels) {
    ifstream answer(path);
    if (!answer.is_open()) {
        cerr << "Cannot open file for answer: " << path << endl;
        return;
    }
    string key, value;
    while (getline(answer, key, '\t')) {
        if (!getline(answer, value, '\n')) { break; }
        if (key.empty() || value.empty()) { continue; }
        name2code[key] = labels.code(value);
    }
}

// Same key derivation gradeGroup uses for query names, minus the GTDB accession regex: the CAMI2
// read name carries a "/1" mate suffix that the answer does not.
static string queryKey(const string &line, const string &testType) {
    string id = line;
    if (!id.empty() && id[0] == '>') { id = id.substr(1); }
    if (testType == "cami" || testType == "cami-long" || testType == "cami-euk") {
        const size_t pos = id.find('/');
        if (pos != string::npos) { id = id.substr(0, pos); }
    } else if (testType == "hiv" || testType == "hiv-ex") {
        const size_t pos = id.find('_');
        if (pos != string::npos) { id = id.substr(0, pos); }
    }
    return id;
}

// purity over groups, recall over labels -- the definitions gradeGroupByCoverage.cpp uses, so the
// tax_id row here has to reproduce gradeGroup's species row on the same run.
static Score score(const unordered_map<int, vector<int32_t>> &group2labels,
                   const unordered_map<int32_t, vector<int>> &label2groups) {
    Score out;

    double major = 0.0, grouped = 0.0;
    unordered_map<int32_t, size_t> freq;
    for (const auto &entry : group2labels) {
        freq.clear();
        for (const int32_t label : entry.second) { freq[label]++; }
        size_t best = 0;
        for (const auto &kv : freq) { best = max(best, kv.second); }
        major += static_cast<double>(best);
        grouped += static_cast<double>(entry.second.size());
    }
    out.purity = (grouped > 0.0) ? major / grouped : 0.0;

    double majorGroup = 0.0, labelReads = 0.0;
    unordered_map<int, size_t> groupFreq;
    for (const auto &entry : label2groups) {
        groupFreq.clear();
        for (const int group : entry.second) { groupFreq[group]++; }
        size_t best = 0;
        for (const auto &kv : groupFreq) { best = max(best, kv.second); }
        majorGroup += static_cast<double>(best);
        labelReads += static_cast<double>(entry.second.size());
    }
    out.recall = (labelReads > 0.0) ? majorGroup / labelReads : 0.0;

    return out;
}

void setGradeStrainDefault(LocalParameters &par) {
    par.testRank = "";
    par.testType = "gtdb";
}

int gradeStrain(int argc, const char **argv, const Command &command) {
    LocalParameters &par = LocalParameters::getLocalInstance();
    setGradeStrainDefault(par);
    par.parseParameters(argc, argv, command, false, Parameters::PARSE_ALLOW_EMPTY, 0);

    const vector<string> groupFiles     = readList(par.filenames[0], "group file list");
    const vector<string> readGroupFiles = readList(par.filenames[1], "read-group file list");
    const vector<string> strainFiles    = readList(par.filenames[2], "strain answer list");
    const vector<string> taxFiles       = readList(par.filenames[3], "tax answer list");
    const vector<string> queryNameFiles = readList(par.filenames[4], "query name file list");

    const size_t numberOfFiles = groupFiles.size();
    if (numberOfFiles == 0) {
        cerr << "No grouping result listed in " << par.filenames[0] << endl;
        return 1;
    }
    if (readGroupFiles.size() != numberOfFiles || strainFiles.size() != numberOfFiles
        || taxFiles.size() != numberOfFiles || queryNameFiles.size() != numberOfFiles) {
        cerr << "The five input lists must have the same number of entries: "
             << numberOfFiles << " groups, " << readGroupFiles.size() << " read-groups, "
             << strainFiles.size() << " strain answers, " << taxFiles.size() << " tax answers, "
             << queryNameFiles.size() << " query names." << endl;
        return 1;
    }

    // Serial on purpose: the run scripts always pass one file, so parallelising over the list
    // would add a shared-state hazard for no gain.
    for (size_t i = 0; i < numberOfFiles; ++i) {
        LabelTable strainLabels, taxLabels;
        unordered_map<string, int32_t> name2strain, name2tax;
        loadAnswer(strainFiles[i], name2strain, strainLabels);
        loadAnswer(taxFiles[i], name2tax, taxLabels);
        cout << "Strain answer loaded: " << strainFiles[i]
             << " (" << name2strain.size() << " reads, " << strainLabels.size() << " genomes)" << endl;
        cout << "Tax answer loaded:    " << taxFiles[i]
             << " (" << name2tax.size() << " reads, " << taxLabels.size() << " tax ids)" << endl;

        // Per query index, 0-based here and 1-based in the grouping files. -1 means the name was
        // not in that answer.
        vector<int32_t> queryStrain, queryTax;
        size_t unmatchedStrain = 0, unmatchedTax = 0;
        {
            ifstream names(queryNameFiles[i]);
            if (!names.is_open()) {
                cerr << "Cannot open file for query names: " << queryNameFiles[i] << endl;
                return 1;
            }
            string line;
            while (getline(names, line)) {
                const string key = queryKey(line, par.testType);
                const auto strainIt = name2strain.find(key);
                const auto taxIt = name2tax.find(key);
                if (strainIt == name2strain.end()) {
                    queryStrain.push_back(-1);
                    unmatchedStrain++;
                } else {
                    queryStrain.push_back(strainIt->second);
                }
                if (taxIt == name2tax.end()) {
                    queryTax.push_back(-1);
                    unmatchedTax++;
                } else {
                    queryTax.push_back(taxIt->second);
                }
            }
        }
        cout << "Query names loaded: " << queryNameFiles[i]
             << " (size=" << queryStrain.size() << ")" << endl;

        // Purity side, from the groups file, exactly as gradeGroup reads it.
        unordered_map<int, vector<int32_t>> group2strain, group2tax;
        size_t numGroup = 0, numReadsInGroup = 0;
        {
            ifstream groups(groupFiles[i]);
            if (!groups.is_open()) {
                cerr << "Cannot open file for groups: " << groupFiles[i] << endl;
                return 1;
            }
            string line;
            while (getline(groups, line)) {
                if (line.empty()) { continue; }
                const vector<string> fields = Util::split(line, "\t");
                if (fields.size() < 2) { continue; }
                const int groupId = stoi(fields[0]);
                for (size_t j = 1; j < fields.size(); ++j) {
                    const long qIdx = stol(fields[j]); // 1-based
                    if (qIdx <= 0 || static_cast<size_t>(qIdx) > queryStrain.size()) {
                        cerr << "Warning: query index " << qIdx << " out of range in "
                             << groupFiles[i] << endl;
                        continue;
                    }
                    const size_t vecIdx = static_cast<size_t>(qIdx) - 1;
                    if (queryStrain[vecIdx] < 0) { continue; }
                    group2strain[groupId].push_back(queryStrain[vecIdx]);
                    group2tax[groupId].push_back(queryTax[vecIdx]);
                    numReadsInGroup++;
                }
                numGroup++;
            }
        }

        // Recall side, from the read-group file, exactly as gradeGroup reads it.
        unordered_map<int32_t, vector<int>> strain2groups, tax2groups;
        size_t numReads = 0;
        {
            ifstream readGroups(readGroupFiles[i]);
            if (!readGroups.is_open()) {
                cerr << "Cannot open file for read-groups: " << readGroupFiles[i] << endl;
                return 1;
            }
            string line;
            while (getline(readGroups, line)) {
                if (line.empty()) { continue; }
                const vector<string> fields = Util::split(line, "\t");
                if (fields.size() < 2) { continue; }
                const long qIdx = stol(fields[0]); // 1-based
                const int groupId = stoi(fields[1]);
                numReads++;
                if (groupId == 0) { continue; }
                if (qIdx <= 0 || static_cast<size_t>(qIdx) > queryStrain.size()) {
                    cerr << "Warning: query index " << qIdx << " out of range in "
                         << readGroupFiles[i] << endl;
                    continue;
                }
                const size_t vecIdx = static_cast<size_t>(qIdx) - 1;
                if (queryStrain[vecIdx] >= 0) {
                    strain2groups[queryStrain[vecIdx]].push_back(groupId);
                }
                if (queryTax[vecIdx] >= 0) {
                    tax2groups[queryTax[vecIdx]].push_back(groupId);
                }
            }
        }

        const Score strainScore = score(group2strain, strain2groups);
        const Score taxScore = score(group2tax, tax2groups);

        // Strain purity restricted to groups that hold exactly one tax id. A group that is impure
        // at strain level but pure at species level is the interesting case: the method separated
        // species and then failed to separate the strains inside one. Mixing the two failures into
        // a single number hides which one the data shows.
        double pureMajor = 0.0, pureReads = 0.0;
        size_t pureGroups = 0;
        unordered_map<int32_t, size_t> freq;
        unordered_map<size_t, size_t> genomesPerGroup;
        size_t singleGenomeGroups = 0;
        for (const auto &entry : group2strain) {
            freq.clear();
            for (const int32_t label : entry.second) { freq[label]++; }
            genomesPerGroup[freq.size()]++;
            if (freq.size() == 1) { singleGenomeGroups++; }

            const auto taxIt = group2tax.find(entry.first);
            if (taxIt == group2tax.end() || taxIt->second.empty()) { continue; }
            bool speciesPure = true;
            const int32_t firstTax = taxIt->second.front();
            for (const int32_t t : taxIt->second) {
                if (t != firstTax) { speciesPure = false; break; }
            }
            if (!speciesPure) { continue; }
            pureGroups++;
            size_t best = 0;
            for (const auto &kv : freq) { best = max(best, kv.second); }
            pureMajor += static_cast<double>(best);
            pureReads += static_cast<double>(entry.second.size());
        }

        cout << readGroupFiles[i] << endl;
        cout << "The number of reads: " << numReads << endl;
        cout << "The number of groups: " << numGroup << endl;
        cout << "The number of reads in groups: " << numReadsInGroup << endl;
        cout << "Unmatched query names (strain answer): " << unmatchedStrain << endl;
        cout << "Unmatched query names (tax answer): " << unmatchedTax << endl;
        if (!queryStrain.empty()
            && unmatchedStrain * 100 > queryStrain.size()) { // more than 1%
            cout << "WARNING: over 1% of query names are not in the strain answer. Check that the "
                    "answer and the reads come from the same dataset, and that --test-type matches."
                 << endl;
        }
        cout << endl;

        cout << "Label\tPurity\tRecall\tF1" << endl;
        cout << "genome_id\t" << strainScore.purity << "\t" << strainScore.recall << "\t"
             << strainScore.f1() << endl;
        cout << "tax_id\t" << taxScore.purity << "\t" << taxScore.recall << "\t"
             << taxScore.f1() << endl;
        cout << endl;

        const double pureStrainPurity = (pureReads > 0.0) ? pureMajor / pureReads : 0.0;
        cout << "Strain purity within species-pure groups: " << pureStrainPurity
             << " (over " << static_cast<size_t>(pureReads) << " reads in " << pureGroups
             << " groups)" << endl;

        cout << "distinct genomes per group:";
        size_t bucket1 = 0, bucket2 = 0, bucket3to5 = 0, bucket6to10 = 0, bucketOver10 = 0;
        for (const auto &kv : genomesPerGroup) {
            if (kv.first == 1) { bucket1 += kv.second; }
            else if (kv.first == 2) { bucket2 += kv.second; }
            else if (kv.first <= 5) { bucket3to5 += kv.second; }
            else if (kv.first <= 10) { bucket6to10 += kv.second; }
            else { bucketOver10 += kv.second; }
        }
        cout << " 1: " << bucket1 << " | 2: " << bucket2 << " | 3-5: " << bucket3to5
             << " | 6-10: " << bucket6to10 << " | >10: " << bucketOver10 << endl;

        const size_t scoredGroups = group2strain.size();
        cout << "single-genome groups: " << singleGenomeGroups << "/" << scoredGroups;
        if (scoredGroups > 0) {
            cout << " (" << (100.0 * static_cast<double>(singleGenomeGroups)
                             / static_cast<double>(scoredGroups)) << "%)";
        }
        cout << endl;
    }

    return 0;
}
