#ifndef METABULI_REPORTER_H
#define METABULI_REPORTER_H
#include "common.h"
#include "iostream"
#include "fstream"
#include <unordered_map>
#include "NcbiTaxonomy.h"
#include "LocalParameters.h"

using namespace std;


class Reporter {
private:
    string outDir;
    string jobId;
    NcbiTaxonomy * taxonomy;

    // Output
    string reportFileName;
    string readClassificationFileName;
    ofstream readClassificationFile;

public:
    Reporter(const LocalParameters &par, NcbiTaxonomy *taxonomy);
    // Write report
    void writeReportFile(int numOfQuery, unordered_map<TaxID, unsigned int> &taxCnt, bool krona = true);
    void writeReport(FILE *FP, const std::unordered_map<TaxID, TaxonCounts> &cladeCounts,
                     unsigned long totalReads, TaxID taxID = 0, int depth = 0);

    // Read by read classification results
    void openReadClassificationFile();
    void writeReadClassification(const vector<Query> & queryList, bool classifiedOnly = false);
    void closeReadClassificationFile();

    unsigned int cladeCountVal(const std::unordered_map<TaxID, TaxonCounts> &map, TaxID key);

    // Setter
    void setReportFileName(const string &reportFileName) {
        Reporter::reportFileName = reportFileName;
    }

    void setReadClassificationFileName(const string &readClassificationFileName) {
        Reporter::readClassificationFileName = readClassificationFileName;
    }
};


#endif //METABULI_REPORTER_H