
#ifndef BINNING_H
#define BINNING_H

#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include "NcbiTaxonomy.h"
#include "LocalParameters.h"
#include "Taxonomer.h"
#include <set>
#include <cassert>

#define BufferSize 16'777'216 //16 * 1024 * 1024 // 16 M
using namespace std;

class Binning {
protected:
    string dbDir;
    
    // Agents
    NcbiTaxonomy *taxonomy;

public:
    Binning(LocalParameters & par);

    ~Binning();

    void startBinning(const LocalParameters & par);

    void makeGroupsFromBinning(const string &groupFileDir, 
                                        unordered_map<uint32_t, unordered_set<string>> &binningReadGroupInfo, 
                                        unordered_map<string, uint32_t> &queryGroupInfo);

    void loadMetabuliResult(const string &resultFileDir, 
                                     vector<pair<string, float>> &binnningMetabuliResult);

    void getRepLabel(const string &groupRepFileDir, 
                              vector<pair<int, float>> &binnningMetabuliResult,
                              unordered_map<uint32_t, unordered_set<string>> &binningReadGroupInfo, 
                              unordered_map<uint32_t, int> &repLabel, 
                              const string &jobId, 
                              int voteMode, 
                              float majorityThr);

    void applyRepLabel(const string &originalResultFileDir, 
                                unordered_map<string, uint32_t> queryGroupInfo, 
                                const unordered_map<uint32_t, int> &repLabel, 
                                const string &jobId);

};

#endif // BINNING_H