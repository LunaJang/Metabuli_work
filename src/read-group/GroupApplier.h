
#ifndef GROUP_APPLIER_H
#define GROUP_APPLIER_H

#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include "NcbiTaxonomy.h"
#include "LocalParameters.h"
#include "Taxonomer.h"
#include "Reporter.h"
#include <set>
#include <cassert>
#include <atomic>

using namespace std;

struct OrgResult {
    int label;
    float score;
    string name;
};

class GroupApplier {
protected:
    const LocalParameters & par;
    string groupFileDir;
    string groupmapFileDir;
    string taxDbDir;
    string orgRes;
    string outDir;
    
    // Agents    
    Reporter *reporter = nullptr;
    TaxonomyWrapper *taxonomy = nullptr;

    // Output
    string updatedResultFileName;
    string updatedReportFileName;

public:
    GroupApplier(LocalParameters & par);

    ~GroupApplier();

    void startGroupApplication(const LocalParameters & par);

    void loadOrgResult(vector<OrgResult>& orgResults, size_t& processedReadCnt);

    void loadGroupsFromFile(unordered_map<uint32_t, unordered_set<uint32_t>> &groupInfo,
                                        vector<uint32_t> &groupMappingInfo,
                                        const string &groupInfoFileName,
                                        const string &groupMappingInfoFileName);
    
    void getRepLabel(vector<OrgResult>& orgResults, 
                     const unordered_map<uint32_t, unordered_set<uint32_t>>& groupInfo, 
                     unordered_map<uint32_t, uint32_t>& repLabel,
                     std::unordered_map<int, int>& external2internalTaxId);
    
    void loadRepLabel(std::unordered_map<uint32_t, uint32_t>& repLabel);

    void applyRepLabel(const vector<OrgResult>& orgResults, 
                       const vector<uint32_t>& queryGroupInfo, 
                       const unordered_map<uint32_t, uint32_t>& repLabel,
                       std::unordered_map<int, int>& external2internalTaxId);
};


#endif // GROUP_APPLIER_H