
#ifndef GROUP_GENERATOR_H
#define GROUP_GENERATOR_H

#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include "IndexCreator.h"
#include "SeqIterator.h"
#include "NcbiTaxonomy.h"
#include "LocalParameters.h"
#include "QueryIndexer.h"
#include "ReducedKmerMatcher.h"
#include "KmerExtractor.h"
#include "Taxonomer.h"
#include "Reporter.h"
#include "KSeqWrapper.h"
#include "DeltaIdxReader.h"
#include <set>
#include <cassert>
#include <thread>
#include <atomic>

#define BufferSize 16'777'216 //16 * 1024 * 1024 // 16 M
using namespace std;

struct Relation {
    uint32_t id1;
    uint32_t id2;
    uint16_t weight;

    Relation(uint32_t a = 0, uint32_t b = 0): id1(a), id2(b), weight(0) {}
    Relation(uint32_t a, uint32_t b, uint16_t w): id1(a), id2(b), weight(w) {}

    static bool compare(const Relation& a, const Relation& b) {
        if (a.id1 != b.id1) return a.id1 < b.id1;
        return a.id2 < b.id2;
    }

    bool operator<(const Relation& other) const {
        if (id1 != other.id1) return id1 < other.id1;
        return id2 < other.id2;
    }

    bool operator==(const Relation& other) const {
        return id1 == other.id1 && id2 == other.id2;
    }
};

class DisjointSet {
public:
    std::vector<uint32_t> parent;
    std::vector<uint32_t> rank;
    std::vector<bool> grouped;

    DisjointSet(size_t numQuery) 
        : parent(numQuery), rank(numQuery), grouped(numQuery) {
        for (size_t i = 0; i < numQuery; ++i) {
            parent[i] = i;
            rank[i] = 0;
            grouped[i] = false;
        }
    }

    uint32_t find(uint32_t element) {
        if (parent[element] != element) {
            parent[element] = find(parent[element]);
        }
        return parent[element];
    }

    void unionSets(uint32_t elem1, uint32_t elem2) {
        uint32_t root1 = find(elem1);
        uint32_t root2 = find(elem2);

        grouped[elem1] = true;
        grouped[elem2] = true;

        if (root1 != root2) {
            if (rank[root1] < rank[root2]) {
                parent[root1] = root2;
            } else if (rank[root1] > rank[root2]) {
                parent[root2] = root1;
            } else {
                if (root1 < root2) {
                    parent[root2] = root1;
                    rank[root1]++;
                }
                else{
                    parent[root1] = root2;
                    rank[root2]++;
                }
            }
        }
    }

    void flatten() {
        for (size_t i = 0; i < parent.size(); ++i) {
            parent[i] = find(i);
        }
    }

    DisjointSet& operator+=(const DisjointSet& rhs) {
        for (size_t i = 0; i < parent.size(); ++i) {
            uint32_t p = rhs.parent[i];
            if (p != i) unionSets(i, p);
        }
        return *this;
    }

    friend DisjointSet operator+(DisjointSet lhs, const DisjointSet& rhs) {
        lhs += rhs;
        return lhs;
    }
};

struct OrgResult {
    int label;
    float score;
    string name;
};

class GroupGenerator {
protected:
    const LocalParameters & par;
    string commonKmerDB;
    string taxDbDir;
    string orgRes;
    string outDir;
    size_t matchPerKmer;
    int kmerFormat;
    
    // Agents    
    GeneticCode * geneticCode = nullptr;
    QueryIndexer *queryIndexer = nullptr;
    KmerExtractor *kmerExtractor = nullptr;
    Reporter *reporter = nullptr;
    TaxonomyWrapper *taxonomy = nullptr;

    // Output
    string updatedResultFileName;
    string updatedReportFileName;
    unordered_map<TaxID, unsigned int> taxCounts;
    
    unordered_map<TaxID, TaxID> taxId2speciesId;
    unordered_map<TaxID, TaxID> taxId2genusId;
    size_t numOfSplits = 0;
    size_t numOfGraph = 0;
    std::vector<uint64_t> kmerBoundaries;
    bool boundariesInitialized = false;
    bool useOnlyTrueRelations = false; // for debug

public:
    GroupGenerator(LocalParameters & par);

    void startGroupGeneration(const LocalParameters & par);
    
    void filterCommonKmers(Buffer<Kmer>& queryKmerBuffer,
                           Buffer<std::pair<uint32_t, uint32_t>> & matchBuffer,
                           const string & db="");

    void writeKmers(Buffer<Kmer>& queryKmerBuffer, 
                    size_t processedReadCnt);

    std::vector<std::pair<size_t, size_t>> getKmerRanges(const Buffer<Kmer>& kmerBuffer, 
                                                         size_t offset);

    void loadOrgResult(vector<OrgResult>& orgResults);

    void makeSubGraph(size_t processedReadCnt);
    
    void saveSubGraphToFile(const unordered_map<uint64_t, uint16_t>& pair2weight,
                            const size_t counter_now);

    void mergeGraph(size_t processedReadCnt);

    void makeGroups(int groupKmerThr,
                    size_t processedReadCnt,
                    unordered_map<uint32_t, unordered_set<uint32_t>>& groupInfo, 
                    vector<uint32_t> &queryGroupInfo);

    void getRepLabel(vector<OrgResult>& orgResults, 
                     const unordered_map<uint32_t, unordered_set<uint32_t>>& groupInfo, 
                     unordered_map<uint32_t, uint32_t>& repLabel,
                     std::unordered_map<int, int>& external2internalTaxId);
    
    void loadRepLabel(std::unordered_map<uint32_t, uint32_t>& repLabel);

    void applyRepLabel(const vector<OrgResult>& orgResults, 
                       const vector<uint32_t>& queryGroupInfo, 
                       const unordered_map<uint32_t, uint32_t>& repLabel,
                       std::unordered_map<int, int>& external2internalTaxId);

    ~GroupGenerator();
};


#endif // GROUP_GENERATOR_H