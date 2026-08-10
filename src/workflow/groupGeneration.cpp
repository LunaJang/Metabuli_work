#include "GroupGenerator.h"
#include "Parameters.h"
#include "LocalParameters.h"
#include "FileUtil.h"
#include "common.h"

void setGroupGenerationDefaults(LocalParameters & par){    
    // Drops k-mers adjacent to a common-k-mer hit, not just the hit itself. On the fixture
    // this removes 17% of k-mers but only 0.7% of edges -- it lowers weights rather than
    // deleting pairs, so it is not a disk-volume lever. --min-overlap-ratio scales with the
    // surviving k-mer count and absorbs the shift; a fixed threshold would not.
    par.neighborKmers = 1;
    par.groupingIter = 10;
    // Phase 2 floor. Sweeping it over 10/5/1 on the species-inclusion benchmark left the useful
    // signal identical (reads landing in their species' dominant group: 0.08602/0.08602/0.08607);
    // it only trades coverage against concentration. 5 is picked for the query reduction it buys
    // (4.1x vs 3.2x at 10) without the purity loss seen at 1.
    par.minEdgeWeight = 5;
    par.convergenceThreshold = 0.01f;
    // Absolute cap is the primary brake on Sum C(m,2): it bounds a single k-mer's edge
    // contribution regardless of read count. The ratio threshold scales with the dataset
    // (readCnt * ratio), so one value cannot fit both a 5k-read and a 62M-read run --
    // 0.0001 skipped everything on the former and 11 k-mers on the latter. Ratio is left
    // off by default and kept only as a secondary net.
    par.maxKmerFreqRatio = 0.0f;
    par.maxKmerReads = 1000; // provisional; C(1000,2) = 499,500 edges per k-mer

    // Phase 1 core threshold as a fraction of k-mers per read. Measured on the species-inclusion
    // benchmark (61.7 M reads, 49.6 k-mers/read -> core threshold 15): 0.3 is the peak, with 0.2,
    // 0.4 and 0.5 all below it. Note the usable range starts above --min-edge / k-mers-per-read;
    // a ratio resolving to <= --min-edge is rejected and knee detection takes over instead.
    par.minOverlapRatio = 0.3f;
    // Phase 1.5 support threshold -- the one knob that moved the useful signal (+62% at a fixed
    // Phase 1 and 2), and it also flattens the sensitivity to the core threshold: without it,
    // core 25 -> 15 gained 73%; with it, 8%. 2 and 3 are within noise of each other; 2 leaves
    // 10% fewer groups.
    par.minSupport = 2;
    par.syncmer = 1;
    par.smerLen = 5;
    par.seqMode = 2;    
    par.verbosity = 3;
    par.ramUsage = 128;
    par.printLog = 0;
    par.maskMode = 0;
    par.maskProb = 0.9;
    par.matchPerKmer = 4; 
}

int groupGeneration(int argc, const char **argv, const Command& command)
{
    LocalParameters & par = LocalParameters::getLocalInstance();
    setGroupGenerationDefaults(par);
    par.parseParameters(argc, argv, command, true, Parameters::PARSE_ALLOW_EMPTY, 0);
    if (par.syncmer == 0) {
        par.kmerFormat = 3;
    } else {
        par.kmerFormat = 5;
    }

    if (par.seqMode == 2) {
        // Check if the second argument is a directory
        if (FileUtil::directoryExists(par.filenames[1].c_str())) {
            cerr << "Error: " << par.filenames[1] << " is a directory. Please specify a query file name." << endl;
            cerr << "       For '--seq-mode 2', please provide two query files." << endl;
            exit(1);
        }

        if (!FileUtil::directoryExists(par.filenames[5].c_str())) {
            FileUtil::makeDir(par.filenames[5].c_str());
        }
    } else {
        // Check if the second argument is file
        if (FileUtil::fileExists(par.filenames[1].c_str()) 
            && !FileUtil::directoryExists(par.filenames[1].c_str())) {
            cerr << "Error: " << par.filenames[1] << " is a file. Please specify a database directory." << endl;
            cerr << "       For '--seq-mode 1' and '--seq-mode 3', please provide one query file." << endl;
            exit(1);
        }

        if (!FileUtil::directoryExists(par.filenames[4].c_str())) {
            FileUtil::makeDir(par.filenames[4].c_str());
        }
    }

#ifdef OPENMP
    omp_set_num_threads(par.threads);
#endif    
    GroupGenerator * groupGenerator = new GroupGenerator(par);
    groupGenerator->startGroupGeneration(par);
    delete groupGenerator;
    return 0;
}