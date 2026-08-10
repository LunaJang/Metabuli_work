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
    par.minEdgeWeight = 10;
    par.convergenceThreshold = 0.01f;
    // Absolute cap is the primary brake on Sum C(m,2): it bounds a single k-mer's edge
    // contribution regardless of read count. The ratio threshold scales with the dataset
    // (readCnt * ratio), so one value cannot fit both a 5k-read and a 62M-read run --
    // 0.0001 skipped everything on the former and 11 k-mers on the latter. Ratio is left
    // off by default and kept only as a secondary net.
    par.maxKmerFreqRatio = 0.0f;
    par.maxKmerReads = 1000; // provisional; C(1000,2) = 499,500 edges per k-mer

    // Phase 1 core threshold as a fraction of k-mers per read. Provisional: on the fixture
    // 0.4 lands near where knee detection used to. Pending gradeGroup measurement.
    par.minOverlapRatio = 0.4f;
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