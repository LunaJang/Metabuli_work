#include "GroupGenerator.h"
#include "Parameters.h"
#include "LocalParameters.h"
#include "FileUtil.h"
#include "common.h"
#include <fstream>
#include <string>

void setGroupGenerationDefaults(LocalParameters & par){
    par.maxKmerReads = 0;
    par.maxKmerQuantile = 0.995f;

    par.minOverlapRatio = 0.5f;
    par.weakBandRatio = 0.3333f;
    par.partitions = 16;
    par.commonKmerSpan = 0;
    par.maxTmpDiskMiB = 0;
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

static int runGroupGeneration(int argc, const char **argv, const Command& command,
                              LocalParameters & par)
{
    par.parseParameters(argc, argv, command, true, Parameters::PARSE_ALLOW_EMPTY, 0);
    if (par.syncmer == 0) {
        par.kmerFormat = 3;
    } else {
        par.kmerFormat = 5;
    }

    if (par.minOverlapRatio <= 0.0f) {
        cerr << "Error: --min-overlap-ratio must be > 0 (given " << par.minOverlapRatio << ")." << endl;
        cerr << "       The Phase 1 core threshold is derived from it as ratio x k-mers per read;" << endl;
        cerr << "       there is no absolute threshold to fall back to." << endl;
        return 1;
    }
    if (par.weakBandRatio <= 0.0f || par.weakBandRatio >= 1.0f) {
        cerr << "Error: --weak-band-ratio must be in (0, 1) (given " << par.weakBandRatio << ")." << endl;
        cerr << "       It is the weak band's lower bound as a fraction of the core threshold." << endl;
        cerr << "       At 0 the band would swallow every edge, including pairs sharing nothing;" << endl;
        cerr << "       at 1 it would be empty and the later passes would have nothing to use." << endl;
        return 1;
    }
    if (par.partitions < 0) {
        cerr << "Error: --partitions must be >= 0 (given " << par.partitions << ")." << endl;
        cerr << "       0 means follow --threads." << endl;
        return 1;
    }
    if (par.maxKmerQuantile < 0.0f || par.maxKmerQuantile > 1.0f) {
        cerr << "Error: --max-kmer-quantile must be in [0, 1] (given " << par.maxKmerQuantile << ")." << endl;
        cerr << "       It is the share of k-mers (counted over those in at least two reads)" << endl;
        cerr << "       that the reads-per-k-mer cap keeps. 0 disables the automatic cap." << endl;
        return 1;
    }

    {
        const std::string dbDir = par.filenames[1 + (par.seqMode == 2)];
        const std::string paramFile = dbDir + "/kmer_params";
        if (FileUtil::fileExists(paramFile.c_str())) {
            int dbSyncmer = -1, dbSmerLen = -1, dbKmerFormat = -1;
            std::ifstream in(paramFile);
            std::string key;
            int value = 0;
            while (in >> key >> value) {
                if (key == "syncmer") { dbSyncmer = value; }
                else if (key == "smer_len") { dbSmerLen = value; }
                else if (key == "kmer_format") { dbKmerFormat = value; }
            }
            if (dbSyncmer != par.syncmer || dbSmerLen != par.smerLen || dbKmerFormat != par.kmerFormat) {
                cerr << "Error: k-mer settings do not match the common k-mer DB at " << dbDir << "." << endl;
                cerr << "       DB:      --syncmer " << dbSyncmer << " --smer-len " << dbSmerLen
                     << " (k-mer format " << dbKmerFormat << ")" << endl;
                cerr << "       Request: --syncmer " << par.syncmer << " --smer-len " << par.smerLen
                     << " (k-mer format " << par.kmerFormat << ")" << endl;
                cerr << "       Rebuild the DB with these settings, or pass the DB's settings." << endl;
                return 1;
            }
        }
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

int groupGeneration(int argc, const char **argv, const Command& command)
{
    LocalParameters & par = LocalParameters::getLocalInstance();
    setGroupGenerationDefaults(par);
    return runGroupGeneration(argc, argv, command, par);
}

void setEasyGroupGenerationDefaults(LocalParameters & par){
    setGroupGenerationDefaults(par);
    par.edgeMode = 1; // EdgeMode::EDGE_MODE_STAR
}

int easyGroupGeneration(int argc, const char **argv, const Command& command)
{
    LocalParameters & par = LocalParameters::getLocalInstance();
    setEasyGroupGenerationDefaults(par);
    return runGroupGeneration(argc, argv, command, par);
}