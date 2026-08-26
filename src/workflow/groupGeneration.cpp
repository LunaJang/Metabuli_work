#include "GroupGenerator.h"
#include "Parameters.h"
#include "LocalParameters.h"
#include "FileUtil.h"
#include "common.h"
#include <fstream>
#include <string>

void setGroupGenerationDefaults(LocalParameters & par){
    // Absolute cap on reads per k-mer -- the primary brake on Sum C(m,2). It bounds a single
    // k-mer's edge contribution regardless of read count, which is why it is kept as a resource
    // knob rather than expressed as a fraction of the read count: a fraction cannot fit both a
    // 5k-read and a 62M-read run (0.0001 skipped everything on the former and 11 k-mers on the
    // latter, which is why --max-kmer-freq-ratio was dropped).
    par.maxKmerReads = 1000; // provisional; C(1000,2) = 499,500 edges per k-mer
    // Automatic selection of the cap above from the reads-per-k-mer distribution. Off for now:
    // --max-kmer-reads is set, and an explicit cap wins, so this changes nothing until the
    // default above is dropped to 0.
    par.maxKmerQuantile = 0.0f;

    // Phase 1 core threshold as a fraction of k-mers per read. Measured on the species-inclusion
    // benchmark (61.7 M reads, 49.6 k-mers/read -> core threshold 15): 0.3 is the peak, with 0.2,
    // 0.4 and 0.5 all below it.
    par.minOverlapRatio = 0.3f;
    // Weak-band lower bound as a fraction of the core threshold, also Phase 3's floor.
    // 5/15 = 0.3333 reproduces the species-inclusion operating point, where the band was (5, 15].
    // The absolute 5 it replaces meant different things per dataset: 5/15 = 0.333 of the core on
    // that benchmark but 5/34 = 0.147 on CAMI2 marine, so marine's band was three times as wide
    // in absolute terms and Phase 2 there absorbed far more chance links.
    // 0.3333 is what every result up to 2026-08-24 was measured with. It was briefly replaced by
    // an absolute --min-edge on the argument that the ratio reproduces an absolute 5 only where
    // the core is 15; measurement reversed that, because an absolute 10 is 0.30 of the core on
    // strain-madness (core 33) but 0.38 on species-exclusion (core 26), and on strain-madness the
    // difference between a bound of 11 and 10 is species purity 0.932 against 0.763.
    par.weakBandRatio = 0.3333f;
    // Intermediate partitions. 0 follows --threads, which is what every run before this parameter
    // existed did; 16 is the value the published numbers were produced at. Kept at 0 until the
    // route-loop and cap work below it is verified, so this build reproduces those runs exactly.
    par.partitions = 0;
    // Phase 2 support as a fraction of the smaller unit's read count. 0 keeps the pre-ratio
    // behaviour (count weak edges, floor 2), which is the measured operating point; the ratio is
    // opt-in until a value is measured on marine. Sweeping the old absolute support over 2/3 was
    // within noise, but it is not scale-free -- chance links between units A and B grow with
    // |A| * |B|, so a fixed count is met automatically once coverage is high.
    // Largest unit Phase 2 may merge, in reads. 0 keeps the pre-gate behaviour, which is what
    // every measurement so far used. Off by default because the bound that fits real metagenomes
    // (where an unbounded Phase 2 chains units into one component holding most of the reads)
    // is not the one that fits the simulated benchmark (where Phase 2 doubled recall at no
    // measurable purity cost). This function overrides LocalParameters' own initialisation, so
    // both places have to agree.
    // K-mers discarded either side of a common-k-mer hit. 0 discards the hit alone; 1 discards
    // three per hit and was the value every result before 2026-08-24 used, on the measurement that
    // it removes 17% of query k-mers but only 0.7% of edges. That measurement counts edges, and
    // edge count is not connectivity. On CAMI2 strain-madness, 1 -> 0 took k-mers per read from
    // 73.6 to 81.4, groups from 2,702,839 to 2,561,814 and species Recall*c from 0.164 to 0.245,
    // for 0.010 of purity. This function overrides LocalParameters' own initialisation, so both
    // places have to agree.
    par.commonKmerSpan = 0;
    // Peak disk the subGraph_* intermediates may hold at once. 0 derives it from the free space
    // at the output directory (80%). A safety ceiling, not a performance target -- a run that
    // fits under it never folds early, so machines with room pay nothing for it. This function
    // overrides LocalParameters' own initialisation, so both places have to agree.
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

    // Both algorithm thresholds are ratios and there is no absolute fallback for either, so a
    // value outside its usable range has to stop the run rather than silently pick something else.
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

    // The query k-mers have to be built the same way as the common k-mer DB, or the two k-mer
    // spaces do not line up and the common-k-mer filter silently removes nothing. DBs built
    // before this record existed cannot be checked; warn instead of failing on those.
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
        } else {
            cerr << "[WARN] " << paramFile << " is missing, so the DB's k-mer settings cannot be"
                 << " checked against --syncmer " << par.syncmer << " --smer-len " << par.smerLen
                 << "." << endl;
            cerr << "[WARN] The DB predates the check. A mismatch would make the common k-mer"
                 << " filter remove nothing, silently." << endl;
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