#include "GroupGenerator.h"
#include "Parameters.h"
#include "LocalParameters.h"
#include "FileUtil.h"
#include "common.h"
#include <fstream>
#include <string>

void setGroupGenerationDefaults(LocalParameters & par){
    // The reads-per-k-mer cap, chosen from the sample rather than fixed. A k-mer shared by m reads
    // contributes C(m,2) edges, so the frequency tail decides the edge volume, and the tail is not
    // the same shape on two datasets: at this quantile the cap resolves to 63 on species-inclusion,
    // 127 on species-exclusion and 255 on CAMI2 strain-madness. An absolute default cannot be
    // right on all three -- 1000, the value this replaces, was above every one of them and so cut
    // nothing where cutting mattered most.
    //
    // maxKmerReads is the absolute form and stays available, but it is off by default: an explicit
    // cap wins over the quantile, so leaving it set would make the quantile dead. Every published
    // measurement passed --max-kmer-reads 0 --max-kmer-quantile 0.995 explicitly for exactly that
    // reason; the defaults now say the same thing.
    par.maxKmerReads = 0;
    par.maxKmerQuantile = 0.995f;

    // Phase 1 core threshold as a fraction of k-mers per read.
    //
    // 0.5 since 2026-09-03, chosen for purity. The earlier 0.3 came from an F1 peak on the
    // species-inclusion benchmark; sweeping 0.30 / 0.40 / 0.50 / 0.60 at 203c83a3 on two datasets
    // says F1 is the wrong thing to read.
    //
    // CAMI2 strain-madness, species purity and Recall x c (c = grouped reads / all reads, because
    // gradeGroup's recall denominator excludes ungrouped reads):
    //
    //   rho    purity     Recall x c   core
    //   0.30   0.670474   0.3379       24
    //   0.40   0.939313   0.1827       33
    //   0.50   0.953477   0.0482       41
    //   0.60   0.958674   0.0274       49
    //
    // Purity saturates: +0.269, then +0.014, then +0.005. It never reaches the 0.98 bar, and past
    // 0.5 the remaining purity is bought at a steep price -- 0.5 -> 0.6 gains 0.005 of purity for
    // 43% of what is left of Recall x c. 0.5 is where that trade stops being worth taking.
    //
    // species-exclusion wants the opposite: purity is 0.999 at every rho, so the bar cannot
    // choose, and Recall x c falls monotonically (0.1679 / 0.1374 / 0.1209 / 0.1092). A simulated
    // benchmark with a per-read key and controlled coverage does not constrain rho; real
    // strain-level data does. The default follows the harder case.
    //
    // This is a default, not a bound: every script passes --min-overlap-ratio explicitly, so the
    // value here only decides what a bare `metabuli grouping` does. It is set to the operating
    // point the paper reports so those two cannot disagree -- the same failure --max-kmer-reads
    // and --max-kmer-quantile had until 2026-08-29.
    par.minOverlapRatio = 0.5f;
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
    // Intermediate partitions. 16 is the value the published numbers were produced at, and it is
    // a constant on purpose: 0 makes the routing follow --threads, so the same input on a machine
    // with a different core count would be split differently and the file and flush counts would
    // move with it. 0 is kept for reproducing runs made before this parameter existed and is not
    // a default any more. This function overrides LocalParameters' own initialisation, so both
    // places have to agree.
    par.partitions = 16;
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
    // Per-read group membership scores for a downstream EM. Off by default: the pass is cheap
    // relative to edge generation but the output is not -- top-k slots at 8 bytes over a
    // hundred-million-read run is tens of gigabytes, and nothing in this command reads it back.
    // This function overrides LocalParameters' own initialisation, so both places have to agree.
    par.score = 0;
    par.scoreTopK = 8;
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

// Everything the two grouping commands share. They differ only in the defaults their caller set,
// so the body lives here once rather than being copied into each entry point -- a copy would let a
// later fix land in one of them.
static int runGroupGeneration(int argc, const char **argv, const Command& command,
                              LocalParameters & par)
{
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

    // --print-log takes a different path through startGroupGeneration: it calls mergeGraph_one and
    // skips the threshold, all three phases and saveGroupsToFile. There are no groups on that
    // path, so there is nothing to score against. Stopping is better than writing an empty file
    // that looks like a result.
    if (par.score != 0 && par.printLog != 0) {
        cerr << "Error: --score and --print-log cannot be used together." << endl;
        cerr << "       --print-log stops after merging the graph; it produces no groups, and" << endl;
        cerr << "       --score scores reads against the groups Phase 1-3 produce." << endl;
        return 1;
    }
    // The score table is scoreTopK slots per read, held for the whole scoring pass. The upper
    // bound is a memory guard, not an algorithmic one: 64 slots over a 200 M-read run is 100 GB.
    if (par.score != 0 && (par.scoreTopK < 1 || par.scoreTopK > 64)) {
        cerr << "Error: --score-top-k must be in [1, 64] (given " << par.scoreTopK << ")." << endl;
        cerr << "       It is how many candidate groups each read keeps, at 8 bytes per slot" << endl;
        cerr << "       held in memory for every read at once." << endl;
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
        }
        // No kmer_params: the DB predates the record and there is nothing to compare against.
        // Silent on purpose. It warned on every run against every DB built before the check, so
        // the line said only "this DB is old" -- never that anything was wrong with this run.
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

// The fast, light variant. Only the edge set differs: one shared k-mer produces a hub and its
// m-1 spokes instead of all C(m,2) pairs, so the edge volume is linear in m rather than quadratic.
// Every other default is grouping's, on purpose -- the cap included. --max-kmer-quantile matters
// less once the quadratic term is gone, but it also raises purity (CAMI2 strain-madness: 0.531
// uncapped against 0.663 at 0.995), and a different value here would add a second variable to
// every comparison between the two commands.
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