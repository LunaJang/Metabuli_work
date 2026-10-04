#include "GroupApplier.h"
#include "Parameters.h"
#include "LocalParameters.h"
#include "FileUtil.h"
#include "common.h"

void setGroupApplicationDefaults(LocalParameters & par){    
    par.seqMode = 2;    
    par.ramUsage = 128;
    par.scoreCol = 5;
    par.readIdCol = 2;
    par.taxidCol = 3;
    par.weightMode = 1; // 0: uniform, 1: score, 2: score^2
    // 0, not 0.15. The threshold does two things at once and only one of them was
    // wanted. It filters the vote that picks a group's label, and it decides who that
    // label is written to: a member below it is overwritten even though it has a label
    // of its own. At 0.15 on Kraken2 that cost 36 M species labels on CAMI2
    // plant-associated and 15 M on strain-madness, because a group's LCA is often
    // coarser than the label it replaced, and the arm came out below both uniform
    // weighting and no propagation at all.
    //
    // At 0 the score still weights the vote and nothing else: only members with no
    // label are written to, so propagation cannot lower the number of classified reads
    // at any rank. A threshold remains available for whoever wants one, and it means
    // what it says rather than two things.
    //
    // 0.15 was never a general value either. It is the short-read threshold of
    // Metabuli-P, derived in that paper from where its own true and false positives
    // separate. Kraken2's manual defines a confidence but recommends no threshold and
    // defaults to 0; Centrifuger's score is an unbounded sum of squared hit lengths
    // whose authors say outright that estimating confidence from it is future work.
    // One number could not have meant the same thing in all three.
    par.minVoteScr = 0.0;
}

int groupApplication(int argc, const char **argv, const Command& command)
{
    LocalParameters & par = LocalParameters::getLocalInstance();
    setGroupApplicationDefaults(par);
    par.parseParameters(argc, argv, command, true, Parameters::PARSE_ALLOW_EMPTY, 0);

    if (par.weightMode != 0) {
        cout << "Warning: --weight-mode " << par.weightMode << " requires classification scores." << endl;
        cout << "         Make sure that score column is correctly set using --score-col." << endl;
    }

    // apply-group takes exactly five arguments -- group result, group mapping result, taxonomy
    // directory, read-by-read result, output directory -- and GroupApplier's constructor reads
    // filenames[0..4]. The output directory is therefore always index 4.
    //
    // This used to branch on --seq-mode and, for the default of 2, index filenames[5] and check
    // filenames[1] for being a directory. Both were wrong: there is no sixth argument, so
    // filenames[5] was an unchecked std::vector::operator[] past the end -- undefined behaviour
    // that happened not to crash -- and filenames[1] is the group mapping FILE, so the check
    // rejected the correct input rather than a wrong one. The branch was copied from
    // groupGeneration, where --seq-mode does shift the argument positions. apply-group never
    // reads par.seqMode at all (nothing in GroupApplier references it).
    if (!FileUtil::directoryExists(par.filenames[4].c_str())) {
        FileUtil::makeDir(par.filenames[4].c_str());
    }

#ifdef OPENMP
    omp_set_num_threads(par.threads);
#endif    
    GroupApplier * groupApplier = new GroupApplier(par);
    groupApplier->startGroupApplication(par);
    delete groupApplier;
    return 0;
}