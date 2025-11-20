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
    par.minVoteScr = 0.15;
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
    GroupApplier * groupApplier = new GroupApplier(par);
    groupApplier->startGroupApplication(par);
    delete groupApplier;
    return 0;
}