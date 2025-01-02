#include "GroupGenerator.h"
#include "Parameters.h"
#include "LocalParameters.h"
#include "FileUtil.h"
#include "common.h"

void setApplyBinningDefaults(LocalParameters & par){
    par.voteMode = 0;
    par.majorityThr = 0.5;
}

int applyBinning(int argc, const char **argv, const Command& command)
{
    LocalParameters & par = LocalParameters::getLocalInstance();
    setApplyBinningDefaults(par);
    par.parseParameters(argc, argv, command, true, Parameters::PARSE_ALLOW_EMPTY, 0);

    // Check if the second argument is file
    if (FileUtil::fileExists(par.filenames[0].c_str()) 
        && !FileUtil::directoryExists(par.filenames[0].c_str())) {
        cerr << "Error: " << par.filenames[0] << " is a file. Please specify a database directory." << endl;
        exit(1);
    }

    if (!FileUtil::directoryExists(par.filenames[2].c_str())) {
        FileUtil::makeDir(par.filenames[2].c_str());
    }


#ifdef OPENMP
    omp_set_num_threads(par.threads);
#endif

    GroupGenerator * groupGenerator = new GroupGenerator(par);
    groupGenerator->startGroupGeneration(par);
    delete groupGenerator;
    return 0;
}