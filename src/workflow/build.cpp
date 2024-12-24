#include "IndexCreator.h"
#include "FileMerger.h"
#include "LocalParameters.h"
#include <Command.h>
#include "FileUtil.h"

void setDefaults_build(LocalParameters & par){
    par.skipRedundancy = 1;
    par.reducedAA = 0;
    par.ramUsage = 128;
    // par.spaceMask = "11111111";
    par.taxonomyPath = "" ;
    par.splitNum = 4096;
    par.maskProb = 0.9;
    par.maskMode = 1;
    par.accessionLevel = 0;
    // Get current date
    time_t now = time(0);
    tm *ltm = localtime(&now);
    par.dbDate = to_string(1900 + ltm->tm_year) + "-" + to_string(1 + ltm->tm_mon) + "-" + to_string(ltm->tm_mday);
    
    // Get random alphanumeric string fore dbName from current time
    srand(time(NULL));
    string randStr = to_string(rand());
    par.dbName = randStr.substr(0, 32);
}

int build(int argc, const char **argv, const Command &command){
    // Load parameters
    LocalParameters &par = LocalParameters::getLocalInstance();
    setDefaults_build(par);
    par.parseParameters(argc, argv, command, true, Parameters::PARSE_ALLOW_EMPTY, 0);
  
    // If dbDirectory does not exist, create it
    if (!FileUtil::directoryExists(par.filenames[0].c_str())) {
        FileUtil::makeDir(par.filenames[0].c_str());
    }

    // Create index
    IndexCreator idxCre(par);
    idxCre.createIndex(par);

    if(idxCre.getNumOfFlush() == 1) {
        cerr << "Index creation completed." << endl;
        return 0;
    }

    // Merge index files
    cout << "Merge reference DB files ... " << endl;
    int numOfSplits = idxCre.getNumOfFlush();
    FileMerger merger(par);
    for (int i = 0; i < numOfSplits; i++) {
        merger.addFilesToMerge(par.filenames[0] + "/" + to_string(i) + "_diffIdx",
                               par.filenames[0] + "/" + to_string(i) + "_info");
    }
    merger.updateTaxId2SpeciesTaxId(par.filenames[0] + "/taxID_list");
    merger.setMergedFileNames(par.filenames[0] + "/diffIdx", par.filenames[0] + "/info", par.filenames[0] + "/split");  
    merger.mergeTargetFiles();
    cerr << "Index creation completed." << endl;
    return 0;
}
