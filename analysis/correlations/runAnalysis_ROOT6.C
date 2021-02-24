#include "CorrelationTask.h"

void runAnalysis_ROOT6()
{
    // set if you want to run the analysis locally (kTRUE), or on grid (kFALSE)
    Bool_t local = kTRUE;
    // if you run on grid, specify test mode (kTRUE) or full grid model (kFALSE)
    Bool_t gridTest = kFALSE;
       TFile::SetCacheFileDir(gSystem->HomeDirectory(), 1, 1);
    // create the analysis manager
    AliAnalysisManager *mgr = new AliAnalysisManager("AnalysisTaskExample");
    AliAODInputHandler *aodH = new AliAODInputHandler();
    mgr->SetInputEventHandler(aodH);

        gInterpreter->LoadMacro("CorrelationTask.cxx++g");
        CorrelationTask *task=reinterpret_cast<CorrelationTask*>(gInterpreter->ExecuteMacro("AddMyTask.C"));
       
    if(!mgr->InitAnalysis()) return;
    mgr->SetDebugLevel(2);
    mgr->PrintStatus();
    mgr->SetUseProgressBar(1, 25);

    if(local) {
        // if you want to run locally, we need to define some input
        TChain* chain = new TChain("aodTree");
        // add a few files to the chain (change this so that your local files are added)

        chain->Add("/var/home/ishaan/Work/CERN/ishaan-ahuja/data_AOD/PbPbLHC11h/0001/AliAOD.root");
        chain->Add("/var/home/ishaan/Work/CERN/ishaan-ahuja/data_AOD/PbPbLHC11h/0002/AliAOD.root");
        chain->Add("/var/home/ishaan/Work/CERN/ishaan-ahuja/data_AOD/PbPbLHC11h/0003/AliAOD.root");
        chain->Add("/var/home/ishaan/Work/CERN/ishaan-ahuja/data_AOD/PbPbLHC11h/0004/AliAOD.root");
        chain->Add("/var/home/ishaan/Work/CERN/ishaan-ahuja/data_AOD/PbPbLHC11h/0005/AliAOD.root");
        chain->Add("/var/home/ishaan/Work/CERN/ishaan-ahuja/data_AOD/PbPbLHC11h/0006/AliAOD.root");
        chain->Add("/var/home/ishaan/Work/CERN/ishaan-ahuja/data_AOD/PbPbLHC11h/0007/AliAOD.root");
        chain->Add("/var/home/ishaan/Work/CERN/ishaan-ahuja/data_AOD/PbPbLHC11h/0008/AliAOD.root");
        chain->Add("/var/home/ishaan/Work/CERN/ishaan-ahuja/data_AOD/PbPbLHC11h/0009/AliAOD.root");
        

        // start the analysis locally, reading the events from the tchain
        mgr->StartAnalysis("local", chain);
    } else {
        // if we want to run on grid, we create and configure the plugin
        AliAnalysisAlien *alienHandler = new AliAnalysisAlien();
        // also specify the include (header) paths on grid
        alienHandler->AddIncludePath("-I. -I$ROOTSYS/include -I$ALICE_ROOT -I$ALICE_ROOT/include -I$ALICE_PHYSICS/include");
        // make sure your source files get copied to grid
        alienHandler->SetAdditionalLibs("CorrelationTask.cxx CorrelationTask.h");
        alienHandler->SetAnalysisSource("CorrelationTask.cxx");
        // select the aliphysics version. all other packages
        // are LOADED AUTOMATICALLY!
        alienHandler->SetAliPhysicsVersion("vAN-20201025_ROOT6-1");
        // set the Alien API version
        alienHandler->SetAPIVersion("V1.1x");
        // select the input data /alice/data/2018/LHC18h/000288804/pass1/AOD234/0030/
        alienHandler->SetGridDataDir("/alice/data/2011/LHC11h_2");
        alienHandler->SetDataPattern("*ESDs/pass2/AOD145/*AOD.root");
        // MC has no prefix, data has prefix 000
        alienHandler->SetRunPrefix("000");
        // runnumber
        alienHandler->AddRunNumber(167813);
        // number of files per subjob
        alienHandler->SetSplitMaxInputFileNumber(40);
        alienHandler->SetExecutable("myTask.sh");
        // specify how many seconds your job may take
        alienHandler->SetTTL(10000);
        alienHandler->SetJDLName("myTask.jdl");

        alienHandler->SetOutputToRunNo(kTRUE);
        alienHandler->SetKeepLogs(kTRUE);
        
		
		
		
        alienHandler->SetMaxMergeStages(1);
        alienHandler->SetMergeViaJDL(kTRUE);
        

        // define the output folders
        alienHandler->SetGridWorkingDir("myWorkingDir");
        alienHandler->SetGridOutputDir("myOutputDir");

        // connect the alien plugin to the manager
        mgr->SetGridHandler(alienHandler);
        if(gridTest) {
            // specify on how many files you want to run
            alienHandler->SetNtestFiles(1);
            // and launch the analysis
            alienHandler->SetRunMode("test");
            mgr->StartAnalysis("grid");
        } else {
            // else launch the full grid analysis
            alienHandler->SetRunMode("full");//
            mgr->StartAnalysis("grid");
        }
    }
}

