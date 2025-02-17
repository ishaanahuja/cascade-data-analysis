#include "AliAnalysisTaskStrangeCascadesRun2.h"

void RunAnalysisTask_StrangeCascadesRun2()
{
    // set if you want to run the analysis locally (kTRUE), or on grid (kFALSE)
    Bool_t local = kFALSE;
    // TString outputFilename = "testLocal.root";
    // if you run on grid, specify test mode (kTRUE) or full grid model (kFALSE)
    Bool_t gridTest = kFALSE;
    // Set Cache
    // if (local)
    //     TFile::SetCacheFileDir(gSystem->HomeDirectory(), 1, 1);
    // create the analysis manager
    AliAnalysisManager *mgr = new AliAnalysisManager("StrangeCascadesAnalysisTask");
    AliAODInputHandler *aodH = new AliAODInputHandler();
    mgr->SetInputEventHandler(aodH);

    // PID task aa
    //  load the macro and add the task
    TMacro PIDadd(gSystem->ExpandPathName("$ALICE_ROOT/ANALYSIS/macros/AddTaskPIDResponse.C"));
    AliAnalysisTaskPIDResponse *PIDresponseTask = reinterpret_cast<AliAnalysisTaskPIDResponse *>(PIDadd.Exec());

    /*gInterpreter->LoadMacro("PIDresponseTask.C++g");
    AliAnalysisTaskMyTask *PIDtask=reinterpret_cast<AliAnalysisTaskMyTask*>(gInterpreter->ExecuteMacro("PIDresponseTask.C"));
    */

    // centrality
    TMacro multSelection(gSystem->ExpandPathName("$ALICE_PHYSICS/OADB/COMMON/MULTIPLICITY/macros/AddTaskMultSelection.C"));
    AliMultSelectionTask *multSelectionTask = reinterpret_cast<AliMultSelectionTask *>(multSelection.Exec());

    gInterpreter->LoadMacro("AliAnalysisTaskStrangeCascadesRun2.cxx++g");

    // AliAnalysisTaskStrangeCascadesRun2 *task = reinterpret_cast<AliAnalysisTaskStrangeCascadesRun2 *>(gInterpreter->ExecuteMacro(TString::Format("AddTaskStrangeCascadesRun2.C'(\"%s\")'", outputFilename.Data())));
    AliAnalysisTaskStrangeCascadesRun2 *task = reinterpret_cast<AliAnalysisTaskStrangeCascadesRun2 *>(gInterpreter->ExecuteMacro("AddTaskStrangeCascadesRun2.C"));

    if (!mgr->InitAnalysis())
        return;
    mgr->SetDebugLevel(2);
    mgr->PrintStatus();
    mgr->SetUseProgressBar(1, 25);

    if (local)
    {
        // if you want to run locally, we need to define some input
        TChain *chain = new TChain("aodTree");
        // add a few files to the chain (change this so that your local files are added)
        // chain->Add("/var/home/ishaan/Work/CERN/ishaan-ahuja/data_AOD/PbPbLHC15o_pass5/AliAOD.root");
        chain->Add("/var/home/ishaan/Work/CERN/ishaan-ahuja/data_AOD/LHC17e1b_cent/265594/AOD202/0004/AliAOD.root"); // MC
        // chain->Add("/var/home/ishaan/alice/sim/2017/LHC17d14_cent/265596/AOD202/0001/AliAOD.root"); // MC
        // chain->Add("/var/home/ishaan/Work/CERN/ishaan-ahuja/data_AOD/alice/sim/2017/LHC17f3a_cent_fix/265594/AOD202/0002/AliAOD.root"); // MC
        // chain->Add("/var/home/ishaan/Work/CERN/ishaan-ahuja/data_AOD/LHC17d14_cent/265594/AOD202/0008/AliAOD.root"); // MC
        // chain->Add("/var/home/ishaan/Work/CERN/ishaan-ahuja/data_AOD/pPbLHC16r/000265596/pass1_CENT_wSDD/AOD190/0001/AliAOD.root"); // good
        // chain->Add("/var/home/ishaan/Work/CERN/ishaan-ahuja/data_AOD/pPbLHC16r/000265596/pass1_CENT_wSDD/AOD190/0035/AliAOD.root"); // good
        // chain->Add("/var/home/ishaan/Work/CERN/ishaan-ahuja/data_AOD/pPbLHC16r/000266076/pass1_CENT_wSDD/AOD190/0060/AliAOD.root"); // pileup
        // chain->Add("/var/home/ishaan/Work/CERN/ishaan-ahuja/data_AOD/pPbLHC16r/000266076/pass1_CENT_wSDD/AOD190/0246/AliAOD.root"); // pileup

        // start the analysis locally, reading the events from the tchain
        mgr->StartAnalysis("local", chain);
    }
    else
    {
        // if we want to run on grid, we create and configure the plugin
        AliAnalysisAlien *alienHandler = new AliAnalysisAlien();
        // also specify the include (header) paths on grid
        alienHandler->AddIncludePath("-I. -I$ROOTSYS/include -I$ALICE_ROOT -I$ALICE_ROOT/include -I$ALICE_PHYSICS/include");
        // make sure your source files get copied to grid
        alienHandler->SetAdditionalLibs("AliAnalysisTaskStrangeCascadesRun2.h AliAnalysisTaskStrangeCascadesRun2.cxx");
        alienHandler->SetAnalysisSource("AliAnalysisTaskStrangeCascadesRun2.cxx");
        // select the aliphysics version. all other packages
        // are LOADED AUTOMATICALLY!
        alienHandler->SetAliPhysicsVersion("vAN-20221009_ROOT6-1");
        // alienHandler->SetAliPhysicsVersion("vAN-20230612_O2-1");
        // alienHandler->SetAliPhysicsVersion("vAN-20191017_ROOT6-1");

        alienHandler->SetAPIVersion("V1.1x");

        // alienHandler->SetGridDataDir("/alice/sim/2018/LHC18f3b_cent_2");          // MC
        // alienHandler->SetGridDataDir("/alice/sim/2017/LHC17e1a_cent");          // MC
        // alienHandler->SetGridDataDir("/alice/sim/2017/LHC17e1b_cent");          // MC
        alienHandler->SetGridDataDir("/alice/sim/2017/LHC17l7a2_cent");          // MC  <<<<<<<<<<<<<< NEXT
        // alienHandler->SetGridDataDir("/alice/sim/2017/LHC17f3a_cent_fix");          // MC 
        // alienHandler->SetGridDataDir("/alice/data/2016/LHC16r");                //Data

        // alienHandler->SetDataPattern("*ESDs/pass2/AOD145/*AOD.root");
        // alienHandler->SetDataPattern("*pass3/AOD252/*AOD.root");
        alienHandler->SetDataPattern("*AOD202/*AliAOD.root");                      // MC
        // alienHandler->SetDataPattern("*pass2_CENT_wSDD/AOD244/*AliAOD.root");   //Data

        // MC has no prefix, data has prefix 000
        // alienHandler->SetRunPrefix("000");                                      //Data
        // runnumber
        /// 6 good runs:
        alienHandler->AddRunNumber(265594);
        alienHandler->AddRunNumber(265596);
        alienHandler->AddRunNumber(265607);
        alienHandler->AddRunNumber(266316);
        alienHandler->AddRunNumber(266317);
        alienHandler->AddRunNumber(266318);

        /// pileup grp1:
        // alienHandler->AddRunNumber(265714);
        // alienHandler->AddRunNumber(265713);
        // alienHandler->AddRunNumber(265756);
        // alienHandler->AddRunNumber(265709);
        // alienHandler->AddRunNumber(265754);
        // alienHandler->AddRunNumber(265705);
        // alienHandler->AddRunNumber(266208);
        // alienHandler->AddRunNumber(265701);
        // alienHandler->AddRunNumber(265797);
        // alienHandler->AddRunNumber(265700);
        // alienHandler->AddRunNumber(265795);

        /// pileup grp3:
        // alienHandler->AddRunNumber(266305);
        // alienHandler->AddRunNumber(266085);
        // alienHandler->AddRunNumber(266296);
        // alienHandler->AddRunNumber(266300);
        // alienHandler->AddRunNumber(265788);
        // alienHandler->AddRunNumber(265744);
        // alienHandler->AddRunNumber(266197);
        // alienHandler->AddRunNumber(266084);
        // alienHandler->AddRunNumber(266083);
        // alienHandler->AddRunNumber(266196);
        // alienHandler->AddRunNumber(266081);
        // alienHandler->AddRunNumber(265742);
        // alienHandler->AddRunNumber(266193);
        // alienHandler->AddRunNumber(266076);
        // alienHandler->AddRunNumber(266190);
        // alienHandler->AddRunNumber(265741);
        // alienHandler->AddRunNumber(266189);
        // alienHandler->AddRunNumber(266074);
        // alienHandler->AddRunNumber(266187);
        // number of files per subjob
        alienHandler->SetSplitMaxInputFileNumber(40);
        alienHandler->SetExecutable("myTask.sh");
        // specify how many seconds your job may take
        alienHandler->SetTTL(35000);
        alienHandler->SetJDLName("myTask.jdl");

        alienHandler->SetOutputToRunNo(kTRUE);
        alienHandler->SetKeepLogs(kTRUE);

        alienHandler->SetMaxMergeStages(1);
        // alienHandler->SetMergeViaJDL(kTRUE); //
        alienHandler->SetMergeViaJDL(kFALSE); //

        // define the output folders
        // alienHandler->SetGridWorkingDir("081023_MC_LHC17f3b_cent_updatedCutsTL_WD");
        alienHandler->SetGridWorkingDir("100824_updatedCutVar_MC_LHC17l7a2cent_WD"); 
        // alienHandler->SetGridWorkingDir("090824_updatedCutVar_MC_LHC17e1bcent_WD");
        // alienHandler->SetGridOutputDir("081023_MC_LHC17f3b_cent_updatedCutsTL_OD");
        alienHandler->SetGridOutputDir("100824_updatedCutVar_MC_LHC17l7a2cent_OD");
        // alienHandler->SetGridOutputDir("090824_updatedCutVar_MC_LHC17e1bcent_OD");
        // connect the alien plugin to the manager
        mgr->SetGridHandler(alienHandler);
        if (gridTest)
        {
            // specify on how many files you want to run
            alienHandler->SetNtestFiles(1);
            // and launch the analysis
            alienHandler->SetRunMode("full");
            mgr->StartAnalysis("grid");
        }
        else
        {
            // else launch the full grid analysis
            // alienHandler->SetRunMode("full"); //
            alienHandler->SetRunMode("terminate"); //
            mgr->StartAnalysis("grid");
        }
    }
}
