#include <TString.h>

AliAnalysisTaskStrangeCascadesRun2 *AddTaskStrangeCascadesRun2(bool Xi = true, bool Om = true, TString suffix = "")
{
  // analysis manager
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr)
  {
    ::Error("AddTaskStrangeCascadesRun2", "No analysis manager to connect to. Please ensure that the AliAnalysisManager is properly initialized.");
    return NULL;
  }
  if (!mgr->GetInputEventHandler())
  {
    ::Error("AddTaskStrangeCascadesRun2", "This task requires an input event handler. Please ensure that an input event handler is set before adding this task.");
    return NULL;
  }

  // Create the task and add it to the manager
  TString tskname = Form("Cascade_analysis%s", suffix.Data());
  AliAnalysisTaskStrangeCascadesRun2 *mytask = new AliAnalysisTaskStrangeCascadesRun2(tskname);
  mgr->AddTask(mytask);
  mytask->SetParticleAnalysisStatus(Xi, Om);

  // output file name
  // The output file name is hardcoded for a specific Monte Carlo/real dataset
  // TString outputFileName = "261124_6Runs_RandomVars.root";
  TString outputFileName = "150125_LHC18f3b_cent_2_RandomVars.root"; // MC1
  // TString outputFileName = "110125_LHC17e1a_cent_RandomVars.root"; //MC2
  // TString outputFileName = "110125_LHC17e1b_cent_RandomVars.root"; //MC3
  // TString outputFileName = "130125_LHC17l7a2_cent_RandomVars.root"; //MC4
  // TString outputFileName = "140125_LHC17f3a_cent_fix_RandomVars.root"; //MC5

  outputFileName += ":Cascade_analysis";
  printf("Set OutputFileName : \n %s\n", outputFileName.Data());

  // create and link only used containers
  Int_t numContainers = 5;
  AliAnalysisDataContainer *coutput[numContainers];
  TString clabels[numContainers] = {"eve", "Xim", "Xip", "Omm", "Omp"};

  mgr->ConnectInput(mytask, 0, mgr->GetCommonInputContainer());
  coutput[0] = mgr->CreateContainer(Form("chists_%s_%s", clabels[0].Data(), suffix.Data()), TList::Class(), AliAnalysisManager::kOutputContainer, outputFileName);
  mgr->ConnectOutput(mytask, 1, coutput[0]);

  int cnumber = 1;
  for (int icont = 1; icont < numContainers; icont++)
  {
    if (mytask->GetParticleAnalysisStatus(icont - 1))
    {
      coutput[cnumber] = mgr->CreateContainer(Form("chists_%s_%s", clabels[icont].Data(), suffix.Data()), TList::Class(), AliAnalysisManager::kOutputContainer, outputFileName);
      mgr->ConnectOutput(mytask, cnumber + 1, coutput[cnumber]);
      cnumber++;
    }
  }

  return mytask;
}
