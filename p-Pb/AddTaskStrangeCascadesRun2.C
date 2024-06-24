AliAnalysisTaskStrangeCascadesRun2 *AddTaskStrangeCascadesRun2(bool Xi = true, bool Om = true, TString suffix = "")
{
  // analysis manager
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr)
  {
    ::Error("AddTaskStrangeCascadesRun2", "No analysis manager to connect to.");
    return NULL;
  }
  if (!mgr->GetInputEventHandler())
  {
    ::Error("AddTaskStrangeCascadesRun2", "This task requires an input event handler");
    return NULL;
  }

  // Create the task and add it to the manager
  TString tskname = Form("Cascade_analysis_testVar%s", suffix.Data());
  AliAnalysisTaskStrangeCascadesRun2 *mytask = new AliAnalysisTaskStrangeCascadesRun2(tskname);
  mgr->AddTask(mytask);
  mytask->SetParticleAnalysisStatus(Xi, Om);

  // output file name
  // TString outputFileName = AliAnalysisManager::GetCommonFileName();
  // TString outputFileName = "test_local_newCutsMC.root";
  // TString outputFileName = "300923_6runs_updatedTL.root";
  // TString outputFileName = "081023_MC_LHC17f3b_cent_updatedCutsTL.root";
  TString outputFileName = "230424_cutVarNew_6Runs.root";

  outputFileName += ":Cascade_analysis";
  printf("Set OutputFileName : \n %s\n", outputFileName.Data());

  //create and link only used containers
  AliAnalysisDataContainer *coutput[8];
  TString clabels[5] = {"eve", "Xim", "Xip", "Omm", "Omp"};

  mgr->ConnectInput(mytask, 0, mgr->GetCommonInputContainer());
  coutput[0] = mgr->CreateContainer(Form("chists_%s_%s", clabels[0].Data(), suffix.Data()), TList::Class(), AliAnalysisManager::kOutputContainer, outputFileName);
  mgr->ConnectOutput(mytask, 1, coutput[0]);

  int cnumber = 1;
  for (int icont = 1; icont < 5; icont++)
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
