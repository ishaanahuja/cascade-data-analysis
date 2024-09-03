#include <THashList.h>
#include <TROOT.h>
#include "TStopwatch.h"

#include "CascadeUtils.h"

int RunVariations(TString inputFilename = "050824_updatedCutVar_6Runs.root", TString outputFilePrefix = "140824", Bool_t fisMC = kFALSE, Int_t verbosity = kInfo)
{

  ROOT::EnableImplicitMT();
  gErrorIgnoreLevel = verbosity;

  TStopwatch totalTimer, macroTimer;
  totalTimer.Start();

  TString histName = "h3_ptmasscent_def";
  TString histSuffix = "";
  Int_t errorCode = 0;
  TString outputFileName;

  Int_t iRun = 0;
  /// Begin input
  TH1::AddDirectory(0);
  TFile *inputFile = OpenFile(inputFilename);

  THashList *list_Omp = (THashList *)inputFile->FindObjectAny("chists_Omp_");
  Int_t nVar = list_Omp->GetEntries();
  Info("RunVariations_chists_Omp_", "Max TH3 histograms per particle = \e[1;32m%d\e[0m", nVar);
  Double_t execTime[nVar][2];

  /// DEFAULT:

  Info("RunVariations_\e[1;32mDEFAULT\e[0m", "Starting analysis with default cuts - \e[1;32m%s\e[0m", histName.Data());

  outputFileName = outputFilePrefix + "_" + histName + ".root";
  macroTimer.Start();
  gROOT->Macro(TString::Format("FitCascades.C+O(\"%s\", \"%s\", \"%s\")", inputFilename.Data(), outputFileName.Data(), histName.Data()), &errorCode);
  macroTimer.Stop();
  Printf("TIMER: <RUN_DEFAULT> Real =%7.3fs, CPU =%7.3fs", macroTimer.RealTime(), macroTimer.CpuTime());
  if (!errorCode)
  {
    execTime[iRun][0] = macroTimer.RealTime();
    execTime[iRun][1] = macroTimer.CpuTime();
    iRun++;

    Info("RunVariations_\e[1;32mDEFAULT\e[0m", "\e[1;32m%s\e[0m - Output successfully saved to <\e[1;94m%s\e[0m>", histName.Data(), outputFileName.Data());
  }
  else
    Error("RunVariations_\e[1;31mDEFAULT\e[0m", "\e[1;31m%s\e[0m - Unexpected Error! TInterpreter::EErrorCode = %d", histName.Data(), errorCode);

  for (Int_t iVar = 0; iVar < nVar; iVar++)
  {
    histName = list_Omp->At(iVar)->GetName(); // Omega analysis contains more no. of variations
    if (histName.Contains("h3Var"))
    {
      Info("RunVariations_Vars", "Analysing cut variation \e[1;93m%d - %s\e[0m ...", iVar, histName.Data());

      histSuffix = histName;
      histSuffix.ReplaceAll("][", "_");
      histSuffix.ReplaceAll("[", "_");
      histSuffix.ReplaceAll("]", "");

      outputFileName = outputFilePrefix + "_" + histSuffix + ".root";

      macroTimer.Start();
      gROOT->Macro(TString::Format("FitCascades.C+O(\"%s\", \"%s\", \"%s\")", inputFilename.Data(), outputFileName.Data(), histName.Data()), &errorCode);
      macroTimer.Stop();
      Printf("TIMER: <RUN_%d> Real =%7.3fs, CPU =%7.3fs", iRun, macroTimer.RealTime(), macroTimer.CpuTime());

      if (!errorCode)
      {
        execTime[iRun][0] = macroTimer.RealTime();
        execTime[iRun][1] = macroTimer.CpuTime();
        iRun++;
        Info("RunVariations_Vars", "Var \e[1;93m%d - %s\e[0m Output successfully saved to <\e[1;94m%s\e[0m>", iVar, histName.Data(), outputFileName.Data());
      }
      else
        Error("RunVariations_Vars", "Var \e[1;31m%d - %s\e[0m Unexpected Error! TInterpreter::EErrorCode = %d", iVar, histName.Data(), errorCode);
    }
  }

  totalTimer.Stop();
  Printf("    ===============================TIMER==================================");
  Printf("    <Run>                           Real (s)                         CPU (s)");
  for (Int_t iTime = 0; iTime < iRun; iTime++)
    Printf("    <Run_%d>                        %7.3f                           %7.3f", iTime, execTime[iTime][0], execTime[iTime][1]);
  Printf("    <Total>                        %7.3f                          %7.3f", totalTimer.RealTime(), totalTimer.CpuTime());

  return 0;
}