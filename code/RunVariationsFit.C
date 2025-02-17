#include <THashList.h>
#include <TROOT.h>
#include "TStopwatch.h"
#include "TDatime.h"
#include "CascadeUtils.h"

/**
 * @brief Runs a series of fit operations over default and variation histograms through compiled macro 'FitCascades.C'.
 *
 * This function enables multi-threading, opens the provided input file, retrieves
 * the relevant histograms, and runs the fit routine either on a default histogram
 * or on each variation histogram if specified. It measures execution times and
 * logs detailed messages based on the verbosity level.
 *
 * @param inputFilename   Path to the input ROOT file containing the histograms.
 * @param outputFolder    Directory where the resulting ROOT files will be saved.
 * @param outputFilePrefix String prefix for output file naming.
 * @param fisMC           Flag indicating whether the analysis is performed on Monte Carlo data.
 * @param fitDiffForVars  Flag indicating whether differential multiplicity fits should be done for variations.
 * @param fitDefault      Flag to perform fit on the default cut histogram as well.
 * @param verbosity       Controls the verbosity level for messages.
 * @return An integer status code (0 if successful, -1 if errors encountered).
 */
int RunVariationsFit(
    TString inputFilename = "/home/alice/iahuja/RandomVars/181224_6Runs_RandomVars.root",
    TString outputFolder = "/home/alice/iahuja/RandomVars/050225_6RunsRandDiff",
    TString outputFilePrefix = "050225_6RunsRandDiff",
    bool fisMC = false,
    bool fitDiffForVars = true,
    bool fitDefault = false,
    Int_t verbosity = kInfo)
{

  ROOT::EnableImplicitMT();
  gErrorIgnoreLevel = verbosity;

  TStopwatch totalTimer, macroTimer;
  totalTimer.Start();

  TString histName = "h3_ptmasscent_def";
  TString histSuffix = "";
  Int_t errorCode = 0;
  TString outputFileName;

  outputFolder = SetOutputFolder(outputFolder);

  Int_t iRun = 0;

  /// Begin input
  TH1::AddDirectory(0);
  TFile *inputFile = OpenFile(inputFilename);
  THashList *list_Omp = (THashList *)inputFile->FindObjectAny("chists_Omp_");
  if (!list_Omp)
  {
    Error("RunVariationsFit", "Failed to find object 'chists_Omp_' in the input file. Exiting ...");
    return -1;
  }
  inputFile->Close();
  delete inputFile;

  Int_t nVar = list_Omp->GetEntries();
  Info("RunVariationsFit_chists_Omp_", "Max TH3 histograms per particle = \e[1;32m%d\e[0m", nVar - 4);
  Double_t execTime[nVar][2];

  /// DEFAULT_CUTS:
  if (fitDefault)
  {
    TDatime now;
    if (fisMC)
      Info("RunVariationsFit_\e[1;32mDEFAULT\e[0m", "%s: Starting MC analysis with default cuts - \e[1;32m%s\e[0m", now.AsSQLString(), histName.Data());
    else
      Info("RunVariationsFit_\e[1;32mDEFAULT\e[0m", "%s: Starting analysis with default cuts - \e[1;32m%s\e[0m", now.AsSQLString(), histName.Data());

    outputFileName = outputFolder + "/" + outputFilePrefix + "_" + histName + ".root";
    macroTimer.Start();
    gROOT->Macro(TString::Format("FitCascades.C+O(\"%s\", \"%s\", \"%s\", %d)", inputFilename.Data(), outputFileName.Data(), histName.Data(), fisMC), &errorCode);
    macroTimer.Stop();
    Printf("TIMER: <RUN_DEFAULT> Real =%7.3fs, CPU =%7.3fs", macroTimer.RealTime(), macroTimer.CpuTime());
    if (!errorCode)
    {
      execTime[iRun][0] = macroTimer.RealTime();
      execTime[iRun][1] = macroTimer.CpuTime();
      iRun++;

      Info("RunVariationsFit_\e[1;32mDEFAULT\e[0m", "\e[1;32m%s\e[0m - Output successfully saved to <\e[1;94m%s\e[0m>", histName.Data(), outputFileName.Data());
    }
    else
      Error("RunVariationsFit_\e[1;31mDEFAULT\e[0m", "\e[1;31m%s\e[0m - Unexpected Error! TInterpreter::EErrorCode = %d", histName.Data(), errorCode);
  }

  /// CUT_VARIATIONS:
  for (Int_t iVar = 0; iVar < nVar; iVar++)
  {
    TDatime now;

    histName = list_Omp->At(iVar)->GetName(); // Omega analysis (used to) contain more variations (obsolete - now it contains the same Rand vars)
    if (histName.Contains("h3Var"))
    {
      Info("RunVariationsFit_Vars", "%s: Analysing cut variation \e[1;93m%d - %s\e[0m ...", now.AsSQLString(), iVar, histName.Data());

      histSuffix = histName;
      outputFileName = outputFolder + "/" + outputFilePrefix + "_" + histSuffix + ".root";

      macroTimer.Start();
      gROOT->Macro(TString::Format("FitCascades.C+O(\"%s\", \"%s\", \"%s\", %d, %d)", inputFilename.Data(), outputFileName.Data(), histName.Data(), fisMC, fitDiffForVars), &errorCode);
      macroTimer.Stop();
      Printf("TIMER: <RUN_%d> Real =%7.3fs, CPU =%7.3fs", iRun, macroTimer.RealTime(), macroTimer.CpuTime());

      if (!errorCode)
      {
        execTime[iRun][0] = macroTimer.RealTime();
        execTime[iRun][1] = macroTimer.CpuTime();
        iRun++;
        Info("RunVariationsFit_Vars", "iVar \e[1;93m%d - %s\e[0m Output successfully saved to <\e[1;94m%s\e[0m>", iVar, histName.Data(), outputFileName.Data());
      }
      else
        Error("RunVariationsFit_Vars", "iVar \e[1;31m%d - %s\e[0m Unexpected Error! TInterpreter::EErrorCode = %d", iVar, histName.Data(), errorCode);
    }
  }

  totalTimer.Stop();
  Printf("    ===============================TIMER==================================");
  Printf("    <Run>                           Real (s)                         CPU (s)");
  for (Int_t iTime = 0; iTime < iRun; iTime++)
    Printf("    <Run_%d>                        %7.3f                           %7.3f", iTime, execTime[iTime][0], execTime[iTime][1]);
  Printf("    <Total>                        %7.3f                          %7.3f", totalTimer.RealTime(), totalTimer.CpuTime());

  delete list_Omp;
  return 0;
}
