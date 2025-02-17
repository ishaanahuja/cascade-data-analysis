#include <TROOT.h>
#include "TStopwatch.h"
#include "TDatime.h"
#include "CascadeUtils.h"

/**
 * @brief Runs efficiency calculations with different variations of cuts
 *
 * This macro performs efficiency calculations for both default cuts and multiple variations.
 * It enables multi-threading, handles output file management, and provides execution timing information.
 *
 * @param inputFilename Path to the input file list (default: "/var/home/ishaan/Work/git/analysis/code/McFileListPrefix.txt")
 * @param outputFolder Directory path for output files (default: "/var/home/ishaan/Work/git/analysis/results/030225_efficiencyVars")
 * @param outputFilePrefix Prefix for output file names (default: "030225_eff")
 * @param calcDefault Flag to calculate efficiency with default cuts (default: kTRUE)
 * @param saveStack Flag to save stack histograms (default: kFALSE)
 * @param verbosity Level of verbosity for error messages (default: kInfo)
 *
 * @return Integer error code (0 for success)
 *
 * The function performs the following operations:
 * - Enables implicit multi-threading
 * - Sets up output directory
 * - If calcDefault is true, calculates efficiency with default cuts
 * - Processes 500 variations of cuts
 * - Records and displays execution times for each run
 * - Saves results to individual ROOT files
 *
 * Output files are named as: outputFolder/outputFilePrefix_histName.root
 */
int RunVariationsEff(
    TString inputFilename = "/var/home/ishaan/Work/git/analysis/code/McFileListPrefix.txt",
    TString outputFolder = "/var/home/ishaan/Work/git/analysis/results/030225_efficiencyVars",
    TString outputFilePrefix = "030225_eff",
    Bool_t calcDefault = kTRUE,
    Bool_t saveStack = kFALSE,
    Int_t verbosity = kInfo)
{

  ROOT::EnableImplicitMT();
  gErrorIgnoreLevel = verbosity;

  TStopwatch totalTimer, macroTimer;
  totalTimer.Start();

  TString histName = "";
  TString histSuffix = "";
  Bool_t analyseDiffs = kFALSE;
  Int_t errorCode = 0;
  TString outputFileName;

  outputFolder = SetOutputFolder(outputFolder);

  Int_t iRun = 0;

  /// Begin input
  Int_t nVar = 500;
  Double_t execTime[nVar][2];

  /// DEFAULT:
  if (calcDefault)
  {
    TDatime now;
    histName = "h3_ptmasscent_def";
    analyseDiffs = kTRUE;
    Info("RunVariationsEff_\e[1;32mDEFAULT\e[0m", "%s: Calculating efficiency for default cuts - \e[1;32m%s\e[0m", now.AsSQLString(), histName.Data());

    outputFileName = outputFolder + "/" + outputFilePrefix + "_" + histName + ".root";
    macroTimer.Start();
    gROOT->Macro(TString::Format("EfficiencyEstimation.C+O(\"%s\", \"%s\", %d,\"%s\", \"%s\", %d)", inputFilename.Data(), histName.Data(), analyseDiffs, outputFileName.Data(), outputFolder.Data(), saveStack), &errorCode);
    macroTimer.Stop();
    Printf("TIMER: <RUN_DEFAULT> Real =%7.3fs, CPU =%7.3fs", macroTimer.RealTime(), macroTimer.CpuTime());
    if (!errorCode)
    {
      execTime[iRun][0] = macroTimer.RealTime();
      execTime[iRun][1] = macroTimer.CpuTime();
      iRun++;

      Info("RunVariationsEff_\e[1;32mDEFAULT\e[0m", "\e[1;32m%s\e[0m - Output successfully saved to <\e[1;94m%s\e[0m>", histName.Data(), outputFileName.Data());
    }
    else
      Error("RunVariationsEff_\e[1;31mDEFAULT\e[0m", "\e[1;31m%s\e[0m - Unexpected Error! TInterpreter::EErrorCode = %d", histName.Data(), errorCode);
  }

  /// VARIATIONS:
  // #pragma omp parallel for
  for (Int_t iVar = 0; iVar < nVar; iVar++)
  {
    TDatime now;
    histName = TString::Format("h3Var_%d", iVar);
    analyseDiffs = kFALSE;
    Info("RunVariationsEff_Vars", "%s: Calculating efficiency for cut variation \e[1;93m%d - %s\e[0m ...", now.AsSQLString(), iVar, histName.Data());

    histSuffix = histName;
    outputFileName = outputFolder + "/" + outputFilePrefix + "_" + histSuffix + ".root";

    macroTimer.Start();
    gROOT->Macro(TString::Format("EfficiencyEstimation.C+O(\"%s\", \"%s\", %d,\"%s\", \"%s\", %d)", inputFilename.Data(), histName.Data(), analyseDiffs, outputFileName.Data(), outputFolder.Data(), saveStack), &errorCode);
    macroTimer.Stop();
    Printf("TIMER: <RUN_%d> Real =%7.3fs, CPU =%7.3fs", iRun, macroTimer.RealTime(), macroTimer.CpuTime());

    if (!errorCode)
    {
      execTime[iRun][0] = macroTimer.RealTime();
      execTime[iRun][1] = macroTimer.CpuTime();
      iRun++;
      Info("RunVariationsEff_Vars", "iVar \e[1;93m%d - %s\e[0m Output successfully saved to <\e[1;94m%s\e[0m>", iVar, histName.Data(), outputFileName.Data());
    }
    else
      Error("RunVariationsEff_Vars", "iVar \e[1;31m%d - %s\e[0m Unexpected Error! TInterpreter::EErrorCode = %d", iVar, histName.Data(), errorCode);
  }

  totalTimer.Stop();
  Printf("    ===============================TIMER==================================");
  Printf("    <Run>                           Real (s)                         CPU (s)");
  for (Int_t iTime = 0; iTime < iRun; iTime++)
    Printf("    <Run_%d>                        %7.3f                           %7.3f", iTime, execTime[iTime][0], execTime[iTime][1]);
  Printf("    <Total>                        %7.3f                          %7.3f", totalTimer.RealTime(), totalTimer.CpuTime());

  return 0;
}
