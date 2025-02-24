#include <TROOT.h>
#include "TStopwatch.h"
#include "TDatime.h"
#include "CascadeUtils.h"


int RunVariationsSigExt(
    TString outputFolder = "/var/home/ishaan/Work/git/analysis/results/RandomVars/190225_SysUncertSigExt",
    TString outputFilePrefix = "190225_SysSigExt",
    TString histSuffix = "DGP2_10bg",
    Int_t fitFunction = kDoubleGausPol2,
    bool fitReal = true,
    Int_t verbosity = kInfo)
{

  ROOT::EnableImplicitMT();
  gErrorIgnoreLevel = verbosity;

  TStopwatch totalTimer, macroTimer;
  totalTimer.Start();

  bool fisMC = false;
  Int_t errorCode = 0;
  Double_t execTime[6][2]; // [Run][Real/CPU]
  TString histName = "h3_ptmasscent_def";

  outputFolder = SetOutputFolder(outputFolder);

  Int_t iRun = 0;

  TString input_realData = "/var/home/ishaan/Work/git/analysis/results/RandomVars/181224_6Runs_RandomVars.root";
  TString output_realData = TString::Format("%s/%s_%s_6Runs.root", outputFolder.Data(), outputFilePrefix.Data(), histSuffix.Data());

  TString input_MC[5] = {"/var/home/ishaan/Work/git/analysis/results/RandomVars/090125_LHC17e1a_cent_RandomVars.root",
                         "/var/home/ishaan/Work/git/analysis/results/RandomVars/110125_LHC17e1b_cent_RandomVars.root",
                         "/var/home/ishaan/Work/git/analysis/results/RandomVars/130125_LHC17l7a2_cent_RandomVars.root",
                         "/var/home/ishaan/Work/git/analysis/results/RandomVars/150125_LHC18f3b_cent_2_RandomVars.root",
                         "/var/home/ishaan/Work/git/analysis/results/RandomVars/140125_LHC17f3a_cent_fix_RandomVars.root"};

  TString output_MC[5] = {TString::Format("%s/%s_%s_OmLHC17e1a.root", outputFolder.Data(), outputFilePrefix.Data(), histSuffix.Data()),
                          TString::Format("%s/%s_%s_XiLHC17e1b.root", outputFolder.Data(), outputFilePrefix.Data(), histSuffix.Data()),
                          TString::Format("%s/%s_%s_LHC17l7a2.root", outputFolder.Data(), outputFilePrefix.Data(), histSuffix.Data()),
                          TString::Format("%s/%s_%s_LHC18f3b.root", outputFolder.Data(), outputFilePrefix.Data(), histSuffix.Data()),
                          TString::Format("%s/%s_%s_LHC17f3a.root", outputFolder.Data(), outputFilePrefix.Data(), histSuffix.Data())};
  TH1::AddDirectory(0);

  /// REAL_DATA:
  if (fitReal)
  {
    fisMC = false;
    TDatime now;

    Info("RunVariationsSigExt_\e[1;32mREAL\e[0m", "%s: Fitting real data \e[1;32m'%s' for '%s'\e[0m", now.AsSQLString(), input_realData.Data(), histSuffix.Data());

    macroTimer.Start();
    gROOT->Macro(TString::Format("FitCascades.C+O(\"%s\", \"%s\", \"%s\", %d, %d, \"%s\",%d)", input_realData.Data(), output_realData.Data(), histName.Data(), fisMC, true, "ye", fitFunction), &errorCode);
    macroTimer.Stop();

    Printf("TIMER: <RUN_REAL> Real =%7.3fs, CPU =%7.3fs", macroTimer.RealTime(), macroTimer.CpuTime());
    if (!errorCode)
    {
      execTime[iRun][0] = macroTimer.RealTime();
      execTime[iRun][1] = macroTimer.CpuTime();
      iRun++;

      Info("RunVariationsSigExt_\e[1;32mREAL\e[0m", "\e[1;32m%s\e[0m - Output successfully saved to <\e[1;94m%s\e[0m>", histSuffix.Data(), output_realData.Data());
    }
    else
      Error("RunVariationsSigExt_\e[1;31mREAL\e[0m", "\e[1;31m%s\e[0m - Unexpected Error! TInterpreter::EErrorCode = %d", histSuffix.Data(), errorCode);
  }

  /// MC_RUNS:
  for (Int_t iVar = 0; iVar < 5; iVar++)
  {
    fisMC = true;
    TDatime now;

    Info("RunVariationsSigExt_MC", "%s: Fitting (%d) MC Data \e[1;93m'%s' for '%s'\e[0m ...", now.AsSQLString(), iVar, input_MC[iVar].Data(), histSuffix.Data());

    macroTimer.Start();
    gROOT->Macro(TString::Format("FitCascades.C+O(\"%s\", \"%s\", \"%s\", %d, %d, \"%s\",%d)", input_MC[iVar].Data(), output_MC[iVar].Data(), histName.Data(), fisMC, true, "ye", fitFunction), &errorCode);
    macroTimer.Stop();
    Printf("TIMER: <RUN_MC_%d> Real =%7.3fs, CPU =%7.3fs", iRun, macroTimer.RealTime(), macroTimer.CpuTime());

    if (!errorCode)
    {
      execTime[iRun][0] = macroTimer.RealTime();
      execTime[iRun][1] = macroTimer.CpuTime();
      iRun++;
      Info("RunVariationsSigExt_MC", "iVar \e[1;93m%d - %s\e[0m Output successfully saved to <\e[1;94m%s\e[0m>", iVar, histSuffix.Data(), output_MC[iVar].Data());
    }
    else
      Error("RunVariationsSigExt_MC", "iVar \e[1;31m%d - %s\e[0m Unexpected Error! TInterpreter::EErrorCode = %d", iVar, histSuffix.Data(), errorCode);
  }

  totalTimer.Stop();
  Printf("    ===============================TIMER==================================");
  Printf("    <Run>                           Real (s)                         CPU (s)");
  for (Int_t iTime = 0; iTime < iRun; iTime++)
    Printf("    <Run_%d>                        %7.3f                           %7.3f", iTime, execTime[iTime][0], execTime[iTime][1]);
  Printf("    <Total>                        %7.3f                          %7.3f", totalTimer.RealTime(), totalTimer.CpuTime());

  return 0;
}
