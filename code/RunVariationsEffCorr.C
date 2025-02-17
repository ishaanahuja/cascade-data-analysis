#include <TROOT.h>
#include "TStopwatch.h"
#include "TDatime.h"
#include "CascadeUtils.h"
#pragma omp parallel for

int RunVariationsEffCorr(
    TString fitInputFilePrefix = "/var/home/ishaan/Work/git/analysis/results/RandomVars/050225_6RunsRandDiff_fit/050225_6RunsRandDiff",
    TString effInputFilePrefix = "/var/home/ishaan/Work/git/analysis/results/RandomVars/030225_efficiencyVars/030225_eff",
    TString outputFilePrefix = "100225_effCorr",
    TString outputFolder = "/var/home/ishaan/Work/git/analysis/results/RandomVars/100225_effCorr",
    Bool_t calcDefault = kTRUE,
    Bool_t saveStack = kTRUE,
    TString imageFormat = "png",
    Int_t verbosity = kInfo)
{
    ROOT::EnableImplicitMT();
    gErrorIgnoreLevel = verbosity;

    TStopwatch totalTimer, macroTimer;
    totalTimer.Start();

    TString histName = "";
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
        Info("RunVariationsEffCorr_\e[1;32mDEFAULT\e[0m", "%s: Calculating efficiency corrected spectra for default cuts - \e[1;32m%s\e[0m", now.AsSQLString(), histName.Data());

        outputFileName = SetOutputFolder(outputFolder + "/" + histName) + "/" + outputFilePrefix + "_" + histName + ".root"; // creates a new folder for each histName if it doesn't exist
        macroTimer.Start();
        gROOT->Macro(TString::Format("EfficiencyCorrection.C+O(\"%s\", \"%s\", \"%s\", \"%s\", \"%s\", %d, \"%s\")", fitInputFilePrefix.Data(), effInputFilePrefix.Data(), histName.Data(), outputFileName.Data(), outputFolder.Data(), saveStack, imageFormat.Data()), &errorCode);
        macroTimer.Stop();
        Printf("TIMER: <RUN_DEFAULT> Real =%7.3fs, CPU =%7.3fs", macroTimer.RealTime(), macroTimer.CpuTime());
        if (!errorCode)
        {
            execTime[iRun][0] = macroTimer.RealTime();
            execTime[iRun][1] = macroTimer.CpuTime();
            iRun++;

            Info("RunVariationsEffCorr_\e[1;32mDEFAULT\e[0m", "\e[1;32m%s\e[0m - Output successfully saved to <\e[1;94m%s\e[0m>", histName.Data(), outputFileName.Data());
        }
        else
            Error("RunVariationsEffCorr_\e[1;31mDEFAULT\e[0m", "\e[1;31m%s\e[0m - Unexpected Error! TInterpreter::EErrorCode = %d", histName.Data(), errorCode);
    }

    /// VARIATIONS:
    for (Int_t iVar = 0; iVar < nVar; iVar++)
    {
        TDatime now;
        histName = TString::Format("h3Var_%d", iVar);
        Info("RunVariationsEffCorr_Vars", "%s: Calculating efficiency corrected spectra for cut variation \e[1;93m%d - %s\e[0m ...", now.AsSQLString(), iVar, histName.Data());

        outputFileName = SetOutputFolder(outputFolder + "/" + histName) + "/" + outputFilePrefix + "_" + histName + ".root"; // creates a new folder for each histName if it doesn't exist

        macroTimer.Start();
        gROOT->Macro(TString::Format("EfficiencyCorrection.C+O(\"%s\", \"%s\", \"%s\", \"%s\", \"%s\", %d, \"%s\")", fitInputFilePrefix.Data(), effInputFilePrefix.Data(), histName.Data(), outputFileName.Data(), outputFolder.Data(), saveStack, imageFormat.Data()), &errorCode);
        macroTimer.Stop();
        Printf("TIMER: <RUN_%d> Real =%7.3fs, CPU =%7.3fs", iRun, macroTimer.RealTime(), macroTimer.CpuTime());

        if (!errorCode)
        {
            execTime[iRun][0] = macroTimer.RealTime();
            execTime[iRun][1] = macroTimer.CpuTime();
            iRun++;
            Info("RunVariationsEffCorr_Vars", "iVar \e[1;93m%d - %s\e[0m Output successfully saved to <\e[1;94m%s\e[0m>", iVar, histName.Data(), outputFileName.Data());
        }
        else
            Error("RunVariationsEffCorr_Vars", "iVar \e[1;31m%d - %s\e[0m Unexpected Error! TInterpreter::EErrorCode = %d", iVar, histName.Data(), errorCode);
    }

    totalTimer.Stop();
    Printf("    ===============================TIMER==================================");
    Printf("    <Run>                           Real (s)                         CPU (s)");
    for (Int_t iTime = 0; iTime < iRun; iTime++)
        Printf("    <Run_%d>                        %7.3f                           %7.3f", iTime, execTime[iTime][0], execTime[iTime][1]);
    Printf("    <Total>                        %7.3f                          %7.3f", totalTimer.RealTime(), totalTimer.CpuTime());

    return 0;
}
