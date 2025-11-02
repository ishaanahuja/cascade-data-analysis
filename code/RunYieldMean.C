#include <TROOT.h>
#include "TStopwatch.h"
#include "TDatime.h"
#include "TFile.h"
#include "TSystem.h"
#include "TError.h"
#include "TString.h"
#include "TH1.h"
#include "TF1.h"
#include "TStyle.h"
#include "TCanvas.h" // Needed for plotting later maybe

/// Dirty, but works for now. In future, consider making a proper library.
#include "YieldMean.C" // The core yield/mean calculation function

//-----------------------------------------------------------------------------
// Main Runner Function
//-----------------------------------------------------------------------------
/**
 * @brief Runs the integrated yield and mean pT analysis for various particles using multiple fit functions.
 *
 * This function processes input ROOT files containing efficiency-corrected spectra histograms,
 * fits the spectra with several physics-motivated functions (e.g., Levy-Tsallis, Boltzmann, Blast-Wave),
 * and computes integrated yields and mean transverse momentum (⟨pT⟩) for each particle species.
 * Results are saved to an output ROOT file, and fit logs are written to a specified directory.
 *
 * The analysis is performed for both multiplicity-integrated and multiplicity-binned spectra.
 * For each particle and multiplicity bin, the function:
 *   - Retrieves statistical and systematic histograms.
 *   - Applies each fit function over a defined pT range.
 *   - Integrates the fit to obtain yields and ⟨pT⟩, including extrapolation outside measured range.
 *   - Stores results and logs for further analysis and plotting.
 *
 * @param inputFileName      Path to the input ROOT file containing spectra histograms.
 * @param outputFileName     Path to the output ROOT file where results will be stored.
 * @param outputLogFolder    Directory where fit log files will be saved.
 * @param verbosity          ROOT message verbosity level (e.g., kInfo, kWarning).
 * @return                   0 on success, non-zero on error (e.g., file open failure).
 *
 * @note
 * - Requires external utilities/macros: YieldMean, ExtraUtils, CascadeUtils, and fit function definitions.
 * - Particle and binning definitions are expected to be set via global constants/arrays.
 * - Fit function parameters and ranges may need tuning for optimal results.
 * - Designed for use in high-energy nuclear/particle physics analyses (ALICE).
 *
 * Example usage (ROOT macro):
 * @code
 * .L RunYieldMean.C+
 * RunYieldMean("input.root", "output.root", "logs/", kInfo);
 * @endcode
 */
int RunYieldMean(
    TString inputFileName = "../RandomVars/260225_sysUncertainty_Total.root",             // Input file with spectra
    TString outputFileName = "../RandomVars/IntegratedYield/030625_IntYieldResults_v5.root", // Output file for results
    TString outputLogFolder = "../RandomVars/IntegratedYield/030625_IntYield_Logs_v5",       // Folder for YieldMean logs
    Int_t verbosity = kInfo)
{
    // --- Setup ---
    // Load and compile YieldMean.C and ExtraUtils.C using ACLiC if running as macro
    // If compiling, ensure they are linked properly.
    // gROOT->LoadMacro("YieldMean.C+"); // Might be needed if running interpreted
    // gROOT->LoadMacro("ExtraUtils.C+"); // Might be needed if running interpreted

    gErrorIgnoreLevel = verbosity;
    TH1::AddDirectory(kFALSE);   // Avoid histograms being added to current file/dir
    TH1::SetDefaultSumw2(kTRUE); // Ensure Sumw2 is set for histograms
    gStyle->SetOptFit(1);        // Show fit parameters on plot (useful for checks)
    gStyle->SetOptStat(0);       // No stats box
    SetCustomColorPalette();     // Apply custom colors from CascadeUtils.h

    TStopwatch totalTimer;
    totalTimer.Start();
    TDatime startTime;
    Info("RunYieldMean", "Analysis started at: %s", startTime.AsString());

    // Create output directory for logs if it doesn't exist
    outputLogFolder = SetOutputFolder(outputLogFolder);

    // --- Input File ---
    TFile *inputFile = OpenFile(inputFileName, "READ");
    if (!inputFile || inputFile->IsZombie())
    {
        Error("RunYieldMean", "Could not open input file: %s", inputFileName.Data());
        return 1;
    }
    Info("RunYieldMean", "Opened input file: %s", inputFileName.Data());

    // --- Output File ---
    TFile *outputFile = OpenFile(outputFileName, "RECREATE");
    if (!outputFile || outputFile->IsZombie())
    {
        Error("RunYieldMean", "Could not create output file: %s", outputFileName.Data());
        inputFile->Close();
        delete inputFile;
        return 1;
    }
    Info("RunYieldMean", "Created output file: %s", outputFileName.Data());

    // --- Define Particles and Fit Functions ---
    const int nParticles = kNumSignedPart; // kXip, kXim, kOmp, kOmm, kXiC, kOmC
    const char *particleName[nParticles] = {"xip", "xim", "omp", "omm", "xiC", "omC"};
    const char *particleTitle[nParticles] = {"#Xi^{+}", "#Xi^{-}", "#Omega^{+}", "#Omega^{-}", "#Xi^{#pm}", "#Omega^{#pm}"};
    Double_t particleMass[nParticles] = {fMass_Xi, fMass_Xi, fMass_Om, fMass_Om, fMass_Xi, fMass_Om};
    const int *nMultBins[nParticles] = {&fNmultbins_Xi, &fNmultbins_Xi, &fNmultbins_Om, &fNmultbins_Om, &fNmultbins_Xi, &fNmultbins_Om};
    const double *multBinEdges[nParticles] = {&fMultbins_Xi[0], &fMultbins_Xi[0], &fMultbins_Om[0], &fMultbins_Om[0], &fMultbins_Xi[0], &fMultbins_Om[0]};
    const int *nPtBins[nParticles] = {&fNptbins_Xi, &fNptbins_Xi, &fNptbins_Om, &fNptbins_Om, &fNptbins_Xi, &fNptbins_Om};
    const double *ptBinEdges[nParticles] = {&fPtbins_Xi[0], &fPtbins_Xi[0], &fPtbins_Om[0], &fPtbins_Om[0], &fPtbins_Xi[0], &fPtbins_Om[0]};

    // Array to hold pointers to the fit functions for each particle type
    TF1 *fitFunctions[nParticles][kNumFitFunctions];

    // Create fit function instances
    for (int iPart = 0; iPart < nParticles; ++iPart)
    {
        Double_t mass = particleMass[iPart];
        Double_t pTmin = ptBinEdges[iPart][0];
        Double_t pTmax = ptBinEdges[iPart][nPtBins[iPart][0]]; // Use actual max edge
        // Use slightly wider range for function definition than just bin edges
        Double_t funcRangeMin = 0.0;
        Double_t funcRangeMax = TMath::Max(8.0, pTmax + 2.0); // Extend beyond data

        // Initial parameters might need tuning per particle/multiplicity
        // For now, use the defaults from ExtraUtils.C
        // fitFunctions[iPart][kLevyTsallis] = LevyTsallis(Form("fLevyTsallis_%s", particleName[iPart]), mass);
        // fitFunctions[iPart][kBoltzmann] = Boltzmann(Form("fBoltzmann_%s", particleName[iPart]), mass);
        // fitFunctions[iPart][kBlastWave] = BGBlastWave(Form("fBlastWave_%s", particleName[iPart]), mass);
        // fitFunctions[iPart][kmTScaling] = mTScaling(Form("fmTScaling_%s", particleName[iPart]), mass);
        // fitFunctions[iPart][kBoseEinstein] = BoseEinstein(Form("fBoseEinstein_%s", particleName[iPart]), mass);
        // fitFunctions[iPart][kFermiDirac] = FermiDirac(Form("fFermiDirac_%s", particleName[iPart]), mass);

        fitFunctions[iPart][kLevyTsallis] = LevyTsallis("fLevyTsallis", mass);
        fitFunctions[iPart][kBoltzmann] = Boltzmann("fBoltzmann", mass);
        fitFunctions[iPart][kBlastWave] = BGBlastWave("fBlastWave", mass);
        // fitFunctions[iPart][kmTScaling] = mTScaling("fmTScaling", mass);
        // fitFunctions[iPart][kBoseEinstein] = BoseEinstein("fBoseEinstein", mass);
        // fitFunctions[iPart][kFermiDirac] = FermiDirac("fFermiDirac", mass);

        // Set range for all functions (can be overridden in YieldMean call if needed)
        for (int iFunc = 0; iFunc < kNumFitFunctions; ++iFunc)
        {
            if (fitFunctions[iPart][iFunc])
            {
                fitFunctions[iPart][iFunc]->SetRange(funcRangeMin, funcRangeMax);
            }
        }
    }

    // --- Main Loop: Particle -> Multiplicity -> Fit Function ---
    for (int iPart = 0; iPart < nParticles; ++iPart)
    {
        Info("RunYieldMean", "Processing particle: %s", particleTitle[iPart]);
        outputFile->cd();
        TDirectory *partDir = outputFile->mkdir(Form("YieldResults_%s", particleName[iPart]));
        if (!partDir)
        {
            Error("RunYieldMean", "Could not create directory for %s", particleName[iPart]);
            continue;
        }
        partDir->cd();

        // Define fit ranges based on particle type
        Double_t pTfitMin = ptBinEdges[iPart][0];               // Use first bin edge
        Double_t pTfitMax = ptBinEdges[iPart][*nPtBins[iPart]]; // Use last bin edge
        // Double_t pTfitMin = 0.0;
        // Double_t pTfitMax = 8.0;
        // Double_t pTfitMax = 5.0; // Or use a fixed value like in the thesis for consistency? Let's try 5.0 first. Check Fig 5.5/5.6 range.
        Double_t pTIntegralMin = 0.0;       // Integrate from 0
        Double_t pTIntegralMax = 8.0;       // Integrate up to 10 (or higher if needed)
        Double_t lowExtrapPrecision = 0.01; // As in original call
        Double_t highExtrapPrecision = 0.1; // As in original call
        TString fitOpt = "QRS";             // Quiet, Range, Minos (for better errors?), Save result

        // --- Process Multiplicity-Integrated ---
        TString histNameStat = Form("stat_effCorrPt_%s", particleName[iPart]);
        TString histNameSys = Form("sys_effCorrPt_%s", particleName[iPart]);

        TH1D *hStatInt = (TH1D *)inputFile->FindObjectAny(histNameStat.Data());
        TH1D *hSysInt = (TH1D *)inputFile->FindObjectAny(histNameSys.Data());

        if (!hStatInt || !hSysInt)
        {
            Warning("RunYieldMean", "Skipping multiplicity-integrated for %s - Histograms not found (%s or %s)",
                    particleName[iPart], histNameStat.Data(), histNameSys.Data());
        }
        else
        {
            TDirectory *intDir = partDir->mkdir("MultiplicityIntegrated");
            // intDir->cd();
            for (int iFunc = 0; iFunc < kNumFitFunctions; ++iFunc)
            {
                if (!fitFunctions[iPart][iFunc])
                {
                    Warning("RunYieldMean", "Fit function %d not defined for particle %s", iFunc, particleName[iPart]);
                    continue;
                }
                TString funcName = fitFunctions[iPart][iFunc]->GetName();
                Info("RunYieldMean", "  Fit Func: %s (Integrated)", funcName.Data());

                TString logFileName = Form("%s/Log_%s.root", outputLogFolder.Data(), funcName.Data());
                // Ensure the log file starts empty for each call if needed, or let YieldMean append
                // TFile::Open(logFileName, "RECREATE")->Close(); // Clear log file before run

                TH1 *yieldInfo = YieldMean(hStatInt, hSysInt, fitFunctions[iPart][iFunc],
                                           pTIntegralMin, pTIntegralMax,
                                           lowExtrapPrecision, highExtrapPrecision,
                                           fitOpt.Data(), logFileName.Data(),
                                           pTfitMin, pTfitMax);

                if (yieldInfo)
                {
                    yieldInfo->SetName(Form("YieldInfo_%s_Integrated_%s", particleName[iPart], funcName.Data()));
                    yieldInfo->SetTitle(Form("%s YieldInfo - Mult: 0-100%%, Fit: %s", particleTitle[iPart], funcName.Data()));
                    intDir->cd();
                    yieldInfo->Write();
                    delete yieldInfo; // Clean up memory if YieldMean allocates 'new'
                }
                else
                {
                    Warning("RunYieldMean", "YieldMean failed for %s, Integrated, Fit: %s", particleName[iPart], funcName.Data());
                }
            } // End loop over fit functions (Integrated)
        } // End multiplicity-integrated block

        // --- Process Multiplicity-Binned ---
        partDir->cd(); // Go back to particle directory
        TDirectory *multDir = partDir->mkdir("MultiplicityBinned");
        // multDir->cd();
        for (int iMult = 0; iMult < *nMultBins[iPart]; ++iMult)
        {
            Double_t multLow = multBinEdges[iPart][iMult];
            Double_t multHigh = multBinEdges[iPart][iMult + 1];
            Info("RunYieldMean", " Processing Multiplicity Bin: %.0f-%.0f %%", multLow, multHigh);

            histNameStat = Form("stat_effCorrPt_%s_mult[%d]", particleName[iPart], iMult);
            histNameSys = Form("sys_effCorrPt_%s_mult[%d]", particleName[iPart], iMult);

            TH1D *hStatMult = (TH1D *)inputFile->FindObjectAny(histNameStat.Data());
            TH1D *hSysMult = (TH1D *)inputFile->FindObjectAny(histNameSys.Data());

            if (!hStatMult || !hSysMult)
            {
                Warning("RunYieldMean", "Skipping bin %d for %s - Histograms not found (%s or %s)",
                        iMult, particleName[iPart], histNameStat.Data(), histNameSys.Data());
                continue;
            }

            // Optional: Adjust fit parameters based on multiplicity if needed
            // e.g., rough scaling of normalization? Or use previous bin's result? For now, keep defaults.

            for (int iFunc = 0; iFunc < kNumFitFunctions; ++iFunc)
            {
                if (!fitFunctions[iPart][iFunc])
                    continue; // Skip if function not defined
                TString funcName = fitFunctions[iPart][iFunc]->GetName();
                Info("RunYieldMean", "  Fit Func: %s (Mult Bin %d)", funcName.Data(), iMult);

                TString logFileName = Form("%s/Log_%s.root", outputLogFolder.Data(), funcName.Data());
                // TFile::Open(logFileName, "RECREATE")->Close(); // Clear log file before run

                // --- !! Crucial: Reset fit parameters before each fit !! ---
                // This is important if TF1 objects are reused across bins.
                // Re-create or reset parameters to defaults before fitting each bin.
                // Let's re-create TF1 for safety within the loop.
                // Note: This creates many TF1 objects. Managing memory or using Clone is an alternative.
                TF1 *currentFitFunc = nullptr;
                Double_t mass = particleMass[iPart];
                // Simple recreation based on type
                if (iFunc == kLevyTsallis)
                    currentFitFunc = LevyTsallis(Form("fLevyTsallis_%s_Mult%d", particleName[iPart], iMult), mass);
                else if (iFunc == kBoltzmann)
                    currentFitFunc = Boltzmann(Form("fBoltzmann_%s_Mult%d", particleName[iPart], iMult), mass);
                else if (iFunc == kBlastWave)
                    currentFitFunc = BGBlastWave(Form("fBlastWave_%s_Mult%d", particleName[iPart], iMult), mass);
                // else if (iFunc == kmTScaling)
                //     currentFitFunc = mTScaling(Form("fmTScaling_%s_Mult%d", particleName[iPart], iMult), mass);
                // else if (iFunc == kBoseEinstein)
                //     currentFitFunc = BoseEinstein(Form("fBoseEinstein_%s_Mult%d", particleName[iPart], iMult), mass);
                // else if (iFunc == kFermiDirac)
                //     currentFitFunc = FermiDirac(Form("fFermiDirac_%s_Mult%d", particleName[iPart], iMult), mass);

                if (!currentFitFunc)
                    continue;
                currentFitFunc->SetRange(pTIntegralMin, pTIntegralMax); // Ensure range is set

                // Adjust fit range for specific functions if needed (e.g., based on Chi2 target from thesis [cite: 59])
                // This logic would go here. For now, use the default pTfitMin, pTfitMax.

                TH1 *yieldInfo = YieldMean(hStatMult, hSysMult, currentFitFunc,
                                           pTIntegralMin, pTIntegralMax,
                                           lowExtrapPrecision, highExtrapPrecision,
                                           fitOpt.Data(), logFileName.Data(),
                                           pTfitMin, pTfitMax);

                if (yieldInfo)
                {
                    yieldInfo->SetName(Form("YieldInfo_%s_Mult%d_%s", particleName[iPart], iMult, funcName.Data()));
                    yieldInfo->SetTitle(Form("%s YieldInfo - Mult: %.0f-%.0f%%, Fit: %s", particleTitle[iPart], multLow, multHigh, funcName.Data()));
                    multDir->cd();
                    yieldInfo->Write();
                    delete yieldInfo;
                }
                else
                {
                    Warning("RunYieldMean", "YieldMean failed for %s, Mult %d, Fit: %s", particleName[iPart], iMult, funcName.Data());
                }
                delete currentFitFunc; // Clean up the temporary TF1
            } // End loop over fit functions (Binned)
        } // End loop over multiplicity bins
    } // End loop over particles

    // --- Cleanup ---
    inputFile->Close();
    outputFile->Close();
    delete inputFile;
    delete outputFile;

    // Clean up TF1 arrays (optional if pointers were managed differently)
    // for (int iPart = 0; iPart < nParticles; ++iPart) {
    //     for (int iFunc = 0; iFunc < kNumFitFunctions; ++iFunc) {
    //         //delete fitFunctions[iPart][iFunc]; // Be careful if these were just pointers to shared objects
    //     }
    // }

    totalTimer.Stop();
    TDatime endTime;
    Info("RunYieldMean", "Analysis finished at: %s", endTime.AsString());
    Info("RunYieldMean", "Total execution time: %.2f seconds (Real Time), %.2f seconds (CPU Time)",
         totalTimer.RealTime(), totalTimer.CpuTime());

    Info("RunYieldMean", "Results saved in: %s", outputFileName.Data());
    Info("RunYieldMean", "Fit logs saved in: %s", outputLogFolder.Data());
    Info("RunYieldMean", "You can now run the IntYieldsMeanPtCompare.C macro to generate the final plots.");

    return 0;
}