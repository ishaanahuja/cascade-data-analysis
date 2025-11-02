/**
 * @file IntYieldsMeanPtCompare.C
 * @author Ishaan Ahuja (ishaanahuja0@gmail.com)
 * @version 1
 * @brief Macro for generating final integrated yield and mean pT plots with fit overlays and ratio plots.
 * @date 13-06-2025
 *
 * @note Dirty and rushed, but works. Should be cleaned up, optimized and modularized.
 *
 * @copyright Copyright (c) 2025
 *
 */

#include <TROOT.h>
#include "TFile.h"
#include "TCanvas.h"
#include "TH1.h"
#include "TF1.h"
#include "TGraphErrors.h"
#include "TLegend.h"
#include "TPad.h"
#include "TStyle.h"
#include "TMath.h"
#include "TError.h"
#include "TString.h"
#include "TLatex.h"
#include "TList.h" // Needed for GetListOfFunctions
#include "TLine.h"
#include "TGraphAsymmErrors.h" // Needed for systematic error boxes
#include "CascadeUtils.h"      // Contains constants, enums, colors, dNch etc.
#include "ExtraUtils.C"        // Contains fit function definitions (enum needed)

// Modified helper function to get individual error components
void GetIndividualErrors(TH1 *hYieldInfo, EValue_t valBin, EValue_t statBin, EValue_t sysLoBin, EValue_t sysHiBin,
                         Double_t &val, Double_t &statError, Double_t &sysErrorLow, Double_t &sysErrorHigh)
{
    if (!hYieldInfo)
    {
        val = 0;
        statError = 0;
        sysErrorLow = 0;
        sysErrorHigh = 0;
        return;
    }
    val = hYieldInfo->GetBinContent(valBin);
    statError = hYieldInfo->GetBinContent(statBin);
    sysErrorLow = hYieldInfo->GetBinContent(sysLoBin);
    sysErrorHigh = hYieldInfo->GetBinContent(sysHiBin);
}

// Helper function to create Data/Fit ratio plot using a retrieved TF1
TH1 *CreateRatioPlot(TH1 *hData, TF1 *fFit, Double_t ratioMin = 0.5, Double_t ratioMax = 1.5)
{
    if (!hData || !fFit)
        return nullptr;

    TString ratioName = TString::Format("%s_Ratio", hData->GetName());
    TH1 *hRatio = (TH1 *)hData->Clone(ratioName.Data());
    hRatio->Reset(); // Clear contents but keep bins
    hRatio->GetYaxis()->SetTitle("Data / Fit");
    hRatio->GetYaxis()->SetRangeUser(ratioMin, ratioMax);
    hRatio->GetYaxis()->SetNdivisions(505);                      // Standard for ratio plots
    hRatio->GetXaxis()->SetTitle(hData->GetXaxis()->GetTitle()); // Keep x-axis title
    hRatio->SetStats(0);                                         // No stats box
    hRatio->SetMarkerStyle(hData->GetMarkerStyle());
    hRatio->SetMarkerSize(hData->GetMarkerSize());
    hRatio->SetMarkerColor(hData->GetMarkerColor());
    hRatio->SetLineColor(hData->GetLineColor());

    for (int iBin = 1; iBin <= hData->GetNbinsX(); ++iBin)
    {
        Double_t dataVal = hData->GetBinContent(iBin);
        Double_t dataErr = hData->GetBinError(iBin); // Use total error for ratio error
        Double_t binCenter = hData->GetBinCenter(iBin);
        Double_t fitVal = fFit->Eval(binCenter);

        if (dataVal != 0 && fitVal != 0 && TMath::Abs(fitVal) > 1e-15) // Avoid division by zero or tiny numbers
        {
            Double_t ratio = dataVal / fitVal;
            Double_t ratioErr = dataErr / TMath::Abs(fitVal); // Basic error propagation
            hRatio->SetBinContent(iBin, ratio);
            hRatio->SetBinError(iBin, ratioErr);
        }
        else
        {
            hRatio->SetBinContent(iBin, 0); // Set to 0 if data or fit is zero
            hRatio->SetBinError(iBin, 0);
        }
    }
    return hRatio;
}
// Modified CreateRatioPlot to take TGraphErrors for data
TH1 *CreateRatioPlotFromGraph(TGraphErrors *gData, TF1 *fFit, TH1 *hDataAxisDef, Double_t ratioMin = 0.6, Double_t ratioMax = 1.4)
{
    if (!gData || !fFit || !hDataAxisDef)
        return nullptr;

    hDataAxisDef->SetTitle(""); // Clear title to avoid confusion
    TString ratioName = TString::Format("%s_Ratio", gData->GetName());
    // Clone axis definition from hDataAxisDef (which has correct binning for pT)
    TH1 *hRatio = (TH1 *)hDataAxisDef->Clone(ratioName.Data());
    hRatio->Reset(); // Clear contents but keep bins and axes
    hRatio->GetYaxis()->SetTitle("Ratio to Fit");
    hRatio->GetYaxis()->SetRangeUser(ratioMin, ratioMax);
    // hRatio->GetYaxis()->SetNdivisions(505);
    hRatio->GetXaxis()->SetTitle(hDataAxisDef->GetXaxis()->GetTitle()); // Keep x-axis title
    hRatio->SetStats(0);

    for (int iPoint = 0; iPoint < gData->GetN(); ++iPoint)
    {
        Double_t pt, dataVal, ptErr, dataErr;
        gData->GetPoint(iPoint, pt, dataVal);
        dataErr = gData->GetErrorY(iPoint); // Get Y stat error
        Double_t fitVal = fFit->Eval(pt);

        Int_t iBin = hRatio->GetXaxis()->FindBin(pt);

        if (dataVal != 0 && fitVal != 0 && TMath::Abs(fitVal) > 1e-15)
        {
            Double_t ratio = dataVal / fitVal;
            Double_t ratioErr = dataErr / TMath::Abs(fitVal);
            hRatio->SetBinContent(iBin, ratio); // Use bin corresponding to pt
            hRatio->SetBinError(iBin, ratioErr);
        }
        else
        {
            hRatio->SetBinContent(iBin, 0);
            hRatio->SetBinError(iBin, 0);
        }
    }
    // Style for points on ratio plot
    hRatio->SetMarkerStyle(gData->GetMarkerStyle());
    hRatio->SetMarkerSize(gData->GetMarkerSize());
    hRatio->SetMarkerColor(gData->GetMarkerColor());
    hRatio->SetLineColor(gData->GetLineColor());
    return hRatio;
}

// Helper to get the standard name for the TF1 based on type, particle, mult
TString GetFitFunctionName(int iFunc, const char *particleName, int iMult)
{
    TString funcName = "";
    if (iFunc == kLevyTsallis)
        funcName = "fLevyTsallis";
    else if (iFunc == kBoltzmann)
        funcName = "fBoltzmann";
    else if (iFunc == kBlastWave)
        funcName = "fBlastWave";
    // else if (iFunc == kmTScaling)
    //     funcName = "fmTScaling";
    // else if (iFunc == kBoseEinstein)
    //     funcName = "fBoseEinstein";
    // else if (iFunc == kFermiDirac)
    //     funcName = "fFermiDirac";
    // return Form("%s_%s_Mult%d", funcName.Data(), particleName, iMult);
    return funcName;
}

// Helper to get the name of the histogram saved in the log file
TString GetLogHistogramName(const char *baseHistName, const char *funcName, const char *partName, int mult)
{
    return Form("%s_%s_%s_Mult%d", baseHistName, funcName, partName, mult);
}

//-----------------------------------------------------------------------------
// Main Plotting Function
//-----------------------------------------------------------------------------
void IntYieldsMeanPtCompare(
    TString inputSpectraFile = "/var/home/ishaan/Work/git/analysis/results/RandomVars/040625_SysUncertainty_Total_final_v2/040625_sysUncertainty_Total_final_v2.root", // Original spectra
    TString resultsFileName = "/var/home/ishaan/Work/git/analysis/results/RandomVars/IntYieldFinal/090625_IntYieldResults.root",                                       // File from RunYieldMean
    TString inputLogFolder = "/var/home/ishaan/Work/git/analysis/results/RandomVars/IntYieldFinal/090625_IntYield_Logs",                                               //// Folder WITH LOGS from RunYieldMean
    TString outputPlotFolder = "/var/home/ishaan/Work/git/analysis/results/RandomVars/IntYieldFinal/120625_Plots_ThisWork",                                            //
    TString imageFormat = "pdf")
{
    // --- Setup ---
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(0);    // Don't show fit parameters box automatically
    SetCustomColorPalette(); // Apply custom colors
    gStyle->SetLineScalePS(1.5);
    gStyle->SetStatFontSize(0.03);
    gStyle->SetPadTickX(1); // Ticks on both top and bottom for X axis
    gStyle->SetPadTickY(1); // Ticks on both left and right for Y axis

    outputPlotFolder = SetOutputFolder(outputPlotFolder); // Create output folder

    // --- Open Files ---
    TFile *resultsFile = OpenFile(resultsFileName, "READ");
    TFile *spectraFile = OpenFile(inputSpectraFile, "READ");

    const int nParticlesPlot = 2; // Plotting combined Xi and Omega
    const char *particleNamePlot[nParticlesPlot] = {"xiC", "omC"};
    const char *particleTitlePlot[nParticlesPlot] = {"#Xi^{+} + #Xi^{-}", "#Omega^{+} + #Omega^{-}"};
    const char *particleDirName[nParticlesPlot] = {"YieldResults_xiC", "YieldResults_omC"}; // Dir names in results file
    const double *dNch[nParticlesPlot] = {dNchXi, dNchOm};
    const double *dNchErr[nParticlesPlot] = {dNchXiErr, dNchOmErr};
    const int *nMultBinsPlot[nParticlesPlot] = {&fNmultbins_Xi, &fNmultbins_Om};
    const double *multBinEdgesPlot[nParticlesPlot] = {&fMultbins_Xi[0], &fMultbins_Om[0]};

    // Store results from all fit functions for all particles and multiplicities
    std::vector<std::vector<double>> yieldsPerFuncAll(nParticlesPlot);
    std::vector<std::vector<double>> meansPerFuncAll(nParticlesPlot);
    for (int iPart = 0; iPart < nParticlesPlot; ++iPart)
    {
        yieldsPerFuncAll[iPart].resize((*nMultBinsPlot[iPart]) * kNumFitFunctions);
        meansPerFuncAll[iPart].resize((*nMultBinsPlot[iPart]) * kNumFitFunctions);
        // Initialize with a value indicating data not found or failed
        std::fill(yieldsPerFuncAll[iPart].begin(), yieldsPerFuncAll[iPart].end(), -999.0);
        std::fill(meansPerFuncAll[iPart].begin(), meansPerFuncAll[iPart].end(), -999.0);
    }
    // Populate yieldsPerFuncAll and meansPerFuncAll first
    Info("IntYieldsMeanPtCompare", "Collecting results from all fit functions for extrapolation uncertainty...");
    for (int iPart = 0; iPart < nParticlesPlot; ++iPart)
    {
        int nMult = *nMultBinsPlot[iPart];
        for (int iMult = 0; iMult < nMult; ++iMult)
        {
            for (int iFunc = 0; iFunc < kNumFitFunctions; ++iFunc)
            {
                TString funcName = GetFitFunctionName(iFunc, particleNamePlot[iPart], iMult);
                TString yieldInfoNameFunc = Form("%s/MultiplicityBinned/YieldInfo_%s_Mult%d_%s",
                                                 particleDirName[iPart], particleNamePlot[iPart], iMult, funcName.Data());
                TH1 *hYieldInfoFunc = (TH1 *)resultsFile->Get(yieldInfoNameFunc.Data());
                int index = iMult * kNumFitFunctions + iFunc;
                if (hYieldInfoFunc)
                {
                    yieldsPerFuncAll[iPart][index] = hYieldInfoFunc->GetBinContent(kYield);
                    meansPerFuncAll[iPart][index] = hYieldInfoFunc->GetBinContent(kMean);
                }
                else
                {
                    // Already initialized to -999, just a warning
                    Warning("IntYieldsMeanPtCompare", "YieldInfo not found (will skip for extrap): %s", yieldInfoNameFunc.Data());
                }
            }
        }
    }

    // Calculate and store extrapolation systematic uncertainties
    std::vector<std::vector<double>> sysExtrapYield(nParticlesPlot);
    std::vector<std::vector<double>> sysExtrapMeanPt(nParticlesPlot);
    for (int iPart = 0; iPart < nParticlesPlot; ++iPart)
    {
        sysExtrapYield[iPart].resize(*nMultBinsPlot[iPart], 0.0);
        sysExtrapMeanPt[iPart].resize(*nMultBinsPlot[iPart], 0.0);
    }

    Info("IntYieldsMeanPtCompare", "Calculating extrapolation systematic uncertainties...");
    for (int iPart = 0; iPart < nParticlesPlot; ++iPart)
    {
        int nMult = *nMultBinsPlot[iPart];
        for (int iMult = 0; iMult < nMult; ++iMult)
        {
            std::vector<double> currentYields;
            std::vector<double> currentMeans;
            double levyYield = -999.0;
            double levyMean = -999.0;

            int levyIndex = iMult * kNumFitFunctions + kLevyTsallis;
            if (levyIndex < yieldsPerFuncAll[iPart].size())
            { // bounds check
                levyYield = yieldsPerFuncAll[iPart][levyIndex];
                levyMean = meansPerFuncAll[iPart][levyIndex];
            }

            if (levyYield <= -998.0)
            { // Check if Levy result was valid
                Warning("IntYieldsMeanPtCompare", "Levy-Tsallis result not found or invalid for %s mult %d. Skipping extrap. uncertainty.", particleNamePlot[iPart], iMult);
                continue;
            }

            double maxYieldDev = 0.0;
            double maxMeanDev = 0.0;

            for (int iFunc = 0; iFunc < kNumFitFunctions; ++iFunc)
            {
                if (iFunc == kLevyTsallis)
                    continue; // Compare alternatives to Levy

                int funcIndex = iMult * kNumFitFunctions + iFunc;
                double funcYield = -999.0;
                double funcMean = -999.0;

                if (funcIndex < yieldsPerFuncAll[iPart].size())
                { // bounds check
                    funcYield = yieldsPerFuncAll[iPart][funcIndex];
                    funcMean = meansPerFuncAll[iPart][funcIndex];
                }

                if (funcYield > -998.0)
                { // Check if alternative function result was valid
                    double devY = TMath::Abs(funcYield - levyYield);
                    if (devY > maxYieldDev)
                        maxYieldDev = devY;
                }
                if (funcMean > -998.0)
                {
                    double devM = TMath::Abs(funcMean - levyMean);
                    if (devM > maxMeanDev)
                        maxMeanDev = devM;
                }
            }
            sysExtrapYield[iPart][iMult] = 0.5 * maxYieldDev;
            sysExtrapMeanPt[iPart][iMult] = 0.5 * maxMeanDev;
        }
    }

    // Output Final Values to Console in a Formatted Table
    Info("IntYieldsMeanPtCompare", "--- Generating Summary Tables for Final Results ---");

    for (int iPart = 0; iPart < nParticlesPlot; ++iPart)
    {
        // Print table header for the current particle
        printf("\n\n");
        printf("========================================================================================================================================\n");
        printf("                                         Summary Table for %s (%s)                                           \n", particleTitlePlot[iPart], particleNamePlot[iPart]);
        printf("========================================================================================================================================\n");
        printf("%-18s | %-22s | %-45s | %-42s | %-15s\n", "Multiplicity Class", "<dNch/deta>", "dN/dy (Yield)", "<pT> (GeV/c)", "Extrap. (%)");
        printf("----------------------------------------------------------------------------------------------------------------------------------------\n");

        int nMult = *nMultBinsPlot[iPart];
        for (int iMult = 0; iMult < nMult; ++iMult)
        {
            // Get the primary results calculated with the default Levy-Tsallis fit
            TString levyFitName = GetFitFunctionName(kLevyTsallis, particleNamePlot[iPart], iMult);
            TString yieldInfoName = Form("%s/MultiplicityBinned/YieldInfo_%s_Mult%d_%s",
                                         particleDirName[iPart], particleNamePlot[iPart], iMult, levyFitName.Data());
            TH1 *hYieldInfo = (TH1 *)resultsFile->Get(yieldInfoName.Data());

            if (!hYieldInfo)
            {
                Warning("SummaryTable", "YieldInfo histogram not found for %s in mult bin %d. Skipping table row.", particleNamePlot[iPart], iMult);
                continue;
            }

            // --- 1. Get Multiplicity Class ---
            TString multClassStr = Form("%.0f-%.0f %%", multBinEdgesPlot[iPart][iMult], multBinEdgesPlot[iPart][iMult + 1]);

            // --- 2. Get <dNch/deta> ---
            TString dNchStr = Form("%.2f #pm %.2f", dNch[iPart][iMult], dNchErr[iPart][iMult]);

            // --- 3. Get Yield and Combine Uncertainties ---
            double yieldVal = hYieldInfo->GetBinContent(kYield);
            double yieldStatErr = hYieldInfo->GetBinContent(kYieldStat);
            // Combine systematic errors: from spectra and from extrapolation choice
            double yieldSysSpectra = (hYieldInfo->GetBinContent(kYieldSysHi) + hYieldInfo->GetBinContent(kYieldSysLo)) / 2.0;
            // CORRECTED: Use the absolute uncertainty from sysExtrapYield directly
            double yieldSysExtrap = sysExtrapYield[iPart][iMult];
            double totalYieldSysErr = TMath::Sqrt(TMath::Power(yieldSysSpectra, 2) + TMath::Power(yieldSysExtrap, 2));
            TString yieldStr = Form("%.4f #pm %.4f (stat) #pm %.4f (sys)", yieldVal, yieldStatErr, totalYieldSysErr);

            // --- 4. Get Mean pT and Combine Uncertainties ---
            double meanPtVal = hYieldInfo->GetBinContent(kMean);
            double meanPtStatErr = hYieldInfo->GetBinContent(kMeanStat);
            // Combine systematic errors similarly
            double meanPtSysSpectra = (hYieldInfo->GetBinContent(kMeanSysHi) + hYieldInfo->GetBinContent(kMeanSysLo)) / 2.0;
            // CORRECTED: Use the absolute uncertainty from sysExtrapMeanPt directly
            double meanPtSysExtrap = sysExtrapMeanPt[iPart][iMult];
            double totalMeanPtSysErr = TMath::Sqrt(TMath::Power(meanPtSysSpectra, 2) + TMath::Power(meanPtSysExtrap, 2));
            TString meanPtStr = Form("%.3f #pm %.3f (stat) #pm %.3f (sys)", meanPtVal, meanPtStatErr, totalMeanPtSysErr);

            // --- 5. Get Extrapolation Percentage ---
            double extrapVal = hYieldInfo->GetBinContent(kExtra);
            double totalYieldForFrac = hYieldInfo->GetBinContent(kYield);
            double extrapPercentage = (totalYieldForFrac > 1e-9) ? (extrapVal / totalYieldForFrac) * 100.0 : 0.0;
            TString extrapStr = Form("%.1f", extrapPercentage);

            // --- Print the fully formatted table row ---
            printf("%-18s | %-22s | %-45s | %-42s | %-15s\n",
                   multClassStr.Data(), dNchStr.Data(), yieldStr.Data(), meanPtStr.Data(), extrapStr.Data());
        }
        printf("========================================================================================================================================\n");
    }
    printf("\n"); // Add a final newline for clean exit
    // --- END OF SUMMARY TABLE SECTION ---

    // --- Generate Yield vs Multiplicity ---
    TCanvas *cYield = new TCanvas("cYieldVsMult", "Yield vs dNch/deta", 1200, 800);
    // cYield->Divide(1, 2);                             // One plot for Xi, one for Omega
    // ** Define TGraphErrors for stat and TGraphAsymmErrors for sys **
    TGraphErrors *gYieldStat[nParticlesPlot] = {nullptr};
    TGraphAsymmErrors *gYieldSys[nParticlesPlot] = {nullptr};
    TF1 *linearYieldFits[nParticlesPlot] = {nullptr}; // To store linear fits to yield
    // TFitResultPtr fitResultYield[nParticlesPlot] = {nullptr}; // Fit results for yield fits

    // TGraphErrors *gYield[nParticlesPlot] = {nullptr}; // Initialize to nullptr
    // TLegend *legYield = new TLegend(0.6, 0.7, 0.88, 0.88);
    // legYield->SetBorderSize(0);
    // legYield->SetFillStyle(0);

    for (int iPart = 0; iPart < nParticlesPlot; ++iPart)
    {
        int nMult = *nMultBinsPlot[iPart];
        std::vector<double> v_dNch, v_dNchErrX, v_Yield;
        std::vector<double> v_YieldStatErr, v_TotalYieldSysLo, v_TotalYieldSysHi; // Total sys
        std::vector<double> v_dNchErrSysXLo, v_dNchErrSysXHi;                     // For x-asymmetry of sys boxes (using dNchErr)

        cYield->Clear();

        for (int iMult = 0; iMult < nMult; ++iMult)
        {
            // Get results for the default Levy-Tsallis fit function name generated by RunYieldMean
            TString funcNameLevy = GetFitFunctionName(kLevyTsallis, particleNamePlot[iPart], iMult);
            TString yieldInfoName = Form("%s/MultiplicityBinned/YieldInfo_%s_Mult%d_%s",
                                         particleDirName[iPart], particleNamePlot[iPart], iMult, funcNameLevy.Data());
            TH1 *hYieldInfo = (TH1 *)resultsFile->Get(yieldInfoName.Data());
            if (!hYieldInfo)
            {
                Warning("IntYieldsMeanPtCompare", "YieldInfo histogram not found for Yield plot: %s", yieldInfoName.Data());
                continue;
            }

            Double_t yieldVal, statErrVal, sysLoSpectra, sysHiSpectra;
            // GetCombinedErrors(hYieldInfo, kYield, kYieldStat, kYieldSysLo, kYieldSysHi, yield, totalErr, statErr, sysLo, sysHi);
            GetIndividualErrors(hYieldInfo, kYield, kYieldStat, kYieldSysLo, kYieldSysHi,
                                yieldVal, statErrVal, sysLoSpectra, sysHiSpectra);

            // Combine sys errors for Yield
            Double_t currentSysExtrapYield = sysExtrapYield[iPart][iMult];
            Double_t totalSysLo = TMath::Sqrt(TMath::Power(sysLoSpectra, 2) + TMath::Power(currentSysExtrapYield, 2));
            Double_t totalSysHi = TMath::Sqrt(TMath::Power(sysHiSpectra, 2) + TMath::Power(currentSysExtrapYield, 2));

            v_dNch.push_back(dNch[iPart][iMult]);
            v_dNchErrX.push_back(dNchErr[iPart][iMult]); // Store dNch error for stat graph x-error
            v_Yield.push_back(yieldVal);
            // Store separate errors
            v_YieldStatErr.push_back(statErrVal);
            // Store combined sys errors
            v_TotalYieldSysLo.push_back(totalSysLo);
            v_TotalYieldSysHi.push_back(totalSysHi);
            v_dNchErrSysXLo.push_back(dNchErr[iPart][iMult]); // Symmetric x-error for boxes
            v_dNchErrSysXHi.push_back(dNchErr[iPart][iMult]);

            // // Store results for all functions for extrapolation uncertainty later
            // for (int iFunc = 0; iFunc < kNumFitFunctions; ++iFunc)
            // {
            //     TString funcName = GetFitFunctionName(iFunc, particleNamePlot[iPart], iMult);
            //     TString yieldInfoNameFunc = Form("%s/MultiplicityBinned/YieldInfo_%s_Mult%d_%s",
            //                                      particleDirName[iPart], particleNamePlot[iPart], iMult, funcName.Data());
            //     TH1 *hYieldInfoFunc = (TH1 *)resultsFile->Get(yieldInfoNameFunc.Data());
            //     int index = iMult * kNumFitFunctions + iFunc;
            //     if (hYieldInfoFunc)
            //     {
            //         yieldsPerFunc[iPart][index] = hYieldInfoFunc->GetBinContent(kYield);
            //         meansPerFunc[iPart][index] = hYieldInfoFunc->GetBinContent(kMean);
            //     }
            //     else
            //     {
            //         yieldsPerFunc[iPart][index] = -1; // Mark as invalid
            //         meansPerFunc[iPart][index] = -1;
            //         Warning("IntYieldsMeanPtCompare", "Missing YieldInfo for extrap uncertainty: %s", yieldInfoNameFunc.Data());
            //     }
            // }
        } // end multiplicity loop

        if (v_dNch.empty())
            continue;

        gYieldStat[iPart] = new TGraphErrors(v_dNch.size(), v_dNch.data(), v_Yield.data(), v_dNchErrX.data(), v_YieldStatErr.data());
        // Use combined sys errors for TGraphAsymmErrors **
        gYieldSys[iPart] = new TGraphAsymmErrors(v_dNch.size(), v_dNch.data(), v_Yield.data(),
                                                 v_dNchErrSysXLo.data(), v_dNchErrSysXHi.data(),
                                                 v_TotalYieldSysLo.data(), v_TotalYieldSysHi.data());
        // gYield[iPart] = new TGraphErrors(v_dNch.size(), v_dNch.data(), v_Yield.data(), v_dNchErr.data(), v_YieldErr.data());
        gYieldStat[iPart]->SetName(Form("gYield_%s", particleNamePlot[iPart]));
        gYieldStat[iPart]->SetTitle(TString::Format("%s : Integrated Yield", particleTitlePlot[iPart]));
        gYieldStat[iPart]->SetMarkerStyle(kFullCircle);              // Use styles from CascadeUtils
        gYieldStat[iPart]->SetMarkerColor(customCascPalette[iPart]); // Use colors from customCascPalette
        gYieldStat[iPart]->SetLineColor(customCascPalette[iPart]);
        // gYieldStat[iPart]->SetMarkerSize(1.5);
        gYieldStat[iPart]->GetYaxis()->SetTitle(Form("%s yield", particleTitlePlot[iPart]));

        gYieldSys[iPart]->SetName(Form("gYieldSys_%s", particleNamePlot[iPart]));
        gYieldSys[iPart]->SetFillColor(customCascPalette[iPart]); // Semi-transparent fill
        gYieldSys[iPart]->SetLineColor(customCascPalette[iPart]); // Box outline color
        gYieldSys[iPart]->SetFillStyle(0);                        // Solid fill, or e.g. 3002 for hatches

        cYield->cd(); // Go to the correct pad (1 for Xi, 2 for Omega)
                      // gPad->SetLeftMargin(0.15);
                      // gPad->SetBottomMargin(0.12);

        // double minDNch = TMath::MinElement(v_dNch.size(), v_dNch.data()) - 5;
        // if (minDNch < 0)
        //     minDNch = 0;
        // double maxDNch = TMath::MaxElement(v_dNch.size(), v_dNch.data()) + 5;
        // double minYield = TMath::MinElement(v_Yield.size(), v_Yield.data());
        // double maxYield = TMath::MaxElement(v_Yield.size(), v_Yield.data());
        // minYield *= (minYield > 0 ? 0.8 : 1.2);
        // if (minYield == 0 && iPart == 0)
        //     minYield = -0.01 * maxYield;
        // if (minYield == 0 && iPart == 1)
        //     minYield = -0.001 * maxYield; // Small offset if min is 0
        // maxYield *= 1.2;

        double minDNch = 0; // Start from 0 for multiplicity
        // double minDNch = TMath::MinElement(v_dNch.size(), v_dNch.data()) - 10;
        // if (minDNch < 0)
        //     minDNch = 0;
        double maxDNch = TMath::MaxElement(v_dNch.size(), v_dNch.data()) + 10;
        double minYieldVal = 0; // Start from 0 for yield
        // double minYieldVal = TMath::MinElement(v_Yield.size(), v_Yield.data());
        double maxYieldVal = TMath::MaxElement(v_Yield.size(), v_Yield.data());
        // Adjust range considering errors
        for (size_t k = 0; k < v_Yield.size(); ++k)
        {
            // if (v_Yield[k] - v_TotalYieldSysLo[k] < minYieldVal)
            // minYieldVal = v_Yield[k] - v_TotalYieldSysLo[k];
            if (v_Yield[k] + v_TotalYieldSysHi[k] > maxYieldVal)
                maxYieldVal = v_Yield[k] + v_TotalYieldSysHi[k];
        }
        // minYieldVal *= (minYieldVal > 0 ? 0.8 : 1.2);
        // if (minYieldVal > -1E-9 && minYieldVal < 1E-9 && iPart == 0)
        //     minYieldVal = -0.01 * maxYieldVal;
        // if (minYieldVal > -1E-9 && minYieldVal < 1E-9 && iPart == 1)
        //     minYieldVal = -0.001 * maxYieldVal;
        maxYieldVal *= 1.2;
        if (maxYieldVal == 0)
            maxYieldVal = 1; // Avoid zero range
        TH1 *hFrameYield = gPad->DrawFrame(minDNch, minYieldVal, maxDNch, maxYieldVal, ";#LTd#it{N}_{ch}/d#it{#eta}#GT_{|#it{#eta}|<0.5};#LTd#it{N}/d#it{y}#GT");
        // hFrameYield->GetYaxis()->SetTitleOffset(1.8);
        hFrameYield->GetXaxis()->SetTitleOffset(1.3);

        // gYield[iPart]->Draw("P SAME"); // Draw points with errors

        // Draw sys boxes first, then stat errors/points
        if (gYieldSys[iPart])
            gYieldSys[iPart]->Draw("E2 SAME"); // "E2" draws a box with fill color
        if (gYieldStat[iPart])
            gYieldStat[iPart]->Draw("PE1 SAME"); // "P" draws points, "PE" would draw stat error bars if needed over box lines

        gPad->SetRightMargin(0.05);
        gPad->SetTopMargin(0.05);
        gPad->SetBottomMargin(0.12);
        TLatex latex;
        latex.SetNDC();
        latex.SetTextSize(0.025);
        // TLatex *latexInfo = new TLatex(0.65, 0.81, "Uncertainties: stat.(bars), sys.(boxes)");
        latex.DrawLatex(0.15, 0.86, "p-Pb #sqrt{s_{NN}} = 8.16 TeV, This work");
        latex.DrawLatex(0.15, 0.82, Form("%s, -0.5 < y < 0", particleTitlePlot[iPart]));
        latex.SetTextSize(0.022);
        latex.SetTextFont(42);
        latex.DrawLatex(0.15, 0.78, "Uncertainties: stat.(bars), sys.(boxes)");
        cYield->Update();
        SaveImage(outputPlotFolder, "YieldVsMult", TString::Format("YieldVsMult_%s", particleNamePlot[iPart]), imageFormat.Data(), cYield);
    } // end particle loop
    delete cYield; // Clean up canvas
    // cYield->cd(); // Draw legend on the first pad
    // // legYield->Draw();
    // TLatex latex;
    // latex.SetNDC();
    // latex.SetTextSize(0.03);
    // latex.DrawLatex(0.2, 0.80, "p-Pb #sqrt{s_{NN}} = 8.16 TeV");
    // latex.DrawLatex(0.2, 0.75, "-0.5 < y < 0");
    // cYield->Update();
    // SaveImage(outputPlotFolder, "YieldVsMult", TString::Format("YieldVsMult_%s", particleNamePlot[iPart]), imageFormat.Data(), cYield);

    // --- Generate Yield vs Multiplicity with Linear Fits ---
    TCanvas *cYieldWithFit = new TCanvas("cYieldVsMultWithFit", "Yield vs dNch/deta with Linear Fit", 1200, 800);

    for (int iPart = 0; iPart < nParticlesPlot; ++iPart)
    {
        cYieldWithFit->Clear(); // Clear canvas for each particle plot

        if (!gYieldStat[iPart] || !gYieldSys[iPart] || gYieldStat[iPart]->GetN() == 0)
        {
            Warning("IntYieldsMeanPtCompare", "gYieldStat/Sys for particle %d not populated. Skipping yield plot with fit.", iPart);
            continue;
        }

        // Dynamic frame calculation
        double minDNch_plot = 0;
        double maxDNch_plot = 60; // Default
        if (gYieldStat[iPart]->GetN() > 0)
        {
            minDNch_plot = TMath::MinElement(gYieldStat[iPart]->GetN(), gYieldStat[iPart]->GetX());
            maxDNch_plot = TMath::MaxElement(gYieldStat[iPart]->GetN(), gYieldStat[iPart]->GetX());
            minDNch_plot = TMath::Max(0.0, minDNch_plot - 5); // Padding
            maxDNch_plot = maxDNch_plot + 5;                  // Padding
        }

        double minYield_plot = 0;
        double maxYield_plot = 0.1; // Default
        if (gYieldStat[iPart]->GetN() > 0)
        {
            minYield_plot = TMath::MinElement(gYieldStat[iPart]->GetN(), gYieldStat[iPart]->GetY());
            maxYield_plot = TMath::MaxElement(gYieldStat[iPart]->GetN(), gYieldStat[iPart]->GetY());
            for (int k = 0; k < gYieldStat[iPart]->GetN(); ++k)
            {
                minYield_plot = TMath::Min(minYield_plot, gYieldStat[iPart]->GetY()[k] - gYieldSys[iPart]->GetErrorYlow(k));
                maxYield_plot = TMath::Max(maxYield_plot, gYieldStat[iPart]->GetY()[k] + gYieldSys[iPart]->GetErrorYhigh(k));
            }
        }
        double yieldPadding_plot = (maxYield_plot - minYield_plot) * 0.15; // 15% padding
        if (yieldPadding_plot == 0)
            yieldPadding_plot = TMath::Max(0.01, TMath::Abs(maxYield_plot * 0.2));

        double final_min_yield = minYield_plot - yieldPadding_plot;
        // Ensure y-axis starts at or below 0 if data is close to 0
        if (minYield_plot >= 0 && final_min_yield > 0 && minYield_plot < 0.1 * (maxYield_plot - minYield_plot))
        {
            final_min_yield = TMath::Min(0.0, final_min_yield);
        }
        if (final_min_yield > 0 && minYield_plot >= 0)
        { // If still positive, make a small negative offset for clarity if near zero
            if (iPart == 0 /*&& particleNamePlot[iPart] == "xiC"*/)
                final_min_yield = -0.005;
            else if (iPart == 1 /*&& particleNamePlot[iPart] == "omC"*/)
                final_min_yield = -0.0005;
        }

        double final_max_yield = maxYield_plot + yieldPadding_plot;
        if (final_max_yield <= final_min_yield)
            final_max_yield = final_min_yield + yieldPadding_plot + 0.01;

        TH1 *hFrameYieldFit = gPad->DrawFrame(minDNch_plot, final_min_yield, maxDNch_plot, final_max_yield, ";#LTd#it{N}_{ch}/d#it{#eta}#GT_{|#it{#eta}|<0.5};#LTd#it{N}/d#it{y}#GT");
        // hFrameYieldFit->GetYaxis()->SetTitleOffset(1.4);
        hFrameYieldFit->GetXaxis()->SetTitleOffset(1.3);
        gPad->SetRightMargin(0.05);
        gPad->SetTopMargin(0.05);
        gPad->SetBottomMargin(0.12);
        // gPad->SetLeftMargin(0.15);

        // Draw existing stat and sys graphs
        if (gYieldSys[iPart])
            gYieldSys[iPart]->Draw("E2 SAME");
        if (gYieldStat[iPart])
            gYieldStat[iPart]->Draw("PE SAME");

        // Fit gYieldStat[iPart] with a linear function
        if (gYieldStat[iPart] && gYieldStat[iPart]->GetN() > 1)
        { // Need at least 2 points
            TString fitName = Form("fitYield_%s", particleNamePlot[iPart]);
            double dNchMin = gYieldStat[iPart]->GetX()[0];
            double dNchMax = gYieldStat[iPart]->GetX()[gYieldStat[iPart]->GetN() - 1];

            linearYieldFits[iPart] = new TF1(fitName.Data(), "pol1", dNchMin - (dNchMax - dNchMin) * 0.1, dNchMax + (dNchMax - dNchMin) * 0.1); // Draw range slightly extended
            linearYieldFits[iPart]->SetLineColor(customCascPalette[iPart]);                                                                     // Use a different color
            linearYieldFits[iPart]->SetLineStyle(kDashed);
            // linearYieldFits[iPart]->SetLineWidth(2);
            linearYieldFits[iPart]->SetNpx(500); // Increase number of points for smoother line

            gYieldStat[iPart]->Fit(linearYieldFits[iPart], "QERNS+rob=0.6", "", dNchMin, dNchMax); // Fit within data range
            linearYieldFits[iPart]->Draw("SAME");
        }

        if (linearYieldFits[iPart])
        {
            // Add fit parameters and chi2/ndf to the plot
            double chi2 = linearYieldFits[iPart]->GetChisquare();
            int ndf = linearYieldFits[iPart]->GetNDF();
            double chi2ndf = (ndf > 0) ? chi2 / ndf : 0.0;              // Avoid division by zero
                                                                        // Create Legend
            TLegend *legYieldFit = new TLegend(0.15, 0.7, 0.27, 0.775); // Adjust as needed
            legYieldFit->SetBorderSize(0);
            legYieldFit->SetFillStyle(0); // Transparent
            legYieldFit->SetTextSize(0.022);
            // legYieldFit->SetHeader(Form("%s, p-Pb #sqrt{s_{NN}} = 8.16 TeV", particleTitlePlot[iPart]), "C");
            // // legYieldFit->AddEntry((TObject*)nullptr, Form("p-Pb #sqrt{s_{NN}} = 8.16 TeV, %s", particleTitlePlot[iPart]), "");
            // legYieldFit->AddEntry((TObject *)nullptr, "-0.5 < y < 0", "");
            // if (gYieldStat[iPart])
            // legYieldFit->AddEntry(gYieldStat[iPart], "Data (stat.)", "pe");
            // if (gYieldSys[iPart])
            // legYieldFit->AddEntry(gYieldSys[iPart], "Data (sys.)", "f");
            // if (linearYieldFits[iPart])
            legYieldFit->AddEntry(linearYieldFits[iPart], Form("Linear Fit  #(){#frac{#chi^{2}}{NDF} = %.2f}", chi2ndf), "l");
            legYieldFit->Draw();
        }

        TLatex latex;
        latex.SetNDC();
        latex.SetTextSize(0.025);
        // TLatex *latexInfo = new TLatex(0.65, 0.81, "Uncertainties: stat.(bars), sys.(boxes)");
        latex.DrawLatex(0.15, 0.86, "p-Pb #sqrt{s_{NN}} = 8.16 TeV, This work");
        latex.DrawLatex(0.15, 0.82, Form("%s, -0.5 < y < 0", particleTitlePlot[iPart]));
        latex.SetTextSize(0.022);
        latex.SetTextFont(42);
        latex.DrawLatex(0.15, 0.78, "Uncertainties: stat.(bars), sys.(boxes)");

        cYieldWithFit->Update();
        SaveImage(outputPlotFolder, "YieldVsMult_WithFit", TString::Format("YieldVsMult_Fit_%s", particleNamePlot[iPart]), imageFormat.Data(), cYieldWithFit);
    }
    delete cYieldWithFit;

    // --- NEW: Plot ratio of linear yield fits (Omega/Xi) ---
    Info("IntYieldsMeanPtCompare", "Generating ratio of linear yield fits (Omega/Xi)...");
    if (nParticlesPlot >= 2 && linearYieldFits[0] && linearYieldFits[1])
    {
        TCanvas *cRatioYieldFits = new TCanvas("cRatioYieldFits", "Ratio of Yield Fits (Omega/Xi) vs dNch/deta", 1200, 800);
        // cRatioYieldFits->SetLeftMargin(0.15);
        cRatioYieldFits->SetBottomMargin(0.12);
        cRatioYieldFits->SetRightMargin(0.05);
        cRatioYieldFits->SetTopMargin(0.05);
        // cRatioYieldFits->SetGrid();
        // cRatioYieldFits->SetGridy();

        // Determine a common dNch range for the ratio plot based on where both original graphs have points.
        double dNch_min_common = 4.47;
        double dNch_max_common = 53.22; // Fallback

        if (gYieldStat[0] && gYieldStat[0]->GetN() > 0 && gYieldStat[1] && gYieldStat[1]->GetN() > 0)
        {
            double xi_min_dNch = gYieldStat[0]->GetX()[gYieldStat[0]->GetN() - 1];
            double xi_max_dNch = gYieldStat[0]->GetX()[0];
            double om_min_dNch = gYieldStat[1]->GetX()[gYieldStat[1]->GetN() - 1];
            double om_max_dNch = gYieldStat[1]->GetX()[0];
            dNch_min_common = TMath::Max(xi_min_dNch, om_min_dNch);
            dNch_max_common = TMath::Min(xi_max_dNch, om_max_dNch);
        }
        // Ensure min < max, provide a sensible default if overlap is problematic
        if (dNch_max_common <= dNch_min_common)
        {
            dNch_min_common = 0;
            dNch_max_common = 55;
            Warning("IntYieldsMeanPtCompare", "Common dNch range for fit ratio is ill-defined, using default [%f, %f]", dNch_min_common, dNch_max_common);
        }

        TF1 *ratioFunc = new TF1("ratio_Omega_Xi_YieldFits_Func", [&](double *x, double *p)
                                 {
            // linearYieldFits[0] is for particleNamePlot[0] (expected Xi)
            // linearYieldFits[1] is for particleNamePlot[1] (expected Omega)
            double val_xi = linearYieldFits[0]->Eval(x[0]);
            double val_om = linearYieldFits[1]->Eval(x[0]);
            if (TMath::Abs(val_xi) > 1e-9) { // Avoid division by zero or very small numbers
                return val_om / val_xi;
            }
            // Return NaN or a very large/small number if denominator is zero to indicate problem
            return TMath::QuietNaN(); }, dNch_min_common, dNch_max_common, 0);

        ratioFunc->SetLineColor(customCascPalette[4]); // Use a distinct color for the ratio
        ratioFunc->SetLineWidth(2);
        ratioFunc->SetNpx(500); // Smooth curve

        // Estimate Y-axis range for ratio
        double y_min_ratio = 0;   // Sensible default
        double y_max_ratio = 0.2; // Sensible default
        bool range_ok = false;
        if (dNch_max_common > dNch_min_common)
        {
            y_min_ratio = ratioFunc->GetMinimum(dNch_min_common, dNch_max_common, 1.e-5, 1000);
            y_max_ratio = ratioFunc->GetMaximum(dNch_min_common, dNch_max_common, 1.e-5, 1000);
            if (!TMath::IsNaN(y_min_ratio) && !TMath::IsNaN(y_max_ratio) && y_max_ratio > y_min_ratio)
            {
                range_ok = true;
            }
        }

        if (!range_ok)
        { // Fallback if GetMinimum/Maximum fails or returns NaN
            y_min_ratio = 0.0;
            y_max_ratio = 0.2; // Default if range finding fails
            Warning("IntYieldsMeanPtCompare", "Could not determine Y-range for ratio plot dynamically, using default [%f, %f]", y_min_ratio, y_max_ratio);
        }

        double y_padding_ratio = (y_max_ratio - y_min_ratio) * 0.15;
        if (y_padding_ratio == 0 && y_max_ratio != 0)
            y_padding_ratio = TMath::Abs(y_max_ratio * 0.2);
        else if (y_padding_ratio == 0)
            y_padding_ratio = 0.02;

        TH1 *hFrameRatioFits = gPad->DrawFrame(dNch_min_common, y_min_ratio - y_padding_ratio,
                                               dNch_max_common, y_max_ratio + y_padding_ratio);

        // hFrameRatioFits->SetTitle("Enhancement: Fit(#Omega) / Fit(#Xi) Yields");
        hFrameRatioFits->GetXaxis()->SetTitle("#LTd#it{N}_{ch}/d#it{#eta}#GT_{|#it{#eta}|<0.5}");
        hFrameRatioFits->GetYaxis()->SetTitle("#Omega^{#pm}/#Xi^{#pm} Fit ratio");
        hFrameRatioFits->GetYaxis()->SetTitleOffset(1.4);
        hFrameRatioFits->GetXaxis()->SetTitleOffset(1.3);

        ratioFunc->Draw("SAME L"); // "L" to ensure line is drawn, default might be points for TF1

        TLatex latex;
        latex.SetNDC();
        latex.SetTextSize(0.025);
        latex.DrawLatex(0.15, 0.86, "p-Pb #sqrt{s_{NN}} = 8.16 TeV, This work");
        // latex.DrawLatex(0.15, 0.82, Form("-0.5 < y < 0"));
        TLegend *legRatioYieldFits = new TLegend(0.15, 0.77, 0.3, 0.815); // Adjust as needed
        legRatioYieldFits->SetBorderSize(0);
        legRatioYieldFits->SetFillStyle(0); // Transparent
        legRatioYieldFits->SetTextSize(0.028);
        legRatioYieldFits->AddEntry(ratioFunc, Form("#frac{#LTd#it{N}/d#it{y}#GT_{%s}^{Fit}}{#LTd#it{N}/d#it{y}#GT_{%s}^{Fit}}", particleTitlePlot[1], particleTitlePlot[0]), "l");
        // latex.DrawLatex(0.15, 0.78, Form("#frac{#LTd#it{N}/d#it{y}#GT_{%s}^{Fit}}{#LTd#it{N}/d#it{y}#GT_{%s}^{Fit}}", particleTitlePlot[1], particleTitlePlot[0]));
        legRatioYieldFits->Draw();
        cRatioYieldFits->Update();
        SaveImage(outputPlotFolder, "YieldFitRatios", "Ratio_YieldFits_Omega_Xi", imageFormat.Data(), cRatioYieldFits);

        delete ratioFunc;
        // delete hFrameRatioFits; // Owned by pad
        delete cRatioYieldFits;
    }
    else
    {
        Warning("IntYieldsMeanPtCompare", "Linear yield fits for Omega and/or Xi not available. Skipping ratio plot.");
        if (nParticlesPlot < 2)
        {
            Warning("IntYieldsMeanPtCompare", "nParticlesPlot is %d, need at least 2 for Omega/Xi ratio.", nParticlesPlot);
        }
        if (!linearYieldFits[0])
            Warning("IntYieldsMeanPtCompare", "Fit for %s (index 0) is missing.", particleNamePlot[0]);
        if (nParticlesPlot >= 2 && !linearYieldFits[1])
            Warning("IntYieldsMeanPtCompare", "Fit for %s (index 1) is missing.", particleNamePlot[1]);
    }
    // --- END OF RATIO PLOT SECTION ---

    // --- FINAL SECTION: Generate Yield vs Multiplicity plots with all data sources ---
    Info("IntYieldsMeanPtCompare", "Generating Yield vs. Multiplicity plots with extensive data comparison from pp, p-Pb, and Pb-Pb...");
    // --- Create transformed dNch values for the 2014 Pb-Pb data by averaging/mapping ---
    double dNch_2014_transformed[Npoints_2014];
    double dNch_err_2014_transformed[Npoints_2014];
    // Bin 0-10%  (average of 0-5 and 5-10)
    dNch_2014_transformed[0] = (dNch_map_2011[0] + dNch_map_2011[1]) / 2.0;
    dNch_err_2014_transformed[0] = TMath::Sqrt(TMath::Power(dNch_err_map_2011[0], 2) + TMath::Power(dNch_err_map_2011[1], 2)) / 2.0;
    // Bin 10-20% (direct map)
    dNch_2014_transformed[1] = dNch_map_2011[2];
    dNch_err_2014_transformed[1] = dNch_err_map_2011[2];
    // Bin 20-40% (average of 20-30 and 30-40)
    dNch_2014_transformed[2] = (dNch_map_2011[3] + dNch_map_2011[4]) / 2.0;
    dNch_err_2014_transformed[2] = TMath::Sqrt(TMath::Power(dNch_err_map_2011[3], 2) + TMath::Power(dNch_err_map_2011[4], 2)) / 2.0;
    // Bin 40-60% (average of 40-50 and 50-60)
    dNch_2014_transformed[3] = (dNch_map_2011[5] + dNch_map_2011[6]) / 2.0;
    dNch_err_2014_transformed[3] = TMath::Sqrt(TMath::Power(dNch_err_map_2011[5], 2) + TMath::Power(dNch_err_map_2011[6], 2)) / 2.0;
    // Bin 60-80% (average of 60-70 and 70-80)
    dNch_2014_transformed[4] = (dNch_map_2011[7] + dNch_map_2011[8]) / 2.0;
    dNch_err_2014_transformed[4] = TMath::Sqrt(TMath::Power(dNch_err_map_2011[7], 2) + TMath::Power(dNch_err_map_2011[8], 2)) / 2.0;

    for (int iPart = 0; iPart < nParticlesPlot; ++iPart)
    {
        TCanvas *cYield_Comp = new TCanvas(Form("cYieldVsMult_FinalComp_%s", particleNamePlot[iPart]),
                                           Form("Yield vs dNch/deta (%s) - All Systems", particleTitlePlot[iPart]),
                                           1200, 900);
        cYield_Comp->cd();
        // gPad->SetLeftMargin(0.15);
        gPad->SetBottomMargin(0.12);
        gPad->SetTopMargin(0.05);
        gPad->SetRightMargin(0.05);
        gPad->SetLogy(); // Use log scale for y-axis to see all yields
        gPad->SetLogx(); // Use log scale for x-axis to see wide range of dNch

        // Step 1: Aggregate data from ALL sources for frame calculation
        std::vector<double> all_dNch_vals;
        std::vector<double> all_Yield_vals_min;
        std::vector<double> all_Yield_vals_max;

        // Source 1: This work's data
        if (gYieldStat[iPart] && gYieldSys[iPart])
        {
            for (int i = 0; i < gYieldStat[iPart]->GetN(); ++i)
            {
                if (gYieldStat[iPart]->GetY()[i] <= 0)
                    continue; // Don't include non-positive values for log scale range
                all_dNch_vals.push_back(gYieldStat[iPart]->GetX()[i]);
                all_Yield_vals_min.push_back(gYieldStat[iPart]->GetY()[i] - gYieldSys[iPart]->GetErrorYlow(i));
                all_Yield_vals_max.push_back(gYieldStat[iPart]->GetY()[i] + gYieldSys[iPart]->GetErrorYhigh(i));
            }
        }

        // Aggregate other data sources based on particle type
        if (iPart == 0)
        { // Xi data sets
            // Source 2: Table 3 (p-Pb 5.02 TeV)
            for (int i = 0; i < nPoints_ppb5_Xi; ++i)
            {
                all_dNch_vals.push_back(dNch_ppb5_Xi[i]);
                all_Yield_vals_min.push_back(yield_ppb5_Xi[i] - sysErr_ppb5_Xi[i]);
                all_Yield_vals_max.push_back(yield_ppb5_Xi[i] + sysErr_ppb5_Xi[i]);
            }
            // Source 3: 2020 (pp 13 TeV)
            for (int i = 0; i < Npoints_Xi_2020; ++i)
            {
                all_dNch_vals.push_back(dNch_Xi_2020[i]);
                all_Yield_vals_min.push_back(XiYields_2020[i] - XiYields_systTot_2020[i]);
                all_Yield_vals_max.push_back(XiYields_2020[i] + XiYields_systTot_2020[i]);
            }
            // Source 4: 2017 (pp 7 TeV)
            for (int i = 0; i < Npoints_Xi2017; ++i)
            {
                all_dNch_vals.push_back(dNch_Xi_2017[i]);
                all_Yield_vals_min.push_back(XiYields2017[i] - XiYields2017_systTot[i]);
                all_Yield_vals_max.push_back(XiYields2017[i] + XiYields2017_systTot[i]);
            }
            // Source 5: 2014 (Pb-Pb 2.76 TeV)
            for (int i = 0; i < Npoints_2014; ++i)
            {
                all_dNch_vals.push_back(dNch_2014_transformed[i]);
                all_Yield_vals_min.push_back(XiCombinedYields2014[i] - XiCombinedYields2014_systTot[i]);
                all_Yield_vals_max.push_back(XiCombinedYields2014[i] + XiCombinedYields2014_systTot[i]);
            }
        }
        else
        { // Omega data sets
            // Source 2: Table 3 (p-Pb 5.02 TeV)
            for (int i = 0; i < nPoints_ppb5_Om; ++i)
            {
                all_dNch_vals.push_back(dNch_ppb5_Om[i]);
                all_Yield_vals_min.push_back(yield_ppb5_Om[i] - sysErr_ppb5_Om[i]);
                all_Yield_vals_max.push_back(yield_ppb5_Om[i] + sysErr_ppb5_Om[i]);
            }
            // Source 3: 2020 (pp 13 TeV)
            for (int i = 0; i < Npoints_Omega_2020; ++i)
            {
                all_dNch_vals.push_back(dNch_Omega_2020[i]);
                all_Yield_vals_min.push_back(OmegaYields_2020[i] - OmegaYields_2020_systTot[i]);
                all_Yield_vals_max.push_back(OmegaYields_2020[i] + OmegaYields_2020_systTot[i]);
            }
            // Source 4: 2017 (pp 7 TeV)
            for (int i = 0; i < Npoints_Omega2017; ++i)
            {
                all_dNch_vals.push_back(dNch_Om_2017[i]);
                all_Yield_vals_min.push_back(OmegaYields2017[i] - OmegaYields2017_systTot[i]);
                all_Yield_vals_max.push_back(OmegaYields2017[i] + OmegaYields2017_systTot[i]);
            }
            // Source 5: 2014 (Pb-Pb 2.76 TeV)
            for (int i = 0; i < Npoints_2014; ++i)
            {
                all_dNch_vals.push_back(dNch_2014_transformed[i]);
                all_Yield_vals_min.push_back(OmegaCombinedYields2014[i] - OmegaCombinedYields2014_systTot[i]);
                all_Yield_vals_max.push_back(OmegaCombinedYields2014[i] + OmegaCombinedYields2014_systTot[i]);
            }
        }

        // Calculate final frame ranges
        double frame_min_dNch = (iPart == 0 ? 1.5 : 1.8);
        double frame_max_dNch = 2000;
        double frame_min_Yield = (iPart == 0 ? 0.003 : 0.0003);
        double frame_max_Yield = (iPart == 0 ? 10 : 2);

        // Step 2: Draw Frame and All Data
        TH1 *hFrameYield_Comp = gPad->DrawFrame(frame_min_dNch, frame_min_Yield, frame_max_dNch, frame_max_Yield, ";#LTd#it{N}_{ch}/d#it{#eta}#GT_{|#it{#eta}|<0.5};#LTd#it{N}/d#it{y}#GT");
        hFrameYield_Comp->GetXaxis()->SetTitleOffset(1.3);

        // --- Draw Data Source 1: This Work (p-Pb, 8.16 TeV) ---
        if (gYieldSys[iPart])
            gYieldSys[iPart]->Draw("E2 SAME");
        if (gYieldStat[iPart])
            gYieldStat[iPart]->Draw("PE SAME");

        // --- Draw other data sources ---
        auto draw_data = [&](int n, double x[], double y[], double xerr[], double ystat[], double ysys[], int marker, int color)
        {
            auto gstat = new TGraphErrors(n, x, y, xerr, ystat);
            auto gsys = new TGraphAsymmErrors(n, x, y, xerr, xerr, ysys, ysys);
            gstat->SetMarkerStyle(marker);
            gstat->SetMarkerColor(color);
            gstat->SetLineColor(color);
            gsys->SetFillStyle(0);
            gsys->SetLineColor(color);
            gsys->Draw("E2 SAME");
            gstat->Draw("PE SAME");
            return std::make_pair(gstat, gsys);
        };

        std::pair<TGraphErrors *, TGraphAsymmErrors *> graphs_table3, graphs_2020, graphs_2017, graphs_2014;

        if (iPart == 0)
        { // Xi
            graphs_table3 = draw_data(nPoints_ppb5_Xi, dNch_ppb5_Xi, yield_ppb5_Xi, dNchErr_ppb5_Xi, statErr_ppb5_Xi, sysErr_ppb5_Xi, markerStyles[2], customCascPalette[2]);
            graphs_2020 = draw_data(Npoints_Xi_2020, dNch_Xi_2020, XiYields_2020, dNch_err_Xi_2020, XiYields_stat_2020, XiYields_systTot_2020, markerStyles[3], customCascPalette[3]);
            graphs_2017 = draw_data(Npoints_Xi2017, dNch_Xi_2017, XiYields2017, dNch_err_Xi_2017, XiYields2017_stat, XiYields2017_systTot, markerStyles[4], customCascPalette[4]);
            graphs_2014 = draw_data(Npoints_2014, dNch_2014_transformed, XiCombinedYields2014, dNch_err_2014_transformed, XiCombinedYields2014_stat, XiCombinedYields2014_systTot, markerStyles[6], customCascPalette[6]);
        }
        else
        { // Omega
            graphs_table3 = draw_data(nPoints_ppb5_Om, dNch_ppb5_Om, yield_ppb5_Om, dNchErr_ppb5_Om, statErr_ppb5_Om, sysErr_ppb5_Om, markerStyles[2], customCascPalette[2]);
            graphs_2020 = draw_data(Npoints_Omega_2020, dNch_Omega_2020, OmegaYields_2020, dNch_err_Omega_2020, OmegaYields_2020_stat, OmegaYields_2020_systTot, markerStyles[3], customCascPalette[3]);
            graphs_2017 = draw_data(Npoints_Omega2017, dNch_Om_2017, OmegaYields2017, dNch_err_Om_2017, OmegaYields2017_stat, OmegaYields2017_systTot, markerStyles[4], customCascPalette[4]);
            graphs_2014 = draw_data(Npoints_2014, dNch_2014_transformed, OmegaCombinedYields2014, dNch_err_2014_transformed, OmegaCombinedYields2014_stat, OmegaCombinedYields2014_systTot, markerStyles[6], customCascPalette[6]);
        }

        // Step 3: Create a comprehensive legend
        TLegend *legYield_Comp = new TLegend(0.12, 0.68, 0.4, 0.87); // Adjusted for more entries
        legYield_Comp->SetBorderSize(0);
        legYield_Comp->SetFillStyle(0); // Transparent background
        legYield_Comp->SetTextSize(0.025);

        if (gYieldStat[iPart])
            legYield_Comp->AddEntry(gYieldStat[iPart], "p-Pb #sqrt{s_{NN}} = 8.16 TeV, This work", "pe");
        if (graphs_table3.first)
            legYield_Comp->AddEntry(graphs_table3.first, "p-Pb #sqrt{s_{NN}} = 5.02 TeV, Phys.Lett.B 758 (2016), 389-401", "pe");
        if (graphs_2020.first)
            legYield_Comp->AddEntry(graphs_2020.first, "pp #sqrt{s} = 13 TeV, Eur.Phys.J.C 80 (2020), 167", "pe");
        if (graphs_2017.first)
            legYield_Comp->AddEntry(graphs_2017.first, "pp #sqrt{s} = 7 TeV, Nature Phys. 13 (2017), 535-539", "pe");
        if (graphs_2014.first)
            legYield_Comp->AddEntry(graphs_2014.first, "Pb-Pb #sqrt{s_{NN}} = 2.76 TeV, Phys.Lett.B 728 (2014) 216-227", "pe");
        legYield_Comp->Draw();

        TLatex latexComp;
        latexComp.SetNDC();
        latexComp.SetTextSize(0.028);
        latexComp.DrawLatex(0.19, 0.88, (Form("%s", particleTitlePlot[iPart])));
        // latexComp.DrawLatex(0.18, 0.84, "-0.5 < y < 0");
        latexComp.SetTextSize(0.023);
        latexComp.SetTextFont(42);
        latexComp.DrawLatex(0.19, 0.65, "Uncertainties: stat.(bars), sys.(boxes)");

        cYield_Comp->Update();
        SaveImage(outputPlotFolder, "YieldVsMult_Comparison", Form("YieldVsMult_FinalComp_%s", particleNamePlot[iPart]), imageFormat.Data(), cYield_Comp);

        // Cleanup
        delete graphs_table3.first;
        delete graphs_table3.second;
        delete graphs_2020.first;
        delete graphs_2020.second;
        delete graphs_2017.first;
        delete graphs_2017.second;
        delete graphs_2014.first;
        delete graphs_2014.second;
        delete cYield_Comp;
    }

    for (int iPart = 0; iPart < nParticlesPlot; ++iPart)
    {
        TCanvas *cYield_Comp = new TCanvas(Form("cYieldVsMult_CompExtensive_%s", particleNamePlot[iPart]),
                                           Form("Yield vs dNch/deta (%s) - Extensive Comparison", particleTitlePlot[iPart]),
                                           1200, 900); // Taller canvas for bigger legend
        cYield_Comp->cd();
        gPad->SetLeftMargin(0.12);
        gPad->SetBottomMargin(0.12);
        gPad->SetTopMargin(0.05);
        gPad->SetRightMargin(0.03);

        // --- Step 1: Aggregate data from ALL sources to calculate plot frame ---
        std::vector<double> all_dNch_vals;
        std::vector<double> all_Yield_vals_min; // yield - total_sys_err
        std::vector<double> all_Yield_vals_max; // yield + total_sys_err

        // Source 1: This work's data (from results file)
        if (gYieldStat[iPart] && gYieldSys[iPart])
        {
            for (int i = 0; i < gYieldStat[iPart]->GetN(); ++i)
            {
                all_dNch_vals.push_back(gYieldStat[iPart]->GetX()[i]);
                all_Yield_vals_min.push_back(gYieldStat[iPart]->GetY()[i] - gYieldSys[iPart]->GetErrorYlow(i));
                all_Yield_vals_max.push_back(gYieldStat[iPart]->GetY()[i] + gYieldSys[iPart]->GetErrorYhigh(i));
            }
        }

        // Add data points for each particle type
        if (iPart == 0)
        { // Xi data sets
            // Source 2: Table 3 Data (p-Pb @ 5.02 TeV)
            for (int i = 0; i < nPoints_ppb5_Xi; ++i)
            {
                all_dNch_vals.push_back(dNch_ppb5_Xi[i]);
                all_Yield_vals_min.push_back(yield_ppb5_Xi[i] - sysErr_ppb5_Xi[i]);
                all_Yield_vals_max.push_back(yield_ppb5_Xi[i] + sysErr_ppb5_Xi[i]);
            }
            // Source 3: 2020 Data (pp @ 13 TeV)
            for (int i = 0; i < Npoints_Xi_2020; ++i)
            {
                all_dNch_vals.push_back(dNch_Xi_2020[i]);
                all_Yield_vals_min.push_back(XiYields_2020[i] - XiYields_systTot_2020[i]);
                all_Yield_vals_max.push_back(XiYields_2020[i] + XiYields_systTot_2020[i]);
            }
            // Source 4: 2017 Data (pp @ 7 TeV)
            // Assumes dNch2017 is renamed to dNch_Xi_2017 in your header
            for (int i = 0; i < Npoints_Xi2017; ++i)
            {
                all_dNch_vals.push_back(dNch_Xi_2017[i]);
                all_Yield_vals_min.push_back(XiYields2017[i] - XiYields2017_systTot[i]);
                all_Yield_vals_max.push_back(XiYields2017[i] + XiYields2017_systTot[i]);
            }
        }
        else
        { // Omega data sets
            // Source 2: Table 3 Data (p-Pb @ 5.02 TeV)
            for (int i = 0; i < nPoints_ppb5_Om; ++i)
            {
                all_dNch_vals.push_back(dNch_ppb5_Om[i]);
                all_Yield_vals_min.push_back(yield_ppb5_Om[i] - sysErr_ppb5_Om[i]);
                all_Yield_vals_max.push_back(yield_ppb5_Om[i] + sysErr_ppb5_Om[i]);
            }
            // Source 3: 2020 Data (pp @ 13 TeV)
            for (int i = 0; i < Npoints_Omega_2020; ++i)
            {
                all_dNch_vals.push_back(dNch_Omega_2020[i]);
                all_Yield_vals_min.push_back(OmegaYields_2020[i] - OmegaYields_2020_systTot[i]);
                all_Yield_vals_max.push_back(OmegaYields_2020[i] + OmegaYields_2020_systTot[i]);
            }
            // Source 4: 2017 Data (pp @ 7 TeV)
            // Assumes dNch2017 is renamed to dNch_Om_2017 in your header
            for (int i = 0; i < Npoints_Omega2017; ++i)
            {
                all_dNch_vals.push_back(dNch_Om_2017[i]);
                all_Yield_vals_min.push_back(OmegaYields2017[i] - OmegaYields2017_systTot[i]);
                all_Yield_vals_max.push_back(OmegaYields2017[i] + OmegaYields2017_systTot[i]);
            }
        }

        // Calculate final frame ranges with padding
        double frame_min_dNch = 0, frame_max_dNch = 60;
        double frame_min_Yield = 0, frame_max_Yield = 0.1;
        if (!all_dNch_vals.empty())
        {
            frame_min_dNch = TMath::MinElement(all_dNch_vals.size(), all_dNch_vals.data());
            frame_max_dNch = TMath::MaxElement(all_dNch_vals.size(), all_dNch_vals.data());
            double dNchPadding = (frame_max_dNch - frame_min_dNch) * 0.1;
            if (dNchPadding == 0)
                dNchPadding = 5;
            frame_min_dNch = TMath::Max(0.0, frame_min_dNch - dNchPadding);
            frame_max_dNch = frame_max_dNch + dNchPadding;
        }
        if (!all_Yield_vals_min.empty())
        {
            frame_min_Yield = TMath::MinElement(all_Yield_vals_min.size(), all_Yield_vals_min.data());
            frame_max_Yield = TMath::MaxElement(all_Yield_vals_max.size(), all_Yield_vals_max.data());
            double yieldPadding = (frame_max_Yield - frame_min_Yield) * 0.1;
            frame_min_Yield -= yieldPadding;
            frame_max_Yield += yieldPadding;
        }

        // --- Step 2: Draw Frame and All Data ---
        TH1 *hFrameYield_Comp = gPad->DrawFrame(frame_min_dNch, frame_min_Yield, frame_max_dNch, frame_max_Yield, ";#LTd#it{N}_{ch}/d#it{#eta}#GT_{|#it{#eta}|<0.5};#LTd#it{N}/d#it{y}#GT");
        // hFrameYield_Comp->GetYaxis()->SetTitleOffset(1.4);
        hFrameYield_Comp->GetXaxis()->SetTitleOffset(1.3);

        // --- Draw Data Source 1: This Work ---
        if (gYieldSys[iPart])
            gYieldSys[iPart]->Draw("E2 SAME");
        if (gYieldStat[iPart])
            gYieldStat[iPart]->Draw("PE SAME");

        // --- Draw Data Source 2: Table 3 (p-Pb @ 5.02 TeV) ---
        TGraphErrors *gYieldStat_Table3 = nullptr;
        TGraphAsymmErrors *gYieldSys_Table3 = nullptr;
        if (iPart == 0)
        { // Xi
            gYieldStat_Table3 = new TGraphErrors(nPoints_ppb5_Xi, dNch_ppb5_Xi, yield_ppb5_Xi, dNchErr_ppb5_Xi, statErr_ppb5_Xi);
            gYieldSys_Table3 = new TGraphAsymmErrors(nPoints_ppb5_Xi, dNch_ppb5_Xi, yield_ppb5_Xi, dNchErr_ppb5_Xi, dNchErr_ppb5_Xi, sysErr_ppb5_Xi, sysErr_ppb5_Xi);
        }
        else
        { // Omega
            gYieldStat_Table3 = new TGraphErrors(nPoints_ppb5_Om, dNch_ppb5_Om, yield_ppb5_Om, dNchErr_ppb5_Om, statErr_ppb5_Om);
            gYieldSys_Table3 = new TGraphAsymmErrors(nPoints_ppb5_Om, dNch_ppb5_Om, yield_ppb5_Om, dNchErr_ppb5_Om, dNchErr_ppb5_Om, sysErr_ppb5_Om, sysErr_ppb5_Om);
        }
        gYieldStat_Table3->SetMarkerStyle(markerStyles[2]);
        gYieldStat_Table3->SetMarkerColor(customCascPalette[2]);
        gYieldStat_Table3->SetLineColor(customCascPalette[2]);
        gYieldSys_Table3->SetLineColor(customCascPalette[2]);
        gYieldSys_Table3->SetFillStyle(0);
        gYieldSys_Table3->Draw("E2 SAME");
        gYieldStat_Table3->Draw("PE SAME");

        // --- Draw Data Source 3: 2020 Data (pp @ 13 TeV) ---
        TGraphErrors *gYieldStat_2020 = nullptr;
        TGraphAsymmErrors *gYieldSys_2020 = nullptr;
        if (iPart == 0)
        { // Xi
            gYieldStat_2020 = new TGraphErrors(Npoints_Xi_2020, dNch_Xi_2020, XiYields_2020, dNch_err_Xi_2020, XiYields_stat_2020);
            gYieldSys_2020 = new TGraphAsymmErrors(Npoints_Xi_2020, dNch_Xi_2020, XiYields_2020, dNch_err_Xi_2020, dNch_err_Xi_2020, XiYields_systTot_2020, XiYields_systTot_2020);
        }
        else
        { // Omega
            gYieldStat_2020 = new TGraphErrors(Npoints_Omega_2020, dNch_Omega_2020, OmegaYields_2020, dNch_err_Omega_2020, OmegaYields_2020_stat);
            gYieldSys_2020 = new TGraphAsymmErrors(Npoints_Omega_2020, dNch_Omega_2020, OmegaYields_2020, dNch_err_Omega_2020, dNch_err_Omega_2020, OmegaYields_2020_systTot, OmegaYields_2020_systTot);
        }
        gYieldStat_2020->SetMarkerStyle(markerStyles[3]);
        gYieldStat_2020->SetMarkerColor(customCascPalette[3]);
        gYieldStat_2020->SetLineColor(customCascPalette[3]);
        gYieldSys_2020->SetLineColor(customCascPalette[3]);
        gYieldSys_2020->SetFillStyle(0);
        gYieldSys_2020->Draw("E2 SAME");
        gYieldStat_2020->Draw("PE SAME");

        // --- Draw Data Source 4: 2017 Data (pp @ 7 TeV) ---
        TGraphErrors *gYieldStat_2017 = nullptr;
        TGraphAsymmErrors *gYieldSys_2017 = nullptr;
        if (iPart == 0)
        { // Xi
            gYieldStat_2017 = new TGraphErrors(Npoints_Xi2017, dNch_Xi_2017, XiYields2017, dNch_err_Xi_2017, XiYields2017_stat);
            gYieldSys_2017 = new TGraphAsymmErrors(Npoints_Xi2017, dNch_Xi_2017, XiYields2017, dNch_err_Xi_2017, dNch_err_Xi_2017, XiYields2017_systTot, XiYields2017_systTot);
        }
        else
        { // Omega
            gYieldStat_2017 = new TGraphErrors(Npoints_Omega2017, dNch_Om_2017, OmegaYields2017, dNch_err_Om_2017, OmegaYields2017_stat);
            gYieldSys_2017 = new TGraphAsymmErrors(Npoints_Omega2017, dNch_Om_2017, OmegaYields2017, dNch_err_Om_2017, dNch_err_Om_2017, OmegaYields2017_systTot, OmegaYields2017_systTot);
        }
        gYieldStat_2017->SetMarkerStyle(markerStyles[4]);
        gYieldStat_2017->SetMarkerColor(customCascPalette[4]);
        gYieldStat_2017->SetLineColor(customCascPalette[4]);
        gYieldSys_2017->SetLineColor(customCascPalette[4]);
        gYieldSys_2017->SetFillStyle(0);
        gYieldSys_2017->Draw("E2 SAME");
        gYieldStat_2017->Draw("PE SAME");

        // --- Step 3: Create a comprehensive legend for all sources ---
        TLegend *legYield_Comp = new TLegend(0.12, 0.7, 0.4, 0.87); // Adjusted for more entries
        legYield_Comp->SetBorderSize(0);
        legYield_Comp->SetFillStyle(0); // Transparent background
        legYield_Comp->SetTextSize(0.025);

        if (gYieldStat[iPart])
            legYield_Comp->AddEntry(gYieldStat[iPart], "p-Pb #sqrt{s_{NN}} = 8.16 TeV, This work", "pe");
        if (gYieldStat_Table3)
            legYield_Comp->AddEntry(gYieldStat_Table3, "p-Pb #sqrt{s_{NN}} = 5.02 TeV, Phys.Lett.B 758 (2016), 389-401", "pe");
        if (gYieldStat_2020)
            legYield_Comp->AddEntry(gYieldStat_2020, "pp #sqrt{s} = 13 TeV, Eur.Phys.J.C 80 (2020), 167", "pe");
        if (gYieldStat_2017)
            legYield_Comp->AddEntry(gYieldStat_2017, "pp #sqrt{s} = 7 TeV, Nature Phys. 13 (2017), 535-539", "pe");
        legYield_Comp->Draw();

        TLatex latexComp;
        latexComp.SetNDC();
        latexComp.SetTextSize(0.028);
        latexComp.DrawLatex(0.19, 0.88, (Form("%s", particleTitlePlot[iPart])));
        //
        // latexComp.DrawLatex(0.18, 0.88, Form("p-Pb #sqrt{s_{NN}} = 8.16 TeV, %s", particleTitlePlot[iPart]));
        // latexComp.DrawLatex(0.18, 0.84, "-0.5 < y < 0");
        latexComp.SetTextSize(0.023);
        latexComp.SetTextFont(42);
        latexComp.DrawLatex(0.19, 0.67, "Uncertainties: stat.(bars), sys.(boxes)");
        // latexComp.DrawLatex(0.18, 0.80, "Data: This work (circles), IntYields PDF (squares)");

        cYield_Comp->Update();
        SaveImage(outputPlotFolder, "YieldVsMult_Comparison",
                  Form("YieldVsMult_CompExtensive_%s", particleNamePlot[iPart]),
                  imageFormat.Data(), cYield_Comp);

        // Cleanup for this loop iteration
        delete gYieldStat_Table3;
        delete gYieldSys_Table3;
        delete gYieldStat_2020;
        delete gYieldSys_2020;
        delete gYieldStat_2017;
        delete gYieldSys_2017;
        delete legYield_Comp;
        delete cYield_Comp;
    }
    // --- END OF EXTENSIVE COMPARISON PLOT SECTION ---

    // --- NEW SECTION: Generate Yield vs Multiplicity plots with Table 3 Data Overlay ---
    Info("IntYieldsMeanPtCompare", "Generating Yield vs. Multiplicity plots with Table 3 Data Comparison...");

    for (int iPart = 0; iPart < nParticlesPlot; ++iPart)
    {
        TCanvas *cYield_Comp = new TCanvas(Form("cYieldVsMult_Compare_%s", particleNamePlot[iPart]),
                                           Form("Yield vs dNch/deta (%s) - Comparison with IntYields from p-Pb at 5.02 TeV", particleTitlePlot[iPart]),
                                           1200, 800); // Match size of original cYield or adjust as needed
        cYield_Comp->cd();
        // gPad->SetLeftMargin(0.15);
        gPad->SetBottomMargin(0.12);
        gPad->SetTopMargin(0.05); // Increased top margin for a comprehensive title
        gPad->SetRightMargin(0.05);

        // Determine combined range for the plot frame
        std::vector<double> all_dNch_vals;
        std::vector<double> all_Yield_vals_min; // yield - total_sys_err
        std::vector<double> all_Yield_vals_max; // yield + total_sys_err

        // Add data from gYieldStat/gYieldSys (original data)
        if (gYieldStat[iPart] && gYieldSys[iPart])
        {
            for (int i = 0; i < gYieldStat[iPart]->GetN(); ++i)
            {
                all_dNch_vals.push_back(gYieldStat[iPart]->GetX()[i]);
                all_Yield_vals_min.push_back(gYieldStat[iPart]->GetY()[i] - gYieldSys[iPart]->GetErrorYlow(i));
                all_Yield_vals_max.push_back(gYieldStat[iPart]->GetY()[i] + gYieldSys[iPart]->GetErrorYhigh(i));
            }
        }

        // Add data from Table 3
        if (iPart == 0)
        { // Xi
            for (int i = 0; i < nPoints_ppb5_Xi; ++i)
            {
                all_dNch_vals.push_back(dNch_ppb5_Xi[i]);
                all_Yield_vals_min.push_back(yield_ppb5_Xi[i] - sysErr_ppb5_Xi[i]);
                all_Yield_vals_max.push_back(yield_ppb5_Xi[i] + sysErr_ppb5_Xi[i]);
            }
        }
        else
        { // Omega
            for (int i = 0; i < nPoints_ppb5_Om; ++i)
            {
                all_dNch_vals.push_back(dNch_ppb5_Om[i]);
                all_Yield_vals_min.push_back(yield_ppb5_Om[i] - sysErr_ppb5_Om[i]);
                all_Yield_vals_max.push_back(yield_ppb5_Om[i] + sysErr_ppb5_Om[i]);
            }
        }

        double frame_min_dNch = 0, frame_max_dNch = 60;    // Defaults
        double frame_min_Yield = 0, frame_max_Yield = 0.1; // Defaults

        if (!all_dNch_vals.empty())
        {
            frame_min_dNch = TMath::MinElement(all_dNch_vals.size(), all_dNch_vals.data());
            frame_max_dNch = TMath::MaxElement(all_dNch_vals.size(), all_dNch_vals.data());
        }
        if (!all_Yield_vals_min.empty())
        {
            frame_min_Yield = TMath::MinElement(all_Yield_vals_min.size(), all_Yield_vals_min.data());
        }
        if (!all_Yield_vals_max.empty())
        {
            frame_max_Yield = TMath::MaxElement(all_Yield_vals_max.size(), all_Yield_vals_max.data());
        }

        // Add padding
        double dNchPadding = (frame_max_dNch - frame_min_dNch) * 0.1;
        if (dNchPadding == 0)
            dNchPadding = 5;
        frame_min_dNch = TMath::Max(0.0, frame_min_dNch - dNchPadding);
        frame_max_dNch = frame_max_dNch + dNchPadding;

        double yieldPadding = (frame_max_Yield - frame_min_Yield) * 0.1;
        if (yieldPadding == 0 && frame_max_Yield != 0)
            yieldPadding = TMath::Abs(frame_max_Yield * 0.2);
        else if (yieldPadding == 0)
            yieldPadding = 0.01;

        double final_min_yield = frame_min_Yield - yieldPadding;
        // Ensure y-axis starts at or below 0 if data is close to 0, for better visualization
        if (frame_min_Yield >= 0 && final_min_yield > 0 && frame_min_Yield < 0.1 * (frame_max_Yield - frame_min_Yield))
        {                                                       // if min is positive but small
            final_min_yield = TMath::Min(0.0, final_min_yield); // Try to include 0
        }
        if (frame_min_Yield >= 0 && final_min_yield > 0 && iPart == 0)
            final_min_yield = -0.005; // Xi specific
        if (frame_min_Yield >= 0 && final_min_yield > 0 && iPart == 1)
            final_min_yield = -0.0005; // Omega specific

        double final_max_yield = frame_max_Yield + yieldPadding;
        if (final_max_yield <= final_min_yield)
            final_max_yield = final_min_yield + 0.01; // Ensure max > min

        TH1 *hFrameYield_Comp = gPad->DrawFrame(frame_min_dNch, final_min_yield, frame_max_dNch, final_max_yield, ";#LTd#it{N}_{ch}/d#it{#eta}#GT_{|#it{#eta}|<0.5};#LTd#it{N}/d#it{y}#GT");
        // hFrameYield_Comp->GetYaxis()->SetTitleOffset(1.4);
        hFrameYield_Comp->GetXaxis()->SetTitleOffset(1.3);

        // Draw original data points (already styled)
        if (gYieldSys[iPart])
            gYieldSys[iPart]->Draw("E2 SAME");
        if (gYieldStat[iPart])
            gYieldStat[iPart]->Draw("PE SAME");

        // Create and style TGraphs for Table 3 data
        TGraphErrors *gYieldStat_Table3 = nullptr;
        TGraphAsymmErrors *gYieldSys_Table3 = nullptr;
        int table3_color = (iPart == 0) ? (kMagenta + 1) : (kGreen + 2); // Pick different colors from palette

        if (iPart == 0)
        { // Xi data
            gYieldStat_Table3 = new TGraphErrors(nPoints_ppb5_Xi, dNch_ppb5_Xi, yield_ppb5_Xi, dNchErr_ppb5_Xi, statErr_ppb5_Xi);
            gYieldSys_Table3 = new TGraphAsymmErrors(nPoints_ppb5_Xi);
            for (int i = 0; i < nPoints_ppb5_Xi; ++i)
            {
                gYieldSys_Table3->SetPoint(i, dNch_ppb5_Xi[i], yield_ppb5_Xi[i]);
                // For TGraphAsymmErrors, x-errors are exl, exh, y-errors are eyl, eyh
                gYieldSys_Table3->SetPointError(i, dNchErr_ppb5_Xi[i], dNchErr_ppb5_Xi[i], sysErr_ppb5_Xi[i], sysErr_ppb5_Xi[i]);
            }
        }
        else
        { // Omega data (iPart == 1)
            gYieldStat_Table3 = new TGraphErrors(nPoints_ppb5_Om, dNch_ppb5_Om, yield_ppb5_Om, dNchErr_ppb5_Om, statErr_ppb5_Om);
            gYieldSys_Table3 = new TGraphAsymmErrors(nPoints_ppb5_Om);
            for (int i = 0; i < nPoints_ppb5_Om; ++i)
            {
                gYieldSys_Table3->SetPoint(i, dNch_ppb5_Om[i], yield_ppb5_Om[i]);
                gYieldSys_Table3->SetPointError(i, dNchErr_ppb5_Om[i], dNchErr_ppb5_Om[i], sysErr_ppb5_Om[i], sysErr_ppb5_Om[i]);
            }
        }

        if (gYieldStat_Table3)
        {
            gYieldStat_Table3->SetMarkerStyle(kFullDiamond); // Different marker style
            gYieldStat_Table3->SetMarkerColor(table3_color);
            gYieldStat_Table3->SetLineColor(table3_color);
            gYieldStat_Table3->SetMarkerSize(1.5);
        }
        if (gYieldSys_Table3)
        {
            // gYieldSys_Table3->SetFillColorAlpha(table3_color, 0.35);
            gYieldSys_Table3->SetLineColor(table3_color);
            gYieldSys_Table3->SetFillStyle(0); // Solid fill, could use a pattern like 3004 for distinction
        }

        // Draw Table 3 data
        if (gYieldSys_Table3)
            gYieldSys_Table3->Draw("E2 SAME"); // Systematics boxes
        if (gYieldStat_Table3)
            gYieldStat_Table3->Draw("PE SAME"); // Statistical error points

        // Add a comprehensive legend
        TLegend *legYield_Comp = new TLegend(0.13, 0.75, 0.4, 0.85); // Adjusted for more entries
        legYield_Comp->SetBorderSize(0);
        legYield_Comp->SetFillStyle(0); // Transparent background
        legYield_Comp->SetTextSize(0.025);
        // legYield_Comp->SetHeader(Form("%s, -0.5 < y < 0", particleTitlePlot[iPart]));

        if (gYieldStat[iPart])
            legYield_Comp->AddEntry(gYieldStat[iPart], "p-Pb #sqrt{s_{NN}} = 8.16 TeV, This work", "pe");
        // if (gYieldSys[iPart]) legYield_Comp->AddEntry(gYieldSys[iPart], "This work (sys.)", "f");
        if (gYieldStat_Table3)
            legYield_Comp->AddEntry(gYieldStat_Table3, "p-Pb #sqrt{s_{NN}} = 5.02 TeV, Phys.Lett.B 758 (2016), 389-401", "pe");
        // if (gYieldSys_Table3) legYield_Comp->AddEntry(gYieldSys_Table3, "IntYields PDF Table 3 (sys.)", "f");
        legYield_Comp->Draw();

        // Add labels (similar to original plot)
        TLatex latexComp;
        latexComp.SetNDC();
        latexComp.SetTextSize(0.028);
        latexComp.DrawLatex(0.2, 0.86, (Form("%s, -0.5 < y < 0", particleTitlePlot[iPart])));

        // latexComp.DrawLatex(0.18, 0.88, Form("p-Pb #sqrt{s_{NN}} = 8.16 TeV, %s", particleTitlePlot[iPart]));
        // latexComp.DrawLatex(0.18, 0.84, "-0.5 < y < 0");
        latexComp.SetTextSize(0.023);
        latexComp.SetTextFont(42);
        latexComp.DrawLatex(0.2, 0.72, "Uncertainties: stat.(bars), sys.(boxes)");
        // latexComp.DrawLatex(0.18, 0.80, "Data: This work (circles), IntYields PDF (squares)");

        cYield_Comp->Update();
        SaveImage(outputPlotFolder, "YieldVsMult_Comparison",
                  Form("YieldVsMult_compare5TeV_%s", particleNamePlot[iPart]),
                  imageFormat.Data(), cYield_Comp);

        // Cleanup for this iteration
        delete gYieldStat_Table3;
        delete gYieldSys_Table3;
        delete legYield_Comp;
        delete cYield_Comp; // Delete the canvas for this particle's comparison plot
    }
    // --- END OF NEW SECTION ---
    // --- Generate Mean pT vs Multiplicity ---
    TCanvas *cMeanPt = new TCanvas("cMeanPtVsMult", "<pT> vs dNch/deta", 1200, 800);
    TGraphErrors *gMeanPtStat[nParticlesPlot] = {nullptr};
    TGraphAsymmErrors *gMeanPtSys[nParticlesPlot] = {nullptr};
    // TGraphErrors *gMeanPt[nParticlesPlot] = {nullptr};
    TLegend *legMeanPt = new TLegend(0.8, 0.2, 0.9, 0.3);
    legMeanPt->SetBorderSize(0);
    legMeanPt->SetFillStyle(0);
    // double globalMinPt = 1e6, globalMaxPt = -1e6;
    // double globalMinDNch = 1e6, globalMaxDNch = -1e6;
    double globalMinPt = 1e6, globalMaxPt = -1e6, globalMinDNch_Pt = 1e6, globalMaxDNch_Pt = -1e6;

    for (int iPart = 0; iPart < nParticlesPlot; ++iPart)
    {
        int nMult = *nMultBinsPlot[iPart];
        std::vector<double> v_dNch, v_dNchErrX, v_MeanPt;
        std::vector<double> v_MeanPtStatErr, v_TotalMeanPtSysLo, v_TotalMeanPtSysHi; // Total sys
        std::vector<double> v_dNchErrSysXLo, v_dNchErrSysXHi;
        for (int iMult = 0; iMult < nMult; ++iMult)
        {
            TString funcNameLevy = GetFitFunctionName(kLevyTsallis, particleNamePlot[iPart], iMult);
            TString yieldInfoName = Form("%s/MultiplicityBinned/YieldInfo_%s_Mult%d_%s",
                                         particleDirName[iPart], particleNamePlot[iPart], iMult, funcNameLevy.Data());
            TH1 *hYieldInfo = (TH1 *)resultsFile->Get(yieldInfoName.Data());
            if (!hYieldInfo)
                continue;

            Double_t meanPtVal, statErrVal, sysLoSpectra, sysHiSpectra;
            GetIndividualErrors(hYieldInfo, kMean, kMeanStat, kMeanSysLo, kMeanSysHi,
                                meanPtVal, statErrVal, sysLoSpectra, sysHiSpectra);

            // Combine sys errors for Mean pT 
            Double_t currentSysExtrapMeanPt = sysExtrapMeanPt[iPart][iMult];
            Double_t totalSysLo = TMath::Sqrt(TMath::Power(sysLoSpectra, 2) + TMath::Power(currentSysExtrapMeanPt, 2));
            Double_t totalSysHi = TMath::Sqrt(TMath::Power(sysHiSpectra, 2) + TMath::Power(currentSysExtrapMeanPt, 2));

            v_dNch.push_back(dNch[iPart][iMult]);
            v_dNchErrX.push_back(dNchErr[iPart][iMult]);
            v_MeanPt.push_back(meanPtVal);
            // v_MeanPtErr.push_back(totalErr);
            v_MeanPtStatErr.push_back(statErrVal);
            // Store combined sys errors
            v_TotalMeanPtSysLo.push_back(totalSysLo);
            v_TotalMeanPtSysHi.push_back(totalSysHi);
            v_dNchErrSysXLo.push_back(dNchErr[iPart][iMult]);
            v_dNchErrSysXHi.push_back(dNchErr[iPart][iMult]);
            // if (dNch[iPart][iMult] < globalMinDNch)
            //     globalMinDNch = dNch[iPart][iMult];
            // if (dNch[iPart][iMult] > globalMaxDNch)
            //     globalMaxDNch = dNch[iPart][iMult];
            // if (meanPt - totalErr < globalMinPt)
            //     globalMinPt = meanPt - totalErr;
            // if (meanPt + totalErr > globalMaxPt)
            //     globalMaxPt = meanPt + totalErr;
            // Update global ranges considering total systematic error
            if (dNch[iPart][iMult] < globalMinDNch_Pt)
                globalMinDNch_Pt = dNch[iPart][iMult];
            if (dNch[iPart][iMult] > globalMaxDNch_Pt)
                globalMaxDNch_Pt = dNch[iPart][iMult];
            if (meanPtVal - totalSysLo < globalMinPt)
                globalMinPt = meanPtVal - totalSysLo;
            if (meanPtVal + totalSysHi > globalMaxPt)
                globalMaxPt = meanPtVal + totalSysHi;
        } // end multiplicity loop

        if (v_dNch.empty())
            continue;
        gMeanPtStat[iPart] = new TGraphErrors(v_dNch.size(), v_dNch.data(), v_MeanPt.data(), v_dNchErrX.data(), v_MeanPtStatErr.data());
        // Use combined sys errors for TGraphAsymmErrors
        gMeanPtSys[iPart] = new TGraphAsymmErrors(v_dNch.size(), v_dNch.data(), v_MeanPt.data(),
                                                  v_dNchErrSysXLo.data(), v_dNchErrSysXHi.data(),
                                                  v_TotalMeanPtSysLo.data(), v_TotalMeanPtSysHi.data());
        // gMeanPt[iPart] = new TGraphErrors(v_dNch.size(), v_dNch.data(), v_MeanPt.data(), v_dNchErr.data(), v_MeanPtErr.data());
        gMeanPtStat[iPart]->SetName(Form("gMeanPtStat_%s", particleNamePlot[iPart]));
        gMeanPtStat[iPart]->SetTitle(TString::Format("%s: Mean #it{p}_{T}", particleTitlePlot[iPart]));
        gMeanPtStat[iPart]->SetMarkerStyle(kFullCircle); // Use styles from CascadeUtils
        // if (iPart == 0)
        // {
        gMeanPtStat[iPart]->SetMarkerColor(customCascPalette[iPart]);
        gMeanPtStat[iPart]->SetLineColor(customCascPalette[iPart]);
        // }
        // else
        // {
        // gMeanPtStat[iPart]->SetMarkerColor(kAzure + 1);
        // gMeanPtStat[iPart]->SetLineColor(kAzure + 1);
        // }
        // gMeanPtStat[iPart]->SetMarkerSize(1.5);

        gMeanPtSys[iPart]->SetName(Form("gMeanPtSys_%s", particleNamePlot[iPart]));
        // if (iPart == 0)
        // {
        // gMeanPtSys[iPart]->SetFillColorAlpha(kMagenta + 2, 0.35);
        gMeanPtSys[iPart]->SetLineColor(customCascPalette[iPart]);
        // }
        // else
        // {
        gMeanPtSys[iPart]->SetFillColor(customCascPalette[iPart]);
        // gMeanPtSys[iPart]->SetLineColor(kAzure + 1);
        // }
        gMeanPtSys[iPart]->SetFillStyle(0);
        legMeanPt->AddEntry(gMeanPtStat[iPart], particleTitlePlot[iPart], "pef");
    } // end particle loop

    cMeanPt->cd();
    // gPad->SetLeftMargin(0.15);
    // gPad->SetBottomMargin(0.12);
    globalMinDNch_Pt -= 5;
    if (globalMinDNch_Pt < 0)
        globalMinDNch_Pt = 0;
    globalMaxDNch_Pt += 5;
    globalMinPt *= 0.9;
    if (globalMinPt == 0)
        globalMinPt = -0.1 * globalMaxPt; // Handle if minPt is 0
    globalMaxPt *= 1.1;
    if (globalMaxPt == 0)
        globalMaxPt = 1;

    TH1 *hFrameMeanPt = gPad->DrawFrame(globalMinDNch_Pt, globalMinPt, globalMaxDNch_Pt, globalMaxPt, ";#LTd#it{N}_{ch}/d#it{#eta}#GT_{|#it{#eta}|<0.5};#LT#it{p}_{T}#GT (GeV/#it{c})");
    // hFrameMeanPt->GetYaxis()->SetTitleOffset(1.8);
    gPad->SetRightMargin(0.05);
    gPad->SetTopMargin(0.05);
    gPad->SetBottomMargin(0.12);
    TLatex latex;
    latex.SetNDC();
    latex.SetTextSize(0.025);
    latex.DrawLatex(0.15, 0.86, "p-Pb #sqrt{s_{NN}} = 8.16 TeV, This work");
    latex.DrawLatex(0.15, 0.82, Form("-0.5 < y < 0"));
    latex.SetTextSize(0.022);
    latex.SetTextFont(42);
    latex.DrawLatex(0.15, 0.78, "Uncertainties: stat.(bars), sys.(boxes)");
    cYield->Update();
    hFrameMeanPt->GetXaxis()->SetTitleOffset(1.3);
    for (int iPart = 0; iPart < nParticlesPlot; ++iPart)
    {
        if (gMeanPtSys[iPart])
            gMeanPtSys[iPart]->Draw("E2 SAME");
        if (gMeanPtStat[iPart])
            gMeanPtStat[iPart]->Draw("PE1 SAME");
    }
    legMeanPt->Draw();
    cMeanPt->Update();
    SaveImage(outputPlotFolder, "MeanPtVsMult", "MeanPtVsMult", imageFormat.Data(), cMeanPt);

    // --- Generate Spectra Fit Plots for ALL bins ---
    // --- Using fits retrieved from log files ---
    Info("IntYieldsMeanPtCompare", "Generating spectra fit plots for all multiplicity bins...");
    const Double_t pT_min_draw = 0.0;
    const Double_t pT_max_draw = 7.0;

    for (int iPart = 0; iPart < nParticlesPlot; ++iPart)
    {
        int nMult = *nMultBinsPlot[iPart];
        for (int iMult = 0; iMult < nMult; ++iMult) // Loop through ALL multiplicity bins
        {
            Double_t multLow = multBinEdgesPlot[iPart][iMult];
            Double_t multHigh = multBinEdgesPlot[iPart][iMult + 1];
            Info("IntYieldsMeanPtCompare", "  Plotting spectra for %s, mult %.0f-%.0f %%", particleNamePlot[iPart], multLow, multHigh);

            TString canvasName = Form("cSpectraFit_%s_Mult%d", particleNamePlot[iPart], iMult);
            TString canvasTitle = Form("%s Spectra Fit (%.0f-%.0f%%)", particleTitlePlot[iPart], multLow, multHigh);
            TCanvas *cSpecFit = new TCanvas(canvasName.Data(), canvasTitle.Data(), 2000, 2000);
            TPad *pad1 = new TPad("pad1", "Main", 0.0, 0.3, 0.99, 0.99);
            TPad *pad2 = new TPad("pad2", "Ratio", 0.0, 0.0, 0.99, 0.3);
            pad1->SetLogy();
            pad1->SetRightMargin(0.02);
            pad1->SetTopMargin(0.02);
            pad1->SetBottomMargin(0.01);
            pad2->SetTopMargin(0.01);
            pad2->SetBottomMargin(0.2);
            pad2->SetRightMargin(0.02);
            pad2->SetGridy();
            pad1->Draw();
            pad2->Draw();

            pad1->cd(); // Activate top pad

            // Get original histograms from spectra file
            TString histNameStat = Form("stat_effCorrPt_%s_mult[%d]", particleNamePlot[iPart], iMult);
            TString histNameSys = Form("sys_effCorrPt_%s_mult[%d]", particleNamePlot[iPart], iMult);
            TH1D *hStat = (TH1D *)spectraFile->FindObjectAny(histNameStat.Data());
            TH1D *hSys = (TH1D *)spectraFile->FindObjectAny(histNameSys.Data());

            if (!hStat || !hSys)
            {
                Error("IntYieldsMeanPtCompare", "Spectra hist not found for %s mult %d. Skipping plot.", particleNamePlot[iPart], iMult);
                delete cSpecFit;
                continue;
            }

            std::vector<double> v_pt, v_ptErr_stat, v_val, v_statErr, v_sysErr; // ptErr_stat is 0 for vertical bars
            // ptErr_sys will be half bin width for sys boxes
            int nDataPoints = 0;
            for (int iBin = 1; iBin <= hStat->GetNbinsX(); ++iBin)
            {
                if (hStat->GetBinContent(iBin) == 0 && hStat->GetBinError(iBin) == 0)
                    continue; // Skip empty bins
                v_pt.push_back(hStat->GetBinCenter(iBin));
                v_ptErr_stat.push_back(hStat->GetBinWidth(iBin) / 2.0); // No horizontal stat error bars, or hStat->GetBinWidth(iBin)/2.0 if desired
                v_val.push_back(hStat->GetBinContent(iBin));
                v_statErr.push_back(hStat->GetBinError(iBin));
                v_sysErr.push_back(hSys->GetBinError(iBin)); // Assuming hSys contains symmetric sys errors
                nDataPoints++;
            }
            if (nDataPoints == 0)
            {
                Error("IntYieldsMeanPtCompare", "No data points in %s", histNameStat.Data());
                delete cSpecFit;
                continue;
            }

            TGraphErrors *gDataPointsStat = new TGraphErrors(nDataPoints, v_pt.data(), v_val.data(), v_ptErr_stat.data(), v_statErr.data());
            TGraphAsymmErrors *gDataBoxesSys = new TGraphAsymmErrors(nDataPoints);
            for (int i = 0; i < nDataPoints; ++i)
            {
                double binWidth = hStat->GetXaxis()->GetBinWidth(hStat->GetXaxis()->FindBin(v_pt[i]));
                gDataBoxesSys->SetPoint(i, v_pt[i], v_val[i]);
                gDataBoxesSys->SetPointError(i, binWidth / 2.0, binWidth / 2.0, v_sysErr[i], v_sysErr[i]); // x-err is half bin width
            }

            // Style data points and boxes
            gDataPointsStat->SetMarkerStyle(kFullCircle);
            gDataPointsStat->SetMarkerSize(1.5);
            gDataPointsStat->SetMarkerColor(customCascPalette[0]); // Use colors from customCascPalette
            gDataPointsStat->SetLineColor(customCascPalette[0]);

            gDataBoxesSys->SetFillColor(customCascPalette[0]); // Light gray for sys boxes
            gDataBoxesSys->SetLineColor(customCascPalette[0]);
            gDataBoxesSys->SetFillStyle(0); // Solid fill

            // Combine stat and sys for plotting errors
            // TH1D *hDataToPlot = (TH1D *)hStat->Clone(Form("hDataToPlot_%s_Mult%d", particleNamePlot[iPart], iMult));
            // hDataToPlot->SetDirectory(0); // Avoid attachment to file
            // for (int iBin = 1; iBin <= hDataToPlot->GetNbinsX(); ++iBin)
            // {
            //     double statErr = hStat->GetBinError(iBin);
            //     double sysErr = hSys->GetBinError(iBin);
            //     hDataToPlot->SetBinError(iBin, TMath::Sqrt(statErr * statErr + sysErr * sysErr));
            // }

            // // Style data points
            // hDataToPlot->SetMarkerStyle(kFullCircle);
            // hDataToPlot->SetMarkerSize(1.0);
            // hDataToPlot->SetMarkerColor(kBlack);
            // hDataToPlot->SetLineColor(kBlack);
            // hDataToPlot->SetTitle(Form("%s (%.0f - %.0f %%)", particleTitlePlot[iPart], multLow, multHigh));
            // hDataToPlot->GetYaxis()->SetTitle("#frac{1}{N_{ev}} #frac{d^{2}N}{dp_{T}dy} (GeV/#it{c})^{-1}"); // Adjust axis titles
            // hDataToPlot->GetYaxis()->SetTitleOffset(1.4);
            // hDataToPlot->GetXaxis()->SetLabelSize(0); // Hide labels on top pad
            // hDataToPlot->GetXaxis()->SetTitleSize(0);
            // hDataToPlot->Draw("PE"); // Draw data points

            // TLegend *legSpecFit = new TLegend(0.55, 0.55, 0.88, 0.88); // Adjust legend position
            // legSpecFit->SetBorderSize(0);
            // legSpecFit->SetFillStyle(0);
            // legSpecFit->AddEntry(hDataToPlot, "Data", "pe");

            // Draw frame with new pT range, then boxes, then stat points
            // Determine Y range dynamically
            double yMin = 1e-5, yMax = 1e2; // Default, will be adjusted
            if (gDataPointsStat->GetN() > 0)
            {
                yMin = TMath::MinElement(gDataPointsStat->GetN(), gDataPointsStat->GetY()) * 0.1;
                yMax = TMath::MaxElement(gDataPointsStat->GetN(), gDataPointsStat->GetY()) * 10.0;
                for (int k = 0; k < gDataPointsStat->GetN(); ++k)
                {
                    if (gDataPointsStat->GetY()[k] + gDataBoxesSys->GetErrorYhigh(k) > yMax)
                        yMax = gDataPointsStat->GetY()[k] + gDataBoxesSys->GetErrorYhigh(k);
                }
                if (yMin <= 0)
                    yMin = 1e-5 * yMax; // Ensure log scale compatibility
                if (yMax <= 0)
                    yMax = 10; // ensure a positive range if all values are zero or negative
            }
            if (yMin >= yMax)
            {
                yMin = 0.001;
                yMax = 100;
            } // Fallback if dynamic range fails

            TH1 *hFrameSpectra = gPad->DrawFrame(pT_min_draw, yMin, pT_max_draw, yMax);
            // hFrameSpectra->SetTitle(Form("%s (%.0f - %.0f %%)", particleTitlePlot[iPart], multLow, multHigh));
            hFrameSpectra->GetYaxis()->SetTitle("#frac{1}{#it{N}_{ev}} #frac{d#it{N}}{d#it{p}_{T}} (GeV/#it{c})^{-1}"); // Example title
            hFrameSpectra->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");                                             // Title will be on ratio plot
            hFrameSpectra->GetXaxis()->SetLabelSize(0);
            hFrameSpectra->GetYaxis()->SetTitleOffset(1.2);

            if (gDataBoxesSys)
                gDataBoxesSys->Draw("E2 SAME"); // Draw sys boxes
            if (gDataPointsStat)
                gDataPointsStat->Draw("PE1 SAME"); // Draw stat points + errors

            TLegend *legSpecFit = new TLegend(0.75, 0.7, 0.9, 0.85);
            legSpecFit->SetBorderSize(0);
            legSpecFit->SetFillStyle(0);
            // Add entry for combined data (stat+sys visually)
            if (gDataPointsStat)
                legSpecFit->AddEntry(gDataPointsStat, "Data (stat #oplus sys)", "pef"); // "f" for fill from sys box

            TF1 *levyFitFunc = nullptr; // To store the Levy fit for the ratio plot

            // Loop through fit functions to retrieve from log files
            for (int iFunc = 0; iFunc < kNumFitFunctions; ++iFunc)
            {
                TString funcName = GetFitFunctionName(iFunc, particleNamePlot[iPart], iMult);
                TString logFileName = Form("%s/Log_%s.root", inputLogFolder.Data(), funcName.Data());

                TFile *logFile = TFile::Open(logFileName.Data());
                if (!logFile || logFile->IsZombie())
                {
                    Warning("IntYieldsMeanPtCompare", "Could not open log file: %s", logFileName.Data());
                    if (logFile)
                        delete logFile;
                    continue;
                }

                // Construct the name of the histogram as saved in the log file by YieldMean
                TString baseHistName = Form("stat_effCorrPt_%s_mult[%d]", particleNamePlot[iPart], iMult);
                TString loggedHistName = GetLogHistogramName(baseHistName.Data(), funcName.Data(), particleNamePlot[iPart], iMult);

                TH1D *hLogged = (TH1D *)logFile->Get(loggedHistName.Data());
                TF1 *fFitted = nullptr;

                if (hLogged)
                {
                    // Retrieve the function attached to the histogram
                    // Info("IntYieldsMeanPtCompare", "Found histogram '%s' in log file '%s'", loggedHistName.Data(), logFileName.Data());
                    fFitted = hLogged->GetFunction(Form("%s_%s_Mult%d", funcName.Data(), particleNamePlot[iPart], iMult));
                }
                else
                {
                    Warning("IntYieldsMeanPtCompare", "Histogram '%s' not found in log file '%s'", loggedHistName.Data(), logFileName.Data());
                }

                if (fFitted)
                {
                    // Successfully retrieved the function
                    // fFitted->SetLineStyle((iFunc == kLevyTsallis) ? kSolid : kDashed + iFunc - 1); // Style lines
                    fFitted->SetLineColor(customCascPalette[iFunc + 1]); // Use colors from customCascPalette
                    fFitted->SetLineWidth(2);
                    // Ensure fit function is drawn across the new extended range **
                    TF1 *fDrawClone = (TF1 *)fFitted->Clone(Form("%s_draw", fFitted->GetName()));
                    fDrawClone->SetRange(pT_min_draw, pT_max_draw); // Explicitly set range for drawing
                    fDrawClone->DrawCopy("SAME");
                    delete fDrawClone;

                    // Colors should be set already by the constructors in ExtraUtils.C
                    // fFitted->DrawCopy("SAME"); // Draw a COPY of the function

                    if (iFunc == kLevyTsallis)
                        levyFitFunc = (TF1 *)fFitted->Clone(); // Clone Levy for ratio plot

                    // Add legend entry
                    TString baseFuncName = ""; // Base name for the TF1 type
                    if (iFunc == kLevyTsallis)
                        baseFuncName = "L#acute{e}vy-Tsallis";
                    else if (iFunc == kBoltzmann)
                        baseFuncName = "Boltzmann";
                    else if (iFunc == kBlastWave)
                        baseFuncName = "BG BlastWave";
                    // else if (iFunc == kmTScaling)
                    //     baseFuncName = "m_{T}-Scaling";
                    // else if (iFunc == kBoseEinstein)
                    //     baseFuncName = "Bose-Einstein";
                    // else if (iFunc == kFermiDirac)
                    //     baseFuncName = "Fermi-Dirac";

                    // double chi2 = fFitted->GetChisquare();
                    // double ndf = fFitted->GetNDF();
                    // double chi2NDF = (ndf > 0) ? chi2 / ndf : 0.0; // Avoid division by zero
                    // if (ndf > 0)
                    // legSpecFit->AddEntry(fFitted, Form("%s (#chi^{2}/NDF = %.2f)", baseFuncName.Data(), chi2NDF), "l");
                    // else
                    legSpecFit->AddEntry(fFitted, Form("%s", baseFuncName.Data()), "l");

                    // double chi2 = fFitted->GetChisquare();
                    // double ndf = fFitted->GetNDF();
                    // if (ndf > 0)
                    // legSpecFit->AddEntry(fFitted, Form("%s (#chi^{2}/NDF = %.2f)", baseFuncName.Data(), chi2 / ndf), "l");
                    // else
                    // legSpecFit->AddEntry(fFitted, Form("%s", baseFuncName.Data()), "l");

                    // Info("IntYieldsMeanPtCompare", "Fit function '%s' found for histogram '%s' in log file '%s'", funcName.Data(), loggedHistName.Data(), logFileName.Data());
                }
                else
                {
                    Warning("IntYieldsMeanPtCompare", "Fit function '%s' not found for histogram '%s' in log file '%s'", funcName.Data(), loggedHistName.Data(), logFileName.Data());
                }

                pad2->cd();
                {
                    TH1 *hRatio = nullptr;

                    // Use temporary TH1 from original hStat for correct bin definitions for ratio
                    TH1D *hAxisDefForRatio = (TH1D *)hStat->Clone("hAxisDefForRatio");
                    hAxisDefForRatio->Reset(); // Keep axes, clear content
                    // Use CreateRatioPlotFromGraph 
                    hRatio = CreateRatioPlotFromGraph(gDataPointsStat, fFitted, hAxisDefForRatio);
                    delete hAxisDefForRatio;

                    if (hRatio)
                    {
                        Info("IntYieldsMeanPtCompare", "Creating ratio plot for %s, mult %d, fit %s", particleNamePlot[iPart], iMult, funcName.Data());
                        // Set X-axis range for ratio plot
                        hRatio->GetXaxis()->SetRangeUser(pT_min_draw, pT_max_draw);
                        hRatio->GetXaxis()->SetLabelSize(0.1);
                        hRatio->GetXaxis()->SetTitleSize(0.1);
                        // hRatio->GetXaxis()->SetTitleOffset(0.9);
                        hRatio->GetYaxis()->SetLabelSize(0.07);
                        hRatio->GetYaxis()->SetTitleSize(0.08);
                        hRatio->GetYaxis()->SetTitleOffset(0.4);
                        // hRatio->GetYaxis()->SetTitle("Data / Levy-Tsallis Fit"); // Title for ratio y-axis
                        hRatio->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");                       // Add title to ratio x-axis
                        hRatio->SetColors(customCascPalette[iFunc + 1], customCascPalette[iFunc + 1]); // Use colors from customCascPalette
                        hRatio->Draw("PE same");                                                       // Draw points on ratio
                        // TLine *line = new TLine(pT_min_draw, 1, pT_max_draw, 1); // Line at 1
                        // line->SetLineStyle(kDashed);
                        // line->SetLineColor(kGray + 2);
                        // line->Draw("same");
                    }

                    // logFile->Close();
                    // delete logFile;
                }
                pad1->cd(); // Back to main pad

            } // End loop over fit functions
            legSpecFit->Draw();

            // // Get Extrapolation fraction from main results file
            // double extrapFrac = 0.0;
            // TString funcNameLevy = GetFitFunctionName(kLevyTsallis, particleNamePlot[iPart], iMult);
            // TString yieldInfoName = Form("%s/MultiplicityBinned/YieldInfo_%s_Mult%d_%s",
            //                              particleDirName[iPart], particleNamePlot[iPart], iMult, funcNameLevy.Data());
            // TH1 *hYieldInfo = (TH1 *)resultsFile->Get(yieldInfoName.Data());
            // if (hYieldInfo)
            // {
            //     extrapFrac = hYieldInfo->GetBinContent(kExtra); // Get from kExtra bin
            // }
            // // else
            // // {
            // //     Warning("IntYieldsMeanPtCompare", "YieldInfo histogram not found for Extrapolation fraction: %s", yieldInfoName.Data());
            // // }

            // // Draw labels and extrapolation fraction
            // latex.SetTextSize(0.030);
            // latex.DrawLatexNDC(0.18, 0.85, "ALICE Preliminary");
            // latex.DrawLatexNDC(0.18, 0.80, "p-Pb #sqrt{s_{NN}} = 8.16 TeV");
            // latex.DrawLatexNDC(0.18, 0.75, "-0.5 < y < 0");
            // if (hYieldInfo)
            //     latex.DrawLatexNDC(0.18, 0.70, Form("Extrap. Frac: %.2f %%", extrapFrac * 100.0 / hYieldInfo->GetBinContent(kYield)));

            // // --- Ratio Plot ---
            // pad2->cd(); // Activate bottom pad
            // TH1 *hRatio = nullptr;
            // if (levyFitFunc)
            // { // Use the cloned Levy fit function
            //     // Info("IntYieldsMeanPtCompare", "Creating ratio plot for %s, mult %d", particleNamePlot[iPart], iMult);
            //     hRatio = CreateRatioPlot(hDataToPlot, levyFitFunc);
            //     // Info("IntYieldsMeanPtCompare", "Ratio plot created for %s, mult %d", particleNamePlot[iPart], iMult);
            //     if (hRatio)
            //     {
            //         // Style ratio plot axes
            //         hRatio->GetXaxis()->SetLabelSize(0.1);
            //         hRatio->GetXaxis()->SetTitleSize(0.12);
            //         hRatio->GetXaxis()->SetTitleOffset(0.9);
            //         hRatio->GetYaxis()->SetLabelSize(0.08);
            //         hRatio->GetYaxis()->SetTitleSize(0.1);
            //         hRatio->GetYaxis()->SetTitleOffset(0.5);
            //         hRatio->Draw("PE");
            //         // Info("IntYieldsMeanPtCompare", "Ratio plot drawn for %s, mult %d", particleNamePlot[iPart], iMult);

            //         // Draw line at 1 for reference
            //         TLine *line = new TLine(hRatio->GetXaxis()->GetXmin(), 1, hRatio->GetXaxis()->GetXmax(), 1);
            //         line->SetLineStyle(kDashed);
            //         line->SetLineColor(kGray + 2);
            //         line->Draw("same");
            //         // Info("IntYieldsMeanPtCompare", "Reference line drawn for %s, mult %d", particleNamePlot[iPart], iMult);
            //     }
            // }
            // else
            // {
            //     Warning("IntYieldsMeanPtCompare", "Levy-Tsallis fit was not found for ratio plot (%s, mult %d).", particleNamePlot[iPart], iMult);
            // }
            double extrapFrac = 0.0, totalYieldForFrac = 1.0;
            TString funcNameLevy = GetFitFunctionName(kLevyTsallis, particleNamePlot[iPart], iMult);
            TString yieldInfoName = Form("%s/MultiplicityBinned/YieldInfo_%s_Mult%d_%s", particleDirName[iPart], particleNamePlot[iPart], iMult, funcNameLevy.Data());
            TH1 *hYieldInfo = (TH1 *)resultsFile->Get(yieldInfoName.Data());
            if (hYieldInfo)
            {
                extrapFrac = hYieldInfo->GetBinContent(kExtra);
                totalYieldForFrac = hYieldInfo->GetBinContent(kYield);
            }
            TLatex latex;
            latex.SetTextSize(0.025);
            latex.DrawLatexNDC(0.15, 0.92, "p-Pb #sqrt{s_{NN}} = 8.16 TeV, This work");
            latex.DrawLatexNDC(0.15, 0.88, Form("%s, -0.5 < y < 0", particleTitlePlot[iPart]));
            latex.DrawLatexNDC(0.15, 0.84, Form("%.0f-%.0f%% V0A", multLow, multHigh));
            latex.SetTextFont(42);
            latex.SetTextSize(0.022);
            latex.DrawLatexNDC(0.15, 0.81, "Uncertainties: stat.(bars), sys.(boxes)");
            if (hYieldInfo && totalYieldForFrac != 0)
                latex.DrawLatexNDC(0.15, 0.78, Form("Extrapolated: %.2f%%", extrapFrac * 100.0 / totalYieldForFrac));

            // Save and Cleanup for this multiplicity bin
            cSpecFit->cd();
            // legSpecFit->Draw();
            cSpecFit->Modified();
            cSpecFit->ForceUpdate();
            // TString figLabel = (iPart == 0) ? "Fig5_5" : "Fig5_6";
            // Info("IntYieldsMeanPtCompare", "Saving spectra fit plot for %s, mult %d", particleNamePlot[iPart], iMult);
            SaveImage(outputPlotFolder, "SpectraFitsPerMult", TString::Format("%s_Mult%d", particleNamePlot[iPart], iMult), imageFormat.Data(), cSpecFit);
            // Info("IntYieldsMeanPtCompare", "Saved spectra fit plot for %s, mult %d", particleNamePlot[iPart], iMult);

            // delete hDataToPlot;
            // if (hRatio)
            //     delete hRatio;
            if (levyFitFunc)
                delete levyFitFunc; // Delete the cloned function
            delete cSpecFit;        // Delete canvas
        } // End loop over multiplicity bins
    } // End loop over particles

    // --- Calculate and Print Extrapolation Uncertainty ---
    printf("\n--- Systematic Uncertainty from Fit Function Choice ---\n");
    printf("Particle | Mult Bin | Yield RMS/Mean (Levy) | Yield MaxDev/Mean (Levy) | MeanPt RMS/Mean (Levy) | MeanPt MaxDev/Mean (Levy)\n");
    printf("---------------------------------------------------------------------------------------------------------------------------\n");
    for (int iPart = 0; iPart < nParticlesPlot; ++iPart)
    {
        int nMult = *nMultBinsPlot[iPart];
        for (int iMult = 0; iMult < nMult; ++iMult)
        {
            std::vector<double> yields;
            std::vector<double> means;
            double levyYield = -1;
            double levyMean = -1;
            for (int iFunc = 0; iFunc < kNumFitFunctions; ++iFunc)
            {
                int index = iMult * kNumFitFunctions + iFunc;
                double currentYield = yieldsPerFuncAll[iPart][index];
                double currentMean = meansPerFuncAll[iPart][index];
                if (currentYield < 0)
                    continue; // Skip if YieldMean failed
                yields.push_back(currentYield);
                means.push_back(currentMean);
                if (iFunc == kLevyTsallis)
                {
                    levyYield = currentYield;
                    levyMean = currentMean;
                }
            }
            if (levyYield <= 0 || yields.size() < 2)
                continue;
            double sumSqYieldDev = 0;
            double maxYieldDev = 0;
            double sumSqMeanDev = 0;
            double maxMeanDev = 0;
            for (double y : yields)
            {
                double dev = y - levyYield;
                sumSqYieldDev += dev * dev;
                if (TMath::Abs(dev) > maxYieldDev)
                    maxYieldDev = TMath::Abs(dev);
            }
            for (double m : means)
            {
                double dev = m - levyMean;
                sumSqMeanDev += dev * dev;
                if (TMath::Abs(dev) > maxMeanDev)
                    maxMeanDev = TMath::Abs(dev);
            }
            double rmsYield = TMath::Sqrt(sumSqYieldDev / yields.size());
            double rmsMean = TMath::Sqrt(sumSqMeanDev / means.size());
            double sysYieldExtrap = 0.5 * maxYieldDev;
            double sysMeanExtrap = 0.5 * maxMeanDev; // Half max deviation
            printf("%-8s | %-8d | %-23.4f | %-26.4f | %-22.4f | %-25.4f\n",
                   particleNamePlot[iPart], iMult, rmsYield / levyYield, sysYieldExtrap / levyYield, rmsMean / levyMean, sysMeanExtrap / levyMean);
        }
    }
    printf("---------------------------------------------------------------------------------------------------------------------------\n");

    // --- Cleanup ---
    resultsFile->Close();
    spectraFile->Close();
    delete resultsFile;
    delete spectraFile;
    // Clean up TGraphErrors
    for (int i = 0; i < nParticlesPlot; ++i)
    {
        delete gYieldStat[i];
        delete gYieldSys[i];
        delete gMeanPtStat[i];
        delete gMeanPtSys[i];
        if (linearYieldFits[i])
        {
            delete linearYieldFits[i];
            linearYieldFits[i] = nullptr;
        }
    }

    Info("IntYieldsMeanPtCompare", "Plots saved in: %s", outputPlotFolder.Data());
    Info("IntYieldsMeanPtCompare", "Log files were read from: %s", inputLogFolder.Data());
}