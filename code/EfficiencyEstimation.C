#include <TROOT.h>
#include <TStyle.h>
#include "CascadeUtils.h"

struct RunInfo
{
    TString runName;
    std::vector<Double_t> effVal;                           // [MC Run][ptbin] : efficiency values in mult integrated (0-100%) Combined (+-) case
    std::vector<Double_t> effWeight;                        // [MC Run][ptbin] : efficiency weights in mult integrated (0-100%) Combined (+-) case
    std::vector<std::vector<Double_t>> effValCommonMult;    // [MC Run][new mult bin][ptbin] : efficiency values in new common multiplicity bins Combined (+-) case
    std::vector<std::vector<Double_t>> effWeightCommonMult; // [MC Run][new mult bin][ptbin] : efficiency weights in new common multiplicity bins Combined (+-) case
};

const char *ParseRunName(TString runName)
{
    // Extract the run name from the filename
    if (runName.Contains("LHC17e1b"))
        return "Trig DPMJET";
    else if (runName.Contains("LHC17e1a"))
        return "Trig DPMJET";
    else if (runName.Contains("LHC18f3b"))
        return "GP DPMJET";
    else if (runName.Contains("LHC17f3a"))
        return "GP EPOS-LHC";
    else if (runName.Contains("LHC17l7a2"))
        return "Inj DPMJET";
    else
        return "Unknown_Run";
};

/**
 * @brief Estimates the efficiency of Xi and Omega particles from Monte Carlo (MC) simulation data.
 *
 * This function processes a list of input files containing histograms of particle data, calculates
 * the efficiency for different particle types (Xi and Omega) across multiple runs, and saves the
 * results to an output file. It also generates and saves various histograms and efficiency plots.
 *
 * @param inputFilename The name of the file containing the list of input files. Default is "McFileListPrefix.txt".
 * @param histName The name of the histogram to be analyzed. Default is "h3_ptmasscent_def".
 * @param analyseDiffs Boolean flag to indicate whether to analyze multiplicity differentiated results. Default is kTRUE.
 * @param outputFilename The name of the output file to save the results. Default is "TEST_effEst.root".
 * @param outputFolder The folder to save the output images. Default is "TEST_effEst_images".
 * @param saveStack Boolean flag to indicate whether to save the histogram stacks as images. Default is kTRUE.
 * @param imageFormat The format of the output images (e.g., "png"). Default is "png".
 * @param verbosity The verbosity level for error messages. Default is kInfo.
 *
 * @return Returns 0 on success, 1 on failure.
 */
int EfficiencyEstimation(
    std::string inputFilename = "McFileListPrefix.txt",
    TString histName = "h3_ptmasscent_def",
    Bool_t analyseDiffs = kTRUE,
    Bool_t analyseCommonMultDiffs = kTRUE,
    TString outputFilename = "/var/home/ishaan/Work/git/thesis/final/images/050625_EffEst/050625_EffEst_finalThesis.root",
    TString outputFolder = "/var/home/ishaan/Work/git/thesis/final/images/050625_EffEst",
    Bool_t saveStack = kTRUE,
    TString imageFormat = "pdf",
    Int_t verbosity = kInfo)
{
    ROOT::EnableImplicitMT();
    gStyle->SetOptFit(1111);
    gErrorIgnoreLevel = verbosity;
    TH1::SetDefaultSumw2(kTRUE);

    // SetCustomColorPalette();
    // extra options for decorating final plots (pdf) in thesis
    gStyle->SetLineScalePS(1.5);
    gStyle->SetStatFontSize(0.03);
    gStyle->SetPadTickX(1); // Ticks on both top and bottom for X axis
    gStyle->SetPadTickY(1); // Ticks on both left and right for Y axis

    // remove ownership of objects from file so we can delete the file ptr
    TH1::AddDirectory(kFALSE);

    std::vector<std::string> inputFileList = GetFileList(inputFilename);
    if (inputFileList.empty())
    {
        Error("EfficiencyEstimation: FileList", "'%s': Could not read the file! Aborting.", inputFilename.data());
        return 1;
    }

    outputFolder = SetOutputFolder(outputFolder);
    // gStyle->SetOptFit(1111);

    Int_t nMcRuns = inputFileList.size() - 1; //-1 because two runs are enriched with only one particle: Xi_LHC17e1b and Omega_LHC17e1a
    TH1 *resultParams_Xip_allInt[nMcRuns];
    TH1 *resultParams_Xim_allInt[nMcRuns];
    TH1 *resultParams_XiC_allInt[nMcRuns];
    TH1 *resultParams_Omp_allInt[nMcRuns];
    TH1 *resultParams_Omm_allInt[nMcRuns];
    TH1 *resultParams_OmC_allInt[nMcRuns];

    TH1 *resultParXip_pt[nMcRuns][fNptbins_Xi];
    TH1 *resultParXim_pt[nMcRuns][fNptbins_Xi];
    TH1 *resultParOmp_pt[nMcRuns][fNptbins_Om];
    TH1 *resultParOmm_pt[nMcRuns][fNptbins_Om];
    TH1 *resultParXiC_pt[nMcRuns][fNptbins_Xi];
    TH1 *resultParOmC_pt[nMcRuns][fNptbins_Om];

    TH1 *resultParXip_pt_mult[nMcRuns][fNptbins_Xi][fNmultbins_Xi];
    TH1 *resultParXim_pt_mult[nMcRuns][fNptbins_Xi][fNmultbins_Xi];
    TH1 *resultParOmp_pt_mult[nMcRuns][fNptbins_Om][fNmultbins_Om];
    TH1 *resultParOmm_pt_mult[nMcRuns][fNptbins_Om][fNmultbins_Om];
    TH1 *resultParXiC_pt_mult[nMcRuns][fNptbins_Xi][fNmultbins_Xi];
    TH1 *resultParOmC_pt_mult[nMcRuns][fNptbins_Om][fNmultbins_Om];

    TH1D *eff_xim_mult[nMcRuns][fNmultbins_Xi]; // mult binned efficiency
    TH1D *eff_xip_mult[nMcRuns][fNmultbins_Xi]; // mult binned efficiency
    TH1D *eff_omm_mult[nMcRuns][fNmultbins_Om]; // mult binned efficiency
    TH1D *eff_omp_mult[nMcRuns][fNmultbins_Om]; // mult binned efficiency
    TH1D *eff_xiC_mult[nMcRuns][fNmultbins_Xi]; // mult binned efficiency
    TH1D *eff_omC_mult[nMcRuns][fNmultbins_Om]; // mult binned efficiency

    TH1D *eff_xiC_ratio_mult[nMcRuns][fNmultbins_Xi]; // mult binned efficiency ratio
    TH1D *eff_omC_ratio_mult[nMcRuns][fNmultbins_Om]; // mult binned efficiency ratio

    TH1D *eff_xim[nMcRuns]; // mult integrated efficiency
    TH1D *eff_xip[nMcRuns]; // mult integrated efficiency
    TH1D *eff_omm[nMcRuns]; // mult integrated efficiency
    TH1D *eff_omp[nMcRuns]; // mult integrated efficiency
    TH1D *eff_xiC[nMcRuns]; // mult integrated efficiency
    TH1D *eff_omC[nMcRuns]; // mult integrated efficiency

    /// EXTRA: Added new multiplicity binning:
    const int nNewCommonMultBins = 3;
    const double newCommonMultEdges[nNewCommonMultBins + 1] = {0, 15, 60, 100};
    // --- Mapping for Xi Bins to New Common Bins ---
    // fMultbins_Xi indices: 0(0-5), 1(5-10), 2(10-15), 3(15-20), 4(20-30), 5(30-40), 6(40-50), 7(50-60), 8(60-80), 9(80-100)
    std::vector<std::vector<int>> originalIndicesMap_Xi(nNewCommonMultBins);
    originalIndicesMap_Xi[0] = {0, 1, 2};       // 0-15%
    originalIndicesMap_Xi[1] = {3, 4, 5, 6, 7}; // 15-60%
    originalIndicesMap_Xi[2] = {8, 9};          // 60-100%

    // --- Mapping for Om Bins to New Common Bins ---
    // fMultbins_Om indices: 0(0-5), 1(5-15), 2(15-30), 3(30-60), 4(60-100)
    std::vector<std::vector<int>> originalIndicesMap_Om(nNewCommonMultBins);
    originalIndicesMap_Om[0] = {0, 1}; // 0-15%
    originalIndicesMap_Om[1] = {2, 3}; // 15-60%
    originalIndicesMap_Om[2] = {4};    // 60-100%

    // const Int_t nMcRuns = 4;
    TH1D *eff_xiC_recalculated_commonMult[nMcRuns][nNewCommonMultBins]; // Use appropriate max runs
    TH1D *eff_omC_recalculated_commonMult[nMcRuns][nNewCommonMultBins]; // Use appropriate max runs

    /// Calculate avg efficiency:
    Double_t effXiAvg[fNptbins_Xi];
    Double_t effXiErrAvg[fNptbins_Xi];
    Double_t effOmAvg[fNptbins_Om];
    Double_t effOmErrAvg[fNptbins_Om];
    TH1D *eff_xiC_avg_ratio[nMcRuns];
    TH1D *eff_omC_avg_ratio[nMcRuns];

    /// Calculate avg efficiency for commonMult:
    Double_t effXiCAvgCommonMult[nNewCommonMultBins][fNptbins_Xi];
    Double_t effXiCAvgCommonMultErr[nNewCommonMultBins][fNptbins_Xi];
    Double_t effOmCAvgCommonMult[nNewCommonMultBins][fNptbins_Om];
    Double_t effOmCAvgCommonMultErr[nNewCommonMultBins][fNptbins_Om];
    TH1D *eff_xiC_avg_commonMult[nNewCommonMultBins];       // mult binned efficiency
    TH1D *eff_omC_avg_commonMult[nNewCommonMultBins];       // mult binned efficiency
    TH1D *eff_xiC_avg_ratio_commonMult[nNewCommonMultBins]; // mult binned efficiency ratio
    TH1D *eff_omC_avg_ratio_commonMult[nNewCommonMultBins]; // mult binned efficiency ratio

    TH1D *eff_xiC_avg = new TH1D("eff_xiC_avg", "#Xi^{#pm} MC Avg", fNptbins_Xi, fPtbins_Xi);
    TH1D *eff_omC_avg = new TH1D("eff_omC_avg", "#Omega^{#pm} MC Avg", fNptbins_Om, fPtbins_Om);

    TH1D *eff_xiC_avg_chi2 = new TH1D("eff_xiC_avg_chi2", "#Xi^{#pm} MC Avg: (#frac{#chi^{2}}{N-1})", fNptbins_Xi, fPtbins_Xi);
    TH1D *eff_omC_avg_chi2 = new TH1D("eff_omC_avg_chi2", "#Omega^{#pm} MC Avg: (#frac{#chi^{2}}{N-1})", fNptbins_Om, fPtbins_Om);

    // THstack efficiency
    auto hs_xip_eff = new THStack("hs_xip_eff", "#Xi^{+}");                                                                                   // mult integrated efficiency + mult binned efficiencies
    auto hs_xim_eff = new THStack("hs_xim_eff", "#Xi^{-}");                                                                                   // mult integrated efficiency + mult binned efficiencies
    auto hs_omp_eff = new THStack("hs_omp_eff", "#Omega^{+}");                                                                                // mult integrated efficiency + mult binned efficiencies
    auto hs_omm_eff = new THStack("hs_omm_eff", "#Omega^{-}");                                                                                // mult integrated efficiency + mult binned efficiencies
    auto hs_xiC_eff = new THStack("hs_xiC_eff", "#Xi^{#pm}");                                                                                 // mult integrated efficiency + mult binned efficiencies
    auto hs_omC_eff = new THStack("hs_omC_eff", "#Omega^{#pm}");                                                                              // mult integrated efficiency + mult binned efficiencies
    auto hs_xiC_effRatio_MCvsAvg = new THStack("hs_xiC_effRatio_MCvsAvg", "Efficiency Ratio #Xi^{+} + #Xi^{-}: MC vs Avg");                   // ratio of (mult integrated / mult binned) efficiencies
    auto hs_omC_effRatio_MCvsAvg = new THStack("hs_omC_effRatio_MCvsAvg", "Efficiency Ratio #Omega^{+} + #Omega^{-}: MC vs Avg");             // ratio of (mult integrated / mult binned) efficiencies
    auto hs_xiC_avg_commonMult = new THStack("hs_xiC_avg_commonMult", "Avg Efficiency #Xi^{+} + #Xi^{-}: Common Mult Eff and Avg Eff");       // mult integrated avg efficiency + mult binned avg efficiencies
    auto hs_omC_avg_commonMult = new THStack("hs_omC_avg_commonMult", "Avg Efficiency #Omega^{+} + #Omega^{-}: Common Mult Eff and Avg Eff"); // mult integrated avg efficiency + mult binned avg efficiencies

    THStack *hs_xiC_eff_ratio[nMcRuns];                                                                                                                                // ratio of (mult integrated / mult binned) efficiencies
    THStack *hs_omC_eff_ratio[nMcRuns];                                                                                                                                // ratio of (mult integrated / mult binned) efficiencies
    auto *hs_xiC_AvgEff_ratio_commonMult = new THStack("hs_xiC_AvgEff_ratio_commonMult", "Efficiency Ratio #Xi^{+} + #Xi^{-}: Avg Eff vs Common Mult Eff vs Avg Eff"); // ratio of (mult binned / mult integrated) efficiencies
    auto *hs_omC_AvgEff_ratio_commonMult = new THStack("hs_omC_AvgEff_ratio_commonMult", "Efficiency Ratio #Omega^{+} + #Omega^{-}: Common Mult Eff vs Avg Eff");      // ratio of (mult binned / mult integrated) efficiencies

    Int_t iRunXi = 0;
    Int_t iRunOm = 0;
    RunInfo XiRunInfo[nMcRuns];
    RunInfo OmegaRunInfo[nMcRuns];

    // iterate over the file list
    for (Int_t iFileList = 0; iFileList < (nMcRuns + 1); iFileList++)
    {
        // get the file name, run name
        TString fileName = inputFileList[iFileList];
        TString McRunName = fileName(fileName.Last('_') + 1, fileName.Length()); // get the last part of the file name - the run name (e.g. fileNamePrefix = "MCfile_LHC17e1b" -> McRunName = "LHC17e1b")
        fileName = fileName + "_" + histName + ".root";                          // append the histName to the file name prefix to get the full file name
        // fileName = fileName + ".root"; // append .root to file name prefix to get the full file name - only for SysUncertainty Signal extraction

        Info("EfficiencyEstimation: histInput", "Getting histograms for '%s': %s", McRunName.Data(), histName.Data());
        TFile *inputFile = OpenFile(fileName);
        if (!inputFile)
        {
            Error("EfficiencyEstimation: inputFile", "Cannot open file '%s' !", fileName.Data());
            return 1;
        }

        if (!McRunName.Contains("Om")) // get Xi hists if filename is NOT flagged for omega only
        {

            XiRunInfo[iRunXi].runName = ParseRunName(McRunName);
            XiRunInfo[iRunXi].effVal.resize(fNptbins_Xi);
            XiRunInfo[iRunXi].effWeight.resize(fNptbins_Xi);
            XiRunInfo[iRunXi].effValCommonMult.resize(nNewCommonMultBins);
            XiRunInfo[iRunXi].effWeightCommonMult.resize(nNewCommonMultBins);
            for (Int_t iNewBin = 0; iNewBin < nNewCommonMultBins; iNewBin++)
            {
                XiRunInfo[iRunXi].effValCommonMult[iNewBin].resize(fNptbins_Xi);
                XiRunInfo[iRunXi].effWeightCommonMult[iNewBin].resize(fNptbins_Xi);
            }

            resultParams_Xip_allInt[iRunXi] = (TH1 *)inputFile->FindObjectAny("resultParams_Xip_allInt");
            resultParams_Xim_allInt[iRunXi] = (TH1 *)inputFile->FindObjectAny("resultParams_Xim_allInt");
            resultParams_XiC_allInt[iRunXi] = (TH1 *)inputFile->FindObjectAny("resultParams_XiC_allInt");

            for (Int_t ptBinXi = 0; ptBinXi < fNptbins_Xi; ptBinXi++)
            {
                resultParXip_pt[iRunXi][ptBinXi] = (TH1 *)inputFile->FindObjectAny(TString::Format(("resultParXip_pt[%d]"), ptBinXi));
                resultParXim_pt[iRunXi][ptBinXi] = (TH1 *)inputFile->FindObjectAny(TString::Format(("resultParXim_pt[%d]"), ptBinXi));
                resultParXiC_pt[iRunXi][ptBinXi] = (TH1 *)inputFile->FindObjectAny(TString::Format(("resultParXiC_pt[%d]"), ptBinXi));

                if (analyseDiffs)
                {
                    for (Int_t multBinXi = 0; multBinXi < fNmultbins_Xi; multBinXi++)
                    {
                        resultParXip_pt_mult[iRunXi][ptBinXi][multBinXi] = (TH1 *)inputFile->FindObjectAny(TString::Format(("resultParXip_pt_mult[%d][%d]"), ptBinXi, multBinXi));
                        resultParXim_pt_mult[iRunXi][ptBinXi][multBinXi] = (TH1 *)inputFile->FindObjectAny(TString::Format(("resultParXim_pt_mult[%d][%d]"), ptBinXi, multBinXi));
                        resultParXiC_pt_mult[iRunXi][ptBinXi][multBinXi] = (TH1 *)inputFile->FindObjectAny(TString::Format(("resultParXiC_pt_mult[%d][%d]"), ptBinXi, multBinXi));
                    }
                }
            }
            iRunXi++;
        }

        if (!McRunName.Contains("Xi")) // get Omega hists if filename is NOT flagged for xi only
        {
            OmegaRunInfo[iRunOm].runName = ParseRunName(McRunName);
            OmegaRunInfo[iRunOm].effVal.resize(fNptbins_Om);
            OmegaRunInfo[iRunOm].effWeight.resize(fNptbins_Om);
            OmegaRunInfo[iRunOm].effValCommonMult.resize(nNewCommonMultBins);
            OmegaRunInfo[iRunOm].effWeightCommonMult.resize(nNewCommonMultBins);
            for (Int_t iNewBin = 0; iNewBin < nNewCommonMultBins; iNewBin++)
            {
                OmegaRunInfo[iRunOm].effValCommonMult[iNewBin].resize(fNptbins_Om);
                OmegaRunInfo[iRunOm].effWeightCommonMult[iNewBin].resize(fNptbins_Om);
            }

            resultParams_Omp_allInt[iRunOm] = (TH1 *)inputFile->FindObjectAny("resultParams_Omp_allInt");
            resultParams_Omm_allInt[iRunOm] = (TH1 *)inputFile->FindObjectAny("resultParams_Omm_allInt");
            resultParams_OmC_allInt[iRunOm] = (TH1 *)inputFile->FindObjectAny("resultParams_OmC_allInt");

            for (Int_t ptBinOm = 0; ptBinOm < fNptbins_Om; ptBinOm++)
            {
                resultParOmp_pt[iRunOm][ptBinOm] = (TH1 *)inputFile->FindObjectAny(TString::Format(("resultParOmp_pt[%d]"), ptBinOm));
                resultParOmm_pt[iRunOm][ptBinOm] = (TH1 *)inputFile->FindObjectAny(TString::Format(("resultParOmm_pt[%d]"), ptBinOm));
                resultParOmC_pt[iRunOm][ptBinOm] = (TH1 *)inputFile->FindObjectAny(TString::Format(("resultParOmC_pt[%d]"), ptBinOm));

                if (analyseDiffs)
                {
                    for (Int_t multBinOm = 0; multBinOm < fNmultbins_Om; multBinOm++)
                    {
                        resultParOmp_pt_mult[iRunOm][ptBinOm][multBinOm] = (TH1 *)inputFile->FindObjectAny(TString::Format(("resultParOmp_pt_mult[%d][%d]"), ptBinOm, multBinOm));
                        resultParOmm_pt_mult[iRunOm][ptBinOm][multBinOm] = (TH1 *)inputFile->FindObjectAny(TString::Format(("resultParOmm_pt_mult[%d][%d]"), ptBinOm, multBinOm));
                        resultParOmC_pt_mult[iRunOm][ptBinOm][multBinOm] = (TH1 *)inputFile->FindObjectAny(TString::Format(("resultParOmC_pt_mult[%d][%d]"), ptBinOm, multBinOm));
                    }
                }
            }
            iRunOm++;
        }
        Info("EfficiencyEstimation: histInput", "Input ended.");
        inputFile->Close();
        delete inputFile;
    } /// Input ended!

    /// Saving output :
    if (outputFilename.IsNull())
    {
        outputFilename = "Efficiency_" + histName + ".root";
    }

    TFile *outputFile = OpenFile(outputFilename, "RECREATE");

    outputFile->mkdir("dirHistStacks");
    outputFile->mkdir("dirEffPt_xim");
    outputFile->mkdir("dirEffPt_xip");
    outputFile->mkdir("dirEffPt_omm");
    outputFile->mkdir("dirEffPt_omp");
    outputFile->mkdir("dirEffPt_xiC");
    outputFile->mkdir("dirEffPt_omC");
    if (analyseDiffs)
    {
        outputFile->mkdir("dirHistStackRatios");
        outputFile->mkdir("dirEffPtMult_xim");
        outputFile->mkdir("dirEffPtMult_xip");
        outputFile->mkdir("dirEffPtMult_omm");
        outputFile->mkdir("dirEffPtMult_omp");
        outputFile->mkdir("dirEffPtMult_xiC");
        outputFile->mkdir("dirEffPtMult_omC");
        outputFile->mkdir("dirEffPtMultRatio_xiC");
        outputFile->mkdir("dirEffPtMultRatio_omC");
    }
    /// Output set.

    Info("EfficiencyEstimation: outputFile", "Output file created: %s", outputFilename.Data());

    /// Begin:

    for (Int_t iRun = 0; iRun < nMcRuns; iRun++)
    {

        /// CASE: Xi
        Info("EfficiencyEstimation", "%d: %s Plotting Xi efficiency ...", iRun, XiRunInfo[iRun].runName.Data());
        /// generate mult integrated efficiency histograms
        eff_xim[iRun] = new TH1D(TString::Format(("eff_xim[%d]"), iRun), TString::Format("%s", XiRunInfo[iRun].runName.Data()), fNptbins_Xi, fPtbins_Xi);
        eff_xip[iRun] = new TH1D(TString::Format(("eff_xip[%d]"), iRun), TString::Format("%s", XiRunInfo[iRun].runName.Data()), fNptbins_Xi, fPtbins_Xi);
        eff_xiC[iRun] = new TH1D(TString::Format(("eff_xiC[%d]"), iRun), TString::Format("%s", XiRunInfo[iRun].runName.Data()), fNptbins_Xi, fPtbins_Xi);

        for (Int_t ptBinXi = 0; ptBinXi < fNptbins_Xi; ptBinXi++)
        {
            /// add XiC mult int. efficiency values and weights to the RunInfo struct
            XiRunInfo[iRun].effVal[ptBinXi] = resultParXiC_pt[iRun][ptBinXi]->GetBinContent(7);
            XiRunInfo[iRun].effWeight[ptBinXi] = 1.0 / (pow(resultParXiC_pt[iRun][ptBinXi]->GetBinError(7), 2));

            /// fill mult integrated efficiency histograms
            eff_xim[iRun]->SetBinContent(ptBinXi + 1, resultParXim_pt[iRun][ptBinXi]->GetBinContent(7));
            eff_xip[iRun]->SetBinContent(ptBinXi + 1, resultParXip_pt[iRun][ptBinXi]->GetBinContent(7));
            eff_xiC[iRun]->SetBinContent(ptBinXi + 1, resultParXiC_pt[iRun][ptBinXi]->GetBinContent(7));

            eff_xim[iRun]->SetBinError(ptBinXi + 1, resultParXim_pt[iRun][ptBinXi]->GetBinError(7));
            eff_xip[iRun]->SetBinError(ptBinXi + 1, resultParXip_pt[iRun][ptBinXi]->GetBinError(7));
            eff_xiC[iRun]->SetBinError(ptBinXi + 1, resultParXiC_pt[iRun][ptBinXi]->GetBinError(7));
        }

        outputFile->cd("dirEffPt_xim");
        eff_xim[iRun]->SetMarkerStyle(openMarkerStyles[iRun]);
        eff_xim[iRun]->Write();
        eff_xim[iRun]->SetName(TString::Format(("%s: Mult: %.0f-%.0f%%"), XiRunInfo[iRun].runName.Data(), fMultbins_Xi[0], fMultbins_Xi[fNmultbins_Xi]));
        hs_xim_eff->Add(eff_xim[iRun]);

        outputFile->cd("dirEffPt_xip");
        eff_xip[iRun]->SetMarkerStyle(openMarkerStyles[iRun]);
        eff_xip[iRun]->Write();
        eff_xip[iRun]->SetName(TString::Format(("%s: Mult: %.0f-%.0f%%"), XiRunInfo[iRun].runName.Data(), fMultbins_Xi[0], fMultbins_Xi[fNmultbins_Xi]));
        hs_xip_eff->Add(eff_xip[iRun]);

        outputFile->cd("dirEffPt_xiC");
        eff_xiC[iRun]->SetMarkerStyle(openMarkerStyles[iRun]);
        eff_xiC[iRun]->Write();
        eff_xiC[iRun]->SetName(TString::Format(("%s: Mult: %.0f-%.0f%%"), XiRunInfo[iRun].runName.Data(), fMultbins_Xi[0], fMultbins_Xi[fNmultbins_Xi]));
        hs_xiC_eff->Add(eff_xiC[iRun]);

        /// mult binned efficiency
        if (analyseDiffs)
        {
            hs_xiC_eff_ratio[iRun] = new THStack(TString::Format("hs_xiC_eff_ratio[%d]", iRun), TString::Format("%s Efficiency #Xi^{+} + #Xi^{-}: Mult_Int/Mult_Diff", XiRunInfo[iRun].runName.Data()));

            for (Int_t multBinXi = 0; multBinXi < fNmultbins_Xi; multBinXi++)
            {
                eff_xim_mult[iRun][multBinXi] = new TH1D(TString::Format(("eff_xim_mult[%d][%d]"), iRun, multBinXi), TString::Format("%s Efficiency: Mult %.0f-%.0f%%", XiRunInfo[iRun].runName.Data(), fMultbins_Xi[multBinXi], fMultbins_Xi[multBinXi + 1]), fNptbins_Xi, fPtbins_Xi);
                eff_xip_mult[iRun][multBinXi] = new TH1D(TString::Format(("eff_xip_mult[%d][%d]"), iRun, multBinXi), TString::Format("%s Efficiency: Mult %.0f-%.0f%%", XiRunInfo[iRun].runName.Data(), fMultbins_Xi[multBinXi], fMultbins_Xi[multBinXi + 1]), fNptbins_Xi, fPtbins_Xi);
                eff_xiC_mult[iRun][multBinXi] = new TH1D(TString::Format(("eff_xiC_mult[%d][%d]"), iRun, multBinXi), TString::Format("%s Efficiency: Mult %.0f-%.0f%%", XiRunInfo[iRun].runName.Data(), fMultbins_Xi[multBinXi], fMultbins_Xi[multBinXi + 1]), fNptbins_Xi, fPtbins_Xi);

                for (Int_t ptBinXi = 0; ptBinXi < fNptbins_Xi; ptBinXi++)
                {
                    eff_xim_mult[iRun][multBinXi]->SetBinContent(ptBinXi + 1, resultParXim_pt_mult[iRun][ptBinXi][multBinXi]->GetBinContent(7));
                    eff_xip_mult[iRun][multBinXi]->SetBinContent(ptBinXi + 1, resultParXip_pt_mult[iRun][ptBinXi][multBinXi]->GetBinContent(7));
                    eff_xiC_mult[iRun][multBinXi]->SetBinContent(ptBinXi + 1, resultParXiC_pt_mult[iRun][ptBinXi][multBinXi]->GetBinContent(7));

                    eff_xim_mult[iRun][multBinXi]->SetBinError(ptBinXi + 1, resultParXim_pt_mult[iRun][ptBinXi][multBinXi]->GetBinError(7));
                    eff_xip_mult[iRun][multBinXi]->SetBinError(ptBinXi + 1, resultParXip_pt_mult[iRun][ptBinXi][multBinXi]->GetBinError(7));
                    eff_xiC_mult[iRun][multBinXi]->SetBinError(ptBinXi + 1, resultParXiC_pt_mult[iRun][ptBinXi][multBinXi]->GetBinError(7));
                }

                /// calculate ratio of mult integrated and mult binned efficiencies
                eff_xiC_ratio_mult[iRun][multBinXi] = (TH1D *)eff_xiC[iRun]->Clone(TString::Format(("eff_xiC_ratio_mult[%d]"), multBinXi));
                eff_xiC_ratio_mult[iRun][multBinXi]->SetTitle(TString::Format("%s Efficiency: Mult #frac{0-100%%}{%.0f-%.0f%%}", XiRunInfo[iRun].runName.Data(), fMultbins_Xi[multBinXi], fMultbins_Xi[multBinXi + 1]));
                eff_xiC_ratio_mult[iRun][multBinXi]->Divide(eff_xiC_mult[iRun][multBinXi]);

                // save the histograms and add them to the stack
                outputFile->cd("dirEffPtMult_xim");
                eff_xim_mult[iRun][multBinXi]->SetMarkerStyle(markerStyles[iRun]);
                eff_xim_mult[iRun][multBinXi]->Write();
                eff_xim_mult[iRun][multBinXi]->SetName(TString::Format(("%s: Mult: %.0f-%.0f%%"), XiRunInfo[iRun].runName.Data(), fMultbins_Xi[multBinXi], fMultbins_Xi[multBinXi + 1]));
                // hs_xim_eff->Add(eff_xim_mult[iRun][multBinXi]);

                outputFile->cd("dirEffPtMult_xip");
                eff_xip_mult[iRun][multBinXi]->SetMarkerStyle(markerStyles[iRun]);
                eff_xip_mult[iRun][multBinXi]->Write();
                eff_xip_mult[iRun][multBinXi]->SetName(TString::Format(("%s: Mult: %.0f-%.0f%%"), XiRunInfo[iRun].runName.Data(), fMultbins_Xi[multBinXi], fMultbins_Xi[multBinXi + 1]));
                // hs_xip_eff->Add(eff_xip_mult[iRun][multBinXi]);

                outputFile->cd("dirEffPtMult_xiC");
                eff_xiC_mult[iRun][multBinXi]->SetMarkerStyle(markerStyles[iRun]);
                eff_xiC_mult[iRun][multBinXi]->Write();
                eff_xiC_mult[iRun][multBinXi]->SetName(TString::Format(("%s: Mult: %.0f-%.0f%%"), XiRunInfo[iRun].runName.Data(), fMultbins_Xi[multBinXi], fMultbins_Xi[multBinXi + 1]));
                // hs_xiC_eff->Add(eff_xiC_mult[iRun][multBinXi]);

                outputFile->cd("dirEffPtMultRatio_xiC");
                eff_xiC_ratio_mult[iRun][multBinXi]->SetMarkerStyle(markerStyles[iRun]);
                eff_xiC_ratio_mult[iRun][multBinXi]->Write();
                eff_xiC_ratio_mult[iRun][multBinXi]->SetName(TString::Format(("%s: Mult: %.0f-%.0f%%"), XiRunInfo[iRun].runName.Data(), fMultbins_Xi[multBinXi], fMultbins_Xi[multBinXi + 1]));
                hs_xiC_eff_ratio[iRun]->Add(eff_xiC_ratio_mult[iRun][multBinXi]);
            }
        }

        /// CASE: Omega
        Info("EfficiencyEstimation", "%d: %s Plotting Omega efficiency ...", iRun, OmegaRunInfo[iRun].runName.Data());
        /// generate mult integrated efficiency histograms
        eff_omm[iRun] = new TH1D(TString::Format(("eff_omm[%d]"), iRun), TString::Format("%s", OmegaRunInfo[iRun].runName.Data()), fNptbins_Om, fPtbins_Om);
        eff_omp[iRun] = new TH1D(TString::Format(("eff_omp[%d]"), iRun), TString::Format("%s", OmegaRunInfo[iRun].runName.Data()), fNptbins_Om, fPtbins_Om);
        eff_omC[iRun] = new TH1D(TString::Format(("eff_omC[%d]"), iRun), TString::Format("%s", OmegaRunInfo[iRun].runName.Data()), fNptbins_Om, fPtbins_Om);

        for (Int_t ptBinOm = 0; ptBinOm < fNptbins_Om; ptBinOm++)
        {
            /// add OmC mult int. efficiency values and weights to the RunInfo struct
            OmegaRunInfo[iRun].effVal[ptBinOm] = resultParOmC_pt[iRun][ptBinOm]->GetBinContent(7);
            OmegaRunInfo[iRun].effWeight[ptBinOm] = 1.0 / (pow(resultParOmC_pt[iRun][ptBinOm]->GetBinError(7), 2));

            /// set the mult integrated efficiency values and errors
            eff_omm[iRun]->SetBinContent(ptBinOm + 1, resultParOmm_pt[iRun][ptBinOm]->GetBinContent(7));
            eff_omp[iRun]->SetBinContent(ptBinOm + 1, resultParOmp_pt[iRun][ptBinOm]->GetBinContent(7));
            eff_omC[iRun]->SetBinContent(ptBinOm + 1, resultParOmC_pt[iRun][ptBinOm]->GetBinContent(7));

            eff_omm[iRun]->SetBinError(ptBinOm + 1, resultParOmm_pt[iRun][ptBinOm]->GetBinError(7));
            eff_omp[iRun]->SetBinError(ptBinOm + 1, resultParOmp_pt[iRun][ptBinOm]->GetBinError(7));
            eff_omC[iRun]->SetBinError(ptBinOm + 1, resultParOmC_pt[iRun][ptBinOm]->GetBinError(7));
        }

        outputFile->cd("dirEffPt_omm");
        eff_omm[iRun]->SetMarkerStyle(openMarkerStyles[iRun]);
        eff_omm[iRun]->Write();
        eff_omm[iRun]->SetName(TString::Format(("%s: Mult: %.0f-%.0f%%"), OmegaRunInfo[iRun].runName.Data(), fMultbins_Om[0], fMultbins_Om[fNmultbins_Om]));
        hs_omm_eff->Add(eff_omm[iRun]);

        outputFile->cd("dirEffPt_omp");
        eff_omp[iRun]->SetMarkerStyle(openMarkerStyles[iRun]);
        eff_omp[iRun]->Write();
        eff_omp[iRun]->SetName(TString::Format(("%s: Mult: %.0f-%.0f%%"), OmegaRunInfo[iRun].runName.Data(), fMultbins_Om[0], fMultbins_Om[fNmultbins_Om]));
        hs_omp_eff->Add(eff_omp[iRun]);

        outputFile->cd("dirEffPt_omC");
        eff_omC[iRun]->SetMarkerStyle(openMarkerStyles[iRun]);
        eff_omC[iRun]->Write();
        eff_omC[iRun]->SetName(TString::Format(("%s: Mult: %.0f-%.0f%%"), OmegaRunInfo[iRun].runName.Data(), fMultbins_Om[0], fMultbins_Om[fNmultbins_Om]));
        hs_omC_eff->Add(eff_omC[iRun]);

        /// mult binned efficiency
        if (analyseDiffs)
        {
            hs_omC_eff_ratio[iRun] = new THStack(TString::Format("hs_omC_eff_ratio[%d]", iRun), TString::Format("%s Efficiency #Omega^{+} + #Omega^{-}: Mult_Int/Mult_Diff", OmegaRunInfo[iRun].runName.Data()));

            for (Int_t multBinOm = 0; multBinOm < fNmultbins_Om; multBinOm++)
            {
                eff_omm_mult[iRun][multBinOm] = new TH1D(TString::Format(("eff_omm_mult[%d][%d]"), iRun, multBinOm), TString::Format("%s Efficiency: Mult %.0f-%.0f%%", OmegaRunInfo[iRun].runName.Data(), fMultbins_Om[multBinOm], fMultbins_Om[multBinOm + 1]), fNptbins_Om, fPtbins_Om);
                eff_omp_mult[iRun][multBinOm] = new TH1D(TString::Format(("eff_omp_mult[%d][%d]"), iRun, multBinOm), TString::Format("%s Efficiency: Mult %.0f-%.0f%%", OmegaRunInfo[iRun].runName.Data(), fMultbins_Om[multBinOm], fMultbins_Om[multBinOm + 1]), fNptbins_Om, fPtbins_Om);
                eff_omC_mult[iRun][multBinOm] = new TH1D(TString::Format(("eff_omC_mult[%d][%d]"), iRun, multBinOm), TString::Format("%s Efficiency: Mult %.0f-%.0f%%", OmegaRunInfo[iRun].runName.Data(), fMultbins_Om[multBinOm], fMultbins_Om[multBinOm + 1]), fNptbins_Om, fPtbins_Om);

                for (Int_t ptBinOm = 0; ptBinOm < fNptbins_Om; ptBinOm++)
                {
                    eff_omm_mult[iRun][multBinOm]->SetBinContent(ptBinOm + 1, resultParOmm_pt_mult[iRun][ptBinOm][multBinOm]->GetBinContent(7));
                    eff_omp_mult[iRun][multBinOm]->SetBinContent(ptBinOm + 1, resultParOmp_pt_mult[iRun][ptBinOm][multBinOm]->GetBinContent(7));
                    eff_omC_mult[iRun][multBinOm]->SetBinContent(ptBinOm + 1, resultParOmC_pt_mult[iRun][ptBinOm][multBinOm]->GetBinContent(7));

                    eff_omm_mult[iRun][multBinOm]->SetBinError(ptBinOm + 1, resultParOmm_pt_mult[iRun][ptBinOm][multBinOm]->GetBinError(7));
                    eff_omp_mult[iRun][multBinOm]->SetBinError(ptBinOm + 1, resultParOmp_pt_mult[iRun][ptBinOm][multBinOm]->GetBinError(7));
                    eff_omC_mult[iRun][multBinOm]->SetBinError(ptBinOm + 1, resultParOmC_pt_mult[iRun][ptBinOm][multBinOm]->GetBinError(7));
                }

                /// calculate ratio of mult integrated and mult binned efficiencies
                eff_omC_ratio_mult[iRun][multBinOm] = (TH1D *)eff_omC[iRun]->Clone(TString::Format(("eff_omC_ratio_mult[%d]"), multBinOm));
                eff_omC_ratio_mult[iRun][multBinOm]->SetTitle(TString::Format("%s Efficiency: Mult #frac{0-100%%}{%.0f-%.0f%%}", OmegaRunInfo[iRun].runName.Data(), fMultbins_Om[multBinOm], fMultbins_Om[multBinOm + 1]));
                eff_omC_ratio_mult[iRun][multBinOm]->Divide(eff_omC_mult[iRun][multBinOm]);

                // save the histograms and add them to the stack
                outputFile->cd("dirEffPtMult_omm");
                eff_omm_mult[iRun][multBinOm]->SetMarkerStyle(markerStyles[iRun]);
                eff_omm_mult[iRun][multBinOm]->Write();
                eff_omm_mult[iRun][multBinOm]->SetName(TString::Format(("%s: Mult: %.0f-%.0f%%"), OmegaRunInfo[iRun].runName.Data(), fMultbins_Om[multBinOm], fMultbins_Om[multBinOm + 1]));
                // hs_omm_eff->Add(eff_omm_mult[iRun][multBinOm]);

                outputFile->cd("dirEffPtMult_omp");
                eff_omp_mult[iRun][multBinOm]->SetMarkerStyle(markerStyles[iRun]);
                eff_omp_mult[iRun][multBinOm]->Write();
                eff_omp_mult[iRun][multBinOm]->SetName(TString::Format(("%s: Mult: %.0f-%.0f%%"), OmegaRunInfo[iRun].runName.Data(), fMultbins_Om[multBinOm], fMultbins_Om[multBinOm + 1]));
                // hs_omp_eff->Add(eff_omp_mult[iRun][multBinOm]);

                outputFile->cd("dirEffPtMult_omC");
                eff_omC_mult[iRun][multBinOm]->SetMarkerStyle(markerStyles[iRun]);
                eff_omC_mult[iRun][multBinOm]->Write();
                eff_omC_mult[iRun][multBinOm]->SetName(TString::Format(("%s: Mult: %.0f-%.0f%%"), OmegaRunInfo[iRun].runName.Data(), fMultbins_Om[multBinOm], fMultbins_Om[multBinOm + 1]));
                // hs_omC_eff->Add(eff_omC_mult[iRun][multBinOm]);

                outputFile->cd("dirEffPtMultRatio_omC");
                eff_omC_ratio_mult[iRun][multBinOm]->SetMarkerStyle(markerStyles[iRun]);
                eff_omC_ratio_mult[iRun][multBinOm]->Write();
                eff_omC_ratio_mult[iRun][multBinOm]->SetName(TString::Format(("%s: Mult: %.0f-%.0f%%"), OmegaRunInfo[iRun].runName.Data(), fMultbins_Om[multBinOm], fMultbins_Om[multBinOm + 1]));
                hs_omC_eff_ratio[iRun]->Add(eff_omC_ratio_mult[iRun][multBinOm]);
            }
        }
    }

    // --- Processing new mult bins for Xi Particles ---
    for (Int_t iRun = 0; iRun < nMcRuns; iRun++)
    {
        for (Int_t iNewBin = 0; iNewBin < nNewCommonMultBins; iNewBin++)
        {
            eff_xiC_recalculated_commonMult[iRun][iNewBin] = new TH1D(
                TString::Format("eff_xiC_recalc_commonMult_run%d_newmult_%.0f_%.0f", iRun, newCommonMultEdges[iNewBin], newCommonMultEdges[iNewBin + 1]),
                TString::Format("%s Xi Recalc. Eff: Mult %.0f-%.0f%%", XiRunInfo[iRun].runName.Data(), newCommonMultEdges[iNewBin], newCommonMultEdges[iNewBin + 1]),
                fNptbins_Xi, fPtbins_Xi);
            eff_xiC_recalculated_commonMult[iRun][iNewBin]->SetXTitle("p_{T} (GeV/c)");
            eff_xiC_recalculated_commonMult[iRun][iNewBin]->SetYTitle("Efficiency");

            for (Int_t ptBin = 0; ptBin < fNptbins_Xi; ptBin++)
            {
                double sum_N_reco_xi = 0.0;
                double sum_N_gen_xi = 0.0;
                double sum_N_reco_err_sq_xi = 0.0;
                double sum_N_gen_err_sq_xi = 0.0;

                for (int orig_xi_mult_idx : originalIndicesMap_Xi[iNewBin])
                {
                    if (orig_xi_mult_idx >= fNmultbins_Xi)
                        continue;
                    TH1 *sourceParamHistXi = resultParXiC_pt_mult[iRun][ptBin][orig_xi_mult_idx];
                    if (sourceParamHistXi != nullptr)
                    {
                        sum_N_reco_xi += sourceParamHistXi->GetBinContent(1);
                        sum_N_reco_err_sq_xi += std::pow(sourceParamHistXi->GetBinError(1), 2);
                        sum_N_gen_xi += sourceParamHistXi->GetBinContent(6);
                        sum_N_gen_err_sq_xi += std::pow(sourceParamHistXi->GetBinError(6), 2);
                    }
                }

                double final_eff_xi = 0.0;
                double final_eff_err_xi = 0.0;
                if (sum_N_gen_xi > 0)
                {
                    final_eff_xi = sum_N_reco_xi / sum_N_gen_xi;
                    double final_err_N_reco_xi = std::sqrt(sum_N_reco_err_sq_xi);
                    double final_err_N_gen_xi = std::sqrt(sum_N_gen_err_sq_xi);
                    final_eff_err_xi = ErrorInRatio(sum_N_reco_xi, final_err_N_reco_xi, sum_N_gen_xi, final_err_N_gen_xi);
                }
                eff_xiC_recalculated_commonMult[iRun][iNewBin]->SetBinContent(ptBin + 1, final_eff_xi);
                eff_xiC_recalculated_commonMult[iRun][iNewBin]->SetBinError(ptBin + 1, final_eff_err_xi);
                /// add XiC mult binned efficiency values and weights to the RunInfo struct
                XiRunInfo[iRun].effValCommonMult[iNewBin][ptBin] = final_eff_xi;
                if (final_eff_err_xi > 0)
                {
                    XiRunInfo[iRun].effWeightCommonMult[iNewBin][ptBin] = 1.0 / (pow(final_eff_err_xi, 2));
                }
                else
                {
                    XiRunInfo[iRun].effWeightCommonMult[iNewBin][ptBin] = 0.0; // Avoid division by zero
                }
            }
        }
    }

    // --- Processing new mult bins for Omega Particles ---
    for (Int_t iRun = 0; iRun < nMcRuns; iRun++)
    {
        for (Int_t iNewBin = 0; iNewBin < nNewCommonMultBins; iNewBin++)
        {
            eff_omC_recalculated_commonMult[iRun][iNewBin] = new TH1D(
                TString::Format("eff_omC_recalc_commonMult_run%d_newmult_%.0f_%.0f", iRun, newCommonMultEdges[iNewBin], newCommonMultEdges[iNewBin + 1]),
                TString::Format("%s Omega Recalc. Eff: Mult %.0f-%.0f%%", OmegaRunInfo[iRun].runName.Data(), newCommonMultEdges[iNewBin], newCommonMultEdges[iNewBin + 1]),
                fNptbins_Om, fPtbins_Om);
            eff_omC_recalculated_commonMult[iRun][iNewBin]->SetXTitle("p_{T} (GeV/c)");
            eff_omC_recalculated_commonMult[iRun][iNewBin]->SetYTitle("Efficiency");

            for (Int_t ptBin = 0; ptBin < fNptbins_Om; ptBin++)
            {
                double sum_N_reco_om = 0.0;
                double sum_N_gen_om = 0.0;
                double sum_N_reco_err_sq_om = 0.0;
                double sum_N_gen_err_sq_om = 0.0;

                for (int orig_om_mult_idx : originalIndicesMap_Om[iNewBin])
                {
                    if (orig_om_mult_idx >= fNmultbins_Om)
                        continue;
                    TH1 *sourceParamHistOm = resultParOmC_pt_mult[iRun][ptBin][orig_om_mult_idx];
                    if (sourceParamHistOm != nullptr)
                    {
                        sum_N_reco_om += sourceParamHistOm->GetBinContent(1);
                        sum_N_reco_err_sq_om += std::pow(sourceParamHistOm->GetBinError(1), 2);
                        sum_N_gen_om += sourceParamHistOm->GetBinContent(6);
                        sum_N_gen_err_sq_om += std::pow(sourceParamHistOm->GetBinError(6), 2);
                    }
                }

                double final_eff_om = 0.0;
                double final_eff_err_om = 0.0;
                if (sum_N_gen_om > 0)
                {
                    final_eff_om = sum_N_reco_om / sum_N_gen_om;
                    double final_err_N_reco_om = std::sqrt(sum_N_reco_err_sq_om);
                    double final_err_N_gen_om = std::sqrt(sum_N_gen_err_sq_om);
                    final_eff_err_om = ErrorInRatio(sum_N_reco_om, final_err_N_reco_om, sum_N_gen_om, final_err_N_gen_om);
                }
                eff_omC_recalculated_commonMult[iRun][iNewBin]->SetBinContent(ptBin + 1, final_eff_om);
                eff_omC_recalculated_commonMult[iRun][iNewBin]->SetBinError(ptBin + 1, final_eff_err_om);

                /// add OmC mult binned efficiency values and weights to the RunInfo struct
                OmegaRunInfo[iRun].effValCommonMult[iNewBin][ptBin] = final_eff_om;
                if (final_eff_err_om > 0)
                {
                    OmegaRunInfo[iRun].effWeightCommonMult[iNewBin][ptBin] = 1.0 / (pow(final_eff_err_om, 2));
                }
                else
                {
                    OmegaRunInfo[iRun].effWeightCommonMult[iNewBin][ptBin] = 0.0; // Avoid division by zero
                }
            }
        }
    }
    /// Calculate avg mult integrated efficiency: Xi
    for (Int_t ptBinXi = 0; ptBinXi < fNptbins_Xi; ptBinXi++)
    {
        Double_t sumEffWeights = 0;
        Double_t sumEffWeightVals = 0;

        for (Int_t iRun = 0; iRun < nMcRuns; iRun++)
        {
            sumEffWeights = sumEffWeights + XiRunInfo[iRun].effWeight[ptBinXi];
            sumEffWeightVals = sumEffWeightVals + (XiRunInfo[iRun].effWeight[ptBinXi] * XiRunInfo[iRun].effVal[ptBinXi]);
        }

        effXiAvg[ptBinXi] = sumEffWeightVals / sumEffWeights;
        effXiErrAvg[ptBinXi] = pow(sumEffWeights, -0.5);

        eff_xiC_avg->SetBinContent(ptBinXi + 1, effXiAvg[ptBinXi]);
        eff_xiC_avg->SetBinError(ptBinXi + 1, effXiErrAvg[ptBinXi]);
    }

    /// Calculate avg mult integrated efficiency: Omega
    for (Int_t ptBinOm = 0; ptBinOm < fNptbins_Om; ptBinOm++)
    {
        Double_t sumEffWeights = 0;
        Double_t sumEffWeightVals = 0;

        for (Int_t iRun = 0; iRun < nMcRuns; iRun++)
        {
            sumEffWeights = sumEffWeights + OmegaRunInfo[iRun].effWeight[ptBinOm];
            sumEffWeightVals = sumEffWeightVals + (OmegaRunInfo[iRun].effWeight[ptBinOm] * OmegaRunInfo[iRun].effVal[ptBinOm]);
        }

        effOmAvg[ptBinOm] = sumEffWeightVals / sumEffWeights;
        effOmErrAvg[ptBinOm] = pow(sumEffWeights, -0.5);

        eff_omC_avg->SetBinContent(ptBinOm + 1, effOmAvg[ptBinOm]);
        eff_omC_avg->SetBinError(ptBinOm + 1, effOmErrAvg[ptBinOm]);
    }

    // Calculate average commonMult binned efficiency: Xi
    for (Int_t iNewBin = 0; iNewBin < nNewCommonMultBins; iNewBin++)
    {
        eff_xiC_avg_commonMult[iNewBin] = new TH1D(
            TString::Format("eff_xiC_avg_commonMult[%d]", iNewBin),
            TString::Format("V0A %.0f-%.0f%%", newCommonMultEdges[iNewBin], newCommonMultEdges[iNewBin + 1]),
            fNptbins_Xi, fPtbins_Xi);

        for (Int_t ptBinXi = 0; ptBinXi < fNptbins_Xi; ptBinXi++)
        {
            Double_t sumEffWeights = 0;
            Double_t sumEffWeightVals = 0;

            for (Int_t iRun = 0; iRun < nMcRuns; iRun++)
            {
                sumEffWeights = sumEffWeights + XiRunInfo[iRun].effWeightCommonMult[iNewBin][ptBinXi];
                sumEffWeightVals = sumEffWeightVals + (XiRunInfo[iRun].effWeightCommonMult[iNewBin][ptBinXi] * XiRunInfo[iRun].effValCommonMult[iNewBin][ptBinXi]);
            }

            effXiCAvgCommonMult[iNewBin][ptBinXi] = sumEffWeightVals / sumEffWeights;
            effXiCAvgCommonMultErr[iNewBin][ptBinXi] = pow(sumEffWeights, -0.5);
            eff_xiC_avg_commonMult[iNewBin]->SetBinContent(ptBinXi + 1, effXiCAvgCommonMult[iNewBin][ptBinXi]);
            eff_xiC_avg_commonMult[iNewBin]->SetBinError(ptBinXi + 1, effXiCAvgCommonMultErr[iNewBin][ptBinXi]);
        }
        outputFile->cd("dirEffPtMultRatio_xiC");
        eff_xiC_avg_commonMult[iNewBin]->SetMarkerStyle(openMarkerStyles[iNewBin]);
        eff_xiC_avg_commonMult[iNewBin]->Write();
        hs_xiC_avg_commonMult->Add(eff_xiC_avg_commonMult[iNewBin]);
    }
    // eff_xiC_avg->SetTitle("V0A 0-100%");
    hs_xiC_avg_commonMult->Add(eff_xiC_avg);

    /// Calculate average commonMult binned efficiency: Omega
    for (Int_t iNewBin = 0; iNewBin < nNewCommonMultBins; iNewBin++)
    {
        eff_omC_avg_commonMult[iNewBin] = new TH1D(
            TString::Format("eff_omC_avg_commonMult[%d]", iNewBin),
            TString::Format("V0A %.0f-%.0f%%", newCommonMultEdges[iNewBin], newCommonMultEdges[iNewBin + 1]),
            fNptbins_Om, fPtbins_Om);
        for (Int_t ptBinOm = 0; ptBinOm < fNptbins_Om; ptBinOm++)
        {
            Double_t sumEffWeights = 0;
            Double_t sumEffWeightVals = 0;

            for (Int_t iRun = 0; iRun < nMcRuns; iRun++)
            {
                sumEffWeights = sumEffWeights + OmegaRunInfo[iRun].effWeightCommonMult[iNewBin][ptBinOm];
                sumEffWeightVals = sumEffWeightVals + (OmegaRunInfo[iRun].effWeightCommonMult[iNewBin][ptBinOm] * OmegaRunInfo[iRun].effValCommonMult[iNewBin][ptBinOm]);
            }

            effOmCAvgCommonMult[iNewBin][ptBinOm] = sumEffWeightVals / sumEffWeights;
            effOmCAvgCommonMultErr[iNewBin][ptBinOm] = pow(sumEffWeights, -0.5);
            eff_omC_avg_commonMult[iNewBin]->SetBinContent(ptBinOm + 1, effOmCAvgCommonMult[iNewBin][ptBinOm]);
            eff_omC_avg_commonMult[iNewBin]->SetBinError(ptBinOm + 1, effOmCAvgCommonMultErr[iNewBin][ptBinOm]);
        }
        outputFile->cd("dirEffPtMultRatio_omC");
        eff_omC_avg_commonMult[iNewBin]->SetMarkerStyle(openMarkerStyles[iNewBin]);
        eff_omC_avg_commonMult[iNewBin]->Write();
        hs_omC_avg_commonMult->Add(eff_omC_avg_commonMult[iNewBin]);
    }
    // eff_omC_avg->SetTitle("V0A 0-100%");
    hs_omC_avg_commonMult->Add(eff_omC_avg);

    if (analyseCommonMultDiffs)
    {

        for (Int_t iNewBin = 0; iNewBin < nNewCommonMultBins; iNewBin++)
        {
            /// calculate ratio of mult binned and mult integrated avg efficiencies
            eff_xiC_avg_ratio_commonMult[iNewBin] = (TH1D *)eff_xiC_avg_commonMult[iNewBin]->Clone(TString::Format(("eff_xiC_avg_ratio_commonMult[%d]"), iNewBin));
            eff_omC_avg_ratio_commonMult[iNewBin] = (TH1D *)eff_omC_avg_commonMult[iNewBin]->Clone(TString::Format(("eff_omC_avg_ratio_commonMult[%d]"), iNewBin));

            eff_xiC_avg_ratio_commonMult[iNewBin]->Divide(eff_xiC_avg);
            eff_omC_avg_ratio_commonMult[iNewBin]->Divide(eff_omC_avg);

            // save the histograms and add them to the stack
            outputFile->cd("dirEffPtMultRatio_xiC");
            eff_xiC_avg_ratio_commonMult[iNewBin]->SetMarkerStyle(openMarkerStyles[iNewBin]);
            eff_xiC_avg_ratio_commonMult[iNewBin]->Write();
            eff_xiC_avg_ratio_commonMult[iNewBin]->SetName(TString::Format(("%s: Common Mult Bin %d"), "Average Efficiency #Xi^{+} + #Xi^{-}", iNewBin));
            hs_xiC_AvgEff_ratio_commonMult->Add(eff_xiC_avg_ratio_commonMult[iNewBin]);

            outputFile->cd("dirEffPtMultRatio_omC");
            eff_omC_avg_ratio_commonMult[iNewBin]->SetMarkerStyle(openMarkerStyles[iNewBin]);
            eff_omC_avg_ratio_commonMult[iNewBin]->Write();
            eff_omC_avg_ratio_commonMult[iNewBin]->SetName(TString::Format(("%s: Common Mult Bin %d"), "Average Efficiency #Omega^{+} + #Omega^{-}", iNewBin));
            hs_omC_AvgEff_ratio_commonMult->Add(eff_omC_avg_ratio_commonMult[iNewBin]);
        }
    }
        /// Create and add empty histogram to hs_***_AvgEff_ratio_commonMult with same binning to skip black (final) colour of palette:
    TH1D *emptyHist_xiC_AvgEff_ratio_commonMult = new TH1D("", "", fNptbins_Xi, fPtbins_Xi);
    TH1D *emptyHist_omC_AvgEff_ratio_commonMult = new TH1D("", "", fNptbins_Om, fPtbins_Om);
    hs_xiC_AvgEff_ratio_commonMult->Add(emptyHist_xiC_AvgEff_ratio_commonMult);
    hs_omC_AvgEff_ratio_commonMult->Add(emptyHist_omC_AvgEff_ratio_commonMult);

    eff_xiC_avg->GetXaxis()->SetTitle("#it{p}_{T} (GeV/c)");
    eff_omC_avg->GetXaxis()->SetTitle("#it{p}_{T} (GeV/c)");
    eff_xiC_avg->GetYaxis()->SetTitle("Efficiency");
    eff_omC_avg->GetYaxis()->SetTitle("Efficiency");
    eff_xiC_avg->SetLineColor(kMagenta);
    eff_omC_avg->SetLineColor(kMagenta);
    eff_xiC_avg->SetLineWidth(2);
    eff_omC_avg->SetLineWidth(2);
    eff_xiC_avg->SetMarkerStyle(kFullDiamond);
    eff_omC_avg->SetMarkerStyle(kFullDiamond);
    eff_xiC_avg->SetMarkerSize(3);
    eff_omC_avg->SetMarkerSize(3);

    outputFile->cd();
    eff_xiC_avg->Write();
    eff_omC_avg->Write();
    hs_xiC_eff->Add(eff_xiC_avg);
    hs_omC_eff->Add(eff_omC_avg);

    /// Compute ratio of MC run vs MC average efficiency and add it to hs_***_effRatio_MCvsAvg:
    for (Int_t iRun = 0; iRun < nMcRuns; iRun++)
    {
        eff_xiC_avg_ratio[iRun] = (TH1D *)eff_xiC[iRun]->Clone(TString::Format(("eff_xiC_avg_ratio[%d]"), iRun));
        eff_omC_avg_ratio[iRun] = (TH1D *)eff_omC[iRun]->Clone(TString::Format(("eff_omC_avg_ratio[%d]"), iRun));

        eff_xiC_avg_ratio[iRun]->Divide(eff_xiC_avg);
        eff_omC_avg_ratio[iRun]->Divide(eff_omC_avg);

        hs_xiC_effRatio_MCvsAvg->Add(eff_xiC_avg_ratio[iRun]);
        hs_omC_effRatio_MCvsAvg->Add(eff_omC_avg_ratio[iRun]);
    }
    /// Create and add empty histogram to hs_***_effRatio_MCvsAvg with same binning to skip black (final) colour of palette:
    TH1D *emptyHist_xiC_effRatio_MCvsAvg = new TH1D("", "", fNptbins_Xi, fPtbins_Xi);
    TH1D *emptyHist_omC_effRatio_MCvsAvg = new TH1D("", "", fNptbins_Om, fPtbins_Om);
    hs_xiC_effRatio_MCvsAvg->Add(emptyHist_xiC_effRatio_MCvsAvg);
    hs_omC_effRatio_MCvsAvg->Add(emptyHist_omC_effRatio_MCvsAvg);

    /// Check if average efficiency is okay:
    // Here xi and δxi are the value and error reported by the ith experiment, and the sums run over the N experiments. We then calculate χ2 and compare it with N − 1, which is the expectation value of χ2 if the measurements are from a Gaussian distribution.
    // If χ2/(N − 1) is less than or equal to 1, and there are no known problems with the data, we accept the results.
    // If χ2/(N − 1) is very large, we may choose not to use the average at all. Alternatively, we may quote the calculated average, but then make an educated guess of the error, a conservative estimate designed to take into account known problems with the data.
    // Finally, if χ2/(N − 1) is greater than 1, but not greatly so, we still average the data, but do some more steps...
    // Full text: https://pdg.lbl.gov/2019/reviews/rpp2019-rev-rpp-intro.pdf (p.15)

    /// check if average is okay
    Double_t chi2Xi[fNptbins_Xi];
    Double_t chi2Om[fNptbins_Om];
    Double_t N_Xi = 0;
    Double_t N_Om = 0;

    for (Int_t ptBinXi = 0; ptBinXi < fNptbins_Xi; ptBinXi++)
    {
        Double_t sumChi2 = 0;
        N_Xi = 0;
        for (Int_t iRun = 0; iRun < nMcRuns; iRun++)
        {
            sumChi2 = sumChi2 + XiRunInfo[iRun].effWeight[ptBinXi] * pow(effXiAvg[ptBinXi] - XiRunInfo[iRun].effVal[ptBinXi], 2);
            N_Xi++;
        }
        chi2Xi[ptBinXi] = sumChi2;
    }

    for (Int_t ptBinOm = 0; ptBinOm < fNptbins_Om; ptBinOm++)
    {
        Double_t sumChi2 = 0;
        N_Om = 0;
        for (Int_t iRun = 0; iRun < nMcRuns; iRun++)
        {
            sumChi2 = sumChi2 + OmegaRunInfo[iRun].effWeight[ptBinOm] * pow(effOmAvg[ptBinOm] - OmegaRunInfo[iRun].effVal[ptBinOm], 2);
            N_Om++;
        }
        chi2Om[ptBinOm] = sumChi2;
    }

    for (Int_t ptBinXi = 0; ptBinXi < fNptbins_Xi; ptBinXi++)
    {
        eff_xiC_avg_chi2->SetBinContent(ptBinXi + 1, (chi2Xi[ptBinXi] / (N_Xi - 1)));
    }

    for (Int_t ptBinOm = 0; ptBinOm < fNptbins_Om; ptBinOm++)
    {
        eff_omC_avg_chi2->SetBinContent(ptBinOm + 1, (chi2Om[ptBinOm] / (N_Om - 1)));
    }

    eff_xiC_avg_chi2->GetXaxis()->SetTitle("#it{p}_{T} (GeV/c)");
    eff_omC_avg_chi2->GetXaxis()->SetTitle("#it{p}_{T} (GeV/c)");
    eff_xiC_avg_chi2->GetYaxis()->SetTitle("#chi^{2}/(N-1)");
    eff_omC_avg_chi2->GetYaxis()->SetTitle("#chi^{2}/(N-1)");
    outputFile->cd();
    eff_xiC_avg_chi2->Write();
    eff_omC_avg_chi2->Write();

    /// Write the histogram stacks to the output file
    outputFile->cd("dirHistStacks");
    hs_xip_eff->Write();
    hs_xim_eff->Write();
    hs_omp_eff->Write();
    hs_omm_eff->Write();
    hs_xiC_eff->Write();
    hs_omC_eff->Write();

    outputFile->cd("dirHistStackRatios");
    hs_xiC_effRatio_MCvsAvg->Write();
    hs_omC_effRatio_MCvsAvg->Write();
    /// Draw the histogram stacks
    gStyle->SetOptStat(0);
    gStyle->SetPalette(kCMYK);
    TString histImageOutFolder = outputFolder + "/" + histName;
    TCanvas *cEff[6];
    for (Int_t iCanvas = 0; iCanvas < 6; iCanvas++)
    {
        cEff[iCanvas] = new TCanvas(TString::Format("cEff%d", iCanvas), TString::Format("cEff%d", iCanvas), 2560, 1440);
    }

    PaintStack(*cEff[0], *hs_xip_eff, kFALSE, "Efficiency #times Acceptance");
    PaintStack(*cEff[1], *hs_omp_eff, kFALSE, "Efficiency #times Acceptance");
    PaintStack(*cEff[2], *hs_xim_eff, kFALSE, "Efficiency #times Acceptance");
    PaintStack(*cEff[3], *hs_omm_eff, kFALSE, "Efficiency #times Acceptance");
    PaintStack(*cEff[4], *hs_xiC_eff, kFALSE, "Efficiency #times Acceptance");
    PaintStack(*cEff[5], *hs_omC_eff, kFALSE, "Efficiency #times Acceptance");

    if (saveStack)
    {
        cEff[0]->cd();
        SaveImage(histImageOutFolder, "eff", "eff_xipN", imageFormat, cEff[0]);

        cEff[1]->cd();
        SaveImage(histImageOutFolder, "eff", "eff_ompN", imageFormat, cEff[1]);

        cEff[2]->cd();
        SaveImage(histImageOutFolder, "eff", "eff_ximN", imageFormat, cEff[2]);

        cEff[3]->cd();
        SaveImage(histImageOutFolder, "eff", "eff_ommN", imageFormat, cEff[3]);

        cEff[4]->cd();
        SaveImage(histImageOutFolder, "eff", "eff_xicN", imageFormat, cEff[4]);

        cEff[5]->cd();
        SaveImage(histImageOutFolder, "eff", "eff_omcN", imageFormat, cEff[5]);
    }
    hs_xiC_effRatio_MCvsAvg->SetMinimum(0.6);
    hs_xiC_effRatio_MCvsAvg->SetMaximum(1.4);
    hs_omC_effRatio_MCvsAvg->SetMinimum(0.6);
    hs_omC_effRatio_MCvsAvg->SetMaximum(1.4);
    hs_xiC_effRatio_MCvsAvg->SetTitle("#Xi^{#pm}");
    hs_omC_effRatio_MCvsAvg->SetTitle("#Omega^{#pm}");
    PaintStack(*cEff[0], *hs_xiC_effRatio_MCvsAvg, kFALSE, "MC Production / MC Average");
    PaintStack(*cEff[1], *hs_omC_effRatio_MCvsAvg, kFALSE, "MC Production / MC Average");
    if (saveStack)
    {
        SaveImage(histImageOutFolder, "effMCvsAvgRatio", hs_xiC_effRatio_MCvsAvg->GetName(), imageFormat, cEff[0]);
        SaveImage(histImageOutFolder, "effMCvsAvgRatio", hs_omC_effRatio_MCvsAvg->GetName(), imageFormat, cEff[1]);
    }
    hs_xiC_AvgEff_ratio_commonMult->SetMinimum(0.8);
    hs_xiC_AvgEff_ratio_commonMult->SetMaximum(1.2);
    hs_omC_AvgEff_ratio_commonMult->SetMinimum(0.8);
    hs_omC_AvgEff_ratio_commonMult->SetMaximum(1.2);
    hs_xiC_AvgEff_ratio_commonMult->SetTitle("#Xi^{#pm}");
    hs_omC_AvgEff_ratio_commonMult->SetTitle("#Omega^{#pm}");
    PaintStack(*cEff[2], *hs_xiC_AvgEff_ratio_commonMult, kFALSE, "Multiplicity class / 0-100%");
    PaintStack(*cEff[3], *hs_omC_AvgEff_ratio_commonMult, kFALSE, "Multiplicity class / 0-100%");
    if (saveStack)
    {
        SaveImage(histImageOutFolder, "effCommonMultRatio", hs_xiC_AvgEff_ratio_commonMult->GetName(), imageFormat, cEff[2]);
        SaveImage(histImageOutFolder, "effCommonMultRatio", hs_omC_AvgEff_ratio_commonMult->GetName(), imageFormat, cEff[3]);
    }

    if (analyseCommonMultDiffs)
    {
        // hs_xiC_avg_commonMult->SetMinimum(0.8);
        // hs_xiC_avg_commonMult->SetMaximum(1.2);
        // hs_omC_avg_commonMult->SetMinimum(0.8);
        // hs_omC_avg_commonMult->SetMaximum(1.2);
        hs_xiC_avg_commonMult->SetTitle("#Xi^{#pm}");
        hs_omC_avg_commonMult->SetTitle("#Omega^{#pm}");
        PaintStack(*cEff[4], *hs_xiC_avg_commonMult, kFALSE, "Efficiency #times Acceptance");
        PaintStack(*cEff[5], *hs_omC_avg_commonMult, kFALSE, "Efficiency #times Acceptance");
        if (saveStack)
        {
            SaveImage(histImageOutFolder, "effAvgCommonMult", hs_xiC_avg_commonMult->GetName(), imageFormat, cEff[4]);
            SaveImage(histImageOutFolder, "effAvgCommonMult", hs_omC_avg_commonMult->GetName(), imageFormat, cEff[5]);
        }
    }
    // delete canvas objects
    for (Int_t iCanvas = 0; iCanvas < 6; iCanvas++)
    {
        delete cEff[iCanvas];
        // delete gROOT->FindObject(TString::Format("cEff%d", iCanvas));
    }

    if (analyseDiffs)
    {
        TCanvas *cDiff[nMcRuns];

        for (Int_t iCanvas = 0; iCanvas < nMcRuns; iCanvas++)
        {
            outputFile->cd("dirHistStackRatios");
            hs_xiC_eff_ratio[iCanvas]->Write();
            hs_omC_eff_ratio[iCanvas]->Write();

            cDiff[iCanvas] = new TCanvas(TString::Format("cDiff%d", iCanvas), TString::Format("cDiff%d", iCanvas), 2560, 1440);

            PaintStack(*cDiff[iCanvas], *hs_xiC_eff_ratio[iCanvas], kFALSE, "Efficiency Ratio");
            if (saveStack)
                SaveImage(histImageOutFolder, "effMultRatio", hs_xiC_eff_ratio[iCanvas]->GetName(), imageFormat, cDiff[iCanvas]);

            PaintStack(*cDiff[iCanvas], *hs_omC_eff_ratio[iCanvas], kFALSE, "Efficiency Ratio");
            if (saveStack)
                SaveImage(histImageOutFolder, "effMultRatio", hs_omC_eff_ratio[iCanvas]->GetName(), imageFormat, cDiff[iCanvas]);
        }

        // delete canvas objects
        for (Int_t iCanvas = 0; iCanvas < nMcRuns; iCanvas++)
        {
            delete cDiff[iCanvas];
            // delete gROOT->FindObject(TString::Format("cDiff%d", iCanvas));
        }
    }

    delete outputFile;

    return 0;
}
