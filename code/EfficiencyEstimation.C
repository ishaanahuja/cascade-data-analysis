#include <TROOT.h>
#include <TStyle.h>
#include "CascadeUtils.h"

struct RunInfo
{
    TString runName;
    std::vector<Double_t> effVal;    // [MC Run][ptbin] : efficiency values in mult integrated (0-100%) Combined (+-) case
    std::vector<Double_t> effWeight; // [MC Run][ptbin] : efficiency weights in mult integrated (0-100%) Combined (+-) case
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
    TString outputFilename = "TEST_effEst.root",
    TString outputFolder = "TEST_effEst_images",
    Bool_t saveStack = kTRUE,
    TString imageFormat = "png",
    Int_t verbosity = kInfo)
{
    ROOT::EnableImplicitMT();
    gStyle->SetOptFit(1111);
    gErrorIgnoreLevel = verbosity;

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

    /// Calculate avg efficiency:

    Double_t effXiAvg[fNptbins_Xi];
    Double_t effXiErrAvg[fNptbins_Xi];
    Double_t effOmAvg[fNptbins_Om];
    Double_t effOmErrAvg[fNptbins_Om];

    TH1D *eff_xiC_avg = new TH1D("eff_xiC_avg", "#Xi^{#pm} Avg MC efficiency", fNptbins_Xi, fPtbins_Xi);
    TH1D *eff_omC_avg = new TH1D("eff_omC_avg", "#Omega^{#pm} Avg MC efficiency", fNptbins_Om, fPtbins_Om);

    TH1D *eff_xiC_avg_chi2 = new TH1D("eff_xiC_avg_chi2", "#Xi^{#pm} Avg MC efficiency: (#frac{#chi^{2}}{N-1})", fNptbins_Xi, fPtbins_Xi);
    TH1D *eff_omC_avg_chi2 = new TH1D("eff_omC_avg_chi2", "#Omega^{#pm} Avg MC efficiency: (#frac{#chi^{2}}{N-1})", fNptbins_Om, fPtbins_Om);

    // THstack efficiency
    auto hs_xip_eff = new THStack("hs_xip_eff", "Efficiency #Xi^{+}");                 // mult integrated efficiency + mult binned efficiencies
    auto hs_xim_eff = new THStack("hs_xim_eff", "Efficiency #Xi^{-}");                 // mult integrated efficiency + mult binned efficiencies
    auto hs_omp_eff = new THStack("hs_omp_eff", "Efficiency #Omega^{+}");              // mult integrated efficiency + mult binned efficiencies
    auto hs_omm_eff = new THStack("hs_omm_eff", "Efficiency #Omega^{-}");              // mult integrated efficiency + mult binned efficiencies
    auto hs_xiC_eff = new THStack("hs_xiC_eff", "Efficiency #Xi^{+} + #Xi^{-}");       // mult integrated efficiency + mult binned efficiencies
    auto hs_omC_eff = new THStack("hs_omC_eff", "Efficiency #Omega^{+} + #Omega^{-}"); // mult integrated efficiency + mult binned efficiencies

    THStack *hs_xiC_eff_ratio[nMcRuns]; // ratio of (mult integrated / mult binned) efficiencies
    THStack *hs_omC_eff_ratio[nMcRuns]; // ratio of (mult integrated / mult binned) efficiencies

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

        Info("EfficiencyEstimation: histInput", "Getting histograms for '%s': %s", McRunName.Data(), histName.Data());
        TFile *inputFile = OpenFile(fileName);
        if (!inputFile)
        {
            Error("EfficiencyEstimation: inputFile", "Cannot open file '%s' !", fileName.Data());
            return 1;
        }

        if (!McRunName.Contains("Om")) // get Xi hists if filename is NOT flagged for omega only
        {

            XiRunInfo[iRunXi].runName = McRunName;
            XiRunInfo[iRunXi].effVal.resize(fNptbins_Xi);
            XiRunInfo[iRunXi].effWeight.resize(fNptbins_Xi);

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
            OmegaRunInfo[iRunOm].runName = McRunName;
            OmegaRunInfo[iRunOm].effVal.resize(fNptbins_Om);
            OmegaRunInfo[iRunOm].effWeight.resize(fNptbins_Om);

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

    Info("DrawCascades: outputFile", "Output file created: %s", outputFilename.Data());

    /// Begin:

    for (Int_t iRun = 0; iRun < nMcRuns; iRun++)
    {

        /// CASE: Xi
        Info("EfficiencyEstimation", "%d: %s Plotting Xi efficiency ...", iRun, XiRunInfo[iRun].runName.Data());
        /// generate mult integrated efficiency histograms
        eff_xim[iRun] = new TH1D(TString::Format(("eff_xim[%d]"), iRun), TString::Format("%s Efficiency: Mult 0-100%%", XiRunInfo[iRun].runName.Data()), fNptbins_Xi, fPtbins_Xi);
        eff_xip[iRun] = new TH1D(TString::Format(("eff_xip[%d]"), iRun), TString::Format("%s Efficiency: Mult 0-100%%", XiRunInfo[iRun].runName.Data()), fNptbins_Xi, fPtbins_Xi);
        eff_xiC[iRun] = new TH1D(TString::Format(("eff_xiC[%d]"), iRun), TString::Format("%s Efficiency: Mult 0-100%%", XiRunInfo[iRun].runName.Data()), fNptbins_Xi, fPtbins_Xi);

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
        eff_xim[iRun]->SetMarkerStyle(markerStyles[iRun]);
        eff_xim[iRun]->Write();
        eff_xim[iRun]->SetName(TString::Format(("%s: Mult: %.0f-%.0f%%"), XiRunInfo[iRun].runName.Data(), fMultbins_Xi[0], fMultbins_Xi[fNmultbins_Xi]));
        hs_xim_eff->Add(eff_xim[iRun]);

        outputFile->cd("dirEffPt_xip");
        eff_xip[iRun]->SetMarkerStyle(markerStyles[iRun]);
        eff_xip[iRun]->Write();
        eff_xip[iRun]->SetName(TString::Format(("%s: Mult: %.0f-%.0f%%"), XiRunInfo[iRun].runName.Data(), fMultbins_Xi[0], fMultbins_Xi[fNmultbins_Xi]));
        hs_xip_eff->Add(eff_xip[iRun]);

        outputFile->cd("dirEffPt_xiC");
        eff_xiC[iRun]->SetMarkerStyle(markerStyles[iRun]);
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
                hs_xim_eff->Add(eff_xim_mult[iRun][multBinXi]);

                outputFile->cd("dirEffPtMult_xip");
                eff_xip_mult[iRun][multBinXi]->SetMarkerStyle(markerStyles[iRun]);
                eff_xip_mult[iRun][multBinXi]->Write();
                eff_xip_mult[iRun][multBinXi]->SetName(TString::Format(("%s: Mult: %.0f-%.0f%%"), XiRunInfo[iRun].runName.Data(), fMultbins_Xi[multBinXi], fMultbins_Xi[multBinXi + 1]));
                hs_xip_eff->Add(eff_xip_mult[iRun][multBinXi]);

                outputFile->cd("dirEffPtMult_xiC");
                eff_xiC_mult[iRun][multBinXi]->SetMarkerStyle(markerStyles[iRun]);
                eff_xiC_mult[iRun][multBinXi]->Write();
                eff_xiC_mult[iRun][multBinXi]->SetName(TString::Format(("%s: Mult: %.0f-%.0f%%"), XiRunInfo[iRun].runName.Data(), fMultbins_Xi[multBinXi], fMultbins_Xi[multBinXi + 1]));
                hs_xiC_eff->Add(eff_xiC_mult[iRun][multBinXi]);

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
        eff_omm[iRun] = new TH1D(TString::Format(("eff_omm[%d]"), iRun), TString::Format("%s Efficiency: Mult 0-100%%", OmegaRunInfo[iRun].runName.Data()), fNptbins_Om, fPtbins_Om);
        eff_omp[iRun] = new TH1D(TString::Format(("eff_omp[%d]"), iRun), TString::Format("%s Efficiency: Mult 0-100%%", OmegaRunInfo[iRun].runName.Data()), fNptbins_Om, fPtbins_Om);
        eff_omC[iRun] = new TH1D(TString::Format(("eff_omC[%d]"), iRun), TString::Format("%s Efficiency: Mult 0-100%%", OmegaRunInfo[iRun].runName.Data()), fNptbins_Om, fPtbins_Om);

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
        eff_omm[iRun]->SetMarkerStyle(markerStyles[iRun]);
        eff_omm[iRun]->Write();
        eff_omm[iRun]->SetName(TString::Format(("%s: Mult: %.0f-%.0f%%"), OmegaRunInfo[iRun].runName.Data(), fMultbins_Om[0], fMultbins_Om[fNmultbins_Om]));
        hs_omm_eff->Add(eff_omm[iRun]);

        outputFile->cd("dirEffPt_omp");
        eff_omp[iRun]->SetMarkerStyle(markerStyles[iRun]);
        eff_omp[iRun]->Write();
        eff_omp[iRun]->SetName(TString::Format(("%s: Mult: %.0f-%.0f%%"), OmegaRunInfo[iRun].runName.Data(), fMultbins_Om[0], fMultbins_Om[fNmultbins_Om]));
        hs_omp_eff->Add(eff_omp[iRun]);

        outputFile->cd("dirEffPt_omC");
        eff_omC[iRun]->SetMarkerStyle(markerStyles[iRun]);
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
                hs_omm_eff->Add(eff_omm_mult[iRun][multBinOm]);

                outputFile->cd("dirEffPtMult_omp");
                eff_omp_mult[iRun][multBinOm]->SetMarkerStyle(markerStyles[iRun]);
                eff_omp_mult[iRun][multBinOm]->Write();
                eff_omp_mult[iRun][multBinOm]->SetName(TString::Format(("%s: Mult: %.0f-%.0f%%"), OmegaRunInfo[iRun].runName.Data(), fMultbins_Om[multBinOm], fMultbins_Om[multBinOm + 1]));
                hs_omp_eff->Add(eff_omp_mult[iRun][multBinOm]);

                outputFile->cd("dirEffPtMult_omC");
                eff_omC_mult[iRun][multBinOm]->SetMarkerStyle(markerStyles[iRun]);
                eff_omC_mult[iRun][multBinOm]->Write();
                eff_omC_mult[iRun][multBinOm]->SetName(TString::Format(("%s: Mult: %.0f-%.0f%%"), OmegaRunInfo[iRun].runName.Data(), fMultbins_Om[multBinOm], fMultbins_Om[multBinOm + 1]));
                hs_omC_eff->Add(eff_omC_mult[iRun][multBinOm]);

                outputFile->cd("dirEffPtMultRatio_omC");
                eff_omC_ratio_mult[iRun][multBinOm]->SetMarkerStyle(markerStyles[iRun]);
                eff_omC_ratio_mult[iRun][multBinOm]->Write();
                eff_omC_ratio_mult[iRun][multBinOm]->SetName(TString::Format(("%s: Mult: %.0f-%.0f%%"), OmegaRunInfo[iRun].runName.Data(), fMultbins_Om[multBinOm], fMultbins_Om[multBinOm + 1]));
                hs_omC_eff_ratio[iRun]->Add(eff_omC_ratio_mult[iRun][multBinOm]);
            }
        }
    }

    /// Calculate avg efficiency: Xi
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

    /// Calculate avg efficiency: Omega
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

    eff_xiC_avg->GetXaxis()->SetTitle("#it{p}_{T} (GeV/c)");
    eff_omC_avg->GetXaxis()->SetTitle("#it{p}_{T} (GeV/c)");
    eff_xiC_avg->GetYaxis()->SetTitle("Efficiency");
    eff_omC_avg->GetYaxis()->SetTitle("Efficiency");
    eff_xiC_avg->SetLineColor(kMagenta);
    eff_omC_avg->SetLineColor(kMagenta);
    eff_xiC_avg->SetLineWidth(3);
    eff_omC_avg->SetLineWidth(3);
    eff_xiC_avg->SetMarkerStyle(kFullCircle);
    eff_omC_avg->SetMarkerStyle(kFullCircle);

    outputFile->cd();
    eff_xiC_avg->Write();
    eff_omC_avg->Write();
    hs_xiC_eff->Add(eff_xiC_avg);
    hs_omC_eff->Add(eff_omC_avg);

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

    /// Draw the histogram stacks
    gStyle->SetOptStat(0);
    gStyle->SetPalette(kVisibleSpectrum);
    TString histImageOutFolder = outputFolder + "/" + histName;
    TCanvas *cEff[6];
    for (Int_t iCanvas = 0; iCanvas < 6; iCanvas++)
    {
        cEff[iCanvas] = new TCanvas(TString::Format("cEff%d", iCanvas), TString::Format("cEff%d", iCanvas), 2560, 1440);
    }

    PaintStack(*cEff[0], *hs_xip_eff, kFALSE, "Efficiency");
    PaintStack(*cEff[1], *hs_omp_eff, kFALSE, "Efficiency");
    PaintStack(*cEff[2], *hs_xim_eff, kFALSE, "Efficiency");
    PaintStack(*cEff[3], *hs_omm_eff, kFALSE, "Efficiency");
    PaintStack(*cEff[4], *hs_xiC_eff, kFALSE, "Efficiency");
    PaintStack(*cEff[5], *hs_omC_eff, kFALSE, "Efficiency");

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
            PaintStack(*cDiff[iCanvas], *hs_omC_eff_ratio[iCanvas], kFALSE, "Efficiency Ratio");
        }

        if (saveStack)
        {
            for (Int_t iCanvas = 0; iCanvas < nMcRuns; iCanvas++)
            {
                cDiff[iCanvas]->cd();
                SaveImage(histImageOutFolder, "effMultRatio", hs_xiC_eff_ratio[iCanvas]->GetName(), imageFormat, cDiff[iCanvas]);

                cDiff[iCanvas]->cd();
                SaveImage(histImageOutFolder, "effMultRatio", hs_omC_eff_ratio[iCanvas]->GetName(), imageFormat, cDiff[iCanvas]);
            }
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
