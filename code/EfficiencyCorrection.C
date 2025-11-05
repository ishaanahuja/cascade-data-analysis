/**
 * @file EfficiencyCorrection.C
 * @author Ishaan Ahuja (ishaanahuja0@gmail.com)
 * Thesis: Multi-strange particle production in p–Pb collisions at √sNN = 8.16 TeV
 * DOI: https://doi.org/10.17181/cwcde-g1z94
 *
 * @brief Perform efficiency correction on particle spectra.
 * @version 1
 * @date 04-06-2025
 *
 * @copyright Copyright (c) 2025
 *
 */

#include <TROOT.h>
#include <TStyle.h>

#include "CascadeUtils.h"

/**
 * @brief Perform efficiency correction on particle spectra.
 *
 * Generates raw pt spectra, gets efficiency values from effEst.root (MC Avg 0-100% xiC/omC),
 * and divides raw spectra by efficiency to generate corrected spectra with stat errors only.
 *
 * This function reads histograms from fit and corresponding efficiency input files,
 * applies efficiency corrections, and saves the corrected histograms to an output
 * file. It also generates and saves images of the histograms if specified.
 *
 * @param fitInputFilePrefix Prefix for the fit input file path.
 * @param effInputFilePrefix Prefix for the efficiency input file path.
 * @param histName Name of the histogram to be processed.
 * @param outputFileName Name of the output file to save the results. If empty, a default name is generated.
 * @param outputFolder Folder to save the output files.
 * @param saveStack Boolean flag to indicate whether to save histogram stack images.
 * @param imageFormat Format of the saved images (e.g., "png").
 * @param verbosity Verbosity level for logging.
 * @return int Status code (0 for success).
 */
int EfficiencyCorrection(
    TString fitInputFilePrefix = "/var/home/ishaan/Work/git/analysis/results/RandomVars/050225_6RunsRandDiff_fit/050225_6RunsRandDiff",                                    //"/var/home/ishaan/Work/git/analysis/results/RandomVars/190225_SysUncertSigExt_Fit/190225_SysSigExt_DGP2_15bg_6Runs.root",
    TString effInputFilePrefix = "/var/home/ishaan/Work/git/analysis/results/RandomVars/030225_EfficiencyEstimation_Vars/030225_eff",                                      //"/var/home/ishaan/Work/git/analysis/results/RandomVars/200225_SysSigExt_effEst/200225_SysSigExt_DGP2_15bg_effEst.root",
    TString histName = "h3_ptmasscent_def",                                                                                                                                //"DGP2_15bg",
    TString outputFileName = "/var/home/ishaan/Work/git/analysis/results/RandomVars/040625_effCorr_multIntCorrected_v2/040625_effCorr_multIntCorrected_def_6runs_v2.root", //"/var/home/ishaan/Work/git/analysis/results/RandomVars/200225_SysSigExt_effCorr/200225_SysSigExt_DGP2_15bg_effCorr.root",
    TString outputFolder = "/var/home/ishaan/Work/git/analysis/results/RandomVars/040625_effCorr_multIntCorrected_v2",                                                     //"/var/home/ishaan/Work/git/analysis/results/RandomVars/200225_SysSigExt_effCorr",
    Bool_t saveStack = kTRUE,
    TString imageFormat = "png",
    Int_t verbosity = kInfo)
{
    ROOT::EnableImplicitMT();
    gStyle->SetOptFit(1111);
    gErrorIgnoreLevel = verbosity;

    // remove ownership of objects from file so we can delete the file ptr
    TH1::AddDirectory(kFALSE);
    TH1::SetDefaultSumw2(kTRUE);
    outputFolder = SetOutputFolder(outputFolder);

    TH1 *h_multBinEntries_Xi;
    TH1 *h_multBinEntries_Om;

    TH1 *resultParXip_pt[fNptbins_Xi];
    TH1 *resultParXim_pt[fNptbins_Xi];
    TH1 *resultParOmp_pt[fNptbins_Om];
    TH1 *resultParOmm_pt[fNptbins_Om];
    TH1 *resultParXiC_pt[fNptbins_Xi];
    TH1 *resultParOmC_pt[fNptbins_Om];

    TH1D *rawPt_xim;
    TH1D *rawPt_xip;
    TH1D *rawPt_omm;
    TH1D *rawPt_omp;
    TH1D *rawPt_xiC;
    TH1D *rawPt_omC;

    TH1D *effCorrPt_xim;
    TH1D *effCorrPt_xip;
    TH1D *effCorrPt_omm;
    TH1D *effCorrPt_omp;
    TH1D *effCorrPt_xiC;
    TH1D *effCorrPt_omC;

    TH1 *resultParXip_pt_mult[fNptbins_Xi][fNmultbins_Xi];
    TH1 *resultParXim_pt_mult[fNptbins_Xi][fNmultbins_Xi];
    TH1 *resultParOmp_pt_mult[fNptbins_Om][fNmultbins_Om];
    TH1 *resultParOmm_pt_mult[fNptbins_Om][fNmultbins_Om];
    TH1 *resultParXiC_pt_mult[fNptbins_Xi][fNmultbins_Xi];
    TH1 *resultParOmC_pt_mult[fNptbins_Om][fNmultbins_Om];

    TH1D *rawPt_xim_mult[fNmultbins_Xi];
    TH1D *rawPt_xip_mult[fNmultbins_Xi];
    TH1D *rawPt_omm_mult[fNmultbins_Om];
    TH1D *rawPt_omp_mult[fNmultbins_Om];
    TH1D *rawPt_xiC_mult[fNmultbins_Xi];
    TH1D *rawPt_omC_mult[fNmultbins_Om];

    TH1D *effCorrPt_xim_mult[fNmultbins_Xi];
    TH1D *effCorrPt_xip_mult[fNmultbins_Xi];
    TH1D *effCorrPt_omm_mult[fNmultbins_Om];
    TH1D *effCorrPt_omp_mult[fNmultbins_Om];
    TH1D *effCorrPt_xiC_mult[fNmultbins_Xi];
    TH1D *effCorrPt_omC_mult[fNmultbins_Om];

    TH1 *eff_xiC_avg;
    TH1 *eff_omC_avg;

    // THStacks:
    auto hs_xip_raw = new THStack("hs_xip_raw", "Raw #it{p}_{T} spectra #Xi^{+}");
    auto hs_xim_raw = new THStack("hs_xim_raw", "Raw #it{p}_{T} spectra #Xi^{-}");
    auto hs_omp_raw = new THStack("hs_omp_raw", "Raw #it{p}_{T} spectra #Omega^{+}");
    auto hs_omm_raw = new THStack("hs_omm_raw", "Raw #it{p}_{T} spectra #Omega^{-}");
    auto hs_xiC_raw = new THStack("hs_xiC_raw", "Raw #it{p}_{T} spectra #Xi^{+} + #Xi^{-}");
    auto hs_omC_raw = new THStack("hs_omC_raw", "Raw #it{p}_{T} spectra #Omega^{+} + #Omega^{-}");
    auto hs_xip_effCorr = new THStack("hs_xip_effCorr", "Efficiency corrected #it{p}_{T} spectra #Xi^{+}");
    auto hs_xim_effCorr = new THStack("hs_xim_effCorr", "Efficiency corrected #it{p}_{T} spectra #Xi^{-}");
    auto hs_omp_effCorr = new THStack("hs_omp_effCorr", "Efficiency corrected #it{p}_{T} spectra #Omega^{+}");
    auto hs_omm_effCorr = new THStack("hs_omm_effCorr", "Efficiency corrected #it{p}_{T} spectra #Omega^{-}");
    auto hs_xiC_effCorr = new THStack("hs_xiC_effCorr", "Efficiency corrected #it{p}_{T} spectra #Xi^{+} + #Xi^{-}");
    auto hs_omC_effCorr = new THStack("hs_omC_effCorr", "Efficiency corrected #it{p}_{T} spectra #Omega^{+} + #Omega^{-}");

    TString fitInputFileName = TString::Format("%s_%s.root", fitInputFilePrefix.Data(), histName.Data());
    TString effInputFileName = TString::Format("%s_%s.root", effInputFilePrefix.Data(), histName.Data());

    // TString fitInputFileName = fitInputFilePrefix; // only for SysSigExt
    // TString effInputFileName = effInputFilePrefix; // only for SysSigExt

    /// Getting histograms from fit input file:
    TFile *fitInputFile = OpenFile(fitInputFileName);
    Info("EfficiencyCorrection: fitInput", "Getting histograms from '%s'", fitInputFileName.Data());
    h_multBinEntries_Xi = (TH1 *)fitInputFile->FindObjectAny("h_multBinEntries_Xi");
    h_multBinEntries_Om = (TH1 *)fitInputFile->FindObjectAny("h_multBinEntries_Om");
    for (Int_t ptBinXi = 0; ptBinXi < fNptbins_Xi; ptBinXi++)
    {
        resultParXip_pt[ptBinXi] = (TH1 *)fitInputFile->FindObjectAny(TString::Format(("resultParXip_pt[%d]"), ptBinXi));
        resultParXim_pt[ptBinXi] = (TH1 *)fitInputFile->FindObjectAny(TString::Format(("resultParXim_pt[%d]"), ptBinXi));
        resultParXiC_pt[ptBinXi] = (TH1 *)fitInputFile->FindObjectAny(TString::Format(("resultParXiC_pt[%d]"), ptBinXi));

        for (Int_t multBinXi = 0; multBinXi < fNmultbins_Xi; multBinXi++)
        {
            resultParXip_pt_mult[ptBinXi][multBinXi] = (TH1 *)fitInputFile->FindObjectAny(TString::Format(("resultParXip_pt_mult[%d][%d]"), ptBinXi, multBinXi));
            resultParXim_pt_mult[ptBinXi][multBinXi] = (TH1 *)fitInputFile->FindObjectAny(TString::Format(("resultParXim_pt_mult[%d][%d]"), ptBinXi, multBinXi));
            resultParXiC_pt_mult[ptBinXi][multBinXi] = (TH1 *)fitInputFile->FindObjectAny(TString::Format(("resultParXiC_pt_mult[%d][%d]"), ptBinXi, multBinXi));
        }
    }

    for (Int_t ptBinOm = 0; ptBinOm < fNptbins_Om; ptBinOm++)
    {
        resultParOmp_pt[ptBinOm] = (TH1 *)fitInputFile->FindObjectAny(TString::Format(("resultParOmp_pt[%d]"), ptBinOm));
        resultParOmm_pt[ptBinOm] = (TH1 *)fitInputFile->FindObjectAny(TString::Format(("resultParOmm_pt[%d]"), ptBinOm));
        resultParOmC_pt[ptBinOm] = (TH1 *)fitInputFile->FindObjectAny(TString::Format(("resultParOmC_pt[%d]"), ptBinOm));

        for (Int_t multBinOm = 0; multBinOm < fNmultbins_Om; multBinOm++)
        {
            resultParOmp_pt_mult[ptBinOm][multBinOm] = (TH1 *)fitInputFile->FindObjectAny(TString::Format(("resultParOmp_pt_mult[%d][%d]"), ptBinOm, multBinOm));
            resultParOmm_pt_mult[ptBinOm][multBinOm] = (TH1 *)fitInputFile->FindObjectAny(TString::Format(("resultParOmm_pt_mult[%d][%d]"), ptBinOm, multBinOm));
            resultParOmC_pt_mult[ptBinOm][multBinOm] = (TH1 *)fitInputFile->FindObjectAny(TString::Format(("resultParOmC_pt_mult[%d][%d]"), ptBinOm, multBinOm));
        }
    }
    // input ended for fitInputFile
    delete fitInputFile;

    /// Getting histograms from efficiency input file:
    TFile *effInputFile = OpenFile(effInputFileName);
    Info("EfficiencyCorrection: effInput", "Getting histograms from '%s'", effInputFileName.Data());
    eff_xiC_avg = (TH1 *)effInputFile->FindObjectAny("eff_xiC_avg");
    eff_omC_avg = (TH1 *)effInputFile->FindObjectAny("eff_omC_avg");
    // input ended for effInputFile
    delete effInputFile;

    /// Saving output :
    if (outputFileName.IsNull())
    {
        outputFileName = TString::Format("%s/EffCorr_%s.root", outputFolder.Data(), histName.Data());
    }

    Info("EfficiencyCorrection: I/O", "Input successful. Creating output file '%s'", outputFileName.Data());

    /// Creating output directory structure:
    TFile *outputFile = OpenFile(outputFileName, "RECREATE");
    outputFile->mkdir("dirRawPt_xim");
    outputFile->mkdir("dirRawPt_xip");
    outputFile->mkdir("dirRawPt_omm");
    outputFile->mkdir("dirRawPt_omp");
    outputFile->mkdir("dirRawPt_xiC");
    outputFile->mkdir("dirRawPt_omC");
    outputFile->mkdir("dirEffCorrPt_xim");
    outputFile->mkdir("dirEffCorrPt_xip");
    outputFile->mkdir("dirEffCorrPt_omm");
    outputFile->mkdir("dirEffCorrPt_omp");
    outputFile->mkdir("dirEffCorrPt_xiC");
    outputFile->mkdir("dirEffCorrPt_omC");
    outputFile->mkdir("dirHistStacks");
    /// Output set.

    /// Begin: write avg {+/-} efficiency to file
    outputFile->cd();
    eff_xiC_avg->Write();
    eff_omC_avg->Write();

    // Calculate total multiplicity entries for 0-100% multiplicity range for Xi and Omega:
    Double_t multEntriesXi = h_multBinEntries_Xi->Integral();
    Double_t multEntriesOm = h_multBinEntries_Om->Integral();

    /// XI:
    {
        /// Mult integrated raw + efficiency corrected spectra:
        rawPt_xim = new TH1D("rawPt_xim", "Raw #Xi^{-}: Mult:0-100%", fNptbins_Xi, fPtbins_Xi);
        rawPt_xip = new TH1D("rawPt_xip", "Raw #Xi^{+}: Mult:0-100%", fNptbins_Xi, fPtbins_Xi);
        rawPt_xiC = new TH1D("rawPt_xiC", "Raw #Xi^{+} + #Xi^{-}: Mult:0-100%", fNptbins_Xi, fPtbins_Xi);

        effCorrPt_xim = new TH1D("effCorrPt_xim", "Eff. corrected #Xi^{-}: Mult:0-100%", fNptbins_Xi, fPtbins_Xi);
        effCorrPt_xip = new TH1D("effCorrPt_xip", "Eff. corrected #Xi^{+}: Mult:0-100%", fNptbins_Xi, fPtbins_Xi);
        effCorrPt_xiC = new TH1D("effCorrPt_xiC", "Eff. corrected #Xi^{+} + #Xi^{-}: Mult:0-100%", fNptbins_Xi, fPtbins_Xi);

        for (Int_t ptBinXi = 0; ptBinXi < fNptbins_Xi; ptBinXi++)
        {
            rawPt_xip->SetBinContent(ptBinXi + 1, (resultParXip_pt[ptBinXi]->GetBinContent(1) * 0.97) / ((rawPt_xip->GetBinWidth(ptBinXi + 1)) * multEntriesXi)); // Bin 1 in resultparams is raw pt's bin counting
            rawPt_xip->SetBinError(ptBinXi + 1, (resultParXip_pt[ptBinXi]->GetBinError(1) * 0.97) / ((rawPt_xip->GetBinWidth(ptBinXi + 1)) * multEntriesXi));
            rawPt_xim->SetBinContent(ptBinXi + 1, (resultParXim_pt[ptBinXi]->GetBinContent(1) * 0.97) / ((rawPt_xim->GetBinWidth(ptBinXi + 1)) * multEntriesXi)); // Bin 1 in resultparams is raw pt's bin counting
            rawPt_xim->SetBinError(ptBinXi + 1, (resultParXim_pt[ptBinXi]->GetBinError(1) * 0.97) / ((rawPt_xim->GetBinWidth(ptBinXi + 1)) * multEntriesXi));
            rawPt_xiC->SetBinContent(ptBinXi + 1, (resultParXiC_pt[ptBinXi]->GetBinContent(1) * 0.97) / ((rawPt_xiC->GetBinWidth(ptBinXi + 1)) * multEntriesXi)); // Bin 1 in resultparams is raw pt's bin counting
            rawPt_xiC->SetBinError(ptBinXi + 1, (resultParXiC_pt[ptBinXi]->GetBinError(1) * 0.97) / ((rawPt_xiC->GetBinWidth(ptBinXi + 1)) * multEntriesXi));
        }

        // Xi+: Mult:0-100%: write to file
        outputFile->cd("dirRawPt_xip");
        rawPt_xip->SetMarkerStyle(markerStyles[11]);
        rawPt_xip->Write();
        outputFile->cd("dirEffCorrPt_xip");
        effCorrPt_xip->Divide(rawPt_xip, eff_xiC_avg);
        effCorrPt_xip->SetMarkerStyle(markerStyles[11]);
        effCorrPt_xip->Write();
        // // Xi+: Mult:0-100%: add to stack
        // rawPt_xip->SetTitle("Mult:0-100%");
        // hs_xip_raw->Add(rawPt_xip);
        // effCorrPt_xip->SetTitle("Mult:0-100%");
        // hs_xip_effCorr->Add(effCorrPt_xip);

        // Xi-: Mult:0-100%: write to file
        outputFile->cd("dirRawPt_xim");
        rawPt_xim->SetMarkerStyle(markerStyles[11]);
        rawPt_xim->Write();
        outputFile->cd("dirEffCorrPt_xim");
        effCorrPt_xim->Divide(rawPt_xim, eff_xiC_avg);
        effCorrPt_xim->SetMarkerStyle(markerStyles[11]);
        effCorrPt_xim->Write();
        // // Xi-: Mult:0-100%: add to stack
        // rawPt_xim->SetTitle("Mult:0-100%");
        // hs_xim_raw->Add(rawPt_xim);
        // effCorrPt_xim->SetTitle("Mult:0-100%");
        // hs_xim_effCorr->Add(effCorrPt_xim);

        // XiC: Mult:0-100%: write to file
        outputFile->cd("dirRawPt_xiC");
        rawPt_xiC->SetMarkerStyle(markerStyles[11]);
        rawPt_xiC->Write();
        outputFile->cd("dirEffCorrPt_xiC");
        effCorrPt_xiC->Divide(rawPt_xiC, eff_xiC_avg);
        effCorrPt_xiC->SetMarkerStyle(markerStyles[11]);
        effCorrPt_xiC->Write();
        // // XiC: Mult:0-100%: add to stack
        // rawPt_xiC->SetTitle("Mult:0-100%");
        // hs_xiC_raw->Add(rawPt_xiC);
        // effCorrPt_xiC->SetTitle("Mult:0-100%");
        // hs_xiC_effCorr->Add(effCorrPt_xiC);

        /// Mult binned raw + efficiency corrected spectra:
        for (Int_t multBinXi = 0; multBinXi < fNmultbins_Xi; multBinXi++)
        {
            // raw pt spectra: mult binned (diff)
            rawPt_xip_mult[multBinXi] = new TH1D(TString::Format(("rawPt_xip_mult[%d]"), multBinXi), TString::Format(("Raw #Xi^{+}: Mult: %.0f-%.0f%%"), fMultbins_Xi[multBinXi], fMultbins_Xi[multBinXi + 1]), fNptbins_Xi, fPtbins_Xi);
            rawPt_xim_mult[multBinXi] = new TH1D(TString::Format(("rawPt_xim_mult[%d]"), multBinXi), TString::Format(("Raw #Xi^{-}: Mult: %.0f-%.0f%%"), fMultbins_Xi[multBinXi], fMultbins_Xi[multBinXi + 1]), fNptbins_Xi, fPtbins_Xi);
            rawPt_xiC_mult[multBinXi] = new TH1D(TString::Format(("rawPt_xiC_mult[%d]"), multBinXi), TString::Format(("Raw #Xi^{+} + #Xi^{-}: Mult: %.0f-%.0f%%"), fMultbins_Xi[multBinXi], fMultbins_Xi[multBinXi + 1]), fNptbins_Xi, fPtbins_Xi);
            // effiency corrected pt spectra: mult binned (diff)
            effCorrPt_xip_mult[multBinXi] = new TH1D(TString::Format(("effCorrPt_xip_mult[%d]"), multBinXi), TString::Format(("Eff. corrected #Xi^{+}: Mult: %.0f-%.0f%%"), fMultbins_Xi[multBinXi], fMultbins_Xi[multBinXi + 1]), fNptbins_Xi, fPtbins_Xi);
            effCorrPt_xim_mult[multBinXi] = new TH1D(TString::Format(("effCorrPt_xim_mult[%d]"), multBinXi), TString::Format(("Eff. corrected #Xi^{-}: Mult: %.0f-%.0f%%"), fMultbins_Xi[multBinXi], fMultbins_Xi[multBinXi + 1]), fNptbins_Xi, fPtbins_Xi);
            effCorrPt_xiC_mult[multBinXi] = new TH1D(TString::Format(("effCorrPt_xiC_mult[%d]"), multBinXi), TString::Format(("Eff. corrected #Xi^{+} + #Xi^{-}: Mult: %.0f-%.0f%%"), fMultbins_Xi[multBinXi], fMultbins_Xi[multBinXi + 1]), fNptbins_Xi, fPtbins_Xi);

            for (Int_t ptBinXi = 0; ptBinXi < fNptbins_Xi; ptBinXi++)
            {
                // Trigger and primary vertex reconstruction efficiency = 0.97
                rawPt_xip_mult[multBinXi]->SetBinContent(ptBinXi + 1, (resultParXip_pt_mult[ptBinXi][multBinXi]->GetBinContent(1) * 0.97) / ((rawPt_xip_mult[multBinXi]->GetBinWidth(ptBinXi + 1)) * h_multBinEntries_Xi->GetBinContent(multBinXi + 1))); // Bin 1 in resultparams is raw pt's bin counting
                rawPt_xip_mult[multBinXi]->SetBinError(ptBinXi + 1, (resultParXip_pt_mult[ptBinXi][multBinXi]->GetBinError(1) * 0.97) / ((rawPt_xip_mult[multBinXi]->GetBinWidth(ptBinXi + 1)) * h_multBinEntries_Xi->GetBinContent(multBinXi + 1)));
                rawPt_xim_mult[multBinXi]->SetBinContent(ptBinXi + 1, (resultParXim_pt_mult[ptBinXi][multBinXi]->GetBinContent(1) * 0.97) / ((rawPt_xim_mult[multBinXi]->GetBinWidth(ptBinXi + 1)) * h_multBinEntries_Xi->GetBinContent(multBinXi + 1))); // Bin 1 in resultparams is raw pt's bin counting
                rawPt_xim_mult[multBinXi]->SetBinError(ptBinXi + 1, (resultParXim_pt_mult[ptBinXi][multBinXi]->GetBinError(1) * 0.97) / ((rawPt_xim_mult[multBinXi]->GetBinWidth(ptBinXi + 1)) * h_multBinEntries_Xi->GetBinContent(multBinXi + 1)));
                rawPt_xiC_mult[multBinXi]->SetBinContent(ptBinXi + 1, (resultParXiC_pt_mult[ptBinXi][multBinXi]->GetBinContent(1) * 0.97) / ((rawPt_xiC_mult[multBinXi]->GetBinWidth(ptBinXi + 1)) * h_multBinEntries_Xi->GetBinContent(multBinXi + 1))); // Bin 1 in resultparams is raw pt's bin counting
                rawPt_xiC_mult[multBinXi]->SetBinError(ptBinXi + 1, (resultParXiC_pt_mult[ptBinXi][multBinXi]->GetBinError(1) * 0.97) / ((rawPt_xiC_mult[multBinXi]->GetBinWidth(ptBinXi + 1)) * h_multBinEntries_Xi->GetBinContent(multBinXi + 1)));
                // rawPt_xip_mult[multBinXi]->SetBinContent(ptBinXi + 1, resultParXip_pt_mult[ptBinXi][multBinXi]->GetBinContent(1) / ((rawPt_xip_mult[multBinXi]->GetBinWidth(ptBinXi + 1)) * (h_multBinEntries_Xi->GetBinContent(multBinXi + 1)))); // Bin 1 in resultparams is raw pt's bin counting
                // rawPt_xip_mult[multBinXi]->SetBinError(ptBinXi + 1, resultParXip_pt_mult[ptBinXi][multBinXi]->GetBinError(1) / ((rawPt_xip_mult[multBinXi]->GetBinWidth(ptBinXi + 1)) * (h_multBinEntries_Xi->GetBinContent(multBinXi + 1))));
                // rawPt_xim_mult[multBinXi]->SetBinContent(ptBinXi + 1, resultParXim_pt_mult[ptBinXi][multBinXi]->GetBinContent(1) / ((rawPt_xim_mult[multBinXi]->GetBinWidth(ptBinXi + 1)) * (h_multBinEntries_Xi->GetBinContent(multBinXi + 1)))); // Bin 1 in resultparams is raw pt's bin counting
                // rawPt_xim_mult[multBinXi]->SetBinError(ptBinXi + 1, resultParXim_pt_mult[ptBinXi][multBinXi]->GetBinError(1) / ((rawPt_xim_mult[multBinXi]->GetBinWidth(ptBinXi + 1)) * (h_multBinEntries_Xi->GetBinContent(multBinXi + 1))));
                // rawPt_xiC_mult[multBinXi]->SetBinContent(ptBinXi + 1, resultParXiC_pt_mult[ptBinXi][multBinXi]->GetBinContent(1) / ((rawPt_xiC_mult[multBinXi]->GetBinWidth(ptBinXi + 1)) * (h_multBinEntries_Xi->GetBinContent(multBinXi + 1)))); // Bin 1 in resultparams is raw pt's bin counting
                // rawPt_xiC_mult[multBinXi]->SetBinError(ptBinXi + 1, resultParXiC_pt_mult[ptBinXi][multBinXi]->GetBinError(1) / ((rawPt_xiC_mult[multBinXi]->GetBinWidth(ptBinXi + 1)) * (h_multBinEntries_Xi->GetBinContent(multBinXi + 1))));
            }

            // Xi+: Mult diff: write to file
            outputFile->cd("dirRawPt_xip");
            rawPt_xip_mult[multBinXi]->SetMarkerStyle(markerStyles[multBinXi]);
            rawPt_xip_mult[multBinXi]->Write();
            outputFile->cd("dirEffCorrPt_xip");
            effCorrPt_xip_mult[multBinXi]->Divide(rawPt_xip_mult[multBinXi], eff_xiC_avg);
            effCorrPt_xip_mult[multBinXi]->SetMarkerStyle(markerStyles[multBinXi]);
            effCorrPt_xip_mult[multBinXi]->Write();
            // Xi+: Mult diff: add to stack
            rawPt_xip_mult[multBinXi]->SetTitle(TString::Format(("Mult: %.0f-%.0f%%, #times2^{%d}"), fMultbins_Xi[multBinXi], fMultbins_Xi[multBinXi + 1], (fNmultbins_Xi - 1) - multBinXi));
            rawPt_xip_mult[multBinXi]->Scale(pow(2, (fNmultbins_Xi - 1) - multBinXi));
            hs_xip_raw->Add(rawPt_xip_mult[multBinXi]);
            effCorrPt_xip_mult[multBinXi]->SetTitle(TString::Format(("Mult: %.0f-%.0f%% #times2^{%d}"), fMultbins_Xi[multBinXi], fMultbins_Xi[multBinXi + 1], (fNmultbins_Xi - 1) - multBinXi));
            effCorrPt_xip_mult[multBinXi]->Scale(pow(2, (fNmultbins_Xi - 1) - multBinXi));
            hs_xip_effCorr->Add(effCorrPt_xip_mult[multBinXi]);

            // Xi-: Mult diff: write to file
            outputFile->cd("dirRawPt_xim");
            rawPt_xim_mult[multBinXi]->SetMarkerStyle(markerStyles[multBinXi]);
            rawPt_xim_mult[multBinXi]->Write();
            outputFile->cd("dirEffCorrPt_xim");
            effCorrPt_xim_mult[multBinXi]->Divide(rawPt_xim_mult[multBinXi], eff_xiC_avg);
            effCorrPt_xim_mult[multBinXi]->SetMarkerStyle(markerStyles[multBinXi]);
            effCorrPt_xim_mult[multBinXi]->Write();
            // Xi-: Mult diff: add to stack
            rawPt_xim_mult[multBinXi]->SetTitle(TString::Format(("Mult: %.0f-%.0f%%, #times2^{%d}"), fMultbins_Xi[multBinXi], fMultbins_Xi[multBinXi + 1], (fNmultbins_Xi - 1) - multBinXi));
            rawPt_xim_mult[multBinXi]->Scale(pow(2, (fNmultbins_Xi - 1) - multBinXi));
            hs_xim_raw->Add(rawPt_xim_mult[multBinXi]);
            effCorrPt_xim_mult[multBinXi]->SetTitle(TString::Format(("Mult: %.0f-%.0f%% #times2^{%d}"), fMultbins_Xi[multBinXi], fMultbins_Xi[multBinXi + 1], (fNmultbins_Xi - 1) - multBinXi));
            effCorrPt_xim_mult[multBinXi]->Scale(pow(2, (fNmultbins_Xi - 1) - multBinXi));
            hs_xim_effCorr->Add(effCorrPt_xim_mult[multBinXi]);

            // XiC: Mult diff: write to file
            outputFile->cd("dirRawPt_xiC");
            rawPt_xiC_mult[multBinXi]->SetMarkerStyle(markerStyles[multBinXi]);
            rawPt_xiC_mult[multBinXi]->Write();
            outputFile->cd("dirEffCorrPt_xiC");
            effCorrPt_xiC_mult[multBinXi]->Divide(rawPt_xiC_mult[multBinXi], eff_xiC_avg);
            effCorrPt_xiC_mult[multBinXi]->SetMarkerStyle(markerStyles[multBinXi]);
            effCorrPt_xiC_mult[multBinXi]->Write();
            // XiC: Mult diff: add to stack
            rawPt_xiC_mult[multBinXi]->SetTitle(TString::Format(("Mult: %.0f-%.0f%%, #times2^{%d}"), fMultbins_Xi[multBinXi], fMultbins_Xi[multBinXi + 1], (fNmultbins_Xi - 1) - multBinXi));
            rawPt_xiC_mult[multBinXi]->Scale(pow(2, (fNmultbins_Xi - 1) - multBinXi));
            hs_xiC_raw->Add(rawPt_xiC_mult[multBinXi]);
            effCorrPt_xiC_mult[multBinXi]->SetTitle(TString::Format(("Mult: %.0f-%.0f%% #times2^{%d}"), fMultbins_Xi[multBinXi], fMultbins_Xi[multBinXi + 1], (fNmultbins_Xi - 1) - multBinXi));
            effCorrPt_xiC_mult[multBinXi]->Scale(pow(2, (fNmultbins_Xi - 1) - multBinXi));
            hs_xiC_effCorr->Add(effCorrPt_xiC_mult[multBinXi]);
        }
    } // End: XI

    /// OMEGA:
    {
        /// Mult integrated raw + efficiency corrected spectra:
        rawPt_omp = new TH1D("rawPt_omp", "Raw #Omega^{+}: Mult:0-100%", fNptbins_Om, fPtbins_Om);
        rawPt_omm = new TH1D("rawPt_omm", "Raw #Omega^{-}: Mult:0-100%", fNptbins_Om, fPtbins_Om);
        rawPt_omC = new TH1D("rawPt_omC", "Raw #Omega^{+} + #Omega^{-}: Mult:0-100%", fNptbins_Om, fPtbins_Om);

        effCorrPt_omp = new TH1D("effCorrPt_omp", "Eff. corrected #Omega^{+}: Mult:0-100%", fNptbins_Om, fPtbins_Om);
        effCorrPt_omm = new TH1D("effCorrPt_omm", "Eff. corrected #Omega^{-}: Mult:0-100%", fNptbins_Om, fPtbins_Om);
        effCorrPt_omC = new TH1D("effCorrPt_omC", "Eff. corrected #Omega^{+} + #Omega^{-}: Mult:0-100%", fNptbins_Om, fPtbins_Om);

        for (Int_t ptBinOm = 0; ptBinOm < fNptbins_Om; ptBinOm++)
        {
            // Trigger and primary vertex reconstruction efficiency = 0.97
            rawPt_omp->SetBinContent(ptBinOm + 1, (resultParOmp_pt[ptBinOm]->GetBinContent(1) * 0.97) / ((rawPt_omp->GetBinWidth(ptBinOm + 1)) * multEntriesOm)); // Bin 1 in resultparams is raw pt's bin counting
            rawPt_omp->SetBinError(ptBinOm + 1, (resultParOmp_pt[ptBinOm]->GetBinError(1) * 0.97) / ((rawPt_omp->GetBinWidth(ptBinOm + 1)) * multEntriesOm));
            rawPt_omm->SetBinContent(ptBinOm + 1, (resultParOmm_pt[ptBinOm]->GetBinContent(1) * 0.97) / ((rawPt_omm->GetBinWidth(ptBinOm + 1)) * multEntriesOm)); // Bin 1 in resultparams is raw pt's bin counting
            rawPt_omm->SetBinError(ptBinOm + 1, (resultParOmm_pt[ptBinOm]->GetBinError(1) * 0.97) / ((rawPt_omm->GetBinWidth(ptBinOm + 1)) * multEntriesOm));
            rawPt_omC->SetBinContent(ptBinOm + 1, (resultParOmC_pt[ptBinOm]->GetBinContent(1) * 0.97) / ((rawPt_omC->GetBinWidth(ptBinOm + 1)) * multEntriesOm)); // Bin 1 in resultparams is raw pt's bin counting
            rawPt_omC->SetBinError(ptBinOm + 1, (resultParOmC_pt[ptBinOm]->GetBinError(1) * 0.97) / ((rawPt_omC->GetBinWidth(ptBinOm + 1)) * multEntriesOm));
            // rawPt_omp->SetBinContent(ptBinOm + 1, resultParOmp_pt[ptBinOm]->GetBinContent(1) / ((rawPt_omp->GetBinWidth(ptBinOm + 1)) * multEntriesOm)); // Bin 1 in resultparams is raw pt's bin counting
            // rawPt_omp->SetBinError(ptBinOm + 1, resultParOmp_pt[ptBinOm]->GetBinError(1) / ((rawPt_omp->GetBinWidth(ptBinOm + 1)) * multEntriesOm));
            // rawPt_omm->SetBinContent(ptBinOm + 1, resultParOmm_pt[ptBinOm]->GetBinContent(1) / ((rawPt_omm->GetBinWidth(ptBinOm + 1)) * multEntriesOm)); // Bin 1 in resultparams is raw pt's bin counting
            // rawPt_omm->SetBinError(ptBinOm + 1, resultParOmm_pt[ptBinOm]->GetBinError(1) / ((rawPt_omm->GetBinWidth(ptBinOm + 1)) * multEntriesOm));
            // rawPt_omC->SetBinContent(ptBinOm + 1, resultParOmC_pt[ptBinOm]->GetBinContent(1) / ((rawPt_omC->GetBinWidth(ptBinOm + 1)) * multEntriesOm)); // Bin 1 in resultparams is raw pt's bin counting
            // rawPt_omC->SetBinError(ptBinOm + 1, resultParOmC_pt[ptBinOm]->GetBinError(1) / ((rawPt_omC->GetBinWidth(ptBinOm + 1)) * multEntriesOm));
        }

        // Om+: Mult:0-100%: write to file
        outputFile->cd("dirRawPt_omp");
        rawPt_omp->SetMarkerStyle(markerStyles[11]);
        rawPt_omp->Write();
        outputFile->cd("dirEffCorrPt_omp");
        effCorrPt_omp->Divide(rawPt_omp, eff_omC_avg);
        effCorrPt_omp->SetMarkerStyle(markerStyles[11]);
        effCorrPt_omp->Write();
        // // Om+: Mult:0-100%: add to stack
        // rawPt_omp->SetTitle("Mult:0-100%");
        // hs_omp_raw->Add(rawPt_omp);
        // effCorrPt_omp->SetTitle("Mult:0-100%");
        // hs_omp_effCorr->Add(effCorrPt_omp);

        // Om-: Mult:0-100%: write to file
        outputFile->cd("dirRawPt_omm");
        rawPt_omm->SetMarkerStyle(markerStyles[11]);
        rawPt_omm->Write();
        outputFile->cd("dirEffCorrPt_omm");
        effCorrPt_omm->Divide(rawPt_omm, eff_omC_avg);
        effCorrPt_omm->SetMarkerStyle(markerStyles[11]);
        effCorrPt_omm->Write();
        // // Om-: Mult:0-100%: add to stack
        // rawPt_omm->SetTitle("Mult:0-100%");
        // hs_omm_raw->Add(rawPt_omm);
        // effCorrPt_omm->SetTitle("Mult:0-100%");
        // hs_omm_effCorr->Add(effCorrPt_omm);

        // OmC: Mult:0-100%: write to file
        outputFile->cd("dirRawPt_omC");
        rawPt_omC->SetMarkerStyle(markerStyles[11]);
        rawPt_omC->Write();
        outputFile->cd("dirEffCorrPt_omC");
        effCorrPt_omC->Divide(rawPt_omC, eff_omC_avg);
        effCorrPt_omC->SetMarkerStyle(markerStyles[11]);
        effCorrPt_omC->Write();
        // // OmC: Mult:0-100%: add to stack
        // rawPt_omC->SetTitle("Mult:0-100%");
        // hs_omC_raw->Add(rawPt_omC);
        // effCorrPt_omC->SetTitle("Mult:0-100%");
        // hs_omC_effCorr->Add(effCorrPt_omC);

        /// Mult binned raw + efficiency corrected spectra:
        for (Int_t multBinOm = 0; multBinOm < fNmultbins_Om; multBinOm++)
        {
            // raw pt spectra: mult binned (diff)
            rawPt_omp_mult[multBinOm] = new TH1D(TString::Format(("rawPt_omp_mult[%d]"), multBinOm), TString::Format(("Raw #Omega^{+}: Mult: %.0f-%.0f%%"), fMultbins_Om[multBinOm], fMultbins_Om[multBinOm + 1]), fNptbins_Om, fPtbins_Om);
            rawPt_omm_mult[multBinOm] = new TH1D(TString::Format(("rawPt_omm_mult[%d]"), multBinOm), TString::Format(("Raw #Omega^{-}: Mult: %.0f-%.0f%%"), fMultbins_Om[multBinOm], fMultbins_Om[multBinOm + 1]), fNptbins_Om, fPtbins_Om);
            rawPt_omC_mult[multBinOm] = new TH1D(TString::Format(("rawPt_omC_mult[%d]"), multBinOm), TString::Format(("Raw #Omega^{+} + #Omega^{-}: Mult: %.0f-%.0f%%"), fMultbins_Om[multBinOm], fMultbins_Om[multBinOm + 1]), fNptbins_Om, fPtbins_Om);
            // effiency corrected pt spectra: mult binned (diff)
            effCorrPt_omp_mult[multBinOm] = new TH1D(TString::Format(("effCorrPt_omp_mult[%d]"), multBinOm), TString::Format(("Eff. corrected #Omega^{+}: Mult: %.0f-%.0f%%"), fMultbins_Om[multBinOm], fMultbins_Om[multBinOm + 1]), fNptbins_Om, fPtbins_Om);
            effCorrPt_omm_mult[multBinOm] = new TH1D(TString::Format(("effCorrPt_omm_mult[%d]"), multBinOm), TString::Format(("Eff. corrected #Omega^{-}: Mult: %.0f-%.0f%%"), fMultbins_Om[multBinOm], fMultbins_Om[multBinOm + 1]), fNptbins_Om, fPtbins_Om);
            effCorrPt_omC_mult[multBinOm] = new TH1D(TString::Format(("effCorrPt_omC_mult[%d]"), multBinOm), TString::Format(("Eff. corrected #Omega^{+} + #Omega^{-}: Mult: %.0f-%.0f%%"), fMultbins_Om[multBinOm], fMultbins_Om[multBinOm + 1]), fNptbins_Om, fPtbins_Om);

            for (Int_t ptBinOm = 0; ptBinOm < fNptbins_Om; ptBinOm++)
            {
                rawPt_omp_mult[multBinOm]->SetBinContent(ptBinOm + 1, (resultParOmp_pt_mult[ptBinOm][multBinOm]->GetBinContent(1) * 0.97) / ((rawPt_omp_mult[multBinOm]->GetBinWidth(ptBinOm + 1)) * h_multBinEntries_Om->GetBinContent(multBinOm + 1))); // Bin 1 in resultparams is raw pt's bin counting
                rawPt_omp_mult[multBinOm]->SetBinError(ptBinOm + 1, (resultParOmp_pt_mult[ptBinOm][multBinOm]->GetBinError(1) * 0.97) / ((rawPt_omp_mult[multBinOm]->GetBinWidth(ptBinOm + 1)) * h_multBinEntries_Om->GetBinContent(multBinOm + 1)));
                rawPt_omm_mult[multBinOm]->SetBinContent(ptBinOm + 1, (resultParOmm_pt_mult[ptBinOm][multBinOm]->GetBinContent(1) * 0.97) / ((rawPt_omm_mult[multBinOm]->GetBinWidth(ptBinOm + 1)) * h_multBinEntries_Om->GetBinContent(multBinOm + 1))); // Bin 1 in resultparams is raw pt's bin counting
                rawPt_omm_mult[multBinOm]->SetBinError(ptBinOm + 1, (resultParOmm_pt_mult[ptBinOm][multBinOm]->GetBinError(1) * 0.97) / ((rawPt_omm_mult[multBinOm]->GetBinWidth(ptBinOm + 1)) * h_multBinEntries_Om->GetBinContent(multBinOm + 1)));
                rawPt_omC_mult[multBinOm]->SetBinContent(ptBinOm + 1, (resultParOmC_pt_mult[ptBinOm][multBinOm]->GetBinContent(1) * 0.97) / ((rawPt_omC_mult[multBinOm]->GetBinWidth(ptBinOm + 1)) * h_multBinEntries_Om->GetBinContent(multBinOm + 1))); // Bin 1 in resultparams is raw pt's bin counting
                rawPt_omC_mult[multBinOm]->SetBinError(ptBinOm + 1, (resultParOmC_pt_mult[ptBinOm][multBinOm]->GetBinError(1) * 0.97) / ((rawPt_omC_mult[multBinOm]->GetBinWidth(ptBinOm + 1)) * h_multBinEntries_Om->GetBinContent(multBinOm + 1)));
                // rawPt_omp_mult[multBinOm]->SetBinContent(ptBinOm + 1, resultParOmp_pt_mult[ptBinOm][multBinOm]->GetBinContent(1) / ((rawPt_omp_mult[multBinOm]->GetBinWidth(ptBinOm + 1)) * (h_multBinEntries_Om->GetBinContent(multBinOm + 1)))); // Bin 1 in resultparams is raw pt's bin counting
                // rawPt_omp_mult[multBinOm]->SetBinError(ptBinOm + 1, resultParOmp_pt_mult[ptBinOm][multBinOm]->GetBinError(1) / ((rawPt_omp_mult[multBinOm]->GetBinWidth(ptBinOm + 1)) * (h_multBinEntries_Om->GetBinContent(multBinOm + 1))));
                // rawPt_omm_mult[multBinOm]->SetBinContent(ptBinOm + 1, resultParOmm_pt_mult[ptBinOm][multBinOm]->GetBinContent(1) / ((rawPt_omm_mult[multBinOm]->GetBinWidth(ptBinOm + 1)) * (h_multBinEntries_Om->GetBinContent(multBinOm + 1)))); // Bin 1 in resultparams is raw pt's bin counting
                // rawPt_omm_mult[multBinOm]->SetBinError(ptBinOm + 1, resultParOmm_pt_mult[ptBinOm][multBinOm]->GetBinError(1) / ((rawPt_omm_mult[multBinOm]->GetBinWidth(ptBinOm + 1)) * (h_multBinEntries_Om->GetBinContent(multBinOm + 1))));
                // rawPt_omC_mult[multBinOm]->SetBinContent(ptBinOm + 1, resultParOmC_pt_mult[ptBinOm][multBinOm]->GetBinContent(1) / ((rawPt_omC_mult[multBinOm]->GetBinWidth(ptBinOm + 1)) * (h_multBinEntries_Om->GetBinContent(multBinOm + 1)))); // Bin 1 in resultparams is raw pt's bin counting
                // rawPt_omC_mult[multBinOm]->SetBinError(ptBinOm + 1, resultParOmC_pt_mult[ptBinOm][multBinOm]->GetBinError(1) / ((rawPt_omC_mult[multBinOm]->GetBinWidth(ptBinOm + 1)) * (h_multBinEntries_Om->GetBinContent(multBinOm + 1))));
            }

            // Om+: Mult diff: write to file
            outputFile->cd("dirRawPt_omp");
            rawPt_omp_mult[multBinOm]->SetMarkerStyle(markerStyles[multBinOm]);
            rawPt_omp_mult[multBinOm]->Write();
            outputFile->cd("dirEffCorrPt_omp");
            effCorrPt_omp_mult[multBinOm]->Divide(rawPt_omp_mult[multBinOm], eff_omC_avg);
            effCorrPt_omp_mult[multBinOm]->SetMarkerStyle(markerStyles[multBinOm]);
            effCorrPt_omp_mult[multBinOm]->Write();
            // Om+: Mult diff: add to stack
            rawPt_omp_mult[multBinOm]->SetTitle(TString::Format(("Mult: %.0f-%.0f%%, #times2^{%d}"), fMultbins_Om[multBinOm], fMultbins_Om[multBinOm + 1], (fNmultbins_Om - 1) - multBinOm));
            rawPt_omp_mult[multBinOm]->Scale(pow(2, (fNmultbins_Om - 1) - multBinOm));
            hs_omp_raw->Add(rawPt_omp_mult[multBinOm]);
            effCorrPt_omp_mult[multBinOm]->SetTitle(TString::Format(("Mult: %.0f-%.0f%% #times2^{%d}"), fMultbins_Om[multBinOm], fMultbins_Om[multBinOm + 1], (fNmultbins_Om - 1) - multBinOm));
            effCorrPt_omp_mult[multBinOm]->Scale(pow(2, (fNmultbins_Om - 1) - multBinOm));
            hs_omp_effCorr->Add(effCorrPt_omp_mult[multBinOm]);

            // Om-: Mult diff: write to file
            outputFile->cd("dirRawPt_omm");
            rawPt_omm_mult[multBinOm]->SetMarkerStyle(markerStyles[multBinOm]);
            rawPt_omm_mult[multBinOm]->Write();
            outputFile->cd("dirEffCorrPt_omm");
            effCorrPt_omm_mult[multBinOm]->Divide(rawPt_omm_mult[multBinOm], eff_omC_avg);
            effCorrPt_omm_mult[multBinOm]->SetMarkerStyle(markerStyles[multBinOm]);
            effCorrPt_omm_mult[multBinOm]->Write();
            // Om-: Mult diff: add to stack
            rawPt_omm_mult[multBinOm]->SetTitle(TString::Format(("Mult: %.0f-%.0f%%, #times2^{%d}"), fMultbins_Om[multBinOm], fMultbins_Om[multBinOm + 1], (fNmultbins_Om - 1) - multBinOm));
            rawPt_omm_mult[multBinOm]->Scale(pow(2, (fNmultbins_Om - 1) - multBinOm));
            hs_omm_raw->Add(rawPt_omm_mult[multBinOm]);
            effCorrPt_omm_mult[multBinOm]->SetTitle(TString::Format(("Mult: %.0f-%.0f%% #times2^{%d}"), fMultbins_Om[multBinOm], fMultbins_Om[multBinOm + 1], (fNmultbins_Om - 1) - multBinOm));
            effCorrPt_omm_mult[multBinOm]->Scale(pow(2, (fNmultbins_Om - 1) - multBinOm));
            hs_omm_effCorr->Add(effCorrPt_omm_mult[multBinOm]);

            // OmC: Mult diff: write to file
            outputFile->cd("dirRawPt_omC");
            rawPt_omC_mult[multBinOm]->SetMarkerStyle(markerStyles[multBinOm]);
            rawPt_omC_mult[multBinOm]->Write();
            outputFile->cd("dirEffCorrPt_omC");
            effCorrPt_omC_mult[multBinOm]->Divide(rawPt_omC_mult[multBinOm], eff_omC_avg);
            effCorrPt_omC_mult[multBinOm]->SetMarkerStyle(markerStyles[multBinOm]);
            effCorrPt_omC_mult[multBinOm]->Write();
            // OmC: Mult diff: add to stack
            rawPt_omC_mult[multBinOm]->SetTitle(TString::Format(("Mult: %.0f-%.0f%%, #times2^{%d}"), fMultbins_Om[multBinOm], fMultbins_Om[multBinOm + 1], (fNmultbins_Om - 1) - multBinOm));
            rawPt_omC_mult[multBinOm]->Scale(pow(2, (fNmultbins_Om - 1) - multBinOm));
            hs_omC_raw->Add(rawPt_omC_mult[multBinOm]);
            effCorrPt_omC_mult[multBinOm]->SetTitle(TString::Format(("Mult: %.0f-%.0f%% #times2^{%d}"), fMultbins_Om[multBinOm], fMultbins_Om[multBinOm + 1], (fNmultbins_Om - 1) - multBinOm));
            effCorrPt_omC_mult[multBinOm]->Scale(pow(2, (fNmultbins_Om - 1) - multBinOm));
            hs_omC_effCorr->Add(effCorrPt_omC_mult[multBinOm]);
        }
    } // End: OMEGA

    /// Draw the histogram stacks
    gStyle->SetOptStat(0);
    gStyle->SetPalette(kVisibleSpectrum);
    TString histImageOutFolder = outputFolder + "/" + histName;
    TCanvas *cPtSpectrum[6];
    for (Int_t iCanvas = 0; iCanvas < 6; iCanvas++)
    {
        cPtSpectrum[iCanvas] = new TCanvas(TString::Format("cPtSpectrum%d", iCanvas), TString::Format("cPtSpectrum%d", iCanvas), 2560, 1440);
    }

    PaintStack(*cPtSpectrum[0], *hs_xip_raw);
    PaintStack(*cPtSpectrum[1], *hs_omp_raw);
    PaintStack(*cPtSpectrum[2], *hs_xim_raw);
    PaintStack(*cPtSpectrum[3], *hs_omm_raw);
    PaintStack(*cPtSpectrum[4], *hs_xiC_raw);
    PaintStack(*cPtSpectrum[5], *hs_omC_raw);

    if (saveStack)
    {
        cPtSpectrum[0]->cd();
        SaveImage(histImageOutFolder, "hs_rawPtSpectra", "hs_xip_raw", imageFormat, cPtSpectrum[0]);

        cPtSpectrum[1]->cd();
        SaveImage(histImageOutFolder, "hs_rawPtSpectra", "hs_omp_raw", imageFormat, cPtSpectrum[1]);

        cPtSpectrum[2]->cd();
        SaveImage(histImageOutFolder, "hs_rawPtSpectra", "hs_xim_raw", imageFormat, cPtSpectrum[2]);

        cPtSpectrum[3]->cd();
        SaveImage(histImageOutFolder, "hs_rawPtSpectra", "hs_omm_raw", imageFormat, cPtSpectrum[3]);

        cPtSpectrum[4]->cd();
        SaveImage(histImageOutFolder, "hs_rawPtSpectra", "hs_xiC_raw", imageFormat, cPtSpectrum[4]);

        cPtSpectrum[5]->cd();
        SaveImage(histImageOutFolder, "hs_rawPtSpectra", "hs_omC_raw", imageFormat, cPtSpectrum[5]);
    }

    PaintStack(*cPtSpectrum[0], *hs_xip_effCorr);
    PaintStack(*cPtSpectrum[1], *hs_omp_effCorr);
    PaintStack(*cPtSpectrum[2], *hs_xim_effCorr);
    PaintStack(*cPtSpectrum[3], *hs_omm_effCorr);
    PaintStack(*cPtSpectrum[4], *hs_xiC_effCorr);
    PaintStack(*cPtSpectrum[5], *hs_omC_effCorr);

    if (saveStack)
    {
        cPtSpectrum[0]->cd();
        SaveImage(histImageOutFolder, "hs_effCorrPtSpectra", "hs_xip_effCorr", imageFormat, cPtSpectrum[0]);

        cPtSpectrum[1]->cd();
        SaveImage(histImageOutFolder, "hs_effCorrPtSpectra", "hs_omp_effCorr", imageFormat, cPtSpectrum[1]);

        cPtSpectrum[2]->cd();
        SaveImage(histImageOutFolder, "hs_effCorrPtSpectra", "hs_xim_effCorr", imageFormat, cPtSpectrum[2]);

        cPtSpectrum[3]->cd();
        SaveImage(histImageOutFolder, "hs_effCorrPtSpectra", "hs_omm_effCorr", imageFormat, cPtSpectrum[3]);

        cPtSpectrum[4]->cd();
        SaveImage(histImageOutFolder, "hs_effCorrPtSpectra", "hs_xiC_effCorr", imageFormat, cPtSpectrum[4]);

        cPtSpectrum[5]->cd();
        SaveImage(histImageOutFolder, "hs_effCorrPtSpectra", "hs_omC_effCorr", imageFormat, cPtSpectrum[5]);
    }

    /// Write histogram stacks to file
    outputFile->cd("dirHistStacks");
    hs_xip_raw->Write();
    hs_xim_raw->Write();
    hs_xiC_raw->Write();
    hs_xip_effCorr->Write();
    hs_xim_effCorr->Write();
    hs_xiC_effCorr->Write();
    hs_omp_raw->Write();
    hs_omm_raw->Write();
    hs_omC_raw->Write();
    hs_omp_effCorr->Write();
    hs_omm_effCorr->Write();
    hs_omC_effCorr->Write();

    /// Clean up
    for (Int_t iCanvas = 0; iCanvas < 6; iCanvas++)
    {
        delete cPtSpectrum[iCanvas];
    }

    outputFile->Close();
    delete outputFile;

    return 0;
}

// End of EfficiencyCorrection function