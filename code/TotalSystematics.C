#include <TROOT.h>
#include <TStyle.h>
#include <TLegend.h>
#include <Math/MinimizerOptions.h>
#include <Math/IntegratorOptions.h>
#include <TLine.h>

#include "CascadeUtils.h"

enum sysTypeSigExt
{
    kFitFunction,
    kSignalRange,
    kBackgroundRange,
    kNumSysTypeSigExt
};

TString gSysTypeSigExt[kNumSysTypeSigExt] = {"FitFunction", "SignalExtractionRegion", "BackgroundFitRange"};

/**
 * @brief Styles a ROOT histogram with specified formatting options
 *
 * @param hist Pointer to the histogram to be styled
 * @param title Optional title for the histogram (default: empty string)
 * @param isHighlighted Optional flag to apply highlighting style (default: false)
 *
 * @details This function applies consistent styling to ROOT histograms for
 * visualization purposes. The styling includes setting colors, line styles,
 * markers, and other visual properties.
 */
void StyleHistogram(TH1 *hist, TString title = "", Bool_t isHighlighted = kFALSE);

/**
 * @brief Styles the provided histogram for systematic analysis.
 *
 * This function modifies the appearance of the input histogram based on
 * certain systematic considerations. It may adjust visual attributes such as
 * color, line width, or marker style, depending on the specified particle type
 * and multiplicity bin.
 *
 * @param hist Pointer to a TH1 histogram that needs styling.
 * @param particle Identifies the type of particle under consideration (-1 indicates a default setting).
 * @param multBin Specifies the multiplicity bin or category (-1 indicates a default or global setting).
 */
void StyleSystematics(TH1 *hist, Int_t particle = -1, Int_t multBin = -1);

/**
 * @brief Applies a custom style configuration to a given histogram.
 *
 * This function modifies the appearance of the provided histogram
 * by adjusting its visual style based on the specified particle
 * type and multiplicity bin. This includes settings such as
 * marker style, color, and other related aesthetic parameters.
 *
 * @param hist Pointer to a histogram object (TH1).
 * @param particle Optional parameter to specify which particle type is being plotted.
 * @param multBin Optional parameter to specify the multiplicity bin for the histogram.
 */
void StyleStatistics(TH1 *hist, Int_t particle = -1, Int_t multBin = -1);

/**
 * @brief Like PaintStack(), but modified for overlapping stat and sys corrected spectra,
 * with a new legend support.
 *
 * Draws a THStack on a specified canvas with optional log-scale on the y-axis.
 * This function overlays the histograms contained in the given THStack onto the
 * provided TCanvas. The user can specify whether to use a logarithmic scale for the
 * y-axis and customize the axis titles.
 *
 * @param c          Reference to the TCanvas on which the THStack is drawn.
 * @param hs         Reference to the THStack containing the histograms to draw.
 * @param setLogY    Boolean controlling whether to set the y-axis to log scale
 *                   (default is kTRUE).
 * @param yAxisTitle String specifying the y-axis title
 *                   (default is "#frac{1}{#it{N}_{inel}} #frac{d#it{N}}{d#it{p}_{T}}").
 * @param xAxisTitle String specifying the x-axis title
 *                   (default is "#it{p}_{T} (GeV/c)").
 */
void PaintStackOverlap(TCanvas &c, THStack &hs, Bool_t setLogY = kTRUE, TString yAxisTitle = "#frac{1}{#it{N}_{inel}} #frac{d#it{N}}{d#it{p}_{T}}", TString xAxisTitle = "#it{p}_{T} (GeV/c)");

/**
 * @brief Computes and combines total systematic uncertainties for Xi and Omega spectra analysis
 *
 * This function:
 * - Takes efficiency corrected spectra with statistical uncertainties as input
 * - Combines systematic uncertainties from:
 *   - Multi-trial analysis (track/topological cuts)
 *   - Signal extraction
 *   - Material budget (constant 4%)
 * - Handles both multiplicity-integrated and multiplicity-dependent results
 * - Creates visualization of uncertainties using THStack objects
 * - Produces output plots and ROOT files containing:
 *   - Individual systematic uncertainty contributions
 *   - Total combined systematic uncertainties
 *   - Spectra with both statistical and systematic uncertainties
 *
 * @param inputEffCorrDefault Path to input ROOT file containing efficiency corrected spectra
 * @param inputSysMultiTrial Path to input ROOT file containing multi-trial systematic uncertainties
 * @param inputSysSigExtraction Path to input ROOT file containing signal extraction systematic uncertainties
 * @param outputFolder Output folder path for plots and results
 * @param outputFileName Output ROOT file name
 * @param saveStack Whether to save THStack plots (default: true)
 * @param imageFormat Image format for saved plots (default: "png")
 * @param verbosity ROOT verbosity level
 *
 * @return 0 on success, error code otherwise
 *
 * The function handles:
 * - Xi minus/plus
 * - Omega minus/plus
 * - Combined charge (Xi/Omega)
 * Both for multiplicity integrated and multiplicity binned cases
 */
int TotalSystematics(
    TString inputEffCorrDefault = "/var/home/ishaan/Work/git/analysis/results/RandomVars/100225_EfficiencyCorrected_Vars/h3_ptmasscent_def/100225_effCorr_h3_ptmasscent_def.root",
    TString inputSysMultiTrial = "/var/home/ishaan/Work/git/analysis/results/RandomVars/240225_SysUncertainty_MultiTrial_212noRBErr/240225_sysUncertainty_noRB_multiTrial.root",
    TString inputSysSigExtractionPrefix = "/var/home/ishaan/Work/git/analysis/results/RandomVars/250225_Systematics_SigExt",
    TString inputSysSigExtractionFileName = "250225_SystematicUncertainty_SigExt.root",
    TString outputFolder = "/var/home/ishaan/Work/git/analysis/results/RandomVars/260225_SysUncertainty_Total",
    TString outputFileName = "/var/home/ishaan/Work/git/analysis/results/RandomVars/260225_SysUncertainty_Total/260225_sysUncertainty_Total.root",
    Bool_t saveStack = kTRUE,
    TString imageFormat = "png",
    Int_t verbosity = kInfo)
{
    ROOT::EnableImplicitMT();

    gStyle->SetPaintTextFormat("1.3f");
    gErrorIgnoreLevel = verbosity;

    // remove ownership of objects from file so we can delete the file ptr
    TH1::AddDirectory(kFALSE);
    TH1::SetDefaultSumw2(kTRUE);

    outputFolder = SetOutputFolder(outputFolder);

    /// Declare output objects:
    TH1D *stat_effCorrPt_xim;                     // mult integrated efficiency corrected spectra for default cuts: errors are statistical uncertainties
    TH1D *stat_effCorrPt_xip;                     // mult integrated efficiency corrected spectra for default cuts: errors are statistical uncertainties
    TH1D *stat_effCorrPt_omm;                     // mult integrated efficiency corrected spectra for default cuts: errors are statistical uncertainties
    TH1D *stat_effCorrPt_omp;                     // mult integrated efficiency corrected spectra for default cuts: errors are statistical uncertainties
    TH1D *stat_effCorrPt_xiC;                     // mult integrated efficiency corrected spectra for default cuts: errors are statistical uncertainties
    TH1D *stat_effCorrPt_omC;                     // mult integrated efficiency corrected spectra for default cuts: errors are statistical uncertainties
    TH1D *stat_effCorrPt_xim_mult[fNmultbins_Xi]; // mult binned efficiency corrected spectra for default cuts: errors are statistical uncertainties
    TH1D *stat_effCorrPt_xip_mult[fNmultbins_Xi]; // mult binned efficiency corrected spectra for default cuts: errors are statistical uncertainties
    TH1D *stat_effCorrPt_omm_mult[fNmultbins_Om]; // mult binned efficiency corrected spectra for default cuts: errors are statistical uncertainties
    TH1D *stat_effCorrPt_omp_mult[fNmultbins_Om]; // mult binned efficiency corrected spectra for default cuts: errors are statistical uncertainties
    TH1D *stat_effCorrPt_xiC_mult[fNmultbins_Xi]; // mult binned efficiency corrected spectra for default cuts: errors are statistical uncertainties
    TH1D *stat_effCorrPt_omC_mult[fNmultbins_Om]; // mult binned efficiency corrected spectra for default cuts: errors are statistical uncertainties

    TH1D *sys_effCorrPt_xim;                     // mult integrated efficiency corrected spectra for default cuts: errors are systematic uncertainties
    TH1D *sys_effCorrPt_xip;                     // mult integrated efficiency corrected spectra for default cuts: errors are systematic uncertainties
    TH1D *sys_effCorrPt_omm;                     // mult integrated efficiency corrected spectra for default cuts: errors are systematic uncertainties
    TH1D *sys_effCorrPt_omp;                     // mult integrated efficiency corrected spectra for default cuts: errors are systematic uncertainties
    TH1D *sys_effCorrPt_xiC;                     // mult integrated efficiency corrected spectra for default cuts: errors are systematic uncertainties
    TH1D *sys_effCorrPt_omC;                     // mult integrated efficiency corrected spectra for default cuts: errors are systematic uncertainties
    TH1D *sys_effCorrPt_xim_mult[fNmultbins_Xi]; // mult binned efficiency corrected spectra for default cuts: errors are systematic uncertainties
    TH1D *sys_effCorrPt_xip_mult[fNmultbins_Xi]; // mult binned efficiency corrected spectra for default cuts: errors are systematic uncertainties
    TH1D *sys_effCorrPt_omm_mult[fNmultbins_Om]; // mult binned efficiency corrected spectra for default cuts: errors are systematic uncertainties
    TH1D *sys_effCorrPt_omp_mult[fNmultbins_Om]; // mult binned efficiency corrected spectra for default cuts: errors are systematic uncertainties
    TH1D *sys_effCorrPt_xiC_mult[fNmultbins_Xi]; // mult binned efficiency corrected spectra for default cuts: errors are systematic uncertainties
    TH1D *sys_effCorrPt_omC_mult[fNmultbins_Om]; // mult binned efficiency corrected spectra for default cuts: errors are systematic uncertainties

    TH1D *sysMultiTrial_xim;                                           // mult integrated systematic uncertainty from MultiTrial
    TH1D *sysMultiTrial_xip;                                           // mult integrated systematic uncertainty from MultiTrial
    TH1D *sysMultiTrial_omm;                                           // mult integrated systematic uncertainty from MultiTrial
    TH1D *sysMultiTrial_omp;                                           // mult integrated systematic uncertainty from MultiTrial
    TH1D *sysMultiTrial_xiC;                                           // mult integrated systematic uncertainty from MultiTrial
    TH1D *sysMultiTrial_omC;                                           // mult integrated systematic uncertainty from MultiTrial
    TH1D *sysMultiTrial_xim_mult[fNmultbins_Xi];                       // mult binned systematic uncertainty from MultiTrial
    TH1D *sysMultiTrial_xip_mult[fNmultbins_Xi];                       // mult binned systematic uncertainty from MultiTrial
    TH1D *sysMultiTrial_omm_mult[fNmultbins_Om];                       // mult binned systematic uncertainty from MultiTrial
    TH1D *sysMultiTrial_omp_mult[fNmultbins_Om];                       // mult binned systematic uncertainty from MultiTrial
    TH1D *sysMultiTrial_xiC_mult[fNmultbins_Xi];                       // mult binned systematic uncertainty from MultiTrial
    TH1D *sysMultiTrial_omC_mult[fNmultbins_Om];                       // mult binned systematic uncertainty from MultiTrial
    TH1D *sysSigExtraction_xim[kNumSysTypeSigExt];                     // mult integrated systematic uncertainty from Signal Extraction
    TH1D *sysSigExtraction_xip[kNumSysTypeSigExt];                     // mult integrated systematic uncertainty from Signal Extraction
    TH1D *sysSigExtraction_omm[kNumSysTypeSigExt];                     // mult integrated systematic uncertainty from Signal Extraction
    TH1D *sysSigExtraction_omp[kNumSysTypeSigExt];                     // mult integrated systematic uncertainty from Signal Extraction
    TH1D *sysSigExtraction_xiC[kNumSysTypeSigExt];                     // mult integrated systematic uncertainty from Signal Extraction
    TH1D *sysSigExtraction_omC[kNumSysTypeSigExt];                     // mult integrated systematic uncertainty from Signal Extraction
    TH1D *sysSigExtraction_xim_mult[kNumSysTypeSigExt][fNmultbins_Xi]; // mult binned systematic uncertainty from Signal Extraction
    TH1D *sysSigExtraction_xip_mult[kNumSysTypeSigExt][fNmultbins_Xi]; // mult binned systematic uncertainty from Signal Extraction
    TH1D *sysSigExtraction_omm_mult[kNumSysTypeSigExt][fNmultbins_Om]; // mult binned systematic uncertainty from Signal Extraction
    TH1D *sysSigExtraction_omp_mult[kNumSysTypeSigExt][fNmultbins_Om]; // mult binned systematic uncertainty from Signal Extraction
    TH1D *sysSigExtraction_xiC_mult[kNumSysTypeSigExt][fNmultbins_Xi]; // mult binned systematic uncertainty from Signal Extraction
    TH1D *sysSigExtraction_omC_mult[kNumSysTypeSigExt][fNmultbins_Om]; // mult binned systematic uncertainty from Signal Extraction

    TH1D *sysMaterialBudget_xi; // material budget systematic uncertainty: constant at 4%
    TH1D *sysMaterialBudget_om; // material budget systematic uncertainty: constant at 4%

    TH1D *sysTotal_xim;                     // mult integrated total systematic uncertainty
    TH1D *sysTotal_xip;                     // mult integrated total systematic uncertainty
    TH1D *sysTotal_omm;                     // mult integrated total systematic uncertainty
    TH1D *sysTotal_omp;                     // mult integrated total systematic uncertainty
    TH1D *sysTotal_xiC;                     // mult integrated total systematic uncertainty
    TH1D *sysTotal_omC;                     // mult integrated total systematic uncertainty
    TH1D *sysTotal_xim_mult[fNmultbins_Xi]; // mult binned total systematic uncertainty
    TH1D *sysTotal_xip_mult[fNmultbins_Xi]; // mult binned total systematic uncertainty
    TH1D *sysTotal_omm_mult[fNmultbins_Om]; // mult binned total systematic uncertainty
    TH1D *sysTotal_omp_mult[fNmultbins_Om]; // mult binned total systematic uncertainty
    TH1D *sysTotal_xiC_mult[fNmultbins_Xi]; // mult binned total systematic uncertainty
    TH1D *sysTotal_omC_mult[fNmultbins_Om]; // mult binned total systematic uncertainty

    /// Create THStack objects for plotting mult integrated statistical, total systematic, multitrial, sig extraction, and material budget systematic uncertainties:
    THStack *hs_uncertainty_xim = new THStack("hs_uncertainty_xim", "#Xi^{-} uncertainties: Mult 0-100%");
    THStack *hs_uncertainty_xip = new THStack("hs_uncertainty_xip", "#Xi^{+} uncertainties: Mult 0-100%");
    THStack *hs_uncertainty_omm = new THStack("hs_uncertainty_omm", "#Omega^{-} uncertainties: Mult 0-100%");
    THStack *hs_uncertainty_omp = new THStack("hs_uncertainty_omp", "#Omega^{+} uncertainties: Mult 0-100%");
    THStack *hs_uncertainty_xiC = new THStack("hs_uncertainty_xiC", "#Xi^{+} + #Xi^{-} uncertainties: Mult 0-100%");
    THStack *hs_uncertainty_omC = new THStack("hs_uncertainty_omC", "#Omega^{+} + #Omega^{-} uncertainties: Mult 0-100%");

    /// Create THStack objects for plotting efficiency corrected spectra with statistical and systematic uncertainties:
    THStack *hs_effCorrUncert_xim = new THStack("hs_effCorrUncert_xim", "Corrected #it{p}_{T} spectra #Xi^{-}");
    THStack *hs_effCorrUncert_xip = new THStack("hs_effCorrUncert_xip", "Corrected #it{p}_{T} spectra #Xi^{+}");
    THStack *hs_effCorrUncert_omm = new THStack("hs_effCorrUncert_omm", "Corrected #it{p}_{T} spectra #Omega^{-}");
    THStack *hs_effCorrUncert_omp = new THStack("hs_effCorrUncert_omp", "Corrected #it{p}_{T} spectra #Omega^{+}");
    THStack *hs_effCorrUncert_xiC = new THStack("hs_effCorrUncert_xiC", "Corrected #it{p}_{T} spectra #Xi^{+} + #Xi^{-}");
    THStack *hs_effCorrUncert_omC = new THStack("hs_effCorrUncert_omC", "Corrected #it{p}_{T} spectra #Omega^{+} + #Omega^{-}");

    /// Create output file and directories:
    TFile *outputFile = new TFile(outputFileName, "RECREATE");
    Info("TotalSystematics: outputFile", "Writing output to '%s'", outputFileName.Data());
    outputFile->mkdir("dirEffCorr_Stat/xim");
    outputFile->mkdir("dirEffCorr_Stat/xip");
    outputFile->mkdir("dirEffCorr_Stat/omm");
    outputFile->mkdir("dirEffCorr_Stat/omp");
    outputFile->mkdir("dirEffCorr_Stat/xiC");
    outputFile->mkdir("dirEffCorr_Stat/omC");
    outputFile->mkdir("dirEffCorr_Sys/xim");
    outputFile->mkdir("dirEffCorr_Sys/xip");
    outputFile->mkdir("dirEffCorr_Sys/omm");
    outputFile->mkdir("dirEffCorr_Sys/omp");
    outputFile->mkdir("dirEffCorr_Sys/xiC");
    outputFile->mkdir("dirEffCorr_Sys/omC");
    outputFile->mkdir("dirSysMultiTrial/xim");
    outputFile->mkdir("dirSysMultiTrial/xip");
    outputFile->mkdir("dirSysMultiTrial/omm");
    outputFile->mkdir("dirSysMultiTrial/omp");
    outputFile->mkdir("dirSysMultiTrial/xiC");
    outputFile->mkdir("dirSysMultiTrial/omC");
    outputFile->mkdir("dirSysTotal/xim");
    outputFile->mkdir("dirSysTotal/xip");
    outputFile->mkdir("dirSysTotal/omm");
    outputFile->mkdir("dirSysTotal/omp");
    outputFile->mkdir("dirSysTotal/xiC");
    outputFile->mkdir("dirSysTotal/omC");
    outputFile->mkdir("dirStack_CompareSystematics");
    outputFile->mkdir("dirStack_EffCorrUncertainty");

    /// Getting default cut histograms:
    TFile *defInputFile = OpenFile(inputEffCorrDefault);
    Info("TotalSystematics: defInput", "Getting default cut histograms from '%s'", inputEffCorrDefault.Data());
    stat_effCorrPt_xim = (TH1D *)defInputFile->FindObjectAny("effCorrPt_xim");
    stat_effCorrPt_xip = (TH1D *)defInputFile->FindObjectAny("effCorrPt_xip");
    stat_effCorrPt_omm = (TH1D *)defInputFile->FindObjectAny("effCorrPt_omm");
    stat_effCorrPt_omp = (TH1D *)defInputFile->FindObjectAny("effCorrPt_omp");
    stat_effCorrPt_xiC = (TH1D *)defInputFile->FindObjectAny("effCorrPt_xiC");
    stat_effCorrPt_omC = (TH1D *)defInputFile->FindObjectAny("effCorrPt_omC");
    stat_effCorrPt_xim->SetNameTitle("stat_effCorrPt_xim", "#Xi^{-}: Eff. corrected #it{p}_{T} spectra w/ stat. uncertainties : Mult 0-100%;#it{p}_{T} (GeV/#it{c});#frac{1}{#it{N}_{inel}} #frac{d#it{N}}{d#it{p}_{T}}");
    stat_effCorrPt_xip->SetNameTitle("stat_effCorrPt_xip", "#Xi^{+}: Eff. corrected #it{p}_{T} spectra w/ stat. uncertainties : Mult 0-100%;#it{p}_{T} (GeV/#it{c});#frac{1}{#it{N}_{inel}} #frac{d#it{N}}{d#it{p}_{T}}");
    stat_effCorrPt_omm->SetNameTitle("stat_effCorrPt_omm", "#Omega^{-}: Eff. corrected #it{p}_{T} spectra w/ stat. uncertainties : Mult 0-100%;#it{p}_{T} (GeV/#it{c});#frac{1}{#it{N}_{inel}} #frac{d#it{N}}{d#it{p}_{T}}");
    stat_effCorrPt_omp->SetNameTitle("stat_effCorrPt_omp", "#Omega^{+}: Eff. corrected #it{p}_{T} spectra w/ stat. uncertainties : Mult 0-100%;#it{p}_{T} (GeV/#it{c});#frac{1}{#it{N}_{inel}} #frac{d#it{N}}{d#it{p}_{T}}");
    stat_effCorrPt_xiC->SetNameTitle("stat_effCorrPt_xiC", "#Xi^{+} + #Xi^{-}: Eff. corrected #it{p}_{T} spectra w/ stat. uncertainties : Mult 0-100%;#it{p}_{T} (GeV/#it{c});#frac{1}{#it{N}_{inel}} #frac{d#it{N}}{d#it{p}_{T}}");
    stat_effCorrPt_omC->SetNameTitle("stat_effCorrPt_omC", "#Omega^{+} + #Omega^{-}: Eff. corrected #it{p}_{T} spectra w/ stat. uncertainties : Mult 0-100%;#it{p}_{T} (GeV/#it{c});#frac{1}{#it{N}_{inel}} #frac{d#it{N}}{d#it{p}_{T}}");

    // Add to output file as no other changes are needed:
    outputFile->cd("dirEffCorr_Stat/xim");
    stat_effCorrPt_xim->Write();
    outputFile->cd("dirEffCorr_Stat/xip");
    stat_effCorrPt_xip->Write();
    outputFile->cd("dirEffCorr_Stat/omm");
    stat_effCorrPt_omm->Write();
    outputFile->cd("dirEffCorr_Stat/omp");
    stat_effCorrPt_omp->Write();
    outputFile->cd("dirEffCorr_Stat/xiC");
    stat_effCorrPt_xiC->Write();
    outputFile->cd("dirEffCorr_Stat/omC");
    stat_effCorrPt_omC->Write();
    for (Int_t multBinXi = 0; multBinXi < fNmultbins_Xi; multBinXi++)
    {
        stat_effCorrPt_xip_mult[multBinXi] = (TH1D *)defInputFile->FindObjectAny(TString::Format(("effCorrPt_xip_mult[%d]"), multBinXi));
        stat_effCorrPt_xim_mult[multBinXi] = (TH1D *)defInputFile->FindObjectAny(TString::Format(("effCorrPt_xim_mult[%d]"), multBinXi));
        stat_effCorrPt_xiC_mult[multBinXi] = (TH1D *)defInputFile->FindObjectAny(TString::Format(("effCorrPt_xiC_mult[%d]"), multBinXi));
        stat_effCorrPt_xip_mult[multBinXi]->SetNameTitle(TString::Format("stat_effCorrPt_xip_mult[%d]", multBinXi), TString::Format("#Xi^{+}: Eff. corrected #it{p}_{T} spectra w/ stat. uncertainties : Mult %.0f-%.0f%%;#it{p}_{T} (GeV/#it{c});#frac{1}{#it{N}_{inel}} #frac{d#it{N}}{d#it{p}_{T}}", fMultbins_Xi[multBinXi], fMultbins_Xi[multBinXi + 1]));
        stat_effCorrPt_xim_mult[multBinXi]->SetNameTitle(TString::Format("stat_effCorrPt_xim_mult[%d]", multBinXi), TString::Format("#Xi^{-}: Eff. corrected #it{p}_{T} spectra w/ stat. uncertainties : Mult %.0f-%.0f%%;#it{p}_{T} (GeV/#it{c});#frac{1}{#it{N}_{inel}} #frac{d#it{N}}{d#it{p}_{T}}", fMultbins_Xi[multBinXi], fMultbins_Xi[multBinXi + 1]));
        stat_effCorrPt_xiC_mult[multBinXi]->SetNameTitle(TString::Format("stat_effCorrPt_xiC_mult[%d]", multBinXi), TString::Format("#Xi^{+} + #Xi^{-}: Eff. corrected #it{p}_{T} spectra w/ stat. uncertainties : Mult %.0f-%.0f%%;#it{p}_{T} (GeV/#it{c});#frac{1}{#it{N}_{inel}} #frac{d#it{N}}{d#it{p}_{T}}", fMultbins_Xi[multBinXi], fMultbins_Xi[multBinXi + 1]));
        // Add to output file as no other changes are needed:
        outputFile->cd("dirEffCorr_Stat/xip");
        stat_effCorrPt_xip_mult[multBinXi]->Write();
        outputFile->cd("dirEffCorr_Stat/xim");
        stat_effCorrPt_xim_mult[multBinXi]->Write();
        outputFile->cd("dirEffCorr_Stat/xiC");
        stat_effCorrPt_xiC_mult[multBinXi]->Write();
    }
    for (Int_t multBinOm = 0; multBinOm < fNmultbins_Om; multBinOm++)
    {
        stat_effCorrPt_omp_mult[multBinOm] = (TH1D *)defInputFile->FindObjectAny(TString::Format(("effCorrPt_omp_mult[%d]"), multBinOm));
        stat_effCorrPt_omm_mult[multBinOm] = (TH1D *)defInputFile->FindObjectAny(TString::Format(("effCorrPt_omm_mult[%d]"), multBinOm));
        stat_effCorrPt_omC_mult[multBinOm] = (TH1D *)defInputFile->FindObjectAny(TString::Format(("effCorrPt_omC_mult[%d]"), multBinOm));
        stat_effCorrPt_omp_mult[multBinOm]->SetNameTitle(TString::Format("stat_effCorrPt_omp_mult[%d]", multBinOm), TString::Format("#Omega^{+}: Eff. corrected #it{p}_{T} spectra w/ stat. uncertainties : Mult %.0f-%.0f%%;#it{p}_{T} (GeV/#it{c});#frac{1}{#it{N}_{inel}} #frac{d#it{N}}{d#it{p}_{T}}", fMultbins_Om[multBinOm], fMultbins_Om[multBinOm + 1]));
        stat_effCorrPt_omm_mult[multBinOm]->SetNameTitle(TString::Format("stat_effCorrPt_omm_mult[%d]", multBinOm), TString::Format("#Omega^{-}: Eff. corrected #it{p}_{T} spectra w/ stat. uncertainties : Mult %.0f-%.0f%%;#it{p}_{T} (GeV/#it{c});#frac{1}{#it{N}_{inel}} #frac{d#it{N}}{d#it{p}_{T}}", fMultbins_Om[multBinOm], fMultbins_Om[multBinOm + 1]));
        stat_effCorrPt_omC_mult[multBinOm]->SetNameTitle(TString::Format("stat_effCorrPt_omC_mult[%d]", multBinOm), TString::Format("#Omega^{+} + #Omega^{-}: Eff. corrected #it{p}_{T} spectra w/ stat. uncertainties : Mult %.0f-%.0f%%;#it{p}_{T} (GeV/#it{c});#frac{1}{#it{N}_{inel}} #frac{d#it{N}}{d#it{p}_{T}}", fMultbins_Om[multBinOm], fMultbins_Om[multBinOm + 1]));
        // Add to output file as no other changes are needed:
        outputFile->cd("dirEffCorr_Stat/omp");
        stat_effCorrPt_omp_mult[multBinOm]->Write();
        outputFile->cd("dirEffCorr_Stat/omm");
        stat_effCorrPt_omm_mult[multBinOm]->Write();
        outputFile->cd("dirEffCorr_Stat/omC");
        stat_effCorrPt_omC_mult[multBinOm]->Write();
    }
    // input ended for defInputFile
    delete defInputFile;

    /// Generate effCorr histograms for holding total systematic uncertainties as error bars:
    sys_effCorrPt_xim = (TH1D *)stat_effCorrPt_xim->Clone();
    sys_effCorrPt_xip = (TH1D *)stat_effCorrPt_xip->Clone();
    sys_effCorrPt_omm = (TH1D *)stat_effCorrPt_omm->Clone();
    sys_effCorrPt_omp = (TH1D *)stat_effCorrPt_omp->Clone();
    sys_effCorrPt_xiC = (TH1D *)stat_effCorrPt_xiC->Clone();
    sys_effCorrPt_omC = (TH1D *)stat_effCorrPt_omC->Clone();
    sys_effCorrPt_xim->SetNameTitle("sys_effCorrPt_xim", "Eff. corrected #it{p}_{T} spectra w/ sys. uncertainties : Mult 0-100%;#it{p}_{T} (GeV/#it{c});#frac{1}{#it{N}_{inel}} #frac{d#it{N}}{d#it{p}_{T}}");
    sys_effCorrPt_xip->SetNameTitle("sys_effCorrPt_xip", "Eff. corrected #it{p}_{T} spectra w/ sys. uncertainties : Mult 0-100%;#it{p}_{T} (GeV/#it{c});#frac{1}{#it{N}_{inel}} #frac{d#it{N}}{d#it{p}_{T}}");
    sys_effCorrPt_omm->SetNameTitle("sys_effCorrPt_omm", "Eff. corrected #it{p}_{T} spectra w/ sys. uncertainties : Mult 0-100%;#it{p}_{T} (GeV/#it{c});#frac{1}{#it{N}_{inel}} #frac{d#it{N}}{d#it{p}_{T}}");
    sys_effCorrPt_omp->SetNameTitle("sys_effCorrPt_omp", "Eff. corrected #it{p}_{T} spectra w/ sys. uncertainties : Mult 0-100%;#it{p}_{T} (GeV/#it{c});#frac{1}{#it{N}_{inel}} #frac{d#it{N}}{d#it{p}_{T}}");
    sys_effCorrPt_xiC->SetNameTitle("sys_effCorrPt_xiC", "Eff. corrected #it{p}_{T} spectra w/ sys. uncertainties : Mult 0-100%;#it{p}_{T} (GeV/#it{c});#frac{1}{#it{N}_{inel}} #frac{d#it{N}}{d#it{p}_{T}}");
    sys_effCorrPt_omC->SetNameTitle("sys_effCorrPt_omC", "Eff. corrected #it{p}_{T} spectra w/ sys. uncertainties : Mult 0-100%;#it{p}_{T} (GeV/#it{c});#frac{1}{#it{N}_{inel}} #frac{d#it{N}}{d#it{p}_{T}}");

    for (Int_t multBinXi = 0; multBinXi < fNmultbins_Xi; multBinXi++)
    {
        sys_effCorrPt_xip_mult[multBinXi] = (TH1D *)stat_effCorrPt_xip_mult[multBinXi]->Clone(TString::Format("sys_effCorrPt_xip_mult[%d]", multBinXi));
        sys_effCorrPt_xim_mult[multBinXi] = (TH1D *)stat_effCorrPt_xim_mult[multBinXi]->Clone(TString::Format("sys_effCorrPt_xim_mult[%d]", multBinXi));
        sys_effCorrPt_xiC_mult[multBinXi] = (TH1D *)stat_effCorrPt_xiC_mult[multBinXi]->Clone(TString::Format("sys_effCorrPt_xiC_mult[%d]", multBinXi));
        sys_effCorrPt_xip_mult[multBinXi]->SetNameTitle(TString::Format("sys_effCorrPt_xip_mult[%d]", multBinXi), TString::Format("#Xi^{+}: Eff. corrected #it{p}_{T} spectra w/ sys. uncertainties : Mult %.0f-%.0f%%;#it{p}_{T} (GeV/#it{c});#frac{1}{#it{N}_{inel}} #frac{d#it{N}}{d#it{p}_{T}}", fMultbins_Xi[multBinXi], fMultbins_Xi[multBinXi + 1]));
        sys_effCorrPt_xim_mult[multBinXi]->SetNameTitle(TString::Format("sys_effCorrPt_xim_mult[%d]", multBinXi), TString::Format("#Xi^{-}: Eff. corrected #it{p}_{T} spectra w/ sys. uncertainties : Mult %.0f-%.0f%%;#it{p}_{T} (GeV/#it{c});#frac{1}{#it{N}_{inel}} #frac{d#it{N}}{d#it{p}_{T}}", fMultbins_Xi[multBinXi], fMultbins_Xi[multBinXi + 1]));
        sys_effCorrPt_xiC_mult[multBinXi]->SetNameTitle(TString::Format("sys_effCorrPt_xiC_mult[%d]", multBinXi), TString::Format("#Xi^{+} + #Xi^{-}: Eff. corrected #it{p}_{T} spectra w/ sys. uncertainties : Mult %.0f-%.0f%%;#it{p}_{T} (GeV/#it{c});#frac{1}{#it{N}_{inel}} #frac{d#it{N}}{d#it{p}_{T}}", fMultbins_Xi[multBinXi], fMultbins_Xi[multBinXi + 1]));
    }
    for (Int_t multBinOm = 0; multBinOm < fNmultbins_Om; multBinOm++)
    {
        sys_effCorrPt_omp_mult[multBinOm] = (TH1D *)stat_effCorrPt_omp_mult[multBinOm]->Clone(TString::Format("sys_effCorrPt_omp_mult[%d]", multBinOm));
        sys_effCorrPt_omm_mult[multBinOm] = (TH1D *)stat_effCorrPt_omm_mult[multBinOm]->Clone(TString::Format("sys_effCorrPt_omm_mult[%d]", multBinOm));
        sys_effCorrPt_omC_mult[multBinOm] = (TH1D *)stat_effCorrPt_omC_mult[multBinOm]->Clone(TString::Format("sys_effCorrPt_omC_mult[%d]", multBinOm));
        sys_effCorrPt_omp_mult[multBinOm]->SetNameTitle(TString::Format("sys_effCorrPt_omp_mult[%d]", multBinOm), TString::Format("#Omega^{+}: Eff. corrected #it{p}_{T} spectra w/ sys. uncertainties : Mult %.0f-%.0f%%;#it{p}_{T} (GeV/#it{c});#frac{1}{#it{N}_{inel}} #frac{d#it{N}}{d#it{p}_{T}}", fMultbins_Om[multBinOm], fMultbins_Om[multBinOm + 1]));
        sys_effCorrPt_omm_mult[multBinOm]->SetNameTitle(TString::Format("sys_effCorrPt_omm_mult[%d]", multBinOm), TString::Format("#Omega^{-}: Eff. corrected #it{p}_{T} spectra w/ sys. uncertainties : Mult %.0f-%.0f%%;#it{p}_{T} (GeV/#it{c});#frac{1}{#it{N}_{inel}} #frac{d#it{N}}{d#it{p}_{T}}", fMultbins_Om[multBinOm], fMultbins_Om[multBinOm + 1]));
        sys_effCorrPt_omC_mult[multBinOm]->SetNameTitle(TString::Format("sys_effCorrPt_omC_mult[%d]", multBinOm), TString::Format("#Omega^{+} + #Omega^{-}: Eff. corrected #it{p}_{T} spectra w/ sys. uncertainties : Mult %.0f-%.0f%%;#it{p}_{T} (GeV/#it{c});#frac{1}{#it{N}_{inel}} #frac{d#it{N}}{d#it{p}_{T}}", fMultbins_Om[multBinOm], fMultbins_Om[multBinOm + 1]));
    }

    /// Getting MultiTrial systematic uncertainty histograms:
    TFile *sysMultiTrialFile = OpenFile(inputSysMultiTrial);
    Info("TotalSystematics: sysMultiTrialInput", "Getting MultiTrial systematic uncertainty histograms from '%s'", inputSysMultiTrial.Data());
    sysMultiTrial_xim = (TH1D *)sysMultiTrialFile->FindObjectAny("sysMultiTrial_xim");
    sysMultiTrial_xip = (TH1D *)sysMultiTrialFile->FindObjectAny("sysMultiTrial_xip");
    sysMultiTrial_omm = (TH1D *)sysMultiTrialFile->FindObjectAny("sysMultiTrial_omm");
    sysMultiTrial_omp = (TH1D *)sysMultiTrialFile->FindObjectAny("sysMultiTrial_omp");
    sysMultiTrial_xiC = (TH1D *)sysMultiTrialFile->FindObjectAny("sysMultiTrial_xiC");
    sysMultiTrial_omC = (TH1D *)sysMultiTrialFile->FindObjectAny("sysMultiTrial_omC");
    outputFile->cd("dirSysMultiTrial/xim");
    sysMultiTrial_xim->Write();
    outputFile->cd("dirSysMultiTrial/xip");
    sysMultiTrial_xip->Write();
    outputFile->cd("dirSysMultiTrial/omm");
    sysMultiTrial_omm->Write();
    outputFile->cd("dirSysMultiTrial/omp");
    sysMultiTrial_omp->Write();
    outputFile->cd("dirSysMultiTrial/xiC");
    sysMultiTrial_xiC->Write();
    outputFile->cd("dirSysMultiTrial/omC");
    sysMultiTrial_omC->Write();

    StyleHistogram(sysMultiTrial_xim, "Multi-trial");
    StyleHistogram(sysMultiTrial_xip, "Multi-trial");
    StyleHistogram(sysMultiTrial_omm, "Multi-trial");
    StyleHistogram(sysMultiTrial_omp, "Multi-trial");
    StyleHistogram(sysMultiTrial_xiC, "Multi-trial");
    StyleHistogram(sysMultiTrial_omC, "Multi-trial");
    // Add MultiTrial systematic uncertainties to uncertainty stack (mult integrated):
    hs_uncertainty_xim->Add(sysMultiTrial_xim, "TEXT00");
    hs_uncertainty_xip->Add(sysMultiTrial_xip, "TEXT00");
    hs_uncertainty_omm->Add(sysMultiTrial_omm, "TEXT00");
    hs_uncertainty_omp->Add(sysMultiTrial_omp, "TEXT00");
    hs_uncertainty_xiC->Add(sysMultiTrial_xiC, "TEXT00");
    hs_uncertainty_omC->Add(sysMultiTrial_omC, "TEXT00");
    for (Int_t multBinXi = 0; multBinXi < fNmultbins_Xi; multBinXi++)
    {
        sysMultiTrial_xim_mult[multBinXi] = (TH1D *)sysMultiTrialFile->FindObjectAny(TString::Format("sysMultiTrial_xim_mult[%d]", multBinXi));
        sysMultiTrial_xip_mult[multBinXi] = (TH1D *)sysMultiTrialFile->FindObjectAny(TString::Format("sysMultiTrial_xip_mult[%d]", multBinXi));
        sysMultiTrial_xiC_mult[multBinXi] = (TH1D *)sysMultiTrialFile->FindObjectAny(TString::Format("sysMultiTrial_xiC_mult[%d]", multBinXi));

        outputFile->cd("dirSysMultiTrial/xim");
        sysMultiTrial_xim_mult[multBinXi]->Write();
        outputFile->cd("dirSysMultiTrial/xip");
        sysMultiTrial_xip_mult[multBinXi]->Write();
        outputFile->cd("dirSysMultiTrial/xiC");
        sysMultiTrial_xiC_mult[multBinXi]->Write();
    }
    for (Int_t multBinOm = 0; multBinOm < fNmultbins_Om; multBinOm++)
    {
        sysMultiTrial_omm_mult[multBinOm] = (TH1D *)sysMultiTrialFile->FindObjectAny(TString::Format("sysMultiTrial_omm_mult[%d]", multBinOm));
        sysMultiTrial_omp_mult[multBinOm] = (TH1D *)sysMultiTrialFile->FindObjectAny(TString::Format("sysMultiTrial_omp_mult[%d]", multBinOm));
        sysMultiTrial_omC_mult[multBinOm] = (TH1D *)sysMultiTrialFile->FindObjectAny(TString::Format("sysMultiTrial_omC_mult[%d]", multBinOm));

        outputFile->cd("dirSysMultiTrial/omm");
        sysMultiTrial_omm_mult[multBinOm]->Write();
        outputFile->cd("dirSysMultiTrial/omp");
        sysMultiTrial_omp_mult[multBinOm]->Write();
        outputFile->cd("dirSysMultiTrial/omC");
        sysMultiTrial_omC_mult[multBinOm]->Write();
    }
    // input ended for sysMultiTrialFile
    delete sysMultiTrialFile;

    /// Getting Signal Extraction systematic uncertainty histograms:
    for (Int_t iType = 0; iType < kNumSysTypeSigExt; iType++)
    {
        TString fullSigExtName = TString::Format("%s/%s/%s", inputSysSigExtractionPrefix.Data(), gSysTypeSigExt[iType].Data(), inputSysSigExtractionFileName.Data());

        TFile *sysSigExtractionFile = OpenFile(fullSigExtName);
        Info("TotalSystematics: sysSigExtractionInput", "Getting Signal Extraction systematic uncertainty histograms from '%s'", fullSigExtName.Data());
        sysSigExtraction_xim[iType] = (TH1D *)sysSigExtractionFile->FindObjectAny("sysSigExt_xim");
        sysSigExtraction_xip[iType] = (TH1D *)sysSigExtractionFile->FindObjectAny("sysSigExt_xip");
        sysSigExtraction_omm[iType] = (TH1D *)sysSigExtractionFile->FindObjectAny("sysSigExt_omm");
        sysSigExtraction_omp[iType] = (TH1D *)sysSigExtractionFile->FindObjectAny("sysSigExt_omp");
        sysSigExtraction_xiC[iType] = (TH1D *)sysSigExtractionFile->FindObjectAny("sysSigExt_xiC");
        sysSigExtraction_omC[iType] = (TH1D *)sysSigExtractionFile->FindObjectAny("sysSigExt_omC");

        // make output directories
        outputFile->mkdir(TString::Format("dirSysSigExtraction/%s/xim", gSysTypeSigExt[iType].Data()));
        outputFile->mkdir(TString::Format("dirSysSigExtraction/%s/xip", gSysTypeSigExt[iType].Data()));
        outputFile->mkdir(TString::Format("dirSysSigExtraction/%s/omm", gSysTypeSigExt[iType].Data()));
        outputFile->mkdir(TString::Format("dirSysSigExtraction/%s/omp", gSysTypeSigExt[iType].Data()));
        outputFile->mkdir(TString::Format("dirSysSigExtraction/%s/xiC", gSysTypeSigExt[iType].Data()));
        outputFile->mkdir(TString::Format("dirSysSigExtraction/%s/omC", gSysTypeSigExt[iType].Data()));

        /// Add to output file as no other changes are needed:
        outputFile->cd(TString::Format("dirSysSigExtraction/%s/xim", gSysTypeSigExt[iType].Data()));
        sysSigExtraction_xim[iType]->Write();
        outputFile->cd(TString::Format("dirSysSigExtraction/%s/xip", gSysTypeSigExt[iType].Data()));
        sysSigExtraction_xip[iType]->Write();
        outputFile->cd(TString::Format("dirSysSigExtraction/%s/omm", gSysTypeSigExt[iType].Data()));
        sysSigExtraction_omm[iType]->Write();
        outputFile->cd(TString::Format("dirSysSigExtraction/%s/omp", gSysTypeSigExt[iType].Data()));
        sysSigExtraction_omp[iType]->Write();
        outputFile->cd(TString::Format("dirSysSigExtraction/%s/xiC", gSysTypeSigExt[iType].Data()));
        sysSigExtraction_xiC[iType]->Write();
        outputFile->cd(TString::Format("dirSysSigExtraction/%s/omC", gSysTypeSigExt[iType].Data()));
        sysSigExtraction_omC[iType]->Write();

        StyleHistogram(sysSigExtraction_xim[iType], gSysTypeSigExt[iType]);
        StyleHistogram(sysSigExtraction_xip[iType], gSysTypeSigExt[iType]);
        StyleHistogram(sysSigExtraction_omm[iType], gSysTypeSigExt[iType]);
        StyleHistogram(sysSigExtraction_omp[iType], gSysTypeSigExt[iType]);
        StyleHistogram(sysSigExtraction_xiC[iType], gSysTypeSigExt[iType]);
        StyleHistogram(sysSigExtraction_omC[iType], gSysTypeSigExt[iType]);

        // Add Signal Extraction systematic uncertainties to uncertainty stack (mult integrated):
        hs_uncertainty_xim->Add(sysSigExtraction_xim[iType], "TEXT00");
        hs_uncertainty_xip->Add(sysSigExtraction_xip[iType], "TEXT00");
        hs_uncertainty_omm->Add(sysSigExtraction_omm[iType], "TEXT00");
        hs_uncertainty_omp->Add(sysSigExtraction_omp[iType], "TEXT00");
        hs_uncertainty_xiC->Add(sysSigExtraction_xiC[iType], "TEXT00");
        hs_uncertainty_omC->Add(sysSigExtraction_omC[iType], "TEXT00");
        for (Int_t multBinXi = 0; multBinXi < fNmultbins_Xi; multBinXi++)
        {
            sysSigExtraction_xim_mult[iType][multBinXi] = (TH1D *)sysSigExtractionFile->FindObjectAny(TString::Format("sysSigExt_xim_mult[%d]", multBinXi));
            sysSigExtraction_xip_mult[iType][multBinXi] = (TH1D *)sysSigExtractionFile->FindObjectAny(TString::Format("sysSigExt_xip_mult[%d]", multBinXi));
            sysSigExtraction_xiC_mult[iType][multBinXi] = (TH1D *)sysSigExtractionFile->FindObjectAny(TString::Format("sysSigExt_xiC_mult[%d]", multBinXi));

            outputFile->cd(TString::Format("dirSysSigExtraction/%s/xim", gSysTypeSigExt[iType].Data()));
            sysSigExtraction_xim_mult[iType][multBinXi]->Write();
            outputFile->cd(TString::Format("dirSysSigExtraction/%s/xip", gSysTypeSigExt[iType].Data()));
            sysSigExtraction_xip_mult[iType][multBinXi]->Write();
            outputFile->cd(TString::Format("dirSysSigExtraction/%s/xiC", gSysTypeSigExt[iType].Data()));
            sysSigExtraction_xiC_mult[iType][multBinXi]->Write();
        }
        for (Int_t multBinOm = 0; multBinOm < fNmultbins_Om; multBinOm++)
        {
            sysSigExtraction_omm_mult[iType][multBinOm] = (TH1D *)sysSigExtractionFile->FindObjectAny(TString::Format("sysSigExt_omm_mult[%d]", multBinOm));
            sysSigExtraction_omp_mult[iType][multBinOm] = (TH1D *)sysSigExtractionFile->FindObjectAny(TString::Format("sysSigExt_omp_mult[%d]", multBinOm));
            sysSigExtraction_omC_mult[iType][multBinOm] = (TH1D *)sysSigExtractionFile->FindObjectAny(TString::Format("sysSigExt_omC_mult[%d]", multBinOm));

            outputFile->cd(TString::Format("dirSysSigExtraction/%s/omm", gSysTypeSigExt[iType].Data()));
            sysSigExtraction_omm_mult[iType][multBinOm]->Write();
            outputFile->cd(TString::Format("dirSysSigExtraction/%s/omp", gSysTypeSigExt[iType].Data()));
            sysSigExtraction_omp_mult[iType][multBinOm]->Write();
            outputFile->cd(TString::Format("dirSysSigExtraction/%s/omC", gSysTypeSigExt[iType].Data()));
            sysSigExtraction_omC_mult[iType][multBinOm]->Write();
        }
        // input ended for sysSigExtractionFile
        delete sysSigExtractionFile;
    }

    /// Compute material budget systematic uncertainty (4% of yield - constant as a function of pT):
    sysMaterialBudget_xi = (TH1D *)sysSigExtraction_xiC[kSignalRange]->Clone("sysMaterialBudget_xi");
    sysMaterialBudget_om = (TH1D *)sysSigExtraction_omC[kSignalRange]->Clone("sysMaterialBudget_om");
    for (Int_t iBin = 1; iBin <= sysMaterialBudget_xi->GetNbinsX(); iBin++)
    {
        sysMaterialBudget_xi->SetBinContent(iBin, 0.04);
        sysMaterialBudget_xi->SetBinError(iBin, 0.0);
    }
    for (Int_t iBin = 1; iBin <= sysMaterialBudget_om->GetNbinsX(); iBin++)
    {
        sysMaterialBudget_om->SetBinContent(iBin, 0.04);
        sysMaterialBudget_om->SetBinError(iBin, 0.0);
    }
    StyleHistogram(sysMaterialBudget_xi, "Material Budget");
    StyleHistogram(sysMaterialBudget_om, "Material Budget");
    // Add Material Budget systematic uncertainties to uncertainty stack:
    hs_uncertainty_xim->Add(sysMaterialBudget_xi, "TEXT00");
    hs_uncertainty_xip->Add(sysMaterialBudget_xi, "TEXT00");
    hs_uncertainty_omm->Add(sysMaterialBudget_om, "TEXT00");
    hs_uncertainty_omp->Add(sysMaterialBudget_om, "TEXT00");
    hs_uncertainty_xiC->Add(sysMaterialBudget_xi, "TEXT00");
    hs_uncertainty_omC->Add(sysMaterialBudget_om, "TEXT00");

    /// Generate TH1D objects for holding total systematic uncertainties:
    sysTotal_xip = new TH1D("sysTotal_xip", "#Xi^{+}: Total systematic uncertainty: Mult 0-100%;#it{p}_{T} (GeV/#it{c});Uncertainty", fNptbins_Xi, fPtbins_Xi);
    sysTotal_xim = new TH1D("sysTotal_xim", "#Xi^{-}: Total systematic uncertainty: Mult 0-100%;#it{p}_{T} (GeV/#it{c});Uncertainty", fNptbins_Xi, fPtbins_Xi);
    sysTotal_omp = new TH1D("sysTotal_omp", "#Omega^{+}: Total systematic uncertainty: Mult 0-100%;#it{p}_{T} (GeV/#it{c});Uncertainty", fNptbins_Om, fPtbins_Om);
    sysTotal_omm = new TH1D("sysTotal_omm", "#Omega^{-}: Total systematic uncertainty: Mult 0-100%;#it{p}_{T} (GeV/#it{c});Uncertainty", fNptbins_Om, fPtbins_Om);
    sysTotal_xiC = new TH1D("sysTotal_xiC", "#Xi^{+} + #Xi^{-}: Total systematic uncertainty: Mult 0-100%;#it{p}_{T} (GeV/#it{c});Uncertainty", fNptbins_Xi, fPtbins_Xi);
    sysTotal_omC = new TH1D("sysTotal_omC", "#Omega^{+} + #Omega^{-}: Total systematic uncertainty: Mult 0-100%;#it{p}_{T} (GeV/#it{c});Uncertainty", fNptbins_Om, fPtbins_Om);
    for (Int_t multBinXi = 0; multBinXi < fNmultbins_Xi; multBinXi++)
    {
        sysTotal_xip_mult[multBinXi] = new TH1D(TString::Format("sysTotal_xip_mult[%d]", multBinXi), TString::Format("#Xi^{+}: Total systematic uncertainty: Mult %.0f-%.0f%%;#it{p}_{T} (GeV/#it{c});Uncertainty", fMultbins_Xi[multBinXi], fMultbins_Xi[multBinXi + 1]), fNptbins_Xi, fPtbins_Xi);
        sysTotal_xim_mult[multBinXi] = new TH1D(TString::Format("sysTotal_xim_mult[%d]", multBinXi), TString::Format("#Xi^{-}: Total systematic uncertainty: Mult %.0f-%.0f%%;#it{p}_{T} (GeV/#it{c});Uncertainty", fMultbins_Xi[multBinXi], fMultbins_Xi[multBinXi + 1]), fNptbins_Xi, fPtbins_Xi);
        sysTotal_xiC_mult[multBinXi] = new TH1D(TString::Format("sysTotal_xiC_mult[%d]", multBinXi), TString::Format("#Xi^{+} + #Xi^{-}: Total systematic uncertainty: Mult %.0f-%.0f%%;#it{p}_{T} (GeV/#it{c});Uncertainty", fMultbins_Xi[multBinXi], fMultbins_Xi[multBinXi + 1]), fNptbins_Xi, fPtbins_Xi);
    }
    for (Int_t multBinOm = 0; multBinOm < fNmultbins_Om; multBinOm++)
    {
        sysTotal_omp_mult[multBinOm] = new TH1D(TString::Format("sysTotal_omp_mult[%d]", multBinOm), TString::Format("#Omega^{+}: Total systematic uncertainty: Mult %.0f-%.0f%%;#it{p}_{T} (GeV/#it{c});Uncertainty", fMultbins_Om[multBinOm], fMultbins_Om[multBinOm + 1]), fNptbins_Om, fPtbins_Om);
        sysTotal_omm_mult[multBinOm] = new TH1D(TString::Format("sysTotal_omm_mult[%d]", multBinOm), TString::Format("#Omega^{-}: Total systematic uncertainty: Mult %.0f-%.0f%%;#it{p}_{T} (GeV/#it{c});Uncertainty", fMultbins_Om[multBinOm], fMultbins_Om[multBinOm + 1]), fNptbins_Om, fPtbins_Om);
        sysTotal_omC_mult[multBinOm] = new TH1D(TString::Format("sysTotal_omC_mult[%d]", multBinOm), TString::Format("#Omega^{+} + #Omega^{-}: Total systematic uncertainty: Mult %.0f-%.0f%%;#it{p}_{T} (GeV/#it{c});Uncertainty", fMultbins_Om[multBinOm], fMultbins_Om[multBinOm + 1]), fNptbins_Om, fPtbins_Om);
    }

    // extract statistical errors from stat_effCorrPt histograms to another histogram for plotting with systematic errors
    TH1D *stat_effCorrPt_xim_err = (TH1D *)stat_effCorrPt_xim->Clone("stat_effCorrPt_xim_err");
    TH1D *stat_effCorrPt_xip_err = (TH1D *)stat_effCorrPt_xip->Clone("stat_effCorrPt_xip_err");
    TH1D *stat_effCorrPt_omm_err = (TH1D *)stat_effCorrPt_omm->Clone("stat_effCorrPt_omm_err");
    TH1D *stat_effCorrPt_omp_err = (TH1D *)stat_effCorrPt_omp->Clone("stat_effCorrPt_omp_err");
    TH1D *stat_effCorrPt_xiC_err = (TH1D *)stat_effCorrPt_xiC->Clone("stat_effCorrPt_xiC_err");
    TH1D *stat_effCorrPt_omC_err = (TH1D *)stat_effCorrPt_omC->Clone("stat_effCorrPt_omC_err");
    stat_effCorrPt_xim_err->Reset();
    stat_effCorrPt_xip_err->Reset();
    stat_effCorrPt_omm_err->Reset();
    stat_effCorrPt_omp_err->Reset();
    stat_effCorrPt_xiC_err->Reset();
    stat_effCorrPt_omC_err->Reset();
    /// Compute total systematic uncertainty:
    Double_t val_sysTotal_xim, val_sysTotal_xip, val_sysTotal_omm, val_sysTotal_omp, val_sysTotal_xiC, val_sysTotal_omC; // temporary variables for storing total systematic uncertainty

    for (Int_t ptBinXi = 0; ptBinXi < fNptbins_Xi; ptBinXi++)
    {
        val_sysTotal_xip = TMath::Sqrt(TMath::Power(sysMultiTrial_xip->GetBinContent(ptBinXi + 1), 2) + TMath::Power(sysSigExtraction_xip[kFitFunction]->GetBinContent(ptBinXi + 1), 2) + TMath::Power(sysSigExtraction_xip[kSignalRange]->GetBinContent(ptBinXi + 1), 2) + TMath::Power(sysSigExtraction_xip[kBackgroundRange]->GetBinContent(ptBinXi + 1), 2) + TMath::Power(0.04, 2) /*const 4% material budget*/);
        val_sysTotal_xim = TMath::Sqrt(TMath::Power(sysMultiTrial_xim->GetBinContent(ptBinXi + 1), 2) + TMath::Power(sysSigExtraction_xim[kFitFunction]->GetBinContent(ptBinXi + 1), 2) + TMath::Power(sysSigExtraction_xim[kSignalRange]->GetBinContent(ptBinXi + 1), 2) + TMath::Power(sysSigExtraction_xim[kBackgroundRange]->GetBinContent(ptBinXi + 1), 2) + TMath::Power(0.04, 2) /*const 4% material budget*/);
        val_sysTotal_xiC = TMath::Sqrt(TMath::Power(sysMultiTrial_xiC->GetBinContent(ptBinXi + 1), 2) + TMath::Power(sysSigExtraction_xiC[kFitFunction]->GetBinContent(ptBinXi + 1), 2) + TMath::Power(sysSigExtraction_xiC[kSignalRange]->GetBinContent(ptBinXi + 1), 2) + TMath::Power(sysSigExtraction_xiC[kBackgroundRange]->GetBinContent(ptBinXi + 1), 2) + TMath::Power(0.04, 2) /*const 4% material budget*/);
        sysTotal_xip->SetBinContent(ptBinXi + 1, val_sysTotal_xip);
        sysTotal_xim->SetBinContent(ptBinXi + 1, val_sysTotal_xim);
        sysTotal_xiC->SetBinContent(ptBinXi + 1, val_sysTotal_xiC);

        val_sysTotal_xip = stat_effCorrPt_xip->GetBinContent(ptBinXi + 1) * val_sysTotal_xip; // multiply by yield to get absolute uncertainty
        val_sysTotal_xim = stat_effCorrPt_xim->GetBinContent(ptBinXi + 1) * val_sysTotal_xim; // multiply by yield to get absolute uncertainty
        val_sysTotal_xiC = stat_effCorrPt_xiC->GetBinContent(ptBinXi + 1) * val_sysTotal_xiC; // multiply by yield to get absolute uncertainty
        sys_effCorrPt_xim->SetBinError(ptBinXi + 1, val_sysTotal_xim);
        sys_effCorrPt_xip->SetBinError(ptBinXi + 1, val_sysTotal_xip);
        sys_effCorrPt_xiC->SetBinError(ptBinXi + 1, val_sysTotal_xiC);

        // Compute relative statistical uncertainty for comparison:
        stat_effCorrPt_xim_err->SetBinContent(ptBinXi + 1, stat_effCorrPt_xim->GetBinError(ptBinXi + 1) / stat_effCorrPt_xim->GetBinContent(ptBinXi + 1));
        stat_effCorrPt_xip_err->SetBinContent(ptBinXi + 1, stat_effCorrPt_xip->GetBinError(ptBinXi + 1) / stat_effCorrPt_xip->GetBinContent(ptBinXi + 1));
        stat_effCorrPt_xiC_err->SetBinContent(ptBinXi + 1, stat_effCorrPt_xiC->GetBinError(ptBinXi + 1) / stat_effCorrPt_xiC->GetBinContent(ptBinXi + 1));

        for (Int_t multBinXi = 0; multBinXi < fNmultbins_Xi; multBinXi++)
        {
            val_sysTotal_xip = TMath::Sqrt(TMath::Power(sysMultiTrial_xip_mult[multBinXi]->GetBinContent(ptBinXi + 1), 2) + TMath::Power(sysSigExtraction_xip_mult[kFitFunction][multBinXi]->GetBinContent(ptBinXi + 1), 2) + TMath::Power(sysSigExtraction_xip_mult[kSignalRange][multBinXi]->GetBinContent(ptBinXi + 1), 2) + TMath::Power(sysSigExtraction_xip_mult[kBackgroundRange][multBinXi]->GetBinContent(ptBinXi + 1), 2) + TMath::Power(0.04, 2) /*const 4% material budget*/);
            val_sysTotal_xim = TMath::Sqrt(TMath::Power(sysMultiTrial_xim_mult[multBinXi]->GetBinContent(ptBinXi + 1), 2) + TMath::Power(sysSigExtraction_xim_mult[kFitFunction][multBinXi]->GetBinContent(ptBinXi + 1), 2) + TMath::Power(sysSigExtraction_xim_mult[kSignalRange][multBinXi]->GetBinContent(ptBinXi + 1), 2) + TMath::Power(sysSigExtraction_xim_mult[kBackgroundRange][multBinXi]->GetBinContent(ptBinXi + 1), 2) + TMath::Power(0.04, 2) /*const 4% material budget*/);
            val_sysTotal_xiC = TMath::Sqrt(TMath::Power(sysMultiTrial_xiC_mult[multBinXi]->GetBinContent(ptBinXi + 1), 2) + TMath::Power(sysSigExtraction_xiC_mult[kFitFunction][multBinXi]->GetBinContent(ptBinXi + 1), 2) + TMath::Power(sysSigExtraction_xiC_mult[kSignalRange][multBinXi]->GetBinContent(ptBinXi + 1), 2) + TMath::Power(sysSigExtraction_xiC_mult[kBackgroundRange][multBinXi]->GetBinContent(ptBinXi + 1), 2) + TMath::Power(0.04, 2) /*const 4% material budget*/);
            sysTotal_xip_mult[multBinXi]->SetBinContent(ptBinXi + 1, val_sysTotal_xip);
            sysTotal_xim_mult[multBinXi]->SetBinContent(ptBinXi + 1, val_sysTotal_xim);
            sysTotal_xiC_mult[multBinXi]->SetBinContent(ptBinXi + 1, val_sysTotal_xiC);

            val_sysTotal_xip = stat_effCorrPt_xip_mult[multBinXi]->GetBinContent(ptBinXi + 1) * val_sysTotal_xip; // multiply by yield to get absolute uncertainty
            val_sysTotal_xim = stat_effCorrPt_xim_mult[multBinXi]->GetBinContent(ptBinXi + 1) * val_sysTotal_xim; // multiply by yield to get absolute uncertainty
            val_sysTotal_xiC = stat_effCorrPt_xiC_mult[multBinXi]->GetBinContent(ptBinXi + 1) * val_sysTotal_xiC; // multiply by yield to get absolute uncertainty

            sys_effCorrPt_xim_mult[multBinXi]->SetBinError(ptBinXi + 1, val_sysTotal_xim);
            sys_effCorrPt_xip_mult[multBinXi]->SetBinError(ptBinXi + 1, val_sysTotal_xip);
            sys_effCorrPt_xiC_mult[multBinXi]->SetBinError(ptBinXi + 1, val_sysTotal_xiC);
        }
    }

    for (Int_t ptBinOm = 0; ptBinOm < fNptbins_Om; ptBinOm++)
    {
        val_sysTotal_omp = TMath::Sqrt(TMath::Power(sysMultiTrial_omp->GetBinContent(ptBinOm + 1), 2) + TMath::Power(sysSigExtraction_omp[kFitFunction]->GetBinContent(ptBinOm + 1), 2) + TMath::Power(sysSigExtraction_omp[kSignalRange]->GetBinContent(ptBinOm + 1), 2) + TMath::Power(sysSigExtraction_omp[kBackgroundRange]->GetBinContent(ptBinOm + 1), 2) + TMath::Power(0.04, 2) /*const 4% material budget*/);
        val_sysTotal_omm = TMath::Sqrt(TMath::Power(sysMultiTrial_omm->GetBinContent(ptBinOm + 1), 2) + TMath::Power(sysSigExtraction_omm[kFitFunction]->GetBinContent(ptBinOm + 1), 2) + TMath::Power(sysSigExtraction_omm[kSignalRange]->GetBinContent(ptBinOm + 1), 2) + TMath::Power(sysSigExtraction_omm[kBackgroundRange]->GetBinContent(ptBinOm + 1), 2) + TMath::Power(0.04, 2) /*const 4% material budget*/);
        val_sysTotal_omC = TMath::Sqrt(TMath::Power(sysMultiTrial_omC->GetBinContent(ptBinOm + 1), 2) + TMath::Power(sysSigExtraction_omC[kFitFunction]->GetBinContent(ptBinOm + 1), 2) + TMath::Power(sysSigExtraction_omC[kSignalRange]->GetBinContent(ptBinOm + 1), 2) + TMath::Power(sysSigExtraction_omC[kBackgroundRange]->GetBinContent(ptBinOm + 1), 2) + TMath::Power(0.04, 2) /*const 4% material budget*/);
        sysTotal_omp->SetBinContent(ptBinOm + 1, val_sysTotal_omp);
        sysTotal_omm->SetBinContent(ptBinOm + 1, val_sysTotal_omm);
        sysTotal_omC->SetBinContent(ptBinOm + 1, val_sysTotal_omC);

        val_sysTotal_omp = stat_effCorrPt_omp->GetBinContent(ptBinOm + 1) * val_sysTotal_omp; // multiply by yield to get absolute uncertainty
        val_sysTotal_omm = stat_effCorrPt_omm->GetBinContent(ptBinOm + 1) * val_sysTotal_omm; // multiply by yield to get absolute uncertainty
        val_sysTotal_omC = stat_effCorrPt_omC->GetBinContent(ptBinOm + 1) * val_sysTotal_omC; // multiply by yield to get absolute uncertainty
        sys_effCorrPt_omp->SetBinError(ptBinOm + 1, val_sysTotal_omp);
        sys_effCorrPt_omm->SetBinError(ptBinOm + 1, val_sysTotal_omm);
        sys_effCorrPt_omC->SetBinError(ptBinOm + 1, val_sysTotal_omC);

        // Compute relative statistical uncertainty for comparison:
        stat_effCorrPt_omm_err->SetBinContent(ptBinOm + 1, stat_effCorrPt_omm->GetBinError(ptBinOm + 1) / stat_effCorrPt_omm->GetBinContent(ptBinOm + 1));
        stat_effCorrPt_omp_err->SetBinContent(ptBinOm + 1, stat_effCorrPt_omp->GetBinError(ptBinOm + 1) / stat_effCorrPt_omp->GetBinContent(ptBinOm + 1));
        stat_effCorrPt_omC_err->SetBinContent(ptBinOm + 1, stat_effCorrPt_omC->GetBinError(ptBinOm + 1) / stat_effCorrPt_omC->GetBinContent(ptBinOm + 1));

        for (Int_t multBinOm = 0; multBinOm < fNmultbins_Om; multBinOm++)
        {
            val_sysTotal_omp = TMath::Sqrt(TMath::Power(sysMultiTrial_omp_mult[multBinOm]->GetBinContent(ptBinOm + 1), 2) + TMath::Power(sysSigExtraction_omp_mult[kFitFunction][multBinOm]->GetBinContent(ptBinOm + 1), 2) + TMath::Power(sysSigExtraction_omp_mult[kSignalRange][multBinOm]->GetBinContent(ptBinOm + 1), 2) + TMath::Power(sysSigExtraction_omp_mult[kBackgroundRange][multBinOm]->GetBinContent(ptBinOm + 1), 2) + TMath::Power(0.04, 2) /*const 4% material budget*/);
            val_sysTotal_omm = TMath::Sqrt(TMath::Power(sysMultiTrial_omm_mult[multBinOm]->GetBinContent(ptBinOm + 1), 2) + TMath::Power(sysSigExtraction_omm_mult[kFitFunction][multBinOm]->GetBinContent(ptBinOm + 1), 2) + TMath::Power(sysSigExtraction_omm_mult[kSignalRange][multBinOm]->GetBinContent(ptBinOm + 1), 2) + TMath::Power(sysSigExtraction_omm_mult[kBackgroundRange][multBinOm]->GetBinContent(ptBinOm + 1), 2) + TMath::Power(0.04, 2) /*const 4% material budget*/);
            val_sysTotal_omC = TMath::Sqrt(TMath::Power(sysMultiTrial_omC_mult[multBinOm]->GetBinContent(ptBinOm + 1), 2) + TMath::Power(sysSigExtraction_omC_mult[kFitFunction][multBinOm]->GetBinContent(ptBinOm + 1), 2) + TMath::Power(sysSigExtraction_omC_mult[kSignalRange][multBinOm]->GetBinContent(ptBinOm + 1), 2) + TMath::Power(sysSigExtraction_omC_mult[kBackgroundRange][multBinOm]->GetBinContent(ptBinOm + 1), 2) + TMath::Power(0.04, 2) /*const 4% material budget*/);
            sysTotal_omp_mult[multBinOm]->SetBinContent(ptBinOm + 1, val_sysTotal_omp);
            sysTotal_omm_mult[multBinOm]->SetBinContent(ptBinOm + 1, val_sysTotal_omm);
            sysTotal_omC_mult[multBinOm]->SetBinContent(ptBinOm + 1, val_sysTotal_omC);

            val_sysTotal_omp = stat_effCorrPt_omp_mult[multBinOm]->GetBinContent(ptBinOm + 1) * val_sysTotal_omp; // multiply by yield to get absolute uncertainty
            val_sysTotal_omm = stat_effCorrPt_omm_mult[multBinOm]->GetBinContent(ptBinOm + 1) * val_sysTotal_omm; // multiply by yield to get absolute uncertainty
            val_sysTotal_omC = stat_effCorrPt_omC_mult[multBinOm]->GetBinContent(ptBinOm + 1) * val_sysTotal_omC; // multiply by yield to get absolute uncertainty
            sys_effCorrPt_omp_mult[multBinOm]->SetBinError(ptBinOm + 1, val_sysTotal_omp);
            sys_effCorrPt_omm_mult[multBinOm]->SetBinError(ptBinOm + 1, val_sysTotal_omm);
            sys_effCorrPt_omC_mult[multBinOm]->SetBinError(ptBinOm + 1, val_sysTotal_omC);
        }
    }
    /// Write total systematic uncertainty histograms to output file before styling:
    outputFile->cd("dirSysTotal/xim");
    sysTotal_xim->Write();
    outputFile->cd("dirSysTotal/xip");
    sysTotal_xip->Write();
    outputFile->cd("dirSysTotal/omm");
    sysTotal_omm->Write();
    outputFile->cd("dirSysTotal/omp");
    sysTotal_omp->Write();
    outputFile->cd("dirSysTotal/xiC");
    sysTotal_xiC->Write();
    outputFile->cd("dirSysTotal/omC");
    sysTotal_omC->Write();
    for (Int_t multBinXi = 0; multBinXi < fNmultbins_Xi; multBinXi++)
    {
        outputFile->cd("dirSysTotal/xim");
        sysTotal_xim_mult[multBinXi]->Write();
        outputFile->cd("dirSysTotal/xip");
        sysTotal_xip_mult[multBinXi]->Write();
        outputFile->cd("dirSysTotal/xiC");
        sysTotal_xiC_mult[multBinXi]->Write();
    }
    for (Int_t multBinOm = 0; multBinOm < fNmultbins_Om; multBinOm++)
    {
        outputFile->cd("dirSysTotal/omm");
        sysTotal_omm_mult[multBinOm]->Write();
        outputFile->cd("dirSysTotal/omp");
        sysTotal_omp_mult[multBinOm]->Write();
        outputFile->cd("dirSysTotal/omC");
        sysTotal_omC_mult[multBinOm]->Write();
    }

    // Style histograms before adding to uncertainty stack:
    StyleHistogram(stat_effCorrPt_xim_err, "Statistical");
    StyleHistogram(stat_effCorrPt_xip_err, "Statistical");
    StyleHistogram(stat_effCorrPt_omm_err, "Statistical");
    StyleHistogram(stat_effCorrPt_omp_err, "Statistical");
    StyleHistogram(stat_effCorrPt_xiC_err, "Statistical");
    StyleHistogram(stat_effCorrPt_omC_err, "Statistical");
    StyleHistogram(sysTotal_xim, "Total Systematic", kTRUE);
    StyleHistogram(sysTotal_xip, "Total Systematic", kTRUE);
    StyleHistogram(sysTotal_omm, "Total Systematic", kTRUE);
    StyleHistogram(sysTotal_omp, "Total Systematic", kTRUE);
    StyleHistogram(sysTotal_xiC, "Total Systematic", kTRUE);
    StyleHistogram(sysTotal_omC, "Total Systematic", kTRUE);
    // Add stat and total systematic uncertainties to uncertainty stack (mult integrated):
    hs_uncertainty_xim->Add(stat_effCorrPt_xim_err, "TEXT00");
    hs_uncertainty_xip->Add(stat_effCorrPt_xip_err, "TEXT00");
    hs_uncertainty_omm->Add(stat_effCorrPt_omm_err, "TEXT00");
    hs_uncertainty_omp->Add(stat_effCorrPt_omp_err, "TEXT00");
    hs_uncertainty_xiC->Add(stat_effCorrPt_xiC_err, "TEXT00");
    hs_uncertainty_omC->Add(stat_effCorrPt_omC_err, "TEXT00");
    hs_uncertainty_xim->Add(sysTotal_xim, "TEXT00");
    hs_uncertainty_xip->Add(sysTotal_xip, "TEXT00");
    hs_uncertainty_omm->Add(sysTotal_omm, "TEXT00");
    hs_uncertainty_omp->Add(sysTotal_omp, "TEXT00");
    hs_uncertainty_xiC->Add(sysTotal_xiC, "TEXT00");
    hs_uncertainty_omC->Add(sysTotal_omC, "TEXT00");

    /// Write effCorr histograms with total systematic uncertainties to output file before styling:
    outputFile->cd("dirEffCorr_Sys/xip");
    sys_effCorrPt_xip->Write();
    outputFile->cd("dirEffCorr_Sys/xim");
    sys_effCorrPt_xim->Write();
    outputFile->cd("dirEffCorr_Sys/omp");
    sys_effCorrPt_omp->Write();
    outputFile->cd("dirEffCorr_Sys/omm");
    sys_effCorrPt_omm->Write();
    outputFile->cd("dirEffCorr_Sys/xiC");
    sys_effCorrPt_xiC->Write();
    outputFile->cd("dirEffCorr_Sys/omC");
    sys_effCorrPt_omC->Write();
    for (Int_t multBinXi = 0; multBinXi < fNmultbins_Xi; multBinXi++)
    {
        outputFile->cd("dirEffCorr_Sys/xip");
        sys_effCorrPt_xip_mult[multBinXi]->Write();
        outputFile->cd("dirEffCorr_Sys/xim");
        sys_effCorrPt_xim_mult[multBinXi]->Write();
        outputFile->cd("dirEffCorr_Sys/xiC");
        sys_effCorrPt_xiC_mult[multBinXi]->Write();
    }
    for (Int_t multBinOm = 0; multBinOm < fNmultbins_Om; multBinOm++)
    {
        outputFile->cd("dirEffCorr_Sys/omp");
        sys_effCorrPt_omp_mult[multBinOm]->Write();
        outputFile->cd("dirEffCorr_Sys/omm");
        sys_effCorrPt_omm_mult[multBinOm]->Write();
        outputFile->cd("dirEffCorr_Sys/omC");
        sys_effCorrPt_omC_mult[multBinOm]->Write();
    }
    // Style and add mult integrated effCorr histograms (sys and stat) to the stack:
    // Set style for systematic uncertainty boxes:
    // StyleSystematics(sys_effCorrPt_xip, kXi, -1);
    // StyleSystematics(sys_effCorrPt_xim, kXi, -1);
    // StyleSystematics(sys_effCorrPt_omp, kXi, -1);
    // StyleSystematics(sys_effCorrPt_omm, kXi, -1);
    // StyleSystematics(sys_effCorrPt_xiC, kXi, -1);
    // StyleSystematics(sys_effCorrPt_omC, kXi, -1);

    // // Set style for statistical error bars:
    // StyleStatistics(stat_effCorrPt_xip, kXi, -1);
    // StyleStatistics(stat_effCorrPt_xim, kXi, -1);
    // StyleStatistics(stat_effCorrPt_omp, kXi, -1);
    // StyleStatistics(stat_effCorrPt_omm, kXi, -1);
    // StyleStatistics(stat_effCorrPt_xiC, kXi, -1);
    // StyleStatistics(stat_effCorrPt_omC, kXi, -1);
    // hs_effCorrUncert_xim->Add(sysTotal_xim, "E2");
    // hs_effCorrUncert_xim->Add(stat_effCorrPt_xim, "P E1");
    // hs_effCorrUncert_xip->Add(sysTotal_xip, "E2");
    // hs_effCorrUncert_xip->Add(stat_effCorrPt_xip, "P E1");
    // hs_effCorrUncert_omm->Add(sysTotal_omm, "E2");
    // hs_effCorrUncert_omm->Add(stat_effCorrPt_omm, "P E1");
    // hs_effCorrUncert_omp->Add(sysTotal_omp, "E2");
    // hs_effCorrUncert_omp->Add(stat_effCorrPt_omp, "P E1");
    // hs_effCorrUncert_xiC->Add(sysTotal_xiC, "E2");
    // hs_effCorrUncert_xiC->Add(stat_effCorrPt_xiC, "P E1");
    // hs_effCorrUncert_omC->Add(sysTotal_omC, "E2");
    // hs_effCorrUncert_omC->Add(stat_effCorrPt_omC, "P E1");
    // Add mult binned effCorr histograms (sys and stat) to the stack:
    for (Int_t multBinXi = 0; multBinXi < fNmultbins_Xi; multBinXi++)
    {
        // Set style for systematic uncertainty boxes:
        StyleSystematics(sys_effCorrPt_xip_mult[multBinXi], kXi, multBinXi);
        StyleSystematics(sys_effCorrPt_xim_mult[multBinXi], kXi, multBinXi);
        StyleSystematics(sys_effCorrPt_xiC_mult[multBinXi], kXi, multBinXi);
        // Set style for statistical error bars:
        StyleStatistics(stat_effCorrPt_xip_mult[multBinXi], kXi, multBinXi);
        StyleStatistics(stat_effCorrPt_xim_mult[multBinXi], kXi, multBinXi);
        StyleStatistics(stat_effCorrPt_xiC_mult[multBinXi], kXi, multBinXi);
        hs_effCorrUncert_xim->Add(sys_effCorrPt_xim_mult[multBinXi], "E2");
        hs_effCorrUncert_xim->Add(stat_effCorrPt_xim_mult[multBinXi], "P E1");
        hs_effCorrUncert_xip->Add(sys_effCorrPt_xip_mult[multBinXi], "E2");
        hs_effCorrUncert_xip->Add(stat_effCorrPt_xip_mult[multBinXi], "P E1");
        hs_effCorrUncert_xiC->Add(sys_effCorrPt_xiC_mult[multBinXi], "E2");
        hs_effCorrUncert_xiC->Add(stat_effCorrPt_xiC_mult[multBinXi], "P E1");
    }
    for (Int_t multBinOm = 0; multBinOm < fNmultbins_Om; multBinOm++)
    {
        // Set style for systematic uncertainty boxes:
        StyleSystematics(sys_effCorrPt_omp_mult[multBinOm], kOm, multBinOm);
        StyleSystematics(sys_effCorrPt_omm_mult[multBinOm], kOm, multBinOm);
        StyleSystematics(sys_effCorrPt_omC_mult[multBinOm], kOm, multBinOm);
        // Set style for statistical error bars:
        StyleStatistics(stat_effCorrPt_omp_mult[multBinOm], kOm, multBinOm);
        StyleStatistics(stat_effCorrPt_omm_mult[multBinOm], kOm, multBinOm);
        StyleStatistics(stat_effCorrPt_omC_mult[multBinOm], kOm, multBinOm);
        hs_effCorrUncert_omm->Add(sys_effCorrPt_omm_mult[multBinOm], "E2");
        hs_effCorrUncert_omm->Add(stat_effCorrPt_omm_mult[multBinOm], "P E1");
        hs_effCorrUncert_omp->Add(sys_effCorrPt_omp_mult[multBinOm], "E2");
        hs_effCorrUncert_omp->Add(stat_effCorrPt_omp_mult[multBinOm], "P E1");
        hs_effCorrUncert_omC->Add(sys_effCorrPt_omC_mult[multBinOm], "E2");
        hs_effCorrUncert_omC->Add(stat_effCorrPt_omC_mult[multBinOm], "P E1");
    }

    // Draw and save stack plots:
    gStyle->SetPalette(kVisibleSpectrum);

    TCanvas *cSpectra[6];
    for (Int_t iCanvas = 0; iCanvas < 6; iCanvas++)
    {
        cSpectra[iCanvas] = new TCanvas(TString::Format("cSpectra%d", iCanvas), TString::Format("cSpectra%d", iCanvas), 3840, 2160);
    }

    PaintStack(*cSpectra[0], *hs_uncertainty_xim, kFALSE, "Relative Uncertainty");
    PaintStack(*cSpectra[1], *hs_uncertainty_xip, kFALSE, "Relative Uncertainty");
    PaintStack(*cSpectra[2], *hs_uncertainty_omm, kFALSE, "Relative Uncertainty");
    PaintStack(*cSpectra[3], *hs_uncertainty_omp, kFALSE, "Relative Uncertainty");
    PaintStack(*cSpectra[4], *hs_uncertainty_xiC, kFALSE, "Relative Uncertainty");
    PaintStack(*cSpectra[5], *hs_uncertainty_omC, kFALSE, "Relative Uncertainty");

    if (saveStack)
    {
        cSpectra[0]->cd();
        SaveImage(outputFolder, "hs_compareSystematics", "hs_uncertainty_xim", imageFormat.Data(), cSpectra[0]);

        cSpectra[1]->cd();
        SaveImage(outputFolder, "hs_compareSystematics", "hs_uncertainty_xip", imageFormat.Data(), cSpectra[1]);

        cSpectra[2]->cd();
        SaveImage(outputFolder, "hs_compareSystematics", "hs_uncertainty_omm", imageFormat.Data(), cSpectra[2]);

        cSpectra[3]->cd();
        SaveImage(outputFolder, "hs_compareSystematics", "hs_uncertainty_omp", imageFormat.Data(), cSpectra[3]);

        cSpectra[4]->cd();
        SaveImage(outputFolder, "hs_compareSystematics", "hs_uncertainty_xiC", imageFormat.Data(), cSpectra[4]);

        cSpectra[5]->cd();
        SaveImage(outputFolder, "hs_compareSystematics", "hs_uncertainty_omC", imageFormat.Data(), cSpectra[5]);
    }

    PaintStackOverlap(*cSpectra[0], *hs_effCorrUncert_xim);
    PaintStackOverlap(*cSpectra[1], *hs_effCorrUncert_xip);
    PaintStackOverlap(*cSpectra[2], *hs_effCorrUncert_omm);
    PaintStackOverlap(*cSpectra[3], *hs_effCorrUncert_omp);
    PaintStackOverlap(*cSpectra[4], *hs_effCorrUncert_xiC);
    PaintStackOverlap(*cSpectra[5], *hs_effCorrUncert_omC);

    if (saveStack)
    {
        cSpectra[0]->cd();
        SaveImage(outputFolder, "hs_effCorrUncert", "hs_effCorrUncert_xim", imageFormat.Data(), cSpectra[0]);

        cSpectra[1]->cd();
        SaveImage(outputFolder, "hs_effCorrUncert", "hs_effCorrUncert_xip", imageFormat.Data(), cSpectra[1]);

        cSpectra[2]->cd();
        SaveImage(outputFolder, "hs_effCorrUncert", "hs_effCorrUncert_omm", imageFormat.Data(), cSpectra[2]);

        cSpectra[3]->cd();
        SaveImage(outputFolder, "hs_effCorrUncert", "hs_effCorrUncert_omp", imageFormat.Data(), cSpectra[3]);

        cSpectra[4]->cd();
        SaveImage(outputFolder, "hs_effCorrUncert", "hs_effCorrUncert_xiC", imageFormat.Data(), cSpectra[4]);

        cSpectra[5]->cd();
        SaveImage(outputFolder, "hs_effCorrUncert", "hs_effCorrUncert_omC", imageFormat.Data(), cSpectra[5]);
    }

    // Write histogram stacks to output file:
    outputFile->cd("dirStack_CompareSystematics");
    hs_uncertainty_xim->Write();
    hs_uncertainty_xip->Write();
    hs_uncertainty_omm->Write();
    hs_uncertainty_omp->Write();
    hs_uncertainty_xiC->Write();
    hs_uncertainty_omC->Write();

    outputFile->cd("dirStack_EffCorrUncertainty");
    hs_effCorrUncert_xim->Write();
    hs_effCorrUncert_xip->Write();
    hs_effCorrUncert_omm->Write();
    hs_effCorrUncert_omp->Write();
    hs_effCorrUncert_xiC->Write();
    hs_effCorrUncert_omC->Write();

    delete outputFile;

    Info("TotalSystematics", "Output written to '%s'", outputFileName.Data());
    return 0;
}

void StyleHistogram(TH1 *hist, TString title, Bool_t isHighlighted)
{

    if (title != "")
    {
        hist->SetTitle(title);
    }
    hist->SetMarkerStyle(kFullCircle);
    // hist->SetMarkerSize(1.5);
    // hist->SetMarkerColor(kBlack);
    // hist->SetLineColor(kBlack);
    if (isHighlighted)
    {
        hist->SetLineStyle(kSolid);
        hist->SetLineWidth(4);
    }
    else
    {
        hist->SetLineWidth(2);
        if (title.Contains("Stat"))
        {
            hist->SetLineStyle(9);
            hist->SetLineColor(kRed + 4);
        }
        else
        {
            hist->SetLineStyle(kDashed);
        }
    }
    hist->SetStats(kFALSE);
}

void StyleSystematics(TH1 *hist, Int_t particle, Int_t multBin)
{
    // Set style for systematic uncertainty boxes:
    // hist->SetFillColor(kGray + 2);
    // hist->SetFillStyle(3004); // a cross-hatched fill
    // hist->SetLineColor(kGray + 2);
    hist->SetFillStyle(0);
    hist->SetLineWidth(1);
    hist->SetDrawOption("E2");
    hist->SetMarkerSize(0);
    hist->SetStats(kFALSE);
    if (multBin >= 0)
    {
        if (particle == kXi)
        {
            hist->SetTitle(TString::Format(("Mult: %.0f-%.0f%% (#times2^{%d})"), fMultbins_Xi[multBin], fMultbins_Xi[multBin + 1], (fNmultbins_Xi - 1) - multBin));
            hist->Scale(pow(2, (fNmultbins_Xi - 1) - multBin));
        }
        else if (particle == kOm)
        {
            hist->SetTitle(TString::Format(("Mult: %.0f-%.0f%% (#times2^{%d})"), fMultbins_Om[multBin], fMultbins_Om[multBin + 1], (fNmultbins_Om - 1) - multBin));
            hist->Scale(pow(2, (fNmultbins_Om - 1) - multBin));
        }
    }
    else // integrated multiplicity case
    {
        hist->SetTitle("Mult: 0-100%");
        hist->Scale(pow(2, -1));
    }
}

void StyleStatistics(TH1 *hist, Int_t particle, Int_t multBin)
{
    // Set style for statistical error bars:
    hist->SetMarkerSize(2);

    hist->SetLineWidth(3);
    hist->SetDrawOption("P E1 SAME");
    // hist->SetStats(kFALSE);
    if (multBin >= 0)
    {
        hist->SetMarkerStyle(markerStyles[multBin]);
        if (particle == kXi)
        {
            hist->SetTitle(TString::Format(("Mult: %.0f-%.0f%% (#times2^{%d})"), fMultbins_Xi[multBin], fMultbins_Xi[multBin + 1], (fNmultbins_Xi - 1) - multBin));
            hist->Scale(pow(2, (fNmultbins_Xi - 1) - multBin));
        }
        else if (particle == kOm)
        {
            hist->SetTitle(TString::Format(("Mult: %.0f-%.0f%% (#times2^{%d})"), fMultbins_Om[multBin], fMultbins_Om[multBin + 1], (fNmultbins_Om - 1) - multBin));
            hist->Scale(pow(2, (fNmultbins_Om - 1) - multBin));
        }
    }
    else // integrated multiplicity case
    {
        hist->SetMarkerStyle(kOpenSquareDiagonal);
        hist->SetTitle("Mult: 0-100%");
        hist->Scale(pow(2, -1));
    }
}

void PaintStackOverlap(TCanvas &c, THStack &hs, Bool_t setLogY, TString yAxisTitle, TString xAxisTitle)
{
    c.Clear();
    c.cd();
    hs.Draw("plc pmc nostack");
    if (setLogY)
        gPad->SetLogy();

    hs.GetXaxis()->SetTitle(xAxisTitle.Data());
    hs.GetYaxis()->SetTitle(yAxisTitle.Data());

    TLegend *legend = new TLegend(0.9, 0.6, 1., 1., "");
    TList *histList = hs.GetHists();
    for (int iHist = 0; iHist < histList->GetSize(); iHist++)
    {
        TH1 *hist = (TH1 *)histList->At(iHist);
        if (TString(hist->GetName()).Contains("stat"))

        {
            legend->AddEntry(hist, hist->GetTitle(), "lpef");
        }
    }
    legend->Draw();
    c.Modified();
    c.ForceUpdate();
}
