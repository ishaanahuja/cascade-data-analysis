#include <TROOT.h>
#include <TStyle.h>
#include <TLegend.h>
#include <Math/MinimizerOptions.h>
#include <Math/IntegratorOptions.h>

#include "CascadeUtils.h"

/**
 * @brief Computes the Roger Barlow criterion for determining significant deviations
 *
 * @param hVar Histogram containing the variation to be tested
 * @param hDef Histogram containing the default/reference values
 * @param hYieldDev Histogram to store yield deviations
 * @param ptBin The pT bin number being analyzed
 * @param iVar The variation number/index
 * @param fDebug Flag to enable debug output
 * @param nSigma Number of sigma for significance criterion (default = 1)
 *
 * @return Double_t Returns the computed Roger Barlow significance value
 *
 * This function implements the Roger Barlow criterion to determine if variations
 * from a default value are statistically significant. It compares histograms
 * of variations against a reference and calculates the statistical significance
 * of any observed differences.
 */
Double_t ComputeRogerBarlow(TH1 *hVar, TH1 *hDef, TH1 *hYieldDev, Int_t ptBin, Int_t iVar, Bool_t fDebug, UInt_t nSigma = 1);

/**
 * @brief Performs a Gaussian fit on a histogram
 *
 * @param hist Pointer to the histogram to be fitted
 * @param fitOptions Options string for the fit procedure (default: "BWLINQ+ MULTITHREAD")
 *        - B: Use better initial values from a Gaussian histogram fit
 *        - W: Set all weights to 1; ignore histogram errors
 *        - L: Use log likelihood method
 *        - I: Use integral of function instead of value in bin center
 *        - N: Do not store the fit results
 *        - Q: Quiet mode (suppress printouts)
 *        - +: Add fit results to existing ones
 *        - MULTITHREAD: Enable multithreading for the fit
 *
 * @return Double_t The sigma value of the fit
 */
Double_t FitGaus(TH1 *hist, TString fitOptions = "BWLINQ+ MULTITHREAD");

/**
 * @brief Draws a histogram and saves it as an image file
 *
 * @param hist Pointer to the histogram to be drawn
 * @param outputFolder Path to the output folder where the image will be saved
 * @param imageFolder Subdirectory name within the output folder
 * @param imageFormat File format for the saved image (e.g., "png", "pdf")
 *
 * @note The function creates the output directory structure if it doesn't exist
 */
void DrawAndSaveImage(TH1 *hist, TString outputFolder, TString imageFolder, TString imageFormat = "png", Int_t failedRbCount = 0);

/**
 * @brief Calculates systematic uncertainties using the multi-trial method for Xi and Omega particles
 *
 * This function processes multiple variations of cut selections to estimate systematic uncertainties
 * in particle yield measurements. It handles both Xi and Omega particles (positive, negative, and combined charges)
 * across different multiplicity and pt bins.
 *
 * @param inputPath Path to the directory containing efficiency corrected input files [default: "/var/home/ishaan/Work/git/analysis/results/RandomVars/100225_effCorr"]
 * @param effCorrInputFilePrefix Prefix for efficiency correction input files [default: "100225_effCorr"]
 * @param outputFileName Path to output ROOT file [default: "/var/home/ishaan/Work/git/analysis/results/RandomVars/120225_sysUncertainty/120225_sysUncertainty_multiTrial.root"]
 * @param outputFolder Path to output directory for results [default: "/var/home/ishaan/Work/git/analysis/results/RandomVars/120225_sysUncertainty"]
 * @param fDebug Enable debug output [default: kFALSE]
 * @param saveImages Save plots as image files [default: kTRUE]
 * @param imageFormat Format for saved images [default: "png"]
 * @param verbosity Level of ROOT verbosity [default: kInfo]
 *
 * @return int Returns 0 on successful completion
 *
 * The function performs the following main operations:
 * 1. Loads default cut histograms for reference
 * 2. Creates output objects for storing results
 * 3. Processes multiple variations (500) of cut selections
 * 4. Calculates yield deviations from default cuts
 * 5. Fits yield deviation distributions with Gaussian functions
 * 6. Extracts systematic uncertainties from fit parameters
 * 7. Saves results to ROOT file and generates plots
 *
 * Supports analysis for:
 * - Xi-minus, Xi-plus, and combined Xi
 * - Omega-minus, Omega-plus, and combined Omega
 * - Both multiplicity integrated and multiplicity dependent results
 */
int SysMultiTrial(
    TString inputPath = "/var/home/ishaan/Work/git/analysis/results/RandomVars/100225_EfficiencyCorrected_Vars",
    TString effCorrInputFilePrefix = "100225_effCorr",
    TString outputFileName = "/var/home/ishaan/Work/git/thesis/final/images/120625_SysUncertainty_MultiTrial_212noRBErr/120625_sysUncertainty_noRB_multiTrial.root",
    TString outputFolder = "/var/home/ishaan/Work/git/thesis/final/images/120625_SysUncertainty_MultiTrial_212noRBErr",
    Bool_t fDebug = kFALSE,
    Bool_t saveImages = kTRUE,
    TString imageFormat = "pdf",
    Int_t verbosity = kInfo)
{
    ROOT::EnableImplicitMT();
    ROOT::Math::MinimizerOptions::SetDefaultMaxFunctionCalls(100000);
    ROOT::Math::MinimizerOptions::SetDefaultMinimizer("Minuit2");
    ROOT::Math::MinimizerOptions::SetDefaultStrategy(2);
    ROOT::Math::MinimizerOptions::SetDefaultPrintLevel(0); // Fit printing: -1 = no printing, 0 (minimal) to 3 (max)

    gROOT->SetBatch(kTRUE);
    gStyle->SetOptFit(1111);
    gStyle->SetPaintTextFormat("1.3f");
    gErrorIgnoreLevel = verbosity;

    // SetCustomColorPalette();
    // extra options for decorating final plots (pdf) in thesis
    gStyle->SetLineScalePS(2);
    gStyle->SetStatFontSize(0.03);
    gStyle->SetPadTickX(1); // Ticks on both top and bottom for X axis
    gStyle->SetPadTickY(1); // Ticks on both left and right for Y axis

    // remove ownership of objects from file so we can delete the file ptr
    TH1::AddDirectory(kFALSE);
    TH1::SetDefaultSumw2(kTRUE);

    outputFolder = SetOutputFolder(outputFolder);

    Int_t binsYieldDev = 212;
    Int_t binsMultiplierYieldDev = 0;
    Double_t minRangeYieldDev = -1.05;
    Double_t maxRangeYieldDev = +1.05;

    TH1 *def_effCorrPt_xim;                     // mult integrated efficiency corrected spectra for default cuts
    TH1 *def_effCorrPt_xip;                     // mult integrated efficiency corrected spectra for default cuts
    TH1 *def_effCorrPt_omm;                     // mult integrated efficiency corrected spectra for default cuts
    TH1 *def_effCorrPt_omp;                     // mult integrated efficiency corrected spectra for default cuts
    TH1 *def_effCorrPt_xiC;                     // mult integrated efficiency corrected spectra for default cuts
    TH1 *def_effCorrPt_omC;                     // mult integrated efficiency corrected spectra for default cuts
    TH1 *def_effCorrPt_xim_mult[fNmultbins_Xi]; // mult binned efficiency corrected spectra for default cuts
    TH1 *def_effCorrPt_xip_mult[fNmultbins_Xi]; // mult binned efficiency corrected spectra for default cuts
    TH1 *def_effCorrPt_omm_mult[fNmultbins_Om]; // mult binned efficiency corrected spectra for default cuts
    TH1 *def_effCorrPt_omp_mult[fNmultbins_Om]; // mult binned efficiency corrected spectra for default cuts
    TH1 *def_effCorrPt_xiC_mult[fNmultbins_Xi]; // mult binned efficiency corrected spectra for default cuts
    TH1 *def_effCorrPt_omC_mult[fNmultbins_Om]; // mult binned efficiency corrected spectra for default cuts

    TH1 *var_effCorrPt_xim;                     // mult integrated efficiency corrected spectra for varied cuts
    TH1 *var_effCorrPt_xip;                     // mult integrated efficiency corrected spectra for varied cuts
    TH1 *var_effCorrPt_omm;                     // mult integrated efficiency corrected spectra for varied cuts
    TH1 *var_effCorrPt_omp;                     // mult integrated efficiency corrected spectra for varied cuts
    TH1 *var_effCorrPt_xiC;                     // mult integrated efficiency corrected spectra for varied cuts
    TH1 *var_effCorrPt_omC;                     // mult integrated efficiency corrected spectra for varied cuts
    TH1 *var_effCorrPt_xim_mult[fNmultbins_Xi]; // mult binned efficiency corrected spectra for varied cuts
    TH1 *var_effCorrPt_xip_mult[fNmultbins_Xi]; // mult binned efficiency corrected spectra for varied cuts
    TH1 *var_effCorrPt_omm_mult[fNmultbins_Om]; // mult binned efficiency corrected spectra for varied cuts
    TH1 *var_effCorrPt_omp_mult[fNmultbins_Om]; // mult binned efficiency corrected spectra for varied cuts
    TH1 *var_effCorrPt_xiC_mult[fNmultbins_Xi]; // mult binned efficiency corrected spectra for varied cuts
    TH1 *var_effCorrPt_omC_mult[fNmultbins_Om]; // mult binned efficiency corrected spectra for varied cuts

    TH1D *varDefYieldDev_xip_pt[fNptbins_Xi];                     // mult integrated yield deviation: var/def
    TH1D *varDefYieldDev_xim_pt[fNptbins_Xi];                     // mult integrated yield deviation: var/def
    TH1D *varDefYieldDev_omp_pt[fNptbins_Om];                     // mult integrated yield deviation: var/def
    TH1D *varDefYieldDev_omm_pt[fNptbins_Om];                     // mult integrated yield deviation: var/def
    TH1D *varDefYieldDev_xiC_pt[fNptbins_Xi];                     // mult integrated yield deviation: var/def
    TH1D *varDefYieldDev_omC_pt[fNptbins_Om];                     // mult integrated yield deviation: var/def
    TH1D *varDefYieldDev_xip_pt_mult[fNptbins_Xi][fNmultbins_Xi]; // mult binned yield deviation: var/def
    TH1D *varDefYieldDev_xim_pt_mult[fNptbins_Xi][fNmultbins_Xi]; // mult binned yield deviation: var/def
    TH1D *varDefYieldDev_omp_pt_mult[fNptbins_Om][fNmultbins_Om]; // mult binned yield deviation: var/def
    TH1D *varDefYieldDev_omm_pt_mult[fNptbins_Om][fNmultbins_Om]; // mult binned yield deviation: var/def
    TH1D *varDefYieldDev_xiC_pt_mult[fNptbins_Xi][fNmultbins_Xi]; // mult binned yield deviation: var/def
    TH1D *varDefYieldDev_omC_pt_mult[fNptbins_Om][fNmultbins_Om]; // mult binned yield deviation: var/def

    TH1D *sysMultiTrial_xim;                     // mult integrated systematic uncertainty from MultiTrial
    TH1D *sysMultiTrial_xip;                     // mult integrated systematic uncertainty from MultiTrial
    TH1D *sysMultiTrial_omm;                     // mult integrated systematic uncertainty from MultiTrial
    TH1D *sysMultiTrial_omp;                     // mult integrated systematic uncertainty from MultiTrial
    TH1D *sysMultiTrial_xiC;                     // mult integrated systematic uncertainty from MultiTrial
    TH1D *sysMultiTrial_omC;                     // mult integrated systematic uncertainty from MultiTrial
    TH1D *sysMultiTrial_xim_mult[fNmultbins_Xi]; // mult binned systematic uncertainty from MultiTrial
    TH1D *sysMultiTrial_xip_mult[fNmultbins_Xi]; // mult binned systematic uncertainty from MultiTrial
    TH1D *sysMultiTrial_omm_mult[fNmultbins_Om]; // mult binned systematic uncertainty from MultiTrial
    TH1D *sysMultiTrial_omp_mult[fNmultbins_Om]; // mult binned systematic uncertainty from MultiTrial
    TH1D *sysMultiTrial_xiC_mult[fNmultbins_Xi]; // mult binned systematic uncertainty from MultiTrial
    TH1D *sysMultiTrial_omC_mult[fNmultbins_Om]; // mult binned systematic uncertainty from MultiTrial

    TH1D *fitMeanYieldDev_xim;                     // mult integrated yield deviation: mean of Gaussian fit
    TH1D *fitMeanYieldDev_xip;                     // mult integrated yield deviation: mean of Gaussian fit
    TH1D *fitMeanYieldDev_omm;                     // mult integrated yield deviation: mean of Gaussian fit
    TH1D *fitMeanYieldDev_omp;                     // mult integrated yield deviation: mean of Gaussian fit
    TH1D *fitMeanYieldDev_xiC;                     // mult integrated yield deviation: mean of Gaussian fit
    TH1D *fitMeanYieldDev_omC;                     // mult integrated yield deviation: mean of Gaussian fit
    TH1D *fitMeanYieldDev_xim_mult[fNmultbins_Xi]; // mult binned yield deviation: mean of Gaussian fit
    TH1D *fitMeanYieldDev_xip_mult[fNmultbins_Xi]; // mult binned yield deviation: mean of Gaussian fit
    TH1D *fitMeanYieldDev_omm_mult[fNmultbins_Om]; // mult binned yield deviation: mean of Gaussian fit
    TH1D *fitMeanYieldDev_omp_mult[fNmultbins_Om]; // mult binned yield deviation: mean of Gaussian fit
    TH1D *fitMeanYieldDev_xiC_mult[fNmultbins_Xi]; // mult binned yield deviation: mean of Gaussian fit
    TH1D *fitMeanYieldDev_omC_mult[fNmultbins_Om]; // mult binned yield deviation: mean of Gaussian fit

    Int_t failedRbTest[kNumSignedPart][fNptbins_Xi] = {0};                     // failed RB test counter
    Int_t failedRbTest_mult[kNumSignedPart][fNptbins_Xi][fNmultbins_Xi] = {0}; // failed RB test counter

    /// Getting default cut histograms:
    TString histName = "h3_ptmasscent_def";
    TString inputFileName = inputPath + "/" + histName + "/" + effCorrInputFilePrefix + "_" + histName + ".root";

    TFile *defInputFile = OpenFile(inputFileName);
    Info("SysMultiTrial: defInput", "Getting default cut histograms from '%s'", inputFileName.Data());
    def_effCorrPt_xim = (TH1 *)defInputFile->FindObjectAny("effCorrPt_xim");
    def_effCorrPt_xip = (TH1 *)defInputFile->FindObjectAny("effCorrPt_xip");
    def_effCorrPt_omm = (TH1 *)defInputFile->FindObjectAny("effCorrPt_omm");
    def_effCorrPt_omp = (TH1 *)defInputFile->FindObjectAny("effCorrPt_omp");
    def_effCorrPt_xiC = (TH1 *)defInputFile->FindObjectAny("effCorrPt_xiC");
    def_effCorrPt_omC = (TH1 *)defInputFile->FindObjectAny("effCorrPt_omC");

    for (Int_t multBinXi = 0; multBinXi < fNmultbins_Xi; multBinXi++)
    {
        def_effCorrPt_xip_mult[multBinXi] = (TH1 *)defInputFile->FindObjectAny(TString::Format(("effCorrPt_xip_mult[%d]"), multBinXi));
        def_effCorrPt_xim_mult[multBinXi] = (TH1 *)defInputFile->FindObjectAny(TString::Format(("effCorrPt_xim_mult[%d]"), multBinXi));
        def_effCorrPt_xiC_mult[multBinXi] = (TH1 *)defInputFile->FindObjectAny(TString::Format(("effCorrPt_xiC_mult[%d]"), multBinXi));
    }

    for (Int_t multBinOm = 0; multBinOm < fNmultbins_Om; multBinOm++)
    {
        def_effCorrPt_omp_mult[multBinOm] = (TH1 *)defInputFile->FindObjectAny(TString::Format(("effCorrPt_omp_mult[%d]"), multBinOm));
        def_effCorrPt_omm_mult[multBinOm] = (TH1 *)defInputFile->FindObjectAny(TString::Format(("effCorrPt_omm_mult[%d]"), multBinOm));
        def_effCorrPt_omC_mult[multBinOm] = (TH1 *)defInputFile->FindObjectAny(TString::Format(("effCorrPt_omC_mult[%d]"), multBinOm));
    }
    // input ended for defInputFile
    delete defInputFile;

    /// Creating output file:
    if (outputFileName.IsNull())
    {
        outputFileName = TString::Format("%s/SysUncertainty_multiTrial.root", outputFolder.Data());
    }

    Info("SysMultiTrial: I/O", "Default input successful. Creating output file '%s'", outputFileName.Data());

    /// Creating output directory structure:
    TFile *outputFile = OpenFile(outputFileName, "RECREATE");
    outputFile->mkdir("dirYieldDev_xip_pt");
    outputFile->mkdir("dirYieldDev_xim_pt");
    outputFile->mkdir("dirYieldDev_omp_pt");
    outputFile->mkdir("dirYieldDev_omm_pt");
    outputFile->mkdir("dirYieldDev_xiC_pt");
    outputFile->mkdir("dirYieldDev_omC_pt");
    outputFile->mkdir("dirYieldDev_xip_pt_mult");
    outputFile->mkdir("dirYieldDev_xim_pt_mult");
    outputFile->mkdir("dirYieldDev_omp_pt_mult");
    outputFile->mkdir("dirYieldDev_omm_pt_mult");
    outputFile->mkdir("dirYieldDev_xiC_pt_mult");
    outputFile->mkdir("dirYieldDev_omC_pt_mult");
    outputFile->mkdir("dirSysMultiTrial_xim_mult");
    outputFile->mkdir("dirSysMultiTrial_xip_mult");
    outputFile->mkdir("dirSysMultiTrial_omm_mult");
    outputFile->mkdir("dirSysMultiTrial_omp_mult");
    outputFile->mkdir("dirSysMultiTrial_xiC_mult");
    outputFile->mkdir("dirSysMultiTrial_omC_mult");
    /// Output set.

    /// Generate output objects
    /// Final output object - Systematic uncertainty from MultiTrial:
    sysMultiTrial_xip = new TH1D("sysMultiTrial_xip", "#Xi^{+}: Systematic uncertainty from MultiTrial: Mult 0-100%;#it{p}_{T} (GeV/#it{c});Relative Uncertainty(#sigma_{Gaus})", fNptbins_Xi, fPtbins_Xi);
    sysMultiTrial_xim = new TH1D("sysMultiTrial_xim", "#Xi^{-}: Systematic uncertainty from MultiTrial: Mult 0-100%;#it{p}_{T} (GeV/#it{c});Relative Uncertainty(#sigma_{Gaus})", fNptbins_Xi, fPtbins_Xi);
    sysMultiTrial_omp = new TH1D("sysMultiTrial_omp", "#Omega^{+}: Systematic uncertainty from MultiTrial: Mult 0-100%;#it{p}_{T} (GeV/#it{c});Relative Uncertainty(#sigma_{Gaus})", fNptbins_Om, fPtbins_Om);
    sysMultiTrial_omm = new TH1D("sysMultiTrial_omm", "#Omega^{-}: Systematic uncertainty from MultiTrial: Mult 0-100%;#it{p}_{T} (GeV/#it{c});Relative Uncertainty(#sigma_{Gaus})", fNptbins_Om, fPtbins_Om);
    sysMultiTrial_xiC = new TH1D("sysMultiTrial_xiC", "#Xi^{+} + #Xi^{-}: Systematic uncertainty from MultiTrial: Mult 0-100%;#it{p}_{T} (GeV/#it{c});Relative Uncertainty(#sigma_{Gaus})", fNptbins_Xi, fPtbins_Xi);
    sysMultiTrial_omC = new TH1D("sysMultiTrial_omC", "#Omega^{+} + #Omega^{-}: Systematic uncertainty from MultiTrial: Mult 0-100%;#it{p}_{T} (GeV/#it{c});Relative Uncertainty(#sigma_{Gaus})", fNptbins_Om, fPtbins_Om);
    for (Int_t multBinXi = 0; multBinXi < fNmultbins_Xi; multBinXi++)
    {
        sysMultiTrial_xip_mult[multBinXi] = new TH1D(TString::Format("sysMultiTrial_xip_mult[%d]", multBinXi), TString::Format("#Xi^{+}: Systematic uncertainty from MultiTrial: Mult %.0f-%.0f%%;#it{p}_{T} (GeV/#it{c});Relative Uncertainty(#sigma_{Gaus})", fMultbins_Xi[multBinXi], fMultbins_Xi[multBinXi + 1]), fNptbins_Xi, fPtbins_Xi);
        sysMultiTrial_xim_mult[multBinXi] = new TH1D(TString::Format("sysMultiTrial_xim_mult[%d]", multBinXi), TString::Format("#Xi^{-}: Systematic uncertainty from MultiTrial: Mult %.0f-%.0f%%;#it{p}_{T} (GeV/#it{c});Relative Uncertainty(#sigma_{Gaus})", fMultbins_Xi[multBinXi], fMultbins_Xi[multBinXi + 1]), fNptbins_Xi, fPtbins_Xi);
        sysMultiTrial_xiC_mult[multBinXi] = new TH1D(TString::Format("sysMultiTrial_xiC_mult[%d]", multBinXi), TString::Format("#Xi^{+} + #Xi^{-}: Systematic uncertainty from MultiTrial: Mult %.0f-%.0f%%;#it{p}_{T} (GeV/#it{c});Relative Uncertainty(#sigma_{Gaus})", fMultbins_Xi[multBinXi], fMultbins_Xi[multBinXi + 1]), fNptbins_Xi, fPtbins_Xi);
    }
    for (Int_t multBinOm = 0; multBinOm < fNmultbins_Om; multBinOm++)
    {
        sysMultiTrial_omp_mult[multBinOm] = new TH1D(TString::Format("sysMultiTrial_omp_mult[%d]", multBinOm), TString::Format("#Omega^{+}: Systematic uncertainty from MultiTrial: Mult %.0f-%.0f%%;#it{p}_{T} (GeV/#it{c});Relative Uncertainty(#sigma_{Gaus})", fMultbins_Om[multBinOm], fMultbins_Om[multBinOm + 1]), fNptbins_Om, fPtbins_Om);
        sysMultiTrial_omm_mult[multBinOm] = new TH1D(TString::Format("sysMultiTrial_omm_mult[%d]", multBinOm), TString::Format("#Omega^{-}: Systematic uncertainty from MultiTrial: Mult %.0f-%.0f%%;#it{p}_{T} (GeV/#it{c});Relative Uncertainty(#sigma_{Gaus})", fMultbins_Om[multBinOm], fMultbins_Om[multBinOm + 1]), fNptbins_Om, fPtbins_Om);
        sysMultiTrial_omC_mult[multBinOm] = new TH1D(TString::Format("sysMultiTrial_omC_mult[%d]", multBinOm), TString::Format("#Omega^{+} + #Omega^{-}: Systematic uncertainty from MultiTrial: Mult %.0f-%.0f%%;#it{p}_{T} (GeV/#it{c});Relative Uncertainty(#sigma_{Gaus})", fMultbins_Om[multBinOm], fMultbins_Om[multBinOm + 1]), fNptbins_Om, fPtbins_Om);
    }

    fitMeanYieldDev_xip = new TH1D("fitMeanYieldDev_xip", "#Xi^{+}: Mean of Gaussian fit: Mult 0-100%;#it{p}_{T} (GeV/#it{c});Mean(#mu_{Gaus})", fNptbins_Xi, fPtbins_Xi);
    fitMeanYieldDev_xim = new TH1D("fitMeanYieldDev_xim", "#Xi^{-}: Mean of Gaussian fit: Mult 0-100%;#it{p}_{T} (GeV/#it{c});Mean(#mu_{Gaus})", fNptbins_Xi, fPtbins_Xi);
    fitMeanYieldDev_omp = new TH1D("fitMeanYieldDev_omp", "#Omega^{+}: Mean of Gaussian fit: Mult 0-100%;#it{p}_{T} (GeV/#it{c});Mean(#mu_{Gaus})", fNptbins_Om, fPtbins_Om);
    fitMeanYieldDev_omm = new TH1D("fitMeanYieldDev_omm", "#Omega^{-}: Mean of Gaussian fit: Mult 0-100%;#it{p}_{T} (GeV/#it{c});Mean(#mu_{Gaus})", fNptbins_Om, fPtbins_Om);
    fitMeanYieldDev_xiC = new TH1D("fitMeanYieldDev_xiC", "#Xi^{+} + #Xi^{-}: Mean of Gaussian fit: Mult 0-100%;#it{p}_{T} (GeV/#it{c});Mean(#mu_{Gaus})", fNptbins_Xi, fPtbins_Xi);
    fitMeanYieldDev_omC = new TH1D("fitMeanYieldDev_omC", "#Omega^{+} + #Omega^{-}: Mean of Gaussian fit: Mult 0-100%;#it{p}_{T} (GeV/#it{c});Mean(#mu_{Gaus})", fNptbins_Om, fPtbins_Om);
    for (Int_t multBinXi = 0; multBinXi < fNmultbins_Xi; multBinXi++)
    {
        fitMeanYieldDev_xip_mult[multBinXi] = new TH1D(TString::Format("fitMeanYieldDev_xip_mult[%d]", multBinXi), TString::Format("#Xi^{+}: Mean of Gaussian fit: Mult %.0f-%.0f%%;#it{p}_{T} (GeV/#it{c});Mean(#mu_{Gaus})", fMultbins_Xi[multBinXi], fMultbins_Xi[multBinXi + 1]), fNptbins_Xi, fPtbins_Xi);
        fitMeanYieldDev_xim_mult[multBinXi] = new TH1D(TString::Format("fitMeanYieldDev_xim_mult[%d]", multBinXi), TString::Format("#Xi^{-}: Mean of Gaussian fit: Mult %.0f-%.0f%%;#it{p}_{T} (GeV/#it{c});Mean(#mu_{Gaus})", fMultbins_Xi[multBinXi], fMultbins_Xi[multBinXi + 1]), fNptbins_Xi, fPtbins_Xi);
        fitMeanYieldDev_xiC_mult[multBinXi] = new TH1D(TString::Format("fitMeanYieldDev_xiC_mult[%d]", multBinXi), TString::Format("#Xi^{+} + #Xi^{-}: Mean of Gaussian fit: Mult %.0f-%.0f%%;#it{p}_{T} (GeV/#it{c});Mean(#mu_{Gaus})", fMultbins_Xi[multBinXi], fMultbins_Xi[multBinXi + 1]), fNptbins_Xi, fPtbins_Xi);
    }
    for (Int_t multBinOm = 0; multBinOm < fNmultbins_Om; multBinOm++)
    {
        fitMeanYieldDev_omp_mult[multBinOm] = new TH1D(TString::Format("fitMeanYieldDev_omp_mult[%d]", multBinOm), TString::Format("#Omega^{+}: Mean of Gaussian fit: Mult %.0f-%.0f%%;#it{p}_{T} (GeV/#it{c});Mean(#mu_{Gaus})", fMultbins_Om[multBinOm], fMultbins_Om[multBinOm + 1]), fNptbins_Om, fPtbins_Om);
        fitMeanYieldDev_omm_mult[multBinOm] = new TH1D(TString::Format("fitMeanYieldDev_omm_mult[%d]", multBinOm), TString::Format("#Omega^{-}: Mean of Gaussian fit: Mult %.0f-%.0f%%;#it{p}_{T} (GeV/#it{c});Mean(#mu_{Gaus})", fMultbins_Om[multBinOm], fMultbins_Om[multBinOm + 1]), fNptbins_Om, fPtbins_Om);
        fitMeanYieldDev_omC_mult[multBinOm] = new TH1D(TString::Format("fitMeanYieldDev_omC_mult[%d]", multBinOm), TString::Format("#Omega^{+} + #Omega^{-}: Mean of Gaussian fit: Mult %.0f-%.0f%%;#it{p}_{T} (GeV/#it{c});Mean(#mu_{Gaus})", fMultbins_Om[multBinOm], fMultbins_Om[multBinOm + 1]), fNptbins_Om, fPtbins_Om);
    }

    /// Reference output objects - yield deviation from default cut, fitted with gaussian:
    for (Int_t ptBinXi = 0; ptBinXi < fNptbins_Xi; ptBinXi++)
    {
        varDefYieldDev_xip_pt[ptBinXi] = new TH1D(TString::Format("varDefYieldDev_xip_pt[%d]", ptBinXi), TString::Format("#Xi^{+}: #it{p}_{T}: <%.1f,%.1f>, Mult: 0-100%%;Y_{sys}/Y_{def} - 1;Counts", fPtbins_Xi[ptBinXi], fPtbins_Xi[ptBinXi + 1]), (binsYieldDev + binsMultiplierYieldDev * ptBinXi), minRangeYieldDev, maxRangeYieldDev);
        varDefYieldDev_xim_pt[ptBinXi] = new TH1D(TString::Format("varDefYieldDev_xim_pt[%d]", ptBinXi), TString::Format("#Xi^{-}: #it{p}_{T}: <%.1f,%.1f>, Mult: 0-100%%;Y_{sys}/Y_{def} - 1;Counts", fPtbins_Xi[ptBinXi], fPtbins_Xi[ptBinXi + 1]), (binsYieldDev + binsMultiplierYieldDev * ptBinXi), minRangeYieldDev, maxRangeYieldDev);
        varDefYieldDev_xiC_pt[ptBinXi] = new TH1D(TString::Format("varDefYieldDev_xiC_pt[%d]", ptBinXi), TString::Format("#Xi^{+} + #Xi^{-}: #it{p}_{T}: <%.1f,%.1f>, Mult: 0-100%%;Y_{sys}/Y_{def} - 1;Counts", fPtbins_Xi[ptBinXi], fPtbins_Xi[ptBinXi + 1]), (binsYieldDev + binsMultiplierYieldDev * ptBinXi), minRangeYieldDev, maxRangeYieldDev);
        for (Int_t multBinXi = 0; multBinXi < fNmultbins_Xi; multBinXi++)
        {
            varDefYieldDev_xip_pt_mult[ptBinXi][multBinXi] = new TH1D(TString::Format("varDefYieldDev_xip_pt_mult[%d][%d]", ptBinXi, multBinXi), TString::Format("#Xi^{+}: #it{p}_{T}: <%.1f,%.1f>, Mult: %.0f-%.0f%%;Y_{sys}/Y_{def} - 1;Counts", fPtbins_Xi[ptBinXi], fPtbins_Xi[ptBinXi + 1], fMultbins_Xi[multBinXi], fMultbins_Xi[multBinXi + 1]), (binsYieldDev + binsMultiplierYieldDev * ptBinXi), minRangeYieldDev, maxRangeYieldDev);
            varDefYieldDev_xim_pt_mult[ptBinXi][multBinXi] = new TH1D(TString::Format("varDefYieldDev_xim_pt_mult[%d][%d]", ptBinXi, multBinXi), TString::Format("#Xi^{-}: #it{p}_{T}: <%.1f,%.1f>, Mult: %.0f-%.0f%%;Y_{sys}/Y_{def} - 1;Counts", fPtbins_Xi[ptBinXi], fPtbins_Xi[ptBinXi + 1], fMultbins_Xi[multBinXi], fMultbins_Xi[multBinXi + 1]), (binsYieldDev + binsMultiplierYieldDev * ptBinXi), minRangeYieldDev, maxRangeYieldDev);
            varDefYieldDev_xiC_pt_mult[ptBinXi][multBinXi] = new TH1D(TString::Format("varDefYieldDev_xiC_pt_mult[%d][%d]", ptBinXi, multBinXi), TString::Format("#Xi^{+} + #Xi^{-}: #it{p}_{T}: <%.1f,%.1f>, Mult: %.0f-%.0f%%;Y_{sys}/Y_{def} - 1;Counts", fPtbins_Xi[ptBinXi], fPtbins_Xi[ptBinXi + 1], fMultbins_Xi[multBinXi], fMultbins_Xi[multBinXi + 1]), (binsYieldDev + binsMultiplierYieldDev * ptBinXi), minRangeYieldDev, maxRangeYieldDev);
        }
    }
    for (Int_t ptBinOm = 0; ptBinOm < fNptbins_Om; ptBinOm++)
    {
        varDefYieldDev_omp_pt[ptBinOm] = new TH1D(TString::Format("varDefYieldDev_omp_pt[%d]", ptBinOm), TString::Format("#Omega^{+}: #it{p}_{T}: <%.1f,%.1f>, Mult: 0-100%%;Y_{sys}/Y_{def} - 1;Counts", fPtbins_Om[ptBinOm], fPtbins_Om[ptBinOm + 1]), (binsYieldDev + binsMultiplierYieldDev * ptBinOm), minRangeYieldDev, maxRangeYieldDev);
        varDefYieldDev_omm_pt[ptBinOm] = new TH1D(TString::Format("varDefYieldDev_omm_pt[%d]", ptBinOm), TString::Format("#Omega^{-}: #it{p}_{T}: <%.1f,%.1f>, Mult: 0-100%%;Y_{sys}/Y_{def} - 1;Counts", fPtbins_Om[ptBinOm], fPtbins_Om[ptBinOm + 1]), (binsYieldDev + binsMultiplierYieldDev * ptBinOm), minRangeYieldDev, maxRangeYieldDev);
        varDefYieldDev_omC_pt[ptBinOm] = new TH1D(TString::Format("varDefYieldDev_omC_pt[%d]", ptBinOm), TString::Format("#Omega^{+} + #Omega^{-}: #it{p}_{T}: <%.1f,%.1f>, Mult: 0-100%%;Y_{sys}/Y_{def} - 1;Counts", fPtbins_Om[ptBinOm], fPtbins_Om[ptBinOm + 1]), (binsYieldDev + binsMultiplierYieldDev * ptBinOm), minRangeYieldDev, maxRangeYieldDev);
        for (Int_t multBinOm = 0; multBinOm < fNmultbins_Om; multBinOm++)
        {
            varDefYieldDev_omp_pt_mult[ptBinOm][multBinOm] = new TH1D(TString::Format("varDefYieldDev_omp_pt_mult[%d][%d]", ptBinOm, multBinOm), TString::Format("#Omega^{+}: #it{p}_{T}: <%.1f,%.1f>, Mult: %.0f-%.0f%%;Y_{sys}/Y_{def} - 1;Counts", fPtbins_Om[ptBinOm], fPtbins_Om[ptBinOm + 1], fMultbins_Om[multBinOm], fMultbins_Om[multBinOm + 1]), (binsYieldDev + binsMultiplierYieldDev * ptBinOm), minRangeYieldDev, maxRangeYieldDev);
            varDefYieldDev_omm_pt_mult[ptBinOm][multBinOm] = new TH1D(TString::Format("varDefYieldDev_omm_pt_mult[%d][%d]", ptBinOm, multBinOm), TString::Format("#Omega^{-}: #it{p}_{T}: <%.1f,%.1f>, Mult: %.0f-%.0f%%;Y_{sys}/Y_{def} - 1;Counts", fPtbins_Om[ptBinOm], fPtbins_Om[ptBinOm + 1], fMultbins_Om[multBinOm], fMultbins_Om[multBinOm + 1]), (binsYieldDev + binsMultiplierYieldDev * ptBinOm), minRangeYieldDev, maxRangeYieldDev);
            varDefYieldDev_omC_pt_mult[ptBinOm][multBinOm] = new TH1D(TString::Format("varDefYieldDev_omC_pt_mult[%d][%d]", ptBinOm, multBinOm), TString::Format("#Omega^{+} + #Omega^{-}: #it{p}_{T}: <%.1f,%.1f>, Mult: %.0f-%.0f%%;Y_{sys}/Y_{def} - 1;Counts", fPtbins_Om[ptBinOm], fPtbins_Om[ptBinOm + 1], fMultbins_Om[multBinOm], fMultbins_Om[multBinOm + 1]), (binsYieldDev + binsMultiplierYieldDev * ptBinOm), minRangeYieldDev, maxRangeYieldDev);
        }
    }
    Info("SysMultiTrial: I/O", "Output objects and output file created successfully.");

    /// Getting varied cut histograms:
    Int_t nVar = 500;
    for (Int_t iVar = 0; iVar < nVar; iVar++)
    {
        histName = TString::Format("h3Var_%d", iVar);
        inputFileName = inputPath + "/" + histName + "/" + effCorrInputFilePrefix + "_" + histName + ".root";

        TFile *varInputFile = OpenFile(inputFileName);
        Info("SysMultiTrial: varInput", "Getting varied cut histograms from '%s'", inputFileName.Data());
        var_effCorrPt_xim = (TH1 *)varInputFile->FindObjectAny("effCorrPt_xim");
        var_effCorrPt_xip = (TH1 *)varInputFile->FindObjectAny("effCorrPt_xip");
        var_effCorrPt_omm = (TH1 *)varInputFile->FindObjectAny("effCorrPt_omm");
        var_effCorrPt_omp = (TH1 *)varInputFile->FindObjectAny("effCorrPt_omp");
        var_effCorrPt_xiC = (TH1 *)varInputFile->FindObjectAny("effCorrPt_xiC");
        var_effCorrPt_omC = (TH1 *)varInputFile->FindObjectAny("effCorrPt_omC");

        for (Int_t multBinXi = 0; multBinXi < fNmultbins_Xi; multBinXi++)
        {
            var_effCorrPt_xip_mult[multBinXi] = (TH1 *)varInputFile->FindObjectAny(TString::Format(("effCorrPt_xip_mult[%d]"), multBinXi));
            var_effCorrPt_xim_mult[multBinXi] = (TH1 *)varInputFile->FindObjectAny(TString::Format(("effCorrPt_xim_mult[%d]"), multBinXi));
            var_effCorrPt_xiC_mult[multBinXi] = (TH1 *)varInputFile->FindObjectAny(TString::Format(("effCorrPt_xiC_mult[%d]"), multBinXi));
        }

        for (Int_t multBinOm = 0; multBinOm < fNmultbins_Om; multBinOm++)
        {
            var_effCorrPt_omp_mult[multBinOm] = (TH1 *)varInputFile->FindObjectAny(TString::Format(("effCorrPt_omp_mult[%d]"), multBinOm));
            var_effCorrPt_omm_mult[multBinOm] = (TH1 *)varInputFile->FindObjectAny(TString::Format(("effCorrPt_omm_mult[%d]"), multBinOm));
            var_effCorrPt_omC_mult[multBinOm] = (TH1 *)varInputFile->FindObjectAny(TString::Format(("effCorrPt_omC_mult[%d]"), multBinOm));
        }
        // input ended for varInputFile
        delete varInputFile;

        /// Compute yield deviation of variation from default cut and pass RB criteria:
        for (Int_t ptBinXi = 0; ptBinXi < fNptbins_Xi; ptBinXi++)
        {
            failedRbTest[kXip][ptBinXi] += ComputeRogerBarlow(var_effCorrPt_xip, def_effCorrPt_xip, varDefYieldDev_xip_pt[ptBinXi], ptBinXi + 1, iVar, fDebug);
            failedRbTest[kXim][ptBinXi] += ComputeRogerBarlow(var_effCorrPt_xim, def_effCorrPt_xim, varDefYieldDev_xim_pt[ptBinXi], ptBinXi + 1, iVar, fDebug);
            failedRbTest[kXiC][ptBinXi] += ComputeRogerBarlow(var_effCorrPt_xiC, def_effCorrPt_xiC, varDefYieldDev_xiC_pt[ptBinXi], ptBinXi + 1, iVar, fDebug);

            for (Int_t multBinXi = 0; multBinXi < fNmultbins_Xi; multBinXi++)
            {
                failedRbTest_mult[kXip][ptBinXi][multBinXi] += ComputeRogerBarlow(var_effCorrPt_xip_mult[multBinXi], def_effCorrPt_xip_mult[multBinXi], varDefYieldDev_xip_pt_mult[ptBinXi][multBinXi], ptBinXi + 1, iVar, fDebug);
                failedRbTest_mult[kXim][ptBinXi][multBinXi] += ComputeRogerBarlow(var_effCorrPt_xim_mult[multBinXi], def_effCorrPt_xim_mult[multBinXi], varDefYieldDev_xim_pt_mult[ptBinXi][multBinXi], ptBinXi + 1, iVar, fDebug);
                failedRbTest_mult[kXiC][ptBinXi][multBinXi] += ComputeRogerBarlow(var_effCorrPt_xiC_mult[multBinXi], def_effCorrPt_xiC_mult[multBinXi], varDefYieldDev_xiC_pt_mult[ptBinXi][multBinXi], ptBinXi + 1, iVar, fDebug);
            }
        }
        for (Int_t ptBinOm = 0; ptBinOm < fNptbins_Om; ptBinOm++)
        {
            failedRbTest[kOmp][ptBinOm] += ComputeRogerBarlow(var_effCorrPt_omp, def_effCorrPt_omp, varDefYieldDev_omp_pt[ptBinOm], ptBinOm + 1, iVar, fDebug);
            failedRbTest[kOmm][ptBinOm] += ComputeRogerBarlow(var_effCorrPt_omm, def_effCorrPt_omm, varDefYieldDev_omm_pt[ptBinOm], ptBinOm + 1, iVar, fDebug);
            failedRbTest[kOmC][ptBinOm] += ComputeRogerBarlow(var_effCorrPt_omC, def_effCorrPt_omC, varDefYieldDev_omC_pt[ptBinOm], ptBinOm + 1, iVar, fDebug);

            for (Int_t multBinOm = 0; multBinOm < fNmultbins_Om; multBinOm++)
            {
                failedRbTest_mult[kOmp][ptBinOm][multBinOm] += ComputeRogerBarlow(var_effCorrPt_omp_mult[multBinOm], def_effCorrPt_omp_mult[multBinOm], varDefYieldDev_omp_pt_mult[ptBinOm][multBinOm], ptBinOm + 1, iVar, fDebug);
                failedRbTest_mult[kOmm][ptBinOm][multBinOm] += ComputeRogerBarlow(var_effCorrPt_omm_mult[multBinOm], def_effCorrPt_omm_mult[multBinOm], varDefYieldDev_omm_pt_mult[ptBinOm][multBinOm], ptBinOm + 1, iVar, fDebug);
                failedRbTest_mult[kOmC][ptBinOm][multBinOm] += ComputeRogerBarlow(var_effCorrPt_omC_mult[multBinOm], def_effCorrPt_omC_mult[multBinOm], varDefYieldDev_omC_pt_mult[ptBinOm][multBinOm], ptBinOm + 1, iVar, fDebug);
            }
        }
    }
    Info("SysMultiTrial: varInput", "Variation input successful. Yield deviation computed.");

    /// Fit the yield deviation histograms with gaussian and save sigmaGaus as systematic uncertainty:
    Info("SysMultiTrial: Fit Xi", "Fitting yield deviation histograms with gaussian and extracting sigmaGaus as systematic uncertainty.");
    for (Int_t ptBinXi = 0; ptBinXi < fNptbins_Xi; ptBinXi++)
    {
        sysMultiTrial_xip->SetBinContent(ptBinXi + 1, FitGaus(varDefYieldDev_xip_pt[ptBinXi]));
        sysMultiTrial_xim->SetBinContent(ptBinXi + 1, FitGaus(varDefYieldDev_xim_pt[ptBinXi]));
        sysMultiTrial_xiC->SetBinContent(ptBinXi + 1, FitGaus(varDefYieldDev_xiC_pt[ptBinXi]));

        // add gaus fit's mean to fitMeanYieldDev histograms
        fitMeanYieldDev_xip->SetBinContent(ptBinXi + 1, ((TF1 *)varDefYieldDev_xip_pt[ptBinXi]->GetListOfFunctions()->At(0))->GetParameter(1));
        fitMeanYieldDev_xim->SetBinContent(ptBinXi + 1, ((TF1 *)varDefYieldDev_xim_pt[ptBinXi]->GetListOfFunctions()->At(0))->GetParameter(1));
        fitMeanYieldDev_xiC->SetBinContent(ptBinXi + 1, ((TF1 *)varDefYieldDev_xiC_pt[ptBinXi]->GetListOfFunctions()->At(0))->GetParameter(1));

        for (Int_t multBinXi = 0; multBinXi < fNmultbins_Xi; multBinXi++)
        {
            sysMultiTrial_xip_mult[multBinXi]->SetBinContent(ptBinXi + 1, FitGaus(varDefYieldDev_xip_pt_mult[ptBinXi][multBinXi]));
            sysMultiTrial_xim_mult[multBinXi]->SetBinContent(ptBinXi + 1, FitGaus(varDefYieldDev_xim_pt_mult[ptBinXi][multBinXi]));
            sysMultiTrial_xiC_mult[multBinXi]->SetBinContent(ptBinXi + 1, FitGaus(varDefYieldDev_xiC_pt_mult[ptBinXi][multBinXi]));

            // add gaus fit's mean to fitMeanYieldDev histograms
            fitMeanYieldDev_xip_mult[multBinXi]->SetBinContent(ptBinXi + 1, ((TF1 *)varDefYieldDev_xip_pt_mult[ptBinXi][multBinXi]->GetListOfFunctions()->At(0))->GetParameter(1));
            fitMeanYieldDev_xim_mult[multBinXi]->SetBinContent(ptBinXi + 1, ((TF1 *)varDefYieldDev_xim_pt_mult[ptBinXi][multBinXi]->GetListOfFunctions()->At(0))->GetParameter(1));
            fitMeanYieldDev_xiC_mult[multBinXi]->SetBinContent(ptBinXi + 1, ((TF1 *)varDefYieldDev_xiC_pt_mult[ptBinXi][multBinXi]->GetListOfFunctions()->At(0))->GetParameter(1));
        }
    }
    Info("SysMultiTrial: Fit Om", "Fitting yield deviation histograms with gaussian and extracting sigmaGaus as systematic uncertainty.");
    for (Int_t ptBinOm = 0; ptBinOm < fNptbins_Om; ptBinOm++)
    {
        sysMultiTrial_omp->SetBinContent(ptBinOm + 1, FitGaus(varDefYieldDev_omp_pt[ptBinOm]));
        sysMultiTrial_omm->SetBinContent(ptBinOm + 1, FitGaus(varDefYieldDev_omm_pt[ptBinOm]));
        sysMultiTrial_omC->SetBinContent(ptBinOm + 1, FitGaus(varDefYieldDev_omC_pt[ptBinOm]));

        // add gaus fit's mean to fitMeanYieldDev histograms
        fitMeanYieldDev_omp->SetBinContent(ptBinOm + 1, ((TF1 *)varDefYieldDev_omp_pt[ptBinOm]->GetListOfFunctions()->At(0))->GetParameter(1));
        fitMeanYieldDev_omm->SetBinContent(ptBinOm + 1, ((TF1 *)varDefYieldDev_omm_pt[ptBinOm]->GetListOfFunctions()->At(0))->GetParameter(1));
        fitMeanYieldDev_omC->SetBinContent(ptBinOm + 1, ((TF1 *)varDefYieldDev_omC_pt[ptBinOm]->GetListOfFunctions()->At(0))->GetParameter(1));

        for (Int_t multBinOm = 0; multBinOm < fNmultbins_Om; multBinOm++)
        {
            sysMultiTrial_omp_mult[multBinOm]->SetBinContent(ptBinOm + 1, FitGaus(varDefYieldDev_omp_pt_mult[ptBinOm][multBinOm]));
            sysMultiTrial_omm_mult[multBinOm]->SetBinContent(ptBinOm + 1, FitGaus(varDefYieldDev_omm_pt_mult[ptBinOm][multBinOm]));
            sysMultiTrial_omC_mult[multBinOm]->SetBinContent(ptBinOm + 1, FitGaus(varDefYieldDev_omC_pt_mult[ptBinOm][multBinOm]));

            // add gaus fit's mean to fitMeanYieldDev histograms
            fitMeanYieldDev_omp_mult[multBinOm]->SetBinContent(ptBinOm + 1, ((TF1 *)varDefYieldDev_omp_pt_mult[ptBinOm][multBinOm]->GetListOfFunctions()->At(0))->GetParameter(1));
            fitMeanYieldDev_omm_mult[multBinOm]->SetBinContent(ptBinOm + 1, ((TF1 *)varDefYieldDev_omm_pt_mult[ptBinOm][multBinOm]->GetListOfFunctions()->At(0))->GetParameter(1));
            fitMeanYieldDev_omC_mult[multBinOm]->SetBinContent(ptBinOm + 1, ((TF1 *)varDefYieldDev_omC_pt_mult[ptBinOm][multBinOm]->GetListOfFunctions()->At(0))->GetParameter(1));
        }
    }

    /// Write output objects to file and save images:
    outputFile->cd("dirSysMultiTrial_xip_mult");
    sysMultiTrial_xip->Write();
    fitMeanYieldDev_xip->Write();
    outputFile->cd("dirSysMultiTrial_xim_mult");
    sysMultiTrial_xim->Write();
    fitMeanYieldDev_xim->Write();
    outputFile->cd("dirSysMultiTrial_xiC_mult");
    sysMultiTrial_xiC->Write();
    fitMeanYieldDev_xiC->Write();
    outputFile->cd("dirSysMultiTrial_omp_mult");
    sysMultiTrial_omp->Write();
    fitMeanYieldDev_omp->Write();
    outputFile->cd("dirSysMultiTrial_omm_mult");
    sysMultiTrial_omm->Write();
    fitMeanYieldDev_omm->Write();
    outputFile->cd("dirSysMultiTrial_omC_mult");
    sysMultiTrial_omC->Write();
    fitMeanYieldDev_omC->Write();
    if (saveImages)
    {
        DrawAndSaveImage(sysMultiTrial_xip, outputFolder, "SysMultiTrial_xip", imageFormat);
        DrawAndSaveImage(sysMultiTrial_xim, outputFolder, "SysMultiTrial_xim", imageFormat);
        DrawAndSaveImage(sysMultiTrial_xiC, outputFolder, "SysMultiTrial_xiC", imageFormat);
        DrawAndSaveImage(sysMultiTrial_omp, outputFolder, "SysMultiTrial_omp", imageFormat);
        DrawAndSaveImage(sysMultiTrial_omm, outputFolder, "SysMultiTrial_omm", imageFormat);
        DrawAndSaveImage(sysMultiTrial_omC, outputFolder, "SysMultiTrial_omC", imageFormat);

        DrawAndSaveImage(fitMeanYieldDev_xip, outputFolder, "FitMeanYieldDev_xip", imageFormat);
        DrawAndSaveImage(fitMeanYieldDev_xim, outputFolder, "FitMeanYieldDev_xim", imageFormat);
        DrawAndSaveImage(fitMeanYieldDev_xiC, outputFolder, "FitMeanYieldDev_xiC", imageFormat);
        DrawAndSaveImage(fitMeanYieldDev_omp, outputFolder, "FitMeanYieldDev_omp", imageFormat);
        DrawAndSaveImage(fitMeanYieldDev_omm, outputFolder, "FitMeanYieldDev_omm", imageFormat);
        DrawAndSaveImage(fitMeanYieldDev_omC, outputFolder, "FitMeanYieldDev_omC", imageFormat);
    }
    for (Int_t multBinXi = 0; multBinXi < fNmultbins_Xi; multBinXi++)
    {
        outputFile->cd("dirSysMultiTrial_xip_mult");
        sysMultiTrial_xip_mult[multBinXi]->Write();
        fitMeanYieldDev_xip_mult[multBinXi]->Write();
        outputFile->cd("dirSysMultiTrial_xim_mult");
        sysMultiTrial_xim_mult[multBinXi]->Write();
        fitMeanYieldDev_xim_mult[multBinXi]->Write();
        outputFile->cd("dirSysMultiTrial_xiC_mult");
        sysMultiTrial_xiC_mult[multBinXi]->Write();
        fitMeanYieldDev_xiC_mult[multBinXi]->Write();
        if (saveImages)
        {
            DrawAndSaveImage(sysMultiTrial_xip_mult[multBinXi], outputFolder, "SysMultiTrial_xip", imageFormat);
            DrawAndSaveImage(sysMultiTrial_xim_mult[multBinXi], outputFolder, "SysMultiTrial_xim", imageFormat);
            DrawAndSaveImage(sysMultiTrial_xiC_mult[multBinXi], outputFolder, "SysMultiTrial_xiC", imageFormat);

            DrawAndSaveImage(fitMeanYieldDev_xip_mult[multBinXi], outputFolder, "FitMeanYieldDev_xip", imageFormat);
            DrawAndSaveImage(fitMeanYieldDev_xim_mult[multBinXi], outputFolder, "FitMeanYieldDev_xim", imageFormat);
            DrawAndSaveImage(fitMeanYieldDev_xiC_mult[multBinXi], outputFolder, "FitMeanYieldDev_xiC", imageFormat);
        }
    }
    for (Int_t multBinOm = 0; multBinOm < fNmultbins_Om; multBinOm++)
    {
        outputFile->cd("dirSysMultiTrial_omp_mult");
        sysMultiTrial_omp_mult[multBinOm]->Write();
        fitMeanYieldDev_omp_mult[multBinOm]->Write();
        outputFile->cd("dirSysMultiTrial_omm_mult");
        sysMultiTrial_omm_mult[multBinOm]->Write();
        fitMeanYieldDev_omm_mult[multBinOm]->Write();
        outputFile->cd("dirSysMultiTrial_omC_mult");
        sysMultiTrial_omC_mult[multBinOm]->Write();
        fitMeanYieldDev_omC_mult[multBinOm]->Write();
        if (saveImages)
        {
            DrawAndSaveImage(sysMultiTrial_omp_mult[multBinOm], outputFolder, "SysMultiTrial_omp", imageFormat);
            DrawAndSaveImage(sysMultiTrial_omm_mult[multBinOm], outputFolder, "SysMultiTrial_omm", imageFormat);
            DrawAndSaveImage(sysMultiTrial_omC_mult[multBinOm], outputFolder, "SysMultiTrial_omC", imageFormat);

            DrawAndSaveImage(fitMeanYieldDev_omp_mult[multBinOm], outputFolder, "FitMeanYieldDev_omp", imageFormat);
            DrawAndSaveImage(fitMeanYieldDev_omm_mult[multBinOm], outputFolder, "FitMeanYieldDev_omm", imageFormat);
            DrawAndSaveImage(fitMeanYieldDev_omC_mult[multBinOm], outputFolder, "FitMeanYieldDev_omC", imageFormat);
        }
    }

    for (Int_t ptBinXi = 0; ptBinXi < fNptbins_Xi; ptBinXi++)
    {
        outputFile->cd("dirYieldDev_xip_pt");
        varDefYieldDev_xip_pt[ptBinXi]->Write();
        outputFile->cd("dirYieldDev_xim_pt");
        varDefYieldDev_xim_pt[ptBinXi]->Write();
        outputFile->cd("dirYieldDev_xiC_pt");
        varDefYieldDev_xiC_pt[ptBinXi]->Write();

        if (saveImages)
        {
            DrawAndSaveImage(varDefYieldDev_xip_pt[ptBinXi], outputFolder, "YieldDev_xip_pt", imageFormat, failedRbTest[kXip][ptBinXi]);
            DrawAndSaveImage(varDefYieldDev_xim_pt[ptBinXi], outputFolder, "YieldDev_xim_pt", imageFormat, failedRbTest[kXim][ptBinXi]);
            DrawAndSaveImage(varDefYieldDev_xiC_pt[ptBinXi], outputFolder, "YieldDev_xiC_pt", imageFormat, failedRbTest[kXiC][ptBinXi]);
        }

        for (Int_t multBinXi = 0; multBinXi < fNmultbins_Xi; multBinXi++)
        {
            outputFile->cd("dirYieldDev_xip_pt_mult");
            varDefYieldDev_xip_pt_mult[ptBinXi][multBinXi]->Write();
            outputFile->cd("dirYieldDev_xim_pt_mult");
            varDefYieldDev_xim_pt_mult[ptBinXi][multBinXi]->Write();
            outputFile->cd("dirYieldDev_xiC_pt_mult");
            varDefYieldDev_xiC_pt_mult[ptBinXi][multBinXi]->Write();

            if (saveImages)
            {
                DrawAndSaveImage(varDefYieldDev_xip_pt_mult[ptBinXi][multBinXi], outputFolder, "YieldDev_xip_pt_mult", imageFormat, failedRbTest_mult[kXip][ptBinXi][multBinXi]);
                DrawAndSaveImage(varDefYieldDev_xim_pt_mult[ptBinXi][multBinXi], outputFolder, "YieldDev_xim_pt_mult", imageFormat, failedRbTest_mult[kXim][ptBinXi][multBinXi]);
                DrawAndSaveImage(varDefYieldDev_xiC_pt_mult[ptBinXi][multBinXi], outputFolder, "YieldDev_xiC_pt_mult", imageFormat, failedRbTest_mult[kXiC][ptBinXi][multBinXi]);
            }
        }
    }
    for (Int_t ptBinOm = 0; ptBinOm < fNptbins_Om; ptBinOm++)
    {
        outputFile->cd("dirYieldDev_omp_pt");
        varDefYieldDev_omp_pt[ptBinOm]->Write();
        outputFile->cd("dirYieldDev_omm_pt");
        varDefYieldDev_omm_pt[ptBinOm]->Write();
        outputFile->cd("dirYieldDev_omC_pt");
        varDefYieldDev_omC_pt[ptBinOm]->Write();

        if (saveImages)
        {
            DrawAndSaveImage(varDefYieldDev_omp_pt[ptBinOm], outputFolder, "YieldDev_omp_pt", imageFormat, failedRbTest[kOmp][ptBinOm]);
            DrawAndSaveImage(varDefYieldDev_omm_pt[ptBinOm], outputFolder, "YieldDev_omm_pt", imageFormat, failedRbTest[kOmm][ptBinOm]);
            DrawAndSaveImage(varDefYieldDev_omC_pt[ptBinOm], outputFolder, "YieldDev_omC_pt", imageFormat, failedRbTest[kOmC][ptBinOm]);
        }
        for (Int_t multBinOm = 0; multBinOm < fNmultbins_Om; multBinOm++)
        {
            outputFile->cd("dirYieldDev_omp_pt_mult");
            varDefYieldDev_omp_pt_mult[ptBinOm][multBinOm]->Write();
            outputFile->cd("dirYieldDev_omm_pt_mult");
            varDefYieldDev_omm_pt_mult[ptBinOm][multBinOm]->Write();
            outputFile->cd("dirYieldDev_omC_pt_mult");
            varDefYieldDev_omC_pt_mult[ptBinOm][multBinOm]->Write();

            if (saveImages)
            {
                DrawAndSaveImage(varDefYieldDev_omp_pt_mult[ptBinOm][multBinOm], outputFolder, "YieldDev_omp_pt_mult", imageFormat, failedRbTest_mult[kOmp][ptBinOm][multBinOm]);
                DrawAndSaveImage(varDefYieldDev_omm_pt_mult[ptBinOm][multBinOm], outputFolder, "YieldDev_omm_pt_mult", imageFormat, failedRbTest_mult[kOmm][ptBinOm][multBinOm]);
                DrawAndSaveImage(varDefYieldDev_omC_pt_mult[ptBinOm][multBinOm], outputFolder, "YieldDev_omC_pt_mult", imageFormat, failedRbTest_mult[kOmC][ptBinOm][multBinOm]);
            }
        }
    }
    Info("SysMultiTrial: I/O", "Output objects written to file '%s' successfully.", outputFileName.Data());

    delete outputFile;

    return 0;
}

Double_t ComputeRogerBarlow(TH1 *hVar, TH1 *hDef, TH1 *hYieldDev, Int_t ptBin, Int_t iVar, Bool_t fDebug, UInt_t nSigma)
{
    Double_t varVal = hVar->GetBinContent(ptBin);
    Double_t varErr = hVar->GetBinError(ptBin);
    Double_t defVal = hDef->GetBinContent(ptBin);
    Double_t defErr = hDef->GetBinError(ptBin);
    Double_t yieldDev = 0.0;
    Double_t sigmaRB = 0.0;
    Double_t nSigmaRB = 0.0;
    Int_t failedRbTest = 0;

    if (defVal > 1e-12)
    {
        yieldDev = varVal / defVal - 1;
        sigmaRB = TMath::Sqrt(TMath::Abs(TMath::Power(varErr, 2) - TMath::Power(defErr, 2))); // Computation of roger barlow sigma_{delta}
        sigmaRB = sigmaRB / defVal;
        nSigmaRB = nSigma * sigmaRB;

        // if (TMath::Abs(yieldDev) > (nSigmaRB))
        // {
        if (fDebug)
            Info("ComputeRogerBarlow", "iVar %d: Filling %s: ptBin: %d, varVal: %f, defVal: %f, yieldDev: %f, nSigmaRB: %f", iVar, hYieldDev->GetName(), ptBin, varVal, defVal, yieldDev, nSigmaRB);
        hYieldDev->Fill(yieldDev);
        // }
        // else
        // {
        // if (fDebug)
        // Info("ComputeRogerBarlow", "iVar %d: Skipping %s: ptBin: %d, varVal: %f, defVal: %f, yieldDev: %f, nSigmaRB: %f. Failed RB criteria.", iVar, hYieldDev->GetName(), ptBin, varVal, defVal, yieldDev, nSigmaRB);
        // }
        if (TMath::Abs(yieldDev) < (nSigmaRB))
            failedRbTest = 1;
    }
    else
    {
        Warning("ComputeRogerBarlow", "iVar %d: Skipping %s: ptBin: %d, varVal: %f, defVal: %f. Division by zero.", iVar, hYieldDev->GetName(), ptBin, varVal, defVal);
    }

    // return nSigmaRB;
    return failedRbTest;
}

Double_t FitGaus(TH1 *hist, TString fitOptions)
{

    Double_t sigmaGaus = 0.0;
    Int_t fitStatus = 0;

    Double_t hAmp = hist->GetMaximum();
    Double_t xAmp = hist->GetXaxis()->GetBinCenter(hist->GetMaximumBin());
    Double_t hSigma = hist->GetStdDev();
    Double_t fitMin = xAmp - 3 * hSigma;
    Double_t fitMax = xAmp + 3 * hSigma;

    TF1 *fTempGaus = new TF1("fTempGaus", "gaus", -1, +1);
    fTempGaus->SetParameter(0, hAmp);
    fTempGaus->SetParameter(1, xAmp);
    for (int iFit = 0; iFit < 5; iFit++)
    {
        hist->Fit(fTempGaus, "WLNQ MULTITHREAD", "", fitMin, fitMax);
        sigmaGaus = TMath::Abs(fTempGaus->GetParameter(2));
        fTempGaus->SetParameter(2, sigmaGaus);
        fitMin = xAmp - 2.5 * sigmaGaus;
        fitMax = xAmp + 2.5 * sigmaGaus;
    }

    TF1 *fGaus = new TF1(TString::Format("fGaus_%s", hist->GetName()), "gaus", -1, 1);
    fGaus->SetNpx(1000);
    fGaus->SetParameter(0, hAmp);
    fGaus->SetParLimits(0, 0.5 * hAmp, 1.1 * hAmp); // Limit amplitude to 50% to 120% of max
    fGaus->SetParameter(1, xAmp);
    fGaus->SetParLimits(1, xAmp - sigmaGaus, xAmp + sigmaGaus);
    fGaus->SetParameter(2, sigmaGaus);

    for (int iFit = 0; iFit < 5; iFit++) // Fit 5 times to get better fit
    {
        fitStatus = hist->Fit(fGaus, fitOptions.Data(), "", fitMin, fitMax);
        sigmaGaus = TMath::Abs(fGaus->GetParameter(2));
        fGaus->SetParameter(2, sigmaGaus);
        fitMin = xAmp - 2.5 * sigmaGaus;
        fitMax = xAmp + 2.5 * sigmaGaus;
    }
    if (fitStatus != 0)
    {
        Warning("FitGaus", "Fit failed for %s. Fit status: %d, Fit mean = %.2f, Fit sigma = %.2f, Fit range: [%.2f, %.2f], Entries: %.1f. Trying again:", hist->GetName(), fitStatus, fGaus->GetParameter(1), sigmaGaus, fitMin, fitMax, hist->GetEntries());
        fitOptions.ReplaceAll("Q", "");
        fitStatus = hist->Fit(fGaus, fitOptions.Data(), "", fitMin, fitMax);
        if (fitStatus != 0)
            Error("FitGaus", "Fit failed for %s. Fit status: %d, Fit mean = %.2f, Fit sigma = %.2f, Fit range: [%.2f, %.2f], Entries: %.1f. Exiting.", hist->GetName(), fitStatus, fGaus->GetParameter(1), sigmaGaus, fitMin, fitMax, hist->GetEntries());
    }

    sigmaGaus = TMath::Abs(fGaus->GetParameter(2));
    hist->GetListOfFunctions()->Add(fGaus);

    return sigmaGaus;
}

void DrawAndSaveImage(TH1 *hist, TString outputFolder, TString imageFolder, TString imageFormat, Int_t failedRbCount)
{
    TString histTitle = hist->GetTitle();
    TCanvas *cDraw = new TCanvas(hist->GetName(), hist->GetTitle(), 1200, 900);
    cDraw->cd();
    // cDraw->SetLeftMargin(0.8);  // New mod for final thesis plots
    cDraw->SetRightMargin(0.03); // New mod for final thesis plots
    cDraw->SetTopMargin(0.03);   // New mod for final thesis plots
    // auto legend = new TLegend(0.1, 0.7, 0.28, 0.9);
    /// Legend for final plots:
    hist->SetStats(0);
    auto legend = new TLegend(0.72, 0.72, 0.92, 0.92);
    legend->SetBorderSize(0);

    TF1 *func = (TF1 *)hist->GetListOfFunctions()->At(0);
    if (func)
    {
        Int_t ndf_sanitized = func->GetNDF();
        if (ndf_sanitized == 0)
            ndf_sanitized = 1;

        // New mod for final thesis plots
        // legend->AddEntry(peak, "p-Pb #sqrt{s_{NN}} = 8.16 TeV", "");
        // legend->AddEntry(peak, TString::Format("%s, %s", ptRange.Data(), multRange.Data()), "");
        legend->SetTextSize(0.025);
        legend->SetHeader("Fit Stats (Gaus)", "C"); // option "C" allows to center the header
        legend->AddEntry(hist, "Random Trial", "lpe");
        legend->AddEntry(func, TString::Format("Fit mean (#mu) = %.3f", func->GetParameter(1)), "l");
        legend->AddEntry(func, TString::Format("Fit sigma (#sigma)= %.3f", func->GetParameter(2)), "l");
        legend->AddEntry(func, TString::Format("#frac{#chi^{2}}{NDF} = %.1f", (func->GetChisquare() / ndf_sanitized)), "l");
        // legend->AddEntry(hist, TString::Format("Failed RB = %d", failedRbCount), "pe");

        TString partTitle = histTitle(0, histTitle.First(":"));
        TString collInfo = "p-Pb #sqrt{s_{NN}} = 8.16 TeV, This work";
        TString ptRange = histTitle(histTitle.First("<") + 1, histTitle.First(">") - histTitle.First("<") - 1) + " GeV/#it{c}";
        TString multRange = histTitle(histTitle.Index("Mult: ") + 6, histTitle.Length()) + " V0A";
        ptRange.ReplaceAll(",", " < #it{p}_{T} < ");
        // multRange.ReplaceAll(",", "-");
        TLatex *latexPartTitle = new TLatex(0.15, 0.88, partTitle.Data());
        TLatex *latexCollInfo = new TLatex(0.15, 0.84, collInfo.Data());
        TLatex *latexPtRange = new TLatex(0.15, 0.80, ptRange.Data());
        TLatex *latexMultRange = new TLatex(0.15, 0.76, multRange.Data());
        latexPartTitle->SetTextSize(0.035);
        latexCollInfo->SetTextSize(0.028);
        latexPtRange->SetTextSize(0.028);
        latexMultRange->SetTextSize(0.028);
        latexPartTitle->SetNDC(kTRUE);
        latexCollInfo->SetNDC(kTRUE);
        latexPtRange->SetNDC(kTRUE);
        latexMultRange->SetNDC(kTRUE);
        hist->SetTitle("");

        hist->GetYaxis()->SetRangeUser(0., hist->GetMaximum() * 1.2);
        hist->SetMarkerStyle(kFullCircle);
        hist->SetColors(kBlue, kBlack);
        // hist->Draw("P LF2 HIST SAME");
        // hist->Draw("HIST");
        hist->Draw("P E1");
        func->Draw("same");

        legend->Draw();
        latexPartTitle->Draw();
        latexCollInfo->Draw();
        latexPtRange->Draw();
        latexMultRange->Draw();
    }
    else if (TString(hist->GetName()).Contains("Mean"))
    {
        hist->SetMarkerStyle(kFullCircle);
        hist->Draw("TEXT00");
    }
    else
    {
        hist->GetYaxis()->SetRangeUser(0., hist->GetMaximum() * 1.4);
        hist->SetMarkerStyle(kFullCircle);
        hist->Draw("TEXT00");
    }

    SaveImage(outputFolder, imageFolder, hist->GetName(), imageFormat, cDraw);

    delete legend;
    delete cDraw;
}