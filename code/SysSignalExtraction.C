#include <TROOT.h>
#include <TStyle.h>
#include <TLegend.h>
#include <TLine.h>
#include <TMath.h>

#include "CascadeUtils.h"
enum stats
{
    kMean,
    kRMS,
    kMax,
    kNumStats
};

/**
 * @brief Computes the yield deviation between a variation and default histogram
 *
 * @param hVar Input variation histogram to compare
 * @param hDef Input default/reference histogram
 * @param hYieldDev Output histogram to store yield deviations
 * @param ptBin pT bin number being processed
 * @param iVar Variation number/index
 * @param fDebug Enable/disable debug output
 *
 * @return Double_t Returns the computed yield deviation value
 */
Double_t ComputeYieldDev(TH1 *hVar, TH1 *hDef, TH1 *hYieldDev, Int_t ptBin, Int_t iVar, Bool_t fDebug);

/**
 * Calculates the uncertainty from yield variations
 *
 * @param hYieldDev Histogram containing yield variations
 * @param stat Statistical method to calculate uncertainty (kMean, kRMS, kMax)
 * @return Double_t Calculated uncertainty value
 *
 * This function processes a histogram of yield variations to compute
 * the associated uncertainty. The calculation method is determined
 * by the stat parameter (enum).
 */
Double_t GetUncertainty(TH1 *hYieldDev, Int_t stat = kMean);

/**
 * @brief Draws a histogram and saves it as an image
 * @param hist The histogram to be drawn
 * @param outputFolder The folder path where the output file will be saved
 * @param imageFolder The subfolder path where the image will be saved
 * @param imageFormat The format of the output image (e.g., "png", "pdf")
 * @param drawLines Boolean flag to determine if yMax, yMean and yRMS lines should be drawn on the histogram
 */
void DrawAndSaveImage(TH1 *hist, TString outputFolder, TString imageFolder, TString imageFormat, Bool_t drawLines);

/**
 * @brief Calculate systematic uncertainty from signal extraction variations
 *
 * This function calculates the systematic uncertainty due to different signal extraction methods.
 * It processes variations of signal extraction methods and compares them to a default method to
 * determine systematic uncertainties for Xi and Omega particle spectra.
 *
 * @param inputPath Path to input files containing signal extraction variations
 * @param effCorrInputFilePrefix Prefix for efficiency corrected input files
 * @param outputFileName Output ROOT file path to store results
 * @param outputFolder Output folder path for plots and results
 * @param fDebug Enable debug output
 * @param saveImages Flag to save plots as image files
 * @param imageFormat Format for saved images (e.g. "png")
 * @param verbosity ROOT verbosity level
 *
 * @return 0 on success, non-zero on error
 *
 * The function:
 * 1. Reads efficiency corrected spectra for default cuts
 * 2. Reads spectra for various signal extraction methods
 * 3. Calculates yield deviations between variations and default
 * 4. Computes systematic uncertainties from the deviations
 * 5. Creates output histograms for:
 *    - Yield deviations per pT and multiplicity bin
 *    - Final systematic uncertainties
 * 6. Saves results to ROOT file and optionally as images
 *
 * Processes both Xi and Omega particles, separately for:
 * - Particles and anti-particles
 * - Combined charge states
 * - Different multiplicity bins
 * - Different pT bins
 */
int SysSignalExtraction(
    TString inputPath = "/var/home/ishaan/Work/git/analysis/results/RandomVars/SysVars_SignalExtraction/200225_SysSigExt_effCorr",
    TString effCorrInputFilePrefix = "200225_SysSigExt",
    TString outputFileName = "/var/home/ishaan/Work/git/analysis/results/RandomVars/250225_Systematics_SigExt/BackgroundFitRange/250225_SystematicUncertainty_SigExt.root",
    TString outputFolder = "/var/home/ishaan/Work/git/analysis/results/RandomVars/250225_Systematics_SigExt/BackgroundFitRange",
    Bool_t fDebug = kFALSE,
    Bool_t saveImages = kTRUE,
    TString imageFormat = "png",
    Int_t verbosity = kInfo)
{

    gErrorIgnoreLevel = verbosity;
    gROOT->SetBatch(kTRUE);
    gStyle->SetOptStat("n");

    // remove ownership of objects from file so we can delete the file ptr
    TH1::AddDirectory(kFALSE);
    TH1::SetDefaultSumw2(kTRUE);

    outputFolder = SetOutputFolder(outputFolder);

    TString varName[] = /*{"GP2_def", "DGP3_def"};*/ /*{"DGP2_3sig", "DGP2_5sig"}*/ {"DGP2_10bg", "DGP2_15bg"};
    UInt_t nVar = sizeof(varName) / sizeof(varName[0]);

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

    TH1D *sysSigExt_xim;                     // mult integrated systematic uncertainty from SigExt
    TH1D *sysSigExt_xip;                     // mult integrated systematic uncertainty from SigExt
    TH1D *sysSigExt_omm;                     // mult integrated systematic uncertainty from SigExt
    TH1D *sysSigExt_omp;                     // mult integrated systematic uncertainty from SigExt
    TH1D *sysSigExt_xiC;                     // mult integrated systematic uncertainty from SigExt
    TH1D *sysSigExt_omC;                     // mult integrated systematic uncertainty from SigExt
    TH1D *sysSigExt_xim_mult[fNmultbins_Xi]; // mult binned systematic uncertainty from SigExt
    TH1D *sysSigExt_xip_mult[fNmultbins_Xi]; // mult binned systematic uncertainty from SigExt
    TH1D *sysSigExt_omm_mult[fNmultbins_Om]; // mult binned systematic uncertainty from SigExt
    TH1D *sysSigExt_omp_mult[fNmultbins_Om]; // mult binned systematic uncertainty from SigExt
    TH1D *sysSigExt_xiC_mult[fNmultbins_Xi]; // mult binned systematic uncertainty from SigExt
    TH1D *sysSigExt_omC_mult[fNmultbins_Om]; // mult binned systematic uncertainty from SigExt

    /// Getting default cut histograms:
    TString varNameDef = "DGP2_def";
    TString inputFileName = TString::Format("%s/%s_%s_effCorr.root", inputPath.Data(), effCorrInputFilePrefix.Data(), varNameDef.Data());

    TFile *defInputFile = OpenFile(inputFileName);
    Info("SysSigExt: defInput", "Getting default cut histograms from '%s'", inputFileName.Data());
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
        outputFileName = TString::Format("%s/SysUncertainty_SigExt.root", outputFolder.Data());
    }

    Info("SysSigExt: I/O", "Default input successful. Creating output file '%s'", outputFileName.Data());

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
    outputFile->mkdir("dirSysSigExt_xim_mult");
    outputFile->mkdir("dirSysSigExt_xip_mult");
    outputFile->mkdir("dirSysSigExt_omm_mult");
    outputFile->mkdir("dirSysSigExt_omp_mult");
    outputFile->mkdir("dirSysSigExt_xiC_mult");
    outputFile->mkdir("dirSysSigExt_omC_mult");
    // Output set.

    // Generate output objects
    // Final output object - Systematic uncertainty from SigExt:
    sysSigExt_xip = new TH1D("sysSigExt_xip", "#Xi^{+}: Systematic uncertainty from SigExt: Mult 0-100%;#it{p}_{T} (GeV/#it{c});Uncertainty", fNptbins_Xi, fPtbins_Xi);
    sysSigExt_xim = new TH1D("sysSigExt_xim", "#Xi^{-}: Systematic uncertainty from SigExt: Mult 0-100%;#it{p}_{T} (GeV/#it{c});Uncertainty", fNptbins_Xi, fPtbins_Xi);
    sysSigExt_omp = new TH1D("sysSigExt_omp", "#Omega^{+}: Systematic uncertainty from SigExt: Mult 0-100%;#it{p}_{T} (GeV/#it{c});Uncertainty", fNptbins_Om, fPtbins_Om);
    sysSigExt_omm = new TH1D("sysSigExt_omm", "#Omega^{-}: Systematic uncertainty from SigExt: Mult 0-100%;#it{p}_{T} (GeV/#it{c});Uncertainty", fNptbins_Om, fPtbins_Om);
    sysSigExt_xiC = new TH1D("sysSigExt_xiC", "#Xi^{+} + #Xi^{-}: Systematic uncertainty from SigExt: Mult 0-100%;#it{p}_{T} (GeV/#it{c});Uncertainty", fNptbins_Xi, fPtbins_Xi);
    sysSigExt_omC = new TH1D("sysSigExt_omC", "#Omega^{+} + #Omega^{-}: Systematic uncertainty from SigExt: Mult 0-100%;#it{p}_{T} (GeV/#it{c});Uncertainty", fNptbins_Om, fPtbins_Om);
    for (Int_t multBinXi = 0; multBinXi < fNmultbins_Xi; multBinXi++)
    {
        sysSigExt_xip_mult[multBinXi] = new TH1D(TString::Format("sysSigExt_xip_mult[%d]", multBinXi), TString::Format("#Xi^{+}: Systematic uncertainty from SigExt: Mult %.0f-%.0f%%;#it{p}_{T} (GeV/#it{c});Uncertainty", fMultbins_Xi[multBinXi], fMultbins_Xi[multBinXi + 1]), fNptbins_Xi, fPtbins_Xi);
        sysSigExt_xim_mult[multBinXi] = new TH1D(TString::Format("sysSigExt_xim_mult[%d]", multBinXi), TString::Format("#Xi^{-}: Systematic uncertainty from SigExt: Mult %.0f-%.0f%%;#it{p}_{T} (GeV/#it{c});Uncertainty", fMultbins_Xi[multBinXi], fMultbins_Xi[multBinXi + 1]), fNptbins_Xi, fPtbins_Xi);
        sysSigExt_xiC_mult[multBinXi] = new TH1D(TString::Format("sysSigExt_xiC_mult[%d]", multBinXi), TString::Format("#Xi^{+} + #Xi^{-}: Systematic uncertainty from SigExt: Mult %.0f-%.0f%%;#it{p}_{T} (GeV/#it{c});Uncertainty", fMultbins_Xi[multBinXi], fMultbins_Xi[multBinXi + 1]), fNptbins_Xi, fPtbins_Xi);
    }
    for (Int_t multBinOm = 0; multBinOm < fNmultbins_Om; multBinOm++)
    {
        sysSigExt_omp_mult[multBinOm] = new TH1D(TString::Format("sysSigExt_omp_mult[%d]", multBinOm), TString::Format("#Omega^{+}: Systematic uncertainty from SigExt: Mult %.0f-%.0f%%;#it{p}_{T} (GeV/#it{c});Uncertainty", fMultbins_Om[multBinOm], fMultbins_Om[multBinOm + 1]), fNptbins_Om, fPtbins_Om);
        sysSigExt_omm_mult[multBinOm] = new TH1D(TString::Format("sysSigExt_omm_mult[%d]", multBinOm), TString::Format("#Omega^{-}: Systematic uncertainty from SigExt: Mult %.0f-%.0f%%;#it{p}_{T} (GeV/#it{c});Uncertainty", fMultbins_Om[multBinOm], fMultbins_Om[multBinOm + 1]), fNptbins_Om, fPtbins_Om);
        sysSigExt_omC_mult[multBinOm] = new TH1D(TString::Format("sysSigExt_omC_mult[%d]", multBinOm), TString::Format("#Omega^{+} + #Omega^{-}: Systematic uncertainty from SigExt: Mult %.0f-%.0f%%;#it{p}_{T} (GeV/#it{c});Uncertainty", fMultbins_Om[multBinOm], fMultbins_Om[multBinOm + 1]), fNptbins_Om, fPtbins_Om);
    }

    /// Reference output objects - yield deviation from default cut:
    for (Int_t ptBinXi = 0; ptBinXi < fNptbins_Xi; ptBinXi++)
    {
        varDefYieldDev_xip_pt[ptBinXi] = new TH1D(TString::Format("varDefYieldDev_xip_pt[%d]", ptBinXi), TString::Format("#Xi^{+}: #it{p}_{T}: <%.1f,%.1f>, Mult: 0-100%%;Vars;|Y_{sys}/Y_{def} - 1|", fPtbins_Xi[ptBinXi], fPtbins_Xi[ptBinXi + 1]), nVar, 0, nVar);
        varDefYieldDev_xim_pt[ptBinXi] = new TH1D(TString::Format("varDefYieldDev_xim_pt[%d]", ptBinXi), TString::Format("#Xi^{-}: #it{p}_{T}: <%.1f,%.1f>, Mult: 0-100%%;Vars;|Y_{sys}/Y_{def} - 1|", fPtbins_Xi[ptBinXi], fPtbins_Xi[ptBinXi + 1]), nVar, 0, nVar);
        varDefYieldDev_xiC_pt[ptBinXi] = new TH1D(TString::Format("varDefYieldDev_xiC_pt[%d]", ptBinXi), TString::Format("#Xi^{+} + #Xi^{-}: #it{p}_{T}: <%.1f,%.1f>, Mult: 0-100%%;Vars;|Y_{sys}/Y_{def} - 1|", fPtbins_Xi[ptBinXi], fPtbins_Xi[ptBinXi + 1]), nVar, 0, nVar);
        for (Int_t multBinXi = 0; multBinXi < fNmultbins_Xi; multBinXi++)
        {
            varDefYieldDev_xip_pt_mult[ptBinXi][multBinXi] = new TH1D(TString::Format("varDefYieldDev_xip_pt_mult[%d][%d]", ptBinXi, multBinXi), TString::Format("#Xi^{+}: #it{p}_{T}: <%.1f,%.1f>, Mult: %.0f-%.0f%%;Vars;|Y_{sys}/Y_{def} - 1|", fPtbins_Xi[ptBinXi], fPtbins_Xi[ptBinXi + 1], fMultbins_Xi[multBinXi], fMultbins_Xi[multBinXi + 1]), nVar, 0, nVar);
            varDefYieldDev_xim_pt_mult[ptBinXi][multBinXi] = new TH1D(TString::Format("varDefYieldDev_xim_pt_mult[%d][%d]", ptBinXi, multBinXi), TString::Format("#Xi^{-}: #it{p}_{T}: <%.1f,%.1f>, Mult: %.0f-%.0f%%;Vars;|Y_{sys}/Y_{def} - 1|", fPtbins_Xi[ptBinXi], fPtbins_Xi[ptBinXi + 1], fMultbins_Xi[multBinXi], fMultbins_Xi[multBinXi + 1]), nVar, 0, nVar);
            varDefYieldDev_xiC_pt_mult[ptBinXi][multBinXi] = new TH1D(TString::Format("varDefYieldDev_xiC_pt_mult[%d][%d]", ptBinXi, multBinXi), TString::Format("#Xi^{+} + #Xi^{-}: #it{p}_{T}: <%.1f,%.1f>, Mult: %.0f-%.0f%%;Vars;|Y_{sys}/Y_{def} - 1|", fPtbins_Xi[ptBinXi], fPtbins_Xi[ptBinXi + 1], fMultbins_Xi[multBinXi], fMultbins_Xi[multBinXi + 1]), nVar, 0, nVar);
        }
    }
    for (Int_t ptBinOm = 0; ptBinOm < fNptbins_Om; ptBinOm++)
    {
        varDefYieldDev_omp_pt[ptBinOm] = new TH1D(TString::Format("varDefYieldDev_omp_pt[%d]", ptBinOm), TString::Format("#Omega^{+}: #it{p}_{T}: <%.1f,%.1f>, Mult: 0-100%%;Vars;|Y_{sys}/Y_{def} - 1|", fPtbins_Om[ptBinOm], fPtbins_Om[ptBinOm + 1]), nVar, 0, nVar);
        varDefYieldDev_omm_pt[ptBinOm] = new TH1D(TString::Format("varDefYieldDev_omm_pt[%d]", ptBinOm), TString::Format("#Omega^{-}: #it{p}_{T}: <%.1f,%.1f>, Mult: 0-100%%;Vars;|Y_{sys}/Y_{def} - 1|", fPtbins_Om[ptBinOm], fPtbins_Om[ptBinOm + 1]), nVar, 0, nVar);
        varDefYieldDev_omC_pt[ptBinOm] = new TH1D(TString::Format("varDefYieldDev_omC_pt[%d]", ptBinOm), TString::Format("#Omega^{+} + #Omega^{-}: #it{p}_{T}: <%.1f,%.1f>, Mult: 0-100%%;Vars;|Y_{sys}/Y_{def} - 1|", fPtbins_Om[ptBinOm], fPtbins_Om[ptBinOm + 1]), nVar, 0, nVar);
        for (Int_t multBinOm = 0; multBinOm < fNmultbins_Om; multBinOm++)
        {
            varDefYieldDev_omp_pt_mult[ptBinOm][multBinOm] = new TH1D(TString::Format("varDefYieldDev_omp_pt_mult[%d][%d]", ptBinOm, multBinOm), TString::Format("#Omega^{+}: #it{p}_{T}: <%.1f,%.1f>, Mult: %.0f-%.0f%%;Vars;|Y_{sys}/Y_{def} - 1|", fPtbins_Om[ptBinOm], fPtbins_Om[ptBinOm + 1], fMultbins_Om[multBinOm], fMultbins_Om[multBinOm + 1]), nVar, 0, nVar);
            varDefYieldDev_omm_pt_mult[ptBinOm][multBinOm] = new TH1D(TString::Format("varDefYieldDev_omm_pt_mult[%d][%d]", ptBinOm, multBinOm), TString::Format("#Omega^{-}: #it{p}_{T}: <%.1f,%.1f>, Mult: %.0f-%.0f%%;Vars;|Y_{sys}/Y_{def} - 1|", fPtbins_Om[ptBinOm], fPtbins_Om[ptBinOm + 1], fMultbins_Om[multBinOm], fMultbins_Om[multBinOm + 1]), nVar, 0, nVar);
            varDefYieldDev_omC_pt_mult[ptBinOm][multBinOm] = new TH1D(TString::Format("varDefYieldDev_omC_pt_mult[%d][%d]", ptBinOm, multBinOm), TString::Format("#Omega^{+} + #Omega^{-}: #it{p}_{T}: <%.1f,%.1f>, Mult: %.0f-%.0f%%;Vars;|Y_{sys}/Y_{def} - 1|", fPtbins_Om[ptBinOm], fPtbins_Om[ptBinOm + 1], fMultbins_Om[multBinOm], fMultbins_Om[multBinOm + 1]), nVar, 0, nVar);
        }
    }
    Info("SysSigExt: I/O", "Output objects and output file created successfully.");

    /// Getting varied cut histograms:
    for (Int_t iVar = 0; iVar < nVar; iVar++)
    {
        inputFileName = TString::Format("%s/%s_%s_effCorr.root", inputPath.Data(), effCorrInputFilePrefix.Data(), varName[iVar].Data());

        TFile *varInputFile = OpenFile(inputFileName);
        Info("SysSigExt: varInput", "Getting varied cut histograms from '%s'", inputFileName.Data());
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
            varDefYieldDev_xip_pt[ptBinXi]->GetXaxis()->SetBinLabel(iVar + 1, varName[iVar].Data());
            varDefYieldDev_xim_pt[ptBinXi]->GetXaxis()->SetBinLabel(iVar + 1, varName[iVar].Data());
            varDefYieldDev_xiC_pt[ptBinXi]->GetXaxis()->SetBinLabel(iVar + 1, varName[iVar].Data());
            ComputeYieldDev(var_effCorrPt_xip, def_effCorrPt_xip, varDefYieldDev_xip_pt[ptBinXi], ptBinXi + 1, iVar, fDebug);
            ComputeYieldDev(var_effCorrPt_xim, def_effCorrPt_xim, varDefYieldDev_xim_pt[ptBinXi], ptBinXi + 1, iVar, fDebug);
            ComputeYieldDev(var_effCorrPt_xiC, def_effCorrPt_xiC, varDefYieldDev_xiC_pt[ptBinXi], ptBinXi + 1, iVar, fDebug);

            for (Int_t multBinXi = 0; multBinXi < fNmultbins_Xi; multBinXi++)
            {
                varDefYieldDev_xip_pt_mult[ptBinXi][multBinXi]->GetXaxis()->SetBinLabel(iVar + 1, varName[iVar].Data());
                varDefYieldDev_xim_pt_mult[ptBinXi][multBinXi]->GetXaxis()->SetBinLabel(iVar + 1, varName[iVar].Data());
                varDefYieldDev_xiC_pt_mult[ptBinXi][multBinXi]->GetXaxis()->SetBinLabel(iVar + 1, varName[iVar].Data());
                ComputeYieldDev(var_effCorrPt_xip_mult[multBinXi], def_effCorrPt_xip_mult[multBinXi], varDefYieldDev_xip_pt_mult[ptBinXi][multBinXi], ptBinXi + 1, iVar, fDebug);
                ComputeYieldDev(var_effCorrPt_xim_mult[multBinXi], def_effCorrPt_xim_mult[multBinXi], varDefYieldDev_xim_pt_mult[ptBinXi][multBinXi], ptBinXi + 1, iVar, fDebug);
                ComputeYieldDev(var_effCorrPt_xiC_mult[multBinXi], def_effCorrPt_xiC_mult[multBinXi], varDefYieldDev_xiC_pt_mult[ptBinXi][multBinXi], ptBinXi + 1, iVar, fDebug);
            }
        }
        for (Int_t ptBinOm = 0; ptBinOm < fNptbins_Om; ptBinOm++)
        {
            varDefYieldDev_omp_pt[ptBinOm]->GetXaxis()->SetBinLabel(iVar + 1, varName[iVar].Data());
            varDefYieldDev_omm_pt[ptBinOm]->GetXaxis()->SetBinLabel(iVar + 1, varName[iVar].Data());
            varDefYieldDev_omC_pt[ptBinOm]->GetXaxis()->SetBinLabel(iVar + 1, varName[iVar].Data());
            ComputeYieldDev(var_effCorrPt_omp, def_effCorrPt_omp, varDefYieldDev_omp_pt[ptBinOm], ptBinOm + 1, iVar, fDebug);
            ComputeYieldDev(var_effCorrPt_omm, def_effCorrPt_omm, varDefYieldDev_omm_pt[ptBinOm], ptBinOm + 1, iVar, fDebug);
            ComputeYieldDev(var_effCorrPt_omC, def_effCorrPt_omC, varDefYieldDev_omC_pt[ptBinOm], ptBinOm + 1, iVar, fDebug);

            for (Int_t multBinOm = 0; multBinOm < fNmultbins_Om; multBinOm++)
            {
                varDefYieldDev_omp_pt_mult[ptBinOm][multBinOm]->GetXaxis()->SetBinLabel(iVar + 1, varName[iVar].Data());
                varDefYieldDev_omm_pt_mult[ptBinOm][multBinOm]->GetXaxis()->SetBinLabel(iVar + 1, varName[iVar].Data());
                varDefYieldDev_omC_pt_mult[ptBinOm][multBinOm]->GetXaxis()->SetBinLabel(iVar + 1, varName[iVar].Data());
                ComputeYieldDev(var_effCorrPt_omp_mult[multBinOm], def_effCorrPt_omp_mult[multBinOm], varDefYieldDev_omp_pt_mult[ptBinOm][multBinOm], ptBinOm + 1, iVar, fDebug);
                ComputeYieldDev(var_effCorrPt_omm_mult[multBinOm], def_effCorrPt_omm_mult[multBinOm], varDefYieldDev_omm_pt_mult[ptBinOm][multBinOm], ptBinOm + 1, iVar, fDebug);
                ComputeYieldDev(var_effCorrPt_omC_mult[multBinOm], def_effCorrPt_omC_mult[multBinOm], varDefYieldDev_omC_pt_mult[ptBinOm][multBinOm], ptBinOm + 1, iVar, fDebug);
            }
        }
    }
    Info("SysSigExt: varInput", "Variation input successful. Yield deviation computed.");

    /// Get systematic uncertainty from yield deviation histograms based on selected statistical method in 'stat' enum:
    Info("SysSigExt: Xi", "Calculating systematic uncertainty from yield deviation histograms");
    for (Int_t ptBinXi = 0; ptBinXi < fNptbins_Xi; ptBinXi++)
    {
        sysSigExt_xip->SetBinContent(ptBinXi + 1, GetUncertainty(varDefYieldDev_xip_pt[ptBinXi]));
        sysSigExt_xim->SetBinContent(ptBinXi + 1, GetUncertainty(varDefYieldDev_xim_pt[ptBinXi]));
        sysSigExt_xiC->SetBinContent(ptBinXi + 1, GetUncertainty(varDefYieldDev_xiC_pt[ptBinXi]));

        for (Int_t multBinXi = 0; multBinXi < fNmultbins_Xi; multBinXi++)
        {
            sysSigExt_xip_mult[multBinXi]->SetBinContent(ptBinXi + 1, GetUncertainty(varDefYieldDev_xip_pt_mult[ptBinXi][multBinXi]));
            sysSigExt_xim_mult[multBinXi]->SetBinContent(ptBinXi + 1, GetUncertainty(varDefYieldDev_xim_pt_mult[ptBinXi][multBinXi]));
            sysSigExt_xiC_mult[multBinXi]->SetBinContent(ptBinXi + 1, GetUncertainty(varDefYieldDev_xiC_pt_mult[ptBinXi][multBinXi]));
        }
    }
    Info("SysSigExt: Omega", "Calculating systematic uncertainty from yield deviation histograms");
    for (Int_t ptBinOm = 0; ptBinOm < fNptbins_Om; ptBinOm++)
    {
        sysSigExt_omp->SetBinContent(ptBinOm + 1, GetUncertainty(varDefYieldDev_omp_pt[ptBinOm]));
        sysSigExt_omm->SetBinContent(ptBinOm + 1, GetUncertainty(varDefYieldDev_omm_pt[ptBinOm]));
        sysSigExt_omC->SetBinContent(ptBinOm + 1, GetUncertainty(varDefYieldDev_omC_pt[ptBinOm]));

        for (Int_t multBinOm = 0; multBinOm < fNmultbins_Om; multBinOm++)
        {
            sysSigExt_omp_mult[multBinOm]->SetBinContent(ptBinOm + 1, GetUncertainty(varDefYieldDev_omp_pt_mult[ptBinOm][multBinOm]));
            sysSigExt_omm_mult[multBinOm]->SetBinContent(ptBinOm + 1, GetUncertainty(varDefYieldDev_omm_pt_mult[ptBinOm][multBinOm]));
            sysSigExt_omC_mult[multBinOm]->SetBinContent(ptBinOm + 1, GetUncertainty(varDefYieldDev_omC_pt_mult[ptBinOm][multBinOm]));
        }
    }

    /// Write output objects to file and save images:
    outputFile->cd("dirSysSigExt_xip_mult");
    sysSigExt_xip->Write();
    outputFile->cd("dirSysSigExt_xim_mult");
    sysSigExt_xim->Write();
    outputFile->cd("dirSysSigExt_xiC_mult");
    sysSigExt_xiC->Write();
    outputFile->cd("dirSysSigExt_omp_mult");
    sysSigExt_omp->Write();
    outputFile->cd("dirSysSigExt_omm_mult");
    sysSigExt_omm->Write();
    outputFile->cd("dirSysSigExt_omC_mult");
    sysSigExt_omC->Write();
    if (saveImages)
    {
        DrawAndSaveImage(sysSigExt_xip, outputFolder, "SysSigExt_xip", imageFormat, kFALSE);
        DrawAndSaveImage(sysSigExt_xim, outputFolder, "SysSigExt_xim", imageFormat, kFALSE);
        DrawAndSaveImage(sysSigExt_xiC, outputFolder, "SysSigExt_xiC", imageFormat, kFALSE);
        DrawAndSaveImage(sysSigExt_omp, outputFolder, "SysSigExt_omp", imageFormat, kFALSE);
        DrawAndSaveImage(sysSigExt_omm, outputFolder, "SysSigExt_omm", imageFormat, kFALSE);
        DrawAndSaveImage(sysSigExt_omC, outputFolder, "SysSigExt_omC", imageFormat, kFALSE);
    }
    for (Int_t multBinXi = 0; multBinXi < fNmultbins_Xi; multBinXi++)
    {
        outputFile->cd("dirSysSigExt_xip_mult");
        sysSigExt_xip_mult[multBinXi]->Write();
        outputFile->cd("dirSysSigExt_xim_mult");
        sysSigExt_xim_mult[multBinXi]->Write();
        outputFile->cd("dirSysSigExt_xiC_mult");
        sysSigExt_xiC_mult[multBinXi]->Write();
        if (saveImages)
        {
            DrawAndSaveImage(sysSigExt_xip_mult[multBinXi], outputFolder, "SysSigExt_xip", imageFormat, kFALSE);
            DrawAndSaveImage(sysSigExt_xim_mult[multBinXi], outputFolder, "SysSigExt_xim", imageFormat, kFALSE);
            DrawAndSaveImage(sysSigExt_xiC_mult[multBinXi], outputFolder, "SysSigExt_xiC", imageFormat, kFALSE);
        }
    }
    for (Int_t multBinOm = 0; multBinOm < fNmultbins_Om; multBinOm++)
    {
        outputFile->cd("dirSysSigExt_omp_mult");
        sysSigExt_omp_mult[multBinOm]->Write();
        outputFile->cd("dirSysSigExt_omm_mult");
        sysSigExt_omm_mult[multBinOm]->Write();
        outputFile->cd("dirSysSigExt_omC_mult");
        sysSigExt_omC_mult[multBinOm]->Write();
        if (saveImages)
        {
            DrawAndSaveImage(sysSigExt_omp_mult[multBinOm], outputFolder, "SysSigExt_omp", imageFormat, kFALSE);
            DrawAndSaveImage(sysSigExt_omm_mult[multBinOm], outputFolder, "SysSigExt_omm", imageFormat, kFALSE);
            DrawAndSaveImage(sysSigExt_omC_mult[multBinOm], outputFolder, "SysSigExt_omC", imageFormat, kFALSE);
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
            DrawAndSaveImage(varDefYieldDev_xip_pt[ptBinXi], outputFolder, "YieldDev_xip_pt", imageFormat, kTRUE);
            DrawAndSaveImage(varDefYieldDev_xim_pt[ptBinXi], outputFolder, "YieldDev_xim_pt", imageFormat, kTRUE);
            DrawAndSaveImage(varDefYieldDev_xiC_pt[ptBinXi], outputFolder, "YieldDev_xiC_pt", imageFormat, kTRUE);
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
                DrawAndSaveImage(varDefYieldDev_xip_pt_mult[ptBinXi][multBinXi], outputFolder, "YieldDev_xip_pt_mult", imageFormat, kTRUE);
                DrawAndSaveImage(varDefYieldDev_xim_pt_mult[ptBinXi][multBinXi], outputFolder, "YieldDev_xim_pt_mult", imageFormat, kTRUE);
                DrawAndSaveImage(varDefYieldDev_xiC_pt_mult[ptBinXi][multBinXi], outputFolder, "YieldDev_xiC_pt_mult", imageFormat, kTRUE);
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
            DrawAndSaveImage(varDefYieldDev_omp_pt[ptBinOm], outputFolder, "YieldDev_omp_pt", imageFormat, kTRUE);
            DrawAndSaveImage(varDefYieldDev_omm_pt[ptBinOm], outputFolder, "YieldDev_omm_pt", imageFormat, kTRUE);
            DrawAndSaveImage(varDefYieldDev_omC_pt[ptBinOm], outputFolder, "YieldDev_omC_pt", imageFormat, kTRUE);
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
                DrawAndSaveImage(varDefYieldDev_omp_pt_mult[ptBinOm][multBinOm], outputFolder, "YieldDev_omp_pt_mult", imageFormat, kTRUE);
                DrawAndSaveImage(varDefYieldDev_omm_pt_mult[ptBinOm][multBinOm], outputFolder, "YieldDev_omm_pt_mult", imageFormat, kTRUE);
                DrawAndSaveImage(varDefYieldDev_omC_pt_mult[ptBinOm][multBinOm], outputFolder, "YieldDev_omC_pt_mult", imageFormat, kTRUE);
            }
        }
    }
    Info("SysSigExt: I/O", "Output objects written to file '%s' successfully.", outputFileName.Data());

    delete outputFile;

    return 0;
}

Double_t GetUncertainty(TH1 *hYieldDev, Int_t stat)
{
    UInt_t nVar = hYieldDev->GetNbinsX();
    Double_t yArray[nVar];
    Double_t uncertainty = 0.0;
    for (int iBin = 0; iBin < nVar; iBin++)
    {
        yArray[iBin] = hYieldDev->GetBinContent(iBin + 1);
    }

    if (stat == kMean)
        uncertainty = TMath::Mean(nVar, yArray);
    else if (stat == kRMS)
        uncertainty = CalculateRMS(nVar, yArray);
    else if (stat == kMax)
        uncertainty = TMath::MaxElement(nVar, yArray);

    return uncertainty;
}

Double_t ComputeYieldDev(TH1 *hVar, TH1 *hDef, TH1 *hYieldDev, Int_t ptBin, Int_t iVar, Bool_t fDebug)
{
    Double_t varVal = hVar->GetBinContent(ptBin);
    Double_t defVal = hDef->GetBinContent(ptBin);
    Double_t yieldDev = 0.0;

    if (defVal > 1e-12)
    {
        yieldDev = (varVal / defVal) - 1;

        if (fDebug)
            Info("ComputeYieldDev", "iVar %d: Filling %s: ptBin: %d, varVal: %f, defVal: %f, yieldDev: %f", iVar, hYieldDev->GetName(), ptBin, varVal, defVal, yieldDev);
        yieldDev = TMath::Abs(yieldDev);
        hYieldDev->SetBinContent(iVar + 1, yieldDev);
    }
    else
    {
        Warning("ComputeYieldDev", "iVar %d: Skipping %s: ptBin: %d, varVal: %f, defVal: %f. Division by zero.", iVar, hYieldDev->GetName(), ptBin, varVal, defVal);
    }

    return yieldDev;
}

void DrawAndSaveImage(TH1 *hist, TString outputFolder, TString imageFolder, TString imageFormat, Bool_t drawLines)
{
    gStyle->SetPaintTextFormat("1.3f");
    hist->SetMarkerStyle(kFullCircle);

    TCanvas *cDraw = new TCanvas(hist->GetName(), hist->GetTitle(), 1920, 1080);
    cDraw->cd();

    if (drawLines)
    {
        UInt_t nVar = hist->GetNbinsX();
        Double_t xMin = hist->GetXaxis()->GetXmin();
        Double_t xMax = hist->GetXaxis()->GetXmax();
        Double_t yArray[nVar];
        for (int iBin = 0; iBin < nVar; iBin++)
        {
            yArray[iBin] = hist->GetBinContent(iBin + 1);
        }
        Double_t yMean = TMath::Mean(nVar, yArray);
        Double_t yRMS = CalculateRMS(nVar, yArray);
        Double_t yMax = TMath::MaxElement(nVar, yArray);

        TLine *lMaxLine = new TLine(xMin, yMax, xMax, yMax);
        lMaxLine->SetLineColor(kRed);
        lMaxLine->SetLineWidth(2);
        lMaxLine->SetLineStyle(2);
        TLine *lMeanLine = new TLine(xMin, yMean, xMax, yMean);
        lMeanLine->SetLineColor(kGreen);
        lMeanLine->SetLineWidth(2);
        lMeanLine->SetLineStyle(9);
        TLine *lRmsLine = new TLine(xMin, yRMS, xMax, yRMS);
        lRmsLine->SetLineColor(kViolet);
        lRmsLine->SetLineWidth(2);
        lRmsLine->SetLineStyle(10);

        auto legend = new TLegend(0.1, 0.79, 0.21, 0.9);
        legend->SetTextSize(0.02);
        legend->AddEntry(lMaxLine, TString::Format("Max = %.5f", yMax), "l");
        legend->AddEntry(lMeanLine, TString::Format("Mean = %.5f", yMean), "l");
        legend->AddEntry(lRmsLine, TString::Format("RMS = %.5f", yRMS), "l");

        hist->SetMaximum(TMath::Max(yMax, yRMS) * 1.2);
        hist->SetMinimum(TMath::Min(hist->GetMinimum(), yRMS) * 0.8);
        hist->Draw("TEXT00");
        lMaxLine->Draw("SAME");
        lMeanLine->Draw("SAME");
        lRmsLine->Draw("SAME");
        legend->Draw("SAME");

        cDraw->Modified();
        cDraw->Update();

        SaveImage(outputFolder, imageFolder, hist->GetName(), imageFormat, cDraw);

        delete lMaxLine;
        delete lMeanLine;
        delete lRmsLine;
        delete legend;
    }
    else
    {
        hist->Draw("TEXT00");
        SaveImage(outputFolder, imageFolder, hist->GetName(), imageFormat, cDraw);
    }

    delete cDraw;
}