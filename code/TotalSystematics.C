#include <TROOT.h>
#include <TStyle.h>
#include <TLegend.h>
#include <Math/MinimizerOptions.h>
#include <Math/IntegratorOptions.h>
#include <TLine.h>

#include "CascadeUtils.h"
#include "ExtraUtils.C"

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
void PaintStackOverlap(TCanvas &c, THStack &hs, Bool_t setLogY = kTRUE, Bool_t fit = kFALSE, TString yAxisTitle = "#frac{1}{2#pi #it{p}_{T} #it{N}_{ev}} #frac{d^{2}#it{N}}{d#it{p}_{T}d#it{y}} [(GeV/#it{c})^{-2}] ", TString xAxisTitle = "#it{p}_{T} (GeV/#it{c})");

// --- Levy-Tsallis function for dN/dpT ---
// par[0] = dN/dy (normalization)
// par[1] = T (temperature)
// par[2] = n (exponent)
// par[3] = mass (particle mass, fixed during fit)
Double_t levyTsallisFuncPt(Double_t *x, Double_t *par)
{
    Double_t pt = x[0];

    Double_t dNdy = par[0];
    Double_t T_par = par[1]; // Use T_par to avoid conflict with TMath::T() if ever used
    Double_t n_par = par[2];
    Double_t mass = par[3];

    // Parameter sanity checks
    if (T_par <= 0 || n_par <= 2.0)
    { // n > 2 is strictly required for this form
        return 0.0;
    }

    Double_t mt = TMath::Sqrt(pt * pt + mass * mass);

    Double_t term_in_denom = n_par * T_par + mass * (n_par - 2.0);
    if (term_in_denom == 0)
    { // Avoid division by zero
        return 0.0;
    }
    Double_t common_factor_coeff = (n_par - 1.0) * (n_par - 2.0) / (n_par * T_par * term_in_denom);

    Double_t base_for_power = 1.0 + (mt - mass) / (n_par * T_par);

    // Handle potential issues with TMath::Power if base_for_power is zero or negative
    // Though with mt >= mass, T_par > 0, n_par > 0, base_for_power should be >= 1.0
    if (base_for_power <= 0 && -n_par != static_cast<int>(-n_par))
    {               // if base is <=0 and exponent is not integer
        return 0.0; // Avoid NaN from pow for non-integer exponent with negative base
    }
    if (base_for_power == 0 && -n_par < 0)
    { // Avoid division by zero if base is 0 and exponent is negative
        return 0.0;
    }

    Double_t tsallis_power_term = TMath::Power(base_for_power, -n_par);

    return pt * dNdy * common_factor_coeff * tsallis_power_term;
};

// --- Helper function to perform the fit on a single histogram ---
void LevyTsallisFit(TH1 *hist, Int_t particle = -1, Int_t multBin = -1)
{
    Double_t xmin = 0.0, xmax = 10.0;
    Double_t mass, fitMin, fitMax = 0.0;
    char *particleName = "";

    if (!hist || hist->GetEntries() == 0)
    {
        printf("FitSkip: %s spectrum is null or empty.\n", particleName);
        return;
    }

    if (particle == kXi)
    {
        mass = fMass_Xi;
        fitMin = fPtbins_Xi[0];           // Use Xi pT bin minimum
        fitMax = fPtbins_Xi[fNptbins_Xi]; // Use Xi pT bin maximum
        particleName = "Xi";
        hist->GetXaxis()->SetRangeUser(fitMin, fitMax); // Set x-axis range for Xi
    }
    else
    {
        mass = fMass_Om;                  // Default mass if not specified
        fitMin = fPtbins_Om[0];           // Use Omega pT bin minimum
        fitMax = fPtbins_Om[fNptbins_Om]; // Use Omega pT bin maximum
        particleName = "Omega";
        hist->GetXaxis()->SetRangeUser(fitMin, fitMax); // Set x-axis range for Omega
    }

    TString fitFuncName = TString::Format("fitLevyTsallis_%s", particleName);

    // Define the function range from 0 up to a bit beyond histogram's pT max
    Double_t funcMinPt = 0.0;
    Double_t funcMaxPt = hist->GetXaxis()->GetXmax() * 1.1; // Extend a bit for drawing

    // Determine actual data range for fitting (from first bin with entries to last)
    Double_t dataFitMinPt = 0.0;
    Double_t dataFitMaxPt = hist->GetXaxis()->GetXmax();

    // Find first and last bins with data to set a sensible fit range
    int firstBin = 0;
    for (int i = 1; i <= hist->GetNbinsX(); ++i)
    {
        if (hist->GetBinContent(i) > 0)
        {
            firstBin = i;
            break;
        }
    }
    if (firstBin > 0)
        dataFitMinPt = hist->GetXaxis()->GetBinLowEdge(firstBin);
    else
        dataFitMinPt = hist->GetXaxis()->GetXmin(); // Fallback

    int lastBin = 0;
    for (int i = hist->GetNbinsX(); i >= 1; --i)
    {
        if (hist->GetBinContent(i) > 0)
        {
            lastBin = i;
            break;
        }
    }
    if (lastBin > 0)
        dataFitMaxPt = hist->GetXaxis()->GetBinUpEdge(lastBin);
    else
        dataFitMaxPt = hist->GetXaxis()->GetXmax(); // Fallback

    if (dataFitMinPt >= dataFitMaxPt && hist->GetEntries() > 0)
    { // If only one bin or range is problematic
        dataFitMinPt = hist->GetXaxis()->GetXmin();
        dataFitMaxPt = hist->GetXaxis()->GetXmax();
    }

    TF1 *fitFunc = new TF1(fitFuncName, levyTsallisFuncPt, funcMinPt, funcMaxPt, 4); // 4 parameters
    fitFunc->SetNpx(1000);                                                           // Increase number of points for better fit resolution
    fitFunc->SetLineStyle(kDashed);                                                  // Dashed line for fits
    if (multBin < 0)
        fitFunc->SetLineColor(kBlack); // Default color for global fits
    else
        fitFunc->SetLineColor(customCascPalette[multBin]); // Use custom color palette for multiplicity bins

    fitFunc->SetParName(0, "dNdy");
    fitFunc->SetParName(1, "T");
    fitFunc->SetParName(2, "n");
    fitFunc->SetParName(3, "mass");

    // Initial parameter estimates
    Double_t initial_dNdy = hist->Integral("width"); // Integral of dN/dpT is dN/dy
    if (initial_dNdy <= 1e-9)
        initial_dNdy = hist->GetMaximum() * hist->GetNbinsX() * hist->GetBinWidth(1) * 0.5; // Rough estimate
    if (initial_dNdy <= 1e-9)
        initial_dNdy = 1.0; // Absolute fallback

    fitFunc->SetParameter(0, initial_dNdy);
    fitFunc->SetParameter(1, 0.200); // T initial: 200 MeV
    fitFunc->SetParameter(2, 10.0);  // n initial
    fitFunc->FixParameter(3, mass);  // Fix particle mass

    // Parameter limits
    fitFunc->SetParLimits(0, initial_dNdy * 0.01, initial_dNdy * 100.0); // dN/dy
    fitFunc->SetParLimits(1, 0.050, 0.800);                              // T (50 MeV to 800 MeV)
    fitFunc->SetParLimits(2, 2.01, 60.0);                                // n > 2

    // printf("\nFitInit: %s spectrum (Title: %s)\n", particleName, hist->GetTitle());
    // printf("  Data fit range: %.2f - %.2f GeV/c\n", dataFitMinPt, dataFitMaxPt);
    // printf("  Initial params: dN/dy=%.3e, T=%.3f, n=%.2f, mass=%.4f (fixed)\n",
    //        fitFunc->GetParameter(0), fitFunc->GetParameter(1), fitFunc->GetParameter(2), fitFunc->GetParameter(3));

    // Perform the fit
    // "S" to return TFitResultPtr
    // "R" to use the range specified in TF1 for drawing (if different from fit range on data)
    // "I" to use integral of function in bin (recommended for binned data, esp. variable width)
    // "Q" for quiet mode (suppresses MINUIT printout during fit) - remove Q for debugging
    // "M" for Minuit "Migrad improved" strategy
    // "E" for Minos error analysis (more robust errors) - can be slow
    TFitResultPtr fitResult;
    for (int i = 0; i < 5; ++i)
    {
        fitResult = hist->Fit(fitFunc, "SRNQE+", "", dataFitMinPt, dataFitMaxPt);
    }

    if (fitResult.Get() && fitResult->IsValid() && fitResult->Status() == 0)
    { // Status = 0 means convergence OK
      // printf("FitOK: %s spectrum\n", particleName);
      // printf("  dN/dy = %.3e +/- %.2e\n", fitResult->Parameter(0), fitResult->ParError(0));
      // printf("  T     = %.4f +/- %.4f GeV\n", fitResult->Parameter(1), fitResult->ParError(1));
      // printf("  n     = %.3f +/- %.3f\n", fitResult->Parameter(2), fitResult->ParError(2));
      // printf("  Mass (fixed) = %.4f GeV/c^2\n", fitResult->Parameter(3));
      // printf("  Chi2/NDF = %.2f / %d = %.2f\n", fitResult->Chi2(), fitResult->Ndf(), (fitResult->Ndf() > 0 ? fitResult->Chi2() / fitResult->Ndf() : 0.0));
      // Fit function is automatically drawn on the histogram if it's displayed in a TCanvas
    }
    else
    {
        printf("FitFAIL: %s spectrum. Status: %d\n", particleName, (fitResult.Get() ? fitResult->Status() : -1));
    }
    // The TF1 object is associated with the histogram and will be deleted when the histogram is deleted,
    // or it can be retrieved using hist->GetFunction(fitFuncName).

    hist->GetXaxis()->SetRangeUser(xmin, xmax); // Set x-axis range for drawing purposes
    fitFunc->SetRange(xmin, xmax);              // Set range for drawing purposes
    hist->GetListOfFunctions()->Add(fitFunc);   // Add the fit function to the histogram's function list
};

/// Fit the stat effCorr spectra with Levy-Tsallis function
/**
 * @brief Fits a histogram with a Levy-Tsallis function and extrapolates the fit to 0.0 GeV/c for drawing purposes.
 *
 * This function fits the provided histogram with a Levy-Tsallis function,
 * which is commonly used in high-energy physics to describe particle spectra.
 * The fit parameters are initialized based on the histogram's content and
 * the specified particle type and multiplicity bin.
 *
 * @param hist Pointer to the histogram to be fitted
 * @param particle Optional parameter to specify which particle type is being fitted (-1 indicates default)
 * @param multBin Optional parameter to specify the multiplicity bin for the fit (-1 indicates default)
 */
// TF1 *FitLevy(TH1 *hist, Int_t particle = -1, Int_t multBin = -1)
// {
//     Double_t mass = 0.0;    // Default mass if not specified
//     Double_t fitMin = 0.0;  // Default fit minimum
//     Double_t fitMax = 10.0; // Default fit maximum

//     if (particle == kXi)
//     {
//         mass = fMass_Xi;
//         fitMin = fPtbins_Xi[0];           // Use Xi pT bin minimum
//         fitMax = fPtbins_Xi[fNptbins_Xi]; // Use Xi pT bin maximum
//     }
//     else
//     {
//         mass = fMass_Om;                  // Default mass if not specified
//         fitMin = fPtbins_Om[0];           // Use Omega pT bin minimum
//         fitMax = fPtbins_Om[fNptbins_Om]; // Use Omega pT bin maximum
//     }
//     TF1 *levyFunc = LevyTsallis("fLevyTsallis", mass);
//     levyFunc->SetNpx(1000); // Increase number of points for better fit resolution
//     if (multBin < 0)
//         levyFunc->SetLineColor(kBlack); // Default color for global fits
//     else
//         levyFunc->SetLineColor(customCascPalette[multBin]); // Use custom color palette for multiplicity bins
//     // Set the fit function to a Levy-Tsallis function
//     // TF1 *levyFit = new TF1("levyFit", "([0] / (TMath::Pi() * [1] * [2])) * (1 + ([3] / [2]) * (x - [4])^2)^(-[1])", 0.0, 10.0);

//     // // Initialize parameters based on the histogram
//     // levyFit->SetParameter(0, hist->GetMaximum() * 1.5); // Normalization
//     // levyFit->SetParameter(1, 0.5);                      // Exponent
//     // levyFit->SetParameter(2, 0.15);                     // Scale parameter
//     // levyFit->SetParameter(3, 0.01);                     // Offset
//     // levyFit->SetParameter(4, 0.0);                      // Mean

//     // Fit the histogram with the Levy-Tsallis function
//     hist->Fit(levyFunc, "QMSLEI MULTITHREAD", "", fitMin, fitMax);

//     // Extrapolate the fit to 0.0 GeV/c for drawing purposes
//     levyFunc->SetRange(0.0, 10.0);
//     hist->GetListOfFunctions()->Add(levyFunc);

//     return levyFunc;
// }

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
    TString inputEffCorrDefault = "/var/home/ishaan/Work/git/analysis/results/RandomVars/040625_effCorr_multIntCorrected_v2/040625_effCorr_multIntCorrected_def_6runs_v2.root",
    TString inputSysMultiTrial = "/var/home/ishaan/Work/git/thesis/final/images/310525_SysUncertainty_MultiTrial_212noRBErr/310525_sysUncertainty_noRB_multiTrial.root",
    TString inputSysSigExtractionPrefix = "/var/home/ishaan/Work/git/analysis/results/RandomVars/250225_Systematics_SigExt",
    TString inputSysSigExtractionFileName = "250225_SystematicUncertainty_SigExt.root",
    TString outputFolder = "/var/home/ishaan/Work/git/analysis/results/RandomVars/120625_SysUncertainty_Total_final",
    TString outputFileName = "/var/home/ishaan/Work/git/analysis/results/RandomVars/120625_SysUncertainty_Total_final/120625_sysUncertainty_Total_final.root",
    Bool_t saveStack = kTRUE,
    TString imageFormat = "pdf",
    Int_t verbosity = kInfo)
{
    ROOT::EnableImplicitMT();
    ROOT::Math::IntegratorOneDimOptions::SetDefaultIntegrator("Adaptive");
    ROOT::Math::IntegratorOneDimOptions::SetDefaultNPoints(6);
    // ROOT::Math::IntegratorOneDimOptions::SetDefaultAbsTolerance(1E-2);
    // ROOT::Math::IntegratorOneDimOptions::SetDefaultRelTolerance(1E-2);
    ROOT::Math::MinimizerOptions::SetDefaultMaxFunctionCalls(100000);
    ROOT::Math::MinimizerOptions::SetDefaultMinimizer("Minuit2");
    ROOT::Math::MinimizerOptions::SetDefaultStrategy(2);
    ROOT::Math::MinimizerOptions::SetDefaultPrintLevel(0); // Fit printing: -1 = no printing, 0 (minimal) to 3 (max)

    gStyle->SetPaintTextFormat("1.3f");
    gErrorIgnoreLevel = verbosity;

    SetCustomColorPalette();
    // gStyle->SetPalette(kRainbow);

    // extra options for decorating final plots (pdf) in thesis
    gStyle->SetLineScalePS(1.5);
    gStyle->SetStatFontSize(0.03);
    gStyle->SetPadTickX(1); // Ticks on both top and bottom for X axis
    gStyle->SetPadTickY(1); // Ticks on both left and right for Y axis
    // gStyle->SetHistTopMargin(10);
    // gStyle->SetHistBottomMargin(0.15);

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
    THStack *hs_uncertainty_xim = new THStack("hs_uncertainty_xim", "#Xi^{-}, 0-100% V0A");
    THStack *hs_uncertainty_xip = new THStack("hs_uncertainty_xip", "#Xi^{+}, 0-100% V0A");
    THStack *hs_uncertainty_omm = new THStack("hs_uncertainty_omm", "#Omega^{-}, 0-100% V0A");
    THStack *hs_uncertainty_omp = new THStack("hs_uncertainty_omp", "#Omega^{+}, 0-100% V0A");
    THStack *hs_uncertainty_xiC = new THStack("hs_uncertainty_xiC", "#Xi^{+} + #Xi^{-}, 0-100% V0A");
    THStack *hs_uncertainty_omC = new THStack("hs_uncertainty_omC", "#Omega^{+} + #Omega^{-}, 0-100% V0A");

    /// Create THStack objects for plotting efficiency corrected spectra with statistical and systematic uncertainties:
    THStack *hs_effCorrUncert_xim = new THStack("hs_effCorrUncert_xim", "#Xi^{-}, -0.5 < #it{y} < 0");
    THStack *hs_effCorrUncert_xip = new THStack("hs_effCorrUncert_xip", "#Xi^{+}, -0.5 < #it{y} < 0");
    THStack *hs_effCorrUncert_omm = new THStack("hs_effCorrUncert_omm", "#Omega^{-}, -0.5 < #it{y} < 0");
    THStack *hs_effCorrUncert_omp = new THStack("hs_effCorrUncert_omp", "#Omega^{+}, -0.5 < #it{y} < 0");
    THStack *hs_effCorrUncert_xiC = new THStack("hs_effCorrUncert_xiC", "#Xi^{+} + #Xi^{-}, -0.5 < #it{y} < 0");
    THStack *hs_effCorrUncert_omC = new THStack("hs_effCorrUncert_omC", "#Omega^{+} + #Omega^{-}, -0.5 < #it{y} < 0");

    /// Create THStack objects for plotting efficiency corrected spectra with statistical and systematic uncertainties and Levy-Tsallis fits:
    THStack *hs_effCorrUncertFit_xim = new THStack("hs_effCorrUncertFit_xim", "#Xi^{-}, -0.5 < #it{y} < 0");
    THStack *hs_effCorrUncertFit_xip = new THStack("hs_effCorrUncertFit_xip", "#Xi^{+}, -0.5 < #it{y} < 0");
    THStack *hs_effCorrUncertFit_omm = new THStack("hs_effCorrUncertFit_omm", "#Omega^{-}, -0.5 < #it{y} < 0");
    THStack *hs_effCorrUncertFit_omp = new THStack("hs_effCorrUncertFit_omp", "#Omega^{+}, -0.5 < #it{y} < 0");
    THStack *hs_effCorrUncertFit_xiC = new THStack("hs_effCorrUncertFit_xiC", "#Xi^{+} + #Xi^{-}, -0.5 < #it{y} < 0");
    THStack *hs_effCorrUncertFit_omC = new THStack("hs_effCorrUncertFit_omC", "#Omega^{+} + #Omega^{-}, -0.5 < #it{y} < 0");

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
    outputFile->mkdir("dirStack_EffCorrUncertaintyFit");

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
    // hs_uncertainty_xim->Add(stat_effCorrPt_xim_err, "TEXT00");
    // hs_uncertainty_xip->Add(stat_effCorrPt_xip_err, "TEXT00");
    // hs_uncertainty_omm->Add(stat_effCorrPt_omm_err, "TEXT00");
    // hs_uncertainty_omp->Add(stat_effCorrPt_omp_err, "TEXT00");
    // hs_uncertainty_xiC->Add(stat_effCorrPt_xiC_err, "TEXT00");
    // hs_uncertainty_omC->Add(stat_effCorrPt_omC_err, "TEXT00");
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
    // Style and add mult integrated effCorr histograms (sys and stat) to the stack:
    // Set style for systematic uncertainty boxes:
    StyleSystematics(sys_effCorrPt_xip, kXi, -1);
    StyleSystematics(sys_effCorrPt_xim, kXi, -1);
    StyleSystematics(sys_effCorrPt_omp, kXi, -1);
    StyleSystematics(sys_effCorrPt_omm, kXi, -1);
    StyleSystematics(sys_effCorrPt_xiC, kXi, -1);
    StyleSystematics(sys_effCorrPt_omC, kXi, -1);

    // // Set style for statistical error bars:
    StyleStatistics(stat_effCorrPt_xip, kXi, -1);
    StyleStatistics(stat_effCorrPt_xim, kXi, -1);
    StyleStatistics(stat_effCorrPt_omp, kXi, -1);
    StyleStatistics(stat_effCorrPt_omm, kXi, -1);
    StyleStatistics(stat_effCorrPt_xiC, kXi, -1);
    StyleStatistics(stat_effCorrPt_omC, kXi, -1);
    hs_effCorrUncert_xim->Add(sys_effCorrPt_xim, "E2");
    hs_effCorrUncert_xim->Add(stat_effCorrPt_xim, "P E1");
    hs_effCorrUncert_xip->Add(sys_effCorrPt_xip, "E2");
    hs_effCorrUncert_xip->Add(stat_effCorrPt_xip, "P E1");
    hs_effCorrUncert_omm->Add(sys_effCorrPt_omm, "E2");
    hs_effCorrUncert_omm->Add(stat_effCorrPt_omm, "P E1");
    hs_effCorrUncert_omp->Add(sys_effCorrPt_omp, "E2");
    hs_effCorrUncert_omp->Add(stat_effCorrPt_omp, "P E1");
    hs_effCorrUncert_xiC->Add(sys_effCorrPt_xiC, "E2");
    hs_effCorrUncert_xiC->Add(stat_effCorrPt_xiC, "P E1");
    hs_effCorrUncert_omC->Add(sys_effCorrPt_omC, "E2");
    hs_effCorrUncert_omC->Add(stat_effCorrPt_omC, "P E1");

    // Draw and save stack plots:
    TCanvas *cSpectra[6];
    for (Int_t iCanvas = 0; iCanvas < 6; iCanvas++)
    {
        cSpectra[iCanvas] = new TCanvas(TString::Format("cSpectra%d", iCanvas), TString::Format("cSpectra%d", iCanvas), 1200, 1000);
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

    // Add mult binned effCorr histograms (sys and stat) with FIT to the stack:
    for (Int_t multBinXi = 0; multBinXi < fNmultbins_Xi; multBinXi++)
    {
        // Style and add mult binned effCorr histograms (sys and stat) with fit:
        LevyTsallisFit(stat_effCorrPt_xip_mult[multBinXi], kXi, multBinXi);
        LevyTsallisFit(stat_effCorrPt_xim_mult[multBinXi], kXi, multBinXi);
        LevyTsallisFit(stat_effCorrPt_xiC_mult[multBinXi], kXi, multBinXi);
        hs_effCorrUncertFit_xim->Add(sys_effCorrPt_xim_mult[multBinXi], "E2");
        hs_effCorrUncertFit_xim->Add(stat_effCorrPt_xim_mult[multBinXi], "P E1");
        hs_effCorrUncertFit_xip->Add(sys_effCorrPt_xip_mult[multBinXi], "E2");
        hs_effCorrUncertFit_xip->Add(stat_effCorrPt_xip_mult[multBinXi], "P E1");
        hs_effCorrUncertFit_xiC->Add(sys_effCorrPt_xiC_mult[multBinXi], "E2");
        hs_effCorrUncertFit_xiC->Add(stat_effCorrPt_xiC_mult[multBinXi], "P E1");
    }
    for (Int_t multBinOm = 0; multBinOm < fNmultbins_Om; multBinOm++)
    {
        // Style and add mult binned effCorr histograms (sys and stat) with fit:
        LevyTsallisFit(stat_effCorrPt_omp_mult[multBinOm], kOm, multBinOm);
        LevyTsallisFit(stat_effCorrPt_omm_mult[multBinOm], kOm, multBinOm);
        LevyTsallisFit(stat_effCorrPt_omC_mult[multBinOm], kOm, multBinOm);
        hs_effCorrUncertFit_omm->Add(sys_effCorrPt_omm_mult[multBinOm], "E2");
        hs_effCorrUncertFit_omm->Add(stat_effCorrPt_omm_mult[multBinOm], "P E1");
        hs_effCorrUncertFit_omp->Add(sys_effCorrPt_omp_mult[multBinOm], "E2");
        hs_effCorrUncertFit_omp->Add(stat_effCorrPt_omp_mult[multBinOm], "P E1");
        hs_effCorrUncertFit_omC->Add(sys_effCorrPt_omC_mult[multBinOm], "E2");
        hs_effCorrUncertFit_omC->Add(stat_effCorrPt_omC_mult[multBinOm], "P E1");
    }

    // Style and add mult integrated effCorr histograms (sys and stat) with fit:
    LevyTsallisFit(stat_effCorrPt_xip, kXi, -1);
    LevyTsallisFit(stat_effCorrPt_xim, kXi, -1);
    LevyTsallisFit(stat_effCorrPt_omp, kOm, -1);
    LevyTsallisFit(stat_effCorrPt_omm, kOm, -1);
    LevyTsallisFit(stat_effCorrPt_xiC, kXi, -1);
    LevyTsallisFit(stat_effCorrPt_omC, kOm, -1);
    hs_effCorrUncertFit_xim->Add(sys_effCorrPt_xim, "E2");
    hs_effCorrUncertFit_xim->Add(stat_effCorrPt_xim, "P E1");
    hs_effCorrUncertFit_xip->Add(sys_effCorrPt_xip, "E2");
    hs_effCorrUncertFit_xip->Add(stat_effCorrPt_xip, "P E1");
    hs_effCorrUncertFit_omm->Add(sys_effCorrPt_omm, "E2");
    hs_effCorrUncertFit_omm->Add(stat_effCorrPt_omm, "P E1");
    hs_effCorrUncertFit_omp->Add(sys_effCorrPt_omp, "E2");
    hs_effCorrUncertFit_omp->Add(stat_effCorrPt_omp, "P E1");
    hs_effCorrUncertFit_xiC->Add(sys_effCorrPt_xiC, "E2");
    hs_effCorrUncertFit_xiC->Add(stat_effCorrPt_xiC, "P E1");
    hs_effCorrUncertFit_omC->Add(sys_effCorrPt_omC, "E2");
    hs_effCorrUncertFit_omC->Add(stat_effCorrPt_omC, "P E1");

    PaintStackOverlap(*cSpectra[0], *hs_effCorrUncertFit_xim, kTRUE, kTRUE);
    PaintStackOverlap(*cSpectra[1], *hs_effCorrUncertFit_xip, kTRUE, kTRUE);
    PaintStackOverlap(*cSpectra[2], *hs_effCorrUncertFit_omm, kTRUE, kTRUE);
    PaintStackOverlap(*cSpectra[3], *hs_effCorrUncertFit_omp, kTRUE, kTRUE);
    PaintStackOverlap(*cSpectra[4], *hs_effCorrUncertFit_xiC, kTRUE, kTRUE);
    PaintStackOverlap(*cSpectra[5], *hs_effCorrUncertFit_omC, kTRUE, kTRUE);

    if (saveStack)
    {
        cSpectra[0]->cd();
        SaveImage(outputFolder, "hs_effCorrUncertFit", "hs_effCorrUncertFit_xim", imageFormat.Data(), cSpectra[0]);

        cSpectra[1]->cd();
        SaveImage(outputFolder, "hs_effCorrUncertFit", "hs_effCorrUncertFit_xip", imageFormat.Data(), cSpectra[1]);

        cSpectra[2]->cd();
        SaveImage(outputFolder, "hs_effCorrUncertFit", "hs_effCorrUncertFit_omm", imageFormat.Data(), cSpectra[2]);

        cSpectra[3]->cd();
        SaveImage(outputFolder, "hs_effCorrUncertFit", "hs_effCorrUncertFit_omp", imageFormat.Data(), cSpectra[3]);

        cSpectra[4]->cd();
        SaveImage(outputFolder, "hs_effCorrUncertFit", "hs_effCorrUncertFit_xiC", imageFormat.Data(), cSpectra[4]);

        cSpectra[5]->cd();
        SaveImage(outputFolder, "hs_effCorrUncertFit", "hs_effCorrUncertFit_omC", imageFormat.Data(), cSpectra[5]);
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

    outputFile->cd("dirStack_EffCorrUncertaintyFit");
    hs_effCorrUncertFit_xim->Write();
    hs_effCorrUncertFit_xip->Write();
    hs_effCorrUncertFit_omm->Write();
    hs_effCorrUncertFit_omp->Write();
    hs_effCorrUncertFit_xiC->Write();
    hs_effCorrUncertFit_omC->Write();

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
    // auto chosenPalette = gStyle->GetColorPalette();

    // Divide the efficiency corrected pt spectrum with 2*pi and each pt bin center for correct normalization
    for (Int_t iBin = 1; iBin <= hist->GetNbinsX(); iBin++)
    {
        Double_t binCenter = hist->GetBinCenter(iBin);
        if (binCenter > 0)
        {
            hist->SetBinContent(iBin, hist->GetBinContent(iBin) / (2 * TMath::Pi() * binCenter));
            hist->SetBinError(iBin, hist->GetBinError(iBin) / (2 * TMath::Pi() * binCenter));
        }
    }

    hist->SetFillStyle(0);
    hist->SetLineWidth(1);
    hist->SetDrawOption("E2");
    hist->SetMarkerSize(0);
    hist->SetStats(kFALSE);
    if (multBin >= 0)
    {
        auto color = TColor::GetColorPalette(multBin); // Get color from the chosen palette based on the multiplicity bin
        hist->SetColors(color, color);                 // Set fill color and line color based on the chosen palette
        // hist->SetLineColor(chosenPalette[multBin]);
        if (particle == kXi)
        {
            hist->SetTitle(TString::Format(("%.0f-%.0f%% (#times2^{%d})"), fMultbins_Xi[multBin], fMultbins_Xi[multBin + 1], (fNmultbins_Xi - 1) - multBin));
            hist->Scale(pow(2, (fNmultbins_Xi - 1) - multBin));
        }
        else if (particle == kOm)
        {
            hist->SetTitle(TString::Format(("%.0f-%.0f%% (#times2^{%d})"), fMultbins_Om[multBin], fMultbins_Om[multBin + 1], (fNmultbins_Om - 1) - multBin));
            hist->Scale(pow(2, (fNmultbins_Om - 1) - multBin));
        }
    }
    else // integrated multiplicity case
    {
        hist->SetTitle("0-100% (#times2^{-5})");
        hist->SetColors(kBlack, kBlack); // Set color to black
        hist->Scale(pow(2, -5));
    }
}

void StyleStatistics(TH1 *hist, Int_t particle, Int_t multBin)
{

    Double_t xmin = 0.0;
    Double_t xmax = 10.0;
    // auto chosenPalette = gStyle->GetColorPalette();

    // Set style for statistical error bars:
    // hist->SetMarkerSize(2);

    // Divide the efficiency corrected pt spectrum with 2*pi and each pt bin center for correct normalization
    for (Int_t iBin = 1; iBin <= hist->GetNbinsX(); iBin++)
    {
        Double_t binCenter = hist->GetBinCenter(iBin);
        if (binCenter > 0)
        {
            hist->SetBinContent(iBin, hist->GetBinContent(iBin) / (2 * TMath::Pi() * binCenter));
            hist->SetBinError(iBin, hist->GetBinError(iBin) / (2 * TMath::Pi() * binCenter));
        }
    }

    // hist->SetLineWidth(3);
    hist->SetDrawOption("P E1 SAME");
    // hist->SetStats(kFALSE);
    if (multBin >= 0)
    {
        auto color = TColor::GetColorPalette(multBin); // Get color from the chosen palette based on the multiplicity bin
        hist->SetColors(color, color);                 // Set fill color and line color based on the chosen palette
        hist->SetMarkerStyle(markerStyles[multBin]);
        if (particle == kXi)
        {
            hist->SetTitle(TString::Format(("%.0f-%.0f%% (#times2^{%d})"), fMultbins_Xi[multBin], fMultbins_Xi[multBin + 1], (fNmultbins_Xi - 1) - multBin));
            hist->Scale(pow(2, (fNmultbins_Xi - 1) - multBin));
            // if (fit)
            //     LevyTsallisFit(hist, particle, multBin); // Fit the histogram with a Levy function
            hist->GetXaxis()->SetRangeUser(xmin, xmax); // Set x-axis range for better visibility
            // levyFit->SetLineColor(color);     // Set the fit line color to match the histogram
        }
        else if (particle == kOm)
        {
            hist->SetTitle(TString::Format(("%.0f-%.0f%% (#times2^{%d})"), fMultbins_Om[multBin], fMultbins_Om[multBin + 1], (fNmultbins_Om - 1) - multBin));
            hist->Scale(pow(2, (fNmultbins_Om - 1) - multBin));
            // if (fit)
            //     LevyTsallisFit(hist, particle, multBin); // Fit the histogram with a Levy function
            hist->GetXaxis()->SetRangeUser(xmin, xmax); // Set x-axis range for better visibility
            // levyFit->SetLineColor(color);     // Set the fit line color to match the histogram
        }
    }
    else // integrated multiplicity case
    {
        hist->SetColors(kBlack, kBlack); // Set color to black
        hist->SetMarkerStyle(kFullCircle);
        hist->SetMarkerSize(1.5);
        hist->SetTitle("0-100% (#times2^{-5})");
        hist->Scale(pow(2, -5));
        // if (fit)
        //     LevyTsallisFit(hist, particle, multBin); // Fit the histogram with a Levy function
        hist->GetXaxis()->SetRangeUser(xmin, xmax); // Set x-axis range for better visibility
    }
}

void PaintStackOverlap(TCanvas &c, THStack &hs, Bool_t setLogY, Bool_t fit, TString yAxisTitle, TString xAxisTitle)
{
    c.Clear();
    c.cd();

    c.SetLeftMargin(0.12);  // New mod for final thesis plots
    c.SetRightMargin(0.05); // New mod for final thesis plots
    c.SetTopMargin(0.05);   // New mod for final thesis plots
    if (setLogY)
        gPad->SetLogy();

    // New mod for final thesis plots
    TString hsTitle = hs.GetTitle();
    TLegend *legend;
    Double_t xmax = 10.0; // Default x-axis maximum
    if (hsTitle.Contains("Xi"))
    {
        legend = new TLegend(0.15, 0.12, 0.8, 0.3, "V0A Multiplicity Percentile");
        legend->SetNColumns(4);         // Set number of columns in the legend
        xmax = fPtbins_Xi[fNptbins_Xi]; // Set x-axis maximum for Xi
    }
    else
    {
        legend = new TLegend(0.17, 0.15, 0.68, 0.3, "V0A Multiplicity Percentile");
        legend->SetNColumns(3);         // Set number of columns in the legend
        xmax = fPtbins_Om[fNptbins_Om]; // Set x-axis maximum for Om
    }

    TH1F *frame = c.DrawFrame(0.0, 1e-8, xmax, hs.GetMaximum() * 5); // New mod for final thesis plots
    frame->SetXTitle(xAxisTitle.Data());
    frame->SetYTitle(yAxisTitle.Data());
    frame->SetTitleOffset(1.8, "Y"); // New mod for final thesis plots
    frame->SetTitleOffset(1.5, "X"); // New mod for final thesis plots
    frame->SetTitleSize(0.03, "Y");  // New mod for final thesis plots
    frame->SetTitleSize(0.03, "X");  // New mod for final thesis plots
    frame->SetLabelSize(0.028, "Y"); // New mod for final thesis plots
    // hs.SetMinimum(1e-7);                // Set minimum to 80% of the lowest bin content
    // hs.SetMaximum(hs.GetMaximum() * 5); // Set maximum to 150% of the highest bin content
    // hs.SetMinimum(hs.GetMinimum() * 0.000001); // Set minimum to 80% of the lowest bin content
    hs.Draw("nostack same");

    // Info("PaintStackOverlap", "xAxisTitle: %s, yAxisTitle: %s", xAxisTitle.Data(), yAxisTitle.Data());
    // hs.GetXaxis()->SetTitle(xAxisTitle.Data());
    // hs.GetYaxis()->SetTitle(yAxisTitle.Data());

    // TString hsTitle = hs.GetTitle();
    // // TLegend *legend = new TLegend(0.9, 0.6, 1., 1., "");
    // // New mod for final thesis plots
    // TLegend *legend;
    // if (hsTitle.Contains("Xi"))
    // {
    //     legend = new TLegend(0.15, 0.12, 0.8, 0.3, "V0A Multiplicity Percentile");
    //     legend->SetNColumns(4); // Set number of columns in the legend
    // }
    // else
    // {
    //     legend = new TLegend(0.17, 0.15, 0.68, 0.3, "V0A Multiplicity Percentile");
    //     legend->SetNColumns(3); // Set number of columns in the legend
    // }
    legend->SetBorderSize(0); // Remove border around the legend
    // legend->SetHeader("Sources", "C"); // Set header for the legend
    legend->SetFillColor(0); // Set legend background color to transparent
    legend->SetFillStyle(0); // Set fill style to transparent
    legend->SetTextSize(0.025);
    TLatex *latexCollTitle = new TLatex(0.63, 0.88, "p-Pb #sqrt{s_{NN}} = 8.16 TeV, This work");
    TLatex *latexPlotTitle = new TLatex(0.63, 0.84, hsTitle.Data());
    TLatex *latexInfo = new TLatex(0.63, 0.81, "Uncertainties: stat.(bars), sys.(boxes)");
    latexCollTitle->SetTextSize(0.025);
    latexPlotTitle->SetTextSize(0.025);
    latexInfo->SetTextSize(0.022);
    latexInfo->SetTextFont(42);
    latexCollTitle->SetNDC(kTRUE);
    latexPlotTitle->SetNDC(kTRUE);
    latexInfo->SetNDC(kTRUE);
    latexCollTitle->Draw();
    latexPlotTitle->Draw();
    latexInfo->Draw();
    hs.SetTitle("");

    TList *histList = hs.GetHists();
    for (int iHist = 0; iHist < histList->GetSize(); iHist++)
    {
        TH1 *hist = (TH1 *)histList->At(iHist);
        if (TString(hist->GetName()).Contains("stat"))
        {
            legend->AddEntry(hist, hist->GetTitle(), "lpef");
        }
    }
    TH1D *h = (TH1D *)histList->At(histList->GetSize() - 1); // Get the last histogram in the list

    // Add Levy fit of 0-100% stat effcorr to legend:
    if (fit)
        legend->AddEntry(h->GetListOfFunctions()->At(0), "L#acute{e}vy-Tsallis fit", "l");

    legend->Draw();
    c.Modified();
    c.ForceUpdate();
}
