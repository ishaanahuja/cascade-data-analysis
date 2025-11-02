#include <TROOT.h>

#include <TVirtualPad.h>
#include <TStyle.h>
#include <TLine.h>
#include <TLegend.h>
#include <TLatex.h>

#include "CascadeUtils.h"

/**
 * @brief Draws and optionally saves histograms with fit results and statistical information
 *
 * This function creates a canvas and draws the peak histogram along with background (if provided),
 * adds vertical lines indicating the signal region (±4σ around the peak position),
 * and displays a legend with fit statistics. The resulting plot can be optionally saved.
 *
 * @param peak Pointer to the peak histogram to be drawn
 * @param bg Pointer to the background histogram to be drawn (can be nullptr)
 * @param resultParams Pointer to histogram containing fit parameters:
 *                    - bin 1: Signal - Background
 *                    - bin 4: Peak position (mean)
 *                    - bin 5: Peak width (sigma)
 * @param saveImages Boolean flag to control whether to save the plot
 * @param outputFolder Directory path where images should be saved
 * @param imageFormat Format of the output image file
 *
 * @note The function automatically handles memory management by deleting created objects
 * @note Signal region is defined as ±4 sigma around the peak position
 */

void DrawAndSave(TH1 *peak, TH1 *bg, TH1 *resultParams, Bool_t saveImages, TString outputFolder, TString imageFormat);

/**
 * @brief Draws a horizontal line at y=1 on the given canvas and histogram stack
 *
 * This function sets the minimum and maximum y-axis values for the histogram stack if provided,
 * and then draws a horizontal line at y=1 on the specified canvas.
 *
 * @param c Reference to the canvas on which to draw the line
 * @param hs Reference to the histogram stack to set y-axis limits and draw the line
 * @param hs_yMin Minimum y-axis value for the histogram stack (default: -100, no change)
 * @param hs_yMax Maximum y-axis value for the histogram stack (default: -100, no change)
 */
void DrawRatioLine(TCanvas &c, THStack &hs, Double_t hs_yMin = -100, Double_t hs_yMax = -100);

/**
 * @brief Draws cascade analysis results and generates various plots
 *
 * This function processes and visualizes cascade particle analysis data, including:
 * - Mass distributions
 * - pT spectra
 * - Efficiencies
 * - Ratios
 * for Xi and Omega particles (both charges) and their combinations
 *
 * @param inputFilename Path to input ROOT file containing fitted analysis histograms
 * @param outputFilename Path for output ROOT file (default: derived from input name)
 * @param outputFolder Path for output images folder
 * @param ptRatioFilename Optional path to reference file (previous *_Draw.root file) for rawPt ratio comparison
 * @param fisMC Flag to indicate if input is from Monte Carlo simulation
 * @param saveImages Flag to save individual histogram images
 * @param saveStack Flag to save stacked histogram images
 * @param imageFormat Output image format (e.g. "png")
 * @param verbosity ROOT verbosity level
 *
 * @return Integer status code (0 for success)
 *
 * The function:
 * 1. Loads input histograms for Xi/Omega particles
 * 2. Processes pT spectra and efficiencies per multiplicity bin
 * 3. Creates stacked visualizations
 * 4. Optionally calculates ratios to reference spectra
 * 5. Saves results to ROOT file and image files
 *
 * Main outputs include:
 * - Mass distribution fits
 * - Raw pT spectra
 * - Efficiency curves (for MC)
 * - Ratio plots (if reference provided)
 * Both integrated and multiplicity-differential results are processed
 *
 * @note This is not mandatory to run anymore, as the functionality is replicated in EfficiencyEstimation.C, but can be used to generate raw pT spectra, fitting and efficiency plots
 */

void PaintStack_Draw(TCanvas &c, THStack &hs, Bool_t setLogY = kTRUE, TString yAxisTitle = "#frac{1}{#it{N}_{inel}} #frac{d#it{N}_{raw}}{d#it{p}_{T}} (GeV/#it{c})^{-1}", TString xAxisTitle = "#it{p}_{T} (GeV/#it{c})")
{
    c.Clear();
    c.cd();

    c.SetLeftMargin(0.12);  // New mod for final thesis plots
    c.SetRightMargin(0.03); // New mod for final thesis plots
    c.SetTopMargin(0.03);   // New mod for final thesis plots
    if (setLogY)
        gPad->SetLogy();

    // New mod for final thesis plots
    TString hsTitle = hs.GetTitle();
    TLegend *legend;
    Double_t xmin, ymin = 0.0;
    Double_t xmax, ymax = 10.0; // Default x-axis maximum
    if (hsTitle.Contains("Xi"))
    {
        legend = new TLegend(0.15, 0.12, 0.8, 0.3, "V0A Multiplicity Percentile");
        legend->SetNColumns(4); // Set number of columns in the legend
        xmin = fPtbins_Xi[0];
        ymin = 8e-7;
        xmax = fPtbins_Xi[fNptbins_Xi]; // Set x-axis maximum for Xi
        ymax =  hs.GetMaximum() * 7;
    }
    else
    {
        legend = new TLegend(0.17, 0.13, 0.68, 0.26, "V0A Multiplicity Percentile");
        legend->SetNColumns(3); // Set number of columns in the legend
        xmin = fPtbins_Om[0];
        ymin = 1e-6;
        xmax = fPtbins_Om[fNptbins_Om]; // Set x-axis maximum for Om
        ymax =  hs.GetMaximum() * 3;
    }

    TH1F *frame = c.DrawFrame(xmin, ymin, xmax, ymax); // New mod for final thesis plots
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
    hs.Draw("plc pmc nostack same");

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
    TLatex *latexCollTitle = new TLatex(0.68, 0.9, "p-Pb #sqrt{s_{NN}} = 8.16 TeV, This work");
    TLatex *latexPlotTitle = new TLatex(0.68, 0.86, hsTitle.Data());
    // TLatex *latexInfo = new TLatex(0.65, 0.81, "Uncertainties: stat.(bars), sys.(boxes)");
    latexCollTitle->SetTextSize(0.025);
    latexPlotTitle->SetTextSize(0.025);
    // latexInfo->SetTextSize(0.022);
    // latexInfo->SetTextFont(42);
    latexCollTitle->SetNDC(kTRUE);
    latexPlotTitle->SetNDC(kTRUE);
    // latexInfo->SetNDC(kTRUE);
    latexCollTitle->Draw();
    latexPlotTitle->Draw();
    // latexInfo->Draw();
    hs.SetTitle("");

    TList *histList = hs.GetHists();
    for (int iHist = 0; iHist < histList->GetSize(); iHist++)
    {
        TH1 *hist = (TH1 *)histList->At(iHist);
        // if (TString(hist->GetName()).Contains("stat"))
        // {
        legend->AddEntry(hist, hist->GetTitle(), "lpe");
        // }
    }
    // TH1D *h = (TH1D *)histList->At(histList->GetSize() - 1); // Get the last histogram in the list

    // Add Levy fit of 0-100% stat effcorr to legend:
    // if (fit)
    //     legend->AddEntry(h->GetListOfFunctions()->At(0), "L#acute{e}vy-Tsallis fit", "l");

    legend->Draw();
    c.Modified();
    c.ForceUpdate();
}

int DrawCascades(
    TString inputFilename = "/var/home/ishaan/Work/git/analysis/results/RandomVars/SysVars_SignalExtraction/190225_SysSigExt_Fit/190225_SysSigExt_DGP2_def_6Runs.root",
    TString outputFilename = "/var/home/ishaan/Work/git/thesis/final/images/120625_InvMassRawPt/120625_InvMassRawPt_def_6Runs_draw.root",
    TString outputFolder = "/var/home/ishaan/Work/git/thesis/final/images/120625_InvMassRawPt",
    TString ptRatioFilename = "",
    Bool_t fisMC = kFALSE,
    Bool_t saveImages = kFALSE,
    Bool_t saveStack = kTRUE,
    TString imageFormat = "pdf",
    Int_t verbosity = kInfo)
{
    // gPrintViaErrorHandler = kTRUE;
    gErrorIgnoreLevel = verbosity;
    gROOT->SetBatch(kTRUE);

    outputFolder = SetOutputFolder(outputFolder);
    gStyle->SetOptFit(1111);

    // SetCustomColorPalette();
    // extra options for decorating final plots (pdf) in thesis
    gStyle->SetLineScalePS(2);
    gStyle->SetStatFontSize(0.03);
    gStyle->SetPadTickX(1); // Ticks on both top and bottom for X axis
    gStyle->SetPadTickY(1); // Ticks on both left and right for Y axis

    TH1 *h_MassXim;
    TH1 *h_MassXip;
    TH1 *h_MassOmm;
    TH1 *h_MassOmp;
    TH1 *h_MassXiC;
    TH1 *h_MassOmC;
    TH1 *h_multBinEntries_Xi;
    TH1 *h_multBinEntries_Om;

    TH1 *resultParams_Xip_allInt, *resultParams_Xim_allInt, *resultParams_XiC_allInt;
    TH1 *resultParams_Omp_allInt, *resultParams_Omm_allInt, *resultParams_OmC_allInt;

    TH1 *h_MassXip_pt[fNptbins_Xi];
    TH1 *h_MassXim_pt[fNptbins_Xi];
    TH1 *h_MassOmp_pt[fNptbins_Om];
    TH1 *h_MassOmm_pt[fNptbins_Om];
    TH1 *h_MassXiC_pt[fNptbins_Xi];
    TH1 *h_MassOmC_pt[fNptbins_Om];

    TH1 *resultParXip_pt[fNptbins_Xi];
    TH1 *resultParXim_pt[fNptbins_Xi];
    TH1 *resultParOmp_pt[fNptbins_Om];
    TH1 *resultParOmm_pt[fNptbins_Om];
    TH1 *resultParXiC_pt[fNptbins_Xi];
    TH1 *resultParOmC_pt[fNptbins_Om];

    TH1 *h_MassXim_pt_mult[fNptbins_Xi][fNmultbins_Xi];
    TH1 *h_MassXip_pt_mult[fNptbins_Xi][fNmultbins_Xi];
    TH1 *h_MassOmm_pt_mult[fNptbins_Om][fNmultbins_Om];
    TH1 *h_MassOmp_pt_mult[fNptbins_Om][fNmultbins_Om];
    TH1 *h_MassXiC_pt_mult[fNptbins_Xi][fNmultbins_Xi];
    TH1 *h_MassOmC_pt_mult[fNptbins_Om][fNmultbins_Om];

    TH1 *resultParXip_pt_mult[fNptbins_Xi][fNmultbins_Xi];
    TH1 *resultParXim_pt_mult[fNptbins_Xi][fNmultbins_Xi];
    TH1 *resultParOmp_pt_mult[fNptbins_Om][fNmultbins_Om];
    TH1 *resultParOmm_pt_mult[fNptbins_Om][fNmultbins_Om];
    TH1 *resultParXiC_pt_mult[fNptbins_Xi][fNmultbins_Xi];
    TH1 *resultParOmC_pt_mult[fNptbins_Om][fNmultbins_Om];

    TH1D *rawPt_xim[fNmultbins_Xi];
    TH1D *rawPt_xip[fNmultbins_Xi];
    TH1D *rawPt_omm[fNmultbins_Om];
    TH1D *rawPt_omp[fNmultbins_Om];
    TH1D *rawPt_xiC[fNmultbins_Xi];
    TH1D *rawPt_omC[fNmultbins_Om];

    TH1D *eff_xim[fNmultbins_Xi];
    TH1D *eff_xip[fNmultbins_Xi];
    TH1D *eff_omm[fNmultbins_Om];
    TH1D *eff_omp[fNmultbins_Om];
    TH1D *eff_xiC[fNmultbins_Xi];
    TH1D *eff_omC[fNmultbins_Om];

    TH1D *eff_xiC_ratio[fNmultbins_Xi];
    TH1D *eff_omC_ratio[fNmultbins_Om];

    TH1D *eff_pt_xim = new TH1D("eff_pt_xim", "V0A 0-100%", fNptbins_Xi, fPtbins_Xi); /// mult integrated efficiency
    TH1D *eff_pt_xip = new TH1D("eff_pt_xip", "V0A 0-100%", fNptbins_Xi, fPtbins_Xi);
    TH1D *eff_pt_omm = new TH1D("eff_pt_omm", "V0A 0-100%", fNptbins_Om, fPtbins_Om);
    TH1D *eff_pt_omp = new TH1D("eff_pt_omp", "V0A 0-100%", fNptbins_Om, fPtbins_Om);
    TH1D *eff_pt_xiC = new TH1D("eff_pt_xiC", "V0A 0-100%", fNptbins_Xi, fPtbins_Xi);
    TH1D *eff_pt_omC = new TH1D("eff_pt_omC", "V0A 0-100%", fNptbins_Om, fPtbins_Om);

    /// background estimation hist through TSpectrum
    TH1 *h_bgXim_pt_mult[fNptbins_Xi][fNmultbins_Xi];
    TH1 *h_bgXip_pt_mult[fNptbins_Xi][fNmultbins_Xi];
    TH1 *h_bgXiC_pt_mult[fNptbins_Xi][fNmultbins_Xi];

    TH1 *h_bgOmm_pt_mult[fNptbins_Om][fNmultbins_Om];
    TH1 *h_bgOmp_pt_mult[fNptbins_Om][fNmultbins_Om];
    TH1 *h_bgOmC_pt_mult[fNptbins_Om][fNmultbins_Om];

    auto hs_xip = new THStack("hs_xip", "");
    auto hs_xim = new THStack("hs_xim", "");
    auto hs_omp = new THStack("hs_omp", "");
    auto hs_omm = new THStack("hs_omm", "");
    auto hs_xiC = new THStack("hs_xiC", "");
    auto hs_omC = new THStack("hs_omC", "");
    if (fisMC)
    {
        hs_xip->SetTitle("Rec. + Ass. #it{p}_{T} spectra #Xi^{+}");
        hs_xim->SetTitle("Rec. + Ass. #it{p}_{T} spectra #Xi^{-}");
        hs_omp->SetTitle("Rec. + Ass. #it{p}_{T} spectra #Omega^{+}");
        hs_omm->SetTitle("Rec. + Ass. #it{p}_{T} spectra #Omega^{-}");
        hs_xiC->SetTitle("Rec. + Ass. #it{p}_{T} spectra #Xi^{+} + #Xi^{-}");
        hs_omC->SetTitle("Rec. + Ass. #it{p}_{T} spectra #Omega^{+} + #Omega^{-}");
    }
    else
    {
        hs_xip->SetTitle("Raw #it{p}_{T} spectra #Xi^{+}");
        hs_xim->SetTitle("Raw #it{p}_{T} spectra #Xi^{-}");
        hs_omp->SetTitle("Raw #it{p}_{T} spectra #Omega^{+}");
        hs_omm->SetTitle("Raw #it{p}_{T} spectra #Omega^{-}");
        hs_xiC->SetTitle("Raw #it{p}_{T} spectra #Xi^{+} + #Xi^{-}");
        hs_omC->SetTitle("Raw #it{p}_{T} spectra #Omega^{+} + #Omega^{-}");
    }

    // THstack efficiency
    auto hs_xip_eff = new THStack("hs_xip_eff", "Efficiency #Xi^{+}");
    auto hs_xim_eff = new THStack("hs_xim_eff", "Efficiency #Xi^{-}");
    auto hs_omp_eff = new THStack("hs_omp_eff", "Efficiency #Omega^{+}");
    auto hs_omm_eff = new THStack("hs_omm_eff", "Efficiency #Omega^{-}");
    auto hs_xiC_eff = new THStack("hs_xiC_eff", "Efficiency #Xi^{+} + #Xi^{-}");
    auto hs_omC_eff = new THStack("hs_omC_eff", "Efficiency #Omega^{+} + #Omega^{-}");

    auto hs_xiC_eff_ratio = new THStack("hs_xiC_eff_ratio", "#Xi^{+} + #Xi^{-} efficiency ratio");
    auto hs_omC_eff_ratio = new THStack("hs_omC_eff_ratio", "#Omega^{+} + #Omega^{-} efficiency ratio");

    /// Getting histograms:

    {
        TFile *inputFile = OpenFile(inputFilename);

        // remove ownership of objects from file so we can delete the file ptr
        TH1::AddDirectory(kFALSE);

        h_MassXim = (TH1 *)inputFile->FindObjectAny("h_MassXim");
        h_MassXip = (TH1 *)inputFile->FindObjectAny("h_MassXip");
        h_MassXiC = (TH1 *)inputFile->FindObjectAny("h_MassXiC");
        h_MassOmm = (TH1 *)inputFile->FindObjectAny("h_MassOmm");
        h_MassOmp = (TH1 *)inputFile->FindObjectAny("h_MassOmp");
        h_MassOmC = (TH1 *)inputFile->FindObjectAny("h_MassOmC");
        h_multBinEntries_Xi = (TH1 *)inputFile->FindObjectAny("h_multBinEntries_Xi");
        h_multBinEntries_Om = (TH1 *)inputFile->FindObjectAny("h_multBinEntries_Om");

        resultParams_Xip_allInt = (TH1 *)inputFile->FindObjectAny("resultParams_Xip_allInt");
        resultParams_Xim_allInt = (TH1 *)inputFile->FindObjectAny("resultParams_Xim_allInt");
        resultParams_XiC_allInt = (TH1 *)inputFile->FindObjectAny("resultParams_XiC_allInt");
        resultParams_Omp_allInt = (TH1 *)inputFile->FindObjectAny("resultParams_Omp_allInt");
        resultParams_Omm_allInt = (TH1 *)inputFile->FindObjectAny("resultParams_Omm_allInt");
        resultParams_OmC_allInt = (TH1 *)inputFile->FindObjectAny("resultParams_OmC_allInt");

        for (Int_t ptBinXi = 0; ptBinXi < fNptbins_Xi; ptBinXi++)
        {
            for (Int_t multBinXi = 0; multBinXi < fNmultbins_Xi; multBinXi++)
            {
                h_MassXim_pt_mult[ptBinXi][multBinXi] = (TH1 *)inputFile->FindObjectAny(TString::Format(("h_MassXim_pt_mult[%d][%d]"), ptBinXi, multBinXi));
                h_MassXip_pt_mult[ptBinXi][multBinXi] = (TH1 *)inputFile->FindObjectAny(TString::Format(("h_MassXip_pt_mult[%d][%d]"), ptBinXi, multBinXi));
                resultParXip_pt_mult[ptBinXi][multBinXi] = (TH1 *)inputFile->FindObjectAny(TString::Format(("resultParXip_pt_mult[%d][%d]"), ptBinXi, multBinXi));
                resultParXim_pt_mult[ptBinXi][multBinXi] = (TH1 *)inputFile->FindObjectAny(TString::Format(("resultParXim_pt_mult[%d][%d]"), ptBinXi, multBinXi));

                h_MassXiC_pt_mult[ptBinXi][multBinXi] = (TH1 *)inputFile->FindObjectAny(TString::Format(("h_MassXiC_pt_mult[%d][%d]"), ptBinXi, multBinXi));
                resultParXiC_pt_mult[ptBinXi][multBinXi] = (TH1 *)inputFile->FindObjectAny(TString::Format(("resultParXiC_pt_mult[%d][%d]"), ptBinXi, multBinXi));

                h_bgXim_pt_mult[ptBinXi][multBinXi] = (TH1 *)inputFile->FindObjectAny(TString::Format(("h_bgXim_pt_mult[%d][%d]"), ptBinXi, multBinXi));
                h_bgXip_pt_mult[ptBinXi][multBinXi] = (TH1 *)inputFile->FindObjectAny(TString::Format(("h_bgXip_pt_mult[%d][%d]"), ptBinXi, multBinXi));
                h_bgXiC_pt_mult[ptBinXi][multBinXi] = (TH1 *)inputFile->FindObjectAny(TString::Format(("h_bgXiC_pt_mult[%d][%d]"), ptBinXi, multBinXi));
            }
            h_MassXip_pt[ptBinXi] = (TH1 *)inputFile->FindObjectAny(TString::Format(("h_MassXip_pt[%d]"), ptBinXi));
            h_MassXim_pt[ptBinXi] = (TH1 *)inputFile->FindObjectAny(TString::Format(("h_MassXim_pt[%d]"), ptBinXi));
            h_MassXiC_pt[ptBinXi] = (TH1 *)inputFile->FindObjectAny(TString::Format(("h_MassXiC_pt[%d]"), ptBinXi));
            resultParXip_pt[ptBinXi] = (TH1 *)inputFile->FindObjectAny(TString::Format(("resultParXip_pt[%d]"), ptBinXi));
            resultParXim_pt[ptBinXi] = (TH1 *)inputFile->FindObjectAny(TString::Format(("resultParXim_pt[%d]"), ptBinXi));
            resultParXiC_pt[ptBinXi] = (TH1 *)inputFile->FindObjectAny(TString::Format(("resultParXiC_pt[%d]"), ptBinXi));
        }
        for (Int_t ptBinOm = 0; ptBinOm < fNptbins_Om; ptBinOm++)
        {
            for (Int_t multBinOm = 0; multBinOm < fNmultbins_Om; multBinOm++)
            {
                h_MassOmm_pt_mult[ptBinOm][multBinOm] = (TH1 *)inputFile->FindObjectAny(TString::Format(("h_MassOmm_pt_mult[%d][%d]"), ptBinOm, multBinOm));
                h_MassOmp_pt_mult[ptBinOm][multBinOm] = (TH1 *)inputFile->FindObjectAny(TString::Format(("h_MassOmp_pt_mult[%d][%d]"), ptBinOm, multBinOm));
                resultParOmp_pt_mult[ptBinOm][multBinOm] = (TH1 *)inputFile->FindObjectAny(TString::Format(("resultParOmp_pt_mult[%d][%d]"), ptBinOm, multBinOm));
                resultParOmm_pt_mult[ptBinOm][multBinOm] = (TH1 *)inputFile->FindObjectAny(TString::Format(("resultParOmm_pt_mult[%d][%d]"), ptBinOm, multBinOm));

                h_MassOmC_pt_mult[ptBinOm][multBinOm] = (TH1 *)inputFile->FindObjectAny(TString::Format(("h_MassOmC_pt_mult[%d][%d]"), ptBinOm, multBinOm));
                resultParOmC_pt_mult[ptBinOm][multBinOm] = (TH1 *)inputFile->FindObjectAny(TString::Format(("resultParOmC_pt_mult[%d][%d]"), ptBinOm, multBinOm));

                h_bgOmm_pt_mult[ptBinOm][multBinOm] = (TH1 *)inputFile->FindObjectAny(TString::Format(("h_bgOmm_pt_mult[%d][%d]"), ptBinOm, multBinOm));
                h_bgOmp_pt_mult[ptBinOm][multBinOm] = (TH1 *)inputFile->FindObjectAny(TString::Format(("h_bgOmp_pt_mult[%d][%d]"), ptBinOm, multBinOm));
                h_bgOmC_pt_mult[ptBinOm][multBinOm] = (TH1 *)inputFile->FindObjectAny(TString::Format(("h_bgOmC_pt_mult[%d][%d]"), ptBinOm, multBinOm));
            }
            h_MassOmp_pt[ptBinOm] = (TH1 *)inputFile->FindObjectAny(TString::Format(("h_MassOmp_pt[%d]"), ptBinOm));
            h_MassOmm_pt[ptBinOm] = (TH1 *)inputFile->FindObjectAny(TString::Format(("h_MassOmm_pt[%d]"), ptBinOm));
            h_MassOmC_pt[ptBinOm] = (TH1 *)inputFile->FindObjectAny(TString::Format(("h_MassOmC_pt[%d]"), ptBinOm));
            resultParOmp_pt[ptBinOm] = (TH1 *)inputFile->FindObjectAny(TString::Format(("resultParOmp_pt[%d]"), ptBinOm));
            resultParOmm_pt[ptBinOm] = (TH1 *)inputFile->FindObjectAny(TString::Format(("resultParOmm_pt[%d]"), ptBinOm));
            resultParOmC_pt[ptBinOm] = (TH1 *)inputFile->FindObjectAny(TString::Format(("resultParOmC_pt[%d]"), ptBinOm));
        }

        Info("DrawCascades: histInput", "Input ended.");
        delete inputFile;
    } /// Input ended!

    /// Saving output :
    if (outputFilename.IsNull())
    {
        outputFilename = inputFilename;
        if (outputFilename.Contains("Fitting"))
            outputFilename.ReplaceAll("Fitting", "Draw");
        else
            outputFilename.ReplaceAll(".root", "_Draw.root");
    }

    TFile *outputFile = OpenFile(outputFilename, "RECREATE");

    outputFile->mkdir("dirRawPt_xim");
    outputFile->mkdir("dirRawPt_xip");
    outputFile->mkdir("dirRawPt_omm");
    outputFile->mkdir("dirRawPt_omp");

    outputFile->mkdir("dirRawPt_xiC");
    outputFile->mkdir("dirRawPt_omC");
    outputFile->mkdir("dirEffPt");
    /// Output set.

    Info("DrawCascades: outputFile", "Output file created.");

    /// Begin
    if (fisMC)
    {
        Info("DrawCascades: fisMC", "Plotting efficiency ...");
        outputFile->cd("dirEffPt");

        for (Int_t ptBinXi = 0; ptBinXi < fNptbins_Xi; ptBinXi++)
        {
            eff_pt_xim->SetBinContent(ptBinXi + 1, resultParXim_pt[ptBinXi]->GetBinContent(7));
            eff_pt_xip->SetBinContent(ptBinXi + 1, resultParXip_pt[ptBinXi]->GetBinContent(7));
            eff_pt_xiC->SetBinContent(ptBinXi + 1, resultParXiC_pt[ptBinXi]->GetBinContent(7));

            eff_pt_xim->SetBinError(ptBinXi + 1, resultParXim_pt[ptBinXi]->GetBinError(7));
            eff_pt_xip->SetBinError(ptBinXi + 1, resultParXip_pt[ptBinXi]->GetBinError(7));
            eff_pt_xiC->SetBinError(ptBinXi + 1, resultParXiC_pt[ptBinXi]->GetBinError(7));
        }

        eff_pt_xim->SetMarkerStyle(markerStyles[10]);
        eff_pt_xim->Write();
        eff_pt_xim->SetName(TString::Format(("V0A %.0f-%.0f%%"), fMultbins_Xi[0], fMultbins_Xi[fNmultbins_Xi]));
        hs_xim_eff->Add(eff_pt_xim);

        eff_pt_xip->SetMarkerStyle(markerStyles[10]);
        eff_pt_xip->Write();
        eff_pt_xip->SetName(TString::Format(("V0A %.0f-%.0f%%"), fMultbins_Xi[0], fMultbins_Xi[fNmultbins_Xi]));
        hs_xip_eff->Add(eff_pt_xip);

        eff_pt_xiC->SetMarkerStyle(markerStyles[10]);
        eff_pt_xiC->Write();
        eff_pt_xiC->SetName(TString::Format(("V0A %.0f-%.0f%%"), fMultbins_Xi[0], fMultbins_Xi[fNmultbins_Xi]));
        hs_xiC_eff->Add(eff_pt_xiC);

        for (Int_t ptBinOm = 0; ptBinOm < fNptbins_Om; ptBinOm++)
        {
            eff_pt_omm->SetBinContent(ptBinOm + 1, resultParOmm_pt[ptBinOm]->GetBinContent(7));
            eff_pt_omp->SetBinContent(ptBinOm + 1, resultParOmp_pt[ptBinOm]->GetBinContent(7));
            eff_pt_omC->SetBinContent(ptBinOm + 1, resultParOmC_pt[ptBinOm]->GetBinContent(7));

            eff_pt_omm->SetBinError(ptBinOm + 1, resultParOmm_pt[ptBinOm]->GetBinError(7));
            eff_pt_omp->SetBinError(ptBinOm + 1, resultParOmp_pt[ptBinOm]->GetBinError(7));
            eff_pt_omC->SetBinError(ptBinOm + 1, resultParOmC_pt[ptBinOm]->GetBinError(7));
        }

        eff_pt_omm->SetMarkerStyle(markerStyles[10]);
        eff_pt_omm->Write();
        eff_pt_omm->SetName(TString::Format(("V0A %.0f-%.0f%%"), fMultbins_Om[0], fMultbins_Om[fNmultbins_Om]));
        hs_omm_eff->Add(eff_pt_omm);

        eff_pt_omp->SetMarkerStyle(markerStyles[10]);
        eff_pt_omp->Write();
        eff_pt_omp->SetName(TString::Format(("V0A %.0f-%.0f%%"), fMultbins_Om[0], fMultbins_Om[fNmultbins_Om]));
        hs_omp_eff->Add(eff_pt_omp);

        eff_pt_omC->SetMarkerStyle(markerStyles[10]);
        eff_pt_omC->Write();
        eff_pt_omC->SetName(TString::Format(("V0A %.0f-%.0f%%"), fMultbins_Om[0], fMultbins_Om[fNmultbins_Om]));
        hs_omC_eff->Add(eff_pt_omC);
    }

    /// XI:
    Info("DrawCascades", "Plotting Xi Data ...");

    h_MassXim->GetXaxis()->SetTitle("M_{#Lambda^{0}#pi^{#pm}} (GeV/ #it{c}^{2})");
    h_MassXip->GetXaxis()->SetTitle("M_{#Lambda^{0}#pi^{#pm}} (GeV/ #it{c}^{2})");
    h_MassXiC->GetXaxis()->SetTitle("M_{#Lambda^{0}#pi^{#pm}} (GeV/ #it{c}^{2})");

    DrawAndSave(h_MassXim, nullptr, resultParams_Xim_allInt, saveImages, outputFolder, imageFormat);
    DrawAndSave(h_MassXip, nullptr, resultParams_Xip_allInt, saveImages, outputFolder, imageFormat);
    DrawAndSave(h_MassXiC, nullptr, resultParams_XiC_allInt, saveImages, outputFolder, imageFormat);

    /// Begin mult integrated part
    for (Int_t ptBinXi = 0; ptBinXi < fNptbins_Xi; ptBinXi++)
    {
        h_MassXim_pt[ptBinXi]->GetXaxis()->SetTitle("M_{#Lambda^{0}#pi^{#pm}} (GeV/ #it{c}^{2})");
        h_MassXip_pt[ptBinXi]->GetXaxis()->SetTitle("M_{#Lambda^{0}#pi^{#pm}} (GeV/ #it{c}^{2})");
        h_MassXiC_pt[ptBinXi]->GetXaxis()->SetTitle("M_{#Lambda^{0}#pi^{#pm}} (GeV/ #it{c}^{2})");
        DrawAndSave(h_MassXim_pt[ptBinXi], nullptr, resultParXim_pt[ptBinXi], saveImages, outputFolder, imageFormat);
        DrawAndSave(h_MassXip_pt[ptBinXi], nullptr, resultParXip_pt[ptBinXi], saveImages, outputFolder, imageFormat);
        DrawAndSave(h_MassXiC_pt[ptBinXi], nullptr, resultParXiC_pt[ptBinXi], saveImages, outputFolder, imageFormat);
    }

    /// Begin differential part
    for (Int_t multBinXi = 0; multBinXi < fNmultbins_Xi; multBinXi++)
    {
        rawPt_xim[multBinXi] = new TH1D(TString::Format(("rawPt_xim[%d]"), multBinXi), TString::Format(("%.0f-%.0f%% #times2^{%d}"), fMultbins_Xi[multBinXi], fMultbins_Xi[multBinXi + 1], (fNmultbins_Xi - 1) - multBinXi), fNptbins_Xi, fPtbins_Xi);
        rawPt_xip[multBinXi] = new TH1D(TString::Format(("rawPt_xip[%d]"), multBinXi), TString::Format(("%.0f-%.0f%% #times2^{%d}"), fMultbins_Xi[multBinXi], fMultbins_Xi[multBinXi + 1], (fNmultbins_Xi - 1) - multBinXi), fNptbins_Xi, fPtbins_Xi);
        rawPt_xiC[multBinXi] = new TH1D(TString::Format(("rawPt_xiC[%d]"), multBinXi), TString::Format(("%.0f-%.0f%% #times2^{%d}"), fMultbins_Xi[multBinXi], fMultbins_Xi[multBinXi + 1], (fNmultbins_Xi - 1) - multBinXi), fNptbins_Xi, fPtbins_Xi);

        if (fisMC)
        {
            eff_xim[multBinXi] = new TH1D(TString::Format(("eff_xim[%d]"), multBinXi), "", fNptbins_Xi, fPtbins_Xi);
            eff_xip[multBinXi] = new TH1D(TString::Format(("eff_xip[%d]"), multBinXi), "", fNptbins_Xi, fPtbins_Xi);

            eff_xiC[multBinXi] = new TH1D(TString::Format(("eff_xiC[%d]"), multBinXi), "", fNptbins_Xi, fPtbins_Xi);
        }
        for (Int_t ptBinXi = 0; ptBinXi < fNptbins_Xi; ptBinXi++)
        {
            rawPt_xim[multBinXi]->SetBinContent(ptBinXi + 1, resultParXim_pt_mult[ptBinXi][multBinXi]->GetBinContent(1) / ((rawPt_xim[multBinXi]->GetBinWidth(ptBinXi + 1)) * (h_multBinEntries_Xi->GetBinContent(multBinXi + 1)))); // Bin 1 in resultparams is raw pt's bin counting
            rawPt_xim[multBinXi]->SetBinError(ptBinXi + 1, resultParXim_pt_mult[ptBinXi][multBinXi]->GetBinError(1) / ((rawPt_xim[multBinXi]->GetBinWidth(ptBinXi + 1)) * (h_multBinEntries_Xi->GetBinContent(multBinXi + 1))));

            rawPt_xip[multBinXi]->SetBinContent(ptBinXi + 1, resultParXip_pt_mult[ptBinXi][multBinXi]->GetBinContent(1) / ((rawPt_xip[multBinXi]->GetBinWidth(ptBinXi + 1)) * (h_multBinEntries_Xi->GetBinContent(multBinXi + 1)))); // Bin 1 in resultparams is raw pt's bin counting
            rawPt_xip[multBinXi]->SetBinError(ptBinXi + 1, resultParXip_pt_mult[ptBinXi][multBinXi]->GetBinError(1) / ((rawPt_xip[multBinXi]->GetBinWidth(ptBinXi + 1)) * (h_multBinEntries_Xi->GetBinContent(multBinXi + 1))));

            rawPt_xiC[multBinXi]->SetBinContent(ptBinXi + 1, resultParXiC_pt_mult[ptBinXi][multBinXi]->GetBinContent(1) / ((rawPt_xiC[multBinXi]->GetBinWidth(ptBinXi + 1)) * (h_multBinEntries_Xi->GetBinContent(multBinXi + 1)))); // Bin 1 in resultparams is raw pt's bin counting
            rawPt_xiC[multBinXi]->SetBinError(ptBinXi + 1, resultParXiC_pt_mult[ptBinXi][multBinXi]->GetBinError(1) / ((rawPt_xiC[multBinXi]->GetBinWidth(ptBinXi + 1)) * (h_multBinEntries_Xi->GetBinContent(multBinXi + 1))));

            if (fisMC)
            {
                eff_xim[multBinXi]->SetBinContent(ptBinXi + 1, resultParXim_pt_mult[ptBinXi][multBinXi]->GetBinContent(7)); // Bin 7 in resultparams is MC efficiency
                eff_xim[multBinXi]->SetBinError(ptBinXi + 1, resultParXim_pt_mult[ptBinXi][multBinXi]->GetBinError(7));

                eff_xip[multBinXi]->SetBinContent(ptBinXi + 1, resultParXip_pt_mult[ptBinXi][multBinXi]->GetBinContent(7)); // Bin 7 in resultparams is MC efficiency
                eff_xip[multBinXi]->SetBinError(ptBinXi + 1, resultParXip_pt_mult[ptBinXi][multBinXi]->GetBinError(7));

                eff_xiC[multBinXi]->SetBinContent(ptBinXi + 1, resultParXiC_pt_mult[ptBinXi][multBinXi]->GetBinContent(7)); // Bin 7 in resultparams is MC efficiency
                eff_xiC[multBinXi]->SetBinError(ptBinXi + 1, resultParXiC_pt_mult[ptBinXi][multBinXi]->GetBinError(7));
            }

            /// Scale for bin width: N->dN/dpt
            // rawPt_xim[multBinXi]->Scale(1, "width");
            // rawPt_xip[multBinXi]->Scale(1, "width");

            h_MassXim_pt_mult[ptBinXi][multBinXi]->GetXaxis()->SetTitle("M_{#Lambda^{0}#pi^{#pm}} (GeV/ #it{c}^{2})");
            h_MassXip_pt_mult[ptBinXi][multBinXi]->GetXaxis()->SetTitle("M_{#Lambda^{0}#pi^{#pm}} (GeV/ #it{c}^{2})");
            h_MassXiC_pt_mult[ptBinXi][multBinXi]->GetXaxis()->SetTitle("M_{#Lambda^{0}#pi^{#pm}} (GeV/ #it{c}^{2})");

            DrawAndSave(h_MassXim_pt_mult[ptBinXi][multBinXi], h_bgXim_pt_mult[ptBinXi][multBinXi], resultParXim_pt_mult[ptBinXi][multBinXi], saveImages, outputFolder, imageFormat);
            DrawAndSave(h_MassXip_pt_mult[ptBinXi][multBinXi], h_bgXip_pt_mult[ptBinXi][multBinXi], resultParXip_pt_mult[ptBinXi][multBinXi], saveImages, outputFolder, imageFormat);
            DrawAndSave(h_MassXiC_pt_mult[ptBinXi][multBinXi], h_bgXiC_pt_mult[ptBinXi][multBinXi], resultParXiC_pt_mult[ptBinXi][multBinXi], saveImages, outputFolder, imageFormat);
        }
        outputFile->cd("dirRawPt_xip");
        rawPt_xip[multBinXi]->SetMarkerStyle(markerStyles[multBinXi]);
        rawPt_xip[multBinXi]->Write();
        rawPt_xip[multBinXi]->SetName(TString::Format(("V0A %.0f-%.0f%% #times2^{%d}"), fMultbins_Xi[multBinXi], fMultbins_Xi[multBinXi + 1], (fNmultbins_Xi - 1) - multBinXi));
        rawPt_xip[multBinXi]->Scale(pow(2, (fNmultbins_Xi - 1) - multBinXi));
        hs_xip->Add(rawPt_xip[multBinXi]);

        if (fisMC)
        {
            eff_xip[multBinXi]->SetMarkerStyle(markerStyles[multBinXi]);
            eff_xip[multBinXi]->Write();
            eff_xip[multBinXi]->SetName(TString::Format(("V0A %.0f-%.0f%%"), fMultbins_Xi[multBinXi], fMultbins_Xi[multBinXi + 1]));
            // eff_xip[multBinXi]->Scale(pow(2, (fNmultbins_Xi - 1) - multBinXi));
            hs_xip_eff->Add(eff_xip[multBinXi]);
        }

        outputFile->cd("dirRawPt_xim");
        rawPt_xim[multBinXi]->SetMarkerStyle(markerStyles[multBinXi]);
        rawPt_xim[multBinXi]->Write();
        rawPt_xim[multBinXi]->SetName(TString::Format(("V0A %.0f-%.0f%% #times2^{%d}"), fMultbins_Xi[multBinXi], fMultbins_Xi[multBinXi + 1], (fNmultbins_Xi - 1) - multBinXi));
        rawPt_xim[multBinXi]->Scale(pow(2, (fNmultbins_Xi - 1) - multBinXi));
        hs_xim->Add(rawPt_xim[multBinXi]);

        if (fisMC)
        {
            eff_xim[multBinXi]->SetMarkerStyle(markerStyles[multBinXi]);
            eff_xim[multBinXi]->Write();
            eff_xim[multBinXi]->SetName(TString::Format(("V0A %.0f-%.0f%%"), fMultbins_Xi[multBinXi], fMultbins_Xi[multBinXi + 1]));
            // eff_xim[multBinXi]->Scale(pow(2, (fNmultbins_Xi - 1) - multBinXi));
            hs_xim_eff->Add(eff_xim[multBinXi]);
        }
        outputFile->cd("dirRawPt_xiC");
        rawPt_xiC[multBinXi]->SetMarkerStyle(markerStyles[multBinXi]);
        rawPt_xiC[multBinXi]->Write();
        rawPt_xiC[multBinXi]->SetName(TString::Format(("V0A %.0f-%.0f%% #times2^{%d}"), fMultbins_Xi[multBinXi], fMultbins_Xi[multBinXi + 1], (fNmultbins_Xi - 1) - multBinXi));
        rawPt_xiC[multBinXi]->Scale(pow(2, (fNmultbins_Xi - 1) - multBinXi));
        hs_xiC->Add(rawPt_xiC[multBinXi]);

        if (fisMC)
        {
            eff_xiC[multBinXi]->SetMarkerStyle(markerStyles[multBinXi]);
            eff_xiC[multBinXi]->Write();
            eff_xiC[multBinXi]->SetName(TString::Format(("V0A %.0f-%.0f%%"), fMultbins_Xi[multBinXi], fMultbins_Xi[multBinXi + 1]));
            // eff_xiC[multBinXi]->Scale(pow(2, (fNmultbins_Xi - 1) - multBinXi));
            hs_xiC_eff->Add(eff_xiC[multBinXi]);

            eff_xiC_ratio[multBinXi] = new TH1D(TString::Format(("eff_xiC_ratio[%d]"), multBinXi), "", fNptbins_Xi, fPtbins_Xi);
            eff_xiC_ratio[multBinXi]->Divide(eff_xiC[multBinXi], eff_pt_xiC);
            eff_xiC_ratio[multBinXi]->SetMarkerStyle(markerStyles[multBinXi]);
            eff_xiC_ratio[multBinXi]->Write();
            eff_xiC_ratio[multBinXi]->SetName(TString::Format(("V0A %.0f-%.0f%%"), fMultbins_Xi[multBinXi], fMultbins_Xi[multBinXi + 1]));
            hs_xiC_eff_ratio->Add(eff_xiC_ratio[multBinXi]);
        }
    }

    /// OMEGA:
    Info("DrawCascades", "Plotting Omega Data ...");

    h_MassOmm->GetXaxis()->SetTitle("M_{#Lambda^{0}K^{#pm}} (GeV/ #it{c}^{2})");
    h_MassOmp->GetXaxis()->SetTitle("M_{#Lambda^{0}K^{#pm}} (GeV/ #it{c}^{2})");
    h_MassOmC->GetXaxis()->SetTitle("M_{#Lambda^{0}K^{#pm}} (GeV/ #it{c}^{2})");

    DrawAndSave(h_MassOmm, nullptr, resultParams_Omm_allInt, saveImages, outputFolder, imageFormat);
    DrawAndSave(h_MassOmp, nullptr, resultParams_Omp_allInt, saveImages, outputFolder, imageFormat);
    DrawAndSave(h_MassOmC, nullptr, resultParams_OmC_allInt, saveImages, outputFolder, imageFormat);

    /// Begin mult integrated part
    for (Int_t ptBinOm = 0; ptBinOm < fNptbins_Om; ptBinOm++)
    {
        h_MassOmm_pt[ptBinOm]->GetXaxis()->SetTitle("M_{#Lambda^{0}K^{#pm}} (GeV/ #it{c}^{2})");
        h_MassOmp_pt[ptBinOm]->GetXaxis()->SetTitle("M_{#Lambda^{0}K^{#pm}} (GeV/ #it{c}^{2})");
        h_MassOmC_pt[ptBinOm]->GetXaxis()->SetTitle("M_{#Lambda^{0}K^{#pm}} (GeV/ #it{c}^{2})");

        DrawAndSave(h_MassOmm_pt[ptBinOm], nullptr, resultParOmm_pt[ptBinOm], saveImages, outputFolder, imageFormat);
        DrawAndSave(h_MassOmp_pt[ptBinOm], nullptr, resultParOmp_pt[ptBinOm], saveImages, outputFolder, imageFormat);
        DrawAndSave(h_MassOmC_pt[ptBinOm], nullptr, resultParOmC_pt[ptBinOm], saveImages, outputFolder, imageFormat);
    }

    /// Begin differential part
    for (Int_t multBinOm = 0; multBinOm < fNmultbins_Om; multBinOm++)
    {
        rawPt_omm[multBinOm] = new TH1D(TString::Format(("rawPt_omm[%d]"), multBinOm), TString::Format(("%.0f-%.0f%% #times2^{%d}"), fMultbins_Om[multBinOm], fMultbins_Om[multBinOm + 1], (fNmultbins_Om - 1) - multBinOm), fNptbins_Om, fPtbins_Om);
        rawPt_omp[multBinOm] = new TH1D(TString::Format(("rawPt_omp[%d]"), multBinOm), TString::Format(("%.0f-%.0f%% #times2^{%d}"), fMultbins_Om[multBinOm], fMultbins_Om[multBinOm + 1], (fNmultbins_Om - 1) - multBinOm), fNptbins_Om, fPtbins_Om);
        rawPt_omC[multBinOm] = new TH1D(TString::Format(("rawPt_omC[%d]"), multBinOm), TString::Format(("%.0f-%.0f%% #times2^{%d}"), fMultbins_Om[multBinOm], fMultbins_Om[multBinOm + 1], (fNmultbins_Om - 1) - multBinOm), fNptbins_Om, fPtbins_Om);

        if (fisMC)
        {
            eff_omm[multBinOm] = new TH1D(TString::Format(("eff_omm[%d]"), multBinOm), TString::Format(("V0A %.0f-%.0f%%"), fMultbins_Om[multBinOm], fMultbins_Om[multBinOm + 1]), fNptbins_Om, fPtbins_Om);
            eff_omp[multBinOm] = new TH1D(TString::Format(("eff_omp[%d]"), multBinOm), TString::Format(("V0A %.0f-%.0f%%"), fMultbins_Om[multBinOm], fMultbins_Om[multBinOm + 1]), fNptbins_Om, fPtbins_Om);

            eff_omC[multBinOm] = new TH1D(TString::Format(("eff_omC[%d]"), multBinOm), TString::Format(("V0A %.0f-%.0f%%"), fMultbins_Om[multBinOm], fMultbins_Om[multBinOm + 1]), fNptbins_Om, fPtbins_Om);
        }
        for (Int_t ptBinOm = 0; ptBinOm < fNptbins_Om; ptBinOm++)
        {
            rawPt_omm[multBinOm]->SetBinContent(ptBinOm + 1, resultParOmm_pt_mult[ptBinOm][multBinOm]->GetBinContent(1) / ((rawPt_omm[multBinOm]->GetBinWidth(ptBinOm + 1)) * (h_multBinEntries_Om->GetBinContent(multBinOm + 1)))); // Bin 1 in resultparams is raw pt's bin counting
            rawPt_omm[multBinOm]->SetBinError(ptBinOm + 1, resultParOmm_pt_mult[ptBinOm][multBinOm]->GetBinError(1) / ((rawPt_omm[multBinOm]->GetBinWidth(ptBinOm + 1)) * (h_multBinEntries_Om->GetBinContent(multBinOm + 1))));
            rawPt_omp[multBinOm]->SetBinContent(ptBinOm + 1, resultParOmp_pt_mult[ptBinOm][multBinOm]->GetBinContent(1) / ((rawPt_omp[multBinOm]->GetBinWidth(ptBinOm + 1)) * (h_multBinEntries_Om->GetBinContent(multBinOm + 1)))); // Bin 1 in resultparams is raw pt's bin counting
            rawPt_omp[multBinOm]->SetBinError(ptBinOm + 1, resultParOmp_pt_mult[ptBinOm][multBinOm]->GetBinError(1) / ((rawPt_omp[multBinOm]->GetBinWidth(ptBinOm + 1)) * (h_multBinEntries_Om->GetBinContent(multBinOm + 1))));

            rawPt_omC[multBinOm]->SetBinContent(ptBinOm + 1, resultParOmC_pt_mult[ptBinOm][multBinOm]->GetBinContent(1) / ((rawPt_omC[multBinOm]->GetBinWidth(ptBinOm + 1)) * (h_multBinEntries_Om->GetBinContent(multBinOm + 1)))); // Bin 1 in resultparams is raw pt's bin counting
            rawPt_omC[multBinOm]->SetBinError(ptBinOm + 1, resultParOmC_pt_mult[ptBinOm][multBinOm]->GetBinError(1) / ((rawPt_omC[multBinOm]->GetBinWidth(ptBinOm + 1)) * (h_multBinEntries_Om->GetBinContent(multBinOm + 1))));

            if (fisMC)
            {
                eff_omm[multBinOm]->SetBinContent(ptBinOm + 1, resultParOmm_pt_mult[ptBinOm][multBinOm]->GetBinContent(7)); // Bin 7 in resultparams is MC efficiency
                eff_omm[multBinOm]->SetBinError(ptBinOm + 1, resultParOmm_pt_mult[ptBinOm][multBinOm]->GetBinError(7));

                eff_omp[multBinOm]->SetBinContent(ptBinOm + 1, resultParOmp_pt_mult[ptBinOm][multBinOm]->GetBinContent(7)); // Bin 7 in resultparams is MC efficiency
                eff_omp[multBinOm]->SetBinError(ptBinOm + 1, resultParOmp_pt_mult[ptBinOm][multBinOm]->GetBinError(7));

                eff_omC[multBinOm]->SetBinContent(ptBinOm + 1, resultParOmC_pt_mult[ptBinOm][multBinOm]->GetBinContent(7)); // Bin 7 in resultparams is MC efficiency
                eff_omC[multBinOm]->SetBinError(ptBinOm + 1, resultParOmC_pt_mult[ptBinOm][multBinOm]->GetBinError(7));
            }
            /// Scale for bin width: N->dN/dpt
            // rawPt_omm[multBinOm]->Scale(1, "width");
            // rawPt_omp[multBinOm]->Scale(1, "width");

            h_MassOmm_pt_mult[ptBinOm][multBinOm]->GetXaxis()->SetTitle("M_{#Lambda^{0}K^{#pm}} (GeV/ #it{c}^{2})");
            h_MassOmp_pt_mult[ptBinOm][multBinOm]->GetXaxis()->SetTitle("M_{#Lambda^{0}K^{#pm}} (GeV/ #it{c}^{2})");
            h_MassOmC_pt_mult[ptBinOm][multBinOm]->GetXaxis()->SetTitle("M_{#Lambda^{0}K^{#pm}} (GeV/ #it{c}^{2})");

            DrawAndSave(h_MassOmm_pt_mult[ptBinOm][multBinOm], h_bgOmm_pt_mult[ptBinOm][multBinOm], resultParOmm_pt_mult[ptBinOm][multBinOm], saveImages, outputFolder, imageFormat);
            DrawAndSave(h_MassOmp_pt_mult[ptBinOm][multBinOm], h_bgOmp_pt_mult[ptBinOm][multBinOm], resultParOmp_pt_mult[ptBinOm][multBinOm], saveImages, outputFolder, imageFormat);
            DrawAndSave(h_MassOmC_pt_mult[ptBinOm][multBinOm], h_bgOmC_pt_mult[ptBinOm][multBinOm], resultParOmC_pt_mult[ptBinOm][multBinOm], saveImages, outputFolder, imageFormat);
        }
        outputFile->cd("dirRawPt_omp");
        rawPt_omp[multBinOm]->SetMarkerStyle(markerStyles[multBinOm]);
        rawPt_omp[multBinOm]->Write();
        rawPt_omp[multBinOm]->SetName(TString::Format(("V0A %.0f-%.0f%% #times2^{%d}"), fMultbins_Om[multBinOm], fMultbins_Om[multBinOm + 1], (fNmultbins_Om - 1) - multBinOm));
        rawPt_omp[multBinOm]->Scale(pow(2, (fNmultbins_Om - 1) - multBinOm));
        hs_omp->Add(rawPt_omp[multBinOm]);

        if (fisMC)
        {
            eff_omp[multBinOm]->SetMarkerStyle(markerStyles[multBinOm]);
            eff_omp[multBinOm]->Write();
            eff_omp[multBinOm]->SetName(TString::Format(("V0A %.0f-%.0f%%"), fMultbins_Om[multBinOm], fMultbins_Om[multBinOm + 1]));
            // eff_omp[multBinOm]->Scale(pow(2, (fNmultbins_Om - 1) - multBinOm));
            hs_omp_eff->Add(eff_omp[multBinOm]);
        }

        outputFile->cd("dirRawPt_omm");
        rawPt_omm[multBinOm]->SetMarkerStyle(markerStyles[multBinOm]);
        rawPt_omm[multBinOm]->Write();
        rawPt_omm[multBinOm]->SetName(TString::Format(("V0A %.0f-%.0f%% #times2^{%d}"), fMultbins_Om[multBinOm], fMultbins_Om[multBinOm + 1], (fNmultbins_Om - 1) - multBinOm));
        rawPt_omm[multBinOm]->Scale(pow(2, (fNmultbins_Om - 1) - multBinOm));
        hs_omm->Add(rawPt_omm[multBinOm]);

        if (fisMC)
        {
            eff_omm[multBinOm]->SetMarkerStyle(markerStyles[multBinOm]);
            eff_omm[multBinOm]->Write();
            eff_omm[multBinOm]->SetName(TString::Format(("V0A %.0f-%.0f%%"), fMultbins_Om[multBinOm], fMultbins_Om[multBinOm + 1]));
            // eff_omm[multBinOm]->Scale(pow(2, (fNmultbins_Om - 1) - multBinOm));
            hs_omm_eff->Add(eff_omm[multBinOm]);
        }

        outputFile->cd("dirRawPt_omC");
        rawPt_omC[multBinOm]->SetMarkerStyle(markerStyles[multBinOm]);
        rawPt_omC[multBinOm]->Write();
        rawPt_omC[multBinOm]->SetName(TString::Format(("V0A %.0f-%.0f%% #times2^{%d}"), fMultbins_Om[multBinOm], fMultbins_Om[multBinOm + 1], (fNmultbins_Om - 1) - multBinOm));
        rawPt_omC[multBinOm]->Scale(pow(2, (fNmultbins_Om - 1) - multBinOm));
        hs_omC->Add(rawPt_omC[multBinOm]);

        if (fisMC)
        {
            eff_omC[multBinOm]->SetMarkerStyle(markerStyles[multBinOm]);
            eff_omC[multBinOm]->Write();
            eff_omC[multBinOm]->SetName(TString::Format(("V0A %.0f-%.0f%%"), fMultbins_Om[multBinOm], fMultbins_Om[multBinOm + 1]));
            // eff_omC[multBinOm]->Scale(pow(2, (fNmultbins_Om - 1) - multBinOm));
            hs_omC_eff->Add(eff_omC[multBinOm]);

            eff_omC_ratio[multBinOm] = new TH1D(TString::Format(("eff_omC_ratio[%d]"), multBinOm), "", fNptbins_Om, fPtbins_Om);
            eff_omC_ratio[multBinOm]->Divide(eff_omC[multBinOm], eff_pt_omC);
            eff_omC_ratio[multBinOm]->SetMarkerStyle(markerStyles[multBinOm]);
            eff_omC_ratio[multBinOm]->Write();
            eff_omC_ratio[multBinOm]->SetName(TString::Format(("V0A %.0f-%.0f%%"), fMultbins_Om[multBinOm], fMultbins_Om[multBinOm + 1]));
            hs_omC_eff_ratio->Add(eff_omC_ratio[multBinOm]);
        }
    }

    gStyle->SetOptStat(0);
    gStyle->SetPalette(kVisibleSpectrum);

    TCanvas *cPtSpectra[8];
    for (Int_t iCanvas = 0; iCanvas < 8; iCanvas++)
    {
        cPtSpectra[iCanvas] = new TCanvas(TString::Format("cPtSpectra%d", iCanvas), TString::Format("cPtSpectra%d", iCanvas), 1200, 900);
    }

    PaintStack_Draw(*cPtSpectra[0], *hs_xip);
    PaintStack_Draw(*cPtSpectra[1], *hs_omp);
    PaintStack_Draw(*cPtSpectra[2], *hs_xim);
    PaintStack_Draw(*cPtSpectra[3], *hs_omm);
    PaintStack_Draw(*cPtSpectra[4], *hs_xiC);
    PaintStack_Draw(*cPtSpectra[5], *hs_omC);

    outputFile->cd();

    hs_xip->Write();
    hs_xim->Write();
    hs_omp->Write();
    hs_omm->Write();
    hs_xiC->Write();
    hs_omC->Write();

    if (saveStack)
    {
        cPtSpectra[0]->cd();
        SaveImage(outputFolder, "rawPt", "xipN", imageFormat.Data(), cPtSpectra[0]);

        cPtSpectra[1]->cd();
        SaveImage(outputFolder, "rawPt", "ompN", imageFormat.Data(), cPtSpectra[1]);

        cPtSpectra[2]->cd();
        SaveImage(outputFolder, "rawPt", "ximN", imageFormat.Data(), cPtSpectra[2]);

        cPtSpectra[3]->cd();
        SaveImage(outputFolder, "rawPt", "ommN", imageFormat.Data(), cPtSpectra[3]);

        cPtSpectra[4]->cd();
        SaveImage(outputFolder, "rawPt", "xicN", imageFormat.Data(), cPtSpectra[4]);

        cPtSpectra[5]->cd();
        SaveImage(outputFolder, "rawPt", "omcN", imageFormat.Data(), cPtSpectra[5]);
    }

    if (fisMC)
    {
        PaintStack_Draw(*cPtSpectra[0], *hs_xip_eff, kFALSE, "Efficiency");
        PaintStack_Draw(*cPtSpectra[1], *hs_omp_eff, kFALSE, "Efficiency");
        PaintStack_Draw(*cPtSpectra[2], *hs_xim_eff, kFALSE, "Efficiency");
        PaintStack_Draw(*cPtSpectra[3], *hs_omm_eff, kFALSE, "Efficiency");
        PaintStack_Draw(*cPtSpectra[4], *hs_xiC_eff, kFALSE, "Efficiency");
        PaintStack_Draw(*cPtSpectra[5], *hs_omC_eff, kFALSE, "Efficiency");

        PaintStack_Draw(*cPtSpectra[6], *hs_xiC_eff_ratio, kFALSE, "Multiplicity classes/0-100%");
        PaintStack_Draw(*cPtSpectra[7], *hs_omC_eff_ratio, kFALSE, "Multiplicity classes/0-100%");

        hs_xip_eff->Write();
        hs_xim_eff->Write();
        hs_omp_eff->Write();
        hs_omm_eff->Write();
        hs_xiC_eff->Write();
        hs_omC_eff->Write();

        hs_xiC_eff_ratio->Write();
        hs_omC_eff_ratio->Write();

        if (saveStack)
        {
            cPtSpectra[0]->cd();
            SaveImage(outputFolder, "eff", "eff_xipN", imageFormat.Data(), cPtSpectra[0]);

            cPtSpectra[1]->cd();
            SaveImage(outputFolder, "eff", "eff_ompN", imageFormat.Data(), cPtSpectra[1]);

            cPtSpectra[2]->cd();
            SaveImage(outputFolder, "eff", "eff_ximN", imageFormat.Data(), cPtSpectra[2]);

            cPtSpectra[3]->cd();
            SaveImage(outputFolder, "eff", "eff_ommN", imageFormat.Data(), cPtSpectra[3]);

            cPtSpectra[4]->cd();
            SaveImage(outputFolder, "eff", "eff_xicN", imageFormat.Data(), cPtSpectra[4]);

            cPtSpectra[5]->cd();
            SaveImage(outputFolder, "eff", "eff_omcN", imageFormat.Data(), cPtSpectra[5]);

            cPtSpectra[6]->cd();
            SaveImage(outputFolder, "effRatio", "effRatio_xicN", imageFormat.Data(), cPtSpectra[6]);

            cPtSpectra[7]->cd();
            SaveImage(outputFolder, "effRatio", "effRatio_omcN", imageFormat.Data(), cPtSpectra[7]);
        }
    }

    if (!ptRatioFilename.IsNull())
    {
        /// get previously calculated pT ratio from file to compute ratio of raw pT spectra

        Printf("\nOpening %s for ratio calculation...", ptRatioFilename.Data());
        TFile *fRatio = OpenFile(ptRatioFilename);

        TString yAxisTitle = outputFilename(outputFilename.Last('/') + 1, outputFilename.Length()) + "/" + ptRatioFilename(ptRatioFilename.Last('/') + 1, ptRatioFilename.Length());
        yAxisTitle.ReplaceAll(".root", "");
        yAxisTitle.ReplaceAll("_Draw", "");

        // needed so we can do file->Close()
        TH1::AddDirectory(0);
        TH1 *rawPt_xim_compare[fNmultbins_Xi];
        TH1 *rawPt_xip_compare[fNmultbins_Xi];
        TH1 *rawPt_omm_compare[fNmultbins_Om];
        TH1 *rawPt_omp_compare[fNmultbins_Om];
        TH1 *rawPt_xiC_compare[fNmultbins_Xi];
        TH1 *rawPt_omC_compare[fNmultbins_Om];

        TH1D *ratioPt_xim[fNmultbins_Xi];
        TH1D *ratioPt_xip[fNmultbins_Xi];
        TH1D *ratioPt_omm[fNmultbins_Om];
        TH1D *ratioPt_omp[fNmultbins_Om];
        TH1D *ratioPt_xiC[fNmultbins_Xi];
        TH1D *ratioPt_omC[fNmultbins_Om];

        auto hs_ratio_xip = new THStack("hs_ratio_xip", "#it{p}_{T} spectra ratio #Xi^{+}");
        auto hs_ratio_xim = new THStack("hs_ratio_xim", "#it{p}_{T} spectra ratio #Xi^{-}");
        auto hs_ratio_omp = new THStack("hs_ratio_omp", "#it{p}_{T} spectra ratio #Omega^{+}");
        auto hs_ratio_omm = new THStack("hs_ratio_omm", "#it{p}_{T} spectra ratio #Omega^{-}");
        auto hs_ratio_xiC = new THStack("hs_ratio_xiC", "#it{p}_{T} spectra ratio #Xi^{+} + #Xi^{-}");
        auto hs_ratio_omC = new THStack("hs_ratio_omC", "#it{p}_{T} spectra ratio #Omega^{+} + #Omega^{-}");

        for (Int_t multBinXi = 0; multBinXi < fNmultbins_Xi; multBinXi++)
        {
            rawPt_xim_compare[multBinXi] = (TH1 *)fRatio->FindObjectAny(TString::Format(("rawPt_xim[%d]"), multBinXi));
            rawPt_xip_compare[multBinXi] = (TH1 *)fRatio->FindObjectAny(TString::Format(("rawPt_xip[%d]"), multBinXi));
            rawPt_xiC_compare[multBinXi] = (TH1 *)fRatio->FindObjectAny(TString::Format(("rawPt_xiC[%d]"), multBinXi));
        }

        for (Int_t multBinOm = 0; multBinOm < fNmultbins_Om; multBinOm++)
        {
            rawPt_omm_compare[multBinOm] = (TH1 *)fRatio->FindObjectAny(TString::Format(("rawPt_omm[%d]"), multBinOm));
            rawPt_omp_compare[multBinOm] = (TH1 *)fRatio->FindObjectAny(TString::Format(("rawPt_omp[%d]"), multBinOm));
            rawPt_omC_compare[multBinOm] = (TH1 *)fRatio->FindObjectAny(TString::Format(("rawPt_omC[%d]"), multBinOm));
        } // input ended for fRatio

        delete fRatio;

        /// calculate ratio of current spectra to previous
        for (Int_t multBinXi = 0; multBinXi < fNmultbins_Xi; multBinXi++)
        {
            /// get current histogram scale back to normal to perform divide
            rawPt_xim[multBinXi]->Scale(pow(2, -((fNmultbins_Xi - 1) - multBinXi)));
            rawPt_xip[multBinXi]->Scale(pow(2, -((fNmultbins_Xi - 1) - multBinXi)));
            rawPt_xiC[multBinXi]->Scale(pow(2, -((fNmultbins_Xi - 1) - multBinXi)));

            /// generate ratio hists
            ratioPt_xim[multBinXi] = new TH1D(TString::Format(("ratioPt_xim[%d]"), multBinXi), "", fNptbins_Xi, fPtbins_Xi);
            ratioPt_xip[multBinXi] = new TH1D(TString::Format(("ratioPt_xip[%d]"), multBinXi), "", fNptbins_Xi, fPtbins_Xi);
            ratioPt_xiC[multBinXi] = new TH1D(TString::Format(("ratioPt_xiC[%d]"), multBinXi), "", fNptbins_Xi, fPtbins_Xi);

            ratioPt_xim[multBinXi]->Divide(rawPt_xim[multBinXi], rawPt_xim_compare[multBinXi]);
            outputFile->cd("dirRawPt_xim");
            ratioPt_xim[multBinXi]->SetMarkerStyle(markerStyles[multBinXi]);
            ratioPt_xim[multBinXi]->Write();
            ratioPt_xim[multBinXi]->SetName(TString::Format(("V0A %.0f-%.0f%%"), fMultbins_Xi[multBinXi], fMultbins_Xi[multBinXi + 1]));
            hs_ratio_xim->Add(ratioPt_xim[multBinXi]);

            ratioPt_xip[multBinXi]->Divide(rawPt_xip[multBinXi], rawPt_xip_compare[multBinXi]);
            outputFile->cd("dirRawPt_xip");
            ratioPt_xip[multBinXi]->SetMarkerStyle(markerStyles[multBinXi]);
            ratioPt_xip[multBinXi]->Write();
            ratioPt_xip[multBinXi]->SetName(TString::Format(("V0A %.0f-%.0f%%"), fMultbins_Xi[multBinXi], fMultbins_Xi[multBinXi + 1]));
            hs_ratio_xip->Add(ratioPt_xip[multBinXi]);

            ratioPt_xiC[multBinXi]->Divide(rawPt_xiC[multBinXi], rawPt_xiC_compare[multBinXi]);
            outputFile->cd("dirRawPt_xiC");
            ratioPt_xiC[multBinXi]->SetMarkerStyle(markerStyles[multBinXi]);
            ratioPt_xiC[multBinXi]->Write();
            ratioPt_xiC[multBinXi]->SetName(TString::Format(("V0A %.0f-%.0f%%"), fMultbins_Xi[multBinXi], fMultbins_Xi[multBinXi + 1]));
            hs_ratio_xiC->Add(ratioPt_xiC[multBinXi]);
        }

        for (Int_t multBinOm = 0; multBinOm < fNmultbins_Om; multBinOm++)
        {
            /// get current histogram scale back to normal to perform divide
            rawPt_omm[multBinOm]->Scale(pow(2, -((fNmultbins_Om - 1) - multBinOm)));
            rawPt_omp[multBinOm]->Scale(pow(2, -((fNmultbins_Om - 1) - multBinOm)));
            rawPt_omC[multBinOm]->Scale(pow(2, -((fNmultbins_Om - 1) - multBinOm)));

            /// generate ratio hists
            ratioPt_omm[multBinOm] = new TH1D(TString::Format(("ratioPt_omm[%d]"), multBinOm), "", fNptbins_Om, fPtbins_Om);
            ratioPt_omp[multBinOm] = new TH1D(TString::Format(("ratioPt_omp[%d]"), multBinOm), "", fNptbins_Om, fPtbins_Om);
            ratioPt_omC[multBinOm] = new TH1D(TString::Format(("ratioPt_omC[%d]"), multBinOm), "", fNptbins_Om, fPtbins_Om);

            ratioPt_omm[multBinOm]->Divide(rawPt_omm[multBinOm], rawPt_omm_compare[multBinOm]);
            outputFile->cd("dirRawPt_omm");
            ratioPt_omm[multBinOm]->SetMarkerStyle(markerStyles[multBinOm]);
            ratioPt_omm[multBinOm]->Write();
            ratioPt_omm[multBinOm]->SetName(TString::Format(("V0A %.0f-%.0f%%"), fMultbins_Om[multBinOm], fMultbins_Om[multBinOm + 1]));
            hs_ratio_omm->Add(ratioPt_omm[multBinOm]);

            ratioPt_omp[multBinOm]->Divide(rawPt_omp[multBinOm], rawPt_omp_compare[multBinOm]);
            outputFile->cd("dirRawPt_omp");
            ratioPt_omp[multBinOm]->SetMarkerStyle(markerStyles[multBinOm]);
            ratioPt_omp[multBinOm]->Write();
            ratioPt_omp[multBinOm]->SetName(TString::Format(("V0A %.0f-%.0f%%"), fMultbins_Om[multBinOm], fMultbins_Om[multBinOm + 1]));
            hs_ratio_omp->Add(ratioPt_omp[multBinOm]);

            ratioPt_omC[multBinOm]->Divide(rawPt_omC[multBinOm], rawPt_omC_compare[multBinOm]);
            outputFile->cd("dirRawPt_omC");
            ratioPt_omC[multBinOm]->SetMarkerStyle(markerStyles[multBinOm]);
            ratioPt_omC[multBinOm]->Write();
            ratioPt_omC[multBinOm]->SetName(TString::Format(("V0A %.0f-%.0f%%"), fMultbins_Om[multBinOm], fMultbins_Om[multBinOm + 1]));
            hs_ratio_omC->Add(ratioPt_omC[multBinOm]);
        }

        PaintStack_Draw(*cPtSpectra[0], *hs_ratio_xip, kFALSE, yAxisTitle);
        PaintStack_Draw(*cPtSpectra[1], *hs_ratio_omp, kFALSE, yAxisTitle);
        PaintStack_Draw(*cPtSpectra[2], *hs_ratio_xim, kFALSE, yAxisTitle);
        PaintStack_Draw(*cPtSpectra[3], *hs_ratio_omm, kFALSE, yAxisTitle);
        PaintStack_Draw(*cPtSpectra[4], *hs_ratio_xiC, kFALSE, yAxisTitle);
        PaintStack_Draw(*cPtSpectra[5], *hs_ratio_omC, kFALSE, yAxisTitle);

        /// Draw a line at y=1 for ratio histStack
        DrawRatioLine(*cPtSpectra[0], *hs_ratio_xip, 0.5, 1.5);
        DrawRatioLine(*cPtSpectra[1], *hs_ratio_omp, 0.5, 1.5);
        DrawRatioLine(*cPtSpectra[2], *hs_ratio_xim, 0.5, 1.5);
        DrawRatioLine(*cPtSpectra[3], *hs_ratio_omm, 0.5, 1.5);
        DrawRatioLine(*cPtSpectra[4], *hs_ratio_xiC, 0.5, 1.5);
        DrawRatioLine(*cPtSpectra[5], *hs_ratio_omC, 0.5, 1.5);

        outputFile->cd();
        hs_ratio_xip->Write();
        hs_ratio_omp->Write();
        hs_ratio_xim->Write();
        hs_ratio_omm->Write();
        hs_ratio_xiC->Write();
        hs_ratio_omC->Write();

        if (saveStack)
        {
            cPtSpectra[0]->cd();
            SaveImage(outputFolder, "ratioPt", "ratio_xipN", imageFormat.Data(), cPtSpectra[0]);

            cPtSpectra[1]->cd();
            SaveImage(outputFolder, "ratioPt", "ratio_ompN", imageFormat.Data(), cPtSpectra[1]);

            cPtSpectra[2]->cd();
            SaveImage(outputFolder, "ratioPt", "ratio_ximN", imageFormat.Data(), cPtSpectra[2]);

            cPtSpectra[3]->cd();
            SaveImage(outputFolder, "ratioPt", "ratio_ommN", imageFormat.Data(), cPtSpectra[3]);

            cPtSpectra[4]->cd();
            SaveImage(outputFolder, "ratioPt", "ratio_xicN", imageFormat.Data(), cPtSpectra[4]);

            cPtSpectra[5]->cd();
            SaveImage(outputFolder, "ratioPt", "ratio_omcN", imageFormat.Data(), cPtSpectra[5]);
        }
    }

    // delete canvas objects
    for (Int_t iCanvas = 0; iCanvas < 8; iCanvas++)
    {
        delete cPtSpectra[iCanvas];
        // delete gROOT->FindObject(TString::Format("c%d", iCanvas));
    }
    delete outputFile;
    return 0;
}

void DrawAndSave(TH1 *peak, TH1 *bg, TH1 *resultParams, Bool_t saveImages, TString outputFolder, TString imageFormat)
{
    TString peakName = peak->GetName();
    TString peakTitle = peak->GetTitle();
    TCanvas *cDraw = new TCanvas(peakName, peakTitle, 1200, 900);
    cDraw->cd();

    Double_t pPosition = resultParams->GetBinContent(4); // avg mean for DGaus fit is stored in resultParams bin 4
    Double_t pWidth = resultParams->GetBinContent(5);    // avg sigma for DGaus fit is stored in resultParams bin 5
    // auto legend = new TLegend(0.1, 0.7, 0.28, 0.9);
    // legend->SetHeader("Fit Stats", "C"); // option "C" allows to center the header
    // legend->AddEntry(peak->GetListOfFunctions()->At(0), "", "l");
    // legend->AddEntry(peak->GetListOfFunctions()->At(0), TString::Format("Fit mean = %.3f +/- %.3f", pPosition, resultParams->GetBinError(4)), "l");
    // legend->AddEntry(peak->GetListOfFunctions()->At(0), TString::Format("Fit sigma = %.3f +/- %.3f", pWidth, resultParams->GetBinError(5)), "l");
    // legend->AddEntry(peak->GetListOfFunctions()->At(1), "", "lpf");
    // legend->AddEntry(peak, TString::Format("Sig - Bg (BC-FF) = %.3f +/- %.3f", resultParams->GetBinContent(1), resultParams->GetBinError(1)), "pe");

    /// Legend for final plots:
    peak->SetStats(0);
    auto legend = new TLegend(0.72, 0.72, 0.88, 0.88);
    legend->SetBorderSize(0);
    legend->SetHeader("Invariant Mass"); // option "C" allows to center the header
    // legend->AddEntry(peak, "p-Pb #sqrt{s_{NN}} = 8.16 TeV", "");
    // legend->AddEntry(peak, TString::Format("%s, %s", ptRange.Data(), multRange.Data()), "");
    legend->AddEntry(peak, "Data", "lpe");
    legend->AddEntry(peak->GetListOfFunctions()->At(1), "Background", "lpf");
    legend->AddEntry(peak->GetListOfFunctions()->At(0), "Signal Fit", "l");

    // TString plotTitle = peakTitle(0, peakTitle.Index(":"));
    TString partTitle = "#bf{" + peakTitle(15, peakTitle.Index(":") - 15);
    if (!partTitle.EndsWith("}"))
        partTitle += "^{#pm}}";
    else
        partTitle += "}";
    TString collInfo = "p-Pb #sqrt{s_{NN}} = 8.16 TeV, This work";
    TString ptRange = peakTitle(peakTitle.Index(":") + 8, peakTitle.Index(", ") - peakTitle.Index(":") - 9) + " GeV/#it{c}";
    TString multRange = peakTitle(peakTitle.Index(", ") + 7, peakTitle.Length() - peakTitle.Index(", ") - 8) + "% V0A";
    ptRange.ReplaceAll(",", " < #it{p}_{T} < ");
    multRange.ReplaceAll(",", "-");
    TLatex *latexPartTitle = new TLatex(0.15, 0.84, partTitle.Data());
    TLatex *latexCollInfo = new TLatex(0.15, 0.80, collInfo.Data());
    TLatex *latexPtRange = new TLatex(0.15, 0.76, ptRange.Data());
    TLatex *latexMultRange = new TLatex(0.15, 0.72, multRange.Data());
    latexPartTitle->SetTextSize(0.03);
    latexCollInfo->SetTextSize(0.025);
    latexPtRange->SetTextSize(0.025);
    latexMultRange->SetTextSize(0.025);
    latexPartTitle->SetNDC(kTRUE);
    latexCollInfo->SetNDC(kTRUE);
    latexPtRange->SetNDC(kTRUE);
    latexMultRange->SetNDC(kTRUE);

    peak->SetTitle("");
    if (peakName.Contains("Xi"))
        peak->GetYaxis()->SetTitle("Counts / (1 MeV/ #it{c}^{2})");
    else
        peak->GetYaxis()->SetTitle("Counts / (2 MeV/ #it{c}^{2})");

    ///  Defining peak limits for signal region (magenta lines):
    ///   par[1] = peak position, par[2] = peak width
    Double_t lPeakLeftLimit = pPosition - 1. * 4 * TMath::Abs(pWidth);
    Double_t lPeakRightLimit = pPosition + 1. * 4 * TMath::Abs(pWidth);

    peak->Draw();
    if (bg)
        bg->Draw("same");
    cDraw->Update();
    TLine *lLineLeft = new TLine(lPeakLeftLimit, gPad->GetUymin(), lPeakLeftLimit, gPad->GetUymax());
    TLine *lLineRight = new TLine(lPeakRightLimit, gPad->GetUymin(), lPeakRightLimit, gPad->GetUymax());
    lLineLeft->SetLineColor(kMagenta);
    lLineLeft->SetLineStyle(kDashed);
    lLineRight->SetLineColor(kMagenta);
    lLineRight->SetLineStyle(kDashed);
    lLineLeft->Draw("same");
    lLineRight->Draw();
    legend->Draw();
    latexPartTitle->Draw();
    latexCollInfo->Draw();
    latexPtRange->Draw();
    latexMultRange->Draw();

    if (saveImages)
    {
        TString imageFolder = peakName;
        if (imageFolder.Contains("["))
            imageFolder.Remove(imageFolder.First('[')); /// getting substring for folder naming purpose
        else
            imageFolder = "h_allInt";
        SaveImage(outputFolder, imageFolder, peakName, imageFormat, cDraw);
    }

    delete legend;
    delete lLineLeft;
    delete lLineRight;
    delete cDraw;
    return;
}

void DrawRatioLine(TCanvas &c, THStack &hs, Double_t hs_yMin, Double_t hs_yMax)
{
    if (hs_yMin != -100)
        hs.SetMinimum(hs_yMin);
    if (hs_yMax != -100)
        hs.SetMaximum(hs_yMax);

    hs.GetYaxis()->SetNdivisions(3, 5, 10);
    TLine *lLineAt1 = new TLine(hs.GetXaxis()->GetXmin(), 1, hs.GetXaxis()->GetXmax(), 1);
    lLineAt1->SetLineColor(kRed);
    c.cd();
    lLineAt1->Draw("same");
    c.Modified();
    c.Update();
    return;
}