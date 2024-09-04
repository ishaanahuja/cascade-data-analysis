#include <TROOT.h>
#include <THashList.h>

#include <TFitter.h>
#include <TFitResultPtr.h>
// #include <TSpectrum.h>
#include <Math/MinimizerOptions.h>
#include <Math/IntegratorOptions.h>
#include <TStyle.h>

#include "CascadeUtils.h"

TH1 *GenerateBg(bool bgTSpectrum, TH1 *h_source);
TF1 *SetFitParametersGaus(TF1 *peakFnc, Double_t *previousPtBinFitParams, TH1 *peak, Double_t mass, Double_t sigma, Bool_t usePrevious);
TF1 *SetFitParametersDG(TF1 *peakFnc, Double_t *previousPtBinFitParams, TH1 *peak, Double_t mass, Double_t sigma, Bool_t usePrevious);

TH1 *FitResults(TH1 *h_bg, TF1 *peakFnc, TH1 *peak, Double_t fitMinSig, Double_t fitMaxSig, Int_t binsMC[], TH2 *MCgen, bool fisMC, Bool_t isGausPol2, TString fitOptSig, TString fitOptBg);

TH1 *FillParams(TH1 *h_bg, TFitResultPtr fFitResult_bg, TF1 *sigBgFnc, TH1 *peak, TFitResultPtr fFitResult, Double_t minInt, Double_t maxInt, TH2 *MCgen, Int_t binsMC[], bool fisMC, bool isGausPol2);

void GetYieldBinCounting(TH1 *h_bg, TH1 *h, TFitResultPtr fFitResult, Double_t minInt, Double_t maxInt, bool fisMC, Double_t &val,
                         Double_t &err, Double_t eps = 1e-2);

void GetYieldFitFunction(TH1 *h, TFitResultPtr fFitResult_sig, TFitResultPtr fFitResult_bg, Double_t minInt, Double_t maxInt, bool fisMC, Double_t &val,
                         Double_t &err, Double_t eps = 1e-2);


int FitCascades(TString inputFilename = "050824_updatedCutVar_6Runs.root", TString outputFilename = "030924_fitUpdatedCuts_def.root", TString h3Name = "h3_ptmasscent_def", TString idAxis = "ye", bool fisMC = false, Int_t fitFunction = 2, bool xi = true, bool om = true, Int_t verbosity = kWarning)
{
    ROOT::EnableImplicitMT();
    ROOT::Math::IntegratorOneDimOptions::SetDefaultIntegrator("Adaptive");
    ROOT::Math::IntegratorOneDimOptions::SetDefaultNPoints(6);
    ROOT::Math::IntegratorOneDimOptions::SetDefaultAbsTolerance(1E-2);
    ROOT::Math::IntegratorOneDimOptions::SetDefaultRelTolerance(1E-2);
    ROOT::Math::MinimizerOptions::SetDefaultMaxFunctionCalls(100000);
    ROOT::Math::MinimizerOptions::SetDefaultMinimizer("Minuit2");
    ROOT::Math::MinimizerOptions::SetDefaultStrategy(2);
    ROOT::Math::MinimizerOptions::SetDefaultPrintLevel(0); // Fit printing: -1 = no printing, 0 (minimal) to 3 (max)
    gStyle->SetOptFit(1111);

    TString fitOptSig = "NMSL+QIE MULTITHREAD";
    TString fitOptBg = "NSB+QLI MULTITHREAD";

    ///special cases where ROOT hangs while fitting variations with option 'E'
    if (h3Name.EqualTo("h3Var_DcaV0Daughters[6][8]")) 
        fitOptSig = "NMSL+QI MULTITHREAD";

    TFitter::SetMaxIterations(100000);
    TVirtualFitter::SetMaxIterations(100000);
    TFitter::SetPrecision(1E-2);
    TH1::AddDirectory(0);

    // gPrintViaErrorHandler = kTRUE;
    gErrorIgnoreLevel = verbosity;

    Bool_t isGausPol2 = kFALSE;
    Bool_t isDoubleGausPol2 = kFALSE;
    // Bool_t bgTSpectrum = kFALSE;
    /// Setting which fit functions to use:
    /// fitFunctions = [GausPol2, DoubleGausPol2, bgTSpectrum]
    /// e.g. fitFunctions = 110 -> fit with GausPol2 and DoubleGausPol2, don't fit bgTSpectrum

    if (fitFunction == 1)
    {
        isGausPol2 = kTRUE;
        Printf("<%s>: Selected fit function: GausPol2", h3Name.Data());
    }
    else if (fitFunction == 2)
    {
        isDoubleGausPol2 = kTRUE;
        Printf("<%s>: Selected fit function: DoubleGausPol2", h3Name.Data());
    }
    // if (fitFunction == '3')
    // {
    //     bgTSpectrum = kTRUE;
    //     Printf("Background fitting with TSpectrum will be used.");
    // }
    // }
    else
    {
        Fatal("FitCascades: <int> fitFunction", "Provided value %d is not valid. fitFunctions = [1:GausPol2, 2:DoubleGausPol2]", fitFunction);
        //  \ne.g. fitFunctions = \"110\" -> fit with GausPol2 and DoubleGausPol2, don't fit bgTSpectrum
        return 3;
    }

    /// declaring histogram objects to read:
    TH3 *Xim, *Xip, *Omm, *Omp;
    TH2 *MCGenXip, *MCGenXim, *MCGenOmp, *MCGenOmm;
    TH2D *MCGenXiC = nullptr;
    TH2D *MCGenOmC = nullptr;

    /// Read input file:
    TFile *inputFile = OpenFile(inputFilename);

    THashList *list_eve = (THashList *)inputFile->FindObjectAny("chists_eve_");
    TH1 *cent = (TH1 *)list_eve->FindObject("hcent");

    if (xi)
    {
        THashList *list_Xim = (THashList *)inputFile->FindObjectAny("chists_Xim_");
        THashList *list_Xip = (THashList *)inputFile->FindObjectAny("chists_Xip_");
        Xim = (TH3 *)list_Xim->FindObject(h3Name.Data());
        Xip = (TH3 *)list_Xip->FindObject(h3Name.Data());
        MCGenXip = (TH2 *)list_Xip->FindObject("h2_gen");
        MCGenXim = (TH2 *)list_Xim->FindObject("h2_gen");

        if (!(Xim && Xip))
        {
            Error("FitCascades: histInput", "XI: Histograms '%s' cannot be found. Skipping XI analysis ...", h3Name.Data());
            xi = false;
        }

        else if ((Xim->GetEntries() < 1) || (Xip->GetEntries() < 1))
        {
            Error("FitCascades: histInput", "XI: Histograms '%s' are empty. Skipping XI analysis ...", h3Name.Data());
            xi = false;
        }
    }

    if (om)
    {
        THashList *list_Omm = (THashList *)inputFile->FindObjectAny("chists_Omm_");
        THashList *list_Omp = (THashList *)inputFile->FindObjectAny("chists_Omp_");
        Omm = (TH3 *)list_Omm->FindObject(h3Name.Data());
        Omp = (TH3 *)list_Omp->FindObject(h3Name.Data());
        MCGenOmp = (TH2 *)list_Omp->FindObject("h2_gen");
        MCGenOmm = (TH2 *)list_Omm->FindObject("h2_gen");

        if (!(Omm && Omp))
        {
            Error("FitCascades: histInput", "OMEGA: Histograms '%s' cannot be found. Skipping OMEGA analysis ...", h3Name.Data());
            om = false;
        }

        else if ((Omm->GetEntries() < 1) || (Omp->GetEntries() < 1))
        {
            Error("FitCascades: histInput", "OMEGA: Histograms '%s' are empty. Skipping OMEGA analysis ...", h3Name.Data());
            om = false;
        }
    }

    if (!(xi || om)) // if none of the analyses can be done
    {
        Error("FitCascades: histInput", "%s: Both particle analysis flags are false. Aborting.", h3Name.Data());
        delete inputFile;
        return 1;
    }

    delete inputFile;

    // if (MCGenXip && MCGenXim)
    // {
    //     MCGenXiC = (TH2D *)MCGenXip->Clone();
    //     MCGenXiC->Add(MCGenXim);
    // }
    // if (MCGenOmp && MCGenOmm)
    // {
    //     MCGenOmC = (TH2D *)MCGenOmp->Clone();
    //     MCGenOmC->Add(MCGenOmm);
    // }

    Int_t binsMC[4] = {0};

    if (outputFilename.IsNull())
    {
        outputFilename = inputFilename;
        outputFilename.ReplaceAll(".root", "_Fitting.root");
    }

    TFile *outputFile = OpenFile(outputFilename, "RECREATE");
    outputFile->mkdir("_allInt");
    if (xi)
    {
        outputFile->mkdir("h_MassXim_pt_mult");
        outputFile->mkdir("h_MassXip_pt_mult");
        outputFile->mkdir("h_MassXiC_pt_mult"); /// for combination of + and -

        outputFile->mkdir("h_MassXim_pt"); // mult: 0-100%
        outputFile->mkdir("h_MassXip_pt");
        outputFile->mkdir("h_MassXiC_pt");
    }

    if (om)
    {
        outputFile->mkdir("h_MassOmm_pt_mult");
        outputFile->mkdir("h_MassOmp_pt_mult");
        outputFile->mkdir("h_MassOmC_pt_mult");

        outputFile->mkdir("h_MassOmm_pt");
        outputFile->mkdir("h_MassOmp_pt");
        outputFile->mkdir("h_MassOmC_pt");
    }

    /// BEGIN:

    Double_t nSigma = 10.;
    Double_t sigmaXi = 0.0018;
    Double_t sigmaOm = 0.0020;

    // Double_t sigma = 0.0025;
    // Double_t sigmaXi = 0.001875;
    // Double_t sigmaOm = 0.0018;

    TH1I *h_multBinEntries_Xi = new TH1I("h_multBinEntries_Xi", "h_multBinEntries_Xi", fNmultbins_Xi, fMultbins_Xi);
    TH1I *h_multBinEntries_Om = new TH1I("h_multBinEntries_Om", "h_multBinEntries_Om", fNmultbins_Om, fMultbins_Om);

    for (Int_t multBinXi = 0; multBinXi < fNmultbins_Xi; multBinXi++)
    {
        h_multBinEntries_Xi->SetBinContent(multBinXi + 1, cent->Integral(cent->FindBin(fMultbins_Xi[multBinXi]), cent->FindBin(fMultbins_Xi[multBinXi + 1])));
    }

    for (Int_t multBinOm = 0; multBinOm < fNmultbins_Om; multBinOm++)
    {
        h_multBinEntries_Om->SetBinContent(multBinOm + 1, cent->Integral(cent->FindBin(fMultbins_Om[multBinOm]), cent->FindBin(fMultbins_Om[multBinOm + 1])));
    }
    outputFile->cd();
    h_multBinEntries_Xi->Write();
    h_multBinEntries_Om->Write();

    TH1 *h_MassXim;
    TH1 *h_MassXip;
    TH1 *h_MassXiC;
    TH1 *h_MassOmm;
    TH1 *h_MassOmp;
    TH1 *h_MassOmC;
    TH1D *h_MassXim_pt_mult[fNptbins_Xi][fNmultbins_Xi];
    TH1D *h_MassXip_pt_mult[fNptbins_Xi][fNmultbins_Xi];
    TH1D *h_MassOmm_pt_mult[fNptbins_Om][fNmultbins_Om];
    TH1D *h_MassOmp_pt_mult[fNptbins_Om][fNmultbins_Om];
    TH1 *resultParXip_pt_mult[fNptbins_Xi][fNmultbins_Xi];
    TH1 *resultParXim_pt_mult[fNptbins_Xi][fNmultbins_Xi];
    TH1 *resultParOmp_pt_mult[fNptbins_Om][fNmultbins_Om];
    TH1 *resultParOmm_pt_mult[fNptbins_Om][fNmultbins_Om];

    /// +/- combination:
    TH1D *h_MassXiC_pt_mult[fNptbins_Xi][fNmultbins_Xi];
    TH1D *h_MassOmC_pt_mult[fNptbins_Om][fNmultbins_Om];
    TH1 *resultParXiC_pt_mult[fNptbins_Xi][fNmultbins_Xi];
    TH1 *resultParOmC_pt_mult[fNptbins_Om][fNmultbins_Om];

    TH1D *h_MassXim_pt[fNptbins_Xi];
    TH1D *h_MassXip_pt[fNptbins_Xi];
    TH1D *h_MassOmm_pt[fNptbins_Om];
    TH1D *h_MassOmp_pt[fNptbins_Om];
    TH1 *resultParXip_pt[fNptbins_Xi];
    TH1 *resultParXim_pt[fNptbins_Xi];
    TH1 *resultParOmp_pt[fNptbins_Om];
    TH1 *resultParOmm_pt[fNptbins_Om];

    /// +/- combination:
    TH1D *h_MassXiC_pt[fNptbins_Xi];
    TH1D *h_MassOmC_pt[fNptbins_Om];
    TH1 *resultParXiC_pt[fNptbins_Xi];
    TH1 *resultParOmC_pt[fNptbins_Om];

    /// fit functions for double gaus fit:
    TF1 *f_Sig_Xip;
    TF1 *f_Sig_Xim;
    TF1 *f_Sig_XiC;
    TF1 *f_Sig_Xim_pt[fNptbins_Xi];
    TF1 *f_Sig_Xip_pt[fNptbins_Xi];
    TF1 *f_Sig_XiC_pt[fNptbins_Xi];
    TF1 *f_Sig_Xim_pt_mult[fNptbins_Xi][fNmultbins_Xi];
    TF1 *f_Sig_Xip_pt_mult[fNptbins_Xi][fNmultbins_Xi];
    TF1 *f_Sig_XiC_pt_mult[fNptbins_Xi][fNmultbins_Xi];

    TF1 *f_Sig_Omp;
    TF1 *f_Sig_Omm;
    TF1 *f_Sig_OmC;
    TF1 *f_Sig_Omm_pt[fNptbins_Om];
    TF1 *f_Sig_Omp_pt[fNptbins_Om];
    TF1 *f_Sig_OmC_pt[fNptbins_Om];
    TF1 *f_Sig_Omm_pt_mult[fNptbins_Om][fNmultbins_Om];
    TF1 *f_Sig_Omp_pt_mult[fNptbins_Om][fNmultbins_Om];
    TF1 *f_Sig_OmC_pt_mult[fNptbins_Om][fNmultbins_Om];

    // background estimation hist through TSpectrum
    // TH1 *h_bgXim_pt[fNptbins_Xi];
    // TH1 *h_bgXip_pt[fNptbins_Xi];
    // TH1 *h_bgXiC_pt[fNptbins_Xi];
    // TH1 *h_bgXim_pt_mult[fNptbins_Xi][fNmultbins_Xi];
    // TH1 *h_bgXip_pt_mult[fNptbins_Xi][fNmultbins_Xi];
    // TH1 *h_bgXiC_pt_mult[fNptbins_Xi][fNmultbins_Xi];

    // TH1 *h_bgOmm_pt[fNptbins_Om];
    // TH1 *h_bgOmp_pt[fNptbins_Om];
    // TH1 *h_bgOmC_pt[fNptbins_Om];
    // TH1 *h_bgOmm_pt_mult[fNptbins_Om][fNmultbins_Om];
    // TH1 *h_bgOmp_pt_mult[fNptbins_Om][fNmultbins_Om];
    // TH1 *h_bgOmC_pt_mult[fNptbins_Om][fNmultbins_Om];

    // Signal+Bg (peak)
    if (xi)
    {
        h_MassXim = (TH1 *)Xim->Project3D(idAxis);
        h_MassXim->SetNameTitle("h_MassXim", "Invariant Mass #Xi^{-}");
        h_MassXip = (TH1 *)Xip->Project3D(idAxis);
        h_MassXip->SetNameTitle("h_MassXip", "Invariant Mass #Xi^{+}");
        h_MassXiC = (TH1 *)h_MassXim->Clone();
        h_MassXiC->Add(h_MassXip);
        h_MassXiC->SetNameTitle("h_MassXiC", "Invariant Mass #Xi^{+} + #Xi^{-}");

        // when fitting variable bin histograms you need first to scale the histogram using the bin width
        h_MassXim->Scale(1, "width");
        h_MassXip->Scale(1, "width");
        h_MassXiC->Scale(1, "width");
    }
    if (om)
    {
        h_MassOmm = (TH1 *)Omm->Project3D(idAxis);
        h_MassOmm->SetNameTitle("h_MassOmm", "Invariant Mass #Omega^{-}");
        h_MassOmp = (TH1 *)Omp->Project3D(idAxis);
        h_MassOmp->SetNameTitle("h_MassOmp", "Invariant Mass #Omega^{+}");
        h_MassOmC = (TH1 *)h_MassOmm->Clone();
        h_MassOmC->Add(h_MassOmp);
        h_MassOmC->SetNameTitle("h_MassOmC", "Invariant Mass #Omega^{+} + #Omega^{-}");

        // when fitting variable bin histograms you need first to scale the histogram using the bin width
        h_MassOmm->Scale(1, "width");
        h_MassOmp->Scale(1, "width");
        h_MassOmC->Scale(1, "width");
    }

    TH1 *resultParams_Xip_allInt, *resultParams_Xim_allInt, *resultParams_XiC_allInt;
    TH1 *resultParams_Omp_allInt, *resultParams_Omm_allInt, *resultParams_OmC_allInt;

    Double_t fitMinSig_Xi = fMass_Xi - nSigma * sigmaXi;
    Double_t fitMaxSig_Xi = fMass_Xi + nSigma * sigmaXi;
    Double_t fitMinSig_Om = fMass_Om - nSigma * sigmaOm;
    Double_t fitMaxSig_Om = fMass_Om + nSigma * sigmaOm;

    if (xi)
    {
        /// pt+mult integrated case
        binsMC[0] = 1;
        binsMC[1] = fNptbins_Xi;
        binsMC[2] = 1;
        binsMC[3] = fNmultbins_Xi;
        Info("FitCascades: Xi_allInt", "Starting XI pt+mult integrated analysis...");

        if (isDoubleGausPol2)
        {
            // Info("FitCascades: Xi_allInt", "Creating FF...");
            f_Sig_Xip = new TF1("DblGausPol2_Xip", DoubleGausPol2, fitMinSig_Xi, fitMaxSig_Xi, 8);
            f_Sig_Xim = new TF1("DblGausPol2_Xim", DoubleGausPol2, fitMinSig_Xi, fitMaxSig_Xi, 8);
            f_Sig_XiC = new TF1("DblGausPol2_XiC", DoubleGausPol2, fitMinSig_Xi, fitMaxSig_Xi, 8);
            f_Sig_Xip->SetNpx(1000);
            f_Sig_Xim->SetNpx(1000);
            f_Sig_XiC->SetNpx(1000);

            // Info("FitCascades: Xi_allInt", "FF created. Setting FF params..");
            f_Sig_Xip = SetFitParametersDG(f_Sig_Xip, nullptr, h_MassXip, fMass_Xi, sigmaXi, kFALSE);
            f_Sig_Xim = SetFitParametersDG(f_Sig_Xim, nullptr, h_MassXim, fMass_Xi, sigmaXi, kFALSE);
            f_Sig_XiC = SetFitParametersDG(f_Sig_XiC, nullptr, h_MassXiC, fMass_Xi, sigmaXi, kFALSE);
            // Info("FitCascades: Xi_allInt", "DG Fit params set.");
        }
        if (isGausPol2)
        {
            f_Sig_Xip = new TF1("GausPol2_Xip", GausPol2, fitMinSig_Xi, fitMaxSig_Xi, 6);
            f_Sig_Xim = new TF1("GausPol2_Xim", GausPol2, fitMinSig_Xi, fitMaxSig_Xi, 6);
            f_Sig_XiC = new TF1("GausPol2_XiC", GausPol2, fitMinSig_Xi, fitMaxSig_Xi, 6);
            f_Sig_Xip->SetNpx(1000);
            f_Sig_Xim->SetNpx(1000);
            f_Sig_XiC->SetNpx(1000);

            f_Sig_Xip = SetFitParametersGaus(f_Sig_Xip, nullptr, h_MassXip, fMass_Xi, sigmaXi, kFALSE);
            f_Sig_Xim = SetFitParametersGaus(f_Sig_Xim, nullptr, h_MassXim, fMass_Xi, sigmaXi, kFALSE);
            f_Sig_XiC = SetFitParametersGaus(f_Sig_XiC, nullptr, h_MassXiC, fMass_Xi, sigmaXi, kFALSE);
        }
        resultParams_Xip_allInt = FitResults(nullptr, f_Sig_Xip, h_MassXip, fitMinSig_Xi, fitMaxSig_Xi, binsMC, MCGenXip, fisMC, isGausPol2, fitOptSig, fitOptBg);
        resultParams_Xim_allInt = FitResults(nullptr, f_Sig_Xim, h_MassXim, fitMinSig_Xi, fitMaxSig_Xi, binsMC, MCGenXim, fisMC, isGausPol2, fitOptSig, fitOptBg);
        resultParams_XiC_allInt = FitResults(nullptr, f_Sig_XiC, h_MassXiC, fitMinSig_Xi, fitMaxSig_Xi, binsMC, MCGenXiC, fisMC, isGausPol2, fitOptSig, fitOptBg);

        outputFile->cd("_allInt");
        h_MassXip->Write();
        h_MassXim->Write();
        h_MassXiC->Write();
        resultParams_Xip_allInt->SetNameTitle("resultParams_Xip_allInt", "Result Parameters #Xi^{+}");
        resultParams_Xim_allInt->SetNameTitle("resultParams_Xim_allInt", "Result Parameters #Xi^{-}");
        resultParams_XiC_allInt->SetNameTitle("resultParams_XiC_allInt", "Result Parameters #Xi^{+} + #Xi^{-}");
        resultParams_Xip_allInt->Write();
        resultParams_Xim_allInt->Write();
        resultParams_XiC_allInt->Write();

        Info("FitCascades: Xi_allInt", "Finished writing allInt hists...");

        for (Int_t ptBinXi = 0; ptBinXi < fNptbins_Xi; ptBinXi++)
        {
            Info("FitCascades: Xi_multInt", "Loop index: pT bin = %d ...", ptBinXi);

            if (fisMC)
            {
                // plotting mult integrated efficiency as a function of pT
                binsMC[0] = ptBinXi + 1;
                binsMC[1] = ptBinXi + 1;
                binsMC[2] = 1;
                binsMC[3] = fNmultbins_Xi;
            }
            h_MassXim_pt[ptBinXi] = (TH1D *)Xim->ProjectionY(TString::Format(("h_MassXim_pt[%d]"), ptBinXi), ptBinXi + 1, ptBinXi + 1, 1, fNmultbins_Xi, "e");
            h_MassXim_pt[ptBinXi]->SetTitle(TString::Format(("Invariant Mass #Xi^{-}: p_{T}<%.1f,%.1f>, Mult<%.1f,%.1f>"), fPtbins_Xi[ptBinXi], fPtbins_Xi[ptBinXi + 1], fMultbins_Xi[0], fMultbins_Xi[fNmultbins_Xi]));
            /// when fitting variable bin histograms you need first to scale the histogram using the bin width
            // h_MassXim_pt[ptBinXi]->Scale(1, "width");
            h_MassXip_pt[ptBinXi] = (TH1D *)Xip->ProjectionY(TString::Format(("h_MassXip_pt[%d]"), ptBinXi), ptBinXi + 1, ptBinXi + 1, 1, fNmultbins_Xi, "e");
            h_MassXip_pt[ptBinXi]->SetTitle(TString::Format(("Invariant Mass #Xi^{+}: p_{T}<%.1f,%.1f>, Mult<%.1f,%.1f>"), fPtbins_Xi[ptBinXi], fPtbins_Xi[ptBinXi + 1], fMultbins_Xi[0], fMultbins_Xi[fNmultbins_Xi]));
            // h_MassXip_pt[ptBinXi]->Scale(1, "width"); // when fitting variable bin histograms you need first to scale the histogram using the bin width
            h_MassXiC_pt[ptBinXi] = (TH1D *)h_MassXip_pt[ptBinXi]->Clone(TString::Format(("h_MassXiC_pt[%d]"), ptBinXi));
            h_MassXiC_pt[ptBinXi]->Add(h_MassXim_pt[ptBinXi]);
            h_MassXiC_pt[ptBinXi]->SetTitle(TString::Format(("Invariant Mass #Xi: p_{T}<%.1f,%.1f>, Mult<%.1f,%.1f>"), fPtbins_Xi[ptBinXi], fPtbins_Xi[ptBinXi + 1], fMultbins_Xi[0], fMultbins_Xi[fNmultbins_Xi]));
            // h_MassXiC_pt[ptBinXi]->Scale(1, "width"); // when fitting variable bin histograms you need first to scale the histogram using the bin width

            if (isDoubleGausPol2)
            {
                f_Sig_Xip_pt[ptBinXi] = new TF1(TString::Format("DblGausPol2_Xip[%d]", ptBinXi), DoubleGausPol2, fitMinSig_Xi, fitMaxSig_Xi, 8);
                f_Sig_Xim_pt[ptBinXi] = new TF1(TString::Format("DblGausPol2_Xim[%d]", ptBinXi), DoubleGausPol2, fitMinSig_Xi, fitMaxSig_Xi, 8);
                f_Sig_XiC_pt[ptBinXi] = new TF1(TString::Format("DblGausPol2_XiC[%d]", ptBinXi), DoubleGausPol2, fitMinSig_Xi, fitMaxSig_Xi, 8);
                f_Sig_Xip_pt[ptBinXi]->SetNpx(1000);
                f_Sig_Xim_pt[ptBinXi]->SetNpx(1000);
                f_Sig_XiC_pt[ptBinXi]->SetNpx(1000);
                if (ptBinXi < 1)
                {
                    f_Sig_Xip_pt[ptBinXi] = SetFitParametersDG(f_Sig_Xip_pt[ptBinXi], f_Sig_Xip->GetParameters(), h_MassXip_pt[ptBinXi], fMass_Xi, sigmaXi, kTRUE);
                    f_Sig_Xim_pt[ptBinXi] = SetFitParametersDG(f_Sig_Xim_pt[ptBinXi], f_Sig_Xim->GetParameters(), h_MassXim_pt[ptBinXi], fMass_Xi, sigmaXi, kTRUE);
                    f_Sig_XiC_pt[ptBinXi] = SetFitParametersDG(f_Sig_XiC_pt[ptBinXi], f_Sig_XiC->GetParameters(), h_MassXiC_pt[ptBinXi], fMass_Xi, sigmaXi, kTRUE);
                }
                else
                {

                    f_Sig_Xip_pt[ptBinXi] = SetFitParametersDG(f_Sig_Xip_pt[ptBinXi], f_Sig_Xip_pt[ptBinXi - 1]->GetParameters(), h_MassXip_pt[ptBinXi], fMass_Xi, sigmaXi, kTRUE);
                    f_Sig_Xim_pt[ptBinXi] = SetFitParametersDG(f_Sig_Xim_pt[ptBinXi], f_Sig_Xim_pt[ptBinXi - 1]->GetParameters(), h_MassXim_pt[ptBinXi], fMass_Xi, sigmaXi, kTRUE);
                    f_Sig_XiC_pt[ptBinXi] = SetFitParametersDG(f_Sig_XiC_pt[ptBinXi], f_Sig_XiC_pt[ptBinXi - 1]->GetParameters(), h_MassXiC_pt[ptBinXi], fMass_Xi, sigmaXi, kTRUE);
                }
            }
            if (isGausPol2)
            {
                f_Sig_Xip_pt[ptBinXi] = new TF1(TString::Format("GausPol2_Xip[%d]", ptBinXi), GausPol2, fitMinSig_Xi, fitMaxSig_Xi, 6);
                f_Sig_Xim_pt[ptBinXi] = new TF1(TString::Format("GausPol2_Xim[%d]", ptBinXi), GausPol2, fitMinSig_Xi, fitMaxSig_Xi, 6);
                f_Sig_XiC_pt[ptBinXi] = new TF1(TString::Format("GausPol2_XiC[%d]", ptBinXi), GausPol2, fitMinSig_Xi, fitMaxSig_Xi, 6);
                f_Sig_Xip_pt[ptBinXi]->SetNpx(1000);
                f_Sig_Xim_pt[ptBinXi]->SetNpx(1000);
                f_Sig_XiC_pt[ptBinXi]->SetNpx(1000);
                if (ptBinXi < 1)
                {
                    f_Sig_Xip_pt[ptBinXi] = SetFitParametersGaus(f_Sig_Xip_pt[ptBinXi], f_Sig_Xip->GetParameters(), h_MassXip_pt[ptBinXi], fMass_Xi, sigmaXi, kTRUE);
                    f_Sig_Xim_pt[ptBinXi] = SetFitParametersGaus(f_Sig_Xim_pt[ptBinXi], f_Sig_Xim->GetParameters(), h_MassXim_pt[ptBinXi], fMass_Xi, sigmaXi, kTRUE);
                    f_Sig_XiC_pt[ptBinXi] = SetFitParametersGaus(f_Sig_XiC_pt[ptBinXi], f_Sig_XiC->GetParameters(), h_MassXiC_pt[ptBinXi], fMass_Xi, sigmaXi, kTRUE);
                }
                else
                {
                    f_Sig_Xip_pt[ptBinXi] = SetFitParametersGaus(f_Sig_Xip_pt[ptBinXi], f_Sig_Xip_pt[ptBinXi - 1]->GetParameters(), h_MassXip_pt[ptBinXi], fMass_Xi, sigmaXi, kTRUE);
                    f_Sig_Xim_pt[ptBinXi] = SetFitParametersGaus(f_Sig_Xim_pt[ptBinXi], f_Sig_Xim_pt[ptBinXi - 1]->GetParameters(), h_MassXim_pt[ptBinXi], fMass_Xi, sigmaXi, kTRUE);
                    f_Sig_XiC_pt[ptBinXi] = SetFitParametersGaus(f_Sig_XiC_pt[ptBinXi], f_Sig_XiC_pt[ptBinXi - 1]->GetParameters(), h_MassXiC_pt[ptBinXi], fMass_Xi, sigmaXi, kTRUE);
                }
            }

            // h_bgXip_pt[ptBinXi] = GenerateBg(bgTSpectrum, h_MassXip_pt[ptBinXi]);
            resultParXip_pt[ptBinXi] = FitResults(nullptr, f_Sig_Xip_pt[ptBinXi], h_MassXip_pt[ptBinXi], fitMinSig_Xi, fitMaxSig_Xi, binsMC, MCGenXip, fisMC, isGausPol2, fitOptSig, fitOptBg);
            resultParXip_pt[ptBinXi]->SetName(TString::Format(("resultParXip_pt[%d]"), ptBinXi));
            resultParXip_pt[ptBinXi]->SetTitle(TString::Format(("Result Parameters #Xi^{+}: p_{T}<%.1f,%.1f>, Mult<%.1f,%.1f>"), fPtbins_Xi[ptBinXi], fPtbins_Xi[ptBinXi + 1], fMultbins_Xi[0], fMultbins_Xi[fNmultbins_Xi]));

            // h_bgXim_pt[ptBinXi] = GenerateBg(bgTSpectrum, h_MassXim_pt[ptBinXi]);
            resultParXim_pt[ptBinXi] = FitResults(nullptr, f_Sig_Xim_pt[ptBinXi], h_MassXim_pt[ptBinXi], fitMinSig_Xi, fitMaxSig_Xi, binsMC, MCGenXim, fisMC, isGausPol2, fitOptSig, fitOptBg);
            resultParXim_pt[ptBinXi]->SetName(TString::Format(("resultParXim_pt[%d]"), ptBinXi));
            resultParXim_pt[ptBinXi]->SetTitle(TString::Format(("Result Parameters #Xi^{-}: p_{T}<%.1f,%.1f>, Mult<%.1f,%.1f>"), fPtbins_Xi[ptBinXi], fPtbins_Xi[ptBinXi + 1], fMultbins_Xi[0], fMultbins_Xi[fNmultbins_Xi]));

            // h_bgXiC_pt[ptBinXi] = GenerateBg(bgTSpectrum, h_MassXiC_pt[ptBinXi]);
            resultParXiC_pt[ptBinXi] = FitResults(nullptr, f_Sig_XiC_pt[ptBinXi], h_MassXiC_pt[ptBinXi], fitMinSig_Xi, fitMaxSig_Xi, binsMC, MCGenXiC, fisMC, isGausPol2, fitOptSig, fitOptBg);
            resultParXiC_pt[ptBinXi]->SetName(TString::Format(("resultParXiC_pt[%d]"), ptBinXi));
            resultParXiC_pt[ptBinXi]->SetTitle(TString::Format(("Result Parameters #Xi: p_{T}<%.1f,%.1f>, Mult<%.1f,%.1f>"), fPtbins_Xi[ptBinXi], fPtbins_Xi[ptBinXi + 1], fMultbins_Xi[0], fMultbins_Xi[fNmultbins_Xi]));

            outputFile->cd("h_MassXim_pt");
            h_MassXim_pt[ptBinXi]->Write();
            resultParXim_pt[ptBinXi]->Write();

            outputFile->cd("h_MassXip_pt");
            h_MassXip_pt[ptBinXi]->Write();
            resultParXip_pt[ptBinXi]->Write();

            outputFile->cd("h_MassXiC_pt");
            h_MassXiC_pt[ptBinXi]->Write();
            resultParXiC_pt[ptBinXi]->Write();

            Info("FitCascades: Xi_multInt", "Successfully written hists for pT bin = %d ...", ptBinXi);
        }

        /// Starting fully differential case
        for (Int_t ptBinXi = 0; ptBinXi < fNptbins_Xi; ptBinXi++)
        {
            if (fisMC)
            {
                binsMC[0] = ptBinXi + 1;
                binsMC[1] = ptBinXi + 1;
            }

            for (Int_t multBinXi = 0; multBinXi < fNmultbins_Xi; multBinXi++)
            {
                Info("FitCascades: Xi_diff", "Loop index: pT bin = %d , Mult bin = %d ...", ptBinXi, multBinXi);

                if (fisMC) /// fully differential case
                {
                    binsMC[2] = multBinXi + 1;
                    binsMC[3] = multBinXi + 1;
                }

                h_MassXim_pt_mult[ptBinXi][multBinXi] = (TH1D *)Xim->ProjectionY(TString::Format(("h_MassXim_pt_mult[%d][%d]"), ptBinXi, multBinXi), ptBinXi + 1, ptBinXi + 1, multBinXi + 1, multBinXi + 1, "e");
                h_MassXim_pt_mult[ptBinXi][multBinXi]->SetTitle(TString::Format(("Invariant Mass #Xi^{-}: p_{T}<%.1f,%.1f>, Mult<%.1f,%.1f>"), fPtbins_Xi[ptBinXi], fPtbins_Xi[ptBinXi + 1], fMultbins_Xi[multBinXi], fMultbins_Xi[multBinXi + 1]));

                h_MassXip_pt_mult[ptBinXi][multBinXi] = (TH1D *)Xip->ProjectionY(TString::Format(("h_MassXip_pt_mult[%d][%d]"), ptBinXi, multBinXi), ptBinXi + 1, ptBinXi + 1, multBinXi + 1, multBinXi + 1, "e");
                h_MassXip_pt_mult[ptBinXi][multBinXi]->SetTitle(TString::Format(("Invariant Mass #Xi^{+}: p_{T}<%.1f,%.1f>, Mult<%.1f,%.1f>"), fPtbins_Xi[ptBinXi], fPtbins_Xi[ptBinXi + 1], fMultbins_Xi[multBinXi], fMultbins_Xi[multBinXi + 1]));

                h_MassXiC_pt_mult[ptBinXi][multBinXi] = (TH1D *)h_MassXip_pt_mult[ptBinXi][multBinXi]->Clone(TString::Format(("h_MassXiC_pt_mult[%d][%d]"), ptBinXi, multBinXi));
                h_MassXiC_pt_mult[ptBinXi][multBinXi]->Add(h_MassXim_pt_mult[ptBinXi][multBinXi]);
                h_MassXiC_pt_mult[ptBinXi][multBinXi]->SetTitle(TString::Format(("Invariant Mass #Xi: p_{T}<%.1f,%.1f>, Mult<%.1f,%.1f>"), fPtbins_Xi[ptBinXi], fPtbins_Xi[ptBinXi + 1], fMultbins_Xi[multBinXi], fMultbins_Xi[multBinXi + 1]));

                if (isDoubleGausPol2)
                {
                    f_Sig_Xip_pt_mult[ptBinXi][multBinXi] = new TF1(TString::Format("DblGausPol2_Xip[%d][%d]", ptBinXi, multBinXi), DoubleGausPol2, fitMinSig_Xi, fitMaxSig_Xi, 8);
                    f_Sig_Xim_pt_mult[ptBinXi][multBinXi] = new TF1(TString::Format("DblGausPol2_Xim[%d][%d]", ptBinXi, multBinXi), DoubleGausPol2, fitMinSig_Xi, fitMaxSig_Xi, 8);
                    f_Sig_XiC_pt_mult[ptBinXi][multBinXi] = new TF1(TString::Format("DblGausPol2_XiC[%d][%d]", ptBinXi, multBinXi), DoubleGausPol2, fitMinSig_Xi, fitMaxSig_Xi, 8);
                    f_Sig_Xip_pt_mult[ptBinXi][multBinXi]->SetNpx(1000);
                    f_Sig_Xim_pt_mult[ptBinXi][multBinXi]->SetNpx(1000);
                    f_Sig_XiC_pt_mult[ptBinXi][multBinXi]->SetNpx(1000);

                    // if (ptBinXi == 0)
                    // {
                    f_Sig_Xip_pt_mult[ptBinXi][multBinXi] = SetFitParametersDG(f_Sig_Xip_pt_mult[ptBinXi][multBinXi], f_Sig_Xip_pt[ptBinXi]->GetParameters(), h_MassXip_pt_mult[ptBinXi][multBinXi], fMass_Xi, resultParXip_pt[ptBinXi]->GetBinContent(4), kTRUE); // ptBin needs to be >1 in this case to ensure correct flow in SetFitParametersDG to not reset fn params
                    f_Sig_Xim_pt_mult[ptBinXi][multBinXi] = SetFitParametersDG(f_Sig_Xim_pt_mult[ptBinXi][multBinXi], f_Sig_Xim_pt[ptBinXi]->GetParameters(), h_MassXim_pt_mult[ptBinXi][multBinXi], fMass_Xi, resultParXim_pt[ptBinXi]->GetBinContent(4), kTRUE);
                    f_Sig_XiC_pt_mult[ptBinXi][multBinXi] = SetFitParametersDG(f_Sig_XiC_pt_mult[ptBinXi][multBinXi], f_Sig_XiC_pt[ptBinXi]->GetParameters(), h_MassXiC_pt_mult[ptBinXi][multBinXi], fMass_Xi, resultParXiC_pt[ptBinXi]->GetBinContent(4), kTRUE);
                    // }
                    // else
                    // {
                    //     f_Sig_Xip_pt_mult[ptBinXi][multBinXi] = SetFitParametersDG(f_Sig_Xip_pt_mult[ptBinXi][multBinXi], f_Sig_Xip_pt_mult[ptBinXi - 1][multBinXi]->GetParameters(), h_MassXip_pt_mult[ptBinXi][multBinXi], fMass_Xi, sigmaXi, ptBinXi + 1, kTRUE); // multInt needs to be true for correct flow in SetFitParametersDG -> TODO: fix it later in SetFitParametersDG
                    //     f_Sig_Xim_pt_mult[ptBinXi][multBinXi] = SetFitParametersDG(f_Sig_Xim_pt_mult[ptBinXi][multBinXi], f_Sig_Xim_pt_mult[ptBinXi - 1][multBinXi]->GetParameters(), h_MassXim_pt_mult[ptBinXi][multBinXi], fMass_Xi, sigmaXi, ptBinXi + 1, kTRUE);
                    //     f_Sig_XiC_pt_mult[ptBinXi][multBinXi] = SetFitParametersDG(f_Sig_XiC_pt_mult[ptBinXi][multBinXi], f_Sig_XiC_pt_mult[ptBinXi - 1][multBinXi]->GetParameters(), h_MassXiC_pt_mult[ptBinXi][multBinXi], fMass_Xi, sigmaXi, ptBinXi + 1, kTRUE);
                    // }
                }

                if (isGausPol2)
                {
                    f_Sig_Xip_pt_mult[ptBinXi][multBinXi] = new TF1(TString::Format("GausPol2_Xip[%d][%d]", ptBinXi, multBinXi), GausPol2, fitMinSig_Xi, fitMaxSig_Xi, 6);
                    f_Sig_Xim_pt_mult[ptBinXi][multBinXi] = new TF1(TString::Format("GausPol2_Xim[%d][%d]", ptBinXi, multBinXi), GausPol2, fitMinSig_Xi, fitMaxSig_Xi, 6);
                    f_Sig_XiC_pt_mult[ptBinXi][multBinXi] = new TF1(TString::Format("GausPol2_XiC[%d][%d]", ptBinXi, multBinXi), GausPol2, fitMinSig_Xi, fitMaxSig_Xi, 6);
                    f_Sig_Xip_pt_mult[ptBinXi][multBinXi]->SetNpx(1000);
                    f_Sig_Xim_pt_mult[ptBinXi][multBinXi]->SetNpx(1000);
                    f_Sig_XiC_pt_mult[ptBinXi][multBinXi]->SetNpx(1000);

                    // if (ptBinXi == 0)
                    // {
                    f_Sig_Xip_pt_mult[ptBinXi][multBinXi] = SetFitParametersGaus(f_Sig_Xip_pt_mult[ptBinXi][multBinXi], f_Sig_Xip_pt[ptBinXi]->GetParameters(), h_MassXip_pt_mult[ptBinXi][multBinXi], fMass_Xi, sigmaXi, kTRUE); // ptBin needs to be >1 in this case to ensure correct flow in SetFitParametersGaus to not reset fn params
                    f_Sig_Xim_pt_mult[ptBinXi][multBinXi] = SetFitParametersGaus(f_Sig_Xim_pt_mult[ptBinXi][multBinXi], f_Sig_Xim_pt[ptBinXi]->GetParameters(), h_MassXim_pt_mult[ptBinXi][multBinXi], fMass_Xi, sigmaXi, kTRUE);
                    f_Sig_XiC_pt_mult[ptBinXi][multBinXi] = SetFitParametersGaus(f_Sig_XiC_pt_mult[ptBinXi][multBinXi], f_Sig_XiC_pt[ptBinXi]->GetParameters(), h_MassXiC_pt_mult[ptBinXi][multBinXi], fMass_Xi, sigmaXi, kTRUE);
                    // }
                    // else
                    // {
                    // f_Sig_Xip_pt_mult[ptBinXi][multBinXi] = SetFitParametersGaus(f_Sig_Xip_pt_mult[ptBinXi][multBinXi], f_Sig_Xip_pt_mult[ptBinXi - 1][multBinXi]->GetParameters(), h_MassXip_pt_mult[ptBinXi][multBinXi], fMass_Xi, sigmaXi, kTRUE); // multInt needs to be true for correct flow in SetFitParametersGaus -> TODO: fix it later in SetFitParametersDG
                    // f_Sig_Xim_pt_mult[ptBinXi][multBinXi] = SetFitParametersGaus(f_Sig_Xim_pt_mult[ptBinXi][multBinXi], f_Sig_Xim_pt_mult[ptBinXi - 1][multBinXi]->GetParameters(), h_MassXim_pt_mult[ptBinXi][multBinXi], fMass_Xi, sigmaXi, kTRUE);
                    // f_Sig_XiC_pt_mult[ptBinXi][multBinXi] = SetFitParametersGaus(f_Sig_XiC_pt_mult[ptBinXi][multBinXi], f_Sig_XiC_pt_mult[ptBinXi - 1][multBinXi]->GetParameters(), h_MassXiC_pt_mult[ptBinXi][multBinXi], fMass_Xi, sigmaXi, kTRUE);
                    // }
                }

                // h_bgXip_pt_mult[ptBinXi][multBinXi] = GenerateBg(bgTSpectrum, h_MassXip_pt_mult[ptBinXi][multBinXi]);
                // resultParXip_pt_mult[ptBinXi][multBinXi] = FitResults(nullptr, f_Sig_Xip_pt[ptBinXi], h_MassXip_pt_mult[ptBinXi][multBinXi], fitMinSig_Xi, fitMaxSig_Xi, &par_allintP[0], &par_allint_errorsP[0], binsMC, MCGenXip, fisMC, kFALSE);
                resultParXip_pt_mult[ptBinXi][multBinXi] = FitResults(nullptr, f_Sig_Xip_pt_mult[ptBinXi][multBinXi], h_MassXip_pt_mult[ptBinXi][multBinXi], fitMinSig_Xi, fitMaxSig_Xi, binsMC, MCGenXip, fisMC, isGausPol2, fitOptSig, fitOptBg);
                resultParXip_pt_mult[ptBinXi][multBinXi]->SetName(TString::Format(("resultParXip_pt_mult[%d][%d]"), ptBinXi, multBinXi));
                resultParXip_pt_mult[ptBinXi][multBinXi]->SetTitle(TString::Format(("Result Parameters #Xi^{+}: p_{T}<%.1f,%.1f>, Mult<%.1f,%.1f>"), fPtbins_Xi[ptBinXi], fPtbins_Xi[ptBinXi + 1], fMultbins_Xi[multBinXi], fMultbins_Xi[multBinXi + 1]));

                // h_bgXim_pt[ptBinXi] = GenerateBg(bgTSpectrum, h_MassXim_pt[ptBinXi]);
                resultParXim_pt_mult[ptBinXi][multBinXi] = FitResults(nullptr, f_Sig_Xim_pt_mult[ptBinXi][multBinXi], h_MassXim_pt_mult[ptBinXi][multBinXi], fitMinSig_Xi, fitMaxSig_Xi, binsMC, MCGenXim, fisMC, isGausPol2, fitOptSig, fitOptBg);
                resultParXim_pt_mult[ptBinXi][multBinXi]->SetName(TString::Format(("resultParXim_pt_mult[%d][%d]"), ptBinXi, multBinXi));
                resultParXim_pt_mult[ptBinXi][multBinXi]->SetTitle(TString::Format(("Result Parameters #Xi^{-}: p_{T}<%.1f,%.1f>, Mult<%.1f,%.1f>"), fPtbins_Xi[ptBinXi], fPtbins_Xi[ptBinXi + 1], fMultbins_Xi[multBinXi], fMultbins_Xi[multBinXi + 1]));

                // h_bgXiC_pt[ptBinXi] = GenerateBg(bgTSpectrum, h_MassXiC_pt_mult[ptBinXi][multBinXi]);
                resultParXiC_pt_mult[ptBinXi][multBinXi] = FitResults(nullptr, f_Sig_XiC_pt_mult[ptBinXi][multBinXi], h_MassXiC_pt_mult[ptBinXi][multBinXi], fitMinSig_Xi, fitMaxSig_Xi, binsMC, MCGenXiC, fisMC, isGausPol2, fitOptSig, fitOptBg);
                resultParXiC_pt_mult[ptBinXi][multBinXi]->SetName(TString::Format(("resultParXiC_pt_mult[%d][%d]"), ptBinXi, multBinXi));
                resultParXiC_pt_mult[ptBinXi][multBinXi]->SetTitle(TString::Format(("Result Parameters #Xi: p_{T}<%.1f,%.1f>, Mult<%.1f,%.1f>"), fPtbins_Xi[ptBinXi], fPtbins_Xi[ptBinXi + 1], fMultbins_Xi[multBinXi], fMultbins_Xi[multBinXi + 1]));

                outputFile->cd("h_MassXim_pt_mult");
                // if (h_bgXim_pt_mult[ptBinXi][multBinXi])
                //     h_bgXim_pt_mult[ptBinXi][multBinXi]->Write();
                h_MassXim_pt_mult[ptBinXi][multBinXi]->Write();
                resultParXim_pt_mult[ptBinXi][multBinXi]->Write();

                outputFile->cd("h_MassXip_pt_mult");
                // if (h_bgXip_pt_mult[ptBinXi][multBinXi])
                //     h_bgXip_pt_mult[ptBinXi][multBinXi]->Write();
                h_MassXip_pt_mult[ptBinXi][multBinXi]->Write();
                resultParXip_pt_mult[ptBinXi][multBinXi]->Write();

                outputFile->cd("h_MassXiC_pt_mult");
                // if (h_bgXiC_pt_mult[ptBinXi][multBinXi])
                //     h_bgXiC_pt_mult[ptBinXi][multBinXi]->Write();
                h_MassXiC_pt_mult[ptBinXi][multBinXi]->Write();
                resultParXiC_pt_mult[ptBinXi][multBinXi]->Write();
                Info("FitCascades: Xi_diff", "Successfully written hists for pT bin = %d , Mult bin = %d ...", ptBinXi, multBinXi);
            }
        }
    }
    if (om)
    {
        /// pt+mult integrated case
        Info("FitCascades: Om_allInt", "Starting Omega pt+mult integrated analysis...");

        binsMC[0] = 1;
        binsMC[1] = fNptbins_Om;
        binsMC[2] = 1;
        binsMC[3] = fNmultbins_Om;
        if (isDoubleGausPol2)
        {
            f_Sig_Omp = new TF1("DblGausPol2_Omp", DoubleGausPol2, fitMinSig_Om, fitMaxSig_Om, 8);
            f_Sig_Omm = new TF1("DblGausPol2_Omm", DoubleGausPol2, fitMinSig_Om, fitMaxSig_Om, 8);
            f_Sig_OmC = new TF1("DblGausPol2_OmC", DoubleGausPol2, fitMinSig_Om, fitMaxSig_Om, 8);
            f_Sig_Omp->SetNpx(1000);
            f_Sig_Omm->SetNpx(1000);
            f_Sig_OmC->SetNpx(1000);

            f_Sig_Omp = SetFitParametersDG(f_Sig_Omp, nullptr, h_MassOmp, fMass_Om, sigmaOm, kFALSE);
            f_Sig_Omm = SetFitParametersDG(f_Sig_Omm, nullptr, h_MassOmm, fMass_Om, sigmaOm, kFALSE);
            f_Sig_OmC = SetFitParametersDG(f_Sig_OmC, nullptr, h_MassOmC, fMass_Om, sigmaOm, kFALSE);
        }
        if (isGausPol2)
        {
            f_Sig_Omp = new TF1("GausPol2_Omp", GausPol2, fitMinSig_Om, fitMaxSig_Om, 6);
            f_Sig_Omm = new TF1("GausPol2_Omm", GausPol2, fitMinSig_Om, fitMaxSig_Om, 6);
            f_Sig_OmC = new TF1("GausPol2_OmC", GausPol2, fitMinSig_Om, fitMaxSig_Om, 6);
            f_Sig_Omp->SetNpx(1000);
            f_Sig_Omm->SetNpx(1000);
            f_Sig_OmC->SetNpx(1000);

            f_Sig_Omp = SetFitParametersGaus(f_Sig_Omp, nullptr, h_MassOmp, fMass_Om, sigmaOm, kFALSE);
            f_Sig_Omm = SetFitParametersGaus(f_Sig_Omm, nullptr, h_MassOmm, fMass_Om, sigmaOm, kFALSE);
            f_Sig_OmC = SetFitParametersGaus(f_Sig_OmC, nullptr, h_MassOmC, fMass_Om, sigmaOm, kFALSE);
        }
        resultParams_Omp_allInt = FitResults(nullptr, f_Sig_Omp, h_MassOmp, fitMinSig_Om, fitMaxSig_Om, binsMC, MCGenOmp, fisMC, isGausPol2, fitOptSig, fitOptBg);
        resultParams_Omm_allInt = FitResults(nullptr, f_Sig_Omm, h_MassOmm, fitMinSig_Om, fitMaxSig_Om, binsMC, MCGenOmm, fisMC, isGausPol2, fitOptSig, fitOptBg);
        resultParams_OmC_allInt = FitResults(nullptr, f_Sig_OmC, h_MassOmC, fitMinSig_Om, fitMaxSig_Om, binsMC, MCGenOmC, fisMC, isGausPol2, fitOptSig, fitOptBg);

        outputFile->cd("_allInt");
        h_MassOmm->Write();
        h_MassOmp->Write();
        h_MassOmC->Write();
        resultParams_Omp_allInt->SetNameTitle("resultParams_Omp_allInt", "Result Parameters #Omega^{+}");
        resultParams_Omm_allInt->SetNameTitle("resultParams_Omm_allInt", "Result Parameters #Omega^{-}");
        resultParams_OmC_allInt->SetNameTitle("resultParams_OmC_allInt", "Result Parameters #Omega^{+} + #Omega^{-}");
        resultParams_Omp_allInt->Write();
        resultParams_Omm_allInt->Write();
        resultParams_OmC_allInt->Write();

        Info("FitCascades: Om_allInt", "Finished writing allInt hists...");

        for (Int_t ptBinOm = 0; ptBinOm < fNptbins_Om; ptBinOm++)
        {
            Info("FitCascades: Om_multInt", "Loop index: pT bin = %d ...", ptBinOm);

            if (fisMC)
            {
                // plotting mult integrated efficiency as a function of pT
                binsMC[0] = ptBinOm + 1;
                binsMC[1] = ptBinOm + 1;
                binsMC[2] = 1;
                binsMC[3] = fNmultbins_Om;
            }

            h_MassOmm_pt[ptBinOm] = (TH1D *)Omm->ProjectionY(TString::Format(("h_MassOmm_pt[%d]"), ptBinOm), ptBinOm + 1, ptBinOm + 1, 1, fNmultbins_Om, "e");
            h_MassOmm_pt[ptBinOm]->SetTitle(TString::Format(("Invariant Mass #Omega^{-}: p_{T}<%.1f,%.1f>, Mult<%.1f,%.1f>"), fPtbins_Om[ptBinOm], fPtbins_Om[ptBinOm + 1], fMultbins_Om[0], fMultbins_Om[fNmultbins_Om]));

            h_MassOmp_pt[ptBinOm] = (TH1D *)Omp->ProjectionY(TString::Format(("h_MassOmp_pt[%d]"), ptBinOm), ptBinOm + 1, ptBinOm + 1, 1, fNmultbins_Om, "e");
            h_MassOmp_pt[ptBinOm]->SetTitle(TString::Format(("Invariant Mass #Omega^{+}: p_{T}<%.1f,%.1f>, Mult<%.1f,%.1f>"), fPtbins_Om[ptBinOm], fPtbins_Om[ptBinOm + 1], fMultbins_Om[0], fMultbins_Om[fNmultbins_Om]));

            h_MassOmC_pt[ptBinOm] = (TH1D *)h_MassOmp_pt[ptBinOm]->Clone(TString::Format(("h_MassOmC_pt[%d]"), ptBinOm));
            h_MassOmC_pt[ptBinOm]->Add(h_MassOmm_pt[ptBinOm]);
            h_MassOmC_pt[ptBinOm]->SetTitle(TString::Format(("Invariant Mass #Omega: p_{T}<%.1f,%.1f>, Mult<%.1f,%.1f>"), fPtbins_Om[ptBinOm], fPtbins_Om[ptBinOm + 1], fMultbins_Om[0], fMultbins_Om[fNmultbins_Om]));

            // omega: merge 2 inv mass bins to 1 -> reducing total no. of bins from 100 to 50, because of low statistics
            h_MassOmm_pt[ptBinOm]->Rebin(2);
            h_MassOmp_pt[ptBinOm]->Rebin(2);
            h_MassOmC_pt[ptBinOm]->Rebin(2);

            if (isDoubleGausPol2)
            {
                f_Sig_Omp_pt[ptBinOm] = new TF1(TString::Format("DblGausPol2_OmP[%d]", ptBinOm), DoubleGausPol2, fitMinSig_Om, fitMaxSig_Om, 8);
                f_Sig_Omm_pt[ptBinOm] = new TF1(TString::Format("DblGausPol2_Omm[%d]", ptBinOm), DoubleGausPol2, fitMinSig_Om, fitMaxSig_Om, 8);
                f_Sig_OmC_pt[ptBinOm] = new TF1(TString::Format("DblGausPol2_OmC[%d]", ptBinOm), DoubleGausPol2, fitMinSig_Om, fitMaxSig_Om, 8);
                f_Sig_Omp_pt[ptBinOm]->SetNpx(1000);
                f_Sig_Omm_pt[ptBinOm]->SetNpx(1000);
                f_Sig_OmC_pt[ptBinOm]->SetNpx(1000);
                if (ptBinOm < 1)
                {
                    f_Sig_Omp_pt[ptBinOm] = SetFitParametersDG(f_Sig_Omp_pt[ptBinOm], f_Sig_Omp->GetParameters(), h_MassOmp_pt[ptBinOm], fMass_Om, sigmaOm, kTRUE);
                    f_Sig_Omm_pt[ptBinOm] = SetFitParametersDG(f_Sig_Omm_pt[ptBinOm], f_Sig_Omm->GetParameters(), h_MassOmm_pt[ptBinOm], fMass_Om, sigmaOm, kTRUE);
                    f_Sig_OmC_pt[ptBinOm] = SetFitParametersDG(f_Sig_OmC_pt[ptBinOm], f_Sig_OmC->GetParameters(), h_MassOmC_pt[ptBinOm], fMass_Om, sigmaOm, kTRUE);
                }
                else
                {

                    f_Sig_Omp_pt[ptBinOm] = SetFitParametersDG(f_Sig_Omp_pt[ptBinOm], f_Sig_Omp_pt[ptBinOm - 1]->GetParameters(), h_MassOmp_pt[ptBinOm], fMass_Om, sigmaOm, kTRUE);
                    f_Sig_Omm_pt[ptBinOm] = SetFitParametersDG(f_Sig_Omm_pt[ptBinOm], f_Sig_Omm_pt[ptBinOm - 1]->GetParameters(), h_MassOmm_pt[ptBinOm], fMass_Om, sigmaOm, kTRUE);
                    f_Sig_OmC_pt[ptBinOm] = SetFitParametersDG(f_Sig_OmC_pt[ptBinOm], f_Sig_OmC_pt[ptBinOm - 1]->GetParameters(), h_MassOmC_pt[ptBinOm], fMass_Om, sigmaOm, kTRUE);
                }
            }
            if (isGausPol2)
            {
                f_Sig_Omp_pt[ptBinOm] = new TF1(TString::Format("GausPol2_OmP[%d]", ptBinOm), GausPol2, fitMinSig_Om, fitMaxSig_Om, 6);
                f_Sig_Omm_pt[ptBinOm] = new TF1(TString::Format("GausPol2_Omm[%d]", ptBinOm), GausPol2, fitMinSig_Om, fitMaxSig_Om, 6);
                f_Sig_OmC_pt[ptBinOm] = new TF1(TString::Format("GausPol2_OmC[%d]", ptBinOm), GausPol2, fitMinSig_Om, fitMaxSig_Om, 6);
                f_Sig_Omp_pt[ptBinOm]->SetNpx(1000);
                f_Sig_Omm_pt[ptBinOm]->SetNpx(1000);
                f_Sig_OmC_pt[ptBinOm]->SetNpx(1000);
                if (ptBinOm < 1)
                {
                    f_Sig_Omp_pt[ptBinOm] = SetFitParametersGaus(f_Sig_Omp_pt[ptBinOm], f_Sig_Omp->GetParameters(), h_MassOmp_pt[ptBinOm], fMass_Om, sigmaOm, kTRUE);
                    f_Sig_Omm_pt[ptBinOm] = SetFitParametersGaus(f_Sig_Omm_pt[ptBinOm], f_Sig_Omm->GetParameters(), h_MassOmm_pt[ptBinOm], fMass_Om, sigmaOm, kTRUE);
                    f_Sig_OmC_pt[ptBinOm] = SetFitParametersGaus(f_Sig_OmC_pt[ptBinOm], f_Sig_OmC->GetParameters(), h_MassOmC_pt[ptBinOm], fMass_Om, sigmaOm, kTRUE);
                }
                else
                {

                    f_Sig_Omp_pt[ptBinOm] = SetFitParametersGaus(f_Sig_Omp_pt[ptBinOm], f_Sig_Omp_pt[ptBinOm - 1]->GetParameters(), h_MassOmp_pt[ptBinOm], fMass_Om, sigmaOm, kTRUE);
                    f_Sig_Omm_pt[ptBinOm] = SetFitParametersGaus(f_Sig_Omm_pt[ptBinOm], f_Sig_Omm_pt[ptBinOm - 1]->GetParameters(), h_MassOmm_pt[ptBinOm], fMass_Om, sigmaOm, kTRUE);
                    f_Sig_OmC_pt[ptBinOm] = SetFitParametersGaus(f_Sig_OmC_pt[ptBinOm], f_Sig_OmC_pt[ptBinOm - 1]->GetParameters(), h_MassOmC_pt[ptBinOm], fMass_Om, sigmaOm, kTRUE);
                }
            }

            // h_bgXip_pt[ptBinXi] = GenerateBg(bgTSpectrum, h_MassXip_pt[ptBinXi]);
            resultParOmm_pt[ptBinOm] = FitResults(nullptr, f_Sig_Omm_pt[ptBinOm], h_MassOmm_pt[ptBinOm], fitMinSig_Om, fitMaxSig_Om, binsMC, MCGenOmm, fisMC, isGausPol2, fitOptSig, fitOptBg);
            resultParOmm_pt[ptBinOm]->SetName(TString::Format(("resultParOmm_pt[%d]"), ptBinOm));
            resultParOmm_pt[ptBinOm]->SetTitle(TString::Format(("Result Parameters #Omega^{-}: p_{T}<%.1f,%.1f>, Mult<%.1f,%.1f>"), fPtbins_Om[ptBinOm], fPtbins_Om[ptBinOm + 1], fMultbins_Om[0], fMultbins_Om[fNmultbins_Om]));

            // h_bgOmm_pt[ptBinOm] = GenerateBg(bgTSpectrum, h_MassOmm_pt[ptBinOm]);
            resultParOmp_pt[ptBinOm] = FitResults(nullptr, f_Sig_Omp_pt[ptBinOm], h_MassOmp_pt[ptBinOm], fitMinSig_Om, fitMaxSig_Om, binsMC, MCGenOmp, fisMC, isGausPol2, fitOptSig, fitOptBg);
            resultParOmp_pt[ptBinOm]->SetName(TString::Format(("resultParOmp_pt[%d]"), ptBinOm));
            resultParOmp_pt[ptBinOm]->SetTitle(TString::Format(("Result Parameters #Omega^{+}: p_{T}<%.1f,%.1f>, Mult<%.1f,%.1f>"), fPtbins_Om[ptBinOm], fPtbins_Om[ptBinOm + 1], fMultbins_Om[0], fMultbins_Om[fNmultbins_Om]));

            // h_bgOmC_pt[ptBinOm] = GenerateBg(bgTSpectrum, h_MassOmC_pt[ptBinOm]);
            resultParOmC_pt[ptBinOm] = FitResults(nullptr, f_Sig_OmC_pt[ptBinOm], h_MassOmC_pt[ptBinOm], fitMinSig_Om, fitMaxSig_Om, binsMC, MCGenOmC, fisMC, isGausPol2, fitOptSig, fitOptBg);
            resultParOmC_pt[ptBinOm]->SetName(TString::Format(("resultParOmC_pt[%d]"), ptBinOm));
            resultParOmC_pt[ptBinOm]->SetTitle(TString::Format(("Result Parameters #Omega: p_{T}<%.1f,%.1f>, Mult<%.1f,%.1f>"), fPtbins_Om[ptBinOm], fPtbins_Om[ptBinOm + 1], fMultbins_Om[0], fMultbins_Om[fNmultbins_Om]));

            outputFile->cd("h_MassOmm_pt");
            h_MassOmm_pt[ptBinOm]->Write();
            // if (h_bgOmm_pt[ptBinOm])
            //     h_bgOmm_pt[ptBinOm]->Write();
            resultParOmm_pt[ptBinOm]->Write();

            outputFile->cd("h_MassOmp_pt");
            h_MassOmp_pt[ptBinOm]->Write();
            // if (h_bgOmp_pt[ptBinOm])
            //     h_bgOmp_pt[ptBinOm]->Write();
            resultParOmp_pt[ptBinOm]->Write();

            outputFile->cd("h_MassOmC_pt");
            h_MassOmC_pt[ptBinOm]->Write();
            // if (h_bgOmC_pt[ptBinOm])
            //     h_bgOmC_pt[ptBinOm]->Write();
            resultParOmC_pt[ptBinOm]->Write();

            Info("FitCascades: Om_multInt", "Successfully written hists for pT bin = %d ...", ptBinOm);
        }

        /// Starting fully differential case
        for (Int_t ptBinOm = 0; ptBinOm < fNptbins_Om; ptBinOm++)
        {
            if (fisMC)
            {
                // plotting mult integrated efficiency as a function of pT
                binsMC[0] = ptBinOm + 1;
                binsMC[1] = ptBinOm + 1;
            }
            for (Int_t multBinOm = 0; multBinOm < fNmultbins_Om; multBinOm++)
            {
                Info("FitCascades: Om_diff", "Loop index: pT bin = %d , Mult bin = %d ...", ptBinOm, multBinOm);

                /// fully differential case
                if (fisMC)
                {
                    binsMC[2] = multBinOm + 1;
                    binsMC[3] = multBinOm + 1;
                }

                h_MassOmm_pt_mult[ptBinOm][multBinOm] = (TH1D *)Omm->ProjectionY(TString::Format(("h_MassOmm_pt_mult[%d][%d]"), ptBinOm, multBinOm), ptBinOm + 1, ptBinOm + 1, multBinOm + 1, multBinOm + 1, "e");
                h_MassOmm_pt_mult[ptBinOm][multBinOm]->SetTitle(TString::Format(("Invariant Mass #Omega^{-}: p_{T}<%.1f,%.1f>, Mult<%.1f,%.1f>"), fPtbins_Om[ptBinOm], fPtbins_Om[ptBinOm + 1], fMultbins_Om[multBinOm], fMultbins_Om[multBinOm + 1]));

                h_MassOmp_pt_mult[ptBinOm][multBinOm] = (TH1D *)Omp->ProjectionY(TString::Format(("h_MassOmp_pt_mult[%d][%d]"), ptBinOm, multBinOm), ptBinOm + 1, ptBinOm + 1, multBinOm + 1, multBinOm + 1, "e");
                h_MassOmp_pt_mult[ptBinOm][multBinOm]->SetTitle(TString::Format(("Invariant Mass #Omega^{+}: p_{T}<%.1f,%.1f>, Mult<%.1f,%.1f>"), fPtbins_Om[ptBinOm], fPtbins_Om[ptBinOm + 1], fMultbins_Om[multBinOm], fMultbins_Om[multBinOm + 1]));

                h_MassOmC_pt_mult[ptBinOm][multBinOm] = (TH1D *)h_MassOmp_pt_mult[ptBinOm][multBinOm]->Clone(TString::Format(("h_MassOmC_pt_mult[%d][%d]"), ptBinOm, multBinOm));
                h_MassOmC_pt_mult[ptBinOm][multBinOm]->Add(h_MassOmm_pt_mult[ptBinOm][multBinOm]);
                h_MassOmC_pt_mult[ptBinOm][multBinOm]->SetTitle(TString::Format(("Invariant Mass #Omega: p_{T}<%.1f,%.1f>, Mult<%.1f,%.1f>"), fPtbins_Om[ptBinOm], fPtbins_Om[ptBinOm + 1], fMultbins_Om[multBinOm], fMultbins_Om[multBinOm + 1]));

                h_MassOmm_pt_mult[ptBinOm][multBinOm]->Rebin(2);
                h_MassOmp_pt_mult[ptBinOm][multBinOm]->Rebin(2);
                h_MassOmC_pt_mult[ptBinOm][multBinOm]->Rebin(2);

                if (isDoubleGausPol2)
                {
                    f_Sig_Omp_pt_mult[ptBinOm][multBinOm] = new TF1(TString::Format("DblGausPol2_Omp[%d][%d]", ptBinOm, multBinOm), DoubleGausPol2, fitMinSig_Om, fitMaxSig_Om, 8);
                    f_Sig_Omm_pt_mult[ptBinOm][multBinOm] = new TF1(TString::Format("DblGausPol2_Omm[%d][%d]", ptBinOm, multBinOm), DoubleGausPol2, fitMinSig_Om, fitMaxSig_Om, 8);
                    f_Sig_OmC_pt_mult[ptBinOm][multBinOm] = new TF1(TString::Format("DblGausPol2_OmC[%d][%d]", ptBinOm, multBinOm), DoubleGausPol2, fitMinSig_Om, fitMaxSig_Om, 8);
                    f_Sig_Omp_pt_mult[ptBinOm][multBinOm]->SetNpx(1000);
                    f_Sig_Omm_pt_mult[ptBinOm][multBinOm]->SetNpx(1000);
                    f_Sig_OmC_pt_mult[ptBinOm][multBinOm]->SetNpx(1000);
                    // if (ptBinOm == 0)
                    // {
                    f_Sig_Omp_pt_mult[ptBinOm][multBinOm] = SetFitParametersDG(f_Sig_Omp_pt_mult[ptBinOm][multBinOm], f_Sig_Omp_pt[ptBinOm]->GetParameters(), h_MassOmp_pt_mult[ptBinOm][multBinOm], fMass_Om, resultParOmp_pt[ptBinOm]->GetBinContent(4), kTRUE);
                    f_Sig_Omm_pt_mult[ptBinOm][multBinOm] = SetFitParametersDG(f_Sig_Omm_pt_mult[ptBinOm][multBinOm], f_Sig_Omm_pt[ptBinOm]->GetParameters(), h_MassOmm_pt_mult[ptBinOm][multBinOm], fMass_Om, resultParOmm_pt[ptBinOm]->GetBinContent(4), kTRUE);
                    f_Sig_OmC_pt_mult[ptBinOm][multBinOm] = SetFitParametersDG(f_Sig_OmC_pt_mult[ptBinOm][multBinOm], f_Sig_OmC_pt[ptBinOm]->GetParameters(), h_MassOmC_pt_mult[ptBinOm][multBinOm], fMass_Om, resultParOmC_pt[ptBinOm]->GetBinContent(4), kTRUE);
                    // }
                    // else
                    // {
                    //     // if (f_Sig_Xip_pt[ptBinXi - 1]->GetParameters())
                    //     // {
                    //     // for (int i = 0; i < 9; i++)
                    //     // {
                    //     //     Printf("%s: Par[%d] = %f", f_Sig_Xip_pt[ptBinXi - 1]->GetName(), i, f_Sig_Xip_pt[ptBinXi - 1]->GetParameter(i));
                    //     // }
                    //     // }
                    //     f_Sig_Omp_pt_mult[ptBinOm][multBinOm] = SetFitParametersDG(f_Sig_Omp_pt_mult[ptBinOm][multBinOm], f_Sig_Omp_pt_mult[ptBinOm - 1][multBinOm]->GetParameters(), h_MassOmp_pt_mult[ptBinOm][multBinOm], fMass_Om, sigmaOm, ptBinOm + 1, kTRUE);
                    //     f_Sig_Omm_pt_mult[ptBinOm][multBinOm] = SetFitParametersDG(f_Sig_Omm_pt_mult[ptBinOm][multBinOm], f_Sig_Omm_pt_mult[ptBinOm - 1][multBinOm]->GetParameters(), h_MassOmm_pt_mult[ptBinOm][multBinOm], fMass_Om, sigmaOm, ptBinOm + 1, kTRUE);
                    //     f_Sig_OmC_pt_mult[ptBinOm][multBinOm] = SetFitParametersDG(f_Sig_OmC_pt_mult[ptBinOm][multBinOm], f_Sig_OmC_pt_mult[ptBinOm - 1][multBinOm]->GetParameters(), h_MassOmp_pt_mult[ptBinOm][multBinOm], fMass_Om, sigmaOm, ptBinOm + 1, kTRUE);
                    // }
                }

                if (isGausPol2)
                {
                    f_Sig_Omp_pt_mult[ptBinOm][multBinOm] = new TF1(TString::Format("GausPol2_Omp[%d][%d]", ptBinOm, multBinOm), GausPol2, fitMinSig_Om, fitMaxSig_Om, 6);
                    f_Sig_Omm_pt_mult[ptBinOm][multBinOm] = new TF1(TString::Format("GausPol2_Omm[%d][%d]", ptBinOm, multBinOm), GausPol2, fitMinSig_Om, fitMaxSig_Om, 6);
                    f_Sig_OmC_pt_mult[ptBinOm][multBinOm] = new TF1(TString::Format("GausPol2_OmC[%d][%d]", ptBinOm, multBinOm), GausPol2, fitMinSig_Om, fitMaxSig_Om, 6);
                    f_Sig_Omp_pt_mult[ptBinOm][multBinOm]->SetNpx(1000);
                    f_Sig_Omm_pt_mult[ptBinOm][multBinOm]->SetNpx(1000);
                    f_Sig_OmC_pt_mult[ptBinOm][multBinOm]->SetNpx(1000);
                    // if (ptBinOm == 0)
                    // {
                    f_Sig_Omp_pt_mult[ptBinOm][multBinOm] = SetFitParametersGaus(f_Sig_Omp_pt_mult[ptBinOm][multBinOm], f_Sig_Omp_pt[ptBinOm]->GetParameters(), h_MassOmp_pt_mult[ptBinOm][multBinOm], fMass_Om, sigmaOm, kTRUE);
                    f_Sig_Omm_pt_mult[ptBinOm][multBinOm] = SetFitParametersGaus(f_Sig_Omm_pt_mult[ptBinOm][multBinOm], f_Sig_Omm_pt[ptBinOm]->GetParameters(), h_MassOmm_pt_mult[ptBinOm][multBinOm], fMass_Om, sigmaOm, kTRUE);
                    f_Sig_OmC_pt_mult[ptBinOm][multBinOm] = SetFitParametersGaus(f_Sig_OmC_pt_mult[ptBinOm][multBinOm], f_Sig_OmC_pt[ptBinOm]->GetParameters(), h_MassOmC_pt_mult[ptBinOm][multBinOm], fMass_Om, sigmaOm, kTRUE);
                    // }
                    // else
                    // {

                    //     f_Sig_Omp_pt_mult[ptBinOm][multBinOm] = SetFitParametersGaus(f_Sig_Omp_pt_mult[ptBinOm][multBinOm], f_Sig_Omp_pt_mult[ptBinOm - 1][multBinOm]->GetParameters(), h_MassOmp_pt_mult[ptBinOm][multBinOm], fMass_Om, sigmaOm, kTRUE);
                    //     f_Sig_Omm_pt_mult[ptBinOm][multBinOm] = SetFitParametersGaus(f_Sig_Omm_pt_mult[ptBinOm][multBinOm], f_Sig_Omm_pt_mult[ptBinOm - 1][multBinOm]->GetParameters(), h_MassOmm_pt_mult[ptBinOm][multBinOm], fMass_Om, sigmaOm, kTRUE);
                    //     f_Sig_OmC_pt_mult[ptBinOm][multBinOm] = SetFitParametersGaus(f_Sig_OmC_pt_mult[ptBinOm][multBinOm], f_Sig_OmC_pt_mult[ptBinOm - 1][multBinOm]->GetParameters(), h_MassOmp_pt_mult[ptBinOm][multBinOm], fMass_Om, sigmaOm, kTRUE);
                    // }
                }

                resultParOmp_pt_mult[ptBinOm][multBinOm] = FitResults(nullptr, f_Sig_Omp_pt_mult[ptBinOm][multBinOm], h_MassOmp_pt_mult[ptBinOm][multBinOm], fitMinSig_Om, fitMaxSig_Om, binsMC, MCGenOmp, fisMC, isGausPol2, fitOptSig, fitOptBg);
                resultParOmp_pt_mult[ptBinOm][multBinOm]->SetName(TString::Format(("resultParOmp_pt_mult[%d][%d]"), ptBinOm, multBinOm));
                resultParOmp_pt_mult[ptBinOm][multBinOm]->SetTitle(TString::Format(("Result Parameters #Omega^{+}: p_{T}<%.1f,%.1f>, Mult<%.1f,%.1f>"), fPtbins_Om[ptBinOm], fPtbins_Om[ptBinOm + 1], fMultbins_Om[multBinOm], fMultbins_Om[multBinOm + 1]));

                // peakFnc_OmM = setFitParameters(peakFnc_OmM, h_MassOmm_pt_mult[ptBinOm][multBinOm], fMass_Om, ptBinOm + 1, kFALSE);
                // h_bgOmm_pt_mult[ptBinOm][multBinOm] = GenerateBg(bgTSpectrum, h_MassOmm_pt_mult[ptBinOm][multBinOm]);
                // resultParOmm_pt_mult[ptBinOm][multBinOm] = FitResults(h_bgOmm_pt_mult[ptBinOm][multBinOm], peakFnc_OmM, h_MassOmm_pt_mult[ptBinOm][multBinOm], fitMinSig_Om, fitMaxSig_Om, &par_allintM[0], &par_allint_errorsM[0], binsMC, MCGenOmm, fisMC, kFALSE);
                resultParOmm_pt_mult[ptBinOm][multBinOm] = FitResults(nullptr, f_Sig_Omm_pt_mult[ptBinOm][multBinOm], h_MassOmm_pt_mult[ptBinOm][multBinOm], fitMinSig_Om, fitMaxSig_Om, binsMC, MCGenOmm, fisMC, isGausPol2, fitOptSig, fitOptBg);
                resultParOmm_pt_mult[ptBinOm][multBinOm]->SetName(TString::Format(("resultParOmm_pt_mult[%d][%d]"), ptBinOm, multBinOm));
                resultParOmm_pt_mult[ptBinOm][multBinOm]->SetTitle(TString::Format(("Result Parameters #Omega^{-}: p_{T}<%.1f,%.1f>, Mult<%.1f,%.1f>"), fPtbins_Om[ptBinOm], fPtbins_Om[ptBinOm + 1], fMultbins_Om[multBinOm], fMultbins_Om[multBinOm + 1]));

                // peakFnc_OmC = setFitParameters(peakFnc_OmC, h_MassOmC_pt_mult[ptBinOm][multBinOm], fMass_Om, ptBinOm + 1, kFALSE);
                // h_bgOmC_pt_mult[ptBinOm][multBinOm] = GenerateBg(bgTSpectrum, h_MassOmC_pt_mult[ptBinOm][multBinOm]);
                // resultParOmC_pt_mult[ptBinOm][multBinOm] = FitResults(h_bgOmC_pt_mult[ptBinOm][multBinOm], peakFnc_OmC, h_MassOmC_pt_mult[ptBinOm][multBinOm], fitMinSig_Om, fitMaxSig_Om, &par_allintC[0], &par_allint_errorsC[0], binsMC, MCGenOmC, fisMC, kFALSE);
                resultParOmC_pt_mult[ptBinOm][multBinOm] = FitResults(nullptr, f_Sig_OmC_pt_mult[ptBinOm][multBinOm], h_MassOmC_pt_mult[ptBinOm][multBinOm], fitMinSig_Om, fitMaxSig_Om, binsMC, MCGenOmC, fisMC, isGausPol2, fitOptSig, fitOptBg);
                resultParOmC_pt_mult[ptBinOm][multBinOm]->SetName(TString::Format(("resultParOmC_pt_mult[%d][%d]"), ptBinOm, multBinOm));
                resultParOmC_pt_mult[ptBinOm][multBinOm]->SetTitle(TString::Format(("Result Parameters #Omega: p_{T}<%.1f,%.1f>, Mult<%.1f,%.1f>"), fPtbins_Om[ptBinOm], fPtbins_Om[ptBinOm + 1], fMultbins_Om[multBinOm], fMultbins_Om[multBinOm + 1]));

                outputFile->cd("h_MassOmm_pt_mult");
                h_MassOmm_pt_mult[ptBinOm][multBinOm]->Write();
                // if (h_bgOmm_pt_mult[ptBinOm][multBinOm])
                //     h_bgOmm_pt_mult[ptBinOm][multBinOm]->Write();
                resultParOmm_pt_mult[ptBinOm][multBinOm]->Write();

                outputFile->cd("h_MassOmp_pt_mult");
                h_MassOmp_pt_mult[ptBinOm][multBinOm]->Write();
                // if (h_bgOmp_pt_mult[ptBinOm][multBinOm])
                //     h_bgOmp_pt_mult[ptBinOm][multBinOm]->Write();
                resultParOmp_pt_mult[ptBinOm][multBinOm]->Write();

                outputFile->cd("h_MassOmC_pt_mult");
                h_MassOmC_pt_mult[ptBinOm][multBinOm]->Write();
                // if (h_bgOmC_pt_mult[ptBinOm][multBinOm])
                //     h_bgOmC_pt_mult[ptBinOm][multBinOm]->Write();
                resultParOmC_pt_mult[ptBinOm][multBinOm]->Write();
                Info("FitCascades: Om_diff", "Successfully written hists for pT bin = %d , Mult bin = %d ...", ptBinOm, multBinOm);
            }
        }
    }

    Printf("FINISHED!!! Output saved to '%s'", outputFile->GetName());
    delete outputFile;
    gErrorIgnoreLevel = kUnset;
    return 0;
}

// TH1 *GenerateBg(bool bgTSpectrum, TH1 *h_source)
// {
//     if (!bgTSpectrum)
//     {
//         return nullptr;
//     }

//     /// testing spectrum
//     Int_t nbinsX = h_source->GetNbinsX();
//     Double_t source[nbinsX];
//     for (int itest = 0; itest < nbinsX; itest++)
//         source[itest] = h_source->GetBinContent(itest + 1);

//     TSpectrum *s = new TSpectrum(1);
//     // int numberIterations[] = {6, 6, 6, 6};
//     // int clipWindow[] = {1, 1};
//     // int backOrder[] = {2, 2, 2, 2};
//     // int backSmoothing[] = {15, 15, 15, 15, 15, 15, 15};
//     for (int i = 0; i < 4; i++)
//     {
//         // s_omC->Background(source, nbinsX, numberIterations[i], clipWindow[i], backOrder[i], kTRUE, backSmoothing[i], kTRUE);
//         // TString hName = TString::Format(("bg_i%d_c%d_o%d_S%d_[%d][%d]"), numberIterations[i], clipWindow[i], backOrder[i], backSmoothing[i], ptBinOm, multBinOm);

//         s->Background(source, nbinsX, 20, TSpectrum::kBackDecreasingWindow, TSpectrum::kBackOrder6, kTRUE, TSpectrum::kBackSmoothing15, kFALSE);

//         // bg[ptBinOm][multBinOm]->GetXaxis()->SetRange(h_MassOmC_pt_mult[ptBinOm][multBinOm]->GetXaxis()->GetXmin(), h_MassOmm_pt_mult[ptBinOm][multBinOm]->GetXaxis()->GetXmax());
//         // Draw the estimated background
//     }
//     // TString hName = TString::Format(("bg_omC[%d][%d]"), ptBinOm, multBinOm);
//     // h_bgOmC_pt_mult[ptBinOm][multBinOm] = new TH1D(hName.Data(), hName.Data(), nbinsX, h_MassOmC_pt_mult[ptBinOm][multBinOm]->GetXaxis()->GetXmin(), h_MassOmC_pt_mult[ptBinOm][multBinOm]->GetXaxis()->GetXmax());
//     TH1D *h_bg = new TH1D(TString::Format(("h_bg%s"), h_source->GetName() + 6), TString::Format(("Estimated Background%s"), h_source->GetTitle() + 14), nbinsX, h_source->GetXaxis()->GetXmin(), h_source->GetXaxis()->GetXmax());
//     h_bg->SetLineColor(kViolet);

//     for (int itest = 0; itest < nbinsX; itest++)
//         h_bg->SetBinContent(itest + 1, source[itest]);

//     return h_bg;
// }

TF1 *SetFitParametersGaus(TF1 *peakFnc, Double_t *previousPtBinFitParams, TH1 *peak, Double_t mass, Double_t sigma, Bool_t usePrevious)
{
    if (usePrevious)
    {
        mass = previousPtBinFitParams[1];
        sigma = previousPtBinFitParams[2];
    }

    TString optBg = "QNFB MULTITHREAD";

    /// gaus fit
    TF1 *bgFnc = new TF1("bgFnc_Pol2", Pol2Exclude, peak->GetXaxis()->GetBinCenter(peak->FindFirstBinAbove(0.)), peak->GetXaxis()->GetBinCenter(peak->FindLastBinAbove(0.)), 6);

    // TF1 *bgFnc = new TF1("bgFnc_Pol2", Pol2Exclude, mass - 15 * sigma, mass + 15 * sigma, 6);
    bgFnc->SetParNames("p[0]", "p[1]", "p[2]", "initMass", "initFourSigma", "rejectBgFlag");
    bgFnc->SetParameter(2, peak->GetMaximum() * 0.1);
    bgFnc->FixParameter(3, mass);        // "mass = mean of fit"
    bgFnc->FixParameter(4, (4 * sigma)); // "4*sigma"
    bgFnc->FixParameter(5, 1);           // "reject" (set 1 or 0)

    for (int iFitBg = 0; iFitBg < 10; iFitBg++)
        peak->Fit(bgFnc, optBg.Data(), "", peak->GetXaxis()->GetBinCenter(peak->FindFirstBinAbove(0.)), peak->GetXaxis()->GetBinCenter(peak->FindLastBinAbove(0.)));

    // peak->Fit(bgFnc, optBg.Data(), "", mass - 15 * sigma, mass + 15 * sigma);
    bgFnc->FixParameter(5, 0); // reject = off (set 1 or 0)
    // peak->GetListOfFunctions()->Add(bgFnc);
    peakFnc->SetParNames("Amplitude", "#mu", "#sigma");
    peakFnc->FixParameter(3, bgFnc->GetParameter(0));
    peakFnc->FixParameter(4, bgFnc->GetParameter(1));
    peakFnc->FixParameter(5, bgFnc->GetParameter(2));

    peakFnc->SetParLimits(0, 0, peak->GetMaximum() * 1.1);
    peakFnc->SetParLimits(1, mass - 2 * sigma, mass + 2 * sigma);
    peakFnc->SetParLimits(2, 0.0009, 0.01);

    if (usePrevious)
    {
        if (previousPtBinFitParams)
        {
            peakFnc->SetParameter(0, previousPtBinFitParams[0]);
            peakFnc->SetParameter(1, previousPtBinFitParams[1]);
            peakFnc->SetParameter(2, previousPtBinFitParams[2]);
        }
        else
        {
            Error("FitCascades: SetFitParametersGaus", "previousPtBinFitParams not passed for %s. Setting defaults...", peak->GetName());
            SetFitParametersGaus(peakFnc, nullptr, peak, mass, sigma, kFALSE);
        }
    }
    else
    {
        peakFnc->SetParameter(0, peak->GetMaximum() * 0.9);
        peakFnc->SetParameter(1, mass);
        peakFnc->SetParameter(2, sigma);
    }

    return peakFnc;
}

TF1 *SetFitParametersDG(TF1 *peakFnc, Double_t *previousPtBinFitParams, TH1 *peak, Double_t mass, Double_t sigma, Bool_t usePrevious)
{
    // if (peak->GetEntries() < 90)
    // {
    //     Error("SetFitParametersDG", "Too low statistics to fit! '%s' entries : %f, skipping ... ", peak->GetName(), peak->GetEntries());
    //     return nullptr;
    // }
    // TString optBg = "QMNR MULTITHREAD";
    TString optBg = "QLNB MULTITHREAD";

    /// double gaus fit
    TF1 *bgFnc = new TF1("bgFnc_Pol2", Pol2Exclude, peak->GetXaxis()->GetBinCenter(peak->FindFirstBinAbove(0.)), peak->GetXaxis()->GetBinCenter(peak->FindLastBinAbove(0.)), 6);

    // TF1 *bgFnc = new TF1("bgFnc_Pol2", Pol2Exclude, mass - 15 * sigma, mass + 15 * sigma, 6);
    bgFnc->SetParNames("p[0]", "p[1]", "p[2]", "initMass", "initFourSigma", "rejectBgFlag");
    bgFnc->SetParameter(2, peak->GetMaximum() * 0.05); 
    bgFnc->FixParameter(3, mass);                      // "mass = mean of fit"
    bgFnc->FixParameter(4, (4 * sigma));               // "4*sigma"

    bgFnc->FixParameter(5, 1); // "reject" (set 1 or 0)

    for (int iFitBg = 0; iFitBg < 3; iFitBg++)
        peak->Fit(bgFnc, optBg.Data(), "", peak->GetXaxis()->GetBinCenter(peak->FindFirstBinAbove(0.)), peak->GetXaxis()->GetBinCenter(peak->FindLastBinAbove(0.)));


    bgFnc->FixParameter(5, 0); // reject = off (set 1 or 0)

    peakFnc->SetParNames("Amp_{1}", "#mu_{1}", "#sigma_{1}", "Amp_{2}", "#sigma_{2}");
    peakFnc->FixParameter(5, bgFnc->GetParameter(2));
    peakFnc->FixParameter(6, bgFnc->GetParameter(1));
    peakFnc->FixParameter(7, bgFnc->GetParameter(0));

    /// set params for mult integrated functions -> these params will be reused for mult differential too
    peakFnc->SetParLimits(0, 0, peak->GetMaximum() * 1.1);
    peakFnc->SetParLimits(1, mass - 2 * sigma, mass + 2 * sigma);
    peakFnc->SetParLimits(2, 0.0009, 0.009);
    peakFnc->SetParLimits(3, 0, peak->GetMaximum() * 0.6);
    peakFnc->SetParLimits(4, 0.0005, 0.025);

    peakFnc->SetParameter(0, peak->GetMaximum() * 0.9);
    peakFnc->SetParameter(3, 0);

    if (usePrevious)
    {

        if (previousPtBinFitParams)
        {
            // peakFnc->SetParameter(0, previousPtBinFitParams[0]);
            peakFnc->SetParameter(1, previousPtBinFitParams[1]);
            peakFnc->SetParameter(2, previousPtBinFitParams[2]);
            // peakFnc->SetParameter(3, previousPtBinFitParams[3]);
            peakFnc->SetParameter(4, previousPtBinFitParams[4]);
            // peakFnc->SetParameter(5, previousPtBinFitParams[5]);
        }
        else
        {
            Error("SetFitParametersDG", "previousPtBinFitParams not passed for %s. Setting defaults ...", peak->GetName());
            SetFitParametersDG(peakFnc, nullptr, peak, mass, sigma, kFALSE);
        }
    }
    else
    {
        peakFnc->SetParameter(1, mass);
        peakFnc->SetParameter(2, sigma);
        peakFnc->SetParameter(4, 0.002);
    }
    peakFnc->Update();

    Info("SetFitParametersDG", "%s: DG Fit params set.", peak->GetName());

    return peakFnc;
}

TH1 *FitResults(TH1 *h_bg, TF1 *peakFnc, TH1 *peak, Double_t fitMinSig, Double_t fitMaxSig, Int_t binsMC[], TH2 *MCgen, bool fisMC, Bool_t isGausPol2, TString fitOptSig, TString fitOptBg)
{
    TH1 *resultParams;

    // if (!peakFnc)
    // {
    //     Error("FitResults", "Cannot fit: function is null, check histogram statistics! '%s' entries : %f", peak->GetName(), peak->GetEntries());
    //     resultParams = new TH1D("emptyParams", "Empty because low statistics!", 1, 0, 1);
    //     return resultParams;
    // }

    Double_t peakFitMass = 0.;
    Double_t peakFitSigma = 0.;
    Double_t peakFitMass_err, peakFitSigma_err = 0.;

    Double_t bgReject = 0.;

    /// use different fit options:
    // TString opt = "QLINES+ MULTITHREAD"; /// best for DG: QLINESR+ MULTITHREAD
    // TString opt = "QMINLSR";

    // TString optBg = "QMNLSR";
    // TString optBg = "QNSFB+ MULTITHREAD"; /// best for DG bkg: QLINSRFB+ MULTITHREAD
    // TString optBg = "QLINSB+ MULTITHREAD"; /// best for DG bkg: QLINSRFB+ MULTITHREAD
    //     TString initFitOptSig = fitOptSig;
    // initFitOptSig.ReplaceAll("E", "");
    //     TString initFitOptBg = fitOptBg;
    // try
    // {
    // Info("FitResults", "%s: Trying out initial peak fit 5 times:", peak->GetName());

    for (int i = 0; i < 5; i++)
        peak->Fit(peakFnc, "QLIN MULTITHREAD", "", fitMinSig, fitMaxSig);

    // Info("FitResults", "%s: Done. Trying out ACTUAL peak fit with options: %s", peak->GetName(), fitOptSig.Data());

    TFitResultPtr fFitResult;
    fFitResult = peak->Fit(peakFnc, fitOptSig.Data(), "", fitMinSig, fitMaxSig);

    if (!fFitResult->IsValid())
    {
        if (fitOptSig.Contains("E "))
        {
            Error("FitResults", "%s: Peak fit not converged. Fit status: %d, Fit Options: '%s'. Trying to fit again without option \"E\" ...", peak->GetName(), fFitResult->Status(), fitOptSig.Data());
            fitOptSig.ReplaceAll("E ", " ");
            fFitResult = peak->Fit(peakFnc, fitOptSig.Data(), "", fitMinSig, fitMaxSig);

            if (!fFitResult->IsValid())
            {
                if (fitOptSig.Contains("QI"))
                {
                    Error("FitResults", "%s: Peak fit STILL not converged. Fit status: %d, Fit Options: '%s'. Trying again a non-quiet fit without option \"I\" : ", peak->GetName(), fFitResult->Status(), fitOptSig.Data());
                    fitOptSig.ReplaceAll("QI", "");
                    fFitResult = peak->Fit(peakFnc, fitOptSig.Data(), "", fitMinSig, fitMaxSig);
                    // fFitResult = peak->Fit(peakFnc, "LNS+ MULTITHREAD", "", fitMinSig, fitMaxSig);

                    if (!fFitResult->IsValid())
                    {
                        Error("FitResults", "\e[1;31m%s\e[0m: Peak fit STILL \e[1;31mNOT\e[0m converged. Fit status: %d, Fit Options: '%s'. Avg mass = %.3f+/-%.3f, sigma = %.3f+/-%.3f. \e[1;31mGiving up\e[0m ...", peak->GetName(), fFitResult->Status(), fitOptSig.Data(), peakFitMass, peakFitMass_err, peakFitSigma, peakFitSigma_err);
                        peakFnc->Print("V");
                    }
                }

                else
                {
                    Error("FitResults", "\e[1;31m%s\e[0m: Peak fit STILL \e[1;31mNOT\e[0m converged. Fit status: %d, Fit Options: '%s'. Avg mass = %.3f+/-%.3f, sigma = %.3f+/-%.3f. \e[1;31mGiving up\e[0m ...", peak->GetName(), fFitResult->Status(), fitOptSig.Data(), peakFitMass, peakFitMass_err, peakFitSigma, peakFitSigma_err);
                    peakFnc->Print("V");
                }
            }
        }
        else
        {
            Error("FitResults", "\e[1;31m%s\e[0m: Peak fit \e[1;31mNOT\e[0m converged. Fit status: %d, Fit Options: '%s'. Avg mass = %.3f+/-%.3f, sigma = %.3f+/-%.3f. \e[1;31mGiving up\e[0m ...", peak->GetName(), fFitResult->Status(), fitOptSig.Data(), peakFitMass, peakFitMass_err, peakFitSigma, peakFitSigma_err);
            // if (gErrorIgnoreLevel < kFatal)
            peakFnc->Print("V");
        }
    }

    Info("FitResults", "%s: Signal fit done. Adding peakFunc '%s' to hist list", peak->GetName(), peakFnc->GetName());

    peak->GetListOfFunctions()->Add(peakFnc);

    if (isGausPol2)
    {
        peakFitMass = peakFnc->GetParameter(1); // getting peak (mass) from peakFnc to be used for excluding bg
        peakFitMass_err = peakFnc->GetParError(1);
        peakFitSigma = TMath::Abs(peakFnc->GetParameter(2));
        peakFitSigma_err = peakFnc->GetParError(2);
    }
    else
        GetMeanSigmaDG(peakFnc, fFitResult, peakFitMass, peakFitMass_err, peakFitSigma, peakFitSigma_err);

    Double_t bgFitMin = peak->GetXaxis()->GetBinCenter(peak->FindFirstBinAbove(0.));
    Double_t bgFitMax = peak->GetXaxis()->GetBinCenter(peak->FindLastBinAbove(0.));
    // Double_t bgFitMin = peakFitMass - 15 * peakFitSigma; // Mean - 15*sigma (approx.)
    // Double_t bgFitMax = peakFitMass + 15 * peakFitSigma; // Mean + 15*sigma (approx.)
    bgReject = 4 * peakFitSigma; // 4*sigma region to be excluded from bg fit

    TFitResultPtr fFitResult_bg;

    if (!fisMC)
    {
        TString bgFncName = peakFnc->GetName();
        bgFncName.ReplaceAll("GausPol2", "_bgPol2Exclude");

        TF1 *bgFnc = new TF1(bgFncName.Data(), Pol2Exclude, bgFitMin, bgFitMax, 6);
        bgFnc->SetParNames("p[0]", "p[1]", "p[2]", "peakFitMass", "bgRejectRegion", "rejectFlag");
        bgFnc->FixParameter(3, peakFitMass); // "mass = mean of fit"
        bgFnc->FixParameter(4, bgReject);    // "4*sigma"
        bgFnc->SetLineColor(kGreen);
        bgFnc->SetFillColor(kYellow);
        bgFnc->SetFillStyle(3009);
        if (isGausPol2)
        {
            bgFnc->SetParameter(0, peakFnc->GetParameter(3));
            bgFnc->SetParameter(1, peakFnc->GetParameter(4));
            bgFnc->SetParameter(2, peakFnc->GetParameter(5));
        }
        else
        {
            bgFnc->SetParameter(0, peakFnc->GetParameter(7));
            bgFnc->SetParameter(1, peakFnc->GetParameter(6));
            bgFnc->SetParameter(2, peakFnc->GetParameter(5));
        }

        // Double_t fitMinBg = peakFitMass - 10 * peakFitSigma;
        // Double_t fitMaxBg = peakFitMass + 10 * peakFitSigma;
        bgFnc->FixParameter(5, 1); // "reject" (set 1 or 0)
                                   // fFitResult_bg = peak->Fit(bgFnc, optBg.Data(), "", fitMinBg, fitMaxBg);

        // Info("FitResults", "%s: Trying out initial background fit 3 times:", peak->GetName());

        for (int i = 0; i < 3; i++)
            peak->Fit(bgFnc, "QLNB MULTITHREAD", "", bgFitMin, bgFitMax);

        // Info("FitResults", "%s: Done. Trying out ACTUAL background fit with options: %s", peak->GetName(), fitOptBg.Data());

        fFitResult_bg = peak->Fit(bgFnc, fitOptBg.Data(), "", bgFitMin, bgFitMax);

        // check if bg fit was successful
        if (!fFitResult_bg->IsValid())
        {
            if (fitOptBg.Contains("I "))
            {
                Error("FitResults", "%s: Bg fit not converged. Fit status: %d, Fit Options: '%s'. Trying again without option \"I\" ... ", peak->GetName(), fFitResult_bg->Status(), fitOptBg.Data());
                fitOptBg.ReplaceAll("I ", " ");
                fFitResult_bg = peak->Fit(bgFnc, fitOptBg.Data(), "", bgFitMin, bgFitMax);

                if (!fFitResult_bg->IsValid())
                {
                    if (fitOptBg.Contains("QL"))
                    {
                        Error("FitResults", "%s: Bg fit STILL not converged. Fit status: %d, Fit Options: '%s'. Trying again a non-quiet fit without option \"L\" :", peak->GetName(), fFitResult_bg->Status(), fitOptBg.Data());
                        fitOptBg.ReplaceAll("QL", "");
                        fFitResult_bg = peak->Fit(bgFnc, fitOptBg.Data(), "", bgFitMin, bgFitMax);

                        if (!fFitResult_bg->IsValid())
                        {
                            Error("FitResults", "%s: Bg fit STILL not converged. Fit status: %d, Fit Options: '%s'. \e[1;31mGiving up\e[0m ...", peak->GetName(), fFitResult_bg->Status(), fitOptBg.Data());

                            if (gErrorIgnoreLevel < kError)
                                bgFnc->Print("V");

                            if (fFitResult_bg->CovMatrixStatus() < 2)
                            {
                                Error("FitResults", "%s: Covariance Matrix not accurate. Status: %d", peak->GetName(), fFitResult_bg->CovMatrixStatus());

                                if (gErrorIgnoreLevel < kError)
                                    fFitResult_bg->GetCovarianceMatrix().Print();
                            }
                        }
                    }
                    else
                    {
                        Error("FitResults", "%s: Bg fit still NOT converged. Fit status: %d, Fit Options: '%s'. \e[1;31mGiving up\e[0m ...", peak->GetName(), fFitResult_bg->Status(), fitOptBg.Data());

                        if (gErrorIgnoreLevel < kError)
                            bgFnc->Print("V");

                        if (fFitResult_bg->CovMatrixStatus() < 2)
                        {
                            Error("FitResults", "%s: Covariance Matrix not accurate. Status: %d", peak->GetName(), fFitResult_bg->CovMatrixStatus());

                            if (gErrorIgnoreLevel < kError)
                                fFitResult_bg->GetCovarianceMatrix().Print();
                        }
                    }
                }
            }
            else
            {
                Error("FitResults", "%s: Bg fit NOT converged. Fit status: %d, Fit Options: '%s'. \e[1;31mGiving up\e[0m ...", peak->GetName(), fFitResult_bg->Status(), fitOptBg.Data());

                if (gErrorIgnoreLevel < kError)
                    bgFnc->Print("V");

                if (fFitResult_bg->CovMatrixStatus() < 2)
                {
                    Error("FitResults", "%s: Covariance Matrix not accurate. Status: %d", peak->GetName(), fFitResult_bg->CovMatrixStatus());

                    if (gErrorIgnoreLevel < kError)
                        fFitResult_bg->GetCovarianceMatrix().Print();
                }
            }
        }
        bgFnc->FixParameter(5, 0); // reject = off (set 1 or 0)

        Info("FitResults", "%s: Background fit done. Adding bg func '%s' to hist list", peak->GetName(), bgFnc->GetName());

        // add fit and background function to histogram so it is automatically drawn
        // with hist
        peak->GetListOfFunctions()->Add(bgFnc);
    }

    /// getting min/max for deciding signal's fit range for integral calculation
    Double_t minInt = peakFitMass - bgReject; /// mean - 4*sigma
    Double_t maxInt = peakFitMass + bgReject; /// mean + 4*sigma
    Info("FitResults", "%s: Fitting successful! Peak region = [%f, %f]. Filling parameters ...", peak->GetName(), minInt, maxInt);

    resultParams = FillParams(h_bg, fFitResult_bg, peakFnc, peak, fFitResult, minInt, maxInt, MCgen, binsMC, fisMC, isGausPol2);
    peak->SetOption("X0E1");
    return resultParams;
    // }
    // catch (...)
    // {
    //     Error("FitResults", "%s: Possible runtime exception during fit. Try again!", peak->GetName());
    //     return nullptr; // Return nullptr on exception
    // }
}

TH1 *FillParams(TH1 *h_bg, TFitResultPtr fFitResult_bg, TF1 *sigBgFnc, TH1 *peak, TFitResultPtr fFitResult, Double_t minInt, Double_t maxInt, TH2 *MCgen, Int_t binsMC[], bool fisMC, bool isGausPol2)
{
    Double_t val, err, gen, genErr;
    Double_t val_intBC, err_intBC, val_intFF, err_intFF, eff, effErr;
    Double_t peakFitMass, peakFitSigma, peakFitMass_err, peakFitSigma_err = 0.;
    if (isGausPol2)
    {
        peakFitMass = sigBgFnc->GetParameter(1); // getting peak (mass) from peakFnc to be used for excluding bg
        peakFitMass_err = sigBgFnc->GetParError(1);
        peakFitSigma = TMath::Abs(sigBgFnc->GetParameter(2));
        peakFitSigma_err = sigBgFnc->GetParError(2);
    }
    else
        GetMeanSigmaDG(sigBgFnc, fFitResult, peakFitMass, peakFitMass_err, peakFitSigma, peakFitSigma_err);

    Int_t nPar = sigBgFnc->GetNpar();
    const Int_t nBins = 10 + nPar;
    Double_t par[nPar];
    Double_t parErr[nPar];
    sigBgFnc->GetParameters(par);
    for (Int_t i = 0; i < nPar; i++)
    {
        parErr[i] = sigBgFnc->GetParError(i);
    }

    TH1 *resultParams = new TH1D("resultParams", peak->GetTitle(), nBins, 0, nBins);
    Int_t iBin = 1;

    GetYieldBinCounting(h_bg, peak, fFitResult_bg, minInt, maxInt, fisMC, val, err);
    val_intBC = val;
    err_intBC = err;
    if (TMath::IsNaN(val_intBC) || TMath::IsNaN(err_intBC))
    {
        Error("FillParams", "%s: GetYieldBinCounting is NaN!", peak->GetName());
        delete resultParams;
        return nullptr;
    }

    resultParams->SetBinContent(iBin, val_intBC);
    resultParams->SetBinError(iBin, err_intBC);
    resultParams->GetXaxis()->SetBinLabel(iBin, "IntBC");
    iBin++;

    val = 0;
    err = 0;

    GetYieldFitFunction(peak, fFitResult, fFitResult_bg, minInt, maxInt, fisMC, val, err);
    val_intFF = val;
    err_intFF = err;
    if (TMath::IsNaN(val_intFF) || TMath::IsNaN(err_intFF))
    {
        Error("FillParams", "%s: GetYieldFitFunction is NaN!", peak->GetName());
        // delete resultParams;
        // return nullptr;
    }
    resultParams->SetBinContent(iBin, val_intFF);
    resultParams->SetBinError(iBin, err_intFF);
    resultParams->GetXaxis()->SetBinLabel(iBin, "IntFF");
    iBin++;

    resultParams->SetBinContent(iBin, peakFitMass);
    resultParams->SetBinError(iBin, peakFitMass_err);
    resultParams->GetXaxis()->SetBinLabel(iBin, "Mass");
    iBin++;

    resultParams->SetBinContent(iBin, peakFitSigma);
    resultParams->SetBinError(iBin, peakFitSigma_err);
    resultParams->GetXaxis()->SetBinLabel(iBin, "Width");
    iBin++;

    if (fisMC)
    {
        gen = GetGeneratedParticles(MCgen, binsMC, genErr);

        resultParams->SetBinContent(iBin, gen);
        resultParams->SetBinError(iBin, genErr);
        resultParams->GetXaxis()->SetBinLabel(iBin, "IntGen");
        iBin++;
        Info("FillParams: fisMC", "%s: IntGen = %f +/- %f", peak->GetName(), gen, genErr);

        if (gen != 0 && genErr != 0)
        {
            eff = val_intBC / gen;
            // effErr = err_intBC / genErr; /// ask if it should be computed like this?

            effErr = ErrorInRatio(val_intBC, err_intBC, gen, genErr);
            Info("FillParams: fisMC", "%s: Efficiency = %f; Error = %f", peak->GetName(), eff, effErr);
            resultParams->SetBinContent(iBin, eff);
            resultParams->SetBinError(iBin, effErr);
            resultParams->GetXaxis()->SetBinLabel(iBin, "Efficiency");
            iBin++;
        }
    }

    resultParams->SetBinContent(iBin, GetChi2(fFitResult));
    resultParams->SetBinError(iBin, 0);
    resultParams->GetXaxis()->SetBinLabel(iBin, "Chi2");
    iBin++;

    resultParams->SetBinContent(iBin, GetNdf(fFitResult));
    resultParams->SetBinError(iBin, 0);
    resultParams->GetXaxis()->SetBinLabel(iBin, "Ndf");
    iBin++;

    resultParams->SetBinContent(iBin, GetReducedChi2(fFitResult));
    resultParams->SetBinError(iBin, 0);
    resultParams->GetXaxis()->SetBinLabel(iBin, "ReducedChi2");
    iBin++;

    resultParams->SetBinContent(iBin, GetProb(fFitResult));
    resultParams->SetBinError(iBin, 0);
    resultParams->GetXaxis()->SetBinLabel(iBin, "Prob");
    iBin++;

    for (Int_t i = 0; i < nPar; i++)
    {
        if (TMath::IsNaN(par[i]) || TMath::IsNaN(parErr[i]))
        {
            delete resultParams;
            return nullptr;
        }
        resultParams->SetBinContent(iBin, par[i]);
        resultParams->SetBinError(iBin, parErr[i]);
        resultParams->GetXaxis()->SetBinLabel(iBin, sigBgFnc->GetParName(i));
        iBin++;
    }

    return resultParams;
}

void GetYieldBinCounting(TH1 *h_bg, TH1 *h, TFitResultPtr fFitResult_bg, Double_t minInt, Double_t maxInt, bool fisMC, Double_t &val,
                         Double_t &err, Double_t eps)
{

    if (!h)
        return;

    Int_t bin_minSig = h->GetXaxis()->FindBin(minInt);
    Int_t bin_maxSig = h->GetXaxis()->FindBin(maxInt);

    /// Histogram integrals always use whole bins. So, for comparisions with function integrals:
    Double_t fnIntMin = h->GetBinLowEdge(bin_minSig);     // left edge of the "bin_minSig" bin
    Double_t fnIntMax = h->GetBinLowEdge(bin_maxSig + 1); // right edge of the "bin_maxSig" bin
    Double_t histWidth = h->GetXaxis()->GetBinWidth(1);

    // if (h->GetEntries() > 50)
    // {
        val = h->IntegralAndError(bin_minSig, bin_maxSig, err);
    // }
    // else
        // Error("GetYieldBinCounting", "Unable to calculate signal IntBC because of low statistics! '%s' entries : %f ", h->GetName(), h->GetEntries());

    if (!fisMC)
    {
        Double_t bgVal = 0;
        Double_t bgErr = 0;
        if (h_bg == nullptr)
        {
            TF1 *bgFnc = (TF1 *)h->GetListOfFunctions()->At(1);

            // // make sub matrix with first three elements -> the last three are fixed user params
            // TMatrixD covMatrix_bg = fFitResult_bg->GetCovarianceMatrix();
            // // (0, 2, fFitResult_bg->GetCovarianceMatrix().GetMatrixArray());
            // // covMatrix_bg.Print();
            // TMatrixD covMatrixSub_bg = covMatrix_bg.GetSub(0, 2, 0, 2);
            // double bgFncPar[3];
            // for (Int_t iPar = 0; iPar < 3; iPar++)
            // {
            //     bgFncPar[iPar] = bgFnc->GetParameter(iPar);
            // }

            // if (h->GetEntries() > 50)
            // {
                bgVal = bgFnc->Integral(fnIntMin, fnIntMax, eps);
                // bgErr = bgFnc->IntegralError(fnIntMin, fnIntMax, bgFncPar, covMatrixSub_bg.GetMatrixArray(), eps);
                bgErr = bgFnc->IntegralError(fnIntMin, fnIntMax, bgFnc->GetParameters(), fFitResult_bg->GetCovarianceMatrix().GetMatrixArray(), eps);

                bgVal /= histWidth;
                bgErr /= histWidth;
            // }
            // else
                // Error("GetYieldBinCounting", "Unable to calculate bg integral because of low statistics! '%s' entries : %f ", h->GetName(), h->GetEntries());

            if (!(bgVal && bgErr && val && err))
                Error("GetYieldBinCounting", "IntBC: %s: [%.2f, %.2f], bins:[%d, %d]: bg(%s) = %.3f +/- %.3f, sig(%s) = %.3f +/- %.3f", h->GetName(), fnIntMin, fnIntMax, bin_minSig, bin_maxSig, bgFnc->GetName(), bgVal, bgErr, h->GetListOfFunctions()->At(0)->GetName(), val, err);
            else
                Info("GetYieldBinCounting", "IntBC: %s: [%.2f, %.2f], bins:[%d, %d]: bg(%s) = %.3f +/- %.3f, sig(%s) = %.3f +/- %.3f", h->GetName(), fnIntMin, fnIntMax, bin_minSig, bin_maxSig, bgFnc->GetName(), bgVal, bgErr, h->GetListOfFunctions()->At(0)->GetName(), val, err);
        }
        else
        {
            // Double_t bgWidth = h_bg->GetXaxis()->GetBinWidth(1);

            bgVal = h_bg->IntegralAndError(bin_minSig, bin_maxSig, bgErr);

            // bgVal /= bgWidth;
            // bgErr /= bgWidth;
        }
        // bgVal /= histWidth;
        // bgErr /= histWidth;
        ////
        // if (TMath::Abs(bgVal) < 1) // test for bg value == 0.00...
        // {
        // bgFnc->Print("V");
        // }
        val -= bgVal;
        err = TMath::Sqrt(TMath::Power(err, 2) + TMath::Power(bgErr, 2));
    }
}

void GetYieldFitFunction(TH1 *h, TFitResultPtr fFitResult_sig, TFitResultPtr fFitResult_bg, Double_t minInt, Double_t maxInt, bool fisMC, Double_t &val,
                         Double_t &err, Double_t eps)
{

    if (!h)
        return;

    // if (h->GetEntries() < 100)
    //     return;

    Double_t histWidth = h->GetXaxis()->GetBinWidth(1);

    TF1 *sigBgFnc = (TF1 *)h->GetListOfFunctions()->At(0);

    val = sigBgFnc->Integral(minInt, maxInt, eps);
    err = sigBgFnc->IntegralError(minInt, maxInt, fFitResult_sig->GetParams(), fFitResult_sig->GetCovarianceMatrix().GetMatrixArray(),
                                  eps);

    val /= histWidth;
    err /= histWidth;
    if (!fisMC)
    {
        TF1 *bgFnc = (TF1 *)h->GetListOfFunctions()->At(1);
        Double_t bg = bgFnc->Integral(minInt, maxInt, eps);

        Double_t bgErr = bgFnc->IntegralError(minInt, maxInt, fFitResult_bg->GetParams(),
                                              fFitResult_bg->GetCovarianceMatrix().GetMatrixArray(), eps);

        bg /= histWidth;
        bgErr /= histWidth;
        // Info("GetYieldFitFunction", "IntFF: %s: [%f, %f], bg = %f +/- %f, sig = %f +/- %f", h->GetName(), minInt, maxInt, bg, bgErr, val, err);

        val -= bg;
        err = TMath::Sqrt(TMath::Power(err, 2) + TMath::Power(bgErr, 2));
    }
}