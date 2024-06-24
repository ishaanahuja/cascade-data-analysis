#include <TFile.h>
#include <TCanvas.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TF1.h>
#include <TString.h>
#include <TFitResult.h>
#include <TFitResultPtr.h>
#include <TError.h>
#include <THashList.h>
#include <TSpectrum.h>
#include <TROOT.h>
#include <Math/MinimizerOptions.h>
#include <Math/IntegratorOptions.h>
#include <TFitter.h>
#include <TMatrixD.h>
#include <TStyle.h>

Double_t DoubleGausPol2(double *x, double *par)
{
    return par[0] * TMath::Gaus(x[0], par[1], par[2]) + par[3] * TMath::Gaus(x[0], par[1], par[4]) + par[5] + par[6] * x[0] + par[7] * x[0] * x[0];
}

Double_t DoubleGausPol3(double *x, double *par)
{
    return par[0] * TMath::Gaus(x[0], par[1], par[2]) + par[3] * TMath::Gaus(x[0], par[4], par[5]) + par[6] + par[7] * x[0] + par[8] * x[0] * x[0] + par[9] * x[0] * x[0] * x[0];
}

Double_t Pol2Exclude(double *x, double *par)
{
    // par[5] = rejectPoint(boolean); par[3] = peakFitMass; par[4]=peakFitSigma;
    if (par[5] && x[0] > (par[3] - par[4]) && x[0] < (par[3] + par[4]))
    {
        TF1::RejectPoint();
        return 0;
    }
    return par[0] * x[0] * x[0] + par[1] * x[0] + par[2];
}
Double_t GausPol2(double *x, double *par)
{
    return par[0] * TMath::Gaus(x[0], par[1], par[2]) + par[3] * x[0] * x[0] + par[4] * x[0] + par[5];
}
Double_t GetChi2(TFitResultPtr fFitResult)
{
    Double_t chi2 = fFitResult->Chi2();
    if (TMath::IsNaN(chi2))
        return 0;
    else
        return chi2;
}

Double_t GetNdf(TFitResultPtr fFitResult)
{
    Double_t ndf = fFitResult->Ndf();
    if (TMath::IsNaN(ndf))
        return 0;
    else
        return ndf;
}

Double_t GetReducedChi2(TFitResultPtr fFitResult)
{
    Double_t reducedChi2 = fFitResult->Chi2() / fFitResult->Ndf();
    if (TMath::IsNaN(reducedChi2))
        return 0;
    else
        return reducedChi2;
}

Double_t GetProb(TFitResultPtr fFitResult)
{
    Double_t prob = fFitResult->Prob();
    if (TMath::IsNaN(prob))
        return 0;
    else
        return prob;
}

Double_t GetGeneratedParticles(TH2 *h2, Int_t binsMC[], Double_t &genErr)
{
    Double_t gen = h2->IntegralAndError(binsMC[0], binsMC[1], binsMC[2], binsMC[3], genErr);
    return gen;
}

Double_t ErrorInRatio(Double_t A, Double_t Aerr, Double_t B, Double_t Berr)
{
    Double_t err = 0.;
    if (B == 0)
        err = -1.;
    else
    {
        Double_t errorfromtop = Aerr * Aerr / (B * B);
        Double_t errorfrombottom = ((A * A) / (B * B * B * B)) * Berr * Berr;
        err = TMath::Sqrt(errorfromtop + errorfrombottom);
    }
    return err;
}
TH1 *generateBg(bool bgTSpectrum, TH1 *h_source);
void GetMeanSigmaDG(TF1 *f_doubleGaus, TFitResultPtr lFitResultPtr, Double_t &mean, Double_t &mean_err, Double_t &sigma, Double_t &sigma_err);
TF1 *setFitParametersGaus(TF1 *peakFnc, Double_t *previousPtBinFitParams, TH1 *peak, Double_t mass, Double_t sigma, Bool_t usePrevious);
TF1 *setFitParametersDG(TF1 *peakFnc, Double_t *previousPtBinFitParams, TH1 *peak, Double_t mass, Double_t sigma, Bool_t usePrevious);

TH1 *fitResults(TH1 *h_bg, TF1 *peakFnc, TH1 *peak, Double_t fitMinSig, Double_t fitMaxSig, Int_t binsMC[], TH2 *MCgen, bool fisMC, Bool_t isGausPol2);

TH1 *fillParams(TH1 *h_bg, TFitResultPtr fFitResult_bg, TF1 *sigBgFnc, TH1 *peak, TFitResultPtr fFitResult, Double_t minInt, Double_t maxInt, TH2 *MCgen, Int_t binsMC[], bool fisMC, bool isGausPol2);

void GetYieldBinCounting(TH1 *h_bg, TH1 *h, TFitResultPtr fFitResult, Double_t minInt, Double_t maxInt, bool fisMC, Double_t &val,
                         Double_t &err, Double_t eps = 1e-3);

void GetYieldFitFunction(TH1 *h, TFitResultPtr fFitResult_sig, TFitResultPtr fFitResult_bg, Double_t minInt, Double_t maxInt, bool fisMC, Double_t &val,
                         Double_t &err, Double_t eps = 1e-3);

int testFitting(TString input = "", TString outputFilename = "", TString idAxis = "ye", bool fisMC = false, Int_t fitFunction = 2, bool xi = true, bool om = true, Int_t verbosity = kWarning)
{
    ROOT::EnableImplicitMT();
    ROOT::Math::IntegratorOneDimOptions::SetDefaultIntegrator("Adaptive");
    ROOT::Math::IntegratorOneDimOptions::SetDefaultAbsTolerance(1.E-3);
    ROOT::Math::IntegratorOneDimOptions::SetDefaultRelTolerance(1.E-3);
    ROOT::Math::MinimizerOptions::SetDefaultMaxFunctionCalls(10000);
    ROOT::Math::MinimizerOptions::SetDefaultMinimizer("Minuit2");
    ROOT::Math::MinimizerOptions::SetDefaultStrategy(2);
    ROOT::Math::MinimizerOptions::SetDefaultPrintLevel(1); // Fit printing: -1 = no printing, 0 (minimal) to 3 (max)
    gStyle->SetOptFit(1111);

    TFitter::SetMaxIterations(100000);
    TVirtualFitter::SetMaxIterations(100000);
    TFitter::SetPrecision(1e-3);
    TH1::AddDirectory(0);

    gPrintViaErrorHandler = kTRUE;
    gErrorIgnoreLevel = verbosity;

    Bool_t isGausPol2 = kFALSE;
    Bool_t isDoubleGausPol2 = kFALSE;
    Bool_t bgTSpectrum = kFALSE;
    /// Setting which fit functions to use:
    /// fitFunctions = [GausPol2, DoubleGausPol2, bgTSpectrum]
    /// e.g. fitFunctions = 110 -> fit with GausPol2 and DoubleGausPol2, don't fit bgTSpectrum

    if (fitFunction == 1)
    {
        isGausPol2 = kTRUE;
        Printf("Fit function GausPol2 will be used.");
    }
    else if (fitFunction == 2)
    {
        isDoubleGausPol2 = kTRUE;
        Printf("Fit function DoubleGausPol2 will be used.");
    }
    // if (fitFunction == '1')
    // {
    //     bgTSpectrum = kTRUE;
    //     Printf("Background fitting with TSpectrum will be used.");
    // }
    // }
    else
    {
        Fatal("arg: fitFunction", "Provided value %d is not valid. fitFunctions = [1:GausPol2, 2:DoubleGausPol2]", fitFunction);
        //  \ne.g. fitFunctions = \"110\" -> fit with GausPol2 and DoubleGausPol2, don't fit bgTSpectrum
        return 999;
    }

    TFile *f = TFile::Open(input);
    if (!f)
    {
        SysError("inputFile", "Error: Cannot open file '%s' !", input.Data());
        return 1;
    }

    THashList *list_eve = (THashList *)f->FindObjectAny("chists_eve_");
    TH1 *cent = (TH1 *)list_eve->FindObject("hcent");

    THashList *list_Xim = (THashList *)f->FindObjectAny("chists_Xim_");
    THashList *list_Xip = (THashList *)f->FindObjectAny("chists_Xip_");
    THashList *list_Omm = (THashList *)f->FindObjectAny("chists_Omm_");
    THashList *list_Omp = (THashList *)f->FindObjectAny("chists_Omp_");

    TH3 *Xim = (TH3 *)list_Xim->FindObject("h3_ptmasscent_def");
    TH3 *Xip = (TH3 *)list_Xip->FindObject("h3_ptmasscent_def");

    TH3 *Omm = (TH3 *)list_Omm->FindObject("h3_ptmasscent_def");
    TH3 *Omp = (TH3 *)list_Omp->FindObject("h3_ptmasscent_def");

    TH2 *MCGenXip = (TH2 *)list_Xip->FindObject("h2_gen");
    TH2 *MCGenXim = (TH2 *)list_Xim->FindObject("h2_gen");
    TH2 *MCGenOmp = (TH2 *)list_Omp->FindObject("h2_gen");
    TH2 *MCGenOmm = (TH2 *)list_Omm->FindObject("h2_gen");
    TH2D *MCGenXiC = nullptr;
    TH2D *MCGenOmC = nullptr;

    if (MCGenXip)
    {
        MCGenXiC = (TH2D *)MCGenXip->Clone();
        MCGenXiC->Add(MCGenXim);
    }
    if (MCGenOmp)
    {
        MCGenOmC = (TH2D *)MCGenOmp->Clone();
        MCGenOmC->Add(MCGenOmm);
    }

    Int_t binsMC[4] = {0};
    // }

    if (!((Xim && Xip) || (Omm && Omp)))
    {
        Fatal("histInput", "Histograms cannot be found. Aborting.");
        return 2;
    }

    TString output = outputFilename;
    if (outputFilename.IsNull())
    {
        output = input;
        output.ReplaceAll("AnalysisResults.root", "RsnProjection.root");
    }

    Info("outputCreation", "Saving output to '%s' ...", output.Data());
    TFile *out = TFile::Open(output.Data(), "RECREATE");
    if (!out)
    {
        SysError("outputCreation", "Error: Cannot open file '%s' !", output.Data());
        return 3;
    }
    out->mkdir("_allInt");
    if (xi)
    {
        out->mkdir("h_MassXim_pt_mult");
        out->mkdir("h_MassXip_pt_mult");
        out->mkdir("h_MassXiC_pt_mult"); /// for combination of + and -

        out->mkdir("h_MassXim_pt"); // mult: 0-100%
        out->mkdir("h_MassXip_pt");
        out->mkdir("h_MassXiC_pt");
    }

    if (om)
    {
        out->mkdir("h_MassOmm_pt_mult");
        out->mkdir("h_MassOmp_pt_mult");
        out->mkdir("h_MassOmC_pt_mult");

        out->mkdir("h_MassOmm_pt");
        out->mkdir("h_MassOmp_pt");
        out->mkdir("h_MassOmC_pt");
    }

    Double_t lMass_Xi = 1.32171;
    Double_t lMass_Om = 1.67245;
    Double_t nSigma = 10.;
    // Double_t sigma = 0.0025;
    // Double_t sigmaXi = 0.001875;
    // Double_t sigmaOm = 0.0018;
    Double_t fourSigXi = 0.0075;
    Double_t fourSigOm = 0.0072;

    // test: Emily's
    Double_t sigmaXi = 0.0018;
    Double_t sigmaOm = 0.0017;
    // Double_t fourSigXi = 0.01;
    // Double_t fourSigOm = 0.01;

    Double_t multbins_Xi[11] = {0, 5, 10, 15, 20, 30, 40, 50, 60, 80, 100}; // V0A
    Double_t multbins_Om[6] = {0, 5, 15, 30, 60, 100};                      // V0A
    Int_t nmultbins_Xi = sizeof(multbins_Xi) / sizeof(Double_t) - 1;
    Int_t nmultbins_Om = sizeof(multbins_Om) / sizeof(Double_t) - 1;

    Double_t ptbins_Xi[] = {0.8, 1.1, 1.3, 1.5, 1.7, 1.9, 2.1, 2.3, 2.5, 2.7, 2.9, 3.1, 3.5, 4, 5.5};
    Double_t ptbins_Om[] = {0.9, 1.6, 2., 2.4, 2.9, 3.5, 5};
    Int_t nptbins_Xi = sizeof(ptbins_Xi) / sizeof(Double_t) - 1;
    Int_t nptbins_Om = sizeof(ptbins_Om) / sizeof(Double_t) - 1;


    TH1I *h_multBinEntries_Xi = new TH1I("h_multBinEntries_Xi", "h_multBinEntries_Xi", nmultbins_Xi, multbins_Xi);
    TH1I *h_multBinEntries_Om = new TH1I("h_multBinEntries_Om", "h_multBinEntries_Om", nmultbins_Om, multbins_Om);

    for (Int_t multBinXi = 0; multBinXi < nmultbins_Xi; multBinXi++)
    {
        h_multBinEntries_Xi->SetBinContent(multBinXi + 1, cent->Integral(cent->FindBin(multbins_Xi[multBinXi]), cent->FindBin(multbins_Xi[multBinXi + 1])));
    }

    for (Int_t multBinOm = 0; multBinOm < nmultbins_Om; multBinOm++)
    {
        h_multBinEntries_Om->SetBinContent(multBinOm + 1, cent->Integral(cent->FindBin(multbins_Om[multBinOm]), cent->FindBin(multbins_Om[multBinOm + 1])));
    }
    out->cd();
    h_multBinEntries_Xi->Write();
    h_multBinEntries_Om->Write();

    TH1D *h_MassXim_pt_mult[nptbins_Xi][nmultbins_Xi];
    TH1D *h_MassXip_pt_mult[nptbins_Xi][nmultbins_Xi];
    TH1D *h_MassOmm_pt_mult[nptbins_Om][nmultbins_Om];
    TH1D *h_MassOmp_pt_mult[nptbins_Om][nmultbins_Om];
    TH1 *resultParXip_pt_mult[nptbins_Xi][nmultbins_Xi];
    TH1 *resultParXim_pt_mult[nptbins_Xi][nmultbins_Xi];
    TH1 *resultParOmp_pt_mult[nptbins_Om][nmultbins_Om];
    TH1 *resultParOmm_pt_mult[nptbins_Om][nmultbins_Om];

    /// +/- combination:
    TH1D *h_MassXiC_pt_mult[nptbins_Xi][nmultbins_Xi];
    TH1D *h_MassOmC_pt_mult[nptbins_Om][nmultbins_Om];
    TH1 *resultParXiC_pt_mult[nptbins_Xi][nmultbins_Xi];
    TH1 *resultParOmC_pt_mult[nptbins_Om][nmultbins_Om];

    TH1D *h_MassXim_pt[nptbins_Xi];
    TH1D *h_MassXip_pt[nptbins_Xi];
    TH1D *h_MassOmm_pt[nptbins_Om];
    TH1D *h_MassOmp_pt[nptbins_Om];
    TH1 *resultParXip_pt[nptbins_Xi];
    TH1 *resultParXim_pt[nptbins_Xi];
    TH1 *resultParOmp_pt[nptbins_Om];
    TH1 *resultParOmm_pt[nptbins_Om];

    /// +/- combination:
    TH1D *h_MassXiC_pt[nptbins_Xi];
    TH1D *h_MassOmC_pt[nptbins_Om];
    TH1 *resultParXiC_pt[nptbins_Xi];
    TH1 *resultParOmC_pt[nptbins_Om];

    /// fit functions for double gaus fit:
    // TF1 *sigFnc_XiP[nptbins_Xi];
    TF1 *f_Sig_Xip;
    TF1 *f_Sig_Xim;
    TF1 *f_Sig_XiC;
    TF1 *f_Sig_Xim_pt[nptbins_Xi];
    TF1 *f_Sig_Xip_pt[nptbins_Xi];
    TF1 *f_Sig_XiC_pt[nptbins_Xi];
    TF1 *f_Sig_Xim_pt_mult[nptbins_Xi][nmultbins_Xi];
    TF1 *f_Sig_Xip_pt_mult[nptbins_Xi][nmultbins_Xi];
    TF1 *f_Sig_XiC_pt_mult[nptbins_Xi][nmultbins_Xi];

    TF1 *f_Sig_Omp;
    TF1 *f_Sig_Omm;
    TF1 *f_Sig_OmC;
    TF1 *f_Sig_Omm_pt[nptbins_Om];
    TF1 *f_Sig_Omp_pt[nptbins_Om];
    TF1 *f_Sig_OmC_pt[nptbins_Om];
    TF1 *f_Sig_Omm_pt_mult[nptbins_Om][nmultbins_Om];
    TF1 *f_Sig_Omp_pt_mult[nptbins_Om][nmultbins_Om];
    TF1 *f_Sig_OmC_pt_mult[nptbins_Om][nmultbins_Om];
    /// background estimation hist through TSpectrum

    TH1 *h_bgXim_pt[nptbins_Xi];
    TH1 *h_bgXip_pt[nptbins_Xi];
    TH1 *h_bgXiC_pt[nptbins_Xi];
    TH1 *h_bgXim_pt_mult[nptbins_Xi][nmultbins_Xi];
    TH1 *h_bgXip_pt_mult[nptbins_Xi][nmultbins_Xi];
    TH1 *h_bgXiC_pt_mult[nptbins_Xi][nmultbins_Xi];

    TH1 *h_bgOmm_pt[nptbins_Om];
    TH1 *h_bgOmp_pt[nptbins_Om];
    TH1 *h_bgOmC_pt[nptbins_Om];
    TH1 *h_bgOmm_pt_mult[nptbins_Om][nmultbins_Om];
    TH1 *h_bgOmp_pt_mult[nptbins_Om][nmultbins_Om];
    TH1 *h_bgOmC_pt_mult[nptbins_Om][nmultbins_Om];

    // Signal+Bg (peak)
    TH1 *h_MassXim = (TH1 *)Xim->Project3D(idAxis);
    h_MassXim->SetNameTitle("h_MassXim", "Invariant Mass #Xi^{-}");
    TH1 *h_MassXip = (TH1 *)Xip->Project3D(idAxis);
    h_MassXip->SetNameTitle("h_MassXip", "Invariant Mass #Xi^{+}");
    TH1 *h_MassXiC = (TH1 *)h_MassXim->Clone();
    h_MassXiC->Add(h_MassXip);
    h_MassXiC->SetNameTitle("h_MassXiC", "Invariant Mass #Xi^{+} + #Xi^{-}");

    TH1 *h_MassOmm = (TH1 *)Omm->Project3D(idAxis);
    h_MassOmm->SetNameTitle("h_MassOmm", "Invariant Mass #Omega^{-}");
    TH1 *h_MassOmp = (TH1 *)Omp->Project3D(idAxis);
    h_MassOmp->SetNameTitle("h_MassOmp", "Invariant Mass #Omega^{+}");
    TH1 *h_MassOmC = (TH1 *)h_MassOmm->Clone();
    h_MassOmC->Add(h_MassOmp);
    h_MassOmC->SetNameTitle("h_MassOmC", "Invariant Mass #Omega^{+} + #Omega^{-}");

    TH1 *resultParams_Xip_allInt, *resultParams_Xim_allInt, *resultParams_XiC_allInt;
    TH1 *resultParams_Omp_allInt, *resultParams_Omm_allInt, *resultParams_OmC_allInt;

    Double_t fitMinSig_Xi = lMass_Xi - nSigma * sigmaXi;
    Double_t fitMaxSig_Xi = lMass_Xi + nSigma * sigmaXi;
    Double_t fitMinSig_Om = lMass_Om - nSigma * sigmaOm;
    Double_t fitMaxSig_Om = lMass_Om + nSigma * sigmaOm;

    /// Fitting gaussian on signal for getting sigma for bg
    // TF1 *gausXi = new TF1("gausXi", "gaus", lMass_Xi - 0.04, lMass_Xi + 0.04);
    // gausXi->SetLineColor(kBlack);
    // h_MassXim->Fit("gausXi", "REM+");
    // Printf("mean = %f, mean error = %f \nsigma = %f, sigma error = %f", gausXi->GetParameter(1), gausXi->GetParError(1), gausXi->GetParameter(2), gausXi->GetParError(2));
    // TF1 *gausOm = new TF1("gausOm", "gaus", lMass_Om - 0.04, lMass_Om + 0.04);
    // gausOm->SetLineColor(kBlack);
    // h_MassOmp->Fit("gausOm", "REM+");
    // Printf("mean = %f, mean error = %f \nsigma = %f, sigma error = %f", gausOm->GetParameter(1), gausOm->GetParError(1), gausOm->GetParameter(2), gausOm->GetParError(2));
    if (xi)
    {

        /// pt+mult integrated case
        binsMC[0] = 1;
        binsMC[1] = nptbins_Xi;
        binsMC[2] = 1;
        binsMC[3] = nmultbins_Xi;
        Info("Xi_allInt", "Starting Xi pt+mult integrated analysis...");

        if (isDoubleGausPol2)
        {
            f_Sig_Xip = new TF1("DblGausPol2_Xip", DoubleGausPol2, fitMinSig_Xi, fitMaxSig_Xi, 8);
            f_Sig_Xim = new TF1("DblGausPol2_Xim", DoubleGausPol2, fitMinSig_Xi, fitMaxSig_Xi, 8);
            f_Sig_XiC = new TF1("DblGausPol2_XiC", DoubleGausPol2, fitMinSig_Xi, fitMaxSig_Xi, 8);
            f_Sig_Xip->SetNpx(1000);
            f_Sig_Xim->SetNpx(1000);
            f_Sig_XiC->SetNpx(1000);

            f_Sig_Xip = setFitParametersDG(f_Sig_Xip, nullptr, h_MassXip, lMass_Xi, sigmaXi, kFALSE);
            f_Sig_Xim = setFitParametersDG(f_Sig_Xim, nullptr, h_MassXim, lMass_Xi, sigmaXi, kFALSE);
            f_Sig_XiC = setFitParametersDG(f_Sig_XiC, nullptr, h_MassXiC, lMass_Xi, sigmaXi, kFALSE);
        }
        if (isGausPol2)
        {
            f_Sig_Xip = new TF1("GausPol2_Xip", GausPol2, fitMinSig_Xi, fitMaxSig_Xi, 6);
            f_Sig_Xim = new TF1("GausPol2_Xim", GausPol2, fitMinSig_Xi, fitMaxSig_Xi, 6);
            f_Sig_XiC = new TF1("GausPol2_XiC", GausPol2, fitMinSig_Xi, fitMaxSig_Xi, 6);
            f_Sig_Xip->SetNpx(1000);
            f_Sig_Xim->SetNpx(1000);
            f_Sig_XiC->SetNpx(1000);

            f_Sig_Xip = setFitParametersGaus(f_Sig_Xip, nullptr, h_MassXip, lMass_Xi, sigmaXi, kFALSE);
            f_Sig_Xim = setFitParametersGaus(f_Sig_Xim, nullptr, h_MassXim, lMass_Xi, sigmaXi, kFALSE);
            f_Sig_XiC = setFitParametersGaus(f_Sig_XiC, nullptr, h_MassXiC, lMass_Xi, sigmaXi, kFALSE);
        }
        resultParams_Xip_allInt = fitResults(nullptr, f_Sig_Xip, h_MassXip, fitMinSig_Xi, fitMaxSig_Xi, binsMC, MCGenXip, fisMC, isGausPol2);
        resultParams_Xim_allInt = fitResults(nullptr, f_Sig_Xim, h_MassXim, fitMinSig_Xi, fitMaxSig_Xi, binsMC, MCGenXim, fisMC, isGausPol2);
        resultParams_XiC_allInt = fitResults(nullptr, f_Sig_XiC, h_MassXiC, fitMinSig_Xi, fitMaxSig_Xi, binsMC, MCGenXiC, fisMC, isGausPol2);

        out->cd("_allInt");
        h_MassXip->Write();
        h_MassXim->Write();
        h_MassXiC->Write();
        resultParams_Xip_allInt->SetNameTitle("resultParams_Xip_allInt", "Result Parameters #Xi^{+}");
        resultParams_Xim_allInt->SetNameTitle("resultParams_Xim_allInt", "Result Parameters #Xi^{-}");
        resultParams_XiC_allInt->SetNameTitle("resultParams_XiC_allInt", "Result Parameters #Xi^{+} + #Xi^{-}");
        resultParams_Xip_allInt->Write();
        resultParams_Xim_allInt->Write();
        resultParams_XiC_allInt->Write();

        Info("Xi_allInt", "Finished writing allInt hists...");

        for (Int_t ptBinXi = 0; ptBinXi < nptbins_Xi; ptBinXi++)
        {
            Info("Xi_multInt", "Loop index: pT bin = %d ...", ptBinXi);

            if (fisMC)
            {
                // plotting mult integrated efficiency as a function of pT
                binsMC[0] = ptBinXi + 1;
                binsMC[1] = ptBinXi + 1;
                binsMC[2] = 1;
                binsMC[3] = nmultbins_Xi;
            }
            h_MassXim_pt[ptBinXi] = (TH1D *)Xim->ProjectionY(TString::Format(("h_MassXim_pt[%d]"), ptBinXi), ptBinXi + 1, ptBinXi + 1, 1, nmultbins_Xi, "e");
            h_MassXim_pt[ptBinXi]->SetTitle(TString::Format(("Invariant Mass #Xi^{-}: p_{T}<%.1f,%.1f>, Mult<%.1f,%.1f>"), ptbins_Xi[ptBinXi], ptbins_Xi[ptBinXi + 1], multbins_Xi[0], multbins_Xi[nmultbins_Xi]));

            h_MassXip_pt[ptBinXi] = (TH1D *)Xip->ProjectionY(TString::Format(("h_MassXip_pt[%d]"), ptBinXi), ptBinXi + 1, ptBinXi + 1, 1, nmultbins_Xi, "e");
            h_MassXip_pt[ptBinXi]->SetTitle(TString::Format(("Invariant Mass #Xi^{+}: p_{T}<%.1f,%.1f>, Mult<%.1f,%.1f>"), ptbins_Xi[ptBinXi], ptbins_Xi[ptBinXi + 1], multbins_Xi[0], multbins_Xi[nmultbins_Xi]));

            h_MassXiC_pt[ptBinXi] = (TH1D *)h_MassXip_pt[ptBinXi]->Clone(TString::Format(("h_MassXiC_pt[%d]"), ptBinXi));
            h_MassXiC_pt[ptBinXi]->Add(h_MassXim_pt[ptBinXi]);
            h_MassXiC_pt[ptBinXi]->SetTitle(TString::Format(("Invariant Mass #Xi: p_{T}<%.1f,%.1f>, Mult<%.1f,%.1f>"), ptbins_Xi[ptBinXi], ptbins_Xi[ptBinXi + 1], multbins_Xi[0], multbins_Xi[nmultbins_Xi]));

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
                    f_Sig_Xip_pt[ptBinXi] = setFitParametersDG(f_Sig_Xip_pt[ptBinXi], f_Sig_Xip->GetParameters(), h_MassXip_pt[ptBinXi], lMass_Xi, sigmaXi, kTRUE);
                    f_Sig_Xim_pt[ptBinXi] = setFitParametersDG(f_Sig_Xim_pt[ptBinXi], f_Sig_Xim->GetParameters(), h_MassXim_pt[ptBinXi], lMass_Xi, sigmaXi, kTRUE);
                    f_Sig_XiC_pt[ptBinXi] = setFitParametersDG(f_Sig_XiC_pt[ptBinXi], f_Sig_XiC->GetParameters(), h_MassXiC_pt[ptBinXi], lMass_Xi, sigmaXi, kTRUE);
                }
                else
                {

                    f_Sig_Xip_pt[ptBinXi] = setFitParametersDG(f_Sig_Xip_pt[ptBinXi], f_Sig_Xip_pt[ptBinXi - 1]->GetParameters(), h_MassXip_pt[ptBinXi], lMass_Xi, sigmaXi, kTRUE);
                    f_Sig_Xim_pt[ptBinXi] = setFitParametersDG(f_Sig_Xim_pt[ptBinXi], f_Sig_Xim_pt[ptBinXi - 1]->GetParameters(), h_MassXim_pt[ptBinXi], lMass_Xi, sigmaXi, kTRUE);
                    f_Sig_XiC_pt[ptBinXi] = setFitParametersDG(f_Sig_XiC_pt[ptBinXi], f_Sig_XiC_pt[ptBinXi - 1]->GetParameters(), h_MassXiC_pt[ptBinXi], lMass_Xi, sigmaXi, kTRUE);
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
                    f_Sig_Xip_pt[ptBinXi] = setFitParametersGaus(f_Sig_Xip_pt[ptBinXi], f_Sig_Xip->GetParameters(), h_MassXip_pt[ptBinXi], lMass_Xi, sigmaXi, kTRUE);
                    f_Sig_Xim_pt[ptBinXi] = setFitParametersGaus(f_Sig_Xim_pt[ptBinXi], f_Sig_Xim->GetParameters(), h_MassXim_pt[ptBinXi], lMass_Xi, sigmaXi, kTRUE);
                    f_Sig_XiC_pt[ptBinXi] = setFitParametersGaus(f_Sig_XiC_pt[ptBinXi], f_Sig_XiC->GetParameters(), h_MassXiC_pt[ptBinXi], lMass_Xi, sigmaXi, kTRUE);
                }
                else
                {
                    f_Sig_Xip_pt[ptBinXi] = setFitParametersGaus(f_Sig_Xip_pt[ptBinXi], f_Sig_Xip_pt[ptBinXi - 1]->GetParameters(), h_MassXip_pt[ptBinXi], lMass_Xi, sigmaXi, kTRUE);
                    f_Sig_Xim_pt[ptBinXi] = setFitParametersGaus(f_Sig_Xim_pt[ptBinXi], f_Sig_Xim_pt[ptBinXi - 1]->GetParameters(), h_MassXim_pt[ptBinXi], lMass_Xi, sigmaXi, kTRUE);
                    f_Sig_XiC_pt[ptBinXi] = setFitParametersGaus(f_Sig_XiC_pt[ptBinXi], f_Sig_XiC_pt[ptBinXi - 1]->GetParameters(), h_MassXiC_pt[ptBinXi], lMass_Xi, sigmaXi, kTRUE);
                }
            }

            // h_bgXip_pt[ptBinXi] = generateBg(bgTSpectrum, h_MassXip_pt[ptBinXi]);
            resultParXip_pt[ptBinXi] = fitResults(nullptr, f_Sig_Xip_pt[ptBinXi], h_MassXip_pt[ptBinXi], fitMinSig_Xi, fitMaxSig_Xi, binsMC, MCGenXip, fisMC, isGausPol2);
            resultParXip_pt[ptBinXi]->SetName(TString::Format(("resultParXip_pt[%d]"), ptBinXi));
            resultParXip_pt[ptBinXi]->SetTitle(TString::Format(("Result Parameters #Xi^{+}: p_{T}<%.1f,%.1f>, Mult<%.1f,%.1f>"), ptbins_Xi[ptBinXi], ptbins_Xi[ptBinXi + 1], multbins_Xi[0], multbins_Xi[nmultbins_Xi]));

            // h_bgXim_pt[ptBinXi] = generateBg(bgTSpectrum, h_MassXim_pt[ptBinXi]);
            resultParXim_pt[ptBinXi] = fitResults(nullptr, f_Sig_Xim_pt[ptBinXi], h_MassXim_pt[ptBinXi], fitMinSig_Xi, fitMaxSig_Xi, binsMC, MCGenXim, fisMC, isGausPol2);
            resultParXim_pt[ptBinXi]->SetName(TString::Format(("resultParXim_pt[%d]"), ptBinXi));
            resultParXim_pt[ptBinXi]->SetTitle(TString::Format(("Result Parameters #Xi^{-}: p_{T}<%.1f,%.1f>, Mult<%.1f,%.1f>"), ptbins_Xi[ptBinXi], ptbins_Xi[ptBinXi + 1], multbins_Xi[0], multbins_Xi[nmultbins_Xi]));

            // h_bgXiC_pt[ptBinXi] = generateBg(bgTSpectrum, h_MassXiC_pt[ptBinXi]);
            resultParXiC_pt[ptBinXi] = fitResults(nullptr, f_Sig_XiC_pt[ptBinXi], h_MassXiC_pt[ptBinXi], fitMinSig_Xi, fitMaxSig_Xi, binsMC, MCGenXiC, fisMC, isGausPol2);
            resultParXiC_pt[ptBinXi]->SetName(TString::Format(("resultParXiC_pt[%d]"), ptBinXi));
            resultParXiC_pt[ptBinXi]->SetTitle(TString::Format(("Result Parameters #Xi: p_{T}<%.1f,%.1f>, Mult<%.1f,%.1f>"), ptbins_Xi[ptBinXi], ptbins_Xi[ptBinXi + 1], multbins_Xi[0], multbins_Xi[nmultbins_Xi]));

            out->cd("h_MassXim_pt");
            h_MassXim_pt[ptBinXi]->Write();
            resultParXim_pt[ptBinXi]->Write();

            out->cd("h_MassXip_pt");
            h_MassXip_pt[ptBinXi]->Write();
            resultParXip_pt[ptBinXi]->Write();

            out->cd("h_MassXiC_pt");
            h_MassXiC_pt[ptBinXi]->Write();
            resultParXiC_pt[ptBinXi]->Write();

            Info("Xi_multInt", "Successfully written hists for pT bin = %d ...", ptBinXi);
        }

        /// Starting fully differential case
        for (Int_t ptBinXi = 0; ptBinXi < nptbins_Xi; ptBinXi++)
        {
            if (fisMC)
            {
                binsMC[0] = ptBinXi + 1;
                binsMC[1] = ptBinXi + 1;
            }

            for (Int_t multBinXi = 0; multBinXi < nmultbins_Xi; multBinXi++)
            {
                Info("Xi_diff", "Loop index: pT bin = %d , Mult bin = %d ...", ptBinXi, multBinXi);

                if (fisMC) /// fully differential case
                {
                    binsMC[2] = multBinXi + 1;
                    binsMC[3] = multBinXi + 1;
                }

                h_MassXim_pt_mult[ptBinXi][multBinXi] = (TH1D *)Xim->ProjectionY(TString::Format(("h_MassXim_pt_mult[%d][%d]"), ptBinXi, multBinXi), ptBinXi + 1, ptBinXi + 1, multBinXi + 1, multBinXi + 1, "e");
                h_MassXim_pt_mult[ptBinXi][multBinXi]->SetTitle(TString::Format(("Invariant Mass #Xi^{-}: p_{T}<%.1f,%.1f>, Mult<%.1f,%.1f>"), ptbins_Xi[ptBinXi], ptbins_Xi[ptBinXi + 1], multbins_Xi[multBinXi], multbins_Xi[multBinXi + 1]));

                h_MassXip_pt_mult[ptBinXi][multBinXi] = (TH1D *)Xip->ProjectionY(TString::Format(("h_MassXip_pt_mult[%d][%d]"), ptBinXi, multBinXi), ptBinXi + 1, ptBinXi + 1, multBinXi + 1, multBinXi + 1, "e");
                h_MassXip_pt_mult[ptBinXi][multBinXi]->SetTitle(TString::Format(("Invariant Mass #Xi^{+}: p_{T}<%.1f,%.1f>, Mult<%.1f,%.1f>"), ptbins_Xi[ptBinXi], ptbins_Xi[ptBinXi + 1], multbins_Xi[multBinXi], multbins_Xi[multBinXi + 1]));

                h_MassXiC_pt_mult[ptBinXi][multBinXi] = (TH1D *)h_MassXip_pt_mult[ptBinXi][multBinXi]->Clone(TString::Format(("h_MassXiC_pt_mult[%d][%d]"), ptBinXi, multBinXi));
                h_MassXiC_pt_mult[ptBinXi][multBinXi]->Add(h_MassXim_pt_mult[ptBinXi][multBinXi]);
                h_MassXiC_pt_mult[ptBinXi][multBinXi]->SetTitle(TString::Format(("Invariant Mass #Xi: p_{T}<%.1f,%.1f>, Mult<%.1f,%.1f>"), ptbins_Xi[ptBinXi], ptbins_Xi[ptBinXi + 1], multbins_Xi[multBinXi], multbins_Xi[multBinXi + 1]));

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
                    f_Sig_Xip_pt_mult[ptBinXi][multBinXi] = setFitParametersDG(f_Sig_Xip_pt_mult[ptBinXi][multBinXi], f_Sig_Xip_pt[ptBinXi]->GetParameters(), h_MassXip_pt_mult[ptBinXi][multBinXi], lMass_Xi, resultParXip_pt[ptBinXi]->GetBinContent(4), kTRUE); // ptBin needs to be >1 in this case to ensure correct flow in setFitParametersDG to not reset fn params
                    f_Sig_Xim_pt_mult[ptBinXi][multBinXi] = setFitParametersDG(f_Sig_Xim_pt_mult[ptBinXi][multBinXi], f_Sig_Xim_pt[ptBinXi]->GetParameters(), h_MassXim_pt_mult[ptBinXi][multBinXi], lMass_Xi, resultParXim_pt[ptBinXi]->GetBinContent(4), kTRUE);
                    f_Sig_XiC_pt_mult[ptBinXi][multBinXi] = setFitParametersDG(f_Sig_XiC_pt_mult[ptBinXi][multBinXi], f_Sig_XiC_pt[ptBinXi]->GetParameters(), h_MassXiC_pt_mult[ptBinXi][multBinXi], lMass_Xi, resultParXiC_pt[ptBinXi]->GetBinContent(4), kTRUE);
                    // }
                    // else
                    // {
                    //     f_Sig_Xip_pt_mult[ptBinXi][multBinXi] = setFitParametersDG(f_Sig_Xip_pt_mult[ptBinXi][multBinXi], f_Sig_Xip_pt_mult[ptBinXi - 1][multBinXi]->GetParameters(), h_MassXip_pt_mult[ptBinXi][multBinXi], lMass_Xi, sigmaXi, ptBinXi + 1, kTRUE); // multInt needs to be true for correct flow in setFitParametersDG -> TODO: fix it later in setFitParametersDG
                    //     f_Sig_Xim_pt_mult[ptBinXi][multBinXi] = setFitParametersDG(f_Sig_Xim_pt_mult[ptBinXi][multBinXi], f_Sig_Xim_pt_mult[ptBinXi - 1][multBinXi]->GetParameters(), h_MassXim_pt_mult[ptBinXi][multBinXi], lMass_Xi, sigmaXi, ptBinXi + 1, kTRUE);
                    //     f_Sig_XiC_pt_mult[ptBinXi][multBinXi] = setFitParametersDG(f_Sig_XiC_pt_mult[ptBinXi][multBinXi], f_Sig_XiC_pt_mult[ptBinXi - 1][multBinXi]->GetParameters(), h_MassXiC_pt_mult[ptBinXi][multBinXi], lMass_Xi, sigmaXi, ptBinXi + 1, kTRUE);
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
                    f_Sig_Xip_pt_mult[ptBinXi][multBinXi] = setFitParametersGaus(f_Sig_Xip_pt_mult[ptBinXi][multBinXi], f_Sig_Xip_pt[ptBinXi]->GetParameters(), h_MassXip_pt_mult[ptBinXi][multBinXi], lMass_Xi, sigmaXi, kTRUE); // ptBin needs to be >1 in this case to ensure correct flow in setFitParametersGaus to not reset fn params
                    f_Sig_Xim_pt_mult[ptBinXi][multBinXi] = setFitParametersGaus(f_Sig_Xim_pt_mult[ptBinXi][multBinXi], f_Sig_Xim_pt[ptBinXi]->GetParameters(), h_MassXim_pt_mult[ptBinXi][multBinXi], lMass_Xi, sigmaXi, kTRUE);
                    f_Sig_XiC_pt_mult[ptBinXi][multBinXi] = setFitParametersGaus(f_Sig_XiC_pt_mult[ptBinXi][multBinXi], f_Sig_XiC_pt[ptBinXi]->GetParameters(), h_MassXiC_pt_mult[ptBinXi][multBinXi], lMass_Xi, sigmaXi, kTRUE);
                    // }
                    // else
                    // {
                    // f_Sig_Xip_pt_mult[ptBinXi][multBinXi] = setFitParametersGaus(f_Sig_Xip_pt_mult[ptBinXi][multBinXi], f_Sig_Xip_pt_mult[ptBinXi - 1][multBinXi]->GetParameters(), h_MassXip_pt_mult[ptBinXi][multBinXi], lMass_Xi, sigmaXi, kTRUE); // multInt needs to be true for correct flow in setFitParametersGaus -> TODO: fix it later in setFitParametersDG
                    // f_Sig_Xim_pt_mult[ptBinXi][multBinXi] = setFitParametersGaus(f_Sig_Xim_pt_mult[ptBinXi][multBinXi], f_Sig_Xim_pt_mult[ptBinXi - 1][multBinXi]->GetParameters(), h_MassXim_pt_mult[ptBinXi][multBinXi], lMass_Xi, sigmaXi, kTRUE);
                    // f_Sig_XiC_pt_mult[ptBinXi][multBinXi] = setFitParametersGaus(f_Sig_XiC_pt_mult[ptBinXi][multBinXi], f_Sig_XiC_pt_mult[ptBinXi - 1][multBinXi]->GetParameters(), h_MassXiC_pt_mult[ptBinXi][multBinXi], lMass_Xi, sigmaXi, kTRUE);
                    // }
                }

                // h_bgXip_pt_mult[ptBinXi][multBinXi] = generateBg(bgTSpectrum, h_MassXip_pt_mult[ptBinXi][multBinXi]);
                // resultParXip_pt_mult[ptBinXi][multBinXi] = fitResults(nullptr, f_Sig_Xip_pt[ptBinXi], h_MassXip_pt_mult[ptBinXi][multBinXi], fitMinSig_Xi, fitMaxSig_Xi, &par_allintP[0], &par_allint_errorsP[0], binsMC, MCGenXip, fisMC, kFALSE);
                resultParXip_pt_mult[ptBinXi][multBinXi] = fitResults(nullptr, f_Sig_Xip_pt_mult[ptBinXi][multBinXi], h_MassXip_pt_mult[ptBinXi][multBinXi], fitMinSig_Xi, fitMaxSig_Xi, binsMC, MCGenXip, fisMC, isGausPol2);
                resultParXip_pt_mult[ptBinXi][multBinXi]->SetName(TString::Format(("resultParXip_pt_mult[%d][%d]"), ptBinXi, multBinXi));
                resultParXip_pt_mult[ptBinXi][multBinXi]->SetTitle(TString::Format(("Result Parameters #Xi^{+}: p_{T}<%.1f,%.1f>, Mult<%.1f,%.1f>"), ptbins_Xi[ptBinXi], ptbins_Xi[ptBinXi + 1], multbins_Xi[multBinXi], multbins_Xi[multBinXi + 1]));
                // Printf(">>>>>>>>>>>>Got %s............", resultParXip_pt_mult[ptBinXi][multBinXi]->GetName());

                // h_bgXim_pt[ptBinXi] = generateBg(bgTSpectrum, h_MassXim_pt[ptBinXi]);
                resultParXim_pt_mult[ptBinXi][multBinXi] = fitResults(nullptr, f_Sig_Xim_pt_mult[ptBinXi][multBinXi], h_MassXim_pt_mult[ptBinXi][multBinXi], fitMinSig_Xi, fitMaxSig_Xi, binsMC, MCGenXim, fisMC, isGausPol2);
                resultParXim_pt_mult[ptBinXi][multBinXi]->SetName(TString::Format(("resultParXim_pt_mult[%d][%d]"), ptBinXi, multBinXi));
                resultParXim_pt_mult[ptBinXi][multBinXi]->SetTitle(TString::Format(("Result Parameters #Xi^{-}: p_{T}<%.1f,%.1f>, Mult<%.1f,%.1f>"), ptbins_Xi[ptBinXi], ptbins_Xi[ptBinXi + 1], multbins_Xi[multBinXi], multbins_Xi[multBinXi + 1]));
                // Printf(">>>>>>>>>>>>Got %s............", resultParXim_pt_mult[ptBinXi][multBinXi]->GetName());

                // h_bgXiC_pt[ptBinXi] = generateBg(bgTSpectrum, h_MassXiC_pt_mult[ptBinXi][multBinXi]);
                resultParXiC_pt_mult[ptBinXi][multBinXi] = fitResults(nullptr, f_Sig_XiC_pt_mult[ptBinXi][multBinXi], h_MassXiC_pt_mult[ptBinXi][multBinXi], fitMinSig_Xi, fitMaxSig_Xi, binsMC, MCGenXiC, fisMC, isGausPol2);
                resultParXiC_pt_mult[ptBinXi][multBinXi]->SetName(TString::Format(("resultParXiC_pt_mult[%d][%d]"), ptBinXi, multBinXi));
                resultParXiC_pt_mult[ptBinXi][multBinXi]->SetTitle(TString::Format(("Result Parameters #Xi: p_{T}<%.1f,%.1f>, Mult<%.1f,%.1f>"), ptbins_Xi[ptBinXi], ptbins_Xi[ptBinXi + 1], multbins_Xi[multBinXi], multbins_Xi[multBinXi + 1]));
                // Printf(">>>>>>>>>>>>Got %s............", resultParXiC_pt_mult[ptBinXi][multBinXi]->GetName());

                out->cd("h_MassXim_pt_mult");
                // if (h_bgXim_pt_mult[ptBinXi][multBinXi])
                //     h_bgXim_pt_mult[ptBinXi][multBinXi]->Write();
                // if (h_MassXim_pt_mult[ptBinXi][multBinXi])
                h_MassXim_pt_mult[ptBinXi][multBinXi]->Write();
                // if (resultParXim_pt_mult[ptBinXi][multBinXi])
                resultParXim_pt_mult[ptBinXi][multBinXi]->Write();

                out->cd("h_MassXip_pt_mult");
                // if (h_bgXip_pt_mult[ptBinXi][multBinXi])
                //     h_bgXip_pt_mult[ptBinXi][multBinXi]->Write();
                // if (h_MassXip_pt_mult[ptBinXi][multBinXi])
                h_MassXip_pt_mult[ptBinXi][multBinXi]->Write();
                // if (resultParXip_pt_mult[ptBinXi][multBinXi])
                resultParXip_pt_mult[ptBinXi][multBinXi]->Write();

                out->cd("h_MassXiC_pt_mult");
                // if (h_bgXiC_pt_mult[ptBinXi][multBinXi])
                //     h_bgXiC_pt_mult[ptBinXi][multBinXi]->Write();
                // if (h_MassXiC_pt_mult[ptBinXi][multBinXi])
                h_MassXiC_pt_mult[ptBinXi][multBinXi]->Write();
                // if (resultParXiC_pt_mult[ptBinXi][multBinXi])
                resultParXiC_pt_mult[ptBinXi][multBinXi]->Write();
                Info("Xi_diff", "Successfully written hists for pT bin = %d , Mult bin = %d ...", ptBinXi, multBinXi);

                // Printf(">>>>>>>>>>>>>>>>>>One cycle complete.......");
            }
        }
    }
    if (om)
    {
        /// pt+mult integrated case
        Info("Om_allInt", "Starting Omega pt+mult integrated analysis...");

        binsMC[0] = 1;
        binsMC[1] = nptbins_Om;
        binsMC[2] = 1;
        binsMC[3] = nmultbins_Om;
        if (isDoubleGausPol2)
        {
            f_Sig_Omp = new TF1("DblGausPol2_Omp", DoubleGausPol2, fitMinSig_Om, fitMaxSig_Om, 8);
            f_Sig_Omm = new TF1("DblGausPol2_Omm", DoubleGausPol2, fitMinSig_Om, fitMaxSig_Om, 8);
            f_Sig_OmC = new TF1("DblGausPol2_OmC", DoubleGausPol2, fitMinSig_Om, fitMaxSig_Om, 8);
            f_Sig_Omp->SetNpx(1000);
            f_Sig_Omm->SetNpx(1000);
            f_Sig_OmC->SetNpx(1000);

            f_Sig_Omp = setFitParametersDG(f_Sig_Omp, nullptr, h_MassOmp, lMass_Om, sigmaOm, kFALSE);
            f_Sig_Omm = setFitParametersDG(f_Sig_Omm, nullptr, h_MassOmm, lMass_Om, sigmaOm, kFALSE);
            f_Sig_OmC = setFitParametersDG(f_Sig_OmC, nullptr, h_MassOmC, lMass_Om, sigmaOm, kFALSE);
        }
        if (isGausPol2)
        {
            f_Sig_Omp = new TF1("GausPol2_Omp", GausPol2, fitMinSig_Om, fitMaxSig_Om, 6);
            f_Sig_Omm = new TF1("GausPol2_Omm", GausPol2, fitMinSig_Om, fitMaxSig_Om, 6);
            f_Sig_OmC = new TF1("GausPol2_OmC", GausPol2, fitMinSig_Om, fitMaxSig_Om, 6);
            f_Sig_Omp->SetNpx(1000);
            f_Sig_Omm->SetNpx(1000);
            f_Sig_OmC->SetNpx(1000);

            f_Sig_Omp = setFitParametersGaus(f_Sig_Omp, nullptr, h_MassOmp, lMass_Om, sigmaOm, kFALSE);
            f_Sig_Omm = setFitParametersGaus(f_Sig_Omm, nullptr, h_MassOmm, lMass_Om, sigmaOm, kFALSE);
            f_Sig_OmC = setFitParametersGaus(f_Sig_OmC, nullptr, h_MassOmC, lMass_Om, sigmaOm, kFALSE);
        }
        resultParams_Omp_allInt = fitResults(nullptr, f_Sig_Omp, h_MassOmp, fitMinSig_Om, fitMaxSig_Om, binsMC, MCGenOmp, fisMC, isGausPol2);
        resultParams_Omm_allInt = fitResults(nullptr, f_Sig_Omm, h_MassOmm, fitMinSig_Om, fitMaxSig_Om, binsMC, MCGenOmm, fisMC, isGausPol2);
        resultParams_OmC_allInt = fitResults(nullptr, f_Sig_OmC, h_MassOmC, fitMinSig_Om, fitMaxSig_Om, binsMC, MCGenOmC, fisMC, isGausPol2);

        out->cd("_allInt");
        h_MassOmm->Write();
        h_MassOmp->Write();
        h_MassOmC->Write();
        resultParams_Omp_allInt->SetNameTitle("resultParams_Omp_allInt", "Result Parameters #Omega^{+}");
        resultParams_Omm_allInt->SetNameTitle("resultParams_Omm_allInt", "Result Parameters #Omega^{-}");
        resultParams_OmC_allInt->SetNameTitle("resultParams_OmC_allInt", "Result Parameters #Omega^{+} + #Omega^{-}");
        resultParams_Omp_allInt->Write();
        resultParams_Omm_allInt->Write();
        resultParams_OmC_allInt->Write();

        Info("Om_allInt", "Finished writing allInt hists...");

        for (Int_t ptBinOm = 0; ptBinOm < nptbins_Om; ptBinOm++)
        {
            Info("Om_multInt", "Loop index: pT bin = %d ...", ptBinOm);

            if (fisMC)
            {
                // plotting mult integrated efficiency as a function of pT
                binsMC[0] = ptBinOm + 1;
                binsMC[1] = ptBinOm + 1;
                binsMC[2] = 1;
                binsMC[3] = nmultbins_Om;
            }

            h_MassOmm_pt[ptBinOm] = (TH1D *)Omm->ProjectionY(TString::Format(("h_MassOmm_pt[%d]"), ptBinOm), ptBinOm + 1, ptBinOm + 1, 1, nmultbins_Om, "e");
            h_MassOmm_pt[ptBinOm]->SetTitle(TString::Format(("Invariant Mass #Omega^{-}: p_{T}<%.1f,%.1f>, Mult<%.1f,%.1f>"), ptbins_Om[ptBinOm], ptbins_Om[ptBinOm + 1], multbins_Om[0], multbins_Om[nmultbins_Om]));

            h_MassOmp_pt[ptBinOm] = (TH1D *)Omp->ProjectionY(TString::Format(("h_MassOmp_pt[%d]"), ptBinOm), ptBinOm + 1, ptBinOm + 1, 1, nmultbins_Om, "e");
            h_MassOmp_pt[ptBinOm]->SetTitle(TString::Format(("Invariant Mass #Omega^{+}: p_{T}<%.1f,%.1f>, Mult<%.1f,%.1f>"), ptbins_Om[ptBinOm], ptbins_Om[ptBinOm + 1], multbins_Om[0], multbins_Om[nmultbins_Om]));

            h_MassOmC_pt[ptBinOm] = (TH1D *)h_MassOmp_pt[ptBinOm]->Clone(TString::Format(("h_MassOmC_pt[%d]"), ptBinOm));
            h_MassOmC_pt[ptBinOm]->Add(h_MassOmm_pt[ptBinOm]);
            h_MassOmC_pt[ptBinOm]->SetTitle(TString::Format(("Invariant Mass #Omega: p_{T}<%.1f,%.1f>, Mult<%.1f,%.1f>"), ptbins_Om[ptBinOm], ptbins_Om[ptBinOm + 1], multbins_Om[0], multbins_Om[nmultbins_Om]));

            // omega: merge 2 inv mass bins to 1 -> reducing total no. of bins from 100 to 50, increasing statistics
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
                    f_Sig_Omp_pt[ptBinOm] = setFitParametersDG(f_Sig_Omp_pt[ptBinOm], f_Sig_Omp->GetParameters(), h_MassOmp_pt[ptBinOm], lMass_Om, sigmaOm, kTRUE);
                    f_Sig_Omm_pt[ptBinOm] = setFitParametersDG(f_Sig_Omm_pt[ptBinOm], f_Sig_Omm->GetParameters(), h_MassOmm_pt[ptBinOm], lMass_Om, sigmaOm, kTRUE);
                    f_Sig_OmC_pt[ptBinOm] = setFitParametersDG(f_Sig_OmC_pt[ptBinOm], f_Sig_OmC->GetParameters(), h_MassOmC_pt[ptBinOm], lMass_Om, sigmaOm, kTRUE);
                }
                else
                {

                    f_Sig_Omp_pt[ptBinOm] = setFitParametersDG(f_Sig_Omp_pt[ptBinOm], f_Sig_Omp_pt[ptBinOm - 1]->GetParameters(), h_MassOmp_pt[ptBinOm], lMass_Om, sigmaOm, kTRUE);
                    f_Sig_Omm_pt[ptBinOm] = setFitParametersDG(f_Sig_Omm_pt[ptBinOm], f_Sig_Omm_pt[ptBinOm - 1]->GetParameters(), h_MassOmm_pt[ptBinOm], lMass_Om, sigmaOm, kTRUE);
                    f_Sig_OmC_pt[ptBinOm] = setFitParametersDG(f_Sig_OmC_pt[ptBinOm], f_Sig_OmC_pt[ptBinOm - 1]->GetParameters(), h_MassOmC_pt[ptBinOm], lMass_Om, sigmaOm, kTRUE);
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
                    f_Sig_Omp_pt[ptBinOm] = setFitParametersGaus(f_Sig_Omp_pt[ptBinOm], f_Sig_Omp->GetParameters(), h_MassOmp_pt[ptBinOm], lMass_Om, sigmaOm, kTRUE);
                    f_Sig_Omm_pt[ptBinOm] = setFitParametersGaus(f_Sig_Omm_pt[ptBinOm], f_Sig_Omm->GetParameters(), h_MassOmm_pt[ptBinOm], lMass_Om, sigmaOm, kTRUE);
                    f_Sig_OmC_pt[ptBinOm] = setFitParametersGaus(f_Sig_OmC_pt[ptBinOm], f_Sig_OmC->GetParameters(), h_MassOmC_pt[ptBinOm], lMass_Om, sigmaOm, kTRUE);
                }
                else
                {

                    f_Sig_Omp_pt[ptBinOm] = setFitParametersGaus(f_Sig_Omp_pt[ptBinOm], f_Sig_Omp_pt[ptBinOm - 1]->GetParameters(), h_MassOmp_pt[ptBinOm], lMass_Om, sigmaOm, kTRUE);
                    f_Sig_Omm_pt[ptBinOm] = setFitParametersGaus(f_Sig_Omm_pt[ptBinOm], f_Sig_Omm_pt[ptBinOm - 1]->GetParameters(), h_MassOmm_pt[ptBinOm], lMass_Om, sigmaOm, kTRUE);
                    f_Sig_OmC_pt[ptBinOm] = setFitParametersGaus(f_Sig_OmC_pt[ptBinOm], f_Sig_OmC_pt[ptBinOm - 1]->GetParameters(), h_MassOmC_pt[ptBinOm], lMass_Om, sigmaOm, kTRUE);
                }
            }

            // h_bgXip_pt[ptBinXi] = generateBg(bgTSpectrum, h_MassXip_pt[ptBinXi]);
            resultParOmm_pt[ptBinOm] = fitResults(nullptr, f_Sig_Omm_pt[ptBinOm], h_MassOmm_pt[ptBinOm], fitMinSig_Om, fitMaxSig_Om, binsMC, MCGenOmm, fisMC, isGausPol2);
            resultParOmm_pt[ptBinOm]->SetName(TString::Format(("resultParOmm_pt[%d]"), ptBinOm));
            resultParOmm_pt[ptBinOm]->SetTitle(TString::Format(("Result Parameters #Omega^{-}: p_{T}<%.1f,%.1f>, Mult<%.1f,%.1f>"), ptbins_Om[ptBinOm], ptbins_Om[ptBinOm + 1], multbins_Om[0], multbins_Om[nmultbins_Om]));

            // h_bgOmm_pt[ptBinOm] = generateBg(bgTSpectrum, h_MassOmm_pt[ptBinOm]);
            resultParOmp_pt[ptBinOm] = fitResults(nullptr, f_Sig_Omp_pt[ptBinOm], h_MassOmp_pt[ptBinOm], fitMinSig_Om, fitMaxSig_Om, binsMC, MCGenOmp, fisMC, isGausPol2);
            resultParOmp_pt[ptBinOm]->SetName(TString::Format(("resultParOmp_pt[%d]"), ptBinOm));
            resultParOmp_pt[ptBinOm]->SetTitle(TString::Format(("Result Parameters #Omega^{+}: p_{T}<%.1f,%.1f>, Mult<%.1f,%.1f>"), ptbins_Om[ptBinOm], ptbins_Om[ptBinOm + 1], multbins_Om[0], multbins_Om[nmultbins_Om]));

            // h_bgOmC_pt[ptBinOm] = generateBg(bgTSpectrum, h_MassOmC_pt[ptBinOm]);
            resultParOmC_pt[ptBinOm] = fitResults(nullptr, f_Sig_OmC_pt[ptBinOm], h_MassOmC_pt[ptBinOm], fitMinSig_Om, fitMaxSig_Om, binsMC, MCGenOmC, fisMC, isGausPol2);
            resultParOmC_pt[ptBinOm]->SetName(TString::Format(("resultParOmC_pt[%d]"), ptBinOm));
            resultParOmC_pt[ptBinOm]->SetTitle(TString::Format(("Result Parameters #Omega: p_{T}<%.1f,%.1f>, Mult<%.1f,%.1f>"), ptbins_Om[ptBinOm], ptbins_Om[ptBinOm + 1], multbins_Om[0], multbins_Om[nmultbins_Om]));

            out->cd("h_MassOmm_pt");
            h_MassOmm_pt[ptBinOm]->Write();
            // if (h_bgOmm_pt[ptBinOm])
            //     h_bgOmm_pt[ptBinOm]->Write();
            resultParOmm_pt[ptBinOm]->Write();

            out->cd("h_MassOmp_pt");
            h_MassOmp_pt[ptBinOm]->Write();
            // if (h_bgOmp_pt[ptBinOm])
            //     h_bgOmp_pt[ptBinOm]->Write();
            resultParOmp_pt[ptBinOm]->Write();

            out->cd("h_MassOmC_pt");
            h_MassOmC_pt[ptBinOm]->Write();
            // if (h_bgOmC_pt[ptBinOm])
            //     h_bgOmC_pt[ptBinOm]->Write();
            resultParOmC_pt[ptBinOm]->Write();

            Info("Om_multInt", "Successfully written hists for pT bin = %d ...", ptBinOm);
        }
        /// Starting fully differential case
        for (Int_t ptBinOm = 0; ptBinOm < nptbins_Om; ptBinOm++)
        {
            if (fisMC)
            {
                // plotting mult integrated efficiency as a function of pT
                binsMC[0] = ptBinOm + 1;
                binsMC[1] = ptBinOm + 1;
            }
            for (Int_t multBinOm = 0; multBinOm < nmultbins_Om; multBinOm++)
            {
                Info("Om_diff", "Loop index: pT bin = %d , Mult bin = %d ...", ptBinOm, multBinOm);

                /// fully differential case
                if (fisMC)
                {
                    binsMC[2] = multBinOm + 1;
                    binsMC[3] = multBinOm + 1;
                }

                h_MassOmm_pt_mult[ptBinOm][multBinOm] = (TH1D *)Omm->ProjectionY(TString::Format(("h_MassOmm_pt_mult[%d][%d]"), ptBinOm, multBinOm), ptBinOm + 1, ptBinOm + 1, multBinOm + 1, multBinOm + 1, "e");
                h_MassOmm_pt_mult[ptBinOm][multBinOm]->SetTitle(TString::Format(("Invariant Mass #Omega^{-}: p_{T}<%.1f,%.1f>, Mult<%.1f,%.1f>"), ptbins_Om[ptBinOm], ptbins_Om[ptBinOm + 1], multbins_Om[multBinOm], multbins_Om[multBinOm + 1]));

                h_MassOmp_pt_mult[ptBinOm][multBinOm] = (TH1D *)Omp->ProjectionY(TString::Format(("h_MassOmp_pt_mult[%d][%d]"), ptBinOm, multBinOm), ptBinOm + 1, ptBinOm + 1, multBinOm + 1, multBinOm + 1, "e");
                h_MassOmp_pt_mult[ptBinOm][multBinOm]->SetTitle(TString::Format(("Invariant Mass #Omega^{+}: p_{T}<%.1f,%.1f>, Mult<%.1f,%.1f>"), ptbins_Om[ptBinOm], ptbins_Om[ptBinOm + 1], multbins_Om[multBinOm], multbins_Om[multBinOm + 1]));

                h_MassOmC_pt_mult[ptBinOm][multBinOm] = (TH1D *)h_MassOmp_pt_mult[ptBinOm][multBinOm]->Clone(TString::Format(("h_MassOmC_pt_mult[%d][%d]"), ptBinOm, multBinOm));
                h_MassOmC_pt_mult[ptBinOm][multBinOm]->Add(h_MassOmm_pt_mult[ptBinOm][multBinOm]);
                h_MassOmC_pt_mult[ptBinOm][multBinOm]->SetTitle(TString::Format(("Invariant Mass #Omega: p_{T}<%.1f,%.1f>, Mult<%.1f,%.1f>"), ptbins_Om[ptBinOm], ptbins_Om[ptBinOm + 1], multbins_Om[multBinOm], multbins_Om[multBinOm + 1]));

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
                    f_Sig_Omp_pt_mult[ptBinOm][multBinOm] = setFitParametersDG(f_Sig_Omp_pt_mult[ptBinOm][multBinOm], f_Sig_Omp_pt[ptBinOm]->GetParameters(), h_MassOmp_pt_mult[ptBinOm][multBinOm], lMass_Om, resultParOmp_pt[ptBinOm]->GetBinContent(4), kTRUE);
                    f_Sig_Omm_pt_mult[ptBinOm][multBinOm] = setFitParametersDG(f_Sig_Omm_pt_mult[ptBinOm][multBinOm], f_Sig_Omm_pt[ptBinOm]->GetParameters(), h_MassOmm_pt_mult[ptBinOm][multBinOm], lMass_Om, resultParOmm_pt[ptBinOm]->GetBinContent(4), kTRUE);
                    f_Sig_OmC_pt_mult[ptBinOm][multBinOm] = setFitParametersDG(f_Sig_OmC_pt_mult[ptBinOm][multBinOm], f_Sig_OmC_pt[ptBinOm]->GetParameters(), h_MassOmC_pt_mult[ptBinOm][multBinOm], lMass_Om, resultParOmC_pt[ptBinOm]->GetBinContent(4), kTRUE);
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
                    //     f_Sig_Omp_pt_mult[ptBinOm][multBinOm] = setFitParametersDG(f_Sig_Omp_pt_mult[ptBinOm][multBinOm], f_Sig_Omp_pt_mult[ptBinOm - 1][multBinOm]->GetParameters(), h_MassOmp_pt_mult[ptBinOm][multBinOm], lMass_Om, sigmaOm, ptBinOm + 1, kTRUE);
                    //     f_Sig_Omm_pt_mult[ptBinOm][multBinOm] = setFitParametersDG(f_Sig_Omm_pt_mult[ptBinOm][multBinOm], f_Sig_Omm_pt_mult[ptBinOm - 1][multBinOm]->GetParameters(), h_MassOmm_pt_mult[ptBinOm][multBinOm], lMass_Om, sigmaOm, ptBinOm + 1, kTRUE);
                    //     f_Sig_OmC_pt_mult[ptBinOm][multBinOm] = setFitParametersDG(f_Sig_OmC_pt_mult[ptBinOm][multBinOm], f_Sig_OmC_pt_mult[ptBinOm - 1][multBinOm]->GetParameters(), h_MassOmp_pt_mult[ptBinOm][multBinOm], lMass_Om, sigmaOm, ptBinOm + 1, kTRUE);
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
                    f_Sig_Omp_pt_mult[ptBinOm][multBinOm] = setFitParametersGaus(f_Sig_Omp_pt_mult[ptBinOm][multBinOm], f_Sig_Omp_pt[ptBinOm]->GetParameters(), h_MassOmp_pt_mult[ptBinOm][multBinOm], lMass_Om, sigmaOm, kTRUE);
                    f_Sig_Omm_pt_mult[ptBinOm][multBinOm] = setFitParametersGaus(f_Sig_Omm_pt_mult[ptBinOm][multBinOm], f_Sig_Omm_pt[ptBinOm]->GetParameters(), h_MassOmm_pt_mult[ptBinOm][multBinOm], lMass_Om, sigmaOm, kTRUE);
                    f_Sig_OmC_pt_mult[ptBinOm][multBinOm] = setFitParametersGaus(f_Sig_OmC_pt_mult[ptBinOm][multBinOm], f_Sig_OmC_pt[ptBinOm]->GetParameters(), h_MassOmC_pt_mult[ptBinOm][multBinOm], lMass_Om, sigmaOm, kTRUE);
                    // }
                    // else
                    // {

                    //     f_Sig_Omp_pt_mult[ptBinOm][multBinOm] = setFitParametersGaus(f_Sig_Omp_pt_mult[ptBinOm][multBinOm], f_Sig_Omp_pt_mult[ptBinOm - 1][multBinOm]->GetParameters(), h_MassOmp_pt_mult[ptBinOm][multBinOm], lMass_Om, sigmaOm, kTRUE);
                    //     f_Sig_Omm_pt_mult[ptBinOm][multBinOm] = setFitParametersGaus(f_Sig_Omm_pt_mult[ptBinOm][multBinOm], f_Sig_Omm_pt_mult[ptBinOm - 1][multBinOm]->GetParameters(), h_MassOmm_pt_mult[ptBinOm][multBinOm], lMass_Om, sigmaOm, kTRUE);
                    //     f_Sig_OmC_pt_mult[ptBinOm][multBinOm] = setFitParametersGaus(f_Sig_OmC_pt_mult[ptBinOm][multBinOm], f_Sig_OmC_pt_mult[ptBinOm - 1][multBinOm]->GetParameters(), h_MassOmp_pt_mult[ptBinOm][multBinOm], lMass_Om, sigmaOm, kTRUE);
                    // }
                }

                resultParOmp_pt_mult[ptBinOm][multBinOm] = fitResults(nullptr, f_Sig_Omp_pt_mult[ptBinOm][multBinOm], h_MassOmp_pt_mult[ptBinOm][multBinOm], fitMinSig_Om, fitMaxSig_Om, binsMC, MCGenOmp, fisMC, isGausPol2);
                resultParOmp_pt_mult[ptBinOm][multBinOm]->SetName(TString::Format(("resultParOmp_pt_mult[%d][%d]"), ptBinOm, multBinOm));
                resultParOmp_pt_mult[ptBinOm][multBinOm]->SetTitle(TString::Format(("Result Parameters #Omega^{+}: p_{T}<%.1f,%.1f>, Mult<%.1f,%.1f>"), ptbins_Om[ptBinOm], ptbins_Om[ptBinOm + 1], multbins_Om[multBinOm], multbins_Om[multBinOm + 1]));

                // peakFnc_OmM = setFitParameters(peakFnc_OmM, h_MassOmm_pt_mult[ptBinOm][multBinOm], lMass_Om, ptBinOm + 1, kFALSE);
                // h_bgOmm_pt_mult[ptBinOm][multBinOm] = generateBg(bgTSpectrum, h_MassOmm_pt_mult[ptBinOm][multBinOm]);
                // resultParOmm_pt_mult[ptBinOm][multBinOm] = fitResults(h_bgOmm_pt_mult[ptBinOm][multBinOm], peakFnc_OmM, h_MassOmm_pt_mult[ptBinOm][multBinOm], fitMinSig_Om, fitMaxSig_Om, &par_allintM[0], &par_allint_errorsM[0], binsMC, MCGenOmm, fisMC, kFALSE);
                resultParOmm_pt_mult[ptBinOm][multBinOm] = fitResults(nullptr, f_Sig_Omm_pt_mult[ptBinOm][multBinOm], h_MassOmm_pt_mult[ptBinOm][multBinOm], fitMinSig_Om, fitMaxSig_Om, binsMC, MCGenOmm, fisMC, isGausPol2);
                resultParOmm_pt_mult[ptBinOm][multBinOm]->SetName(TString::Format(("resultParOmm_pt_mult[%d][%d]"), ptBinOm, multBinOm));
                resultParOmm_pt_mult[ptBinOm][multBinOm]->SetTitle(TString::Format(("Result Parameters #Omega^{-}: p_{T}<%.1f,%.1f>, Mult<%.1f,%.1f>"), ptbins_Om[ptBinOm], ptbins_Om[ptBinOm + 1], multbins_Om[multBinOm], multbins_Om[multBinOm + 1]));

                // peakFnc_OmC = setFitParameters(peakFnc_OmC, h_MassOmC_pt_mult[ptBinOm][multBinOm], lMass_Om, ptBinOm + 1, kFALSE);
                // h_bgOmC_pt_mult[ptBinOm][multBinOm] = generateBg(bgTSpectrum, h_MassOmC_pt_mult[ptBinOm][multBinOm]);
                // resultParOmC_pt_mult[ptBinOm][multBinOm] = fitResults(h_bgOmC_pt_mult[ptBinOm][multBinOm], peakFnc_OmC, h_MassOmC_pt_mult[ptBinOm][multBinOm], fitMinSig_Om, fitMaxSig_Om, &par_allintC[0], &par_allint_errorsC[0], binsMC, MCGenOmC, fisMC, kFALSE);
                resultParOmC_pt_mult[ptBinOm][multBinOm] = fitResults(nullptr, f_Sig_OmC_pt_mult[ptBinOm][multBinOm], h_MassOmC_pt_mult[ptBinOm][multBinOm], fitMinSig_Om, fitMaxSig_Om, binsMC, MCGenOmC, fisMC, isGausPol2);
                resultParOmC_pt_mult[ptBinOm][multBinOm]->SetName(TString::Format(("resultParOmC_pt_mult[%d][%d]"), ptBinOm, multBinOm));
                resultParOmC_pt_mult[ptBinOm][multBinOm]->SetTitle(TString::Format(("Result Parameters #Omega: p_{T}<%.1f,%.1f>, Mult<%.1f,%.1f>"), ptbins_Om[ptBinOm], ptbins_Om[ptBinOm + 1], multbins_Om[multBinOm], multbins_Om[multBinOm + 1]));

                out->cd("h_MassOmm_pt_mult");
                h_MassOmm_pt_mult[ptBinOm][multBinOm]->Write();
                // if (h_bgOmm_pt_mult[ptBinOm][multBinOm])
                //     h_bgOmm_pt_mult[ptBinOm][multBinOm]->Write();
                resultParOmm_pt_mult[ptBinOm][multBinOm]->Write();

                out->cd("h_MassOmp_pt_mult");
                h_MassOmp_pt_mult[ptBinOm][multBinOm]->Write();
                // if (h_bgOmp_pt_mult[ptBinOm][multBinOm])
                //     h_bgOmp_pt_mult[ptBinOm][multBinOm]->Write();
                resultParOmp_pt_mult[ptBinOm][multBinOm]->Write();

                out->cd("h_MassOmC_pt_mult");
                h_MassOmC_pt_mult[ptBinOm][multBinOm]->Write();
                // if (h_bgOmC_pt_mult[ptBinOm][multBinOm])
                //     h_bgOmC_pt_mult[ptBinOm][multBinOm]->Write();
                resultParOmC_pt_mult[ptBinOm][multBinOm]->Write();
                Info("Om_diff", "Successfully written hists for pT bin = %d , Mult bin = %d ...", ptBinOm, multBinOm);
            }
        }
    }

    out->Close();
    f->Close();
    Printf("FINISHED!!! Output saved to '%s'", output.Data());

    // delete f, out;
    return 0;
}

TH1 *generateBg(bool bgTSpectrum, TH1 *h_source)
{
    if (!bgTSpectrum)
    {
        return nullptr;
    }

    /// testing spectrum
    Int_t nbinsX = h_source->GetNbinsX();
    Double_t source[nbinsX];
    for (int itest = 0; itest < nbinsX; itest++)
        source[itest] = h_source->GetBinContent(itest + 1);

    TSpectrum *s = new TSpectrum(1);
    // int numberIterations[] = {6, 6, 6, 6};
    // int clipWindow[] = {1, 1};
    // int backOrder[] = {2, 2, 2, 2};
    // int backSmoothing[] = {15, 15, 15, 15, 15, 15, 15};
    for (int i = 0; i < 4; i++)
    {
        // s_omC->Background(source, nbinsX, numberIterations[i], clipWindow[i], backOrder[i], kTRUE, backSmoothing[i], kTRUE);
        // TString hName = TString::Format(("bg_i%d_c%d_o%d_S%d_[%d][%d]"), numberIterations[i], clipWindow[i], backOrder[i], backSmoothing[i], ptBinOm, multBinOm);

        s->Background(source, nbinsX, 20, TSpectrum::kBackDecreasingWindow, TSpectrum::kBackOrder6, kTRUE, TSpectrum::kBackSmoothing15, kFALSE);

        // bg[ptBinOm][multBinOm]->GetXaxis()->SetRange(h_MassOmC_pt_mult[ptBinOm][multBinOm]->GetXaxis()->GetXmin(), h_MassOmm_pt_mult[ptBinOm][multBinOm]->GetXaxis()->GetXmax());
        // Draw the estimated background
    }
    // TString hName = TString::Format(("bg_omC[%d][%d]"), ptBinOm, multBinOm);
    // h_bgOmC_pt_mult[ptBinOm][multBinOm] = new TH1D(hName.Data(), hName.Data(), nbinsX, h_MassOmC_pt_mult[ptBinOm][multBinOm]->GetXaxis()->GetXmin(), h_MassOmC_pt_mult[ptBinOm][multBinOm]->GetXaxis()->GetXmax());
    TH1D *h_bg = new TH1D(TString::Format(("h_bg%s"), h_source->GetName() + 6), TString::Format(("Estimated Background%s"), h_source->GetTitle() + 14), nbinsX, h_source->GetXaxis()->GetXmin(), h_source->GetXaxis()->GetXmax());
    h_bg->SetLineColor(kViolet);

    for (int itest = 0; itest < nbinsX; itest++)
        h_bg->SetBinContent(itest + 1, source[itest]);

    return h_bg;
}

TF1 *setFitParametersGaus(TF1 *peakFnc, Double_t *previousPtBinFitParams, TH1 *peak, Double_t mass, Double_t sigma, Bool_t usePrevious)
{
    if (usePrevious)
    {
        mass = previousPtBinFitParams[1];
        sigma = previousPtBinFitParams[2];
    }

    // TString optBg = "QMNR MULTITHREAD";
    TString optBg = "QNFB MULTITHREAD";

    Double_t fourSigma = 4 * sigma;

    // Double_t zeroes[9] = {0.};
    // peakFnc->SetParameters(zeroes);
    // peakFnc->SetParErrors(zeroes);
    // for (Int_t i = 0; i < 9; i++)
    // {
    //     peakFnc->ReleaseParameter(i);
    // }
    /// gaus fit
    TF1 *bgFnc = new TF1("bgFnc_Pol2", Pol2Exclude, peak->GetXaxis()->GetBinCenter(peak->FindFirstBinAbove(0.)), peak->GetXaxis()->GetBinCenter(peak->FindLastBinAbove(0.)), 6);
    // TString bgFncName = peakFnc->GetName();
    // bgFncName.ReplaceAll("DblGausPol2", "_DG_bgPol2");

    // TF1 *bgFnc = new TF1("bgFnc_Pol2", Pol2Exclude, mass - 15 * sigma, mass + 15 * sigma, 6);
    bgFnc->SetParNames("p[0]", "p[1]", "p[2]", "peakFitMass", "peakFitSigma", "rejectBgPoint");
    bgFnc->SetParameter(2, peak->GetMaximum() * 0.1);
    bgFnc->FixParameter(3, mass);      // "mass = mean of fit"
    bgFnc->FixParameter(4, fourSigma); // "4*sigma"
    // bgFnc->SetLineColor(kOrange);
    // bgFnc->SetFillColor(kBlue);
    // bgFnc->SetFillStyle(3016);
    bgFnc->FixParameter(5, 1); // "reject" (set 1 or 0)
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
    // peakFnc->SetParLimits(4, mass - 2 * sigma, mass + 2 * sigma);
    // peakFnc->SetParLimits(3, 0, peak->GetMaximum() * 1.1);
    peakFnc->SetParLimits(2, 0.0009, 0.01);
    // peakFnc->SetParLimits(5, 0.001, 0.01);

    if (usePrevious)
    {
        if (previousPtBinFitParams)
        {
            peakFnc->SetParameter(0, previousPtBinFitParams[0]);
            peakFnc->SetParameter(1, previousPtBinFitParams[1]);
            // peakFnc->SetParameter(4, previousPtBinFitParams[4]);
            // peakFnc->SetParameter(5, previousPtBinFitParams[5]);
            peakFnc->SetParameter(2, previousPtBinFitParams[2]);
            // peakFnc->SetParameter(3, previousPtBinFitParams[3]);
        }
        else
        {
            Error("setFitParametersGaus", "previousPtBinFitParams not passed for %s. Setting defaults...", peak->GetName());
            setFitParametersGaus(peakFnc, nullptr, peak, mass, sigma, kFALSE);
        }
    }
    else
    {
        peakFnc->SetParameter(0, peak->GetMaximum() * 0.9);
        peakFnc->SetParameter(1, mass);
        // peakFnc->SetParameter(4, mass);
        // peakFnc->SetParameter(5, 0.001);
        peakFnc->SetParameter(2, sigma);
        // peakFnc->SetParameter(3, peak->GetMaximum() * 0.3);
    }

    return peakFnc;
}

TF1 *setFitParametersDG(TF1 *peakFnc, Double_t *previousPtBinFitParams, TH1 *peak, Double_t mass, Double_t sigma, Bool_t usePrevious)
{
    // if (ptBin > 1)
    // {
    //     mass = previousPtBinFitParams[1];
    //     sigma = previousPtBinFitParams[2];
    // }

    // TString optBg = "QMNR MULTITHREAD";
    TString optBg = "QNFB MULTITHREAD";

    Double_t fourSigma = 4 * sigma;

    // Double_t zeroes[9] = {0.};
    // peakFnc->SetParameters(zeroes);
    // peakFnc->SetParErrors(zeroes);
    // for (Int_t i = 0; i < 9; i++)
    // {
    //     peakFnc->ReleaseParameter(i);
    // }
    /// double gaus fit
    TF1 *bgFnc = new TF1("bgFnc_Pol2", Pol2Exclude, peak->GetXaxis()->GetBinCenter(peak->FindFirstBinAbove(0.)), peak->GetXaxis()->GetBinCenter(peak->FindLastBinAbove(0.)), 6);
    // TString bgFncName = peakFnc->GetName();
    // bgFncName.ReplaceAll("DblGausPol2", "_DG_bgPol2");

    // TF1 *bgFnc = new TF1("bgFnc_Pol2", Pol2Exclude, mass - 15 * sigma, mass + 15 * sigma, 6);
    bgFnc->SetParNames("p[0]", "p[1]", "p[2]", "peakFitMass", "peakFitSigma", "rejectBgPoint");
    bgFnc->SetParameter(2, peak->GetMaximum() * 0.1);
    bgFnc->FixParameter(3, mass);      // "mass = mean of fit"
    bgFnc->FixParameter(4, fourSigma); // "4*sigma"
    // bgFnc->SetLineColor(kOrange);
    // bgFnc->SetFillColor(kBlue);
    // bgFnc->SetFillStyle(3016);
    bgFnc->FixParameter(5, 1); // "reject" (set 1 or 0)
    for (int iFitBg = 0; iFitBg < 3; iFitBg++)
        peak->Fit(bgFnc, optBg.Data(), "", peak->GetXaxis()->GetBinCenter(peak->FindFirstBinAbove(0.)), peak->GetXaxis()->GetBinCenter(peak->FindLastBinAbove(0.)));
    // peak->Fit(bgFnc, optBg.Data(), "", mass - 15 * sigma, mass + 15 * sigma);
    bgFnc->FixParameter(5, 0); // reject = off (set 1 or 0)
    // peak->GetListOfFunctions()->Add(bgFnc);

    peakFnc->SetParNames("Amp_{1}", "#mu_{1}", "#sigma_{1}", "Amp_{2}", "#sigma_{2}");
    peakFnc->FixParameter(5, bgFnc->GetParameter(2));
    peakFnc->FixParameter(6, bgFnc->GetParameter(1));
    peakFnc->FixParameter(7, bgFnc->GetParameter(0));

    /// set params for mult integrated functions -> these params will be reused for mult differential too
    peakFnc->SetParLimits(0, 0, peak->GetMaximum() * 1.1);
    peakFnc->SetParLimits(1, mass - 2 * sigma, mass + 2 * sigma);
    // peakFnc->FixParameter(1, mass);
    peakFnc->SetParLimits(2, 0.0009, 0.009);
    peakFnc->SetParLimits(3, 0, peak->GetMaximum() * 0.6);
    // peakFnc->SetParLimits(4, mass - 2 * sigma, mass + 2 * sigma);
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
            Error("setFitParametersDG", "previousPtBinFitParams not passed for %s. Setting defaults...", peak->GetName());
            setFitParametersDG(peakFnc, nullptr, peak, mass, sigma, kFALSE);
        }
    }
    else
    {
        peakFnc->SetParameter(1, mass);
        peakFnc->SetParameter(2, sigma);
        // peakFnc->SetParameter(4, mass);
        peakFnc->SetParameter(4, 0.002);
    }
    peakFnc->Update();
    return peakFnc;
}

TH1 *fitResults(TH1 *h_bg, TF1 *peakFnc, TH1 *peak, Double_t fitMinSig, Double_t fitMaxSig, Int_t binsMC[], TH2 *MCgen, bool fisMC, Bool_t isGausPol2)
{
    // ROOT::Math::IntegratorOneDimOptions::SetDefaultAbsTolerance(1.E-6);
    // ROOT::Math::IntegratorOneDimOptions::SetDefaultRelTolerance(1.E-6);

    Double_t peakFitMass = 0.;
    Double_t peakFitSigma = 0.;
    Double_t peakFitMass_err, peakFitSigma_err = 0.;

    Double_t bgReject = 0.;

    if (TMath::IsNaN(peak->GetEntries()))
    {
        return 0;
    }

    /// use different fit options:
    TString opt = "QLINES+ MULTITHREAD"; /// best for DG: QLINESR+ MULTITHREAD
    // TString opt = "QMINLSR";

    // TString optBg = "QMNLSR";
    TString optBg = "QNSFB+ MULTITHREAD"; /// best for DG bkg: QLINSRFB+ MULTITHREAD

    for (int i = 0; i < 10; i++)
        peak->Fit(peakFnc, "QLIN MULTITHREAD", "", fitMinSig, fitMaxSig);

    TFitResultPtr fFitResult;
    fFitResult = peak->Fit(peakFnc, opt.Data(), "", fitMinSig, fitMaxSig);

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

    if (!fFitResult->IsValid())
    {
        // Printf(">>>>>>>>>>>>%s: Peak fit not converged. Fit status: %d, Minimise: %d, Minos: %d, Hesse: %d, Avg mass = %f+/-%f, sigma = %f+/-%f<<<<<<<<<<<<<<<<", peak->GetName(), fitStatus, ROOT::Minuit2::Minuit2Minimizer::Minimize(), ROOT::Minuit2::Minuit2Minimizer::MinosStatus(), ROOT::Minuit2::Minuit2Minimizer::Hesse(), peakFitMass, peakFitMass_err, peakFitSigma, peakFitSigma_err);
        Error("fitResult", "%s: Peak fit not converged. Fit status: %d, Avg mass = %f+/-%f, sigma = %f+/-%f", peak->GetName(), fFitResult->Status(), peakFitMass, peakFitMass_err, peakFitSigma, peakFitSigma_err);
        // fFitResult->Print();
        //  fFitResult_bg->Print();
        // fFitResult->GetCovarianceMatrix().Print();
    }
    // Double_t bgFitMin = peakFitMass - 0.03; // Mean - 15*sigma (approx.)
    // Double_t bgFitMax = peakFitMass + 0.03; // Mean + 15*sigma (approx.)
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
        bgFnc->SetParNames("p[0]", "p[1]", "p[2]", "peakFitMass", "peakFitSigma", "rejectBgPoint");
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
            bgFnc->SetParError(0, peakFnc->GetParError(3));
            bgFnc->SetParError(1, peakFnc->GetParError(4));
            bgFnc->SetParError(2, peakFnc->GetParError(5));
        }
        else
        {
            bgFnc->SetParameter(0, peakFnc->GetParameter(7));
            bgFnc->SetParameter(1, peakFnc->GetParameter(6));
            bgFnc->SetParameter(2, peakFnc->GetParameter(5));
            bgFnc->SetParError(0, peakFnc->GetParError(7));
            bgFnc->SetParError(1, peakFnc->GetParError(6));
            bgFnc->SetParError(2, peakFnc->GetParError(5));
        }

        // Double_t fitMinBg = peakFitMass - 10 * peakFitSigma;
        // Double_t fitMaxBg = peakFitMass + 10 * peakFitSigma;
        bgFnc->FixParameter(5, 1); // "reject" (set 1 or 0)
        // fFitResult_bg = peak->Fit(bgFnc, optBg.Data(), "", fitMinBg, fitMaxBg);
        for (int i = 0; i < 10; i++)
            peak->Fit(bgFnc, "QNFB MULTITHREAD", "", bgFitMin, bgFitMax);
        // TMinuit::SetPrintLevel(1);
        fFitResult_bg = peak->Fit(bgFnc, optBg.Data(), "", bgFitMin, bgFitMax);
        bgFnc->FixParameter(5, 0); // reject = off (set 1 or 0)
                                   // Printf("bg fit done!");
        // fitStatus = fFitResult_bg;
        if (!fFitResult_bg->IsValid())
        {
            Error("fitResult_bg", "%s: Bg fit not converged. Fit status: %d", peak->GetName(), fFitResult_bg->Status());
            if (gErrorIgnoreLevel < kFatal)
                fFitResult_bg->GetCovarianceMatrix().Print();
        }

        if (fFitResult_bg->CovMatrixStatus() < 2)
        {
            Warning("fitResult_bg", "%s: Covariance Matrix not accurate. Status: %d", peak->GetName(), fFitResult_bg->CovMatrixStatus());
        }

        // add fit and background function to histogram so it is automatically drawn
        // with hist
        peak->GetListOfFunctions()->Add(bgFnc);
    }

    /// getting min/max for deciding signal's fit range for integral calculation
    Double_t minInt = peakFitMass - bgReject; /// mean - 4*sigma
    Double_t maxInt = peakFitMass + bgReject; /// mean + 4*sigma
    Info("fitResults", "%s: MinInt, MaxInt = %f, %f. Filling parameters ...", peak->GetName(), minInt, maxInt);

    TH1 *resultParams = fillParams(h_bg, fFitResult_bg, peakFnc, peak, fFitResult, minInt, maxInt, MCgen, binsMC, fisMC, isGausPol2);
    peak->SetOption("X0E1");
    return resultParams;
}

TH1 *fillParams(TH1 *h_bg, TFitResultPtr fFitResult_bg, TF1 *sigBgFnc, TH1 *peak, TFitResultPtr fFitResult, Double_t minInt, Double_t maxInt, TH2 *MCgen, Int_t binsMC[], bool fisMC, bool isGausPol2)
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
        Error("fillParams", "%s: GetYieldBinCounting is NaN!", peak->GetName());
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
        Error("fillParams", "%s: GetYieldFitFunction is NaN!", peak->GetName());
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
        Info("fillParams: fisMC", "%s: IntGen = %f +/- %f", peak->GetName(), gen, genErr);

        if (gen != 0 && genErr != 0)
        {
            eff = val_intBC / gen;
            // effErr = err_intBC / genErr; /// ask if it should be computed like this?

            effErr = ErrorInRatio(val_intBC, err_intBC, gen, genErr);
            Info("fillParams: fisMC", "%s: Efficiency = %f; Error = %f", peak->GetName(), eff, effErr);
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

    for (Int_t i = 0; i < sigBgFnc->GetNpar(); i++)
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

    // Info("fillParams", "Parameters filled. Returning %s ...", resultParams->GetName());
    return resultParams;
}

void GetYieldBinCounting(TH1 *h_bg, TH1 *h, TFitResultPtr fFitResult_bg, Double_t minInt, Double_t maxInt, bool fisMC, Double_t &val,
                         Double_t &err, Double_t eps = 1e-3)
{

    if (!h)
        return;

    Int_t bin_minSig = h->GetXaxis()->FindBin(minInt);
    Int_t bin_maxSig = h->GetXaxis()->FindBin(maxInt);

    /// Histogram integrals always use whole bins. So, for comparisions with function integrals:
    Double_t fnIntMin = h->GetBinLowEdge(bin_minSig);     // left edge of the "bin_minSig" bin
    Double_t fnIntMax = h->GetBinLowEdge(bin_maxSig + 1); // right edge of the "bin_maxSig" bin
    Double_t histWidth = h->GetXaxis()->GetBinWidth(1);

    val = h->IntegralAndError(bin_minSig, bin_maxSig, err);

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
            // Printf("Bgfnc Param: [0] = %f , [1] = %f , [2] = %f, [3] = %f , [4] = %f , [5] = %f", bgFnc->GetParameter(0), bgFnc->GetParameter(1), bgFnc->GetParameter(2), bgFnc->GetParameter(3), bgFnc->GetParameter(4), bgFnc->GetParameter(5));
            // Printf("FitRs Param: [0] = %f , [1] = %f , [2] = %f, [3] = %f , [4] = %f , [5] = %f", params[0], params[1], params[2], params[3], params[4], params[5]);
            // Printf("BgPar error: [0] = %f , [1] = %f , [2] = %f, [3] = %f , [4] = %f , [5] = %f", bgFnc->GetParError(0), bgFnc->GetParError(1), bgFnc->GetParError(2), bgFnc->GetParError(3), bgFnc->GetParError(4), bgFnc->GetParError(5));
            bgVal = bgFnc->Integral(fnIntMin, fnIntMax, eps);
            bgErr = bgFnc->IntegralError(fnIntMin, fnIntMax, bgFnc->GetParameters(), fFitResult_bg->GetCovarianceMatrix().GetMatrixArray(), eps);
            bgVal /= histWidth;
            bgErr /= histWidth;
            if (!(bgVal && bgErr && val && err))
                Error("GetYieldBinCounting", "IntBC: %s: [%f, %f], bins:[%d, %d]: bg(%s) = %f +/- %f, sig(%s) = %f +/- %f", h->GetName(), fnIntMin, fnIntMax, bin_minSig, bin_maxSig, bgFnc->GetName(), bgVal, bgErr, h->GetListOfFunctions()->At(0)->GetName(), val, err);
            else
                Info("GetYieldBinCounting", "IntBC: %s: [%f, %f], bins:[%d, %d]: bg(%s) = %f +/- %f, sig(%s) = %f +/- %f", h->GetName(), fnIntMin, fnIntMax, bin_minSig, bin_maxSig, bgFnc->GetName(), bgVal, bgErr, h->GetListOfFunctions()->At(0)->GetName(), val, err);
        }
        else
        {
            Double_t bgWidth = h_bg->GetXaxis()->GetBinWidth(1);

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
                         Double_t &err, Double_t eps = 1e-3)
{

    if (!h)
        return;

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

void GetMeanSigmaDG(TF1 *f_doubleGaus, TFitResultPtr lFitResultPtr, Double_t &mean, Double_t &mean_err, Double_t &sigma, Double_t &sigma_err)
{
    //[6]+[7]*x+[8]*x*x+[0]*TMath::Gaus(x,[1],[2])+[3]*TMath::Gaus(x,[4],[5])",1.305,1.34
    Double_t N1a = f_doubleGaus->GetParameter(0);         // N1
    Double_t N2a = f_doubleGaus->GetParameter(3);         // N2
    Double_t mu1a = f_doubleGaus->GetParameter(1);        // mu1
    Double_t mu2a = f_doubleGaus->GetParameter(1);        // mu2
    Double_t sigma1a = f_doubleGaus->GetParameter(2);     // sigma1
    Double_t sigma2a = f_doubleGaus->GetParameter(4);     // sigma2
    TMatrixD cova = lFitResultPtr->GetCovarianceMatrix(); // lFitResultPtr being the FitResult ptr

    Double_t mu_wa = (N1a * mu1a + N2a * mu2a) / (N1a + N2a);
    // cout << "func   " << mu1a << "  " << mu2a << endl;
    // Double_t mu_wa = (mu1a + mu2a)/2;
    Double_t sigma_wa = (N1a * sigma1a + N2a * sigma2a) / (N1a + N2a);

    Double_t sa = N1a + N2a;
    Double_t wa_mu = N1a * mu1a + N2a * mu2a;
    Double_t wa_sigma = N1a * sigma1a + N2a * sigma2a;

    Double_t mu_wa_step = pow((mu1a - mu2a), 2) * (pow(N1a, 2) * cova(3, 3) + pow(N2a, 2) * cova(0, 0)) + 2 * cova(0, 3) * (wa_mu - sa * mu1a) * (wa_mu - sa * mu2a) + pow(sa, 2) * (pow(N1a, 2) * cova(1, 1) + pow(N2a, 2) * cova(4, 4) + 2 * N1a * N2a * cova(1, 4)) - 2 * sa * (N1a * (cova(0, 1) * (wa_mu - sa * mu1a) + cova(3, 1) * (wa_mu - sa * mu2a)) + N2a * (cova(0, 4) * (wa_mu - sa * mu1a) + cova(3, 4) * (wa_mu - sa * mu2a)));
    Double_t mu_wa_err = sqrt(mu_wa_step / pow(sa, 4));

    Double_t sigma_wa_step = pow((sigma1a - sigma2a), 2) * (pow(N1a, 2) * cova(3, 3) + pow(N2a, 2) * cova(0, 0)) + 2 * cova(0, 3) * (wa_sigma - sa * sigma1a) * (wa_sigma - sa * sigma2a) + pow(sa, 2) * (pow(N1a, 2) * cova(2, 2) + pow(N2a, 2) * cova(5, 5) + 2 * N1a * N2a * cova(2, 5)) - 2 * sa * (N1a * (cova(0, 2) * (wa_sigma - sa * sigma1a) + cova(3, 2) * (wa_sigma - sa * sigma2a)) + N2a * (cova(0, 5) * (wa_sigma - sa * sigma1a) + cova(3, 5) * (wa_sigma - sa * sigma2a)));
    Double_t sigma_wa_err = sqrt(sigma_wa_step / pow(sa, 4));

    mean = mu_wa;
    mean_err = mu_wa_err;
    sigma = sigma_wa;
    sigma_err = sigma_wa_err;
    Info("GetMeanSigmaDG", "Avg Fit Mass = %f +/- %f; Avg Fit Sigma = %f +/- %f", mean, mean_err, sigma, sigma_err);
    // Info("GetMeanSigmaDG", "Avg Fit Sigma = %f +/- %f", sigma, sigma_err);
}
