#include <TString.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TH1.h>
#include <TF1.h>
#include <TFitResult.h>
#include <TFitResultPtr.h>
#include "TLine.h"

double bgReject_Xi = 0.0075; // 4*sigma=0.01, 6*sigma=0.015
double bgReject_Om = 0.0072; // 4*sigma=0.01, 6*sigma=0.015
// double bgReject_Xi = 0.01125; //4*sigma=0.01, 6*sigma=0.015
// double bgReject_Om = 0.0108;  //4*sigma=0.01, 6*sigma=0.015
bool reject; /// rejecting the region +/-4sigma for fitting background

double mean_Xi[282] = {0};
double mean_Om[62] = {0};
double fourSigma_Xi[282] = {0};
double fourSigma_Om[62] = {0};
int count_Xi = 0;
int count_Om = 0;
// Double_t genErr = 0;

Double_t Pol1(double *x, double *par)
{
    return par[0] * x[0] + par[1];
}
Double_t Pol1Xi(double *x, double *par)
{
    if (reject && x[0] > (1.32171 - bgReject_Xi) && x[0] < (1.32171 + bgReject_Xi))
    {
        TF1::RejectPoint();
        return 0;
    }
    return par[0] * x[0] + par[1];
}
Double_t Pol1Om(double *x, double *par)
{
    if (reject && x[0] > (1.67245 - bgReject_Om) && x[0] < (1.67245 + bgReject_Om))
    {
        TF1::RejectPoint();
        return 0;
    }
    return par[0] * x[0] + par[1];
}
Double_t Pol2(double *x, double *par)
{
    return par[0] + x[0] * par[1] + x[0] * x[0] * par[2];
}
Double_t Pol2Xi(double *x, double *par)
{
    if (reject && x[0] > (1.32171 - bgReject_Xi) && x[0] < (1.32171 + bgReject_Xi))
    {
        TF1::RejectPoint();
        return 0;
    }
    return par[0] + x[0] * par[1] + x[0] * x[0] * par[2];
}
Double_t Pol2Om(double *x, double *par)
{
    if (reject && x[0] > (1.67245 - bgReject_Om) && x[0] < (1.67245 + bgReject_Om))
    {
        TF1::RejectPoint();
        return 0;
    }
    return par[0] + x[0] * par[1] + x[0] * x[0] * par[2];
}
Double_t GausPol2(double *x, double *par)
{
    return par[0] * TMath::Gaus(x[0], par[1], par[2]) + Pol2(x, &par[3]);
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

TF1 *setFitParameters(TF1 *peakFnc, TH1 *peak, Double_t mass, Int_t ptBin);

TH1 *fitResults(TF1 *peakFnc, TF1 *bgFnc, TH1 *peak, Double_t fitMinSig, Double_t fitMaxSig, Double_t fitMinBg, Double_t fitMaxBg, TH2 *MCgen, Int_t binsMC[], bool options = kFALSE);

TH1 *fillParams(TF1 *sigBgFnc, TH1 *peak, TFitResultPtr fFitResult, Double_t minInt, Double_t maxInt, TH2 *MCgen, Int_t binsMC[]);

void GetYieldBinCounting(TH1 *h, TFitResultPtr fFitResult, Double_t minInt, Double_t maxInt, Double_t &val,
                         Double_t &err, Double_t eps = 1e-6);
void GetYieldFitFunction(TH1 *h, TFitResultPtr fFitResult, Double_t minInt, Double_t maxInt, Double_t &val,
                         Double_t &err, Double_t eps = 1e-6);

const bool fisMC = true;

int AliCascadeAnalysisMC_Ishaan_Fitting(TString input = "AnalysisResults.root", TString idAxis = "y",
                                        TString outputFilename = "AliCascadeAnalysisMC_Ishaan_Fitting.root", bool xi = true, bool om = true)
{
    TH1::AddDirectory(0);

    TFile *f = TFile::Open(input);
    if (!f)
    {
        Printf("Error: Cannot open file '%s' !", input.Data());
        return 1;
    }

    THashList *list_eve = (THashList *)f->FindObjectAny("chists_eve_");
    TH1 *cent = (TH1 *)list_eve->FindObject("hcent");

    // if (xi)
    // {
    THashList *list_Xim = (THashList *)f->FindObjectAny("chists_Xim_");
    THashList *list_Xip = (THashList *)f->FindObjectAny("chists_Xip_");
    // }

    // if (om)
    // {
    THashList *list_Omm = (THashList *)f->FindObjectAny("chists_Omm_");
    THashList *list_Omp = (THashList *)f->FindObjectAny("chists_Omp_");
    // }

    TH3 *Xim = (TH3 *)list_Xim->FindObject("h3_ptmasscent_def");
    TH3 *Xip = (TH3 *)list_Xip->FindObject("h3_ptmasscent_def");
    TH1 *resultParams_Xip, *resultParams_Xim;

    TH3 *Omm = (TH3 *)list_Omm->FindObject("h3_ptmasscent_def");
    TH3 *Omp = (TH3 *)list_Omp->FindObject("h3_ptmasscent_def");
    TH1 *resultParams_Omp, *resultParams_Omm;

    TH2 *MCGenXip = (TH2 *)list_Xip->FindObject("h2_gen");
    TH2 *MCGenXim = (TH2 *)list_Xim->FindObject("h2_gen");
    TH2 *MCGenOmp = (TH2 *)list_Omp->FindObject("h2_gen");
    TH2 *MCGenOmm = (TH2 *)list_Omm->FindObject("h2_gen");

    TH2D *MCGenXiC = (TH2D *)MCGenXip->Clone();
    MCGenXiC->Add(MCGenXim);

    TH2D *MCGenOmC = (TH2D *)MCGenOmp->Clone();
    MCGenOmC->Add(MCGenOmm);

    if (!((Xim && Xip) || (Omm && Omp)))
    {
        Printf("Histograms cannot be found. Aborting.");
        return 2;
    }

    TString output = outputFilename;
    if (outputFilename.IsNull())
    {
        output = input;
        output.ReplaceAll("AnalysisResults.root", "RsnProjection.root");
    }
    Printf("Saving output to '%s' ...", output.Data());
    TFile *out = TFile::Open(output.Data(), "RECREATE");
    if (!out)
    {
        Printf("Error: Cannot open file '%s' !", output.Data());
        return 3;
    }

    if (xi)
    {
        out->mkdir("h_MassXim_pt_mult");
        out->mkdir("h_MassXip_pt_mult");
        out->mkdir("h_MassXiC_pt_mult"); /// for combination of + and -

        out->mkdir("h_MassXim_pt");
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
    Double_t nSigma = 4.;
    // Double_t sigmaExp = 0.0025;
    Double_t sigma = 0.0025;
    Double_t sigmaXi = 0.001875;
    Double_t sigmaOm = 0.0018;
    Double_t fourSigXi = 0.0075;
    Double_t fourSigOm = 0.0072;

    Double_t multbins_Xi[11] = {0, 5, 10, 15, 20, 30, 40, 50, 60, 80, 100}; // V0A
    Double_t multbins_Om[6] = {0, 5, 15, 30, 60, 100};                      // V0A
    Int_t nmultbins_Xi = sizeof(multbins_Xi) / sizeof(Double_t) - 1;
    Int_t nmultbins_Om = sizeof(multbins_Om) / sizeof(Double_t) - 1;

    Double_t ptbins_Xi[] = {0.8, 1.1, 1.3, 1.5, 1.7, 1.9, 2.1, 2.3, 2.5, 2.7, 2.9, 3.1, 3.5, 4, 5.5};
    Double_t ptbins_Om[] = {0.9, 1.6, 2., 2.4, 2.9, 3.5, 5};
    Int_t nptbins_Xi = sizeof(ptbins_Xi) / sizeof(Double_t) - 1;
    Int_t nptbins_Om = sizeof(ptbins_Om) / sizeof(Double_t) - 1;

    Int_t binsMC[4] = {0};

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

    TH1D *h_MassXim_pt_mult[nptbins_Xi][nmultbins_Xi];
    TH1D *h_MassXip_pt_mult[nptbins_Xi][nmultbins_Xi];
    TH1D *h_MassOmm_pt_mult[nptbins_Om][nmultbins_Om];
    TH1D *h_MassOmp_pt_mult[nptbins_Om][nmultbins_Om];
    TH1 *resultParams_xip[nptbins_Xi][nmultbins_Xi];
    TH1 *resultParams_xim[nptbins_Xi][nmultbins_Xi];
    TH1 *resultParams_omp[nptbins_Om][nmultbins_Om];
    TH1 *resultParams_omm[nptbins_Om][nmultbins_Om];

    /// +/- combination:
    TH1D *h_MassXiC_pt_mult[nptbins_Xi][nmultbins_Xi];
    TH1D *h_MassOmC_pt_mult[nptbins_Om][nmultbins_Om];
    TH1 *resultParams_xiC[nptbins_Xi][nmultbins_Xi];
    TH1 *resultParams_omC[nptbins_Om][nmultbins_Om];

    // mult integrated:
    //  TH1D *effXim_pt = new TH1D("effXim_pt", "Efficiency #Xi^{-}", nptbins_Xi, ptbins_Xi);
    //  TH1D *effXip_pt = new TH1D("effXip_pt", "Efficiency #Xi^{+}", nptbins_Xi, ptbins_Xi);
    //  TH1D *effXiC_pt = new TH1D("effXiC_pt", "Efficiency #Xi^{+} + #Xi^{-}", nptbins_Xi, ptbins_Xi);
    //  TH1D *effOmm_pt = new TH1D("effOmm_pt", "Efficiency #Omega^{-}", nptbins_Om, ptbins_Om);
    //  TH1D *effOmp_pt = new TH1D("effOmp_pt", "Efficiency #Omega^{+}", nptbins_Om, ptbins_Om);
    //  TH1D *effOmC_pt = new TH1D("effOmC_pt", "Efficiency #Omega^{+} + #Omega^{-}", nptbins_Om, ptbins_Om);

    TH1D *h_MassXim_pt[nptbins_Xi];
    TH1D *h_MassXip_pt[nptbins_Xi];
    TH1D *h_MassOmm_pt[nptbins_Om];
    TH1D *h_MassOmp_pt[nptbins_Om];
    TH1 *resultParams_pt_xip[nptbins_Xi];
    TH1 *resultParams_pt_xim[nptbins_Xi];
    TH1 *resultParams_pt_omp[nptbins_Om];
    TH1 *resultParams_pt_omm[nptbins_Om];

    /// +/- combination:
    TH1D *h_MassXiC_pt[nptbins_Xi];
    TH1D *h_MassOmC_pt[nptbins_Om];
    TH1 *resultParams_pt_xiC[nptbins_Xi];
    TH1 *resultParams_pt_omC[nptbins_Om];

    // Signal+Bg (peak)
    // if (xi)
    // {
    TH1 *InvMass_Xim = (TH1 *)Xim->Project3D(idAxis);
    InvMass_Xim->SetNameTitle("InvMass_Xim", "Invariant Mass #Xi^{-}");
    TH1 *InvMass_Xip = (TH1 *)Xip->Project3D(idAxis);
    InvMass_Xip->SetNameTitle("InvMass_Xip", "Invariant Mass #Xi^{+}");
    // }

    // if (om)
    // {
    TH1 *InvMass_Omm = (TH1 *)Omm->Project3D(idAxis);
    InvMass_Omm->SetNameTitle("InvMass_Omm", "Invariant Mass #Omega^{-}");
    TH1 *InvMass_Omp = (TH1 *)Omp->Project3D(idAxis);
    InvMass_Omp->SetNameTitle("InvMass_Omp", "Invariant Mass #Omega^{+}");
    // }

    // Setting fitting range for signal
    // Double_t fitMinSig_Xi = lMass_Xi - nSigma * sigmaXi;
    // Double_t fitMaxSig_Xi = lMass_Xi + nSigma * sigmaXi;
    // Double_t fitMinSig_Om = lMass_Om - nSigma * sigmaOm;
    // Double_t fitMaxSig_Om = lMass_Om + nSigma * sigmaOm;
    Double_t fitMinSig_Xi = lMass_Xi - fourSigXi;
    Double_t fitMaxSig_Xi = lMass_Xi + fourSigXi;
    Double_t fitMinSig_Om = lMass_Om - fourSigOm;
    Double_t fitMaxSig_Om = lMass_Om + fourSigOm;

    // Setting fitting range for bg
    // Double_t fitMinBg_Xi = lMass_Xi - 2 * fourSigXi;
    // Double_t fitMaxBg_Xi = lMass_Xi + 2 * fourSigXi;
    // Double_t fitMinBg_Om = lMass_Om - 2 * fourSigOm;
    // Double_t fitMaxBg_Om = lMass_Om + 2 * fourSigOm;
    Double_t fitMinBg_Xi = lMass_Xi - 3 * fourSigXi;
    Double_t fitMaxBg_Xi = lMass_Xi + 3 * fourSigXi;
    Double_t fitMinBg_Om = lMass_Om - 3 * fourSigOm;
    Double_t fitMaxBg_Om = lMass_Om + 3 * fourSigOm;
    // Double_t fitMinBg_Xi = lMass_Xi - 2.5 * nSigma * sigmaXi;
    // Double_t fitMaxBg_Xi = lMass_Xi + 2.5 * nSigma * sigmaXi;
    // Double_t fitMinBg_Om = lMass_Om - 2.5 * nSigma * sigmaOm;
    // Double_t fitMaxBg_Om = lMass_Om + 2.5 * nSigma * sigmaOm;

    // Double_t effXip, effXim, effOmp, effOmm = 0;
    Double_t genXip, genXim, genXiC, genOmp, genOmm, genOmC = 0;

    /// Fitting gaussian on signal for getting sigma for bg
    // TF1 *gausXi = new TF1("gausXi", "gaus", lMass_Xi - 0.04, lMass_Xi + 0.04);
    // gausXi->SetLineColor(kBlack);
    // InvMass_Xim->Fit("gausXi", "REM+");
    // Printf("mean = %f, mean error = %f \nsigma = %f, sigma error = %f", gausXi->GetParameter(1), gausXi->GetParError(1), gausXi->GetParameter(2), gausXi->GetParError(2));
    // TF1 *gausOm = new TF1("gausOm", "gaus", lMass_Om - 0.04, lMass_Om + 0.04);
    // gausOm->SetLineColor(kBlack);
    // InvMass_Omp->Fit("gausOm", "REM+");
    // Printf("mean = %f, mean error = %f \nsigma = %f, sigma error = %f", gausOm->GetParameter(1), gausOm->GetParError(1), gausOm->GetParameter(2), gausOm->GetParError(2));

    // for (Int_t iFit = 0; iFit < 10; iFit++)
    // {
    if (xi)
    {
        for (Int_t ptBinXi = 0; ptBinXi < nptbins_Xi; ptBinXi++)
        {
            for (Int_t multBinXi = 0; multBinXi < nmultbins_Xi; multBinXi++)
            {

                TF1 *peakFnc_Xip = new TF1("GausPol2_Xip", GausPol2, fitMinSig_Xi, fitMaxSig_Xi, 6);
                TF1 *peakFnc_Xim = new TF1("GausPol2_Xim", GausPol2, fitMinSig_Xi, fitMaxSig_Xi, 6);
                TF1 *peakFnc_XiC = new TF1("GausPol2_XiC", GausPol2, fitMinSig_Xi, fitMaxSig_Xi, 6);
                // reject = true;
                TF1 *bgFnc_Xip = new TF1("Pol2Xip", Pol2Xi, fitMinBg_Xi, fitMaxBg_Xi, 3);
                TF1 *bgFnc_Xim = new TF1("Pol2Xim", Pol2Xi, fitMinBg_Xi, fitMaxBg_Xi, 3);
                TF1 *bgFnc_XiC = new TF1("Pol2XiC", Pol2Xi, fitMinBg_Xi, fitMaxBg_Xi, 3);

                bgFnc_Xip->SetLineColor(kGreen);
                bgFnc_Xip->SetFillColor(kYellow);
                bgFnc_Xip->SetFillStyle(3009);

                bgFnc_Xim->SetLineColor(kGreen);
                bgFnc_Xim->SetFillColor(kYellow);
                bgFnc_Xim->SetFillStyle(3009);

                bgFnc_XiC->SetLineColor(kGreen);
                bgFnc_XiC->SetFillColor(kYellow);
                bgFnc_XiC->SetFillStyle(3009);

                // reject = false;
                h_MassXim_pt_mult[ptBinXi][multBinXi] = (TH1D *)Xim->ProjectionY(TString::Format(("h_MassXim_pt_mult[%d][%d]"), ptBinXi, multBinXi), ptBinXi + 1, ptBinXi + 1, multBinXi + 1, multBinXi + 1);
                h_MassXip_pt_mult[ptBinXi][multBinXi] = (TH1D *)Xip->ProjectionY(TString::Format(("h_MassXip_pt_mult[%d][%d]"), ptBinXi, multBinXi), ptBinXi + 1, ptBinXi + 1, multBinXi + 1, multBinXi + 1);
                // h_MassXim_pt_mult[ptBinXi][multBinXi]->Sumw2();
                // h_MassXip_pt_mult[ptBinXi][multBinXi]->Sumw2();
                h_MassXiC_pt_mult[ptBinXi][multBinXi] = (TH1D *)h_MassXip_pt_mult[ptBinXi][multBinXi]->Clone(TString::Format(("h_MassXiC_pt_mult[%d][%d]"), ptBinXi, multBinXi));
                h_MassXiC_pt_mult[ptBinXi][multBinXi]->Add(h_MassXim_pt_mult[ptBinXi][multBinXi]);

                binsMC[0] = ptBinXi + 1;
                binsMC[1] = ptBinXi + 1;
                binsMC[2] = multBinXi + 1;
                binsMC[3] = multBinXi + 1;

                // if (fisMC)
                // {
                //     genXip = GetGeneratedParticles(MCGenXip, ptBinXi + 1, ptBinXi + 2, multBinXi + 1, multBinXi + 2, genErr);
                //     genXim = GetGeneratedParticles(MCGenXim, ptBinXi + 1, ptBinXi + 2, multBinXi + 1, multBinXi + 2, genErr);
                //     genXiC = GetGeneratedParticles(MCGenXiC, ptBinXi + 1, ptBinXi + 2, multBinXi + 1, multBinXi + 2, genErr);
                // }
                peakFnc_Xip = setFitParameters(peakFnc_Xip, h_MassXip_pt_mult[ptBinXi][multBinXi], lMass_Xi, ptBinXi + 1);
                resultParams_xip[ptBinXi][multBinXi] = fitResults(peakFnc_Xip, bgFnc_Xip, h_MassXip_pt_mult[ptBinXi][multBinXi], fitMinSig_Xi, fitMaxSig_Xi, fitMinBg_Xi, fitMaxBg_Xi, MCGenXip, binsMC, kFALSE);
                resultParams_xip[ptBinXi][multBinXi]->SetName(TString::Format(("resultParams_xip[%d][%d]"), ptBinXi, multBinXi));

                peakFnc_Xim = setFitParameters(peakFnc_Xim, h_MassXim_pt_mult[ptBinXi][multBinXi], lMass_Xi, ptBinXi + 1);
                resultParams_xim[ptBinXi][multBinXi] = fitResults(peakFnc_Xim, bgFnc_Xim, h_MassXim_pt_mult[ptBinXi][multBinXi], fitMinSig_Xi, fitMaxSig_Xi, fitMinBg_Xi, fitMaxBg_Xi, MCGenXim, binsMC, kFALSE);
                resultParams_xim[ptBinXi][multBinXi]->SetName(TString::Format(("resultParams_xim[%d][%d]"), ptBinXi, multBinXi));

                peakFnc_XiC = setFitParameters(peakFnc_XiC, h_MassXiC_pt_mult[ptBinXi][multBinXi], lMass_Xi, ptBinXi + 1);
                resultParams_xiC[ptBinXi][multBinXi] = fitResults(peakFnc_XiC, bgFnc_XiC, h_MassXiC_pt_mult[ptBinXi][multBinXi], fitMinSig_Xi, fitMaxSig_Xi, fitMinBg_Xi, fitMaxBg_Xi, MCGenXiC, binsMC, kFALSE);
                resultParams_xiC[ptBinXi][multBinXi]->SetName(TString::Format(("resultParams_xiC[%d][%d]"), ptBinXi, multBinXi));

                out->cd("h_MassXim_pt_mult");
                h_MassXim_pt_mult[ptBinXi][multBinXi]->SetTitle(TString::Format(("Invariant Mass #Xi^{-}: p_{T}<%.1f,%.1f>, Mult<%.1f,%.1f>"), ptbins_Xi[ptBinXi], ptbins_Xi[ptBinXi + 1], multbins_Xi[multBinXi], multbins_Xi[multBinXi + 1]));
                h_MassXim_pt_mult[ptBinXi][multBinXi]->Write();
                resultParams_xim[ptBinXi][multBinXi]->Write();
                out->cd("h_MassXip_pt_mult");
                h_MassXip_pt_mult[ptBinXi][multBinXi]->SetTitle(TString::Format(("Invariant Mass #Xi^{+}: p_{T}<%.1f,%.1f>, Mult<%.1f,%.1f>"), ptbins_Xi[ptBinXi], ptbins_Xi[ptBinXi + 1], multbins_Xi[multBinXi], multbins_Xi[multBinXi + 1]));
                h_MassXip_pt_mult[ptBinXi][multBinXi]->Write();
                resultParams_xip[ptBinXi][multBinXi]->Write();

                out->cd("h_MassXiC_pt_mult");
                h_MassXiC_pt_mult[ptBinXi][multBinXi]->SetTitle(TString::Format(("Invariant Mass #Xi: p_{T}<%.1f,%.1f>, Mult<%.1f,%.1f>"), ptbins_Xi[ptBinXi], ptbins_Xi[ptBinXi + 1], multbins_Xi[multBinXi], multbins_Xi[multBinXi + 1]));
                h_MassXiC_pt_mult[ptBinXi][multBinXi]->Write();
                resultParams_xiC[ptBinXi][multBinXi]->Write();

                // Printf("Xip: IntBC = %f, Xim: IntBC = %f", resultParams_xip[ptBinXi][multBinXi]->GetBinContent(1));
            }

            // plotting mult integrated efficiency as a function of pT
            binsMC[0] = ptBinXi + 1;
            binsMC[1] = ptBinXi + 1;
            binsMC[2] = 1;
            binsMC[3] = nmultbins_Xi;

            TF1 *peakFnc_Xip = new TF1("GausPol2_Xip", GausPol2, fitMinSig_Xi, fitMaxSig_Xi, 6);
            TF1 *peakFnc_Xim = new TF1("GausPol2_Xim", GausPol2, fitMinSig_Xi, fitMaxSig_Xi, 6);
            TF1 *peakFnc_XiC = new TF1("GausPol2_XiC", GausPol2, fitMinSig_Xi, fitMaxSig_Xi, 6);
            // reject = true;
            TF1 *bgFnc_Xip = new TF1("Pol2Xip", Pol2Xi, fitMinBg_Xi, fitMaxBg_Xi, 3);
            TF1 *bgFnc_Xim = new TF1("Pol2Xim", Pol2Xi, fitMinBg_Xi, fitMaxBg_Xi, 3);
            TF1 *bgFnc_XiC = new TF1("Pol2XiC", Pol2Xi, fitMinBg_Xi, fitMaxBg_Xi, 3);

            bgFnc_Xip->SetLineColor(kGreen);
            bgFnc_Xip->SetFillColor(kYellow);
            bgFnc_Xip->SetFillStyle(3009);

            bgFnc_Xim->SetLineColor(kGreen);
            bgFnc_Xim->SetFillColor(kYellow);
            bgFnc_Xim->SetFillStyle(3009);

            bgFnc_XiC->SetLineColor(kGreen);
            bgFnc_XiC->SetFillColor(kYellow);
            bgFnc_XiC->SetFillStyle(3009);

            // reject = false;
            h_MassXim_pt[ptBinXi] = (TH1D *)Xim->ProjectionY(TString::Format(("h_MassXim_pt[%d]"), ptBinXi), ptBinXi + 1, ptBinXi + 1, 1, nmultbins_Xi);
            h_MassXip_pt[ptBinXi] = (TH1D *)Xip->ProjectionY(TString::Format(("h_MassXip_pt_mult[%d]"), ptBinXi), ptBinXi + 1, ptBinXi + 1, 1, nmultbins_Xi);
            h_MassXiC_pt[ptBinXi] = (TH1D *)h_MassXip_pt[ptBinXi]->Clone(TString::Format(("h_MassXiC_pt_mult[%d]"), ptBinXi));
            h_MassXiC_pt[ptBinXi]->Add(h_MassXim_pt[ptBinXi]);

            peakFnc_Xip = setFitParameters(peakFnc_Xip, h_MassXip_pt[ptBinXi], lMass_Xi, ptBinXi + 1);
            resultParams_pt_xip[ptBinXi] = fitResults(peakFnc_Xip, bgFnc_Xip, h_MassXip_pt[ptBinXi], fitMinSig_Xi, fitMaxSig_Xi, fitMinBg_Xi, fitMaxBg_Xi, MCGenXip, binsMC, kFALSE);
            resultParams_pt_xip[ptBinXi]->SetName(TString::Format(("resultParams_pt_xip[%d]"), ptBinXi));

            peakFnc_Xim = setFitParameters(peakFnc_Xim, h_MassXim_pt[ptBinXi], lMass_Xi, ptBinXi + 1);
            resultParams_pt_xim[ptBinXi] = fitResults(peakFnc_Xim, bgFnc_Xim, h_MassXim_pt[ptBinXi], fitMinSig_Xi, fitMaxSig_Xi, fitMinBg_Xi, fitMaxBg_Xi, MCGenXim, binsMC, kFALSE);
            resultParams_pt_xim[ptBinXi]->SetName(TString::Format(("resultParams_pt_xim[%d]"), ptBinXi));

            peakFnc_XiC = setFitParameters(peakFnc_XiC, h_MassXiC_pt[ptBinXi], lMass_Xi, ptBinXi + 1);
            resultParams_pt_xiC[ptBinXi] = fitResults(peakFnc_XiC, bgFnc_XiC, h_MassXiC_pt[ptBinXi], fitMinSig_Xi, fitMaxSig_Xi, fitMinBg_Xi, fitMaxBg_Xi, MCGenXiC, binsMC, kFALSE);
            resultParams_pt_xiC[ptBinXi]->SetName(TString::Format(("resultParams_pt_xiC[%d]"), ptBinXi));

            out->cd("h_MassXim_pt");
            h_MassXim_pt[ptBinXi]->SetTitle(TString::Format(("Invariant Mass #Xi^{-}: p_{T}<%.1f,%.1f>, Mult<%.1f,%.1f>"), ptbins_Xi[ptBinXi], ptbins_Xi[ptBinXi + 1], multbins_Xi[0], multbins_Xi[nmultbins_Xi]));
            h_MassXim_pt[ptBinXi]->Write();
            resultParams_pt_xim[ptBinXi]->Write();
            out->cd("h_MassXip_pt");
            h_MassXip_pt[ptBinXi]->SetTitle(TString::Format(("Invariant Mass #Xi^{+}: p_{T}<%.1f,%.1f>, Mult<%.1f,%.1f>"), ptbins_Xi[ptBinXi], ptbins_Xi[ptBinXi + 1], multbins_Xi[0], multbins_Xi[nmultbins_Xi]));
            h_MassXip_pt[ptBinXi]->Write();
            resultParams_pt_xip[ptBinXi]->Write();

            out->cd("h_MassXiC_pt");
            h_MassXiC_pt[ptBinXi]->SetTitle(TString::Format(("Invariant Mass #Xi: p_{T}<%.1f,%.1f>, Mult<%.1f,%.1f>"), ptbins_Xi[ptBinXi], ptbins_Xi[ptBinXi + 1], multbins_Xi[0], multbins_Xi[nmultbins_Xi]));
            h_MassXiC_pt[ptBinXi]->Write();
            resultParams_pt_xiC[ptBinXi]->Write();
        }

        TF1 *peakFnc_Xip = new TF1("GausPol2_Xip", GausPol2, fitMinSig_Xi, fitMaxSig_Xi, 6);
        TF1 *peakFnc_Xim = new TF1("GausPol2_Xim", GausPol2, fitMinSig_Xi, fitMaxSig_Xi, 6);
        TF1 *peakFnc_XiC = new TF1("GausPol2_XiC", GausPol2, fitMinSig_Xi, fitMaxSig_Xi, 6);
        // reject = true;
        TF1 *bgFnc_Xip = new TF1("Pol2Xip", Pol2Xi, fitMinBg_Xi, fitMaxBg_Xi, 3);
        TF1 *bgFnc_Xim = new TF1("Pol2Xim", Pol2Xi, fitMinBg_Xi, fitMaxBg_Xi, 3);
        TF1 *bgFnc_XiC = new TF1("Pol2XiC", Pol2Xi, fitMinBg_Xi, fitMaxBg_Xi, 3);

        bgFnc_Xip->SetLineColor(kGreen);
        bgFnc_Xip->SetFillColor(kYellow);
        bgFnc_Xip->SetFillStyle(3009);

        bgFnc_Xim->SetLineColor(kGreen);
        bgFnc_Xim->SetFillColor(kYellow);
        bgFnc_Xim->SetFillStyle(3009);

        bgFnc_XiC->SetLineColor(kGreen);
        bgFnc_XiC->SetFillColor(kYellow);
        bgFnc_XiC->SetFillStyle(3009);

        // if (fisMC)
        // {
        //     genXip = GetGeneratedParticles(MCGenXip, 1, nptbins_Xi, 1, nmultbins_Xi, genErr);
        //     genXim = GetGeneratedParticles(MCGenXim, 1, nptbins_Xi, 1, nmultbins_Xi, genErr);
        //     genXiC = GetGeneratedParticles(MCGenXiC, 1, nptbins_Xi, 1, nmultbins_Xi, genErr);
        // }

        binsMC[0] = 1;
        binsMC[1] = nptbins_Xi;
        binsMC[2] = 1;
        binsMC[3] = nmultbins_Xi;

        peakFnc_Xip = setFitParameters(peakFnc_Xip, InvMass_Xip, lMass_Xi, -1);
        resultParams_Xip = fitResults(peakFnc_Xip, bgFnc_Xip, InvMass_Xip, fitMinSig_Xi, fitMaxSig_Xi, fitMinBg_Xi, fitMaxBg_Xi, MCGenXip, binsMC, kFALSE);
        peakFnc_Xim = setFitParameters(peakFnc_Xim, InvMass_Xim, lMass_Xi, -1);
        resultParams_Xim = fitResults(peakFnc_Xim, bgFnc_Xim, InvMass_Xim, fitMinSig_Xi, fitMaxSig_Xi, fitMinBg_Xi, fitMaxBg_Xi, MCGenXim, binsMC, kFALSE);

        out->cd();
        h_multBinEntries_Xi->Write();
        InvMass_Xip->Write();
        InvMass_Xim->Write();
        if (resultParams_Xip)
            resultParams_Xip->Write();
        if (resultParams_Xim)
            resultParams_Xim->Write();
    }
    if (om)
    {
        for (Int_t ptBinOm = 0; ptBinOm < nptbins_Om; ptBinOm++)
        {
            for (Int_t multBinOm = 0; multBinOm < nmultbins_Om; multBinOm++)
            {
                TF1 *peakFnc_Omp = new TF1("GausPol2_Omp", GausPol2, fitMinSig_Om, fitMaxSig_Om, 6);
                TF1 *peakFnc_Omm = new TF1("GausPol2_Omm", GausPol2, fitMinSig_Om, fitMaxSig_Om, 6);
                TF1 *peakFnc_OmC = new TF1("GausPol2_OmC", GausPol2, fitMinSig_Om, fitMaxSig_Om, 6);
                TF1 *bgFnc_Omp = new TF1("Pol2Omp", Pol2Om, fitMinBg_Om, fitMaxBg_Om, 3);
                TF1 *bgFnc_Omm = new TF1("Pol2Omm", Pol2Om, fitMinBg_Om, fitMaxBg_Om, 3);
                TF1 *bgFnc_OmC = new TF1("Pol2OmC", Pol2Om, fitMinBg_Om, fitMaxBg_Om, 3);
                bgFnc_Omp->SetLineColor(kGreen);
                bgFnc_Omp->SetFillColor(kYellow);
                bgFnc_Omp->SetFillStyle(3009);
                bgFnc_Omm->SetLineColor(kGreen);
                bgFnc_Omm->SetFillColor(kYellow);
                bgFnc_Omm->SetFillStyle(3009);
                bgFnc_OmC->SetLineColor(kGreen);
                bgFnc_OmC->SetFillColor(kYellow);
                bgFnc_OmC->SetFillStyle(3009);

                binsMC[0] = ptBinOm + 1;
                binsMC[1] = ptBinOm + 1;
                binsMC[2] = multBinOm + 1;
                binsMC[3] = multBinOm + 1;

                // if (fisMC)
                // {
                //     genOmp = GetGeneratedParticles(MCGenOmp, ptBinOm + 1, ptBinOm + 2, multBinOm + 1, multBinOm + 2, genErr);
                //     genOmm = GetGeneratedParticles(MCGenOmm, ptBinOm + 1, ptBinOm + 2, multBinOm + 1, multBinOm + 2, genErr);
                //     genOmC = GetGeneratedParticles(MCGenOmC, ptBinOm + 1, ptBinOm + 2, multBinOm + 1, multBinOm + 2, genErr);
                // }

                h_MassOmm_pt_mult[ptBinOm][multBinOm] = (TH1D *)Omm->ProjectionY(TString::Format(("h_MassOmm_pt_mult[%d][%d]"), ptBinOm, multBinOm), ptBinOm + 1, ptBinOm + 1, multBinOm + 1, multBinOm + 1);
                h_MassOmp_pt_mult[ptBinOm][multBinOm] = (TH1D *)Omp->ProjectionY(TString::Format(("h_MassOmp_pt_mult[%d][%d]"), ptBinOm, multBinOm), ptBinOm + 1, ptBinOm + 1, multBinOm + 1, multBinOm + 1);
                // h_MassOmm_pt_mult[ptBinOm][multBinOm]->Sumw2();
                // h_MassOmp_pt_mult[ptBinOm][multBinOm]->Sumw2();

                h_MassOmC_pt_mult[ptBinOm][multBinOm] = (TH1D *)h_MassOmp_pt_mult[ptBinOm][multBinOm]->Clone(TString::Format(("h_MassOmC_pt_mult[%d][%d]"), ptBinOm, multBinOm));
                h_MassOmC_pt_mult[ptBinOm][multBinOm]->Add(h_MassOmm_pt_mult[ptBinOm][multBinOm]);

                peakFnc_Omp = setFitParameters(peakFnc_Omp, h_MassOmp_pt_mult[ptBinOm][multBinOm], lMass_Om, ptBinOm + 1);
                resultParams_omp[ptBinOm][multBinOm] = fitResults(peakFnc_Omp, bgFnc_Omp, h_MassOmp_pt_mult[ptBinOm][multBinOm], fitMinSig_Om, fitMaxSig_Om, fitMinBg_Om, fitMaxBg_Om, MCGenOmp, binsMC, kFALSE);
                resultParams_omp[ptBinOm][multBinOm]->SetName(TString::Format(("resultParams_omp[%d][%d]"), ptBinOm, multBinOm));
                peakFnc_Omm = setFitParameters(peakFnc_Omm, h_MassOmm_pt_mult[ptBinOm][multBinOm], lMass_Om, ptBinOm + 1);
                resultParams_omm[ptBinOm][multBinOm] = fitResults(peakFnc_Omm, bgFnc_Omm, h_MassOmm_pt_mult[ptBinOm][multBinOm], fitMinSig_Om, fitMaxSig_Om, fitMinBg_Om, fitMaxBg_Om, MCGenOmm, binsMC, kFALSE);
                resultParams_omm[ptBinOm][multBinOm]->SetName(TString::Format(("resultParams_omm[%d][%d]"), ptBinOm, multBinOm));

                peakFnc_OmC = setFitParameters(peakFnc_OmC, h_MassOmC_pt_mult[ptBinOm][multBinOm], lMass_Om, ptBinOm + 1);
                resultParams_omC[ptBinOm][multBinOm] = fitResults(peakFnc_OmC, bgFnc_OmC, h_MassOmC_pt_mult[ptBinOm][multBinOm], fitMinSig_Om, fitMaxSig_Om, fitMinBg_Om, fitMaxBg_Om, MCGenOmC, binsMC, kFALSE);
                resultParams_omC[ptBinOm][multBinOm]->SetName(TString::Format(("resultParams_omC[%d][%d]"), ptBinOm, multBinOm));

                out->cd("h_MassOmm_pt_mult");
                h_MassOmm_pt_mult[ptBinOm][multBinOm]->SetTitle(TString::Format(("Invariant Mass #Omega^{-}: p_{T}<%.1f,%.1f>, Mult<%.1f,%.1f>"), ptbins_Om[ptBinOm], ptbins_Om[ptBinOm + 1], multbins_Om[multBinOm], multbins_Om[multBinOm + 1]));
                h_MassOmm_pt_mult[ptBinOm][multBinOm]->Write();
                resultParams_omm[ptBinOm][multBinOm]->Write();
                out->cd("h_MassOmp_pt_mult");
                h_MassOmp_pt_mult[ptBinOm][multBinOm]->SetTitle(TString::Format(("Invariant Mass #Omega^{+}: p_{T}<%.1f,%.1f>, Mult<%.1f,%.1f>"), ptbins_Om[ptBinOm], ptbins_Om[ptBinOm + 1], multbins_Om[multBinOm], multbins_Om[multBinOm + 1]));
                h_MassOmp_pt_mult[ptBinOm][multBinOm]->Write();
                resultParams_omp[ptBinOm][multBinOm]->Write();

                out->cd("h_MassOmC_pt_mult");
                h_MassOmC_pt_mult[ptBinOm][multBinOm]->SetTitle(TString::Format(("Invariant Mass #Omega: p_{T}<%.1f,%.1f>, Mult<%.1f,%.1f>"), ptbins_Om[ptBinOm], ptbins_Om[ptBinOm + 1], multbins_Om[multBinOm], multbins_Om[multBinOm + 1]));
                h_MassOmC_pt_mult[ptBinOm][multBinOm]->Write();
                resultParams_omC[ptBinOm][multBinOm]->Write();
            }

            // plotting mult integrated efficiency as a function of pT
            binsMC[0] = ptBinOm + 1;
            binsMC[1] = ptBinOm + 1;
            binsMC[2] = 1;
            binsMC[3] = nmultbins_Om;

            TF1 *peakFnc_Omp = new TF1("GausPol2_Omp", GausPol2, fitMinSig_Om, fitMaxSig_Om, 6);
            TF1 *peakFnc_Omm = new TF1("GausPol2_Omm", GausPol2, fitMinSig_Om, fitMaxSig_Om, 6);
            TF1 *peakFnc_OmC = new TF1("GausPol2_OmC", GausPol2, fitMinSig_Om, fitMaxSig_Om, 6);
            // reject = true;
            TF1 *bgFnc_Omp = new TF1("Pol2Omp", Pol2Om, fitMinBg_Om, fitMaxBg_Om, 3);
            TF1 *bgFnc_Omm = new TF1("Pol2Omm", Pol2Om, fitMinBg_Om, fitMaxBg_Om, 3);
            TF1 *bgFnc_OmC = new TF1("Pol2OmC", Pol2Om, fitMinBg_Om, fitMaxBg_Om, 3);

            bgFnc_Omp->SetLineColor(kGreen);
            bgFnc_Omp->SetFillColor(kYellow);
            bgFnc_Omp->SetFillStyle(3009);

            bgFnc_Omm->SetLineColor(kGreen);
            bgFnc_Omm->SetFillColor(kYellow);
            bgFnc_Omm->SetFillStyle(3009);

            bgFnc_OmC->SetLineColor(kGreen);
            bgFnc_OmC->SetFillColor(kYellow);
            bgFnc_OmC->SetFillStyle(3009);

            // reject = false;
            h_MassOmm_pt[ptBinOm] = (TH1D *)Omm->ProjectionY(TString::Format(("h_MassOmm_pt[%d]"), ptBinOm), ptBinOm + 1, ptBinOm + 1, 1, nmultbins_Om);
            h_MassOmp_pt[ptBinOm] = (TH1D *)Omp->ProjectionY(TString::Format(("h_MassOmp_pt_mult[%d]"), ptBinOm), ptBinOm + 1, ptBinOm + 1, 1, nmultbins_Om);
            h_MassOmC_pt[ptBinOm] = (TH1D *)h_MassOmp_pt[ptBinOm]->Clone(TString::Format(("h_MassOmC_pt_mult[%d]"), ptBinOm));
            h_MassOmC_pt[ptBinOm]->Add(h_MassOmm_pt[ptBinOm]);

            peakFnc_Omp = setFitParameters(peakFnc_Omp, h_MassOmp_pt[ptBinOm], lMass_Om, ptBinOm + 1);
            resultParams_pt_omp[ptBinOm] = fitResults(peakFnc_Omp, bgFnc_Omp, h_MassOmp_pt[ptBinOm], fitMinSig_Om, fitMaxSig_Om, fitMinBg_Om, fitMaxBg_Om, MCGenOmp, binsMC, kFALSE);
            resultParams_pt_omp[ptBinOm]->SetName(TString::Format(("resultParams_pt_omp[%d]"), ptBinOm));

            peakFnc_Omm = setFitParameters(peakFnc_Omm, h_MassOmm_pt[ptBinOm], lMass_Om, ptBinOm + 1);
            resultParams_pt_omm[ptBinOm] = fitResults(peakFnc_Omm, bgFnc_Omm, h_MassOmm_pt[ptBinOm], fitMinSig_Om, fitMaxSig_Om, fitMinBg_Om, fitMaxBg_Om, MCGenOmm, binsMC, kFALSE);
            resultParams_pt_omm[ptBinOm]->SetName(TString::Format(("resultParams_pt_omm[%d]"), ptBinOm));

            peakFnc_OmC = setFitParameters(peakFnc_OmC, h_MassOmC_pt[ptBinOm], lMass_Om, ptBinOm + 1);
            resultParams_pt_omC[ptBinOm] = fitResults(peakFnc_OmC, bgFnc_OmC, h_MassOmC_pt[ptBinOm], fitMinSig_Om, fitMaxSig_Om, fitMinBg_Om, fitMaxBg_Om, MCGenOmC, binsMC, kFALSE);
            resultParams_pt_omC[ptBinOm]->SetName(TString::Format(("resultParams_pt_omC[%d]"), ptBinOm));

            out->cd("h_MassOmm_pt");
            h_MassOmm_pt[ptBinOm]->SetTitle(TString::Format(("Invariant Mass #Omega^{-}: p_{T}<%.1f,%.1f>, Mult<%.1f,%.1f>"), ptbins_Om[ptBinOm], ptbins_Om[ptBinOm + 1], multbins_Om[0], multbins_Om[nmultbins_Om]));
            h_MassOmm_pt[ptBinOm]->Write();
            resultParams_pt_omm[ptBinOm]->Write();
            out->cd("h_MassOmp_pt");
            h_MassOmp_pt[ptBinOm]->SetTitle(TString::Format(("Invariant Mass #Omega^{+}: p_{T}<%.1f,%.1f>, Mult<%.1f,%.1f>"), ptbins_Om[ptBinOm], ptbins_Om[ptBinOm + 1], multbins_Om[0], multbins_Om[nmultbins_Om]));
            h_MassOmp_pt[ptBinOm]->Write();
            resultParams_pt_omp[ptBinOm]->Write();

            out->cd("h_MassOmC_pt");
            h_MassOmC_pt[ptBinOm]->SetTitle(TString::Format(("Invariant Mass #Omega: p_{T}<%.1f,%.1f>, Mult<%.1f,%.1f>"), ptbins_Om[ptBinOm], ptbins_Om[ptBinOm + 1], multbins_Om[0], multbins_Om[nmultbins_Om]));
            h_MassOmC_pt[ptBinOm]->Write();
            resultParams_pt_omC[ptBinOm]->Write();
        }

        // InvMass_Xim->Fit(peakFnc_Xi, "QN MFCLES R", "", fitMinSig_Xi, fitMaxSig_Xi);
        // InvMass_Xip->Fit(peakFnc_Xi, "QN MFCLES R", "", fitMinSig_Xi, fitMaxSig_Xi);
        // InvMass_Omm->Fit(peakFnc_Om, "QN MFCLES R", "", fitMinSig_Om, fitMaxSig_Om);
        // InvMass_Omp->Fit(peakFnc_Om, "QN MFCLES R", "", fitMinSig_Om, fitMaxSig_Om);

        TF1 *peakFnc_Omp = new TF1("GausPol2_Omp", GausPol2, fitMinSig_Om, fitMaxSig_Om, 6);
        TF1 *peakFnc_Omm = new TF1("GausPol2_Omm", GausPol2, fitMinSig_Om, fitMaxSig_Om, 6);
        TF1 *peakFnc_OmC = new TF1("GausPol2_OmC", GausPol2, fitMinSig_Om, fitMaxSig_Om, 6);
        TF1 *bgFnc_Omp = new TF1("Pol2Omp", Pol2Om, fitMinBg_Om, fitMaxBg_Om, 3);
        TF1 *bgFnc_Omm = new TF1("Pol2Omm", Pol2Om, fitMinBg_Om, fitMaxBg_Om, 3);
        TF1 *bgFnc_OmC = new TF1("Pol2OmC", Pol2Om, fitMinBg_Om, fitMaxBg_Om, 3);
        bgFnc_Omp->SetLineColor(kGreen);
        bgFnc_Omp->SetFillColor(kYellow);
        bgFnc_Omp->SetFillStyle(3009);
        bgFnc_Omm->SetLineColor(kGreen);
        bgFnc_Omm->SetFillColor(kYellow);
        bgFnc_Omm->SetFillStyle(3009);
        bgFnc_OmC->SetLineColor(kGreen);
        bgFnc_OmC->SetFillColor(kYellow);
        bgFnc_OmC->SetFillStyle(3009);

        // if (fisMC)
        // {
        //     genOmp = GetGeneratedParticles(MCGenOmp, 1, nptbins_Om, 1, nmultbins_Om, genErr);
        //     genOmm = GetGeneratedParticles(MCGenOmm, 1, nptbins_Om, 1, nmultbins_Om, genErr);
        //     genOmC = GetGeneratedParticles(MCGenOmC, 1, nptbins_Om, 1, nmultbins_Om, genErr);
        // }

        binsMC[0] = 1;
        binsMC[1] = nptbins_Om;
        binsMC[2] = 1;
        binsMC[3] = nmultbins_Om;

        peakFnc_Omp = setFitParameters(peakFnc_Omp, InvMass_Omp, lMass_Om, -1);
        resultParams_Omp = fitResults(peakFnc_Omp, bgFnc_Omp, InvMass_Omp, fitMinSig_Om, fitMaxSig_Om, fitMinBg_Om, fitMaxBg_Om, MCGenOmp, binsMC, kFALSE);
        peakFnc_Omm = setFitParameters(peakFnc_Omm, InvMass_Omm, lMass_Om, -1);
        resultParams_Omm = fitResults(peakFnc_Omm, bgFnc_Omm, InvMass_Omm, fitMinSig_Om, fitMaxSig_Om, fitMinBg_Om, fitMaxBg_Om, MCGenOmm, binsMC, kFALSE);

        out->cd();
        InvMass_Omm->Write();
        InvMass_Omp->Write();
        h_multBinEntries_Om->Write();

        if (resultParams_Omp)
            resultParams_Omp->Write();
        if (resultParams_Omm)
            resultParams_Omm->Write();
    }

    // if (xi)
    // {
    // }
    // if (om)
    // {
    //     // }
    // }

    /// temporarily commented out:
    /*if (!gROOT->IsBatch())
    {
        TCanvas *c = new TCanvas("XiOmegaFitting", "Cascade Analysis: Xi and Omega");
        c->Divide(2, 2);
        c->cd(1);
        InvMass_Xip->DrawCopy();
        c->cd(2);
        InvMass_Xim->DrawCopy();
        c->cd(3);
        InvMass_Omp->DrawCopy();
        c->cd(4);
        InvMass_Omm->DrawCopy();
    }*/

    /// end of temporarily commented out

    /// testing
    // int nXi = count_Xi;
    // int nOm = count_Om;
    // double avg_mean_Xi = 0;
    // double avg_fourSigma_Xi = 0;
    // double avg_mean_Om = 0;
    // double avg_fourSigma_Om = 0;

    // // TH1D *testing_Xi = new TH1D("testing_Xi", "fourSigma_Xi", 400, 0.0065, 0.0085);
    // // Printf("count_Xi = %d , count_Om = %d", nXi, nOm);

    // for (int iXi = 0; iXi < nXi; iXi++)
    // {
    //     if (mean_Xi[iXi] != 0 && !(TMath::IsNaN(mean_Xi[iXi])))
    //     {
    //         avg_mean_Xi = avg_mean_Xi + mean_Xi[iXi];
    //     }

    //     if (fourSigma_Xi[iXi] != 0 && !(TMath::IsNaN(fourSigma_Xi[iXi])))
    //     {
    //         avg_fourSigma_Xi = avg_fourSigma_Xi + fourSigma_Xi[iXi];
    //         // testing_Xi->Fill(fourSigma_Xi[iXi]);
    //     }
    //     Printf("mean_Xi[%d] = %f, fourSigma_Xi[%d] = %f", iXi, mean_Xi[iXi], iXi, fourSigma_Xi[iXi]);
    // }

    // avg_mean_Xi = avg_mean_Xi / nXi;
    // avg_fourSigma_Xi = avg_fourSigma_Xi / nXi;

    // //Omega testing
    // for (int iOm = 0; iOm < nOm; iOm++)
    // {
    //     if (mean_Om[iOm] != 0 && !(TMath::IsNaN(mean_Om[iOm])))
    //     {
    //         avg_mean_Om = avg_mean_Om + mean_Om[iOm];
    //     }

    //     if (fourSigma_Om[iOm] != 0 && !(TMath::IsNaN(fourSigma_Om[iOm])))
    //     {
    //         avg_fourSigma_Om = avg_fourSigma_Om + fourSigma_Om[iOm];
    //     }
    //     Printf("mean_Om[%d] = %f, fourSigma_Om[%d] = %f", iOm, mean_Om[iOm], iOm, fourSigma_Om[iOm]);
    // }

    // avg_mean_Om = avg_mean_Om / nOm;
    // avg_fourSigma_Om = avg_fourSigma_Om / nOm;

    // Printf("avg_mean_Xi (Xi Mass) = %f \navg_fourSigma_Xi (signal interval Xi) = %f \navg_mean_Om (Omega Mass) = %f \navg_fourSigma_Om (signal interval Omega) = %f", avg_mean_Xi, avg_fourSigma_Xi, avg_mean_Om, avg_fourSigma_Om);

    // TF1 *testFit = new TF1("testFit", "gaus", 0.0065, 0.0085);
    // testing_Xi->Fit(testFit);
    // TCanvas *c1 = new TCanvas("Canvas_testing_Xi", "fourSigma_Xi testing");
    // c1->cd();
    // testing_Xi->Draw();

    out->Close();
    f->Close();

    // delete f, out;
    return 0;
}

TF1 *setFitParameters(TF1 *peakFnc, TH1 *peak, Double_t mass, Int_t ptBin)
{
    peakFnc->SetParNames("Amplitude", "Mean", "#sigma");
    peakFnc->SetParameter(1, mass);
    peakFnc->SetParameter(2, 0.0018);
    peakFnc->SetParLimits(2, 0.001, 0.01);

    if (ptBin == -1)
    { /// case when all bins are integrated
        peakFnc->SetParameter(0, peak->GetMaximum() * 0.9);
        peakFnc->SetParameter(3, 0);
        peakFnc->SetParameter(4, 0);
        peakFnc->SetParameter(5, peak->GetMaximum() * 0.1);
        return peakFnc;
    }
    else
    {
        peakFnc->SetParameter(1, mass);
        peakFnc->SetParLimits(1, mass - 0.003, mass + 0.003);
        peakFnc->SetParameter(2, 0.0018);

        if (ptBin <= 4)
            peakFnc->SetParLimits(2, 0.0005, 0.0032);
        else
            peakFnc->SetParLimits(2, 0.0005, 0.0038); // 35
        peakFnc->SetParameter(0, peak->GetMaximum() * 0.9);
        peakFnc->SetParameter(3, 0);
        peakFnc->SetParameter(4, 0);
        peakFnc->SetParameter(5, peak->GetMaximum() * 0.1);
        return peakFnc;
    }
}

TH1 *fitResults(TF1 *peakFnc, TF1 *bgFnc, TH1 *peak, Double_t fitMinSig, Double_t fitMaxSig, Double_t fitMinBg, Double_t fitMaxBg, TH2 *MCgen, Int_t binsMC[], bool options = kFALSE)
{
    // reject = true;

    if (TMath::IsNaN(peak->GetEntries()))
    {
        return 0;
    }

    /// use different fit options:
    TString opt = "QMINLSR";

    TString optBg = "QMNLSR";
    if (options)
        opt = "QNLS";

    TFitResultPtr fFitResult;

    for (Int_t iFit = 0; iFit < 1; iFit++)
    {
        fFitResult = peak->Fit(peakFnc, opt.Data(), "", fitMinSig, fitMaxSig);

        // Double_t par[6];
        // peakFnc->GetParameters(par);
        // const Double_t *parErr = peakFnc->GetParErrors();

        /// Defining peak limits for signal region (green lines):
        ///  par[1] = peak position, par[2] = peak width
        // Double_t lPeakLeftLimit = par[1] - 1. * 4 * TMath::Abs(par[2]);
        // Double_t lPeakRightLimit = par[1] + 1. * 4 * TMath::Abs(par[2]);
        // TLine *lLineLeft = new TLine(lPeakLeftLimit, 0, lPeakLeftLimit, peak->GetMaximum() * 0.95);
        // TLine *lLineRight = new TLine(lPeakRightLimit, 0, lPeakRightLimit, peak->GetMaximum() * 0.95);
        // lLineLeft->SetLineColor(kGreen + 2);
        // lLineRight->SetLineColor(kGreen + 2);
        /// TODO: need to draw the line with each histogram...can't save it on top of Hist?

        // Getting pol2 part from fit and set it to bgFnc for drawing and then
        // subtracting
        reject = true;

        // bgFnc->SetParameters(&par[3]);

        // bgFnc->SetParameter(0, par[3]);
        // bgFnc->SetParameter(1, par[4]);
        // bgFnc->SetParameter(2, par[5]);
        // bgFnc->SetParError(0, parErr[3]);
        // bgFnc->SetParError(1, parErr[4]);
        // bgFnc->SetParError(2, parErr[5]);
        peak->Fit(bgFnc, optBg.Data(), "", fitMinBg, fitMaxBg);

        // TF1 *bgLeft = new TF1("bgLeft", Pol2, fitMinBg, (par[1] - bgReject), 3);
        // peak->Fit(bgLeft, opt.Data(), "", fitMinBg, (par[1] - bgReject));
        // bgLeft->SetParameters(bgFnc->GetParameters());
        // peak->GetListOfFunctions()->Add(bgLeft);

        // TF1 *bgRight = new TF1("bgRight", Pol2, (par[1] + bgReject), fitMaxBg, 3);
        // bgRight->SetParameters(bgFnc->GetParameters());
        // peak->Fit(bgRight, opt.Data(), "", (par[1] + bgReject), fitMaxBg);
        // peak->GetListOfFunctions()->Add(bgRight);

        reject = false;
    }

    // add fit and background function to histogram so it is automatically drawn
    // with hist

    peak->GetListOfFunctions()->Add(peakFnc);
    peak->GetListOfFunctions()->Add(bgFnc);

    // store 2 separate functions for visualization
    //  TF1 *fleft = new TF1("fleft", fline, 0, 2.5, 2);
    //  fleft->SetParameters(fl->GetParameters());
    //  h->GetListOfFunctions()->Add(fleft);
    //  gROOT->GetListOfFunctions()->Remove(fleft);
    //  TF1 *fright = new TF1("fright", fline, 3.5, 5, 2);
    //  fright->SetParameters(fl->GetParameters());
    //  h->GetListOfFunctions()->Add(fright);
    //  gROOT->GetListOfFunctions()->Remove(fright);

    /// getting min/max for deciding signal's fit range for integral calculation
    Double_t minInt = peakFnc->GetParameter(1) - 4 * peakFnc->GetParameter(2); /// mean - 4*sigma
    Double_t maxInt = peakFnc->GetParameter(1) + 4 * peakFnc->GetParameter(2); /// mean + 4*sigma
    // Printf("mean = %f, 4*sigma = %f, %s", peakFnc->GetParameter(1), 4 * peakFnc->GetParameter(2), peakFnc->GetName());
    // if (strcmp(peakFnc->GetName(), "GausPol2_Xi") == 0)
    // {
    //     mean_Xi[count_Xi] = peakFnc->GetParameter(1);
    //     fourSigma_Xi[count_Xi] = 4 * peakFnc->GetParameter(2);
    //     count_Xi++;
    // }
    // if (strcmp(peakFnc->GetName(), "GausPol2_Om") == 0)
    // {
    //     mean_Om[count_Om] = peakFnc->GetParameter(1);
    //     fourSigma_Om[count_Om] = 4 * peakFnc->GetParameter(2);
    //     count_Om++;
    // }
    // Int_t bin_minSig = peak->GetXaxis()->FindBin(minInt);
    // Int_t bin_maxSig = peak->GetXaxis()->FindBin(maxInt);

    // Int_t bin_minBg = peak->GetXaxis()->FindBin(fitMinBg);
    // Int_t bin_maxBg = peak->GetXaxis()->FindBin(fitMaxBg);

    TH1 *resultParams = fillParams(peakFnc, peak, fFitResult, minInt, maxInt, MCgen, binsMC);
    peak->SetOption("X0E1");
    return resultParams;
}
TH1 *fillParams(TF1 *sigBgFnc, TH1 *peak, TFitResultPtr fFitResult, Double_t minInt, Double_t maxInt, TH2 *MCgen, Int_t binsMC[])
{
    Double_t val, err, gen, genErr;
    Double_t val_intBC, err_intBC, eff, effErr;
    const Int_t nBins = 8 + sigBgFnc->GetNpar();
    Double_t par[6];
    Double_t parErr[6];
    TH1 *resultParams = new TH1D("resultParams", "Result parameters", nBins, 0, nBins);
    Int_t iBin = 1;

    GetYieldBinCounting(peak, fFitResult, minInt, maxInt, val, err);
    if (TMath::IsNaN(val) || TMath::IsNaN(err))
    {
        Printf("GetYieldBinCounting is nan");
        delete resultParams;
        return nullptr;
    }
    val_intBC = val;
    err_intBC = err;

    resultParams->SetBinContent(iBin, val_intBC);
    resultParams->SetBinError(iBin, err_intBC);
    resultParams->GetXaxis()->SetBinLabel(iBin, "IntBC");
    iBin++;

    GetYieldFitFunction(peak, fFitResult, minInt, maxInt, val, err);
    if (TMath::IsNaN(val) || TMath::IsNaN(err))
    {
        Printf("GetYieldFitFunction is nan.");
        delete resultParams;
        return nullptr;
    }
    resultParams->SetBinContent(iBin, val);
    resultParams->SetBinError(iBin, err);
    resultParams->GetXaxis()->SetBinLabel(iBin, "IntFF");
    iBin++;

    // if (fisMC)
    // {
    gen = GetGeneratedParticles(MCgen, binsMC, genErr);

    resultParams->SetBinContent(iBin, gen);
    resultParams->SetBinError(iBin, genErr);
    resultParams->GetXaxis()->SetBinLabel(iBin, "IntGen");
    iBin++;

    // Printf("[%f, %f]: genErr = %f", minInt, maxInt, genErr);
    if (gen != 0 && genErr != 0)
    {
        eff = val_intBC / gen;
        // effErr = err_intBC / genErr; /// ask if it should be computed like this?

        effErr = ErrorInRatio(val_intBC, err_intBC, gen, genErr);
        resultParams->SetBinContent(iBin, eff);
        // resultParams->SetBinError(iBin, err_intBC);
        resultParams->SetBinError(iBin, effErr);
        resultParams->GetXaxis()->SetBinLabel(iBin, "Efficiency");
        iBin++;
    }
    // }
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
    return resultParams;
}

void GetYieldBinCounting(TH1 *h, TFitResultPtr fFitResult, Double_t minInt, Double_t maxInt, Double_t &val,
                         Double_t &err, Double_t eps = 1e-6)
{

    if (!h)
        return;

    Int_t bin_minSig = h->GetXaxis()->FindBin(minInt);
    Int_t bin_maxSig = h->GetXaxis()->FindBin(maxInt);

    Double_t histWidth = h->GetXaxis()->GetBinWidth(1);

    val = h->IntegralAndError(bin_minSig, bin_maxSig, err);

    TF1 *bgFnc = (TF1 *)h->GetListOfFunctions()->At(1);

    Double_t bg = bgFnc->Integral(minInt, maxInt, eps);

    // TODO Verify it
    Double_t bgErr = bgFnc->IntegralError(minInt, maxInt, fFitResult->GetParams(),
                                          fFitResult->GetCovarianceMatrix().GetMatrixArray(), eps);

    bg /= histWidth;
    bgErr /= histWidth;
    ////
    Printf("%s: [%f, %f], bins:[%d, %d]: bg = %f, val(signal) = %f", h->GetName(), minInt, maxInt, bin_minSig, bin_maxSig, bg, val);
    if (val > 30)
    {
        val -= bg;
        err = TMath::Sqrt(TMath::Power(err, 2) + TMath::Power(bgErr, 2));
    }
}

void GetYieldFitFunction(TH1 *h, TFitResultPtr fFitResult, Double_t minInt, Double_t maxInt, Double_t &val,
                         Double_t &err, Double_t eps = 1e-6)
{

    if (!h)
        return;

    Double_t histWidth = h->GetXaxis()->GetBinWidth(1);

    TF1 *sigBgFnc = (TF1 *)h->GetListOfFunctions()->At(0);
    TF1 *bgFnc = (TF1 *)h->GetListOfFunctions()->At(1);

    val = sigBgFnc->Integral(minInt, maxInt, eps);
    err = sigBgFnc->IntegralError(minInt, maxInt, fFitResult->GetParams(), fFitResult->GetCovarianceMatrix().GetMatrixArray(),
                                  eps);

    val /= histWidth;
    err /= histWidth;

    Double_t bg = bgFnc->Integral(minInt, maxInt, eps);

    // TODO Verify it
    Double_t bgErr = bgFnc->IntegralError(minInt, maxInt, fFitResult->GetParams(),
                                          fFitResult->GetCovarianceMatrix().GetMatrixArray(), eps);

    bg /= histWidth;
    bgErr /= histWidth;

    if (val > 30)
    {
        val -= bg;
        err = TMath::Sqrt(TMath::Power(err, 2) + TMath::Power(bgErr, 2));
    }
}
