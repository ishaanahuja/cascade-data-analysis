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
TF1 *setFitParameters(TF1 *peakFnc, TH1 *peak, Double_t mass, Int_t ptBin);

TH1 *fitResults(TF1 *peakFnc, TF1 *bgFnc, TH1 *peak, Double_t fitMinSig, Double_t fitMaxSig, Double_t fitMinBg, Double_t fitMaxBg, bool options = kFALSE);

TH1 *fillParams(TF1 *sigBgFnc, TH1 *peak, TFitResultPtr fFitResult, Double_t minInt, Double_t maxInt);

void GetYieldBinCounting(TH1 *h, TFitResultPtr fFitResult, Double_t minInt, Double_t maxInt, Double_t &val,
                         Double_t &err, Double_t eps = 1e-6);
void GetYieldFitFunction(TH1 *h, TFitResultPtr fFitResult, Double_t minInt, Double_t maxInt, Double_t &val,
                         Double_t &err, Double_t eps = 1e-6);

int AliCascadeAnalysis_Ishaan_Fitting(TString input = "AnalysisResults.root", TString idAxis = "y",
                                      TString outputFilename = "AliCascadeAnalysis_Ishaan_Fitting.root")
{
    TH1::AddDirectory(0);

    TFile *f = TFile::Open(input);
    if (!f)
    {
        Printf("Error: Cannot open file '%s' !", input.Data());
        return 1;
    }

    THashList *list_eve = (THashList *)f->FindObjectAny("chists_eve_");
    THashList *list_Xim = (THashList *)f->FindObjectAny("chists_Xim_");
    THashList *list_Xip = (THashList *)f->FindObjectAny("chists_Xip_");
    THashList *list_Omm = (THashList *)f->FindObjectAny("chists_Omm_");
    THashList *list_Omp = (THashList *)f->FindObjectAny("chists_Omp_");

    TH1 *cent = (TH1 *)list_eve->FindObject("hcent");
    TH3 *Xim = (TH3 *)list_Xim->FindObject("h3_ptmasscent_def");
    TH3 *Xip = (TH3 *)list_Xip->FindObject("h3_ptmasscent_def");
    TH3 *Omm = (TH3 *)list_Omm->FindObject("h3_ptmasscent_def");
    TH3 *Omp = (TH3 *)list_Omp->FindObject("h3_ptmasscent_def");
    if (!(Xim && Xip && Omm && Omp))
    {
        Printf("Histograms cannot be found. Aborting.");
        return 2;
    }

    TH1 *resultParams_Xip, *resultParams_Xim, *resultParams_Omp, *resultParams_Omm;

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
    out->mkdir("h_MassXim_pt_mult");
    out->mkdir("h_MassXip_pt_mult");
    out->mkdir("h_MassOmm_pt_mult");
    out->mkdir("h_MassOmp_pt_mult");

    /// for combination of + and -
    out->mkdir("h_MassXiC_pt_mult");
    out->mkdir("h_MassOmC_pt_mult");

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

    // Signal+Bg (peak)
    TH1 *InvMass_Xim = (TH1 *)Xim->Project3D(idAxis);
    InvMass_Xim->SetNameTitle("InvMass_Xim", "Invariant Mass #Xi^{-}");
    TH1 *InvMass_Xip = (TH1 *)Xip->Project3D(idAxis);
    InvMass_Xip->SetNameTitle("InvMass_Xip", "Invariant Mass #Xi^{+}");
    TH1 *InvMass_Omm = (TH1 *)Omm->Project3D(idAxis);
    InvMass_Omm->SetNameTitle("InvMass_Omm", "Invariant Mass #Omega^{-}");
    TH1 *InvMass_Omp = (TH1 *)Omp->Project3D(idAxis);
    InvMass_Omp->SetNameTitle("InvMass_Omp", "Invariant Mass #Omega^{+}");

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
    Double_t fitMinBg_Xi = lMass_Xi - 2 * fourSigXi;
    Double_t fitMaxBg_Xi = lMass_Xi + 2 * fourSigXi;
    Double_t fitMinBg_Om = lMass_Om - 2 * fourSigOm;
    Double_t fitMaxBg_Om = lMass_Om + 2 * fourSigOm;
    // Double_t fitMinBg_Xi = lMass_Xi - 2.5 * nSigma * sigmaXi;
    // Double_t fitMaxBg_Xi = lMass_Xi + 2.5 * nSigma * sigmaXi;
    // Double_t fitMinBg_Om = lMass_Om - 2.5 * nSigma * sigmaOm;
    // Double_t fitMaxBg_Om = lMass_Om + 2.5 * nSigma * sigmaOm;

    /// Fitting gaussian on signal for getting sigma for bg
    // TF1 *gausXi = new TF1("gausXi", "gaus", lMass_Xi - 0.04, lMass_Xi + 0.04);
    // gausXi->SetLineColor(kBlack);
    // InvMass_Xim->Fit("gausXi", "REM+");
    // Printf("mean = %f, mean error = %f \nsigma = %f, sigma error = %f", gausXi->GetParameter(1), gausXi->GetParError(1), gausXi->GetParameter(2), gausXi->GetParError(2));
    // TF1 *gausOm = new TF1("gausOm", "gaus", lMass_Om - 0.04, lMass_Om + 0.04);
    // gausOm->SetLineColor(kBlack);
    // InvMass_Omp->Fit("gausOm", "REM+");
    // Printf("mean = %f, mean error = %f \nsigma = %f, sigma error = %f", gausOm->GetParameter(1), gausOm->GetParError(1), gausOm->GetParameter(2), gausOm->GetParError(2));

    for (Int_t iFit = 0; iFit < 1; iFit++)
    {
        for (Int_t ptBinXi = 0; ptBinXi < nptbins_Xi; ptBinXi++)
        {
            for (Int_t multBinXi = 0; multBinXi < nmultbins_Xi; multBinXi++)
            {

                TF1 *peakFnc_Xip = new TF1("GausPol2_Xip", GausPol2, fitMinSig_Xi, fitMaxSig_Xi, 6);
                TF1 *peakFnc_Xim = new TF1("GausPol2_Xim", GausPol2, fitMinSig_Xi, fitMaxSig_Xi, 6);
                TF1 *peakFnc_XiC = new TF1("GausPol2_XiC", GausPol2, fitMinSig_Xi, fitMaxSig_Xi, 6);
                // reject = true;
                TF1 *bgFnc_Xip = new TF1("Pol1Xip", Pol1Xi, fitMinBg_Xi, fitMaxBg_Xi, 2);
                TF1 *bgFnc_Xim = new TF1("Pol1Xim", Pol1Xi, fitMinBg_Xi, fitMaxBg_Xi, 2);
                TF1 *bgFnc_XiC = new TF1("Pol1XiC", Pol1Xi, fitMinBg_Xi, fitMaxBg_Xi, 2);

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

                peakFnc_Xip = setFitParameters(peakFnc_Xip, h_MassXip_pt_mult[ptBinXi][multBinXi], lMass_Xi, ptBinXi + 1);
                resultParams_xip[ptBinXi][multBinXi] = fitResults(peakFnc_Xip, bgFnc_Xip, h_MassXip_pt_mult[ptBinXi][multBinXi], fitMinSig_Xi, fitMaxSig_Xi, fitMinBg_Xi, fitMaxBg_Xi, kFALSE);
                resultParams_xip[ptBinXi][multBinXi]->SetName(TString::Format(("resultParams_xip[%d][%d]"), ptBinXi, multBinXi));

                peakFnc_Xim = setFitParameters(peakFnc_Xim, h_MassXim_pt_mult[ptBinXi][multBinXi], lMass_Xi, ptBinXi + 1);
                resultParams_xim[ptBinXi][multBinXi] = fitResults(peakFnc_Xim, bgFnc_Xim, h_MassXim_pt_mult[ptBinXi][multBinXi], fitMinSig_Xi, fitMaxSig_Xi, fitMinBg_Xi, fitMaxBg_Xi, kFALSE);
                resultParams_xim[ptBinXi][multBinXi]->SetName(TString::Format(("resultParams_xim[%d][%d]"), ptBinXi, multBinXi));

                peakFnc_XiC = setFitParameters(peakFnc_XiC, h_MassXiC_pt_mult[ptBinXi][multBinXi], lMass_Xi, ptBinXi + 1);
                resultParams_xiC[ptBinXi][multBinXi] = fitResults(peakFnc_XiC, bgFnc_XiC, h_MassXiC_pt_mult[ptBinXi][multBinXi], fitMinSig_Xi, fitMaxSig_Xi, fitMinBg_Xi, fitMaxBg_Xi, kFALSE);
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
        }
        for (Int_t ptBinOm = 0; ptBinOm < nptbins_Om; ptBinOm++)
        {
            for (Int_t multBinOm = 0; multBinOm < nmultbins_Om; multBinOm++)
            {
                TF1 *peakFnc_Omp = new TF1("GausPol2_Omp", GausPol2, fitMinSig_Om, fitMaxSig_Om, 6);
                TF1 *peakFnc_Omm = new TF1("GausPol2_Omm", GausPol2, fitMinSig_Om, fitMaxSig_Om, 6);
                TF1 *peakFnc_OmC = new TF1("GausPol2_OmC", GausPol2, fitMinSig_Om, fitMaxSig_Om, 6);
                TF1 *bgFnc_Omp = new TF1("Pol1Omp", Pol1Om, fitMinBg_Om, fitMaxBg_Om, 2);
                TF1 *bgFnc_Omm = new TF1("Pol1Omm", Pol1Om, fitMinBg_Om, fitMaxBg_Om, 2);
                TF1 *bgFnc_OmC = new TF1("Pol1OmC", Pol1Om, fitMinBg_Om, fitMaxBg_Om, 2);
                bgFnc_Omp->SetLineColor(kGreen);
                bgFnc_Omp->SetFillColor(kYellow);
                bgFnc_Omp->SetFillStyle(3009);
                bgFnc_Omm->SetLineColor(kGreen);
                bgFnc_Omm->SetFillColor(kYellow);
                bgFnc_Omm->SetFillStyle(3009);
                bgFnc_OmC->SetLineColor(kGreen);
                bgFnc_OmC->SetFillColor(kYellow);
                bgFnc_OmC->SetFillStyle(3009);

                h_MassOmm_pt_mult[ptBinOm][multBinOm] = (TH1D *)Omm->ProjectionY(TString::Format(("h_MassOmm_pt_mult[%d][%d]"), ptBinOm, multBinOm), ptBinOm + 1, ptBinOm + 1, multBinOm + 1, multBinOm + 1);
                h_MassOmp_pt_mult[ptBinOm][multBinOm] = (TH1D *)Omp->ProjectionY(TString::Format(("h_MassOmp_pt_mult[%d][%d]"), ptBinOm, multBinOm), ptBinOm + 1, ptBinOm + 1, multBinOm + 1, multBinOm + 1);
                // h_MassOmm_pt_mult[ptBinOm][multBinOm]->Sumw2();
                // h_MassOmp_pt_mult[ptBinOm][multBinOm]->Sumw2();

                h_MassOmC_pt_mult[ptBinOm][multBinOm] = (TH1D *)h_MassOmp_pt_mult[ptBinOm][multBinOm]->Clone(TString::Format(("h_MassOmC_pt_mult[%d][%d]"), ptBinOm, multBinOm));
                h_MassOmC_pt_mult[ptBinOm][multBinOm]->Add(h_MassOmm_pt_mult[ptBinOm][multBinOm]);

                peakFnc_Omp = setFitParameters(peakFnc_Omp, h_MassOmp_pt_mult[ptBinOm][multBinOm], lMass_Om, ptBinOm + 1);
                resultParams_omp[ptBinOm][multBinOm] = fitResults(peakFnc_Omp, bgFnc_Omp, h_MassOmp_pt_mult[ptBinOm][multBinOm], fitMinSig_Om, fitMaxSig_Om, fitMinBg_Om, fitMaxBg_Om, kFALSE);
                resultParams_omp[ptBinOm][multBinOm]->SetName(TString::Format(("resultParams_omp[%d][%d]"), ptBinOm, multBinOm));
                peakFnc_Omm = setFitParameters(peakFnc_Omm, h_MassOmm_pt_mult[ptBinOm][multBinOm], lMass_Om, ptBinOm + 1);
                resultParams_omm[ptBinOm][multBinOm] = fitResults(peakFnc_Omm, bgFnc_Omm, h_MassOmm_pt_mult[ptBinOm][multBinOm], fitMinSig_Om, fitMaxSig_Om, fitMinBg_Om, fitMaxBg_Om, kFALSE);
                resultParams_omm[ptBinOm][multBinOm]->SetName(TString::Format(("resultParams_omm[%d][%d]"), ptBinOm, multBinOm));

                peakFnc_OmC = setFitParameters(peakFnc_OmC, h_MassOmC_pt_mult[ptBinOm][multBinOm], lMass_Om, ptBinOm + 1);
                resultParams_omC[ptBinOm][multBinOm] = fitResults(peakFnc_OmC, bgFnc_OmC, h_MassOmC_pt_mult[ptBinOm][multBinOm], fitMinSig_Om, fitMaxSig_Om, fitMinBg_Om, fitMaxBg_Om, kFALSE);
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
        }

        // InvMass_Xim->Fit(peakFnc_Xi, "QN MFCLES R", "", fitMinSig_Xi, fitMaxSig_Xi);
        // InvMass_Xip->Fit(peakFnc_Xi, "QN MFCLES R", "", fitMinSig_Xi, fitMaxSig_Xi);
        // InvMass_Omm->Fit(peakFnc_Om, "QN MFCLES R", "", fitMinSig_Om, fitMaxSig_Om);
        // InvMass_Omp->Fit(peakFnc_Om, "QN MFCLES R", "", fitMinSig_Om, fitMaxSig_Om);
        TF1 *peakFnc_Xip = new TF1("GausPol2_Xip", GausPol2, fitMinSig_Xi, fitMaxSig_Xi, 6);
        TF1 *peakFnc_Xim = new TF1("GausPol2_Xim", GausPol2, fitMinSig_Xi, fitMaxSig_Xi, 6);
        TF1 *peakFnc_XiC = new TF1("GausPol2_XiC", GausPol2, fitMinSig_Xi, fitMaxSig_Xi, 6);
        // reject = true;
        TF1 *bgFnc_Xip = new TF1("Pol1Xip", Pol1Xi, fitMinBg_Xi, fitMaxBg_Xi, 2);
        TF1 *bgFnc_Xim = new TF1("Pol1Xim", Pol1Xi, fitMinBg_Xi, fitMaxBg_Xi, 2);
        TF1 *bgFnc_XiC = new TF1("Pol1XiC", Pol1Xi, fitMinBg_Xi, fitMaxBg_Xi, 2);

        bgFnc_Xip->SetLineColor(kGreen);
        bgFnc_Xip->SetFillColor(kYellow);
        bgFnc_Xip->SetFillStyle(3009);

        bgFnc_Xim->SetLineColor(kGreen);
        bgFnc_Xim->SetFillColor(kYellow);
        bgFnc_Xim->SetFillStyle(3009);

        bgFnc_XiC->SetLineColor(kGreen);
        bgFnc_XiC->SetFillColor(kYellow);
        bgFnc_XiC->SetFillStyle(3009);

        TF1 *peakFnc_Omp = new TF1("GausPol2_Omp", GausPol2, fitMinSig_Om, fitMaxSig_Om, 6);
        TF1 *peakFnc_Omm = new TF1("GausPol2_Omm", GausPol2, fitMinSig_Om, fitMaxSig_Om, 6);
        TF1 *peakFnc_OmC = new TF1("GausPol2_OmC", GausPol2, fitMinSig_Om, fitMaxSig_Om, 6);
        TF1 *bgFnc_Omp = new TF1("Pol1Omp", Pol1Om, fitMinBg_Om, fitMaxBg_Om, 2);
        TF1 *bgFnc_Omm = new TF1("Pol1Omm", Pol1Om, fitMinBg_Om, fitMaxBg_Om, 2);
        TF1 *bgFnc_OmC = new TF1("Pol1OmC", Pol1Om, fitMinBg_Om, fitMaxBg_Om, 2);
        bgFnc_Omp->SetLineColor(kGreen);
        bgFnc_Omp->SetFillColor(kYellow);
        bgFnc_Omp->SetFillStyle(3009);
        bgFnc_Omm->SetLineColor(kGreen);
        bgFnc_Omm->SetFillColor(kYellow);
        bgFnc_Omm->SetFillStyle(3009);
        bgFnc_OmC->SetLineColor(kGreen);
        bgFnc_OmC->SetFillColor(kYellow);
        bgFnc_OmC->SetFillStyle(3009);

        peakFnc_Xip = setFitParameters(peakFnc_Xip, InvMass_Xip, lMass_Xi, -1);
        resultParams_Xip = fitResults(peakFnc_Xip, bgFnc_Xip, InvMass_Xip, fitMinSig_Xi, fitMaxSig_Xi, fitMinBg_Xi, fitMaxBg_Xi, kFALSE);
        peakFnc_Xim = setFitParameters(peakFnc_Xim, InvMass_Xim, lMass_Xi, -1);
        resultParams_Xim = fitResults(peakFnc_Xim, bgFnc_Xim, InvMass_Xim, fitMinSig_Xi, fitMaxSig_Xi, fitMinBg_Xi, fitMaxBg_Xi, kFALSE);
        peakFnc_Omp = setFitParameters(peakFnc_Omp, InvMass_Omp, lMass_Om, -1);
        resultParams_Omp = fitResults(peakFnc_Omp, bgFnc_Omp, InvMass_Omp, fitMinSig_Om, fitMaxSig_Om, fitMinBg_Om, fitMaxBg_Om, kFALSE);
        peakFnc_Omm = setFitParameters(peakFnc_Omm, InvMass_Omm, lMass_Om, -1);
        resultParams_Omm = fitResults(peakFnc_Omm, bgFnc_Omm, InvMass_Omm, fitMinSig_Om, fitMaxSig_Om, fitMinBg_Om, fitMaxBg_Om, kFALSE);
    }

    if (!gROOT->IsBatch())
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
    }

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

    out->cd();
    h_multBinEntries_Xi->Write();
    h_multBinEntries_Om->Write();
    InvMass_Xip->Write();
    InvMass_Xim->Write();
    InvMass_Omm->Write();
    InvMass_Omp->Write();
    if (resultParams_Xip)
        resultParams_Xip->Write();
    if (resultParams_Xim)
        resultParams_Xim->Write();
    if (resultParams_Omp)
        resultParams_Omp->Write();
    if (resultParams_Omm)
        resultParams_Omm->Write();
    out->Close();
    f->Close();

    delete f, out;
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

TH1 *fitResults(TF1 *peakFnc, TF1 *bgFnc, TH1 *peak, Double_t fitMinSig, Double_t fitMaxSig, Double_t fitMinBg, Double_t fitMaxBg, bool options = kFALSE)
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
        opt = "QNLSR";
    TFitResultPtr fFitResult = peak->Fit(peakFnc, opt.Data(), "", fitMinSig, fitMaxSig);

    Double_t par[6];
    peakFnc->GetParameters(par);
    const Double_t *parErr = peakFnc->GetParErrors();

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

    TH1 *resultParams = fillParams(peakFnc, peak, fFitResult, minInt, maxInt);
    peak->SetOption("X0E1");
    return resultParams;
}
TH1 *fillParams(TF1 *sigBgFnc, TH1 *peak, TFitResultPtr fFitResult, Double_t minInt, Double_t maxInt)
{
    Double_t val, err;
    const Int_t nBins = 6 + sigBgFnc->GetNpar();
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
    resultParams->SetBinContent(iBin, val);
    resultParams->SetBinError(iBin, err);
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
    val -= bg;
    err = TMath::Sqrt(TMath::Power(err, 2) + TMath::Power(bgErr, 2));
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

    val -= bg;
    err = TMath::Sqrt(TMath::Power(err, 2) + TMath::Power(bgErr, 2));
}
