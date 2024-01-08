#include <TString.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TF1.h>
#include <TFitResult.h>
#include <TFitResultPtr.h>
#include <TError.h>
#include <THashList.h>

// double bgReject_Xi = 0.0075; // 4*sigma=0.01, 6*sigma=0.015
// double bgReject_Om = 0.0072; // 4*sigma=0.01, 6*sigma=0.015
// double bgReject_Xi = 0.01125; //4*sigma=0.01, 6*sigma=0.015
// double bgReject_Om = 0.0108;  //4*sigma=0.01, 6*sigma=0.015
// bool reject = false; /// rejecting the region +/-4sigma for fitting background
// double bgReject = 0.00;
// double peakFitMass = 0.00;

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

TF1 *setFitParameters(TF1 *peakFnc, TH1 *peak, Double_t mass, Int_t ptBin, Bool_t multInt = kFALSE);

TH1 *fitResults(TF1 *peakFnc, TH1 *peak, Double_t fitMinSig, Double_t fitMaxSig, Double_t *par_allint, Double_t *par_allint_errors, Int_t binsMC[], TH2 *MCgen, bool fisMC = false, Bool_t allInt = kFALSE, Bool_t options = kFALSE);

TH1 *fillParams(TF1 *sigBgFnc, TH1 *peak, TFitResultPtr fFitResult, Double_t minInt, Double_t maxInt, TH2 *MCgen, Int_t binsMC[], bool fisMC);

void GetYieldBinCounting(TH1 *h, TFitResultPtr fFitResult, Double_t minInt, Double_t maxInt, bool fisMC, Double_t &val,
                         Double_t &err, Double_t eps = 1e-6);
void GetYieldFitFunction(TH1 *h, TFitResultPtr fFitResult, Double_t minInt, Double_t maxInt, bool fisMC, Double_t &val,
                         Double_t &err, Double_t eps = 1e-6);

int AliCascadeAnalysisMC_Ishaan_Fitting(TString input = "AnalysisResults.root", TString idAxis = "ye",
                                        TString outputFilename = "AliCascadeAnalysisMC_Ishaan_Fitting.root", bool fisMC = true, bool xi = true, bool om = true)
{
    TH1::AddDirectory(0);
    // TVirtualFitter::SetMaxIterations(1000000);

    // gPrintViaErrorHandler = kTRUE; // set if using option "E" in fitting (disables messages by TMinuit) -> redirects all printing via ROOT's error handler
    // gErrorIgnoreLevel = kWarning;   // same as before

    TFile *f = TFile::Open(input);
    if (!f)
    {
        Printf("Error: Cannot open file '%s' !", input.Data());
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
    Double_t nSigma = 4.;
    Double_t sigma = 0.0025;
    Double_t sigmaXi = 0.001875;
    Double_t sigmaOm = 0.0018;
    Double_t fourSigXi = 0.0075;
    Double_t fourSigOm = 0.0072;

    // test: Emily's
    // Double_t sigmaXi = 0.0025;
    // Double_t sigmaOm = 0.0025;
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

    Double_t par_allintP[6];
    Double_t par_allint_errorsP[6];
    Double_t par_allintM[6];
    Double_t par_allint_errorsM[6];
    Double_t par_allintC[6];
    Double_t par_allint_errorsC[6];

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
        TF1 *peakFnc_XiP = new TF1("GausPol2_XiP", GausPol2, fitMinSig_Xi, fitMaxSig_Xi, 6);
        TF1 *peakFnc_XiM = new TF1("GausPol2_XiM", GausPol2, fitMinSig_Xi, fitMaxSig_Xi, 6);
        TF1 *peakFnc_XiC = new TF1("GausPol2_XiC", GausPol2, fitMinSig_Xi, fitMaxSig_Xi, 6);
        /// pt+mult integrated case
        binsMC[0] = 1;
        binsMC[1] = nptbins_Xi;
        binsMC[2] = 1;
        binsMC[3] = nmultbins_Xi;

        peakFnc_XiP = setFitParameters(peakFnc_XiP, h_MassXip, lMass_Xi, -1, kTRUE);
        resultParams_Xip_allInt = fitResults(peakFnc_XiP, h_MassXip, fitMinSig_Xi, fitMaxSig_Xi, &par_allintP[0], &par_allint_errorsP[0], binsMC, MCGenXip, fisMC, kTRUE);
        peakFnc_XiM = setFitParameters(peakFnc_XiM, h_MassXim, lMass_Xi, -1, kTRUE);
        resultParams_Xim_allInt = fitResults(peakFnc_XiM, h_MassXim, fitMinSig_Xi, fitMaxSig_Xi, &par_allintM[0], &par_allint_errorsM[0], binsMC, MCGenXim, fisMC, kTRUE);
        peakFnc_XiC = setFitParameters(peakFnc_XiC, h_MassXiC, lMass_Xi, -1, kTRUE);
        resultParams_XiC_allInt = fitResults(peakFnc_XiC, h_MassXiC, fitMinSig_Xi, fitMaxSig_Xi, &par_allintC[0], &par_allint_errorsC[0], binsMC, MCGenXiC, fisMC, kTRUE);

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

        for (Int_t ptBinXi = 0; ptBinXi < nptbins_Xi; ptBinXi++)
        {
            if (fisMC)
            {
                // plotting mult integrated efficiency as a function of pT
                binsMC[0] = ptBinXi + 1;
                binsMC[1] = ptBinXi + 1;
                binsMC[2] = 1;
                binsMC[3] = nmultbins_Xi;
            }

            h_MassXim_pt[ptBinXi] = (TH1D *)Xim->ProjectionY(TString::Format(("h_MassXim_pt[%d]"), ptBinXi), ptBinXi + 1, ptBinXi + 1, 1, nmultbins_Xi, "e");

            h_MassXip_pt[ptBinXi] = (TH1D *)Xip->ProjectionY(TString::Format(("h_MassXip_pt[%d]"), ptBinXi), ptBinXi + 1, ptBinXi + 1, 1, nmultbins_Xi, "e");

            h_MassXiC_pt[ptBinXi] = (TH1D *)h_MassXip_pt[ptBinXi]->Clone(TString::Format(("h_MassXiC_pt[%d]"), ptBinXi));
            h_MassXiC_pt[ptBinXi]->Add(h_MassXim_pt[ptBinXi]);

            peakFnc_XiP = setFitParameters(peakFnc_XiP, h_MassXip_pt[ptBinXi], lMass_Xi, ptBinXi + 1, kTRUE);
            resultParXip_pt[ptBinXi] = fitResults(peakFnc_XiP, h_MassXip_pt[ptBinXi], fitMinSig_Xi, fitMaxSig_Xi, &par_allintP[0], &par_allint_errorsP[0], binsMC, MCGenXip, fisMC, kFALSE);
            resultParXip_pt[ptBinXi]->SetName(TString::Format(("resultParXip_pt[%d]"), ptBinXi));

            peakFnc_XiM = setFitParameters(peakFnc_XiM, h_MassXim_pt[ptBinXi], lMass_Xi, ptBinXi + 1, kTRUE);
            resultParXim_pt[ptBinXi] = fitResults(peakFnc_XiM, h_MassXim_pt[ptBinXi], fitMinSig_Xi, fitMaxSig_Xi, &par_allintM[0], &par_allint_errorsM[0], binsMC, MCGenXim, fisMC, kFALSE);
            resultParXim_pt[ptBinXi]->SetName(TString::Format(("resultParXim_pt[%d]"), ptBinXi));

            peakFnc_XiC = setFitParameters(peakFnc_XiC, h_MassXiC_pt[ptBinXi], lMass_Xi, ptBinXi + 1, kTRUE);
            resultParXiC_pt[ptBinXi] = fitResults(peakFnc_XiC, h_MassXiC_pt[ptBinXi], fitMinSig_Xi, fitMaxSig_Xi, &par_allintC[0], &par_allint_errorsC[0], binsMC, MCGenXiC, fisMC, kFALSE);
            resultParXiC_pt[ptBinXi]->SetName(TString::Format(("resultParXiC_pt[%d]"), ptBinXi));

            out->cd("h_MassXim_pt");
            h_MassXim_pt[ptBinXi]->SetTitle(TString::Format(("Invariant Mass #Xi^{-}: p_{T}<%.1f,%.1f>, Mult<%.1f,%.1f>"), ptbins_Xi[ptBinXi], ptbins_Xi[ptBinXi + 1], multbins_Xi[0], multbins_Xi[nmultbins_Xi]));
            resultParXim_pt[ptBinXi]->SetTitle(TString::Format(("Result Parameters #Xi^{-}: p_{T}<%.1f,%.1f>, Mult<%.1f,%.1f>"), ptbins_Xi[ptBinXi], ptbins_Xi[ptBinXi + 1], multbins_Xi[0], multbins_Xi[nmultbins_Xi]));
            h_MassXim_pt[ptBinXi]->Write();
            resultParXim_pt[ptBinXi]->Write();
            out->cd("h_MassXip_pt");
            h_MassXip_pt[ptBinXi]->SetTitle(TString::Format(("Invariant Mass #Xi^{+}: p_{T}<%.1f,%.1f>, Mult<%.1f,%.1f>"), ptbins_Xi[ptBinXi], ptbins_Xi[ptBinXi + 1], multbins_Xi[0], multbins_Xi[nmultbins_Xi]));
            resultParXip_pt[ptBinXi]->SetTitle(TString::Format(("Result Parameters #Xi^{+}: p_{T}<%.1f,%.1f>, Mult<%.1f,%.1f>"), ptbins_Xi[ptBinXi], ptbins_Xi[ptBinXi + 1], multbins_Xi[0], multbins_Xi[nmultbins_Xi]));
            h_MassXip_pt[ptBinXi]->Write();
            resultParXip_pt[ptBinXi]->Write();

            out->cd("h_MassXiC_pt");
            h_MassXiC_pt[ptBinXi]->SetTitle(TString::Format(("Invariant Mass #Xi: p_{T}<%.1f,%.1f>, Mult<%.1f,%.1f>"), ptbins_Xi[ptBinXi], ptbins_Xi[ptBinXi + 1], multbins_Xi[0], multbins_Xi[nmultbins_Xi]));
            resultParXiC_pt[ptBinXi]->SetTitle(TString::Format(("Result Parameters #Xi: p_{T}<%.1f,%.1f>, Mult<%.1f,%.1f>"), ptbins_Xi[ptBinXi], ptbins_Xi[ptBinXi + 1], multbins_Xi[0], multbins_Xi[nmultbins_Xi]));
            h_MassXiC_pt[ptBinXi]->Write();
            resultParXiC_pt[ptBinXi]->Write();

            for (Int_t multBinXi = 0; multBinXi < nmultbins_Xi; multBinXi++)
            {

                if (fisMC) /// fully differential case
                {
                    binsMC[2] = multBinXi + 1;
                    binsMC[3] = multBinXi + 1;
                }

                // reject = false;
                h_MassXim_pt_mult[ptBinXi][multBinXi] = (TH1D *)Xim->ProjectionY(TString::Format(("h_MassXim_pt_mult[%d][%d]"), ptBinXi, multBinXi), ptBinXi + 1, ptBinXi + 1, multBinXi + 1, multBinXi + 1, "e");

                h_MassXip_pt_mult[ptBinXi][multBinXi] = (TH1D *)Xip->ProjectionY(TString::Format(("h_MassXip_pt_mult[%d][%d]"), ptBinXi, multBinXi), ptBinXi + 1, ptBinXi + 1, multBinXi + 1, multBinXi + 1, "e");

                h_MassXiC_pt_mult[ptBinXi][multBinXi] = (TH1D *)h_MassXip_pt_mult[ptBinXi][multBinXi]->Clone(TString::Format(("h_MassXiC_pt_mult[%d][%d]"), ptBinXi, multBinXi));
                h_MassXiC_pt_mult[ptBinXi][multBinXi]->Add(h_MassXim_pt_mult[ptBinXi][multBinXi]);

                peakFnc_XiP = setFitParameters(peakFnc_XiP, h_MassXip_pt_mult[ptBinXi][multBinXi], lMass_Xi, ptBinXi + 1, kFALSE);
                resultParXip_pt_mult[ptBinXi][multBinXi] = fitResults(peakFnc_XiP, h_MassXip_pt_mult[ptBinXi][multBinXi], fitMinSig_Xi, fitMaxSig_Xi, &par_allintP[0], &par_allint_errorsP[0], binsMC, MCGenXip, fisMC, kFALSE);
                resultParXip_pt_mult[ptBinXi][multBinXi]->SetName(TString::Format(("resultParXip_pt_mult[%d][%d]"), ptBinXi, multBinXi));

                peakFnc_XiM = setFitParameters(peakFnc_XiM, h_MassXim_pt_mult[ptBinXi][multBinXi], lMass_Xi, ptBinXi + 1, kFALSE);
                resultParXim_pt_mult[ptBinXi][multBinXi] = fitResults(peakFnc_XiM, h_MassXim_pt_mult[ptBinXi][multBinXi], fitMinSig_Xi, fitMaxSig_Xi, &par_allintM[0], &par_allint_errorsM[0], binsMC, MCGenXim, fisMC, kFALSE);
                resultParXim_pt_mult[ptBinXi][multBinXi]->SetName(TString::Format(("resultParXim_pt_mult[%d][%d]"), ptBinXi, multBinXi));

                peakFnc_XiC = setFitParameters(peakFnc_XiC, h_MassXiC_pt_mult[ptBinXi][multBinXi], lMass_Xi, ptBinXi + 1, kFALSE);
                resultParXiC_pt_mult[ptBinXi][multBinXi] = fitResults(peakFnc_XiC, h_MassXiC_pt_mult[ptBinXi][multBinXi], fitMinSig_Xi, fitMaxSig_Xi, &par_allintC[0], &par_allint_errorsC[0], binsMC, MCGenXiC, fisMC, kFALSE);
                resultParXiC_pt_mult[ptBinXi][multBinXi]->SetName(TString::Format(("resultParXiC_pt_mult[%d][%d]"), ptBinXi, multBinXi));

                out->cd("h_MassXim_pt_mult");
                h_MassXim_pt_mult[ptBinXi][multBinXi]->SetTitle(TString::Format(("Invariant Mass #Xi^{-}: p_{T}<%.1f,%.1f>, Mult<%.1f,%.1f>"), ptbins_Xi[ptBinXi], ptbins_Xi[ptBinXi + 1], multbins_Xi[multBinXi], multbins_Xi[multBinXi + 1]));
                resultParXim_pt_mult[ptBinXi][multBinXi]->SetTitle(TString::Format(("Result Parameters #Xi^{-}: p_{T}<%.1f,%.1f>, Mult<%.1f,%.1f>"), ptbins_Xi[ptBinXi], ptbins_Xi[ptBinXi + 1], multbins_Xi[multBinXi], multbins_Xi[multBinXi + 1]));
                h_MassXim_pt_mult[ptBinXi][multBinXi]->Write();
                resultParXim_pt_mult[ptBinXi][multBinXi]->Write();
                out->cd("h_MassXip_pt_mult");
                h_MassXip_pt_mult[ptBinXi][multBinXi]->SetTitle(TString::Format(("Invariant Mass #Xi^{+}: p_{T}<%.1f,%.1f>, Mult<%.1f,%.1f>"), ptbins_Xi[ptBinXi], ptbins_Xi[ptBinXi + 1], multbins_Xi[multBinXi], multbins_Xi[multBinXi + 1]));
                resultParXip_pt_mult[ptBinXi][multBinXi]->SetTitle(TString::Format(("Result Parameters #Xi^{+}: p_{T}<%.1f,%.1f>, Mult<%.1f,%.1f>"), ptbins_Xi[ptBinXi], ptbins_Xi[ptBinXi + 1], multbins_Xi[multBinXi], multbins_Xi[multBinXi + 1]));
                h_MassXip_pt_mult[ptBinXi][multBinXi]->Write();
                resultParXip_pt_mult[ptBinXi][multBinXi]->Write();

                out->cd("h_MassXiC_pt_mult");
                h_MassXiC_pt_mult[ptBinXi][multBinXi]->SetTitle(TString::Format(("Invariant Mass #Xi: p_{T}<%.1f,%.1f>, Mult<%.1f,%.1f>"), ptbins_Xi[ptBinXi], ptbins_Xi[ptBinXi + 1], multbins_Xi[multBinXi], multbins_Xi[multBinXi + 1]));
                resultParXiC_pt_mult[ptBinXi][multBinXi]->SetTitle(TString::Format(("Result Parameters #Xi: p_{T}<%.1f,%.1f>, Mult<%.1f,%.1f>"), ptbins_Xi[ptBinXi], ptbins_Xi[ptBinXi + 1], multbins_Xi[multBinXi], multbins_Xi[multBinXi + 1]));
                h_MassXiC_pt_mult[ptBinXi][multBinXi]->Write();
                resultParXiC_pt_mult[ptBinXi][multBinXi]->Write();
            }
        }
    }

    if (om)
    {
        TF1 *peakFnc_OmP = new TF1("GausPol2_OmP", GausPol2, fitMinSig_Om, fitMaxSig_Om, 6);
        TF1 *peakFnc_OmM = new TF1("GausPol2_OmM", GausPol2, fitMinSig_Om, fitMaxSig_Om, 6);
        TF1 *peakFnc_OmC = new TF1("GausPol2_OmC", GausPol2, fitMinSig_Om, fitMaxSig_Om, 6);

        /// pt+mult integrated case

        binsMC[0] = 1;
        binsMC[1] = nptbins_Om;
        binsMC[2] = 1;
        binsMC[3] = nmultbins_Om;

        peakFnc_OmP = setFitParameters(peakFnc_OmP, h_MassOmp, lMass_Om, -1, kTRUE);
        resultParams_Omp_allInt = fitResults(peakFnc_OmP, h_MassOmp, fitMinSig_Om, fitMaxSig_Om, &par_allintP[0], &par_allint_errorsP[0], binsMC, MCGenOmp, fisMC, kTRUE);
        peakFnc_OmM = setFitParameters(peakFnc_OmM, h_MassOmm, lMass_Om, -1, kTRUE);
        resultParams_Omm_allInt = fitResults(peakFnc_OmM, h_MassOmm, fitMinSig_Om, fitMaxSig_Om, &par_allintM[0], &par_allint_errorsM[0], binsMC, MCGenOmm, fisMC, kTRUE);
        peakFnc_OmC = setFitParameters(peakFnc_OmC, h_MassOmC, lMass_Om, -1, kTRUE);
        resultParams_OmC_allInt = fitResults(peakFnc_OmC, h_MassOmC, fitMinSig_Om, fitMaxSig_Om, &par_allintC[0], &par_allint_errorsC[0], binsMC, MCGenOmC, fisMC, kTRUE);

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

        for (Int_t ptBinOm = 0; ptBinOm < nptbins_Om; ptBinOm++)
        {
            if (fisMC)
            {
                // plotting mult integrated efficiency as a function of pT
                binsMC[0] = ptBinOm + 1;
                binsMC[1] = ptBinOm + 1;
                binsMC[2] = 1;
                binsMC[3] = nmultbins_Om;
            }

            h_MassOmm_pt[ptBinOm] = (TH1D *)Omm->ProjectionY(TString::Format(("h_MassOmm_pt[%d]"), ptBinOm), ptBinOm + 1, ptBinOm + 1, 1, nmultbins_Om, "e");

            h_MassOmp_pt[ptBinOm] = (TH1D *)Omp->ProjectionY(TString::Format(("h_MassOmp_pt[%d]"), ptBinOm), ptBinOm + 1, ptBinOm + 1, 1, nmultbins_Om, "e");

            h_MassOmC_pt[ptBinOm] = (TH1D *)h_MassOmp_pt[ptBinOm]->Clone(TString::Format(("h_MassOmC_pt[%d]"), ptBinOm));
            h_MassOmC_pt[ptBinOm]->Add(h_MassOmm_pt[ptBinOm]);

            peakFnc_OmP = setFitParameters(peakFnc_OmP, h_MassOmp_pt[ptBinOm], lMass_Om, ptBinOm + 1, kTRUE);
            resultParOmp_pt[ptBinOm] = fitResults(peakFnc_OmP, h_MassOmp_pt[ptBinOm], fitMinSig_Om, fitMaxSig_Om, &par_allintP[0], &par_allint_errorsP[0], binsMC, MCGenOmp, fisMC, kFALSE);
            resultParOmp_pt[ptBinOm]->SetName(TString::Format(("resultParOmp_pt[%d]"), ptBinOm));

            peakFnc_OmM = setFitParameters(peakFnc_OmM, h_MassOmm_pt[ptBinOm], lMass_Om, ptBinOm + 1, kTRUE);
            resultParOmm_pt[ptBinOm] = fitResults(peakFnc_OmM, h_MassOmm_pt[ptBinOm], fitMinSig_Om, fitMaxSig_Om, &par_allintM[0], &par_allint_errorsM[0], binsMC, MCGenOmm, fisMC, kFALSE);
            resultParOmm_pt[ptBinOm]->SetName(TString::Format(("resultParOmm_pt[%d]"), ptBinOm));

            peakFnc_OmC = setFitParameters(peakFnc_OmC, h_MassOmC_pt[ptBinOm], lMass_Om, ptBinOm + 1, kTRUE);
            resultParOmC_pt[ptBinOm] = fitResults(peakFnc_OmC, h_MassOmC_pt[ptBinOm], fitMinSig_Om, fitMaxSig_Om, &par_allintC[0], &par_allint_errorsC[0], binsMC, MCGenOmC, fisMC, kFALSE);
            resultParOmC_pt[ptBinOm]->SetName(TString::Format(("resultParOmC_pt[%d]"), ptBinOm));

            out->cd("h_MassOmm_pt");
            h_MassOmm_pt[ptBinOm]->SetTitle(TString::Format(("Invariant Mass #Omega^{-}: p_{T}<%.1f,%.1f>, Mult<%.1f,%.1f>"), ptbins_Om[ptBinOm], ptbins_Om[ptBinOm + 1], multbins_Om[0], multbins_Om[nmultbins_Om]));
            resultParOmm_pt[ptBinOm]->SetTitle(TString::Format(("Result Parameters #Omega^{-}: p_{T}<%.1f,%.1f>, Mult<%.1f,%.1f>"), ptbins_Om[ptBinOm], ptbins_Om[ptBinOm + 1], multbins_Om[0], multbins_Om[nmultbins_Om]));
            h_MassOmm_pt[ptBinOm]->Write();
            resultParOmm_pt[ptBinOm]->Write();

            out->cd("h_MassOmp_pt");
            h_MassOmp_pt[ptBinOm]->SetTitle(TString::Format(("Invariant Mass #Omega^{+}: p_{T}<%.1f,%.1f>, Mult<%.1f,%.1f>"), ptbins_Om[ptBinOm], ptbins_Om[ptBinOm + 1], multbins_Om[0], multbins_Om[nmultbins_Om]));
            resultParOmp_pt[ptBinOm]->SetTitle(TString::Format(("Result Parameters #Omega^{+}: p_{T}<%.1f,%.1f>, Mult<%.1f,%.1f>"), ptbins_Om[ptBinOm], ptbins_Om[ptBinOm + 1], multbins_Om[0], multbins_Om[nmultbins_Om]));
            h_MassOmp_pt[ptBinOm]->Write();
            resultParOmp_pt[ptBinOm]->Write();

            out->cd("h_MassOmC_pt");
            h_MassOmC_pt[ptBinOm]->SetTitle(TString::Format(("Invariant Mass #Omega: p_{T}<%.1f,%.1f>, Mult<%.1f,%.1f>"), ptbins_Om[ptBinOm], ptbins_Om[ptBinOm + 1], multbins_Om[0], multbins_Om[nmultbins_Om]));
            resultParOmC_pt[ptBinOm]->SetTitle(TString::Format(("Result Parameters #Omega: p_{T}<%.1f,%.1f>, Mult<%.1f,%.1f>"), ptbins_Om[ptBinOm], ptbins_Om[ptBinOm + 1], multbins_Om[0], multbins_Om[nmultbins_Om]));
            h_MassOmC_pt[ptBinOm]->Write();
            resultParOmC_pt[ptBinOm]->Write();

            for (Int_t multBinOm = 0; multBinOm < nmultbins_Om; multBinOm++)
            {
                /// fully differential case
                if (fisMC)
                {
                    binsMC[2] = multBinOm + 1;
                    binsMC[3] = multBinOm + 1;
                }

                h_MassOmm_pt_mult[ptBinOm][multBinOm] = (TH1D *)Omm->ProjectionY(TString::Format(("h_MassOmm_pt_mult[%d][%d]"), ptBinOm, multBinOm), ptBinOm + 1, ptBinOm + 1, multBinOm + 1, multBinOm + 1, "e");

                h_MassOmp_pt_mult[ptBinOm][multBinOm] = (TH1D *)Omp->ProjectionY(TString::Format(("h_MassOmp_pt_mult[%d][%d]"), ptBinOm, multBinOm), ptBinOm + 1, ptBinOm + 1, multBinOm + 1, multBinOm + 1, "e");

                h_MassOmC_pt_mult[ptBinOm][multBinOm] = (TH1D *)h_MassOmp_pt_mult[ptBinOm][multBinOm]->Clone(TString::Format(("h_MassOmC_pt_mult[%d][%d]"), ptBinOm, multBinOm));
                h_MassOmC_pt_mult[ptBinOm][multBinOm]->Add(h_MassOmm_pt_mult[ptBinOm][multBinOm]);

                peakFnc_OmP = setFitParameters(peakFnc_OmP, h_MassOmp_pt_mult[ptBinOm][multBinOm], lMass_Om, ptBinOm + 1, kFALSE);
                resultParOmp_pt_mult[ptBinOm][multBinOm] = fitResults(peakFnc_OmP, h_MassOmp_pt_mult[ptBinOm][multBinOm], fitMinSig_Om, fitMaxSig_Om, &par_allintP[0], &par_allint_errorsP[0], binsMC, MCGenOmp, fisMC, kFALSE);
                resultParOmp_pt_mult[ptBinOm][multBinOm]->SetName(TString::Format(("resultParOmp_pt_mult[%d][%d]"), ptBinOm, multBinOm));

                peakFnc_OmM = setFitParameters(peakFnc_OmM, h_MassOmm_pt_mult[ptBinOm][multBinOm], lMass_Om, ptBinOm + 1, kFALSE);
                resultParOmm_pt_mult[ptBinOm][multBinOm] = fitResults(peakFnc_OmM, h_MassOmm_pt_mult[ptBinOm][multBinOm], fitMinSig_Om, fitMaxSig_Om, &par_allintM[0], &par_allint_errorsM[0], binsMC, MCGenOmm, fisMC, kFALSE);
                resultParOmm_pt_mult[ptBinOm][multBinOm]->SetName(TString::Format(("resultParOmm_pt_mult[%d][%d]"), ptBinOm, multBinOm));

                peakFnc_OmC = setFitParameters(peakFnc_OmC, h_MassOmC_pt_mult[ptBinOm][multBinOm], lMass_Om, ptBinOm + 1, kFALSE);
                resultParOmC_pt_mult[ptBinOm][multBinOm] = fitResults(peakFnc_OmC, h_MassOmC_pt_mult[ptBinOm][multBinOm], fitMinSig_Om, fitMaxSig_Om, &par_allintC[0], &par_allint_errorsC[0], binsMC, MCGenOmC, fisMC, kFALSE);
                resultParOmC_pt_mult[ptBinOm][multBinOm]->SetName(TString::Format(("resultParOmC_pt_mult[%d][%d]"), ptBinOm, multBinOm));

                out->cd("h_MassOmm_pt_mult");
                h_MassOmm_pt_mult[ptBinOm][multBinOm]->SetTitle(TString::Format(("Invariant Mass #Omega^{-}: p_{T}<%.1f,%.1f>, Mult<%.1f,%.1f>"), ptbins_Om[ptBinOm], ptbins_Om[ptBinOm + 1], multbins_Om[multBinOm], multbins_Om[multBinOm + 1]));
                resultParOmm_pt_mult[ptBinOm][multBinOm]->SetTitle(TString::Format(("Result Parameters #Omega^{-}: p_{T}<%.1f,%.1f>, Mult<%.1f,%.1f>"), ptbins_Om[ptBinOm], ptbins_Om[ptBinOm + 1], multbins_Om[multBinOm], multbins_Om[multBinOm + 1]));
                h_MassOmm_pt_mult[ptBinOm][multBinOm]->Write();
                resultParOmm_pt_mult[ptBinOm][multBinOm]->Write();

                out->cd("h_MassOmp_pt_mult");
                h_MassOmp_pt_mult[ptBinOm][multBinOm]->SetTitle(TString::Format(("Invariant Mass #Omega^{+}: p_{T}<%.1f,%.1f>, Mult<%.1f,%.1f>"), ptbins_Om[ptBinOm], ptbins_Om[ptBinOm + 1], multbins_Om[multBinOm], multbins_Om[multBinOm + 1]));
                resultParOmp_pt_mult[ptBinOm][multBinOm]->SetTitle(TString::Format(("Result Parameters #Omega^{+}: p_{T}<%.1f,%.1f>, Mult<%.1f,%.1f>"), ptbins_Om[ptBinOm], ptbins_Om[ptBinOm + 1], multbins_Om[multBinOm], multbins_Om[multBinOm + 1]));
                h_MassOmp_pt_mult[ptBinOm][multBinOm]->Write();
                resultParOmp_pt_mult[ptBinOm][multBinOm]->Write();

                out->cd("h_MassOmC_pt_mult");
                h_MassOmC_pt_mult[ptBinOm][multBinOm]->SetTitle(TString::Format(("Invariant Mass #Omega: p_{T}<%.1f,%.1f>, Mult<%.1f,%.1f>"), ptbins_Om[ptBinOm], ptbins_Om[ptBinOm + 1], multbins_Om[multBinOm], multbins_Om[multBinOm + 1]));
                resultParOmC_pt_mult[ptBinOm][multBinOm]->SetTitle(TString::Format(("Result Parameters #Omega: p_{T}<%.1f,%.1f>, Mult<%.1f,%.1f>"), ptbins_Om[ptBinOm], ptbins_Om[ptBinOm + 1], multbins_Om[multBinOm], multbins_Om[multBinOm + 1]));
                h_MassOmC_pt_mult[ptBinOm][multBinOm]->Write();
                resultParOmC_pt_mult[ptBinOm][multBinOm]->Write();
            }
        }
    }

    /// temporarily commented out:
    // if (!gROOT->IsBatch())
    // {
    //     TCanvas *c = new TCanvas("XiOmegaFitting", "Cascade Analysis: Xi and Omega");
    //     c->Divide(2, 3);
    //     c->cd(1);
    //     h_MassXip->DrawCopy();
    //     c->cd(2);
    //     h_MassXim->DrawCopy();
    //     c->cd(3);
    //     h_MassOmp->DrawCopy();
    //     c->cd(4);
    //     h_MassOmm->DrawCopy();
    //     c->cd(5);
    //     h_MassXiC->DrawCopy();
    //     c->cd(6);
    //     h_MassOmC->DrawCopy();
    // }

    /// end of temporarily commented out

    out->Close();
    f->Close();

    // delete f, out;
    return 0;
}

TF1 *setFitParameters(TF1 *peakFnc, TH1 *peak, Double_t mass, Int_t ptBin, Bool_t multInt)
{
    Double_t zeroes[6] = {0.};
    peakFnc->SetParameters(zeroes);
    peakFnc->SetParErrors(zeroes);
    for (Int_t i = 0; i < 6; i++)
    {
        peakFnc->ReleaseParameter(i);
    }

    peakFnc->SetParNames("Amplitude", "Mean", "#sigma");
    peakFnc->SetParameter(1, mass);
    // peakFnc->SetParameter(2, 0.0018);
    // peakFnc->SetParLimits(2, 0.001, 0.01);
    if (multInt && ptBin == -1)
    {
        /// case when all bins are integrated: pt+mult integrated
        peakFnc->SetParameter(0, peak->GetMaximum() * 0.9);
        // peakFnc->SetParameter(2, 0.0025); // used by Emily - not so good...
        peakFnc->SetParameter(2, 0.0018); // better...
        peakFnc->SetParameter(3, 0);
        peakFnc->SetParameter(4, 0);
        peakFnc->SetParameter(5, peak->GetMaximum() * 0.1);
        return peakFnc;
    }

    if (multInt && ptBin > 0)
    {
        /// mult integrated case
        peakFnc->SetParameter(1, mass);
        peakFnc->SetParLimits(1, mass - 0.003, mass + 0.003);
        peakFnc->SetParameter(2, 0.0018);
        // peakFnc->SetParameter(2, 0.0025);      // used by Emily - not so good...

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
    else
    {
        /// Fully differential case
        peakFnc->SetParameter(1, mass);
        peakFnc->SetParLimits(1, mass - 0.003, mass + 0.003);
        peakFnc->SetParameter(2, 0.0018);
        // peakFnc->SetParameter(2, 0.0025); // used by Emily - not so good...

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

TH1 *fitResults(TF1 *peakFnc, TH1 *peak, Double_t fitMinSig, Double_t fitMaxSig, Double_t *par_allint, Double_t *par_allint_errors, Int_t binsMC[], TH2 *MCgen, bool fisMC, Bool_t allInt, Bool_t options)
{
    // ROOT::Math::IntegratorOneDimOptions::SetDefaultAbsTolerance(1.E-2);
    // ROOT::Math::IntegratorOneDimOptions::SetDefaultRelTolerance(1.E-2);

    Double_t peakFitMass = 0.;
    Double_t peakFitSigma = 0.;
    Double_t bgReject = 0.;

    if (TMath::IsNaN(peak->GetEntries()))
    {
        return 0;
    }

    /// use different fit options:
    TString opt = "QMINLSR+ MULTITHREAD";
    // TString opt = "QMINLSR";

    // TString optBg = "QMNLSR";
    TString optBg = "QMINSR+ MULTITHREAD";
    if (options)
        opt = "QNLS";

    TFitResultPtr fFitResult;

    fFitResult = peak->Fit(peakFnc, opt.Data(), "", fitMinSig, fitMaxSig);

    peak->GetListOfFunctions()->Add(peakFnc);

    // Getting pol2 part from fit and set it to bgFnc for drawing and then
    // subtracting

    if (allInt)
    {
        // copy par_allInt from allInt case to another array to set bg params.
        peakFnc->GetParameters(par_allint);
        // par_allint_errors = peakFnc->GetParErrors();
        for (Int_t i = 0; i < 6; i++)
        {
            par_allint_errors[i] = peakFnc->GetParError(i);
        }
        Printf("%s: par_allint[] = %.2f, %f, %.2f, %.2f, %.2f, %.2f", peak->GetName(), par_allint[0], par_allint[1], par_allint[2], par_allint[3], par_allint[4], par_allint[5]);
        Printf("%s: par_allint_errors[] = %.2f, %.2f, %.2f, %.2f, %.2f, %.2f", peak->GetName(), par_allint_errors[0], par_allint_errors[1], par_allint_errors[2], par_allint_errors[3], par_allint_errors[4], par_allint_errors[5]);
    }

    // Printf("%s: par_allint[] = %.2f, %.2f, %.2f, %.2f, %.2f, %.2f", peak->GetName(), par_allint[0], par_allint[1], par_allint[2], par_allint[3], par_allint[4], par_allint[5]);
    // Printf("%s: par_allint_errors[] = %.2f, %.2f, %.2f, %.2f, %.2f, %.2f", peak->GetName(), par_allint_errors[0], par_allint_errors[1], par_allint_errors[2], par_allint_errors[3], par_allint_errors[4], par_allint_errors[5]);

    peakFitMass = peakFnc->GetParameter(1); // getting peak (mass) from peakFnc to be used for excluding bg
    peakFitSigma = TMath::Abs(peakFnc->GetParameter(2));
    bgReject = 4 * peakFitSigma; // 4*sigma region to be excluded from bg fit

    if (!fisMC)
    {
        TF1 *bgFnc = new TF1("Pol2Exclude", Pol2Exclude, peak->GetXaxis()->GetBinCenter(peak->FindFirstBinAbove(0.)), peak->GetXaxis()->GetBinCenter(peak->FindLastBinAbove(0.)), 6);
        bgFnc->SetParNames("p[0]", "p[1]", "p[2]", "peakFitMass", "peakFitSigma", "rejectBgPoint");
        bgFnc->FixParameter(3, peakFitMass); // "mass = mean of fit"
        bgFnc->FixParameter(4, bgReject);    // "4*sigma"
        bgFnc->SetLineColor(kGreen);
        bgFnc->SetFillColor(kYellow);
        bgFnc->SetFillStyle(3009);

        bgFnc->SetParameter(0, par_allint[3]);
        bgFnc->SetParameter(1, par_allint[4]);
        bgFnc->SetParameter(2, par_allint[5]);
        bgFnc->SetParError(0, par_allint_errors[3]);
        bgFnc->SetParError(1, par_allint_errors[4]);
        bgFnc->SetParError(2, par_allint_errors[5]);

        // fitMinBg = peakFitMass - 10 * peakFitSigma;
        // fitMaxBg = peakFitMass + 10 * peakFitSigma;
        // peak->Fit(bgFnc, optBg.Data(), "", fitMinBg, fitMaxBg);
        bgFnc->FixParameter(5, 1); // "reject" (set 1 or 0)
        peak->Fit(bgFnc, optBg.Data(), "", peak->GetXaxis()->GetBinCenter(peak->FindFirstBinAbove(0.)), peak->GetXaxis()->GetBinCenter(peak->FindLastBinAbove(0.)));
        bgFnc->FixParameter(5, 0); // reject = off (set 1 or 0)

        // add fit and background function to histogram so it is automatically drawn
        // with hist
        peak->GetListOfFunctions()->Add(bgFnc);
    }

    /// getting min/max for deciding signal's fit range for integral calculation
    Double_t minInt = peakFitMass - bgReject; /// mean - 4*sigma
    Double_t maxInt = peakFitMass + bgReject; /// mean + 4*sigma
    // Printf("%s: MinInt, MaxInt = %f, %f ", peak->GetName(), minInt, maxInt);

    TH1 *resultParams = fillParams(peakFnc, peak, fFitResult, minInt, maxInt, MCgen, binsMC, fisMC);
    peak->SetOption("X0E1");
    return resultParams;
}
TH1 *fillParams(TF1 *sigBgFnc, TH1 *peak, TFitResultPtr fFitResult, Double_t minInt, Double_t maxInt, TH2 *MCgen, Int_t binsMC[], bool fisMC)
{
    Double_t val, err, gen, genErr;
    Double_t val_intBC, err_intBC, val_intFF, err_intFF, eff, effErr;
    Int_t nPar = sigBgFnc->GetNpar();
    const Int_t nBins = 8 + nPar;
    Double_t par[nPar];
    Double_t parErr[nPar];
    sigBgFnc->GetParameters(par);
    for (Int_t i = 0; i < nPar; i++)
    {
        parErr[i] = sigBgFnc->GetParError(i);
    }

    TH1 *resultParams = new TH1D("resultParams", "Result parameters", nBins, 0, nBins);
    Int_t iBin = 1;

    GetYieldBinCounting(peak, fFitResult, minInt, maxInt, fisMC, val, err);
    val_intBC = val;
    err_intBC = err;
    if (TMath::IsNaN(val_intBC) || TMath::IsNaN(err_intBC))
    {
        Printf("%s: GetYieldBinCounting is NaN!", peak->GetName());
        delete resultParams;
        return nullptr;
    }

    resultParams->SetBinContent(iBin, val_intBC);
    resultParams->SetBinError(iBin, err_intBC);
    resultParams->GetXaxis()->SetBinLabel(iBin, "IntBC");
    iBin++;

    val = 0;
    err = 0;

    GetYieldFitFunction(peak, fFitResult, minInt, maxInt, fisMC, val, err);
    val_intFF = val;
    err_intFF = err;
    if (TMath::IsNaN(val_intFF) || TMath::IsNaN(err_intFF))
    {
        Printf("%s: GetYieldFitFunction is NaN!", peak->GetName());
        delete resultParams;
        return nullptr;
    }
    resultParams->SetBinContent(iBin, val_intFF);
    resultParams->SetBinError(iBin, err_intFF);
    resultParams->GetXaxis()->SetBinLabel(iBin, "IntFF");
    iBin++;

    if (fisMC)
    {
        gen = GetGeneratedParticles(MCgen, binsMC, genErr);

        resultParams->SetBinContent(iBin, gen);
        resultParams->SetBinError(iBin, genErr);
        resultParams->GetXaxis()->SetBinLabel(iBin, "IntGen");
        iBin++;
        // Printf("%s: IntGen = %f; IntGenErr = %f", peak->GetName(), gen, genErr);

        // Printf("[%f, %f]: genErr = %f", minInt, maxInt, genErr);
        if (gen != 0 && genErr != 0)
        {
            eff = val_intBC / gen;
            // effErr = err_intBC / genErr; /// ask if it should be computed like this?

            effErr = ErrorInRatio(val_intBC, err_intBC, gen, genErr);
            // Printf("%s: Efficiency = %f; Error = %f", peak->GetName(), eff, effErr);
            resultParams->SetBinContent(iBin, eff);
            // resultParams->SetBinError(iBin, err_intBC);
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
    return resultParams;
}

void GetYieldBinCounting(TH1 *h, TFitResultPtr fFitResult, Double_t minInt, Double_t maxInt, bool fisMC, Double_t &val,
                         Double_t &err, Double_t eps = 1e-6)
{

    if (!h)
        return;

    Int_t bin_minSig = h->GetXaxis()->FindBin(minInt);
    Int_t bin_maxSig = h->GetXaxis()->FindBin(maxInt);

    Double_t histWidth = h->GetXaxis()->GetBinWidth(1); /// different width for each bin -> bin widths are the same!!!

    val = h->IntegralAndError(bin_minSig, bin_maxSig, err);

    if (!fisMC)
    {
        TF1 *bgFnc = (TF1 *)h->GetListOfFunctions()->At(1);

        Double_t bg = bgFnc->Integral(minInt, maxInt, eps);

        // TODO Verify it
        Double_t bgErr = bgFnc->IntegralError(minInt, maxInt, fFitResult->GetParams(),
                                              fFitResult->GetCovarianceMatrix().GetMatrixArray(), eps);

        bg /= histWidth;
        bgErr /= histWidth;
        ////
        if (TMath::Abs(bg) < 1) // test for bg value == 0.00...
        {
            Printf("IntBC: %s: [%f, %f], bins:[%d, %d]: bg = %f, val(signal) = %f", h->GetName(), minInt, maxInt, bin_minSig, bin_maxSig, bg, val);
            // bgFnc->Print("V");
        }
        if (val > 30) /// why?
        {
            val -= bg;
            err = TMath::Sqrt(TMath::Power(err, 2) + TMath::Power(bgErr, 2));
        }
    }
}

void GetYieldFitFunction(TH1 *h, TFitResultPtr fFitResult, Double_t minInt, Double_t maxInt, bool fisMC, Double_t &val,
                         Double_t &err, Double_t eps = 1e-6)
{

    if (!h)
        return;

    Double_t histWidth = h->GetXaxis()->GetBinWidth(1);

    TF1 *sigBgFnc = (TF1 *)h->GetListOfFunctions()->At(0);

    val = sigBgFnc->Integral(minInt, maxInt, eps);
    err = sigBgFnc->IntegralError(minInt, maxInt, fFitResult->GetParams(), fFitResult->GetCovarianceMatrix().GetMatrixArray(),
                                  eps);

    val /= histWidth;
    err /= histWidth;
    if (!fisMC)
    {
        TF1 *bgFnc = (TF1 *)h->GetListOfFunctions()->At(1);
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
}
