#if !defined(__CINT__) || (defined(__MAKECINT__))
#include <iostream>
#include "TH1.h"
#include "TVirtualFitter.h"
#include "TMath.h"
#include "TFile.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TROOT.h"
#include "TRandom.h"
#include "TH3.h"
#include "TH2.h"
#include "TGraphErrors.h"
#include "TLatex.h"

#endif
using namespace std;
#include <TString.h>
/* definition of the fields in the histogram returned */
enum EValue_t
{
    kYield = 1,
    kYieldStat,
    kYieldSysHi,
    kYieldSysLo,
    kMean,
    kMeanStat,
    kMeanSysHi,
    kMeanSysLo,
    kExtra,
    kChi

};

TH1 *YieldMean(TH1 *hstat, TH1 *hsys, TF1 *f = NULL, Double_t minfit = 0.4, Double_t maxfit = 8., TString logfilename = "log.root", Double_t min = 0.8, Double_t max = 8., Double_t loprecision = 0.01, Double_t hiprecision = 0.1, Option_t *opt = "0q");
void YieldMean_IntegralMean(TH1 *hdata, TH1 *hlo, TH1 *hhi, Double_t &integral, Double_t &mean, Double_t &extra, Bool_t printinfo = kFALSE);
TH1 *YieldMean_LowExtrapolationHisto(TH1 *h, TF1 *f, Double_t min, Double_t binwidth = 0.01);
TH1 *YieldMean_HighExtrapolationHisto(TH1 *h, TF1 *f, Double_t max, Double_t binwidth = 0.1);
TH1 *YieldMean_ReturnRandom(TH1 *hin);
TH1 *YieldMean_ReturnCoherentRandom(TH1 *hin);
TH1 *YieldMean_ReturnExtremeHisto(TH1 *hin, Float_t sign = 1.);
TH1 *YieldMean_ReturnExtremeHardHisto(TH1 *hin);
TH1 *YieldMean_ReturnExtremeSoftHisto(TH1 *hin);
TH1 *YieldMean_ReturnExtremeLowHisto(TH1 *hin);
TH1 *YieldMean_ReturnExtremeHighHisto(TH1 *hin);
TF1 *LevyTsallis(const Char_t *name, Double_t mass, Double_t n = 5., Double_t C = 0.1, Double_t norm = 1.);
// void Init(Int_t lClassCode,
//           TString lWhichEnergyEstimator,
//           TString lWhichMultEstimator,
//           Bool_t lDoMB,
//           TString fWhichParticle,
//           vector<Double_t> &percentileMult_low,
//           vector<Double_t> &percentileMult_high,
//           vector<Double_t> &percentileEnergy_low,
//           vector<Double_t> &percentileEnergy_high,
//           vector<Double_t> &nch,
//           vector<Double_t> &nchErr);

void CalculateYields()
{

    const int ntrials = 500;
    TH1D *hout[ntrials];
    TF1 *LevyTsallisfunc[ntrials];
    double minfit = 0.8;
    double maxfit = 8.0;
    double mass = 1.67245;
    TFile *fdef = TFile::Open("CorrSpectra/YieldEffCorrOmega_Mult0100_Default.root");
    TH1D *hdef = (TH1D *)fdef->Get("histoYieldCorr");
    TFile *fvaried[ntrials];
    TH1D *hvaried[ntrials];
    double YieldVar[ntrials];
    double YieldVarStat[ntrials];

    TH1D *hdist = new TH1D("hdist", "1 - (dNdy_{var}/dNdy_{def})", 50, -0.5, 0.5);

    TF1 *LevyTsallisfuncdef = LevyTsallis(Form("LevyTsallisfuncdef"), mass);
    TH1D *houtdef = (TH1D *)YieldMean(hdef, hdef, LevyTsallisfuncdef, minfit, maxfit, Form("logfiles/logfiledef.root"));
    // Get yield
    double YieldDef = houtdef->GetBinContent(kYield);
    double YieldVarStatDef = houtdef->GetBinContent(kYieldStat);

    for (int i = 0; i < ntrials; i++)
    {
        // Fit function
        LevyTsallisfunc[i] = LevyTsallis(Form("LevyTsallisfunc%i", i), mass);

        fvaried[i] = TFile::Open(Form("CorrSpectra/YieldEffCorrOmega_Mult0100_%i.root", i));
        hvaried[i] = (TH1D *)fvaried[i]->Get("histoYieldCorr");

        hout[i] = (TH1D *)YieldMean(hvaried[i], hvaried[i], LevyTsallisfunc[i], minfit, maxfit, Form("logfiles/logfile%i.root", i));

        // Get yield
        YieldVar[i] = hout[i]->GetBinContent(kYield);
        YieldVarStat[i] = hout[i]->GetBinContent(kYieldStat);

        hdist->Fill(1 - YieldVar[i] / YieldDef);

    } // end loop over classes

    TFile *write = new TFile("pTuncorr_noextrap.root", "RECREATE");
    hdist->Write();

    TF1 *fgaus = new TF1(Form("fgaus"), "gaus", -1, +1);
    hdist->Fit(fgaus, "q");

    cout << "Mean: " << fgaus->GetParameter(1) << endl;
    cout << "Sigma: " << fgaus->GetParameter(2) << endl;

    TLatex *ltx = new TLatex();
    ltx->SetTextSize(0.04);
    ltx->SetTextColor(kRed);

    TCanvas *c = new TCanvas("c", "c", 1000, 900);
    c->cd(1);
    c->SetTopMargin(0.05);
    c->SetBottomMargin(0.15);
    c->SetLeftMargin(0.15);
    c->SetRightMargin(0.05);
    hdist->SetMarkerStyle(20);
    hdist->SetMarkerColor(kBlack);
    hdist->SetMarkerSize(1);
    hdist->SetStats(0);
    hdist->SetLineColor(kBlack);
    hdist->GetYaxis()->SetTitle("Counts");
    hdist->GetXaxis()->SetTitle("1 - (dN/dy_{var}/dN/dy_{def})");
    hdist->GetXaxis()->SetTitleSize(0.05);
    hdist->GetXaxis()->SetTitleOffset(1.1);
    hdist->GetYaxis()->SetTitleSize(0.05);
    hdist->GetYaxis()->SetTitleOffset(1.1);
    hdist->SetTitle("");
    hdist->Draw("EP");
    ltx->DrawLatexNDC(0.2, 0.8, Form("#mu = %.3f", fgaus->GetParameter(1)));
    ltx->DrawLatexNDC(0.2, 0.7, Form("#sigma = %.3f", fgaus->GetParameter(2)));
    ltx->DrawLatexNDC(0.2, 0.6, Form("#chi^{2}/ndf = %.1f", fgaus->GetChisquare() / fgaus->GetNDF()));
    c->SaveAs("images/pTuncorr_noextrap.png");
    c->SaveAs("images/pTuncorr_noextrap.pdf");
}

TH1 *YieldMean(TH1 *hstat, TH1 *hsys, TF1 *f = NULL, Double_t minfit = .4, Double_t maxfit = 8., TString logfilename = "log.root", Double_t min = 0.8, Double_t max = 8., Double_t loprecision = 0.01, Double_t hiprecision = 0.1, Option_t *opt = "0q")
{
    if (maxfit > max)
        max = maxfit;
    if (minfit < min)
        min = minfit;

    /* set many iterations when fitting the data so we don't
       stop minimization with MAX_CALLS */
    TVirtualFitter::SetMaxIterations(1000000);

    /* create output histo */
    Double_t integral, mean, extra;
    TH1 *hout = new TH1D("hout", "", 9, 0, 9);
    TH1 *hlo, *hhi;

    /* create histo with stat+sys errors */
    TH1 *htot = (TH1 *)hstat->Clone(Form("%sfittedwith%s", hstat->GetName(), "Func"));
    for (Int_t ibin = 0; ibin < htot->GetNbinsX(); ibin++)
    {
        // htot->SetBinError(ibin + 1, TMath::Sqrt(hsys->GetBinError(ibin + 1) * hsys->GetBinError(ibin + 1) + hstat->GetBinError(ibin + 1) * hstat->GetBinError(ibin + 1)));
    }

    /*
     *   measure the central value
     */
    Int_t fitres;
    Int_t trials = 0;
    trials = 0;
    do
    {
        fitres = htot->Fit(f, opt, "", minfit, maxfit);
        Printf("Trial: %d", trials++);
        if (trials > 10)
        {
            Printf("FIT DOES NOT CONVERGE IN LINE %d", __LINE__);
            break;
        }
    } while (fitres != 0);
    TFile *filewithfits = TFile::Open(logfilename.Data(), "UPDATE");
    htot->Write();
    f->SetName(Form("Levyfitto%s", hstat->GetName()));
    f->Write();
    filewithfits->Close();
    delete filewithfits;

    cout << " Fit sys+stat for " << f->GetName() << endl;
    cout << "NDF=" << f->GetNDF() << " Chi^2=" << f->GetChisquare() << " Chi^2/NDF=" << f->GetChisquare() / f->GetNDF() << endl;
    hout->SetBinContent(kChi, f->GetChisquare() / f->GetNDF());

    hlo = YieldMean_LowExtrapolationHisto(htot, f, min, loprecision);
    hhi = YieldMean_HighExtrapolationHisto(htot, f, max, hiprecision);
    YieldMean_IntegralMean(htot, hlo, hhi, integral, mean, extra, kTRUE);
    hout->SetBinContent(kYield, integral);
    hout->SetBinContent(kMean, mean);
    hout->SetBinContent(kExtra, extra);

    /*
     * STATISTICS
     */

    TCanvas *cCanvasStat = new TCanvas("cCanvasStat");
    cCanvasStat->Divide(2, 1);

    /*
     * measure statistical error
     */

    /* fit with stat error */
    trials = 0;
    do
    {
        fitres = hstat->Fit(f, opt, "", minfit, maxfit);
        Printf("Trial: %d", trials++);
        if (trials > 10)
        {
            Printf("FIT DOES NOT CONVERGE IN LINE %d", __LINE__);
            break;
        }
    } while (fitres != 0);
    hlo = YieldMean_LowExtrapolationHisto(hstat, f, min, loprecision);
    hhi = YieldMean_HighExtrapolationHisto(hstat, f, max, hiprecision);

    /* random generation with integration (coarse) */
    TH1 *hIntegral_tmp = new TH1F("hIntegral_tmp", "", 1000, 0.75 * integral, 1.25 * integral);
    TH1 *hMean_tmp = new TH1F("hMean_tmp", "", 1000, 0.75 * mean, 1.25 * mean);
    for (Int_t irnd = 0; irnd < 100; irnd++)
    {
        /* get random histogram */
        TH1 *hrnd = YieldMean_ReturnRandom(hstat);
        /* fit */
        TH1 *hrndlo = YieldMean_ReturnCoherentRandom(hlo);
        TH1 *hrndhi = YieldMean_ReturnCoherentRandom(hhi);
        /* integrate */
        YieldMean_IntegralMean(hrnd, hrndlo, hrndhi, integral, mean, extra);
        hIntegral_tmp->Fill(integral);
        hMean_tmp->Fill(mean);
        delete hrnd;
        delete hrndlo;
        delete hrndhi;
    }
    /* random generation with integration (fine) */
    TH1 *hIntegral = new TH1F("hIntegral", "", 100,
                              hIntegral_tmp->GetMean() - 10. * hIntegral_tmp->GetRMS(),
                              hIntegral_tmp->GetMean() + 10. * hIntegral_tmp->GetRMS());
    TH1 *hMean = new TH1F("hMean", "", 100,
                          hMean_tmp->GetMean() - 10. * hMean_tmp->GetRMS(),
                          hMean_tmp->GetMean() + 10. * hMean_tmp->GetRMS());
    for (Int_t irnd = 0; irnd < 1000; irnd++)
    {
        /* get random histogram */
        TH1 *hrnd = YieldMean_ReturnRandom(hstat);
        /* fit */
        TH1 *hrndlo = YieldMean_ReturnCoherentRandom(hlo);
        TH1 *hrndhi = YieldMean_ReturnCoherentRandom(hhi);
        /* integrate */
        YieldMean_IntegralMean(hrnd, hrndlo, hrndhi, integral, mean, extra);
        hIntegral->Fill(integral);
        hMean->Fill(mean);
        delete hrnd;
        delete hrndlo;
        delete hrndhi;
    }
    TF1 *gaus = (TF1 *)gROOT->GetFunction("gaus");

    cCanvasStat->cd(1);
    hIntegral->Fit(gaus, "q");
    integral = hout->GetBinContent(kYield) * gaus->GetParameter(2) / gaus->GetParameter(1);
    hout->SetBinContent(kYieldStat, integral);

    cCanvasStat->cd(2);
    hMean->Fit(gaus, "q");
    mean = hout->GetBinContent(kMean) * gaus->GetParameter(2) / gaus->GetParameter(1);
    hout->SetBinContent(kMeanStat, mean);

    /*
     * SYSTEMATICS
     */

    TCanvas *cCanvasSys = new TCanvas("cCanvasYieldSys");
    cCanvasSys->Divide(2, 1);
    cCanvasSys->cd(1)->DrawFrame(min, 1.e-3, max, 1.e3);
    hsys->SetMarkerStyle(20);
    hsys->SetMarkerColor(1);
    hsys->SetMarkerSize(1);
    hsys->Draw("same");
    cCanvasSys->cd(2)->DrawFrame(min, 1.e-3, max, 1.e3);
    hsys->Draw("same");

    /*
     * systematic error high
     */

    TH1 *hhigh = YieldMean_ReturnExtremeHighHisto(hsys);
    trials = 0;
    do
    {
        fitres = hhigh->Fit(f, opt, "", minfit, maxfit);
        Printf("Trial: %d", trials++);
        if (trials > 10)
        {
            Printf("FIT DOES NOT CONVERGE IN LINE %d", __LINE__);
            break;
        }
    } while (fitres != 0);
    hlo = YieldMean_LowExtrapolationHisto(hhigh, f, min, loprecision);
    hhi = YieldMean_HighExtrapolationHisto(hhigh, f, max, hiprecision);
    YieldMean_IntegralMean(hhigh, hlo, hhi, integral, mean, extra);
    integral = TMath::Abs(integral - hout->GetBinContent(kYield));
    hout->SetBinContent(kYieldSysHi, integral);

    cCanvasSys->cd(1);
    f->SetLineColor(2);
    f->DrawCopy("same");

    /*
     * systematic error hard
     */

    TH1 *hhard = YieldMean_ReturnExtremeHardHisto(hsys);
    trials = 0;
    do
    {
        fitres = hhard->Fit(f, opt, "", minfit, maxfit);
        Printf("Trial: %d", trials++);
        if (trials > 10)
        {
            Printf("FIT DOES NOT CONVERGE IN LINE %d", __LINE__);
            break;
        }
    } while (fitres != 0);
    hlo = YieldMean_LowExtrapolationHisto(hhard, f, min, loprecision);
    hhi = YieldMean_HighExtrapolationHisto(hhard, f, max, hiprecision);
    YieldMean_IntegralMean(hhard, hlo, hhi, integral, mean, extra);
    mean = TMath::Abs(mean - hout->GetBinContent(kMean));
    hout->SetBinContent(kMeanSysHi, mean);

    cCanvasSys->cd(2);
    f->SetLineColor(2);
    f->DrawCopy("same");

    /*
     * systematic error low
     */

    TH1 *hlow = YieldMean_ReturnExtremeLowHisto(hsys);
    trials = 0;
    do
    {
        fitres = hlow->Fit(f, opt, "", minfit, maxfit);
        Printf("Trial: %d", trials++);
        if (trials > 10)
        {
            Printf("FIT DOES NOT CONVERGE IN LINE %d", __LINE__);
            break;
        }
    } while (fitres != 0);
    hlo = YieldMean_LowExtrapolationHisto(hlow, f, min, loprecision);
    hhi = YieldMean_HighExtrapolationHisto(hlow, f, max, hiprecision);
    YieldMean_IntegralMean(hlow, hlo, hhi, integral, mean, extra);
    integral = TMath::Abs(integral - hout->GetBinContent(kYield));
    hout->SetBinContent(kYieldSysLo, integral);

    cCanvasSys->cd(1);
    f->SetLineColor(4);
    f->DrawCopy("same");

    /*
     * systematic error soft
     */

    TH1 *hsoft = YieldMean_ReturnExtremeSoftHisto(hsys);
    trials = 0;
    do
    {
        fitres = hsoft->Fit(f, opt, "", minfit, maxfit);
        Printf("Trial: %d", trials++);
        if (trials > 10)
        {
            Printf("FIT DOES NOT CONVERGE IN LINE %d", __LINE__);
            break;
        }
    } while (fitres != 0);
    hlo = YieldMean_LowExtrapolationHisto(hsoft, f, min, loprecision);
    hhi = YieldMean_HighExtrapolationHisto(hsoft, f, max, hiprecision);
    YieldMean_IntegralMean(hsoft, hlo, hhi, integral, mean, extra);
    mean = TMath::Abs(mean - hout->GetBinContent(kMean));
    hout->SetBinContent(kMeanSysLo, mean);

    cCanvasSys->cd(2);
    f->SetLineColor(4);
    f->DrawCopy("same");

    return hout;
}

TH1 *YieldMean_LowExtrapolationHisto(TH1 *h, TF1 *f, Double_t min, Double_t binwidth)
{
    /* find lowest edge in histo */
    Int_t binlo;
    Double_t lo;
    for (Int_t ibin = 1; ibin < h->GetNbinsX() + 1; ibin++)
    {
        if (h->GetBinContent(ibin) != 0.)
        {
            binlo = ibin;
            lo = h->GetBinLowEdge(ibin);
            break;
        }
    }

    Int_t nbins = (lo - min) / binwidth;
    if (nbins < 1)
        return 0x0;
    TH1 *hlo = new TH1F("hlo", "", nbins, min, lo);

    /* integrate function in histogram bins */
    Double_t cont, err, width;
    for (Int_t ibin = 0; ibin < hlo->GetNbinsX(); ibin++)
    {
        width = hlo->GetBinWidth(ibin + 1);
        cont = f->Integral(hlo->GetBinLowEdge(ibin + 1), hlo->GetBinLowEdge(ibin + 2)); //, (Double_t *)0, 1.e-6);
        err = f->IntegralError(hlo->GetBinLowEdge(ibin + 1), hlo->GetBinLowEdge(ibin + 2), (Double_t *)0, (Double_t *)0, 1.e-6);
        hlo->SetBinContent(ibin + 1, cont / width);
        hlo->SetBinError(ibin + 1, err / width);
    }

    return hlo;
}

TH1 *YieldMean_HighExtrapolationHisto(TH1 *h, TF1 *f, Double_t max, Double_t binwidth)
{
    /* find highest edge in histo */
    Int_t binhi;
    Double_t hi;
    for (Int_t ibin = h->GetNbinsX(); ibin > 0; ibin--)
    {
        if (h->GetBinContent(ibin) != 0.)
        {
            binhi = ibin + 1;
            hi = h->GetBinLowEdge(ibin + 1);
            break;
        }
    }
    if (max < hi)
    {
        Printf("Warning! You should probably set a higher max value (Max = %f, hi = %f)", max, hi);
        return 0x0;
    }
    Int_t nbins = (max - hi) / binwidth;
    if (nbins < 1)
        return 0x0;
    TH1 *hhi = new TH1F("hhi", "", nbins, hi, max);

    /* integrate function in histogram bins */
    Double_t cont, err, width;
    for (Int_t ibin = 0; ibin < hhi->GetNbinsX(); ibin++)
    {
        width = hhi->GetBinWidth(ibin + 1);
        cont = f->Integral(hhi->GetBinLowEdge(ibin + 1), hhi->GetBinLowEdge(ibin + 2)); //, (Double_t *)0, 1.e-6);
        err = f->IntegralError(hhi->GetBinLowEdge(ibin + 1), hhi->GetBinLowEdge(ibin + 2), (Double_t *)0, (Double_t *)0, 1.e-6);
        hhi->SetBinContent(ibin + 1, cont / width);
        hhi->SetBinError(ibin + 1, err / width);
    }

    return hhi;
}

TH1 *YieldMean_ReturnRandom(TH1 *hin)
{
    TH1 *hout = (TH1 *)hin->Clone("hout");
    hout->Reset();
    Double_t cont, err;
    for (Int_t ibin = 0; ibin < hin->GetNbinsX(); ibin++)
    {
        if (hin->GetBinError(ibin + 1) <= 0.)
            continue;
        cont = hin->GetBinContent(ibin + 1);
        err = hin->GetBinError(ibin + 1);
        hout->SetBinContent(ibin + 1, gRandom->Gaus(cont, err));
        hout->SetBinError(ibin + 1, err);
    }
    return hout;
}

TH1 *YieldMean_ReturnCoherentRandom(TH1 *hin)
{
    if (!hin)
        return 0x0;
    TH1 *hout = (TH1 *)hin->Clone("hout");
    hout->Reset();
    Double_t cont, err, cohe;
    cohe = gRandom->Gaus(0., 1.);
    for (Int_t ibin = 0; ibin < hin->GetNbinsX(); ibin++)
    {
        if (hin->GetBinError(ibin + 1) <= 0.)
            continue;
        cont = hin->GetBinContent(ibin + 1);
        err = hin->GetBinError(ibin + 1);
        hout->SetBinContent(ibin + 1, cont + cohe * err);
        hout->SetBinError(ibin + 1, err);
    }
    return hout;
}

TH1 *YieldMean_ReturnExtremeHighHisto(TH1 *hin)
{
    TH1 *hout = (TH1 *)hin->Clone(Form("%s_extremehigh", hin->GetName()));
    for (Int_t ibin = 0; ibin < hin->GetNbinsX(); ibin++)
    {
        if (hin->GetBinError(ibin + 1) <= 0.)
            continue;
        Double_t val = hin->GetBinContent(ibin + 1);
        Double_t err = hin->GetBinError(ibin + 1);
        hout->SetBinContent(ibin + 1, val + err);
    }
    return hout;
}

TH1 *YieldMean_ReturnExtremeLowHisto(TH1 *hin)
{
    TH1 *hout = (TH1 *)hin->Clone(Form("%s_extremelow", hin->GetName()));
    for (Int_t ibin = 0; ibin < hin->GetNbinsX(); ibin++)
    {
        if (hin->GetBinError(ibin + 1) <= 0.)
            continue;
        Double_t val = hin->GetBinContent(ibin + 1);
        Double_t err = hin->GetBinError(ibin + 1);
        hout->SetBinContent(ibin + 1, val - err);
    }
    return hout;
}

TH1 *YieldMean_ReturnExtremeSoftHisto(TH1 *hin)
{
    return YieldMean_ReturnExtremeHisto(hin, -1.);
}

TH1 *YieldMean_ReturnExtremeHardHisto(TH1 *hin)
{
    return YieldMean_ReturnExtremeHisto(hin, 1.);
}

TH1 *YieldMean_ReturnExtremeHisto(TH1 *hin, Float_t sign)
{
    Double_t ptlow, pthigh;
    for (Int_t ibin = 0; ibin < hin->GetNbinsX(); ibin++)
    {
        if (hin->GetBinError(ibin + 1) <= 0.)
            continue;
        ptlow = hin->GetBinLowEdge(ibin + 1);
        break;
    }
    for (Int_t ibin = hin->GetNbinsX(); ibin >= 0; ibin--)
    {
        if (hin->GetBinError(ibin + 1) <= 0.)
            continue;
        pthigh = hin->GetBinLowEdge(ibin + 2);
        break;
    }

    Double_t mean = hin->GetMean();
    Double_t maxdiff = 0.;
    TH1 *hmax = NULL;
    for (Int_t inode = 0; inode < hin->GetNbinsX(); inode++)
    {

        Double_t ptnode = hin->GetBinCenter(inode + 1);
        TH1 *hout = (TH1 *)hin->Clone(Form("%s_extremehard", hin->GetName()));

        for (Int_t ibin = 0; ibin < hin->GetNbinsX(); ibin++)
        {
            if (hin->GetBinError(ibin + 1) <= 0.)
                continue;
            Double_t val = hin->GetBinContent(ibin + 1);
            Double_t err = hin->GetBinError(ibin + 1);
            Double_t cen = hin->GetBinCenter(ibin + 1);
            if (cen < ptnode)
                err *= -1. + (cen - ptlow) / (ptnode - ptlow);
            else
                err *= (cen - ptnode) / (pthigh - ptnode);

            hout->SetBinContent(ibin + 1, val + sign * err);
        }

        Double_t diff = TMath::Abs(mean - hout->GetMean());
        if (diff > maxdiff)
        {
            //      printf("found max at %f\n", ptnode);
            if (hmax)
                delete hmax;
            hmax = (TH1 *)hout->Clone("hmax");
            maxdiff = diff;
        }
        delete hout;
    }
    return hmax;
}

void YieldMean_IntegralMean(TH1 *hdata, TH1 *hlo, TH1 *hhi, Double_t &integral, Double_t &mean, Double_t &extra, Bool_t printinfo)
{

    /*
     * compute integrals
     */

    Double_t cont, err, width, cent;
    Double_t I = 0., IX = 0., Ierr = 0., IXerr = 0., Ilerr = 0., IXlerr = 0.;
    Double_t M = 0., Merr = 0., Mlerr = 0., C;
    Double_t E = 0;
    Double_t dataonly = 0.0;

    /* integrate the data */
    for (Int_t ibin = 0; ibin < hdata->GetNbinsX(); ibin++)
    {
        cent = hdata->GetBinCenter(ibin + 1);
        width = hdata->GetBinWidth(ibin + 1);
        cont = width * hdata->GetBinContent(ibin + 1);
        err = width * hdata->GetBinError(ibin + 1);
        if (err <= 0.)
            continue;
        I += cont;
        IX += cont * cent;
    }

    dataonly = I;
    /* integrate low */
    if (hlo)
    {
        for (Int_t ibin = 0; ibin < hlo->GetNbinsX(); ibin++)
        {
            cent = hlo->GetBinCenter(ibin + 1);
            width = hlo->GetBinWidth(ibin + 1);
            cont = width * hlo->GetBinContent(ibin + 1);
            err = width * hlo->GetBinError(ibin + 1);
            if (err <= 0.)
                continue;
            I += cont;
            IX += cont * cent;
            E += cont;
        }
    }
    /* integrate high */
    if (printinfo)
        cout << "low part data only = " << dataonly << " total = " << I << " ratio= " << dataonly / I << endl;
    if (hhi)
    {
        for (Int_t ibin = 0; ibin < hhi->GetNbinsX(); ibin++)
        {
            cent = hhi->GetBinCenter(ibin + 1);
            width = hhi->GetBinWidth(ibin + 1);
            cont = width * hhi->GetBinContent(ibin + 1);
            err = width * hhi->GetBinError(ibin + 1);
            if (err <= 0.)
                continue;
            I += cont;
            IX += cont * cent;
            E += cont;
        }
    }
    /* set values */
    integral = I;
    mean = IX / I;
    extra = E;
    if (printinfo)
        cout << "low+high data only = " << dataonly << " total = " << I << " ratio= " << dataonly / I << endl;
}

//---------------------------------------------------------------------
/**
 * @brief Calculates the Levy-Tsallis function for given parameters.
 *
 * This function computes the value of the Levy-Tsallis distribution, which is 
 * often used in high-energy physics to describe the transverse momentum (pt) 
 * distribution of particles.
 *
 * @param x Pointer to an array containing the transverse momentum (pt) value.
 * @param p Pointer to an array containing the parameters of the function:
 *          - p[0]: mass of the particle
 *          - p[1]: parameter n (shape parameter)
 *          - p[2]: parameter C (scale parameter)
 *          - p[3]: normalization factor
 * @return The value of the Levy-Tsallis function for the given pt and parameters.
 */
Double_t LevyTsallis_Func(const Double_t *x, const Double_t *p)
{
    /* dN/dpt */
    Double_t pt = x[0];
    Double_t mass = p[0];
    Double_t mt = TMath::Sqrt(pt * pt + mass * mass);
    Double_t n = p[1];
    Double_t C = p[2];
    Double_t norm = p[3];

    Double_t part1 = (n - 1.) * (n - 2.);
    Double_t part2 = n * C * (n * C + mass * (n - 2.));
    Double_t part3 = part1 / part2;
    Double_t part4 = 1. + (mt - mass) / n / C;
    Double_t part5 = TMath::Power(part4, -n);
    return pt * norm * part3 * part5;
}

TF1 *LevyTsallis(const Char_t *name, Double_t mass, Double_t n = 5., Double_t C = 0.1, Double_t norm = 1.)
{

    TF1 *fLevyTsallis = new TF1(name, LevyTsallis_Func, 0., 10., 4);
    fLevyTsallis->SetParameters(mass, n, C, norm);
    fLevyTsallis->SetParNames("mass", "n", "C", "norm");
    fLevyTsallis->FixParameter(0, mass);
    fLevyTsallis->SetParLimits(1, 1.e-3, 1.e3);
    fLevyTsallis->SetParLimits(2, 1.e-4, 1.e4); // 4
    fLevyTsallis->SetParLimits(3, 1.e-5, 1.e5); // 5
    return fLevyTsallis;
}
//------------------------------------------------------
