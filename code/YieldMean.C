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
#include "ExtraUtils.C"
#include <TFitter.h>
#include <TFitResultPtr.h>
#include <Math/MinimizerOptions.h>
#include <Math/IntegratorOptions.h>

#endif
using namespace std;
/* definition of the fields in the histogram returned -- moved to ExtraUtils */
// enum EValue_t
// {
//   kYield = 1,
//   kYieldStat,
//   kYieldSysHi,
//   kYieldSysLo,
//   kMean,
//   kMeanStat,
//   kMeanSysHi,
//   kMeanSysLo,
//   kExtra,
//   kChi
// };

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

TF1 *SetFitParams(TF1 *func, TH1 *hist, Bool_t usePrevious = kFALSE, Double_t *previousFitParams = nullptr);

/**
 * @brief Calculates yield, mean, and their statistical and systematic uncertainties
 *        from input histograms using fitting and extrapolation methods.
 *
 * This function fits a histogram with statistical and systematic errors, calculates
 * the yield (integral) and mean values, and estimates their statistical and systematic
 * uncertainties through multiple approaches including random sampling and variations
 * of the input histograms.
 *
 * The function performs:
 * 1. Fitting of the combined statistical and systematic data
 * 2. Computation of central values (yield, mean, extrapolation)
 * 3. Statistical error estimation through Monte Carlo random sampling
 * 4. Systematic error estimation by using extreme variations of input histograms
 *    (high, low, hard, and soft variations)
 * 5. Creation of diagnostic canvases showing the fits
 *
 * @param hstat Input histogram containing statistical data
 * @param hsys Input histogram containing systematic uncertainties
 * @param f Function to fit the histograms (must be pre-defined)
 * @param min Minimum range for integration/extrapolation
 * @param max Maximum range for integration/extrapolation
 * @param loprecision Precision for low extrapolation
 * @param hiprecision Precision for high extrapolation
 * @param opt ROOT fitting options
 * @param logfilename Output ROOT file name to store the fits
 * @param minfit Minimum range for fitting
 * @param maxfit Maximum range for fitting
 *
 * @return TH1 histogram with 9 bins containing:
 *         - bin 1: Yield (integral)
 *         - bin 2: Mean value
 *         - bin 3: Extra (extrapolated fraction)
 *         - bin 4: Statistical error on yield
 *         - bin 5: Statistical error on mean
 *         - bin 6: High systematic error on yield
 *         - bin 7: High systematic error on mean
 *         - bin 8: Low systematic error on yield
 *         - bin 9: Low systematic error on mean
 */
TH1 *YieldMean(TH1 *hstat, TH1 *hsys, TF1 *f = NULL, Double_t min = 0., Double_t max = 10., Double_t loprecision = 0.01, Double_t hiprecision = 0.1, Option_t *opt = "0q", TString logfilename = "log.root", Double_t minfit = 0.0, Double_t maxfit = 10.0)
{

  if (maxfit > max)
    max = maxfit;
  if (minfit < min)
    min = minfit;

  ////
  ROOT::Math::IntegratorOneDimOptions::SetDefaultIntegrator("Adaptive");
  ROOT::Math::IntegratorOneDimOptions::SetDefaultNPoints(6);
  ROOT::Math::IntegratorOneDimOptions::SetDefaultAbsTolerance(1E-2);
  ROOT::Math::IntegratorOneDimOptions::SetDefaultRelTolerance(1E-2);
  ROOT::Math::MinimizerOptions::SetDefaultMaxFunctionCalls(100000);
  ROOT::Math::MinimizerOptions::SetDefaultMinimizer("Minuit2");
  ROOT::Math::MinimizerOptions::SetDefaultStrategy(2);
  ROOT::Math::MinimizerOptions::SetDefaultPrintLevel(0); // Fit printing: -1 = no printing, 0 (minimal) to 3 (max)
  gStyle->SetOptFit(1111);
  ////

  /* set many iterations when fitting the data so we don't
     stop minimization with MAX_CALLS */
  TVirtualFitter::SetMaxIterations(1000000);

  /* create output histo */
  Double_t integral, mean, extra;
  TH1 *hout = new TH1D("hout", "", 10, 0, 10);
  hout->GetXaxis()->SetBinLabel(kYield, "Yield");
  hout->GetXaxis()->SetBinLabel(kYieldStat, "YieldStat");
  hout->GetXaxis()->SetBinLabel(kYieldSysHi, "YieldSysHi");
  hout->GetXaxis()->SetBinLabel(kYieldSysLo, "YieldSysLo");
  hout->GetXaxis()->SetBinLabel(kMean, "Mean");
  hout->GetXaxis()->SetBinLabel(kMeanStat, "MeanStat");
  hout->GetXaxis()->SetBinLabel(kMeanSysHi, "MeanSysHi");
  hout->GetXaxis()->SetBinLabel(kMeanSysLo, "MeanSysLo");
  hout->GetXaxis()->SetBinLabel(kExtra, "Extra");
  hout->GetXaxis()->SetBinLabel(kChi, "Chi2/NDF");
  TH1 *hlo, *hhi;
  TF1 *func = new TF1();

  /* create histo with stat+sys errors */
  TH1 *htot = (TH1 *)hstat->Clone(Form("%s_%s", hstat->GetName(), f->GetName()));
  for (Int_t ibin = 0; ibin < htot->GetNbinsX(); ibin++)
  {
    htot->SetBinError(ibin + 1, TMath::Sqrt(hsys->GetBinError(ibin + 1) * hsys->GetBinError(ibin + 1) + hstat->GetBinError(ibin + 1) * hstat->GetBinError(ibin + 1)));
  }

  /*
   *   measure the central value
   */
  Int_t fitres;
  Int_t trials = 0;
  trials = 0;
  f->Copy(*func);
  func = SetFitParams(func, htot);
  do
  {
    fitres = htot->Fit(func, opt, "", minfit, maxfit);
    // Printf("Trial: %d", trials++);
    if (trials > 20)
    {
      Printf("htot %s: %s -- FIT DOES NOT CONVERGE IN LINE %d", htot->GetName(), func->GetName(), __LINE__);
      break;
    }
  } while (fitres != 0);

  // Add the fitted function to the histogram's list before writing
  // This ensures the function is saved with the histogram in the log file.
  htot->GetListOfFunctions()->Clear();   // Remove any pre-existing functions first
  htot->GetListOfFunctions()->Add(func); // Add the current fit function 'func'
  htot->SetDirectory(0);                 // Decouple histogram from memory directory

  // Save the histogram (which now contains the fit function) to the log file
  TFile *filewithfits = TFile::Open(logfilename.Data(), "UPDATE");
  if (filewithfits && !filewithfits->IsZombie())
  {
    // Write htot using its current name, overwrite if it exists in the file
    htot->Write(htot->GetName(), TObject::kOverwrite);
    filewithfits->Close();
    delete filewithfits; // Close and delete the file pointer
  }
  else
  {
    Error("YieldMean", "Could not open log file %s to write fitted histogram.", logfilename.Data());
    if (filewithfits)
      delete filewithfits; // Clean up if zombie
  }

  cout << " Fit sys+stat for " << func->GetName() << endl;
  cout << "NDF=" << func->GetNDF() << " Chi^2=" << func->GetChisquare() << " Chi^2/NDF=" << func->GetChisquare() / func->GetNDF() << endl;

  hlo = YieldMean_LowExtrapolationHisto(htot, func, min, loprecision);
  hhi = YieldMean_HighExtrapolationHisto(htot, func, max, hiprecision);
  YieldMean_IntegralMean(htot, hlo, hhi, integral, mean, extra, kTRUE);
  hout->SetBinContent(kYield, integral);
  hout->SetBinContent(kMean, mean);
  hout->SetBinContent(kExtra, extra);
  hout->SetBinContent(kChi, func->GetChisquare() / func->GetNDF());

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
  func = nullptr;
  func = new TF1();
  f->Copy(*func);
  func = SetFitParams(func, hstat);
  do
  {
    fitres = hstat->Fit(func, opt, "", minfit, maxfit);
    // Printf("Trial: %d", trials++);
    if (trials > 20)
    {
      Printf("hstat %s: %s -- FIT DOES NOT CONVERGE IN LINE %d", hstat->GetName(), func->GetName(), __LINE__);
      break;
    }
  } while (fitres != 0);
  hlo = YieldMean_LowExtrapolationHisto(hstat, func, min, loprecision);
  hhi = YieldMean_HighExtrapolationHisto(hstat, func, max, hiprecision);

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
  func = nullptr;
  func = new TF1();
  f->Copy(*func);
  func = SetFitParams(func, hhigh);
  do
  {
    fitres = hhigh->Fit(func, opt, "", minfit, maxfit);
    // Printf("Trial: %d", trials++);
    if (trials > 20)
    {
      Printf("hhigh %s: %s -- FIT DOES NOT CONVERGE IN LINE %d", hhigh->GetName(), func->GetName(), __LINE__);
      break;
    }
  } while (fitres != 0);
  hlo = YieldMean_LowExtrapolationHisto(hhigh, func, min, loprecision);
  hhi = YieldMean_HighExtrapolationHisto(hhigh, func, max, hiprecision);
  YieldMean_IntegralMean(hhigh, hlo, hhi, integral, mean, extra);
  integral = TMath::Abs(integral - hout->GetBinContent(kYield));
  hout->SetBinContent(kYieldSysHi, integral);

  cCanvasSys->cd(1);
  func->SetLineColor(2);
  func->DrawCopy("same");

  /*
   * systematic error hard
   */

  TH1 *hhard = YieldMean_ReturnExtremeHardHisto(hsys);
  trials = 0;
  func = nullptr;
  func = new TF1();
  f->Copy(*func);
  func = SetFitParams(func, hhard);
  do
  {
    fitres = hhard->Fit(func, opt, "", minfit, maxfit);
    // Printf("Trial: %d", trials++);
    if (trials > 20)
    {
      Printf("hhard %s: %s -- FIT DOES NOT CONVERGE IN LINE %d", hhard->GetName(), func->GetName(), __LINE__);
      break;
    }
  } while (fitres != 0);
  hlo = YieldMean_LowExtrapolationHisto(hhard, func, min, loprecision);
  hhi = YieldMean_HighExtrapolationHisto(hhard, func, max, hiprecision);
  YieldMean_IntegralMean(hhard, hlo, hhi, integral, mean, extra);
  mean = TMath::Abs(mean - hout->GetBinContent(kMean));
  hout->SetBinContent(kMeanSysHi, mean);

  cCanvasSys->cd(2);
  func->SetLineColor(2);
  func->DrawCopy("same");

  /*
   * systematic error low
   */

  TH1 *hlow = YieldMean_ReturnExtremeLowHisto(hsys);
  trials = 0;
  func = nullptr;
  func = new TF1();
  f->Copy(*func);
  func = SetFitParams(func, hlow);
  do
  {
    fitres = hlow->Fit(func, opt, "", minfit, maxfit);
    // Printf("Trial: %d", trials++);
    if (trials > 20)
    {
      Printf("hlow %s: %s -- FIT DOES NOT CONVERGE IN LINE %d", hlow->GetName(), func->GetName(), __LINE__);
      break;
    }
  } while (fitres != 0);
  hlo = YieldMean_LowExtrapolationHisto(hlow, func, min, loprecision);
  hhi = YieldMean_HighExtrapolationHisto(hlow, func, max, hiprecision);
  YieldMean_IntegralMean(hlow, hlo, hhi, integral, mean, extra);
  integral = TMath::Abs(integral - hout->GetBinContent(kYield));
  hout->SetBinContent(kYieldSysLo, integral);

  cCanvasSys->cd(1);
  func->SetLineColor(4);
  func->DrawCopy("same");

  /*
   * systematic error soft
   */

  TH1 *hsoft = YieldMean_ReturnExtremeSoftHisto(hsys);
  trials = 0;
  func = nullptr;
  func = new TF1();
  f->Copy(*func);
  func = SetFitParams(func, hsoft);
  do
  {
    fitres = hsoft->Fit(func, opt, "", minfit, maxfit);
    // Printf("Trial: %d", trials++);
    if (trials > 20)
    {
      Printf("hsoft %s: %s -- FIT DOES NOT CONVERGE IN LINE %d", hsoft->GetName(), func->GetName(), __LINE__);
      break;
    }
  } while (fitres != 0);
  hlo = YieldMean_LowExtrapolationHisto(hsoft, func, min, loprecision);
  hhi = YieldMean_HighExtrapolationHisto(hsoft, func, max, hiprecision);
  YieldMean_IntegralMean(hsoft, hlo, hhi, integral, mean, extra);
  mean = TMath::Abs(mean - hout->GetBinContent(kMean));
  hout->SetBinContent(kMeanSysLo, mean);

  cCanvasSys->cd(2);
  func->SetLineColor(4);
  func->DrawCopy("same");

  // hout->Write();
  // filewithfits->Close();
  // delete filewithfits;
  return hout;
}

TF1 *SetFitParams(TF1 *func, TH1 *hist, Bool_t usePrevious, Double_t *previousFitParams)
{

  if (usePrevious)
  {
    if (previousFitParams)
      func->SetParameters(previousFitParams);

    else
      Error("SetFitParameters", "previousFitParams not passed for %s - %s.", func->GetName(), hist->GetName());
  }
  else
  {
    TString opt = "QN";
    Double_t fitMin = hist->GetXaxis()->GetBinCenter(hist->FindFirstBinAbove(0.));
    Double_t fitMax = hist->GetXaxis()->GetBinCenter(hist->FindLastBinAbove(0.));

    for (int iFit = 0; iFit < 3; iFit++)
      hist->Fit(func, opt.Data(), "", fitMin, fitMax);
  }
  func->Update();

  // Info("SetFitParameters", "%s: Fit params set.", hist->GetName());

  return func;
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
    cont = f->Integral(hlo->GetBinLowEdge(ibin + 1), hlo->GetBinLowEdge(ibin + 2), 1.e-6);
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
    cont = f->Integral(hhi->GetBinLowEdge(ibin + 1), hhi->GetBinLowEdge(ibin + 2), 1.e-6);
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