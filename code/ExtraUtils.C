/**
 * @file ExtraUtils.C
 * @author Ishaan Ahuja (ishaanahuja0@gmail.com)
 * @brief Utils and functions for extrapolation and integrated yield calculations
 * @version 1
 * @date 03-06-2025
 *
 * @copyright Copyright (c) 2025
 *
 */

#include "CascadeUtils.h"
#include "TF1.h"
#include "TMath.h"

/* definition of the fields in the histogram returned from YieldMean*/
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

enum fitFunction
{
  kLevyTsallis = 0,
  kBoltzmann,
  kBlastWave,
  // kmTScaling,
  // kBoseEinstein, //<< dont use
  // kFermiDirac,    //<< dont use
  kNumFitFunctions
};

/*****************************************************************/
// BOLTZMANN
/*****************************************************************/

Double_t
Boltzmann_Func(const Double_t *x, const Double_t *p)
{
  /* dN/dpt */

  Double_t pt = x[0];
  Double_t mass = p[0];
  Double_t mt = TMath::Sqrt(pt * pt + mass * mass);
  Double_t T = p[1];
  Double_t norm = p[2];

  return pt * norm * mt * TMath::Exp(-mt / T);
}

inline TF1 *
Boltzmann(const Char_t *name, Double_t mass, Double_t T = 0.1, Double_t norm = 1.)
{

  TF1 *fBoltzmann = new TF1(name, Boltzmann_Func, 0., 10., 3);
  fBoltzmann->SetParameters(mass, T, norm);
  fBoltzmann->SetParNames("mass", "T", "norm");
  fBoltzmann->FixParameter(0, mass);
  fBoltzmann->SetLineColor(customCascPalette[kBoltzmann]);
  return fBoltzmann;
}

/*****************************************************************/
/* LEVY-TSALLIS */
/*****************************************************************/

Double_t
LevyTsallis_Func(const Double_t *x, const Double_t *p)
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

inline TF1 *
LevyTsallis(const Char_t *name, Double_t mass, Double_t n = 5., Double_t C = 0.1, Double_t norm = 1.)
{

  TF1 *fLevyTsallis = new TF1(name, LevyTsallis_Func, 0., 10., 4);
  fLevyTsallis->SetParameters(mass, n, C, norm);
  fLevyTsallis->SetParNames("mass", "n", "C", "norm");
  fLevyTsallis->FixParameter(0, mass);
  fLevyTsallis->SetParLimits(1, 1.e-3, 1.e3); // n > 0
  fLevyTsallis->SetParLimits(2, 1.e-3, 1.e3); // C > 0
  fLevyTsallis->SetParLimits(3, 1.e-6, 1.e6); // norm > 0
  fLevyTsallis->SetLineColor(customCascPalette[kLevyTsallis]);
  return fLevyTsallis;
}
/*****************************************************************/
/* BOLTZMANN-GIBBS BLAST-WAVE */
/*****************************************************************/

static TF1 *fBGBlastWave_Integrand = NULL;
static TF1 *fBGBlastWave_Integrand_num = NULL;
static TF1 *fBGBlastWave_Integrand_den = NULL;
Double_t
BGBlastWave_Integrand(const Double_t *x, const Double_t *p)
{

  /*
     x[0] -> r (radius)
     p[0] -> mT (transverse mass)
     p[1] -> pT (transverse momentum)
     p[2] -> beta_max (surface velocity)
     p[3] -> T (freezout temperature)
     p[4] -> n (velocity profile)
  */

  Double_t r = x[0];
  Double_t mt = p[0];
  Double_t pt = p[1];
  Double_t beta_max = p[2];
  Double_t temp_1 = 1. / p[3];
  Double_t n = p[4];

  Double_t beta = beta_max * TMath::Power(r, n);
  if (beta > 0.9999999999999999)
    beta = 0.9999999999999999;
  Double_t rho = TMath::ATanH(beta);
  Double_t argI0 = pt * TMath::SinH(rho) * temp_1;
  if (argI0 > 700.)
    argI0 = 700.;
  Double_t argK1 = mt * TMath::CosH(rho) * temp_1;
  //  if (argI0 > 100 || argI0 < -100)
  //    printf("r=%f, pt=%f, beta_max=%f, temp=%f, n=%f, mt=%f, beta=%f, rho=%f, argI0=%f, argK1=%f\n", r, pt, beta_max, 1. / temp_1, n, mt, beta, rho, argI0, argK1);
  return r * mt * TMath::BesselI0(argI0) * TMath::BesselK1(argK1);
}

Double_t
BGBlastWave_Func(const Double_t *x, const Double_t *p)
{
  /* dN/dpt */

  Double_t pt = x[0];
  Double_t mass = p[0];
  Double_t mt = TMath::Sqrt(pt * pt + mass * mass);
  Double_t beta_max = p[1];
  Double_t temp = p[2];
  Double_t n = p[3];
  Double_t norm = p[4];

  if (!fBGBlastWave_Integrand)
    fBGBlastWave_Integrand = new TF1("fBGBlastWave_Integrand", BGBlastWave_Integrand, 0., 1., 5);
  fBGBlastWave_Integrand->SetParameters(mt, pt, beta_max, temp, n);
  Double_t integral = fBGBlastWave_Integrand->Integral(0., 1., (Double_t)0 /*, 1.e-6*/);
  return norm * pt * integral;
}
inline TF1 *
BGBlastWave(const Char_t *name, Double_t mass, Double_t beta_max = 0.9, Double_t temp = 0.1, Double_t n = 1., Double_t norm = 1.e6)
{

  TF1 *fBGBlastWave = new TF1(name, BGBlastWave_Func, 0., 10., 5);
  fBGBlastWave->SetParameters(mass, beta_max, temp, n, norm);
  fBGBlastWave->SetParNames("mass", "beta_max", "T", "n", "norm");
  fBGBlastWave->FixParameter(0, mass);
  fBGBlastWave->SetParLimits(1, 0.01, 0.99);
  fBGBlastWave->SetParLimits(2, 0.01, 1.);
  fBGBlastWave->SetParLimits(3, 0.01, 50.);
  fBGBlastWave->SetLineColor(customCascPalette[kBlastWave]);
  return fBGBlastWave;
}

/*****************************************************************/
// mT SCALING (EXPONENTIAL)
/*****************************************************************/
Double_t mTScaling_Func(Double_t *x, Double_t *par)
{
  // par[0] = norm
  // par[1] = T (slope parameter)
  // par[2] = mass (fixed)
  double pT = x[0];
  double mass = par[2];
  double mT = TMath::Sqrt(pT * pT + mass * mass);
  double norm = par[0];
  double T = par[1];
  // Function form: N * pT * exp(-mT / T) - Note: The pT factor is conventional for dN/dpT
  // Check if the form should be just N*exp(-mT/T) based on context,
  // but typically dN/dpT includes a pT factor from phase space (pT*dpT*dphi*dy)
  // The GetMTExpdNdptTimesPt helper function in RunYieldMean.C had norm*x*exp(-sqrt(x^2+m^2)/T)
  // which corresponds to pT * dN/(pT dpT dy) = dN/(dpT dy). Let's use that form.
  if (T < 1e-9)
    return 0; // Avoid division by zero
  return norm * pT * TMath::Exp(-mT / T);
}

inline TF1 *mTScaling(const char *name, Double_t mass, Double_t T = 0.2, Double_t norm = 1.0)
{
  TF1 *fmTScaling = new TF1(name, mTScaling_Func, 0., 10., 3);
  fmTScaling->SetParameters(norm, T, mass);
  fmTScaling->SetParNames("norm", "T", "mass");
  fmTScaling->FixParameter(2, mass);       // Fix mass
  fmTScaling->SetParLimits(0, 1e-9, 1e9);  // Norm limits
  fmTScaling->SetParLimits(1, 0.01, 10.0); // T limits
  // fmTScaling->SetLineColor(customCascPalette[kmTScaling]); // Use customCascPalette from CascadeUtils.h
  fmTScaling->SetLineWidth(2); // Consistent line width
  return fmTScaling;
}

/*****************************************************************/
// BOSE-EINSTEIN
/*****************************************************************/
Double_t BoseEinstein_Func(Double_t *x, Double_t *par)
{
  // par[0] = norm
  // par[1] = T (temperature)
  // par[2] = mass (fixed)
  // par[3] = mu (chemical potential, often assumed 0 for photons/pions, maybe for strange baryons too?)
  // Let's assume mu=0 based on typical use unless specified otherwise.
  double pT = x[0];
  double mass = par[2];
  double mT = TMath::Sqrt(pT * pT + mass * mass);
  double norm = par[0];
  double T = par[1];
  double mu = 0.0; // Assuming chemical potential is zero

  if (T < 1e-9)
    return 0; // Avoid division by zero
  double exp_arg = (mT - mu) / T;
  // Avoid overflow for large arguments
  if (exp_arg > 700)
    exp_arg = 700;

  double denominator = TMath::Exp(exp_arg) - 1.0;
  if (TMath::Abs(denominator) < 1e-9)
    return 0; // Avoid division by zero

  // Form: pT * norm / (exp((mT-mu)/T) - 1)
  // The pT factor is for dN/dpT
  return norm * pT / denominator;
}

inline TF1 *BoseEinstein(const char *name, Double_t mass, Double_t T = 0.15, Double_t norm = 1.0)
{
  TF1 *fBoseEinstein = new TF1(name, BoseEinstein_Func, 0., 10., 3); // 3 parameters (norm, T, mass) if mu=0 fixed
  fBoseEinstein->SetParameters(norm, T, mass);
  fBoseEinstein->SetParNames("norm", "T", "mass");
  fBoseEinstein->FixParameter(2, mass);       // Fix mass
  fBoseEinstein->SetParLimits(0, 1e-9, 1e9);  // Norm limits
  fBoseEinstein->SetParLimits(1, 0.01, 10.0); // T limits
  // fBoseEinstein->SetLineColor(customCascPalette[kBoseEinstein]);
  fBoseEinstein->SetLineWidth(2);
  // Warning: Bose-Einstein might diverge at low pT depending on parameters
  // May need careful fit range setting. The thesis suggests not using it? Let's keep it for now.
  return fBoseEinstein;
}

/*****************************************************************/
// FERMI-DIRAC
/*****************************************************************/
Double_t FermiDirac_Func(Double_t *x, Double_t *par)
{
  // par[0] = norm
  // par[1] = T (temperature)
  // par[2] = mass (fixed)
  // par[3] = mu (chemical potential, assumed 0)
  double pT = x[0];
  double mass = par[2];
  double mT = TMath::Sqrt(pT * pT + mass * mass);
  double norm = par[0];
  double T = par[1];
  double mu = 0.0; // Assuming chemical potential is zero

  if (T < 1e-9)
    return 0; // Avoid division by zero
  double exp_arg = (mT - mu) / T;
  // Avoid overflow for large arguments
  if (exp_arg > 700)
    exp_arg = 700;

  double denominator = TMath::Exp(exp_arg) + 1.0;
  if (TMath::Abs(denominator) < 1e-9)
    return 0; // Avoid division by zero? (Should always be > 1)

  // Form: pT * norm / (exp((mT-mu)/T) + 1)
  return norm * pT / denominator;
}

inline TF1 *FermiDirac(const char *name, Double_t mass, Double_t T = 0.15, Double_t norm = 1.0)
{
  TF1 *fFermiDirac = new TF1(name, FermiDirac_Func, 0., 10., 3); // 3 parameters if mu=0 fixed
  fFermiDirac->SetParameters(norm, T, mass);
  fFermiDirac->SetParNames("norm", "T", "mass");
  fFermiDirac->FixParameter(2, mass);       // Fix mass
  fFermiDirac->SetParLimits(0, 1e-9, 1e9);  // Norm limits
  fFermiDirac->SetParLimits(1, 0.01, 10.0); // T limits
  // fFermiDirac->SetLineColor(customCascPalette[kFermiDirac]);
  fFermiDirac->SetLineWidth(2);
  // Warning: The thesis suggests not using it? Let's keep it for now.
  return fFermiDirac;
}
/*****************************************************************/
