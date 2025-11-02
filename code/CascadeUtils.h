#ifndef CASCADEUTILS_H
#define CASCADEUTILS_H

#include <fstream>
#include <TFile.h>
#include <TSystem.h>
#include <TMath.h>
#include <TError.h>
#include <TString.h>

#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TF1.h>

#include <TMatrixD.h>
#include <TFitResult.h>

#include <TCanvas.h>
#include <THStack.h>
#include "TColor.h"
#include "TStyle.h"
#include <TLatex.h>
#include <TLegend.h>

enum particles
{
    kXi,
    kOm,
    kNumPart
};

enum signedParticles
{
    kXip,
    kXim,
    kOmp,
    kOmm,
    kXiC,
    kOmC,
    kNumSignedPart
};

enum fitFunctions
{
    kGausPol2 = 1,
    kDoubleGausPol2 = 2,
    kDoubleGausPol3 = 3
    // bgTSpectrum = 4
};

const Double_t fMass_Xi = 1.32171;
const Double_t fMass_Om = 1.67245;

const double fMultbins_Xi[11] = {0, 5, 10, 15, 20, 30, 40, 50, 60, 80, 100}; // V0A
const double fMultbins_Om[6] = {0, 5, 15, 30, 60, 100};                      // V0A
const Int_t fNmultbins_Xi = sizeof(fMultbins_Xi) / sizeof(double) - 1;
const Int_t fNmultbins_Om = sizeof(fMultbins_Om) / sizeof(double) - 1;

const double fPtbins_Xi[] = {0.8, 1.1, 1.3, 1.5, 1.7, 1.9, 2.1, 2.3, 2.5, 2.7, 2.9, 3.1, 3.5, 4, 5.5};
const double fPtbins_Om[] = {0.9, 1.6, 2., 2.4, 2.9, 3.5, 5};
const Int_t fNptbins_Xi = sizeof(fPtbins_Xi) / sizeof(double) - 1;
const Int_t fNptbins_Om = sizeof(fPtbins_Om) / sizeof(double) - 1;

double dNchXi[10] = {53.22, 42.40, 37.30, 33.64, 29.30, 24.49, 20.34, 16.46, 10.97, 4.47};
double dNchXiErr[10] = {1.38, 1.10, 0.97, 0.87, 0.76, 0.66, 0.53, 0.43, 0.29, 0.12};
double dNchOm[5] = {53.22, 39.86, 30.77, 20.42, 7.63};
double dNchOmErr[5] = {1.38, 1.10, 0.86, 0.56, 0.20};

// Define Table 3 data from "Multi-strange baryon production in p–Pb collisions at 5TeV.pdf"
// For Xi (particleNamePlot[0] == "xiC")
const int nPoints_ppb5_Xi = 7;
double dNch_ppb5_Xi[] = {45.0, 36.2, 30.5, 23.2, 16.1, 9.8, 4.3};
double dNchErr_ppb5_Xi[] = {1.0, 0.8, 0.7, 0.5, 0.4, 0.2, 0.1};                    // Used for x-error of points and boxes
double yield_ppb5_Xi[] = {0.2354, 0.1861, 0.1500, 0.1100, 0.0726, 0.0398, 0.0143}; // Assuming 0.1100 for "0.1~100"
double statErr_ppb5_Xi[] = {0.0020, 0.0016, 0.0010, 0.0006, 0.0006, 0.0004, 0.0003};
double sysErr_ppb5_Xi[] = {0.0161, 0.0138, 0.0112, 0.0085, 0.0065, 0.0031, 0.0015}; // Symmetric systematic error

// For Omega (particleNamePlot[1] == "omC")
const int nPoints_ppb5_Om = 7;
double dNch_ppb5_Om[] = {45.0, 36.2, 30.5, 23.2, 16.1, 9.8, 4.3};
double dNchErr_ppb5_Om[] = {1.0, 0.8, 0.7, 0.5, 0.4, 0.2, 0.1}; // Used for x-error of points and boxes
double yield_ppb5_Om[] = {0.0260, 0.0215, 0.0167, 0.0120, 0.0072, 0.0042, 0.0013};
double statErr_ppb5_Om[] = {0.0011, 0.0008, 0.0006, 0.0005, 0.0003, 0.0002, 0.0003};
double sysErr_ppb5_Om[] = {0.0034, 0.0029, 0.0022, 0.0016, 0.0010, 0.0006, 0.0003}; // Symmetric systematic error

/// 2020: Multiplicity dependence of (multi-)strange hadron production in proton-proton collisions at snn = 13 TeV
// ───(A) Number of multiplicity classes / points ──────────────────────────
const Int_t Npoints_Xi_2020 = 10;

// ───(B) Xi± integrated yields (dN/dy) and their errors ───────────────────
// Extracted from HEPData Table 8c (see file “HEPData-ins1748157-v2-Table_8c.yaml”)
// for sqrt(s)=13 TeV, |y|<0.5, in ten V0M multiplicity classes.
Double_t XiYields_2020[Npoints_Xi_2020] = {
    0.1284456594,
    0.0924905535,
    0.0750700817,
    0.0609966459,
    0.0550614327,
    0.0447087709,
    0.0333413242,
    0.0261796735,
    0.0158834917,
    0.0061615164};

Double_t XiYields_stat_2020[Npoints_Xi_2020] = {
    0.0029034064,
    0.0012247650,
    0.0010058078,
    0.0008704892,
    0.0008879130,
    0.0005695243,
    0.0005306301,
    0.0004594446,
    0.0002620426,
    0.0001610817};

Double_t XiYields_systTot_2020[Npoints_Xi_2020] = {
    0.0083107249,
    0.0059671232,
    0.0053410856,
    0.0041266370,
    0.0042970402,
    0.0036712523,
    0.0029716380,
    0.0024463140,
    0.0015786774,
    0.0010696392};

// ───(C) Mean charged-particle multiplicity ⟨dNch/dη⟩ and “error” ────────
// The “value” is the central ⟨dNch/dη⟩, and we take (high−low)/2 as its bin‐width.
Double_t dNch_Xi_2020[Npoints_Xi_2020] = {
    25.75,
    19.83,
    16.12,
    13.76,
    12.06,
    10.11,
    8.07,
    6.48,
    4.64,
    2.52};

Double_t dNch_err_Xi_2020[Npoints_Xi_2020] = {
    0.40,
    0.30,
    0.24,
    0.21,
    0.18,
    0.15,
    0.12,
    0.10,
    0.07,
    0.04};

// ───(A) Number of multiplicity classes / points ──────────────────────────
const Int_t Npoints_Omega_2020 = 5;

// ───(B) Ω± integrated yields (dN/dy) and their errors ────────────────────
// Extracted from “HEPData-ins1748157-v2-Table_8cO.yaml” for sqrt(s)=13 TeV, |y|<0.5,
// in five V0M multiplicity classes:
Double_t OmegaYields_2020[Npoints_Omega_2020] = {
    0.0103214156,
    0.0070094986,
    0.0043486953,
    0.0022854405,
    0.0006284247};

Double_t OmegaYields_2020_stat[Npoints_Omega_2020] = {
    0.0005386190,
    0.0002650247,
    0.0002384983,
    0.0001724946,
    0.0000597499};

Double_t OmegaYields_2020_systTot[Npoints_Omega_2020] = {
    0.0009359367,
    0.0006899374,
    0.0004698625,
    0.0003108356,
    0.0000987166};

// ───(C) Mean charged-particle multiplicity ⟨dNch/dη⟩ and “error” ────────
// For each V0M class: “value” = central ⟨dNch/dη⟩, and “error” = (high−low)/2.
Double_t dNch_Omega_2020[Npoints_Omega_2020] = {
    21.01,
    14.94,
    10.76,
    7.28,
    3.37};

Double_t dNch_err_Omega_2020[Npoints_Omega_2020] = {
    0.32,
    0.225,
    0.16,
    0.11,
    0.05};

/// 2017: Enhanced production of multi-strange hadrons in high-multiplicity proton-proton collisions

// ───(A) Number of multiplicity classes / points ──────────────────────────
const Int_t Npoints_Xi2017 = 10;

// ───(B) Ξ± integrated yields (dN/dy) for 2017 and their errors ───────────
// Extracted from “HEPData-ins1471838-v1-Table_44.yaml” for √s = 7 TeV, |y|<0.5,
// in ten multiplicity bins.
Double_t XiYields2017[Npoints_Xi2017] = {
    0.096160, // ⟨dN/dy⟩ in highest-multiplicity bin
    0.071690,
    0.059310,
    0.048320,
    0.040680,
    0.033330,
    0.024850,
    0.018890,
    0.012050,
    0.005199 // ⟨dN/dy⟩ in lowest-multiplicity bin
};

Double_t XiYields2017_stat[Npoints_Xi2017] = {
    0.001598, // stat error, 1st bin
    0.0007437,
    0.0006245,
    0.0005393,
    0.0004974,
    0.0003545,
    0.0002936,
    0.0002831,
    0.0001746,
    0.0001069 // stat error, last bin
};

Double_t XiYields2017_systTot[Npoints_Xi2017] = {
    0.006081, // syst,total error, 1st bin
    0.004629,
    0.003992,
    0.003352,
    0.003042,
    0.002549,
    0.002045,
    0.001698,
    0.001275,
    0.0007461 // syst,total error, last bin
};

// ───(C) Mean charged-particle multiplicity ⟨dNch/dη⟩ and “error” ────────
// From the same table: “value” = central ⟨dNch/dη⟩, and error = (high−low)/2.
Double_t dNch_Xi_2017[Npoints_Xi2017] = {
    21.29, // bin: [20.65, 21.93]
    16.51, // bin: [16.01, 17.01]
    13.46, // bin: [13.06, 13.86]
    11.51, // bin: [11.16, 11.86]
    10.08, // bin: [9.78, 10.38]
    8.45,  // bin: [8.20, 8.70]
    6.72,  // bin: [6.51, 6.93]
    5.40,  // bin: [5.23, 5.57]
    3.90,  // bin: [3.76, 4.04]
    2.26   // bin: [2.14, 2.38]
};

Double_t dNch_err_Xi_2017[Npoints_Xi2017] = {
    0.64, // (21.93−20.65)/2
    0.50, // (17.01−16.01)/2
    0.40, // (13.86−13.06)/2
    0.35, // (11.86−11.16)/2
    0.30, // (10.38−9.78)/2
    0.25, // (8.70−8.20)/2
    0.21, // (6.93−6.51)/2
    0.17, // (5.57−5.23)/2
    0.14, // (4.04−3.76)/2
    0.12  // (2.38−2.14)/2
};

// ───(A) Number of multiplicity classes / points ──────────────────────────
const Int_t Npoints_Omega2017 = 5;

// ───(B) Ω± integrated yields (dN/dy) for 2017 and their errors ───────────
// Extracted from the provided YAML for √s = 7 TeV, |y|<0.5,
// in five multiplicity bins.
Double_t OmegaYields2017[Npoints_Omega2017] = {
    0.008274, // 1st bin
    0.005082, // 2nd bin
    0.002976, // 3rd bin
    0.001863, // 4th bin
    0.0005767 // 5th bin
};

Double_t OmegaYields2017_stat[Npoints_Omega2017] = {
    0.0002640, // stat error, 1st bin
    0.0001951, // stat error, 2nd bin
    0.0000951, // stat error, 3rd bin (9.51e-05)
    0.0001179, // stat error, 4th bin
    0.00006029 // stat error, 5th bin (6.029e-05)
};

Double_t OmegaYields2017_systTot[Npoints_Omega2017] = {
    0.0007787, // syst,total, 1st bin
    0.0005235, // syst,total, 2nd bin
    0.0002974, // syst,total, 3rd bin
    0.0002373, // syst,total, 4th bin
    0.0001050  // syst,total, 5th bin
};

// ───(C) Mean charged-particle multiplicity ⟨dNch/dη⟩ and “error” ────────
// From the independent-variable section: “value” = central ⟨dNch/dη⟩,
// and “error” = (high – low) / 2 for each bin.
Double_t dNch_Om_2017[Npoints_Omega2017] = {
    17.47, // bin: [16.95, 17.99]
    12.48, // bin: [12.11, 12.85]
    8.99,  // bin: [8.72,  9.26]
    6.06,  // bin: [5.87,  6.25]
    2.89   // bin: [2.75,  3.03]
};

Double_t dNch_err_Om_2017[Npoints_Omega2017] = {
    0.52, // (17.99 – 16.95) / 2 = 1.04 / 2
    0.37, // (12.85 – 12.11) / 2 = 0.74 / 2
    0.27, // (9.26  –  8.72) / 2 = 0.54 / 2
    0.19, // (6.25  –  5.87) / 2 = 0.38 / 2
    0.14  // (3.03  –  2.75) / 2 = 0.28 / 2
};

// 2014: Multi-strange baryon production at mid-rapidity in Pb-Pb collisions at sqrt(s_NN) = 2.76 TeV
//  --- Mapping Data: Centrality to dNch/deta for Pb-Pb @ 2.76 TeV (from PRL 106, 032301) ---
const int Nmap_pbpb_2011 = 9;
double dNch_map_2011[Nmap_pbpb_2011] = {1601.0, 1294.0, 966.0, 649.0, 426.0, 261.0, 149.0, 76.0, 35.0};
double dNch_err_map_2011[Nmap_pbpb_2011] = {60.0, 49.0, 37.0, 23.0, 15.0, 9.0, 6.0, 4.0, 2.0};

// ───(A) Number of centrality bins ─────────────────────────────────────────
const Int_t Npoints_2014 = 5;

// ───(B) Ξ– + Ξ¯+ combined integrated yields (dN/dy) and their errors ─────
// Taken from the “(Ξ– + Ξ¯+)” block in the YAML (√sNN = 2.76 TeV, |y|<0.5).
Double_t XiCombinedYields2014[Npoints_2014] = {
    6.67, // 0–10%
    5.14, // 10–20%
    3.03, // 20–40%
    1.07, // 40–60%
    0.24  // 60–80%
};

Double_t XiCombinedYields2014_stat[Npoints_2014] = {
    0.08, // stat error, 0–10%
    0.06, // stat error, 10–20%
    0.03, // stat error, 20–40%
    0.01, // stat error, 40–60%
    0.006 // stat error, 60–80%
};

Double_t XiCombinedYields2014_systTot[Npoints_2014] = {
    0.47, // syst,total, 0–10%
    0.36, // syst,total, 10–20%
    0.22, // syst,total, 20–40%
    0.08, // syst,total, 40–60%
    0.019 // syst,total, 60–80%
};

// ───(C) Ω– + Ω¯+ combined integrated yields (dN/dy) and their errors ────
// Taken from the “(Ω– + Ω¯+)” block in the YAML (√sNN = 2.76 TeV, |y|<0.5).
Double_t OmegaCombinedYields2014[Npoints_2014] = {
    1.19, // 0–10%
    0.78, // 10–20%
    0.48, // 20–40%
    0.17, // 40–60%
    0.032 // 60–80%
};

Double_t OmegaCombinedYields2014_stat[Npoints_2014] = {
    0.06,  // stat error, 0–10%
    0.04,  // stat error, 10–20%
    0.02,  // stat error, 20–40%
    0.007, // stat error, 40–60%
    0.003  // stat error, 60–80%
};

Double_t OmegaCombinedYields2014_systTot[Npoints_2014] = {
    0.19,  // syst,total, 0–10%
    0.15,  // syst,total, 10–20%
    0.08,  // syst,total, 20–40%
    0.029, // syst,total, 40–60%
    0.005  // syst,total, 60–80%
};

// From the “CENTRALITY” independent variable: bins [0–10], [10–20], [20–40], [40–60], [60–80].
Double_t centBinLow2014[Npoints_2014] = {
    0.0,  // (0 + 10) / 2
    10.0, // (10 + 20) / 2
    20.0, // (20 + 40) / 2
    40.0, // (40 + 60) / 2
    60.0  // (60 + 80) / 2
};
Double_t centBinHigh2014[Npoints_2014] = {
    10.0, // (0 + 10) / 2
    20.0, // (10 + 20) / 2
    40.0, // (20 + 40) / 2
    60.0, // (40 + 60) / 2
    80.0  // (60 + 80) / 2
};

// chosen marker style palette
Int_t markerStyles[] = {kCircle, kFullSquare, kFullTriangleUp, kStar, kFullTriangleDown, kOpenSquare, kFullDiamond, kFullCross, kOpenDiamond, kFullFourTrianglesX, kFullThreeTriangles, kFullCircle};
Int_t openMarkerStyles[] = {kOpenCircle, kOpenTriangleUp, kOpenDiamond, kOpenTriangleDown, kOpenDiamond, kOpenStar, kOpenCross, kOpenSquare, kOpenFourTrianglesX, kOpenThreeTriangles, kOpenSquareDiagonal};

/**
 * @brief Sets a custom color palette for ROOT graphics
 *
 * This function defines a palette of 15 visually distinct colors and
 * sets it as the default palette in ROOT's global style. The colors
 * are carefully chosen to be distinguishable from each other, making
 * them suitable for plots and visualizations that need to differentiate
 * between multiple data series.
 *
 * Colors included: red, blue, green, orange, purple, yellow, brown,
 * pink, teal, coral, light blue, light magenta, light green, gold,
 * and gray.
 *
 * Note: The yellow color (#ffff33) should be used with caution on
 * white backgrounds due to its high saturation.
 *
 * @see TColor::GetColor
 * @see TStyle::SetPalette
 */
int customCascPalette[15] = {
    // TColor::GetColor("#377eb8"), // Blue (Azure?)
    TColor::GetColor("#4169e1"), // Royal Blue (for extensive yield comparison only)
    TColor::GetColor("#e41a1c"), // Red
    TColor::GetColor("#4daf4a"), // Green
    TColor::GetColor("#984ea3"), // Purple
    TColor::GetColor("#ff7f00"), // Orange
    // TColor::GetColor("#ffff33"), // Yellow (saturated, but use with caution on white backgrounds)
    TColor::GetColor("#999999"), // Gray
    TColor::GetColor("#a65628"), // Brown
    TColor::GetColor("#f781bf"), // Pink
    TColor::GetColor("#66c2a5"), // Teal
    TColor::GetColor("#8da0cb"), // Light Blue
    TColor::GetColor("#e78ac3"), // Light Magenta
    TColor::GetColor("#a6d854"), // Light Green
    TColor::GetColor("#fc8d62"), // Coral
    // TColor::GetColor(kGreen-3), // Light Green
    TColor::GetColor("#ffd92f"), // Gold
    // TColor::GetColor(kViolet + 10), // Gold
    TColor::GetColor("#FFDB58") // Mustard Yellow (saturated, but use with caution on white backgrounds)

};
inline void SetCustomColorPalette()
{
    // Set the palette in gStyle so that it is used by default when drawing
    gStyle->SetPalette(15, customCascPalette);
}

/**
 * @brief Evaluates a double Gaussian function with a polynomial of degree 2
 *
 * This function returns the sum of two Gaussian distributions with the same mean
 * but different widths, plus a quadratic polynomial background.
 *
 * @param x[0] The x coordinate where the function is evaluated
 * @param par[0] Amplitude of the first Gaussian
 * @param par[1] Mean of both Gaussians
 * @param par[2] Sigma (width) of the first Gaussian
 * @param par[3] Amplitude of the second Gaussian
 * @param par[4] Sigma (width) of the second Gaussian
 * @param par[5] Constant term of the polynomial (p0)
 * @param par[6] Linear term coefficient of the polynomial (p1)
 * @param par[7] Quadratic term coefficient of the polynomial (p2)
 *
 * @return Value of the combined function: (Gaussian1 + Gaussian2 + p0 + p1*x + p2*x^2)
 */
inline Double_t DoubleGausPol2(double *x, double *par) // with same mean parameter for both gaussians
{
    return par[0] * TMath::Gaus(x[0], par[1], par[2]) + par[3] * TMath::Gaus(x[0], par[1], par[4]) + par[5] + par[6] * x[0] + par[7] * x[0] * x[0];
}

/**
 * Function that combines two Gaussian distributions with a polynomial background of order 3.
 * Both Gaussian distributions share the same mean parameter.
 *
 * @param x[0] Input value at which to evaluate the function
 * @param par[0] Amplitude of first Gaussian
 * @param par[1] Mean value (shared between both Gaussians)
 * @param par[2] Sigma (width) of first Gaussian
 * @param par[3] Amplitude of second Gaussian
 * @param par[4] Sigma (width) of second Gaussian
 * @param par[5] Constant term of polynomial (p0)
 * @param par[6] Linear term coefficient of polynomial (p1)
 * @param par[7] Quadratic term coefficient of polynomial (p2)
 * @param par[8] Cubic term coefficient of polynomial (p3)
 *
 * @return Value of the combined function: (Gaussian1 + Gaussian2 + p0 + p1*x + p2*x^2 + p3*x^3)
 */
inline Double_t DoubleGausPol3(double *x, double *par) // with same mean parameter for both gaussians
{
    return par[0] * TMath::Gaus(x[0], par[1], par[2]) + par[3] * TMath::Gaus(x[0], par[1], par[4]) + par[5] + par[6] * x[0] + par[7] * x[0] * x[0] + par[8] * x[0] * x[0] * x[0];
}

/**
 * @brief Second degree polynomial function with exclusion region capability
 *
 * This function implements a second degree polynomial (ax^2 + bx + c) with the ability
 * to exclude points within a specified region around a peak. Used primarily for background
 * fitting in mass spectra analysis.
 *
 * @param x Pointer to the x coordinate [input]
 * @param par Array of parameters:
 *        par[0] = quadratic term coefficient (a)
 *        par[1] = linear term coefficient (b)
 *        par[2] = constant term (c)
 *        par[3] = peak mass position
 *        par[4] = peak width (sigma)
 *        par[5] = rejection flag (boolean)
 *
 * @return Double_t The function value at x, or 0 if point is rejected
 *
 * @note When par[5] is true, points within ±sigma around peak mass are rejected using TF1::RejectPoint()
 */
inline Double_t Pol2Exclude(double *x, double *par)
{
    // par[5] = rejectPoint(boolean); par[3] = peakFitMass; par[4]=peakFitSigma;
    if (par[5] && x[0] > (par[3] - par[4]) && x[0] < (par[3] + par[4]))
    {
        TF1::RejectPoint();
        return 0;
    }
    return par[0] * x[0] * x[0] + par[1] * x[0] + par[2];
}

/**
 * @brief A third-order polynomial function with optional point exclusion around a peak.
 *
 * This function evaluates a third-order polynomial and can optionally exclude points
 * within a specified range around a peak. Used typically for background fitting in
 * particle physics analysis.
 *
 * @param x    Pointer to the x value where the function is evaluated
 * @param par  Array of parameters:
 *             - par[0]: Coefficient of x³ term
 *             - par[1]: Coefficient of x² term
 *             - par[2]: Coefficient of x term
 *             - par[3]: Peak position (mass)
 *             - par[4]: Peak width (sigma)
 *             - par[5]: Boolean flag for point rejection
 *             - par[6]: Constant term
 *
 * @return The value of the polynomial at x[0], or 0 if the point is rejected
 */
inline Double_t Pol3Exclude(double *x, double *par)
{
    // par[5] = rejectPoint(boolean); par[3] = peakFitMass; par[4]=peakFitSigma;
    if (par[5] && x[0] > (par[3] - par[4]) && x[0] < (par[3] + par[4]))
    {
        TF1::RejectPoint();
        return 0;
    }
    return par[0] * x[0] * x[0] * x[0] + par[1] * x[0] * x[0] + par[2] * x[0] + par[6];
}

/**
 * @brief Evaluates a combination of a Gaussian function and a quadratic polynomial
 *
 * @param x Array containing the x coordinate [x[0]]
 * @param par Array of parameters:
 *        par[0]: Amplitude of Gaussian
 *        par[1]: Mean of Gaussian
 *        par[2]: Sigma of Gaussian
 *        par[3]: Coefficient of x²
 *        par[4]: Coefficient of x
 *        par[5]: Constant term
 *
 * @return Double_t Result of f(x) = A*Gaus(x,μ,σ) + ax² + bx + c
 *                  where Gaus is the normalized Gaussian function
 */
inline Double_t GausPol2(double *x, double *par)
{
    return par[0] * TMath::Gaus(x[0], par[1], par[2]) + par[3] * x[0] * x[0] + par[4] * x[0] + par[5];
}

/**
 * @brief Reads a list of file paths from a text file and returns them as a vector of strings
 *
 * This function opens the specified input file and reads it line by line,
 * storing each non-empty line as a string in a vector. Each line is expected
 * to contain a file path.
 *
 * @param inputFilePath The path to the text file containing the list of files
 * @return std::vector<std::string> A vector containing all non-empty lines from the input file.
 *         Returns an empty vector if the file cannot be opened or read.
 *
 * @note The function automatically skips empty lines in the input file
 * @note The file is automatically closed after reading
 */
std::vector<std::string> GetFileList(const std::string &inputFilePath)
{
    std::vector<std::string> fileList;
    std::ifstream inputFile(inputFilePath);
    std::string line;

    // Check if the file is open and can be read
    if (!inputFile.is_open())
    {
        Error("Utils: GetFileList", "Cannot open file '%s' !", inputFilePath.c_str());
        return fileList; // Return an empty vector in case of error
    }

    // Read the file line by line
    while (std::getline(inputFile, line))
    {
        // Avoid adding empty lines (if any)
        if (!line.empty())
        {
            fileList.push_back(line);
        }
    }

    // Close the file
    inputFile.close();

    return fileList;
}

/**
 * @brief Sets up and returns the output folder path
 *
 * If no folder name is provided, defaults to "Output" directory in current working directory.
 * Expands any environment variables or special characters in the provided path.
 * Creates the folder if it does not exist and sets appropriate permissions (755).
 *
 * @param folderName Target folder path (optional). If empty, defaults to "./Output"
 * @return TString The full path to the output folder
 *
 * @warning Prints warning if no folder name is provided
 * @note Creates intermediate directories as needed with permissions set to 755
 */
inline TString SetOutputFolder(TString folderName = "")
{
    if (folderName.IsNull())
    {
        folderName = gSystem->GetWorkingDirectory();
        folderName.Append("/Output");
        Warning("Utils: SetOutputFolder", "Folder Name is empty! Defaulting to '%s' ...", folderName.Data());
    }
    else
        gSystem->ExpandPathName(folderName);

    if (gSystem->AccessPathName(folderName.Data())) /// returns true if folder path does NOT exist
    {
        Info("Utils: SetOutputFolder", "Creating folder '%s' ...", folderName.Data());
        gSystem->mkdir(folderName.Data(), kTRUE); // makes the path if it doesn't exist
        gSystem->Chmod(folderName.Data(), 0755);
    }
    return folderName;
}

/**
 * @brief Opens a ROOT file with specified options and performs validity checks
 *
 * @param fileName Name of the file to open. If no extension is provided, .root will be appended
 * @param options File access mode ("READ", "NEW", "CREATE", "RECREATE"). Defaults to "READ"
 *
 * @return TFile* Pointer to the opened file, nullptr if operation fails
 *
 * @details This function:
 * - Checks if the provided filename is not empty
 * - Appends .root extension if filename has no extension
 * - Creates directory structure if file is being created/recreated
 * - Opens the file and performs zombie check
 * - Prints informative messages about the operation
 *
 * @note The caller is responsible for closing and deleting the returned TFile pointer
 */
inline TFile *OpenFile(TString fileName, TString options = "READ")
{
    if (fileName.IsNull())
    {
        SysError("Utils: OpenFile", "Error: Provided file name is empty for %sING!", options.Data());
        return nullptr;
    }

    // append .root at the end if filename is without extension
    if (!fileName.EndsWith(".root"))
        fileName.Append(".root");

    // create file path if it doesn't exist
    if (options == "NEW" || options == "CREATE" || options == "RECREATE")
    {
        TString folderPath = fileName(0, fileName.Last('/'));
        if (folderPath != "")
        {
            SetOutputFolder(folderPath);
        }
    }

    Info("Utils: OpenFile", "Opening file '%s' for %sING ...", fileName.Data(), options.Data());

    TFile *fileHandle = TFile::Open(fileName.Data(), options.Data());

    // validity check
    if ((!fileHandle) || fileHandle->IsZombie())
    {
        delete fileHandle;
        SysError("Utils: OpenFile", "Error: Cannot open file '%s' for %sING!", fileName.Data(), options.Data());
        return nullptr;
    }
    else
        return fileHandle;
}

/**
 * @brief Paints a THStack object on a TCanvas with specified settings
 *
 * @param c Reference to the TCanvas where the stack will be painted
 * @param hs Reference to the THStack object to be painted
 * @param setLogY Boolean flag to set logarithmic Y axis (default: true)
 * @param yAxisTitle Title for Y axis (default: "#frac{1}{#it{N}_{inel}} #frac{d#it{N}}{d#it{p}_{T}}")
 * @param xAxisTitle Title for X axis (default: "#it{p}_{T} (GeV/c)")
 *
 * This function:
 * - Clears the canvas
 * - Draws the stack without stacking ("nostack" option)
 * - Sets logarithmic Y axis if specified
 * - Builds a legend
 * - Sets X and Y axis titles
 * - Updates the canvas
 * 
 *  --- Comment/Uncomment/Change Values for different plotting styles ---
 */
inline void PaintStack(TCanvas &c, THStack &hs, Bool_t setLogY = kTRUE, TString yAxisTitle = "#frac{1}{#it{N}_{inel}} #frac{d#it{N}_{raw}}{d#it{p}_{T}} (GeV/#it{c})^{-1}", TString xAxisTitle = "#it{p}_{T} (GeV/#it{c})")
{
    c.Clear();
    c.cd();
    // c.SetLeftMargin(0.12);  // New mod for final thesis plots
    c.SetRightMargin(0.05); // New mod for final thesis plots
    c.SetTopMargin(0.05);   // New mod for final thesis plots
    // c.SetGridy();

    // hs.SetMaximum(hs.GetMinimum()-0.03); // New mod for final thesis plots

    hs.Draw("plc pmc nostack");
    // gPad->SetGridy(1); // Enable grid on y-axis

    if (setLogY)
        gPad->SetLogy();
    // gPad->BuildLegend(0.85, 0.5, 0.995, 0.9, ""); <<bigger legend size
    // gPad->BuildLegend(0.85, 0.6, 0.995, 0.9, "Multiplicity Classes");


    // TLegend *legend = new TLegend(0.85, 0.7, 0.995, 0.9, "", "brNDC");
    // legend->SetHeader("Multiplicity Classes", "C"); // Set header for the legend
    // for (int i = 0; i < hs.GetHists()->GetEntries()-1; ++i)
    // {
    //     legend->AddEntry(hs.GetHists()->At(i), hs.GetHists()->At(i)->GetTitle(), "lp"); // Add remaining histograms to legend
    // }

    hs.GetXaxis()->SetTitle(xAxisTitle.Data());
    hs.GetYaxis()->SetTitle(yAxisTitle.Data());

    /// New mod for final thesis plots
    TLegend *legend = c.BuildLegend(0.48, 0.74, 0.92, 0.90);
    legend->SetNColumns(2);   // Set number of columns in the legend
    legend->SetBorderSize(0); // Remove border around the legend
    // legend->SetHeader("Sources", "C"); // Set header for the legend
    legend->SetFillColor(0); // Set legend background color to transparent
    legend->SetFillStyle(0); // Set fill style to transparent
    legend->SetTextSize(0.025);
    TLatex *latexCollTitle = new TLatex(0.15, 0.84, "p-Pb #sqrt{s_{NN}} = 8.16 TeV");
    TLatex *latexPlotTitle = new TLatex(0.15, 0.80, hs.GetTitle());
    latexCollTitle->SetTextSize(0.025);
    latexPlotTitle->SetTextSize(0.025);
    latexCollTitle->SetNDC(kTRUE);
    latexPlotTitle->SetNDC(kTRUE);
    latexCollTitle->Draw();
    latexPlotTitle->Draw();
    hs.SetTitle("");
    legend->Draw();

    ///

    c.Modified();
    c.ForceUpdate();
}

/**
 * @brief Saves a ROOT canvas as an image file
 *
 * This function saves either a specified TCanvas or the current pad (gPad) as an image file
 * in the specified output location. After saving, it sets the file permissions to 0755.
 *
 * @param outputFolder Base output directory path
 * @param imageFolder Subdirectory within output folder where image will be saved
 * @param imageName Name of the output image file (without extension)
 * @param imageFormat Format/extension of the output image (e.g., "png", "pdf")
 * @param c Pointer to TCanvas to be saved. If nullptr, current pad (gPad) will be used
 *
 * @note The function uses SetOutputFolder to ensure the directory exists before saving
 * @note File permissions are set to 0755 (rwxr-xr-x) after saving
 */
inline void SaveImage(TString outputFolder, TString imageFolder, TString imageName, TString imageFormat, TCanvas *c = nullptr)
{
    TString imagePath = SetOutputFolder(outputFolder + "/" + imageFolder);
    if (c)
        c->Print(TString::Format("%s/%s.%s", imagePath.Data(), imageName.Data(), imageFormat.Data()), imageFormat.Data());
    else
        gPad->Print(TString::Format("%s/%s.%s", imagePath.Data(), imageName.Data(), imageFormat.Data()), imageFormat.Data());
    gSystem->Chmod(TString::Format("%s/%s.%s", imagePath.Data(), imageName.Data(), imageFormat.Data()), 0755);
}

inline Double_t GetChi2(TFitResultPtr fFitResult)
{
    if (!fFitResult)
        return 0;
    Double_t chi2 = fFitResult->Chi2();
    if (TMath::IsNaN(chi2))
        return 0;
    else
        return chi2;
}

inline Double_t GetNdf(TFitResultPtr fFitResult)
{
    if (!fFitResult)
        return 0;
    Double_t ndf = fFitResult->Ndf();
    if (TMath::IsNaN(ndf))
        return 0;
    else
        return ndf;
}

inline Double_t GetReducedChi2(TFitResultPtr fFitResult)
{
    if (!fFitResult)
        return 0;
    Double_t reducedChi2 = fFitResult->Chi2() / fFitResult->Ndf();
    if (TMath::IsNaN(reducedChi2))
        return 0;
    else
        return reducedChi2;
}

inline Double_t GetProb(TFitResultPtr fFitResult)
{
    if (!fFitResult)
        return 0;
    Double_t prob = fFitResult->Prob();
    if (TMath::IsNaN(prob))
        return 0;
    else
        return prob;
}

/**
 * @brief Calculates the Root Mean Square (RMS) value of an array of doubles
 *
 * This function computes the RMS by:
 * 1. Squaring each element
 * 2. Finding the mean of squared values
 * 3. Taking the square root of the mean
 *
 * @param arraySize The number of elements in the input array
 * @param array Pointer to the array of double values
 * @return double The RMS value of the array elements. Returns 0.0 if array is nullptr or arraySize <= 0
 */
double CalculateRMS(int arraySize, double *array)
{
    if (!array || arraySize <= 0)
        return 0.0;

    double square = 0.0, mean = 0.0, root = 0.0;

    // Calculate square
    for (int iArr = 0; iArr < arraySize; iArr++)
    {
        square += array[iArr] * array[iArr];
    }

    // Calculate Mean
    mean = square / arraySize;

    // Calculate Root
    root = sqrt(mean);

    return root;
}

/**
 * Calculates the integral and its error for a specified region in a 2D histogram
 *
 * @param h2 Pointer to the 2D histogram
 * @param binsMC Array containing bin range [xmin, xmax, ymin, ymax]
 * @param genErr Reference to store the error of the integral
 * @return Double_t The integral value for the specified region
 */
inline Double_t GetGeneratedParticles(TH2 *h2, Int_t binsMC[], Double_t &genErr)
{
    Double_t gen = h2->IntegralAndError(binsMC[0], binsMC[1], binsMC[2], binsMC[3], genErr);
    return gen;
}

/**
 * @brief Calculates the error propagation in a ratio A/B
 *
 * This function computes the uncertainty in a ratio of two quantities using error propagation formula.
 * If the denominator (B) is zero, returns -1 to indicate invalid ratio.
 * The error is calculated using the formula: sqrt(|σ²ᴀ/B² - (A²/B⁴)σ²ʙ|)
 * where σᴀ and σʙ are the uncertainties in A and B respectively.
 *
 * @param A Numerator value
 * @param Aerr Error/uncertainty in numerator (σᴀ)
 * @param B Denominator value
 * @param Berr Error/uncertainty in denominator (σʙ)
 *
 * @return Double_t Error in the ratio A/B, or -1 if B is zero
 */
inline Double_t ErrorInRatio(Double_t A, Double_t Aerr, Double_t B, Double_t Berr)
{
    Double_t err = 0.;
    if (B == 0)
        err = -1.;
    else
    {
        Double_t errorfromtop = Aerr * Aerr / (B * B);
        Double_t errorfrombottom = ((A * A) / (B * B * B * B)) * Berr * Berr;
        // err = TMath::Sqrt(errorfromtop + errorfrombottom);
        err = TMath::Sqrt(TMath::Abs(errorfromtop - errorfrombottom));
    }
    return err;
}

/// @brief Override ROOT's TF1 print function - for clarity
/// @details In case of printing function parameters, skips the list of saved points when using verbose ("V") option. See commented part below.
/// @param option "V": prints TF1 GetParameters()
void TF1::Print(Option_t *option) const
{
    if (fType == EFType::kFormula)
    {
        printf("Formula based function:     %s \n", GetName());
        assert(fFormula);
        fFormula->Print(option);
    }
    else if (fType > 0)
    {
        if (fType == EFType::kInterpreted)
            printf("Interpreted based function: %s(double *x, double *p).  Ndim = %d, Npar = %d  \n", GetName(), GetNdim(),
                   GetNpar());
        else if (fType == EFType::kCompositionFcn)
        {
            printf("Composition based function: %s. Ndim = %d, Npar = %d \n", GetName(), GetNdim(), GetNpar());
            if (!fComposition)
                printf("fComposition not found!\n"); // this would be bad
        }
        else
        {
            if (fFunctor)
                printf("Compiled based function: %s  based on a functor object.  Ndim = %d, Npar = %d\n", GetName(),
                       GetNdim(), GetNpar());
            else
            {
                printf("Function based on a list of points from a compiled based function: %s.  Ndim = %d, Npar = %d, Npx "
                       "= %zu\n",
                       GetName(), GetNdim(), GetNpar(), fSave.size());
                if (fSave.empty())
                    Warning("Print", "Function %s is based on a list of points but list is empty", GetName());
            }
        }
        TString opt(option);
        opt.ToUpper();
        if (opt.Contains("V"))
        {
            // print list of parameters
            if (fNpar > 0)
            {
                printf("List of  Parameters: \n");
                for (int i = 0; i < fNpar; ++i)
                    printf(" %20s =  %10f \n", GetParName(i), GetParameter(i));
            }
            //  if (!fSave.empty()) {                                               /// skip printing list of saved points: 1000 points flood the output!
            //     // print list of saved points
            //     printf("List of  Saved points (N=%d): \n", int(fSave.size()));
            //     for (auto &x : fSave)
            //        printf("( %10f )  ", x);
            //     printf("\n");
            //  }
        }
    }
    if (fHistogram)
    {
        printf("Contained histogram\n");
        fHistogram->Print(option);
    }
}

/// Get avg mean and sigma from both gaussians - adapted for both gaussians with same mean
/**
 * @brief Calculates weighted mean and sigma from a double Gaussian fit with background
 *
 * Extracts parameters from a double Gaussian fit function of the form:
 * [6]+[7]*x+[8]*x*x+[0]*TMath::Gaus(x,[1],[2])+[3]*TMath::Gaus(x,[4],[5])
 * and calculates the weighted mean and sigma values with their uncertainties.
 *
 * @param f_doubleGaus Pointer to the double Gaussian fit function
 * @param lFitResultPtr Pointer to the fit result containing covariance matrix
 * @param mean Output parameter for the calculated weighted mean
 * @param mean_err Output parameter for the error on weighted mean
 * @param sigma Output parameter for the calculated weighted sigma
 * @param sigma_err Output parameter for the error on weighted sigma
 *
 * The function calculates weighted averages using:
 * mean = (N1*mu1 + N2*mu2)/(N1 + N2)
 * sigma = (N1*sigma1 + N2*sigma2)/(N1 + N2)
 *
 * Errors are calculated using either:
 * - Full error propagation with covariance matrix if available
 * - Simplified error estimation if covariance matrix is not available
 */
inline void GetMeanSigmaDG(TF1 *f_doubleGaus, TFitResultPtr lFitResultPtr, Double_t &mean, Double_t &mean_err, Double_t &sigma, Double_t &sigma_err)
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

    Double_t mu_wa_err = 0.;
    Double_t sigma_wa_err = 0.;

    Double_t sa = N1a + N2a;
    Double_t wa_mu = N1a * mu1a + N2a * mu2a;
    Double_t wa_sigma = N1a * sigma1a + N2a * sigma2a;

    if (cova.NonZeros() > 0)
    {
        Double_t mu_wa_step = pow((mu1a - mu2a), 2) * (pow(N1a, 2) * cova(3, 3) + pow(N2a, 2) * cova(0, 0)) + 2 * cova(0, 3) * (wa_mu - sa * mu1a) * (wa_mu - sa * mu2a) + pow(sa, 2) * (pow(N1a, 2) * cova(1, 1) + pow(N2a, 2) * cova(4, 4) + 2 * N1a * N2a * cova(1, 4)) - 2 * sa * (N1a * (cova(0, 1) * (wa_mu - sa * mu1a) + cova(3, 1) * (wa_mu - sa * mu2a)) + N2a * (cova(0, 4) * (wa_mu - sa * mu1a) + cova(3, 4) * (wa_mu - sa * mu2a)));
        mu_wa_err = sqrt(mu_wa_step / pow(sa, 4));

        Double_t sigma_wa_step = pow((sigma1a - sigma2a), 2) * (pow(N1a, 2) * cova(3, 3) + pow(N2a, 2) * cova(0, 0)) + 2 * cova(0, 3) * (wa_sigma - sa * sigma1a) * (wa_sigma - sa * sigma2a) + pow(sa, 2) * (pow(N1a, 2) * cova(2, 2) + pow(N2a, 2) * cova(5, 5) + 2 * N1a * N2a * cova(2, 5)) - 2 * sa * (N1a * (cova(0, 2) * (wa_sigma - sa * sigma1a) + cova(3, 2) * (wa_sigma - sa * sigma2a)) + N2a * (cova(0, 5) * (wa_sigma - sa * sigma1a) + cova(3, 5) * (wa_sigma - sa * sigma2a)));
        sigma_wa_err = sqrt(sigma_wa_step / pow(sa, 4));
        Info("GetMeanSigmaDG", "%s: Avg Fit Mass = %.3f +/- %f; Avg Fit Sigma = %f +/- %f", f_doubleGaus->GetName(), mu_wa, mu_wa_err, sigma_wa, sigma_wa_err);
    }
    else
    {
        Double_t mu1e = f_doubleGaus->GetParError(1);
        Double_t mu2e = f_doubleGaus->GetParError(1);
        Double_t sigma1e = f_doubleGaus->GetParError(2);
        Double_t sigma2e = f_doubleGaus->GetParError(4);
        mu_wa_err = sqrt(pow((N1a * mu1e), 2) + pow((N2a * mu2e), 2)) / sa;
        sigma_wa_err = sqrt(pow((N1a * sigma1e), 2) + pow((N2a * sigma2e), 2)) / sa;
        Warning("GetMeanSigmaDG", "%s: Covariance matrix not available! Approximating errors. Avg Fit Mass = %.3f +/- %f; Avg Fit Sigma = %f +/- %f", f_doubleGaus->GetName(), mu_wa, mu_wa_err, sigma_wa, sigma_wa_err);
    }
    mean = mu_wa;
    mean_err = mu_wa_err;
    sigma = sigma_wa;
    sigma_err = sigma_wa_err;
}

#endif // CASCADEUTILS_H