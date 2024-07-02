#ifndef CASCADEUTILS_H
#define CASCADEUTILS_H

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

// chosen marker style palette
Int_t markerStyles[] = {4, 21, 22, 23, 29, 33, 34, 43, 47, 41, 20};

/// @brief Open a root file and return the TFile's handle.
/// @param fileName Full path to the file, appends .root if absent
/// @param options Open a file in Read, Write, Recreate mode etc. Check options in TFile::Open
/// @return Returns the pointer to file.
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

/// @brief Accepts path where output should be stored. Creates folder if it doesn't exist.
/// @param folderName NOTE: if folderName is empty, creates a folder called "Output" in CWD.
/// @return Full path of output folder
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

/// @brief Draw a THStack on the current canvas, with legend and axes titles
/// @param c current canvas to use (*)
/// @param hs THStack object to paint (*)
/// @param setLogY Set Y-axis as Logarithmic
/// @param yAxisTitle Title for y-axis
/// @param xAxisTitle Title for x-axis
inline void PaintStack(TCanvas &c, THStack &hs, Bool_t setLogY = kTRUE, TString yAxisTitle = "#frac{1}{#it{N}_{inel}} #frac{d#it{N}}{d#it{p}_{T}}", TString xAxisTitle = "#it{p}_{T} (GeV/c)")
{
    c.Clear();
    c.cd();
    hs.Draw("plc pmc nostack");
    if (setLogY)
        gPad->SetLogy();
    gPad->BuildLegend(0.9, 0.6, 1., 1., "");
    hs.GetXaxis()->SetTitle(xAxisTitle.Data());
    hs.GetYaxis()->SetTitle(yAxisTitle.Data());
    c.Modified();
    c.ForceUpdate();
}

/// @brief Save image to directory.
/// @param outputFolder Parent directory of images
/// @param imageFolder Sub-directory for image groups
/// @param imageName name of the image file
/// @param imageFormat image format, e.g. png, pdf, jpg ...
inline void SaveImage(TString outputFolder, TString imageFolder, TString imageName, TString imageFormat)
{
    TString imagePath = SetOutputFolder(outputFolder + "/" + imageFolder);
    gPad->Print(TString::Format("%s/%s.%s", imagePath.Data(), imageName.Data(), imageFormat.Data()), imageFormat.Data());
    gSystem->Chmod(TString::Format("%s/%s.%s", imagePath.Data(), imageName.Data(), imageFormat.Data()), 0755);
}

inline Double_t DoubleGausPol2(double *x, double *par) // with same mean parameter for both gaussians
{
    return par[0] * TMath::Gaus(x[0], par[1], par[2]) + par[3] * TMath::Gaus(x[0], par[1], par[4]) + par[5] + par[6] * x[0] + par[7] * x[0] * x[0];
}

inline Double_t DoubleGausPol3(double *x, double *par)
{
    return par[0] * TMath::Gaus(x[0], par[1], par[2]) + par[3] * TMath::Gaus(x[0], par[4], par[5]) + par[6] + par[7] * x[0] + par[8] * x[0] * x[0] + par[9] * x[0] * x[0] * x[0];
}

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

inline Double_t GausPol2(double *x, double *par)
{
    return par[0] * TMath::Gaus(x[0], par[1], par[2]) + par[3] * x[0] * x[0] + par[4] * x[0] + par[5];
}

inline Double_t GetChi2(TFitResultPtr fFitResult)
{
    Double_t chi2 = fFitResult->Chi2();
    if (TMath::IsNaN(chi2))
        return 0;
    else
        return chi2;
}

inline Double_t GetNdf(TFitResultPtr fFitResult)
{
    Double_t ndf = fFitResult->Ndf();
    if (TMath::IsNaN(ndf))
        return 0;
    else
        return ndf;
}

inline Double_t GetReducedChi2(TFitResultPtr fFitResult)
{
    Double_t reducedChi2 = fFitResult->Chi2() / fFitResult->Ndf();
    if (TMath::IsNaN(reducedChi2))
        return 0;
    else
        return reducedChi2;
}

inline Double_t GetProb(TFitResultPtr fFitResult)
{
    Double_t prob = fFitResult->Prob();
    if (TMath::IsNaN(prob))
        return 0;
    else
        return prob;
}

inline Double_t GetGeneratedParticles(TH2 *h2, Int_t binsMC[], Double_t &genErr)
{
    Double_t gen = h2->IntegralAndError(binsMC[0], binsMC[1], binsMC[2], binsMC[3], genErr);
    return gen;
}

inline Double_t ErrorInRatio(Double_t A, Double_t Aerr, Double_t B, Double_t Berr)
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

/// Get avg mean and sigma from both gaussians - adapted for both gaussians with same mean
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
}

#endif // CASCADEUTILS_H