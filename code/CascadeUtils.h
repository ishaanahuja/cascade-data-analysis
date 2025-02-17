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

enum particles
{
    kXi,
    kOm,
    kNumPart
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

// chosen marker style palette
Int_t markerStyles[] = {kCircle, kFullSquare, kFullTriangleUp, kFullTriangleDown, kFullStar, kFullDiamond, kFullCross, kFullDoubleDiamond, kFullCrossX, kFullFourTrianglesX, kFullThreeTriangles, kFullCircle};

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
 * @brief Opens a ROOT file with specified options
 *
 * This function opens a ROOT file and performs basic validity checks.
 * If the file name is provided without the .root extension, it will be automatically added.
 *
 * @param fileName The name/path of the ROOT file to open
 * @param options File access options (default: "READ")
 * @return TFile* Pointer to the opened file, nullptr if opening fails
 *
 * @note The function performs the following checks:
 *       - Empty file name
 *       - File opening success
 *       - Zombie file status
 *
 * @warning The caller is responsible for closing and deleting the returned TFile pointer
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
 */
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

inline Double_t DoubleGausPol2(double *x, double *par) // with same mean parameter for both gaussians
{
    return par[0] * TMath::Gaus(x[0], par[1], par[2]) + par[3] * TMath::Gaus(x[0], par[1], par[4]) + par[5] + par[6] * x[0] + par[7] * x[0] * x[0];
}

inline Double_t DoubleGausPol3(double *x, double *par)
{
    return par[0] * TMath::Gaus(x[0], par[1], par[2]) + par[3] * TMath::Gaus(x[0], par[4], par[5]) + par[6] + par[7] * x[0] + par[8] * x[0] * x[0] + par[9] * x[0] * x[0] * x[0];
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