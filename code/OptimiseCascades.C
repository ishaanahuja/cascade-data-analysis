#include <TROOT.h>

#include <TStyle.h>
#include <TLine.h>
#include <TLegend.h>
#include <TGraphErrors.h>

#include "CascadeUtils.h"

enum particles
{
    kXi,
    kOm,
    kNumPart
};
enum ptInterval
{
    kLow,
    kMid,
    kHigh,
    kNumPtInterval
};
enum cuts
{
    kCascTransDecayRadius,
    kV0TransDecayRadius,
    kDcaBachToPv,
    kDcaV0ToPv,
    kDcaMesV0ToPv,
    kDcaBarV0ToPv,
    kDcaV0Daughters,
    kDcaBachToV0,
    kCascCosPa,
    kV0CosPa,
    kDefaultVals,
    kV0InvMassWindow,
    kTpcDedxPidSigma,
    kDeviationPropLifetime,
    kLeastTpcClusters,
    kCompetingCascRejectOm,
    kLeastCRows,
    kLeastCRowsOvF,
    kTrackLengthCut,
    kNumCuts
};

struct cutInfo
{
    Int_t cutId = 0;
    Bool_t isTopo = kFALSE;
    TString cutName = "";
    Double_t defVal[kNumPart][kNumPtInterval] = {0.}; // default value for each cut [particle] [pt_interval]
    Double_t optVal[kNumPart][kNumPtInterval] = {0.}; // optimum value for each cut [particle] [pt_interval]
    Int_t nVars = 0;
    std::vector<Double_t> cutVals;                           // cutVals [var]
    std::vector<std::vector<std::vector<Int_t>>> isValDefAt; // isValDefAt[particle][pt_interval[var]: true if cutVals[var] == defVal[var]
    std::vector<std::vector<std::vector<Double_t>>> sigVal;  // sigVal [particle] [var] [pt_interval]
    std::vector<std::vector<std::vector<Double_t>>> sigErr;  // sigErr [particle] [var] [pt_interval]
    std::vector<std::vector<std::vector<Double_t>>> bgVal;   // bgVal [particle] [var] [pt_interval]
    std::vector<std::vector<std::vector<Double_t>>> bgErr;   // bgErr [particle] [var] [pt_interval]

    std::vector<std::vector<std::vector<Double_t>>> significance;    // significance [particle] [pt_interval] [var]
    std::vector<std::vector<std::vector<Double_t>>> significanceErr; // significanceErr [particle] [pt_interval] [var]

    std::vector<std::vector<Int_t>> loosestVar;                 // index of which variation is loosest (max sigVal) [particle] [pt_interval]
    std::vector<std::vector<Double_t>> maxSigVal;               // maxSigVal [particle] [pt_interval]
    std::vector<std::vector<Double_t>> maxSigErr;               // maxSigErr [particle] [pt_interval]
    std::vector<std::vector<std::vector<Double_t>>> sigLoss;    // sigLoss [particle] [pt_interval] [var]
    std::vector<std::vector<std::vector<Double_t>>> sigLossErr; // sigLossErr [particle] [pt_interval] [var]

    TGraphErrors *g_significance[kNumPart][kNumPtInterval]; // significance graph [particle] [pt_interval]
    TGraphErrors *g_sigLoss[kNumPart][kNumPtInterval];      // signal loss graph [particle] [pt_interval]
};

double fPtInterval[2][4] = {{0.8, 1.9, 3.1, 5.5}, {0.9, 2., 2.9, 5}}; // values for low, mid, high pt intervals for Xi and Omega
int fPtIntervalBins[2][4] = {{0, 5, 11, 13}, {0, 2, 4, 5}};           // binning for low, mid, high pt intervals for Xi and Omega

inline Bool_t CheckFileName(TString fileName, Int_t &fileCutNum, Int_t &fileVarNum);
void SetCutInfo(cutInfo &cut, Int_t cutId, Bool_t isTopo, TString cutName, Int_t nVars, Double_t lowVal, Double_t highVal);
void SetDefCutValue(cutInfo &cut, double cutVal, int particle = -1, int ptInterval = -1);
void GetSignalBg(cutInfo &cut, TH1 &h, Int_t varNum, Int_t particle, Int_t ptBin, Double_t sig[], Double_t bg[], Double_t errSigSq[], Double_t errBgSq[]);
void CalculateSignificance(cutInfo &cut, cutInfo &defCut);
void CalculateSigLoss(cutInfo &cut);
void PrintCutInfo(cutInfo &cut);
void Initialise(std::vector<cutInfo> &cuts);
void DrawAndSave(TGraphErrors *graph1, TGraphErrors *graph2, Double_t xDefLine, Double_t xOptLine, TString outputFolder, TString imageFolder, TString imageName, TString imageFormat);

void FindOptimumDefVal(cutInfo &cut, cutInfo &defCut);

int OptimiseCascades(std::string inputFilename = "/var/home/ishaan/Work/git/analysis/results/230924_results_6Runs/230924_results_6Runs_fileList.txt", TString inputFileDefaultCuts = "/var/home/ishaan/Work/git/analysis/results/230924_results_6Runs/230924_6Runs_h3_ptmasscent_def.root", TString outputFilename = "significance_test2.root", TString outputFolder = "031024_images_significance", Bool_t saveImages = kTRUE, TString imageFormat = "png", Int_t verbosity = kInfo)
{
    gErrorIgnoreLevel = verbosity;

    outputFolder = SetOutputFolder(outputFolder);

    // remove ownership of objects from files -> we can delete the file ptr later
    TH1::AddDirectory(kFALSE);

    std::vector<std::string> inputFileList = GetFileList(inputFilename);
    if (inputFileList.empty())
    {
        Error("OptimiseCascades: FileList", "'%s': Could not read the file! Aborting.", inputFilename.data());
        return 3;
    }

    std::vector<cutInfo> cuts(kNumCuts);
    Initialise(cuts);

    TH1 *resultParXiC_pt[fNptbins_Xi];
    TH1 *resultParOmC_pt[fNptbins_Om];
    // TH1 *h_multBinEntries_Xi;
    // TH1 *h_multBinEntries_Om;

    // Double_t numVars[kNumCuts] = {6, 4, 3, 3, 4, 4, 4, 4, 10, 10, 9, 9, 9, 13, 10, 10, 15, 15};
    // Double_t sigVal[kNumCuts][kNumPtInterval];
    TString fileName = "";
    Int_t fileCutNum = -1;
    Int_t fileVarNum = -1;
    TFile *currentFile = nullptr;

    /// temporary arrays for GetSignalBg
    /// handle case for default cuts file - cutId = 10
    Double_t sigDef[kNumPtInterval] = {0.};
    Double_t bgDef[kNumPtInterval] = {0.};
    Double_t errSigDefSq[kNumPtInterval] = {0.};
    Double_t errBgDefSq[kNumPtInterval] = {0.};

    currentFile = OpenFile(inputFileDefaultCuts);

    for (Int_t ptBinXi = 0; ptBinXi < fNptbins_Xi; ptBinXi++)
    {
        resultParXiC_pt[ptBinXi] = (TH1 *)currentFile->FindObjectAny(TString::Format(("resultParXiC_pt[%d]"), ptBinXi));

        if ((!resultParXiC_pt[ptBinXi]) || resultParXiC_pt[ptBinXi]->IsZombie())
        {
            Error("OptimiseCascades", "Default_cuts: resultParXiC_pt[%d] not found. Skipping Xi analysis for file '%s'", ptBinXi, fileName.Data());
            break;
        }
        GetSignalBg(cuts[kDefaultVals], *resultParXiC_pt[ptBinXi], 0, kXi, ptBinXi, sigDef, bgDef, errSigDefSq, errBgDefSq);
    }

    for (Int_t ptBinOm = 0; ptBinOm < fNptbins_Om; ptBinOm++)
    {

        resultParOmC_pt[ptBinOm] = (TH1 *)currentFile->FindObjectAny(TString::Format(("resultParOmC_pt[%d]"), ptBinOm));

        if ((!resultParOmC_pt[ptBinOm]) || resultParOmC_pt[ptBinOm]->IsZombie())
        {
            Error("OptimiseCascades", "Default_cuts: resultParOmC_pt[%d] not found. Skipping Omega analysis for file '%s'", ptBinOm, fileName.Data());
            break;
        }
        GetSignalBg(cuts[kDefaultVals], *resultParOmC_pt[ptBinOm], 0, kOm, ptBinOm, sigDef, bgDef, errSigDefSq, errBgDefSq);
    }

    currentFile->Close(); // input ended for this file

    /// Loop over all files in the variation list
    for (size_t iFile = 0; iFile < inputFileList.size(); iFile++)
    {
        fileName = inputFileList.at(iFile);
        fileCutNum = -1;
        fileVarNum = -1;

        Bool_t isFileNameValid = CheckFileName(fileName, fileCutNum, fileVarNum);
        if (!isFileNameValid)
            continue;

        currentFile = OpenFile(fileName);

        /// temporary arrays for GetSignalBg
        Double_t sig[kNumPtInterval] = {0.};
        Double_t bg[kNumPtInterval] = {0.};
        Double_t errSigSq[kNumPtInterval] = {0.};
        Double_t errBgSq[kNumPtInterval] = {0.};

        /// XI input
        for (Int_t ptBinXi = 0; ptBinXi < fNptbins_Xi; ptBinXi++)
        {

            resultParXiC_pt[ptBinXi] = (TH1 *)currentFile->FindObjectAny(TString::Format(("resultParXiC_pt[%d]"), ptBinXi));

            if ((!resultParXiC_pt[ptBinXi]) || resultParXiC_pt[ptBinXi]->IsZombie())
            {
                Error("OptimiseCascades", "%d_%d: resultParXiC_pt[%d] not found. Skipping Xi analysis for file '%s'", fileCutNum, fileVarNum, ptBinXi, fileName.Data());
                break;
            }
            GetSignalBg(cuts[fileCutNum], *resultParXiC_pt[ptBinXi], fileVarNum, kXi, ptBinXi, sig, bg, errSigSq, errBgSq);
        }

        /// OMEGA input
        for (Int_t ptBinOm = 0; ptBinOm < fNptbins_Om; ptBinOm++)
        {

            resultParOmC_pt[ptBinOm] = (TH1 *)currentFile->FindObjectAny(TString::Format(("resultParOmC_pt[%d]"), ptBinOm));

            if ((!resultParOmC_pt[ptBinOm]) || resultParOmC_pt[ptBinOm]->IsZombie())
            {
                Error("OptimiseCascades", "%d_%d: resultParOmC_pt[%d] not found. Skipping Omega analysis for file '%s'", fileCutNum, fileVarNum, ptBinOm, fileName.Data());
                break;
            }
            GetSignalBg(cuts[fileCutNum], *resultParOmC_pt[ptBinOm], fileVarNum, kOm, ptBinOm, sig, bg, errSigSq, errBgSq);
        }

        currentFile->Close(); // input ended for this file
    }
    delete currentFile;

    Info("OptimiseCascades", "Input ended.");

    /// Saving output :
    TFile *outputFile = OpenFile(outputFilename, "RECREATE");
    outputFile->mkdir("sig_XiC");
    outputFile->mkdir("sig_OmC");
    outputFile->mkdir("sLoss_XiC");
    outputFile->mkdir("sLoss_OmC");
    TString particleTitle;
    TString particleName;

    for (Int_t iCut = 0; iCut < kNumCuts; iCut++)
    {
        CalculateSignificance(cuts[iCut], cuts[kDefaultVals]);
        CalculateSigLoss(cuts[iCut]);
        FindOptimumDefVal(cuts[iCut], cuts[kDefaultVals]);
        PrintCutInfo(cuts[iCut]);
    }

    /// Create graphs
    for (Int_t iCut = 0; iCut < kNumCuts; iCut++)
    {
        if (iCut == kDefaultVals) // skip generating graphs for default cut object - values will be added in all graphs
            continue;

        for (Int_t iPart = 0; iPart < kNumPart; iPart++)
        {
            if (iPart == kXi)
            {
                particleName = "XiC";
                particleTitle = "#Xi^{+} + #Xi^{-}";
            }
            else
            {
                particleName = "OmC";
                particleTitle = "#Omega^{+} + #Omega^{-}";
            }
            for (Int_t iPt = 0; iPt < kNumPtInterval; iPt++)
            {
                outputFile->cd("sig_" + particleName);

                /// Create significance TGraphErrors for each cut, particle, pt interval and populate with significance values on y-axis and cut values on x-axis
                cuts[iCut].g_significance[iPart][iPt] = new TGraphErrors(cuts[iCut].cutVals.size(), cuts[iCut].cutVals.data(), cuts[iCut].significance[iPart][iPt].data(), nullptr, cuts[iCut].significanceErr[iPart][iPt].data());

                bool defValFound = false;

                for (Int_t iPoint = 0; iPoint < cuts[iCut].g_significance[iPart][iPt]->GetN(); iPoint++)
                {
                    if (TMath::Abs(cuts[iCut].g_significance[iPart][iPt]->GetPointX(iPoint) - cuts[iCut].defVal[iPart][iPt]) < 1e-9)
                    {
                        Printf("Default value (%.4f) for cut %s found at Var(point): %d, Part: %s, Pt: %.1f, previous significance = %.2f+/-%.2f", cuts[iCut].defVal[iPart][iPt], cuts[iCut].cutName.Data(), iPoint, particleName.Data(), fPtInterval[iPart][iPt], cuts[iCut].g_significance[iPart][iPt]->GetPointY(iPoint), cuts[iCut].g_significance[iPart][iPt]->GetErrorY(iPoint));

                        /// Update the significance for the default X value with the actual value for default cuts
                        cuts[iCut].g_significance[iPart][iPt]->SetPointY(iPoint, cuts[kDefaultVals].significance[iPart][iPt][0]);
                        cuts[iCut].g_significance[iPart][iPt]->SetPointError(iPoint, 0.0, cuts[kDefaultVals].significanceErr[iPart][iPt][0]);
                        defValFound = true;
                        break;
                    }
                }

                if (!defValFound)
                {
                    Printf("Default value (%.4f) for cut %s, Part: %s, Pt: %.1f not found in the graph, adding manually", cuts[iCut].defVal[iPart][iPt], cuts[iCut].cutName.Data(), particleName.Data(), fPtInterval[iPart][iPt]);
                    cuts[iCut].g_significance[iPart][iPt]->AddPoint(cuts[iCut].defVal[iPart][iPt], cuts[kDefaultVals].significance[iPart][iPt][0]);
                    cuts[iCut].g_significance[iPart][iPt]->SetPointError(cuts[iCut].g_significance[iPart][iPt]->GetN() - 1, 0.0, cuts[kDefaultVals].significanceErr[iPart][iPt][0]);
                }

                /// Remove points with significance = 0 from the graph: default values will be here instead
                // Int_t iPoint = 0;
                // while (iPoint < cuts[iCut].g_significance[iPart][iPt]->GetN())
                // {
                //     if (cuts[iCut].g_significance[iPart][iPt]->GetPointY(iPoint) == 0.0)
                //         cuts[iCut].g_significance[iPart][iPt]->RemovePoint(iPoint);
                //     else
                //         iPoint++;
                // }

                cuts[iCut].g_significance[iPart][iPt]->Sort();

                cuts[iCut].g_significance[iPart][iPt]->SetNameTitle(TString::Format("g_significance_%d_%d_%d", iCut, iPart, iPt), TString::Format("%s: %.1f < #it{p}_{T} < %.1f", particleTitle.Data(), fPtInterval[iPart][iPt], fPtInterval[iPart][iPt + 1]));
                cuts[iCut].g_significance[iPart][iPt]->GetXaxis()->SetTitle(cuts[iCut].cutName);
                cuts[iCut].g_significance[iPart][iPt]->GetYaxis()->SetTitle("Signal Significance (#frac{S}{#sqrt{S+B}})");
                cuts[iCut].g_significance[iPart][iPt]->Write();

                outputFile->cd("sLoss_" + particleName);

                /// Create signal loss TGraphErrors for each cut, particle, pt interval and populate with signal loss values on y-axis and cut values on x-axis
                cuts[iCut].g_sigLoss[iPart][iPt] = new TGraphErrors(cuts[iCut].cutVals.size(), cuts[iCut].cutVals.data(), cuts[iCut].sigLoss[iPart][iPt].data(), nullptr, nullptr); // cuts[iCut].sigLossErr[iPart][iPt].data());

                defValFound = false;

                for (Int_t iPoint = 0; iPoint < cuts[iCut].g_sigLoss[iPart][iPt]->GetN(); iPoint++)
                {
                    if (TMath::Abs(cuts[iCut].g_sigLoss[iPart][iPt]->GetPointX(iPoint) - cuts[iCut].defVal[iPart][iPt]) < 1e-9)
                    {
                        Printf("Default value (%.4f) for cut %s found at Var(point): %d, Part: %s, Pt: %.1f, previous sigLoss = %.2f+/-%.2f, actual = %.2f+/-%.2f", cuts[iCut].defVal[iPart][iPt], cuts[iCut].cutName.Data(), iPoint, particleName.Data(), fPtInterval[iPart][iPt], cuts[iCut].g_sigLoss[iPart][iPt]->GetPointY(iPoint), cuts[iCut].g_sigLoss[iPart][iPt]->GetErrorY(iPoint), cuts[kDefaultVals].sigLoss[iPart][iPt][0], cuts[kDefaultVals].sigLossErr[iPart][iPt][0]);

                        /// Update the signal loss for the default X value with the actual value for default cuts
                        cuts[iCut].g_sigLoss[iPart][iPt]->SetPointY(iPoint, cuts[kDefaultVals].sigLoss[iPart][iPt][0]);
                        cuts[iCut].g_sigLoss[iPart][iPt]->SetPointError(iPoint, 0.0, cuts[kDefaultVals].sigLossErr[iPart][iPt][0]);
                        defValFound = true;
                        break;
                    }
                }

                if (!defValFound)
                {
                    Printf("Default value (%.4f) for cut %s, Part: %s, Pt: %.1f not found in the graph, adding manually", cuts[iCut].defVal[iPart][iPt], cuts[iCut].cutName.Data(), particleName.Data(), fPtInterval[iPart][iPt]);
                    cuts[iCut].g_sigLoss[iPart][iPt]->AddPoint(cuts[iCut].defVal[iPart][iPt], cuts[kDefaultVals].sigLoss[iPart][iPt][0]);
                    cuts[iCut].g_sigLoss[iPart][iPt]->SetPointError(cuts[iCut].g_sigLoss[iPart][iPt]->GetN() - 1, 0.0, cuts[kDefaultVals].sigLossErr[iPart][iPt][0]);
                }

                // /// Remove points with signal loss = 100 from the graph
                // iPoint = 0;
                // while (iPoint < cuts[iCut].g_sigLoss[iPart][iPt]->GetN())
                // {
                //     if (cuts[iCut].g_sigLoss[iPart][iPt]->GetPointY(iPoint) == 100.0)
                //         cuts[iCut].g_sigLoss[iPart][iPt]->RemovePoint(iPoint);
                //     else
                //         iPoint++;
                // }

                cuts[iCut].g_sigLoss[iPart][iPt]->Sort();

                cuts[iCut].g_sigLoss[iPart][iPt]->SetNameTitle(TString::Format("g_sigLoss_%d_%d_%d", iCut, iPart, iPt), TString::Format("%s: %.1f < #it{p}_{T} < %.1f", particleTitle.Data(), fPtInterval[iPart][iPt], fPtInterval[iPart][iPt + 1]));
                cuts[iCut].g_sigLoss[iPart][iPt]->GetXaxis()->SetTitle(cuts[iCut].cutName);
                cuts[iCut].g_sigLoss[iPart][iPt]->GetYaxis()->SetTitle("Signal Loss (%)");
                if (saveImages)
                    DrawAndSave(cuts[iCut].g_significance[iPart][iPt], cuts[iCut].g_sigLoss[iPart][iPt], cuts[iCut].defVal[iPart][iPt], cuts[iCut].optVal[iPart][iPt], outputFolder, "multi_" + particleName, cuts[iCut].cutName + TString::Format("_%d_%d", iPart, iPt), imageFormat);
                cuts[iCut].g_sigLoss[iPart][iPt]->Write();
            }
        }
    }

    Printf("FINISHED!!! Output saved to '%s', images saved to '%s'.", outputFile->GetName(), outputFolder.Data());
    outputFile->Close();
    delete outputFile;
    return 0;
}

inline void Initialise(std::vector<cutInfo> &cuts)
{
    Int_t all = -1;

    /// Set cut information - cutId, isTopo, cutName, nVars, cutValues (x-axis), initialise sigVal, sigErr, bgVal, bgErr vectors
    SetCutInfo(cuts[0], kCascTransDecayRadius, kTRUE, "CascTransDecayRadius", 9, 0.4, 2.);
    SetCutInfo(cuts[1], kV0TransDecayRadius, kTRUE, "V0TransDecayRadius", 13, 1., 7.);
    SetCutInfo(cuts[2], kDcaBachToPv, kTRUE, "DcaBachToPv", 10, 0.02, 0.2);
    SetCutInfo(cuts[3], kDcaV0ToPv, kTRUE, "DcaV0ToPv", 10, 0.02, 0.2);
    SetCutInfo(cuts[4], kDcaMesV0ToPv, kTRUE, "DcaMesV0ToPv", 15, 0.02, 0.3);
    SetCutInfo(cuts[5], kDcaBarV0ToPv, kTRUE, "DcaBarV0ToPv", 15, 0.02, 0.3);
    SetCutInfo(cuts[6], kDcaV0Daughters, kTRUE, "DcaV0Daughters", 9, 0.4, 2.);
    SetCutInfo(cuts[7], kDcaBachToV0, kTRUE, "DcaBachToV0", 9, 0.4, 2.);
    SetCutInfo(cuts[8], kCascCosPa, kTRUE, "CascCosPa", 10, 0.95, 0.995);
    SetCutInfo(cuts[9], kV0CosPa, kTRUE, "V0CosPa", 10, 0.95, 0.995);
    SetCutInfo(cuts[10], kDefaultVals, kFALSE, "Default", 1, 0.0, 0.0);
    SetCutInfo(cuts[11], kV0InvMassWindow, kFALSE, "V0InvMassWindow", 4, 0.003, 0.012);
    SetCutInfo(cuts[12], kTpcDedxPidSigma, kFALSE, "TpcDedxPidSigma", 6, 2, 7);
    SetCutInfo(cuts[13], kDeviationPropLifetime, kFALSE, "DeviationPropLifetime", 4, 2, 5);
    SetCutInfo(cuts[14], kLeastTpcClusters, kFALSE, "LeastTpcClusters", 3, 75, 85);
    SetCutInfo(cuts[15], kCompetingCascRejectOm, kFALSE, "CompetingCascRejectOm", 3, 0.0, 0.010);
    SetCutInfo(cuts[16], kLeastCRows, kFALSE, "LeastCRows", 4, 60, 90);
    SetCutInfo(cuts[17], kLeastCRowsOvF, kFALSE, "LeastCRowsOvF", 4, 0.75, 0.9);
    SetCutInfo(cuts[18], kTrackLengthCut, kFALSE, "TrackLengthCut", 4, 0, 3);

    /// Set default cut values for each cut and particle-pt interval
    SetDefCutValue(cuts[kTpcDedxPidSigma], 4.0, all, all);
    SetDefCutValue(cuts[kDeviationPropLifetime], 3.0, all, all);
    SetDefCutValue(cuts[kLeastTpcClusters], 70.0, all, all);
    SetDefCutValue(cuts[kCompetingCascRejectOm], 0.008, all, all);
    SetDefCutValue(cuts[kLeastCRows], 70.0, all, all);
    SetDefCutValue(cuts[kLeastCRowsOvF], 0.8, all, all);
    SetDefCutValue(cuts[kTrackLengthCut], 1.0, all, all);
    SetDefCutValue(cuts[kV0InvMassWindow], 0.008, all, all);
    SetDefCutValue(cuts[kCascTransDecayRadius], 0.6, kXi, all);
    SetDefCutValue(cuts[kCascTransDecayRadius], 0.5, kOm, all);
    SetDefCutValue(cuts[kV0TransDecayRadius], 1.2, kXi, all);
    SetDefCutValue(cuts[kV0TransDecayRadius], 1.1, kOm, all);
    SetDefCutValue(cuts[kDcaBachToPv], 0.04, all, all);
    SetDefCutValue(cuts[kDcaV0ToPv], 0.06, all, all);
    SetDefCutValue(cuts[kDcaMesV0ToPv], 0.04, all, all);
    SetDefCutValue(cuts[kDcaBarV0ToPv], 0.03, all, all);
    SetDefCutValue(cuts[kDcaV0Daughters], 1.7, all, kLow);
    SetDefCutValue(cuts[kDcaV0Daughters], 1.6, all, kMid);
    SetDefCutValue(cuts[kDcaV0Daughters], 1.5, all, kHigh);
    SetDefCutValue(cuts[kDcaBachToV0], 1.6, all, kLow);
    SetDefCutValue(cuts[kDcaBachToV0], 1.5, all, kMid);
    SetDefCutValue(cuts[kDcaBachToV0], 1.4, all, kHigh);
    SetDefCutValue(cuts[kCascCosPa], 0.96, all, kLow);
    SetDefCutValue(cuts[kCascCosPa], 0.97, all, kMid);
    SetDefCutValue(cuts[kCascCosPa], 0.97, all, kHigh);
    SetDefCutValue(cuts[kV0CosPa], 0.98, kXi, all);
    SetDefCutValue(cuts[kV0CosPa], 0.97, kOm, kLow);
    SetDefCutValue(cuts[kV0CosPa], 0.98, kOm, kMid);
    SetDefCutValue(cuts[kV0CosPa], 0.98, kOm, kHigh);
}

inline Bool_t CheckFileName(TString fileName, Int_t &fileCutNum, Int_t &fileVarNum)
{

    TString token = "";

    /// Check File extension before proceeding
    if (!fileName.EndsWith(".root"))
    {
        Warning("OptimiseCascades: CheckFileName", "'%s': Not a '.root' file. Skipping ...", fileName.Data());
        return kFALSE;
    }

    TObjArray *fileNameTokens = (TObjArray *)fileName.Tokenize("_.");
    if (fileNameTokens->GetEntries() < 4)
    {
        Error("OptimiseCascades: CheckFileName", "'%s': Not enough tokens in the file name. Skipping ...", fileName.Data());
        fileNameTokens->Delete();
        return kFALSE;
    }

    token = fileNameTokens->At(fileNameTokens->GetLast() - 2)->GetName(); /// third last token is the cut number
    if (token.IsDigit())
        fileCutNum = token.Atoi();

    token = fileNameTokens->At(fileNameTokens->GetLast() - 1)->GetName(); /// second last token is the var number
    if (token.IsDigit())
        fileVarNum = token.Atoi();

    fileNameTokens->Delete(); // done with tokens - cleanup

    if (fileCutNum < 0 || fileVarNum < 0)
    {
        Error("OptimiseCascades: CheckFileName", "'%s': Invalid cut or var number. Skipping ...", fileName.Data());
        return kFALSE;
    }

    return kTRUE;
}

void SetCutInfo(cutInfo &cut, Int_t cutId, Bool_t isTopo, TString cutName, Int_t nVars, Double_t lowVal, Double_t highVal)
{
    // Set cut information
    cut.cutId = cutId;
    cut.isTopo = isTopo;
    cut.cutName = cutName;
    cut.nVars = nVars;

    // Set cut values for x-axis
    cut.cutVals.resize(nVars);
    if (nVars > 1)
    {
        for (int iVar = 0; iVar < nVars; iVar++)
        {
            cut.cutVals[iVar] = lowVal + iVar * ((highVal - lowVal) / (nVars - 1));
        }
    }

    // Initialise and resize vectors to 0
    cut.isValDefAt.resize(kNumPart, std::vector<std::vector<Int_t>>(kNumPtInterval, std::vector<Int_t>(nVars, 0)));
    cut.sigVal.resize(kNumPart, std::vector<std::vector<Double_t>>(nVars, std::vector<Double_t>(kNumPtInterval, 0.0)));
    cut.sigErr.resize(kNumPart, std::vector<std::vector<Double_t>>(nVars, std::vector<Double_t>(kNumPtInterval, 0.0)));
    cut.bgVal.resize(kNumPart, std::vector<std::vector<Double_t>>(nVars, std::vector<Double_t>(kNumPtInterval, 0.0)));
    cut.bgErr.resize(kNumPart, std::vector<std::vector<Double_t>>(nVars, std::vector<Double_t>(kNumPtInterval, 0.0)));
    cut.significance.resize(kNumPart, std::vector<std::vector<Double_t>>(kNumPtInterval, std::vector<Double_t>(nVars, 0.0)));
    cut.significanceErr.resize(kNumPart, std::vector<std::vector<Double_t>>(kNumPtInterval, std::vector<Double_t>(nVars, 0.0)));
    cut.loosestVar.resize(kNumPart, std::vector<Int_t>(kNumPtInterval, 0));
    cut.maxSigVal.resize(kNumPart, std::vector<Double_t>(kNumPtInterval, 0.0));
    cut.maxSigErr.resize(kNumPart, std::vector<Double_t>(kNumPtInterval, 0.0));
    cut.sigLoss.resize(kNumPart, std::vector<std::vector<Double_t>>(kNumPtInterval, std::vector<Double_t>(nVars, 100.0)));
    cut.sigLossErr.resize(kNumPart, std::vector<std::vector<Double_t>>(kNumPtInterval, std::vector<Double_t>(nVars, 0.0)));
}

void SetDefCutValue(cutInfo &cut, double cutValDefault, int particle, int ptInterval)
{
    int defFoundAt = -1;

    auto found = find(cut.cutVals.begin(), cut.cutVals.end(), cutValDefault);

    // if default value is in the cut values, set isValDefAt to true
    if (found != cut.cutVals.end())
    {
        defFoundAt = distance(cut.cutVals.begin(), found);
    }

    if (particle == -1)
    {
        for (int iPart = kXi; iPart < kNumPart; iPart++)
        {
            if (ptInterval == -1)
            {
                for (int iPt = kLow; iPt < kNumPtInterval; iPt++)
                {
                    cut.defVal[iPart][iPt] = cutValDefault;
                    if (defFoundAt > -1)
                        cut.isValDefAt[iPart][iPt][defFoundAt] = defFoundAt;
                }
            }
            else
            {
                cut.defVal[iPart][ptInterval] = cutValDefault;
                if (defFoundAt > -1)
                    cut.isValDefAt[iPart][ptInterval][defFoundAt] = defFoundAt;
            }
        }
    }
    else if (ptInterval == -1)
    {
        for (int iPt = kLow; iPt < kNumPtInterval; iPt++)
        {
            cut.defVal[particle][iPt] = cutValDefault;
            if (defFoundAt > -1)
                cut.isValDefAt[particle][iPt][defFoundAt] = defFoundAt;
        }
    }
    else
    {
        cut.defVal[particle][ptInterval] = cutValDefault;
        if (defFoundAt > -1)
            cut.isValDefAt[particle][ptInterval][defFoundAt] = defFoundAt;
    }
}

void GetSignalBg(cutInfo &cut, TH1 &h, Int_t varNum, Int_t particle, Int_t ptBin, Double_t sig[], Double_t bg[], Double_t errSigSq[], Double_t errBgSq[])
{
    /// low pt
    if (ptBin < fPtIntervalBins[particle][kMid])
    {
        sig[kLow] += h.GetBinContent(1);
        bg[kLow] += h.GetBinContent(2);
        errSigSq[kLow] += pow(h.GetBinError(1), 2);
        errBgSq[kLow] += pow(h.GetBinError(2), 2);
    }
    /// mid pt
    else if (ptBin < fPtIntervalBins[particle][kHigh])
    {
        sig[kMid] += h.GetBinContent(1);
        bg[kMid] += h.GetBinContent(2);
        errSigSq[kMid] += pow(h.GetBinError(1), 2);
        errBgSq[kMid] += pow(h.GetBinError(2), 2);
    }
    /// high pt
    else
    {
        sig[kHigh] += h.GetBinContent(1);
        bg[kHigh] += h.GetBinContent(2);
        errSigSq[kHigh] += pow(h.GetBinError(1), 2);
        errBgSq[kHigh] += pow(h.GetBinError(2), 2);
    }

    for (Int_t iPt = 0; iPt < kNumPtInterval; iPt++)
    {
        if (cut.isValDefAt[particle][iPt][varNum]) // skip default values
        {
            Warning("GetSignalBg", "%s_%d: Part = %s, PtInterval = %d : Default value (%.3f). Skipping ...", cut.cutName.Data(), varNum, (particle == kXi) ? "Xi" : "Om", iPt, cut.cutVals[varNum]);
            continue;
        }
        cut.sigVal[particle][varNum][iPt] = sig[iPt] / (fPtInterval[particle][iPt + 1] - fPtInterval[particle][iPt]);
        cut.bgVal[particle][varNum][iPt] = bg[iPt] / (fPtInterval[particle][iPt + 1] - fPtInterval[particle][iPt]);
        cut.sigErr[particle][varNum][iPt] = sqrt(errSigSq[iPt]) / (fPtInterval[particle][iPt + 1] - fPtInterval[particle][iPt]);
        cut.bgErr[particle][varNum][iPt] = sqrt(errBgSq[iPt]) / (fPtInterval[particle][iPt + 1] - fPtInterval[particle][iPt]);
    }
    // cut.sigVal[particle][varNum][kLow] = sig[kLow] / (fPtInterval[particle][kMid] - fPtInterval[particle][kLow]);
    // cut.bgVal[particle][varNum][kLow] = bg[kLow] / (fPtInterval[particle][kMid] - fPtInterval[particle][kLow]);
    // cut.sigVal[particle][varNum][kMid] = sig[kMid] / (fPtInterval[particle][kHigh] - fPtInterval[particle][kMid]);
    // cut.bgVal[particle][varNum][kMid] = bg[kMid] / (fPtInterval[particle][kHigh] - fPtInterval[particle][kMid]);
    // cut.sigVal[particle][varNum][kHigh] = sig[kHigh] / (fPtInterval[particle][kNumPtInterval] - fPtInterval[particle][kHigh]);
    // cut.bgVal[particle][varNum][kHigh] = bg[kHigh] / (fPtInterval[particle][kNumPtInterval] - fPtInterval[particle][kHigh]);

    // cut.sigErr[particle][varNum][kLow] = sqrt(errSigSq[kLow]) / (fPtInterval[particle][kMid] - fPtInterval[particle][kLow]);
    // cut.bgErr[particle][varNum][kLow] = sqrt(errBgSq[kLow]) / (fPtInterval[particle][kMid] - fPtInterval[particle][kLow]);
    // cut.sigErr[particle][varNum][kMid] = sqrt(errSigSq[kMid]) / (fPtInterval[particle][kHigh] - fPtInterval[particle][kMid]);
    // cut.bgErr[particle][varNum][kMid] = sqrt(errBgSq[kMid]) / (fPtInterval[particle][kHigh] - fPtInterval[particle][kMid]);
    // cut.sigErr[particle][varNum][kHigh] = sqrt(errSigSq[kHigh]) / (fPtInterval[particle][kNumPtInterval] - fPtInterval[particle][kHigh]);
    // cut.bgErr[particle][varNum][kHigh] = sqrt(errBgSq[kHigh]) / (fPtInterval[particle][kNumPtInterval] - fPtInterval[particle][kHigh]);

    Info("OptimiseCascades: GetSignalBg", "%s_%d : %s Low {%f+/-%f}, Mid {%f+/-%f}, High {%f+/-%f}", cut.cutName.Data(), varNum, h.GetName(), cut.sigVal[particle][varNum][kLow], cut.sigErr[particle][varNum][kLow], cut.sigVal[particle][varNum][kMid], cut.sigErr[particle][varNum][kMid], cut.sigVal[particle][varNum][kHigh], cut.sigErr[particle][varNum][kHigh]);
}

void CalculateSignificance(cutInfo &cut, cutInfo &defCut)
{
    for (Int_t iPart = kXi; iPart < kNumPart; iPart++)
    {
        for (Int_t iPt = 0; iPt < kNumPtInterval; iPt++)
        {
            Double_t maxSigVal = defCut.sigVal[iPart][0][iPt]; // default value
            Double_t maxSigErr = defCut.sigErr[iPart][0][iPt];
            Int_t varNumMaxSigVal = 0; // index of the loosest variation (max sigVal) for sig loss calculation

            for (Int_t iVar = 0; iVar < cut.nVars; iVar++)
            {
                Double_t N_signal = cut.sigVal[iPart][iVar][iPt];
                Double_t N_background = cut.bgVal[iPart][iVar][iPt];
                Double_t err_signal = cut.sigErr[iPart][iVar][iPt];
                Double_t err_background = cut.bgErr[iPart][iVar][iPt];

                if (N_signal < 1e-6) // if (N_signal + N_background) == 0)
                {
                    Warning("CalculateSignificance", "%s_%d: Part = %d, PtInterval = %d : Signal = %f+/-%f, Bg = %f+/-%f. Skipping ...", cut.cutName.Data(), iVar, iPart, iPt, N_signal, err_signal, N_background, err_background);
                    continue;
                }

                if (cut.isValDefAt[iPart][iPt][iVar]) /// should not happen
                {
                    Warning("CalculateSignificance", "%s_%d: Part = %d, PtInterval = %d : Default value (%.3f). Skipping ...", cut.cutName.Data(), iVar, iPart, iPt, cut.cutVals[iVar]);
                    continue;
                }

                /// Find the loosest variation (max sigval) for sig loss calculation
                if (N_signal > maxSigVal)
                {
                    maxSigVal = N_signal;
                    maxSigErr = err_signal;
                    varNumMaxSigVal = iVar;
                }

                // Step 1: Calculate the significance: N_signal / sqrt(N_signal + N_background)
                Double_t significance = N_signal / sqrt(N_signal + N_background);

                // Step 2: Calculate the error on significance using error propagation: (1 / sqrt(N_signal + N_background)) * sqrt(err_signal^2 + (N_signal / (N_signal + N_background))^2 * err_background^2)
                Double_t denom = sqrt(N_signal + N_background);
                Double_t signal_fraction = N_signal / (N_signal + N_background);

                Double_t sigma_significance = (1.0 / denom) * sqrt(pow(err_signal, 2) +
                                                                   pow(signal_fraction, 2) * pow(err_background, 2));

                // Step 3: Store the significance and error on significance
                cut.significance[iPart][iPt][iVar] = significance;
                cut.significanceErr[iPart][iPt][iVar] = sigma_significance;
            }

            /// Store the loosest variation (max sigval) for sig loss calculation
            cut.loosestVar[iPart][iPt] = varNumMaxSigVal;
            cut.maxSigVal[iPart][iPt] = maxSigVal;
            cut.maxSigErr[iPart][iPt] = maxSigErr;
        }
    }
}

void CalculateSigLoss(cutInfo &cut)
{
    // if (cut.nVars == 1)
    // {
    //     Warning("CalculateSigLoss", "%s: Only one variation. Skipping ...", cut.cutName.Data());
    //     return;
    // }

    for (Int_t iPart = kXi; iPart < kNumPart; iPart++)
    {
        for (Int_t iPt = 0; iPt < kNumPtInterval; iPt++)
        {
            Double_t maxSigVal = cut.maxSigVal[iPart][iPt];
            Double_t maxSigErr = cut.maxSigErr[iPart][iPt];

            for (Int_t iVar = 0; iVar < cut.nVars; iVar++)
            {
                Double_t N_signal = cut.sigVal[iPart][iVar][iPt];
                Double_t err_signal = cut.sigErr[iPart][iVar][iPt];

                if (N_signal < 1e-6) // if (N_signal + N_background) == 0)
                {
                    Warning("CalculateSigLoss", "%s_%d: Part = %d, PtInterval = %d : Signal = %f+/-%f. Skipping ...", cut.cutName.Data(), iVar, iPart, iPt, N_signal, err_signal);
                    continue;
                }

                if (cut.isValDefAt[iPart][iPt][iVar]) /// should not happen
                {
                    Warning("CalculateSigLoss", "%s_%d: Part = %d, PtInterval = %d : Default value (%.3f). Skipping ...", cut.cutName.Data(), iVar, iPart, iPt, cut.cutVals[iVar]);
                    continue;
                }

                /// Calculate the signal loss: (maxSigVal - N_signal) / maxSigVal
                Double_t sigLoss = (maxSigVal - N_signal) / maxSigVal;
                Double_t sigma_sigLoss = ErrorInRatio(N_signal, err_signal, maxSigVal, maxSigErr);

                /// calculate the percentage error
                if (sigLoss == 0)
                    Info("CalculateSigLoss", "sigLoss = %f+/-%f for %s_%d, Part = %d, PtInterval = %d", sigLoss, sigma_sigLoss, cut.cutName.Data(), iVar, iPart, iPt);
                else
                    sigma_sigLoss = (sigma_sigLoss / sigLoss) * 100;

                /// Convert to percentage
                sigLoss *= 100;
                /// Calculate the error on signal loss using error propagation: sqrt((1 / maxSigVal)^2 * err_signal^2 + (N_signal / maxSigVal)^2 * err_maxSig^2)
                // Double_t sigma_sigLoss = sqrt(pow(1.0 / maxSigVal, 2) * pow(err_signal, 2) +            pow(N_signal / maxSigVal, 2) * pow(maxSigErr, 2));

                /// Store the signal loss and error on signal loss
                cut.sigLoss[iPart][iPt][iVar] = sigLoss;
                cut.sigLossErr[iPart][iPt][iVar] = sigma_sigLoss;
            }
        }
    }

    if (cut.cutName == "DcaMesV0ToPv")
    {
        for (Int_t iPart = kXi; iPart < kNumPart; iPart++)
        {
            for (Int_t iPt = 0; iPt < kNumPtInterval; iPt++)
            {
                for (Int_t iVar = 0; iVar < cut.nVars; iVar++)
                {
                    Info("CalculateSigLoss", "%s_%d: Part = %d, PtInterval = %d : Signal Loss = %.3f +/- %.3f", cut.cutName.Data(), iVar, iPart, iPt, cut.sigLoss[iPart][iPt][iVar], cut.sigLossErr[iPart][iPt][iVar]);
                }
            }
        }
    }
}

void FindOptimumDefVal(cutInfo &cut, cutInfo &defCut)
{
    for (Int_t iPart = kXi; iPart < kNumPart; iPart++)
    {
        for (Int_t iPt = 0; iPt < kNumPtInterval; iPt++)
        {
            Double_t maxSignificance = 0.0;

            if (defCut.sigLoss[iPart][iPt][0] < 2)
            {
                maxSignificance = defCut.significance[iPart][iPt][0]; // initialise with default value
                cut.optVal[iPart][iPt] = cut.defVal[iPart][iPt];
            }
            else
                Warning("FindOptimumDefVal", "%s: Default value (%.4f) has sigLoss > 2 (%.1f%%), Part = %s, PtInterval = %d", cut.cutName.Data(), cut.defVal[iPart][iPt], defCut.sigLoss[iPart][iPt][0], iPart == kXi ? "Xi" : "Om", iPt);

            for (Int_t iVar = 0; iVar < cut.nVars; iVar++)
            {
                if (cut.sigLoss[iPart][iPt][iVar] < 2)
                {

                    Double_t significance = cut.significance[iPart][iPt][iVar];

                    if (significance > maxSignificance)
                    {
                        maxSignificance = significance;
                        cut.optVal[iPart][iPt] = cut.cutVals[iVar];
                    }
                }
            }
        }
    }
}

void PrintCutInfo(cutInfo &cut)
{
    Info("PrintCutInfo", "Cut: %d, isTopo: %d, cutName: %s, nVars: %d", cut.cutId, cut.isTopo, cut.cutName.Data(), cut.nVars);

    for (int iPart = kXi; iPart < kNumPart; iPart++)
    {
        for (int iVar = 0; iVar < cut.nVars; iVar++)
        {
            for (int iPt = 0; iPt < kNumPtInterval; iPt++)
            {
                Info("PrintCutInfo", "%s%s_%s_%d, cutVal %d: %.4f, Signal[%d][%d][%d]: %.2f +/- %.2f, Bg: %.2f +/- %.2f, Significance = %.2f +/- %.2f, SigLoss = %.3f +/- %.3f", (cut.isValDefAt[iPart][iPt][iVar]) ? "DEF_" : "", cut.cutName.Data(), iPart == kXi ? "Xi" : "Om", iPt, iVar, cut.cutVals[iVar], iPart, iVar, iPt, cut.sigVal[iPart][iVar][iPt], cut.sigErr[iPart][iVar][iPt], cut.bgVal[iPart][iVar][iPt], cut.bgErr[iPart][iVar][iPt], cut.significance[iPart][iPt][iVar], cut.significanceErr[iPart][iPt][iVar], cut.sigLoss[iPart][iPt][iVar], cut.sigLossErr[iPart][iPt][iVar]);
            }
        }
    }
}

void DrawAndSave(TGraphErrors *graph1, TGraphErrors *graph2, Double_t xDefLine, Double_t xOptLine, TString outputFolder, TString imageFolder, TString imageName, TString imageFormat)
{

    gROOT->SetBatch(kTRUE);
    gROOT->SetStyle("Plain");
    gStyle->SetOptStat(0);

    graph1->SetMarkerStyle(kFullCircle);
    graph1->SetMarkerSize(1.5);
    graph2->SetMarkerStyle(kStar);
    graph2->SetMarkerSize(1.5);
    graph2->GetYaxis()->SetRangeUser(0, graph2->GetYaxis()->GetXmax());

    TLine *defLine1 = new TLine(xDefLine, graph1->GetYaxis()->GetXmin(), xDefLine, graph1->GetYaxis()->GetXmax());
    defLine1->SetLineColor(kRed);
    defLine1->SetLineWidth(2);
    defLine1->SetLineStyle(2);
    TLine *defLine2 = new TLine(xDefLine, graph2->GetYaxis()->GetXmin(), xDefLine, graph2->GetYaxis()->GetXmax());
    defLine2->SetLineColor(kRed);
    defLine2->SetLineWidth(2);
    defLine2->SetLineStyle(2);

    TLine *sigLoss2 = new TLine(graph2->GetXaxis()->GetXmin(), 2, graph2->GetXaxis()->GetXmax(), 2);
    sigLoss2->SetLineColor(kGreen);
    sigLoss2->SetLineWidth(3);
    sigLoss2->SetLineStyle(6);
    TLine *sigLoss5 = new TLine(graph2->GetXaxis()->GetXmin(), 5, graph2->GetXaxis()->GetXmax(), 5);
    sigLoss5->SetLineColor(kMagenta);
    sigLoss5->SetLineWidth(2);
    sigLoss5->SetLineStyle(9);
    TLine *sigLoss10 = new TLine(graph2->GetXaxis()->GetXmin(), 10, graph2->GetXaxis()->GetXmax(), 10);
    sigLoss10->SetLineColor(kMagenta);
    sigLoss10->SetLineWidth(2);
    sigLoss10->SetLineStyle(10);

    TLine *optLine1 = new TLine(xOptLine, graph1->GetYaxis()->GetXmin(), xOptLine, graph1->GetYaxis()->GetXmax());
    optLine1->SetLineColor(kBlue);
    optLine1->SetLineWidth(3);
    optLine1->SetLineStyle(6);
    TLine *optLine2 = new TLine(xOptLine, graph2->GetYaxis()->GetXmin(), xOptLine, graph2->GetYaxis()->GetXmax());
    optLine2->SetLineColor(kBlue);
    optLine2->SetLineWidth(3);
    optLine2->SetLineStyle(6);

    TCanvas *cDraw = new TCanvas("cDraw_multipads", "Significance and Signal Loss", 1920, 2160);
    gStyle->SetOptTitle(0);
    // cDraw->Divide(1, 2, 0.0, 0.0);
    TPad *p1 = new TPad("p1", "p1", 0., 0., 1., 0.49, 0, 0, 0);
    p1->SetTopMargin(0);
    p1->Draw();

    TPad *p2 = new TPad("p2", "p2", 0., 0.49, 1., 0.98, 0, 0, 0);
    p2->SetTopMargin(0);
    p2->SetBottomMargin(0);
    p2->Draw();

    p1->cd();
    gPad->SetTicks(2, 2);
    gPad->SetGrid();
    graph2->Draw("APLX");
    defLine2->Draw("same");
    sigLoss2->Draw("same");
    sigLoss5->Draw("same");
    sigLoss10->Draw("same");
    optLine2->Draw("same");

    p2->cd();
    gPad->SetTicks(2, 2);
    graph1->GetXaxis()->SetLabelSize(0);
    graph1->Draw("APL");
    defLine1->Draw("same");
    optLine1->Draw("same");

    // cDraw->BuildLegend();

    auto legend = new TLegend(0.74, 0.80, 0.9, 0.98);
    legend->SetHeader(graph1->GetTitle(), "C"); // option "C" allows to center the header
    legend->AddEntry(graph1, "Signal Significance", "lp");
    legend->AddEntry(defLine1, "Chosen Default", "l");
    legend->AddEntry(optLine1, "Optimum Default", "l");
    legend->AddEntry(graph2, "Signal Loss", "lp");
    legend->Draw();

    SaveImage(outputFolder, imageFolder, imageName, imageFormat, cDraw);

    p1->Delete();
    p2->Delete();
    defLine1->Delete();
    defLine2->Delete();
    sigLoss2->Delete();
    sigLoss5->Delete();
    sigLoss10->Delete();
    optLine1->Delete();
    optLine2->Delete();
    legend->Delete();
    cDraw->Close();
    delete cDraw;
}