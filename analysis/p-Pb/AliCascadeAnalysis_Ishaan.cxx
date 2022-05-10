class TTree;
class TParticle;
class TVector3;

class AliESDVertex;
class AliESDv0;
/// class AliESDcascade;
class AliAODVertex;
class AliAODv0;
class AliAODcascade;

#include <Riostream.h>
#include "TH3.h"
#include "TFile.h"
#include "TCanvas.h"
#include "THistManager.h"
#include "AliESDEvent.h"
#include "AliESDtrack.h"
#include "AliAODEvent.h"
#include "AliAODTrack.h"
#include "AliAODVertex.h"
#include "AliPID.h"
#include "AliInputEventHandler.h"
#include "AliAnalysisManager.h"
#include "AliMultSelection.h"
#include "AliESDcascade.h"
#include "AliAODMCParticle.h"
#include "AliMCParticle.h"
#include "AliMCEvent.h"
#include "AliAnalysisTaskESDfilter.h"
#include "AliAnalysisUtils.h"
#include "AliAODMCHeader.h"

#include "AliCascadeAnalysis_Ishaan.h"

#include "TMath.h"

ClassImp(AliCascadeAnalysis_Ishaan)

    AliCascadeAnalysis_Ishaan::AliCascadeAnalysis_Ishaan() : AliAnalysisTaskSE(),
                                                             /// outputs
                                                             fHistos_eve(nullptr),
                                                             //    fHistos_K0S(nullptr),
                                                             //  fHistos_Lam(nullptr),
                                                             //  fHistos_ALam(nullptr),
                                                             fHistos_XiMin(nullptr),
                                                             fHistos_XiPlu(nullptr),
                                                             fHistos_OmMin(nullptr),
                                                             fHistos_OmPlu(nullptr),
                                                             /// objects from the manager
                                                             fPIDResponse(0),
                                                             fTriggerMask(0),

                                                             // MC-related variables
                                                             //   fisMC(kFALSE),
                                                             //   fisMCassoc(kFALSE),

                                                             /// default cuts configuration
                                                             fDefOnly(kTRUE),
                                                             //  fV0_Cuts{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
                                                             fCasc_Cuts{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
                                                             /// particle to be analysed
                                                             fParticleAnalysisStatus{true, true, true, true, true, true, true},
                                                             /// variables for V0 cuts
                                                             //  fV0_DcaV0Daught(0),
                                                             //  fV0_DcaPosToPV(0),
                                                             //  fV0_DcaNegToPV(0),
                                                             //  fV0_V0CosPA(0),
                                                             //  fV0_V0Rad(0),
                                                             //  fV0_Pt(0),
                                                             //  fV0_yK0S(0),
                                                             //  fV0_yLam(0),
                                                             //  fV0_etaPos(0),
                                                             //  fV0_etaNeg(0),
                                                             //  fV0_InvMassK0s(0),
                                                             //  fV0_InvMassLam(0),
                                                             //  fV0_InvMassALam(0),
                                                             //  fV0_LeastCRaws(0),
                                                             //  fV0_LeastCRawsOvF(0),
                                                             //  fV0_LeastTPCcls(0),
                                                             //  fV0_NSigPosProton(0),
                                                             //  fV0_NSigPosPion(0),
                                                             //  fV0_NSigNegProton(0),
                                                             //  fV0_NSigNegPion(0),
                                                             //  fV0_DistOverTotP(0),
                                                             //  fV0_ITSTOFtracks(0),
                                                             //  fV0_NegTrackStatus(0),
                                                             //  fV0_PosTrackStatus(0),
                                                             //  fV0_kinkidx(0),

                                                             /// variables for Cascade analysis
                                                             fCasc_DcaCascDaught(0),
                                                             fCasc_CascCosPA(0),
                                                             fCasc_CascRad(0),
                                                             fCasc_etaPos(0),
                                                             fCasc_etaNeg(0),
                                                             fCasc_etaBac(0),
                                                             fCasc_kinkidx(0),
                                                             fCasc_NSigPosProton(0),
                                                             fCasc_NSigPosPion(0),
                                                             fCasc_NSigNegProton(0),
                                                             fCasc_NSigNegPion(0),
                                                             fCasc_NSigBacPion(0),
                                                             fCasc_NSigBacKaon(0),
                                                             fCasc_LeastCRaws(0),
                                                             fCasc_LeastCRawsOvF(0),
                                                             fCasc_LeastTPCcls(0),
                                                             fCasc_InvMassLam(0),
                                                             fCasc_DcaV0Daught(0),
                                                             fCasc_V0CosPA(0),
                                                             fCasc_DcaV0ToPV(0),
                                                             fCasc_DcaBachToPV(0),
                                                             fCasc_ITSTOFtracks(0),
                                                             fCasc_yXi(0),
                                                             fCasc_yOm(0),
                                                             fCasc_charge(0),
                                                             fCasc_Pt(0),
                                                             fCasc_DistOverTotP(0),
                                                             fCasc_InvMassXiMin(0),
                                                             fCasc_InvMassXiPlu(0),
                                                             fCasc_InvMassOmMin(0),
                                                             fCasc_InvMassOmPlu(0),
                                                             fCasc_V0Rad(0),
                                                             fCasc_DcaPosToPV(0),
                                                             fCasc_DcaNegToPV(0),
                                                             fCasc_NegTrackStatus(0),
                                                             fCasc_PosTrackStatus(0),
                                                             fCasc_BacTrackStatus(0),
                                                             fCasc_BacBarCosPA(0),
                                                             fisParametricBacBarCosPA(kFALSE),
                                                             fHist_PtBacBarCosPA(0),
                                                             ptXiBoundary_LowMid(0),
                                                             ptXiBoundary_MidHigh(0),
                                                             ptOmBoundary_LowMid(0),
                                                             ptOmBoundary_MidHigh(0),
                                                             fCentLimit_BacBarCosPA(0)

{
    ////default constructor
}

AliCascadeAnalysis_Ishaan::AliCascadeAnalysis_Ishaan(const char *name, TString lExtraOptions) : AliAnalysisTaskSE(name),
                                                                                                ////outputs
                                                                                                fHistos_eve(nullptr),
                                                                                                // fHistos_K0S(nullptr),
                                                                                                // fHistos_Lam(nullptr),
                                                                                                // fHistos_ALam(nullptr),
                                                                                                fHistos_XiMin(nullptr),
                                                                                                fHistos_XiPlu(nullptr),
                                                                                                fHistos_OmMin(nullptr),
                                                                                                fHistos_OmPlu(nullptr),
                                                                                                ////objects from the manager
                                                                                                fPIDResponse(0),
                                                                                                fTriggerMask(0),

                                                                                                ////MC-related variables
                                                                                                // fisMC(kFALSE),
                                                                                                // fisMCassoc(kFALSE),

                                                                                                ////default cuts configuration
                                                                                                fDefOnly(kTRUE),
                                                                                                // fV0_Cuts{1., 0.11, 0.11, 0.97, 1., 0.5, 0.8, 70., 0.8, 50., 5., 20., 30., 1.},
                                                                                                fCasc_Cuts{1.5, 0.96, 0.6, 4., 70., 0.8, 70., 0.008, 1.6, 0.98, 0.06, 0.04, 1., -0.5, 0.8, 3., 3., 1.2, 0.04, 0.03, 1., 0.5, 1.1, 1.6, 1.4, 0.97, 0.97, 1.7, 1.5, 0.97, 0.98, 0.98},
                                                                                                ////particle to be analysed
                                                                                                fParticleAnalysisStatus{true, true, true, true, true, true, true},
                                                                                                ////variables for V0 cuts
                                                                                                // fV0_DcaV0Daught(0),
                                                                                                // fV0_DcaPosToPV(0),
                                                                                                // fV0_DcaNegToPV(0),
                                                                                                // fV0_V0CosPA(0),
                                                                                                // fV0_V0Rad(0),
                                                                                                // fV0_Pt(0),
                                                                                                // fV0_yK0S(0),
                                                                                                // fV0_yLam(0),
                                                                                                // fV0_etaPos(0),
                                                                                                // fV0_etaNeg(0),
                                                                                                // fV0_InvMassK0s(0),
                                                                                                // fV0_InvMassLam(0),
                                                                                                // fV0_InvMassALam(0),
                                                                                                // fV0_LeastCRaws(0),
                                                                                                // fV0_LeastCRawsOvF(0),
                                                                                                // fV0_LeastTPCcls(0),
                                                                                                // fV0_NSigPosProton(0),
                                                                                                // fV0_NSigPosPion(0),
                                                                                                // fV0_NSigNegProton(0),
                                                                                                // fV0_NSigNegPion(0),
                                                                                                // fV0_DistOverTotP(0),
                                                                                                // fV0_ITSTOFtracks(0),
                                                                                                // fV0_NegTrackStatus(0),
                                                                                                // fV0_PosTrackStatus(0),
                                                                                                // fV0_kinkidx(0),
                                                                                                ////variables for Cascade analysis
                                                                                                fCasc_DcaCascDaught(0),
                                                                                                fCasc_CascCosPA(0),
                                                                                                fCasc_CascRad(0),
                                                                                                fCasc_etaPos(0),
                                                                                                fCasc_etaNeg(0),
                                                                                                fCasc_etaBac(0),
                                                                                                fCasc_kinkidx(0),
                                                                                                fCasc_NSigPosProton(0),
                                                                                                fCasc_NSigPosPion(0),
                                                                                                fCasc_NSigNegProton(0),
                                                                                                fCasc_NSigNegPion(0),
                                                                                                fCasc_NSigBacPion(0),
                                                                                                fCasc_NSigBacKaon(0),
                                                                                                fCasc_LeastCRaws(0),
                                                                                                fCasc_LeastCRawsOvF(0),
                                                                                                fCasc_LeastTPCcls(0),
                                                                                                fCasc_InvMassLam(0),
                                                                                                fCasc_DcaV0Daught(0),
                                                                                                fCasc_V0CosPA(0),
                                                                                                fCasc_DcaV0ToPV(0),
                                                                                                fCasc_DcaBachToPV(0),
                                                                                                fCasc_ITSTOFtracks(0),
                                                                                                fCasc_yXi(0),
                                                                                                fCasc_yOm(0),
                                                                                                fCasc_charge(0),
                                                                                                fCasc_Pt(0),
                                                                                                fCasc_DistOverTotP(0),
                                                                                                fCasc_InvMassXiMin(0),
                                                                                                fCasc_InvMassXiPlu(0),
                                                                                                fCasc_InvMassOmMin(0),
                                                                                                fCasc_InvMassOmPlu(0),
                                                                                                fCasc_V0Rad(0),
                                                                                                fCasc_DcaPosToPV(0),
                                                                                                fCasc_DcaNegToPV(0),
                                                                                                fCasc_NegTrackStatus(0),
                                                                                                fCasc_PosTrackStatus(0),
                                                                                                fCasc_BacTrackStatus(0),
                                                                                                fCasc_BacBarCosPA(0),
                                                                                                fisParametricBacBarCosPA(kFALSE),
                                                                                                fHist_PtBacBarCosPA(0),
                                                                                                ptXiBoundary_LowMid(0),
                                                                                                ptXiBoundary_MidHigh(0),
                                                                                                ptOmBoundary_LowMid(0),
                                                                                                ptOmBoundary_MidHigh(0),
                                                                                                fCentLimit_BacBarCosPA(0)

{
    ////setting default cuts
    SetDefCutVals();
    // SetDefCutVariations();
    ////setting default centrality binning
    Double_t centbins_Xi[11] = {0, 5, 10, 15, 20, 30, 40, 50, 60, 80, 100}; // V0A
    Double_t centbins_Om[6] = {0, 5, 15, 30, 60, 100};                      // V0A
    Int_t ncentbins_Xi = sizeof(centbins_Xi) / sizeof(Double_t) - 1;
    Int_t ncentbins_Om = sizeof(centbins_Om) / sizeof(Double_t) - 1;

    // for (int ipart = 0; ipart < knumpart; ipart++)
    // {
    //     SetCentbinning(ipart, 11, centbins);
    // }
    SetCentbinning(kXi, ncentbins_Xi, centbins_Xi);
    SetCentbinning(kOm, ncentbins_Om, centbins_Om);

    ////setting default mass binning
    int massbins[4] = {90, 50, 100, 100};
    double minmass[4] = {0.41, 1.07, 1.272, 1.622};
    double maxmass[4] = {0.59, 1.17, 1.372, 1.722};
    for (int ipart = 0; ipart < knumpart; ipart++)
    {
        SetMassbinning(ipart, massbins[ipart], minmass[ipart], maxmass[ipart]);
    }
    ////setting default pt binning
    // double ptbins[4][251];
    // int nptbins[4] = {250, 250, 75, 75};

    // double maxpt[4] = {25., 25., 15., 15.};
    // for (int ipart = 0; ipart < knumpart; ipart++)
    // {
    //     for (int ipt = 0; ipt < nptbins[ipart] + 1; ipt++)
    //     {
    //         ptbins[ipart][ipt] = ipt * maxpt[ipart] / nptbins[ipart];
    //     }
    //     SetPtbinning(ipart, nptbins[ipart], ptbins[ipart]);
    // }

    Double_t ptbins_Xi[] = {0.8, 1.1, 1.3, 1.5, 1.7, 1.9, 2.1, 2.3, 2.5, 2.7, 2.9, 3.1, 3.5, 4, 5.5};
    Double_t ptbins_Om[] = {0.9, 1.6, 2., 2.4, 2.9, 3.5, 5};
    Int_t nptbins_Xi = sizeof(ptbins_Xi) / sizeof(Double_t) - 1;
    Int_t nptbins_Om = sizeof(ptbins_Om) / sizeof(Double_t) - 1;

    SetPtbinning(kXi, nptbins_Xi, ptbins_Xi);
    SetPtbinning(kOm, nptbins_Om, ptbins_Om);

    ////Standard output
    DefineOutput(1, TList::Class()); //// Event Histograms
    // DefineOutput(2, TList::Class()); //// K0S Histograms
    // DefineOutput(3, TList::Class()); //// Lambda Histograms
    // DefineOutput(4, TList::Class()); //// AntiLambda Histograms
    DefineOutput(2, TList::Class()); //// XiMinus Histograms
    DefineOutput(3, TList::Class()); //// XiPlus Histograms
    DefineOutput(4, TList::Class()); //// OmegaMinus Histograms
    DefineOutput(5, TList::Class()); //// OmegaPlus Histograms
}

AliCascadeAnalysis_Ishaan::~AliCascadeAnalysis_Ishaan()
{
    ////------------------------------------------------
    //// DESTRUCTOR
    ////------------------------------------------------

    ////Destroy output objects if present
    ///     if (fListHist) {
    ///         delete fListHist;
    ///         fListHist = 0x0;
    ///     }
}

//_____________________________________________________________________________
void AliCascadeAnalysis_Ishaan::UserCreateOutputObjects()
{

    ////histograms for event variables
    fHistos_eve = new THistManager("histos_eve");

    fHistos_eve->CreateTH1("hcent", "Multiplicity Distribution", 100, 0, 100, "s"); ////storing #events in bins of centrality
    fHistos_eve->CreateTH1("henum", "", 1, 0, 1);          ////storing total #events
    fHistos_eve->CreateTH1("fCuts", "Accepted events after event selection", 5, 0, 5, "s"); ////storing impact of cuts on no. of events
    fHistos_eve->CreateTH1("fHistPtXiP", "fHistPtXiP", fnptbins[kXi], fptbinning[kXi], "s");
    fHistos_eve->CreateTH1("fHistPtXiM", "fHistPtXiM", fnptbins[kXi], fptbinning[kXi], "s");
    fHistos_eve->CreateTH1("fHistPtOmP", "fHistPtOmP", fnptbins[kOm], fptbinning[kOm], "s");
    fHistos_eve->CreateTH1("fHistPtOmM", "fHistPtOmM", fnptbins[kOm], fptbinning[kOm], "s");
    fHistos_eve->CreateTH1("fCasc_InvMassXiMin", "fCasc_InvMassXiMin", 100, 1.272, 1.372, "s");
    fHistos_eve->CreateTH1("fCasc_InvMassXiPlu", "fCasc_InvMassXiPlu", 100, 1.272, 1.372, "s");
    fHistos_eve->CreateTH1("fCasc_InvMassOmMin", "fCasc_InvMassOmMin", 100, 1.622, 1.722, "s");
    fHistos_eve->CreateTH1("fCasc_InvMassOmPlu", "fCasc_InvMassOmPlu", 100, 1.622, 1.722, "s");
    fHistos_eve->CreateTH1("fHistCentXiP", "fHistCentXiP", fncentbins[kXi], fcentbinning[kXi], "s");
    fHistos_eve->CreateTH1("fHistCentXiM", "fHistCentXiM", fncentbins[kXi], fcentbinning[kXi], "s");
    fHistos_eve->CreateTH1("fHistCentOmP", "fHistCentOmP", fncentbins[kOm], fcentbinning[kOm], "s");
    fHistos_eve->CreateTH1("fHistCentOmM", "fHistCentOmM", fncentbins[kOm], fcentbinning[kOm], "s");
    
    TH1 *fCuts = (TH1 *)fHistos_eve->FindObject("fCuts");
    fCuts->GetXaxis()->SetBinLabel(1, "Total Events");
    fCuts->GetXaxis()->SetBinLabel(2, "V0A_MultSel");
    fCuts->GetXaxis()->SetBinLabel(3, "InelGT0");
    fCuts->GetXaxis()->SetBinLabel(4, "|PVz|<10");
    fCuts->GetXaxis()->SetBinLabel(5, "Mult<0,100> && kINT7");

    // ////histograms for V0 variables
    // if (fParticleAnalysisStatus[kk0s])
    // {
    //     fHistos_K0S = new THistManager("histos_K0S");
    //     if (fisMC)
    //         fHistos_K0S->CreateTH2("h2_gen", "", fnptbins[kK0s], fptbinning[kK0s], fncentbins[kK0s], fcentbinning[kK0s]);
    //     fHistos_K0S->CreateTH3("h3_ptmasscent_def", "", fnptbins[kK0s], fptbinning[kK0s], fnmassbins[kK0s], fmassbinning[kK0s], fncentbins[kK0s], fcentbinning[kK0s]);
    // }
    // if (fParticleAnalysisStatus[klam])
    // {
    //     fHistos_Lam = new THistManager("histos_Lam");
    //     fHistos_ALam = new THistManager("histos_ALam");
    //     fHistos_Lam->CreateTH3("h3_ptmasscent_def", "", fnptbins[kLam], fptbinning[kLam], fnmassbins[kLam], fmassbinning[kLam], fncentbins[kLam], fcentbinning[kLam]);
    //     fHistos_ALam->CreateTH3("h3_ptmasscent_def", "", fnptbins[kLam], fptbinning[kLam], fnmassbins[kLam], fmassbinning[kLam], fncentbins[kLam], fcentbinning[kLam]);
    //     if (fisMC)
    //     {
    //         fHistos_Lam->CreateTH2("h2_gen", "", fnptbins[kLam], fptbinning[kLam], fncentbins[kLam], fcentbinning[kLam]);
    //         fHistos_Lam->CreateTH3("h3_FDmtxNUM_def", "", fnptbins[kLam], fptbinning[kLam], fnptbins[kXi], fptbinning[kXi], fncentbins[kLam], fcentbinning[kLam]);
    //         fHistos_Lam->CreateTH2("h2_FDmtxDEN_def", "", fnptbins[kXi], fptbinning[kXi], fncentbins[kLam], fcentbinning[kLam]);
    //         fHistos_ALam->CreateTH2("h2_gen", "", fnptbins[kLam], fptbinning[kLam], fncentbins[kLam], fcentbinning[kLam]);
    //         fHistos_ALam->CreateTH3("h3_FDmtxNUM_def", "", fnptbins[kLam], fptbinning[kLam], fnptbins[kXi], fptbinning[kXi], fncentbins[kLam], fcentbinning[kLam]);
    //         fHistos_ALam->CreateTH2("h2_FDmtxDEN_def", "", fnptbins[kXi], fptbinning[kXi], fncentbins[kLam], fcentbinning[kLam]);
    //     }
    // }
    // if (!fDefOnly && (fParticleAnalysisStatus[kk0s] || fParticleAnalysisStatus[klam]))
    // {
    //     for (int icut = 0; icut < kV0cutsnum; icut++)
    //     {
    //         if (icut == kV0_y || icut == kV0_etaDaugh)
    //             continue;
    //         for (int ivar = 0; ivar < nvarcut_V0[icut]; ivar++)
    //         {
    //             if (fParticleAnalysisStatus[kk0s])
    //             {
    //                 if (icut != kV0_PropLifetLam)
    //                     fHistos_K0S->CreateTH3(Form("h3_ptmasscent[%d][%d]", icut, ivar), "", fnptbins[kK0s], fptbinning[kK0s], fnmassbins[kK0s], fmassbinning[kK0s], fncentbins[kK0s], fcentbinning[kK0s]);
    //             }
    //             if (fParticleAnalysisStatus[klam])
    //             {
    //                 if (icut != kV0_PropLifetK0s)
    //                 {
    //                     fHistos_Lam->CreateTH3(Form("h3_ptmasscent[%d][%d]", icut, ivar), "", fnptbins[kLam], fptbinning[kLam], fnmassbins[kLam], fmassbinning[kLam], fncentbins[kLam], fcentbinning[kLam]);
    //                     fHistos_ALam->CreateTH3(Form("h3_ptmasscent[%d][%d]", icut, ivar), "", fnptbins[kLam], fptbinning[kLam], fnmassbins[kLam], fmassbinning[kLam], fncentbins[kLam], fcentbinning[kLam]);
    //                     if (fisMC)
    //                     {
    //                         fHistos_Lam->CreateTH3(Form("h3_FDmtxNUM[%d][%d]", icut, ivar), "", fnptbins[kLam], fptbinning[kLam], fnptbins[kXi], fptbinning[kXi], fncentbins[kLam], fcentbinning[kLam]);
    //                         fHistos_Lam->CreateTH2(Form("h2_FDmtxDEN[%d][%d]", icut, ivar), "", fnptbins[kXi], fptbinning[kXi], fncentbins[kLam], fcentbinning[kLam]);
    //                         fHistos_ALam->CreateTH3(Form("h3_FDmtxNUM[%d][%d]", icut, ivar), "", fnptbins[kLam], fptbinning[kLam], fnptbins[kXi], fptbinning[kXi], fncentbins[kLam], fcentbinning[kLam]);
    //                         fHistos_ALam->CreateTH2(Form("h2_FDmtxDEN[%d][%d]", icut, ivar), "", fnptbins[kXi], fptbinning[kXi], fncentbins[kLam], fcentbinning[kLam]);
    //                     }
    //                 }
    //             }
    //         }
    //     }
    // }
    ////histograms for Cascade variables
    if (fParticleAnalysisStatus[kxip])
    {
        fHistos_XiMin = new THistManager("histos_XiMin");
        fHistos_XiPlu = new THistManager("histos_XiPlu");
        fHistos_XiMin->CreateTH3("h3_ptmasscent_def", "", fnptbins[kXi], fptbinning[kXi], fnmassbins[kXi], fmassbinning[kXi], fncentbins[kXi], fcentbinning[kXi]);
        fHistos_XiPlu->CreateTH3("h3_ptmasscent_def", "", fnptbins[kXi], fptbinning[kXi], fnmassbins[kXi], fmassbinning[kXi], fncentbins[kXi], fcentbinning[kXi]);
        // if (fisMC)
        //     fHistos_XiMin->CreateTH2("h2_gen", "", fnptbins[kXi], fptbinning[kXi], fncentbins[kXi], fcentbinning[kXi]);
        // if (fisMC)
        //     fHistos_XiPlu->CreateTH2("h2_gen", "", fnptbins[kXi], fptbinning[kXi], fncentbins[kXi], fcentbinning[kXi]);
    }
    if (fParticleAnalysisStatus[komp])
    {
        fHistos_OmMin = new THistManager("histos_OmMin");
        fHistos_OmPlu = new THistManager("histos_OmPlu");
        fHistos_OmMin->CreateTH3("h3_ptmasscent_def", "", fnptbins[kOm], fptbinning[kOm], fnmassbins[kOm], fmassbinning[kOm], fncentbins[kOm], fcentbinning[kOm]);
        fHistos_OmPlu->CreateTH3("h3_ptmasscent_def", "", fnptbins[kOm], fptbinning[kOm], fnmassbins[kOm], fmassbinning[kOm], fncentbins[kOm], fcentbinning[kOm]);
        // if (fisMC)
        //     fHistos_OmMin->CreateTH2("h2_gen", "", fnptbins[kOm], fptbinning[kOm], fncentbins[kOm], fcentbinning[kOm]);
        // if (fisMC)
        //     fHistos_OmPlu->CreateTH2("h2_gen", "", fnptbins[kOm], fptbinning[kOm], fncentbins[kOm], fcentbinning[kOm]);
    }
    if (!fDefOnly && (fParticleAnalysisStatus[kxip] || fParticleAnalysisStatus[komp]))
    {
        for (int icut = 0; icut < kCasccutsnum; icut++)
        {
            if (icut == kCasc_y || icut == kCasc_etaDaugh)
                continue;
            for (int ivar = 0; ivar < nvarcut_Casc[icut]; ivar++)
            {
                if (fParticleAnalysisStatus[kxip])
                {
                    if (icut != kCasc_PropLifetOm)
                        fHistos_XiMin->CreateTH3(Form("h3_ptmasscent[%d][%d]", icut, ivar), "", fnptbins[kXi], fptbinning[kXi], fnmassbins[kXi], fmassbinning[kXi], fncentbins[kXi], fcentbinning[kXi]);
                    if (icut != kCasc_PropLifetOm)
                        fHistos_XiPlu->CreateTH3(Form("h3_ptmasscent[%d][%d]", icut, ivar), "", fnptbins[kXi], fptbinning[kXi], fnmassbins[kXi], fmassbinning[kXi], fncentbins[kXi], fcentbinning[kXi]);
                }
                if (fParticleAnalysisStatus[komp])
                {
                    if (icut != kCasc_PropLifetXi)
                        fHistos_OmMin->CreateTH3(Form("h3_ptmasscent[%d][%d]", icut, ivar), "", fnptbins[kOm], fptbinning[kOm], fnmassbins[kOm], fmassbinning[kOm], fncentbins[kOm], fcentbinning[kOm]);
                    if (icut != kCasc_PropLifetXi)
                        fHistos_OmPlu->CreateTH3(Form("h3_ptmasscent[%d][%d]", icut, ivar), "", fnptbins[kOm], fptbinning[kOm], fnmassbins[kOm], fmassbinning[kOm], fncentbins[kOm], fcentbinning[kOm]);
                }
            }
        }
    }
    //// PID Setup
    AliAnalysisManager *man = AliAnalysisManager::GetAnalysisManager();
    AliInputEventHandler *inputHandler = (AliInputEventHandler *)(man->GetInputEventHandler());
    fPIDResponse = inputHandler->GetPIDResponse();
    inputHandler->SetNeedField();

    ////Output posting
    DataPosting();
}
//_____________________________________________________________________________
void AliCascadeAnalysis_Ishaan::UserExec(Option_t *)
{

    AliAODEvent *lAODevent = dynamic_cast<AliAODEvent *>(InputEvent());
    if (!lAODevent)
    {
        AliWarning("ERROR: event not available \n");
        DataPosting();
        return;
    }

    ////dumb histo for checking
    fHistos_eve->FillTH1("henum", 0.5);
    fHistos_eve->FillTH1("fCuts", 0.5);

    ////get trigger information
    fTriggerMask = ((AliInputEventHandler *)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected();

    //// Multiplicity Information
    double lPercentile = -666;
    int lEvSelCode = -666;
    AliMultSelection *MultSelection = (AliMultSelection *)lAODevent->FindListObject("MultSelection");
    if (!MultSelection)
    {
        AliWarning("AliMultSelection object not found!");
        DataPosting();
        return;
    }
    else
    {
        lPercentile = MultSelection->GetMultiplicityPercentile("V0A"); /// changed V0M to V0A
        // fHistos_eve->FillTH1("htest", lPercentile);
        lEvSelCode = MultSelection->GetEvSelCode(); ////==0 means event is good. Set by AliMultSelectionTask
    }

    ////skip everything if event selection code !=0
    if (lEvSelCode != 0)
    {
        DataPosting();
        return;
    }

    fHistos_eve->FillTH1("fCuts", 1.5);

    /// Inel>0 in the |eta|<1 range
    const AliMultiplicity *mult = (AliMultiplicity *)lAODevent->GetMultiplicity();
    Bool_t inelgt0 = kFALSE;
    Int_t counterTracklets = 0;
    for (Int_t i = 0; i < mult->GetNumberOfTracklets(); ++i)
    {
        if (TMath::Abs(mult->GetEta(i)) < 1.)
        {
            inelgt0 = kTRUE;
            counterTracklets++;
        }
    }

    if (counterTracklets < 1)
    {
        DataPosting();
        return;
    }

    fHistos_eve->FillTH1("fCuts", 2.5);

    Double_t cutPrimVertex = 10.0;

    ////acquire best PV
    const AliVVertex *lBestPrimVtx = lAODevent->GetPrimaryVertex();

    if (!((TMath::Abs(lBestPrimVtx->GetZ())) < cutPrimVertex)) /// vertex_z<10cm
    {
        DataPosting();
        return;
    }

    fHistos_eve->FillTH1("fCuts", 3.5);

    double lBestPV[3] = {-666., -666., -666.};
    lBestPrimVtx->GetXYZ(lBestPV);
    float lfBestPV[3];
    for (int i = 0; i < 3; i++)
        lfBestPV[i] = (float)lBestPV[i];

    ////acquire magnetic field
    double lMagField = -666;
    lMagField = lAODevent->GetMagneticField();

    // bool assFlag[ksignednumpart] = {1, 1, 1, 1, 1, 1, 1}; ////MC ass flags are true if associated and, by construction, always for data. They can be false only if fisMC && notassociated

    /// centrality selection and triggermask, kINT7 (MB)
    if ((lPercentile > 0. && lPercentile < 100.) && fTriggerMask & AliVEvent::kINT7)
    {
        fHistos_eve->FillTH1("hcent", lPercentile);
    }
    else
    {
        DataPosting();
        return;
    }

    fHistos_eve->FillTH1("fCuts", 4.5);

    /// set boundaries for low, mid , high pt ranges
    ptXiBoundary_LowMid = 1.9;
    ptXiBoundary_MidHigh = 3.1;
    ptOmBoundary_LowMid = 2.;
    ptOmBoundary_MidHigh = 2.9;
    /// start cascade part if xi or omega analysis is requested
    if (fParticleAnalysisStatus[kxip] || fParticleAnalysisStatus[komp])
    {
        int ncasc = 0;
        ncasc = lAODevent->GetNumberOfCascades();
        for (int i_casc = 0; i_casc < ncasc; i_casc++)
        {
            /// start AOD part
            AliAODcascade *casc = lAODevent->GetCascade(i_casc);
            if (!casc)
                continue;

            /// cascade and V0 2D radii
            double lVtxCasc[3];
            lVtxCasc[0] = casc->DecayVertexXiX();
            lVtxCasc[1] = casc->DecayVertexXiY();
            lVtxCasc[2] = casc->DecayVertexXiZ();
            fCasc_CascRad = TMath::Sqrt(lVtxCasc[0] * lVtxCasc[0] + lVtxCasc[1] * lVtxCasc[1]);
            fCasc_V0Rad = casc->RadiusSecVtx();

            /// get daughter tracks (positive, negative and bachelor)
            AliAODTrack *pTrackCasc = dynamic_cast<AliAODTrack *>(casc->GetDaughter(0));
            AliAODTrack *nTrackCasc = dynamic_cast<AliAODTrack *>(casc->GetDaughter(1));
            AliAODTrack *bTrackCasc = dynamic_cast<AliAODTrack *>(casc->GetDecayVertexXi()->GetDaughter(0));
            if (!pTrackCasc || !nTrackCasc || !bTrackCasc)
            {
                AliWarning("ERROR: Could not retrieve one of the 3 AOD daughter tracks of the cascade ...\n");
                continue;
            }

            /// preliminary check
            if (pTrackCasc->GetTPCNclsF() <= 0 || nTrackCasc->GetTPCNclsF() <= 0 || bTrackCasc->GetTPCNclsF() <= 0)
                continue; /// check here to avoid division by zero later
            /// daughters' etas
            fCasc_etaPos = pTrackCasc->Eta();
            fCasc_etaNeg = nTrackCasc->Eta();
            fCasc_etaBac = bTrackCasc->Eta();

            /// PID
            fCasc_NSigPosProton = fPIDResponse->NumberOfSigmasTPC(pTrackCasc, AliPID::kProton);
            fCasc_NSigPosPion = fPIDResponse->NumberOfSigmasTPC(pTrackCasc, AliPID::kPion);
            fCasc_NSigNegProton = fPIDResponse->NumberOfSigmasTPC(nTrackCasc, AliPID::kProton);
            fCasc_NSigNegPion = fPIDResponse->NumberOfSigmasTPC(nTrackCasc, AliPID::kPion);
            fCasc_NSigBacPion = fPIDResponse->NumberOfSigmasTPC(bTrackCasc, AliPID::kPion);
            fCasc_NSigBacKaon = fPIDResponse->NumberOfSigmasTPC(bTrackCasc, AliPID::kKaon);

            /// crossed raws
            double lCrosRawsPos = pTrackCasc->GetTPCClusterInfo(2, 1);
            double lCrosRawsNeg = nTrackCasc->GetTPCClusterInfo(2, 1);
            double lCrosRawsBac = bTrackCasc->GetTPCClusterInfo(2, 1);
            fCasc_LeastCRaws = (int)(lCrosRawsPos < lCrosRawsNeg ? std::min(lCrosRawsPos, lCrosRawsBac) : std::min(lCrosRawsNeg, lCrosRawsBac));
            /// crossed raws / Findable clusters
            double lCrosRawsOvFPos = lCrosRawsPos / ((double)(pTrackCasc->GetTPCNclsF()));
            double lCrosRawsOvFNeg = lCrosRawsNeg / ((double)(nTrackCasc->GetTPCNclsF()));
            double lCrosRawsOvFBac = lCrosRawsBac / ((double)(bTrackCasc->GetTPCNclsF()));
            fCasc_LeastCRawsOvF = lCrosRawsOvFPos < lCrosRawsOvFNeg ? std::min(lCrosRawsOvFPos, lCrosRawsOvFBac) : std::min(lCrosRawsOvFNeg, lCrosRawsOvFBac);

            /// clusters for TPC PID
            double_t lTPCclsPos = pTrackCasc->GetTPCsignalN();
            double_t lTPCclsNeg = nTrackCasc->GetTPCsignalN();
            double_t lTPCclsBac = bTrackCasc->GetTPCsignalN();
            fCasc_LeastTPCcls = (int)(lTPCclsPos < lTPCclsNeg ? std::min(lTPCclsPos, lTPCclsBac) : std::min(lTPCclsNeg, lTPCclsBac));

            /// DCA info
            fCasc_DcaCascDaught = casc->DcaXiDaughters();
            fCasc_DcaBachToPV = casc->DcaBachToPrimVertex();
            fCasc_DcaPosToPV = casc->DcaPosToPrimVertex();
            fCasc_DcaNegToPV = casc->DcaNegToPrimVertex();
            fCasc_DcaV0Daught = casc->DcaV0Daughters();
            fCasc_DcaV0ToPV = casc->DcaV0ToPrimVertex();

            /// cascade and V0 cosine of pointing angle
            fCasc_CascCosPA = casc->CosPointingAngleXi((const Double_t &)lBestPV[0], (const Double_t &)lBestPV[1], (const Double_t &)lBestPV[2]);
            fCasc_V0CosPA = casc->CosPointingAngle(lBestPV);

            /// track status: ( fCasc_NegTrackStatus & AliESDtrack::kITSrefit ) is the codition to check kITSrefit
            fCasc_NegTrackStatus = nTrackCasc->GetStatus();
            fCasc_PosTrackStatus = pTrackCasc->GetStatus();
            fCasc_BacTrackStatus = bTrackCasc->GetStatus();

            /// check if at least one of candidate's daughter has a hit in the TOF or has ITSrefit flag (removes Out Of Bunch Pileup)
            fCasc_ITSTOFtracks = 0;
            if ((fCasc_NegTrackStatus & AliESDtrack::kITSrefit) || (nTrackCasc->GetTOFBunchCrossing(lMagField) > -95.))
                fCasc_ITSTOFtracks++;
            if ((fCasc_PosTrackStatus & AliESDtrack::kITSrefit) || (pTrackCasc->GetTOFBunchCrossing(lMagField) > -95.))
                fCasc_ITSTOFtracks++;
            if ((fCasc_BacTrackStatus & AliESDtrack::kITSrefit) || (bTrackCasc->GetTOFBunchCrossing(lMagField) > -95.))
                fCasc_ITSTOFtracks++;

            /// candidate's rapidity (mass hypothesis dependent)
            fCasc_yXi = casc->RapXi();
            fCasc_yOm = casc->RapOmega();

            /// charge
            fCasc_charge = (int)casc->ChargeXi();

            /// V0 daughter mass (later to be checked against nominal)
            if (fCasc_charge < 0)
            {
                fCasc_InvMassLam = casc->MassLambda();
            }
            else
            {
                fCasc_InvMassLam = casc->MassAntiLambda();
            }

            /// transverse momentum
            fCasc_Pt = TMath::Sqrt(casc->Pt2Xi());

            /// distance over total momentum
            fCasc_DistOverTotP = (TMath::Sqrt(TMath::Power(lVtxCasc[0] - lBestPV[0], 2) + TMath::Power(lVtxCasc[1] - lBestPV[1], 2) + TMath::Power(lVtxCasc[2] - lBestPV[2], 2))) / (TMath::Sqrt(casc->Ptot2Xi()) + 1e-10);

            /// candidate's invariant mass
            fCasc_InvMassXiMin = casc->MassXi();
            fCasc_InvMassXiPlu = casc->MassXi();
            fCasc_InvMassOmMin = casc->MassOmega();
            fCasc_InvMassOmPlu = casc->MassOmega();

            /// calculate DCA Bachelor-Baryon to remove "bump" structure in InvMass
            fCasc_BacBarCosPA = casc->BachBaryonCosPA();
        }

        /// Apply parametric BacBarCosPA cut, if requested
        if (fisParametricBacBarCosPA)
        {
            if (fCasc_Pt >= fHist_PtBacBarCosPA->GetXaxis()->GetXmin() && fCasc_Pt <= fHist_PtBacBarCosPA->GetXaxis()->GetXmax() && lPercentile < fCentLimit_BacBarCosPA)
            {
                SetCutVal(kFALSE, kTRUE, kCasc_BacBarCosPA, fHist_PtBacBarCosPA->GetBinContent(fHist_PtBacBarCosPA->GetXaxis()->FindBin(fCasc_Pt)));
            }
            else
            {
                if ((fParticleAnalysisStatus[kxip] || fParticleAnalysisStatus[komp]))
                    SetCutVal(kFALSE, kTRUE, kCasc_BacBarCosPA, fCasc_Cuts[kCasc_BacBarCosPA]);
            }
        }

        /// default cuts

        if (fParticleAnalysisStatus[kxip])
        {
            if (ApplyCuts(kxip))
            {
                fHistos_eve->FillTH1("fHistPtXiP", fCasc_Pt);
                fHistos_eve->FillTH1("fCasc_InvMassXiPlu", fCasc_InvMassXiPlu);
                fHistos_eve->FillTH1("fHistCentXiP", lPercentile);

                fHistos_XiPlu->FillTH3("h3_ptmasscent_def", fCasc_Pt, fCasc_InvMassXiPlu, lPercentile);
            }
        }
        if (fParticleAnalysisStatus[kxim])
        {
            if (ApplyCuts(kxim))
            {
                fHistos_eve->FillTH1("fHistPtXiM", fCasc_Pt);
                fHistos_eve->FillTH1("fCasc_InvMassXiMin", fCasc_InvMassXiMin);
                fHistos_eve->FillTH1("fHistCentXiM", lPercentile);

                fHistos_XiMin->FillTH3("h3_ptmasscent_def", fCasc_Pt, fCasc_InvMassXiMin, lPercentile);
            }
        }
        if (fParticleAnalysisStatus[komp])
        {
            if (ApplyCuts(komp))
            {
                fHistos_eve->FillTH1("fHistPtOmP", fCasc_Pt);
                fHistos_eve->FillTH1("fCasc_InvMassOmPlu", fCasc_InvMassOmPlu);
                fHistos_eve->FillTH1("fHistCentOmP", lPercentile);

                fHistos_OmPlu->FillTH3("h3_ptmasscent_def", fCasc_Pt, fCasc_InvMassOmPlu, lPercentile);
            }
        }
        if (fParticleAnalysisStatus[komm])
        {
            if (ApplyCuts(komm))
            {
                fHistos_eve->FillTH1("fHistPtOmM", fCasc_Pt);
                fHistos_eve->FillTH1("fCasc_InvMassOmMin", fCasc_InvMassOmMin);
                fHistos_eve->FillTH1("fHistCentOmM", lPercentile);

                fHistos_OmMin->FillTH3("h3_ptmasscent_def", fCasc_Pt, fCasc_InvMassOmMin, lPercentile);
            }
        }

        /// filling 3D histograms
        // if (!fDefOnly)
        //     FillHistCutVariations(kTRUE, lPercentile, 1, assFlag, 0);
    }

    DataPosting();
}

void AliCascadeAnalysis_Ishaan::SetCutVal(bool defchange, bool iscasc, int cutnum, double cval)
{

    cutval_Casc[cutnum] = cval;
    if (defchange)
    {
        fCasc_Cuts[cutnum] = cval;
    }
}

///________________________________________________________________________
void AliCascadeAnalysis_Ishaan::SetParametricBacBarCosPA(int nbins, float *ptbins, float *values, int cent_limit)
{
    fisParametricBacBarCosPA = kTRUE;
    fCentLimit_BacBarCosPA = cent_limit;
    fHist_PtBacBarCosPA = new TH1F("", "", nbins, ptbins);
    for (int iBin = 1; iBin <= nbins; iBin++)
        fHist_PtBacBarCosPA->SetBinContent(iBin, values[iBin - 1]);
}

///________________________________________________________________________
void AliCascadeAnalysis_Ishaan::SetDefCutVals()
{

    /// Cascade part
    if (fParticleAnalysisStatus[kxip] || fParticleAnalysisStatus[komp])
    {
        for (int iCasc_cut = 0; iCasc_cut < kCasccutsnum; iCasc_cut++)
        {
            SetCutVal(kFALSE, kTRUE, iCasc_cut, fCasc_Cuts[iCasc_cut]);
        }
    }
}

///________________________________________________________________________
void AliCascadeAnalysis_Ishaan::SetCutVariation(bool iscasc, int cutnum, int nvar, double lowval, double highval)
{

    {
        nvarcut_Casc[cutnum] = nvar;
        varlowcut_Casc[cutnum] = lowval;
        varhighcut_Casc[cutnum] = highval;
    }
}

///________________________________________________________________________
// void AliCascadeAnalysis_Ishaan::SetDefCutVariations()
// {

//     SetCutVariation(kTRUE, kCasc_DcaCascDaught, 3, 1.4, 1.6);
//     SetCutVariation(kTRUE, kCasc_CascCosPA, 2, 0.96, 0.97);
//     SetCutVariation(kTRUE, kCasc_CascRadXi, 2, 0.5, 0.6);
//     SetCutVariation(kTRUE, kCasc_CascRadOm, 2, 0.5, 0.6);
//     SetCutVariation(kTRUE, kCasc_NSigPID, 2, 4, 4);
//     // SetCutVariation(kTRUE, kCasc_LeastCRaws, 11, 70, 90);
//     // SetCutVariation(kTRUE, kCasc_LeastCRawsOvF, 11, 0.75, 0.9);
//     SetCutVariation(kTRUE, kCasc_LeastTPCcls, 2, 70, 70);
//     SetCutVariation(kTRUE, kCasc_InvMassLam, 2, 0.008, 0.008);
//     SetCutVariation(kTRUE, kCasc_DcaV0Daught, 3, 1.5, 1.7);
//     SetCutVariation(kTRUE, kCasc_V0CosPA, 2, 0.98, 0.98);
//     SetCutVariation(kTRUE, kCasc_DcaV0ToPV, 2, 0.06, 0.06);
//     SetCutVariation(kTRUE, kCasc_DcaBachToPV, 2, 0.04, 0.04);
//     // SetCutVariation(kTRUE, kCasc_ITSTOFtracks, 3, 1, 3);
//     SetCutVariation(kTRUE, kCasc_PropLifetXi, 7, 2, 5);
//     SetCutVariation(kTRUE, kCasc_PropLifetOm, 7, 2, 5);
//     SetCutVariation(kTRUE, kCasc_V0RadXi, 2, 1.1, 1.2);
//     SetCutVariation(kTRUE, kCasc_V0RadOm, 2, 1.1, 1.2);
//     SetCutVariation(kTRUE, kCasc_DcaMesToPV, 2, 0.04, 0.04);
//     SetCutVariation(kTRUE, kCasc_DcaBarToPV, 2, 0.03, 0.03);
//     SetCutVariation(kTRUE, kCasc_BacBarCosPA, 10, 0.999, 0.99999);
// }

///________________________________________________________________________
void AliCascadeAnalysis_Ishaan::SetDefOnly(bool isdefonly)
{
    fDefOnly = isdefonly;
}

///________________________________________________________________________
bool AliCascadeAnalysis_Ishaan::ApplyCuts(int part)
{
    /// we are checking cascades

    /// check candidate's charge
    if ((part == kxim || part == komm) && fCasc_charge > 0)
        return kFALSE;

    if ((part == kxip || part == komp) && fCasc_charge < 0)
        return kFALSE;

    /// check candidate daughters' pseudo-rapidity
    if (TMath::Abs(fCasc_etaPos) > cutval_Casc[kCasc_etaDaugh] || TMath::Abs(fCasc_etaNeg) > cutval_Casc[kCasc_etaDaugh] || TMath::Abs(fCasc_etaBac) > cutval_Casc[kCasc_etaDaugh])
        return kFALSE;
    /// check candidate daughters' crossed TPC raws (note that the checked value is the lowest among the daughters)
    if (fCasc_LeastCRaws < cutval_Casc[kCasc_LeastCRaws])
        return kFALSE;
    /// check candidate daughters' crossed TPC raws over findable
    if (fCasc_LeastCRawsOvF < cutval_Casc[kCasc_LeastCRawsOvF])
        return kFALSE;
    /// check candidate daughters' TPC clusters
    if (fCasc_LeastTPCcls < cutval_Casc[kCasc_LeastTPCcls])
        return kFALSE;

    /// check candidate V0 daughter's 2D decay distance from PV (if it is too small, then it's not a weak decay)
    if ((part == kxip || part == kxim) && fCasc_V0Rad < cutval_Casc[kCasc_V0RadXi])
        return kFALSE;
    if ((part == komp || part == komm) && fCasc_V0Rad < cutval_Casc[kCasc_V0RadOm])
        return kFALSE;

    /// check candidate daughters' DCA to Primary Vertex (needs to be large because decay is far from the Primary Vertex)
    if (fCasc_DcaBachToPV < cutval_Casc[kCasc_DcaBachToPV])
        return kFALSE;
    if (fCasc_DcaV0ToPV < cutval_Casc[kCasc_DcaV0ToPV])
        return kFALSE;
    /// check V0 daughters' DCA to Primary Vertex. Different cut for meson and baryon daughters, so different conditions for + and - candidates
    if ((part == kxim || part == komm) && (fCasc_DcaPosToPV < cutval_Casc[kCasc_DcaBarToPV] || fCasc_DcaNegToPV < cutval_Casc[kCasc_DcaMesToPV]))
        return kFALSE; /// in this case pos=p, neg=pi-, bach=pi- or K-
    if ((part == kxip || part == komp) && (fCasc_DcaPosToPV < cutval_Casc[kCasc_DcaMesToPV] || fCasc_DcaNegToPV < cutval_Casc[kCasc_DcaBarToPV]))
        return kFALSE; /// in this case pos=pi+, neg=anti-p, bach=pi+ or K+

    /// Xi: Particle dependent cuts
    if ((part == kxip || part == kxim))
    {

        /// check candidate's rapidity (particle hypothesis' dependent)
        if (fCasc_yXi < cutval_Casc[kCasc_y] && fCasc_yXi > 0.)
            return kFALSE;

        /// check candidate's 2D decay distance from PV (if it is too small, then it's not a weak decay)
        if (fCasc_CascRad < cutval_Casc[kCasc_CascRadXi])
            return kFALSE;

        if (fCasc_V0CosPA < cutval_Casc[kCasc_V0CosPAXi])
            return kFALSE;

        /// check candidate's proper lifetime (particle hypothesis' dependent). Remember: c*tau = L*m/p
        if (((1.32171 * fCasc_DistOverTotP) > (4.91 * cutval_Casc[kCasc_PropLifetXi])))
            return kFALSE; /// 4.91 is the ctau of xi in cm

        if (fCasc_Pt > 0. && fCasc_Pt < ptXiBoundary_LowMid) /// low pt range
        {
            /// check the cosine of the Pointing Angle for both cascade and V0 (angle between candidate's momentum and vector connecting Primary and secondary vertices)
            if (fCasc_CascCosPA < cutval_Casc[kCasc_CascCosPALow])
                return kFALSE;

            /// check V0 daughter's daughters DCA between them (needs to be small because they have to come from the same secondary vertex)
            if (fCasc_DcaV0Daught > cutval_Casc[kCasc_DcaV0DaughtLow])
                return kFALSE;

            if (fCasc_DcaCascDaught > cutval_Casc[kCasc_DcaCascDaughtLow]) /// check candidate daughter's DCA between them (needs to be small because they have to come from the same secondary vertex)
                return kFALSE;
        }
        if (fCasc_Pt >= ptXiBoundary_LowMid && fCasc_Pt < ptXiBoundary_MidHigh) /// mid pt range
        {
            if (fCasc_CascCosPA < cutval_Casc[kCasc_CascCosPAMid])
                return kFALSE;

            if (fCasc_DcaV0Daught > cutval_Casc[kCasc_DcaV0DaughtMid])
                return kFALSE;
            if (fCasc_DcaCascDaught > cutval_Casc[kCasc_DcaCascDaughtMid]) /// check candidate daughter's DCA between them (needs to be small because they have to come from the same secondary vertex)
                return kFALSE;
        }
        if (fCasc_Pt >= ptXiBoundary_MidHigh) /// high pt range
        {
            if (fCasc_CascCosPA < cutval_Casc[kCasc_CascCosPAHigh])
                return kFALSE;

            if (fCasc_DcaV0Daught > cutval_Casc[kCasc_DcaV0DaughtHigh])
                return kFALSE;
            if (fCasc_DcaCascDaught > cutval_Casc[kCasc_DcaCascDaughtHigh]) /// check candidate daughter's DCA between them (needs to be small because they have to come from the same secondary vertex)
                return kFALSE;
        }
    }

    /// Omega: Particle dependent cuts
    if ((part == komp || part == komm))
    {
        /// check candidate's rapidity (particle hypothesis' dependent)
        if ((part == komp || part == komm) && fCasc_yOm < cutval_Casc[kCasc_y] && fCasc_yOm > 0.)
            return kFALSE;

        /// check candidate's 2D decay distance from PV (if it is too small, then it's not a weak decay)
        if (fCasc_CascRad < cutval_Casc[kCasc_CascRadOm])
            return kFALSE;

        if (((1.67245 * fCasc_DistOverTotP) > (2.461 * cutval_Casc[kCasc_PropLifetOm])))
            return kFALSE; /// 2.461 is the ctau of om in cm

        if (fCasc_Pt > 0. && fCasc_Pt < ptOmBoundary_LowMid) /// low pt range
        {
            if (fCasc_CascCosPA < cutval_Casc[kCasc_CascCosPALow])
                return kFALSE;

            if (fCasc_V0CosPA < cutval_Casc[kCasc_V0CosPAOmLow])
                return kFALSE;

            /// check V0 daughter's daughters DCA between them (needs to be small because they have to come from the same secondary vertex)
            if (fCasc_DcaV0Daught > cutval_Casc[kCasc_DcaV0DaughtLow])
                return kFALSE;

            if (fCasc_DcaCascDaught > cutval_Casc[kCasc_DcaCascDaughtLow]) /// check candidate daughter's DCA between them (needs to be small because they have to come from the same secondary vertex)
                return kFALSE;
        }
        if (fCasc_Pt >= ptOmBoundary_LowMid && fCasc_Pt < ptOmBoundary_MidHigh) /// mid pt range
        {
            if (fCasc_V0CosPA < cutval_Casc[kCasc_V0CosPAOmMid])
                return kFALSE;

            if (fCasc_CascCosPA < cutval_Casc[kCasc_CascCosPAMid])
                return kFALSE;

            if (fCasc_DcaV0Daught > cutval_Casc[kCasc_DcaV0DaughtMid])
                return kFALSE;

            if (fCasc_DcaCascDaught > cutval_Casc[kCasc_DcaCascDaughtMid]) /// check candidate daughter's DCA between them (needs to be small because they have to come from the same secondary vertex)
                return kFALSE;
        }
        if (fCasc_Pt >= ptOmBoundary_MidHigh) /// high pt range
        {
            if (fCasc_V0CosPA < cutval_Casc[kCasc_V0CosPAOmHigh])
                return kFALSE;

            if (fCasc_CascCosPA < cutval_Casc[kCasc_CascCosPAHigh])
                return kFALSE;

            if (fCasc_DcaV0Daught > cutval_Casc[kCasc_DcaV0DaughtHigh])
                return kFALSE;

            if (fCasc_DcaCascDaught > cutval_Casc[kCasc_DcaCascDaughtHigh]) /// check candidate daughter's DCA between them (needs to be small because they have to come from the same secondary vertex)
                return kFALSE;
        }
    }

    /// check candidate V0 daughter's mass difference from nominal Lambda mass
    if (TMath::Abs(fCasc_InvMassLam - 1.115683) > cutval_Casc[kCasc_InvMassLam])
        return kFALSE;
    /// check that none of daughters is a kink
    // if (fCasc_kinkidx > 0)
    //   return kFALSE;
    /// check if at least one of candidate's daughter has a hit in the TOF or has ITSrefit flag (removes Out Of Bunch Pileup)
    // if (fCasc_ITSTOFtracks < cutval_Casc[kCasc_ITSTOFtracks])
    //   return kFALSE;
    /// TPC refit, should be already verified for Offline V0s
    if (!(fCasc_PosTrackStatus & AliESDtrack::kTPCrefit) ||
        !(fCasc_NegTrackStatus & AliESDtrack::kTPCrefit) ||
        !(fCasc_BacTrackStatus & AliESDtrack::kTPCrefit))
        return kFALSE;

    /// check DCA bachelor-baryon. If it is too small --> bump structure in Inv Mass
    // if (fCasc_BacBarCosPA > cutval_Casc[kCasc_BacBarCosPA])
    //     return kFALSE;
    /// check PID for all daughters (particle hypothesis' dependent)
    if ((part == kxip) && (TMath::Abs(fCasc_NSigPosPion) > cutval_Casc[kCasc_NSigPID] || TMath::Abs(fCasc_NSigNegProton) > cutval_Casc[kCasc_NSigPID] || TMath::Abs(fCasc_NSigBacPion) > cutval_Casc[kCasc_NSigPID]))
        return kFALSE;
    if ((part == kxim) && (TMath::Abs(fCasc_NSigNegPion) > cutval_Casc[kCasc_NSigPID] || TMath::Abs(fCasc_NSigPosProton) > cutval_Casc[kCasc_NSigPID] || TMath::Abs(fCasc_NSigBacPion) > cutval_Casc[kCasc_NSigPID]))
        return kFALSE;
    if ((part == komp) && (TMath::Abs(fCasc_NSigPosPion) > cutval_Casc[kCasc_NSigPID] || TMath::Abs(fCasc_NSigNegProton) > cutval_Casc[kCasc_NSigPID] || TMath::Abs(fCasc_NSigBacKaon) > cutval_Casc[kCasc_NSigPID]))
        return kFALSE;
    if ((part == komm) && (TMath::Abs(fCasc_NSigNegPion) > cutval_Casc[kCasc_NSigPID] || TMath::Abs(fCasc_NSigPosProton) > cutval_Casc[kCasc_NSigPID] || TMath::Abs(fCasc_NSigBacKaon) > cutval_Casc[kCasc_NSigPID]))
        return kFALSE;

    AliInfo("All cuts passed!\n");
    return kTRUE; /// survived!
}

void AliCascadeAnalysis_Ishaan::SetParticleAnalysisStatus(bool k0s, bool lambda, bool xi, bool omega)
{
    fParticleAnalysisStatus[kk0s] = k0s;
    fParticleAnalysisStatus[klam] = lambda;
    fParticleAnalysisStatus[kalam] = lambda;
    fParticleAnalysisStatus[kxip] = xi;
    fParticleAnalysisStatus[kxim] = xi;
    fParticleAnalysisStatus[komp] = omega;
    fParticleAnalysisStatus[komm] = omega;
}

bool AliCascadeAnalysis_Ishaan::GetParticleAnalysisStatus(int part)
{
    if (part < 0 || part >= ksignednumpart)
    {
        ::Error("AliCascadeAnalysis_Ishaan::GetParticleAnalysisStatus", "Wrong particle selected: accepted values from 0 to 3");
        return false;
    }
    return fParticleAnalysisStatus[part];
}

void AliCascadeAnalysis_Ishaan::SetCentbinning(int ipart, int numcentbins, double *centbins)
{
    fncentbins[ipart] = numcentbins;
    for (int i = 0; i < fncentbins[ipart] + 1; i++)
    {
        fcentbinning[ipart][i] = centbins[i];
    }
}

void AliCascadeAnalysis_Ishaan::SetPtbinning(int ipart, int numptbins, double *ptbins)
{
    fnptbins[ipart] = numptbins;
    for (int i = 0; i < fnptbins[ipart] + 1; i++)
    {
        fptbinning[ipart][i] = ptbins[i];
    }
}

void AliCascadeAnalysis_Ishaan::SetMassbinning(int ipart, int nummassbins, double valminmass, double valmaxmass)
{
    fnmassbins[ipart] = nummassbins;
    for (int i = 0; i < fnmassbins[ipart] + 1; i++)
    {
        fmassbinning[ipart][i] = valminmass + i * (valmaxmass - valminmass) / fnmassbins[ipart];
    }
}

///________________________________________________________________________
void AliCascadeAnalysis_Ishaan::DataPosting()
{

    PostData(1, fHistos_eve->GetListOfHistograms());

    /// Histograms for analysed particle specie
    int histnumber = 1;
    // if (fParticleAnalysisStatus[kk0s])
    // {
    //     histnumber++;
    //     PostData(histnumber, fHistos_K0S->GetListOfHistograms());
    // }
    // if (fParticleAnalysisStatus[klam])
    // {
    //     histnumber = histnumber + 2;
    //     PostData(histnumber - 1, fHistos_Lam->GetListOfHistograms());
    //     PostData(histnumber, fHistos_ALam->GetListOfHistograms());
    // }
    if (fParticleAnalysisStatus[kxip])
    {
        histnumber = histnumber + 2;
        PostData(histnumber - 1, fHistos_XiMin->GetListOfHistograms());
        PostData(histnumber, fHistos_XiPlu->GetListOfHistograms());
    }
    if (fParticleAnalysisStatus[komp])
    {
        histnumber = histnumber + 2;
        PostData(histnumber - 1, fHistos_OmMin->GetListOfHistograms());
        PostData(histnumber, fHistos_OmPlu->GetListOfHistograms());
    }
}

///________________________________________________________________________
// void AliCascadeAnalysis_Ishaan::FillHistCutVariations(bool iscasc, double perc, bool phypri, bool *associFlag, double ptassxi)
// {

//     for (int i_cut = 0; i_cut < kCasccutsnum; i_cut++)
//     {
//         if (i_cut == kCasc_y || i_cut == kCasc_etaDaugh)
//             continue;
//         for (int i_var = 0; i_var < nvarcut_Casc[i_cut]; i_var++)
//         {
//             if (fisParametricBacBarCosPA && i_cut != kCasc_BacBarCosPA && perc < fCentLimit_BacBarCosPA)
//             {
//                 if (fCasc_Pt >= fHist_PtBacBarCosPA->GetXaxis()->GetXmin() && fCasc_Pt <= fHist_PtBacBarCosPA->GetXaxis()->GetXmax())
//                 {
//                     SetCutVal(kFALSE, kTRUE, kCasc_BacBarCosPA, fHist_PtBacBarCosPA->GetBinContent(fHist_PtBacBarCosPA->GetXaxis()->FindBin(fCasc_Pt)));
//                 }
//                 else
//                 {
//                     SetCutVal(kFALSE, kTRUE, kCasc_BacBarCosPA, fCasc_Cuts[kCasc_BacBarCosPA]);
//                 }
//             }
//             ///Xi filling
//             if (i_cut != kCasc_PropLifetOm)
//             {
//                 if (fParticleAnalysisStatus[kxip])
//                 {
//                     SetCutVal(kFALSE, kTRUE, i_cut, varlowcut_Casc[i_cut] + i_var * (varhighcut_Casc[i_cut] - varlowcut_Casc[i_cut]) / (nvarcut_Casc[i_cut] - 1));
//                     if (phypri && associFlag[kxim] && ApplyCuts(kxim))
//                         fHistos_XiMin->FillTH3(Form("h3_ptmasscent[%d][%d]", i_cut, i_var), fCasc_Pt, fCasc_InvMassXiMin, perc);
//                     if (phypri && associFlag[kxip] && ApplyCuts(kxip))
//                         fHistos_XiPlu->FillTH3(Form("h3_ptmasscent[%d][%d]", i_cut, i_var), fCasc_Pt, fCasc_InvMassXiPlu, perc);
//                 }
//             }
//             ///Om filling
//             if (i_cut != kCasc_PropLifetXi)
//             {
//                 if (fParticleAnalysisStatus[komp])
//                 {
//                     SetCutVal(kFALSE, kTRUE, i_cut, varlowcut_Casc[i_cut] + i_var * (varhighcut_Casc[i_cut] - varlowcut_Casc[i_cut]) / (nvarcut_Casc[i_cut] - 1));
//                     if (phypri && associFlag[komm] && ApplyCuts(komm))
//                         fHistos_OmMin->FillTH3(Form("h3_ptmasscent[%d][%d]", i_cut, i_var), fCasc_Pt, fCasc_InvMassOmMin, perc);
//                     if (phypri && associFlag[komp] && ApplyCuts(komp))
//                         fHistos_OmPlu->FillTH3(Form("h3_ptmasscent[%d][%d]", i_cut, i_var), fCasc_Pt, fCasc_InvMassOmPlu, perc);
//                 }
//             }
//         }
//         SetDefCutVals(); ///reset defaults
//     }
// }

void AliCascadeAnalysis_Ishaan::Terminate(Option_t *)
{
}
