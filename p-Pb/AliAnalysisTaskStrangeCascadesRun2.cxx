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
#include "AliEventCuts.h"
#include "AliESDtrackCuts.h"
#include "TError.h"

#include "AliAnalysisTaskStrangeCascadesRun2.h"

#include "TMath.h"

ClassImp(AliAnalysisTaskStrangeCascadesRun2)

    AliAnalysisTaskStrangeCascadesRun2::AliAnalysisTaskStrangeCascadesRun2() : AliAnalysisTaskSE(),
                                                                 /// outputs
                                                                 fHistos_eve(nullptr),
                                                                 fHistos_XiMin(nullptr),
                                                                 fHistos_XiPlu(nullptr),
                                                                 fHistos_OmMin(nullptr),
                                                                 fHistos_OmPlu(nullptr),
                                                                 /// objects from the manager
                                                                 fPIDResponse(0),
                                                                 fTriggerMask(0),

                                                                 // AliEventCuts object
                                                                 fEventCuts(0),
                                                                 // pile-up rejection flag
                                                                 fPileupCut(0),

                                                                 // MC-related variables
                                                                 fisMC(kFALSE),
                                                                 fisMCassoc(kFALSE),

                                                                 /// default cuts configuration
                                                                 fDefOnly(kFALSE),
                                                                 //  fCasc_Cuts{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
                                                                 /// particle to be analysed
                                                                 fParticleAnalysisStatus{true, true, true, true},
                                                                 // geometrical cut usage
                                                                 fESDTrackCuts(0),
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
                                                                 fCasc_LeastCRows(0),
                                                                 fCasc_LeastCRowsOvF(0),
                                                                 fCasc_LeastTPCcls(0),
                                                                 fCasc_TrackLengthCut(0),
                                                                 fCasc_MaxChi2perCls(0),

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
                                                                 fCentLimit_BacBarCosPA(0),
                                                                 fisParametricTrackLengthCut(kFALSE),
                                                                 fHist_CentTrackLengthCut(0),
                                                                 //  fHistCutsEv(0),
                                                                 //  fHistCutsTopo(0),
                                                                 fDeadZoneWidth_GeoCut(0),
                                                                 fNcrNclLength_GeoCut(0),
                                                                 fTPCsignalNCut(0),
                                                                 fCasc_TrackLength(0) /*,
                                                                  fncentbins[](0)*/

{
    ////default constructor
}

AliAnalysisTaskStrangeCascadesRun2::AliAnalysisTaskStrangeCascadesRun2(const char *name, TString lExtraOptions) : AliAnalysisTaskSE(name),
                                                                                                    ////outputs
                                                                                                    fHistos_eve(nullptr),

                                                                                                    fHistos_XiMin(nullptr),
                                                                                                    fHistos_XiPlu(nullptr),
                                                                                                    fHistos_OmMin(nullptr),
                                                                                                    fHistos_OmPlu(nullptr),
                                                                                                    ////objects from the manager
                                                                                                    fPIDResponse(0),
                                                                                                    fTriggerMask(0),

                                                                                                    // AliEventCuts object
                                                                                                    fEventCuts(0),
                                                                                                    // pile-up rejection flag
                                                                                                    fPileupCut(1),

                                                                                                    ////MC-related variables - make true for MC only
                                                                                                    fisMC(kFALSE),
                                                                                                    fisMCassoc(kFALSE),

                                                                                                    ////default cuts configuration
                                                                                                    fDefOnly(kFALSE),
                                                                                                    // fCasc_Cuts{1.5, 0.96, 0.6, 4., 70., 0.8, 70., 1., 2.5, 0.008, 1.6, 0.98, 0.06, 0.04, 1., -0.5, 0.8, 3., 3., 1.2, 0.04, 0.03, 1., 0.5, 1.1, 1.6, 1.4, 0.97, 0.97, 1.7, 1.5, 0.97, 0.98, 0.98, 0.008, 0.98, 0.98, 0.6, 0.6, 1.2, 1.2},
                                                                                                    ////particle to be analysed
                                                                                                    fParticleAnalysisStatus{true, true, true, true},

                                                                                                    // geometrical cut usage
                                                                                                    fESDTrackCuts(0),
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
                                                                                                    fCasc_LeastCRows(0),
                                                                                                    fCasc_LeastCRowsOvF(0),
                                                                                                    fCasc_LeastTPCcls(0),
                                                                                                    fCasc_TrackLengthCut(0),
                                                                                                    fCasc_MaxChi2perCls(0),

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
                                                                                                    fCentLimit_BacBarCosPA(0),
                                                                                                    fisParametricTrackLengthCut(kFALSE),
                                                                                                    fHist_CentTrackLengthCut(0),
                                                                                                    // fHistCutsEv(0),
                                                                                                    // fHistCutsTopo(0),
                                                                                                    fDeadZoneWidth_GeoCut(0),
                                                                                                    fNcrNclLength_GeoCut(0),
                                                                                                    fTPCsignalNCut(50),
                                                                                                    fCasc_TrackLength(90) /*,
                                                                                                    fncentbins[](0)*/

{

    ////setting default cuts
    SetDefCuts();
    SetDefCutVariations();
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
    int massbins[2] = {100, 100};
    double minmass[2] = {1.272, 1.622};
    double maxmass[2] = {1.372, 1.722};
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
    DefineOutput(2, TList::Class()); //// XiMinus Histograms
    DefineOutput(3, TList::Class()); //// XiPlus Histograms
    DefineOutput(4, TList::Class()); //// OmegaMinus Histograms
    DefineOutput(5, TList::Class()); //// OmegaPlus Histograms
}

AliAnalysisTaskStrangeCascadesRun2::~AliAnalysisTaskStrangeCascadesRun2()
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
void AliAnalysisTaskStrangeCascadesRun2::UserCreateOutputObjects()
{

    ////histograms for event variables
    fHistos_eve = new THistManager("histos_eve");

    fHistos_eve->CreateTH1("hcent", "Multiplicity Distribution", 100, 0, 100, "s"); ////storing #events in bins of centrality
    fHistos_eve->CreateTH1("henum", "", 4, -0.5, 3.5);                              ////storing total #events
    // fHistos_eve->CreateTH1("hCheckCuts", "", 200, -1, 4);                          ////storing total #events

    /// test
    fHistos_eve->CreateTH1("hTrackLengthP", "pTrack Length (cm)", 200, 0, 200, "s"); ////check track length for cut
    fHistos_eve->CreateTH1("hTrackLengthN", "nTrack Length (cm)", 200, 0, 200, "s"); ////check track length for cut
    fHistos_eve->CreateTH1("hTrackLengthB", "bTrack Length (cm)", 200, 0, 200, "s"); ////check track length for cut

    const char *labels[4] = {"Total", "MultSelection", "AliEventCuts", "Pile-up rejection"};
    for (int iLab = 1; iLab <= 4; iLab++)
        ((TH1 *)fHistos_eve->FindObject("henum"))->GetXaxis()->SetBinLabel(iLab, labels[iLab - 1]);
    fHistos_eve->CreateTH1("fCuts", "Accepted events after event selection", 6, 0, 6, "s"); ////storing impact of cuts on no. of events

    TH1 *fCuts = (TH1 *)fHistos_eve->FindObject("fCuts");
    fCuts->GetXaxis()->SetBinLabel(1, "Total Events");
    fCuts->GetXaxis()->SetBinLabel(2, "V0A_MultSel");
    fCuts->GetXaxis()->SetBinLabel(3, "InelGT0");
    fCuts->GetXaxis()->SetBinLabel(4, "|PVz|<10");
    fCuts->GetXaxis()->SetBinLabel(5, "Mult<0,100> && kINT7");
    fCuts->GetXaxis()->SetBinLabel(6, "Track Length>90cm");

    ////storing impact of evCuts on no. of events
    fHistos_eve->CreateTH1("fHistCutsEv", "Event Selection cuts", kNumCascEvCuts + 1, 0, kNumCascEvCuts + 1, "s");                                                               // last bin (entry 12.5) for kTPCrefit
    fHistos_eve->CreateTH3("fHistCutsTopo", "Topological Selection cuts", knumpart, 0, knumpart, kNumCascTopoCuts, 0, kNumCascTopoCuts, kNumPtInterval, 0, kNumPtInterval, "s"); ////storing impact of topoCuts on no. of events

    TH1 *fHistCutsEv = (TH1 *)fHistos_eve->FindObject("fHistCutsEv");
    for (int iCutEv = 0; iCutEv < kNumCascEvCuts; iCutEv++)
    {
        fHistCutsEv->GetXaxis()->SetBinLabel(iCutEv + 1, cutNamesEv[iCutEv].Data());
    }
    fHistCutsEv->GetXaxis()->SetBinLabel(13, "kTPCrefit");

    TH3 *fHistCutsTopo = (TH3 *)fHistos_eve->FindObject("fHistCutsTopo");
    fHistCutsTopo->GetXaxis()->SetBinLabel(1, "#Xi");
    fHistCutsTopo->GetXaxis()->SetBinLabel(2, "#Omega");
    fHistCutsTopo->GetZaxis()->SetBinLabel(1, "Low p_{T}");
    fHistCutsTopo->GetZaxis()->SetBinLabel(2, "Mid p_{T}");
    fHistCutsTopo->GetZaxis()->SetBinLabel(3, "High p_{T}");
    for (int iCutTopo = 0; iCutTopo < kNumCascTopoCuts; iCutTopo++)
    {
        fHistCutsTopo->GetYaxis()->SetBinLabel(iCutTopo + 1, cutNamesTopo[iCutTopo].Data());
    }

    ////histograms for Cascade variables
    if (fParticleAnalysisStatus[kxip])
    {
        fHistos_XiMin = new THistManager("histos_XiMin");
        fHistos_XiPlu = new THistManager("histos_XiPlu");

        fHistos_XiPlu->CreateTH1("fHistPtXiP", "fHistPtXiP", fnptbins[kXi], fptbinning[kXi], "s");
        fHistos_XiPlu->CreateTH1("fCasc_InvMassXiPlu", "fCasc_InvMassXiPlu", 100, 1.272, 1.372, "s");
        fHistos_XiPlu->CreateTH1("fHistCentXiP", "fHistCentXiP", fncentbins[kXi], fcentbinning[kXi], "s");
        fHistos_XiPlu->CreateTH3("h3_ptmasscent_def", "", fnptbins[kXi], fptbinning[kXi], fnmassbins[kXi], fmassbinning[kXi], fncentbins[kXi], fcentbinning[kXi]);

        fHistos_XiMin->CreateTH1("fHistPtXiM", "fHistPtXiM", fnptbins[kXi], fptbinning[kXi], "s");
        fHistos_XiMin->CreateTH1("fCasc_InvMassXiMin", "fCasc_InvMassXiMin", 100, 1.272, 1.372, "s");
        fHistos_XiMin->CreateTH1("fHistCentXiM", "fHistCentXiM", fncentbins[kXi], fcentbinning[kXi], "s");
        fHistos_XiMin->CreateTH3("h3_ptmasscent_def", "", fnptbins[kXi], fptbinning[kXi], fnmassbins[kXi], fmassbinning[kXi], fncentbins[kXi], fcentbinning[kXi]);
        if (fisMC)
        {
            fHistos_XiMin->CreateTH2("h2_gen", "", fnptbins[kXi], fptbinning[kXi], fncentbins[kXi], fcentbinning[kXi]);
            fHistos_XiPlu->CreateTH2("h2_gen", "", fnptbins[kXi], fptbinning[kXi], fncentbins[kXi], fcentbinning[kXi]);
        }
    }
    if (fParticleAnalysisStatus[komp])
    {
        fHistos_OmMin = new THistManager("histos_OmMin");
        fHistos_OmPlu = new THistManager("histos_OmPlu");

        fHistos_OmPlu->CreateTH1("fHistPtOmP", "fHistPtOmP", fnptbins[kOm], fptbinning[kOm], "s");
        fHistos_OmPlu->CreateTH1("fCasc_InvMassOmPlu", "fCasc_InvMassOmPlu", 100, 1.622, 1.722, "s");
        fHistos_OmPlu->CreateTH1("fHistCentOmP", "fHistCentOmP", fncentbins[kOm], fcentbinning[kOm], "s");
        fHistos_OmPlu->CreateTH3("h3_ptmasscent_def", "", fnptbins[kOm], fptbinning[kOm], fnmassbins[kOm], fmassbinning[kOm], fncentbins[kOm], fcentbinning[kOm]);

        fHistos_OmMin->CreateTH1("fHistPtOmM", "fHistPtOmM", fnptbins[kOm], fptbinning[kOm], "s");
        fHistos_OmMin->CreateTH1("fCasc_InvMassOmMin", "fCasc_InvMassOmMin", 100, 1.622, 1.722, "s");
        fHistos_OmMin->CreateTH1("fHistCentOmM", "fHistCentOmM", fncentbins[kOm], fcentbinning[kOm], "s");
        fHistos_OmMin->CreateTH3("h3_ptmasscent_def", "", fnptbins[kOm], fptbinning[kOm], fnmassbins[kOm], fmassbinning[kOm], fncentbins[kOm], fcentbinning[kOm]);
        if (fisMC)
        {
            fHistos_OmMin->CreateTH2("h2_gen", "", fnptbins[kOm], fptbinning[kOm], fncentbins[kOm], fcentbinning[kOm]);
            fHistos_OmPlu->CreateTH2("h2_gen", "", fnptbins[kOm], fptbinning[kOm], fncentbins[kOm], fcentbinning[kOm]);
        }
    }

    /// creating histograms for cut variations
    if (!fDefOnly && (fParticleAnalysisStatus[kxip] || fParticleAnalysisStatus[komp]))
    {

        for (int iCutEv = 0; iCutEv < kNumCascEvCuts; iCutEv++)
        {
            if (iCutEv == kRapidityIntervalMin || iCutEv == kRapidityIntervalMax || iCutEv == kLeastCRows || iCutEv == kLeastCRowsOvF || iCutEv == kTrackLengthCut || iCutEv == kEtaDaughter || iCutEv == kBacBarCosPa)
                continue;
            if (nvarcut_Ev[iCutEv] == -1) // skip the cut if value == -1 -> unused variations
                continue;
            for (int iVarEv = 0; iVarEv < nvarcut_Ev[iCutEv]; iVarEv++)
            {
                if (iCutEv != kCompetingCascRejectOm)
                {
                    if (fParticleAnalysisStatus[kxip])
                    {
                        // if (icut != kCasc_PropLifetOm)
                        // if (icut != kCasc_PropLifetOm && icut != kCompetingCascRejectOm)
                        fHistos_XiMin->CreateTH3(TString::Format("h3Var_%s[%d][%d]", cutNamesEv[iCutEv].Data(), iCutEv, iVarEv), TString::Format("h3Var_%s[%d][%d]", cutNamesEv[iCutEv].Data(), iCutEv, iVarEv), fnptbins[kXi], fptbinning[kXi], fnmassbins[kXi], fmassbinning[kXi], fncentbins[kXi], fcentbinning[kXi]);
                        // if (icut != kCasc_PropLifetOm)
                        // if (icut != kCasc_PropLifetOm && icut != kCompetingCascRejectOm)
                        fHistos_XiPlu->CreateTH3(TString::Format("h3Var_%s[%d][%d]", cutNamesEv[iCutEv].Data(), iCutEv, iVarEv), TString::Format("h3Var_%s[%d][%d]", cutNamesEv[iCutEv].Data(), iCutEv, iVarEv), fnptbins[kXi], fptbinning[kXi], fnmassbins[kXi], fmassbinning[kXi], fncentbins[kXi], fcentbinning[kXi]);
                    }
                }
                if (fParticleAnalysisStatus[komp])
                {
                    // if (icut != kCasc_PropLifetXi)
                    fHistos_OmMin->CreateTH3(TString::Format("h3Var_%s[%d][%d]", cutNamesEv[iCutEv].Data(), iCutEv, iVarEv), TString::Format("h3Var_%s[%d][%d]", cutNamesEv[iCutEv].Data(), iCutEv, iVarEv), fnptbins[kOm], fptbinning[kOm], fnmassbins[kOm], fmassbinning[kOm], fncentbins[kOm], fcentbinning[kOm]);
                    // if (icut != kCasc_PropLifetXi)
                    fHistos_OmPlu->CreateTH3(TString::Format("h3Var_%s[%d][%d]", cutNamesEv[iCutEv].Data(), iCutEv, iVarEv), TString::Format("h3Var_%s[%d][%d]", cutNamesEv[iCutEv].Data(), iCutEv, iVarEv), fnptbins[kOm], fptbinning[kOm], fnmassbins[kOm], fmassbinning[kOm], fncentbins[kOm], fcentbinning[kOm]);
                }
            }
        }

        for (int iCutTopo = 0; iCutTopo < kNumCascTopoCuts; iCutTopo++)
        {
            if (nvarcut_Topo[iCutTopo] == -1) // skip the cut if value == -1 -> unused variations
                continue;
            for (int iVarTopo = 0; iVarTopo < nvarcut_Topo[iCutTopo]; iVarTopo++)
            {

                // if (fParticleAnalysisStatus[kxip] && (var_cutValTopo[kXi][iCutTopo][kMid][iVarTopo] != -1))
                // {
                if (fParticleAnalysisStatus[kxip])
                {
                    // if (icut != kCasc_PropLifetOm)
                    // if (icut != kCasc_PropLifetOm && icut != kCompetingCascRejectOm)
                    fHistos_XiMin->CreateTH3(TString::Format("h3Var_%s[%d][%d]", cutNamesTopo[iCutTopo].Data(), iCutTopo, iVarTopo), TString::Format("h3Var_%s[%d][%d]", cutNamesTopo[iCutTopo].Data(), iCutTopo, iVarTopo), fnptbins[kXi], fptbinning[kXi], fnmassbins[kXi], fmassbinning[kXi], fncentbins[kXi], fcentbinning[kXi]);
                    // if (icut != kCasc_PropLifetOm)
                    // if (icut != kCasc_PropLifetOm && icut != kCompetingCascRejectOm)
                    fHistos_XiPlu->CreateTH3(TString::Format("h3Var_%s[%d][%d]", cutNamesTopo[iCutTopo].Data(), iCutTopo, iVarTopo), TString::Format("h3Var_%s[%d][%d]", cutNamesTopo[iCutTopo].Data(), iCutTopo, iVarTopo), fnptbins[kXi], fptbinning[kXi], fnmassbins[kXi], fmassbinning[kXi], fncentbins[kXi], fcentbinning[kXi]);
                }
                //
                if (fParticleAnalysisStatus[komp])
                {
                    // if (icut != kCasc_PropLifetXi)
                    fHistos_OmMin->CreateTH3(TString::Format("h3Var_%s[%d][%d]", cutNamesTopo[iCutTopo].Data(), iCutTopo, iVarTopo), TString::Format("h3Var_%s[%d][%d]", cutNamesTopo[iCutTopo].Data(), iCutTopo, iVarTopo), fnptbins[kOm], fptbinning[kOm], fnmassbins[kOm], fmassbinning[kOm], fncentbins[kOm], fcentbinning[kOm]);
                    // if (icut != kCasc_PropLifetXi)
                    fHistos_OmPlu->CreateTH3(TString::Format("h3Var_%s[%d][%d]", cutNamesTopo[iCutTopo].Data(), iCutTopo, iVarTopo), TString::Format("h3Var_%s[%d][%d]", cutNamesTopo[iCutTopo].Data(), iCutTopo, iVarTopo), fnptbins[kOm], fptbinning[kOm], fnmassbins[kOm], fmassbinning[kOm], fncentbins[kOm], fcentbinning[kOm]);
                }
            }
        }
    }
    //// PID Setup
    AliAnalysisManager *man = AliAnalysisManager::GetAnalysisManager();
    AliInputEventHandler *inputHandler = (AliInputEventHandler *)(man->GetInputEventHandler());
    fPIDResponse = inputHandler->GetPIDResponse();
    inputHandler->SetNeedField();

    // fEventCuts Setup
    if (fPileupCut == 1)
        fEventCuts.SetRejectTPCPileupWithITSTPCnCluCorr(kTRUE);
    if (fPileupCut == 2)
        fEventCuts.SetRejectTPCPileupWithV0CentTPCnTracksCorr(kTRUE);

    ////Output posting
    DataPosting();
}
//_____________________________________________________________________________
void AliAnalysisTaskStrangeCascadesRun2::UserExec(Option_t *)
{
    // ignore all warnings
    gErrorIgnoreLevel = kError;

    // get event from the input handler and cast it into the desired type of event
    AliVEvent *lVevent = dynamic_cast<AliVEvent *>(InputEvent());
    if (!lVevent)
    {
        AliWarning("ERROR: event not available \n");
        DataPosting();
        return;
    }

    AliAODEvent *lAODevent = 0x0;

    if (lVevent->InheritsFrom("AliAODEvent"))
    {
        lAODevent = dynamic_cast<AliAODEvent *>(InputEvent());
    }
    else if (!lAODevent)
    {
        AliWarning("ERROR: event does not inherit from AliAODEvent! \n");
        DataPosting();
        return;
    }

    // AliAODEvent *lAODevent = dynamic_cast<AliAODEvent *>(InputEvent());
    // if (!lAODevent)
    // {
    //     AliWarning("ERROR: event not available \n");
    //     DataPosting();
    //     return;
    // }

    // get MC event
    AliMCEvent *lMCev = 0x0;
    if (fisMC)
    {
        lMCev = MCEvent();
        if (!lMCev)
        {
            Printf("ERROR: Could not retrieve MC event in file %s\n", fInputHandler->GetTree()->GetCurrentFile()->GetName());
            DataPosting();
            return;
        }
    }

    // get MC header and MC array
    AliAODMCHeader *header = 0x0;
    TClonesArray *MCTrackArray = 0x0;
    if (fisMC)
    {
        header = static_cast<AliAODMCHeader *>(lAODevent->FindListObject(AliAODMCHeader::StdBranchName()));
        if (!header)
        {
            AliWarning("No MC header found.");
            DataPosting();
            return;
        }
        MCTrackArray = dynamic_cast<TClonesArray *>(lAODevent->FindListObject(AliAODMCParticle::StdBranchName()));
        if (MCTrackArray == NULL)
        {
            AliWarning("No MC track array found.");
            DataPosting();
            return;
        }
    }

    ////dumb histo for checking
    fHistos_eve->FillTH1("henum", 0.);
    fHistos_eve->FillTH1("fCuts", 0.5);

    ////get trigger information
    fTriggerMask = ((AliInputEventHandler *)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected();

    //// Multiplicity Information
    double lPercentile = -666;
    int lEvSelCode = -666;
    AliMultSelection *MultSelection = (AliMultSelection *)lVevent->FindListObject("MultSelection");
    if (!MultSelection)
    {
        AliWarning("AliMultSelection object not found!");
        DataPosting();
        return;
    }
    else
    {
        lPercentile = MultSelection->GetMultiplicityPercentile("V0A"); /// changed V0M to V0A
        lEvSelCode = MultSelection->GetEvSelCode();                    ////==0 means event is good. Set by AliMultSelectionTask
    }

    ////skip everything if event selection code !=0
    if (lEvSelCode != 0)
    {
        DataPosting();
        return;
    }

    // fill number of events after AliMultSelection
    fHistos_eve->FillTH1("henum", 1.);
    fHistos_eve->FillTH1("fCuts", 1.5);

    bool isEvtAccepted = fEventCuts.AcceptEvent(lVevent);

    if (!isEvtAccepted && fEventCuts.PassedCut(AliEventCuts::kTPCPileUp))
    {
        DataPosting();
        return;
    }

    // fill number of events after event cuts
    fHistos_eve->FillTH1("henum", 2.);

    if (fPileupCut)
    {
        if (!fEventCuts.PassedCut(AliEventCuts::kTPCPileUp))
        {
            DataPosting();
            return;
        }
        if (fisMC)
        {
            if (AliAnalysisUtils::IsPileupInGeneratedEvent(header, "ijing")) // reject HIJING (Heavy Ion Jet INteraction Generator)
            {
                DataPosting();
                return;
            }
        }
    }

    // fill number of events after pile-up rejection
    fHistos_eve->FillTH1("henum", 3.);

    /// Inel>0 in the |eta|<1 range
    // check validity of this step (i.e. if it is needed)?
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

    // MC truth
    if (fisMC)
    {
        for (int i_MCtrk = 0; i_MCtrk < lMCev->GetNumberOfTracks(); i_MCtrk++)
        {
            AliVParticle *lPart;
            lPart = (AliAODMCParticle *)lMCev->GetTrack(i_MCtrk);
            bool isOOBpileup = kFALSE;
            isOOBpileup = AliAnalysisUtils::IsParticleFromOutOfBunchPileupCollision(i_MCtrk, header, MCTrackArray);
            if (!lPart || lPart->Y() < -0.5 || lPart->Y() > 0.5 || !lPart->IsPhysicalPrimary() || isOOBpileup)
                continue;

            if (fParticleAnalysisStatus[kxim] && lPart->PdgCode() == 3312)
                fHistos_XiMin->FillTH2("h2_gen", lPart->Pt(), lPercentile);
            if (fParticleAnalysisStatus[kxip] && lPart->PdgCode() == -3312)
                fHistos_XiPlu->FillTH2("h2_gen", lPart->Pt(), lPercentile);
            if (fParticleAnalysisStatus[komm] && lPart->PdgCode() == 3334)
                fHistos_OmMin->FillTH2("h2_gen", lPart->Pt(), lPercentile);
            if (fParticleAnalysisStatus[komp] && lPart->PdgCode() == -3334)
                fHistos_OmPlu->FillTH2("h2_gen", lPart->Pt(), lPercentile);
        }
    }

    Double_t cutPrimVertex = 10.0;

    ////acquire best PV
    // const AliVVertex *lBestPrimVtx = lAODevent->GetPrimaryVertex();
    const AliVVertex *lBestPrimVtx = lVevent->GetPrimaryVertex();
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
    lMagField = lVevent->GetMagneticField();
    // MC association
    int pdgPosDaught = 0;
    int pdgNegDaught = 0;
    int pdgBachelor = 0;
    int pdgV0 = 0;
    int pdgCasc = 0;
    //  bool assFlag[ksignednumpart] = {1, 1, 1, 1, 1, 1, 1}; ////MC ass flags are true if associated and, by construction, always for data. They can be false only if fisMC && notassociated
    bool assFlag[ksignednumpart]; // MC ass flags are true if associated and, by construction, always for data. They can be false only if fisMC && notassociated
    // double fdmtx_ptxi = 0;        // 0 if no secondary lambda, value corresponding to generated-xi pT, - for lambda (from xim) and + for anti-lambda (from xip)

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
        ncasc = lVevent->GetNumberOfCascades();
        for (int i_casc = 0; i_casc < ncasc; i_casc++)
        {
            // MC ass flag reset
            for (int iflag = kxip; iflag < ksignednumpart; iflag++)
                assFlag[iflag] = kTRUE;
            bool physprim = kTRUE;

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

            /// crossed Rows
            double lCrosRowsPos = pTrackCasc->GetTPCClusterInfo(2, 1);
            double lCrosRowsNeg = nTrackCasc->GetTPCClusterInfo(2, 1);
            double lCrosRowsBac = bTrackCasc->GetTPCClusterInfo(2, 1);
            fCasc_LeastCRows = (int)(lCrosRowsPos < lCrosRowsNeg ? std::min(lCrosRowsPos, lCrosRowsBac) : std::min(lCrosRowsNeg, lCrosRowsBac));
            /// crossed Rows / Findable clusters
            double lCrosRowsOvFPos = lCrosRowsPos / ((double)(pTrackCasc->GetTPCNclsF()));
            double lCrosRowsOvFNeg = lCrosRowsNeg / ((double)(nTrackCasc->GetTPCNclsF()));
            double lCrosRowsOvFBac = lCrosRowsBac / ((double)(bTrackCasc->GetTPCNclsF()));
            fCasc_LeastCRowsOvF = lCrosRowsOvFPos < lCrosRowsOvFNeg ? std::min(lCrosRowsOvFPos, lCrosRowsOvFBac) : std::min(lCrosRowsOvFNeg, lCrosRowsOvFBac);

            /// clusters for TPC PID
            double_t lTPCclsPos = pTrackCasc->GetTPCsignalN();
            double_t lTPCclsNeg = nTrackCasc->GetTPCsignalN();
            double_t lTPCclsBac = bTrackCasc->GetTPCsignalN();
            fCasc_LeastTPCcls = (int)(lTPCclsPos < lTPCclsNeg ? std::min(lTPCclsPos, lTPCclsBac) : std::min(lTPCclsNeg, lTPCclsBac));

            /// New: check effectiveness?
            // track length cut
            fCasc_TrackLengthCut = (pTrackCasc->GetTPCsignalN() > fTPCsignalNCut && nTrackCasc->GetTPCsignalN() > fTPCsignalNCut && bTrackCasc->GetTPCsignalN() > fTPCsignalNCut) ? 1 : 0;
            // if (fESDTrackCuts.AcceptVTrack(pTrackCasc) && fESDTrackCuts.AcceptVTrack(nTrackCasc) && fESDTrackCuts.AcceptVTrack(bTrackCasc))
            //     fCasc_TrackLengthCut = fCasc_TrackLengthCut + 2;
            Float_t lTrackLengthPos, lTrackLengthNeg, lTrackLengthBac = 0;
            lTrackLengthPos = GetLengthInActiveZone(pTrackCasc, 2.0, 220.0, lMagField);
            lTrackLengthNeg = GetLengthInActiveZone(nTrackCasc, 2.0, 220.0, lMagField);
            lTrackLengthBac = GetLengthInActiveZone(bTrackCasc, 2.0, 220.0, lMagField);

            /// test:
            fHistos_eve->FillTH1("hTrackLengthP", lTrackLengthPos);
            fHistos_eve->FillTH1("hTrackLengthN", lTrackLengthNeg);
            fHistos_eve->FillTH1("hTrackLengthB", lTrackLengthBac);

            if ((lTrackLengthPos > fCasc_TrackLength) && (lTrackLengthNeg > fCasc_TrackLength))
            {
                fCasc_TrackLengthCut = fCasc_TrackLengthCut + 2;
                fHistos_eve->FillTH1("fCuts", 5.5);
            }

            // chi^2 per TPC cluster
            // double_t lChi2perTPCclsPos = pTrackCasc->GetTPCchi2() / pTrackCasc->GetNcls(1);
            // double_t lChi2perTPCclsNeg = nTrackCasc->GetTPCchi2() / nTrackCasc->GetNcls(1);
            // double_t lChi2perTPCclsBac = bTrackCasc->GetTPCchi2() / bTrackCasc->GetNcls(1);
            // fCasc_MaxChi2perCls = (lChi2perTPCclsPos > lChi2perTPCclsNeg ? TMath::Max(lChi2perTPCclsPos, lChi2perTPCclsBac) : TMath::Max(lChi2perTPCclsNeg, lChi2perTPCclsBac));

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

            // /// check if at least one of candidate's daughter has a hit in the TOF or has ITSrefit flag (removes Out Of Bunch Pileup)
            // fCasc_ITSTOFtracks = 0;
            // if ((fCasc_NegTrackStatus & AliESDtrack::kITSrefit) || (nTrackCasc->GetTOFBunchCrossing(lMagField) > -95.))
            //     fCasc_ITSTOFtracks++;
            // if ((fCasc_PosTrackStatus & AliESDtrack::kITSrefit) || (pTrackCasc->GetTOFBunchCrossing(lMagField) > -95.))
            //     fCasc_ITSTOFtracks++;
            // if ((fCasc_BacTrackStatus & AliESDtrack::kITSrefit) || (bTrackCasc->GetTOFBunchCrossing(lMagField) > -95.))
            //     fCasc_ITSTOFtracks++;

            // check if at least one of candidate's daughter has a hit in the TOF or has ITSrefit flag (removes Out Of Bunch Pileup)
            fCasc_ITSTOFtracks = ((fCasc_NegTrackStatus & AliESDtrack::kITSrefit) || (nTrackCasc->GetTOFBunchCrossing(lMagField) > -95.)) ? 1 : 0;
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

            // fHistos_eve->FillTH1("hCheckCuts", cutValEv[kDeviationPropLifetime]);
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

            if (fisMC)
            {
                pdgPosDaught = ((AliAODMCParticle *)lMCev->GetTrack((int)TMath::Abs(pTrackCasc->GetLabel())))->PdgCode();
                pdgNegDaught = ((AliAODMCParticle *)lMCev->GetTrack((int)TMath::Abs(nTrackCasc->GetLabel())))->PdgCode();
                pdgBachelor = ((AliAODMCParticle *)lMCev->GetTrack((int)TMath::Abs(bTrackCasc->GetLabel())))->PdgCode();
                int labMothPosDaught = ((AliAODMCParticle *)lMCev->GetTrack((int)TMath::Abs(pTrackCasc->GetLabel())))->GetMother();
                int labMothNegDaught = ((AliAODMCParticle *)lMCev->GetTrack((int)TMath::Abs(nTrackCasc->GetLabel())))->GetMother();
                int labMothV0 = ((AliAODMCParticle *)lMCev->GetTrack((int)TMath::Abs(labMothPosDaught)))->GetMother();
                int labMothBach = ((AliAODMCParticle *)lMCev->GetTrack((int)TMath::Abs(bTrackCasc->GetLabel())))->GetMother();
                pdgV0 = ((AliAODMCParticle *)lMCev->GetTrack((int)TMath::Abs(labMothPosDaught)))->PdgCode();
                pdgCasc = ((AliAODMCParticle *)lMCev->GetTrack((int)TMath::Abs(labMothBach)))->PdgCode();
                physprim = ((AliMCParticle *)lMCev->GetTrack((int)TMath::Abs(labMothBach)))->IsPhysicalPrimary();
                if (fisMCassoc && (pdgPosDaught != 2212 || pdgNegDaught != -211 || pdgBachelor != -211 || labMothPosDaught != labMothNegDaught || pdgV0 != 3122 || labMothV0 != labMothBach || pdgCasc != 3312))
                    assFlag[kxim] = kFALSE;
                if (fisMCassoc && (pdgPosDaught != 211 || pdgNegDaught != -2212 || pdgBachelor != 211 || labMothPosDaught != labMothNegDaught || pdgV0 != -3122 || labMothV0 != labMothBach || pdgCasc != -3312))
                    assFlag[kxip] = kFALSE;
                if (fisMCassoc && (pdgPosDaught != 2212 || pdgNegDaught != -211 || pdgBachelor != -321 || labMothPosDaught != labMothNegDaught || pdgV0 != 3122 || labMothV0 != labMothBach || pdgCasc != 3334))
                    assFlag[komm] = kFALSE;
                if (fisMCassoc && (pdgPosDaught != 211 || pdgNegDaught != -2212 || pdgBachelor != 321 || labMothPosDaught != labMothNegDaught || pdgV0 != -3122 || labMothV0 != labMothBach || pdgCasc != -3334))
                    assFlag[komp] = kFALSE;
            }

            /// Apply parametric BacBarCosPA cut, if requested
            if (fisParametricBacBarCosPA)
            {
                if (fCasc_Pt >= fHist_PtBacBarCosPA->GetXaxis()->GetXmin() && fCasc_Pt <= fHist_PtBacBarCosPA->GetXaxis()->GetXmax() && lPercentile < fCentLimit_BacBarCosPA)
                {
                    SetEvCutVal(kFALSE, kTRUE, kBacBarCosPa, fHist_PtBacBarCosPA->GetBinContent(fHist_PtBacBarCosPA->GetXaxis()->FindBin(fCasc_Pt)));
                }
                else
                {
                    if ((fParticleAnalysisStatus[kxip] || fParticleAnalysisStatus[komp]))
                        SetEvCutVal(kFALSE, kTRUE, kBacBarCosPa, def_cutValEv[kBacBarCosPa]);
                }
            }

            // Apply parametric TrackLengthCut cut, if requested
            if (fisParametricTrackLengthCut)
            {
                if (lPercentile >= fHist_CentTrackLengthCut->GetXaxis()->GetXmin() && lPercentile <= fHist_CentTrackLengthCut->GetXaxis()->GetXmax())
                {
                    SetEvCutVal(kFALSE, kTRUE, kTrackLengthCut, fHist_CentTrackLengthCut->GetBinContent(fHist_CentTrackLengthCut->GetXaxis()->FindBin(lPercentile)));
                }
                else
                {
                    SetEvCutVal(kFALSE, kTRUE, kTrackLengthCut, cutValEv[kTrackLengthCut]);
                }
            }

            // fills TH3 with default cuts
            if (fParticleAnalysisStatus[kxip])
            {
                if (physprim && assFlag[kxim] && ApplyCuts(kxim))
                {
                    fHistos_XiMin->FillTH1("fHistPtXiM", fCasc_Pt);
                    fHistos_XiMin->FillTH1("fCasc_InvMassXiMin", fCasc_InvMassXiMin);
                    fHistos_XiMin->FillTH1("fHistCentXiM", lPercentile);
                    fHistos_XiMin->FillTH3("h3_ptmasscent_def", fCasc_Pt, fCasc_InvMassXiMin, lPercentile);
                }

                if (physprim && assFlag[kxip] && ApplyCuts(kxip))
                {
                    fHistos_XiPlu->FillTH1("fHistPtXiP", fCasc_Pt);
                    fHistos_XiPlu->FillTH1("fCasc_InvMassXiPlu", fCasc_InvMassXiPlu);
                    fHistos_XiPlu->FillTH1("fHistCentXiP", lPercentile);
                    fHistos_XiPlu->FillTH3("h3_ptmasscent_def", fCasc_Pt, fCasc_InvMassXiPlu, lPercentile);
                }
            }
            if (fParticleAnalysisStatus[komp])
            {
                if (physprim && assFlag[komm] && ApplyCuts(komm))
                {
                    fHistos_OmMin->FillTH1("fHistPtOmM", fCasc_Pt);
                    fHistos_OmMin->FillTH1("fCasc_InvMassOmMin", fCasc_InvMassOmMin);
                    fHistos_OmMin->FillTH1("fHistCentOmM", lPercentile);
                    fHistos_OmMin->FillTH3("h3_ptmasscent_def", fCasc_Pt, fCasc_InvMassOmMin, lPercentile);
                }
                if (physprim && assFlag[komp] && ApplyCuts(komp))
                {
                    fHistos_OmPlu->FillTH1("fHistPtOmP", fCasc_Pt);
                    fHistos_OmPlu->FillTH1("fCasc_InvMassOmPlu", fCasc_InvMassOmPlu);
                    fHistos_OmPlu->FillTH1("fHistCentOmP", lPercentile);
                    fHistos_OmPlu->FillTH3("h3_ptmasscent_def", fCasc_Pt, fCasc_InvMassOmPlu, lPercentile);
                }
            }

            /// filling 3D histograms
            // if (!fDefOnly)
            // FillHistCutVariations(kTRUE, lPercentile, 1, assFlag, 0);

            // filling 3D histograms
            if (!fDefOnly)
                FillHistCutVariations(lPercentile, physprim, assFlag);
        }
    }
    DataPosting();
}

Float_t AliAnalysisTaskStrangeCascadesRun2::GetLengthInActiveZone(AliAODTrack *gt, Float_t deltaY, Float_t deltaZ, Float_t b)
{
    // Input parameters:
    //   deltaY - user defined "dead region" in cm
    //   deltaZ - user defined "active region" in cm (250 cm drift lenght - 14 cm L1 delay
    //   b     - magnetic field
    AliESDtrack esdTrack(gt);
    esdTrack.SetESDEvent((AliESDEvent *)gt->GetEvent());
    AliExternalTrackParam etp;
    etp.CopyFromVTrack(gt);
    esdTrack.ResetTrackParamIp(&etp);
    return esdTrack.GetLengthInActiveZone(1, deltaY, deltaZ, b);
}

void AliAnalysisTaskStrangeCascadesRun2::SetCutValue(bool isTopo, int cutName, double cutVal, int particle = -1, int ptInterval = -1)
{
    if (!isTopo)
        cutValEv[cutName] = cutVal;

    else
    {
        if (particle == -1)
        {
            for (int iPart = kXi; iPart < knumpart; iPart++)
            {
                if (ptInterval == -1)
                {
                    for (int iPt = kLow; iPt < kNumPtInterval; iPt++)
                    {
                        cutValTopo[iPart][cutName][iPt] = cutVal;
                    }
                }
                else
                {
                    cutValTopo[iPart][cutName][ptInterval] = cutVal;
                }
            }
        }
        else if (ptInterval == -1)
        {
            for (int iPt = kLow; iPt < kNumPtInterval; iPt++)
            {
                cutValTopo[particle][cutName][iPt] = cutVal;
            }
        }
        else
        {
            cutValTopo[particle][cutName][ptInterval] = cutVal;
        }
    }
}

void AliAnalysisTaskStrangeCascadesRun2::SetDefCutValue(bool isTopo, int cutName, double cutVal, int particle = -1, int ptInterval = -1)
{
    if (!isTopo)
    {
        // SetCutValue(def_cutValEv, isTopo, cutName, cutVal, particle, ptInterval);
        // SetCutValue(cutValEv, isTopo, cutName, cutVal, particle, ptInterval);

        def_cutValEv[cutName] = cutVal;
        cutValEv[cutName] = cutVal;
    }
    else
    {

        // SetCutValue(def_cutValTopo, isTopo, cutName, cutVal, particle, ptInterval);
        // SetCutValue(cutValTopo, isTopo, cutName, cutVal, particle, ptInterval);
        if (particle == -1)
        {
            for (int iPart = kXi; iPart < knumpart; iPart++)
            {
                if (ptInterval == -1)
                {
                    for (int iPt = kLow; iPt < kNumPtInterval; iPt++)
                    {
                        def_cutValTopo[iPart][cutName][iPt] = cutVal;
                        cutValTopo[iPart][cutName][iPt] = cutVal;
                    }
                }
                else
                {
                    def_cutValTopo[iPart][cutName][ptInterval] = cutVal;
                    cutValTopo[iPart][cutName][ptInterval] = cutVal;
                }
            }
        }
        else if (ptInterval == -1)
        {
            for (int iPt = kLow; iPt < kNumPtInterval; iPt++)
            {
                def_cutValTopo[particle][cutName][iPt] = cutVal;
                cutValTopo[particle][cutName][iPt] = cutVal;
            }
        }
        else
        {
            def_cutValTopo[particle][cutName][ptInterval] = cutVal;
            cutValTopo[particle][cutName][ptInterval] = cutVal;
        }
    }
}

void AliAnalysisTaskStrangeCascadesRun2::SetVarCutValue(bool isTopo, int cutName, int cutVar, double cutVal, int particle = -1, int ptInterval = -1)
{
    if (!isTopo)
    {
        var_cutValEv[cutName][cutVar] = cutVal;
    }
    else
    {
        if (particle == -1)
        {
            for (int iPart = kXi; iPart < knumpart; iPart++)
            {
                if (ptInterval == -1)
                {
                    for (int iPt = kLow; iPt < kNumPtInterval; iPt++)
                    {
                        var_cutValTopo[iPart][cutName][iPt][cutVar] = cutVal;
                    }
                }
                else
                {
                    var_cutValTopo[iPart][cutName][ptInterval][cutVar] = cutVal;
                }
            }
        }
        else if (ptInterval == -1)
        {
            for (int iPt = kLow; iPt < kNumPtInterval; iPt++)
            {
                var_cutValTopo[particle][cutName][iPt][cutVar] = cutVal;
            }
        }
        else
        {
            var_cutValTopo[particle][cutName][ptInterval][cutVar] = cutVal;
        }
    }
}

///________________________________________________________________________
void AliAnalysisTaskStrangeCascadesRun2::SetParametricBacBarCosPA(int nbins, float *ptbins, float *values, int cent_limit)
{
    fisParametricBacBarCosPA = kTRUE;
    fCentLimit_BacBarCosPA = cent_limit;
    fHist_PtBacBarCosPA = new TH1F("", "", nbins, ptbins);
    for (int iBin = 1; iBin <= nbins; iBin++)
        fHist_PtBacBarCosPA->SetBinContent(iBin, values[iBin - 1]);
}

///________________________________________________________________________
void AliAnalysisTaskStrangeCascadesRun2::SetDefCuts()
{
    Int_t all = -1;

    /// event selection cuts
    SetDefCutValue(kFALSE, kRapidityIntervalMin, -0.5);    // kCasc_y
    SetDefCutValue(kFALSE, kRapidityIntervalMax, 0.0);     // == 0
    SetDefCutValue(kFALSE, kTpcDedxPidSigma, 4);           // kCasc_NSigPID
    SetDefCutValue(kFALSE, kDeviationPropLifetime, 3.0);   // kCasc_PropLifetXi
    SetDefCutValue(kFALSE, kLeastTpcClusters, 70);         // kCasc_LeastTPCcls,    // --> DON'T USE
    SetDefCutValue(kFALSE, kCompetingCascRejectOm, 0.008); // kCasc_CompetingXiMass
    SetDefCutValue(kFALSE, kLeastCRows, 70);               // not used --> USE -> 70 (emily), 80 (marek's suggestion - from Michal)
    SetDefCutValue(kFALSE, kLeastCRowsOvF, 0.8);           // not used --> USE -> (0.8)
    SetDefCutValue(kFALSE, kTrackLengthCut, 1.0);          // not used --> USE
    SetDefCutValue(kFALSE, kEtaDaughter, 0.8);             // kCasc_etaDaugh
    SetDefCutValue(kFALSE, kBacBarCosPa, 1.0);             // kCasc_BacBarCosPA :: not used --> USE
    SetDefCutValue(kFALSE, kV0InvMassWindow, 0.008);       // kCasc_InvMassLam,     // set to 0.008 instead of 0.005 (lambda == v0)

    /// topological cuts
    SetDefCutValue(kTRUE, kCascTransDecayRadius, 0.6, kXi, all); // kCasc_CascRad(XiMid)
    SetDefCutValue(kTRUE, kCascTransDecayRadius, 0.5, kOm, all); // kCasc_CascRad(XiMid)
    SetDefCutValue(kTRUE, kV0TransDecayRadius, 1.2, kXi, all);   // kCasc_V0RadXi(Mid)
    SetDefCutValue(kTRUE, kV0TransDecayRadius, 1.1, kOm, all);   // kCasc_V0RadXi(Mid)
    SetDefCutValue(kTRUE, kDcaBachToPv, 0.04, all, all);         // kCasc_DcaBachToPV
    SetDefCutValue(kTRUE, kDcaV0ToPv, 0.06, all, all);           // kCasc_DcaV0ToPV
    SetDefCutValue(kTRUE, kDcaMesV0ToPv, 0.04, all, all);        // kCasc_DcaMesToPV
    SetDefCutValue(kTRUE, kDcaBarV0ToPv, 0.03, all, all);        // kCasc_DcaBarToPV
    SetDefCutValue(kTRUE, kDcaV0Daughters, 1.7, all, kLow);      // kCasc_DcaV0Daught(Mid), // DCA V0 daughters (sigma)
    SetDefCutValue(kTRUE, kDcaV0Daughters, 1.6, all, kMid);
    SetDefCutValue(kTRUE, kDcaV0Daughters, 1.5, all, kHigh);
    SetDefCutValue(kTRUE, kDcaBachToV0, 1.6, all, kLow); // kCasc_DcaCascDaught(Mid)
    SetDefCutValue(kTRUE, kDcaBachToV0, 1.5, all, kMid);
    SetDefCutValue(kTRUE, kDcaBachToV0, 1.4, all, kHigh);
    SetDefCutValue(kTRUE, kCascCosPa, 0.96, all, kLow); // kCasc_CascCosPA(Low)
    SetDefCutValue(kTRUE, kCascCosPa, 0.97, all, kMid);
    SetDefCutValue(kTRUE, kCascCosPa, 0.97, all, kHigh);
    SetDefCutValue(kTRUE, kV0CosPa, 0.98, kXi, all);
    SetDefCutValue(kTRUE, kV0CosPa, 0.97, kOm, kLow);
    SetDefCutValue(kTRUE, kV0CosPa, 0.98, kOm, kMid); // kCasc_V0CosPAXi(Mid)
    SetDefCutValue(kTRUE, kV0CosPa, 0.98, kOm, kHigh);
}

///________________________________________________________________________
void AliAnalysisTaskStrangeCascadesRun2::SetDefCutVariations()
// {

//     /// Initialise all elements of var array to -1 to skip unused indices
//     // CAUTION: make sure none of the var values are -1: this value is ignored for applying cuts
//     std::fill_n(*var_cutValEv, kNumCascEvCuts * kNumCutVars, -1.);
//     std::fill_n(***var_cutValTopo, knumpart * kNumCascTopoCuts * kNumPtInterval * kNumCutVars, -1.);

//     Int_t all = -1;

//     /// event selection cuts
//     // SetVarCutValue(kFALSE, kRapidityIntervalMin, -0.5);    // kCasc_y
//     // SetVarCutValue(kFALSE, kRapidityIntervalMax, 0.0);     // == 0
//     SetVarCutValue(kFALSE, kTpcDedxPidSigma, kLoose, 5);         // kCasc_NSigPID
//     SetVarCutValue(kFALSE, kTpcDedxPidSigma, kTight, 3.5);       // kCasc_NSigPID
//     SetVarCutValue(kFALSE, kDeviationPropLifetime, kLoose, 4);   // kCasc_PropLifetXi
//     SetVarCutValue(kFALSE, kDeviationPropLifetime, kTight, 2);   // kCasc_PropLifetXi
//     SetVarCutValue(kFALSE, kLeastTpcClusters, kTight, 75);       // kCasc_LeastTPCcls,    // --> DON'T USE
//     SetVarCutValue(kFALSE, kLeastTpcClusters, kVeryTight, 80);   // kCasc_LeastTPCcls,    // --> DON'T USE
//     SetVarCutValue(kFALSE, kCompetingCascRejectOm, kLoose, 0.0); // kCasc_CompetingXiMass

//     // SetVarCutValue(kFALSE, kLeastCRows, 70);               // not used --> USE -> 70 (emily), 80 (marek's suggestion - from Michal)
//     // SetVarCutValue(kFALSE, kLeastCRowsOvF, 0.8);           // not used --> USE -> (0.8)
//     // SetVarCutValue(kFALSE, kTrackLengthCut, 1.0);          // not used --> USE
//     // SetVarCutValue(kFALSE, kEtaDaughter, 0.8);             // kCasc_etaDaugh
//     // SetVarCutValue(kFALSE, kBacBarCosPa, 1.0);             // kCasc_BacBarCosPA :: not used --> USE
//     SetVarCutValue(kFALSE, kV0InvMassWindow, kLoose, 0.009); // kCasc_InvMassLam,     // set to 0.008 instead of 0.005 (lambda == v0)
//     SetVarCutValue(kFALSE, kV0InvMassWindow, kTight, 0.006); // kCasc_InvMassLam,     // set to 0.008 instead of 0.005 (lambda == v0)

//     /// topological cuts
//     SetVarCutValue(kTRUE, kV0CosPa, kVeryLoose, 0.95, all, all);
//     SetVarCutValue(kTRUE, kV0CosPa, kLoose, 0.96, all, all);
//     SetVarCutValue(kTRUE, kV0CosPa, kTight, 0.985, kXi, kLow);
//     SetVarCutValue(kTRUE, kV0CosPa, kTight, 0.996, kXi, kMid);
//     SetVarCutValue(kTRUE, kV0CosPa, kTight, 0.998, kXi, kHigh);
//     SetVarCutValue(kTRUE, kV0CosPa, kTight, 0.975, kOm, kLow);
//     SetVarCutValue(kTRUE, kV0CosPa, kTight, 0.994, kOm, kMid);
//     SetVarCutValue(kTRUE, kV0CosPa, kTight, 0.996, kOm, kHigh);
//     SetVarCutValue(kTRUE, kV0CosPa, kVeryTight, 0.99, kXi, kLow);
//     SetVarCutValue(kTRUE, kV0CosPa, kVeryTight, 0.998, kXi, kMid);
//     SetVarCutValue(kTRUE, kV0CosPa, kVeryTight, 0.999, kXi, kHigh);
//     SetVarCutValue(kTRUE, kV0CosPa, kVeryTight, 0.98, kOm, kLow);
//     SetVarCutValue(kTRUE, kV0CosPa, kVeryTight, 0.994, kOm, kMid);
//     SetVarCutValue(kTRUE, kV0CosPa, kVeryTight, 0.996, kOm, kHigh);

//     SetVarCutValue(kTRUE, kCascCosPa, kVeryLoose, 0.95, all, all);
//     SetVarCutValue(kTRUE, kCascCosPa, kLoose, 0.95, all, kLow);
//     SetVarCutValue(kTRUE, kCascCosPa, kLoose, 0.96, all, kMid);
//     SetVarCutValue(kTRUE, kCascCosPa, kLoose, 0.96, all, kHigh);
//     SetVarCutValue(kTRUE, kCascCosPa, kTight, 0.965, all, kLow);
//     SetVarCutValue(kTRUE, kCascCosPa, kTight, 0.985, all, kMid);
//     SetVarCutValue(kTRUE, kCascCosPa, kTight, 0.985, kXi, kHigh);
//     SetVarCutValue(kTRUE, kCascCosPa, kTight, 0.98, kOm, kHigh);
//     SetVarCutValue(kTRUE, kCascCosPa, kVeryTight, 0.97, all, kLow);
//     SetVarCutValue(kTRUE, kCascCosPa, kVeryTight, 0.992, all, kMid);
//     SetVarCutValue(kTRUE, kCascCosPa, kVeryTight, 0.992, kXi, kHigh);
//     SetVarCutValue(kTRUE, kCascCosPa, kVeryTight, 0.985, kOm, kHigh);

//     SetVarCutValue(kTRUE, kDcaBachToV0, kVeryLoose, 2.0, all, all);
//     SetVarCutValue(kTRUE, kDcaBachToV0, kLoose, 1.8, all, all);
//     SetVarCutValue(kTRUE, kDcaBachToV0, kTight, 1.4, all, kLow);
//     SetVarCutValue(kTRUE, kDcaBachToV0, kTight, 1.0, all, kMid);
//     SetVarCutValue(kTRUE, kDcaBachToV0, kTight, 0.9, kXi, kHigh);
//     SetVarCutValue(kTRUE, kDcaBachToV0, kTight, 1.0, kOm, kHigh);
//     SetVarCutValue(kTRUE, kDcaBachToV0, kVeryTight, 1.2, all, kLow);
//     SetVarCutValue(kTRUE, kDcaBachToV0, kVeryTight, 0.8, all, kMid);
//     SetVarCutValue(kTRUE, kDcaBachToV0, kVeryTight, 0.7, kXi, kHigh);
//     SetVarCutValue(kTRUE, kDcaBachToV0, kVeryTight, 0.8, kOm, kHigh);

//     SetVarCutValue(kTRUE, kDcaV0Daughters, kVeryLoose, 2.0, all, all);
//     SetVarCutValue(kTRUE, kDcaV0Daughters, kLoose, 1.9, all, kLow);
//     SetVarCutValue(kTRUE, kDcaV0Daughters, kLoose, 1.8, all, kMid);
//     SetVarCutValue(kTRUE, kDcaV0Daughters, kLoose, 1.8, all, kHigh);
//     SetVarCutValue(kTRUE, kDcaV0Daughters, kTight, 1.5, all, kLow);
//     SetVarCutValue(kTRUE, kDcaV0Daughters, kTight, 0.9, kXi, kMid);
//     SetVarCutValue(kTRUE, kDcaV0Daughters, kTight, 0.6, kXi, kHigh);
//     SetVarCutValue(kTRUE, kDcaV0Daughters, kTight, 1.4, kOm, kMid);
//     SetVarCutValue(kTRUE, kDcaV0Daughters, kTight, 1.4, kOm, kHigh);
//     SetVarCutValue(kTRUE, kDcaV0Daughters, kVeryTight, 1.3, all, kLow);
//     SetVarCutValue(kTRUE, kDcaV0Daughters, kVeryTight, 0.8, kXi, kMid);
//     SetVarCutValue(kTRUE, kDcaV0Daughters, kVeryTight, 0.4, kXi, kHigh);
//     SetVarCutValue(kTRUE, kDcaV0Daughters, kVeryTight, 1.2, kOm, kMid);
//     SetVarCutValue(kTRUE, kDcaV0Daughters, kVeryTight, 1.2, kOm, kHigh);

//     SetVarCutValue(kTRUE, kCascTransDecayRadius, kLoose, 0.5, all, all);
//     SetVarCutValue(kTRUE, kCascTransDecayRadius, kTight, 1.0, kXi, kLow);
//     SetVarCutValue(kTRUE, kCascTransDecayRadius, kTight, 1.4, kXi, kMid);
//     SetVarCutValue(kTRUE, kCascTransDecayRadius, kTight, 2.0, kXi, kHigh);
//     SetVarCutValue(kTRUE, kCascTransDecayRadius, kTight, 0.6, kOm, kLow);
//     SetVarCutValue(kTRUE, kCascTransDecayRadius, kTight, 0.9, kOm, kMid);
//     SetVarCutValue(kTRUE, kCascTransDecayRadius, kTight, 1.0, kOm, kHigh);
//     SetVarCutValue(kTRUE, kCascTransDecayRadius, kVeryTight, 1.4, kXi, kLow);
//     SetVarCutValue(kTRUE, kCascTransDecayRadius, kVeryTight, 2, kXi, kMid);
//     SetVarCutValue(kTRUE, kCascTransDecayRadius, kVeryTight, 3, kXi, kHigh);

//     SetVarCutValue(kTRUE, kV0TransDecayRadius, kLoose, 1.1, all, all);
//     SetVarCutValue(kTRUE, kV0TransDecayRadius, kTight, 3.0, kXi, kLow);
//     SetVarCutValue(kTRUE, kV0TransDecayRadius, kTight, 5.5, kXi, kMid);
//     SetVarCutValue(kTRUE, kV0TransDecayRadius, kTight, 7.0, kXi, kHigh);
//     SetVarCutValue(kTRUE, kV0TransDecayRadius, kTight, 2.0, kOm, kLow);
//     SetVarCutValue(kTRUE, kV0TransDecayRadius, kTight, 2.5, kOm, kMid);
//     SetVarCutValue(kTRUE, kV0TransDecayRadius, kTight, 5.0, kOm, kHigh);

//     SetVarCutValue(kTRUE, kDcaBachToPv, kLoose, 0.03, all, all);
//     SetVarCutValue(kTRUE, kDcaBachToPv, kTight, 0.11, kXi, all);
//     SetVarCutValue(kTRUE, kDcaBachToPv, kTight, 0.06, kOm, all);
//     SetVarCutValue(kTRUE, kDcaBachToPv, kVeryTight, 0.19, kXi, all);
//     SetVarCutValue(kTRUE, kDcaBachToPv, kVeryTight, 0.09, kOm, all);

//     SetVarCutValue(kTRUE, kDcaV0ToPv, kLoose, 0.05, all, all);
//     SetVarCutValue(kTRUE, kDcaV0ToPv, kTight, 0.09, kXi, all);
//     SetVarCutValue(kTRUE, kDcaV0ToPv, kTight, 0.08, kOm, all);
//     SetVarCutValue(kTRUE, kDcaV0ToPv, kVeryTight, 0.14, kXi, kLow);
//     SetVarCutValue(kTRUE, kDcaV0ToPv, kVeryTight, 0.13, kXi, kMid);
//     SetVarCutValue(kTRUE, kDcaV0ToPv, kVeryTight, 0.11, kXi, kHigh);
//     SetVarCutValue(kTRUE, kDcaV0ToPv, kVeryTight, 0.1, kOm, all);

//     SetVarCutValue(kTRUE, kDcaMesV0ToPv, kVeryLoose, 0.02, all, all);
//     SetVarCutValue(kTRUE, kDcaMesV0ToPv, kLoose, 0.03, all, all);
//     SetVarCutValue(kTRUE, kDcaMesV0ToPv, kTight, 0.24, kXi, all);
//     SetVarCutValue(kTRUE, kDcaMesV0ToPv, kTight, 0.11, kOm, kLow);
//     SetVarCutValue(kTRUE, kDcaMesV0ToPv, kTight, 0.06, kOm, kMid);
//     SetVarCutValue(kTRUE, kDcaMesV0ToPv, kTight, 0.15, kOm, kHigh);
//     SetVarCutValue(kTRUE, kDcaMesV0ToPv, kVeryTight, 0.42, kXi, kLow);
//     SetVarCutValue(kTRUE, kDcaMesV0ToPv, kVeryTight, 0.5, kXi, kMid);
//     SetVarCutValue(kTRUE, kDcaMesV0ToPv, kVeryTight, 0.5, kXi, kHigh);

//     SetVarCutValue(kTRUE, kDcaBarV0ToPv, kLoose, 0.02, all, all);
//     SetVarCutValue(kTRUE, kDcaBarV0ToPv, kTight, 0.07, kXi, all);
//     SetVarCutValue(kTRUE, kDcaBarV0ToPv, kTight, 0.06, kOm, all);
//     SetVarCutValue(kTRUE, kDcaBarV0ToPv, kVeryTight, 0.13, kXi, kLow);
//     SetVarCutValue(kTRUE, kDcaBarV0ToPv, kVeryTight, 0.11, kXi, kMid);
//     SetVarCutValue(kTRUE, kDcaBarV0ToPv, kVeryTight, 0.11, kXi, kHigh);
//     SetVarCutValue(kTRUE, kDcaBarV0ToPv, kVeryTight, 0.1, kOm, all);
// }
{
    /// Initialise all elements of nvarcut array to -1 to skip unused indices
    std::fill_n(nvarcut_Ev, kNumCascEvCuts, -1);
    std::fill_n(nvarcut_Topo, kNumCascTopoCuts, -1);

    /// event selection cuts
    // SetVarCutValue(kFALSE, kRapidityIntervalMin, -0.5);    // kCasc_y
    // SetVarCutValue(kFALSE, kRapidityIntervalMax, 0.0);     // == 0
    SetCutVariation(kFALSE, kTpcDedxPidSigma, 6, 2, 7);             // kCasc_NSigPID
    SetCutVariation(kFALSE, kDeviationPropLifetime, 7, 2, 5);       // kCasc_PropLifetXi
    SetCutVariation(kFALSE, kLeastTpcClusters, 5, 65, 85);          // kCasc_LeastTPCcls,    // --> DON'T USE
    SetCutVariation(kFALSE, kCompetingCascRejectOm, 6, 0.0, 0.010); // kCasc_CompetingXiMass

    // SetVarCutValue(kFALSE, kLeastCRows, 70);               // not used --> USE -> 70 (emily), 80 (marek's suggestion - from Michal)
    // SetVarCutValue(kFALSE, kLeastCRowsOvF, 0.8);           // not used --> USE -> (0.8)
    // SetVarCutValue(kFALSE, kTrackLengthCut, 1.0);          // not used --> USE
    // SetVarCutValue(kFALSE, kEtaDaughter, 0.8);             // kCasc_etaDaugh
    // SetVarCutValue(kFALSE, kBacBarCosPa, 1.0);             // kCasc_BacBarCosPA :: not used --> USE
    SetCutVariation(kFALSE, kV0InvMassWindow, 4, 0.003, 0.012); // kCasc_InvMassLam,     // set to 0.008 instead of 0.005 (lambda == v0)

    /// topological cuts
    SetCutVariation(kTRUE, kV0CosPa, 21, 0.94, 0.999);
    SetCutVariation(kTRUE, kCascCosPa, 21, 0.94, 0.999);
    SetCutVariation(kTRUE, kDcaBachToV0, 19, 0.3, 2.1);
    SetCutVariation(kTRUE, kDcaV0Daughters, 19, 0.3, 2.1);
    SetCutVariation(kTRUE, kCascTransDecayRadius, 11, 0.5, 1.5);
    SetCutVariation(kTRUE, kV0TransDecayRadius, 13, 1., 7.);
    SetCutVariation(kTRUE, kDcaBachToPv, 19, 0.02, 0.2);
    SetCutVariation(kTRUE, kDcaV0ToPv, 19, 0.02, 0.2);
    SetCutVariation(kTRUE, kDcaMesV0ToPv, 24, 0.02, 0.5);
    SetCutVariation(kTRUE, kDcaBarV0ToPv, 15, 0.01, 0.15);
}

//________________________________________________________________________
void AliAnalysisTaskStrangeCascadesRun2::SetCutVariation(bool isTopo, int cutnum, int nvar, double lowval, double highval)
{
    if (!isTopo)
    {
        nvarcut_Ev[cutnum] = nvar;
        varlowcut_Ev[cutnum] = lowval;
        varhighcut_Ev[cutnum] = highval;
    }
    else
    {
        nvarcut_Topo[cutnum] = nvar;
        varlowcut_Topo[cutnum] = lowval;
        varhighcut_Topo[cutnum] = highval;
    }
}

///________________________________________________________________________
void AliAnalysisTaskStrangeCascadesRun2::SetDefOnly(bool isdefonly)
{
    fDefOnly = isdefonly;
}

///________________________________________________________________________
bool AliAnalysisTaskStrangeCascadesRun2::ApplyCuts(int part)
{
    /// Event Selection Cuts

    /// check candidate's charge
    if ((part == kxim || part == komm) && fCasc_charge > 0)
        return kFALSE;

    if ((part == kxip || part == komp) && fCasc_charge < 0)
        return kFALSE;

    /// check DCA bachelor-baryon. If it is too small --> bump structure in Inv Mass
    if (fCasc_BacBarCosPA > cutValEv[kBacBarCosPa])
        return kFALSE;
    else
        fHistos_eve->FillTH1("fHistCutsEv", kBacBarCosPa);

    /// check PID for all daughters (particle hypothesis' dependent)
    if ((part == kxip) && (TMath::Abs(fCasc_NSigPosPion) > cutValEv[kTpcDedxPidSigma] || TMath::Abs(fCasc_NSigNegProton) > cutValEv[kTpcDedxPidSigma] || TMath::Abs(fCasc_NSigBacPion) > cutValEv[kTpcDedxPidSigma]))
        return kFALSE;
    else
        fHistos_eve->FillTH1("fHistCutsEv", kTpcDedxPidSigma);
    if ((part == kxim) && (TMath::Abs(fCasc_NSigNegPion) > cutValEv[kTpcDedxPidSigma] || TMath::Abs(fCasc_NSigPosProton) > cutValEv[kTpcDedxPidSigma] || TMath::Abs(fCasc_NSigBacPion) > cutValEv[kTpcDedxPidSigma]))
        return kFALSE;
    else
        fHistos_eve->FillTH1("fHistCutsEv", kTpcDedxPidSigma);
    if ((part == komp) && (TMath::Abs(fCasc_NSigPosPion) > cutValEv[kTpcDedxPidSigma] || TMath::Abs(fCasc_NSigNegProton) > cutValEv[kTpcDedxPidSigma] || TMath::Abs(fCasc_NSigBacKaon) > cutValEv[kTpcDedxPidSigma]))
        return kFALSE;
    else
        fHistos_eve->FillTH1("fHistCutsEv", kTpcDedxPidSigma);
    if ((part == komm) && (TMath::Abs(fCasc_NSigNegPion) > cutValEv[kTpcDedxPidSigma] || TMath::Abs(fCasc_NSigPosProton) > cutValEv[kTpcDedxPidSigma] || TMath::Abs(fCasc_NSigBacKaon) > cutValEv[kTpcDedxPidSigma]))
        return kFALSE;
    else
        fHistos_eve->FillTH1("fHistCutsEv", kTpcDedxPidSigma);

    /// check candidate daughters' pseudo-rapidity
    if (TMath::Abs(fCasc_etaPos) > cutValEv[kEtaDaughter] || TMath::Abs(fCasc_etaNeg) > cutValEv[kEtaDaughter] || TMath::Abs(fCasc_etaBac) > cutValEv[kEtaDaughter])
        return kFALSE;
    else
        fHistos_eve->FillTH1("fHistCutsEv", kEtaDaughter);

    /// check candidate daughters' crossed TPC Rows (note that the checked value is the lowest among the daughters)
    if (fCasc_LeastCRows < cutValEv[kLeastCRows])
        return kFALSE;
    else
        fHistos_eve->FillTH1("fHistCutsEv", kLeastCRows);
    /// check candidate daughters' crossed TPC Rows over findable
    if (fCasc_LeastCRowsOvF < cutValEv[kLeastCRowsOvF])
        return kFALSE;
    else
        fHistos_eve->FillTH1("fHistCutsEv", kLeastCRowsOvF);

    /// check candidate daughters' TPC clusters
    // if (fCasc_LeastTPCcls < cutValEv[kCasc_LeastTPCcls])
    //     return kFALSE;

    // check candidate daughters' TPC clusters or/and apply geometrical cut
    if (fCasc_TrackLengthCut < (cutValEv[kTrackLengthCut] - 0.1))
        return kFALSE;
    if (TMath::Abs(cutValEv[kTrackLengthCut] - 1) < 0.1 && fCasc_TrackLengthCut == 2)
        return kFALSE;
    else
        fHistos_eve->FillTH1("fHistCutsEv", kTrackLengthCut);

    /// check candidate V0 daughter's mass difference from nominal Lambda mass
    if (TMath::Abs(fCasc_InvMassLam - 1.115683) > cutValEv[kV0InvMassWindow])
        return kFALSE;
    else
        fHistos_eve->FillTH1("fHistCutsEv", kV0InvMassWindow);
    // TPC refit, should be already verified for Offline V0s
    if (!(fCasc_PosTrackStatus & AliESDtrack::kTPCrefit) ||
        !(fCasc_NegTrackStatus & AliESDtrack::kTPCrefit) ||
        !(fCasc_BacTrackStatus & AliESDtrack::kTPCrefit))
        return kFALSE;
    else
        fHistos_eve->FillTH1("fHistCutsEv", 12.5);
    // check candidate daughters' Chi^2 per TPC cluster
    // if (fCasc_MaxChi2perCls > cutValEv[kCasc_MaxChi2perCls])
    //     return kFALSE;

    /// Xi: Particle dependent cuts
    if ((part == kxip || part == kxim))
    {

        /// check candidate's rapidity (particle hypothesis' dependent)
        if (fCasc_yXi < cutValEv[kRapidityIntervalMin] && fCasc_yXi > cutValEv[kRapidityIntervalMax])
            return kFALSE;
        else
        {
            fHistos_eve->FillTH1("fHistCutsEv", kRapidityIntervalMin);
            fHistos_eve->FillTH1("fHistCutsEv", kRapidityIntervalMax);
        }

        /// check candidate's proper lifetime (particle hypothesis' dependent). Remember: c*tau = L*m/p
        if (((1.32171 * fCasc_DistOverTotP) > (4.91 * cutValEv[kDeviationPropLifetime])))
            return kFALSE; /// 4.91 is the ctau of xi in cm
        else
            fHistos_eve->FillTH1("fHistCutsEv", kDeviationPropLifetime);

        if (fCasc_Pt > 0. && fCasc_Pt < ptXiBoundary_LowMid) /// low pt range
        {
            /// check candidate's 2D decay distance from PV (if it is too small, then it's not a weak decay)
            if (fCasc_CascRad < cutValTopo[kXi][kCascTransDecayRadius][kLow])
                return kFALSE;
            else
                fHistos_eve->FillTH3("fHistCutsTopo", kXi, kCascTransDecayRadius, kLow);
            /// check candidate V0 daughter's 2D decay distance from PV (if it is too small, then it's not a weak decay)
            if (fCasc_V0Rad < cutValTopo[kXi][kV0TransDecayRadius][kLow])
                return kFALSE;
            else
                fHistos_eve->FillTH3("fHistCutsTopo", kXi, kV0TransDecayRadius, kLow);
            /// check candidate daughters' DCA to Primary Vertex (needs to be large because decay is far from the Primary Vertex)
            if (fCasc_DcaBachToPV < cutValTopo[kXi][kDcaBachToPv][kLow])
                return kFALSE;
            else
                fHistos_eve->FillTH3("fHistCutsTopo", kXi, kDcaBachToPv, kLow);
            if (fCasc_DcaV0ToPV < cutValTopo[kXi][kDcaV0ToPv][kLow])
                return kFALSE;
            else
                fHistos_eve->FillTH3("fHistCutsTopo", kXi, kDcaV0ToPv, kLow);

            /// check V0 daughters' DCA to Primary Vertex. Different cut for meson and baryon daughters, so different conditions for + and - candidates
            if ((part == kxim) && (fCasc_DcaPosToPV < cutValTopo[kXi][kDcaBarV0ToPv][kLow] || fCasc_DcaNegToPV < cutValTopo[kXi][kDcaMesV0ToPv][kLow]))
                return kFALSE; // in this case pos=p, neg=pi-, bach=pi- or K-
            if ((part == kxip) && (fCasc_DcaPosToPV < cutValTopo[kXi][kDcaMesV0ToPv][kLow] || fCasc_DcaNegToPV < cutValTopo[kXi][kDcaBarV0ToPv][kLow]))
                return kFALSE; // in this case pos=pi+, neg=anti-p, bach=pi+ or K+
            else
            {
                fHistos_eve->FillTH3("fHistCutsTopo", kXi, kDcaBarV0ToPv, kLow);
                fHistos_eve->FillTH3("fHistCutsTopo", kXi, kDcaMesV0ToPv, kLow);
            }
            /// check V0 daughter's daughters DCA between them (needs to be small because they have to come from the same secondary vertex)
            if (fCasc_DcaV0Daught > cutValTopo[kXi][kDcaV0Daughters][kLow])
                return kFALSE;
            else
                fHistos_eve->FillTH3("fHistCutsTopo", kXi, kDcaV0Daughters, kLow);
            if (fCasc_DcaCascDaught > cutValTopo[kXi][kDcaBachToV0][kLow]) /// check candidate daughter's DCA between them (needs to be small because they have to come from the same secondary vertex)
                return kFALSE;
            else
                fHistos_eve->FillTH3("fHistCutsTopo", kXi, kDcaBachToV0, kLow);
            /// check the cosine of the Pointing Angle for both cascade and V0 (angle between candidate's momentum and vector connecting Primary and secondary vertices)
            if (fCasc_CascCosPA < cutValTopo[kXi][kCascCosPa][kLow])
                return kFALSE;
            else
                fHistos_eve->FillTH3("fHistCutsTopo", kXi, kCascCosPa, kLow);
            if (fCasc_V0CosPA < cutValTopo[kXi][kV0CosPa][kLow])
                return kFALSE;
            else
                fHistos_eve->FillTH3("fHistCutsTopo", kXi, kV0CosPa, kLow);
        }

        if (fCasc_Pt >= ptXiBoundary_LowMid && fCasc_Pt < ptXiBoundary_MidHigh) /// mid pt range
        {
            /// check candidate's 2D decay distance from PV (if it is too small, then it's not a weak decay)
            if (fCasc_CascRad < cutValTopo[kXi][kCascTransDecayRadius][kMid])
                return kFALSE;

            /// check candidate V0 daughter's 2D decay distance from PV (if it is too small, then it's not a weak decay)
            if (fCasc_V0Rad < cutValTopo[kXi][kV0TransDecayRadius][kMid])
                return kFALSE;

            /// check candidate daughters' DCA to Primary Vertex (needs to be large because decay is far from the Primary Vertex)
            if (fCasc_DcaBachToPV < cutValTopo[kXi][kDcaBachToPv][kMid])
                return kFALSE;

            if (fCasc_DcaV0ToPV < cutValTopo[kXi][kDcaV0ToPv][kMid])
                return kFALSE;

            /// check V0 daughters' DCA to Primary Vertex. Different cut for meson and baryon daughters, so different conditions for + and - candidates
            if ((part == kxim) && (fCasc_DcaPosToPV < cutValTopo[kXi][kDcaBarV0ToPv][kMid] || fCasc_DcaNegToPV < cutValTopo[kXi][kDcaMesV0ToPv][kMid]))
                return kFALSE; // in this case pos=p, neg=pi-, bach=pi- or K-

            if ((part == kxip) && (fCasc_DcaPosToPV < cutValTopo[kXi][kDcaMesV0ToPv][kMid] || fCasc_DcaNegToPV < cutValTopo[kXi][kDcaBarV0ToPv][kMid]))
                return kFALSE; // in this case pos=pi+, neg=anti-p, bach=pi+ or K+

            /// check V0 daughter's daughters DCA between them (needs to be small because they have to come from the same secondary vertex)
            if (fCasc_DcaV0Daught > cutValTopo[kXi][kDcaV0Daughters][kMid])
                return kFALSE;

            if (fCasc_DcaCascDaught > cutValTopo[kXi][kDcaBachToV0][kMid]) /// check candidate daughter's DCA between them (needs to be small because they have to come from the same secondary vertex)
                return kFALSE;

            /// check the cosine of the Pointing Angle for both cascade and V0 (angle between candidate's momentum and vector connecting Primary and secondary vertices)
            if (fCasc_CascCosPA < cutValTopo[kXi][kCascCosPa][kMid])
                return kFALSE;

            if (fCasc_V0CosPA < cutValTopo[kXi][kV0CosPa][kMid])
                return kFALSE;
        }

        if (fCasc_Pt >= ptXiBoundary_MidHigh) /// high pt range
        {                                     // check candidate V0 daughter's 2D decay distance from PV (if it is too small, then it's not a weak decay)

            /// check candidate's 2D decay distance from PV (if it is too small, then it's not a weak decay)
            if (fCasc_CascRad < cutValTopo[kXi][kCascTransDecayRadius][kHigh])
                return kFALSE;

            /// check candidate V0 daughter's 2D decay distance from PV (if it is too small, then it's not a weak decay)
            if (fCasc_V0Rad < cutValTopo[kXi][kV0TransDecayRadius][kHigh])
                return kFALSE;

            /// check candidate daughters' DCA to Primary Vertex (needs to be large because decay is far from the Primary Vertex)
            if (fCasc_DcaBachToPV < cutValTopo[kXi][kDcaBachToPv][kHigh])
                return kFALSE;

            if (fCasc_DcaV0ToPV < cutValTopo[kXi][kDcaV0ToPv][kHigh])
                return kFALSE;

            /// check V0 daughters' DCA to Primary Vertex. Different cut for meson and baryon daughters, so different conditions for + and - candidates
            if ((part == kxim) && (fCasc_DcaPosToPV < cutValTopo[kXi][kDcaBarV0ToPv][kHigh] || fCasc_DcaNegToPV < cutValTopo[kXi][kDcaMesV0ToPv][kHigh]))
                return kFALSE; // in this case pos=p, neg=pi-, bach=pi- or K-

            if ((part == kxip) && (fCasc_DcaPosToPV < cutValTopo[kXi][kDcaMesV0ToPv][kHigh] || fCasc_DcaNegToPV < cutValTopo[kXi][kDcaBarV0ToPv][kHigh]))
                return kFALSE; // in this case pos=pi+, neg=anti-p, bach=pi+ or K+

            /// check V0 daughter's daughters DCA between them (needs to be small because they have to come from the same secondary vertex)
            if (fCasc_DcaV0Daught > cutValTopo[kXi][kDcaV0Daughters][kHigh])
                return kFALSE;

            if (fCasc_DcaCascDaught > cutValTopo[kXi][kDcaBachToV0][kHigh]) /// check candidate daughter's DCA between them (needs to be small because they have to come from the same secondary vertex)
                return kFALSE;

            /// check the cosine of the Pointing Angle for both cascade and V0 (angle between candidate's momentum and vector connecting Primary and secondary vertices)
            if (fCasc_CascCosPA < cutValTopo[kXi][kCascCosPa][kHigh])
                return kFALSE;

            if (fCasc_V0CosPA < cutValTopo[kXi][kV0CosPa][kHigh])
                return kFALSE;
        }
    }

    /// Omega: Particle dependent cuts
    if ((part == komp || part == komm))
    {
        /// check candidate's rapidity (particle hypothesis' dependent)
        if (fCasc_yOm < cutValEv[kRapidityIntervalMin] && fCasc_yOm > cutValEv[kRapidityIntervalMax])
            return kFALSE;

        // /// check candidate's 2D decay distance from PV (if it is too small, then it's not a weak decay)
        // if (fCasc_CascRad < cutValTopo[kCasc_CascRadOm])
        //     return kFALSE;

        if (((1.67245 * fCasc_DistOverTotP) > (2.461 * cutValEv[kDeviationPropLifetime])))
            return kFALSE; /// 2.461 is the ctau of om in cm

        /// competing cascade rejection (only for omegas)
        if ((part == komm) && ((TMath::Abs(fCasc_InvMassXiMin - 1.32171)) < cutValEv[kCompetingCascRejectOm]))
            return kFALSE;
        if ((part == komp) && ((TMath::Abs(fCasc_InvMassXiPlu - 1.32171)) < cutValEv[kCompetingCascRejectOm]))
            return kFALSE;

        if (fCasc_Pt > 0. && fCasc_Pt < ptOmBoundary_LowMid) /// low pt range
        {
            /// check candidate's 2D decay distance from PV (if it is too small, then it's not a weak decay)
            if (fCasc_CascRad < cutValTopo[kOm][kCascTransDecayRadius][kLow])
                return kFALSE;

            /// check candidate V0 daughter's 2D decay distance from PV (if it is too small, then it's not a weak decay)
            if (fCasc_V0Rad < cutValTopo[kOm][kV0TransDecayRadius][kLow])
                return kFALSE;

            /// check candidate daughters' DCA to Primary Vertex (needs to be large because decay is far from the Primary Vertex)
            if (fCasc_DcaBachToPV < cutValTopo[kOm][kDcaBachToPv][kLow])
                return kFALSE;

            if (fCasc_DcaV0ToPV < cutValTopo[kOm][kDcaV0ToPv][kLow])
                return kFALSE;

            /// check V0 daughters' DCA to Primary Vertex. Different cut for meson and baryon daughters, so different conditions for + and - candidates
            if ((part == komm) && (fCasc_DcaPosToPV < cutValTopo[kOm][kDcaBarV0ToPv][kLow] || fCasc_DcaNegToPV < cutValTopo[kOm][kDcaMesV0ToPv][kLow]))
                return kFALSE; // in this case pos=p, neg=pi-, bach=pi- or K-

            if ((part == komp) && (fCasc_DcaPosToPV < cutValTopo[kOm][kDcaMesV0ToPv][kLow] || fCasc_DcaNegToPV < cutValTopo[kOm][kDcaBarV0ToPv][kLow]))
                return kFALSE; // in this case pos=pi+, neg=anti-p, bach=pi+ or K+

            /// check V0 daughter's daughters DCA between them (needs to be small because they have to come from the same secondary vertex)
            if (fCasc_DcaV0Daught > cutValTopo[kOm][kDcaV0Daughters][kLow])
                return kFALSE;

            if (fCasc_DcaCascDaught > cutValTopo[kOm][kDcaBachToV0][kLow]) /// check candidate daughter's DCA between them (needs to be small because they have to come from the same secondary vertex)
                return kFALSE;

            /// check the cosine of the Pointing Angle for both cascade and V0 (angle between candidate's momentum and vector connecting Primary and secondary vertices)
            if (fCasc_CascCosPA < cutValTopo[kOm][kCascCosPa][kLow])
                return kFALSE;

            if (fCasc_V0CosPA < cutValTopo[kOm][kV0CosPa][kLow])
                return kFALSE;
        }

        if (fCasc_Pt >= ptOmBoundary_LowMid && fCasc_Pt < ptOmBoundary_MidHigh) /// mid pt range
        {
            /// check candidate's 2D decay distance from PV (if it is too small, then it's not a weak decay)
            if (fCasc_CascRad < cutValTopo[kOm][kCascTransDecayRadius][kMid])
                return kFALSE;

            /// check candidate V0 daughter's 2D decay distance from PV (if it is too small, then it's not a weak decay)
            if (fCasc_V0Rad < cutValTopo[kOm][kV0TransDecayRadius][kMid])
                return kFALSE;

            /// check candidate daughters' DCA to Primary Vertex (needs to be large because decay is far from the Primary Vertex)
            if (fCasc_DcaBachToPV < cutValTopo[kOm][kDcaBachToPv][kMid])
                return kFALSE;

            if (fCasc_DcaV0ToPV < cutValTopo[kOm][kDcaV0ToPv][kMid])
                return kFALSE;

            /// check V0 daughters' DCA to Primary Vertex. Different cut for meson and baryon daughters, so different conditions for + and - candidates
            if ((part == komm) && (fCasc_DcaPosToPV < cutValTopo[kOm][kDcaBarV0ToPv][kMid] || fCasc_DcaNegToPV < cutValTopo[kOm][kDcaMesV0ToPv][kMid]))
                return kFALSE; // in this case pos=p, neg=pi-, bach=pi- or K-

            if ((part == komp) && (fCasc_DcaPosToPV < cutValTopo[kOm][kDcaMesV0ToPv][kMid] || fCasc_DcaNegToPV < cutValTopo[kOm][kDcaBarV0ToPv][kMid]))
                return kFALSE; // in this case pos=pi+, neg=anti-p, bach=pi+ or K+

            /// check V0 daughter's daughters DCA between them (needs to be small because they have to come from the same secondary vertex)
            if (fCasc_DcaV0Daught > cutValTopo[kOm][kDcaV0Daughters][kMid])
                return kFALSE;

            if (fCasc_DcaCascDaught > cutValTopo[kOm][kDcaBachToV0][kMid]) /// check candidate daughter's DCA between them (needs to be small because they have to come from the same secondary vertex)
                return kFALSE;

            /// check the cosine of the Pointing Angle for both cascade and V0 (angle between candidate's momentum and vector connecting Primary and secondary vertices)
            if (fCasc_CascCosPA < cutValTopo[kOm][kCascCosPa][kMid])
                return kFALSE;

            if (fCasc_V0CosPA < cutValTopo[kOm][kV0CosPa][kMid])
                return kFALSE;
        }

        if (fCasc_Pt >= ptOmBoundary_MidHigh) /// high pt range
        {
            /// check candidate's 2D decay distance from PV (if it is too small, then it's not a weak decay)
            if (fCasc_CascRad < cutValTopo[kOm][kCascTransDecayRadius][kHigh])
                return kFALSE;

            /// check candidate V0 daughter's 2D decay distance from PV (if it is too small, then it's not a weak decay)
            if (fCasc_V0Rad < cutValTopo[kOm][kV0TransDecayRadius][kHigh])
                return kFALSE;

            /// check candidate daughters' DCA to Primary Vertex (needs to be large because decay is far from the Primary Vertex)
            if (fCasc_DcaBachToPV < cutValTopo[kOm][kDcaBachToPv][kHigh])
                return kFALSE;

            if (fCasc_DcaV0ToPV < cutValTopo[kOm][kDcaV0ToPv][kHigh])
                return kFALSE;

            /// check V0 daughters' DCA to Primary Vertex. Different cut for meson and baryon daughters, so different conditions for + and - candidates
            if ((part == komm) && (fCasc_DcaPosToPV < cutValTopo[kOm][kDcaBarV0ToPv][kHigh] || fCasc_DcaNegToPV < cutValTopo[kOm][kDcaMesV0ToPv][kHigh]))
                return kFALSE; // in this case pos=p, neg=pi-, bach=pi- or K-

            if ((part == komp) && (fCasc_DcaPosToPV < cutValTopo[kOm][kDcaMesV0ToPv][kHigh] || fCasc_DcaNegToPV < cutValTopo[kOm][kDcaBarV0ToPv][kHigh]))
                return kFALSE; // in this case pos=pi+, neg=anti-p, bach=pi+ or K+

            /// check V0 daughter's daughters DCA between them (needs to be small because they have to come from the same secondary vertex)
            if (fCasc_DcaV0Daught > cutValTopo[kOm][kDcaV0Daughters][kHigh])
                return kFALSE;

            if (fCasc_DcaCascDaught > cutValTopo[kOm][kDcaBachToV0][kHigh]) /// check candidate daughter's DCA between them (needs to be small because they have to come from the same secondary vertex)
                return kFALSE;

            /// check the cosine of the Pointing Angle for both cascade and V0 (angle between candidate's momentum and vector connecting Primary and secondary vertices)
            if (fCasc_CascCosPA < cutValTopo[kOm][kCascCosPa][kHigh])
                return kFALSE;

            if (fCasc_V0CosPA < cutValTopo[kOm][kV0CosPa][kHigh])
                return kFALSE;
        }
    }

    /// check that none of daughters is a kink
    // if (fCasc_kinkidx > 0)
    //   return kFALSE;
    /// older, new version below - check if at least one of candidate's daughter has a hit in the TOF or has ITSrefit flag (removes Out Of Bunch Pileup)
    // if (fCasc_ITSTOFtracks < cutval_Casc[kCasc_ITSTOFtracks])
    //   return kFALSE;

    // check if at least one of candidate's daughter has a hit in the TOF or has ITSrefit flag (removes Out Of Bunch Pileup)
    // if (fCasc_ITSTOFtracks < cutval_Casc[kCasc_ITSTOFtracks] - 0.1)
    //     return kFALSE;

    AliInfo("All cuts passed!\n");
    return kTRUE; /// survived!
}

void AliAnalysisTaskStrangeCascadesRun2::SetParticleAnalysisStatus(bool xi, bool omega)
{
    fParticleAnalysisStatus[kxip] = xi;
    fParticleAnalysisStatus[kxim] = xi;
    fParticleAnalysisStatus[komp] = omega;
    fParticleAnalysisStatus[komm] = omega;
}

bool AliAnalysisTaskStrangeCascadesRun2::GetParticleAnalysisStatus(int part)
{
    if (part < 0 || part >= ksignednumpart)
    {
        ::Error("AliAnalysisTaskStrangeCascadesRun2::GetParticleAnalysisStatus", "Wrong particle selected: accepted values from 0 to 1");
        return false;
    }
    return fParticleAnalysisStatus[part];
}

void AliAnalysisTaskStrangeCascadesRun2::SetCentbinning(int ipart, int numcentbins, double *centbins)
{
    fncentbins[ipart] = numcentbins;
    for (int i = 0; i < fncentbins[ipart] + 1; i++)
    {
        fcentbinning[ipart][i] = centbins[i];
    }
}

void AliAnalysisTaskStrangeCascadesRun2::SetPtbinning(int ipart, int numptbins, double *ptbins)
{
    fnptbins[ipart] = numptbins;
    for (int i = 0; i < fnptbins[ipart] + 1; i++)
    {
        fptbinning[ipart][i] = ptbins[i];
    }
}

void AliAnalysisTaskStrangeCascadesRun2::SetMassbinning(int ipart, int nummassbins, double valminmass, double valmaxmass)
{
    fnmassbins[ipart] = nummassbins;
    for (int i = 0; i < fnmassbins[ipart] + 1; i++)
    {
        fmassbinning[ipart][i] = valminmass + i * (valmaxmass - valminmass) / fnmassbins[ipart];
    }
}

///________________________________________________________________________
void AliAnalysisTaskStrangeCascadesRun2::DataPosting()
{

    PostData(1, fHistos_eve->GetListOfHistograms());

    /// Histograms for analysed particle specie
    int histnumber = 1;

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

/*
void AliAnalysisTaskStrangeCascadesRun2::FillHistCutVariations(double perc, bool phypri, bool *associFlag)
{
    for (int iCutEv = 0; iCutEv < kNumCascEvCuts; iCutEv++)
    {
        if (iCutEv == kRapidityIntervalMin || iCutEv == kRapidityIntervalMax || iCutEv == kLeastCRows || iCutEv == kLeastCRowsOvF || iCutEv == kTrackLengthCut || iCutEv == kEtaDaughter || iCutEv == kBacBarCosPa)
            continue;
        for (int iVarEv = 0; iVarEv < kNumCutVars; iVarEv++)
        {
            if (var_cutValEv[iCutEv][iVarEv] == -1) // skip the variation if value == -1 -> unused variation
                continue;
            //             if (fisParametricBacBarCosPA && iCutEv != kCasc_BacBarCosPA && perc < fCentLimit_BacBarCosPA)
            //             {
            //                 if (fCasc_Pt >= fHist_PtBacBarCosPA->GetXaxis()->GetXmin() && fCasc_Pt <= fHist_PtBacBarCosPA->GetXaxis()->GetXmax())
            //                 {
            //                     SetEvCutVal(kFALSE, kTRUE, kCasc_BacBarCosPA, fHist_PtBacBarCosPA->GetBinContent(fHist_PtBacBarCosPA->GetXaxis()->FindBin(fCasc_Pt)));
            //                 }
            //                 else
            //                 {
            //                     SetEvCutVal(kFALSE, kTRUE, kCasc_BacBarCosPA, fCasc_Cuts[kCasc_BacBarCosPA]);
            //                 }
            //             }
            //             if (fisParametricTrackLengthCut && iCutEv != kCasc_TrackLengthCut)
            //             {
            //                 if (perc >= fHist_CentTrackLengthCut->GetXaxis()->GetXmin() && perc <= fHist_CentTrackLengthCut->GetXaxis()->GetXmax())
            //                 {
            //                     SetEvCutVal(kFALSE, kTRUE, kCasc_TrackLengthCut, fHist_CentTrackLengthCut->GetBinContent(fHist_CentTrackLengthCut->GetXaxis()->FindBin(perc)));
            //                 }
            //                 else
            //                 {
            //                     SetEvCutVal(kFALSE, kTRUE, kCasc_TrackLengthCut, fCasc_Cuts[kCasc_TrackLengthCut]);
            //                 }
            //             }
            // Xi filling
            // if (iCutEv != kCasc_PropLifetOm && iCutEv != kCompetingCascRejectOm)
            //             {
            // if (fParticleAnalysisStatus[kxip])
            // {
            // SetEvCutVal(kFALSE, kTRUE, iCutEv, varlowcut_Casc[iCutEv] + iVarEv * (varhighcut_Casc[iCutEv] - varlowcut_Casc[iCutEv]) / (nvarcut_Casc[iCutEv] - 1));
            SetCutValue(false, iCutEv, var_cutValEv[iCutEv][iVarEv]);

            if (fParticleAnalysisStatus[kxip])
            {
                if (phypri && associFlag[kxim] && ApplyCuts(kxim))
                    fHistos_XiMin->FillTH3(TString::Format("h3Var_%s_%s[%d][%d]", cutNamesEv[iCutEv].Data(), cutVarNames[iVarEv].Data(), iCutEv, iVarEv), fCasc_Pt, fCasc_InvMassXiMin, perc);
                if (phypri && associFlag[kxip] && ApplyCuts(kxip))
                    fHistos_XiPlu->FillTH3(TString::Format("h3Var_%s_%s[%d][%d]", cutNamesEv[iCutEv].Data(), cutVarNames[iVarEv].Data(), iCutEv, iVarEv), fCasc_Pt, fCasc_InvMassXiPlu, perc);
            }
            if (fParticleAnalysisStatus[komp])
            {
                if (phypri && associFlag[komm] && ApplyCuts(komm))
                    fHistos_OmMin->FillTH3(TString::Format("h3Var_%s_%s[%d][%d]", cutNamesEv[iCutEv].Data(), cutVarNames[iVarEv].Data(), iCutEv, iVarEv), fCasc_Pt, fCasc_InvMassOmMin, perc);
                if (phypri && associFlag[komp] && ApplyCuts(komp))
                    fHistos_OmPlu->FillTH3(TString::Format("h3Var_%s_%s[%d][%d]", cutNamesEv[iCutEv].Data(), cutVarNames[iVarEv].Data(), iCutEv, iVarEv), fCasc_Pt, fCasc_InvMassOmPlu, perc);
            }
        }
        SetDefCuts(); // reset defaults
    }

    for (int iCutTopo = 0; iCutTopo < kNumCascTopoCuts; iCutTopo++)
    {
        for (int iVarTopo = 0; iVarTopo < kNumCutVars; iVarTopo++)
        {
            // if (var_cutValTopo[kXi][iCutTopo][kMid][iVarTopo] == -1) // skip the variation if ptInterval values == -1 -> unused variation (including all pT intervals)
            //     continue;

            for (int iPtTopo = 0; iPtTopo < kNumPtInterval; iPtTopo++)
            {
                SetCutValue(true, iCutTopo, var_cutValTopo[kXi][iCutTopo][iPtTopo][iVarTopo], kXi, iPtTopo);
                SetCutValue(true, iCutTopo, var_cutValTopo[kOm][iCutTopo][iPtTopo][iVarTopo], kOm, iPtTopo);
            }

            if (fParticleAnalysisStatus[kxip] && (var_cutValTopo[kXi][iCutTopo][kMid][iVarTopo] != -1))
            {
                // SetEvCutVal(kFALSE, kTRUE, iCutEv, varlowcut_Casc[iCutEv] + iVarEv * (varhighcut_Casc[iCutEv] - varlowcut_Casc[iCutEv]) / (nvarcut_Casc[iCutEv] - 1));
                if (phypri && associFlag[kxim] && ApplyCuts(kxim))
                    fHistos_XiMin->FillTH3(TString::Format("h3Var_%s_%s[%d][%d]", cutNamesTopo[iCutTopo].Data(), cutVarNames[iVarTopo].Data(), iCutTopo, iVarTopo), fCasc_Pt, fCasc_InvMassXiMin, perc);
                if (phypri && associFlag[kxip] && ApplyCuts(kxip))
                    fHistos_XiPlu->FillTH3(TString::Format("h3Var_%s_%s[%d][%d]", cutNamesTopo[iCutTopo].Data(), cutVarNames[iVarTopo].Data(), iCutTopo, iVarTopo), fCasc_Pt, fCasc_InvMassXiPlu, perc);
            }

            if (fParticleAnalysisStatus[komp] && (var_cutValTopo[kOm][iCutTopo][kMid][iVarTopo] != -1))
            {
                // SetEvCutVal(kFALSE, kTRUE, iCutEv, varlowcut_Casc[iCutEv] + iVarEv * (varhighcut_Casc[iCutEv] - varlowcut_Casc[iCutEv]) / (nvarcut_Casc[iCutEv] - 1));
                if (phypri && associFlag[komm] && ApplyCuts(komm))
                    fHistos_OmMin->FillTH3(TString::Format("h3Var_%s_%s[%d][%d]", cutNamesTopo[iCutTopo].Data(), cutVarNames[iVarTopo].Data(), iCutTopo, iVarTopo), fCasc_Pt, fCasc_InvMassOmMin, perc);
                if (phypri && associFlag[komp] && ApplyCuts(komp))
                    fHistos_OmPlu->FillTH3(TString::Format("h3Var_%s_%s[%d][%d]", cutNamesTopo[iCutTopo].Data(), cutVarNames[iVarTopo].Data(), iCutTopo, iVarTopo), fCasc_Pt, fCasc_InvMassOmPlu, perc);
            }
        }
        SetDefCuts(); // reset defaults
    }
}
*/

//________________________________________________________________________
void AliAnalysisTaskStrangeCascadesRun2::FillHistCutVariations(double perc, bool phypri, bool *associFlag)
{
    int all = -1;

    for (int iCutEv = 0; iCutEv < kNumCascEvCuts; iCutEv++)
    {
        if (iCutEv == kRapidityIntervalMin || iCutEv == kRapidityIntervalMax || iCutEv == kLeastCRows || iCutEv == kLeastCRowsOvF || iCutEv == kTrackLengthCut || iCutEv == kEtaDaughter || iCutEv == kBacBarCosPa)
            continue;
        if (nvarcut_Ev[iCutEv] == -1) // skip the cut if value == -1 -> unused variations
            continue;
        for (int iVarEv = 0; iVarEv < nvarcut_Ev[iCutEv]; iVarEv++)
        {
            //             if (fisParametricBacBarCosPA && iCutEv != kCasc_BacBarCosPA && perc < fCentLimit_BacBarCosPA)
            //             {
            //                 if (fCasc_Pt >= fHist_PtBacBarCosPA->GetXaxis()->GetXmin() && fCasc_Pt <= fHist_PtBacBarCosPA->GetXaxis()->GetXmax())
            //                 {
            //                     SetEvCutVal(kFALSE, kTRUE, kCasc_BacBarCosPA, fHist_PtBacBarCosPA->GetBinContent(fHist_PtBacBarCosPA->GetXaxis()->FindBin(fCasc_Pt)));
            //                 }
            //                 else
            //                 {
            //                     SetEvCutVal(kFALSE, kTRUE, kCasc_BacBarCosPA, fCasc_Cuts[kCasc_BacBarCosPA]);
            //                 }
            //             }
            //             if (fisParametricTrackLengthCut && iCutEv != kCasc_TrackLengthCut)
            //             {
            //                 if (perc >= fHist_CentTrackLengthCut->GetXaxis()->GetXmin() && perc <= fHist_CentTrackLengthCut->GetXaxis()->GetXmax())
            //                 {
            //                     SetEvCutVal(kFALSE, kTRUE, kCasc_TrackLengthCut, fHist_CentTrackLengthCut->GetBinContent(fHist_CentTrackLengthCut->GetXaxis()->FindBin(perc)));
            //                 }
            //                 else
            //                 {
            //                     SetEvCutVal(kFALSE, kTRUE, kCasc_TrackLengthCut, fCasc_Cuts[kCasc_TrackLengthCut]);
            //                 }
            //             }

            // SetEvCutVal(kFALSE, kTRUE, iCutEv, varlowcut_Casc[iCutEv] + iVarEv * (varhighcut_Casc[iCutEv] - varlowcut_Casc[iCutEv]) / (nvarcut_Casc[iCutEv] - 1));
            SetCutValue(kFALSE, iCutEv, (varlowcut_Ev[iCutEv] + iVarEv * (varhighcut_Ev[iCutEv] - varlowcut_Ev[iCutEv]) / (nvarcut_Ev[iCutEv] - 1)), all, all);

            // Xi filling
            if (iCutEv != kCompetingCascRejectOm)
            {
                if (fParticleAnalysisStatus[kxip])
                {
                    if (phypri && associFlag[kxim] && ApplyCuts(kxim))
                        fHistos_XiMin->FillTH3(TString::Format("h3Var_%s[%d][%d]", cutNamesEv[iCutEv].Data(), iCutEv, iVarEv), fCasc_Pt, fCasc_InvMassXiMin, perc);
                    if (phypri && associFlag[kxip] && ApplyCuts(kxip))
                        fHistos_XiPlu->FillTH3(TString::Format("h3Var_%s[%d][%d]", cutNamesEv[iCutEv].Data(), iCutEv, iVarEv), fCasc_Pt, fCasc_InvMassXiPlu, perc);
                }
            }

            // Omega filling
            if (fParticleAnalysisStatus[komp])
            {
                if (phypri && associFlag[komm] && ApplyCuts(komm))
                    fHistos_OmMin->FillTH3(TString::Format("h3Var_%s[%d][%d]", cutNamesEv[iCutEv].Data(), iCutEv, iVarEv), fCasc_Pt, fCasc_InvMassOmMin, perc);
                if (phypri && associFlag[komp] && ApplyCuts(komp))
                    fHistos_OmPlu->FillTH3(TString::Format("h3Var_%s[%d][%d]", cutNamesEv[iCutEv].Data(), iCutEv, iVarEv), fCasc_Pt, fCasc_InvMassOmPlu, perc);
            }
        }
        SetDefCuts(); // reset defaults
    }

    for (int iCutTopo = 0; iCutTopo < kNumCascTopoCuts; iCutTopo++)
    {
        if (nvarcut_Topo[iCutTopo] == -1) // skip the cut if value == -1 -> unused variations
            continue;

        for (int iVarTopo = 0; iVarTopo < nvarcut_Topo[iCutTopo]; iVarTopo++)
        {
            // if (var_cutValTopo[kXi][iCutTopo][kMid][iVarTopo] == -1) // skip the variation if ptInterval values == -1 -> unused variation (including all pT intervals)
            //     continue;

            // for (int iPtTopo = 0; iPtTopo < kNumPtInterval; iPtTopo++)
            // {
            //     SetCutValue(true, iCutTopo, var_cutValTopo[kXi][iCutTopo][iPtTopo][iVarTopo], kXi, iPtTopo);
            //     SetCutValue(true, iCutTopo, var_cutValTopo[kOm][iCutTopo][iPtTopo][iVarTopo], kOm, iPtTopo);
            // }

            // SetEvCutVal(kFALSE, kTRUE, iCutEv, varlowcut_Casc[iCutEv] + iVarEv * (varhighcut_Casc[iCutEv] - varlowcut_Casc[iCutEv]) / (nvarcut_Casc[iCutEv] - 1));
            SetCutValue(kTRUE, iCutTopo, (varlowcut_Topo[iCutTopo] + iVarTopo * (varhighcut_Topo[iCutTopo] - varlowcut_Topo[iCutTopo]) / (nvarcut_Topo[iCutTopo] - 1)), all, all);

            // if (fParticleAnalysisStatus[kxip] && (var_cutValTopo[kXi][iCutTopo][kMid][iVarTopo] != -1))
            // {
            if (fParticleAnalysisStatus[kxip])
            {
                // SetEvCutVal(kFALSE, kTRUE, iCutEv, varlowcut_Casc[iCutEv] + iVarEv * (varhighcut_Casc[iCutEv] - varlowcut_Casc[iCutEv]) / (nvarcut_Casc[iCutEv] - 1));
                if (phypri && associFlag[kxim] && ApplyCuts(kxim))
                    fHistos_XiMin->FillTH3(TString::Format("h3Var_%s[%d][%d]", cutNamesTopo[iCutTopo].Data(), iCutTopo, iVarTopo), fCasc_Pt, fCasc_InvMassXiMin, perc);
                if (phypri && associFlag[kxip] && ApplyCuts(kxip))
                    fHistos_XiPlu->FillTH3(TString::Format("h3Var_%s[%d][%d]", cutNamesTopo[iCutTopo].Data(), iCutTopo, iVarTopo), fCasc_Pt, fCasc_InvMassXiPlu, perc);
            }

            // if (fParticleAnalysisStatus[komp] && (var_cutValTopo[kOm][iCutTopo][kMid][iVarTopo] != -1))
            // {
            if (fParticleAnalysisStatus[komp])
            {
                // SetEvCutVal(kFALSE, kTRUE, iCutEv, varlowcut_Casc[iCutEv] + iVarEv * (varhighcut_Casc[iCutEv] - varlowcut_Casc[iCutEv]) / (nvarcut_Casc[iCutEv] - 1));
                if (phypri && associFlag[komm] && ApplyCuts(komm))
                    fHistos_OmMin->FillTH3(TString::Format("h3Var_%s[%d][%d]", cutNamesTopo[iCutTopo].Data(), iCutTopo, iVarTopo), fCasc_Pt, fCasc_InvMassOmMin, perc);
                if (phypri && associFlag[komp] && ApplyCuts(komp))
                    fHistos_OmPlu->FillTH3(TString::Format("h3Var_%s[%d][%d]", cutNamesTopo[iCutTopo].Data(), iCutTopo, iVarTopo), fCasc_Pt, fCasc_InvMassOmPlu, perc);
            }
        }
        SetDefCuts(); // reset defaults
    }
}

void AliAnalysisTaskStrangeCascadesRun2::Terminate(Option_t *)
{
}
