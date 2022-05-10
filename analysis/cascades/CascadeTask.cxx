#include "TChain.h"
#include "TH1D.h"
#include "TList.h"
#include "TMath.h"
#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"
#include "AliAODEvent.h"
#include "AliAODInputHandler.h"
#include "AliAODcascade.h"
#include "AliPIDResponse.h"
#include "AliPID.h"
#include "CascadeTask.h"
// #include "TopologicalCuts.h"

// #include "AliCentrality.h" /// enable for AODs < 2015
// #include "AliMultSelection.h"

class CascadeTask;

// Double_t const kPI = TMath::Pi();

// Double_t calcPhi(Double_t phi)
// {
//     if (phi > (1.5 * kPI))
//         phi -= 2.0 * kPI;
//     if (phi < (-0.5 * kPI))
//         phi += 2.0 * kPI;
//     return phi;
// }

/// Declaring constants
const Double_t M_LAMBDA = 1.115683;
const Double_t M_XI = 1.32171;
const Double_t M_OMEGA = 1.67245;

Bool_t Min(Double_t var1, Double_t var2)
{
    if ((var1 - var2) < 1e-8)
        return kTRUE;
    else
        return kFALSE;
}

Bool_t Max(Double_t var1, Double_t var2)
{
    if ((var1 - var2) > 1e-8)
        return kTRUE;
    else
        return kFALSE;
}

Bool_t IsSelected(TObject *cascade, TObject *pVertex, TObject *fCuts, /* test*/ TObject *testSparse /* */)
{
    /// Declaring cut variables
    Double_t fCascRadius = 1.2;           // Cascade transverse decay radius R2D
    Double_t fV0Radius = 3.0;             // V0 transverse decay radius
    Double_t fCascDcaBachToPrimVtx = 0.1; // DCA (bachelor - PV)
    Double_t fCascDcaV0ToPrimVtx = 0.1;   // DCA (V0 - PV)
    Double_t fCascDcaPosToPrimVtx = 0.2;  // DCA (meson V0 track - PV)
    Double_t fCascDcaNegToPrimVtx = 0.2;  // DCA (baryon V0 track - PV)
    Double_t fV0DcaDaughters = 1.0;       // DCA (V0 tracks)
    Double_t fCascDcaBachToV0 = 1.0;      // DCA (bachelor - V0)
    Double_t fCascCosPA = 0.95;           // Cascade cos(PA)
    Double_t fV0CosPA = 0.95;             // V0 cos(PA)
    Double_t fV0InvMassWindow = 0.005;    // V0 invariant mass window
    Double_t fCascRapInterval = 0.5;      // Rapidity Interval
    Double_t fTPCNSigma = 4;              // TPC dE/dxSelection
    Double_t fCascPropTime = 15;          // Proper Lifetime (mL/p)
    Double_t fV0PropTime = 30;            // V0 Lifetime (mL/p)
    Double_t fCascOmegaRejection = 0.008; // Competing cascade rejection (onlyΩ)
    Double_t fNCrossedRows = 80;          // Daughter TrackNcrossedrows
    Double_t fTrackLength = 90;           // Daughter Track Length,L
    Double_t fNcrOverLength = 0.8;        // Daughter TrackNcrossed rows/L

    TH1D *cuts = (TH1D *)fCuts;                ////////////
    THnSparse *test = (THnSparse *)testSparse; ////////////
    AliAODcascade *Cascade = (AliAODcascade *)cascade;
    AliAODVertex *PrimVtx = (AliAODVertex *)pVertex;
    AliAODTrack *pTrack = dynamic_cast<AliAODTrack *>(Cascade->GetDaughter(0));
    AliAODTrack *nTrack = dynamic_cast<AliAODTrack *>(Cascade->GetDaughter(1));
    AliAODTrack *bachTrack = dynamic_cast<AliAODTrack *>(Cascade->GetDecayVertexXi()->GetDaughter(0));

    if (!pTrack || !nTrack || !bachTrack)
    {
        return kFALSE;
    }

    // AliAnalysisManager *man = AliAnalysisManager::GetAnalysisManager();
    // if (man)
    // {
    //     AliAODInputHandler *inputHandler = (AliAODInputHandler *)(man->GetInputEventHandler());
    //     if (inputHandler)
    //         PIDResponse = inputHandler->GetPIDResponse();
    // }

    /// Calculating physical quantities:
    Double_t CascP = TMath::Sqrt(TMath::Sq(Cascade->MomXiX()) + TMath::Sq(Cascade->MomXiY()) + TMath::Sq(Cascade->MomXiZ()));                                                                                                                            // Cascade Total Momentum
    Double_t V0P = TMath::Sqrt(TMath::Sq(Cascade->MomV0X()) + TMath::Sq(Cascade->MomV0Y()) + TMath::Sq(Cascade->MomV0Z()));                                                                                                                              // V0 Total Momentum
    Double_t CascLength = TMath::Sqrt(TMath::Sq((Cascade->DecayVertexXiX()) - (PrimVtx->GetX())) + TMath::Sq((Cascade->DecayVertexXiZ()) - (PrimVtx->GetZ())) + TMath::Sq((Cascade->DecayVertexXiY()) - (PrimVtx->GetY())));                             // Distance between cascade decay point and the primary vertex
    Double_t V0Length = TMath::Sqrt(TMath::Sq((Cascade->DecayVertexV0X()) - (Cascade->DecayVertexXiX())) + TMath::Sq((Cascade->DecayVertexV0Y()) - (Cascade->DecayVertexXiY())) + TMath::Sq((Cascade->DecayVertexV0Z()) - (Cascade->DecayVertexXiZ()))); // Distance between cascade decay point and the V0 decay point

    /// Getting cascade variables:
    Double_t CascRadius = TMath::Sqrt(TMath::Sq(Cascade->DecayVertexXiX()) + TMath::Sq(Cascade->DecayVertexXiY())); // Cascade transverse decay radius R2D
    Double_t V0Radius = TMath::Sqrt(TMath::Sq(Cascade->DecayVertexV0X()) + TMath::Sq(Cascade->DecayVertexV0Y()));   // V0 transverse decay radius
    Double_t CascDcaBachToPrimVtx = Cascade->DcaBachToPrimVertex();                                                 // DCA (bachelor - PV)
    Double_t CascDcaV0ToPrimVtx = Cascade->DcaV0ToPrimVertex();                                                     // DCA (V0 - PV)
    Double_t CascDcaPosToPrimVtx = Cascade->DcaPosToPrimVertex();                                                   // DCA (meson V0 track - PV)
    Double_t CascDcaNegToPrimVtx = Cascade->DcaNegToPrimVertex();                                                   // DCA (baryon V0 track - PV)
    Double_t V0DcaDaughters = Cascade->DcaV0Daughters();                                                            // DCA (V0 tracks)
    Double_t CascDcaBachToV0 = Cascade->DcaXiDaughters();                                                           // DCA (bachelor - V0)
    Double_t CascCosPA = Cascade->CosPointingAngleXi(PrimVtx->GetX(), PrimVtx->GetY(), PrimVtx->GetZ());            // Cascade cos(PA)
    Double_t V0CosPA = Cascade->CosPointingAngle(Cascade->GetDecayVertexXi());                                      // V0 cos(PA)
    Double_t V0InvMassWindow = TMath::Abs((Cascade->MassLambda()) - M_LAMBDA);                                      // V0 invariant mass window
    Double_t RapXi = TMath::Abs(Cascade->RapXi());                                                                  // Rapidity Interval
    Double_t RapOmega = TMath::Abs(Cascade->RapOmega());                                                            // Rapidity Interval
    // Double_t CascBachNSigmaKaon = TMath::Abs(PIDResponse->NumberOfSigmasTPC(bachTrack, AliPID::kKaon));             // NSigma Bachelor Kaon (Omega) for TPC dE/dx Selection
    // Double_t CascBachNSigmaPion = TMath::Abs(PIDResponse->NumberOfSigmasTPC(bachTrack, AliPID::kPion));             // NSigma Bachelor Pion (Xi) for TPC dE/dx Selection
    // Double_t CascNegNSigmaPion = TMath::Abs(PIDResponse->NumberOfSigmasTPC(nTrack, AliPID::kPion));                 // NSigma Negative Track Pion for TPC dE/dx Selection
    // Double_t CascNegNSigmaProton = TMath::Abs(PIDResponse->NumberOfSigmasTPC(nTrack, AliPID::kProton));             // NSigma Negative Track Proton for TPC dE/dx Selection
    // Double_t CascPosNSigmaPion = TMath::Abs(PIDResponse->NumberOfSigmasTPC(pTrack, AliPID::kPion));                 // NSigma Positive Track Pion for TPC dE/dx Selection
    // Double_t CascPosNSigmaProton = TMath::Abs(PIDResponse->NumberOfSigmasTPC(pTrack, AliPID::kProton));             // NSigma Positive Track Proton for TPC dE/dx Selection
    Double_t CascPropTime = -10.0; // Cascade Proper Lifetime (mL/p)
    Double_t V0PropTime = -10.0;   // V0 Lifetime (mL/p)

    if (CascP != 0)
    {
        CascPropTime = (M_XI * CascLength) / CascP;
    }
    if (V0P != 0)
    {
        V0PropTime = (M_LAMBDA * V0Length) / V0P;
    }

    Double_t CascOmegaRejection = TMath::Abs((Cascade->MassXi()) - M_XI); // Competing Cascade rejection (only Ω)

    Double_t pTrackCrossedRows = pTrack->GetTPCClusterInfo(2, 1); // Daughter Track TPC Crossed Rows: Ncrossedrows
    Double_t nTrackCrossedRows = nTrack->GetTPCClusterInfo(2, 1); // Daughter Track TPC Crossed Rows: Ncrossedrows
    /// TODO: check GetIntegratedLength in AN Code: should be 0-2m (~1m)
    Double_t pTrackLength = pTrack->GetIntegratedLength(); // Daughter Track Length: L
    Double_t nTrackLength = nTrack->GetIntegratedLength(); // Daughter Track Length: L
    Double_t pTrackCrOverLength = 0;                       // Crossed rows over track length: Ncrossed rows/L
    Double_t nTrackCrOverLength = 0;                       // Crossed rows over track length: Ncrossed rows/L

    if (pTrackLength > 0 && nTrackLength > 0)
    {
        pTrackCrOverLength = pTrackCrossedRows / pTrackLength;

        nTrackCrOverLength = nTrackCrossedRows / nTrackLength;
    }

    ///testing//////////////////
    Double_t testing[5] = {CascPropTime, V0PropTime, pTrackCrossedRows, pTrackLength, pTrackCrOverLength};
    test->Fill(testing);
    /////////////

    /// Applying cuts:
    if (Min(CascRadius, fCascRadius))
        return kFALSE;
    else
        cuts->AddBinContent(1);
    if (Min(V0Radius, fV0Radius))
        return kFALSE;
    else
        cuts->AddBinContent(2);
    if (Min(CascDcaBachToPrimVtx, fCascDcaBachToPrimVtx))
        return kFALSE;
    else
        cuts->AddBinContent(3);
    if (Min(CascDcaV0ToPrimVtx, fCascDcaV0ToPrimVtx))
        return kFALSE;
    else
        cuts->AddBinContent(4);
    if (Min(CascDcaPosToPrimVtx, fCascDcaPosToPrimVtx))
        return kFALSE;
    else
        cuts->AddBinContent(5);
    if (Min(CascDcaNegToPrimVtx, fCascDcaNegToPrimVtx))
        return kFALSE;
    else
        cuts->AddBinContent(6);
    if (Max(V0DcaDaughters, fV0DcaDaughters))
        return kFALSE;
    else
        cuts->AddBinContent(7);
    if (Min(CascDcaBachToV0, fCascDcaBachToV0))
        return kFALSE;
    else
        cuts->AddBinContent(8);
    if (Min(CascCosPA, fCascCosPA))
        return kFALSE;
    else
        cuts->AddBinContent(9);
    if (Min(V0CosPA, fV0CosPA))
        return kFALSE;
    else
        cuts->AddBinContent(10);
    if (Min(V0InvMassWindow, fV0InvMassWindow))
        return kFALSE;
    else
        cuts->AddBinContent(11);
    if (Min(RapXi, fCascRapInterval))
        return kFALSE;
    else
        cuts->AddBinContent(12);
    if (Min(RapOmega, fCascRapInterval))
        return kFALSE;
    else
        cuts->AddBinContent(13);
    // if (Min(CascBachNSigmaKaon, fTPCNSigma))
    //     return kFALSE;
    // if (Min(CascBachNSigmaPion, fTPCNSigma))
    //     return kFALSE;
    // if (Min(CascNegNSigmaPion, fTPCNSigma))
    //     return kFALSE;
    // if (Min(CascNegNSigmaProton, fTPCNSigma))
    //     return kFALSE;
    // if (Min(CascPosNSigmaPion, fTPCNSigma))
    //     return kFALSE;
    // if (Min(CascPosNSigmaProton, fTPCNSigma))
    //     return kFALSE;
    if (Min(CascPropTime, fCascPropTime))
        return kFALSE;
    else
        cuts->AddBinContent(14);
    if (V0PropTime < fV0PropTime)
        // return kFALSE;
        // else
        cuts->AddBinContent(15);
    if (Max(CascOmegaRejection, fCascOmegaRejection))
        return kFALSE;
    else
        cuts->AddBinContent(16);
    if (!(pTrack->IsOn(AliAODTrack::kTPCrefit)))
        return kFALSE;
    else
        cuts->AddBinContent(17);
    if (!(nTrack->IsOn(AliAODTrack::kTPCrefit)))
        return kFALSE;
    else
        cuts->AddBinContent(18);
    if (!(bachTrack->IsOn(AliAODTrack::kTPCrefit)))
        return kFALSE;
    else
        cuts->AddBinContent(19);
    if (Min(pTrackCrossedRows, fNCrossedRows))
        return kFALSE;
    else
        cuts->AddBinContent(20);
    if (Min(nTrackCrossedRows, fNCrossedRows))
        return kFALSE;
    else
        cuts->AddBinContent(21);
    if (Min(pTrackLength, fTrackLength))
        return kFALSE;
    else
        cuts->AddBinContent(22);
    if (Min(nTrackLength, fTrackLength))
        return kFALSE;
    else
        cuts->AddBinContent(23);
    if (Min(pTrackCrOverLength, fNcrOverLength))
        return kFALSE;
    else
        cuts->AddBinContent(24);
    if (Min(nTrackCrOverLength, fNcrOverLength))
        return kFALSE;
    else
        cuts->AddBinContent(25);

    return kTRUE;
}

using namespace std;

ClassImp(CascadeTask)

    //_____________________________________________________________________________
    CascadeTask::CascadeTask(const char *name)

    : AliAnalysisTaskSE(name),
      fAOD(0),
      fOutputList(0),
      fPIDResponse(0),
      fHistXiPt3(0),
      fHistXiPt4(0),
      fHistXiPt5(0),
      fHistXiPt6(0),
      fHistXiPt7(0),
      fHistOmegaPt3(0),
      fHistOmegaPt4(0),
      fHistOmegaPt5(0),
      fHistOmegaPt6(0),
      fHistOmegaPt7(0),
      fHistEtaXi(0),
      fHistEtaOmega(0),
      test(0), ///////////////////////////////
      fCuts(0) ///////////////////////////////////////////////
{
    // constructor

    DefineInput(0, TChain::Class());
    DefineOutput(1, TList::Class());
}
//_____________________________________________________________________________
CascadeTask::~CascadeTask()
{
    // destructor
    if (fOutputList)
    {
        delete fOutputList;
    }
}
//_____________________________________________________________________________
void CascadeTask::UserCreateOutputObjects()
{

    fOutputList = new TList();
    fOutputList->SetOwner(kTRUE);

    fHistXiPt3 = new TH1D("fHistXiPt3", "Invariant Mass Xi: pT(3,4) GeV/c", 150, 1.255, 1.405);
    fHistXiPt4 = new TH1D("fHistXiPt4", "Invariant Mass Xi: pT(4,5) GeV/c", 150, 1.255, 1.405);
    fHistXiPt5 = new TH1D("fHistXiPt5", "Invariant Mass Xi: pT(5,6) GeV/c", 150, 1.255, 1.405);
    fHistXiPt6 = new TH1D("fHistXiPt6", "Invariant Mass Xi: pT(6,7) GeV/c", 150, 1.255, 1.405);
    fHistXiPt7 = new TH1D("fHistXiPt7", "Invariant Mass Xi: pT(7,8) GeV/c", 150, 1.255, 1.405);

    fHistOmegaPt3 = new TH1D("fHistOmegaPt3", "Invariant Mass Omega: pT(3,4) GeV/c", 150, 1.60, 1.90);
    fHistOmegaPt4 = new TH1D("fHistOmegaPt4", "Invariant Mass Omega: pT(4,5) GeV/c", 150, 1.60, 1.90);
    fHistOmegaPt5 = new TH1D("fHistOmegaPt5", "Invariant Mass Omega: pT(5,6) GeV/c", 150, 1.60, 1.90);
    fHistOmegaPt6 = new TH1D("fHistOmegaPt6", "Invariant Mass Omega: pT(6,7) GeV/c", 150, 1.60, 1.90);
    fHistOmegaPt7 = new TH1D("fHistOmegaPt7", "Invariant Mass Omega: pT(7,8) GeV/c", 150, 1.60, 1.90);

    fHistEtaXi = new TH1D("fHistEtaXi", "Rapidity (eta) Xi", 100, -1, 1);
    fHistEtaOmega = new TH1D("fHistEtaOmega", "Rapidity (eta) Omega", 100, -1, 1);

    //////////////////////
    fCuts = new TH1D("fCuts", "Topological cuts", 25, 0, 25); ///////////
    // const Double_t cascProp[20] = {-50.0, -10.0, -5.0, 0.0, 1.0, 2.0, 4.0, 6.0, 8.0, 10.0, 12.0, 14.0, 16.0, 20.0, 25.0, 30.0, 35.0, 40.0, 50.0, 100.0};
    // const Double_t V0Prop[22] = {-50.0, -10.0, -5.0, 0.0, 1.0, 2.0, 4.0, 6.0, 8.0, 10.0, 12.0, 14.0, 16.0, 20.0, 25.0, 30.0, 35.0, 40.0, 50.0, 60.0, 70.0, 100.0};
    // const Double_t pCrossedRows[13] = {0.0, 10.0, 20.0, 40.0, 60.0, 80.0, 100.0, 120.0, 140.0, 160.0, 180.0, 200.0, 400.0};
    // const Double_t pLength[13] = {0, 10, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 400};
    // const Double_t pCrLength[17] = {0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 2.0};

    const Int_t bins[5] = {20, 22, 13, 13, 17};
    // const Double_t min[5] = {cascProp[0], V0Prop[0], pCrossedRows[0], pLength[0], pCrLength[0]};
    // const Double_t max[5] = {cascProp[20], V0Prop[22], pCrossedRows[13], pLength[13], pCrLength[17]};
    const Double_t min[5] = {-10.0, -10.0, 60.0, 360.0, 0.0};
    const Double_t max[5] = {50.0, 80.0, 190.0, 450.0, 1.0};
    test = new THnSparseD("test", "test cuts", 5, bins, min, max);
    test->GetAxis(0)->SetNameTitle("cascProp", "cascProp");
    test->GetAxis(1)->SetNameTitle("V0Prop", "V0Prop");
    test->GetAxis(2)->SetNameTitle("pCrossedRows", "pCrossedRows");
    test->GetAxis(3)->SetNameTitle("pLength", "pLength");
    test->GetAxis(4)->SetNameTitle("pCrLength", "pCrLength");

    fOutputList->Add(fCuts); /////////////
    fOutputList->Add(test);  /////////////
    fOutputList->Add(fHistXiPt3);
    fOutputList->Add(fHistXiPt4);
    fOutputList->Add(fHistXiPt5);
    fOutputList->Add(fHistXiPt6);
    fOutputList->Add(fHistXiPt7);

    fOutputList->Add(fHistOmegaPt3);
    fOutputList->Add(fHistOmegaPt4);
    fOutputList->Add(fHistOmegaPt5);
    fOutputList->Add(fHistOmegaPt6);
    fOutputList->Add(fHistOmegaPt7);

    fOutputList->Add(fHistEtaXi);
    fOutputList->Add(fHistEtaOmega);

    /// Set custom histogram drawing options
    fHistXiPt3->SetOption("EP");
    fHistXiPt4->SetOption("EP");
    fHistXiPt5->SetOption("EP");
    fHistXiPt6->SetOption("EP");
    fHistXiPt7->SetOption("EP");

    fHistOmegaPt3->SetOption("EP");
    fHistOmegaPt4->SetOption("EP");
    fHistOmegaPt5->SetOption("EP");
    fHistOmegaPt6->SetOption("EP");
    fHistOmegaPt7->SetOption("EP");

    fHistEtaXi->SetOption("EP");
    fHistEtaOmega->SetOption("EP");

    PostData(1, fOutputList);
}
//_____________________________________________________________________________
void CascadeTask::UserExec(Option_t *)
{
    fAOD = dynamic_cast<AliAODEvent *>(InputEvent());
    if (!fAOD)
        return;
    Double_t fPtMin = 3.;
    Double_t fPtMax = 8.;
    AliAODVertex *PrimVtx = fAOD->GetPrimaryVertex();
    Int_t nCascades = fAOD->GetNumberOfCascades();

    TObjArray *selectedCascades = new TObjArray;
    selectedCascades->SetOwner(kTRUE);

    // Double_t fCascRadius = -1.0;
    // Double_t fV0Radius = -1.0;
    // Double_t fCascDcaBachToPrimVtx = -1.0;
    // Double_t fCascDcaV0ToPrimVtx = -1.0;
    // Double_t fCascDcaPosToPrimVtx = -1.0;
    // Double_t fCascDcaNegToPrimVtx = -1.0;
    // Double_t fCascDcaBachToV0 = -1.0;
    // Double_t fCascCosPA = -10.0;
    // Double_t fV0CosPA = -10.0;
    // Double_t fCascPropTime = -10.0;
    // Double_t fV0PropTime = -10.0;
    // Double_t fCascP;
    // Double_t fV0P;
    // // Double_t dist10cm;
    // Double_t fCascBachNSigmaKaon = 0.0, fCascBachNSigmaPion = 0.0, fCascPosNSigmaProton = 0.0, fCascNegNSigmaProton = 0.0, fCascPosNSigmaPion = 0.0, fCascNegNSigmaPion = 0.0;
    // Double_t fV0DcaDaughters = -1.0;
    // Double_t fMassLambda=1.115683; ///only .h
    // Double_t fV0InvMassWindow = -1.0;
    // // Double_t deltaM_Lbar = -1.0;
    // Double_t fRapOmega = -1.0;
    // Double_t fRapXi = -1.0;
    // Bool_t fCascBachTracking = kFALSE, fV0PosTracking = kFALSE, fV0NegTracking = kFALSE;
    // // Bool_t bachNcls = kFALSE, V0posNcls = kFALSE, V0negNcls = kFALSE;
    // Float_t pTrackCrossedRows = -1.0;
    // Float_t nTrackCrossedRows = -1.0;
    // Double_t pTrackLength = -1.0;
    // Double_t nTrackLength = -1.0;
    // Double_t pTrackCrOverLength = -1.0;
    // Double_t nTrackCrOverLength = -1.0;

    for (Int_t i = 0; i < nCascades; i++)
    {
        AliAODcascade *cascade = static_cast<AliAODcascade *>(fAOD->GetCascade(i));
        // TopologicalCuts *cuts = new TopologicalCuts("cuts", "cuts");
        /// Filtering Cascades
        if (!cascade)
            continue;

        if (cascade->Pt() < fPtMin || cascade->Pt() > fPtMax)
            continue;
        if (!(IsSelected(cascade, PrimVtx, fCuts, test)))
            continue;
        selectedCascades->Add(cascade);
    }

    Int_t nselectedCascades = selectedCascades->GetEntries();

    for (Int_t j = 0; j < nselectedCascades; j++)
    {
        AliAODcascade *selectedCascade = (AliAODcascade *)selectedCascades->At(j);

        fHistEtaXi->Fill(TMath::Abs(selectedCascade->RapXi()));
        fHistEtaOmega->Fill(TMath::Abs(selectedCascade->RapOmega()));

        Double_t PtCascade = TMath::Sqrt(selectedCascade->Pt2Xi());

        if (PtCascade >= 3. && PtCascade < 4.)
        {
            fHistXiPt3->Fill(selectedCascade->MassXi());
            fHistOmegaPt3->Fill(selectedCascade->MassOmega());
        }

        else if (PtCascade >= 4. && PtCascade < 5.)
        {
            fHistXiPt4->Fill(selectedCascade->MassXi());
            fHistOmegaPt4->Fill(selectedCascade->MassOmega());
        }

        else if (PtCascade >= 5. && PtCascade < 6.)
        {
            fHistXiPt5->Fill(selectedCascade->MassXi());
            fHistOmegaPt5->Fill(selectedCascade->MassOmega());
        }
        else if (PtCascade >= 6. && PtCascade < 7.)
        {
            fHistXiPt6->Fill(selectedCascade->MassXi());
            fHistOmegaPt6->Fill(selectedCascade->MassOmega());
        }
        else if (PtCascade >= 7. && PtCascade < 8.)
        {
            fHistXiPt7->Fill(selectedCascade->MassXi());
            fHistOmegaPt7->Fill(selectedCascade->MassOmega());
        }
    }

    /// Write objects to output list
    PostData(1, fOutputList);
}

void CascadeTask::Terminate(Option_t *)
{
}
