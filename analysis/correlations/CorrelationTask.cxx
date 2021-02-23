#include "TChain.h"
#include "TH1D.h"
#include "TList.h"
#include "TMath.h"
#include "TH2D.h"
#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"
#include "AliAODEvent.h"
#include "AliAODInputHandler.h"
#include "CorrelationTask.h"

/// Pool Manager for mixing events
#include "AliEventPoolManager.h"
#include <AliMultiInputEventHandler.h>
#include <AliMixInputEventHandler.h>

#include "AliCentrality.h"

class CorrelationTask;

Double_t const kPI = TMath::Pi();

Double_t calcPhi(Double_t phi)
{
    if (phi > (1.5 * kPI))
        phi -= 2.0 * kPI;
    if (phi < (-0.5 * kPI))
        phi += 2.0 * kPI;
    return phi;
}

using namespace std;

ClassImp(CorrelationTask)

    //_____________________________________________________________________________
    CorrelationTask::CorrelationTask(const char *name)

    : AliAnalysisTaskSE(name),
      fAOD(0),
      fOutputList(0),
      fHistdEta(0),
      fHistEtaTrig(0),
      fHistEtaAssoc(0),
      fHistdPhi(0),
      fHistPhiTrig(0),
      fHistPhiAssoc(0),
      fHistdPhidEta(0),
      fFillMixed(kTRUE),
      fMixingTracks(500),
      fPoolMgr(0x0),
      fHistMixC1(0),
      fHistMixC2(0)

{
    // constructor

    DefineInput(0, TChain::Class());
    DefineOutput(1, TList::Class());
}
//_____________________________________________________________________________
CorrelationTask::~CorrelationTask()
{
    // destructor
    if (fOutputList)
    {
        delete fOutputList;
    }
}
//_____________________________________________________________________________
void CorrelationTask::UserCreateOutputObjects()
{

    fOutputList = new TList();
    fOutputList->SetOwner(kTRUE);

    Int_t nCentralityBins = 9;
    Double_t centBins[] = {0., 10., 20., 30., 40., 50., 60., 70., 80., 90.};
    // Double_t centBins[2] = {0., 100.};
    const Double_t *centralityBins = centBins;
    // defining bins for Z vertex
    Int_t nZvtxBins = 7;
    Double_t vertexBins[] = {-7., -5., -3., -1., 1., 3., 5., 7.};
    // Double_t vertexBins[2] = {-7., 7.};
    const Double_t *zvtxBins = vertexBins;
    // pt bins of associated particles for the analysis
    Int_t nPtBins = 7;
    const Double_t PtBins[8] = {2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0};
    // pt bins of trigger particles for the analysis
    Int_t nPtBinsCh = 11;
    const Double_t PtBinsCh[12] = {4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0};
    // defining bins for dPhi distributions
    const Int_t nbPhiBins = 72;
    Double_t PhiMin = -1.57;
    Double_t PhiMax = 4.71;
    Double_t PhiBins[nbPhiBins + 1] = {0.};
    PhiBins[0] = PhiMin;
    for (Int_t i = 0; i < nbPhiBins; i++)
    {
        PhiBins[i + 1] = PhiBins[i] + (PhiMax - PhiMin) / nbPhiBins;
    }

    // defining bins for dEta distributions
    const Int_t nbEtaBins = 40;
    Double_t EtaMin = -1.6;
    Double_t EtaMax = 1.6;
    Double_t EtaBins[nbEtaBins + 1] = {0.};
    EtaBins[0] = EtaMin;
    for (Int_t i = 0; i < nbEtaBins; i++)
    {
        EtaBins[i + 1] = EtaBins[i] + (EtaMax - EtaMin) / nbEtaBins;
    }

    fHistdEta = new TH1D("fHistdEta", "dEta (triggers-associated); dEta; n", nbEtaBins, EtaMin, EtaMax);
    fHistEtaTrig = new TH1D("fHistEtaTrig", "Eta (Triggers); Eta; n", nbEtaBins, EtaMin, EtaMax);
    fHistEtaAssoc = new TH1D("fHistEtaAssoc", "Eta (Associated); Eta; n", nbEtaBins, EtaMin, EtaMax);

    fHistdPhi = new TH1D("fHistdPhi", "dPhi (triggers-associated); dPhi; n", nbPhiBins, PhiMin, PhiMax);
    fHistPhiTrig = new TH1D("fHistPhiTrig", "Phi (Triggers); Phi; n", nbPhiBins, PhiMin, PhiMax);
    fHistPhiAssoc = new TH1D("fHistPhiAssoc", "Phi (Associated); Phi; n", nbPhiBins, PhiMin, PhiMax);

    fHistdPhidEta = new TH2D("fHistdPhidEta", "dPhi vs. dEta; dPhi; dEta", nbPhiBins, PhiMin, PhiMax, nbEtaBins, EtaMin, EtaMax);
    const Int_t corBins[6] = {nbPhiBins, nbEtaBins, nPtBinsCh, nPtBins, nCentralityBins, nZvtxBins};
    const Double_t corMin[6] = {PhiBins[0], EtaBins[0], PtBinsCh[0], PtBins[0], centralityBins[0], zvtxBins[0]};
    const Double_t corMax[6] = {PhiBins[nbPhiBins], EtaBins[nbEtaBins], PtBinsCh[nPtBinsCh], PtBins[nPtBins], centralityBins[nCentralityBins], zvtxBins[nZvtxBins]};

    fHistMixC1 = new THnSparseD("fHistMixC1", "dPhi vs. dEta mixed", 6, corBins, corMin, corMax);
    fHistMixC1->GetAxis(0)->SetTitle("dPhiMix");
    fHistMixC1->GetAxis(1)->SetTitle("dEtaMix");
    fHistMixC1->GetAxis(2)->SetTitle("chTrigPt");
    fHistMixC1->GetAxis(3)->SetTitle("assocPt");
    fHistMixC1->GetAxis(4)->SetTitle("lCent");
    fHistMixC1->GetAxis(5)->SetTitle("lPVz");

    fHistMixC2 = (THnSparseD *)fHistMixC1->Clone("fHistMixC2");

    // Settings for event mixing
    Int_t trackDepth = fMixingTracks;
    // Int_t trackDepth = 5;
    Int_t poolSize = 200; // Maximum number of events, ignored in the present implemented of AliEventPoolManager
    // Int_t poolSize = 100;
    fPoolMgr = new AliEventPoolManager(poolSize, trackDepth, nCentralityBins, centBins, nZvtxBins, vertexBins);

    fOutputList->Add(fHistdEta);
    fOutputList->Add(fHistEtaTrig);
    fOutputList->Add(fHistEtaAssoc);

    fOutputList->Add(fHistdPhi);
    fOutputList->Add(fHistPhiTrig);
    fOutputList->Add(fHistPhiAssoc);

    fOutputList->Add(fHistdPhidEta);

    fOutputList->Add(fHistMixC1);
    fOutputList->Add(fHistMixC2);

    PostData(1, fOutputList);
}
//_____________________________________________________________________________
void CorrelationTask::UserExec(Option_t *)
{
    fAOD = dynamic_cast<AliAODEvent *>(InputEvent());
    if (!fAOD)
        return;
    Double_t PtAssocMin = 3;

    Int_t nTracks = fAOD->GetNumberOfTracks();

    TObjArray *selectedTracks = new TObjArray;
    selectedTracks->SetOwner(kTRUE);

    TObjArray *selectedChargedAssoc = new TObjArray;
    selectedChargedAssoc->SetOwner(kTRUE);

    TObjArray *selectedChargedTriggers = new TObjArray;
    selectedChargedTriggers->SetOwner(kTRUE);
    for (Int_t i = 0; i < nTracks; i++)
    {
        AliAODTrack *tr = dynamic_cast<AliAODTrack *>(fAOD->GetTrack(i));

        /// Filtering tracks
        if (!tr)
            continue;
        if (!tr->TestFilterBit(768))
            continue;
        if (tr->Pt() < PtAssocMin)
            continue;
        if (TMath::Abs(tr->Eta()) > 0.9)
            continue;

        /// saving selected tracks
        selectedTracks->Add(tr);

        /// saving associated tracks/particles
        if (tr->Pt() < 8.)
            selectedChargedAssoc->Add(tr);
        /// saving the Charged trigger particles
        if ((tr->Pt() >= 8.) && (tr->Pt() < 15.))
        {
            selectedChargedTriggers->Add(tr);
        }
    }

    Int_t nSelectedChargedTriggers = selectedChargedTriggers->GetEntries();
    Int_t nselectedChargedAssoc = selectedChargedAssoc->GetEntries();
    Double_t dEta = 0, dPhi = 0;
    Int_t atrCount = 0;

    for (Int_t j = 0; j < nSelectedChargedTriggers; j++)
    {
        AliAODTrack *chTrig = (AliAODTrack *)selectedChargedTriggers->At(j);
        fHistEtaTrig->Fill(chTrig->Eta());
        fHistPhiTrig->Fill(calcPhi(chTrig->Phi()));

        for (Int_t m = 0; m < nselectedChargedAssoc; m++)
        {
            AliAODTrack *atr = (AliAODTrack *)selectedChargedAssoc->At(m);

            dEta = (chTrig->Eta()) - (atr->Eta());
            dPhi = calcPhi((chTrig->Phi()) - (atr->Phi()));

            fHistdEta->Fill(dEta);
            fHistdPhi->Fill(dPhi);
            fHistdPhidEta->Fill(dPhi, dEta);
            if (atrCount < nselectedChargedAssoc)
            {
                fHistEtaAssoc->Fill(atr->Eta());
                fHistPhiAssoc->Fill(calcPhi(atr->Phi()));
                atrCount++;
            }
        }
    }

    // ________________Mixing_______________________

    // Vertex cut
    Double_t cutPrimVertex = 7.0;
    AliAODVertex *myPrimVertex = fAOD->GetPrimaryVertex();
    if (!myPrimVertex)
    {
        return;
    }
    if ((TMath::Abs(myPrimVertex->GetZ())) >= cutPrimVertex)
        return;

    Double_t lPVx = myPrimVertex->GetX();
    Double_t lPVy = myPrimVertex->GetY();
    Double_t lPVz = myPrimVertex->GetZ();

    if (TMath::Abs(lPVx) < 10e-5 && TMath::Abs(lPVy) < 10e-5 && TMath::Abs(lPVz) < 10e-5)
        return;
    // Centrality definition
    Double_t lCent = 0.0;
    AliCentrality *centralityObj = 0;
    centralityObj = ((AliVAODHeader *)fAOD->GetHeader())->GetCentralityP();
    lCent = centralityObj->GetCentralityPercentile("V0M");
    if ((lCent < 0.) || (lCent > 90.)) /// Centrality ranges for strangeness
        return;
    if ((lCent > 10.) && (lCent < 60.)) /// Centrality range for this particular case - ignore for (10 < lCent < 60)
        return;

    fHistMixC1->Sumw2();
    fHistMixC2->Sumw2();
    AliEventPool *pool = fPoolMgr->GetEventPool(lCent, lPVz);
    if (!pool)
        AliFatal(Form("No pool found for centrality = %f, zVtx = %f", lCent, lPVz));
    // Int_t tracks = 0; ///////////////////

    if (pool->IsReady() || pool->NTracksInPool() > fMixingTracks / 10 || pool->GetCurrentNEvents() >= 5)
    {

        Int_t nMix = pool->GetCurrentNEvents();
        for (Int_t jMix = 0; jMix < nMix; jMix++)
        { // loop through mixing events

            TObjArray *bgTracks = pool->GetEvent(jMix);
            for (Int_t i = 0; i < selectedChargedTriggers->GetEntriesFast(); i++)    /// instead of selected V0
            {                                                                        /// loop through selected charged trigger particles
                AliAODTrack *chTrig = (AliAODTrack *)selectedChargedTriggers->At(i); /// instead of AliV0ChBasicParticle
                for (Int_t j = 0; j < bgTracks->GetEntriesFast(); j++)
                { // mixing tracks loop
                    AliVParticle *assoc = (AliVParticle *)bgTracks->At(j);
                    // be careful tracks may have bigger pt than v0s.
                    if (((assoc->Pt()) >= chTrig->Pt()) || ((assoc->Pt()) < PtAssocMin))
                        continue;
                    Double_t dEtaMix = assoc->Eta() - chTrig->Eta();
                    Double_t dPhiMix = assoc->Phi() - chTrig->Phi();
                    if (dPhiMix > (1.5 * kPI))
                        dPhiMix -= 2.0 * kPI;
                    if (dPhiMix < (-0.5 * kPI))
                        dPhiMix += 2.0 * kPI;
                    Double_t spMix[6] = {dPhiMix, dEtaMix, chTrig->Pt(), assoc->Pt(), lCent, lPVz};
                    if (lCent < 10.)
                        fHistMixC1->Fill(spMix); /// fill for centrality range (0-10)
                    else
                        fHistMixC2->Fill(spMix); /// fill for remaining centrality ranges (60-90)
                }                                // end of mixing track loop
            }                                    // end of loop through selected charged trigger particles
        }                                        // end of loop of mixing events
    }

    // if (pool->NTracksInPool() != 0)
    // {
    //     tracks = tracks + pool->NTracksInPool();
    //     Printf("n tracks=%d", tracks); ///////////////
    //     // pool->PrintInfo();
    // }
    TObjArray *tracksClone = (TObjArray *)selectedTracks->Clone();
    tracksClone->SetOwner(kTRUE);
    pool->UpdatePool(tracksClone);
    /// Set custom histogram drawing options
    fHistEtaTrig->SetOption("EP");
    fHistEtaAssoc->SetOption("EP");
    fHistdEta->SetOption("EP");

    fHistPhiTrig->SetOption("EP");
    fHistPhiAssoc->SetOption("EP");
    fHistdPhi->SetOption("EP");

    fHistdPhidEta->SetOption("SURF1");

    /// Write objects to output list
    PostData(1, fOutputList);
}

void CorrelationTask::Terminate(Option_t *)
{
}
