#include "TChain.h"
#include "TH1D.h"
#include "TList.h"
#include "TMath.h"
#include "TH2F.h"
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
      fHistdEtadPhi(0),
      fFillMixed(kTRUE),
      fMixingTracks(50000),
      fPoolMgr(0x0),
      fHistdPhidEtaMix(0)

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

    fHistdEta = new TH1D("fHistdEta", "dEta (triggers-associated); dEta; n", 200, -2, 2);
    fHistEtaTrig = new TH1D("fHistEtaTrig", "Eta (Triggers); Eta; n", 100, -1, 1);
    fHistEtaAssoc = new TH1D("fHistEtaAssoc", "Eta (Associated); Eta; n", 100, -1, 1);

    fHistdPhi = new TH1D("fHistdPhi", "dPhi (triggers-associated); dPhi; n", 144, -1.57, 4.71);
    fHistPhiTrig = new TH1D("fHistPhiTrig", "Phi (Triggers); Phi; n", 144, -1.57, 4.71);
    fHistPhiAssoc = new TH1D("fHistPhiAssoc", "Phi (Associated); Phi; n", 144, -1.57, 4.71);

    fHistdEtadPhi = new TH2F("fHistdEtadPhi", "dEta vs. dPhi; dEta; dPhi", 100, -2, 2, 72, -1.57, 4.71);

    Int_t nCentralityBins = 9;
    Double_t centBins[] = {0., 10., 20., 30., 40., 50., 60., 70., 80., 90.};
    // Int_t nCentralityBins = 1;
    // Double_t centBins[] = {0., 90.};
    const Double_t *centralityBins = centBins;
    // defining bins for Z vertex
    Int_t nZvtxBins = 7;
    Double_t vertexBins[] = {-7., -5., -3., -1., 1., 3., 5., 7.};
    // Int_t nZvtxBins = 1;
    // Double_t vertexBins[] = {-10., 10.};
    const Double_t *zvtxBins = vertexBins;
    // pt bins of associated particles for the analysis
    Int_t nPtBins = 7;
    const Double_t PtBins[8] = {2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0};
    //const Double_t PtBins[2] = {3.0,15.0};
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
    Double_t EtaMin = -2.0;
    Double_t EtaMax = 2.0;
    Double_t EtaBins[nbEtaBins + 1] = {0.};
    EtaBins[0] = EtaMin;
    for (Int_t i = 0; i < nbEtaBins; i++)
    {
        EtaBins[i + 1] = EtaBins[i] + (EtaMax - EtaMin) / nbEtaBins;
    }

    const Int_t corBins[6] = {nbPhiBins, nbEtaBins, nPtBinsCh, nPtBins, nCentralityBins, nZvtxBins};
    const Double_t corMin[6] = {PhiBins[0], EtaBins[0], PtBinsCh[0], PtBins[0], centralityBins[0], zvtxBins[0]};
    const Double_t corMax[6] = {PhiBins[72], EtaBins[40], PtBinsCh[11], PtBins[7], centralityBins[1], zvtxBins[1]};

    fHistdPhidEtaMix = new THnSparseF("fHistdPhidEtaMix", "dPhi vs. dEta mixed", 6, corBins, corMin, corMax);
    fHistdPhidEtaMix->GetAxis(0)->SetTitle("dPhiMix");
    fHistdPhidEtaMix->GetAxis(1)->SetTitle("dEtaMix");
    fHistdPhidEtaMix->GetAxis(2)->SetTitle("chTrigPt");
    fHistdPhidEtaMix->GetAxis(3)->SetTitle("assocPt");
    fHistdPhidEtaMix->GetAxis(4)->SetTitle("lCent");
    fHistdPhidEtaMix->GetAxis(5)->SetTitle("lPVz");

    // Settings for event mixing
    Int_t trackDepth = fMixingTracks;
    Int_t poolSize = 200; // Maximum number of events, ignored in the present implemented of AliEventPoolManager

    fPoolMgr = new AliEventPoolManager(poolSize, trackDepth, nCentralityBins, centBins, nZvtxBins, vertexBins);
    // ncentralitybins=1, centbins[2]=0..90,nZvtxBins=1, vertexBins=[-10,10]

    fOutputList->Add(fHistdEta);
    fOutputList->Add(fHistEtaTrig);
    fOutputList->Add(fHistEtaAssoc);

    fOutputList->Add(fHistdPhi);
    fOutputList->Add(fHistPhiTrig);
    fOutputList->Add(fHistPhiAssoc);

    fOutputList->Add(fHistdEtadPhi);
    fOutputList->Add(fHistdPhidEtaMix);

    PostData(1, fOutputList);
}
//_____________________________________________________________________________
void CorrelationTask::UserExec(Option_t *)
{
    fAOD = dynamic_cast<AliAODEvent *>(InputEvent());
    if (!fAOD)
        return;
    Double_t PtAssocMin = 1;

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
        if (tr->Pt() < 2.)
            selectedChargedAssoc->Add(tr);
        /// saving the Charged trigger particles
        if ((tr->Pt() >= 2.) && (tr->Pt() < 15.))
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
        // instead of AliAODTrack try TParticle, maybe it will work
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
            fHistdEtadPhi->Fill(dEta, dPhi);
            if (atrCount < nselectedChargedAssoc)
            {
                fHistEtaAssoc->Fill(atr->Eta());
                fHistPhiAssoc->Fill(calcPhi(atr->Phi()));
                atrCount++;
            }
        }
    }

    // ________________Mixing stuff_______________________

    // Vertex cut
    Double_t cutPrimVertex = 7.0;
    AliAODVertex *myPrimVertex = fAOD->GetPrimaryVertex();
    if (!myPrimVertex)
    {
        // Printf("AliAODVertex myPrimVertex not found. Skipping...\n"); ///////////////
        return;
    }
    if ((TMath::Abs(myPrimVertex->GetZ())) >= cutPrimVertex)
        return;
    // Printf("PV GetZ>primary vertex cut\n"); //////////////////////

    Double_t lPVx = myPrimVertex->GetX();
    Double_t lPVy = myPrimVertex->GetY();
    Double_t lPVz = myPrimVertex->GetZ();

    if (TMath::Abs(lPVx) < 10e-5 && TMath::Abs(lPVy) < 10e-5 && TMath::Abs(lPVz) < 10e-5)
        return;
    // Printf("PV out of bounds!!\n");/////////////////////
    // Centrality definition
    Double_t lCent = 0.0;
    AliCentrality *centralityObj = 0;
    centralityObj = ((AliVAODHeader *)fAOD->GetHeader())->GetCentralityP();
    lCent = centralityObj->GetCentralityPercentile("V0M");
    if ((lCent < 0.) || (lCent > 90.))
    {
        // Printf("Centrality out of bounds!!\n lCent=%f\n", lCent);
        return;
    }

    // Printf("TEST before mixing\n");

    // Mixing ==============================================

    fHistdPhidEtaMix->Sumw2();
    AliEventPool *pool = fPoolMgr->GetEventPool(lCent, lPVz);
    if (!pool)
        AliFatal(Form("No pool found for centrality = %f, zVtx = %f", lCent, lPVz));
    //pool->SetDebug(1);
    // Printf("test1\n");
    // pool->PrintInfo();
    if (pool->IsReady() || pool->NTracksInPool() > fMixingTracks / 10 || pool->GetCurrentNEvents() >= 5)
    {

        Int_t nMix = pool->GetCurrentNEvents();
        // Printf("test2\n nMix=%d", nMix);

        for (Int_t jMix = 0; jMix < nMix; jMix++)
        { // loop through mixing events
            // Printf("test3\n");

            TObjArray *bgTracks = pool->GetEvent(jMix);
            for (Int_t i = 0; i < selectedChargedTriggers->GetEntriesFast(); i++)    ///instead of selected V0
            {                                                                        // loop through selected charged trigger particles
                AliAODTrack *chTrig = (AliAODTrack *)selectedChargedTriggers->At(i); /// instead of AliV0ChBasicParticle
                Double_t chTrigPhi = chTrig->Phi();
                Double_t chTrigEta = chTrig->Eta();
                Double_t chTrigPt = chTrig->Pt();
                // Short_t trigC = trig->WhichCandidate();
                for (Int_t j = 0; j < bgTracks->GetEntriesFast(); j++)
                { // mixing tracks loop
                    AliVParticle *assoc = (AliVParticle *)bgTracks->At(j);
                    // be careful tracks may have bigger pt than v0s.
                    if (((assoc->Pt()) >= chTrigPt) || ((assoc->Pt()) < PtAssocMin))
                        continue;
                    Double_t dEtaMix = assoc->Eta() - chTrigEta;
                    Double_t dPhiMix = assoc->Phi() - chTrigPhi;
                    if (dPhiMix > (1.5 * kPI))
                        dPhiMix -= 2.0 * kPI;
                    if (dPhiMix < (-0.5 * kPI))
                        dPhiMix += 2.0 * kPI;

                    Double_t spMix[6] = {dPhiMix, dEtaMix, chTrigPt, assoc->Pt(), lCent, lPVz};
                    fHistdPhidEtaMix->Fill(spMix);
                    // Printf("test4\n %f\n", spMix[0]);
                } // end of mixing track loop
            }     // end of loop through selected charged trigger particles
        }         // end of loop of mixing events
    }

    TObjArray *tracksClone = (TObjArray *)selectedTracks->Clone();
    tracksClone->SetOwner(kTRUE);
    pool->UpdatePool(tracksClone);

    /// Set custom histogram drawing options
    // fHistEtaTrig->SetLineColor(kRED);
    fHistEtaTrig->SetOption("CP*");
    fHistEtaAssoc->SetOption("CP*");
    fHistdEta->SetOption("CP*");

    // fHistPhiTrig->SetLineColor(kRED);
    fHistPhiTrig->SetOption("CP*");
    fHistPhiAssoc->SetOption("CP*");
    fHistdPhi->SetOption("CP*");

    fHistdEtadPhi->SetOption("SURF1");

    /// Write objects to output list
    PostData(1, fOutputList);
}

void CorrelationTask::Terminate(Option_t *)
{
}
