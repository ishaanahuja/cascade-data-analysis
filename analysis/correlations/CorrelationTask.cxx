#include "TChain.h"
#include "TH1F.h"
#include "TList.h"
#include "TMath.h"
#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"
#include "AliAODEvent.h"
#include "AliAODInputHandler.h"
#include "CorrelationTask.h"

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

    CorrelationTask::CorrelationTask() : AliAnalysisTaskSE(),
                                         fAOD(0), fOutputList(0), fHistdEta(0), fHistEtaTrig(0), fHistEtaAssoc(0), fHistdPhi(0), fHistPhiTrig(0), fHistPhiAssoc(0), fHistdEtadPhi(0)
{
}
//_____________________________________________________________________________
CorrelationTask::CorrelationTask(const char *name) : AliAnalysisTaskSE(name),
                                                     fAOD(0), fOutputList(0), fHistdEta(0), fHistEtaTrig(0), fHistEtaAssoc(0), fHistdPhi(0), fHistPhiTrig(0), fHistPhiAssoc(0), fHistdEtadPhi(0)
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

    fHistdEtadPhi = new TH2F("fHistdEtadPhi", "dEta vs. dPhi; dEta; dPhi", 200, -2, 2, 144, -1.57, 4.71);

    fOutputList->Add(fHistdEta);
    fOutputList->Add(fHistEtaTrig);
    fOutputList->Add(fHistEtaAssoc);

    fOutputList->Add(fHistdPhi);
    fOutputList->Add(fHistPhiTrig);
    fOutputList->Add(fHistPhiAssoc);

    fOutputList->Add(fHistdEtadPhi);

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
        // instead of AliAODTrack you can try TParticle, maybe it will work
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

    /// Set custom histogram drawing options
    fHistEtaTrig->SetLineColor(kRED);
    fHistEtaTrig->SetOption("CP*");
    fHistEtaAssoc->SetOption("CP*");
    fHistdEta->SetOption("CP*");

    fHistPhiTrig->SetLineColor(kRED);
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
