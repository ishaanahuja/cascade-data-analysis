#include "TChain.h"
#include "TH1D.h"
#include "TList.h"
#include "TMath.h"
#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"
#include "AliAODEvent.h"
#include "AliAODInputHandler.h"
#include "AliAODcascade.h"

#include "CascadeTask.h"

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

using namespace std;

ClassImp(CascadeTask)
    // ClassImp(AliMixBasicParticle)

    //_____________________________________________________________________________
    CascadeTask::CascadeTask(const char *name)

    : AliAnalysisTaskSE(name),
      fAOD(0),
      fOutputList(0),
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
      fHistEtaOmega(0)
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
    Double_t PtMin = 3.;
    Double_t PtMax = 8.;

    Int_t nCascades = fAOD->GetNumberOfCascades();

    TObjArray *selectedCascades = new TObjArray;
    selectedCascades->SetOwner(kTRUE);

    for (Int_t i = 0; i < nCascades; i++)
    {
        AliAODcascade *Cascade = static_cast<AliAODcascade *>(fAOD->GetCascade(i));

        /// Filtering Cascades
        if (!Cascade)
            continue;

        if (Cascade->Pt() < PtMin || Cascade->Pt() > PtMax)
            continue;
        // if (TMath::Abs(Cascade->Eta()) > 0.9)
        //     continue;

        /// saving selected Cascades
        selectedCascades->Add(Cascade);
    }

    Int_t nselectedCascades = selectedCascades->GetEntries();

    for (Int_t j = 0; j < nselectedCascades; j++)
    {
        AliAODcascade *selectedCascade = (AliAODcascade *)selectedCascades->At(j);

        fHistEtaXi->Fill(selectedCascade->RapXi());
        fHistEtaOmega->Fill(selectedCascade->RapOmega());

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
