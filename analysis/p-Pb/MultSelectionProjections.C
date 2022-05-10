#include <TString.h>

int MultSelectionProjections(TString outputFilename = "MultSelectionProjections.root")
{
    // TDirectory::AddDirectory(0);

    TString inputRunFiles[] = {"/var/home/ishaan/Work/CERN/ishaan-ahuja/analysis/Analysis_Results/pPb6runsAfterEventSel/pPb6runsAfterEventSel.root", "/var/home/ishaan/Work/CERN/ishaan-ahuja/analysis/Analysis_Results/pPbGrp1PileUp2/pPbGrp1PileUp2.root", "/var/home/ishaan/Work/CERN/ishaan-ahuja/analysis/Analysis_Results/pPbGrp2PileUp3/pPbGrp2PileUp3.root", "/var/home/ishaan/Work/CERN/ishaan-ahuja/analysis/Analysis_Results/pPbGrp3PileUp5/pPbGrp3PileUp5.root"};
    Int_t totalRunFiles = 4;

    TH2 *eventSelections[totalRunFiles];
    TH1D *eventSelectionProj[totalRunFiles];
    TString histName = "fHistEventSelections";

    for (int inputIndex = 0; inputIndex < totalRunFiles; inputIndex++)
    {
        Printf("Opening file '%s'...", inputRunFiles[inputIndex].Data());
        TFile *f = TFile::Open(inputRunFiles[inputIndex].Data());
        if (!f)
        {
            Printf("Error: Cannot open file '%s' !", inputRunFiles[inputIndex].Data());
            return 1;
        }

        // needed so we can do file->Close()
        TH1::AddDirectory(0);

        TList *listMultSel = (TList *)f->FindObjectAny("cListMultSelection");

        eventSelections[inputIndex] = (TH2 *)listMultSel->FindObject(histName.Data());
        if (!eventSelections[inputIndex])
        {
            Printf("%s not found in %s. Aborting...", histName.Data(), inputRunFiles[inputIndex].Data());
            return 99;
        }
        f->Close();
    }
    /// Input ended!

    /// Saving output :
    TString output = outputFilename;

    Printf("Saving output to '%s' ...", output.Data());
    TFile *out = TFile::Open(output.Data(), "RECREATE");
    if (!out)
    {
        Printf("Error: Cannot open file '%s' !", output.Data());
        return 3;
    }

    /// Begin

    for (int inputIndex = 0; inputIndex < totalRunFiles; inputIndex++)
    {
        eventSelectionProj[inputIndex] = (TH1D *)eventSelections[inputIndex]->ProjectionY(TString::Format(("eventSelectionProj[%d]"), inputIndex));
        eventSelectionProj[inputIndex]->SetTitle(TString::Format(("Dataset [%d]"), inputIndex));
        eventSelectionProj[inputIndex]->Write();
    }

    out->Close();
    delete out;

    TCanvas *c = new TCanvas("c", "c", 1920, 1080);
    c->Divide(2, 2);

    for (int inputIndex = 0; inputIndex < totalRunFiles; inputIndex++)
    {
        c->cd(inputIndex + 1);
        eventSelectionProj[inputIndex]->Draw();
    }
    return 0;
}