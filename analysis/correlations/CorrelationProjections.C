int CorrelationProjections(TString input = "CorrelationTask.root", TString output = "CorrelationProjections.root", TString sparseMix = "fHistMixC1", TString listName = "MyOutputContainer")
{
        TH1::AddDirectory(0);
        TFile *f = TFile::Open(input);
        if (!f)
        {
                Printf("Error: Cannot open file '%s' !", input.Data());
                return 1;
        }
        f->cd("dEta_dPhi");
        TList *list = (TList *)f->FindObjectAny(listName.Data());
        if (!list)
        {
                Printf("Error: Cannot open TList '%s' !", listName.Data());
                return 2;
        }

        TH2D *histdPhidEta = (TH2D *)list->FindObject("fHistdPhidEta");
        THnSparseD *histMix = (THnSparseD *)list->FindObject(sparseMix.Data());
        TH1D *histPhiTrig = (TH1D *)list->FindObject("fHistPhiTrig");

        if (!(histdPhidEta && histMix && histPhiTrig))
        {
                Printf("Histograms cannot be found. Aborting.");
                return 3;
        }

        TH2D *histMixProj = (TH2D *)histMix->Projection(1, 0, "E");
        TH2D *histRatio = new TH2D("histRatio", "dPhi vs. dEta : Sibling/Mixing; dPhi; dEta", histdPhidEta->GetXaxis()->GetNbins(), histdPhidEta->GetXaxis()->GetXmin(), histdPhidEta->GetXaxis()->GetXmax(), histdPhidEta->GetYaxis()->GetNbins(), histdPhidEta->GetYaxis()->GetXmin(), histdPhidEta->GetYaxis()->GetXmax());

        // Obtaining Ratio histogram
        Bool_t IsDivide = histRatio->Divide(histdPhidEta, histMixProj);
        if (!IsDivide)
        {
                Printf("Cannot divide histograms. Aborting.");
                return 4;
        }

        // Normalizing histogram
        histRatio->Scale(1 / (histPhiTrig->GetEntries()));

        /// TODO: normalize mixing 2D hist with highest bin (average)- dEta should be 1 at max
        /// You'll get (sort of) pair efficiency (2-particle)

        TH1D *histRatiodPhi = histRatio->ProjectionX();
        histRatiodPhi->SetNameTitle("histRatiodPhi", "dPhi : Sibling/Mixing");
        TH1D *histRatiodEta = histRatio->ProjectionY();
        histRatiodEta->SetNameTitle("histRatiodEta", "dEta : Sibling/Mixing");

        if (!(histRatiodPhi && histRatiodEta))
        {
                Printf("Projection failed. Aborting.");
                return 5;
        }

        TFile *out = TFile::Open(output.Data(), "RECREATE");
        if (!out)
        {
                Printf("Error: Cannot open file '%s' !", output.Data());
                return 6;
        }
        histMixProj->SetOption("SURF1");
        // histdPhidEta->SetOption("SURF1");
        histRatio->SetOption("SURF1");
        histRatiodPhi->SetOption("EP");
        histRatiodEta->SetOption("EP");

        Printf("Saving output to '%s' ...", output.Data());
        histMixProj->Write();
        // histdPhidEta->Write();
        histRatio->Write();
        histRatiodPhi->Write();
        histRatiodEta->Write();
        out->Close();
        f->Close();
        delete out, histMixProj, histdPhidEta, histRatio, histRatiodPhi, histRatiodEta, histMix, histPhiTrig, list, f;

        return 0;
}