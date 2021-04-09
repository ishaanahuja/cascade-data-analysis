int ProjectionsMC(TString output = "ProjectionsMC.root")
{
    Double_t const kPI = TMath::Pi();

    /// Defining random number generator parameters
    TRandom3 rndgen;
    // Double_t rndMin = -1.;
    // Double_t rndMax = 1.;

    /// Defining histogram parameters
    Int_t const dim = 6;
    Int_t nEta = 80;
    Int_t nPhi = 144;
    Int_t nbins[dim] = {nEta, nEta, 2 * nEta, nPhi, nPhi, nPhi};
    Double_t min[dim] = {-0.9, -0.9, -1.8, -kPI, -kPI, -1.57};
    Double_t max[dim] = {0.9, 0.9, 1.8, kPI, kPI, 4.71};
    Double_t values[dim];

    /// Simulate tracks using TVector3
    // TVector3 *track1 = new TVector3();
    // TVector3 *track2 = new TVector3();

    /// Creating histogram THnSparse
    // THnSparseF *sparsePhiEta = new THnSparseF("sparsePhiEta", "Uniform Distribution", dim, nbins, min, max);
    // sparsePhiEta->GetAxis(0)->SetTitle("Eta Track1");
    // sparsePhiEta->GetAxis(1)->SetTitle("Eta Track2");
    // sparsePhiEta->GetAxis(2)->SetTitle("dEta");
    // sparsePhiEta->GetAxis(3)->SetTitle("Phi Track1");
    // sparsePhiEta->GetAxis(4)->SetTitle("Phi Track2");
    // sparsePhiEta->GetAxis(5)->SetTitle("dPhi");
    // sparsePhiEta->Sumw2();

    /// Method 1: Generate uniform random track position distribution
    // for (int i = 0; i < 10000; ++i)
    // {
    //     /// filling tracks with random cartesian coordinates
    //     track1->SetXYZ(rndgen.Uniform(rndMin, rndMax), rndgen.Uniform(rndMin, rndMax), rndgen.Uniform(rndMin, rndMax));
    //     track2->SetXYZ(rndgen.Uniform(rndMin, rndMax), rndgen.Uniform(rndMin, rndMax), rndgen.Uniform(rndMin, rndMax));
    //     values[0] = track1->Eta();                 // eta1
    //     values[1] = track2->Eta();                 // eta2
    //     values[2] = track1->Eta() - track2->Eta(); // dEta
    //     values[3] = track1->Phi();                 // Phi1
    //     values[4] = track2->Phi();                 // Phi2
    //     values[5] = track1->Phi() - track2->Phi(); // dPhi

    //     sparsePhiEta->Fill(values);
    // }

    /// Method 2: Generate random Eta and Phi distributions
    TH1D *etaTrack1 = new TH1D("etaTrack1", "Track 1 Eta", nEta, min[0], max[0]);
    TH1D *etaTrack2 = new TH1D("etaTrack2", "Track 2 Eta", nEta, min[1], max[1]);
    TH1D *phiTrack1 = new TH1D("phiTrack1", "Track 1 Phi", nPhi, min[3], max[3]);
    TH1D *phiTrack2 = new TH1D("phiTrack2", "Track 2 Phi", nPhi, min[4], max[4]);
    TH1D *dPhi = new TH1D("dPhi", "dPhi", nPhi, min[5], max[5]);
    TH1D *dEta = new TH1D("dEta", "dEta", nEta, min[2], max[2]);

    TH2D *track1PhiEta = new TH2D("track1PhiEta", "track1PhiEta", nPhi / 2, min[3], max[3], nEta / 2, min[0], max[0]);
    TH2D *track2PhiEta = new TH2D("track2PhiEta", "track2PhiEta", nPhi / 2, min[4], max[4], nEta / 2, min[1], max[1]);

    TH2D *dPhidEta = new TH2D("dPhidEta", "dPhidEta", nPhi / 2, min[5], max[5], nEta / 2, min[2], max[2]);
    Double_t dPhiVal = 0.;

    for (int i = 0; i < 1000000; ++i)
    {
        values[0] = rndgen.Uniform(min[0], max[0]); // Eta1
        values[1] = rndgen.Uniform(min[1], max[1]); // Eta2
        values[2] = values[1] - values[0];          // dEta
        values[3] = rndgen.Uniform(min[3], max[3]); // Phi1
        values[4] = rndgen.Uniform(min[4], max[4]); // Phi1
        values[5] = values[4] - values[3];          // dPhi
        dPhiVal = values[5];
        if (dPhiVal > (1.5 * kPI))
            dPhiVal -= 2.0 * kPI;
        if (dPhiVal < (-0.5 * kPI))
            dPhiVal += 2.0 * kPI;
        etaTrack1->Fill(values[0]);
        etaTrack2->Fill(values[1]);
        phiTrack1->Fill(values[3]);
        phiTrack2->Fill(values[4]);
        dPhi->Fill(dPhiVal);
        dEta->Fill(values[2]);
        track1PhiEta->Fill(values[3], values[0]);
        track2PhiEta->Fill(values[4], values[1]);
        dPhidEta->Fill(dPhiVal, values[2]);
    }

    // sparsePhiEta->Print();
    track1PhiEta->SetOption("SURF1");
    track2PhiEta->SetOption("SURF1");
    dPhidEta->SetOption("SURF1");

    // Writing histogram "sparsePhiEta" to file
    TFile *out = TFile::Open(output.Data(), "RECREATE");
    if (!out)
    {
        Printf("Error: Cannot open file '%s' !", output.Data());
        return 1;
    }
    Printf("Saving output to '%s' ...", output.Data());

    // sparsePhiEta->Write();
    etaTrack1->Write();
    etaTrack2->Write();
    phiTrack1->Write();
    phiTrack2->Write();
    dPhi->Write();
    dEta->Write();
    track1PhiEta->Write();
    track2PhiEta->Write();
    dPhidEta->Write();
    out->Close();
    delete out;
    delete track1PhiEta, track2PhiEta, dPhidEta, dEta, dPhi, phiTrack2, phiTrack1, etaTrack2, etaTrack1;
    return 0;
}
