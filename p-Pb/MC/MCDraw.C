#include <TString.h>
#include <TSystem.h>

void DrawAndSave(TH1 *peak, Double_t pPosition, Double_t pWidth, Bool_t saveImages, TString outputFolder, TString imageFormat);

inline void SaveImage(TString imagePath, TString imageName, TString imageFormat)
{
    if (gSystem->AccessPathName(imagePath.Data())) /// returns true if folder path does NOT exist
    {
        gSystem->mkdir(imagePath.Data(), kTRUE); // makes the path if it doesn't exist
        gSystem->Chmod(imagePath.Data(), 0755);
    }
    gPad->Print(Form("%s/%s.%s", imagePath.Data(), imageName.Data(), imageFormat.Data()), imageFormat.Data());
    gSystem->Chmod(Form("%s/%s.%s", imagePath.Data(), imageName.Data(), imageFormat.Data()), 0755);
}

int MCDraw(std::string input = "AliCascadeAnalysisMC_Ishaan_Fitting.root", TString outputFilename = "MCDraw.root", TString outputFolder = ".", Bool_t fisMC = kTRUE, Bool_t saveImages = kFALSE, Bool_t saveStack = kTRUE, TString imageFormat = "png")
{
    // TDirectory::AddDirectory(0);
    outputFolder = gSystem->ExpandPathName(outputFolder.Data());
    if (gSystem->AccessPathName(outputFolder.Data())) /// returns true if folder path does NOT exist
    {
        gSystem->mkdir(outputFolder.Data(), kTRUE); // makes the path if it doesn't exist
        gSystem->Chmod(outputFolder.Data(), 0755);
    }

    gStyle->SetOptFit(1111);
    Double_t lMass_Xi = 1.32171;
    Double_t lMass_Om = 1.67245;

    double multbins_Xi[11] = {0, 5, 10, 15, 20, 30, 40, 50, 60, 80, 100}; // V0A
    double multbins_Om[6] = {0, 5, 15, 30, 60, 100};                      // V0A
    Int_t nmultbins_Xi = sizeof(multbins_Xi) / sizeof(double) - 1;
    Int_t nmultbins_Om = sizeof(multbins_Om) / sizeof(double) - 1;

    double ptbins_Xi[] = {0.8, 1.1, 1.3, 1.5, 1.7, 1.9, 2.1, 2.3, 2.5, 2.7, 2.9, 3.1, 3.5, 4, 5.5};
    double ptbins_Om[] = {0.9, 1.6, 2., 2.4, 2.9, 3.5, 5};
    Int_t nptbins_Xi = sizeof(ptbins_Xi) / sizeof(double) - 1;
    Int_t nptbins_Om = sizeof(ptbins_Om) / sizeof(double) - 1;

    TH1 *h_MassXim_pt_mult[nptbins_Xi][nmultbins_Xi];
    TH1 *h_MassXip_pt_mult[nptbins_Xi][nmultbins_Xi];
    TH1 *h_MassOmm_pt_mult[nptbins_Om][nmultbins_Om];
    TH1 *h_MassOmp_pt_mult[nptbins_Om][nmultbins_Om];
    TH1 *resultParXip_pt_mult[nptbins_Xi][nmultbins_Xi];
    TH1 *resultParXim_pt_mult[nptbins_Xi][nmultbins_Xi];
    TH1 *resultParOmp_pt_mult[nptbins_Om][nmultbins_Om];
    TH1 *resultParOmm_pt_mult[nptbins_Om][nmultbins_Om];
    TH1D *rawPt_xim[nmultbins_Xi];
    TH1D *rawPt_xip[nmultbins_Xi];
    TH1D *rawPt_omm[nmultbins_Om];
    TH1D *rawPt_omp[nmultbins_Om];

    TH1 *h_MassXiC_pt_mult[nptbins_Xi][nmultbins_Xi];
    TH1 *h_MassOmC_pt_mult[nptbins_Om][nmultbins_Om];
    TH1 *resultParXiC_pt_mult[nptbins_Xi][nmultbins_Xi];
    TH1 *resultParOmC_pt_mult[nptbins_Om][nmultbins_Om];
    TH1D *rawPt_xiC[nmultbins_Xi];
    TH1D *rawPt_omC[nmultbins_Om];

    TH1D *eff_xim[nmultbins_Xi];
    TH1D *eff_xip[nmultbins_Xi];
    TH1D *eff_omm[nmultbins_Om];
    TH1D *eff_omp[nmultbins_Om];
    TH1D *eff_xiC[nmultbins_Xi];
    TH1D *eff_omC[nmultbins_Om];

    TH1D *eff_xiC_ratio[nmultbins_Xi];
    TH1D *eff_omC_ratio[nmultbins_Om];

    TH1 *resultParXip_pt[nptbins_Xi];
    TH1 *resultParXim_pt[nptbins_Xi];
    TH1 *resultParOmp_pt[nptbins_Om];
    TH1 *resultParOmm_pt[nptbins_Om];
    TH1 *resultParXiC_pt[nptbins_Xi];
    TH1 *resultParOmC_pt[nptbins_Om];

    TH1D *eff_pt_xim = new TH1D("eff_pt_xim", "Mult: 0-100%", nptbins_Xi, ptbins_Xi);
    TH1D *eff_pt_xip = new TH1D("eff_pt_xip", "Mult: 0-100%", nptbins_Xi, ptbins_Xi);
    TH1D *eff_pt_omm = new TH1D("eff_pt_omm", "Mult: 0-100%", nptbins_Om, ptbins_Om);
    TH1D *eff_pt_omp = new TH1D("eff_pt_omp", "Mult: 0-100%", nptbins_Om, ptbins_Om);
    TH1D *eff_pt_xiC = new TH1D("eff_pt_xiC", "Mult: 0-100%", nptbins_Xi, ptbins_Xi);
    TH1D *eff_pt_omC = new TH1D("eff_pt_omC", "Mult: 0-100%", nptbins_Om, ptbins_Om);

    TH1 *h_MassXim;
    TH1 *h_MassXip;
    TH1 *h_MassOmm;
    TH1 *h_MassOmp;
    TH1 *h_MassXiC;
    TH1 *h_MassOmC;
    TH1 *h_multBinEntries_Xi;
    TH1 *h_multBinEntries_Om;

    auto hs_xip = new THStack("xip", "");
    auto hs_xim = new THStack("xim", "");
    auto hs_omp = new THStack("omp", "");
    auto hs_omm = new THStack("omm", "");
    auto hs_xiC = new THStack("xiC", "");
    auto hs_omC = new THStack("omC", "");
    if (fisMC)
    {
        hs_xip->SetTitle("Rec. + Ass. #it{p}_{T} spectra #Xi^{+}");
        hs_xim->SetTitle("Rec. + Ass. #it{p}_{T} spectra #Xi^{-}");
        hs_omp->SetTitle("Rec. + Ass. #it{p}_{T} spectra #Omega^{+}");
        hs_omm->SetTitle("Rec. + Ass. #it{p}_{T} spectra #Omega^{-}");
        hs_xiC->SetTitle("Rec. + Ass. #it{p}_{T} spectra #Xi^{+} + #Xi^{-}");
        hs_omC->SetTitle("Rec. + Ass. #it{p}_{T} spectra #Omega^{+} + #Omega^{-}");
    }
    else
    {
        hs_xip->SetTitle("Raw #it{p}_{T} spectra #Xi^{+}");
        hs_xim->SetTitle("Raw #it{p}_{T} spectra #Xi^{-}");
        hs_omp->SetTitle("Raw #it{p}_{T} spectra #Omega^{+}");
        hs_omm->SetTitle("Raw #it{p}_{T} spectra #Omega^{-}");
        hs_xiC->SetTitle("Raw #it{p}_{T} spectra #Xi^{+} + #Xi^{-}");
        hs_omC->SetTitle("Raw #it{p}_{T} spectra #Omega^{+} + #Omega^{-}");
    }

    // THstack efficiency
    auto hs_xip_eff = new THStack("xip_eff", "Efficiency #Xi^{+}");
    auto hs_xim_eff = new THStack("xim_eff", "Efficiency #Xi^{-}");
    auto hs_omp_eff = new THStack("omp_eff", "Efficiency #Omega^{+}");
    auto hs_omm_eff = new THStack("omm_eff", "Efficiency #Omega^{-}");
    auto hs_xiC_eff = new THStack("xiC_eff", "Efficiency #Xi^{+} + #Xi^{-}");
    auto hs_omC_eff = new THStack("omC_eff", "Efficiency #Omega^{+} + #Omega^{-}");

    auto hs_xiC_eff_ratio = new THStack("xiC_eff_ratio", "#Xi^{+} + #Xi^{-} efficiency ratio");
    auto hs_omC_eff_ratio = new THStack("omC_eff_ratio", "#Omega^{+} + #Omega^{-} efficiency ratio");

    /// Getting histograms:

    {
        Printf("Opening file ...");
        TFile *f = TFile::Open(input.data());
        if (!f)
        {
            Printf("Error: Cannot open file '%s' !", input.data());
            return 1;
        }

        // needed so we can do file->Close()
        TH1::AddDirectory(0);

        h_MassXim = (TH1 *)f->Get("h_MassXim");
        h_MassXip = (TH1 *)f->Get("h_MassXip");
        h_MassXiC = (TH1 *)f->Get("h_MassXiC");
        h_MassOmm = (TH1 *)f->Get("h_MassOmm");
        h_MassOmp = (TH1 *)f->Get("h_MassOmp");
        h_MassOmC = (TH1 *)f->Get("h_MassOmC");
        h_multBinEntries_Xi = (TH1 *)f->Get("h_multBinEntries_Xi");
        h_multBinEntries_Om = (TH1 *)f->Get("h_multBinEntries_Om");
        for (Int_t ptBinXi = 0; ptBinXi < nptbins_Xi; ptBinXi++)
        {
            for (Int_t multBinXi = 0; multBinXi < nmultbins_Xi; multBinXi++)
            {
                h_MassXim_pt_mult[ptBinXi][multBinXi] = (TH1 *)f->FindObjectAny(TString::Format(("h_MassXim_pt_mult[%d][%d]"), ptBinXi, multBinXi));
                h_MassXip_pt_mult[ptBinXi][multBinXi] = (TH1 *)f->FindObjectAny(TString::Format(("h_MassXip_pt_mult[%d][%d]"), ptBinXi, multBinXi));
                resultParXip_pt_mult[ptBinXi][multBinXi] = (TH1 *)f->FindObjectAny(TString::Format(("resultParXip_pt_mult[%d][%d]"), ptBinXi, multBinXi));
                resultParXim_pt_mult[ptBinXi][multBinXi] = (TH1 *)f->FindObjectAny(TString::Format(("resultParXim_pt_mult[%d][%d]"), ptBinXi, multBinXi));

                h_MassXiC_pt_mult[ptBinXi][multBinXi] = (TH1 *)f->FindObjectAny(TString::Format(("h_MassXiC_pt_mult[%d][%d]"), ptBinXi, multBinXi));
                resultParXiC_pt_mult[ptBinXi][multBinXi] = (TH1 *)f->FindObjectAny(TString::Format(("resultParXiC_pt_mult[%d][%d]"), ptBinXi, multBinXi));

                resultParXip_pt[ptBinXi] = (TH1 *)f->FindObjectAny(TString::Format(("resultParXip_pt[%d]"), ptBinXi));
                resultParXim_pt[ptBinXi] = (TH1 *)f->FindObjectAny(TString::Format(("resultParXim_pt[%d]"), ptBinXi));
                resultParXiC_pt[ptBinXi] = (TH1 *)f->FindObjectAny(TString::Format(("resultParXiC_pt[%d]"), ptBinXi));
            }
        }
        for (Int_t ptBinOm = 0; ptBinOm < nptbins_Om; ptBinOm++)
        {
            for (Int_t multBinOm = 0; multBinOm < nmultbins_Om; multBinOm++)
            {
                h_MassOmm_pt_mult[ptBinOm][multBinOm] = (TH1 *)f->FindObjectAny(TString::Format(("h_MassOmm_pt_mult[%d][%d]"), ptBinOm, multBinOm));
                h_MassOmp_pt_mult[ptBinOm][multBinOm] = (TH1 *)f->FindObjectAny(TString::Format(("h_MassOmp_pt_mult[%d][%d]"), ptBinOm, multBinOm));
                resultParOmp_pt_mult[ptBinOm][multBinOm] = (TH1 *)f->FindObjectAny(TString::Format(("resultParOmp_pt_mult[%d][%d]"), ptBinOm, multBinOm));
                resultParOmm_pt_mult[ptBinOm][multBinOm] = (TH1 *)f->FindObjectAny(TString::Format(("resultParOmm_pt_mult[%d][%d]"), ptBinOm, multBinOm));

                h_MassOmC_pt_mult[ptBinOm][multBinOm] = (TH1 *)f->FindObjectAny(TString::Format(("h_MassOmC_pt_mult[%d][%d]"), ptBinOm, multBinOm));
                resultParOmC_pt_mult[ptBinOm][multBinOm] = (TH1 *)f->FindObjectAny(TString::Format(("resultParOmC_pt_mult[%d][%d]"), ptBinOm, multBinOm));

                resultParOmp_pt[ptBinOm] = (TH1 *)f->FindObjectAny(TString::Format(("resultParOmp_pt[%d]"), ptBinOm));
                resultParOmm_pt[ptBinOm] = (TH1 *)f->FindObjectAny(TString::Format(("resultParOmm_pt[%d]"), ptBinOm));
                resultParOmC_pt[ptBinOm] = (TH1 *)f->FindObjectAny(TString::Format(("resultParOmC_pt[%d]"), ptBinOm));
            }
        }

        f->Close();
    } /// Input ended!

    /// Saving output :
    TString output = outputFilename;
    if (outputFilename.IsNull())
    {
        output = input;
        output.ReplaceAll("AliCascadeAnalysisMC_Ishaan_Fitting.root", "MCDraw.root");
    }
    Printf("Saving output to '%s' ...", output.Data());
    TFile *out = TFile::Open(output.Data(), "RECREATE");
    if (!out)
    {
        Printf("Error: Cannot open file '%s' !", output.Data());
        return 3;
    }

    out->mkdir("dirRawPt_xim");
    out->mkdir("dirRawPt_xip");
    out->mkdir("dirRawPt_omm");
    out->mkdir("dirRawPt_omp");

    out->mkdir("dirRawPt_xiC");
    out->mkdir("dirRawPt_omC");
    out->mkdir("dirEffPt");
    /// Output set.

    /// Begin
    Int_t markerStyles[] = {4, 21, 22, 23, 29, 33, 34, 43, 47, 41, 20}; // chosen marker style palette
    if (fisMC)
    {
        out->cd("dirEffPt");

        for (Int_t ptBinXi = 0; ptBinXi < nptbins_Xi; ptBinXi++)
        {
            eff_pt_xim->SetBinContent(ptBinXi + 1, resultParXim_pt[ptBinXi]->GetBinContent(4));
            eff_pt_xip->SetBinContent(ptBinXi + 1, resultParXip_pt[ptBinXi]->GetBinContent(4));
            eff_pt_xiC->SetBinContent(ptBinXi + 1, resultParXiC_pt[ptBinXi]->GetBinContent(4));

            eff_pt_xim->SetBinError(ptBinXi + 1, resultParXim_pt[ptBinXi]->GetBinError(4));
            eff_pt_xip->SetBinError(ptBinXi + 1, resultParXip_pt[ptBinXi]->GetBinError(4));
            eff_pt_xiC->SetBinError(ptBinXi + 1, resultParXiC_pt[ptBinXi]->GetBinError(4));
        }

        eff_pt_xim->SetMarkerStyle(markerStyles[10]);
        eff_pt_xim->Write();
        eff_pt_xim->SetName(TString::Format(("Mult: %.0f-%.0f%%"), multbins_Xi[0], multbins_Xi[nmultbins_Xi]));
        hs_xim_eff->Add(eff_pt_xim);

        eff_pt_xip->SetMarkerStyle(markerStyles[10]);
        eff_pt_xip->Write();
        eff_pt_xip->SetName(TString::Format(("Mult: %.0f-%.0f%%"), multbins_Xi[0], multbins_Xi[nmultbins_Xi]));
        hs_xip_eff->Add(eff_pt_xip);

        eff_pt_xiC->SetMarkerStyle(markerStyles[10]);
        eff_pt_xiC->Write();
        eff_pt_xiC->SetName(TString::Format(("Mult: %.0f-%.0f%%"), multbins_Xi[0], multbins_Xi[nmultbins_Xi]));
        hs_xiC_eff->Add(eff_pt_xiC);

        for (Int_t ptBinOm = 0; ptBinOm < nptbins_Om; ptBinOm++)
        {
            eff_pt_omm->SetBinContent(ptBinOm + 1, resultParOmm_pt[ptBinOm]->GetBinContent(4));
            eff_pt_omp->SetBinContent(ptBinOm + 1, resultParOmp_pt[ptBinOm]->GetBinContent(4));
            eff_pt_omC->SetBinContent(ptBinOm + 1, resultParOmC_pt[ptBinOm]->GetBinContent(4));

            eff_pt_omm->SetBinError(ptBinOm + 1, resultParOmm_pt[ptBinOm]->GetBinError(4));
            eff_pt_omp->SetBinError(ptBinOm + 1, resultParOmp_pt[ptBinOm]->GetBinError(4));
            eff_pt_omC->SetBinError(ptBinOm + 1, resultParOmC_pt[ptBinOm]->GetBinError(4));
        }

        eff_pt_omm->SetMarkerStyle(markerStyles[10]);
        eff_pt_omm->Write();
        eff_pt_omm->SetName(TString::Format(("Mult: %.0f-%.0f%%"), multbins_Om[0], multbins_Om[nmultbins_Om]));
        hs_omm_eff->Add(eff_pt_omm);

        eff_pt_omp->SetMarkerStyle(markerStyles[10]);
        eff_pt_omp->Write();
        eff_pt_omp->SetName(TString::Format(("Mult: %.0f-%.0f%%"), multbins_Om[0], multbins_Om[nmultbins_Om]));
        hs_omp_eff->Add(eff_pt_omp);

        eff_pt_omC->SetMarkerStyle(markerStyles[10]);
        eff_pt_omC->Write();
        eff_pt_omC->SetName(TString::Format(("Mult: %.0f-%.0f%%"), multbins_Om[0], multbins_Om[nmultbins_Om]));
        hs_omC_eff->Add(eff_pt_omC);
    }
    /// Begin differential part

    for (Int_t multBinXi = 0; multBinXi < nmultbins_Xi; multBinXi++)
    {
        rawPt_xim[multBinXi] = new TH1D(TString::Format(("rawPt_xim[%d]"), multBinXi), "", nptbins_Xi, ptbins_Xi);
        rawPt_xip[multBinXi] = new TH1D(TString::Format(("rawPt_xip[%d]"), multBinXi), "", nptbins_Xi, ptbins_Xi);

        rawPt_xiC[multBinXi] = new TH1D(TString::Format(("rawPt_xiC[%d]"), multBinXi), "", nptbins_Xi, ptbins_Xi);

        if (fisMC)
        {
            eff_xim[multBinXi] = new TH1D(TString::Format(("eff_xim[%d]"), multBinXi), "", nptbins_Xi, ptbins_Xi);
            eff_xip[multBinXi] = new TH1D(TString::Format(("eff_xip[%d]"), multBinXi), "", nptbins_Xi, ptbins_Xi);

            eff_xiC[multBinXi] = new TH1D(TString::Format(("eff_xiC[%d]"), multBinXi), "", nptbins_Xi, ptbins_Xi);
        }
        // redChi2_xim[multBinXi] = new TH1D(TString::Format(("redChi2_xim[%d]"), multBinXi), TString::Format(("Mult: %.0f-%.0f%%"), multbins_Xi[multBinXi ], multbins_Xi[multBinXi+1]), nptbins_Xi, ptbins_Xi);
        // redChi2_xip[multBinXi] = new TH1D(TString::Format(("redChi2_xip[%d]"), multBinXi), TString::Format(("Mult: %.0f-%.0f%%"), multbins_Xi[multBinXi ], multbins_Xi[multBinXi+1]), nptbins_Xi, ptbins_Xi);
        for (Int_t ptBinXi = 0; ptBinXi < nptbins_Xi; ptBinXi++)
        {
            rawPt_xim[multBinXi]->SetBinContent(ptBinXi + 1, resultParXim_pt_mult[ptBinXi][multBinXi]->GetBinContent(1) / ((rawPt_xim[multBinXi]->GetBinWidth(ptBinXi + 1)) * (h_multBinEntries_Xi->GetBinContent(multBinXi + 1)))); // Bin 1 in resultparams is raw pt's bin counting
            rawPt_xim[multBinXi]->SetBinError(ptBinXi + 1, resultParXim_pt_mult[ptBinXi][multBinXi]->GetBinError(1) / ((rawPt_xim[multBinXi]->GetBinWidth(ptBinXi + 1)) * (h_multBinEntries_Xi->GetBinContent(multBinXi + 1))));

            rawPt_xip[multBinXi]->SetBinContent(ptBinXi + 1, resultParXip_pt_mult[ptBinXi][multBinXi]->GetBinContent(1) / ((rawPt_xip[multBinXi]->GetBinWidth(ptBinXi + 1)) * (h_multBinEntries_Xi->GetBinContent(multBinXi + 1)))); // Bin 1 in resultparams is raw pt's bin counting
            rawPt_xip[multBinXi]->SetBinError(ptBinXi + 1, resultParXip_pt_mult[ptBinXi][multBinXi]->GetBinError(1) / ((rawPt_xip[multBinXi]->GetBinWidth(ptBinXi + 1)) * (h_multBinEntries_Xi->GetBinContent(multBinXi + 1))));

            rawPt_xiC[multBinXi]->SetBinContent(ptBinXi + 1, resultParXiC_pt_mult[ptBinXi][multBinXi]->GetBinContent(1) / ((rawPt_xiC[multBinXi]->GetBinWidth(ptBinXi + 1)) * (h_multBinEntries_Xi->GetBinContent(multBinXi + 1)))); // Bin 1 in resultparams is raw pt's bin counting
            rawPt_xiC[multBinXi]->SetBinError(ptBinXi + 1, resultParXiC_pt_mult[ptBinXi][multBinXi]->GetBinError(1) / ((rawPt_xiC[multBinXi]->GetBinWidth(ptBinXi + 1)) * (h_multBinEntries_Xi->GetBinContent(multBinXi + 1))));

            if (fisMC)
            {
                eff_xim[multBinXi]->SetBinContent(ptBinXi + 1, resultParXim_pt_mult[ptBinXi][multBinXi]->GetBinContent(4)); // Bin 4 in resultparams is MC efficiency
                eff_xim[multBinXi]->SetBinError(ptBinXi + 1, resultParXim_pt_mult[ptBinXi][multBinXi]->GetBinError(4));

                eff_xip[multBinXi]->SetBinContent(ptBinXi + 1, resultParXip_pt_mult[ptBinXi][multBinXi]->GetBinContent(4)); // Bin 4 in resultparams is MC efficiency
                eff_xip[multBinXi]->SetBinError(ptBinXi + 1, resultParXip_pt_mult[ptBinXi][multBinXi]->GetBinError(4));

                eff_xiC[multBinXi]->SetBinContent(ptBinXi + 1, resultParXiC_pt_mult[ptBinXi][multBinXi]->GetBinContent(4)); // Bin 4 in resultparams is MC efficiency
                eff_xiC[multBinXi]->SetBinError(ptBinXi + 1, resultParXiC_pt_mult[ptBinXi][multBinXi]->GetBinError(4));
            }

            /// Scale for bin width: N->dN/dpt
            // rawPt_xim[multBinXi]->Scale(1, "width");
            // rawPt_xip[multBinXi]->Scale(1, "width");

            h_MassXim_pt_mult[ptBinXi][multBinXi]->GetXaxis()->SetTitle("#Lambda^{0}-#pi^{#pm} Inv. Mass (GeV/c^{2})");
            h_MassXip_pt_mult[ptBinXi][multBinXi]->GetXaxis()->SetTitle("#Lambda^{0}-#pi^{#pm} Inv. Mass (GeV/c^{2})");
            h_MassXiC_pt_mult[ptBinXi][multBinXi]->GetXaxis()->SetTitle("#Lambda^{0}-#pi^{#pm} Inv. Mass (GeV/c^{2})");

            TF1 *f1 = (TF1 *)h_MassXim_pt_mult[ptBinXi][multBinXi]->GetListOfFunctions()->At(0); /// 0 position = peak Function (GausPol2)
            DrawAndSave(h_MassXim_pt_mult[ptBinXi][multBinXi], f1->GetParameter(1), f1->GetParameter(2), saveImages, outputFolder, imageFormat);
            TF1 *f2 = (TF1 *)h_MassXip_pt_mult[ptBinXi][multBinXi]->GetListOfFunctions()->At(0); /// peak Function (GausPol2)
            DrawAndSave(h_MassXip_pt_mult[ptBinXi][multBinXi], f2->GetParameter(1), f2->GetParameter(2), saveImages, outputFolder, imageFormat);
            TF1 *f3 = (TF1 *)h_MassXiC_pt_mult[ptBinXi][multBinXi]->GetListOfFunctions()->At(0); /// peak Function (GausPol2)
            DrawAndSave(h_MassXiC_pt_mult[ptBinXi][multBinXi], f3->GetParameter(1), f3->GetParameter(2), saveImages, outputFolder, imageFormat);
        }
        out->cd("dirRawPt_xip");
        rawPt_xip[multBinXi]->SetMarkerStyle(markerStyles[multBinXi]);
        rawPt_xip[multBinXi]->Write();
        rawPt_xip[multBinXi]->SetName(TString::Format(("Mult: %.0f-%.0f%%, #times2^{%d}"), multbins_Xi[multBinXi], multbins_Xi[multBinXi + 1], (nmultbins_Xi - 1) - multBinXi));
        rawPt_xip[multBinXi]->Scale(pow(2, (nmultbins_Xi - 1) - multBinXi));
        hs_xip->Add(rawPt_xip[multBinXi]);

        if (fisMC)
        {
            eff_xip[multBinXi]->SetMarkerStyle(markerStyles[multBinXi]);
            eff_xip[multBinXi]->Write();
            eff_xip[multBinXi]->SetName(TString::Format(("Mult: %.0f-%.0f%%"), multbins_Xi[multBinXi], multbins_Xi[multBinXi + 1]));
            // eff_xip[multBinXi]->Scale(pow(2, (nmultbins_Xi - 1) - multBinXi));
            hs_xip_eff->Add(eff_xip[multBinXi]);
        }

        out->cd("dirRawPt_xim");
        rawPt_xim[multBinXi]->SetMarkerStyle(markerStyles[multBinXi]);
        rawPt_xim[multBinXi]->Write();
        rawPt_xim[multBinXi]->SetName(TString::Format(("Mult: %.0f-%.0f%%, #times2^{%d}"), multbins_Xi[multBinXi], multbins_Xi[multBinXi + 1], (nmultbins_Xi - 1) - multBinXi));
        rawPt_xim[multBinXi]->Scale(pow(2, (nmultbins_Xi - 1) - multBinXi));
        hs_xim->Add(rawPt_xim[multBinXi]);

        if (fisMC)
        {
            eff_xim[multBinXi]->SetMarkerStyle(markerStyles[multBinXi]);
            eff_xim[multBinXi]->Write();
            eff_xim[multBinXi]->SetName(TString::Format(("Mult: %.0f-%.0f%%"), multbins_Xi[multBinXi], multbins_Xi[multBinXi + 1]));
            // eff_xim[multBinXi]->Scale(pow(2, (nmultbins_Xi - 1) - multBinXi));
            hs_xim_eff->Add(eff_xim[multBinXi]);
        }
        out->cd("dirRawPt_xiC");
        rawPt_xiC[multBinXi]->SetMarkerStyle(markerStyles[multBinXi]);
        rawPt_xiC[multBinXi]->Write();
        rawPt_xiC[multBinXi]->SetName(TString::Format(("Mult: %.0f-%.0f%%, #times2^{%d}"), multbins_Xi[multBinXi], multbins_Xi[multBinXi + 1], (nmultbins_Xi - 1) - multBinXi));
        rawPt_xiC[multBinXi]->Scale(pow(2, (nmultbins_Xi - 1) - multBinXi));
        hs_xiC->Add(rawPt_xiC[multBinXi]);

        if (fisMC)
        {
            eff_xiC[multBinXi]->SetMarkerStyle(markerStyles[multBinXi]);
            eff_xiC[multBinXi]->Write();
            eff_xiC[multBinXi]->SetName(TString::Format(("Mult: %.0f-%.0f%%"), multbins_Xi[multBinXi], multbins_Xi[multBinXi + 1]));
            // eff_xiC[multBinXi]->Scale(pow(2, (nmultbins_Xi - 1) - multBinXi));
            hs_xiC_eff->Add(eff_xiC[multBinXi]);

            eff_xiC_ratio[multBinXi] = new TH1D(TString::Format(("eff_xiC_ratio[%d]"), multBinXi), "", nptbins_Xi, ptbins_Xi);
            eff_xiC_ratio[multBinXi]->Divide(eff_xiC[multBinXi], eff_pt_xiC);
            eff_xiC_ratio[multBinXi]->SetMarkerStyle(markerStyles[multBinXi]);
            eff_xiC_ratio[multBinXi]->Write();
            eff_xiC_ratio[multBinXi]->SetName(TString::Format(("Mult: %.0f-%.0f%%"), multbins_Xi[multBinXi], multbins_Xi[multBinXi + 1]));
            hs_xiC_eff_ratio->Add(eff_xiC_ratio[multBinXi]);
        }
    }
    for (Int_t multBinOm = 0; multBinOm < nmultbins_Om; multBinOm++)
    {
        rawPt_omm[multBinOm] = new TH1D(TString::Format(("rawPt_omm[%d]"), multBinOm), TString::Format(("Mult: %.0f-%.0f%%, #times2^{%d}"), multbins_Om[multBinOm], multbins_Om[multBinOm + 1], (nmultbins_Om - 1) - multBinOm), nptbins_Om, ptbins_Om);
        rawPt_omp[multBinOm] = new TH1D(TString::Format(("rawPt_omp[%d]"), multBinOm), TString::Format(("Mult: %.0f-%.0f%%, #times2^{%d}"), multbins_Om[multBinOm], multbins_Om[multBinOm + 1], (nmultbins_Om - 1) - multBinOm), nptbins_Om, ptbins_Om);

        rawPt_omC[multBinOm] = new TH1D(TString::Format(("rawPt_omC[%d]"), multBinOm), TString::Format(("Mult: %.0f-%.0f%%, #times2^{%d}"), multbins_Om[multBinOm], multbins_Om[multBinOm + 1], (nmultbins_Om - 1) - multBinOm), nptbins_Om, ptbins_Om);

        if (fisMC)
        {
            eff_omm[multBinOm] = new TH1D(TString::Format(("eff_omm[%d]"), multBinOm), TString::Format(("Mult: %.0f-%.0f%%"), multbins_Om[multBinOm], multbins_Om[multBinOm + 1]), nptbins_Om, ptbins_Om);
            eff_omp[multBinOm] = new TH1D(TString::Format(("eff_omp[%d]"), multBinOm), TString::Format(("Mult: %.0f-%.0f%%"), multbins_Om[multBinOm], multbins_Om[multBinOm + 1]), nptbins_Om, ptbins_Om);

            eff_omC[multBinOm] = new TH1D(TString::Format(("eff_omC[%d]"), multBinOm), TString::Format(("Mult: %.0f-%.0f%%"), multbins_Om[multBinOm], multbins_Om[multBinOm + 1]), nptbins_Om, ptbins_Om);
        }
        for (Int_t ptBinOm = 0; ptBinOm < nptbins_Om; ptBinOm++)
        {
            rawPt_omm[multBinOm]->SetBinContent(ptBinOm + 1, resultParOmm_pt_mult[ptBinOm][multBinOm]->GetBinContent(1) / ((rawPt_omm[multBinOm]->GetBinWidth(ptBinOm + 1)) * (h_multBinEntries_Om->GetBinContent(multBinOm + 1)))); // Bin 1 in resultparams is raw pt's bin counting
            rawPt_omm[multBinOm]->SetBinError(ptBinOm + 1, resultParOmm_pt_mult[ptBinOm][multBinOm]->GetBinError(1) / ((rawPt_omm[multBinOm]->GetBinWidth(ptBinOm + 1)) * (h_multBinEntries_Om->GetBinContent(multBinOm + 1))));
            rawPt_omp[multBinOm]->SetBinContent(ptBinOm + 1, resultParOmp_pt_mult[ptBinOm][multBinOm]->GetBinContent(1) / ((rawPt_omp[multBinOm]->GetBinWidth(ptBinOm + 1)) * (h_multBinEntries_Om->GetBinContent(multBinOm + 1)))); // Bin 1 in resultparams is raw pt's bin counting
            rawPt_omp[multBinOm]->SetBinError(ptBinOm + 1, resultParOmp_pt_mult[ptBinOm][multBinOm]->GetBinError(1) / ((rawPt_omp[multBinOm]->GetBinWidth(ptBinOm + 1)) * (h_multBinEntries_Om->GetBinContent(multBinOm + 1))));

            rawPt_omC[multBinOm]->SetBinContent(ptBinOm + 1, resultParOmC_pt_mult[ptBinOm][multBinOm]->GetBinContent(1) / ((rawPt_omC[multBinOm]->GetBinWidth(ptBinOm + 1)) * (h_multBinEntries_Om->GetBinContent(multBinOm + 1)))); // Bin 1 in resultparams is raw pt's bin counting
            rawPt_omC[multBinOm]->SetBinError(ptBinOm + 1, resultParOmC_pt_mult[ptBinOm][multBinOm]->GetBinError(1) / ((rawPt_omC[multBinOm]->GetBinWidth(ptBinOm + 1)) * (h_multBinEntries_Om->GetBinContent(multBinOm + 1))));

            if (fisMC)
            {
                eff_omm[multBinOm]->SetBinContent(ptBinOm + 1, resultParOmm_pt_mult[ptBinOm][multBinOm]->GetBinContent(4)); // Bin 4 in resultparams is MC efficiency
                eff_omm[multBinOm]->SetBinError(ptBinOm + 1, resultParOmm_pt_mult[ptBinOm][multBinOm]->GetBinError(4));

                eff_omp[multBinOm]->SetBinContent(ptBinOm + 1, resultParOmp_pt_mult[ptBinOm][multBinOm]->GetBinContent(4)); // Bin 4 in resultparams is MC efficiency
                eff_omp[multBinOm]->SetBinError(ptBinOm + 1, resultParOmp_pt_mult[ptBinOm][multBinOm]->GetBinError(4));

                eff_omC[multBinOm]->SetBinContent(ptBinOm + 1, resultParOmC_pt_mult[ptBinOm][multBinOm]->GetBinContent(4)); // Bin 4 in resultparams is MC efficiency
                eff_omC[multBinOm]->SetBinError(ptBinOm + 1, resultParOmC_pt_mult[ptBinOm][multBinOm]->GetBinError(4));
            }
            /// Scale for bin width: N->dN/dpt
            // rawPt_omm[multBinOm]->Scale(1, "width");
            // rawPt_omp[multBinOm]->Scale(1, "width");

            h_MassOmm_pt_mult[ptBinOm][multBinOm]->GetXaxis()->SetTitle("#Lambda^{0}-K^{#pm} Inv. Mass (GeV/c^{2})");
            h_MassOmp_pt_mult[ptBinOm][multBinOm]->GetXaxis()->SetTitle("#Lambda^{0}-K^{#pm} Inv. Mass (GeV/c^{2})");
            h_MassOmC_pt_mult[ptBinOm][multBinOm]->GetXaxis()->SetTitle("#Lambda^{0}-K^{#pm} Inv. Mass (GeV/c^{2})");

            TF1 *f1 = (TF1 *)h_MassOmm_pt_mult[ptBinOm][multBinOm]->GetListOfFunctions()->At(0); /// peak Function (GausPol2)
            DrawAndSave(h_MassOmm_pt_mult[ptBinOm][multBinOm], f1->GetParameter(1), f1->GetParameter(2), saveImages, outputFolder, imageFormat);
            TF1 *f2 = (TF1 *)h_MassOmp_pt_mult[ptBinOm][multBinOm]->GetListOfFunctions()->At(0); /// peak Function (GausPol2)
            DrawAndSave(h_MassOmp_pt_mult[ptBinOm][multBinOm], f2->GetParameter(1), f2->GetParameter(2), saveImages, outputFolder, imageFormat);
            TF1 *f3 = (TF1 *)h_MassOmC_pt_mult[ptBinOm][multBinOm]->GetListOfFunctions()->At(0); /// peak Function (GausPol2)
            DrawAndSave(h_MassOmC_pt_mult[ptBinOm][multBinOm], f3->GetParameter(1), f3->GetParameter(2), saveImages, outputFolder, imageFormat);
        }
        out->cd("dirRawPt_omp");
        rawPt_omp[multBinOm]->SetMarkerStyle(markerStyles[multBinOm]);
        rawPt_omp[multBinOm]->Write();
        rawPt_omp[multBinOm]->SetName(TString::Format(("Mult: %.0f-%.0f%%, #times2^{%d}"), multbins_Om[multBinOm], multbins_Om[multBinOm + 1], (nmultbins_Om - 1) - multBinOm));
        rawPt_omp[multBinOm]->Scale(pow(2, (nmultbins_Om - 1) - multBinOm));
        hs_omp->Add(rawPt_omp[multBinOm]);

        if (fisMC)
        {
            eff_omp[multBinOm]->SetMarkerStyle(markerStyles[multBinOm]);
            eff_omp[multBinOm]->Write();
            eff_omp[multBinOm]->SetName(TString::Format(("Mult: %.0f-%.0f%%"), multbins_Om[multBinOm], multbins_Om[multBinOm + 1]));
            // eff_omp[multBinOm]->Scale(pow(2, (nmultbins_Om - 1) - multBinOm));
            hs_omp_eff->Add(eff_omp[multBinOm]);
        }

        out->cd("dirRawPt_omm");
        rawPt_omm[multBinOm]->SetMarkerStyle(markerStyles[multBinOm]);
        rawPt_omm[multBinOm]->Write();
        rawPt_omm[multBinOm]->SetName(TString::Format(("Mult: %.0f-%.0f%%, #times2^{%d}"), multbins_Om[multBinOm], multbins_Om[multBinOm + 1], (nmultbins_Om - 1) - multBinOm));
        rawPt_omm[multBinOm]->Scale(pow(2, (nmultbins_Om - 1) - multBinOm));
        hs_omm->Add(rawPt_omm[multBinOm]);

        if (fisMC)
        {
            eff_omm[multBinOm]->SetMarkerStyle(markerStyles[multBinOm]);
            eff_omm[multBinOm]->Write();
            eff_omm[multBinOm]->SetName(TString::Format(("Mult: %.0f-%.0f%%"), multbins_Om[multBinOm], multbins_Om[multBinOm + 1]));
            // eff_omm[multBinOm]->Scale(pow(2, (nmultbins_Om - 1) - multBinOm));
            hs_omm_eff->Add(eff_omm[multBinOm]);
        }

        out->cd("dirRawPt_omC");
        rawPt_omC[multBinOm]->SetMarkerStyle(markerStyles[multBinOm]);
        rawPt_omC[multBinOm]->Write();
        rawPt_omC[multBinOm]->SetName(TString::Format(("Mult: %.0f-%.0f%%, #times2^{%d}"), multbins_Om[multBinOm], multbins_Om[multBinOm + 1], (nmultbins_Om - 1) - multBinOm));
        rawPt_omC[multBinOm]->Scale(pow(2, (nmultbins_Om - 1) - multBinOm));
        hs_omC->Add(rawPt_omC[multBinOm]);

        if (fisMC)
        {
            eff_omC[multBinOm]->SetMarkerStyle(markerStyles[multBinOm]);
            eff_omC[multBinOm]->Write();
            eff_omC[multBinOm]->SetName(TString::Format(("Mult: %.0f-%.0f%%"), multbins_Om[multBinOm], multbins_Om[multBinOm + 1]));
            // eff_omC[multBinOm]->Scale(pow(2, (nmultbins_Om - 1) - multBinOm));
            hs_omC_eff->Add(eff_omC[multBinOm]);

            eff_omC_ratio[multBinOm] = new TH1D(TString::Format(("eff_omC_ratio[%d]"), multBinOm), "", nptbins_Om, ptbins_Om);
            eff_omC_ratio[multBinOm]->Divide(eff_omC[multBinOm], eff_pt_omC);
            eff_omC_ratio[multBinOm]->SetMarkerStyle(markerStyles[multBinOm]);
            eff_omC_ratio[multBinOm]->Write();
            eff_omC_ratio[multBinOm]->SetName(TString::Format(("Mult: %.0f-%.0f%%"), multbins_Om[multBinOm], multbins_Om[multBinOm + 1]));
            hs_omC_eff_ratio->Add(eff_omC_ratio[multBinOm]);
        }
    }

    gStyle->SetOptStat(0);
    gStyle->SetPalette(kVisibleSpectrum);

    TCanvas *c1 = new TCanvas("c1", "c1", 1920, 1080);
    TCanvas *c2 = new TCanvas("c2", "c2", 1920, 1080);
    TCanvas *c3 = new TCanvas("c3", "c3", 1920, 1080);
    TCanvas *c4 = new TCanvas("c4", "c4", 1920, 1080);
    TCanvas *c5 = new TCanvas("c5", "c5", 1920, 1080);
    TCanvas *c6 = new TCanvas("c6", "c6", 1920, 1080);
    c1->cd();
    hs_xip->Draw("plc pmc nostack");
    gPad->SetLogy();
    gPad->BuildLegend(0.9, 0.6, 1., 1., "");
    hs_xip->GetXaxis()->SetTitle("#it{p}_{T} (GeV/c)");
    hs_xip->GetYaxis()->SetTitle("#frac{1}{#it{N}_{inel}} #frac{d#it{N}}{d#it{p}_{T}}");
    c1->Modified();
    c1->ForceUpdate();

    c2->cd();
    hs_omp->Draw("plc pmc nostack");
    gPad->SetLogy();
    gPad->BuildLegend(0.9, 0.6, 1., 1., "");
    hs_omp->GetXaxis()->SetTitle("#it{p}_{T} (GeV/c)");
    hs_omp->GetYaxis()->SetTitle("#frac{1}{#it{N}_{inel}} #frac{d#it{N}}{d#it{p}_{T}}");
    c2->Modified();
    c2->ForceUpdate();

    c3->cd();
    hs_xim->Draw("plc pmc nostack");
    gPad->SetLogy();
    gPad->BuildLegend(0.9, 0.6, 1., 1., "");
    hs_xim->GetXaxis()->SetTitle("#it{p}_{T} (GeV/c)");
    hs_xim->GetYaxis()->SetTitle("#frac{1}{#it{N}_{inel}} #frac{d#it{N}}{d#it{p}_{T}}");
    c3->Modified();
    c3->ForceUpdate();

    c4->cd();
    hs_omm->Draw("plc pmc nostack");
    gPad->SetLogy();
    gPad->BuildLegend(0.9, 0.6, 1., 1., "");
    hs_omm->GetXaxis()->SetTitle("#it{p}_{T} (GeV/c)");
    hs_omm->GetYaxis()->SetTitle("#frac{1}{#it{N}_{inel}} #frac{d#it{N}}{d#it{p}_{T}}");
    c4->Modified();
    c4->ForceUpdate();

    c5->cd();
    hs_xiC->Draw("plc pmc nostack");
    gPad->SetLogy();
    gPad->BuildLegend(0.9, 0.6, 1., 1., "");
    hs_xiC->GetXaxis()->SetTitle("#it{p}_{T} (GeV/c)");
    hs_xiC->GetYaxis()->SetTitle("#frac{1}{#it{N}_{inel}} #frac{d#it{N}}{d#it{p}_{T}}");
    c5->Modified();
    c5->ForceUpdate();

    c6->cd();
    hs_omC->Draw("plc pmc nostack");
    gPad->SetLogy();
    gPad->BuildLegend(0.9, 0.6, 1., 1., "");
    hs_omC->GetXaxis()->SetTitle("#it{p}_{T} (GeV/c)");
    hs_omC->GetYaxis()->SetTitle("#frac{1}{#it{N}_{inel}} #frac{d#it{N}}{d#it{p}_{T}}");
    c6->Modified();
    c6->ForceUpdate();

    c1->Modified();
    c1->ForceUpdate();
    c2->Modified();
    c2->ForceUpdate();
    c3->Modified();
    c3->ForceUpdate();
    c4->Modified();
    c4->ForceUpdate();
    c5->Modified();
    c5->ForceUpdate();
    c6->Modified();
    c6->ForceUpdate();

    out->cd();

    hs_xip->Write();
    hs_xim->Write();
    hs_omp->Write();
    hs_omm->Write();
    hs_xiC->Write();
    hs_omC->Write();

    if (saveStack)
    {
        c1->cd();
        SaveImage(Form("%s/images/rawPt", outputFolder.Data()), "xipN", imageFormat.Data());

        c2->cd();
        SaveImage(Form("%s/images/rawPt", outputFolder.Data()), "ompN", imageFormat.Data());

        c3->cd();
        SaveImage(Form("%s/images/rawPt", outputFolder.Data()), "ximN", imageFormat.Data());

        c4->cd();
        SaveImage(Form("%s/images/rawPt", outputFolder.Data()), "ommN", imageFormat.Data());

        c5->cd();
        SaveImage(Form("%s/images/rawPt", outputFolder.Data()), "xicN", imageFormat.Data());

        c6->cd();
        SaveImage(Form("%s/images/rawPt", outputFolder.Data()), "omcN", imageFormat.Data());
    }

    if (fisMC)
    {
        TCanvas *e1 = new TCanvas("e1", "e1", 1920, 1080);
        TCanvas *e2 = new TCanvas("e2", "e2", 1920, 1080);
        TCanvas *e3 = new TCanvas("e3", "e3", 1920, 1080);
        TCanvas *e4 = new TCanvas("e4", "e4", 1920, 1080);
        TCanvas *e5 = new TCanvas("e5", "e5", 1920, 1080);
        TCanvas *e6 = new TCanvas("e6", "e6", 1920, 1080);

        TCanvas *r1 = new TCanvas("r1", "r1", 1920, 1080);
        TCanvas *r2 = new TCanvas("r2", "r2", 1920, 1080);

        e1->cd();
        hs_xip_eff->Draw("plc pmc nostack");
        gPad->BuildLegend(0.9, 0.6, 1., 1., "");
        hs_xip_eff->GetXaxis()->SetTitle("#it{p}_{T} (GeV/c)");
        hs_xip_eff->GetYaxis()->SetTitle("Efficiency");
        e1->Modified();
        e1->ForceUpdate();

        e2->cd();
        hs_omp_eff->Draw("plc pmc nostack");
        gPad->BuildLegend(0.9, 0.6, 1., 1., "");
        hs_omp_eff->GetXaxis()->SetTitle("#it{p}_{T} (GeV/c)");
        hs_omp_eff->GetYaxis()->SetTitle("Efficiency");
        e2->Modified();
        e2->ForceUpdate();

        e3->cd();
        hs_xim_eff->Draw("plc pmc nostack");
        gPad->BuildLegend(0.9, 0.6, 1., 1., "");
        hs_xim_eff->GetXaxis()->SetTitle("#it{p}_{T} (GeV/c)");
        hs_xim_eff->GetYaxis()->SetTitle("Efficiency");
        e3->Modified();
        e3->ForceUpdate();

        e4->cd();
        hs_omm_eff->Draw("plc pmc nostack");
        gPad->BuildLegend(0.9, 0.6, 1., 1., "");
        hs_omm_eff->GetXaxis()->SetTitle("#it{p}_{T} (GeV/c)");
        hs_omm_eff->GetYaxis()->SetTitle("Efficiency");
        e4->Modified();
        e4->ForceUpdate();

        e5->cd();
        hs_xiC_eff->Draw("plc pmc nostack");
        gPad->BuildLegend(0.9, 0.6, 1., 1., "");
        hs_xiC_eff->GetXaxis()->SetTitle("#it{p}_{T} (GeV/c)");
        hs_xiC_eff->GetYaxis()->SetTitle("Efficiency");
        e5->Modified();
        e5->ForceUpdate();

        e6->cd();
        hs_omC_eff->Draw("plc pmc nostack");
        gPad->BuildLegend(0.9, 0.6, 1., 1., "");
        hs_omC_eff->GetXaxis()->SetTitle("#it{p}_{T} (GeV/c)");
        hs_omC_eff->GetYaxis()->SetTitle("Efficiency");
        e6->Modified();
        e6->ForceUpdate();

        r1->cd();
        hs_xiC_eff_ratio->Draw("plc pmc nostack");
        gPad->BuildLegend(0.9, 0.6, 1., 1., "");
        hs_xiC_eff_ratio->GetXaxis()->SetTitle("#it{p}_{T} (GeV/c)");
        hs_xiC_eff_ratio->GetYaxis()->SetTitle("Multiplicity classes/0-100%");
        r1->Modified();
        r1->ForceUpdate();

        r2->cd();
        hs_omC_eff_ratio->Draw("plc pmc nostack");
        gPad->BuildLegend(0.9, 0.6, 1., 1., "");
        hs_omC_eff_ratio->GetXaxis()->SetTitle("#it{p}_{T} (GeV/c)");
        hs_omC_eff_ratio->GetYaxis()->SetTitle("Multiplicity classes/0-100%");
        r2->Modified();
        r2->ForceUpdate();

        hs_xip->GetXaxis()->SetTitle("#it{p}_{T} (GeV/c)");
        hs_xim->GetXaxis()->SetTitle("#it{p}_{T} (GeV/c)");
        hs_omp->GetXaxis()->SetTitle("#it{p}_{T} (GeV/c)");
        hs_omm->GetXaxis()->SetTitle("#it{p}_{T} (GeV/c)");
        hs_xiC->GetXaxis()->SetTitle("#it{p}_{T} (GeV/c)");
        hs_omC->GetXaxis()->SetTitle("#it{p}_{T} (GeV/c)");

        hs_xip->GetYaxis()->SetTitle("#frac{1}{#it{N}_{inel}} #frac{d#it{N}}{d#it{p}_{T}}");
        hs_xim->GetYaxis()->SetTitle("#frac{1}{#it{N}_{inel}} #frac{d#it{N}}{d#it{p}_{T}}");
        hs_omp->GetYaxis()->SetTitle("#frac{1}{#it{N}_{inel}} #frac{d#it{N}}{d#it{p}_{T}}");
        hs_omm->GetYaxis()->SetTitle("#frac{1}{#it{N}_{inel}} #frac{d#it{N}}{d#it{p}_{T}}");
        hs_xiC->GetYaxis()->SetTitle("#frac{1}{#it{N}_{inel}} #frac{d#it{N}}{d#it{p}_{T}}");
        hs_omC->GetYaxis()->SetTitle("#frac{1}{#it{N}_{inel}} #frac{d#it{N}}{d#it{p}_{T}}");

        hs_xip_eff->GetXaxis()->SetTitle("#it{p}_{T} (GeV/c)");
        hs_xim_eff->GetXaxis()->SetTitle("#it{p}_{T} (GeV/c)");
        hs_omp_eff->GetXaxis()->SetTitle("#it{p}_{T} (GeV/c)");
        hs_omm_eff->GetXaxis()->SetTitle("#it{p}_{T} (GeV/c)");
        hs_xiC_eff->GetXaxis()->SetTitle("#it{p}_{T} (GeV/c)");
        hs_omC_eff->GetXaxis()->SetTitle("#it{p}_{T} (GeV/c)");

        hs_xip_eff->GetYaxis()->SetTitle("Efficiency");
        hs_xim_eff->GetYaxis()->SetTitle("Efficiency");
        hs_omp_eff->GetYaxis()->SetTitle("Efficiency");
        hs_omm_eff->GetYaxis()->SetTitle("Efficiency");
        hs_xiC_eff->GetYaxis()->SetTitle("Efficiency");
        hs_omC_eff->GetYaxis()->SetTitle("Efficiency");

        hs_xiC_eff_ratio->GetXaxis()->SetTitle("#it{p}_{T} (GeV/c)");
        hs_omC_eff_ratio->GetXaxis()->SetTitle("#it{p}_{T} (GeV/c)");
        hs_xiC_eff_ratio->GetYaxis()->SetTitle("Multiplicity classes/0-100%");
        hs_omC_eff_ratio->GetYaxis()->SetTitle("Multiplicity classes/0-100%");

        e1->Modified();
        e1->ForceUpdate();
        e2->Modified();
        e2->ForceUpdate();
        e3->Modified();
        e3->ForceUpdate();
        e4->Modified();
        e4->ForceUpdate();
        e5->Modified();
        e5->ForceUpdate();
        e6->Modified();
        e6->ForceUpdate();

        r1->Modified();
        r1->ForceUpdate();
        r2->Modified();
        r2->ForceUpdate();

        hs_xip_eff->Write();
        hs_xim_eff->Write();
        hs_omp_eff->Write();
        hs_omm_eff->Write();
        hs_xiC_eff->Write();
        hs_omC_eff->Write();

        hs_xiC_eff_ratio->Write();
        hs_omC_eff_ratio->Write();

        if (saveStack)
        {
            e1->cd();
            SaveImage(Form("%s/images/eff", outputFolder.Data()), "eff_xipN", imageFormat.Data());
            e2->cd();
            SaveImage(Form("%s/images/eff", outputFolder.Data()), "eff_ompN", imageFormat.Data());

            e3->cd();
            SaveImage(Form("%s/images/eff", outputFolder.Data()), "eff_ximN", imageFormat.Data());

            e4->cd();
            SaveImage(Form("%s/images/eff", outputFolder.Data()), "eff_ommN", imageFormat.Data());

            e5->cd();
            SaveImage(Form("%s/images/eff", outputFolder.Data()), "eff_xicN", imageFormat.Data());

            e6->cd();
            SaveImage(Form("%s/images/eff", outputFolder.Data()), "eff_omcN", imageFormat.Data());

            r1->cd();
            SaveImage(Form("%s/images/effRatio", outputFolder.Data()), "effRatio_xicN", imageFormat.Data());

            r2->cd();
            SaveImage(Form("%s/images/effRatio", outputFolder.Data()), "effRatio_omcN", imageFormat.Data());
        }
    }

    out->Close();
    // delete out;
    return 0;
}
void DrawAndSave(TH1 *peak, Double_t pPosition, Double_t pWidth, Bool_t saveImages, TString outputFolder, TString imageFormat)
{

    gROOT->SetBatch(kTRUE);
    TCanvas *c1 = new TCanvas(peak->GetName(), peak->GetTitle(), 1920, 1080);

    /// Defining peak limits for signal region (green lines):
    ///  par[1] = peak position, par[2] = peak width
    Double_t lPeakLeftLimit = pPosition - 1. * 4 * TMath::Abs(pWidth);
    Double_t lPeakRightLimit = pPosition + 1. * 4 * TMath::Abs(pWidth);
    TLine *lLineLeft = new TLine(lPeakLeftLimit, 0, lPeakLeftLimit, peak->GetMaximum());
    TLine *lLineRight = new TLine(lPeakRightLimit, 0, lPeakRightLimit, peak->GetMaximum());
    lLineLeft->SetLineColor(kMagenta);
    lLineRight->SetLineColor(kMagenta);
    peak->Draw();
    lLineLeft->Draw("same");
    lLineRight->Draw("same");

    if (saveImages)
    {
        TString imageFolder = peak->GetName();
        imageFolder = imageFolder(0, 17); /// getting substring for naming purpose
        SaveImage(Form("%s/images/%s", outputFolder.Data(), imageFolder.Data()), peak->GetName(), imageFormat.Data());
    }
}
