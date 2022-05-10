#include <TString.h>
#include <TSystem.h>

void DrawAndSave(TH1 *peak, Double_t pPosition, Double_t pWidth, Bool_t saveImages, TString outputFolder);

int AliCascadeAnalysis_Ishaan_Draw(std::string input = "AliCascadeAnalysis_Ishaan_Fitting.root", TString outputFilename = "AliCascadeAnalysis_Ishaan_Draw.root", TString outputFolder = ".", Bool_t saveImages = kFALSE, Bool_t saveStack = kFALSE)
{
    // TDirectory::AddDirectory(0);
    gStyle->SetOptFit(1111);
    // gSystem->Chmod(outputFolder, 0777);
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
    TH1 *resultParams_xip[nptbins_Xi][nmultbins_Xi];
    TH1 *resultParams_xim[nptbins_Xi][nmultbins_Xi];
    TH1 *resultParams_omp[nptbins_Om][nmultbins_Om];
    TH1 *resultParams_omm[nptbins_Om][nmultbins_Om];
    TH1D *rawPt_xim[nmultbins_Xi];
    TH1D *rawPt_xip[nmultbins_Xi];
    TH1D *rawPt_omm[nmultbins_Om];
    TH1D *rawPt_omp[nmultbins_Om];

    TH1 *h_MassXiC_pt_mult[nptbins_Xi][nmultbins_Xi];
    TH1 *h_MassOmC_pt_mult[nptbins_Om][nmultbins_Om];
    TH1 *resultParams_xiC[nptbins_Xi][nmultbins_Xi];
    TH1 *resultParams_omC[nptbins_Om][nmultbins_Om];
    TH1D *rawPt_xiC[nmultbins_Xi];
    TH1D *rawPt_omC[nmultbins_Om];

    TH1 *InvMass_Xim;
    TH1 *InvMass_Xip;
    TH1 *InvMass_Omm;
    TH1 *InvMass_Omp;
    TH1 *h_multBinEntries_Xi;
    TH1 *h_multBinEntries_Om;

    auto hs_xip = new THStack("xip", "Raw #it{p}_{T} spectra #Xi^{+}");
    auto hs_xim = new THStack("xim", "Raw #it{p}_{T} spectra #Xi^{-}");
    auto hs_omp = new THStack("omp", "Raw #it{p}_{T} spectra #Omega^{+}");
    auto hs_omm = new THStack("omm", "Raw #it{p}_{T} spectra #Omega^{-}");

    auto hs_xiC = new THStack("xiC", "Raw #it{p}_{T} spectra #Xi^{+} + #Xi^{-}");
    auto hs_omC = new THStack("omC", "Raw #it{p}_{T} spectra #Omega^{+} + #Omega^{-}");
    // auto hs_redChi2_xip = new THStack("chi2_xip", "#chi^{2} #Xi^{+}");
    // auto hs_redChi2_xim = new THStack("chi2_xim", "#chi^{2} #Xi^{-}");

    // TH1D *dummy_xip = new TH1D("dummy_xip", "Raw #it{p}_{T} spectra #Xi^{+}", nptbins_Xi, ptbins_Xi);
    // TH1D *dummy_xim = new TH1D("dummy_xim", "Raw #it{p}_{T} spectra #Xi^{-}", nptbins_Xi, ptbins_Xi);
    // TH1D *dummy_omp = new TH1D("dummy_omp", "Raw #it{p}_{T} spectra #Omega^{+}", nptbins_Om, ptbins_Om);
    // TH1D *dummy_omm = new TH1D("dummy_omm", "Raw #it{p}_{T} spectra #Omega^{-}", nptbins_Om, ptbins_Om);

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

        InvMass_Xim = (TH1 *)f->Get("InvMass_Xim");
        InvMass_Xip = (TH1 *)f->Get("InvMass_Xip");
        InvMass_Omm = (TH1 *)f->Get("InvMass_Omm");
        InvMass_Omp = (TH1 *)f->Get("InvMass_Omp");
        h_multBinEntries_Xi = (TH1 *)f->Get("h_multBinEntries_Xi");
        h_multBinEntries_Om = (TH1 *)f->Get("h_multBinEntries_Om");
        for (Int_t ptBinXi = 0; ptBinXi < nptbins_Xi; ptBinXi++)
        {
            for (Int_t multBinXi = 0; multBinXi < nmultbins_Xi; multBinXi++)
            {
                h_MassXim_pt_mult[ptBinXi][multBinXi] = (TH1 *)f->FindObjectAny(TString::Format(("h_MassXim_pt_mult[%d][%d]"), ptBinXi, multBinXi));
                h_MassXip_pt_mult[ptBinXi][multBinXi] = (TH1 *)f->FindObjectAny(TString::Format(("h_MassXip_pt_mult[%d][%d]"), ptBinXi, multBinXi));
                resultParams_xip[ptBinXi][multBinXi] = (TH1 *)f->FindObjectAny(TString::Format(("resultParams_xip[%d][%d]"), ptBinXi, multBinXi));
                resultParams_xim[ptBinXi][multBinXi] = (TH1 *)f->FindObjectAny(TString::Format(("resultParams_xim[%d][%d]"), ptBinXi, multBinXi));

                h_MassXiC_pt_mult[ptBinXi][multBinXi] = (TH1 *)f->FindObjectAny(TString::Format(("h_MassXiC_pt_mult[%d][%d]"), ptBinXi, multBinXi));
                resultParams_xiC[ptBinXi][multBinXi] = (TH1 *)f->FindObjectAny(TString::Format(("resultParams_xiC[%d][%d]"), ptBinXi, multBinXi));
            }
        }
        for (Int_t ptBinOm = 0; ptBinOm < nptbins_Om; ptBinOm++)
        {
            for (Int_t multBinOm = 0; multBinOm < nmultbins_Om; multBinOm++)
            {
                h_MassOmm_pt_mult[ptBinOm][multBinOm] = (TH1 *)f->FindObjectAny(TString::Format(("h_MassOmm_pt_mult[%d][%d]"), ptBinOm, multBinOm));
                h_MassOmp_pt_mult[ptBinOm][multBinOm] = (TH1 *)f->FindObjectAny(TString::Format(("h_MassOmp_pt_mult[%d][%d]"), ptBinOm, multBinOm));
                resultParams_omp[ptBinOm][multBinOm] = (TH1 *)f->FindObjectAny(TString::Format(("resultParams_omp[%d][%d]"), ptBinOm, multBinOm));
                resultParams_omm[ptBinOm][multBinOm] = (TH1 *)f->FindObjectAny(TString::Format(("resultParams_omm[%d][%d]"), ptBinOm, multBinOm));

                h_MassOmC_pt_mult[ptBinOm][multBinOm] = (TH1 *)f->FindObjectAny(TString::Format(("h_MassOmC_pt_mult[%d][%d]"), ptBinOm, multBinOm));
                resultParams_omC[ptBinOm][multBinOm] = (TH1 *)f->FindObjectAny(TString::Format(("resultParams_omC[%d][%d]"), ptBinOm, multBinOm));
            }
        }

        f->Close();
    } /// Input ended!

    /// Saving output :
    TString output = outputFilename;
    if (outputFilename.IsNull())
    {
        output = input;
        output.ReplaceAll("AliCascadeAnalysis_Ishaan_Fitting.root", "AliCascadeAnalysis_Ishaan_Draw.root");
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
    /// Output set.

    /// Begin
    Int_t markerStyles[] = {20, 21, 22, 23, 29, 33, 34, 43, 47, 41}; // chosen marker style palette

    for (Int_t multBinXi = 0; multBinXi < nmultbins_Xi; multBinXi++)
    {
        rawPt_xim[multBinXi] = new TH1D(TString::Format(("rawPt_xim[%d]"), multBinXi), TString::Format(("Mult: %.0f-%.0f%%, #times2^{%d}"), multbins_Xi[multBinXi], multbins_Xi[multBinXi + 1], (nmultbins_Xi - 1) - multBinXi), nptbins_Xi, ptbins_Xi);
        rawPt_xip[multBinXi] = new TH1D(TString::Format(("rawPt_xip[%d]"), multBinXi), TString::Format(("Mult: %.0f-%.0f%%, #times2^{%d}"), multbins_Xi[multBinXi], multbins_Xi[multBinXi + 1], (nmultbins_Xi - 1) - multBinXi), nptbins_Xi, ptbins_Xi);

        rawPt_xiC[multBinXi] = new TH1D(TString::Format(("rawPt_xiC[%d]"), multBinXi), TString::Format(("Mult: %.0f-%.0f%%, #times2^{%d}"), multbins_Xi[multBinXi], multbins_Xi[multBinXi + 1], (nmultbins_Xi - 1) - multBinXi), nptbins_Xi, ptbins_Xi);
        // redChi2_xim[multBinXi] = new TH1D(TString::Format(("redChi2_xim[%d]"), multBinXi), TString::Format(("Mult: %.0f-%.0f%%"), multbins_Xi[multBinXi ], multbins_Xi[multBinXi+1]), nptbins_Xi, ptbins_Xi);
        // redChi2_xip[multBinXi] = new TH1D(TString::Format(("redChi2_xip[%d]"), multBinXi), TString::Format(("Mult: %.0f-%.0f%%"), multbins_Xi[multBinXi ], multbins_Xi[multBinXi+1]), nptbins_Xi, ptbins_Xi);
        for (Int_t ptBinXi = 0; ptBinXi < nptbins_Xi; ptBinXi++)
        {
            rawPt_xim[multBinXi]->SetBinContent(ptBinXi + 1, resultParams_xim[ptBinXi][multBinXi]->GetBinContent(1) / (rawPt_xim[multBinXi]->GetBinWidth(ptBinXi + 1) * (h_multBinEntries_Xi->GetBinContent(multBinXi + 1)))); // Bin 1 in resultparams is raw pt's bin counting
            rawPt_xim[multBinXi]->SetBinError(ptBinXi + 1, resultParams_xim[ptBinXi][multBinXi]->GetBinError(1) / (rawPt_xim[multBinXi]->GetBinWidth(ptBinXi + 1) * (h_multBinEntries_Xi->GetBinContent(multBinXi + 1))));

            rawPt_xip[multBinXi]->SetBinContent(ptBinXi + 1, resultParams_xip[ptBinXi][multBinXi]->GetBinContent(1) / (rawPt_xip[multBinXi]->GetBinWidth(ptBinXi + 1) * (h_multBinEntries_Xi->GetBinContent(multBinXi + 1)))); // Bin 1 in resultparams is raw pt's bin counting
            rawPt_xip[multBinXi]->SetBinError(ptBinXi + 1, resultParams_xip[ptBinXi][multBinXi]->GetBinError(1) / (rawPt_xip[multBinXi]->GetBinWidth(ptBinXi + 1) * (h_multBinEntries_Xi->GetBinContent(multBinXi + 1))));

            rawPt_xiC[multBinXi]->SetBinContent(ptBinXi + 1, resultParams_xiC[ptBinXi][multBinXi]->GetBinContent(1) / (rawPt_xiC[multBinXi]->GetBinWidth(ptBinXi + 1) * (h_multBinEntries_Xi->GetBinContent(multBinXi + 1)))); // Bin 1 in resultparams is raw pt's bin counting
            rawPt_xiC[multBinXi]->SetBinError(ptBinXi + 1, resultParams_xiC[ptBinXi][multBinXi]->GetBinError(1) / (rawPt_xiC[multBinXi]->GetBinWidth(ptBinXi + 1) * (h_multBinEntries_Xi->GetBinContent(multBinXi + 1))));

            // redChi2_xim[multBinXi]->SetBinContent(ptBinXi, resultParams_xim[ptBinXi][multBinXi]->GetBinContent(2) / (redChi2_xim[multBinXi]->GetBinWidth(ptBinXi) * multbinEntries_Xim[multBinXi])); //Bin 1 in resultparams is raw pt's bin counting
            // redChi2_xim[multBinXi]->SetBinError(ptBinXi, resultParams_xim[ptBinXi][multBinXi]->GetBinError(2) / (redChi2_xim[multBinXi]->GetBinWidth(ptBinXi) * multbinEntries_Xim[multBinXi]));
            // redChi2_xip[multBinXi]->SetBinContent(ptBinXi, resultParams_xip[ptBinXi][multBinXi]->GetBinContent(2) / (redChi2_xip[multBinXi]->GetBinWidth(ptBinXi) * multbinEntries_Xip[multBinXi])); //Bin 1 in resultparams is raw pt's bin counting
            // redChi2_xip[multBinXi]->SetBinError(ptBinXi, resultParams_xip[ptBinXi][multBinXi]->GetBinError(2) / (redChi2_xip[multBinXi]->GetBinWidth(ptBinXi) * multbinEntries_Xip[multBinXi]));

            /// Scale for bin width: N->dN/dpt
            // rawPt_xim[multBinXi]->Scale(1, "width");
            // rawPt_xip[multBinXi]->Scale(1, "width");

            h_MassXim_pt_mult[ptBinXi][multBinXi]->GetXaxis()->SetTitle("#Lambda^{0}-#pi^{#pm} Inv. Mass (GeV/c^{2})");
            h_MassXip_pt_mult[ptBinXi][multBinXi]->GetXaxis()->SetTitle("#Lambda^{0}-#pi^{#pm} Inv. Mass (GeV/c^{2})");
            h_MassXiC_pt_mult[ptBinXi][multBinXi]->GetXaxis()->SetTitle("#Lambda^{0}-#pi^{#pm} Inv. Mass (GeV/c^{2})");

            TF1 *f1 = (TF1 *)h_MassXim_pt_mult[ptBinXi][multBinXi]->GetListOfFunctions()->At(0); /// 0 position = peak Function (GausPol2)
            DrawAndSave(h_MassXim_pt_mult[ptBinXi][multBinXi], f1->GetParameter(1), f1->GetParameter(2), saveImages, outputFolder);
            TF1 *f2 = (TF1 *)h_MassXip_pt_mult[ptBinXi][multBinXi]->GetListOfFunctions()->At(0); /// peak Function (GausPol2)
            DrawAndSave(h_MassXip_pt_mult[ptBinXi][multBinXi], f2->GetParameter(1), f2->GetParameter(2), saveImages, outputFolder);
            TF1 *f3 = (TF1 *)h_MassXiC_pt_mult[ptBinXi][multBinXi]->GetListOfFunctions()->At(0); /// peak Function (GausPol2)
            DrawAndSave(h_MassXiC_pt_mult[ptBinXi][multBinXi], f3->GetParameter(1), f3->GetParameter(2), saveImages, outputFolder);
        }
        out->cd("dirRawPt_xip");
        rawPt_xip[multBinXi]->SetMarkerStyle(markerStyles[multBinXi]);
        rawPt_xip[multBinXi]->Write();
        rawPt_xip[multBinXi]->Scale(pow(2, (nmultbins_Xi - 1) - multBinXi));
        hs_xip->Add(rawPt_xip[multBinXi]);

        // //
        // redChi2_xip[multBinXi]->SetMarkerStyle(markerStyles[multBinXi]);
        // hs_redChi2_xip->Add(redChi2_xip[multBinXi]);
        ///
        out->cd("dirRawPt_xim");
        rawPt_xim[multBinXi]->SetMarkerStyle(markerStyles[multBinXi]);
        rawPt_xim[multBinXi]->Write();
        rawPt_xim[multBinXi]->Scale(pow(2, (nmultbins_Xi - 1) - multBinXi));
        hs_xim->Add(rawPt_xim[multBinXi]);

        out->cd("dirRawPt_xiC");
        rawPt_xiC[multBinXi]->SetMarkerStyle(markerStyles[multBinXi]);
        rawPt_xiC[multBinXi]->Write();
        rawPt_xiC[multBinXi]->Scale(pow(2, (nmultbins_Xi - 1) - multBinXi));
        hs_xiC->Add(rawPt_xiC[multBinXi]);

        //
        // redChi2_xim[multBinXi]->SetMarkerStyle(markerStyles[multBinXi]);
        // hs_redChi2_xim->Add(redChi2_xim[multBinXi]);
        ///
    }
    for (Int_t multBinOm = 0; multBinOm < nmultbins_Om; multBinOm++)
    {
        rawPt_omm[multBinOm] = new TH1D(TString::Format(("rawPt_omm[%d]"), multBinOm), TString::Format(("Mult: %.0f-%.0f%%, #times2^{%d}"), multbins_Om[multBinOm], multbins_Om[multBinOm + 1], (nmultbins_Om - 1) - multBinOm), nptbins_Om, ptbins_Om);
        rawPt_omp[multBinOm] = new TH1D(TString::Format(("rawPt_omp[%d]"), multBinOm), TString::Format(("Mult: %.0f-%.0f%%, #times2^{%d}"), multbins_Om[multBinOm], multbins_Om[multBinOm + 1], (nmultbins_Om - 1) - multBinOm), nptbins_Om, ptbins_Om);

        rawPt_omC[multBinOm] = new TH1D(TString::Format(("rawPt_omC[%d]"), multBinOm), TString::Format(("Mult: %.0f-%.0f%%, #times2^{%d}"), multbins_Om[multBinOm], multbins_Om[multBinOm + 1], (nmultbins_Om - 1) - multBinOm), nptbins_Om, ptbins_Om);

        for (Int_t ptBinOm = 0; ptBinOm < nptbins_Om; ptBinOm++)
        {
            rawPt_omm[multBinOm]->SetBinContent(ptBinOm + 1, resultParams_omm[ptBinOm][multBinOm]->GetBinContent(1) / (rawPt_omm[multBinOm]->GetBinWidth(ptBinOm + 1) * (h_multBinEntries_Om->GetBinContent(multBinOm + 1)))); // Bin 1 in resultparams is raw pt's bin counting
            rawPt_omm[multBinOm]->SetBinError(ptBinOm + 1, resultParams_omm[ptBinOm][multBinOm]->GetBinError(1) / (rawPt_omm[multBinOm]->GetBinWidth(ptBinOm + 1) * (h_multBinEntries_Om->GetBinContent(multBinOm + 1))));
            rawPt_omp[multBinOm]->SetBinContent(ptBinOm + 1, resultParams_omp[ptBinOm][multBinOm]->GetBinContent(1) / (rawPt_omp[multBinOm]->GetBinWidth(ptBinOm + 1) * (h_multBinEntries_Om->GetBinContent(multBinOm + 1)))); // Bin 1 in resultparams is raw pt's bin counting
            rawPt_omp[multBinOm]->SetBinError(ptBinOm + 1, resultParams_omp[ptBinOm][multBinOm]->GetBinError(1) / (rawPt_omp[multBinOm]->GetBinWidth(ptBinOm + 1) * (h_multBinEntries_Om->GetBinContent(multBinOm + 1))));

            rawPt_omC[multBinOm]->SetBinContent(ptBinOm + 1, resultParams_omC[ptBinOm][multBinOm]->GetBinContent(1) / (rawPt_omC[multBinOm]->GetBinWidth(ptBinOm + 1) * (h_multBinEntries_Om->GetBinContent(multBinOm + 1)))); // Bin 1 in resultparams is raw pt's bin counting
            rawPt_omC[multBinOm]->SetBinError(ptBinOm + 1, resultParams_omC[ptBinOm][multBinOm]->GetBinError(1) / (rawPt_omC[multBinOm]->GetBinWidth(ptBinOm + 1) * (h_multBinEntries_Om->GetBinContent(multBinOm + 1))));

            /// Scale for bin width: N->dN/dpt
            // rawPt_omm[multBinOm]->Scale(1, "width");
            // rawPt_omp[multBinOm]->Scale(1, "width");

            h_MassOmm_pt_mult[ptBinOm][multBinOm]->GetXaxis()->SetTitle("#Lambda^{0}-K^{#pm} Inv. Mass (GeV/c^{2})");
            h_MassOmp_pt_mult[ptBinOm][multBinOm]->GetXaxis()->SetTitle("#Lambda^{0}-K^{#pm} Inv. Mass (GeV/c^{2})");
            h_MassOmC_pt_mult[ptBinOm][multBinOm]->GetXaxis()->SetTitle("#Lambda^{0}-K^{#pm} Inv. Mass (GeV/c^{2})");

            TF1 *f1 = (TF1 *)h_MassOmm_pt_mult[ptBinOm][multBinOm]->GetListOfFunctions()->At(0); /// peak Function (GausPol2)
            DrawAndSave(h_MassOmm_pt_mult[ptBinOm][multBinOm], f1->GetParameter(1), f1->GetParameter(2), saveImages, outputFolder);
            TF1 *f2 = (TF1 *)h_MassOmp_pt_mult[ptBinOm][multBinOm]->GetListOfFunctions()->At(0); /// peak Function (GausPol2)
            DrawAndSave(h_MassOmp_pt_mult[ptBinOm][multBinOm], f2->GetParameter(1), f2->GetParameter(2), saveImages, outputFolder);
            TF1 *f3 = (TF1 *)h_MassOmC_pt_mult[ptBinOm][multBinOm]->GetListOfFunctions()->At(0); /// peak Function (GausPol2)
            DrawAndSave(h_MassOmC_pt_mult[ptBinOm][multBinOm], f3->GetParameter(1), f3->GetParameter(2), saveImages, outputFolder);
        }
        out->cd("dirRawPt_omp");
        rawPt_omp[multBinOm]->SetMarkerStyle(markerStyles[multBinOm]);
        rawPt_omp[multBinOm]->Write();
        rawPt_omp[multBinOm]->Scale(pow(2, (nmultbins_Om - 1) - multBinOm));
        hs_omp->Add(rawPt_omp[multBinOm]);

        out->cd("dirRawPt_omm");
        rawPt_omm[multBinOm]->SetMarkerStyle(markerStyles[multBinOm]);
        rawPt_omm[multBinOm]->Write();
        rawPt_omm[multBinOm]->Scale(pow(2, (nmultbins_Om - 1) - multBinOm));
        hs_omm->Add(rawPt_omm[multBinOm]);

        out->cd("dirRawPt_omC");
        rawPt_omC[multBinOm]->SetMarkerStyle(markerStyles[multBinOm]);
        rawPt_omC[multBinOm]->Write();
        rawPt_omC[multBinOm]->Scale(pow(2, (nmultbins_Om - 1) - multBinOm));
        hs_omC->Add(rawPt_omC[multBinOm]);
    }
    out->Close();

    gStyle->SetOptStat(0);
    gStyle->SetPalette(kVisibleSpectrum);

    TCanvas *c = new TCanvas("c", "c", 2560, 1440);
    c->Divide(2, 3);

    c->cd(1);
    hs_xip->Draw("plc pmc nostack");
    gPad->SetLogy();
    gPad->BuildLegend(0.8, 0.55, 0.99, 0.99, "");
    // if (saveStack)
    // {
    //     gPad->Print(Form("%s/images/rawPt/xipN.png", outputFolder.Data()), "png");
    // }

    c->cd(2);
    hs_omp->Draw("plc pmc nostack");
    gPad->SetLogy();
    gPad->BuildLegend(0.8, 0.75, 0.99, 0.99, "");
    // if (saveStack)
    // {
    //     gPad->Print(Form("%s/images/rawPt/ompN.png", outputFolder.Data()), "png");
    // }
    c->cd(3);
    hs_xim->Draw("plc pmc nostack");
    gPad->SetLogy();
    gPad->BuildLegend(0.8, 0.55, 0.99, 0.99, "");
    // if (saveStack)
    // {
    //     gPad->Print(Form("%s/images/rawPt/ximN.png", outputFolder.Data()), "png");
    // }
    c->cd(4);
    hs_omm->Draw("plc pmc nostack");
    gPad->SetLogy();
    gPad->BuildLegend(0.8, 0.75, 0.99, 0.99, "");
    // if (saveStack)
    // {
    //     gPad->Print(Form("%s/images/rawPt/ommN.png", outputFolder.Data()), "png");
    // }
    c->cd(5);
    hs_xiC->Draw("plc pmc nostack");
    gPad->SetLogy();
    gPad->BuildLegend(0.8, 0.55, 0.99, 0.99, "");
    // if (saveStack)
    // {
    //     gPad->Print(Form("%s/images/rawPt/xicN.png", outputFolder.Data()), "png");
    // }
    c->cd(6);
    hs_omC->Draw("plc pmc nostack");
    gPad->SetLogy();
    gPad->BuildLegend(0.8, 0.75, 0.99, 0.99, "");
    // if (saveStack)
    // {
    //     gPad->Print(Form("%s/images/rawPt/omcN.png", outputFolder.Data()), "png");
    // }
    // c->cd(5);
    // hs_redChi2_xip->Draw("plc pmc nostack");
    // gPad->SetLogy();
    // gPad->BuildLegend(0.8, 0.55, 0.99, 0.99, "");

    // c->cd(6);
    // hs_redChi2_xim->Draw("plc pmc nostack");
    // gPad->SetLogy();
    // gPad->BuildLegend(0.8, 0.55, 0.99, 0.99, "");

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

    c->Modified();
    c->ForceUpdate();

    if (saveStack)
    {
        c->cd(1);
        gPad->Print(Form("%s/images/rawPt/xipN.png", outputFolder.Data()), "png");

        c->cd(2);
        gPad->Print(Form("%s/images/rawPt/ximN.png", outputFolder.Data()), "png");

        c->cd(3);
        gPad->Print(Form("%s/images/rawPt/ompN.png", outputFolder.Data()), "png");

        c->cd(4);
        gPad->Print(Form("%s/images/rawPt/ommN.png", outputFolder.Data()), "png");

        c->cd(5);
        gPad->Print(Form("%s/images/rawPt/xicN.png", outputFolder.Data()), "png");

        c->cd(6);
        gPad->Print(Form("%s/images/rawPt/omcN.png", outputFolder.Data()), "png");
    }


    // delete out;
    return 0;
}
void DrawAndSave(TH1 *peak, Double_t pPosition, Double_t pWidth, Bool_t saveImages, TString outputFolder)
{
    if (saveImages)
    {
        gROOT->SetBatch(kTRUE);
        TCanvas *c1 = new TCanvas(peak->GetName(), peak->GetTitle(), 900, 600);

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
        TString imageFolder = peak->GetName();
        imageFolder = imageFolder(0, 17); /// getting substring for naming purpose
        gPad->Print(Form("%s/images/%s/%s.png", outputFolder.Data(), imageFolder.Data(), peak->GetName()), "png");
        // delete c1, lLineLeft, lLineRight;
    }
}
