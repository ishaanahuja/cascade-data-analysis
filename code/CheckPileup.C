/// DEPRECATED: Needs rewrite. Use at your own risk.

/**
 * @file CheckPileup.C
 * @brief Deprecated.
 * 
 * @note Was supposed to check for pileup - initial version. Don't use. For refernce only. Obsolete and probably incorrect.
 * @version 0.1
 * @date 2025-11-02
*/

#include <TString.h>

int CheckPileup(TString outputFilename = "AliCascadeAnalysis_Ishaan_Pileup.root")
{
    TDirectory::AddDirectory(0);
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

    Double_t Np[3][10] = {{257843, 256216, 256389, 255536, 469123, 468926, 467426, 467782, 891895, 796079}, {371959, 368884, 369641, 368505, 673370, 672042, 672240, 672469, 1280758, 1143463}, {362815, 359509, 358759, 356060, 650986, 654554, 651482, 652266, 1240716, 1099142}};

    TString inputRunFiles[] = {"/var/home/ishaan/Work/CERN/ishaan-ahuja/analysis/Analysis_Results/pPb6runsAfterEventSel/pPb6runsAfterEventSel_Draw.root", "/var/home/ishaan/Work/CERN/ishaan-ahuja/analysis/Analysis_Results/pPbGrp1PileUp2/pPbGrp1PileUp2_Draw.root", "/var/home/ishaan/Work/CERN/ishaan-ahuja/analysis/Analysis_Results/pPbGrp2PileUp3/pPbGrp2PileUp3_Draw.root", "/var/home/ishaan/Work/CERN/ishaan-ahuja/analysis/Analysis_Results/pPbGrp3PileUp5/pPbGrp3PileUp5_Draw.root"};
    Int_t totalRunFiles = 4;

    TH1 *rawPt_xim[totalRunFiles][nmultbins_Xi];
    TH1 *rawPt_xip[totalRunFiles][nmultbins_Xi];
    TH1 *rawPt_omm[totalRunFiles][nmultbins_Om];
    TH1 *rawPt_omp[totalRunFiles][nmultbins_Om];
    TH1 *rawPt_xiC[totalRunFiles][nmultbins_Xi];
    TH1 *rawPt_omC[totalRunFiles][nmultbins_Om];

    TH1D *rawPtRatio_xim[totalRunFiles][nmultbins_Xi];
    TH1D *rawPtRatio_xip[totalRunFiles][nmultbins_Xi];
    TH1D *rawPtRatio_omm[totalRunFiles][nmultbins_Om];
    TH1D *rawPtRatio_omp[totalRunFiles][nmultbins_Om];
    TH1D *rawPtRatio_xiC[totalRunFiles][nmultbins_Xi];
    TH1D *rawPtRatio_omC[totalRunFiles][nmultbins_Om];

    ///
    TH1D *rawPtNp_xim[totalRunFiles][nmultbins_Xi];
    TH1D *rawPtNp_xip[totalRunFiles][nmultbins_Xi];
    TH1D *rawPtNp_xiC[totalRunFiles][nmultbins_Xi];
    ///

    auto hs_xip = new THStack("xip", "Ratio of raw #it{p}_{T} spectra #Xi^{+}");
    auto hs_xim = new THStack("xim", "Ratio of raw #it{p}_{T} spectra #Xi^{-}");
    auto hs_omp = new THStack("omp", "Ratio of raw #it{p}_{T} spectra #Omega^{+}");
    auto hs_omm = new THStack("omm", "Ratio of raw #it{p}_{T} spectra #Omega^{-}");

    auto hs_xiC = new THStack("xiC", "Ratio of raw #it{p}_{T} spectra #Xi^{+} + #Xi^{-}");
    auto hs_omC = new THStack("omC", "Ratio of raw #it{p}_{T} spectra #Omega^{+} + #Omega^{-}");

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

        for (Int_t multBinXi = 0; multBinXi < nmultbins_Xi; multBinXi++)
        {
            rawPt_xim[inputIndex][multBinXi] = (TH1 *)f->FindObjectAny(TString::Format(("rawPt_xim[%d]"), multBinXi));
            rawPt_xip[inputIndex][multBinXi] = (TH1 *)f->FindObjectAny(TString::Format(("rawPt_xip[%d]"), multBinXi));

            rawPt_xiC[inputIndex][multBinXi] = (TH1 *)f->FindObjectAny(TString::Format(("rawPt_xiC[%d]"), multBinXi));
        }

        for (Int_t multBinOm = 0; multBinOm < nmultbins_Om; multBinOm++)
        {
            rawPt_omm[inputIndex][multBinOm] = (TH1 *)f->FindObjectAny(TString::Format(("rawPt_omm[%d]"), multBinOm));
            rawPt_omp[inputIndex][multBinOm] = (TH1 *)f->FindObjectAny(TString::Format(("rawPt_omp[%d]"), multBinOm));
            rawPt_omC[inputIndex][multBinOm] = (TH1 *)f->FindObjectAny(TString::Format(("rawPt_omC[%d]"), multBinOm));
        }
        f->Close();
    }

    /// Saving output :
    TString output = outputFilename;
    Printf("Saving output to '%s' ...", output.Data());
    TFile *out = TFile::Open(output.Data(), "RECREATE");
    if (!out)
    {
        Printf("Error: Cannot open file '%s' !", output.Data());
        return 3;
    }

    ////
    TFile *outNp = TFile::Open("Pileup_Np.root", "RECREATE");
    if (!outNp)
    {
        Printf("Error: Cannot open file Pileup_Np.root!");
        return 4;
    }
    ////

    out->mkdir("dirRawPtRatio_xim");
    out->mkdir("dirRawPtRatio_xip");
    out->mkdir("dirRawPtRatio_omm");
    out->mkdir("dirRawPtRatio_omp");

    out->mkdir("dirRawPtRatio_xiC");
    out->mkdir("dirRawPtRatio_omC");
    /// Output set.

    /// Begin
    // Int_t markerStyles[] = {20, 21, 22, 23, 29, 33, 34, 43, 47, 41}; // chosen marker style palette

    for (int inputIndex = 1; inputIndex < totalRunFiles; inputIndex++) /// go over pile-ups only
    {
        for (Int_t multBinXi = 0; multBinXi < nmultbins_Xi; multBinXi++)
        {
            ////
            rawPtNp_xim[inputIndex][multBinXi] = (TH1D *)rawPt_xim[0][multBinXi]->Clone(TString::Format(("rawPtNp_xim[%d][%d]"), inputIndex, multBinXi));
            rawPtNp_xip[inputIndex][multBinXi] = (TH1D *)rawPt_xip[0][multBinXi]->Clone(TString::Format(("rawPtNp_xip[%d][%d]"), inputIndex, multBinXi));
            rawPtNp_xiC[inputIndex][multBinXi] = (TH1D *)rawPt_xiC[0][multBinXi]->Clone(TString::Format(("rawPtNp_xiC[%d][%d]"), inputIndex, multBinXi));

            rawPtNp_xim[inputIndex][multBinXi]->SetTitle(TString::Format(("#Xi^{-}: Pileup Grp: %d, Mult: %.0f-%.0f%%"), inputIndex, multbins_Xi[multBinXi], multbins_Xi[multBinXi + 1]));
            rawPtNp_xip[inputIndex][multBinXi]->SetTitle(TString::Format(("#Xi^{+}: Pileup Grp: %d, Mult: %.0f-%.0f%%"), inputIndex, multbins_Xi[multBinXi], multbins_Xi[multBinXi + 1]));
            rawPtNp_xiC[inputIndex][multBinXi]->SetTitle(TString::Format(("#Xi^{+} + #Xi^{-}: Pileup Grp: %d, Mult: %.0f-%.0f%%"), inputIndex, multbins_Xi[multBinXi], multbins_Xi[multBinXi + 1]));

            rawPtNp_xim[inputIndex][multBinXi]->GetXaxis()->SetTitle("#it{p}_{T} (GeV/c)");
            rawPtNp_xip[inputIndex][multBinXi]->GetXaxis()->SetTitle("#it{p}_{T} (GeV/c)");
            rawPtNp_xiC[inputIndex][multBinXi]->GetXaxis()->SetTitle("#it{p}_{T} (GeV/c)");

            rawPtNp_xim[inputIndex][multBinXi]->GetYaxis()->SetTitle("N_{p} #times #frac{#it{p}_{T}}{N_{g}}");
            rawPtNp_xip[inputIndex][multBinXi]->GetYaxis()->SetTitle("N_{p} #times #frac{#it{p}_{T}}{N_{g}}");
            rawPtNp_xiC[inputIndex][multBinXi]->GetYaxis()->SetTitle("N_{p} #times #frac{#it{p}_{T}}{N_{g}}");

            for (Int_t ptBinXi = 0; ptBinXi < nptbins_Xi; ptBinXi++)
            {
                rawPtNp_xim[inputIndex][multBinXi]->SetBinContent(ptBinXi + 1, (rawPt_xim[0][multBinXi]->GetBinContent(ptBinXi + 1) * Np[inputIndex - 1][multBinXi])); // for each pT bin: (pt/Ng)*Np
                rawPtNp_xip[inputIndex][multBinXi]->SetBinContent(ptBinXi + 1, (rawPt_xip[0][multBinXi]->GetBinContent(ptBinXi + 1) * Np[inputIndex - 1][multBinXi])); // for each pT bin: (pt/Ng)*Np
                rawPtNp_xiC[inputIndex][multBinXi]->SetBinContent(ptBinXi + 1, (rawPt_xiC[0][multBinXi]->GetBinContent(ptBinXi + 1) * Np[inputIndex - 1][multBinXi])); // for each pT bin: (pt/Ng)*Np
            }
            ////

            rawPtRatio_xim[inputIndex][multBinXi] = (TH1D *)rawPt_xim[inputIndex][multBinXi]->Clone(TString::Format(("rawPtRatio_xim[%d][%d]"), inputIndex, multBinXi));
            rawPtRatio_xip[inputIndex][multBinXi] = (TH1D *)rawPt_xip[inputIndex][multBinXi]->Clone(TString::Format(("rawPtRatio_xip[%d][%d]"), inputIndex, multBinXi));
            rawPtRatio_xiC[inputIndex][multBinXi] = (TH1D *)rawPt_xiC[inputIndex][multBinXi]->Clone(TString::Format(("rawPtRatio_xiC[%d][%d]"), inputIndex, multBinXi));

            rawPtRatio_xim[inputIndex][multBinXi]->SetTitle(TString::Format(("Pileup Grp: %d, Mult: %.0f-%.0f%%"), inputIndex, multbins_Xi[multBinXi], multbins_Xi[multBinXi + 1]));
            rawPtRatio_xip[inputIndex][multBinXi]->SetTitle(TString::Format(("Pileup Grp: %d, Mult: %.0f-%.0f%%"), inputIndex, multbins_Xi[multBinXi], multbins_Xi[multBinXi + 1]));
            rawPtRatio_xiC[inputIndex][multBinXi]->SetTitle(TString::Format(("Pileup Grp: %d, Mult: %.0f-%.0f%%"), inputIndex, multbins_Xi[multBinXi], multbins_Xi[multBinXi + 1]));

            Bool_t divide_xim = rawPtRatio_xim[inputIndex][multBinXi]->Divide(rawPt_xim[0][multBinXi], rawPt_xim[inputIndex][multBinXi]);
            Bool_t divide_xip = rawPtRatio_xip[inputIndex][multBinXi]->Divide(rawPt_xip[0][multBinXi], rawPt_xip[inputIndex][multBinXi]);
            Bool_t divide_xiC = rawPtRatio_xiC[inputIndex][multBinXi]->Divide(rawPt_xiC[0][multBinXi], rawPt_xiC[inputIndex][multBinXi]);

            if (!(divide_xim && divide_xip && divide_xiC))
            {
                Printf("Problems with divide in input index = %d and multbinXi = %d", inputIndex, multBinXi);
                return 9;
            }

            rawPtRatio_xim[inputIndex][multBinXi]->GetXaxis()->SetTitle("#it{p}_{T} (GeV/c)");
            rawPtRatio_xip[inputIndex][multBinXi]->GetXaxis()->SetTitle("#it{p}_{T} (GeV/c)");
            rawPtRatio_xiC[inputIndex][multBinXi]->GetXaxis()->SetTitle("#it{p}_{T} (GeV/c)");

            rawPtRatio_xim[inputIndex][multBinXi]->GetYaxis()->SetTitle("Ratio");
            rawPtRatio_xip[inputIndex][multBinXi]->GetYaxis()->SetTitle("Ratio");
            rawPtRatio_xiC[inputIndex][multBinXi]->GetYaxis()->SetTitle("Ratio");

            out->cd("dirRawPtRatio_xip");
            rawPtRatio_xip[inputIndex][multBinXi]->Write();
            // rawPtRatio_xip[inputIndex][multBinXi]->Scale(pow(2, (nmultbins_Xi - 1) - multBinXi));
            hs_xip->Add(rawPtRatio_xip[inputIndex][multBinXi]);

            out->cd("dirRawPtRatio_xim");
            rawPtRatio_xim[inputIndex][multBinXi]->Write();
            // rawPtRatio_xim[inputIndex][multBinXi]->Scale(pow(2, (nmultbins_Xi - 1) - multBinXi));
            hs_xim->Add(rawPtRatio_xim[inputIndex][multBinXi]);

            out->cd("dirRawPtRatio_xiC");
            rawPtRatio_xiC[inputIndex][multBinXi]->Write();
            // rawPtRatio_xiC[inputIndex][multBinXi]->Scale(pow(2, (nmultbins_Xi - 1) - multBinXi));
            hs_xiC->Add(rawPtRatio_xiC[inputIndex][multBinXi]);

            ////
            outNp->cd();
            rawPtNp_xim[inputIndex][multBinXi]->Write();
            rawPtNp_xip[inputIndex][multBinXi]->Write();
            rawPtNp_xiC[inputIndex][multBinXi]->Write();
            ////
        }

        for (Int_t multBinOm = 0; multBinOm < nmultbins_Om; multBinOm++)
        {

            rawPtRatio_omm[inputIndex][multBinOm] = (TH1D *)rawPt_omm[inputIndex][multBinOm]->Clone(TString::Format(("rawPtRatio_omm[%d][%d]"), inputIndex, multBinOm));
            rawPtRatio_omp[inputIndex][multBinOm] = (TH1D *)rawPt_omp[inputIndex][multBinOm]->Clone(TString::Format(("rawPtRatio_omp[%d][%d]"), inputIndex, multBinOm));
            rawPtRatio_omC[inputIndex][multBinOm] = (TH1D *)rawPt_omC[inputIndex][multBinOm]->Clone(TString::Format(("rawPtRatio_omC[%d][%d]"), inputIndex, multBinOm));

            rawPtRatio_omm[inputIndex][multBinOm]->SetTitle(TString::Format(("Pileup Grp: %d, Mult: %.0f-%.0f%%"), inputIndex, multbins_Om[multBinOm], multbins_Om[multBinOm + 1]));
            rawPtRatio_omp[inputIndex][multBinOm]->SetTitle(TString::Format(("Pileup Grp: %d, Mult: %.0f-%.0f%%"), inputIndex, multbins_Om[multBinOm], multbins_Om[multBinOm + 1]));
            rawPtRatio_omC[inputIndex][multBinOm]->SetTitle(TString::Format(("Pileup Grp: %d, Mult: %.0f-%.0f%%"), inputIndex, multbins_Om[multBinOm], multbins_Om[multBinOm + 1]));

            Bool_t divide_omm = rawPtRatio_omm[inputIndex][multBinOm]->Divide(rawPt_omm[0][multBinOm], rawPt_omm[inputIndex][multBinOm]);
            Bool_t divide_omp = rawPtRatio_omp[inputIndex][multBinOm]->Divide(rawPt_omp[0][multBinOm], rawPt_omp[inputIndex][multBinOm]);
            Bool_t divide_omC = rawPtRatio_omC[inputIndex][multBinOm]->Divide(rawPt_omC[0][multBinOm], rawPt_omC[inputIndex][multBinOm]);

            if (divide_omm && divide_omp && divide_omC == kFALSE)
            {
                Printf("Problems with divide in input index = %d and multBinOm = %d", inputIndex, multBinOm);
                return 99;
            }

            rawPtRatio_omm[inputIndex][multBinOm]->GetXaxis()->SetTitle("#it{p}_{T} (GeV/c)");
            rawPtRatio_omp[inputIndex][multBinOm]->GetXaxis()->SetTitle("#it{p}_{T} (GeV/c)");
            rawPtRatio_omC[inputIndex][multBinOm]->GetXaxis()->SetTitle("#it{p}_{T} (GeV/c)");

            rawPtRatio_omm[inputIndex][multBinOm]->GetYaxis()->SetTitle("Ratio");
            rawPtRatio_omp[inputIndex][multBinOm]->GetYaxis()->SetTitle("Ratio");
            rawPtRatio_omC[inputIndex][multBinOm]->GetYaxis()->SetTitle("Ratio");

            out->cd("dirRawPtRatio_omp");
            rawPtRatio_omp[inputIndex][multBinOm]->Write();
            // rawPtRatio_omp[inputIndex][multBinOm]->Scale(pow(2, (nmultbins_Om - 1) - multBinOm));
            hs_omp->Add(rawPtRatio_omp[inputIndex][multBinOm]);

            out->cd("dirRawPtRatio_omm");
            rawPtRatio_omm[inputIndex][multBinOm]->Write();
            // rawPtRatio_omm[inputIndex][multBinOm]->Scale(pow(2, (nmultbins_Om - 1) - multBinOm));
            hs_omm->Add(rawPtRatio_omm[inputIndex][multBinOm]);

            out->cd("dirRawPtRatio_omC");
            rawPtRatio_omC[inputIndex][multBinOm]->Write();
            // rawPtRatio_omC[inputIndex][multBinOm]->Scale(pow(2, (nmultbins_Om - 1) - multBinOm));
            hs_omC->Add(rawPtRatio_omC[inputIndex][multBinOm]);
        }
    }
    out->Close();
    ////
    outNp->Close();
    ////
    gStyle->SetOptStat(0);
    gStyle->SetPalette(kVisibleSpectrum);

    TCanvas *c = new TCanvas("c", "c", 2560, 1440);
    c->Divide(2, 3);

    c->cd(1);
    hs_xip->Draw("plc pmc nostack");
    // gPad->SetLogy();
    gPad->BuildLegend(0.8, 0.55, 0.99, 0.99, "");

    c->cd(2);
    hs_xim->Draw("plc pmc nostack");
    // gPad->SetLogy();
    gPad->BuildLegend(0.8, 0.55, 0.99, 0.99, "");

    c->cd(3);
    hs_omp->Draw("plc pmc nostack");
    // gPad->SetLogy();
    gPad->BuildLegend(0.8, 0.75, 0.99, 0.99, "");

    c->cd(4);
    hs_omm->Draw("plc pmc nostack");
    // gPad->SetLogy();
    gPad->BuildLegend(0.8, 0.75, 0.99, 0.99, "");

    c->cd(5);
    hs_xiC->Draw("plc pmc nostack");
    // gPad->SetLogy();
    gPad->BuildLegend(0.8, 0.55, 0.99, 0.99, "");

    c->cd(6);
    hs_omC->Draw("plc pmc nostack");
    // gPad->SetLogy();
    gPad->BuildLegend(0.8, 0.75, 0.99, 0.99, "");

    hs_xip->GetXaxis()->SetTitle("#it{p}_{T} (GeV/c)");
    hs_xim->GetXaxis()->SetTitle("#it{p}_{T} (GeV/c)");
    hs_omp->GetXaxis()->SetTitle("#it{p}_{T} (GeV/c)");
    hs_omm->GetXaxis()->SetTitle("#it{p}_{T} (GeV/c)");
    hs_xiC->GetXaxis()->SetTitle("#it{p}_{T} (GeV/c)");
    hs_omC->GetXaxis()->SetTitle("#it{p}_{T} (GeV/c)");

    hs_xip->GetYaxis()->SetTitle("Ratio");
    hs_xim->GetYaxis()->SetTitle("Ratio");
    hs_omp->GetYaxis()->SetTitle("Ratio");
    hs_omm->GetYaxis()->SetTitle("Ratio");
    hs_xiC->GetYaxis()->SetTitle("Ratio");
    hs_omC->GetYaxis()->SetTitle("Ratio");

    c->Modified();
    c->ForceUpdate();

    // delete out;
    return 0;
}
