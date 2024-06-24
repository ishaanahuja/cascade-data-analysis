#include <TString.h>
#include <TSystem.h>
#include <TH1.h>
#include <TF1.h>
#include <TVirtualPad.h>
#include <TStyle.h>
#include <THStack.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TLine.h>
#include <TROOT.h>
#include <TLegend.h>

void DrawAndSave(TH1 *peak, TH1 *bg, TH1 *resultParams, Bool_t saveImages, TString outputFolder, TString imageFormat);

inline void PaintStack(TCanvas &c, THStack &hs, Bool_t setLogY = kTRUE, TString yAxisTitle = "#frac{1}{#it{N}_{inel}} #frac{d#it{N}}{d#it{p}_{T}}", TString xAxisTitle = "#it{p}_{T} (GeV/c)")
{
    c.cd();
    hs.Draw("plc pmc nostack");
    if (setLogY)
        gPad->SetLogy();
    gPad->BuildLegend(0.9, 0.6, 1., 1., "");
    hs.GetXaxis()->SetTitle(xAxisTitle.Data());
    hs.GetYaxis()->SetTitle(yAxisTitle.Data());
    c.Modified();
    c.ForceUpdate();
}

void SaveImage(TString imagePath, TString imageName, TString imageFormat)
{
    if (gSystem->AccessPathName(imagePath.Data())) /// returns true if folder path does NOT exist
    {
        gSystem->mkdir(imagePath.Data(), kTRUE); // makes the path if it doesn't exist
        gSystem->Chmod(imagePath.Data(), 0755);
    }
    gPad->Print(Form("%s/%s.%s", imagePath.Data(), imageName.Data(), imageFormat.Data()), imageFormat.Data());
    gSystem->Chmod(Form("%s/%s.%s", imagePath.Data(), imageName.Data(), imageFormat.Data()), 0755);
}

int CascadesDraw(std::string input = "", TString outputFilename = "", TString outputFolder = ".", TString ptRatioFilename = "", Bool_t fisMC = kFALSE, Bool_t saveImages = kFALSE, Bool_t saveStack = kTRUE, TString imageFormat = "png")
{
    gErrorIgnoreLevel = kWarning; /// suppresses printing of Info messages

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

    TH1 *h_MassXip_pt[nptbins_Xi];
    TH1 *h_MassXim_pt[nptbins_Xi];
    TH1 *h_MassOmp_pt[nptbins_Om];
    TH1 *h_MassOmm_pt[nptbins_Om];
    TH1 *h_MassXiC_pt[nptbins_Xi];
    TH1 *h_MassOmC_pt[nptbins_Om];

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

    TH1 *resultParams_Xip_allInt, *resultParams_Xim_allInt, *resultParams_XiC_allInt;
    TH1 *resultParams_Omp_allInt, *resultParams_Omm_allInt, *resultParams_OmC_allInt;
    /// background estimation hist through TSpectrum

    TH1 *h_bgXim_pt_mult[nptbins_Xi][nmultbins_Xi];
    TH1 *h_bgXip_pt_mult[nptbins_Xi][nmultbins_Xi];
    TH1 *h_bgXiC_pt_mult[nptbins_Xi][nmultbins_Xi];

    TH1 *h_bgOmm_pt_mult[nptbins_Om][nmultbins_Om];
    TH1 *h_bgOmp_pt_mult[nptbins_Om][nmultbins_Om];
    TH1 *h_bgOmC_pt_mult[nptbins_Om][nmultbins_Om];

    auto hs_xip = new THStack("hs_xip", "");
    auto hs_xim = new THStack("hs_xim", "");
    auto hs_omp = new THStack("hs_omp", "");
    auto hs_omm = new THStack("hs_omm", "");
    auto hs_xiC = new THStack("hs_xiC", "");
    auto hs_omC = new THStack("hs_omC", "");
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
    auto hs_xip_eff = new THStack("hs_xip_eff", "Efficiency #Xi^{+}");
    auto hs_xim_eff = new THStack("hs_xim_eff", "Efficiency #Xi^{-}");
    auto hs_omp_eff = new THStack("hs_omp_eff", "Efficiency #Omega^{+}");
    auto hs_omm_eff = new THStack("hs_omm_eff", "Efficiency #Omega^{-}");
    auto hs_xiC_eff = new THStack("hs_xiC_eff", "Efficiency #Xi^{+} + #Xi^{-}");
    auto hs_omC_eff = new THStack("hs_omC_eff", "Efficiency #Omega^{+} + #Omega^{-}");

    auto hs_xiC_eff_ratio = new THStack("hs_xiC_eff_ratio", "#Xi^{+} + #Xi^{-} efficiency ratio");
    auto hs_omC_eff_ratio = new THStack("hs_omC_eff_ratio", "#Omega^{+} + #Omega^{-} efficiency ratio");

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

        h_MassXim = (TH1 *)f->FindObjectAny("h_MassXim");
        h_MassXip = (TH1 *)f->FindObjectAny("h_MassXip");
        h_MassXiC = (TH1 *)f->FindObjectAny("h_MassXiC");
        h_MassOmm = (TH1 *)f->FindObjectAny("h_MassOmm");
        h_MassOmp = (TH1 *)f->FindObjectAny("h_MassOmp");
        h_MassOmC = (TH1 *)f->FindObjectAny("h_MassOmC");
        h_multBinEntries_Xi = (TH1 *)f->FindObjectAny("h_multBinEntries_Xi");
        h_multBinEntries_Om = (TH1 *)f->FindObjectAny("h_multBinEntries_Om");

        resultParams_Xip_allInt = (TH1 *)f->FindObjectAny("resultParams_Xip_allInt");
        resultParams_Xim_allInt = (TH1 *)f->FindObjectAny("resultParams_Xim_allInt");
        resultParams_XiC_allInt = (TH1 *)f->FindObjectAny("resultParams_XiC_allInt");
        resultParams_Omp_allInt = (TH1 *)f->FindObjectAny("resultParams_Omp_allInt");
        resultParams_Omm_allInt = (TH1 *)f->FindObjectAny("resultParams_Omm_allInt");
        resultParams_OmC_allInt = (TH1 *)f->FindObjectAny("resultParams_OmC_allInt");

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

                h_bgXim_pt_mult[ptBinXi][multBinXi] = (TH1 *)f->FindObjectAny(TString::Format(("h_bgXim_pt_mult[%d][%d]"), ptBinXi, multBinXi));
                h_bgXip_pt_mult[ptBinXi][multBinXi] = (TH1 *)f->FindObjectAny(TString::Format(("h_bgXip_pt_mult[%d][%d]"), ptBinXi, multBinXi));
                h_bgXiC_pt_mult[ptBinXi][multBinXi] = (TH1 *)f->FindObjectAny(TString::Format(("h_bgXiC_pt_mult[%d][%d]"), ptBinXi, multBinXi));
            }
            h_MassXip_pt[ptBinXi] = (TH1 *)f->FindObjectAny(TString::Format(("h_MassXip_pt[%d]"), ptBinXi));
            h_MassXim_pt[ptBinXi] = (TH1 *)f->FindObjectAny(TString::Format(("h_MassXim_pt[%d]"), ptBinXi));
            h_MassXiC_pt[ptBinXi] = (TH1 *)f->FindObjectAny(TString::Format(("h_MassXiC_pt[%d]"), ptBinXi));
            resultParXip_pt[ptBinXi] = (TH1 *)f->FindObjectAny(TString::Format(("resultParXip_pt[%d]"), ptBinXi));
            resultParXim_pt[ptBinXi] = (TH1 *)f->FindObjectAny(TString::Format(("resultParXim_pt[%d]"), ptBinXi));
            resultParXiC_pt[ptBinXi] = (TH1 *)f->FindObjectAny(TString::Format(("resultParXiC_pt[%d]"), ptBinXi));
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

                h_bgOmm_pt_mult[ptBinOm][multBinOm] = (TH1 *)f->FindObjectAny(TString::Format(("h_bgOmm_pt_mult[%d][%d]"), ptBinOm, multBinOm));
                h_bgOmp_pt_mult[ptBinOm][multBinOm] = (TH1 *)f->FindObjectAny(TString::Format(("h_bgOmp_pt_mult[%d][%d]"), ptBinOm, multBinOm));
                h_bgOmC_pt_mult[ptBinOm][multBinOm] = (TH1 *)f->FindObjectAny(TString::Format(("h_bgOmC_pt_mult[%d][%d]"), ptBinOm, multBinOm));
            }
            h_MassOmp_pt[ptBinOm] = (TH1 *)f->FindObjectAny(TString::Format(("h_MassOmp_pt[%d]"), ptBinOm));
            h_MassOmm_pt[ptBinOm] = (TH1 *)f->FindObjectAny(TString::Format(("h_MassOmm_pt[%d]"), ptBinOm));
            h_MassOmC_pt[ptBinOm] = (TH1 *)f->FindObjectAny(TString::Format(("h_MassOmC_pt[%d]"), ptBinOm));
            resultParOmp_pt[ptBinOm] = (TH1 *)f->FindObjectAny(TString::Format(("resultParOmp_pt[%d]"), ptBinOm));
            resultParOmm_pt[ptBinOm] = (TH1 *)f->FindObjectAny(TString::Format(("resultParOmm_pt[%d]"), ptBinOm));
            resultParOmC_pt[ptBinOm] = (TH1 *)f->FindObjectAny(TString::Format(("resultParOmC_pt[%d]"), ptBinOm));
        }

        f->Close();
    } /// Input ended!

    /// Saving output :
    TString output = outputFilename;
    if (outputFilename.IsNull())
    {
        // output = input;
        // output.ReplaceAll("Fit", "Draw");
        Printf("Empty output filename, aborting...");
        return 2;
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
            eff_pt_xim->SetBinContent(ptBinXi + 1, resultParXim_pt[ptBinXi]->GetBinContent(6));
            eff_pt_xip->SetBinContent(ptBinXi + 1, resultParXip_pt[ptBinXi]->GetBinContent(6));
            eff_pt_xiC->SetBinContent(ptBinXi + 1, resultParXiC_pt[ptBinXi]->GetBinContent(6));

            eff_pt_xim->SetBinError(ptBinXi + 1, resultParXim_pt[ptBinXi]->GetBinError(6));
            eff_pt_xip->SetBinError(ptBinXi + 1, resultParXip_pt[ptBinXi]->GetBinError(6));
            eff_pt_xiC->SetBinError(ptBinXi + 1, resultParXiC_pt[ptBinXi]->GetBinError(6));
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
            eff_pt_omm->SetBinContent(ptBinOm + 1, resultParOmm_pt[ptBinOm]->GetBinContent(6));
            eff_pt_omp->SetBinContent(ptBinOm + 1, resultParOmp_pt[ptBinOm]->GetBinContent(6));
            eff_pt_omC->SetBinContent(ptBinOm + 1, resultParOmC_pt[ptBinOm]->GetBinContent(6));

            eff_pt_omm->SetBinError(ptBinOm + 1, resultParOmm_pt[ptBinOm]->GetBinError(6));
            eff_pt_omp->SetBinError(ptBinOm + 1, resultParOmp_pt[ptBinOm]->GetBinError(6));
            eff_pt_omC->SetBinError(ptBinOm + 1, resultParOmC_pt[ptBinOm]->GetBinError(6));
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

    /// XI:
    h_MassXim->GetXaxis()->SetTitle("#Lambda^{0}-#pi^{#pm} Inv. Mass (GeV/c^{2})");
    h_MassXip->GetXaxis()->SetTitle("#Lambda^{0}-#pi^{#pm} Inv. Mass (GeV/c^{2})");
    h_MassXiC->GetXaxis()->SetTitle("#Lambda^{0}-#pi^{#pm} Inv. Mass (GeV/c^{2})");

    DrawAndSave(h_MassXim, nullptr, resultParams_Xim_allInt, saveImages, outputFolder, imageFormat);
    DrawAndSave(h_MassXip, nullptr, resultParams_Xip_allInt, saveImages, outputFolder, imageFormat);
    DrawAndSave(h_MassXiC, nullptr, resultParams_XiC_allInt, saveImages, outputFolder, imageFormat);

    /// Begin mult integrated part
    for (Int_t ptBinXi = 0; ptBinXi < nptbins_Xi; ptBinXi++)
    {
        h_MassXim_pt[ptBinXi]->GetXaxis()->SetTitle("#Lambda^{0}-#pi^{#pm} Inv. Mass (GeV/c^{2})");
        h_MassXip_pt[ptBinXi]->GetXaxis()->SetTitle("#Lambda^{0}-#pi^{#pm} Inv. Mass (GeV/c^{2})");
        h_MassXiC_pt[ptBinXi]->GetXaxis()->SetTitle("#Lambda^{0}-#pi^{#pm} Inv. Mass (GeV/c^{2})");
        DrawAndSave(h_MassXim_pt[ptBinXi], nullptr, resultParXim_pt[ptBinXi], saveImages, outputFolder, imageFormat);
        DrawAndSave(h_MassXip_pt[ptBinXi], nullptr, resultParXip_pt[ptBinXi], saveImages, outputFolder, imageFormat);
        DrawAndSave(h_MassXiC_pt[ptBinXi], nullptr, resultParXiC_pt[ptBinXi], saveImages, outputFolder, imageFormat);
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
                eff_xim[multBinXi]->SetBinContent(ptBinXi + 1, resultParXim_pt_mult[ptBinXi][multBinXi]->GetBinContent(6)); // Bin 6 in resultparams is MC efficiency
                eff_xim[multBinXi]->SetBinError(ptBinXi + 1, resultParXim_pt_mult[ptBinXi][multBinXi]->GetBinError(6));

                eff_xip[multBinXi]->SetBinContent(ptBinXi + 1, resultParXip_pt_mult[ptBinXi][multBinXi]->GetBinContent(6)); // Bin 6 in resultparams is MC efficiency
                eff_xip[multBinXi]->SetBinError(ptBinXi + 1, resultParXip_pt_mult[ptBinXi][multBinXi]->GetBinError(6));

                eff_xiC[multBinXi]->SetBinContent(ptBinXi + 1, resultParXiC_pt_mult[ptBinXi][multBinXi]->GetBinContent(6)); // Bin 6 in resultparams is MC efficiency
                eff_xiC[multBinXi]->SetBinError(ptBinXi + 1, resultParXiC_pt_mult[ptBinXi][multBinXi]->GetBinError(6));
            }

            /// Scale for bin width: N->dN/dpt
            // rawPt_xim[multBinXi]->Scale(1, "width");
            // rawPt_xip[multBinXi]->Scale(1, "width");

            h_MassXim_pt_mult[ptBinXi][multBinXi]->GetXaxis()->SetTitle("#Lambda^{0}-#pi^{#pm} Inv. Mass (GeV/c^{2})");
            h_MassXip_pt_mult[ptBinXi][multBinXi]->GetXaxis()->SetTitle("#Lambda^{0}-#pi^{#pm} Inv. Mass (GeV/c^{2})");
            h_MassXiC_pt_mult[ptBinXi][multBinXi]->GetXaxis()->SetTitle("#Lambda^{0}-#pi^{#pm} Inv. Mass (GeV/c^{2})");

            // TF1 *f1 = (TF1 *)h_MassXim_pt_mult[ptBinXi][multBinXi]->GetListOfFunctions()->At(0); /// 0 position = peak Function (GausPol2)
            DrawAndSave(h_MassXim_pt_mult[ptBinXi][multBinXi], h_bgXim_pt_mult[ptBinXi][multBinXi], resultParXim_pt_mult[ptBinXi][multBinXi], saveImages, outputFolder, imageFormat);
            // TF1 *f2 = (TF1 *)h_MassXip_pt_mult[ptBinXi][multBinXi]->GetListOfFunctions()->At(0); /// peak Function (GausPol2)
            DrawAndSave(h_MassXip_pt_mult[ptBinXi][multBinXi], h_bgXip_pt_mult[ptBinXi][multBinXi], resultParXip_pt_mult[ptBinXi][multBinXi], saveImages, outputFolder, imageFormat);
            // TF1 *f3 = (TF1 *)h_MassXiC_pt_mult[ptBinXi][multBinXi]->GetListOfFunctions()->At(0); /// peak Function (GausPol2)
            DrawAndSave(h_MassXiC_pt_mult[ptBinXi][multBinXi], h_bgXiC_pt_mult[ptBinXi][multBinXi], resultParXiC_pt_mult[ptBinXi][multBinXi], saveImages, outputFolder, imageFormat);
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

    /// OMEGA:
    h_MassOmm->GetXaxis()->SetTitle("#Lambda^{0}-K^{#pm} Inv. Mass (GeV/c^{2})");
    h_MassOmp->GetXaxis()->SetTitle("#Lambda^{0}-K^{#pm} Inv. Mass (GeV/c^{2})");
    h_MassOmC->GetXaxis()->SetTitle("#Lambda^{0}-K^{#pm} Inv. Mass (GeV/c^{2})");
    DrawAndSave(h_MassOmm, nullptr, resultParams_Omm_allInt, saveImages, outputFolder, imageFormat);
    DrawAndSave(h_MassOmp, nullptr, resultParams_Omp_allInt, saveImages, outputFolder, imageFormat);
    DrawAndSave(h_MassOmC, nullptr, resultParams_OmC_allInt, saveImages, outputFolder, imageFormat);

    for (Int_t ptBinOm = 0; ptBinOm < nptbins_Om; ptBinOm++)
    {
        h_MassOmm_pt[ptBinOm]->GetXaxis()->SetTitle("#Lambda^{0}-K^{#pm} Inv. Mass (GeV/c^{2})");
        h_MassOmp_pt[ptBinOm]->GetXaxis()->SetTitle("#Lambda^{0}-K^{#pm} Inv. Mass (GeV/c^{2})");
        h_MassOmC_pt[ptBinOm]->GetXaxis()->SetTitle("#Lambda^{0}-K^{#pm} Inv. Mass (GeV/c^{2})");

        DrawAndSave(h_MassOmm_pt[ptBinOm], nullptr, resultParOmm_pt[ptBinOm], saveImages, outputFolder, imageFormat);
        DrawAndSave(h_MassOmp_pt[ptBinOm], nullptr, resultParOmp_pt[ptBinOm], saveImages, outputFolder, imageFormat);
        DrawAndSave(h_MassOmC_pt[ptBinOm], nullptr, resultParOmC_pt[ptBinOm], saveImages, outputFolder, imageFormat);
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
                eff_omm[multBinOm]->SetBinContent(ptBinOm + 1, resultParOmm_pt_mult[ptBinOm][multBinOm]->GetBinContent(6)); // Bin 6 in resultparams is MC efficiency
                eff_omm[multBinOm]->SetBinError(ptBinOm + 1, resultParOmm_pt_mult[ptBinOm][multBinOm]->GetBinError(6));

                eff_omp[multBinOm]->SetBinContent(ptBinOm + 1, resultParOmp_pt_mult[ptBinOm][multBinOm]->GetBinContent(6)); // Bin 6 in resultparams is MC efficiency
                eff_omp[multBinOm]->SetBinError(ptBinOm + 1, resultParOmp_pt_mult[ptBinOm][multBinOm]->GetBinError(6));

                eff_omC[multBinOm]->SetBinContent(ptBinOm + 1, resultParOmC_pt_mult[ptBinOm][multBinOm]->GetBinContent(6)); // Bin 6 in resultparams is MC efficiency
                eff_omC[multBinOm]->SetBinError(ptBinOm + 1, resultParOmC_pt_mult[ptBinOm][multBinOm]->GetBinError(6));
            }
            /// Scale for bin width: N->dN/dpt
            // rawPt_omm[multBinOm]->Scale(1, "width");
            // rawPt_omp[multBinOm]->Scale(1, "width");

            h_MassOmm_pt_mult[ptBinOm][multBinOm]->GetXaxis()->SetTitle("#Lambda^{0}-K^{#pm} Inv. Mass (GeV/c^{2})");
            h_MassOmp_pt_mult[ptBinOm][multBinOm]->GetXaxis()->SetTitle("#Lambda^{0}-K^{#pm} Inv. Mass (GeV/c^{2})");
            h_MassOmC_pt_mult[ptBinOm][multBinOm]->GetXaxis()->SetTitle("#Lambda^{0}-K^{#pm} Inv. Mass (GeV/c^{2})");

            // TF1 *f1 = (TF1 *)h_MassOmm_pt_mult[ptBinOm][multBinOm]->GetListOfFunctions()->At(0); /// peak Function (GausPol2)
            DrawAndSave(h_MassOmm_pt_mult[ptBinOm][multBinOm], h_bgOmm_pt_mult[ptBinOm][multBinOm], resultParOmm_pt_mult[ptBinOm][multBinOm], saveImages, outputFolder, imageFormat);
            // TF1 *f2 = (TF1 *)h_MassOmp_pt_mult[ptBinOm][multBinOm]->GetListOfFunctions()->At(0); /// peak Function (GausPol2)
            DrawAndSave(h_MassOmp_pt_mult[ptBinOm][multBinOm], h_bgOmp_pt_mult[ptBinOm][multBinOm], resultParOmp_pt_mult[ptBinOm][multBinOm], saveImages, outputFolder, imageFormat);
            // TF1 *f3 = (TF1 *)h_MassOmC_pt_mult[ptBinOm][multBinOm]->GetListOfFunctions()->At(0); /// peak Function (GausPol2)
            DrawAndSave(h_MassOmC_pt_mult[ptBinOm][multBinOm], h_bgOmC_pt_mult[ptBinOm][multBinOm], resultParOmC_pt_mult[ptBinOm][multBinOm], saveImages, outputFolder, imageFormat);
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

    PaintStack(*c1, *hs_xip);
    PaintStack(*c2, *hs_omp);
    PaintStack(*c3, *hs_xim);
    PaintStack(*c4, *hs_omm);
    PaintStack(*c5, *hs_xiC);
    PaintStack(*c6, *hs_omC);

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

        PaintStack(*e1, *hs_xip_eff, kFALSE, "Efficiency");
        PaintStack(*e2, *hs_omp_eff, kFALSE, "Efficiency");
        PaintStack(*e3, *hs_xim_eff, kFALSE, "Efficiency");
        PaintStack(*e4, *hs_omm_eff, kFALSE, "Efficiency");
        PaintStack(*e5, *hs_xiC_eff, kFALSE, "Efficiency");
        PaintStack(*e6, *hs_omC_eff, kFALSE, "Efficiency");

        PaintStack(*r1, *hs_xiC_eff_ratio, kFALSE, "Multiplicity classes/0-100%");
        PaintStack(*r2, *hs_omC_eff_ratio, kFALSE, "Multiplicity classes/0-100%");

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

    if (!ptRatioFilename.IsNull())
    {
        /// get previusly calculated pT ratio from file to compute ratio of pT spectra

        Printf("\nOpening %s for ratio calculation...", ptRatioFilename.Data());
        TFile *fRatio = TFile::Open(ptRatioFilename.Data());
        if (!fRatio)
        {
            Printf("Error: Cannot open file '%s' !", ptRatioFilename.Data());
            return 1;
        }

        TString yAxisTitle = outputFilename(outputFilename.Last('/') + 1, outputFilename.Length()) + "/" + ptRatioFilename(ptRatioFilename.Last('/') + 1, ptRatioFilename.Length());
        yAxisTitle.ReplaceAll(".root", "");
        yAxisTitle.ReplaceAll("_draw", "");

        // needed so we can do file->Close()
        TH1::AddDirectory(0);
        TH1 *rawPt_xim_compare[nmultbins_Xi];
        TH1 *rawPt_xip_compare[nmultbins_Xi];
        TH1 *rawPt_omm_compare[nmultbins_Om];
        TH1 *rawPt_omp_compare[nmultbins_Om];
        TH1 *rawPt_xiC_compare[nmultbins_Xi];
        TH1 *rawPt_omC_compare[nmultbins_Om];

        TH1D *ratioPt_xim[nmultbins_Xi];
        TH1D *ratioPt_xip[nmultbins_Xi];
        TH1D *ratioPt_omm[nmultbins_Om];
        TH1D *ratioPt_omp[nmultbins_Om];
        TH1D *ratioPt_xiC[nmultbins_Xi];
        TH1D *ratioPt_omC[nmultbins_Om];

        auto hs_ratio_xip = new THStack("hs_ratio_xip", "#it{p}_{T} spectra ratio #Xi^{+}");
        auto hs_ratio_xim = new THStack("hs_ratio_xim", "#it{p}_{T} spectra ratio #Xi^{-}");
        auto hs_ratio_omp = new THStack("hs_ratio_omp", "#it{p}_{T} spectra ratio #Omega^{+}");
        auto hs_ratio_omm = new THStack("hs_ratio_omm", "#it{p}_{T} spectra ratio #Omega^{-}");
        auto hs_ratio_xiC = new THStack("hs_ratio_xiC", "#it{p}_{T} spectra ratio #Xi^{+} + #Xi^{-}");
        auto hs_ratio_omC = new THStack("hs_ratio_omC", "#it{p}_{T} spectra ratio #Omega^{+} + #Omega^{-}");

        for (Int_t multBinXi = 0; multBinXi < nmultbins_Xi; multBinXi++)
        {
            rawPt_xim_compare[multBinXi] = (TH1 *)fRatio->FindObjectAny(TString::Format(("rawPt_xim[%d]"), multBinXi));
            rawPt_xip_compare[multBinXi] = (TH1 *)fRatio->FindObjectAny(TString::Format(("rawPt_xip[%d]"), multBinXi));
            rawPt_xiC_compare[multBinXi] = (TH1 *)fRatio->FindObjectAny(TString::Format(("rawPt_xiC[%d]"), multBinXi));
        }

        for (Int_t multBinOm = 0; multBinOm < nmultbins_Om; multBinOm++)
        {
            rawPt_omm_compare[multBinOm] = (TH1 *)fRatio->FindObjectAny(TString::Format(("rawPt_omm[%d]"), multBinOm));
            rawPt_omp_compare[multBinOm] = (TH1 *)fRatio->FindObjectAny(TString::Format(("rawPt_omp[%d]"), multBinOm));
            rawPt_omC_compare[multBinOm] = (TH1 *)fRatio->FindObjectAny(TString::Format(("rawPt_omC[%d]"), multBinOm));
        }
        fRatio->Close(); /// close file for input

        /// calculate ratio of current spectra to previous
        for (Int_t multBinXi = 0; multBinXi < nmultbins_Xi; multBinXi++)
        {
            /// get current histogram scale back to normal to perform divide
            rawPt_xim[multBinXi]->Scale(pow(2, -((nmultbins_Xi - 1) - multBinXi)));
            rawPt_xip[multBinXi]->Scale(pow(2, -((nmultbins_Xi - 1) - multBinXi)));
            rawPt_xiC[multBinXi]->Scale(pow(2, -((nmultbins_Xi - 1) - multBinXi)));

            /// generate ratio hists
            ratioPt_xim[multBinXi] = new TH1D(TString::Format(("ratioPt_xim[%d]"), multBinXi), "", nptbins_Xi, ptbins_Xi);
            ratioPt_xip[multBinXi] = new TH1D(TString::Format(("ratioPt_xip[%d]"), multBinXi), "", nptbins_Xi, ptbins_Xi);
            ratioPt_xiC[multBinXi] = new TH1D(TString::Format(("ratioPt_xiC[%d]"), multBinXi), "", nptbins_Xi, ptbins_Xi);

            ratioPt_xim[multBinXi]->Divide(rawPt_xim[multBinXi], rawPt_xim_compare[multBinXi]);
            out->cd("dirRawPt_xim");
            ratioPt_xim[multBinXi]->SetMarkerStyle(markerStyles[multBinXi]);
            ratioPt_xim[multBinXi]->Write();
            ratioPt_xim[multBinXi]->SetName(TString::Format(("Mult: %.0f-%.0f%%"), multbins_Xi[multBinXi], multbins_Xi[multBinXi + 1]));
            hs_ratio_xim->Add(ratioPt_xim[multBinXi]);

            ratioPt_xip[multBinXi]->Divide(rawPt_xip[multBinXi], rawPt_xip_compare[multBinXi]);
            out->cd("dirRawPt_xip");
            ratioPt_xip[multBinXi]->SetMarkerStyle(markerStyles[multBinXi]);
            ratioPt_xip[multBinXi]->Write();
            ratioPt_xip[multBinXi]->SetName(TString::Format(("Mult: %.0f-%.0f%%"), multbins_Xi[multBinXi], multbins_Xi[multBinXi + 1]));
            hs_ratio_xip->Add(ratioPt_xip[multBinXi]);

            ratioPt_xiC[multBinXi]->Divide(rawPt_xiC[multBinXi], rawPt_xiC_compare[multBinXi]);
            out->cd("dirRawPt_xiC");
            ratioPt_xiC[multBinXi]->SetMarkerStyle(markerStyles[multBinXi]);
            ratioPt_xiC[multBinXi]->Write();
            ratioPt_xiC[multBinXi]->SetName(TString::Format(("Mult: %.0f-%.0f%%"), multbins_Xi[multBinXi], multbins_Xi[multBinXi + 1]));
            hs_ratio_xiC->Add(ratioPt_xiC[multBinXi]);
        }

        for (Int_t multBinOm = 0; multBinOm < nmultbins_Om; multBinOm++)
        {
            /// get current histogram scale back to normal to perform divide
            rawPt_omm[multBinOm]->Scale(pow(2, -((nmultbins_Om - 1) - multBinOm)));
            rawPt_omp[multBinOm]->Scale(pow(2, -((nmultbins_Om - 1) - multBinOm)));
            rawPt_omC[multBinOm]->Scale(pow(2, -((nmultbins_Om - 1) - multBinOm)));

            /// generate ratio hists
            ratioPt_omm[multBinOm] = new TH1D(TString::Format(("ratioPt_omm[%d]"), multBinOm), "", nptbins_Om, ptbins_Om);
            ratioPt_omp[multBinOm] = new TH1D(TString::Format(("ratioPt_omm[%d]"), multBinOm), "", nptbins_Om, ptbins_Om);
            ratioPt_omC[multBinOm] = new TH1D(TString::Format(("ratioPt_omm[%d]"), multBinOm), "", nptbins_Om, ptbins_Om);

            ratioPt_omm[multBinOm]->Divide(rawPt_omm[multBinOm], rawPt_omm_compare[multBinOm]);
            out->cd("dirRawPt_omm");
            ratioPt_omm[multBinOm]->SetMarkerStyle(markerStyles[multBinOm]);
            ratioPt_omm[multBinOm]->Write();
            ratioPt_omm[multBinOm]->SetName(TString::Format(("Mult: %.0f-%.0f%%"), multbins_Om[multBinOm], multbins_Om[multBinOm + 1]));
            hs_ratio_omm->Add(ratioPt_omm[multBinOm]);

            ratioPt_omp[multBinOm]->Divide(rawPt_omp[multBinOm], rawPt_omp_compare[multBinOm]);
            out->cd("dirRawPt_omp");
            ratioPt_omp[multBinOm]->SetMarkerStyle(markerStyles[multBinOm]);
            ratioPt_omp[multBinOm]->Write();
            ratioPt_omp[multBinOm]->SetName(TString::Format(("Mult: %.0f-%.0f%%"), multbins_Om[multBinOm], multbins_Om[multBinOm + 1]));
            hs_ratio_omp->Add(ratioPt_omp[multBinOm]);

            ratioPt_omC[multBinOm]->Divide(rawPt_omC[multBinOm], rawPt_omC_compare[multBinOm]);
            out->cd("dirRawPt_omC");
            ratioPt_omC[multBinOm]->SetMarkerStyle(markerStyles[multBinOm]);
            ratioPt_omC[multBinOm]->Write();
            ratioPt_omC[multBinOm]->SetName(TString::Format(("Mult: %.0f-%.0f%%"), multbins_Om[multBinOm], multbins_Om[multBinOm + 1]));
            hs_ratio_omC->Add(ratioPt_omC[multBinOm]);
        }

        TCanvas *cRatio1 = new TCanvas("cRatio1", "cRatio1", 1920, 1080);
        TCanvas *cRatio2 = new TCanvas("cRatio2", "cRatio2", 1920, 1080);
        TCanvas *cRatio3 = new TCanvas("cRatio3", "cRatio3", 1920, 1080);
        TCanvas *cRatio4 = new TCanvas("cRatio4", "cRatio4", 1920, 1080);
        TCanvas *cRatio5 = new TCanvas("cRatio5", "cRatio5", 1920, 1080);
        TCanvas *cRatio6 = new TCanvas("cRatio6", "cRatio6", 1920, 1080);

        PaintStack(*cRatio1, *hs_ratio_xip, kFALSE, yAxisTitle);
        PaintStack(*cRatio2, *hs_ratio_omp, kFALSE, yAxisTitle);
        PaintStack(*cRatio3, *hs_ratio_xim, kFALSE, yAxisTitle);
        PaintStack(*cRatio4, *hs_ratio_omm, kFALSE, yAxisTitle);
        PaintStack(*cRatio5, *hs_ratio_xiC, kFALSE, yAxisTitle);
        PaintStack(*cRatio6, *hs_ratio_omC, kFALSE, yAxisTitle);
        {
            /// Draw a line at y=1 for ratio histStack
            TLine *lLineAt1 = new TLine(0.8, 1, 5.3, 1);
            lLineAt1->SetLineColor(kRed);

            cRatio1->cd();
            lLineAt1->Draw("same");
            cRatio2->cd();
            lLineAt1->Draw("same");
            cRatio3->cd();
            lLineAt1->Draw("same");
            cRatio4->cd();
            lLineAt1->Draw("same");
            cRatio5->cd();
            lLineAt1->Draw("same");
            cRatio6->cd();
            lLineAt1->Draw("same");
        }
        hs_ratio_xim->SetMinimum(0.7);
        hs_ratio_xip->SetMinimum(0.7);
        hs_ratio_xiC->SetMinimum(0.7);
        hs_ratio_omm->SetMinimum(0.5);
        hs_ratio_omp->SetMinimum(0.5);
        hs_ratio_omC->SetMinimum(0.5);

        hs_ratio_xim->SetMaximum(1.1);
        hs_ratio_xip->SetMaximum(1.1);
        hs_ratio_xiC->SetMaximum(1.1);
        hs_ratio_omm->SetMaximum(1.2);
        hs_ratio_omp->SetMaximum(1.2);
        hs_ratio_omC->SetMaximum(1.2);

        out->cd();
        hs_ratio_xip->Write();
        hs_ratio_omp->Write();
        hs_ratio_xim->Write();
        hs_ratio_omm->Write();
        hs_ratio_xiC->Write();
        hs_ratio_omC->Write();

        if (saveStack)
        {
            cRatio1->cd();
            SaveImage(Form("%s/images/ratioPt", outputFolder.Data()), "ratio_xipN", imageFormat.Data());

            cRatio2->cd();
            SaveImage(Form("%s/images/ratioPt", outputFolder.Data()), "ratio_ompN", imageFormat.Data());

            cRatio3->cd();
            SaveImage(Form("%s/images/ratioPt", outputFolder.Data()), "ratio_ximN", imageFormat.Data());

            cRatio4->cd();
            SaveImage(Form("%s/images/ratioPt", outputFolder.Data()), "ratio_ommN", imageFormat.Data());

            cRatio5->cd();
            SaveImage(Form("%s/images/ratioPt", outputFolder.Data()), "ratio_xicN", imageFormat.Data());

            cRatio6->cd();
            SaveImage(Form("%s/images/ratioPt", outputFolder.Data()), "ratio_omcN", imageFormat.Data());
        }
    }

    out->Close();
    // delete out;
    return 0;
}
void DrawAndSave(TH1 *peak, TH1 *bg, TH1 *resultParams, Bool_t saveImages, TString outputFolder, TString imageFormat)
{

    gROOT->SetBatch(kTRUE);
    TCanvas *c1 = new TCanvas(peak->GetName(), peak->GetTitle(), 1920, 1080);

    /// TODO: CREATE CUSTOM LEGEND SHOWING AVG MASS AND SIGMA ON GRAPH
    Double_t pPosition = resultParams->GetBinContent(3); // avg mean for DGaus fit is stored in resultParams bin 3
    Double_t pWidth = resultParams->GetBinContent(4);    // avg sigma for DGaus fit is stored in resultParams bin 4
    auto legend = new TLegend(0.1, 0.7, 0.28, 0.9);
    legend->SetHeader("Fit Stats", "C"); // option "C" allows to center the header
    legend->AddEntry(peak->GetListOfFunctions()->At(0), "", "l");
    legend->AddEntry(peak->GetListOfFunctions()->At(0), TString::Format("Fit mean = %f +/- %f", pPosition, resultParams->GetBinError(3)), "l");
    legend->AddEntry(peak->GetListOfFunctions()->At(0), TString::Format("Fit sigma = %f +/- %f", pWidth, resultParams->GetBinError(4)), "l");
    // legend->AddEntry(peak->GetListOfFunctions()->At(1), "", "lpf");
    legend->AddEntry(peak, TString::Format("Sig - Bg (BC-FF) = %f +/- %f", resultParams->GetBinContent(1), resultParams->GetBinError(1)), "pe");

    ///  Defining peak limits for signal region (green lines):
    ///   par[1] = peak position, par[2] = peak width
    Double_t lPeakLeftLimit = pPosition - 1. * 4 * TMath::Abs(pWidth);
    Double_t lPeakRightLimit = pPosition + 1. * 4 * TMath::Abs(pWidth);
    TLine *lLineLeft = new TLine(lPeakLeftLimit, 0, lPeakLeftLimit, peak->GetMaximum());
    TLine *lLineRight = new TLine(lPeakRightLimit, 0, lPeakRightLimit, peak->GetMaximum());
    lLineLeft->SetLineColor(kMagenta);
    lLineRight->SetLineColor(kMagenta);
    peak->Draw();
    if (bg)
        bg->Draw("same");
    lLineLeft->Draw("same");
    lLineRight->Draw("same");
    legend->Draw();

    if (saveImages)
    {
        TString imageFolder = peak->GetName();
        if (imageFolder.Contains("["))
            imageFolder.Remove(imageFolder.First('[')); /// getting substring for folder naming purpose
        else
            imageFolder = "h_allInt";
        SaveImage(Form("%s/images/%s", outputFolder.Data(), imageFolder.Data()), peak->GetName(), imageFormat.Data());
    }
}
