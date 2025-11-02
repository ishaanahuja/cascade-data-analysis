/// DEPRECATED AND OBSOLETE. Read note in function definition.

#include <fstream>
#include <vector>

#include <TCanvas.h>
#include <TDirectory.h>
#include <TFile.h>
#include <TF1.h>
#include <TH1.h>
#include <THStack.h>
#include <TLine.h>
#include <TROOT.h>
#include <TString.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TVirtualPad.h>

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

/**
 * @brief This function calculates and corrects the efficiency of Xi and Omega particles from Monte Carlo (MC) simulations.
 *
 * @param inputFileList The name of the file containing the list of MC files to be processed. Default is "McFileList.txt".
 * @param outputFilename The name of the output ROOT file where results will be saved. Default is "test_CorrectedSpectra.root".
 * @param outputFolder The folder where output images and results will be saved. Default is "rootResults/tstCorrected".
 * @param fGenerateCorrected A flag to indicate whether to generate corrected spectra. Default is kTRUE.
 * @param input_RawPt The path to the input ROOT file containing raw Pt spectra. Default is "/var/home/ishaan/Work/git/analysis/p-Pb/MC/rootResults/6runs/testDraw_all.root".
 * @param saveStack A flag to indicate whether to save the stack plots as images. Default is kTRUE.
 * @param imageFormat The format of the output images. Default is "png".
 *
 * @return int Returns 0 on success, 1 if there is an error opening a file or finding a histogram, and 3 if there is an error opening the output file.
 *
 * The function performs the following steps:
 * 1. Initializes histograms and directories.
 * 2. Reads the list of MC files and extracts efficiency histograms for Xi and Omega particles.
 * 3. Calculates the average efficiency for Xi and Omega particles.
 * 4. Saves the average efficiency histograms to the output file.
 * 5. Draws and saves stack plots of the efficiencies.
 * 6. If fGenerateCorrected is true, applies efficiency corrections to raw Pt spectra and generates corrected spectra.
 * 7. Saves the corrected spectra to the output file and optionally saves the stack plots as images.
 *
 * The SaveImage function saves the current pad as an image file in the specified format.
 *
 * @warning Obsolete and replaced by EfficiencyEstimation.C and EfficiencyCorrection.C
 * @note Might work, not supported anymore, use the above mentioned scripts for efficiency estimation and correction.
 *
 */
int CorrectedSpectra(
    std::string inputFileList = "McFileList.txt",
    TString outputFilename = "test_CorrectedSpectra.root",
    TString outputFolder = "rootResults/tstCorrected",
    Bool_t fGenerateCorrected = kTRUE,
    TString input_RawPt = "/var/home/ishaan/Work/git/analysis/p-Pb/MC/rootResults/6runs/testDraw_all.root",
    Bool_t saveStack = kTRUE,
    TString imageFormat = "png")
{
    // keep histograms independent from their directory
    TH1::AddDirectory(0);
    TDirectory::AddDirectory(0);

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

    TH1 *eff_xiC;
    TH1 *eff_omC;
    std::vector<TH1 *> hvec_eff_xiC;
    std::vector<TH1 *> hvec_eff_omC;

    std::vector<std::string> MCname;

    auto hs_xiC_eff = new THStack("hs_xiC_eff", "Efficiency #Xi^{+} + #Xi^{-}");
    auto hs_omC_eff = new THStack("hs_omC_eff", "Efficiency #Omega^{+} + #Omega^{-}");

    /// Getting histograms:
    int totalFilesInFileList = 0;
    std::string fileName;

    Printf("Opening MC file list ...");
    std::ifstream fileList(inputFileList.data());

    // read MC filelist
    while (std::getline(fileList, fileName))
    {
        // find name of MC from parent folder naming scheme
        std::string delimiter = "/";
        size_t pos1, pos2 = 0;
        std::string MCrun = fileName;
        pos1 = MCrun.rfind(delimiter);
        MCrun.erase(pos1, std::string::npos);
        pos2 = MCrun.rfind(delimiter);
        MCrun = MCrun.substr(pos2 + delimiter.length(), std::string::npos);
        MCname.push_back(MCrun);

        Printf("%d: %s", totalFilesInFileList, MCrun.data());

        TFile *f = TFile::Open(fileName.data());
        if (!f)
        {
            Printf("Error: Cannot open file '%s' !", fileName.data());
            // return 1;
            continue;
        }

        /// getting efficiency histogram
        eff_xiC = (TH1 *)f->FindObjectAny("eff_pt_xiC");
        if (!eff_xiC)
        {
            Printf("Error: Cannot find histogram 'eff_pt_xiC' in '%s'!", fileName.data());
            return 1;
        }
        eff_xiC->SetNameTitle(TString::Format(("eff_xiC_[%d]"), totalFilesInFileList), TString::Format(("%s"), MCrun.data()));
        hvec_eff_xiC.push_back(eff_xiC);

        eff_omC = (TH1 *)f->FindObjectAny("eff_pt_omC");
        if (!eff_omC)
        {
            Printf("Error: Cannot find histogram 'eff_pt_omC' in '%s'!", fileName.data());
            return 1;
        }
        eff_omC->SetNameTitle((TString::Format(("eff_omC_[%d]"), totalFilesInFileList)), (TString::Format(("%s"), MCrun.data())));
        hvec_eff_omC.push_back(eff_omC);

        totalFilesInFileList++;

        /// Input ended!
        f->Close();
    }

    /// Calculate avg efficiency:
    Double_t weightsXi[totalFilesInFileList][nptbins_Xi];
    Double_t weightsOm[totalFilesInFileList][nptbins_Om];
    Double_t valsXi[totalFilesInFileList][nptbins_Xi];
    Double_t valsOm[totalFilesInFileList][nptbins_Om];

    Double_t effXiAvg[nptbins_Xi];
    Double_t effXiErrAvg[nptbins_Xi];
    Double_t effOmAvg[nptbins_Om];
    Double_t effOmErrAvg[nptbins_Om];

    TH1D *eff_xiC_avg = new TH1D("eff_xiC_avg", "#Xi^{#pm} Avg MC efficiency", nptbins_Xi, ptbins_Xi);
    TH1D *eff_omC_avg = new TH1D("eff_omC_avg", "#Omega^{#pm} Avg MC efficiency", nptbins_Om, ptbins_Om);

    /// Saving output :
    TString output = outputFilename;
    if (outputFilename.IsNull())
        output = "CorrectedSpectra.root";

    Printf("Saving output to '%s' ...", output.Data());
    TFile *out = TFile::Open(output.Data(), "RECREATE");
    if (!out)
    {
        Printf("Error: Cannot open file '%s' !", output.Data());
        return 3;
    }

    out->mkdir("eff_xiC");
    out->mkdir("eff_omC");
    /// Output set.

    /// Begin
    Int_t markerStyles[] = {31, 47, 21, 33, 22, 23, 29, 33, 34, 43, 47, 41, 20}; // chosen marker style palette

    for (Int_t i = 0; i < totalFilesInFileList; i++)
    {
        out->cd("eff_xiC");
        hvec_eff_xiC.at(i)->Write();
        hvec_eff_xiC.at(i)->SetMarkerStyle(markerStyles[i]);
        if (MCname.at(i).find("Om") == std::string::npos) // don't add it to histStack if MC run name contains omega
            hs_xiC_eff->Add(hvec_eff_xiC.at(i));

        out->cd("eff_omC");
        hvec_eff_omC.at(i)->Write();
        hvec_eff_omC.at(i)->SetMarkerStyle(markerStyles[i]);
        if (MCname.at(i).find("Xi") == std::string::npos) // don't add it to histStack if MC run name contains Xi
            hs_omC_eff->Add(hvec_eff_omC.at(i));
    }

    for (Int_t j = 0; j < totalFilesInFileList; j++)
    {
        if (MCname.at(j).find("Om") == std::string::npos) // only take Xi MC runs into consideration
        {
            for (Int_t ptBin = 1; ptBin <= nptbins_Xi; ptBin++)
            {
                valsXi[j][ptBin] = hvec_eff_xiC.at(j)->GetBinContent(ptBin);
                weightsXi[j][ptBin] = 1 / pow(hvec_eff_xiC.at(j)->GetBinError(ptBin), 2);
            }
        }
        if (MCname.at(j).find("Xi") == std::string::npos) // only take Omega MC runs into consideration
        {
            for (Int_t ptBin = 1; ptBin <= nptbins_Om; ptBin++)
            {
                valsOm[j][ptBin] = hvec_eff_omC.at(j)->GetBinContent(ptBin);
                weightsOm[j][ptBin] = 1 / pow(hvec_eff_omC.at(j)->GetBinError(ptBin), 2);
            }
        }
    }

    for (Int_t ptBin = 1; ptBin <= nptbins_Xi; ptBin++)
    {

        Double_t sumWeights = 0;
        Double_t sumWeightVals = 0;

        for (Int_t j = 0; j < totalFilesInFileList; j++)
        {
            if (MCname.at(j).find("Om") == std::string::npos) // only take Xi MC runs into consideration
            {
                sumWeights = sumWeights + weightsXi[j][ptBin];
                sumWeightVals = sumWeightVals + (weightsXi[j][ptBin] * valsXi[j][ptBin]);
            }
        }
        effXiAvg[ptBin] = sumWeightVals / sumWeights;
        effXiErrAvg[ptBin] = pow(sumWeights, -0.5);

        eff_xiC_avg->SetBinContent(ptBin, effXiAvg[ptBin]);
        eff_xiC_avg->SetBinError(ptBin, effXiErrAvg[ptBin]);
    }

    for (Int_t ptBin = 1; ptBin <= nptbins_Om; ptBin++)
    {

        Double_t sumWeights = 0;
        Double_t sumWeightVals = 0;

        for (Int_t j = 0; j < totalFilesInFileList; j++)
        {
            if (MCname.at(j).find("Xi") == std::string::npos) // only take Omega MC runs into consideration
            {
                sumWeights = sumWeights + weightsOm[j][ptBin];
                sumWeightVals = sumWeightVals + (weightsOm[j][ptBin] * valsOm[j][ptBin]);
            }
        }
        effOmAvg[ptBin] = sumWeightVals / sumWeights;
        effOmErrAvg[ptBin] = pow(sumWeights, -0.5);
        eff_omC_avg->SetBinContent(ptBin, effOmAvg[ptBin]);
        eff_omC_avg->SetBinError(ptBin, effOmErrAvg[ptBin]);
    }

    /// check if average is okay
    Double_t chi2Xi[nptbins_Xi];
    Double_t chi2Om[nptbins_Om];
    Double_t N_Xi = 0;
    Double_t N_Om = 0;

    for (Int_t ptBin = 1; ptBin <= nptbins_Xi; ptBin++)
    {
        Double_t sumChi2 = 0;
        N_Xi = 0;
        for (Int_t j = 0; j < totalFilesInFileList; j++)
        {
            if (MCname.at(j).find("Om") == std::string::npos) // only take Xi MC runs into consideration
            {
                sumChi2 = sumChi2 + weightsXi[j][ptBin] * pow(effXiAvg[ptBin] - valsXi[j][ptBin], 2);
                N_Xi++;
            }
        }
        chi2Xi[ptBin] = sumChi2;
    }
    for (Int_t ptBin = 1; ptBin <= nptbins_Om; ptBin++)
    {
        Double_t sumChi2 = 0;
        N_Om = 0;

        for (Int_t j = 0; j < totalFilesInFileList; j++)
        {
            if (MCname.at(j).find("Xi") == std::string::npos) // only take Omega MC runs into consideration
            {
                sumChi2 = sumChi2 + weightsOm[j][ptBin] * pow(effOmAvg[ptBin] - valsOm[j][ptBin], 2);
                N_Om++;
            }
        }
        chi2Om[ptBin] = sumChi2;
    }
    for (Int_t ptBin = 1; ptBin <= nptbins_Xi; ptBin++)
    {
        Printf("PtBin: %d chi2Xi = %f chi2Xi/N-1 = %f", ptBin, chi2Xi[ptBin], (chi2Xi[ptBin] / (N_Xi - 1)));
    }
    for (Int_t ptBin = 1; ptBin <= nptbins_Om; ptBin++)
    {
        Printf("PtBin: %d chi2Om = %f chi2Om/N-1 = %f", ptBin, chi2Om[ptBin], (chi2Om[ptBin] / (N_Om - 1)));
    }

    eff_xiC_avg->GetXaxis()->SetTitle("#it{p}_{T} (GeV/c)");
    eff_omC_avg->GetXaxis()->SetTitle("#it{p}_{T} (GeV/c)");
    eff_xiC_avg->GetYaxis()->SetTitle("Efficiency");
    eff_omC_avg->GetYaxis()->SetTitle("Efficiency");
    eff_xiC_avg->SetLineColor(kRed);
    eff_omC_avg->SetLineColor(kRed);
    eff_xiC_avg->SetLineWidth(2);
    eff_omC_avg->SetLineWidth(2);
    eff_xiC_avg->SetMarkerStyle(kFullCircle);
    eff_omC_avg->SetMarkerStyle(kFullCircle);

    out->cd();
    eff_xiC_avg->Write();
    eff_omC_avg->Write();
    hs_xiC_eff->Add(eff_xiC_avg);
    hs_omC_eff->Add(eff_omC_avg);

    /// Begin drawing
    gStyle->SetOptStat("i");
    // gStyle->SetPalette(kVisibleSpectrum);
    gStyle->SetPalette(kRainBow);

    TCanvas *c1 = new TCanvas("c1", "#Xi efficiency comparison", 1920, 1080);
    TCanvas *c2 = new TCanvas("c2", "#Omega efficiency comparison", 1920, 1080);

    c1->cd();
    hs_xiC_eff->Draw("PLC PMC NOSTACK");
    // gPad->SetLogy();
    gPad->BuildLegend(0.9, 0.74, 1., 0.9, "MC Run");
    hs_xiC_eff->GetXaxis()->SetTitle("#it{p}_{T} (GeV/c)");
    hs_xiC_eff->GetYaxis()->SetTitle("Efficiency");
    c1->Modified();
    c1->ForceUpdate();

    c2->cd();
    hs_omC_eff->Draw("PLC PMC NOSTACK");
    // gPad->SetLogy();
    gPad->BuildLegend(0.9, 0.74, 1., 0.9, "MC Run");
    hs_omC_eff->GetXaxis()->SetTitle("#it{p}_{T} (GeV/c)");
    hs_omC_eff->GetYaxis()->SetTitle("Efficiency");
    c2->Modified();
    c2->ForceUpdate();

    c1->Modified();
    c1->ForceUpdate();
    c2->Modified();
    c2->ForceUpdate();

    hs_xiC_eff->Write();
    hs_omC_eff->Write();

    if (saveStack)
    {
        c1->cd();
        SaveImage(Form("%s/effComp", outputFolder.Data()), "xiC_eff_final", imageFormat.Data());

        c2->cd();
        SaveImage(Form("%s/effComp", outputFolder.Data()), "omC_eff_final", imageFormat.Data());
    }

    if (fGenerateCorrected)
    {
        Printf("Applying efficiency corrections and generating corrected spectra...");

        TH1D *rawPt_xim[nmultbins_Xi];
        TH1D *rawPt_xip[nmultbins_Xi];
        TH1D *rawPt_omm[nmultbins_Om];
        TH1D *rawPt_omp[nmultbins_Om];
        TH1D *rawPt_xiC[nmultbins_Xi];
        TH1D *rawPt_omC[nmultbins_Om];

        TH1D *corrPt_xim[nmultbins_Xi];
        TH1D *corrPt_xip[nmultbins_Xi];
        TH1D *corrPt_omm[nmultbins_Om];
        TH1D *corrPt_omp[nmultbins_Om];
        TH1D *corrPt_xiC[nmultbins_Xi];
        TH1D *corrPt_omC[nmultbins_Om];

        auto hs_xip_corrected = new THStack("hs_xip_corrected", "Efficiency corrected #it{p}_{T} spectra #Xi^{+}");
        auto hs_xim_corrected = new THStack("hs_xim_corrected", "Efficiency corrected #it{p}_{T} spectra #Xi^{-}");
        auto hs_omp_corrected = new THStack("hs_omp_corrected", "Efficiency corrected #it{p}_{T} spectra #Omega^{+}");
        auto hs_omm_corrected = new THStack("hs_omm_corrected", "Efficiency corrected #it{p}_{T} spectra #Omega^{-}");
        auto hs_xiC_corrected = new THStack("hs_xiC_corrected", "Efficiency corrected #it{p}_{T} spectra #Xi^{+} + #Xi^{-}");
        auto hs_omC_corrected = new THStack("hs_omC_corrected", "Efficiency corrected #it{p}_{T} spectra #Omega^{+} + #Omega^{-}");

        Printf("Opening file '%s'... ", input_RawPt.Data());
        TFile *f = TFile::Open(input_RawPt.Data());
        if (!f)
        {
            Printf("Error: Cannot open file '%s' !", input_RawPt.Data());
            return 1;
        }

        for (Int_t multBinXi = 0; multBinXi < nmultbins_Xi; multBinXi++)
        {
            rawPt_xim[multBinXi] = (TH1D *)f->FindObjectAny(TString::Format(("rawPt_xim[%d]"), multBinXi));
            rawPt_xip[multBinXi] = (TH1D *)f->FindObjectAny(TString::Format(("rawPt_xip[%d]"), multBinXi));
            rawPt_xiC[multBinXi] = (TH1D *)f->FindObjectAny(TString::Format(("rawPt_xiC[%d]"), multBinXi));
        }
        for (Int_t multBinOm = 0; multBinOm < nmultbins_Om; multBinOm++)
        {
            rawPt_omm[multBinOm] = (TH1D *)f->FindObjectAny(TString::Format(("rawPt_omm[%d]"), multBinOm));
            rawPt_omp[multBinOm] = (TH1D *)f->FindObjectAny(TString::Format(("rawPt_omp[%d]"), multBinOm));
            rawPt_omC[multBinOm] = (TH1D *)f->FindObjectAny(TString::Format(("rawPt_omC[%d]"), multBinOm));
        }

        f->Close();

        out->mkdir("CorrectedSpectra");
        out->cd("CorrectedSpectra");
        gDirectory->mkdir("Xim");
        gDirectory->mkdir("Xip");
        gDirectory->mkdir("Omm");
        gDirectory->mkdir("Omp");
        gDirectory->mkdir("XiC");
        gDirectory->mkdir("OmC");

        for (Int_t multBinXi = 0; multBinXi < nmultbins_Xi; multBinXi++)
        {
            out->cd("CorrectedSpectra");
            gDirectory->cd("Xim");
            corrPt_xim[multBinXi] = new TH1D(TString::Format(("corrPt_xim[%d]"), multBinXi), TString::Format(("Eff. corrected #Xi^{-}: Mult: %.0f-%.0f%%"), multbins_Xi[multBinXi], multbins_Xi[multBinXi + 1]), nptbins_Xi, ptbins_Xi);
            // corrPt_xim[multBinXi]->Sumw2();
            corrPt_xim[multBinXi]->Divide(rawPt_xim[multBinXi], eff_xiC_avg);
            corrPt_xim[multBinXi]->SetMarkerStyle(markerStyles[multBinXi]);
            corrPt_xim[multBinXi]->Write();

            corrPt_xim[multBinXi]->SetTitle(TString::Format(("#Xi^{-}: Mult: %.0f-%.0f%%"), multbins_Xi[multBinXi], multbins_Xi[multBinXi + 1]));
            corrPt_xim[multBinXi]->Scale(pow(2, (nmultbins_Xi - 1) - multBinXi));
            hs_xim_corrected->Add(corrPt_xim[multBinXi]);

            out->cd("CorrectedSpectra");
            gDirectory->cd("Xip");
            corrPt_xip[multBinXi] = new TH1D(TString::Format(("corrPt_xip[%d]"), multBinXi), TString::Format(("Eff. corrected #Xi^{+}: Mult: %.0f-%.0f%%"), multbins_Xi[multBinXi], multbins_Xi[multBinXi + 1]), nptbins_Xi, ptbins_Xi);
            // corrPt_xip[multBinXi]->Sumw2();
            corrPt_xip[multBinXi]->Divide(rawPt_xip[multBinXi], eff_xiC_avg);
            corrPt_xip[multBinXi]->SetMarkerStyle(markerStyles[multBinXi]);
            corrPt_xip[multBinXi]->Write();

            corrPt_xip[multBinXi]->SetTitle(TString::Format(("#Xi^{+}: Mult: %.0f-%.0f%%"), multbins_Xi[multBinXi], multbins_Xi[multBinXi + 1]));
            corrPt_xip[multBinXi]->Scale(pow(2, (nmultbins_Xi - 1) - multBinXi));
            hs_xip_corrected->Add(corrPt_xip[multBinXi]);

            out->cd("CorrectedSpectra");
            gDirectory->cd("XiC");
            corrPt_xiC[multBinXi] = new TH1D(TString::Format(("corrPt_xiC[%d]"), multBinXi), TString::Format(("Eff. corrected #Xi^{#pm}: Mult: %.0f-%.0f%%"), multbins_Xi[multBinXi], multbins_Xi[multBinXi + 1]), nptbins_Xi, ptbins_Xi);
            // corrPt_xiC[multBinXi]->Sumw2();
            corrPt_xiC[multBinXi]->Divide(rawPt_xiC[multBinXi], eff_xiC_avg);
            corrPt_xiC[multBinXi]->SetMarkerStyle(markerStyles[multBinXi]);
            corrPt_xiC[multBinXi]->Write();

            corrPt_xiC[multBinXi]->SetTitle(TString::Format(("#Xi^{#pm}: Mult: %.0f-%.0f%%"), multbins_Xi[multBinXi], multbins_Xi[multBinXi + 1]));
            corrPt_xiC[multBinXi]->Scale(pow(2, (nmultbins_Xi - 1) - multBinXi));
            hs_xiC_corrected->Add(corrPt_xiC[multBinXi]);
        }

        for (Int_t multBinOm = 0; multBinOm < nmultbins_Om; multBinOm++)
        {
            out->cd("CorrectedSpectra");
            gDirectory->cd("Omm");
            corrPt_omm[multBinOm] = new TH1D(TString::Format(("corrPt_omm[%d]"), multBinOm), TString::Format(("Eff. corrected #Omega^{-}: Mult: %.0f-%.0f%%"), multbins_Om[multBinOm], multbins_Om[multBinOm + 1]), nptbins_Om, ptbins_Om);
            // corrPt_omm[multBinOm]->Sumw2();
            corrPt_omm[multBinOm]->Divide(rawPt_omm[multBinOm], eff_omC_avg);
            corrPt_omm[multBinOm]->SetMarkerStyle(markerStyles[multBinOm]);
            corrPt_omm[multBinOm]->Write();

            corrPt_omm[multBinOm]->SetTitle(TString::Format(("#Omega^{-}: Mult: %.0f-%.0f%%"), multbins_Om[multBinOm], multbins_Om[multBinOm + 1]));
            corrPt_omm[multBinOm]->Scale(pow(2, (nmultbins_Om - 1) - multBinOm));
            hs_omm_corrected->Add(corrPt_omm[multBinOm]);

            out->cd("CorrectedSpectra");
            gDirectory->cd("Omp");
            corrPt_omp[multBinOm] = new TH1D(TString::Format(("corrPt_omp[%d]"), multBinOm), TString::Format(("Eff. corrected #Omega^{+}: Mult: %.0f-%.0f%%"), multbins_Om[multBinOm], multbins_Om[multBinOm + 1]), nptbins_Om, ptbins_Om);
            // corrPt_omp[multBinOm]->Sumw2();
            corrPt_omp[multBinOm]->Divide(rawPt_omp[multBinOm], eff_omC_avg);
            corrPt_omp[multBinOm]->SetMarkerStyle(markerStyles[multBinOm]);
            corrPt_omp[multBinOm]->Write();

            corrPt_omp[multBinOm]->SetTitle(TString::Format(("#Omega^{+}: Mult: %.0f-%.0f%%"), multbins_Om[multBinOm], multbins_Om[multBinOm + 1]));
            corrPt_omp[multBinOm]->Scale(pow(2, (nmultbins_Om - 1) - multBinOm));
            hs_omp_corrected->Add(corrPt_omp[multBinOm]);

            out->cd("CorrectedSpectra");
            gDirectory->cd("OmC");
            corrPt_omC[multBinOm] = new TH1D(TString::Format(("corrPt_omC[%d]"), multBinOm), TString::Format(("Eff. corrected #Omega^{#pm}: Mult: %.0f-%.0f%%"), multbins_Om[multBinOm], multbins_Om[multBinOm + 1]), nptbins_Om, ptbins_Om);
            // corrPt_omC[multBinOm]->Sumw2();
            corrPt_omC[multBinOm]->Divide(rawPt_omC[multBinOm], eff_omC_avg);
            corrPt_omC[multBinOm]->SetMarkerStyle(markerStyles[multBinOm]);
            corrPt_omC[multBinOm]->Write();

            corrPt_omC[multBinOm]->SetTitle(TString::Format(("#Omega^{#pm}: Mult: %.0f-%.0f%%"), multbins_Om[multBinOm], multbins_Om[multBinOm + 1]));
            corrPt_omC[multBinOm]->Scale(pow(2, (nmultbins_Om - 1) - multBinOm));
            hs_omC_corrected->Add(corrPt_omC[multBinOm]);
        }

        TCanvas *e1 = new TCanvas("e1", "e1", 1920, 1080);
        TCanvas *e2 = new TCanvas("e2", "e2", 1920, 1080);
        TCanvas *e3 = new TCanvas("e3", "e3", 1920, 1080);
        TCanvas *e4 = new TCanvas("e4", "e4", 1920, 1080);
        TCanvas *e5 = new TCanvas("e5", "e5", 1920, 1080);
        TCanvas *e6 = new TCanvas("e6", "e6", 1920, 1080);

        e1->cd();
        hs_xip_corrected->Draw("PLC PMC NOSTACK");
        gPad->SetLogy();
        gPad->BuildLegend(0.9, 0.6, 1., 1., "");
        hs_xip_corrected->GetXaxis()->SetTitle("#it{p}_{T} (GeV/c)");
        hs_xip_corrected->GetYaxis()->SetTitle("#frac{1}{#it{N}_{inel}} #frac{d#it{N}}{d#it{p}_{T}}");
        e1->Modified();
        e1->ForceUpdate();

        e2->cd();
        hs_omp_corrected->Draw("PLC PMC NOSTACK");
        gPad->SetLogy();
        gPad->BuildLegend(0.9, 0.7, 1., 0.9, "");
        hs_omp_corrected->GetXaxis()->SetTitle("#it{p}_{T} (GeV/c)");
        hs_omp_corrected->GetYaxis()->SetTitle("#frac{1}{#it{N}_{inel}} #frac{d#it{N}}{d#it{p}_{T}}");
        e2->Modified();
        e2->ForceUpdate();

        e3->cd();
        hs_xim_corrected->Draw("PLC PMC NOSTACK");
        gPad->SetLogy();
        gPad->BuildLegend(0.9, 0.6, 1., 1., "");
        hs_xim_corrected->GetXaxis()->SetTitle("#it{p}_{T} (GeV/c)");
        hs_xim_corrected->GetYaxis()->SetTitle("#frac{1}{#it{N}_{inel}} #frac{d#it{N}}{d#it{p}_{T}}");
        e3->Modified();
        e3->ForceUpdate();

        e4->cd();
        hs_omm_corrected->Draw("PLC PMC NOSTACK");
        gPad->SetLogy();
        gPad->BuildLegend(0.9, 0.7, 1., 0.9, "");
        hs_omm_corrected->GetXaxis()->SetTitle("#it{p}_{T} (GeV/c)");
        hs_omm_corrected->GetYaxis()->SetTitle("#frac{1}{#it{N}_{inel}} #frac{d#it{N}}{d#it{p}_{T}}");
        e4->Modified();
        e4->ForceUpdate();

        e5->cd();
        hs_xiC_corrected->Draw("PLC PMC NOSTACK");
        gPad->SetLogy();
        gPad->BuildLegend(0.9, 0.6, 1., 1., "");
        hs_xiC_corrected->GetXaxis()->SetTitle("#it{p}_{T} (GeV/c)");
        hs_xiC_corrected->GetYaxis()->SetTitle("#frac{1}{#it{N}_{inel}} #frac{d#it{N}}{d#it{p}_{T}}");
        e5->Modified();
        e5->ForceUpdate();

        e6->cd();
        hs_omC_corrected->Draw("PLC PMC NOSTACK");
        gPad->SetLogy();
        gPad->BuildLegend(0.9, 0.7, 1., 0.9, "");
        hs_omC_corrected->GetXaxis()->SetTitle("#it{p}_{T} (GeV/c)");
        hs_omC_corrected->GetYaxis()->SetTitle("#frac{1}{#it{N}_{inel}} #frac{d#it{N}}{d#it{p}_{T}}");
        e6->Modified();
        e6->ForceUpdate();

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

        out->cd("CorrectedSpectra");
        hs_xip_corrected->Write();
        hs_xim_corrected->Write();
        hs_omp_corrected->Write();
        hs_omm_corrected->Write();
        hs_xiC_corrected->Write();
        hs_omC_corrected->Write();

        if (saveStack)
        {
            e1->cd();
            SaveImage(Form("%s/correctedPt", outputFolder.Data()), "corr_xipN", imageFormat.Data());

            e2->cd();
            SaveImage(Form("%s/correctedPt", outputFolder.Data()), "corr_ompN", imageFormat.Data());

            e3->cd();
            SaveImage(Form("%s/correctedPt", outputFolder.Data()), "corr_ximN", imageFormat.Data());

            e4->cd();
            SaveImage(Form("%s/correctedPt", outputFolder.Data()), "corr_ommN", imageFormat.Data());

            e5->cd();
            SaveImage(Form("%s/correctedPt", outputFolder.Data()), "corr_xicN", imageFormat.Data());

            e6->cd();
            SaveImage(Form("%s/correctedPt", outputFolder.Data()), "corr_omcN", imageFormat.Data());
        }
    }

    out->Close();
    // delete out;

    // compare averages

    // TFile *outTemp = TFile::Open("/var/home/ishaan/Work/git/analysis/p-Pb/MC/rootResults/effComp/compareAvgs.root", "UPDATE");
    // outTemp->cd();
    // eff_xiC_avg->SetName("eff_xiC_avg_4");
    // eff_xiC_avg->Write();
    // eff_omC_avg->SetName("eff_omC_avg_4");
    // eff_omC_avg->Write();
    // outTemp->Close();

    return 0;
}
