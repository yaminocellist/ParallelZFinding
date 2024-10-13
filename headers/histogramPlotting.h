#include "globalDefinitions.h"
#include "TLine.h"
#include "TLegend.h"
#include "TStyle.h"

void angularPlot1D (TH1D* const histo, std::vector<std::string> method, const std::string &fileTitle) {
    std::cout << histo->GetXaxis()->GetXmin() << std::endl;
    int maxBin = histo->GetMaximumBin();
    double maxBinCenter = histo->GetBinCenter(maxBin);
    int maxEntry = histo -> GetBinContent(maxBin);
    TCanvas *can1 = new TCanvas("c1d","c1d",0,50,1920,1056);
    histo -> Draw();
    histo -> SetFillColor(kYellow - 7);
    histo -> SetLineWidth(1);
    histo -> SetFillStyle(1001);
    histo -> GetXaxis() -> SetTitleSize(.05);
    histo -> GetXaxis() -> SetLabelSize(.03);
    histo -> GetXaxis() -> CenterTitle(true);
    histo -> GetXaxis() -> SetNdivisions(31, 5, 0);
    histo -> GetXaxis() -> SetTitleOffset(.8);
    histo -> GetYaxis() -> SetTitleSize(.05);
    histo -> GetYaxis() -> SetLabelSize(.03);
    histo -> GetYaxis() -> SetTitleOffset(.8);
    histo -> GetYaxis() -> CenterTitle(true);

    if (histo->GetXaxis()->GetXmin() < -M_PI) {
        double two_pi = TMath::Pi()*2;
        int bin_min = 1;                   // The first bin
        int bin_max = histo->GetNbinsX();  // The last bin

        // Calculate bin positions for each label
        int bin_pi   = bin_max;
        int bin_0    = bin_min + (bin_max - bin_min)/2;
        int binPi_2  = bin_min + 3*(bin_max - bin_min)/4;
        int bin_pi_2 = bin_min + (bin_max - bin_min)/4;

        // Set the labels at the calculated positions
        histo->GetXaxis()->SetBinLabel(bin_0, "0");
        histo->GetXaxis()->SetBinLabel(bin_pi_2, "-#pi");
        histo->GetXaxis()->SetBinLabel(binPi_2, "#pi");
        histo->GetXaxis()->SetBinLabel(bin_pi, "2#pi");
        histo->GetXaxis()->SetBinLabel(bin_min, "-2#pi");
        histo->GetXaxis()->LabelsOption("h"); // Draw the labels vertically

        // Ensure the custom labels are displayed by setting the number of divisions
        histo->GetXaxis()->SetNdivisions(9, 0, 0, kFALSE);
        histo->GetXaxis()->SetLabelSize(0.04);
    }
    else {
        double pi = TMath::Pi();
        int bin_min = 1;                   // The first bin
        int bin_max = histo->GetNbinsX();  // The last bin

        // Calculate bin positions for each label
        int bin_pi   = bin_max;
        int bin_0    = bin_min + (bin_max - bin_min)/2;
        int binPi_2  = bin_min + 3*(bin_max - bin_min)/4;
        int bin_pi_2 = bin_min + (bin_max - bin_min)/4;

        // Set the labels at the calculated positions
        histo->GetXaxis()->SetBinLabel(bin_0, "0");
        histo->GetXaxis()->SetBinLabel(bin_pi_2, "#frac{-#pi}{2}");
        histo->GetXaxis()->SetBinLabel(binPi_2, "#frac{#pi}{2}");
        // histo->GetXaxis()->SetBinLabel(bin_3pi_4, "#frac{3#pi}{4}");
        histo->GetXaxis()->SetBinLabel(bin_pi, "#pi");
        histo->GetXaxis()->SetBinLabel(bin_min, "-#pi");
        histo->GetXaxis()->LabelsOption("h"); // Draw the labels vertically
        // histo->GetXaxis()->SetBinLabel(bin_max - 3*(bin_max - bin_min)/4, "-#frac{3#pi}{4}");
        // histo->GetXaxis()->SetBinLabel(bin_max - (bin_max - bin_min)/2, "-#frac{#pi}{2}");
        // histo->GetXaxis()->SetBinLabel(bin_max - (bin_max - bin_min)/4, "-#frac{#pi}{4}");

        // Ensure the custom labels are displayed by setting the number of divisions
        histo->GetXaxis()->SetNdivisions(9, 0, 0, kFALSE);
        histo->GetXaxis()->SetLabelSize(0.04);
    }
    

    TLine *l = new TLine(0, 0, 0, maxEntry);
	l -> Draw("same"); 
    l -> SetLineColor(kRed);
    TLegend *lg = new TLegend(0.12, 0.82, 0.46, 0.9);
    lg -> AddEntry(histo, "-250 <= MBD_z_vtx <= -50 mm, MBD_centrality <= 0.70", "f");
    gStyle -> SetLegendTextSize(.022);
    lg->Draw("same");
    // gPad -> SetLogy();
    can1 -> SaveAs(Form("../../External/zFindingPlots/%s.png", fileTitle.c_str()));
}

void angularPlot2D(TH2D *const h2D, std::vector<std::string> method, const std::string &fileTitle) {
    TCanvas *can = new TCanvas("c2d","c2d",0,50,1920,1056);
    // h2D -> SetTitle(Form("All dPhi values (no event mixed) with %d bins of Z vtx", h2D->GetNbinsY()));
    // h2D -> Draw("lego2");

    if (h2D->GetXaxis()->GetXmin() < -M_PI) {
        double two_pi = TMath::Pi()*2;
        int bin_min = 1;                   // The first bin
        int bin_max = h2D->GetNbinsX();  // The last bin

        // Calculate bin positions for each label
        int bin_pi   = bin_max;
        int bin_0    = bin_min + (bin_max - bin_min)/2;
        int binPi_2  = bin_min + 3*(bin_max - bin_min)/4;
        int bin_pi_2 = bin_min + (bin_max - bin_min)/4;

        // Set the labels at the calculated positions
        h2D->GetXaxis()->SetBinLabel(bin_0, "0");
        h2D->GetXaxis()->SetBinLabel(bin_pi_2, "-#pi");
        h2D->GetXaxis()->SetBinLabel(binPi_2, "#pi");
        h2D->GetXaxis()->SetBinLabel(bin_pi, "2#pi");
        h2D->GetXaxis()->SetBinLabel(bin_min, "-2#pi");
        h2D->GetXaxis()->LabelsOption("h"); // Draw the labels vertically

        // Ensure the custom labels are displayed by setting the number of divisions
        h2D->GetXaxis()->SetNdivisions(9, 0, 0, kFALSE);
        h2D->GetXaxis()->SetLabelSize(0.04);
    }
    else {
        double pi = TMath::Pi();
        int bin_min = 1;                   // The first bin
        int bin_max = h2D->GetNbinsX();  // The last bin

        // Calculate bin positions for each label
        int bin_pi   = bin_max;
        int bin_0    = bin_min + (bin_max - bin_min)/2;
        int binPi_2  = bin_min + 3*(bin_max - bin_min)/4;
        int bin_pi_2 = bin_min + (bin_max - bin_min)/4;

        // Set the labels at the calculated positions
        h2D->GetXaxis()->SetBinLabel(bin_0, "0");
        h2D->GetXaxis()->SetBinLabel(bin_pi_2, "#frac{-#pi}{2}");
        h2D->GetXaxis()->SetBinLabel(binPi_2, "#frac{#pi}{2}");
        // h2D->GetXaxis()->SetBinLabel(bin_3pi_4, "#frac{3#pi}{4}");
        h2D->GetXaxis()->SetBinLabel(bin_pi, "#pi");
        h2D->GetXaxis()->SetBinLabel(bin_min, "-#pi");
        h2D->GetXaxis()->LabelsOption("h"); // Draw the labels vertically
        // h2D->GetXaxis()->SetBinLabel(bin_max - 3*(bin_max - bin_min)/4, "-#frac{3#pi}{4}");
        // h2D->GetXaxis()->SetBinLabel(bin_max - (bin_max - bin_min)/2, "-#frac{#pi}{2}");
        // h2D->GetXaxis()->SetBinLabel(bin_max - (bin_max - bin_min)/4, "-#frac{#pi}{4}");

        // Ensure the custom labels are displayed by setting the number of divisions
        h2D->GetXaxis()->SetNdivisions(9, 0, 0, kFALSE);
        h2D->GetXaxis()->SetLabelSize(0.04);
    }

    h2D -> Draw("colz");
    h2D -> GetXaxis() -> CenterTitle(true);
    h2D -> GetYaxis() -> CenterTitle(true);
    h2D -> GetXaxis() -> SetTitleOffset(1.4);   h2D -> GetYaxis() -> SetTitleOffset(1.4);
    can -> Modified();
    can -> Update();
    can -> SaveAs(Form("../../External/zFindingPlots/%s.png", fileTitle.c_str()));
}

void angularPlot3D(TH2D * h_dPhi_Z, std::vector<std::string> method, const std::string &fileTitle) {
    TCanvas *can1 = new TCanvas("c3d","c3d",0,50,1920,1056);
    h_dPhi_Z -> SetTitle(Form("All dPhi values (no event mixed) with %d bins of Z vtx", h_dPhi_Z->GetNbinsY()));
    h_dPhi_Z -> Draw("lego2");
    h_dPhi_Z -> GetXaxis() -> CenterTitle(true);
    h_dPhi_Z -> GetYaxis() -> CenterTitle(true);
    h_dPhi_Z -> GetXaxis() -> SetTitleOffset(1.8);   h_dPhi_Z -> GetYaxis() -> SetTitleOffset(1.8);
    can1 -> Modified();
    can1 -> Update();

    can1 -> SaveAs(Form("../External/zFindingPlots/%s.png", fileTitle.c_str()));
}

void ZResolutionSinglePlot (TH1D* const histo, std::vector<std::string> method, const std::string &fileTitle) {
    if (isInteger(method[1]))   int lowerRange = stoi(method[1]);
    if (isInteger(method[2]))   int upperRange = stoi(method[2]);
    int maxBin = histo->GetMaximumBin();
    double maxBinCenter = histo->GetBinCenter(maxBin);
    int maxEntry = histo -> GetBinContent(maxBin);
    TCanvas *can1 = new TCanvas("c1d","c1d",0,50,1920,1056);
    histo -> Draw();
    histo -> SetFillColor(kYellow - 7);
    histo -> SetLineWidth(1);
    histo -> SetFillStyle(1001);
    histo -> GetXaxis() -> SetTitleSize(.05);
    histo -> GetXaxis() -> SetLabelSize(.03);
    histo -> GetXaxis() -> CenterTitle(true);
    histo -> GetXaxis() -> SetNdivisions(31, 5, 0);
    histo -> GetXaxis() -> SetTitleOffset(.8);
    histo -> GetYaxis() -> SetTitleSize(.05);
    histo -> GetYaxis() -> SetLabelSize(.03);
    histo -> GetYaxis() -> SetTitleOffset(.8);
    histo -> GetYaxis() -> CenterTitle(true);

    TLine *l = new TLine(0, 0, 0, maxEntry);
	l -> Draw("same"); 
    l -> SetLineColor(kRed);
    TLegend *lg = new TLegend(0.12, 0.8, 0.42, 0.9);
    lg -> AddEntry(histo, "Fiducial cut as -250 < MBD_z_vtx < -50 mm", "f");
    gStyle -> SetLegendTextSize(.02);
    lg->Draw("same");
    // gPad -> SetLogy();
    can1 -> SaveAs(Form("../External/zFindingPlots/%s.png", fileTitle.c_str()));
}

void TGraphSinglePlot (TGraph* g0, const char *title, const char *Xtitle, const std::string &fileTitle) {
    TCanvas *can1 = new TCanvas("cg","cg",0,50,2100,1200);
    g0 -> SetMarkerStyle(29);
    g0 -> SetMarkerSize(1.1);
    g0 -> SetMarkerColor(kBlue - 7);
    g0 -> SetLineWidth(3);
    g0 -> SetLineColor(kWhite);
    gStyle -> SetTitleW(0.7);  //per cent of the pad width
    gStyle -> SetTitleH(0.08); //per cent of the pad height
    g0 -> SetTitle(title);
    // g0 -> GetXaxis() -> SetTitle("# of hits");
    g0 -> GetXaxis() -> SetTitle(Xtitle);
    g0 -> GetXaxis() -> SetTitleSize(0.05);
    g0 -> GetXaxis() -> SetLabelSize(0.04);
    g0 -> GetXaxis() -> CenterTitle(true);
    g0 -> GetYaxis() -> SetTitle("found z resolution [mm]");
    g0 -> GetYaxis() -> SetTitleSize(0.05);
    g0 -> GetYaxis() -> SetLabelSize(0.025);
    g0 -> GetYaxis() -> CenterTitle(true);
    // g0 -> SetMinimum(-200); // Setting y range;
    // g0 -> SetMaximum(200);  // Setting y range;
    g0 -> GetYaxis() -> SetTitleOffset(0.8); 
    g0 -> GetXaxis() -> SetTitleOffset(0.8); 
    // g0 -> GetXaxis() -> SetLimits(0, 6000); // Setting x range;
    g0 -> Draw("AP SAME");
    gPad->SetGrid(5, 2); gPad->Update();

    TLine *l = new TLine(-250, 0, 8000, 0);
	l -> Draw("same"); 
    l -> SetLineColor(kRed);
    TLegend *lg = new TLegend(0.55, 0.85, 0.9, 0.9);
    lg -> AddEntry(g0, "Fiducial cut as -250 < MBD_z_vtx < -50 mm, centrality <= 0.7", "f");
    gStyle -> SetLegendTextSize(.02);
    lg->Draw("same");
    can1 -> SaveAs(Form("../External/zFindingPlots/%s.png", fileTitle.c_str()));
}

void TGraphSinglePlot_Squared (TGraph* g0, const char *title, const char *Xtitle, const char *Ytitle, const std::string &fileTitle) {
    TCanvas *can1 = new TCanvas("cgs","cgs",0,50,1200,1200);
    g0 -> SetMarkerStyle(29);
    g0 -> SetMarkerSize(1.1);
    g0 -> SetMarkerColor(kBlue - 7);
    g0 -> SetLineWidth(3);
    g0 -> SetLineColor(kWhite);
    gStyle -> SetTitleW(0.7);  //per cent of the pad width
    gStyle -> SetTitleH(0.08); //per cent of the pad height
    g0 -> SetTitle(title);
    g0 -> GetXaxis() -> SetTitle(Xtitle);   g0 -> GetYaxis() -> SetTitle(Ytitle);
    g0 -> GetXaxis() -> SetTitleSize(0.05);
    g0 -> GetXaxis() -> SetLabelSize(0.04);
    g0 -> GetXaxis() -> CenterTitle(true);
    g0 -> GetYaxis() -> SetTitleSize(0.05);
    g0 -> GetYaxis() -> SetLabelSize(0.025);
    g0 -> GetYaxis() -> CenterTitle(true);
    // g0 -> SetMinimum(-200); // Setting y range;
    // g0 -> SetMaximum(200);  // Setting y range;
    g0 -> GetYaxis() -> SetTitleOffset(0.8); 
    g0 -> GetXaxis() -> SetTitleOffset(0.8); 
    // g0 -> GetXaxis() -> SetLimits(0, 6000); // Setting x range;
    g0 -> Draw("AP SAME");
    gPad->SetGrid(5, 2); gPad->Update();

    TLine *l = new TLine(-250, -250, -50, -50);
	l -> Draw("same"); 
    l -> SetLineColor(kRed);
    TLegend *lg = new TLegend(0.36, 0.85, 0.9, 0.9);
    lg -> AddEntry(g0, "Fiducial cut as -250 < MBD_z_vtx < -50 mm, centrality <= 0.7", "p");
    gStyle -> SetLegendTextSize(.02);
    lg->Draw("same");
    can1 -> SaveAs(Form("../External/zFindingPlots/%s.png", fileTitle.c_str()));
}

void TGraphMultiPlot (TMultiGraph *mg, const char *title, const char *Xtitle, const char *Ytitle, const std::string &fileTitle) {
    TList *graphs = mg->GetListOfGraphs();
    int nGraphs = graphs->GetSize();    // number of graphs;
    TIter next(graphs);
    TGraph *graph;
    for (int i = 0; i < nGraphs; i++) {
        graph = (TGraph *)next();
        graph -> SetMarkerStyle(29);
        graph -> SetMarkerSize(1.1);
        graph -> SetMarkerColor(kBlue - 7);
        if (i == 1)     graph -> SetMarkerColor(2);
        graph -> SetLineWidth(3);
        graph -> SetLineColor(kWhite);
        graph -> GetXaxis() -> SetTitleSize(0.05);
        graph -> GetXaxis() -> SetLabelSize(0.04);
        graph -> GetXaxis() -> CenterTitle(true);
        graph -> GetYaxis() -> SetTitleSize(0.05);
        graph -> GetYaxis() -> SetLabelSize(0.025);
        graph -> GetYaxis() -> CenterTitle(true);
        // graph -> SetMinimum(-200); // Setting y range;
        // graph -> SetMaximum(200);  // Setting y range;
        graph -> GetYaxis() -> SetTitleOffset(0.8); 
        graph -> GetXaxis() -> SetTitleOffset(0.8); 
    }

    // Use "xdpyinfo | grep dimensions" in terminal to get your display's max size;
    TCanvas *can1 = new TCanvas("cgm","cgm",0,50,1056,1056);
    mg->SetTitle(title);
    mg->GetXaxis()->SetTitle(Xtitle);
    mg->GetYaxis()->SetTitle(Ytitle);
    mg -> Draw("AP SAME");
    gStyle -> SetTitleW(0.7);  //per cent of the pad width
    gStyle -> SetTitleH(0.08); //per cent of the pad height
    gPad->SetGrid(5, 2); gPad->Update();

    // TLine *l = new TLine(-250, -250, -50, -50);
	// l -> Draw("same"); 
    // l -> SetLineColor(kRed);
        // TLegend *lg = new TLegend(0.36, 0.85, 0.9, 0.9);
        // lg -> AddEntry(graphs[0], "Fiducial cut as -250 < MBD_z_vtx < -50 mm, centrality <= 0.7", "p");
        // gStyle -> SetLegendTextSize(.02);
        // lg->Draw("same");
    can1 -> SaveAs(Form("../External/zFindingPlots/%s.png", fileTitle.c_str()));
}

void EtaPhiSinglePlot (TH1D* const histo, std::vector<std::string> method, Int_t const & target) {
    if (isInteger(method[1]))   int lowerRange = stoi(method[1]);
    if (isInteger(method[2]))   int upperRange = stoi(method[2]);
    int maxBin = histo->GetMaximumBin();
    double maxBinCenter = histo->GetBinCenter(maxBin);
    TCanvas *can1 = new TCanvas("c1","c1",0,50,1920,1056);
    histo -> Draw();
    histo -> SetFillColor(kYellow - 7);
    histo -> SetLineWidth(1);
    histo -> SetFillStyle(1001);
    histo -> GetXaxis() -> SetTitleSize(.05);
    histo -> GetXaxis() -> SetLabelSize(.03);
    histo -> GetXaxis() -> CenterTitle(true);
    histo -> GetXaxis() -> SetNdivisions(31, 5, 0);
    histo -> GetXaxis() -> SetTitleOffset(.8);
    histo -> GetYaxis() -> SetTitleSize(.05);
    histo -> GetYaxis() -> SetLabelSize(.03);
    histo -> GetYaxis() -> SetTitleOffset(.8);
    histo -> GetYaxis() -> CenterTitle(true);

    double pi = TMath::Pi();
    int bin_min = 1;  // The first bin
    int bin_max = histo->GetNbinsX();  // The last bin

    // Calculate bin positions for each label
    int bin_pi = bin_max;
    int bin_0 = bin_min + (bin_max - bin_min)/2;
    int binPi_2 = bin_min + 3*(bin_max - bin_min)/4;
    // int bin_pi_2 = bin_min + (bin_max - bin_min)/2;
    int bin_pi_2 = bin_min + (bin_max - bin_min)/4;

    // Set the labels at the calculated positions
    histo->GetXaxis()->SetBinLabel(bin_0, "0");
    histo->GetXaxis()->SetBinLabel(bin_pi_2, "#frac{-#pi}{2}");
    histo->GetXaxis()->SetBinLabel(binPi_2, "#frac{#pi}{2}");
    // histo->GetXaxis()->SetBinLabel(bin_3pi_4, "#frac{3#pi}{4}");
    histo->GetXaxis()->SetBinLabel(bin_pi, "#pi");
    histo->GetXaxis()->SetBinLabel(bin_min, "-#pi");
    // histo->GetXaxis()->SetBinLabel(bin_max - 3*(bin_max - bin_min)/4, "-#frac{3#pi}{4}");
    // histo->GetXaxis()->SetBinLabel(bin_max - (bin_max - bin_min)/2, "-#frac{#pi}{2}");
    // histo->GetXaxis()->SetBinLabel(bin_max - (bin_max - bin_min)/4, "-#frac{#pi}{4}");

    // Ensure the custom labels are displayed by setting the number of divisions
    histo->GetXaxis()->SetNdivisions(9, 0, 0, kFALSE);
    histo->GetXaxis()->SetLabelSize(0.04);
    // Update histogram to refresh the axis
    histo->Draw("HIST");
    histo->GetXaxis()->LabelsOption("h"); // Draw the labels vertically
    // histo -> GetYaxis() -> SetRangeUser(5300e3, 6130e3);
    // histo -> SetTitle(Form("dPhi data of all events whose found z vtx is ~ [-%d, -%d] mm, centered at %0.4f", lowerRange, upperRange, maxBinCenter));
    // gPad -> SetLogy();
    // can1 -> SaveAs(Form("../External/xyFindingPlots/dPhi_all_%d_%d.png", lowerRange, upperRange));
    // histo-> SetTitle(Form("dPhi of %d events mixed up in range of [-209.375, -207.5] mm", target));
    // can1 -> SaveAs(Form("../External/xyFindingPlots/dPhi_mixed_%d.png", target));
}

void doublePlot (TH1D* const hBackground, TH1D* const hSignal, std::vector<std::string> method, Int_t const & target) {
    int upperRange = stoi(method[2]);
    int lowerRange = stoi(method[1]);
    int maxBin = hBackground->GetMaximumBin();
    double maxBinCenter = hBackground->GetBinCenter(maxBin);
    TCanvas *can1 = new TCanvas("c1","c1",0,50,1920,1056);

    // Draw the first histogram
    hBackground -> Draw();
    hBackground -> SetFillColor(kYellow - 7);
    hBackground -> SetLineWidth(1);
    hBackground -> SetFillStyle(1001);
    hBackground -> GetXaxis() -> SetTitle("Phi value");
    hBackground -> GetXaxis() -> SetTitleSize(.05);
    hBackground -> GetXaxis() -> SetLabelSize(.03);
    hBackground -> GetXaxis() -> CenterTitle(true);
    hBackground -> GetXaxis() -> SetNdivisions(31, 5, 0);
    hBackground -> GetXaxis() -> SetTitleOffset(.8);
    hBackground -> GetYaxis() -> SetTitle("# of Counts");
    hBackground -> GetYaxis() -> SetTitleSize(.05);
    hBackground -> GetYaxis() -> SetLabelSize(.03);
    hBackground -> GetYaxis() -> SetTitleOffset(.8);
    hBackground -> GetYaxis() -> CenterTitle(true);
    hBackground -> GetYaxis() -> SetRangeUser(4000e3, 6200e3); // Setting x range;

    // Draw the second histogram on the same canvas
    hSignal -> Draw("SAME");
    hSignal -> SetFillColor(kBlue - 7);
    hSignal -> SetLineWidth(1);
    hSignal -> SetFillStyle(1001);

    // hBackground -> SetTitle(Form("dPhi of %d events mix-up/non-mix-up in range of [-209.375, -207.5] mm", target));
    // hBackground-> SetTitle(Form("dPhi of %d events mixed up in range of [-209.375, -207.5] mm", target));
    // can1 -> SaveAs(Form("../External/xyFindingPlots/dPhi_mixed_%d.png", target));
}

void backgroundCancelling (TH1D* const hBackground, TH1D* const hSignal, std::vector<std::string> method, Int_t const & target) {
    TCanvas *can1 = new TCanvas("csub","csub",0,50,1920,1056);
    can1 -> Divide(1, 2);
    can1 -> cd(1);
    hBackground -> Draw();
    hBackground -> SetFillColor(kYellow - 7);
    hBackground -> SetLineWidth(1);
    hBackground -> SetFillStyle(1001);
    hBackground -> GetXaxis() -> SetTitle("dPhi");
    hBackground -> GetXaxis() -> SetTitleSize(.05);
    hBackground -> GetXaxis() -> SetLabelSize(.03);
    hBackground -> GetXaxis() -> CenterTitle(true);
    hBackground -> GetXaxis() -> SetNdivisions(31, 5, 0);
    hBackground -> GetXaxis() -> SetTitleOffset(.8);
    hBackground -> GetYaxis() -> SetTitle("# of Counts");
    hBackground -> GetYaxis() -> SetTitleSize(.05);
    hBackground -> GetYaxis() -> SetLabelSize(.03);
    hBackground -> GetYaxis() -> SetTitleOffset(.8);
    hBackground -> GetYaxis() -> CenterTitle(true);
    hBackground -> GetYaxis() -> SetRangeUser(4300e3, 6500e3); // Setting x range;

    // Draw the second histogram on the same canvas
    hSignal -> Draw("SAME");
    hSignal -> SetFillColor(kBlue - 7);
    hSignal -> SetLineWidth(1);
    hSignal -> SetFillStyle(1001);

    hBackground -> SetTitle(Form("dPhi of %d events mix-up/non-mix-up in range of [-209.375, -207.5] mm", target));

    can1 -> cd(2);
    double centralPeak = 0.05;
    double N1 = hBackground -> Integral(hBackground->FindFixBin(-M_PI), hBackground->FindFixBin(-centralPeak), "") + 
                hBackground -> Integral(hBackground->FindFixBin(centralPeak), hBackground->FindFixBin(M_PI), "");
    double N2 = hSignal -> Integral(hSignal->FindFixBin(-M_PI), hSignal->FindFixBin(-centralPeak), "") + 
                hSignal -> Integral(hSignal->FindFixBin(centralPeak), hSignal->FindFixBin(M_PI), "");
    double N  = N2/N1;

    TH1D* hNormalized = (TH1D*) hBackground->Clone("Normalized hBackground");
    hNormalized -> Add(hBackground, N-1);
    TH1D* hDiff = (TH1D*) hSignal->Clone("Background Subtracted Signal");
    hDiff -> Add(hNormalized, -1);
    hDiff -> Sumw2();
    hDiff -> Draw("HIST SAME");
    hDiff -> Draw("e1psame");
    hDiff -> GetXaxis() -> SetTitle("dPhi");
    hDiff -> GetXaxis() -> CenterTitle(true);
    hDiff -> GetXaxis() -> SetTitleSize(.05);
    hDiff -> SetFillColor(kYellow - 7);
    hDiff -> SetFillStyle(1001);
    //hDiff -> GetXaxis() -> SetRange(hDiff -> FindFixBin(-0.5),hDiff -> FindFixBin(0.5));
    double peak  = hDiff -> Integral(hDiff->FindFixBin(-centralPeak), hDiff->FindFixBin(centralPeak), "");
    // double Ratio = peak/events;
    
    // hDiff -> SetTitle(Form("Subtracted Signal for %5.0f Events", events));
    hDiff -> SetTitle("Background Subtracted dPhi");
    gStyle->SetEndErrorSize(6);
    gStyle->SetErrorX(0.5);
    
    gStyle -> SetTitleFont(100,"t");
    gStyle -> SetTitleSize(0.065,"t");

    //gPad -> SetLogy(1);
    TLine *l = new TLine(-2,0,2,0);
	l -> Draw("same"); 
    l -> SetLineColor(kGreen);

    can1 -> SaveAs("../../External/zFindingPlots/dPhi_subtraction.png");
}

void backgroundCancelling_dPhi (TH1D* const hBackground, TH1D* const hSignal, std::vector<std::string> method, Int_t const & target) {
    // TCanvas *can1 = new TCanvas("csub","csub",0,50,1920,1056);
    TCanvas *can1 = new TCanvas("csub","csub",0,50,2560,1440);
    can1 -> Divide(1, 2);
    can1 -> cd(1);
    double phi_range_low = -2.4, phi_range_high = -1.8;
    int bin_range_low = hBackground->FindBin(phi_range_low), bin_range_high = hBackground->FindBin(phi_range_high);
    double max_unmixed = -1, max_mixed = -1, current_binContent;
    for (int bin = bin_range_low; bin <= bin_range_high; bin++) {
        current_binContent = hSignal->GetBinContent(bin);
        if (max_unmixed < current_binContent)  max_unmixed = current_binContent;
        current_binContent = hBackground->GetBinContent(bin);
        if (max_mixed < current_binContent)  max_mixed = current_binContent;
    }
    hSignal -> Add(hSignal, max_mixed/max_unmixed - 1);
    hSignal -> Draw("SAME");
    hSignal -> GetXaxis()->SetLabelSize(0.06);
    hSignal -> SetTitle(Form("%d events, -%2.2fcm < z vtx < -%2.2fcm, %1.2f < centrality < %1.2f", target, std::stod(method[4]), std::stod(method[5]), std::stod(method[2]), std::stod(method[3])));

    max_mixed   = hBackground->GetBinContent(hBackground->GetMaximumBin());
    max_unmixed = hSignal->GetBinContent(hSignal->GetMaximumBin());

    if (max_unmixed > max_mixed) {
        hSignal -> GetYaxis() -> SetRangeUser(max_unmixed*0.8, max_unmixed*1.1);
        hBackground -> GetYaxis() -> SetRangeUser(max_unmixed*0.8, max_unmixed*1.1);
    }   
    else {
        hSignal -> GetYaxis() -> SetRangeUser(max_mixed*0.8, max_mixed*1.1);
        hBackground -> GetYaxis() -> SetRangeUser(max_mixed*0.8, max_mixed*1.1);
    }                      
    hBackground -> Draw("SAME");
    hBackground -> SetLineColor(2);

    TLegend *lg = new TLegend(0.12, 0.8, 0.33, 0.9);
    lg -> AddEntry(hSignal, "Unmixed Events' dPhi", "l");
    lg -> AddEntry(hBackground, "Mixed Events' dPhi", "l");
    gStyle -> SetLegendTextSize(.043);
    lg->Draw("same");
    can1 -> cd(2);
    double centralPeak = 0.05;
    double N1 = hBackground -> Integral(hBackground->FindFixBin(-M_PI), hBackground->FindFixBin(-centralPeak), "") + 
                hBackground -> Integral(hBackground->FindFixBin(centralPeak), hBackground->FindFixBin(M_PI), "");
    double N2 = hSignal -> Integral(hSignal->FindFixBin(-M_PI), hSignal->FindFixBin(-centralPeak), "") + 
                hSignal -> Integral(hSignal->FindFixBin(centralPeak), hSignal->FindFixBin(M_PI), "");
    double N  = N2/N1;

    TH1D* hNormalized = (TH1D*) hBackground->Clone("Normalized hBackground");
    hNormalized -> Add(hBackground, N-1);
    TH1D* hDiff = (TH1D*) hSignal->Clone("Background Subtracted Signal");
    hDiff -> Add(hNormalized, -1);
    hDiff -> Sumw2();
    hDiff -> Draw("HIST SAME");
    hDiff -> Draw("e1psame");
    hDiff -> GetXaxis() -> SetTitle("dPhi");
    hDiff -> GetXaxis() -> CenterTitle(true);
    hDiff -> GetXaxis() -> SetTitleSize(.05);
    hDiff->GetXaxis()->SetLabelSize(0.06);
    hDiff -> SetFillColor(kYellow - 7);
    hDiff -> SetFillStyle(1001);
    //hDiff -> GetXaxis() -> SetRange(hDiff -> FindFixBin(-0.5),hDiff -> FindFixBin(0.5));
    double peak  = hDiff -> Integral(hDiff->FindFixBin(-centralPeak), hDiff->FindFixBin(centralPeak), "");
    // double Ratio = peak/events;
    
    // hDiff -> SetTitle(Form("Subtracted Signal for %5.0f Events", events));
    hDiff -> SetTitle("Background Subtracted dPhi");
    gStyle->SetEndErrorSize(6);
    gStyle->SetErrorX(0.5);
    
    gStyle -> SetTitleFont(100,"t");
    gStyle -> SetTitleSize(0.065,"t");

    //gPad -> SetLogy(1);
    TLine *l = new TLine(-2,0,2,0);
	l -> Draw("same"); 
    l -> SetLineColor(kGreen);

    can1 -> SaveAs(Form("../../External/zFindingPlots/dPhi_mixedsubtract_%2.2f_%2.2f_%1.2f_%1.2f.png", std::stod(method[4]), std::stod(method[5]), std::stod(method[2]), std::stod(method[3])));
}

void ArrayPlot1D_Logy (const std::vector<TH1D*>& h, std::vector<std::string> method, const std::string &fileTitle) {
    double entries = 0;
    for (int i = 0; i < h.size(); i++) {
        entries += h[i] -> GetEntries();
    }
    entries /= 14;
    // h[0] -> Scale(1.0 / entries);
    // h[0]->SetFillStyle(3025);

    TCanvas *can1 = new TCanvas("c1d","c1d",0,50,2100,1200);
    h[0] -> SetLineWidth(3);
    h[0] -> SetLineColor(2);
    h[0] -> Draw();
    for (int d = 1; d < h.size(); d++) {
        // entries = h[d] -> GetEntries();
        // h[d] -> Scale(1.0 / entries);
        // h[d]->SetLineColor(d + 1); // ROOT colors start from 1 (0 is white)
        // h[d]->SetFillColor(d + 1);
        // h[d]->SetFillStyle(3001 + d); 
        h[d] -> SetLineWidth(3);
        h[d] -> Draw("SAME");
    }
    h[0] -> GetYaxis() -> SetRangeUser(5e4, h[1]->GetEntries()/2e2); // Setting x range;
    h[0] -> GetXaxis() -> CenterTitle(true);    h[0] -> GetYaxis() -> CenterTitle(true);

    double pi = TMath::Pi();
    int bin_min  = 1;  // The first bin
    int bin_max  = h[0]->GetNbinsX();  // The last bin
    // Calculate bin positions for each label
    int bin_pi   = bin_max;
    int bin_0    = bin_min + (bin_max - bin_min)/2;
    int binPi_2  = bin_min + 3*(bin_max - bin_min)/4;
    int bin_pi_2 = bin_min + (bin_max - bin_min)/4;
    // Set the labels at the calculated positions
    h[0]->GetXaxis()->SetBinLabel(bin_0, "0");
    h[0]->GetXaxis()->SetBinLabel(bin_pi_2, "#frac{-#pi}{2}");
    h[0]->GetXaxis()->SetBinLabel(binPi_2, "#frac{#pi}{2}");
    h[0]->GetXaxis()->SetBinLabel(bin_pi, "#pi");
    h[0]->GetXaxis()->SetBinLabel(bin_min, "-#pi");
    // Ensure the custom labels are displayed by setting the number of divisions
    h[0]->GetXaxis()->SetNdivisions(9, 0, 0, kFALSE);
    h[0]->GetXaxis()->SetLabelSize(0.04);
    // Update histogram to refresh the axis
    // h[0]->Draw("HIST");
    h[0]->GetXaxis()->LabelsOption("h"); // Draw the labels vertically

    gPad -> SetLogy(1);
    can1 -> SaveAs(Form("../External/zFindingPlots/%s.png", fileTitle.c_str()));
}

void ArrayPlot1D_Logy_ver2 (const std::vector<TH1D*>& h, std::vector<std::string> method, const std::string &fileTitle) {
    TCanvas *can1 = new TCanvas("c1d","c1d",0,50,2100,1200);
    // h[0] -> Draw();
    // h[1] -> Draw("SAME");
    // h[2] -> Draw("SAME");
    // h[3] -> Draw("SAME");
    // h[4] -> Draw("SAME");
    // h[5] -> Draw("SAME");
    // h[6] -> Draw("SAME");
    // h[7] -> Draw("SAME");
    // h[8] -> Draw("SAME");
    // h[9] -> Draw("SAME");
    // h[10] -> Draw("SAME");
    // h[11] -> Draw("SAME");
    // h[12] -> Draw("SAME");
    // h[13] -> Draw("SAME");
    // h[14] -> Draw("SAME");
    // h[15] -> Draw("SAME");
    // h[16] -> Draw("SAME");
    // h[17] -> Draw("SAME");
    // h[18] -> Draw("SAME");
    // h[19] -> Draw("SAME");

    double pi = TMath::Pi();
    int bin_min  = 1;  // The first bin
    int bin_max  = h[16]->GetNbinsX();  // The last bin
    // Calculate bin positions for each label
    int bin_pi   = bin_max;
    int bin_0    = bin_min + (bin_max - bin_min)/2;
    int binPi_2  = bin_min + 3*(bin_max - bin_min)/4;
    int bin_pi_2 = bin_min + (bin_max - bin_min)/4;
    // Set the labels at the calculated positions
    h[16]->GetXaxis()->SetBinLabel(bin_0, "0");
    h[16]->GetXaxis()->SetBinLabel(bin_pi_2, "#frac{-#pi}{2}");
    h[16]->GetXaxis()->SetBinLabel(binPi_2, "#frac{#pi}{2}");
    h[16]->GetXaxis()->SetBinLabel(bin_pi, "#pi");
    h[16]->GetXaxis()->SetBinLabel(bin_min, "-#pi");
    // Ensure the custom labels are displayed by setting the number of divisions
    h[16]->GetXaxis()->SetNdivisions(9, 0, 0, kFALSE);
    h[16]->GetXaxis()->SetLabelSize(0.04);
    // Update histogram to refresh the axis
    h[16]->Draw("HIST");
    h[16]->GetXaxis()->LabelsOption("h"); // Draw the labels vertically

    for (int i = 0; i < h.size(); i++) {
        // std::cout << h[i] -> GetEntries() << std::endl;
        h[i] -> Draw("SAME");
        h[i] -> SetLineWidth(3);
        h[i] -> GetYaxis() -> SetRangeUser(1e5, h[16]->GetEntries()/5e2);
        h[i] -> GetXaxis() -> CenterTitle(true);    h[i] -> GetYaxis() -> CenterTitle(true);
    }
    gPad -> SetLogy(1);
}

void ArrayPlot1D_Rescale (const std::vector<TH1D*>& h, std::vector<std::string> method, const std::string &fileTitle) {
    std::vector<double> max_entries;
    double min_y = std::numeric_limits<double>::max();
    for (int i = 0; i < h.size(); i++) {
        max_entries.push_back(h[i]->GetBinContent(h[i]->GetMaximumBin()));
        h[i] -> Add(h[i], 1/(max_entries[i]+1) - 1);

        // Find minimum bin content across all histograms
        for (int bin = 1; bin <= h[i]->GetNbinsX(); ++bin) {
            double bin_content = h[i]->GetBinContent(bin);
            if (bin_content > 0 && bin_content < min_y) {
                min_y = bin_content;
            }
        }
        // h[i] -> GetYaxis() -> SetRangeUser(0.6, 1.1);
    }
    TLegend *lg = new TLegend(0.12, 0.73, 0.25, 0.9);
    gStyle -> SetLegendTextSize(.015);
    for (int i = 0; i < h.size(); i++) {
        h[i] -> GetYaxis() -> SetRangeUser(min_y*0.9, 1.1);
        h[i] -> SetLineColor(32 + i);
        lg -> AddEntry(h[i], Form("Centrality between %.2f and %.2f", static_cast<double>(i)*0.05, static_cast<double>(i + 1)*0.05), "l");
        h[i] -> SetTitle(fileTitle.c_str());
    }
    TCanvas *can1 = new TCanvas("c1d","c1d",0,50,2560,1600);
    h[0] -> SetLineWidth(2);
    h[0] -> SetLineColor(2);
    h[0] -> Draw();
    for (int d = 1; d < h.size(); d++) {
        // entries = h[d] -> GetEntries();
        // h[d] -> Scale(1.0 / entries);
        // h[d]->SetLineColor(d + 1); // ROOT colors start from 1 (0 is white)
        // h[d]->SetFillColor(d + 1);
        // h[d]->SetFillStyle(3001 + d); 
        h[d] -> SetLineWidth(2);
        h[d] -> Draw("SAME");
    }
    // h[0] -> GetYaxis() -> SetRangeUser(5e4, h[1]->GetEntries()/2e2); // Setting x range;
    h[0] -> GetXaxis() -> CenterTitle(true);    h[0] -> GetYaxis() -> CenterTitle(true);

    double pi = TMath::Pi();
    int bin_min  = 1;  // The first bin
    int bin_max  = h[0]->GetNbinsX();  // The last bin
    // Calculate bin positions for each label
    int bin_pi   = bin_max;
    int bin_0    = bin_min + (bin_max - bin_min)/2;
    int binPi_2  = bin_min + 3*(bin_max - bin_min)/4;
    int bin_pi_2 = bin_min + (bin_max - bin_min)/4;
    // Set the labels at the calculated positions
    h[0]->GetXaxis()->SetBinLabel(bin_0, "0");
    h[0]->GetXaxis()->SetBinLabel(bin_pi_2, "#frac{-#pi}{2}");
    h[0]->GetXaxis()->SetBinLabel(binPi_2, "#frac{#pi}{2}");
    h[0]->GetXaxis()->SetBinLabel(bin_pi, "#pi");
    h[0]->GetXaxis()->SetBinLabel(bin_min, "-#pi");
    // Ensure the custom labels are displayed by setting the number of divisions
    h[0]->GetXaxis()->SetNdivisions(9, 0, 0, kFALSE);
    h[0]->GetXaxis()->SetLabelSize(0.04);
    // Update histogram to refresh the axis
    // h[0]->Draw("HIST");
    h[0]->GetXaxis()->LabelsOption("h"); // Draw the labels vertically

    lg->Draw("same");
    can1 -> Update();

    can1 -> SaveAs(Form("../../External/zFindingPlots/%s.png", fileTitle.c_str()));
}

void ArrayPlot1D_Rescale_ver2 (const std::vector<TH1D*>& h, std::vector<std::string> method, const std::string &fileTitle) {
    TCanvas *can1 = new TCanvas("c1d","c1d",0,50,2100,1200);
    std::vector<double> max_entries;
    double min_y = std::numeric_limits<double>::max();
    for (int i = 0; i < h.size(); i++) {
        max_entries.push_back(h[i]->GetBinContent(h[i]->GetMaximumBin()));
        h[i] -> Add(h[i], 1/(max_entries[i]+1) - 1);
        for (int bin = 1; bin <= h[i]->GetNbinsX(); ++bin) {
            double bin_content = h[i]->GetBinContent(bin);
            if (bin_content > 0 && bin_content < min_y) {
                min_y = bin_content;
            }
        }
    }
    h[16] -> Draw();
    TLegend *lg = new TLegend(0.12, 0.73, 0.28, 0.9);
    TLegend *lg2 = new TLegend(0.28 , 0.73, 0.44, 0.9);
    gStyle -> SetLegendTextSize(.015);
    for (int i = 0; i < h.size(); i++) {
        if (i < 10) {
            lg -> AddEntry(h[i], Form("Z VTX between %1.2f and %1.2f", -6 - static_cast<double>(i), -5 - static_cast<double>(i)), "l");
        } else {
            lg2 -> AddEntry(h[i], Form("Z VTX between %.2f and %.2f", -6 - static_cast<double>(i), -5 - static_cast<double>(i)), "l");
        }
        h[i] -> GetYaxis() -> SetRangeUser(min_y*0.95, 1.05);
        h[i] -> GetXaxis() -> CenterTitle(true);    h[i] -> GetYaxis() -> CenterTitle(true);
        h[i] -> SetLineWidth(2);
        h[i] -> SetLineColor(29 + i);
        h[i] -> Draw("SAME");
        h[i] -> SetTitle(fileTitle.c_str());
    }

    h[16] -> SetLineColor(2);
    double pi = TMath::Pi();
    int bin_min  = 1;  // The first bin
    int bin_max  = h[0]->GetNbinsX();  // The last bin
    // Calculate bin positions for each label
    int bin_pi   = bin_max;
    int bin_0    = bin_min + (bin_max - bin_min)/2;
    int binPi_2  = bin_min + 3*(bin_max - bin_min)/4;
    int bin_pi_2 = bin_min + (bin_max - bin_min)/4;
    // Set the labels at the calculated positions
    h[16]->GetXaxis()->SetBinLabel(bin_0, "0");
    h[16]->GetXaxis()->SetBinLabel(bin_pi_2, "#frac{-#pi}{2}");
    h[16]->GetXaxis()->SetBinLabel(binPi_2, "#frac{#pi}{2}");
    h[16]->GetXaxis()->SetBinLabel(bin_pi, "#pi");
    h[16]->GetXaxis()->SetBinLabel(bin_min, "-#pi");
    // Ensure the custom labels are displayed by setting the number of divisions
    h[16]->GetXaxis()->SetNdivisions(9, 0, 0, kFALSE);
    h[16]->GetXaxis()->SetLabelSize(0.04);
    // Update histogram to refresh the axis
    // h[0]->Draw("HIST");
    h[16]->GetXaxis()->LabelsOption("h"); // Draw the labels vertically

    lg->Draw("same");   lg2->Draw("same");
    can1 -> Update();

    can1 -> SaveAs(Form("../../External/zFindingPlots/%s.png", fileTitle.c_str()));
}

void ArrayPlot1D_Rescale_dEta (const std::vector<TH1D*>& h, std::vector<std::string> method, const std::string &fileTitle) {
    double minRange = h[0]->GetXaxis()->GetXmin();
    double maxRange = h[0]->GetXaxis()->GetXmax();
    int index = 0;
    double max = 0;
    TCanvas *can1 = new TCanvas("c1d","c1d",0,50,2100,1200);
    for (int i = 0; i < h.size(); i++) {
        double N = h[i] -> Integral(h[i]->FindFixBin(minRange), h[i]->FindFixBin(maxRange), "");
        h[i] -> Add(h[i], 1/(N/1000+1) - 1);
        h[i] -> GetXaxis() -> SetRangeUser(minRange*5/8, maxRange*5/8);
        h[i] -> GetXaxis() -> CenterTitle(true);    h[i] -> GetYaxis() -> CenterTitle(true);
        int maxEntry = h[i] -> GetBinContent(h[i]->GetMaximumBin());
        if (maxEntry > max) {
            max = maxEntry;
            index = i;
        } 
        // h[i] -> SetLineColor(30+i);
        h[i] -> SetLineWidth(3);
        // h[i] -> Draw("SAME");
    }
    std::cout << index << ", " << max << std::endl;
    h[index] -> Draw("SAME");
    h[index] -> GetYaxis() -> SetRangeUser(0, max*1.05);
    for (int i = 0; i < h.size(); i++) {
        if (i != index)     h[i] -> Draw("SAME");
    }
    h[0] -> SetLineColor(2);
    h[1] -> SetLineColor(4);
    h[2] -> SetLineColor(6);
    h[3] -> SetLineColor(8);
    h[4] -> SetLineColor(9);
    h[5] -> SetLineColor(12);
    h[6] -> SetLineColor(28);
    h[7] -> SetLineColor(30);
    h[8] -> SetLineColor(31);
    h[9] -> SetLineColor(32);
    h[10] -> SetLineColor(38);
    h[11] -> SetLineColor(40);
    h[12] -> SetLineColor(42);
    h[13] -> SetLineColor(46);

    can1 -> SaveAs(Form("../../External/zFindingPlots/%s.png", fileTitle.c_str()));
}

void ArrayPlot_dEta_1D_Logy (
    const std::vector<TH1D*> &h,
    const std::vector<std::string> &method,
    const std::string &fileName
) {
    TCanvas *c1 = new TCanvas("c1dEta", "dPhi Histogram", 1920, 1056);
    // h[0]->Draw("same");
    // h[1]->Draw("same");   h[1]->SetLineColor(kGreen+2);
    // h[2]->Draw("same");
    // h[3]->Draw("same");
    // h[4]->Draw("same");   h[4]->SetLineColor(2);
    // h[5]->Draw("same");
    // h[6]->Draw("same");
    // h[7]->Draw("same");
    // h[8]->Draw("same");
    // h[9]->Draw("same");
    // h[10]->Draw("same");  h[10]->SetLineColor(2);
    // h[11]->Draw("same");
    // h[12]->Draw("same");  h[12]->SetLineColor(kMagenta);
    // h[13]->Draw("same");

    int maxEntry = 0;
    for (int i = 0; i < h.size(); i++) {
            h[i]->Draw("same");
            h[i]->GetXaxis()->CenterTitle(true);    h[i]->GetYaxis()->CenterTitle(true);
            h[i]->SetLineWidth(3);
            h[i]->SetLineColor(20+i);
            int temp = h[i]->GetBinContent(h[i]->GetMaximumBin());
            if (temp >= maxEntry)   maxEntry = temp;
    }
    for (int i = 0; i < h.size(); i++)
        h[i]->GetYaxis()->SetRangeUser(1, maxEntry*1.2);
    TLine *l = new TLine(0, 0, 0, maxEntry);
	l -> Draw("same"); 
    l -> SetLineColor(kRed);
    l -> SetLineWidth(3);
    gPad->SetLogy();
    c1->Update();    c1->Modified();

    c1 -> SaveAs(Form("../../External/zFindingPlots/%s.png", fileName.c_str()));
}