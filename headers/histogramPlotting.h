#include "globalDefinitions.h"

void angularPlot1D (TH1D* const histo, std::vector<std::string> method, const std::string &fileTitle) {
    if (isInteger(method[1]))   int lowerRange = stoi(method[1]);
    if (isInteger(method[2]))   int upperRange = stoi(method[2]);
    int maxBin = histo->GetMaximumBin();
    double maxBinCenter = histo->GetBinCenter(maxBin);
    int maxEntry = histo -> GetBinContent(maxBin);
    TCanvas *can1 = new TCanvas("c1d","c1d",0,50,1800,1200);
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

    TLine *l = new TLine(0, 0, 0, maxEntry);
	l -> Draw("same"); 
    l -> SetLineColor(kRed);
    TLegend *lg = new TLegend(0.12, 0.8, 0.43, 0.9);
    lg -> AddEntry(histo, "Fiducial cut as -250 <= MBD_z_vtx <= -50 mm, MBD_centrality <= 0.70", "f");
    gStyle -> SetLegendTextSize(.015);
    lg->Draw("same");
    // gPad -> SetLogy();
    can1 -> SaveAs(Form("../External/zFindingPlots/%s.png", fileTitle.c_str()));
}

void angularPlot2D(TH2D *const h_dPhi_Z, std::vector<std::string> method, const std::string &fileTitle) {
    TCanvas *can = new TCanvas("c2d","c2d",0,50,1800,1200);
    h_dPhi_Z -> Draw("colz");
    h_dPhi_Z -> SetTitle(Form("All dPhi values (no event mixed) with %d bins of Z vtx", h_dPhi_Z->GetNbinsY()));
    // h_dPhi_Z -> GetXaxis() -> SetTitle("dPhi"); 
    h_dPhi_Z -> GetXaxis() -> CenterTitle(true);
    // h_dPhi_Z -> GetYaxis() -> SetTitle("found z vtx position [mm]");   
    h_dPhi_Z -> GetYaxis() -> CenterTitle(true);
    h_dPhi_Z -> GetXaxis() -> SetTitleOffset(1.8);   h_dPhi_Z -> GetYaxis() -> SetTitleOffset(1.8);
    can -> SaveAs(Form("../External/zFindingPlots/%s.png", fileTitle.c_str()));
}

void angularPlot3D(TH2D *const h_dPhi_Z, std::vector<std::string> method, const std::string &fileTitle) {
    TCanvas *can = new TCanvas("c3d","c3d",0,50,1800,1200);
    h_dPhi_Z -> Draw("lego2");
    h_dPhi_Z -> SetTitle(Form("All dPhi values (no event mixed) with %d bins of Z vtx", h_dPhi_Z->GetNbinsY()));
    // h_dPhi_Z -> GetXaxis() -> SetTitle("dPhi"); 
    h_dPhi_Z -> GetXaxis() -> CenterTitle(true);
    // h_dPhi_Z -> GetYaxis() -> SetTitle("found z vtx position [mm]");   
    h_dPhi_Z -> GetYaxis() -> CenterTitle(true);
    h_dPhi_Z -> GetXaxis() -> SetTitleOffset(1.8);   h_dPhi_Z -> GetYaxis() -> SetTitleOffset(1.8);
    can -> SaveAs(Form("../External/zFindingPlots/%s.png", fileTitle.c_str()));
}

void ZResolutionSinglePlot (TH1D* const histo, std::vector<std::string> method, const std::string &fileTitle) {
    if (isInteger(method[1]))   int lowerRange = stoi(method[1]);
    if (isInteger(method[2]))   int upperRange = stoi(method[2]);
    int maxBin = histo->GetMaximumBin();
    double maxBinCenter = histo->GetBinCenter(maxBin);
    int maxEntry = histo -> GetBinContent(maxBin);
    TCanvas *can1 = new TCanvas("c1d","c1d",0,50,1800,1200);
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
    TCanvas *can1 = new TCanvas("cg","cg",0,50,1800,1200);
    g0 -> SetMarkerStyle(29);
    g0 -> SetMarkerSize(2);
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
    g0 -> SetMinimum(-200); // Setting y range;
    g0 -> SetMaximum(200);  // Setting y range;
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

void EtaPhiSinglePlot (TH1D* const histo, std::vector<std::string> method, Int_t const & target) {
    if (isInteger(method[1]))   int lowerRange = stoi(method[1]);
    if (isInteger(method[2]))   int upperRange = stoi(method[2]);
    int maxBin = histo->GetMaximumBin();
    double maxBinCenter = histo->GetBinCenter(maxBin);
    TCanvas *can1 = new TCanvas("c1","c1",0,50,1800,1200);
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
    TCanvas *can1 = new TCanvas("c1","c1",0,50,1800,1200);

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
    TCanvas *can1 = new TCanvas("c1","c1",0,50,1800,1200);
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
    hDiff -> SetTitle("Background Subtraced dPhi");
    gStyle->SetEndErrorSize(6);
    gStyle->SetErrorX(0.5);
    
    gStyle -> SetTitleFont(100,"t");
    gStyle -> SetTitleSize(0.065,"t");

    //gPad -> SetLogy(1);
    TLine *l = new TLine(-2,0,2,0);
	l -> Draw("same"); 
    l -> SetLineColor(kGreen);
}