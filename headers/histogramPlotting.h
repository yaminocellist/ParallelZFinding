#include "globalDefinitions.h"

void singlePlot (TH1D* const histo, std::vector<std::string> method, Int_t const & target) {
    int upperRange = stoi(method[2]);
    int lowerRange = stoi(method[1]);
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
    // histo -> GetXaxis() -> SetRangeUser(-M_PI, M_PI); // Setting x range;
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