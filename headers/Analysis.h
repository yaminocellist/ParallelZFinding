#include "globalDefinitions.h"
#include <ROOT/RDataFrame.hxx>


void drawHistogram (const int &evt, TH1D *h, TH1D *h2) {
    TCanvas *can1 = new TCanvas("c1","c1",0,50,1800,550);
    can1 -> Divide(1, 2);
    can1->cd(1);
    h -> Draw();
    h -> SetFillColor(kYellow - 7);
    h -> SetLineWidth(1);
    h -> SetFillStyle(1001);
    h -> GetXaxis() -> SetTitle("Eta value");
    h -> GetXaxis() -> SetTitleSize(.05);
    h -> GetXaxis() -> SetLabelSize(.04);
    h -> GetXaxis() -> CenterTitle(true);
    h -> GetXaxis() -> SetNdivisions(31, 5, 0);
    h -> GetXaxis() -> SetTitleOffset(.8);
    h -> GetYaxis() -> SetTitle("# of Counts");
    h -> GetYaxis() -> SetTitleSize(.05);
    h -> GetYaxis() -> SetLabelSize(.04);
    h -> GetYaxis() -> SetTitleOffset(.62);
    h -> GetYaxis() -> CenterTitle(true);
    h -> SetTitle(Form("Eta data from AncG4P of Event %d", evt));
    can1 -> cd(2);
    h2 -> Draw();
    h2 -> SetFillColor(kBlue - 7);
    h2 -> SetLineWidth(1);
    h2 -> SetFillStyle(1001);
    h2 -> GetXaxis() -> SetTitle("Eta value");
    h2 -> GetXaxis() -> SetTitleSize(.05);
    h2 -> GetXaxis() -> SetLabelSize(.04);
    h2 -> GetXaxis() -> CenterTitle(true);
    h2 -> GetXaxis() -> SetNdivisions(31, 5, 0);
    h2 -> GetXaxis() -> SetTitleOffset(.8);
    h2 -> GetYaxis() -> SetTitle("# of Counts");
    h2 -> GetYaxis() -> SetTitleSize(.05);
    h2 -> GetYaxis() -> SetLabelSize(.04);
    h2 -> GetYaxis() -> SetTitleOffset(.62);
    h2 -> GetYaxis() -> CenterTitle(true);
    h2 -> SetTitle(Form("Eta data from Cluster of Event %d", evt));
}

void dEtaCheck (const int &evt, std::vector<myTrackletMember> t0, std::vector<myTrackletMember> t1) {
    TH1D *h_dEta = new TH1D("", "", 1601, -4 - .0025, 4 + .0025);
    for (int i = 0; i < t0.size(); i++) {
        for (int j = 0; j < t1.size(); j++) {
            h_dEta -> Fill(t0[i].eta - t1[j].eta);
        }
    }
    int maxBin = h_dEta->GetMaximumBin();
    double maxBinCenter = h_dEta->GetBinCenter(maxBin);
    TCanvas *can1 = new TCanvas("c1","c1",0,50,1800,1200);
    h_dEta -> Draw();
    h_dEta -> SetFillColor(kYellow - 7);
    h_dEta -> SetLineWidth(1);
    h_dEta -> SetFillStyle(1001);
    h_dEta -> GetXaxis() -> SetTitle("Eta value");
    h_dEta -> GetXaxis() -> SetTitleSize(.05);
    h_dEta -> GetXaxis() -> SetLabelSize(.03);
    h_dEta -> GetXaxis() -> CenterTitle(true);
    h_dEta -> GetXaxis() -> SetNdivisions(31, 5, 0);
    h_dEta -> GetXaxis() -> SetTitleOffset(.8);
    h_dEta -> GetYaxis() -> SetTitle("# of Counts");
    h_dEta -> GetYaxis() -> SetTitleSize(.05);
    h_dEta -> GetYaxis() -> SetLabelSize(.03);
    h_dEta -> GetYaxis() -> SetTitleOffset(.8);
    h_dEta -> GetYaxis() -> CenterTitle(true);
    h_dEta -> GetXaxis() -> SetRangeUser(-3., +3.); // Setting x range;
    static TLine *l1 = new TLine();
    l1 -> SetLineColor(kRed);
    l1 -> SetLineStyle(0);
    l1 -> SetLineWidth(1);
    l1 -> SetX1(maxBinCenter);   l1 -> SetY1(0);
    l1 -> SetX2(maxBinCenter);   l1 -> SetY2(h_dEta->GetMaximum());
    // l1 -> Draw("same");

    h_dEta -> SetTitle(Form("dEta data of Event %d, centered at %0.4f", evt, maxBinCenter));
    // gPad -> SetLogy();
    can1 -> SaveAs(Form("../External/xyFindingPlots/dEta_single_%d.png", evt));
}

void dEtaCheckAll (TH1D* const h_dEta) {
    int maxBin = h_dEta->GetMaximumBin();
    double maxBinCenter = h_dEta->GetBinCenter(maxBin);
    TCanvas *can1 = new TCanvas("c1","c1",0,50,1800,1200);
    h_dEta -> Draw();
    h_dEta -> SetFillColor(kYellow - 7);
    h_dEta -> SetLineWidth(1);
    h_dEta -> SetFillStyle(1001);
    h_dEta -> GetXaxis() -> SetTitle("Eta value");
    h_dEta -> GetXaxis() -> SetTitleSize(.05);
    h_dEta -> GetXaxis() -> SetLabelSize(.03);
    h_dEta -> GetXaxis() -> CenterTitle(true);
    h_dEta -> GetXaxis() -> SetNdivisions(31, 5, 0);
    h_dEta -> GetXaxis() -> SetTitleOffset(.8);
    h_dEta -> GetYaxis() -> SetTitle("# of Counts");
    h_dEta -> GetYaxis() -> SetTitleSize(.05);
    h_dEta -> GetYaxis() -> SetLabelSize(.03);
    h_dEta -> GetYaxis() -> SetTitleOffset(.8);
    h_dEta -> GetYaxis() -> CenterTitle(true);
    h_dEta -> GetXaxis() -> SetRangeUser(-3., +3.); // Setting x range;
    static TLine *l1 = new TLine();
    l1 -> SetLineColor(kRed);
    l1 -> SetLineStyle(0);
    l1 -> SetLineWidth(1);
    l1 -> SetX1(maxBinCenter);   l1 -> SetY1(0);
    l1 -> SetX2(maxBinCenter);   l1 -> SetY2(h_dEta->GetMaximum());
    // l1 -> Draw("same");

    h_dEta -> SetTitle(Form("dEta data of all events, centered at %0.4f", maxBinCenter));
    // gPad -> SetLogy();

    can1 -> SaveAs("../External/xyFindingPlots/dEta_all.png");
}

void dPhiCheck (const int &evt, std::vector<myTrackletMember> t0, std::vector<myTrackletMember> t1) {
    // TH1D *h_dphi = new TH1D("", "", 1601, -4 - .0025, 4 + .0025);
    int N = 1000;  // Choose an odd number of bins
    double range_min = -M_PI;
    double range_max = M_PI;
    double bin_width = (range_max - range_min) / N;
    TH1D *h_dphi = new TH1D("", "", N, range_min, range_max);
    for (int i = 0; i < t0.size(); i++) {
        for (int j = 0; j < t1.size(); j++) {
            if (t0[i].phi - t1[j].phi > M_PI) {
                h_dphi -> Fill(t0[i].phi - t1[j].phi - 2*M_PI);
            }
            else if (t0[i].phi - t1[j].phi < -M_PI) {
                h_dphi -> Fill(t0[i].phi - t1[j].phi + 2*M_PI);
            }
            else {
                h_dphi -> Fill(t0[i].phi - t1[j].phi);
            }
                // h_dphi -> Fill(t0[i].phi - t1[j].phi);
        }
    }
    int maxBin = h_dphi->GetMaximumBin();
    double maxBinCenter = h_dphi->GetBinCenter(maxBin);
    TCanvas *can1 = new TCanvas("c1","c1",0,50,1800,1200);
    h_dphi -> Draw();
    h_dphi -> SetFillColor(kYellow - 7);
    h_dphi -> SetLineWidth(1);
    h_dphi -> SetFillStyle(1001);
    h_dphi -> GetXaxis() -> SetTitle("Phi value");
    h_dphi -> GetXaxis() -> SetTitleSize(.05);
    h_dphi -> GetXaxis() -> SetLabelSize(.03);
    h_dphi -> GetXaxis() -> CenterTitle(true);
    h_dphi -> GetXaxis() -> SetNdivisions(31, 5, 0);
    h_dphi -> GetXaxis() -> SetTitleOffset(.8);
    h_dphi -> GetYaxis() -> SetTitle("# of Counts");
    h_dphi -> GetYaxis() -> SetTitleSize(.05);
    h_dphi -> GetYaxis() -> SetLabelSize(.03);
    h_dphi -> GetYaxis() -> SetTitleOffset(.8);
    h_dphi -> GetYaxis() -> CenterTitle(true);
    // h_dphi -> GetXaxis() -> SetRangeUser(-3., 3.); // Setting x range;
    h_dphi -> SetTitle(Form("dPhi data of Event %d, centered at %0.4f", evt, maxBinCenter));
    gPad -> SetLogy();
    double pi = TMath::Pi();
    int bin_min = 1;  // The first bin
    int bin_max = h_dphi->GetNbinsX();  // The last bin

// Calculate bin positions for each label
int bin_pi = bin_max;
int bin_0 = bin_min + (bin_max - bin_min)/2;
int binPi_2 = bin_min + 3*(bin_max - bin_min)/4;
// int bin_pi_2 = bin_min + (bin_max - bin_min)/2;
int bin_pi_2 = bin_min + (bin_max - bin_min)/4;



// Set the labels at the calculated positions
h_dphi->GetXaxis()->SetBinLabel(bin_0, "0");
h_dphi->GetXaxis()->SetBinLabel(bin_pi_2, "#frac{-#pi}{2}");
h_dphi->GetXaxis()->SetBinLabel(binPi_2, "#frac{#pi}{2}");
// h_dphi->GetXaxis()->SetBinLabel(bin_3pi_4, "#frac{3#pi}{4}");
h_dphi->GetXaxis()->SetBinLabel(bin_pi, "#pi");
h_dphi->GetXaxis()->SetBinLabel(bin_min, "-#pi");
// h_dphi->GetXaxis()->SetBinLabel(bin_max - 3*(bin_max - bin_min)/4, "-#frac{3#pi}{4}");
// h_dphi->GetXaxis()->SetBinLabel(bin_max - (bin_max - bin_min)/2, "-#frac{#pi}{2}");
// h_dphi->GetXaxis()->SetBinLabel(bin_max - (bin_max - bin_min)/4, "-#frac{#pi}{4}");

// Ensure the custom labels are displayed by setting the number of divisions
h_dphi->GetXaxis()->SetNdivisions(9, 0, 0, kFALSE);

// Update histogram to refresh the axis
h_dphi->Draw("HIST");
h_dphi->GetXaxis()->LabelsOption("h"); // Draw the labels vertically
    can1 -> SaveAs(Form("../External/xyFindingPlots/dPhi_single_%d.png", evt));
}

void dPhiCheckAll (TH1D* const h_dphi, std::vector<std::string> method, Int_t const & target) {
    int upperRange = stoi(method[2]);
    int lowerRange = stoi(method[1]);
    int maxBin = h_dphi->GetMaximumBin();
    double maxBinCenter = h_dphi->GetBinCenter(maxBin);
    TCanvas *can1 = new TCanvas("c1","c1",0,50,1800,1200);
    h_dphi -> Draw();
    h_dphi -> SetFillColor(kYellow - 7);
    h_dphi -> SetLineWidth(1);
    h_dphi -> SetFillStyle(1001);
    h_dphi -> GetXaxis() -> SetTitle("Phi value");
    h_dphi -> GetXaxis() -> SetTitleSize(.05);
    h_dphi -> GetXaxis() -> SetLabelSize(.03);
    h_dphi -> GetXaxis() -> CenterTitle(true);
    h_dphi -> GetXaxis() -> SetNdivisions(31, 5, 0);
    h_dphi -> GetXaxis() -> SetTitleOffset(.8);
    h_dphi -> GetYaxis() -> SetTitle("# of Counts");
    h_dphi -> GetYaxis() -> SetTitleSize(.05);
    h_dphi -> GetYaxis() -> SetLabelSize(.03);
    h_dphi -> GetYaxis() -> SetTitleOffset(.8);
    h_dphi -> GetYaxis() -> CenterTitle(true);
    // h_dphi -> GetXaxis() -> SetRangeUser(-M_PI, M_PI); // Setting x range;
    // h_dphi -> GetYaxis() -> SetRangeUser(5300e3, 6130e3);
    // h_dphi -> SetTitle(Form("dPhi data of all events whose found z vtx is ~ [-%d, -%d] mm, centered at %0.4f", lowerRange, upperRange, maxBinCenter));
    // gPad -> SetLogy();
    // can1 -> SaveAs(Form("../External/xyFindingPlots/dPhi_all_%d_%d.png", lowerRange, upperRange));
    h_dphi-> SetTitle(Form("dPhi of %d events mixed up in range of [-209.375, -207.5] mm", target));
    can1 -> SaveAs(Form("../External/xyFindingPlots/dPhi_mixed_%d.png", target));
}

void dPhiCheckDouble (TH1D* const h_dphi, TH1D* const h_phi, std::vector<std::string> method, Int_t const & target) {
    int upperRange = stoi(method[2]);
    int lowerRange = stoi(method[1]);
    int maxBin = h_dphi->GetMaximumBin();
    double maxBinCenter = h_dphi->GetBinCenter(maxBin);
    TCanvas *can1 = new TCanvas("c1","c1",0,50,1800,1200);

    

    // Draw the first histogram
    h_dphi -> Draw();
    h_dphi -> SetFillColor(kYellow - 7);
    h_dphi -> SetLineWidth(1);
    h_dphi -> SetFillStyle(1001);
    h_dphi -> GetXaxis() -> SetTitle("Phi value");
    h_dphi -> GetXaxis() -> SetTitleSize(.05);
    h_dphi -> GetXaxis() -> SetLabelSize(.03);
    h_dphi -> GetXaxis() -> CenterTitle(true);
    h_dphi -> GetXaxis() -> SetNdivisions(31, 5, 0);
    h_dphi -> GetXaxis() -> SetTitleOffset(.8);
    h_dphi -> GetYaxis() -> SetTitle("# of Counts");
    h_dphi -> GetYaxis() -> SetTitleSize(.05);
    h_dphi -> GetYaxis() -> SetLabelSize(.03);
    h_dphi -> GetYaxis() -> SetTitleOffset(.8);
    h_dphi -> GetYaxis() -> CenterTitle(true);
    h_dphi -> GetYaxis() -> SetRangeUser(0, 6200e3); // Setting x range;

    // Draw the second histogram on the same canvas
    h_phi -> Draw("SAME");
    h_phi -> SetFillColor(kBlue - 7);
    h_phi -> SetLineWidth(1);
    h_phi -> SetFillStyle(1001);

    // h_dphi-> SetTitle(Form("dPhi of %d events mixed up in range of [-209.375, -207.5] mm", target));
    // can1 -> SaveAs(Form("../External/xyFindingPlots/dPhi_mixed_%d.png", target));
}

void dRCheckAll (TH1D* const h) {
    int bin1 = h->FindBin(0.01);  // find the bin number corresponding to 0.0
    int bin2 = h->FindBin(0.02);  // find the bin number corresponding to 0.02
    double minContent = h->GetBinContent(bin1);
    int minBin = bin1;  
    for (int i = bin1 + 1; i <= bin2; ++i) {
        double currentContent = h -> GetBinContent(i);
        if (currentContent < minContent) {
            minContent = currentContent;
            minBin = i;
        }
    }

    TCanvas *can2 = new TCanvas("c2","c2",0,50,1800,550);
    h -> Draw();
    h -> SetFillColor(kYellow - 7);
    h -> SetLineWidth(1);
    h -> SetFillStyle(1001);
    h -> GetXaxis() -> SetTitle("Angular Distance");
    h -> GetXaxis() -> SetTitleSize(.05);
    h -> GetXaxis() -> SetLabelSize(.04);
    h -> GetXaxis() -> CenterTitle(true);
    h -> GetXaxis() -> SetNdivisions(31, 5, 0);
    h -> GetXaxis() -> SetTitleOffset(.8);
    h -> GetYaxis() -> SetTitle("# of Counts");
    h -> GetYaxis() -> SetTitleSize(.05);
    h -> GetYaxis() -> SetLabelSize(.04);
    h -> GetYaxis() -> SetTitleOffset(.62);
    h -> GetYaxis() -> CenterTitle(true);
    h -> GetXaxis() -> SetRangeUser(0, 0.2); // Setting x range;
    h -> GetYaxis()->SetRangeUser(0, 5e6);
    gPad -> SetGrid(1,1); gPad -> Update();
    // gPad -> SetLogx();
    can2 -> SaveAs("../External/xyFindingPlots/angular_distance_all.png");
}

void DCCheck (const int &evt, std::vector<myTrackletMember> t0, std::vector<myTrackletMember> t1) {
    Double_t binEdges[401];
    double dx_1 = 0.0002;
    binEdges[0] = 0.;
    for (int i = 1; i <= 100; i++) {
        binEdges[i] = binEdges[0] + i*dx_1;
    }
    double dx_2 = 0.001;
    for (int j = 1; j <= 100; j++) {
        binEdges[j+100] = binEdges[100] + j*dx_2;
    }
    double dx_3 = 0.01;
    for (int l = 1; l <= 200; l++) {
        binEdges[l + 200] = binEdges[200] + l*dx_3;
    }
    // TH1D *h_CD = new TH1D("", "", 400, binEdges);
    TH1D *h_CD = new TH1D("", "", 400, 0., 7.);

    for (int i = 0; i < t0.size(); i++) {
        for (int j = 0; j < t1.size(); j++) {
            myPoint3D p0 = {t0[i].x, t0[i].y, t0[i].z};
            myPoint3D p1 = {t1[j].x, t1[j].y, t1[j].z};
            double closest_distance = nearestZ(p0, p1).second;
            h_CD -> Fill(closest_distance);
        }
    }

    TCanvas *can1 = new TCanvas("c1","c1",0,50,1800,550);
    h_CD -> Draw();
    h_CD -> SetFillColor(kYellow - 7);
    h_CD -> SetLineWidth(1);
    h_CD -> SetFillStyle(1001);
    h_CD -> GetXaxis() -> SetTitle("Closest Distance values");
    h_CD -> GetXaxis() -> SetTitleSize(.05);
    h_CD -> GetXaxis() -> SetLabelSize(.04);
    h_CD -> GetXaxis() -> CenterTitle(true);
    h_CD -> GetXaxis() -> SetNdivisions(31, 5, 0);
    h_CD -> GetXaxis() -> SetTitleOffset(.8);
    h_CD -> GetYaxis() -> SetTitle("# of Counts");
    h_CD -> GetYaxis() -> SetTitleSize(.05);
    h_CD -> GetYaxis() -> SetLabelSize(.04);
    h_CD -> GetYaxis() -> SetTitleOffset(.62);
    h_CD -> GetYaxis() -> CenterTitle(true);
    // h_CD -> GetXaxis() -> SetRangeUser(0, 0.1); // Setting x range;
    h_CD -> SetTitle(Form("Closest Distances of Event %d", evt));
    gPad -> SetLogy();  gPad -> SetLogx();
    can1 -> SaveAs(Form("../External/zFindingPlot/cd_single_%d.png", evt));
}

void PhiCheck (const int &evt, std::vector<myTrackletMember> t0, std::vector<myTrackletMember> t1,
                const double &eta_cut_low, const double &eta_cut_high,
                const double &phi_cut_low, const double &phi_cut_high, const double &trueZ) {
    TH1D *h_phi = new TH1D("", "", 70, -3.5 - .05, 3.5 + .05);

    for (int i = 0; i < t0.size(); i++) {
        h_phi -> Fill(t0[i].phi);
    }
    for (int i = 0; i < t1.size(); i++) {
        h_phi -> Fill(t1[i].phi);
    }

    TCanvas *can1 = new TCanvas("c1","c1",0,50,1800,550);
    h_phi -> Draw();
    h_phi -> SetFillColor(kYellow - 7);
    h_phi -> SetLineWidth(1);
    h_phi -> SetFillStyle(1001);
    h_phi -> GetXaxis() -> SetTitle("Phi value");
    h_phi -> GetXaxis() -> SetTitleSize(.05);
    h_phi -> GetXaxis() -> SetLabelSize(.04);
    h_phi -> GetXaxis() -> CenterTitle(true);
    h_phi -> GetXaxis() -> SetNdivisions(31, 5, 0);
    h_phi -> GetXaxis() -> SetTitleOffset(.8);
    h_phi -> GetYaxis() -> SetTitle("# of Counts");
    h_phi -> GetYaxis() -> SetTitleSize(.05);
    h_phi -> GetYaxis() -> SetLabelSize(.04);
    h_phi -> GetYaxis() -> SetTitleOffset(.62);
    h_phi -> GetYaxis() -> CenterTitle(true);
    h_phi -> SetTitle(Form("Phi data of Event %d", evt));
    gPad -> SetLogy();
    can1 -> SaveAs(Form("../External/zFindingPlots/Phi_single_%d.png", evt));
}

double EtaCheck (const int &evt, std::vector<myTrackletMember> t0, std::vector<myTrackletMember> t1,
                const double &eta_cut_low, const double &eta_cut_high,
                const double &phi_cut_low, const double &phi_cut_high, const double &trueZ) {
    TH1D *h_eta = new TH1D("", "", 80, -4 - .05, 4 + .05);

    for (int i = 0; i < t0.size(); i++) {
        h_eta -> Fill(t0[i].eta);
    }
    for (int i = 0; i < t1.size(); i++) {
        h_eta -> Fill(t1[i].eta);
    }

    double firstNonZeroBinEdge = -4.05; // Start with the minimum possible value
    double lastNonZeroBinEdge = 4.05; // Start with the maximum possible value
    bool foundFirstNonZeroBin = false;

    for (int bin = 1; bin <= h_eta->GetNbinsX(); bin++) { // Bins are numbered from 1 to N
        if (h_eta->GetBinContent(bin) > 0) {
            if (!foundFirstNonZeroBin) {
                firstNonZeroBinEdge = h_eta->GetBinLowEdge(bin);
                foundFirstNonZeroBin = true;
            }
            // For the last non-zero bin, keep updating this value until the last iteration with non-zero content
            lastNonZeroBinEdge = h_eta->GetBinLowEdge(bin) + h_eta->GetBinWidth(bin);
        }
    }

    /*
    TCanvas *can1 = new TCanvas("c1","c1",0,50,1800,550);
    h_eta -> Draw();
    h_eta -> SetFillColor(kYellow - 7);
    h_eta -> SetLineWidth(1);
    h_eta -> SetFillStyle(1001);
    h_eta -> GetXaxis() -> SetTitle("Eta value");
    h_eta -> GetXaxis() -> SetTitleSize(.05);
    h_eta -> GetXaxis() -> SetLabelSize(.04);
    h_eta -> GetXaxis() -> CenterTitle(true);
    h_eta -> GetXaxis() -> SetNdivisions(31, 5, 0);
    h_eta -> GetYaxis() -> SetTitle("# of Counts");
    h_eta -> GetYaxis() -> SetTitleSize(.05);
    h_eta -> GetYaxis() -> SetLabelSize(.04);
    h_eta -> GetYaxis() -> SetTitleOffset(.62);
    h_eta -> GetYaxis() -> CenterTitle(true);
    h_eta -> SetTitle(Form("All #eta data of Event %d", evt));
    gPad -> SetLogy();
    can1 -> SaveAs(Form("../External/zFindingPlots/Eta_single_%d.png", evt));
    */
    return lastNonZeroBinEdge - firstNonZeroBinEdge;
}

void foundZAnalysisLite (string filePath, string select = "zscan", string dim = "2D") {
    // if (select == "zscanE") {
    //     filePath = "foundZ_debug_zscan_with_error.txt";
    // } else if (select == "zscan") {
    //     filePath = "foundZ_debug_zscan.txt";
    // } else if (select == "e") {
    //     filePath = "foundZ_debug_nearest_z_with_error.txt";
    // } else {
        // filePath = "/Users/yaminocellist/MIT_mentorship/3rd_semester/Meeting5/foundZ_debug2_nearest_z.txt";
    // }
    ifstream myfile(filePath);
    if (!myfile.is_open()){
		  cout << "Unable to open linelabel" << endl;
		  system("read -n 1 -s -p \"Press any key to continue...\" echo");
		  exit(1);
 	}

    TH1D *h = new TH1D("", "", 151, -1e2, 1e2);

    int evt, Nparticles;
    string line, substr;
    double found_z, true_z;
    std::vector<double> foundz, dZ;
    std::vector<double> truez, uppersize, lowersize, totalsize;

    while (getline(myfile, line)) {
        stringstream str(line);
        getline(str, substr, ',');
        getline(str, substr, ',');
        evt = stoi(substr);
        getline(str, substr, ',');
        Nparticles = stoi(substr);
        getline(str, substr, ',');
        found_z = stod(substr)*1e1;
        getline(str, substr, ',');
        getline(str, substr, ',');
        getline(str, substr, ',');
        true_z = stod(substr)*1e1;

        foundz.push_back(found_z);
        truez.push_back(true_z);
        totalsize.push_back(Nparticles);
        dZ.push_back(found_z - true_z);
        h -> Fill(found_z - true_z);
    }
    
    TGraph *g0 = new TGraph(dZ.size(), totalsize.data(), dZ.data());
    TGraph *g1 = new TGraph(dZ.size(), truez.data(), dZ.data());
    TCanvas *canvas = new TCanvas("c2","c2", 0, 0,1350,800);
    
    if (dim == "1D") {
        h -> Draw();
        h -> SetFillColor(kYellow - 7);
        h -> SetLineWidth(1);
        h -> SetFillStyle(1001);
        h -> GetXaxis() -> SetTitle("foundZ - trueZ [mm]");
        h -> GetXaxis() -> SetTitleSize(.05);
        h -> GetXaxis() -> SetLabelSize(.03);
        h -> GetXaxis() -> CenterTitle(true);
        h -> GetXaxis() -> SetNdivisions(31, 5, 0);
        h -> GetXaxis() -> SetTitleOffset(.8);
        h -> GetXaxis() -> SetRangeUser(-50, 50); // Setting x range;
        h -> GetYaxis() -> SetTitle("# of Counts");
        h -> GetYaxis() -> SetTitleSize(.05);
        h -> GetYaxis() -> SetLabelSize(.04);
        h -> GetYaxis() -> SetTitleOffset(.8);
        h -> GetYaxis() -> CenterTitle(true);

        static TLine *l1 = new TLine();
        l1 -> SetLineColor(kRed);
        l1 -> SetLineStyle(0);
        l1 -> SetLineWidth(2);
        l1 -> SetX1(0);   l1 -> SetY1(0);
        l1 -> SetX2(0);   l1 -> SetY2(h->GetMaximum());
        l1 -> Draw("same");

        static TLine *l2 = new TLine();
        l2 -> SetLineColor(kRed);
        l2 -> SetLineStyle(2);
        l2 -> SetLineWidth(4);
        l2 -> SetX1(-16);   l2 -> SetY1(0);
        l2 -> SetX2(-16);   l2 -> SetY2(h->GetMaximum());
        l2 -> Draw("same");

        static TLine *l3 = new TLine();
        l3 -> SetLineColor(kRed);
        l3 -> SetLineStyle(2);
        l3 -> SetLineWidth(4);
        l3 -> SetX1(+16);   l3 -> SetY1(0);
        l3 -> SetX2(+16);   l3 -> SetY2(h->GetMaximum());
        l3 -> Draw("same");
        TLegend *legend = new TLegend(0.5,0.8,0.9,0.9);
        // legend->SetHeader("Legend","C"); // option "C" allows to center the header
        legend->AddEntry(h,"Segments expanded between [-8, 8]mm","f");
        legend->Draw("same");
        legend->SetTextSize(0.03);

        if (select == "zscanE") {
            h -> SetTitle("Z Scan with Errors, |dEta| < 0.01 (for small N) or 0.001 (for large N), |dPhi| < 0.1");
        } else if (select == "zscan") {
            h -> SetTitle("Traditional Z Scan method, |dEta| < 0.01 (for small N) or 0.001 (for large N), |dPhi| < 0.1");
        } else if (select == "e") {
            h -> SetTitle("Nearest Z with Errors, |dPhi| < 0.1");
        } else if (select == "zscanN") {
            h -> SetTitle("Z Scan with Background Normalized");
        }
        else {
            h -> SetTitle(Form("DCA indicated by greatest bin, %lu events, |dPhi| < 0.01,  closest distance cut = %1.2f", foundz.size(), DCA_cut));
            // h -> SetTitle(Form("DCA indicated with Gaussian Fit, |dPhi| < 0.01, closest distance cut = %1.2f", DCA_cut));
        }
    }
    else if (dim == "2D") {
        g0 -> SetMarkerStyle(29);
        g0 -> SetMarkerSize(2);
        g0 -> SetMarkerColor(kBlue);
        g0 -> SetLineWidth(3);
        g0 -> SetLineColor(kWhite);
        if (select == "zscanE") {
            g0 -> SetTitle("Z Scan with Errors, |dEta| < 0.01 (for small N) or 0.001 (for large N), |dPhi| < 0.1");
        } else if (select == "zscan") {
            g0 -> SetTitle("Traditional Z Scan method, |dEta| < 0.01 (for small N) or 0.001 (for large N), |dPhi| < 0.1");
        } else if (select == "e") {
            g0 -> SetTitle("Nearest Z with Errors, |dPhi| < 0.1");
        } 
        else if (select == "zscanN") {
            g0 -> SetTitle("Z Scan with Background Normalized");
        }
        else {
            g0 -> SetTitle(Form("DCA indicated by greatest bin, %lu events, |dPhi| < 0.01, closest distance cut = %1.2f", foundz.size(), DCA_cut));
            // g0 -> SetTitle(Form("DCA indicated with Gaussian Fit, |dPhi| < 0.01, closest distance cut = %1.2f", DCA_cut));
        }
        
        gStyle -> SetTitleW(0.7);  //per cent of the pad width
        gStyle -> SetTitleH(0.08); //per cent of the pad height
        g0 -> GetXaxis() -> SetTitle("Number of particles");
        g0 -> GetXaxis() -> SetTitleSize(0.06);
        g0 -> GetXaxis() -> SetLabelSize(0.04);
        g0 -> GetXaxis() -> CenterTitle(true);
        g0 -> GetYaxis() -> SetTitle("Found_z - True_z [mm]");
        g0 -> GetYaxis() -> SetTitleSize(0.06);
        g0 -> GetYaxis() -> SetLabelSize(0.025);
        g0 -> GetYaxis() -> CenterTitle(true);
        g0 -> SetMinimum(-300); // Setting y range;
        g0 -> SetMaximum(300);  // Setting y range;
        g0 -> GetYaxis() -> SetTitleOffset(0.8); 
        g0 -> GetXaxis() -> SetTitleOffset(0.8); 
        // g0 -> GetXaxis() -> SetLimits(0, 6000); // Setting x range;
        g0 -> Draw("AP SAME");

        gPad->SetGrid(5, 2); gPad->Update();

        TLine *line1 = new TLine(0, -16, 10000, -16);
        line1->SetLineColor(kRed);
        line1->SetLineWidth(2);
        line1->Draw("same");  // Draw line on the same canvas
        TLine *line2 = new TLine(0, 16, 10000, 16);
        line2->SetLineColor(kRed);
        line2->SetLineWidth(2);
        line2->Draw("same");  // Draw line on the same canvas
        TLegend *legend = new TLegend(0.4,0.8,0.9,0.9);
        // legend->SetHeader("Legend","C"); // option "C" allows to center the header
        legend->AddEntry(g0,"Segments expanded between [-8, 8]mm","f");
        legend->SetTextSize(0.04);
        legend->Draw("same");
    }
    else if (dim == "zVSdz") {
        g1 -> SetMarkerStyle(29);
        g1 -> SetMarkerSize(1);
        g1 -> SetMarkerColor(kBlue - 7);
        g1 -> SetLineWidth(3);
        g1 -> SetLineColor(kWhite);
        if (select == "zscanE") {
            g1 -> SetTitle("Z Scan with Errors, |dEta| < 0.01 (for small N) or 0.001 (for large N), |dPhi| < 0.1");
        } else if (select == "zscan") {
            g1 -> SetTitle("Traditional Z Scan method, |dEta| < 0.01 (for small N) or 0.001 (for large N), |dPhi| < 0.1");
        } else if (select == "e") {
            g1 -> SetTitle("Nearest Z with Errors, |dPhi| < 0.1");
        } 
        else if (select == "zscanN") {
            g1 -> SetTitle("Z Scan with Background Normalized");
        }
        else if (select == "DCA") {
            g1 -> SetTitle(Form("DCA indicated by greatest bin, %lu events, |dPhi| < 0.01, closest distance cut = %1.2f", foundz.size(), DCA_cut));
        }
        else {
            g1 -> SetTitle(Form("DCA indicated with Gaussian Fit, |dPhi| < 0.01, closest distance cut = %1.2f", DCA_cut));
        }
        
        gStyle -> SetTitleW(0.7);  //per cent of the pad width
        gStyle -> SetTitleH(0.08); //per cent of the pad height
        g1 -> GetXaxis() -> SetTitle("True Z Vertex [mm]");
        g1 -> GetXaxis() -> SetTitleSize(0.06);
        g1 -> GetXaxis() -> SetLabelSize(0.04);
        g1 -> GetXaxis() -> CenterTitle(true);
        g1 -> GetYaxis() -> SetTitle("Found_z - True_z [mm]");
        g1 -> GetYaxis() -> SetTitleSize(0.06);
        g1 -> GetYaxis() -> SetLabelSize(0.025);
        g1 -> GetYaxis() -> CenterTitle(true);
        g1 -> SetMinimum(-100); // Setting y range;
        g1 -> SetMaximum(100);  // Setting y range;
        g1 -> GetYaxis() -> SetTitleOffset(0.8); 
        g1 -> GetXaxis() -> SetTitleOffset(0.8); 
        // g1 -> GetXaxis() -> SetLimits(0, 6000); // Setting x range;
        g1 -> Draw("AP SAME");

        gPad->SetGrid(5, 2); gPad->Update();

        TLine *line1 = new TLine(-500, 0, +500, 0);
        line1->SetLineColor(kRed);
        line1->SetLineWidth(2);
        line1->Draw("same");  // Draw line on the same canvas
        TLegend *legend = new TLegend(0.4,0.8,0.9,0.9);
        // legend->SetHeader("Legend","C"); // option "C" allows to center the header
        legend->AddEntry(g1,"Segments expanded between [0, 16] mm","f");
        legend->SetTextSize(0.04);
        legend->Draw("same");
    }
    else {
        double ymin = g1->GetHistogram()->GetMinimum();
        double ymax = g1->GetHistogram()->GetMaximum();
        int nbins = 1000;
        TH1F *h1 = new TH1F("ProjectionY", "Histogram of vtx_z_resolution vs vtx_z 1D projection;Found_z - True_z [mm];Entries", nbins, ymin, ymax);
        int npoints = g1->GetN();
        double x, y;
        for (int i = 0; i < npoints; i++) {
            g1->GetPoint(i, x, y);
            h1->Fill(y);
        }
        h1->Draw();
        h1 -> GetXaxis() -> CenterTitle(true);
        h1 -> GetYaxis() -> CenterTitle(true);
        h1 -> GetXaxis() -> SetRangeUser(-60, 60); // Setting x range;
        h1 -> SetFillColor(kYellow - 7);
        h1 -> SetLineWidth(1);
        h1 -> SetFillStyle(1001);
        static TLine *l1 = new TLine();
        l1 -> SetLineColor(kRed);
        l1 -> SetLineStyle(0);
        l1 -> SetLineWidth(2);
        l1 -> SetX1(0);   l1 -> SetY1(0);
        l1 -> SetX2(0);   l1 -> SetY2(h1->GetMaximum());
        l1 -> Draw("same");
        gPad -> SetLogy();
    }
}

void foundZAnalysis (string filePath, string select = "zscan", string dim = "2D") {
    // if (select == "zscanE") {
    //     filePath = "foundZ_debug_zscan_with_error.txt";
    // } else if (select == "zscan") {
    //     filePath = "foundZ_debug_zscan.txt";
    // } else if (select == "e") {
    //     filePath = "foundZ_debug_nearest_z_with_error.txt";
    // } else {
        // filePath = "/Users/yaminocellist/MIT_mentorship/3rd_semester/Meeting5/foundZ_debug2_nearest_z.txt";
    // }
    ifstream myfile(filePath);
    if (!myfile.is_open()){
		  cout << "Unable to open linelabel" << endl;
		  system("read -n 1 -s -p \"Press any key to continue...\" echo");
		  exit(1);
 	}

    TH1D *h = new TH1D("", "", 151, -1e2, 1e2);

    int evt, Nparticles;
    string line, substr;
    double found_z, true_z;
    std::vector<double> foundz, dZ;
    std::vector<double> truez, uppersize, lowersize, totalsize;
    getline(myfile, line);

    while (getline(myfile, line)) {
        stringstream str(line);
        getline(str, substr, ',');
        evt = stoi(substr);
        getline(str, substr, ',');  // index;
        getline(str, substr, ',');
        Nparticles = stoi(substr);
        getline(str, substr, ',');
        found_z = stod(substr)*1e1;
        getline(str, substr, ',');  // True_x;
        getline(str, substr, ',');  // True_y;
        getline(str, substr, ',');
        true_z = stod(substr)*1e1;

        foundz.push_back(found_z);
        truez.push_back(true_z);
        totalsize.push_back(Nparticles);
        dZ.push_back(found_z - true_z);
        h -> Fill(found_z - true_z);
    }
    
    TGraph *g0 = new TGraph(dZ.size(), totalsize.data(), dZ.data());
    TGraph *g1 = new TGraph(dZ.size(), truez.data(), dZ.data());
    TCanvas *canvas = new TCanvas("c2","c2", 0, 0,1350,800);
    
    if (dim == "1D") {
        h -> Draw();
        h -> SetFillColor(kYellow - 7);
        h -> SetLineWidth(1);
        h -> SetFillStyle(1001);
        h -> GetXaxis() -> SetTitle("foundZ - trueZ [mm]");
        h -> GetXaxis() -> SetTitleSize(.05);
        h -> GetXaxis() -> SetLabelSize(.03);
        h -> GetXaxis() -> CenterTitle(true);
        h -> GetXaxis() -> SetNdivisions(31, 5, 0);
        h -> GetXaxis() -> SetTitleOffset(.8);
        h -> GetXaxis() -> SetRangeUser(-50, 50); // Setting x range;
        h -> GetYaxis() -> SetTitle("# of Counts");
        h -> GetYaxis() -> SetTitleSize(.05);
        h -> GetYaxis() -> SetLabelSize(.04);
        h -> GetYaxis() -> SetTitleOffset(.8);
        h -> GetYaxis() -> CenterTitle(true);

        static TLine *l1 = new TLine();
        l1 -> SetLineColor(kRed);
        l1 -> SetLineStyle(0);
        l1 -> SetLineWidth(2);
        l1 -> SetX1(0);   l1 -> SetY1(0);
        l1 -> SetX2(0);   l1 -> SetY2(h->GetMaximum());
        l1 -> Draw("same");

        static TLine *l2 = new TLine();
        l2 -> SetLineColor(kRed);
        l2 -> SetLineStyle(2);
        l2 -> SetLineWidth(4);
        l2 -> SetX1(-16);   l2 -> SetY1(0);
        l2 -> SetX2(-16);   l2 -> SetY2(h->GetMaximum());
        l2 -> Draw("same");

        static TLine *l3 = new TLine();
        l3 -> SetLineColor(kRed);
        l3 -> SetLineStyle(2);
        l3 -> SetLineWidth(4);
        l3 -> SetX1(+16);   l3 -> SetY1(0);
        l3 -> SetX2(+16);   l3 -> SetY2(h->GetMaximum());
        l3 -> Draw("same");
        TLegend *legend = new TLegend(0.5,0.8,0.9,0.9);
        // legend->SetHeader("Legend","C"); // option "C" allows to center the header
        legend->AddEntry(h,"Segments expanded between [-8, 8]mm","f");
        legend->Draw("same");
        legend->SetTextSize(0.03);

        if (select == "zscanE") {
            h -> SetTitle("Z Scan with Errors, |dEta| < 0.01 (for small N) or 0.001 (for large N), |dPhi| < 0.1");
        } else if (select == "zscan") {
            h -> SetTitle("Traditional Z Scan method, |dEta| < 0.01 (for small N) or 0.001 (for large N), |dPhi| < 0.1");
        } else if (select == "e") {
            h -> SetTitle("Nearest Z with Errors, |dPhi| < 0.1");
        } else if (select == "zscanN") {
            h -> SetTitle("Z Scan with Background Normalized");
        }
        else {
            h -> SetTitle(Form("DCA indicated by greatest bin, %lu events, |dPhi| < 0.01,  closest distance cut = %1.2f", foundz.size(), DCA_cut));
            // h -> SetTitle(Form("DCA indicated with Gaussian Fit, |dPhi| < 0.01, closest distance cut = %1.2f", DCA_cut));
        }
    }
    else if (dim == "2D") {
        g0 -> SetMarkerStyle(29);
        g0 -> SetMarkerSize(2);
        g0 -> SetMarkerColor(kBlue);
        g0 -> SetLineWidth(3);
        g0 -> SetLineColor(kWhite);
        if (select == "zscanE") {
            g0 -> SetTitle("Z Scan with Errors, |dEta| < 0.01 (for small N) or 0.001 (for large N), |dPhi| < 0.1");
        } else if (select == "zscan") {
            g0 -> SetTitle("Traditional Z Scan method, |dEta| < 0.01 (for small N) or 0.001 (for large N), |dPhi| < 0.1");
        } else if (select == "e") {
            g0 -> SetTitle("Nearest Z with Errors, |dPhi| < 0.1");
        } 
        else if (select == "zscanN") {
            g0 -> SetTitle("Z Scan with Background Normalized");
        }
        else {
            g0 -> SetTitle(Form("DCA indicated by greatest bin, %lu events, |dPhi| < 0.01, closest distance cut = %1.2f", foundz.size(), DCA_cut));
            // g0 -> SetTitle(Form("DCA indicated with Gaussian Fit, |dPhi| < 0.01, closest distance cut = %1.2f", DCA_cut));
        }
        
        gStyle -> SetTitleW(0.7);  //per cent of the pad width
        gStyle -> SetTitleH(0.08); //per cent of the pad height
        g0 -> GetXaxis() -> SetTitle("Number of particles");
        g0 -> GetXaxis() -> SetTitleSize(0.06);
        g0 -> GetXaxis() -> SetLabelSize(0.04);
        g0 -> GetXaxis() -> CenterTitle(true);
        g0 -> GetYaxis() -> SetTitle("Found_z - True_z [mm]");
        g0 -> GetYaxis() -> SetTitleSize(0.06);
        g0 -> GetYaxis() -> SetLabelSize(0.025);
        g0 -> GetYaxis() -> CenterTitle(true);
        g0 -> SetMinimum(-300); // Setting y range;
        g0 -> SetMaximum(300);  // Setting y range;
        g0 -> GetYaxis() -> SetTitleOffset(0.8); 
        g0 -> GetXaxis() -> SetTitleOffset(0.8); 
        // g0 -> GetXaxis() -> SetLimits(0, 6000); // Setting x range;
        g0 -> Draw("AP SAME");

        gPad->SetGrid(5, 2); gPad->Update();

        TLine *line1 = new TLine(0, -16, 10000, -16);
        line1->SetLineColor(kRed);
        line1->SetLineWidth(2);
        line1->Draw("same");  // Draw line on the same canvas
        TLine *line2 = new TLine(0, 16, 10000, 16);
        line2->SetLineColor(kRed);
        line2->SetLineWidth(2);
        line2->Draw("same");  // Draw line on the same canvas
        TLegend *legend = new TLegend(0.4,0.8,0.9,0.9);
        // legend->SetHeader("Legend","C"); // option "C" allows to center the header
        legend->AddEntry(g0,"Segments expanded between [-8, 8]mm","f");
        legend->SetTextSize(0.04);
        legend->Draw("same");
    }
    else if (dim == "zVSdz") {
        g1 -> SetMarkerStyle(29);
        g1 -> SetMarkerSize(1);
        g1 -> SetMarkerColor(kBlue - 7);
        g1 -> SetLineWidth(3);
        g1 -> SetLineColor(kWhite);
        if (select == "zscanE") {
            g1 -> SetTitle("Z Scan with Errors, |dEta| < 0.01 (for small N) or 0.001 (for large N), |dPhi| < 0.1");
        } else if (select == "zscan") {
            g1 -> SetTitle("Traditional Z Scan method, |dEta| < 0.01 (for small N) or 0.001 (for large N), |dPhi| < 0.1");
        } else if (select == "e") {
            g1 -> SetTitle("Nearest Z with Errors, |dPhi| < 0.1");
        } 
        else if (select == "zscanN") {
            g1 -> SetTitle("Z Scan with Background Normalized");
        }
        else if (select == "DCA") {
            g1 -> SetTitle(Form("DCA indicated by greatest bin, %lu events, |dPhi| < 0.01, closest distance cut = %1.2f", foundz.size(), DCA_cut));
        }
        else {
            g1 -> SetTitle(Form("DCA indicated with Gaussian Fit, |dPhi| < 0.01, closest distance cut = %1.2f", DCA_cut));
        }
        
        gStyle -> SetTitleW(0.7);  //per cent of the pad width
        gStyle -> SetTitleH(0.08); //per cent of the pad height
        g1 -> GetXaxis() -> SetTitle("True Z Vertex [mm]");
        g1 -> GetXaxis() -> SetTitleSize(0.06);
        g1 -> GetXaxis() -> SetLabelSize(0.04);
        g1 -> GetXaxis() -> CenterTitle(true);
        g1 -> GetYaxis() -> SetTitle("Found_z - True_z [mm]");
        g1 -> GetYaxis() -> SetTitleSize(0.06);
        g1 -> GetYaxis() -> SetLabelSize(0.025);
        g1 -> GetYaxis() -> CenterTitle(true);
        g1 -> SetMinimum(-100); // Setting y range;
        g1 -> SetMaximum(100);  // Setting y range;
        g1 -> GetYaxis() -> SetTitleOffset(0.8); 
        g1 -> GetXaxis() -> SetTitleOffset(0.8); 
        // g1 -> GetXaxis() -> SetLimits(0, 6000); // Setting x range;
        g1 -> Draw("AP SAME");

        gPad->SetGrid(5, 2); gPad->Update();

        TLine *line1 = new TLine(-500, 0, +500, 0);
        line1->SetLineColor(kRed);
        line1->SetLineWidth(2);
        line1->Draw("same");  // Draw line on the same canvas
        TLegend *legend = new TLegend(0.4,0.8,0.9,0.9);
        // legend->SetHeader("Legend","C"); // option "C" allows to center the header
        legend->AddEntry(g1,"Segments expanded between [0, 16] mm","f");
        legend->SetTextSize(0.04);
        legend->Draw("same");
    }
    else {
        double ymin = g1->GetHistogram()->GetMinimum();
        double ymax = g1->GetHistogram()->GetMaximum();
        int nbins = 1000;
        TH1F *h1 = new TH1F("ProjectionY", "Histogram of vtx_z_resolution vs vtx_z 1D projection;Found_z - True_z [mm];Entries", nbins, ymin, ymax);
        int npoints = g1->GetN();
        double x, y;
        for (int i = 0; i < npoints; i++) {
            g1->GetPoint(i, x, y);
            h1->Fill(y);
        }
        h1->Draw();
        h1 -> GetXaxis() -> CenterTitle(true);
        h1 -> GetYaxis() -> CenterTitle(true);
        h1 -> GetXaxis() -> SetRangeUser(-60, 60); // Setting x range;
        h1 -> SetFillColor(kYellow - 7);
        h1 -> SetLineWidth(1);
        h1 -> SetFillStyle(1001);
        static TLine *l1 = new TLine();
        l1 -> SetLineColor(kRed);
        l1 -> SetLineStyle(0);
        l1 -> SetLineWidth(2);
        l1 -> SetX1(0);   l1 -> SetY1(0);
        l1 -> SetX2(0);   l1 -> SetY2(h1->GetMaximum());
        l1 -> Draw("same");
        gPad -> SetLogy();
    }
}

void foundXYAnalysis (TTree *EventTree, string savePath, string select = "zscan", string dim = "2D") {
    ifstream myfile(savePath);
    if (!myfile.is_open()){
		  cout << "Unable to open linelabel" << endl;
		  system("read -n 1 -s -p \"Press any key to continue...\" echo");
		  exit(1);
 	}

    int evt, Nparticles;
    string line, substr;
    double found_x, found_y, found_z, true_x, true_y, true_z;
    std::vector<double> foundx, foundy, foundz, dX, dY, dZ;
    std::vector<double> truex, truey, truez, totalsize;

    while (getline(myfile, line)) {
        stringstream str(line);
        getline(str, substr, ',');  evt        = stoi(substr);
        getline(str, substr, ',');  // index;
        getline(str, substr, ',');  Nparticles = stoi(substr);
        getline(str, substr, ',');  found_x    = stod(substr)*1e1;
        getline(str, substr, ',');  found_y    = stod(substr)*1e1;
        getline(str, substr, ',');  found_z    = stod(substr)*1e1;
        getline(str, substr, ',');  true_x     = stod(substr)*1e1;
        getline(str, substr, ',');  true_y     = stod(substr)*1e1;
        getline(str, substr, ',');  true_z     = stod(substr)*1e1;
        
        foundx.push_back(found_x);      foundy.push_back(found_y);      foundz.push_back(found_z);
        truex.push_back(true_x);        truey.push_back(true_y);        truez.push_back(true_z);
        dX.push_back(found_x - true_x); dY.push_back(found_y - true_y); dZ.push_back(found_z - true_z);
        totalsize.push_back(Nparticles);
    }
    
    if (select == "dX") {
        TCanvas* canvas = new TCanvas("canvas", "foundX resolution as a function of z vtx", 800, 600);
        TGraph* graph = new TGraph(foundz.size(), &foundz[0], &dX[0]);
        graph->SetTitle("foundX resolution as a function of z vtx, MBD < 10");
        graph->GetXaxis()->SetTitle("foundz [mm]");
        graph->GetYaxis()->SetTitle("dX = foundX - trueX [mm]");
        graph->GetXaxis()->CenterTitle();
        graph->GetYaxis()->CenterTitle();
        graph->GetXaxis()->SetRangeUser(-260, -140);
        graph->GetYaxis()->SetRangeUser(-2., +2.);
        graph->SetMarkerStyle(20);
        graph->SetMarkerSize(0.7);
        canvas->SetGrid();
        graph->Draw("AP");
    }
    else if (select == "dY") {
        TCanvas* canvas = new TCanvas("canvas", "Scatter Plot", 800, 600);
        TGraph* graph = new TGraph(foundz.size(), &foundz[0], &dY[0]);
        graph->SetTitle("foundY resolution as a function of z vtx, MBD < 10");
        graph->GetXaxis()->SetTitle("foundz [mm]");
        graph->GetYaxis()->SetTitle("dY = foundY - trueY [mm]");
        graph->GetXaxis()->CenterTitle();
        graph->GetYaxis()->CenterTitle();
        graph->GetXaxis()->SetRangeUser(-260, -140);
        graph->GetYaxis()->SetRangeUser(-2., +2.);
        graph->SetMarkerStyle(20);
        graph->SetMarkerSize(0.7);
        canvas->SetGrid();
        graph->Draw("AP");
    }
}


void foundEtaPhiAnalysis (TTree *EventTree, string savePath, Int_t target, std::vector<std::string> method) {
    ifstream myfile(savePath);
    if (!myfile.is_open()){
		  cout << "Unable to open linelabel" << endl;
		  system("read -n 1 -s -p \"Press any key to continue...\" echo");
		  exit(1);
 	}
    int Nparticles;
    string line, substr;
    double found_x, found_y, found_z, true_x, true_y, true_z;
    vector<int> evt, idx;
    std::vector<double> foundx, foundy, foundz, dX, dY, dZ;
    std::vector<double> truex, truey, truez, totalsize;
    while (getline(myfile, line)) {
        stringstream str(line);
        getline(str, substr, ',');  int e      = stoi(substr);
        getline(str, substr, ',');  int index  = stoi(substr);
        getline(str, substr, ',');  Nparticles = stoi(substr);
        getline(str, substr, ',');  found_x    = stod(substr);
        getline(str, substr, ',');  found_y    = stod(substr);
        getline(str, substr, ',');  found_z    = stod(substr);
        getline(str, substr, ',');  true_x     = stod(substr);
        getline(str, substr, ',');  true_y     = stod(substr);
        getline(str, substr, ',');  true_z     = stod(substr);
        
        evt.push_back(e);   idx.push_back(index);
        foundx.push_back(found_x);      foundy.push_back(found_y);      foundz.push_back(found_z);
        truex.push_back(true_x);        truey.push_back(true_y);        truez.push_back(true_z);
        dX.push_back(found_x - true_x); dY.push_back(found_y - true_y); dZ.push_back(found_z - true_z);
        totalsize.push_back(Nparticles);
    }
    // Set up variables to hold the data from .root file
    std::vector<float> *ClusX = nullptr;
    std::vector<float> *ClusY = nullptr;
    std::vector<float> *ClusZ = nullptr;
    std::vector<int>   *ClusLayer = nullptr;
    std::vector<int>   *ClusLadderPhiId = nullptr;
    std::vector<int>   *UniqueAncG4P_TrackID = nullptr;
    std::vector<float> *TruthPV_x = nullptr;
    std::vector<float> *TruthPV_y = nullptr;
    std::vector<float> *TruthPV_z = nullptr;
    std::vector<float> *TruthPV_Npart = nullptr;
    int event, NTruthVtx;
    float centrality_mbd;
    EventTree -> SetBranchAddress("event", &event);
    EventTree -> SetBranchAddress("ClusX", &ClusX);
    EventTree -> SetBranchAddress("ClusY", &ClusY);
    EventTree -> SetBranchAddress("ClusZ", &ClusZ);
    EventTree -> SetBranchAddress("ClusLayer", &ClusLayer);
    EventTree -> SetBranchAddress("ClusLadderPhiId", &ClusLadderPhiId);
    EventTree -> SetBranchAddress("UniqueAncG4P_TrackID", &UniqueAncG4P_TrackID);
    EventTree -> SetBranchAddress("TruthPV_x", &TruthPV_x);
    EventTree -> SetBranchAddress("TruthPV_y", &TruthPV_y);
    EventTree -> SetBranchAddress("TruthPV_z", &TruthPV_z);
    EventTree -> SetBranchAddress("NTruthVtx", &NTruthVtx);
    EventTree -> SetBranchAddress("TruthPV_Npart", &TruthPV_Npart);
    EventTree -> SetBranchAddress("centrality_mbd", &centrality_mbd);

    std::vector<myTrackletMember> tracklet_layer_0, tracklet_layer_1;
    double r, currentZ, theta, eta, phi;   // intermediate variables;
    double dx, dy, dz, dPhi;
    if (method[1] == "single") {
        for (int i = 0; i < idx.size(); i++) {
            if (idx[i] == target) {
                target = evt[i];
                EventTree->GetEntry(idx[i]);
                for (int k = 0; k < ClusX->size(); k++) {
                    r   = std::sqrt((ClusX->at(k) - foundx[i])*(ClusX->at(k) - foundx[i]) + (ClusY->at(k) - foundy[i])*(ClusY->at(k) - foundy[i]));
                    phi = std::atan2(ClusY->at(k) - foundy[i], ClusX->at(k) - foundx[i]);
                    currentZ         = ClusZ->at(k);
                    dz               = currentZ - foundz[i];
                    theta            = std::atan2(r, dz);
                    if (dz >= 0) {
                        eta = -log(tan(theta/2));
                    } else {
                        eta = log(tan((M_PI - theta)/2));
                    }         
                    if (ClusLayer->at(k) == 3 || ClusLayer->at(k) == 4) {
                        tracklet_layer_0.push_back({ClusX->at(k), ClusY->at(k), currentZ, r,
                                                    eta, phi, 
                                                    0});
                    } else {
                        tracklet_layer_1.push_back({ClusX->at(k), ClusY->at(k), currentZ, r,
                                                    eta, phi, 
                                                    1});
                    }
                }
                break;
            }
            if (idx[i] > target)    break;
        }
        std::cout << tracklet_layer_0.size() << "," << tracklet_layer_1.size() << std::endl;
        if (method[2] == "deta") {
            dEtaCheck(target, tracklet_layer_0, tracklet_layer_1);
        }
        if (method[2] == "dphi") {
            dPhiCheck(target, tracklet_layer_0, tracklet_layer_1);
        }
    } else if (method[1] == "all") {
        TH1D *h_dEta = new TH1D("", "", 1601, -4 - .0025, 4 + .0025);
        TH2D *h_eta_phi = new TH2D("", "", 80, -4 - .05, 4 + .05, 70, -3.5 - .05, 3.5 + .05);
        int N = 1000;
        double range_min = -M_PI;
        double range_max = M_PI;
        double bin_width = (range_max - range_min) / N;
        // TH1D *h_dPhi = new TH1D("", "", 1601, -4 - .0025, 4 + .0025);
        TH1D *h_dPhi = new TH1D("", "", N, range_min, range_max);
        // TH2D *h_dPhi_Z = new TH2D("", "", N, range_min, range_max, 120, -25 - 0.5, -15. + 0.5);
        TH2D *h_dPhi_Z = new TH2D("", "", N, range_min, range_max, 7, -25.6*10., -14.4*10.);
        const int n1 = 100;         // number of bins in [0, 0.1]
        const int n2 = 10;          // number of bins in [0.1, 1.1], [1.1, 2.1], ...
        const int n3 = 8 * n2 + n1; // total number of bins
        float unequal_bins[n3 + 1]; // array of bin edges
        // Set bin edges
        for (int i = 0; i < n1; ++i) {
            unequal_bins[i] = 0.1 * static_cast<double>(i) / n1;  // divide [0, 0.1] into n1 bins
        }
        for (int i = 0; i < 8; ++i) {
            for (int j = 0; j < n2; ++j) {
                unequal_bins[n1 + i * n2 + j] = 0.1 + static_cast<double>(i) + static_cast<double>(j) / n2;  // divide [0.1 + i, 1.1 + i] into n2 bins
            }
        }
        unequal_bins[n3] = 8.1;  // set the last bin edge
        TH1D *h_ad = new TH1D("", "Angular Distance Histogram of all Events", n3, unequal_bins);

        if (method[2] == "deta") {
            for (int i = 0; i < idx.size(); i++) {
                EventTree->GetEntry(idx[i]);
                std::cout << idx[i] << std::endl;
                for (int k = 0; k < ClusX->size(); k++) {
                    r   = std::sqrt((ClusX->at(k) - foundx[i])*(ClusX->at(k) - foundx[i]) + (ClusY->at(k) - foundy[i])*(ClusY->at(k) - foundy[i]));
                    phi = std::atan2(ClusY->at(k) - foundy[i], ClusX->at(k) - foundx[i]);
                    currentZ         = ClusZ->at(k);
                    dz               = currentZ - foundz[i];
                    theta            = std::atan2(r, dz);
                    if (dz >= 0) {
                        eta = -log(tan(theta/2));
                    } else {
                        eta = log(tan((M_PI - theta)/2));
                    }         
                    if (ClusLayer->at(k) == 3 || ClusLayer->at(k) == 4) {
                        tracklet_layer_0.push_back({ClusX->at(k), ClusY->at(k), currentZ, r,
                                                    eta, phi, 
                                                    0});
                    } else {
                        tracklet_layer_1.push_back({ClusX->at(k), ClusY->at(k), currentZ, r,
                                                    eta, phi, 
                                                    1});
                    }
                }
                for (int i = 0; i < tracklet_layer_0.size(); i++) {
                    for (int j = 0; j < tracklet_layer_1.size(); j++) {
                        // dEta.push_back(tracklet_layer_0[i].eta - tracklet_layer_1[j].eta);
                        // dPhi.push_back(tracklet_layer_0[i].phi - tracklet_layer_1[j].phi);
                        double dEta = tracklet_layer_0[i].eta - tracklet_layer_1[j].eta;
                        // double dPhi = tracklet_layer_0[i].phi - tracklet_layer_1[j].phi;
                        // double dR   = std::sqrt(dEta*dEta + dPhi*dPhi);
                        h_dEta    -> Fill(dEta);
                        // h_dPhi    -> Fill(dPhi);
                        // h_eta_phi -> Fill(dEta, dPhi);
                        // h_ad      -> Fill(dR);
                    }
                }
                tracklet_layer_0.clear();   tracklet_layer_1.clear();
            }
            dEtaCheckAll(h_dEta);
        }
        else if (method[2] == "dphi") {
            for (int i = 0; i < idx.size(); i++) {
                EventTree->GetEntry(idx[i]);
                // std::cout << idx[i] << "," << event << "," << evt[i] << std::endl;
                if (event != evt[i])    break;
                for (int k = 0; k < ClusX->size(); k++) {
                    r        = std::sqrt((ClusX->at(k) - foundx[i])*(ClusX->at(k) - foundx[i]) + (ClusY->at(k) - foundy[i])*(ClusY->at(k) - foundy[i]));
                    phi      = std::atan2(ClusY->at(k) - foundy[i], ClusX->at(k) - foundx[i]);
                    currentZ = ClusZ->at(k);
                    dz       = currentZ - foundz[i];
                    theta    = std::atan2(r, dz);
                    if (dz >= 0) {
                        eta = -log(tan(theta/2));
                    } else {
                        eta = log(tan((M_PI - theta)/2));
                    }         
                    if (ClusLayer->at(k) == 3 || ClusLayer->at(k) == 4) {
                        tracklet_layer_0.push_back({ClusX->at(k), ClusY->at(k), currentZ, r,
                                                    eta, phi, 
                                                    0});
                    } else {
                        tracklet_layer_1.push_back({ClusX->at(k), ClusY->at(k), currentZ, r,
                                                    eta, phi, 
                                                    1});
                    }
                }
                for (int i = 0; i < tracklet_layer_0.size(); i++) {
                    for (int j = 0; j < tracklet_layer_1.size(); j++) {
                        double dPhi = tracklet_layer_0[i].phi - tracklet_layer_1[j].phi;
                        if (dPhi > M_PI) {
                            h_dPhi -> Fill(dPhi - 2*M_PI);
                        }
                        else if (dPhi < -M_PI) {
                            h_dPhi -> Fill(dPhi + 2*M_PI);
                        }
                        else {
                            h_dPhi -> Fill(dPhi);
                        }
                        // h_dPhi    -> Fill(dPhi);
                    }
                }
                tracklet_layer_0.clear();   tracklet_layer_1.clear();
            }
            dPhiCheckAll(h_dPhi, method, target);
        }
        else if (method[2] == "dR") {
            for (int i = 0; i < idx.size(); i++) {
                EventTree->GetEntry(idx[i]);
                std::cout << idx[i] << std::endl;
                for (int k = 0; k < ClusX->size(); k++) {
                    r   = std::sqrt((ClusX->at(k) - foundx[i])*(ClusX->at(k) - foundx[i]) + (ClusY->at(k) - foundy[i])*(ClusY->at(k) - foundy[i]));
                    phi = std::atan2(ClusY->at(k) - foundy[i], ClusX->at(k) - foundx[i]);
                    currentZ         = ClusZ->at(k);
                    dz               = currentZ - foundz[i];
                    theta            = std::atan2(r, dz);
                    if (dz >= 0) {
                        eta = -log(tan(theta/2));
                    } else {
                        eta = log(tan((M_PI - theta)/2));
                    }         
                    if (ClusLayer->at(k) == 3 || ClusLayer->at(k) == 4) {
                        tracklet_layer_0.push_back({ClusX->at(k), ClusY->at(k), currentZ, r,
                                                    eta, phi, 
                                                    0});
                    } else {
                        tracklet_layer_1.push_back({ClusX->at(k), ClusY->at(k), currentZ, r,
                                                    eta, phi, 
                                                    1});
                    }
                }
                for (int i = 0; i < tracklet_layer_0.size(); i++) {
                    for (int j = 0; j < tracklet_layer_1.size(); j++) {
                        // dEta.push_back(tracklet_layer_0[i].eta - tracklet_layer_1[j].eta);
                        // dPhi.push_back(tracklet_layer_0[i].phi - tracklet_layer_1[j].phi);
                        double dEta = tracklet_layer_0[i].eta - tracklet_layer_1[j].eta;
                        double dPhi = tracklet_layer_0[i].phi - tracklet_layer_1[j].phi;
                        double dR   = std::sqrt(dEta*dEta + dPhi*dPhi);
                        // h_dEta    -> Fill(dEta);
                        // h_dPhi    -> Fill(dPhi);
                        // h_eta_phi -> Fill(dEta, dPhi);
                        h_ad      -> Fill(dR);
                    }
                }
                tracklet_layer_0.clear();   tracklet_layer_1.clear();
            }
            dRCheckAll(h_ad);
        }
        else if (method[2] == "dphiZ") {
            for (int i = 0; i < idx.size(); i++) {
                EventTree->GetEntry(idx[i]);
                // std::cout << idx[i] << "," << event << "," << evt[i] << std::endl;
                // if (event != evt[i])    break;
                double fx = foundx[i];
                double fy = foundy[i];
                double fz = foundz[i];
                for (int k = 0; k < ClusX->size(); k++) {
                    dx       = ClusX->at(k) - fx;
                    dy       = ClusY->at(k) - fy;
                    r        = std::sqrt(dx * dx + dy * dy);
                    phi      = std::atan2(dy, dx);
                    // r        = std::sqrt((ClusX->at(k) - fx)*(ClusX->at(k) - fx) + (ClusY->at(k) - fy)*(ClusY->at(k) - fy));
                    // phi      = std::atan2(ClusY->at(k) - foundy[i], ClusX->at(k) - foundx[i]);
                    currentZ = ClusZ->at(k);
                    dz       = currentZ - fz;
                    theta    = std::atan2(r, dz);
                    if (dz >= 0) {
                        eta = -log(tan(theta/2));
                    } else {
                        eta = log(tan((M_PI - theta)/2));
                    }         
                    if (ClusLayer->at(k) == 3 || ClusLayer->at(k) == 4) {
                        tracklet_layer_0.push_back({ClusX->at(k), ClusY->at(k), currentZ, r,
                                                    eta, phi, 
                                                    0});
                    } else {
                        tracklet_layer_1.push_back({ClusX->at(k), ClusY->at(k), currentZ, r,
                                                    eta, phi, 
                                                    1});
                    }
                }
                for (int i = 0; i < tracklet_layer_0.size(); i++) {
                    for (int j = 0; j < tracklet_layer_1.size(); j++) {
                        dPhi = tracklet_layer_0[i].phi - tracklet_layer_1[j].phi;
                        if (dPhi > M_PI) {
                            h_dPhi_Z -> Fill(dPhi - 2*M_PI, fz*10);
                        }
                        else if (dPhi < -M_PI) {
                            h_dPhi_Z -> Fill(dPhi + 2*M_PI, fz*10);
                        }
                        else {
                            h_dPhi_Z -> Fill(dPhi, fz*10);
                        }
                        // h_dPhi    -> Fill(dPhi);
                    }
                }
                tracklet_layer_0.clear();   tracklet_layer_1.clear();
            }
            // dPhiCheckAll(h_dPhi_Z);
            h_dPhi_Z -> Draw("lego2");
            h_dPhi_Z -> SetTitle("All dPhi values (no event mix) with 10 bins of Z vtx");
            h_dPhi_Z -> GetXaxis() -> SetTitle("dPhi"); h_dPhi_Z -> GetXaxis() -> CenterTitle(true);
            h_dPhi_Z -> GetYaxis() -> SetTitle("found z vtx position [mm]");   h_dPhi_Z -> GetYaxis() -> CenterTitle(true);
            h_dPhi_Z -> GetXaxis() -> SetTitleOffset(1.8);   h_dPhi_Z -> GetYaxis() -> SetTitleOffset(1.8);
        }
    }
}

void dPhiInZVtx (TTree *EventTree, string savePath, Int_t target, std::vector<std::string> method) {
    double upperRange = stod(method[2])/10.;
    double lowerRange = stod(method[1])/10.;
    ifstream myfile(savePath);
    if (!myfile.is_open()){
		  cout << "Unable to open linelabel" << endl;
		  system("read -n 1 -s -p \"Press any key to continue...\" echo");
		  exit(1);
 	}
    int Nparticles;
    string line, substr;
    double found_x, found_y, found_z, true_x, true_y, true_z;
    std::vector<int> evt, idx;
    std::vector<double> foundx, foundy, foundz, dX, dY, dZ;
    std::vector<double> truex, truey, truez, totalsize;
    std::vector<int> selected_evt, selected_idx;
    std::ofstream outFile("test.txt");
    std::ofstream outFile2("test2.txt");
    while (getline(myfile, line)) {
        stringstream str(line);
        getline(str, substr, ',');  int e      = stoi(substr);
        getline(str, substr, ',');  int index  = stoi(substr);
        getline(str, substr, ',');  Nparticles = stoi(substr);
        getline(str, substr, ',');  found_x    = stod(substr);
        getline(str, substr, ',');  found_y    = stod(substr);
        getline(str, substr, ',');  found_z    = stod(substr);
        getline(str, substr, ',');  true_x     = stod(substr);
        getline(str, substr, ',');  true_y     = stod(substr);
        getline(str, substr, ',');  true_z     = stod(substr);
        
        // if (found_z >= -lowerRange && found_z <= -upperRange) {
        if (found_z >= -20.9375 && found_z <= -20.75) {
            evt.push_back(e);   idx.push_back(index);
            foundx.push_back(found_x);      foundy.push_back(found_y);      foundz.push_back(found_z);
            truex.push_back(true_x);        truey.push_back(true_y);        truez.push_back(true_z);
            dX.push_back(found_x - true_x); dY.push_back(found_y - true_y); dZ.push_back(found_z - true_z);
            totalsize.push_back(Nparticles);
        }
    }

    // Set up variables to hold the data from .root file
    std::vector<float> *ClusX = nullptr;
    std::vector<float> *ClusY = nullptr;
    std::vector<float> *ClusZ = nullptr;
    std::vector<int>   *ClusLayer = nullptr;
    std::vector<int>   *ClusLadderPhiId = nullptr;
    std::vector<int>   *UniqueAncG4P_TrackID = nullptr;
    std::vector<float> *TruthPV_x = nullptr;
    std::vector<float> *TruthPV_y = nullptr;
    std::vector<float> *TruthPV_z = nullptr;
    std::vector<float> *TruthPV_Npart = nullptr;
    int event, NTruthVtx;
    float centrality_mbd;
    EventTree -> SetBranchAddress("event", &event);
    EventTree -> SetBranchAddress("ClusX", &ClusX);
    EventTree -> SetBranchAddress("ClusY", &ClusY);
    EventTree -> SetBranchAddress("ClusZ", &ClusZ);
    EventTree -> SetBranchAddress("ClusLayer", &ClusLayer);
    EventTree -> SetBranchAddress("ClusLadderPhiId", &ClusLadderPhiId);
    EventTree -> SetBranchAddress("UniqueAncG4P_TrackID", &UniqueAncG4P_TrackID);
    EventTree -> SetBranchAddress("TruthPV_x", &TruthPV_x);
    EventTree -> SetBranchAddress("TruthPV_y", &TruthPV_y);
    EventTree -> SetBranchAddress("TruthPV_z", &TruthPV_z);
    EventTree -> SetBranchAddress("NTruthVtx", &NTruthVtx);
    EventTree -> SetBranchAddress("TruthPV_Npart", &TruthPV_Npart);
    EventTree -> SetBranchAddress("centrality_mbd", &centrality_mbd);
    // setup ends

    int N = 1000;
    double range_min = -M_PI;
    double range_max = M_PI;
    double bin_width = (range_max - range_min) / N;
    // TH1D *h_dPhi = new TH1D("", "", 1601, -4 - .0025, 4 + .0025);
    TH1D *h_dPhi = new TH1D("", "", N, range_min, range_max);
    TH1D *h_Phi  = new TH1D("", "", N, range_min, range_max);
    double r, currentZ, theta, eta, phi;   // intermediate variables;
    double dx, dy, dz, dPhi;
    std::vector<std::vector<double>> all_Phi_0, all_Phi_1;
    for (int i = 0; i < idx.size(); i++) {
        EventTree->GetEntry(idx[i]);
        // std::cout << idx[i] << "," << event << "," << evt[i] << std::endl;
        // if (event != evt[i])    break;
        double fx = foundx[i];
        double fy = foundy[i];
        double fz = foundz[i];
        all_Phi_0.push_back(std::vector<double>());   all_Phi_1.push_back(std::vector<double>());
        for (int k = 0; k < ClusX->size(); k++) {
            dx       = ClusX->at(k) - fx;
            dy       = ClusY->at(k) - fy;
            phi      = std::atan2(dy, dx);
            // currentZ = ClusZ->at(k);
            // dz       = currentZ - fz;
            // r        = std::sqrt(dx * dx + dy * dy);
            // theta    = std::atan2(r, dz);
            // if (dz >= 0) {
            //     eta = -log(tan(theta/2));
            // } else {
            //     eta = log(tan((M_PI - theta)/2));
            // }         
            if (ClusLayer->at(k) == 3 || ClusLayer->at(k) == 4) {
                all_Phi_0[i].push_back(phi);
            } else {
                all_Phi_1[i].push_back(phi);
            }
        }
        // for (int i = 0; i < tracklet_layer_0.size(); i++) {
        //     for (int j = 0; j < tracklet_layer_1.size(); j++) {
        //         dPhi = tracklet_layer_0[i].phi - tracklet_layer_1[j].phi;
        //         if (dPhi > M_PI) {
        //             h_dPhi -> Fill(dPhi - 2*M_PI);
        //         }
        //         else if (dPhi < -M_PI) {
        //             h_dPhi -> Fill(dPhi + 2*M_PI);
        //         }
        //         else {
        //             h_dPhi -> Fill(dPhi);
        //         }
        //     }
        // }
        // tracklet_layer_0.clear();   tracklet_layer_1.clear();
    }
    std::cout << all_Phi_0.size() << "," << all_Phi_1.size() << std::endl;
    std::cout << all_Phi_0[0].size() << "," << all_Phi_1[0].size() << std::endl;
    // /*
    for (int i = 0; i < all_Phi_0.size(); i++) {
        if (i >= target)    break;
        // std::cout << i << ". Event: " << evt[i] << ", index: " << idx[i] << std::endl;
        for (int j = i; j < all_Phi_1.size(); j++) {
            if (j >= target)    break;
            // std::cout << " ---- interacting with Event: " << evt[j] << ", index: " << idx[j] << std::endl;
            for (int k = 0; k < all_Phi_0[i].size(); k++) {
                for (int l = 0; l < all_Phi_1[j].size(); l++) {
                    dPhi = all_Phi_0[i][k] - all_Phi_1[j][l];
                    if (dPhi > M_PI) {
                        h_dPhi -> Fill(dPhi - 2*M_PI);
                    }
                    else if (dPhi < -M_PI) {
                        h_dPhi -> Fill(dPhi + 2*M_PI);
                    }
                    else {
                        h_dPhi -> Fill(dPhi);
                    }
                    // h_dPhi -> Fill(dPhi);
                }
            }
        }
    }
    // */
    // dPhiCheckAll(h_dPhi, method, target);
    for (int n = 0; n < 20; n++) {
    for (int i = 0; i < all_Phi_0.size(); i++) {
        if (i >= target)    break;
            for (int k = 0; k < all_Phi_0[i].size(); k++) {
                for (int l = 0; l < all_Phi_1[i].size(); l++) {
                    dPhi = all_Phi_0[i][k] - all_Phi_1[i][l];
                    if (dPhi > M_PI) {
                        h_Phi -> Fill(dPhi - 2*M_PI);
                    }
                    else if (dPhi < -M_PI) {
                        h_Phi -> Fill(dPhi + 2*M_PI);
                    }
                    else {
                        h_Phi -> Fill(dPhi);
                    }
                    // h_dPhi -> Fill(dPhi);
                }
            }
    }
    }
    // for (int k = 0; k < 1e4; k++) {
    // for (int i = 0; i < all_Phi_0.size(); i++) {
    //     if (i >= target) break;
    //     for (int j = 0; j < all_Phi_0[i].size(); j++) {
    //         h_Phi -> Fill(all_Phi_0[i][j]);
    //     }
    // }
    // for (int i = 0; i < all_Phi_1.size(); i++) {
    //     if (i >= target) break;
    //     for (int j = 0; j < all_Phi_1[i].size(); j++) {
    //         h_Phi -> Fill(all_Phi_1[i][j]);
    //     }
    // }
    // }
    // dPhiCheckAll(h_Phi, method, target);
    dPhiCheckDouble(h_dPhi, h_Phi, method, target);
    /*
    TH1D *h_foundz = new TH1D("", "", 161, (-25. - 0.03125)*10, (-15. + 0.03125)*10);
    for (int i = 0; i < foundz.size(); i++) {
        h_foundz -> Fill(foundz[i]*10);
    }
    h_foundz -> Draw();
    h_foundz -> SetFillStyle(1001);
    h_foundz -> SetFillColor(kYellow - 7);
    h_foundz -> SetLineWidth(1);
    h_foundz -> GetXaxis() -> SetTitle("found z vtx position [mm]");
    h_foundz -> GetXaxis() -> CenterTitle(true);


    int binIndex = 1; // Example bin index
    double binWidth = h_foundz->GetBinWidth(binIndex);
    std::cout << "Width of bin " << binIndex << " is " << binWidth << std::endl;

    for (int i = 1; i <= h_foundz->GetNbinsX(); i++) {
        if (h_foundz->GetBinContent(i) == 0) {
            // Get the center of the bin
            double binCenter = h_foundz->GetBinCenter(i);
            std::cout << "Empty bin at x = " << std::setprecision(8) << binCenter/10. << std::endl;
        }
    }
    */
}

// void dPhiInZVtx (TTree *EventTree, string savePath, Int_t target, std::vector<std::string> method) {
//     ifstream myfile(savePath);
//     if (!myfile.is_open()){
// 		  cout << "Unable to open linelabel" << endl;
// 		  system("read -n 1 -s -p \"Press any key to continue...\" echo");
// 		  exit(1);
//  	}
//     int Nparticles;
//     string line, substr;
//     double found_x, found_y, found_z, true_x, true_y, true_z;
//     std::vector<int> evt, idx;
//     std::vector<double> foundx, foundy, foundz, dX, dY, dZ;
//     std::vector<double> truex, truey, truez, totalsize;
//     std::vector<int> selected_evt, selected_idx;
//     std::ofstream outFile("test.txt");
//     std::ofstream outFile2("test2.txt");
//     while (getline(myfile, line)) {
//         stringstream str(line);
//         getline(str, substr, ',');  int e      = stoi(substr);
//         getline(str, substr, ',');  int index  = stoi(substr);
//         getline(str, substr, ',');  Nparticles = stoi(substr);
//         getline(str, substr, ',');  found_x    = stod(substr);
//         getline(str, substr, ',');  found_y    = stod(substr);
//         getline(str, substr, ',');  found_z    = stod(substr);
//         getline(str, substr, ',');  true_x     = stod(substr);
//         getline(str, substr, ',');  true_y     = stod(substr);
//         getline(str, substr, ',');  true_z     = stod(substr);
        
//         evt.push_back(e);   idx.push_back(index);
//         foundx.push_back(found_x);      foundy.push_back(found_y);      foundz.push_back(found_z);
//         truex.push_back(true_x);        truey.push_back(true_y);        truez.push_back(true_z);
//         dX.push_back(found_x - true_x); dY.push_back(found_y - true_y); dZ.push_back(found_z - true_z);
//         totalsize.push_back(Nparticles);
//         if (found_z >= -21.4375 && found_z <= -20.5625) {
//             outFile << e << "," << index << std::endl;
//             selected_evt.push_back(e);
//             selected_idx.push_back(index);
//         }
//     }

//     // Set up variables to hold the data from .root file
//     std::vector<float> *ClusX = nullptr;
//     std::vector<float> *ClusY = nullptr;
//     std::vector<float> *ClusZ = nullptr;
//     std::vector<int>   *ClusLayer = nullptr;
//     std::vector<int>   *ClusLadderPhiId = nullptr;
//     std::vector<int>   *UniqueAncG4P_TrackID = nullptr;
//     std::vector<float> *TruthPV_x = nullptr;
//     std::vector<float> *TruthPV_y = nullptr;
//     std::vector<float> *TruthPV_z = nullptr;
//     std::vector<float> *TruthPV_Npart = nullptr;
//     int event, NTruthVtx;
//     float centrality_mbd;
//     EventTree -> SetBranchAddress("event", &event);
//     EventTree -> SetBranchAddress("ClusX", &ClusX);
//     EventTree -> SetBranchAddress("ClusY", &ClusY);
//     EventTree -> SetBranchAddress("ClusZ", &ClusZ);
//     EventTree -> SetBranchAddress("ClusLayer", &ClusLayer);
//     EventTree -> SetBranchAddress("ClusLadderPhiId", &ClusLadderPhiId);
//     EventTree -> SetBranchAddress("UniqueAncG4P_TrackID", &UniqueAncG4P_TrackID);
//     EventTree -> SetBranchAddress("TruthPV_x", &TruthPV_x);
//     EventTree -> SetBranchAddress("TruthPV_y", &TruthPV_y);
//     EventTree -> SetBranchAddress("TruthPV_z", &TruthPV_z);
//     EventTree -> SetBranchAddress("NTruthVtx", &NTruthVtx);
//     EventTree -> SetBranchAddress("TruthPV_Npart", &TruthPV_Npart);
//     EventTree -> SetBranchAddress("centrality_mbd", &centrality_mbd);
//     // setup ends

//     std::vector<std::vector<double>> all_Phi;
//     for (int i = 0; i < selected_evt.size(); i++) {
//         EventTree -> GetEntry(selected_idx[i]);
//         // std::cout << event << "," << selected_idx[i] << std::endl;
//         all_Phi.push_back(std::vector<double>());
//         for (int j = 0; j < ClusX->size(); j++) {
//             all_Phi[i].push_back(ClusX->at(j));
//         }
//         // for (int j = 0; j < selected_evt.size(); j++) {
//         //     EventTree2 -> GetEntry(selected_idx[j]);
//         //     // std::cout << event2 << "," << selected_idx[j] << "              ";
//         // }
//     }
//     // int size_i, size_j;
//     // TH1F *h = new TH1F("h", "Histogram of Phi Differences", 100, -10., 10.);
//     // for (int i = 0; i < all_Phi.size(); i++) {
//     //     // std::cout << i << std::endl;
//     //     for (int j = i; j < all_Phi.size(); j++) {
            
//     //         size_i = all_Phi[i].size();
//     //         size_j = all_Phi[j].size();
//     //         for (int k = 0; k < size_i; k++) {
//     //             for (int l = 0; l < size_j; l++) {
//     //                 h -> Fill(all_Phi[i][k] - all_Phi[j][l]);
//     //             }
//     //         }
//     //     }
//     // }
//     #include <vector>

// std::vector<float> flattened_Phi;
// std::vector<int> indices;

// int index = 0;
// indices.push_back(index);

// // Flatten the vector and store the starting indices of each original sub-vector
// for (const auto& subVec : all_Phi) {
//     flattened_Phi.insert(flattened_Phi.end(), subVec.begin(), subVec.end());
//     index += subVec.size();
//     indices.push_back(index); // This marks the end of a sub-vector
// }
// TH1F *h = new TH1F("h", "Histogram of Phi Differences", 100, -10., 10.);

// // Adjust the loop to use the flattened structure
// for (int i = 0; i < all_Phi.size(); i++) {
//     for (int j = i; j < all_Phi.size(); j++) {
//         std::cout << i << "," << j << std::endl;
//         int start_i = indices[i];
//         int end_i = indices[i+1];
//         int start_j = indices[j];
//         int end_j = indices[j+1];

//         for (int k = start_i; k < end_i; k++) {
//             for (int l = start_j; l < end_j; l++) {
//                 h->Fill(flattened_Phi[k] - flattened_Phi[l]);
//             }
//         }
//     }
// }


// }