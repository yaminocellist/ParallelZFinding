#include "globalDefinitions.h"

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

void dEtaCheck (const int &evt, std::vector<myTrackletMember> t0, std::vector<myTrackletMember> t1,
                const double &eta_cut_low, const double &eta_cut_high,
                const double &phi_cut_low, const double &phi_cut_high, const double &trueZ) {
    TH1D *h_dEta = new TH1D("", "", 160, -4 - .025, 4 + .025);

    for (int i = 0; i < t0.size(); i++) {
        for (int j = 0; j < t1.size(); j++) {
            h_dEta -> Fill(t0[i].eta - t1[j].eta);
        }
    }

    TCanvas *can1 = new TCanvas("c1","c1",0,50,1800,550);
    h_dEta -> Draw();
    h_dEta -> SetFillColor(kYellow - 7);
    h_dEta -> SetLineWidth(1);
    h_dEta -> SetFillStyle(1001);
    h_dEta -> GetXaxis() -> SetTitle("Eta value");
    h_dEta -> GetXaxis() -> SetTitleSize(.05);
    h_dEta -> GetXaxis() -> SetLabelSize(.04);
    h_dEta -> GetXaxis() -> CenterTitle(true);
    h_dEta -> GetXaxis() -> SetNdivisions(31, 5, 0);
    h_dEta -> GetXaxis() -> SetTitleOffset(.8);
    h_dEta -> GetYaxis() -> SetTitle("# of Counts");
    h_dEta -> GetYaxis() -> SetTitleSize(.05);
    h_dEta -> GetYaxis() -> SetLabelSize(.04);
    h_dEta -> GetYaxis() -> SetTitleOffset(.62);
    h_dEta -> GetYaxis() -> CenterTitle(true);
    h_dEta -> GetXaxis() -> SetRangeUser(-0.9, 1.1); // Setting x range;
    static TLine *l1 = new TLine();
    l1 -> SetLineColor(kRed);
    l1 -> SetLineStyle(0);
    l1 -> SetLineWidth(2);
    l1 -> SetX1(0.1);   l1 -> SetY1(0);
    l1 -> SetX2(0.1);   l1 -> SetY2(h_dEta->GetMaximum());
    l1 -> Draw("same");

    h_dEta -> SetTitle(Form("dEta data of Event %d", evt));
    // gPad -> SetLogy();
    can1 -> SaveAs(Form("../External/zFindingPlots/dEta_single_%d.png", evt));
}

void dPhiCheck (const int &evt, std::vector<myTrackletMember> t0, std::vector<myTrackletMember> t1) {

    TH1D *h_dphi = new TH1D("", "", 2801, -7 - .0025, 7 + .0025);

    for (int i = 0; i < t0.size(); i++) {
        for (int j = 0; j < t1.size(); j++) {
            // if (t0[i].trackID == t1[i].trackID)
                h_dphi -> Fill(t0[i].phi - t1[j].phi);
        }
    }

    TCanvas *can1 = new TCanvas("c1","c1",0,50,1800,550);
    h_dphi -> Draw();
    h_dphi -> SetFillColor(kYellow - 7);
    h_dphi -> SetLineWidth(1);
    h_dphi -> SetFillStyle(1001);
    h_dphi -> GetXaxis() -> SetTitle("Phi value");
    h_dphi -> GetXaxis() -> SetTitleSize(.05);
    h_dphi -> GetXaxis() -> SetLabelSize(.04);
    h_dphi -> GetXaxis() -> CenterTitle(true);
    h_dphi -> GetXaxis() -> SetNdivisions(31, 5, 0);
    h_dphi -> GetXaxis() -> SetTitleOffset(.8);
    h_dphi -> GetYaxis() -> SetTitle("# of Counts");
    h_dphi -> GetYaxis() -> SetTitleSize(.05);
    h_dphi -> GetYaxis() -> SetLabelSize(.04);
    h_dphi -> GetYaxis() -> SetTitleOffset(.62);
    h_dphi -> GetYaxis() -> CenterTitle(true);
    h_dphi -> GetXaxis() -> SetRangeUser(-0.1, 0.1); // Setting x range;
    h_dphi -> SetTitle(Form("dPhi data of Event %d", evt));
    gPad -> SetLogy();
    can1 -> SaveAs(Form("../External/zFindingPlots/dPhi_single_%d.png", evt));
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
    getline(myfile, line);

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
        TCanvas* canvas = new TCanvas("canvas", "Scatter Plot", 800, 600);
        TGraph* graph = new TGraph(foundz.size(), &foundz[0], &dX[0]);
        graph->SetTitle("Scatter Plot");
        graph->GetXaxis()->SetTitle("foundz");
        graph->GetYaxis()->SetTitle("dX");
        graph->GetXaxis()->SetRangeUser(-260, -140);
        graph->SetMarkerStyle(20);
        graph->SetMarkerSize(1);
        graph->Draw("AP");
    }
    else if (select == "dY") {
        TCanvas* canvas = new TCanvas("canvas", "Scatter Plot", 800, 600);
        TGraph* graph = new TGraph(foundz.size(), &foundz[0], &dY[0]);
        graph->SetTitle("Scatter Plot");
        graph->GetXaxis()->SetTitle("foundz");
        graph->GetYaxis()->SetTitle("dX");
        graph->GetXaxis()->SetRangeUser(-260, -140);
        graph->SetMarkerStyle(20);
        graph->SetMarkerSize(1);
        graph->Draw("AP");
    }
}