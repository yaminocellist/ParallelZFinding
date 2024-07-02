#include "headers/Analysis.h"

void INTTZAnalysisLite (string filePath, std::vector<string> options) {
    ifstream myfile(filePath);
    if (!myfile.is_open()){
		std::cout << "Unable to open linelabel" << std::endl;
		system("read -n 1 -s -p \"Press any key to continue...\" echo");
		exit(1);
 	}

    string line, value;
    std::vector<int> index, event, NClus;
    std::vector<double> foundZ, MBD_z_vtx, MBD_centrality;

    getline(myfile, line);
    while (getline(myfile, line)) {
        stringstream data(line);
        getline(data, value, ',');  int i     = std::stoi(value);
        getline(data, value, ',');  int e     = std::stoi(value);
        getline(data, value, ',');  int N     = std::stoi(value);
        getline(data, value, ',');  double f  = std::stod(value)*10;
        getline(data, value, ',');  double Mz = std::stod(value)*10;
        getline(data, value, ',');  double Mc = std::stod(value);

        if (f >= -250. && f <= -50. && Mc <= 0.7) {
            index.push_back(i);     event.push_back(e);         NClus.push_back(N);
            foundZ.push_back(f);    MBD_z_vtx.push_back(Mz);    MBD_centrality.push_back(Mc);
        }
    }
    
    TH1D *h_resolution = new TH1D("found z resolution", Form("INTT found z resolution of %lu events;dZ [mm];# of counts", index.size()), 151, -1e2, 1e2);
    TGraph *g_NFunction = new TGraph();
    TGraph *g_ZFunction = new TGraph();
    TGraph *g_CFunction = new TGraph();
    double z_resolution;
    for (int i = 0; i < index.size(); i++) {
        z_resolution = foundZ[i] - MBD_z_vtx[i];
        h_resolution -> Fill(z_resolution);
        g_NFunction  -> SetPoint(g_NFunction->GetN(), NClus[i], z_resolution);
        g_ZFunction  -> SetPoint(g_ZFunction->GetN(), MBD_z_vtx[i], z_resolution);
        g_CFunction  -> SetPoint(g_CFunction->GetN(), MBD_centrality[i], z_resolution);
    }
    ZResolutionSinglePlot(h_resolution, options, "foundZResolution_1D");
    TGraphSinglePlot(g_NFunction, Form("INTT found z resolution of %d events", g_NFunction->GetN()), "# of Hits", "foundZResolution_vs_NClus");
    TGraphSinglePlot(g_ZFunction, Form("INTT found z resolution of %d events", g_ZFunction->GetN()), "MBD z vtx [mm]", "foundZResolution_vs_MBDZ");
    TGraphSinglePlot(g_CFunction, Form("INTT found z resolution of %d events", g_CFunction->GetN()), "MBD centrality", "foundZResolution_vs_MBC");

    // TCanvas *canvas = new TCanvas("c2","c2", 0, 0,1350,800);
    
    // if (dim == "1D") {
    //     h -> Draw();
    //     h -> SetFillColor(kYellow - 7);
    //     h -> SetLineWidth(1);
    //     h -> SetFillStyle(1001);
    //     h -> GetXaxis() -> SetTitle("foundZ - trueZ [mm]");
    //     h -> GetXaxis() -> SetTitleSize(.05);
    //     h -> GetXaxis() -> SetLabelSize(.03);
    //     h -> GetXaxis() -> CenterTitle(true);
    //     h -> GetXaxis() -> SetNdivisions(31, 5, 0);
    //     h -> GetXaxis() -> SetTitleOffset(.8);
    //     h -> GetXaxis() -> SetRangeUser(-50, 50); // Setting x range;
    //     h -> GetYaxis() -> SetTitle("# of Counts");
    //     h -> GetYaxis() -> SetTitleSize(.05);
    //     h -> GetYaxis() -> SetLabelSize(.04);
    //     h -> GetYaxis() -> SetTitleOffset(.8);
    //     h -> GetYaxis() -> CenterTitle(true);

    //     static TLine *l1 = new TLine();
    //     l1 -> SetLineColor(kRed);
    //     l1 -> SetLineStyle(0);
    //     l1 -> SetLineWidth(2);
    //     l1 -> SetX1(0);   l1 -> SetY1(0);
    //     l1 -> SetX2(0);   l1 -> SetY2(h->GetMaximum());
    //     l1 -> Draw("same");

    //     static TLine *l2 = new TLine();
    //     l2 -> SetLineColor(kRed);
    //     l2 -> SetLineStyle(2);
    //     l2 -> SetLineWidth(4);
    //     l2 -> SetX1(-16);   l2 -> SetY1(0);
    //     l2 -> SetX2(-16);   l2 -> SetY2(h->GetMaximum());
    //     l2 -> Draw("same");

    //     static TLine *l3 = new TLine();
    //     l3 -> SetLineColor(kRed);
    //     l3 -> SetLineStyle(2);
    //     l3 -> SetLineWidth(4);
    //     l3 -> SetX1(+16);   l3 -> SetY1(0);
    //     l3 -> SetX2(+16);   l3 -> SetY2(h->GetMaximum());
    //     l3 -> Draw("same");
    //     TLegend *legend = new TLegend(0.5,0.8,0.9,0.9);
    //     // legend->SetHeader("Legend","C"); // option "C" allows to center the header
    //     legend->AddEntry(h,"Segments expanded between [-8, 8]mm","f");
    //     legend->Draw("same");
    //     legend->SetTextSize(0.03);

    //     if (select == "zscanE") {
    //         h -> SetTitle("Z Scan with Errors, |dEta| < 0.01 (for small N) or 0.001 (for large N), |dPhi| < 0.1");
    //     } else if (select == "zscan") {
    //         h -> SetTitle("Traditional Z Scan method, |dEta| < 0.01 (for small N) or 0.001 (for large N), |dPhi| < 0.1");
    //     } else if (select == "e") {
    //         h -> SetTitle("Nearest Z with Errors, |dPhi| < 0.1");
    //     } else if (select == "zscanN") {
    //         h -> SetTitle("Z Scan with Background Normalized");
    //     }
    //     else {
    //         h -> SetTitle(Form("DCA indicated by greatest bin, %lu events, |dPhi| < 0.01,  closest distance cut = %1.2f", foundz.size(), DCA_cut));
    //         // h -> SetTitle(Form("DCA indicated with Gaussian Fit, |dPhi| < 0.01, closest distance cut = %1.2f", DCA_cut));
    //     }
    // }
    // else if (dim == "2D") {
    //     g0 -> SetMarkerStyle(29);
    //     g0 -> SetMarkerSize(2);
    //     g0 -> SetMarkerColor(kBlue);
    //     g0 -> SetLineWidth(3);
    //     g0 -> SetLineColor(kWhite);
    //     if (select == "zscanE") {
    //         g0 -> SetTitle("Z Scan with Errors, |dEta| < 0.01 (for small N) or 0.001 (for large N), |dPhi| < 0.1");
    //     } else if (select == "zscan") {
    //         g0 -> SetTitle("Traditional Z Scan method, |dEta| < 0.01 (for small N) or 0.001 (for large N), |dPhi| < 0.1");
    //     } else if (select == "e") {
    //         g0 -> SetTitle("Nearest Z with Errors, |dPhi| < 0.1");
    //     } 
    //     else if (select == "zscanN") {
    //         g0 -> SetTitle("Z Scan with Background Normalized");
    //     }
    //     else {
    //         g0 -> SetTitle(Form("DCA indicated by greatest bin, %lu events, |dPhi| < 0.01, closest distance cut = %1.2f", foundz.size(), DCA_cut));
    //         // g0 -> SetTitle(Form("DCA indicated with Gaussian Fit, |dPhi| < 0.01, closest distance cut = %1.2f", DCA_cut));
    //     }
        
    //     gStyle -> SetTitleW(0.7);  //per cent of the pad width
    //     gStyle -> SetTitleH(0.08); //per cent of the pad height
    //     g0 -> GetXaxis() -> SetTitle("Number of particles");
    //     g0 -> GetXaxis() -> SetTitleSize(0.06);
    //     g0 -> GetXaxis() -> SetLabelSize(0.04);
    //     g0 -> GetXaxis() -> CenterTitle(true);
    //     g0 -> GetYaxis() -> SetTitle("Found_z - True_z [mm]");
    //     g0 -> GetYaxis() -> SetTitleSize(0.06);
    //     g0 -> GetYaxis() -> SetLabelSize(0.025);
    //     g0 -> GetYaxis() -> CenterTitle(true);
    //     g0 -> SetMinimum(-300); // Setting y range;
    //     g0 -> SetMaximum(300);  // Setting y range;
    //     g0 -> GetYaxis() -> SetTitleOffset(0.8); 
    //     g0 -> GetXaxis() -> SetTitleOffset(0.8); 
    //     // g0 -> GetXaxis() -> SetLimits(0, 6000); // Setting x range;
    //     g0 -> Draw("AP SAME");

    //     gPad->SetGrid(5, 2); gPad->Update();

    //     TLine *line1 = new TLine(0, -16, 10000, -16);
    //     line1->SetLineColor(kRed);
    //     line1->SetLineWidth(2);
    //     line1->Draw("same");  // Draw line on the same canvas
    //     TLine *line2 = new TLine(0, 16, 10000, 16);
    //     line2->SetLineColor(kRed);
    //     line2->SetLineWidth(2);
    //     line2->Draw("same");  // Draw line on the same canvas
    //     TLegend *legend = new TLegend(0.4,0.8,0.9,0.9);
    //     // legend->SetHeader("Legend","C"); // option "C" allows to center the header
    //     legend->AddEntry(g0,"Segments expanded between [-8, 8]mm","f");
    //     legend->SetTextSize(0.04);
    //     legend->Draw("same");
    // }
    // else if (dim == "zVSdz") {
    //     g1 -> SetMarkerStyle(29);
    //     g1 -> SetMarkerSize(1);
    //     g1 -> SetMarkerColor(kBlue - 7);
    //     g1 -> SetLineWidth(3);
    //     g1 -> SetLineColor(kWhite);
    //     if (select == "zscanE") {
    //         g1 -> SetTitle("Z Scan with Errors, |dEta| < 0.01 (for small N) or 0.001 (for large N), |dPhi| < 0.1");
    //     } else if (select == "zscan") {
    //         g1 -> SetTitle("Traditional Z Scan method, |dEta| < 0.01 (for small N) or 0.001 (for large N), |dPhi| < 0.1");
    //     } else if (select == "e") {
    //         g1 -> SetTitle("Nearest Z with Errors, |dPhi| < 0.1");
    //     } 
    //     else if (select == "zscanN") {
    //         g1 -> SetTitle("Z Scan with Background Normalized");
    //     }
    //     else if (select == "DCA") {
    //         g1 -> SetTitle(Form("DCA indicated by greatest bin, %lu events, |dPhi| < 0.01, closest distance cut = %1.2f", foundz.size(), DCA_cut));
    //     }
    //     else {
    //         g1 -> SetTitle(Form("DCA indicated with Gaussian Fit, |dPhi| < 0.01, closest distance cut = %1.2f", DCA_cut));
    //     }
        
    //     gStyle -> SetTitleW(0.7);  //per cent of the pad width
    //     gStyle -> SetTitleH(0.08); //per cent of the pad height
    //     g1 -> GetXaxis() -> SetTitle("True Z Vertex [mm]");
    //     g1 -> GetXaxis() -> SetTitleSize(0.06);
    //     g1 -> GetXaxis() -> SetLabelSize(0.04);
    //     g1 -> GetXaxis() -> CenterTitle(true);
    //     g1 -> GetYaxis() -> SetTitle("Found_z - True_z [mm]");
    //     g1 -> GetYaxis() -> SetTitleSize(0.06);
    //     g1 -> GetYaxis() -> SetLabelSize(0.025);
    //     g1 -> GetYaxis() -> CenterTitle(true);
    //     g1 -> SetMinimum(-100); // Setting y range;
    //     g1 -> SetMaximum(100);  // Setting y range;
    //     g1 -> GetYaxis() -> SetTitleOffset(0.8); 
    //     g1 -> GetXaxis() -> SetTitleOffset(0.8); 
    //     // g1 -> GetXaxis() -> SetLimits(0, 6000); // Setting x range;
    //     g1 -> Draw("AP SAME");

    //     gPad->SetGrid(5, 2); gPad->Update();

    //     TLine *line1 = new TLine(-500, 0, +500, 0);
    //     line1->SetLineColor(kRed);
    //     line1->SetLineWidth(2);
    //     line1->Draw("same");  // Draw line on the same canvas
    //     TLegend *legend = new TLegend(0.4,0.8,0.9,0.9);
    //     // legend->SetHeader("Legend","C"); // option "C" allows to center the header
    //     legend->AddEntry(g1,"Segments expanded between [0, 16] mm","f");
    //     legend->SetTextSize(0.04);
    //     legend->Draw("same");
    // }
    // else {
    //     double ymin = g1->GetHistogram()->GetMinimum();
    //     double ymax = g1->GetHistogram()->GetMaximum();
    //     int nbins = 1000;
    //     TH1F *h1 = new TH1F("ProjectionY", "Histogram of vtx_z_resolution vs vtx_z 1D projection;Found_z - True_z [mm];Entries", nbins, ymin, ymax);
    //     int npoints = g1->GetN();
    //     double x, y;
    //     for (int i = 0; i < npoints; i++) {
    //         g1->GetPoint(i, x, y);
    //         h1->Fill(y);
    //     }
    //     h1->Draw();
    //     h1 -> GetXaxis() -> CenterTitle(true);
    //     h1 -> GetYaxis() -> CenterTitle(true);
    //     h1 -> GetXaxis() -> SetRangeUser(-60, 60); // Setting x range;
    //     h1 -> SetFillColor(kYellow - 7);
    //     h1 -> SetLineWidth(1);
    //     h1 -> SetFillStyle(1001);
    //     static TLine *l1 = new TLine();
    //     l1 -> SetLineColor(kRed);
    //     l1 -> SetLineStyle(0);
    //     l1 -> SetLineWidth(2);
    //     l1 -> SetX1(0);   l1 -> SetY1(0);
    //     l1 -> SetX2(0);   l1 -> SetY2(h1->GetMaximum());
    //     l1 -> Draw("same");
    //     gPad -> SetLogy();
    // }
}

void INTTdPhiAnalysis (TTree *EventTree, string filePath, std::vector<string> options) {
    ifstream myfile(filePath);
    if (!myfile.is_open()){
		std::cout << "Unable to open linelabel" << std::endl;
		system("read -n 1 -s -p \"Press any key to continue...\" echo");
		exit(1);
 	}

    string line, value;
    std::vector<int> index, event, NHits;
    std::vector<double> foundZ, MBD_true_z, MBD_cen;

    getline(myfile, line);
    while (getline(myfile, line)) {
        stringstream data(line);
        getline(data, value, ',');  int i     = std::stoi(value);
        getline(data, value, ',');  int e     = std::stoi(value);
        getline(data, value, ',');  int N     = std::stoi(value);
        getline(data, value, ',');  double f  = std::stod(value)*10;
        getline(data, value, ',');  double Mz = std::stod(value)*10;
        getline(data, value, ',');  double Mc = std::stod(value);

        if (f >= -250. && f <= -50. && Mc <= 0.7) {
            index.push_back(i);     event.push_back(e);         NHits.push_back(N);
            foundZ.push_back(f);    MBD_true_z.push_back(Mz);   MBD_cen.push_back(Mc);
        }
    }

    Long64_t nEntries = EventTree -> GetEntries();
    int event25, NClus;
    float MBD_centrality, MBD_z_vtx;
    std::vector<float> *ClusX     = nullptr;
    std::vector<float> *ClusY     = nullptr;
    std::vector<float> *ClusZ     = nullptr;
    std::vector<int>   *ClusLayer = nullptr;
    std::vector<float> *ClusR     = nullptr;    // FYI only;
    std::vector<float> *ClusPhi   = nullptr;    // FYI only;
    std::vector<float> *ClusEta   = nullptr;    // FYI only;

    TBranch *branch25 = (TBranch*)EventTree->GetListOfBranches()->At(25);    // event25;
    TBranch *branch10 = (TBranch*)EventTree->GetListOfBranches()->At(10);    // NClus;
    TBranch *branch11 = (TBranch*)EventTree->GetListOfBranches()->At(11);    // ClusLayer;
    TBranch *branch12 = (TBranch*)EventTree->GetListOfBranches()->At(12);    // ClusX;
    TBranch *branch13 = (TBranch*)EventTree->GetListOfBranches()->At(13);    // ClusY;
    TBranch *branch14 = (TBranch*)EventTree->GetListOfBranches()->At(14);    // ClusZ;
    TBranch *branch30 = (TBranch*)EventTree->GetListOfBranches()->At(30);    // MBD_centrality;
    TBranch *branch31 = (TBranch*)EventTree->GetListOfBranches()->At(31);    // MBD_z_vtx;
    TBranch *branch15 = (TBranch*)EventTree->GetListOfBranches()->At(15);    // ClusR;
    TBranch *branch16 = (TBranch*)EventTree->GetListOfBranches()->At(16);    // ClusPhi;
    TBranch *branch17 = (TBranch*)EventTree->GetListOfBranches()->At(17);    // ClusEta;
    
    branch25->SetAddress(&event25);
    branch10->SetAddress(&NClus);
    branch11->SetAddress(&ClusLayer);
    branch12->SetAddress(&ClusX);
    branch13->SetAddress(&ClusY);
    branch14->SetAddress(&ClusZ);
    branch15->SetAddress(&ClusR);
    branch16->SetAddress(&ClusPhi);
    branch17->SetAddress(&ClusEta);
    branch30->SetAddress(&MBD_centrality);
    branch31->SetAddress(&MBD_z_vtx);

    double dPhi, phi_0, phi_1;
    std::vector<double> Phi0, Phi1;
    int N = 1000;
    double range_min = -M_PI;
    double range_max = M_PI;
    double bin_width = (range_max - range_min) / N;
    TH1D *h_unmixed_dPhi  = new TH1D("dPhi of unmixed", Form("dPhi of unmixed %lu events; dPhi value; # of counts", event.size()), N, range_min, range_max);
    TH2D *h_dPhi_Z = new TH2D("dPhi of unmixed 2D", Form("dPhi of unmixed %lu events; dPhi value; # of counts", event.size()), N, range_min, range_max, 20, (-25 - 0.5)*10, (-5. + 0.5)*10);
    for (int i = 0; i < event.size(); i++) {
        std::cout << i << std::endl;
        branch25->GetEntry(index[i]);  branch10->GetEntry(index[i]);  branch11->GetEntry(index[i]);
        branch12->GetEntry(index[i]);  branch13->GetEntry(index[i]);  branch14->GetEntry(index[i]);
        branch15->GetEntry(index[i]);  branch16->GetEntry(index[i]);  branch17->GetEntry(index[i]);
        branch30->GetEntry(index[i]);  branch31->GetEntry(index[i]);

        // if (event25 != event[i] || NHits[i] != NClus || std::abs(MBD_cen[i] - MBD_centrality) > 3e-08) {
        //     std::cout << red << "# of hit insanity at " << event[i] << ',' << event25
        //     << ", " << NHits[i] << "," << NClus << ", " << MBD_cen[i] << "," << MBD_centrality << color_reset << std::endl;
        // }

        for (int j = 0; j < ClusX->size(); j++) {
            double phi = std::atan2(ClusY->at(j), ClusX->at(j));
            if (ClusLayer->at(j) == 3 || ClusLayer->at(j) == 4) {
                Phi0.push_back(phi);
            }
            else {
                Phi1.push_back(phi);
            }
        }
        double fz = foundZ[i];
        for (int k = 0; k < Phi0.size(); k++) {
            for (int l = 0; l < Phi1.size(); l++) {
                dPhi = Phi0[k] - Phi1[l];
                if (dPhi > M_PI)    dPhi = dPhi - 2*M_PI;
                if (dPhi < -M_PI)   dPhi = dPhi + 2*M_PI;
                h_unmixed_dPhi -> Fill(dPhi);
                h_dPhi_Z -> Fill(dPhi, fz);
            }
        }
        Phi0.clear();   Phi1.clear();
    }
    angularPlot1D(h_unmixed_dPhi, options, "dPhi of unmixed");
    angularPlot2D(h_dPhi_Z, options, "dPhi of unmixed vs Z");
}

void ResultsAnalysis (std::string method) {
    // Open the ROOT file
    TFile *file = TFile::Open("../External/Data_CombinedNtuple_Run20869_HotDead_BCO_ADC_Survey.root");
    if (!file || file->IsZombie()) {
        std::cerr << "Error opening file!" << std::endl;
        return;
    }
    // Get the TTree
    TTree *tree = (TTree*)file->Get("EventTree");
    if (!tree) {
        std::cerr << "Error getting TTree!" << std::endl;
        return;
    }

    std::string filePath = "zFindingResults/new_foundZ_DCAfit_0_16_-001_001_step2e-1_z_25_5.txt";
    std::vector<std::string> sub_options = splitString(method, '_');
    if (sub_options.size() == 2) {
        sub_options.push_back("");
    }
    else if (sub_options.size() == 1) {
        sub_options.push_back("");
        sub_options.push_back("");
    }

    if (sub_options[0] == "foundz") INTTZAnalysisLite(filePath, sub_options);
    if (sub_options[0] == "dPhi")   INTTdPhiAnalysis(tree, filePath, sub_options);
}
