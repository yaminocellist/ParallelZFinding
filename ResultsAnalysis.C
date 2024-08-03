#include "headers/Analysis.h"
#include <TStopwatch.h>

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

        if (Mz >= -250. && Mz <= -50. && f >= -250. && f <= -50. && Mc <= 0.7) {
        // if ((Mz >= -250. && Mz <= -180. && Mc <= 0.7) || (Mz >= -160. && Mz <= -50. && Mc <= 0.7)) {
            index.push_back(i);     event.push_back(e);         NClus.push_back(N);
            foundZ.push_back(f);    MBD_z_vtx.push_back(Mz);    MBD_centrality.push_back(Mc);
        }
    }
    
    TH1D *h_resolution = new TH1D("found z resolution", Form("INTT found z resolution of %lu events;dZ [mm];# of counts", index.size()), 151, -1e2, 1e2);
    TH2D *h = new TH2D("", ";MBD Z;INTT found Z", 100, -250, -50, 100, -250, -50);
    TGraph *g_NFunction    = new TGraph();
    TGraph *g_ZFunction    = new TGraph();
    TGraph *g_CFunction    = new TGraph();
    TGraph *g_ZCorrelation = new TGraph();
    TGraph *g_Calibration  = new TGraph();
    double z_resolution;
    for (int i = 0; i < index.size(); i++) {
        h -> Fill(MBD_z_vtx[i], foundZ[i]);
        z_resolution = foundZ[i] - MBD_z_vtx[i];
        h_resolution -> Fill(z_resolution);
        g_NFunction  -> SetPoint(g_NFunction->GetN(), NClus[i], z_resolution);
        g_ZFunction  -> SetPoint(g_ZFunction->GetN(), MBD_z_vtx[i], z_resolution);
        g_CFunction  -> SetPoint(g_CFunction->GetN(), MBD_centrality[i], z_resolution);
        g_ZCorrelation -> SetPoint(g_ZCorrelation->GetN(), foundZ[i], MBD_z_vtx[i]);
        g_Calibration -> SetPoint(g_Calibration->GetN(), MBD_z_vtx[i], MBD_z_vtx[i]);
    }
    // h->Draw("lego2");
    TH1D *h_projX = h->ProjectionX("projX", h->GetXaxis()->FindBin(-250), h->GetXaxis()->FindBin(-50));
    h_projX -> Draw();
    // TH1D *h_projY = h->ProjectionY("projY", h->GetYaxis()->FindBin(-250), h->GetYaxis()->FindBin(-50));
    // h_projY -> Draw("same");
    // h_projY -> SetLineColor(2);
    // ZResolutionSinglePlot(h_resolution, options, "foundZResolution_1D");
    // TGraphSinglePlot(g_NFunction, Form("INTT found z resolution of %d events", g_NFunction->GetN()), "# of Hits", "foundZResolution_vs_NClus");
    // TGraphSinglePlot(g_ZFunction, Form("INTT found z resolution of %d events", g_ZFunction->GetN()), "MBD z vtx [mm]", "foundZResolution_vs_MBDZ");
    // TGraphSinglePlot(g_CFunction, Form("INTT found z resolution of %d events", g_CFunction->GetN()), "MBD centrality", "foundZResolution_vs_MBC");
    // TGraphSinglePlot_Squared(g_ZCorrelation, Form("INTT found z resolution of %d events", g_ZCorrelation->GetN()), "INTT found z [mm]", "MBD z vtx [mm]", "foundZ_vs_MBZ");
    TMultiGraph *mg = new TMultiGraph();
    mg -> Add(g_ZCorrelation);
    mg -> Add(g_Calibration);
    // TGraphMultiPlot(mg, Form("INTT found z resolution of %d events", g_ZCorrelation->GetN()), "INTT found z [mm]", "MBD z vtx [mm]", "foundZ_vs_MBZ_with_Cal");

    // TCanvas *canvas = new TCanvas("c2","c2", 0, 0,1350,800);
}

void INTTMixingEvent (TTree *EventTree, string filePath, std::vector<string> options) {
    ifstream myfile(filePath);
    if (!myfile.is_open()){
		std::cout << "Unable to open linelabel" << std::endl;
		system("read -n 1 -s -p \"Press any key to continue...\" echo");
		exit(1);
 	}
    double z_lower_range = 24.05;
    double z_upper_range = 23.85;
    double c_lower_range = 24.05;
    double c_upper_range = 23.85;
    string line, value;
    std::vector<int>    index, event, NHits;
    std::vector<double> foundZ, MBD_true_z, MBD_cen;
    getline(myfile, line);
    while (getline(myfile, line)) {
        stringstream data(line);
        getline(data, value, ',');  int i     = std::stoi(value);
        getline(data, value, ',');  int e     = std::stoi(value);
        getline(data, value, ',');  int N     = std::stoi(value);
        getline(data, value, ',');  double f  = std::stod(value);
        getline(data, value, ',');  double Mz = std::stod(value);
        getline(data, value, ',');  double Mc = std::stod(value);

        if (options[1] == "Z") {
            // if (f >= -21.05 && f <= -20.85 && Mc <= 0.7) {
            if (f >= -z_lower_range && f <= -z_upper_range && Mc <= 0.7) {
                index.push_back(i);     event.push_back(e);         NHits.push_back(N);
                foundZ.push_back(f);    MBD_true_z.push_back(Mz);   MBD_cen.push_back(Mc);
            }
        }
        else if (options[1] == "cen") {
            if (f >= -15.05 && f <= -14.85 && Mc <= 0.1 && Mc >= 0.05) {
                index.push_back(i);     event.push_back(e);         NHits.push_back(N);
                foundZ.push_back(f);    MBD_true_z.push_back(Mz);   MBD_cen.push_back(Mc);
            }
        }
        
    }
    // std::cout << event.size() << std::endl;

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

    int target = std::stoi(options[2]) > event.size() ? event.size() : std::stoi(options[2]);
    double phi, dPhi, phi_0, phi_1;
    std::vector<double> Phi0, Phi1, mixed_Phi0, mixed_Phi1;
    std::vector<std::vector <double>> event_Phi0, event_Phi1;
    int N = 1000, clus_layer;
    double range_min = -M_PI;
    double range_max = M_PI;
    double bin_width = (range_max - range_min) / N;
    TH1D *h_unmixed_dPhi = new TH1D("dPhi of unmixed", Form("dPhi of unmixed %d events;dPhi value;# of counts", target), N, range_min, range_max);
    TH1D *h_mixed_dPhi   = new TH1D("dPhi of mixed", Form("dPhi of %d mixed events;dPhi value;# of counts", target), N, range_min, range_max);

    // dPhi of unmixed events:
    
    for (int i = 0; i < target; i++) {
        branch16->GetEntry(index[i]);   branch11->GetEntry(index[i]);
        event_Phi0.push_back(std::vector <double>());   event_Phi1.push_back(std::vector <double>());
        for (int j = 0; j < ClusPhi->size(); j++) {
            phi        = ClusPhi->at(j);
            clus_layer = ClusLayer->at(j);
            if (clus_layer == 3 || clus_layer == 4) {
                Phi0.push_back(phi);
                event_Phi0[i].push_back(phi);
            }
            else {
                Phi1.push_back(phi);
                event_Phi1[i].push_back(phi);
            }
        }
        for (int k = 0; k < Phi0.size(); k++) {
            for (int l = 0; l < Phi1.size(); l++) {
                dPhi = Phi0[k] - Phi1[l];
                if (dPhi > M_PI)    dPhi = dPhi - 2*M_PI;
                if (dPhi < -M_PI)   dPhi = dPhi + 2*M_PI;
                h_unmixed_dPhi -> Fill(dPhi);
            }
        }
        Phi0.clear();   Phi1.clear();
    }

    // std::cout << event_Phi0.size() << ", " << event_Phi1.size() << std::endl;
    // dPhi of all events mixed up (one pair another, not mixing them all up):
    // for (int i = 0; i < event_Phi0.size(); i++) {
    //     std::cout << i << std::endl;
    //     for (int j = i; j < event_Phi1.size(); j++) {
    //         mixed_Phi0.insert(mixed_Phi0.end(), event_Phi0[i].begin(), event_Phi0[i].end());
    //         mixed_Phi0.insert(mixed_Phi0.end(), event_Phi0[j].begin(), event_Phi0[j].end());
    //         mixed_Phi1.insert(mixed_Phi1.end(), event_Phi1[i].begin(), event_Phi1[i].end());
    //         mixed_Phi1.insert(mixed_Phi1.end(), event_Phi1[j].begin(), event_Phi1[j].end());
    //         // std::cout << "  ---  " << event_Phi0[i].size() << ", " << event_Phi0[j].size() << ", " << event_Phi1[i].size() << ", " << event_Phi1[j].size() << ", " << mixed_Phi0.size() << ", " << mixed_Phi1.size() << std::endl;
    //         for (int k = 0; k < mixed_Phi0.size(); k++) {
    //             for (int l = 0; l < mixed_Phi1.size(); l++) {
    //                 dPhi = mixed_Phi0[k] - mixed_Phi1[l];
    //                 if (dPhi > M_PI)    dPhi = dPhi - 2*M_PI;
    //                 if (dPhi < -M_PI)   dPhi = dPhi + 2*M_PI;
    //                 h_mixed_dPhi -> Fill(dPhi);
    //             }
    //         }
    //         mixed_Phi0.clear(); mixed_Phi1.clear();
    //     }
    // }

    for (int i = 0; i < event_Phi0.size(); i++) {
        std::cout << i << std::endl;
        for (int j = i; j < event_Phi1.size(); j++) {
            // Process pairs from event_Phi0[i] and event_Phi1[j]
            for (double phi0 : event_Phi0[i]) {
                for (double phi1 : event_Phi1[j]) {
                    dPhi = phi0 - phi1;
                    if (dPhi > M_PI) dPhi -= 2 * M_PI;
                    if (dPhi < -M_PI) dPhi += 2 * M_PI;
                    h_mixed_dPhi->Fill(dPhi);
                }
            }
            // Process pairs from event_Phi0[j] and event_Phi1[i] to ensure symmetry
            for (double phi0 : event_Phi0[j]) {
                for (double phi1 : event_Phi1[i]) {
                    dPhi = phi0 - phi1;
                    if (dPhi > M_PI) dPhi -= 2 * M_PI;
                    if (dPhi < -M_PI) dPhi += 2 * M_PI;
                    h_mixed_dPhi->Fill(dPhi);
                }
            }
        }
    }
    
    TCanvas *can1 = new TCanvas("c1d","c1d",0,50,2100,1200);
    double phi_range_low = -2.4, phi_range_high = -1.8;
    int bin_range_low = h_mixed_dPhi->FindBin(phi_range_low), bin_range_high = h_mixed_dPhi->FindBin(phi_range_high);
    double max_unmixed = -1, max_mixed = -1, current_binContent;
    for (int bin = bin_range_low; bin <= bin_range_high; bin++) {
        current_binContent = h_unmixed_dPhi->GetBinContent(bin);
        if (max_unmixed < current_binContent)  max_unmixed = current_binContent;
        current_binContent = h_mixed_dPhi->GetBinContent(bin);
        if (max_mixed < current_binContent)  max_mixed = current_binContent;
    }
    h_unmixed_dPhi -> Add(h_unmixed_dPhi, max_mixed/max_unmixed - 1);
    h_unmixed_dPhi -> Draw("SAME");

    max_mixed   = h_mixed_dPhi->GetBinContent(h_mixed_dPhi->GetMaximumBin());
    max_unmixed = h_unmixed_dPhi->GetBinContent(h_unmixed_dPhi->GetMaximumBin());
    std::cout << h_unmixed_dPhi->GetBinContent(h_unmixed_dPhi->GetMaximumBin()) << ", " << h_mixed_dPhi->GetBinContent(h_mixed_dPhi->GetMaximumBin())  << std::endl;
    if (max_unmixed > max_mixed) {
        h_unmixed_dPhi -> GetYaxis() -> SetRangeUser(1e7, max_unmixed*1.2);
        h_mixed_dPhi -> GetYaxis() -> SetRangeUser(1e7, max_unmixed*1.2);
    }   
    else {
        h_unmixed_dPhi -> GetYaxis() -> SetRangeUser(1e7, max_mixed*1.2);
        h_mixed_dPhi -> GetYaxis() -> SetRangeUser(1e7, max_mixed*1.2);
    }                      
    h_mixed_dPhi -> Draw("SAME");
    h_mixed_dPhi -> SetLineColor(2);

    double pi = TMath::Pi();
    int bin_min  = 1;  // The first bin
    int bin_max  = h_unmixed_dPhi->GetNbinsX();  // The last bin
    // Calculate bin positions for each label
    int bin_pi   = bin_max;
    int bin_0    = bin_min + (bin_max - bin_min)/2;
    int binPi_2  = bin_min + 3*(bin_max - bin_min)/4;
    int bin_pi_2 = bin_min + (bin_max - bin_min)/4;
    // Set the labels at the calculated positions
    h_mixed_dPhi->GetXaxis()->SetBinLabel(bin_0, "0");
    h_mixed_dPhi->GetXaxis()->SetBinLabel(bin_pi_2, "#frac{-#pi}{2}");
    h_mixed_dPhi->GetXaxis()->SetBinLabel(binPi_2, "#frac{#pi}{2}");
    h_mixed_dPhi->GetXaxis()->SetBinLabel(bin_pi, "#pi");
    h_mixed_dPhi->GetXaxis()->SetBinLabel(bin_min, "-#pi");
    // Ensure the custom labels are displayed by setting the number of divisions
    h_mixed_dPhi->GetXaxis()->SetNdivisions(9, 0, 0, kFALSE);
    h_mixed_dPhi->GetXaxis()->SetLabelSize(0.04);
    // Update histogram to refresh the axis
    // h[0]->Draw("HIST");
    h_mixed_dPhi->GetXaxis()->LabelsOption("h"); // Draw the labels vertically

    h_unmixed_dPhi->GetXaxis()->SetBinLabel(bin_0, "0");
    h_unmixed_dPhi->GetXaxis()->SetBinLabel(bin_pi_2, "#frac{-#pi}{2}");
    h_unmixed_dPhi->GetXaxis()->SetBinLabel(binPi_2, "#frac{#pi}{2}");
    h_unmixed_dPhi->GetXaxis()->SetBinLabel(bin_pi, "#pi");
    h_unmixed_dPhi->GetXaxis()->SetBinLabel(bin_min, "-#pi");
    // Ensure the custom labels are displayed by setting the number of divisions
    h_unmixed_dPhi->GetXaxis()->SetNdivisions(9, 0, 0, kFALSE);
    h_unmixed_dPhi->GetXaxis()->SetLabelSize(0.04);
    // Update histogram to refresh the axis
    // h[0]->Draw("HIST");
    h_unmixed_dPhi->GetXaxis()->LabelsOption("h"); // Draw the labels vertically

    h_unmixed_dPhi -> GetXaxis() -> CenterTitle(true);
    h_unmixed_dPhi -> GetYaxis() -> CenterTitle(true);
    h_mixed_dPhi -> GetXaxis() -> CenterTitle(true);
    h_mixed_dPhi -> GetYaxis() -> CenterTitle(true);

    TLegend *lg = new TLegend(0.12, 0.8, 0.33, 0.9);
    lg -> AddEntry(h_mixed_dPhi, Form("z_vtx between %2.2f and %2.2f", z_lower_range, z_upper_range), "l");
    gStyle -> SetLegendTextSize(.023);
    lg->Draw("same");
    
    if (options[1] == "Z")     can1 -> SaveAs("../External/zFindingPlots/dPhi_mixed_Z.png");
    if (options[1] == "cen")   can1 -> SaveAs("../External/zFindingPlots/dPhi_mixed_centrality.png");
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

        if (Mz >= -250. && Mz <= -50. && Mc <= 0.7) {
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
    TH1D *h_unmixed_dPhi = new TH1D("dPhi of unmixed", Form("dPhi of unmixed %lu events; dPhi value; # of counts", event.size()), N, range_min, range_max);
    TH2D *h_dPhi_Z       = new TH2D("dPhi of unmixed v.s. Z 2D", Form("dPhi of unmixed %lu events; dPhi value; foundZ vtx [mm]", event.size()), N, range_min, range_max, 20, (-25 - 0.5)*10, (-5. + 0.5)*10);
    TH2D *h_dPhi_cen     = new TH2D("dPhi of unmixed v.s. MBD centrality 2D", Form("dPhi of unmixed %lu events; dPhi value; MBD centrality", event.size()), N, range_min, range_max, 14, 0.0, 0.7);

    // std::cout << options[1] << std::endl;
    int target = std::stoi(options[2]) > event.size() ? event.size() : std::stoi(options[2]);
    if (options[1] == "fzStrips") {
        TCanvas *can1 = new TCanvas("c1d","c1d",0,50,2100,1200);
        TH1D *h_all_foundz = new TH1D("", "INTT found z v.s. MBD Z;Z Vertex Position [mm];# of counts", 200, -250., -50.);
        TH1D *h_all_MBD_Z  = new TH1D("", "INTT found z v.s. MBD Z;Z Vertex Position [mm];# of counts", 100, -250., -50.);
        for (int i = 0; i < event.size(); i++) {
            h_all_foundz -> Fill(foundZ[i]);
            h_all_MBD_Z  -> Fill(MBD_true_z[i]);
        }
        h_all_foundz -> Draw();
        h_all_foundz -> SetLineWidth(3);
        h_all_foundz -> GetXaxis() -> CenterTitle(true);   h_all_foundz -> GetYaxis() -> CenterTitle(true);
        h_all_MBD_Z  -> GetXaxis() -> CenterTitle(true);   h_all_MBD_Z  -> GetYaxis() -> CenterTitle(true);
        h_all_MBD_Z  -> Draw("SAME");
        h_all_MBD_Z  -> SetLineColor(2);
        h_all_MBD_Z  -> SetLineWidth(5);
        TLegend *lg = new TLegend(0.62, 0.8, 0.73, 0.9);
        lg -> AddEntry(h_all_foundz, "INTT found Z", "l");
        lg -> AddEntry(h_all_MBD_Z, "MBD Z vertex", "l");
        gStyle -> SetLegendTextSize(.028);
        lg->Draw("same");
        can1 -> SaveAs("../External/zFindingPlots/INTT_MBD_Z_comparison.png");

        foundZSlicer(h_all_foundz);
    }
    if (options[1] == "fCenStrips") {
        TCanvas *can1 = new TCanvas("c1d","c1d",0,50,2100,1200);
        TH2D *h_Z_cen = new TH2D("", "Found Z vertex in bins of centrality;Z Vertex Position [mm];MBD centrality", 200, -250., -50., 14, 0.0, 0.7);
        for (int i = 0; i < event.size(); i++) {
            h_Z_cen -> Fill(foundZ[i], MBD_cen[i]);
        }
        h_Z_cen -> Draw("lego2");
        // h_Z_cen -> SetLineWidth(3);
        h_Z_cen -> GetXaxis() -> CenterTitle(true);   h_Z_cen -> GetYaxis() -> CenterTitle(true);
        h_Z_cen -> GetXaxis() -> SetTitleOffset(2);
        h_Z_cen -> GetYaxis() -> SetTitleOffset(2);
        // TLegend *lg = new TLegend(0.62, 0.8, 0.73, 0.9);
        // lg -> AddEntry(h_all_foundz, "INTT found Z", "l");
        // lg -> AddEntry(h_all_MBD_Z, "MBD Z vertex", "l");
        // gStyle -> SetLegendTextSize(.028);
        // lg->Draw("same");
        can1 -> SaveAs("../External/zFindingPlots/foundZ_vs_MBD_centrality.png");
    }
    if (options[1] == "nomix") {
        for (int i = 0; i < target; i++) {
            // std::cout << i << std::endl;
            // branch25->GetEntry(index[i]);  branch10->GetEntry(index[i]);  branch11->GetEntry(index[i]);
            // branch12->GetEntry(index[i]);  branch13->GetEntry(index[i]);  branch14->GetEntry(index[i]);
            // branch15->GetEntry(index[i]);  branch16->GetEntry(index[i]);  branch17->GetEntry(index[i]);
            // branch30->GetEntry(index[i]);  branch31->GetEntry(index[i]);
            branch16->GetEntry(index[i]);   branch11->GetEntry(index[i]);

            // if (event25 != event[i] || NHits[i] != NClus || std::abs(MBD_cen[i] - MBD_centrality) > 3e-08) {
            //     std::cout << red << "# of hit insanity at " << event[i] << ',' << event25
            //     << ", " << NHits[i] << "," << NClus << ", " << MBD_cen[i] << "," << MBD_centrality << color_reset << std::endl;
            // }

            double fz = foundZ[i], cen = MBD_cen[i];
            for (int j = 0; j < ClusPhi->size(); j++) {
                // double phi = std::atan2(ClusY->at(j), ClusX->at(j));
                double phi = ClusPhi->at(j);
                if (ClusLayer->at(j) == 3 || ClusLayer->at(j) == 4) {
                    Phi0.push_back(phi);
                }
                else {
                    Phi1.push_back(phi);
                }
            }
            for (int k = 0; k < Phi0.size(); k++) {
                for (int l = 0; l < Phi1.size(); l++) {
                    dPhi = Phi0[k] - Phi1[l];
                    if (dPhi > M_PI)    dPhi = dPhi - 2*M_PI;
                    if (dPhi < -M_PI)   dPhi = dPhi + 2*M_PI;
                    h_unmixed_dPhi -> Fill(dPhi);
                    // h_dPhi_Z       -> Fill(dPhi, fz);
                    // h_dPhi_cen     -> Fill(dPhi, cen);
                }
            }
            Phi0.clear();   Phi1.clear();
        }
        angularPlot1D(h_unmixed_dPhi, options, "dPhi of unmixed");
        // angularPlot2D(h_dPhi_Z, options, "dPhi of unmixed vs Z");
        // angularPlot2D(h_dPhi_cen, options, "dPhi of unmixed vs MBD centrality");
        // angularPlot3D(h_dPhi_cen, options, "dPhi of unmixed vs MBD centrality 3D");
        // angularPlot3D(h_dPhi_Z, options, "dPhi of unmixed vs Z");
    }
    else if (options[1] == "perCenOnOne") {
        std::vector<double> Phi0, Phi1;
        TH1D *h_onOne[14];
        h_onOne[0] = new TH1D("dPhi of different centralities", "dPhi of different centralities;dPhi;# of counts", N, range_min, range_max);
        for (int i = 1; i < 14; i++) {
            h_onOne[i] = new TH1D(Form("dPhi of %f to %f", static_cast<double>(i)*0.05, static_cast<double>(i + 1)*0.05), Form("dPhi of %f to %f;dPhi;# of counts", static_cast<double>(i)*0.05, static_cast<double>(i + 1)*0.05), N, range_min, range_max);
        }
        for (int i = 0; i < target; i++) {
            // std::cout << i << std::endl;
            branch16->GetEntry(index[i]);   // ClusPhi;
            branch11->GetEntry(index[i]);   // ClusLayer;
            // branch30->GetEntry(index[i]);   // MBD_centrality;
            // if (std::abs(MBD_centrality - MBD_cen[i]) > 1e-07) {
            //     std::cout << MBD_centrality - MBD_cen[i] << ", " << MBD_cen[i] << std::endl;
            // }
            for (int n = 0; n < 14; n++) {
                if (MBD_cen[i] >= static_cast<double>(n)*0.05 && MBD_cen[i] <= static_cast<double>(n + 1)*0.05) {
                    for (int j = 0; j < ClusPhi->size(); j++) {
                        // double phi = std::atan2(ClusY->at(j), ClusX->at(j));
                        double phi = ClusPhi->at(j);
                        if (ClusLayer->at(j) == 3 || ClusLayer->at(j) == 4) {
                            Phi0.push_back(phi);
                        }
                        else {
                            Phi1.push_back(phi);
                        }
                    }
                    for (int k = 0; k < Phi0.size(); k++) {
                        for (int l = 0; l < Phi1.size(); l++) {
                            dPhi = Phi0[k] - Phi1[l];
                            if (dPhi > M_PI)    dPhi = dPhi - 2*M_PI;
                            if (dPhi < -M_PI)   dPhi = dPhi + 2*M_PI;
                            h_onOne[n] -> Fill(dPhi);
                        }
                    }
                    Phi0.clear();   Phi1.clear();
                }
            }
            
        }
        std::vector<TH1D*> h(h_onOne, h_onOne + 14);
        // ArrayPlot1D_Logy(h, options, "dPhi_per_centralities");
        ArrayPlot1D_Rescale(h, options, "dPhi_per_centralities_rescale");
    }
    else if (options[1] == "perZOnOne") {
        std::vector<double> Phi0, Phi1;
        TH1D *h_onOne[20];
        for (int i = 0; i < 16; i++) {
            h_onOne[i] = new TH1D(Form("dPhi of %1.0f to %1.0f", -5 - static_cast<double>(i), -6 - static_cast<double>(i)), Form("dPhi of %1.0f to %1.0f;dPhi;# of counts", -5 - static_cast<double>(i), -6 - static_cast<double>(i)), N, range_min, range_max);
        }
        h_onOne[16] = new TH1D("dPhis of different Z vertex", "dPhis of different Z vertex;dPhi;# of counts", N, range_min, range_max);
        for (int i = 17; i < 20; i++) {
            h_onOne[i] = new TH1D(Form("dPhi of %1.0f to %1.0f", -5 - static_cast<double>(i), -6 - static_cast<double>(i)), Form("dPhi of %1.0f to %1.0f;dPhi;# of counts", -5 - static_cast<double>(i), -6 - static_cast<double>(i)), N, range_min, range_max);
        }
        for (int i = 0; i < target; i++) {
            // std::cout << i << std::endl;
            branch16->GetEntry(index[i]);   // ClusPhi;
            branch11->GetEntry(index[i]);   // ClusLayer;

            for (int n = 0; n < 20; n++) {
                if (foundZ[i] >= -60. - static_cast<double>(n)*10 && foundZ[i] <= -50. - static_cast<double>(n)*10) {
                    for (int j = 0; j < ClusPhi->size(); j++) {
                        // double phi = std::atan2(ClusY->at(j), ClusX->at(j));
                        double phi = ClusPhi->at(j);
                        if (ClusLayer->at(j) == 3 || ClusLayer->at(j) == 4) {
                            Phi0.push_back(phi);
                        }
                        else {
                            Phi1.push_back(phi);
                        }
                    }
                    for (int k = 0; k < Phi0.size(); k++) {
                        for (int l = 0; l < Phi1.size(); l++) {
                            dPhi = Phi0[k] - Phi1[l];
                            if (dPhi > M_PI)    dPhi = dPhi - 2*M_PI;
                            if (dPhi < -M_PI)   dPhi = dPhi + 2*M_PI;
                            h_onOne[n] -> Fill(dPhi);
                        }
                    }
                    Phi0.clear();   Phi1.clear();
                }
            }
        }
        std::vector<TH1D*> h(h_onOne, h_onOne + 20);
        // ArrayPlot1D_Logy_ver2(h, options, "dPhi_per_Zvtx");
        ArrayPlot1D_Rescale_ver2(h, options, "dPhi_per_Zvtx_rescale");
    }

}



void INTTdEtaAnalysis (TTree *EventTree, string filePath, std::vector<string> options) {
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
        getline(data, value, ',');  double f  = std::stod(value);
        getline(data, value, ',');  double Mz = std::stod(value);
        getline(data, value, ',');  double Mc = std::stod(value);

        if (f >= -25. && f <= -5. && Mc <= 0.7) {
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

    double dZ, R, RSquared, theta, eta, dEta;
    std::vector<double> Eta0, Eta1;
    int N = 2000;
    double range_min = -2*M_PI;
    double range_max = 2*M_PI;
    double bin_width = (range_max - range_min) / N;
    TH1D *h_unmixed_dEta = new TH1D("dEta of unmixed", Form("dEta of unmixed %lu events; dEta value; # of counts", event.size()), N, range_min, range_max);
    // TH2D *h_dPhi_Z       = new TH2D("dPhi of unmixed 2D", Form("dPhi of unmixed %lu events; dPhi value; foundZ vtx [mm]", event.size()), N, range_min, range_max, 20, (-25 - 0.5)*10, (-5. + 0.5)*10);
    // TH2D *h_dPhi_cen     = new TH2D("dPhi of unmixed 2D", Form("dPhi of unmixed %lu events; dPhi value; MBD centrality", event.size()), N, range_min, range_max, 14, 0.0, 0.7);
    for (int i = 0; i < event.size(); i++) {
        // std::cout << i << std::endl;
        branch25->GetEntry(index[i]);  branch10->GetEntry(index[i]);  branch11->GetEntry(index[i]);
        branch12->GetEntry(index[i]);  branch13->GetEntry(index[i]);  branch14->GetEntry(index[i]);
        branch15->GetEntry(index[i]);  branch16->GetEntry(index[i]);  branch17->GetEntry(index[i]);
        branch30->GetEntry(index[i]);  branch31->GetEntry(index[i]);

        double fz = foundZ[i], cen = MBD_cen[i];
        for (int j = 0; j < ClusX->size(); j++) {
            dZ       = ClusZ->at(j) - fz;
            // RSquared = (ClusX->at(j))*(ClusX->at(j)) + (ClusY->at(j))*(ClusY->at(j));
            // R        = std::sqrt(RSquared);
            theta    = std::atan2(ClusR->at(j), dZ);
            // theta    = std::atan2(R, dZ);
            if (dZ >= 0)    eta = -std::log(std::tan(theta/2));
            if (dZ <  0)    eta = std::log(std::tan((M_PI - theta)/2));
            if (ClusLayer->at(j) == 3 || ClusLayer->at(j) == 4) {
                Eta0.push_back(eta);
            }
            else {
                Eta1.push_back(eta);
            }
        }
        for (int k = 0; k < Eta0.size(); k++) {
            for (int l = 0; l < Eta1.size(); l++) {
                dEta = Eta0[k] - Eta1[l];
                h_unmixed_dEta -> Fill(dEta);
                // h_dPhi_Z       -> Fill(dPhi, fz);
                // h_dPhi_cen     -> Fill(dPhi, cen);
            }
        }

        // for (int k = 0; k < Eta0.size(); k++) {
        //     if (Eta0[k] < 2)    h_unmixed_dEta -> Fill(Eta0[k]);
        // }
        // for (int l = 0; l < Eta1.size(); l++) {
        //     if (Eta1[l] < 2)    h_unmixed_dEta -> Fill(Eta1[l]);
        // }
        Eta0.clear();   Eta1.clear();
    }
    angularPlot1D(h_unmixed_dEta, options, "dEta of unmixed");
    // angularPlot2D(h_dPhi_Z, options, "dPhi of unmixed vs Z");
    // angularPlot2D(h_dPhi_cen, options, "dPhi of unmixed vs MBD centrality");
    // angularPlot3D(h_dPhi_cen, options, "dPhi of unmixed vs MBD centrality 3D");
}

void ResultsAnalysis (std::string method) {
    // Start the stopwatch
    TStopwatch timer;
    timer.Start();

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

    if (sub_options[0] == "foundz")      INTTZAnalysisLite(filePath, sub_options);
    if (sub_options[0] == "dPhi")        INTTdPhiAnalysis(tree, filePath, sub_options);
    if (sub_options[0] == "dPhiMix")     INTTMixingEvent(tree, filePath, sub_options);
    if (sub_options[0] == "dEta")        INTTdEtaAnalysis(tree, filePath, sub_options);

    // Stop the stopwatch and print the runtime
    timer.Stop();
    timer.Print();
}
