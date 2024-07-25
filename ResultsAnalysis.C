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
    TH1D *h_unmixed_dPhi = new TH1D("dPhi of unmixed", Form("dPhi of unmixed %lu events; dPhi value; # of counts", event.size()), N, range_min, range_max);
    TH2D *h_dPhi_Z       = new TH2D("dPhi of unmixed v.s. Z 2D", Form("dPhi of unmixed %lu events; dPhi value; foundZ vtx [mm]", event.size()), N, range_min, range_max, 20, (-25 - 0.5)*10, (-5. + 0.5)*10);
    TH2D *h_dPhi_cen     = new TH2D("dPhi of unmixed v.s. MBD centrality 2D", Form("dPhi of unmixed %lu events; dPhi value; MBD centrality", event.size()), N, range_min, range_max, 14, 0.0, 0.7);

    // std::cout << options[1] << std::endl;
    if (options[1] == "nomix") {
        for (int i = 0; i < 10; i++) {
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
                    // h_unmixed_dPhi -> Fill(dPhi);
                    h_dPhi_Z       -> Fill(dPhi, fz);
                    // h_dPhi_cen     -> Fill(dPhi, cen);
                }
            }
            Phi0.clear();   Phi1.clear();
        }
        // angularPlot1D(h_unmixed_dPhi, options, "dPhi of unmixed");
        angularPlot2D(h_dPhi_Z, options, "dPhi of unmixed vs Z");
        // angularPlot2D(h_dPhi_cen, options, "dPhi of unmixed vs MBD centrality");
        // angularPlot3D(h_dPhi_cen, options, "dPhi of unmixed vs MBD centrality 3D");
        angularPlot3D(h_dPhi_Z, options, "dPhi of unmixed vs Z");
    }
    else if (options[1] == "perCenOnOne") {
        std::vector<double> Phi0, Phi1;
        TH1D *h_onOne[14];
        for (int i = 0; i < 14; i++) {
            h_onOne[i] = new TH1D(Form("dPhi of %f to %f", static_cast<double>(i)*0.05, static_cast<double>(i + 1)*0.05), Form("dPhi of %f to %f", static_cast<double>(i)*0.05, static_cast<double>(i + 1)*0.05), N, range_min, range_max);
        }
        for (int i = 0; i < 100; i++) {
            // std::cout << i << std::endl;
            branch16->GetEntry(index[i]);   // ClusPhi;
            branch11->GetEntry(index[i]);   // ClusLayer;
            branch30->GetEntry(index[i]);   // MBD_centrality;
            // if (std::abs(MBD_centrality - MBD_cen[i]) > 1e-07) {
            //     std::cout << MBD_centrality - MBD_cen[i] << ", " << MBD_cen[i] << std::endl;
            // }
            for (int n = 0; n < 14; n++) {
                if (MBD_centrality >= static_cast<double>(n)*0.05 && MBD_centrality <= static_cast<double>(n + 1)*0.05) {
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
        ArrayPlot1D(h, options, "Hi");
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

    if (sub_options[0] == "foundz") INTTZAnalysisLite(filePath, sub_options);
    if (sub_options[0] == "dPhi")   INTTdPhiAnalysis(tree, filePath, sub_options);
    if (sub_options[0] == "dEta")   INTTdEtaAnalysis(tree, filePath, sub_options);

    // Stop the stopwatch and print the runtime
    timer.Stop();
    timer.Print();
}
