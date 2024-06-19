#include "./headers/Analysis.h"
#include "./headers/xyFinding.h"

void single_xyFinding (TTree *EventTree, Int_t target, std::vector<std::string> method, const std::string &filePath) {
    std::ifstream myfile(filePath);
    if (!myfile.is_open()) {
		  std::cout << "Unable to open text file" << std::endl;
		  system("read -n 1 -s -p \"Press any key to continue...\" echo");
		  exit(1);
 	}
    
    std::vector<double> foundz, truez;
    std::vector<double> totalsize;
    std::vector<int>    evt, idx;
    std::string line, substr;
    getline(myfile, line);
    while (getline(myfile, line)) {
        stringstream str(line);
        getline(str, substr, ',');  int    e          = stoi(substr);
        getline(str, substr, ',');  int    index      = stoi(substr);
        getline(str, substr, ',');  int    Nparticles = stoi(substr);
        getline(str, substr, ',');  double found_z    = stod(substr);
        getline(str, substr, ',');  double true_x     = stod(substr);
        getline(str, substr, ',');  double true_y     = stod(substr);
        getline(str, substr, ',');  double true_z     = stod(substr);

        evt.push_back(e);   idx.push_back(index);
        foundz.push_back(found_z);
        truez.push_back(true_z);
        totalsize.push_back(Nparticles);
    }

    // Set up event variables to inherit the data
    std::vector<float> *ClusX         = nullptr;
    std::vector<float> *ClusY         = nullptr;
    std::vector<float> *ClusZ         = nullptr;
    std::vector<int>   *ClusLayer     = nullptr;
    std::vector<float> *TruthPV_x     = nullptr;
    std::vector<float> *TruthPV_y     = nullptr;
    std::vector<float> *TruthPV_z     = nullptr;
    std::vector<int>   *TruthPV_Npart = nullptr;
    int   event, NTruthVtx;
    float centrality_bimp, centrality_impactparam, centrality_mbdquantity, centrality_mbd;
    EventTree -> SetBranchAddress("event", &event);
    EventTree -> SetBranchAddress("ClusX", &ClusX);
    EventTree -> SetBranchAddress("ClusY", &ClusY);
    EventTree -> SetBranchAddress("ClusZ", &ClusZ);
    EventTree -> SetBranchAddress("NTruthVtx", &NTruthVtx);
    EventTree -> SetBranchAddress("centrality_bimp", &centrality_bimp);
    EventTree -> SetBranchAddress("centrality_impactparam", &centrality_impactparam);
    EventTree -> SetBranchAddress("centrality_mbdquantity", &centrality_mbdquantity);
    EventTree -> SetBranchAddress("centrality_mbd", &centrality_mbd);
    EventTree -> SetBranchAddress("ClusLayer", &ClusLayer);
    EventTree -> SetBranchAddress("TruthPV_Npart", &TruthPV_Npart);
    EventTree -> SetBranchAddress("TruthPV_x", &TruthPV_x);
    EventTree -> SetBranchAddress("TruthPV_y", &TruthPV_y);
    EventTree -> SetBranchAddress("TruthPV_z", &TruthPV_z);

    Long64_t j = 0;
    // "ad" stand for Angular distance;
    // https://particle.wiki/wiki/Angular_distance
    std::vector<myTrackletMember> t0_for_ad, t1_for_ad, tracklet_layer_0, tracklet_layer_1;
    double r, currentZ, dZ, theta, eta, phi;   // intermediate variables;
    for (Long64_t i = 0; i < EventTree->GetEntries(); ++i) {
        // EventTree->GetEntry(i);
        EventTree -> GetEntry(idx[i]);
        // if (NTruthVtx == 1 && TruthPV_Npart->at(0) > 500 && j < evt.size() &&
        //     centrality_mbd <= 70 && TruthPV_z->at(0) >= -25. && TruthPV_z->at(0) <= -15.) {
            // while (evt[j] != event || i != idx[j] || TruthPV_Npart->at(0) != totalsize[j]) {
            //     i++;    EventTree->GetEntry(i);
            // }
            // if (i == target) {
            if (idx[i] == target && NTruthVtx == 1 && TruthPV_Npart->at(0) > 500 && j < evt.size() &&
                centrality_mbd <= 70 && TruthPV_z->at(0) >= -25. && TruthPV_z->at(0) <= -15.) {
                std::cout << evt[i] << ", " << event << std::endl;
                std::cout << truez[i] << ", " << TruthPV_z->at(0) << std::endl;
                std::cout << totalsize[i] << ", " << TruthPV_Npart->at(0) << std::endl;
                for (int k = 0; k < ClusX->size(); k++) {
                    // r   = std::sqrt(std::pow(ClusX->at(k), 2) + std::pow(ClusY->at(k), 2));
                    r   = std::sqrt((ClusX->at(k) - TruthPV_x->at(0))*(ClusX->at(k) - TruthPV_x->at(0)) + (ClusY->at(k) - TruthPV_y->at(0))*(ClusY->at(k) - TruthPV_y->at(0)));
                    phi = std::atan2(ClusY->at(k), ClusX->at(k));
                    currentZ         = ClusZ->at(k);
                    dZ               = currentZ - foundz[j];
                    // dZ               = currentZ - TruthPV_z->at(0);
                    theta            = std::atan2(r, dZ);
                    if (dZ >= 0) {
                        eta = -log(tan(theta/2));
                    } else {
                        eta = log(tan((M_PI - theta)/2));
                    }         
                    if (ClusLayer->at(k) == 3 || ClusLayer->at(k) == 4) {
                        t0_for_ad.push_back({ClusX->at(k), ClusY->at(k), currentZ, r,
                                                    eta, phi, 
                                                    0});
                        tracklet_layer_0.push_back({ClusX->at(k), ClusY->at(k), currentZ, r,
                                                    eta, phi, 
                                                    0});
                    } else {
                        t1_for_ad.push_back({ClusX->at(k), ClusY->at(k), currentZ, r,
                                                    eta, phi, 
                                                    1});
                        tracklet_layer_1.push_back({ClusX->at(k), ClusY->at(k), currentZ, r,
                                                    eta, phi, 
                                                    1});
                    }
                    //
                    for (int l = 1; l <= 32; l++) {
                        currentZ         = ClusZ->at(k) + l*0.05;
                        dZ               = currentZ - foundz[j];
                        // dZ               = currentZ - TruthPV_z->at(0);
                        theta            = std::atan2(r, dZ);
                        if (dZ >= 0) {
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
                }
                std::cout << evt[j] << "," << event << "," << idx[j] << "," << i << "," 
                          << foundz[j] << "," << TruthPV_z->at(0) << "," << TruthPV_x->at(0) << "," << TruthPV_y->at(0)
                          << std::endl;
                double adCut = angularDistance(evt[j], t0_for_ad, t1_for_ad);
                std::pair<double, double> found_xy = findXY(evt[j], tracklet_layer_0, tracklet_layer_1, adCut*adCut, TruthPV_z->at(0));
                t0_for_ad.clear(); t1_for_ad.clear(); tracklet_layer_0.clear(); tracklet_layer_1.clear();
                break;
            }
            j++;
        }
    // }
}

void all_xyFinding (TTree *EventTree, Int_t target, std::vector<std::string> method, const std::string &filePath, const std::string &saveFile) {
    std::ifstream myfile(filePath);
    if (!myfile.is_open()) {
		std::cout << "Unable to open text file" << std::endl;
		system("read -n 1 -s -p \"Press any key to continue...\" echo");
		exit(1);
 	}
    std::ofstream outputFile(saveFile, std::ios_base::app);
    if (!outputFile.is_open()) {
		std::cout << "Unable to open the file to be written." << std::endl;
		system("read -n 1 -s -p \"Press any key to continue...\" echo");
		exit(1);
	}
    
    std::vector<double> foundz, truez;
    std::vector<double> totalsize;
    std::vector<int>    evt, idx;
    string line, substr;
    getline(myfile, line);
    while (getline(myfile, line)) {
        stringstream str(line);
        getline(str, substr, ',');  int    e          = stoi(substr);
        getline(str, substr, ',');  int    index      = stoi(substr);
        getline(str, substr, ',');  int    Nparticles = stoi(substr);
        getline(str, substr, ',');  double found_z    = stod(substr);
        getline(str, substr, ',');  double true_x     = stod(substr);
        getline(str, substr, ',');  double true_y     = stod(substr);
        getline(str, substr, ',');  double true_z     = stod(substr);

        evt.push_back(e);   idx.push_back(index);
        foundz.push_back(found_z);
        truez.push_back(true_z);
        totalsize.push_back(Nparticles);
    }

    // Set up event variables to inherit the data
    std::vector<float> *ClusX = nullptr;
    std::vector<float> *ClusY = nullptr;
    std::vector<float> *ClusZ = nullptr;
    std::vector<int>   *ClusLayer = nullptr;
    std::vector<float> *TruthPV_x = nullptr;
    std::vector<float> *TruthPV_y = nullptr;
    std::vector<float> *TruthPV_z = nullptr;
    std::vector<int>   *TruthPV_Npart = nullptr;
    int   event, NTruthVtx;
    float centrality_bimp, centrality_impactparam, centrality_mbdquantity, centrality_mbd;
    EventTree -> SetBranchAddress("event", &event);
    EventTree -> SetBranchAddress("ClusX", &ClusX);
    EventTree -> SetBranchAddress("ClusY", &ClusY);
    EventTree -> SetBranchAddress("ClusZ", &ClusZ);
    EventTree -> SetBranchAddress("NTruthVtx", &NTruthVtx);
    EventTree -> SetBranchAddress("centrality_bimp", &centrality_bimp);
    EventTree -> SetBranchAddress("centrality_impactparam", &centrality_impactparam);
    EventTree -> SetBranchAddress("centrality_mbdquantity", &centrality_mbdquantity);
    EventTree -> SetBranchAddress("centrality_mbd", &centrality_mbd);
    EventTree -> SetBranchAddress("ClusLayer", &ClusLayer);
    EventTree -> SetBranchAddress("TruthPV_Npart", &TruthPV_Npart);
    EventTree -> SetBranchAddress("TruthPV_x", &TruthPV_x);
    EventTree -> SetBranchAddress("TruthPV_y", &TruthPV_y);
    EventTree -> SetBranchAddress("TruthPV_z", &TruthPV_z);

    Long64_t j = 0;
    std::vector<myTrackletMember> t0_for_ad, t1_for_ad, tracklet_layer_0, tracklet_layer_1;
    double r, currentZ, dZ, theta, eta, phi;   // intermediate variables;
    for (Long64_t i = 0; i < EventTree->GetEntries(); ++i) {
        EventTree->GetEntry(i);
        if (NTruthVtx == 1 && TruthPV_Npart->at(0) > 500 && j < evt.size() &&
            centrality_mbd <= 70 && TruthPV_z->at(0) >= -25. && TruthPV_z->at(0) <= -15.) {
            while (evt[j] != event || i != idx[j] || TruthPV_Npart->at(0) != totalsize[j]) {
                i++;    EventTree->GetEntry(i);
            }
            if (i >= target) {
                
                for (int k = 0; k < ClusX->size(); k++) {
                    // r   = std::sqrt(std::pow(ClusX->at(k), 2) + std::pow(ClusY->at(k), 2));
                    r   = std::sqrt((ClusX->at(k) - TruthPV_x->at(0))*(ClusX->at(k) - TruthPV_x->at(0)) + (ClusY->at(k) - TruthPV_y->at(0))*(ClusY->at(k) - TruthPV_y->at(0)));
                    phi = std::atan2(ClusY->at(k), ClusX->at(k));
                    // 
                        currentZ         = ClusZ->at(k);
                        // dZ               = currentZ - foundz[j];
                        dZ               = currentZ - TruthPV_z->at(0);
                        theta            = std::atan2(r, dZ);
                        if (dZ >= 0) {
                            eta = -log(tan(theta/2));
                        } else {
                            eta = log(tan((M_PI - theta)/2));
                        }         
                        if (ClusLayer->at(k) == 3 || ClusLayer->at(k) == 4) {
                            t0_for_ad.push_back({ClusX->at(k), ClusY->at(k), currentZ, r,
                                                        eta, phi, 
                                                        0});
                            tracklet_layer_0.push_back({ClusX->at(k), ClusY->at(k), currentZ, r,
                                                        eta, phi, 
                                                        0});
                        } else {
                            t1_for_ad.push_back({ClusX->at(k), ClusY->at(k), currentZ, r,
                                                        eta, phi, 
                                                        1});
                            tracklet_layer_1.push_back({ClusX->at(k), ClusY->at(k), currentZ, r,
                                                        eta, phi, 
                                                        1});
                        }
                    //
                    for (int l = 1; l <= 32; l++) {
                        currentZ         = ClusZ->at(k) + l*0.05;
                        // dZ               = currentZ - foundz[j];
                        dZ               = currentZ - TruthPV_z->at(0);
                        theta            = std::atan2(r, dZ);
                        if (dZ >= 0) {
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
                }
                // /*
                double adCut = angularDistance(evt[j], t0_for_ad, t1_for_ad);
                std::pair<double, double> found_xy = findXY(evt[j], tracklet_layer_0, tracklet_layer_1, adCut*adCut, foundz[j]);
                outputFile << event << "," << i << "," << TruthPV_Npart->at(0) << "," << 
                            std::fixed << std::setprecision(6) << found_xy.first << "," << found_xy.second << "," << foundz[j] << std::defaultfloat << "," << 
                            TruthPV_x->at(0) << "," << TruthPV_y->at(0) << "," << TruthPV_z->at(0) << "," 
                            << centrality_bimp << "," << centrality_impactparam << "," << centrality_mbdquantity << "," << centrality_mbd << std::endl;

                // */
                t0_for_ad.clear(); t1_for_ad.clear(); tracklet_layer_0.clear(); tracklet_layer_1.clear();
                
            //    std::cout << evt[j] << "," << event << "," << idx[j] << "," << i << "," 
            //               << foundz[j] << "," << TruthPV_z->at(0) << "," << TruthPV_x->at(0) << "," << TruthPV_y->at(0)
            //               << std::endl;
            }
            j++;
        }
    }
}

void INTT_xyFinding (std::string method = "single", Int_t target = 0) {
    // foundZ data:
    std::string filePath = "./zFindingResults/foundZ_MBD0_10_DCAfit_0_16_-001_001_2e-1.txt";
    // foundXY results, whose name is alighend with foundZ file:
    std::string saveFile = "../External/xyFindingResults/foundXY_MBD0_10_DCAfit_0_16_-001_001_2e-1.txt";

    /*   Data from original .root file:   */
    TFile *f = TFile::Open("../External/INTTRecoClusters_sim_ana382_zvtx-20cm_Bfield0T.root");
    if (!f || f->IsZombie()) {
        std::cerr << "Error opening file" << std::endl;
        return;
    }
    TTree *EventTree;
    f->GetObject("EventTree", EventTree);
    if (!EventTree) {
        std::cerr << "Tree not found" << std::endl;
        return;
    }
    
    std::vector<std::string> substrings = splitString(method, '_');
    if (substrings.size() == 2) {
        substrings.push_back("");
    }
    else if (substrings.size() == 1) {
        substrings.push_back("");
        substrings.push_back("");
    }

    if        (substrings[0] == "single") {
        single_xyFinding(EventTree, target, substrings, filePath);
    } else if (substrings[0] == "all") {
        all_xyFinding(EventTree, target, substrings, filePath, saveFile);
    } else if (substrings[0] == "vtxAnal") {
        foundXYAnalysis(EventTree, saveFile, substrings[1], substrings[2]);
    } else if (substrings[0] == "etaPhiAnal") {
        foundEtaPhiAnalysis(EventTree, saveFile, target, substrings);
    }
    else if   (substrings[0] == "dPhiInZ") {
        dPhiInZVtx(EventTree, saveFile, target, substrings);
    }
    else {
        // etaPhiAnalysis(EventTree, target, substrings, lower_bound, upper_bound);
    }
}