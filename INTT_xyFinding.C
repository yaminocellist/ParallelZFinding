#include <TFile.h>
#include <TTree.h>
#include <vector>
#include <cmath>
void single_xyFinding (TTree *EventTree, Int_t target, std::vector<std::string> method, const std::string &filePath) {
    std::ifstream myfile(filePath);
    if (!myfile.is_open()) {
		  cout << "Unable to open text file" << endl;
		  system("read -n 1 -s -p \"Press any key to continue...\" echo");
		  exit(1);
 	}
    
    vector<double> foundz, truez;
    vector<double> totalsize;
    vector<int>    evt, idx;
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
    // std::cout << evt.size() << std::endl;
    std::vector<myTrackletMember> tracklet_layer_0, tracklet_layer_1;
    for (Long64_t i = 0; i < EventTree->GetEntries(); ++i) {
        EventTree->GetEntry(i);
        if (NTruthVtx == 1 && TruthPV_Npart->at(0) > 500 && j < evt.size()) {
            while (evt[j] != event || i != idx[j] || TruthPV_Npart->at(0) != totalsize[j]) {
                // std::cout << i << "," << idx[j] << std::endl;
                i++;    EventTree->GetEntry(i);
            }
            if (i == target) {
                for (int k = 0; k < ClusX->size(); k++) {
                if (ClusLayer->at(k) == 3 || ClusLayer->at(k) == 4) {
                    for (int l = 0; l <= 32; l++) {
                        tracklet_layer_0.push_back({ClusX->at(k), ClusY->at(k), ClusZ->at(k) + l*0.05, std::sqrt(std::pow(ClusX->at(k), 2) + std::pow(ClusY->at(k), 2)), INFINITY, std::atan2(ClusY->at(k), ClusX->at(k)), 0});
                    }
                } else {
                    for (int l = 0; l <= 32; l++) {
                        tracklet_layer_1.push_back({ClusX->at(k), ClusY->at(k), ClusZ->at(k) + l*0.05, std::sqrt(std::pow(ClusX->at(k), 2) + std::pow(ClusY->at(k), 2)), INFINITY, std::atan2(ClusY->at(k), ClusX->at(k)), 1});
                    }
                }
            }
            std::cout << evt[j] << "," << event << "," << idx[j] << "," << i << "," 
                      << foundz[j] << "," << TruthPV_z->at(0) << "," << tracklet_layer_0[0].layer << "," << tracklet_layer_1[0].layer
                      << std::endl;
            }


            

            j++;
        }
    }
}

void INTT_xyFinding (std::string method = "single", Int_t target = 0) {
    // foundZ data:
    std::string filePath = "./zFindingResults/foundZ_allInfo_DCA_0_16_-001_001_2e-1.txt"；

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
}