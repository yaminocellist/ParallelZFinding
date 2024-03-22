#include <TFile.h>
#include <TTree.h>
#include <vector>
#include <cmath>
#include "./headers/Analysis.h"
void MultiplicityAnalysis () {}

void trueSingle (Int_t target) {
    // Open the ROOT file
    TFile *f = TFile::Open("../External/INTTRecoClusters_sim_ana382_zvtx-20cm_Bfield0T.root");
    if (!f || f->IsZombie()) {
        std::cerr << "Error opening file" << std::endl;
        return;
    }
    // Get the tree
    TTree *EventTree;
    f->GetObject("EventTree", EventTree);
    if (!EventTree) {
        std::cerr << "Tree not found" << std::endl;
        return;
    }
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
    int NTruthVtx;
    float centrality_mbd;
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
    std::cout << target << std::endl;
    for (int i = 0; i < EventTree -> GetEntries(); ++i) {
        EventTree -> GetEntry(i);
        if (i == target && NTruthVtx == 1 && TruthPV_Npart->at(0) > 500
            && centrality_mbd <= 70 && TruthPV_z->at(0) >= -25. && TruthPV_z->at(0) <= -15.) {
            std::cout << ClusX->size() << "," << UniqueAncG4P_TrackID->size() << std::endl;
            for (int j = 0; j < ClusX->size(); j++) {
                // std::cout << "Hi" << j << std::endl;
                double z = ClusZ->at(j) - TruthPV_z->at(0);
                double r = std::sqrt(std::pow(ClusX->at(j), 2) + std::pow(ClusY->at(j), 2));
                double theta = std::atan2(r, z);
                double eta;
                if (z >= 0) {
                    eta = -std::log(std::tan(theta/2));
                } else {
                    eta = std::log(std::tan((M_PI - theta)/2));
                }
                if (ClusLayer->at(j) == 3 || ClusLayer->at(j) == 4) {
                    tracklet_layer_0.push_back({ClusX->at(j), ClusY->at(j), ClusZ->at(j), r, eta, std::atan2(ClusY->at(j), ClusX->at(j)), 0});
                } else {
                    tracklet_layer_1.push_back({ClusX->at(j), ClusY->at(j), ClusZ->at(j), r, eta, std::atan2(ClusY->at(j), ClusX->at(j)), 1});
                }
            }
            break;
        }
    }
    //  First, calculate the dη:
    double dEta = EtaCheck(target, tracklet_layer_0, tracklet_layer_1, -0.001, 0.001, -0.01, 0.01, TruthPV_z->at(0));
    std::cout << dEta << std::endl;
}