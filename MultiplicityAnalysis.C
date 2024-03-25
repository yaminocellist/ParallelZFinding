#include <TFile.h>
#include <TTree.h>
#include <vector>
#include <cmath>
#include "./headers/Analysis.h"
#include "./headers/zFinding.h"
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
    for (int i = 0; i < EventTree -> GetEntries(); ++i) {
        EventTree -> GetEntry(i);
        if (i == target && NTruthVtx == 1 && TruthPV_Npart->at(0) > 500
            && centrality_mbd <= 70 && TruthPV_z->at(0) >= -25. && TruthPV_z->at(0) <= -15.) {
            std::cout << ClusX->size() << "," << UniqueAncG4P_TrackID->size() << std::endl;
            for (int j = 0; j < ClusX->size(); j++) {
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
    //  1st, calculate the dη:
    double dEta = EtaCheck(target, tracklet_layer_0, tracklet_layer_1, -0.001, 0.001, -0.01, 0.01, TruthPV_z->at(0));
    std::cout << dEta << std::endl;

    //  2nd, calculate the dN_ch ('ch' means charged particle):
        //  Use DC (distance closest) to distinguish signal/background;
        //  the same value as doing z-finding;
    // std::cout << tracklet_layer_0.size() << "," << tracklet_layer_1.size() << std::endl;
    TH1D *HSignal     = new TH1D("HSignal","HSignal; Delta_Eta; Hits",400,-2,2);
    TH1D *HBackground = new TH1D("HBackground","HBackground;Delta_Eta; Hits",400,-2,2);
    for (int i = 0; i < tracklet_layer_0.size(); i++) {
        for (int j = 0; j < tracklet_layer_1.size(); j++) {
            myPoint3D P1 = {tracklet_layer_0[i].x, tracklet_layer_0[i].y, tracklet_layer_0[i].z};
            myPoint3D P2 = {tracklet_layer_1[j].x, tracklet_layer_1[j].y, tracklet_layer_1[j].z};
            double DC = nearestZ(P1, P2).second;
            double d_phi = tracklet_layer_0[i].phi - tracklet_layer_1[j].phi;
            if (DC <= DCA_cut) {
                HSignal -> Fill(d_phi);
            } else {
                HBackground -> Fill(d_phi);
            }
        }
    }
    double N1 = HBackground -> Integral(HBackground->FindFixBin(-2), HBackground->FindFixBin(-0.1), "") + 
                HBackground -> Integral(HBackground->FindFixBin(0.1), HBackground->FindFixBin(2), "");
    double N2 = HSignal -> Integral(HSignal->FindFixBin(-2), HSignal->FindFixBin(-0.1), "") + 
                HSignal -> Integral(HSignal->FindFixBin(0.1), HSignal->FindFixBin(2), "");
    double N  = N2/N1;
    
    TCanvas *c1 = new TCanvas("c1", "Analysis Results", 10, 10, 900, 600);
    c1 -> Divide(1, 2); // This divides the canvas into two parts
    c1 -> cd(1);
    HSignal -> Sumw2();
    HSignal -> Draw("HIST SAME");
    HSignal -> Draw("e1psame");
    HSignal -> SetFillColor(kRed-10);
    HSignal -> SetFillStyle(1001);
    // HSignal -> SetTitle(Form("HSignal & Normalized Background for %5.0f Events", events));
    // HSignal -> GetXaxis() -> SetRange(HSignal -> FindFixBin(-0.5),HSignal -> FindFixBin(0.5));
    HSignal -> SetFillColor(kCyan - 7);
    HSignal -> SetMarkerSize(0.55);
    HSignal -> SetFillStyle(1001);
    // gPad -> SetLogy(1);

    TH1D* HNormalized = (TH1D*) HBackground->Clone("Normalized HBackground");
    HNormalized -> Add(HBackground, N-1);
    HNormalized -> Sumw2();
    // HNormalized -> Draw("HIST SAME");
    // HNormalized -> Draw("e1psame");

    HNormalized -> SetFillColor(kRed-7);
    HNormalized -> SetFillStyle(1001);
    //HNormalized -> GetXaxis() -> SetRange(HNormalized -> FindFixBin(-0.5),HNormalized -> FindFixBin(0.5));
    // HNormalized -> SetTitle(Form("Normalized HBackground Between 0.1 and 0.2 for %5.0f Events", events));
    c1 -> cd(2);
    TH1D* hDiff = (TH1D*) HSignal->Clone("Subtracted Signal");
    hDiff -> Add(HNormalized, -1);
    hDiff -> Sumw2();
    hDiff -> Draw("HIST SAME");
    hDiff -> Draw("e1psame");

    hDiff -> SetFillColor(kYellow - 7);
    hDiff -> SetFillStyle(1001);
    //hDiff -> GetXaxis() -> SetRange(hDiff -> FindFixBin(-0.5),hDiff -> FindFixBin(0.5));
    double peak = hDiff -> Integral(hDiff->FindFixBin(-0.1), hDiff->FindFixBin(0.1), "");
    // double Ratio = peak/events;
   
    
    // hDiff -> SetTitle(Form("Subtracted Signal for %5.0f Events", events));

    gStyle->SetEndErrorSize(6);
    gStyle->SetErrorX(0.5);
    
    gStyle -> SetTitleFont(100,"t");
    gStyle -> SetTitleSize(0.065,"t");

    //gPad -> SetLogy(1);
    TLine *l = new TLine(-2,0,2,0);
	l ->Draw("same"); 
    l -> SetLineColor(kGreen);
}