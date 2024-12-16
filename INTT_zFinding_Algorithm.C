#include "headers/zFinding.h"
#include "headers/Analysis.h"

double singleEventZFinding (TTree *EventTree, Int_t target, std::vector<string> options, double lower_dPhi_cut, double upper_dPhi_cut) {
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

    branch25->GetEntry(target);  branch10->GetEntry(target);
    branch11->GetEntry(target);  branch12->GetEntry(target);  branch13->GetEntry(target);  branch14->GetEntry(target);
    branch15->GetEntry(target);  branch16->GetEntry(target);  branch17->GetEntry(target);
    branch30->GetEntry(target);  branch31->GetEntry(target);

    double found_z;
    if (MBD_z_vtx >= -25. & MBD_z_vtx <= -5.) {
        std::vector<myTrackletMemberLite> tracklet_layer_0, tracklet_layer_1;
        for (int i = 0; i < ClusX->size(); i++) {
            if (ClusLayer->at(i) == 3 || ClusLayer->at(i) == 4) {
                // std::cout << std::atan2(ClusY->at(i), ClusX->at(i)) - ClusPhi->at(i) << std::endl;
                for (int k = 0; k <= 32; k++) {
                    tracklet_layer_0.push_back({ClusX->at(i), ClusY->at(i), ClusZ->at(i) + k*0.05,
                                                ClusPhi->at(i)});
                }
            }
            else {
                for (int k = 0; k <= 32; k++) {
                    tracklet_layer_1.push_back({ClusX->at(i), ClusY->at(i), ClusZ->at(i) + k*0.05,
                                                ClusPhi->at(i)});
                }
            }
        }
        found_z = DCA_npeaks_fit(event25, tracklet_layer_0, tracklet_layer_1, -0.001, 0.001, lower_dPhi_cut, upper_dPhi_cut, MBD_z_vtx);
        tracklet_layer_0.clear();   tracklet_layer_1.clear();
    }
    else {
        std::cout << "Event " << event25 << " is NOT the event we need." << std::endl;
        exit(1);
    }
    std::cout << found_z << ",   " << MBD_z_vtx << std::endl;
    return found_z;
}

void allEventZFinding (TTree *EventTree, Int_t target, std::vector<string> options, double lower_dPhi_cut, double upper_dPhi_cut, const std::string &filePath) {
    std::ofstream outputFile(filePath, std::ios_base::app);
    if (!outputFile.is_open()) {
		std::cout << "Unable to open the file to be written." << std::endl;
		system("read -n 1 -s -p \"Press any key to continue...\" echo");
		exit(1);
	}
    std::ifstream checkFile(filePath, std::ios::ate); // Check if the text file is empty;
    if (!checkFile.is_open()) {
        std::cerr << "Failed to open file for reading.\n";
        system("read -n 1 -s -p \"Press any key to continue...\" echo");
        exit(1);
    }
    bool isEmpty = (checkFile.tellg() == 0);
    if(isEmpty) outputFile << "index,event,NClus,foundZ,trueZ,centrality" << std::endl;
    
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

    double found_z;
    std::vector<myTrackletMemberLite> tracklet_layer_0, tracklet_layer_1;
    for (Int_t i = target; i < nEntries; i++) {
        branch25->GetEntry(i);  branch10->GetEntry(i);
        branch11->GetEntry(i);  branch12->GetEntry(i);  branch13->GetEntry(i);  branch14->GetEntry(i);
        branch15->GetEntry(i);  branch16->GetEntry(i);  branch17->GetEntry(i);
        branch30->GetEntry(i);  branch31->GetEntry(i);
        if (MBD_z_vtx >= -25. & MBD_z_vtx <= -5. && MBD_centrality <= 0.7) {
            for (int j = 0; j < ClusX->size(); j++) {
                if (ClusLayer->at(j) == 3 || ClusLayer->at(j) == 4) {
                    // std::cout << std::atan2(ClusY->at(i), ClusX->at(i)) - ClusPhi->at(i) << std::endl;
                    for (int k = 0; k <= 32; k++) {
                        tracklet_layer_0.push_back({ClusX->at(j), ClusY->at(j), ClusZ->at(j) + k*0.05,
                                                    ClusPhi->at(j)});
                    }
                } 
                else {
                    for (int k = 0; k <= 32; k++) {
                        tracklet_layer_1.push_back({ClusX->at(j), ClusY->at(j), ClusZ->at(j) + k*0.05,
                                                    ClusPhi->at(j)});
                    }
                }
            }
            found_z = DCA_npeaks_fitLite(event25, tracklet_layer_0, tracklet_layer_1, -0.001, 0.001, lower_dPhi_cut, upper_dPhi_cut, MBD_z_vtx);
            tracklet_layer_0.clear();   tracklet_layer_1.clear();
            outputFile << i << "," << event25 << "," << NClus << "," << found_z << "," << MBD_z_vtx << "," << MBD_centrality << std::endl;
        }
    }
}

void INTT_zFinding_Algorithm (std::string method = "single", double lower_dPhi_cut = -0.1, double upper_dPhi_cut = 0.1, Int_t target = 0) {
    std::string filePath = "zFindingResults/new_foundZ_DCAfit_0_16_-001_001_step2e-1_z_25_5.txt";
    
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

    std::vector<std::string> sub_options = splitString(method, '_');
    if (sub_options.size() == 2) {
        sub_options.push_back("");
    }
    else if (sub_options.size() == 1) {
        sub_options.push_back("");
        sub_options.push_back("");
    }
    if (sub_options[0] == "single")    double rst = singleEventZFinding(tree, target, sub_options, lower_dPhi_cut, upper_dPhi_cut);
    if (sub_options[0] == "all")       allEventZFinding(tree, target, sub_options, lower_dPhi_cut, upper_dPhi_cut, filePath);
}