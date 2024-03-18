{
    #include <TFile.h>
    #include <TTree.h>
    #include <vector>
    #include <cmath>

    /*   retreive data from foundz results text file:   */
    ifstream myfile("../Meeting5/foundZ_debug2_DCA_biased.txt");
    // ifstream myfile("../Meeting5/foundZ_debug2_DCA_-8_8_-010_010_2e-1.txt");
    if (!myfile.is_open()) {
		  cout << "Unable to open text file" << endl;
		  system("read -n 1 -s -p \"Press any key to continue...\" echo");
		  exit(1);
 	}
    vector<double> foundz, truez;
    vector<double> totalsize;
    vector<int>    evt;
    string line, substr;
    while (getline(myfile, line)) {
        stringstream str(line);
        // getline(str, substr, ',');
        getline(str, substr, ',');
        int e = stoi(substr);
        getline(str, substr, ',');
        int Nparticles = stoi(substr);
        getline(str, substr, ',');
        double found_z = stod(substr);
        getline(str, substr, ',');
        double true_z = stod(substr);

        evt.push_back(e);
        foundz.push_back(found_z);
        truez.push_back(true_z);
        totalsize.push_back(Nparticles);
    }

    /*   Read data from .root file:   */
    TFile *f = TFile::Open("/Users/yaminocellist/MIT_mentorship/3rd_semester/INTTRecoClusters_sim_ana382_zvtx-20cm_Bfield0T.root");
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
    // Set up event variables to inherit the data
    std::vector<int>   *ClusLayer = nullptr;
    std::vector<int>   *ClusLadderPhiId = nullptr;
    std::vector<float> *TruthPV_x = nullptr;
    std::vector<float> *TruthPV_y = nullptr;
    std::vector<float> *TruthPV_z = nullptr;
    std::vector<int>   *TruthPV_Npart = nullptr;
    int   event, NTruthVtx;
    float centrality_bimp, centrality_impactparam, centrality_mbdquantity, centrality_mbd;
    EventTree -> SetBranchAddress("event", &event);
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
    EventTree -> SetBranchAddress("ClusLadderPhiId", &ClusLadderPhiId);
    int k = 0;

    /*    Create and write a metric .root file:    */
    // Create a new ROOT file
    TFile *newFile = new TFile("metrics.root", "RECREATE");
    // Create a tree within the file
    TTree *MetricTree = new TTree("MetricTree", "MetricTree");
    // Define variables to hold data for each branch
    Int_t Event, Index, Npart, NpartFromSource;
    Float_t trueZ_vtx, foundZ_vtx, trueXfromSource, trueYfromSource, trueZfromSource;
    Float_t centralityBimp, centralityImpactparam, centralityMbdquantity, centralityMbd;
    // Create branches in the tree:
    MetricTree -> Branch("Event", &Event, "Event/I");
    MetricTree -> Branch("Index", &Index, "Index/I");
    MetricTree -> Branch("foundZ_vtx", &foundZ_vtx, "foundZ_vtx/F");
    MetricTree -> Branch("Npart", &Npart, "Npart/I");
    MetricTree -> Branch("trueXfromSource", &trueXfromSource, "trueXfromSource/F");
    MetricTree -> Branch("trueYfromSource", &trueYfromSource, "trueYfromSource/F");
    MetricTree -> Branch("trueZfromSource", &trueZfromSource, "trueZfromSource/F");
    MetricTree -> Branch("NpartFromSource", &NpartFromSource, "NpartFromSource/I");
    MetricTree -> Branch("centralityBimp", &centralityBimp, "centralityBimp/F");
    MetricTree -> Branch("centralityImpactparam", &centralityImpactparam, "centralityImpactparam/F");
    MetricTree -> Branch("centralityMbdquantity", &centralityMbdquantity, "centralityMbdquantity/F");
    MetricTree -> Branch("centralityMbd", &centralityMbd, "centralityMbd/F");
    // Fill the tree with data
    Long64_t j = 0;

    for (Long64_t i = 0; i < EventTree->GetEntries(); ++i) {
        
        EventTree->GetEntry(i);
        if (NTruthVtx == 1 && TruthPV_Npart->at(0) > 500) {
            while (i != evt[j] || TruthPV_Npart->at(0) != totalsize[j]) {
                // std::cout << i << "," << evt[k] << std::endl;
                i++;    EventTree->GetEntry(i);
            }
            Event = event;
            Index = i;
            // std::cout << i << "," << j << "," << evt[j] << "," << Event << std::endl;
            foundZ_vtx = foundz[j];
            trueXfromSource = TruthPV_x->at(0);
            trueYfromSource = TruthPV_y->at(0);
            trueZfromSource = TruthPV_z->at(0);
            Npart = totalsize[j];
            NpartFromSource = TruthPV_Npart->at(0);
            centralityBimp = centrality_bimp;
            centralityImpactparam = centrality_impactparam;
            centralityMbdquantity = centrality_mbdquantity;
            centralityMbd = centrality_mbd;

            // Fill the tree with the current event data
            MetricTree->Fill();
            // if (std::abs(TruthPV_z->at(0)*10 - truez[k]) > 1e-3 )   std::cout << std::abs(TruthPV_z->at(0)*10 - truez[k]) << std::endl;
            // std::cout << TruthPV_Npart->at(0) << "," << totalsize[k] << "," << TruthPV_z->at(0)*10 << "," << truez[k] << std::endl; 
            j++;
        }
    }

    // Write the tree to the file and close it
    MetricTree->Write();
    newFile->Close();
}