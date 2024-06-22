#include <TFile.h>
#include <TTree.h>
#include <iostream>

void readTree() {
    // Open the ROOT file
    TFile *file = TFile::Open("path/to/your/file.root");
    if (!file || file->IsZombie()) {
        std::cerr << "Error opening file!" << std::endl;
        return;
    }

    // Get the TTree
    TTree *tree = (TTree*)file->Get("treeName"); // Replace "treeName" with the actual name of your TTree
    if (!tree) {
        std::cerr << "Error getting TTree!" << std::endl;
        return;
    }

    // Define the variables to hold the data
    int event0, event25;

    // Get the branches directly by their indices
    TBranch *branch0 = (TBranch*)tree->GetListOfBranches()->At(0);
    TBranch *branch25 = (TBranch*)tree->GetListOfBranches()->At(25);

    if (!branch0 || !branch25) {
        std::cerr << "Error getting branches!" << std::endl;
        return;
    }

    // Set the branch addresses manually
    branch0->SetAddress(&event0);
    branch25->SetAddress(&event25);

    // Loop over the entries in the TTree and print the data
    Long64_t nEntries = tree->GetEntries();
    for (Long64_t i = 0; i < nEntries; ++i) {
        branch0->GetEntry(i);
        branch25->GetEntry(i);

        std::cout << "Entry " << i << ":" << std::endl;
        std::cout << "event0: " << event0 << std::endl;
        std::cout << "event25: " << event25 << std::endl;
    }

    // Close the file
    file->Close();
}

void RootFileChecker () {
    
}