#include "globalDefinitions.h"

void save_histo_toRootFile (TH1D const* h, const std::vector<std::string> method, const std::string fileName) {
    TFile* file = TFile::Open(Form("%s.root", fileName.c_str()), "UPDATE");
    if (!file || file -> IsZombie()) {
        std::cerr << "You messed up!" << std::endl;
        exit(1);
    }
    h -> Write(Form("%s", method[0].c_str()), TObject::kOverwrite);
}

void histogramOverflowCheck (TH1D const* h) {
    int underflowBin = 0;
    int overflowBin  = h->GetNbinsX() + 1;
    // Check underflow bin
    double underflow = h->GetBinContent(underflowBin);
    // Check overflow bin
    double overflow  = h->GetBinContent(overflowBin);
    if (underflow > 0) {
        std::cout << "There are " << underflow << " entries in the underflow bin." << std::endl;
    }
    if (overflow > 0) {
        std::cout << "There are " << overflow << " entries in the overflow bin." << std::endl;
    }
    if (underflow == 0 && overflow == 0) {
        std::cout << "No entries outside the range." << std::endl;
    }
}