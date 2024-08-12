#include "globalDefinitions.h"


void current_PC_time () {
    auto now = std::chrono::system_clock::now();
    auto now_c = std::chrono::system_clock::to_time_t(now);
    std::cout << "Run started at: " << std::put_time(std::localtime(&now_c), "%Y-%m-%d %H:%M:%S") << std::endl;
}

void save_histo_toRootFile (TH1D const* h, const std::vector<std::string> method, const std::string fileName) {
    TFile* file = TFile::Open(Form("%s.root", fileName.c_str()), "UPDATE");
    if (!file || file -> IsZombie()) {
        std::cerr << "You messed up!" << std::endl;
        exit(1);
    }
    h -> Write(Form("%s", method[0].c_str()), TObject::kOverwrite);
}