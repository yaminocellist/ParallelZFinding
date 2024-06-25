#include "headers/Analysis.h"

void ResultsAnalysis (std::string method) {
    std::string filePath = "zFindingResults/new_foundZ_DCAfit_0_16_-001_001_step2e-1_z_25_5.txt";
    std::vector<std::string> sub_options = splitString(method, '_');
    if (sub_options.size() == 2) {
        sub_options.push_back("");
    }
    else if (sub_options.size() == 1) {
        sub_options.push_back("");
        sub_options.push_back("");
    }

    if (sub_options[0] == "foundz") foundZAnalysisLite(filePath, sub_options);   
}
