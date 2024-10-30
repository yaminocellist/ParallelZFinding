#ifndef GLOBAL_DEFINITIONS_H
#define GLOBAL_DEFINITIONS_H
#include <vector>
#include <cmath>
#include <string>
#include <memory>
#include <iostream>
#include <fstream>
#include <functional>
#include <chrono>

#include <TH1.h>
#include <TH2F.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TPolyMarker.h>
#include <TAxis.h>
#include <TFile.h>
#include <TTree.h>
#include <TMultiGraph.h>


/*********************************************************************
 *                      GLOBAL VARIABLES;
 * ******************************************************************/
double zmin = -45;
double zmax = 5;        // From sPHENIX paper, the stave's length is aroung 27.12 cm;
double scanstep = 0.2;  // unit: cm; 
int bins = (zmax - zmin)/scanstep + 1;
double dPhi_cut       = 0.01;
double dEta_cut       = 0.2;
double DCA_cut        = 0.2;    // unit: cm;
double DCA_cutSQUARED = 0.04;   // unit: cm;
double MBD_lower = 0., MBD_upper = 10.;
double abs_fit_range = 0.2;

 // ANSI escape code for red text
const std::string COLOR_GREEN="\033[0;32m";
const std::string COLOR_RED="\033[0;31m";
const std::string COLOR_YELLOW="\033[0;33m";
const std::string COLOR_OCHRE="\033[38;5;95m";
const std::string COLOR_BLUE="\033[0;34m";
const std::string COLOR_WHITE="\033[0;37m";
const std::string COLOR_RESET="\033[0m";

/*  Genearting unequal bin ranges for a historgram  */
// Double_t binEdges[481];
// double dx_1 = 0.01;
// binEdges[0] = -2;
// for (int i = 1; i <= 190; i++) {
//     binEdges[i] = binEdges[0] + i*dx_1;
// }
// double dx_2 = 0.002;
// for (int j = 1; j <= 100; j++) {
//     binEdges[j+190] = binEdges[190] + j*dx_2;
// }
// for (int l = 1; l <= 190; l++) {
//     binEdges[l + 290] = binEdges[290] + l*dx_1;
// }

//////////////////////////////////////////////////////////////////////

struct myPoint3D {
    double x, y, z;
};

struct myTrackletMemberLite {
    Double_t x, y, z;
    Double_t phi;
};

struct myTrackletMember {
    Double_t x, y, z, r;
    Double_t eta, phi;
    Int_t layer;
};

struct myTrackletMemberExtended : myTrackletMember {
    Int_t trackID;
};

struct EtaWithPhi {
    double eta_value;
    double phi_value;

    // You can add constructors, methods, etc., to enhance functionality
    EtaWithPhi(double e, double p) : eta_value(e), phi_value(p) {}
};

std::pair<double, double> nearestZ (const myPoint3D &p1, const myPoint3D &p2) {
    double numeratorZ = (p1.y - p2.y)*(p1.y*p2.z - p2.y*p1.z) - (p1.x - p2.x)*(p2.x*p1.z - p1.x*p2.z);
    double denominatorZ = (p1.x - p2.x)*(p1.x - p2.x) + (p1.y - p2.y)*(p1.y - p2.y);
    double nearest_z = numeratorZ/denominatorZ;

    double part123 = (p1.x*p1.x + p2.x*p2.x + p1.y*p1.y + p2.y*p2.y - 2*p1.x*p2.x - 2*p1.y*p2.y)*nearest_z*nearest_z +
                    2*((p1.x*p2.x + p1.y*p2.y - p2.x*p2.x - p2.y*p2.y)*p1.z + (p1.x*p2.x + p1.y*p2.y - p1.x*p1.x - p1.y*p1.y)*p2.z)*nearest_z +
                    p1.y*p1.y*p2.z*p2.z + p2.y*p2.y*p1.z*p1.z + p1.x*p1.x*p2.z*p2.z + p2.x*p2.x*p1.z*p1.z + p1.x*p1.x*p2.y*p2.y + p2.x*p2.x*p1.y*p1.y -
                    2*p1.x*p2.x*p1.y*p2.y - 2*p1.y*p2.y*p1.z*p2.z - 2*p1.x*p2.x*p1.z*p2.z;

    double denominator = sqrt(p1.x*p1.x + p2.x*p2.x + p1.y*p1.y + p2.y*p2.y + p1.z*p1.z + p2.z*p2.z - 2*p1.x*p2.x - 2*p1.y*p2.y - 2*p1.z*p2.z);

    return std::make_pair(nearest_z, sqrt(part123)/denominator);
}

std::pair<double, double> nearestZ_revised (const myPoint3D &p1, const myPoint3D &p2) {
    double x1 = p1.x, y1 = p1.y, z1 = p1.z;
    double dx = p2.x - x1, dy = p2.y - y1, dz = p2.z - z1;
    double bPlusde            = 2*(x1*dx + y1*dy);
    double fourAMinusdSQUARED = 4*(dx*dx + dy*dy); 

    double z_m = z1 + (-2)*dz*bPlusde/fourAMinusdSQUARED;
    // double dist_m = std::sqrt(x1*x1 + y1*y1 - bPlusde*bPlusde/fourAMinusdSQUARED);
    double dist_mSQUARED = x1*x1 + y1*y1 - bPlusde*bPlusde/fourAMinusdSQUARED;

    return std::make_pair(z_m, dist_mSQUARED);
}

std::vector<std::string> splitString(const std::string &str, char delim) {
    std::vector<std::string> tokens;
    std::string token;
    std::istringstream tokenStream(str);

    while (std::getline(tokenStream, token, delim)) {
        tokens.push_back(token);
    }

    return tokens;
}

bool isInteger(const std::string& s) {
    std::istringstream iss(s);
    int x;
    return (iss >> x) && (iss.eof());
}

void current_PC_time () {
    auto now   = std::chrono::system_clock::now();
    auto now_c = std::chrono::system_clock::to_time_t(now);
    std::cout << "Run started at: " << std::put_time(std::localtime(&now_c), "%Y-%m-%d %H:%M:%S") << std::endl;
}

std::vector<int> readCsvToVector(const std::string& filename) {
    std::vector<int> result;
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Failed to open file: " << filename << std::endl;
        return result;
    }
    
    std::string line;
    if (std::getline(file, line)) {
        std::stringstream ss(line);
        std::string value;
        
        while (std::getline(ss, value, ',')) {
            result.push_back(std::stoi(value));
        }
    }
    
    file.close();
    return result;
}

std::string findRootFile (const std::string &directory, const std::string &prefix) {
    for (const auto &entry : std::filesystem::directory_iterator(directory)) {
        if (entry.is_regular_file()) {
            std::string fileName = entry.path().filename().string();
            if (fileName.find(prefix) == 0 && fileName.find(".root") != std::string::npos) {
                return entry.path().string();
            }
        }
    }
    return "";
}

std::vector<std::string> findRootFiles (const std::string &directory, const std::string &prefix) {
    std::vector<std::string> matchingFiles;
    for (const auto &entry : std::filesystem::directory_iterator(directory)) {
        if (entry.is_regular_file()) {
            std::string fileName = entry.path().filename().string();
            if (fileName.find(prefix) == 0 && fileName.find(".root") != std::string::npos) {
                matchingFiles.push_back(entry.path().string());
            }
        }
    }
    std::sort(matchingFiles.begin(), matchingFiles.end());
    return matchingFiles;
}

template <typename T>
void printRed(const T &content) {
    std::cout << COLOR_RED << content << COLOR_RESET << std::endl;
}

template <typename T>
void printWhite(const T &content) {
    std::cout << COLOR_WHITE << content << COLOR_RESET << std::endl;
}

template <typename T>
void printBlue(const T &content) {
    std::cout << COLOR_BLUE << content << COLOR_RESET << std::endl;
}

void printSeparation () {
    printBlue("=======================================================================");
}

#endif