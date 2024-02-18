#ifndef GLOBAL_DEFINITIONS_H
#define GLOBAL_DEFINITIONS_H
#include <vector>
#include <cmath>
#include <string>
#include <memory>
#include <TH1.h>
#include <TH2F.h>
/*********************************************************************
 *                      GLOBAL VARIABLES;
 * ******************************************************************/
double zmin = -40.6;
double zmax = 3.2;        // From sPHENIX paper, the stave's length is aroung 27.12 cm;
double scanstep = 0.2;  // unit: cm; 
int bins = (zmax - zmin)/scanstep + 1;
double DCA_cut = 0.2;   // unit: cm;
double MBD_lower = 60., MBD_upper = 70;

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

std::vector<std::string> splitString(const std::string &str, char delim) {
    std::vector<std::string> tokens;
    std::string token;
    std::istringstream tokenStream(str);

    while (std::getline(tokenStream, token, delim)) {
        tokens.push_back(token);
    }

    return tokens;
}

struct myPoint3D {
    double x, y, z;
};

struct myTrackletMember {
    Double_t x, y, z, r;
    Double_t eta, phi;
    Int_t layer, trackID;
};

#endif