#ifndef GLOBAL_DEFINITIONS_H
#define GLOBAL_DEFINITIONS_H
#include <vector>
#include <cmath>
#include <string>
#include <memory>
#include <TH1.h>
#include <TH2F.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TPolyMarker.h>
#include <TAxis.h>
/*********************************************************************
 *                      GLOBAL VARIABLES;
 * ******************************************************************/
double zmin = -40.6;
double zmax = 3.2;        // From sPHENIX paper, the stave's length is aroung 27.12 cm;
double scanstep = 0.1;  // unit: cm; 
int bins = (zmax - zmin)/scanstep + 1;
double DCA_cut = 0.2;   // unit: cm;
double MBD_lower = 0., MBD_upper = 10.;

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

struct myTrackletMember {
    Double_t x, y, z, r;
    Double_t eta, phi;
    Int_t layer;
};

struct myTrackletMemberExtended : myTrackletMember {
    Int_t trackID;
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

std::vector<std::string> splitString(const std::string &str, char delim) {
    std::vector<std::string> tokens;
    std::string token;
    std::istringstream tokenStream(str);

    while (std::getline(tokenStream, token, delim)) {
        tokens.push_back(token);
    }

    return tokens;
}



#endif