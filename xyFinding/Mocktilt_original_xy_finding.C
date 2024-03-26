#include <iostream>
#include "TH1.h"
#include <limits>
#include <fstream>
#include <string>
#include <vector>
#include <sstream>
#include <math.h>
#include "TMath.h"
#include "TArrayD.h"
#include "TMinuit.h"
#include <cmath>
#include <algorithm>

using namespace std;
//////////////////////////////////////////////////////////////////////
vector<float> upper_gx, upper_gy, upper_gz, lower_gx, lower_gy, lower_gz;
vector<float> upper_eta_value, upper_phi_value;
vector<float> lower_eta_value, lower_phi_value;
vector<int> upper_layer, lower_layer;
//////////////////////////////////////////////////////////////////////
float upper_angular_distance_cut, lower_angular_distance_cut;
float upper_found_x, upper_found_y, upper_found_z, lower_found_x, lower_found_y, lower_found_z;
//////////////////////////////////////////////////////////////////////


double gMedian(const TH1D * h1) {
   int n = h1->GetXaxis()->GetNbins();
   std::vector<double>  x(n);
   h1->GetXaxis()->GetCenter( &x[0] );
   const double * y = h1->GetArray();
   // exclude underflow/overflows from bin content array y
   return TMath::Median(n, &x[0], &y[1]);
}

std::tuple<vector<float>, vector<float>, vector<int>> find_eta_phi_layer_label(const vector<float> &gx, const vector<float> &gy, const vector<float> &gz, const float &found_z) {
    vector<float> eta_value, phi_value;
    vector<int> layer;
    for (int i = 0; i < gx.size(); i++) {
        float phi_i = atan2(gy[i], gx[i]);
        phi_value.push_back(phi_i);
        float r_i = sqrt(pow(gx[i], 2) + pow(gy[i], 2));  float z_i = gz[i] - found_z;
        if (r_i >= 2.37 && r_i <= 2.80) {
            layer.push_back(0);
	    }
	    else if (r_i >= 3.14 && r_i <= 3.59) {
            layer.push_back(1);
	    }
	    else if (r_i >= 3.91 && r_i <= 4.3925) {
            layer.push_back(2);
	    }
        else {
            layer.push_back(3);
        }
        float theta_i = atan2(r_i, z_i);
        float eta_i;
        if (z_i >= 0) {
            eta_i = -log(tan(theta_i/2));
        }
        else {
            eta_i = log(tan((M_PI - theta_i)/2));
        }
        eta_value.push_back(eta_i);
    }
    return std::make_tuple(eta_value, phi_value, layer);
}

/*******************************************
 * @brief find a cut for angular distance;
 * 
 * @return ** double 
 ******************************************/
float aggregate_angular_distance (const int & evt, const vector<float> &eta_value, const vector<float> &phi_value, const vector<int> &layer) {
    const int n1 = 10;  // number of bins in [0, 0.1]
    const int n2 = 100;   // number of bins in [0.1, 1.1], [1.1, 2.1], ...
    const int n3 = 8 * n2 + n1; // total number of bins
    float unequal_bins[n3 + 1]; // array of bin edges
    // Set bin edges
    for (int i = 0; i < n1; ++i) {
        unequal_bins[i] = 0.1 * static_cast<float>(i) / n1;  // divide [0, 0.1] into n1 bins
    }
    for (int i = 0; i < 8; ++i) {
        for (int j = 0; j < n2; ++j) {
            unequal_bins[n1 + i * n2 + j] = 0.1 + static_cast<float>(i) + static_cast<float>(j) / n2;  // divide [0.1 + i, 1.1 + i] into n2 bins
        }
    }
    unequal_bins[n3] = 8.1;  // set the last bin edge
    TH1F *h = new TH1F(Form("h%d", evt), Form("Histogram Title %d", evt), n3, unequal_bins); // 

    for (int i = 0; i < eta_value.size(); i++) {
        for (int j = i + 1; j < eta_value.size(); j++) {
            bool different_layers = layer[i] != layer[j] && layer[i] != 3 && layer[j] != 3;
            if (different_layers) {
                float angular_distance = sqrt(pow(eta_value[i] - eta_value[j], 2) + pow(phi_value[i] - phi_value[j], 2));
                h -> Fill(angular_distance);
            }        
        }
    }

    int bin1 = h->FindBin(0.0);  // find the bin number corresponding to 0.0
    int bin2 = h->FindBin(0.1);  // find the bin number corresponding to 0.1
    float minContent = h->GetBinContent(bin1);
    int minBin = bin1;  
    for (int i = bin1 + 1; i <= bin2; ++i) {
        float currentContent = h -> GetBinContent(i);
        if (currentContent < minContent) {
            minContent = currentContent;
            minBin = i;
        }
    }
    float rst = (h -> GetBinCenter(minBin)) - 0.1/(n1*2);
    h -> Reset("ICESM");
    delete h;   h = 0;
    return rst;
}

std::pair<float, float> find_x_y (const char &flag, const int & evt, const vector<float> &gx, const vector<float> &gy, const vector<float> &gz, const vector<float> &eta_value, const vector<float> &phi_value, const vector<int> &layer, const float &angular_distance_cut, const float &found_z) {
    TH2F *h_intersection = new TH2F("", "", 3600, -4.5, 4.5, 3600, -4.5, 4.5);
    float pixel_size = 9.0/3600;
    for (int i = 0; i < gx.size(); i++) {
        for (int j = i + 1; j < gx.size(); j++) {
            float dx = eta_value[i] - eta_value[j];
            float dy = phi_value[i] - phi_value[j];
            float delta_r = sqrt(dx*dx + dy*dy);
            bool a = delta_r <= angular_distance_cut && layer[i] != layer[j] && layer[i] != 3 && layer[j] != 3;
            if (a) {
                float t = (found_z - gz[i]) / (gz[j] - gz[i]);
                float x_at_found_z = gx[i] + t * (gx[j] - gx[i]);
                float y_at_found_z = gy[i] + t * (gy[j] - gy[i]);
                h_intersection -> Fill(x_at_found_z, y_at_found_z);
            }       
        }
    }

    TH1D *h_proj_x = h_intersection -> ProjectionX("h_proj_x");
    TH1D *h_proj_y = h_intersection -> ProjectionY("h_proj_y");
    float median_x = gMedian(h_proj_x);
    float median_y = gMedian(h_proj_y);
    // TCanvas *canvas = new TCanvas("canvas", "Classifier output", 800, 600);
    // h_intersection -> Draw("colz");
    // h_intersection -> GetYaxis() -> SetRangeUser(-0.1, 0.1);
    // h_intersection -> GetXaxis() -> SetRangeUser(-0.1, 0.1);
    // h_intersection -> GetXaxis() -> SetTitle("X (cm)");
    // h_intersection -> GetYaxis() -> SetTitle("Y (cm)");
    // h_intersection -> GetYaxis() -> SetTitleOffset(1);
    // h_intersection -> GetXaxis() -> SetTitleSize(0.05);
    // h_intersection -> GetYaxis() -> SetTitleSize(0.05); 
    // h_intersection -> SetTitle(Form("#splitline{Data from event #%d of MVTXRecoClusters_HijingMBwoPU0T_MockSim.root}{      Angular Separation cut at %1.3f}", evt, angular_distance_cut));
    // gStyle -> SetTitleSize(0.032, "t");
    // h_intersection -> GetXaxis() -> CenterTitle(true);
    // h_intersection -> GetYaxis() -> CenterTitle(true);
    // Color_t colorS = kTeal - 7;
    // TMarker *markerS = new TMarker(0, 0, 21); // 8 is the default marker style for TH2D
    // markerS->SetMarkerColor(colorS);
    // h_intersection -> SetMarkerColor(colorS);
    // TLegend *lg = new TLegend(0.12, 0.8, 0.7, 0.9);
    // lg -> AddEntry(markerS, Form("Found x and y vertex at: %1.5f, %1.5f (cm)", median_x, median_y), "P");
    // gStyle -> SetLegendTextSize(.03);
    // lg -> Draw("same");
    //canvas -> SaveAs(Form("Found_xy/field_off_%d_%c.png", evt, flag));
    //delete canvas;
    h_intersection -> Reset("ICESM");
    h_proj_y -> Reset("ICESM");
    h_proj_x -> Reset("ICESM");
    delete h_intersection;  delete h_proj_x;    delete h_proj_y;
    h_intersection = 0; h_proj_x = 0;   h_proj_y = 0;
    return std::make_pair(median_x, median_y);
}

void Mocktilt_original_xy_finding (Int_t target = 0) {
    ifstream found_information_file("original_foundz.txt");
    if (!found_information_file.is_open()) {
        cout << "Unable to open the text file!" << endl;
		system("read -n 1 -s -p \"Press any key to continue...\" echo");
		exit(1);
    }

    ofstream outfile2("original_lower_foundxy.txt", std::ios_base::app);
    if (!outfile2.is_open()) {
        cout << "Unable to open outfile" << endl;
		system("read -n 1 -s -p \"Press any key to continue...\" echo");
		exit(1);
    }
    ofstream outfile("original_all_foundxy.txt");
    if (!outfile.is_open()) {
        cout << "Unable to open outfile" << endl;
		system("read -n 1 -s -p \"Press any key to continue...\" echo");
		exit(1);
    }
    
    TFile *f = new TFile("/Users/yaminocellist/MIT_mentorship/2nd_semester/week7/MVTXRecoClusters_HijingMBwoPU0T_MockSim.root");
    TIter next(f->GetListOfKeys());
    TKey *key = (TKey*)next();
    TTree *T = (TTree*)key->ReadObj();

    int evt = -1;
    T -> SetBranchAddress("event", &evt);
    vector<float> *hitfromG4P_X = nullptr;
    T -> SetBranchAddress("HitfromG4P_X", &hitfromG4P_X);
    vector<float> *hitfromG4P_Y = nullptr;
    T -> SetBranchAddress("HitfromG4P_Y", &hitfromG4P_Y);
    vector<float> *hitfromG4P_Z = nullptr;
    T -> SetBranchAddress("HitfromG4P_Z", &hitfromG4P_Z);
    vector<float> *truthPV_x = nullptr;
    T -> SetBranchAddress("TruthPV_x", &truthPV_x);
    vector<float> *truthPV_y = nullptr;
    T -> SetBranchAddress("TruthPV_y", &truthPV_y);
    vector<float> *truthPV_z = nullptr;
    T -> SetBranchAddress("TruthPV_z", &truthPV_z);

    Long64_t nentries = T->GetEntries();

    for(Long64_t i = 0; i < nentries; i++)
    {
        T -> GetEntry(i);  // load data from file
        string line, frag;
        getline(found_information_file, line);
        stringstream str(line);
        getline(str, frag, ',');
        getline(str, frag, ',');
        upper_found_z = stof(frag);
        getline(str, frag, ',');
        lower_found_z = stof(frag);
        if (i >= target) {
            // if (size_of_event > 1100) {
            for (int i = 0; i < hitfromG4P_X -> size(); i++) {
                if (hitfromG4P_X -> at(i) >= 0) {
                    upper_gx.push_back(hitfromG4P_X -> at(i));
                    upper_gy.push_back(hitfromG4P_Y -> at(i));
                    upper_gz.push_back(hitfromG4P_Z -> at(i));
                }
                else {
                    lower_gx.push_back(hitfromG4P_X -> at(i));
                    lower_gy.push_back(hitfromG4P_Y -> at(i));
                    lower_gz.push_back(hitfromG4P_Z -> at(i));
                }
            }
            tie(upper_eta_value, upper_phi_value, upper_layer) = find_eta_phi_layer_label(upper_gx, upper_gy, upper_gz, upper_found_z);
            tie(lower_eta_value, lower_phi_value, lower_layer) = find_eta_phi_layer_label(lower_gx, lower_gy, lower_gz, lower_found_z);
            upper_angular_distance_cut = aggregate_angular_distance(evt, upper_eta_value, upper_phi_value, upper_layer);
            lower_angular_distance_cut = aggregate_angular_distance(evt, lower_eta_value, lower_phi_value, lower_layer);
            std::pair<float, float> upper_found_xy = find_x_y('u', evt, upper_gx, upper_gy, upper_gz, upper_eta_value, upper_phi_value, upper_layer, upper_angular_distance_cut, upper_found_z);
            std::pair<float, float> lower_found_xy = find_x_y('l', evt, lower_gx, lower_gy, lower_gz, lower_eta_value, lower_phi_value, lower_layer, lower_angular_distance_cut, lower_found_z);
            outfile << evt << "," << upper_found_xy.first << "," << upper_found_xy.second << "," << upper_found_z << "," << lower_found_xy.first << "," << lower_found_xy.second << "," << lower_found_z << "," << truthPV_x -> at(0) << "," << truthPV_y -> at(0) << "," << truthPV_z -> at(0) << "," << upper_gx.size() << "," << lower_gx.size() << "," << hitfromG4P_X -> size() << endl;
            outfile2 << evt << "," << lower_found_xy.first << "," << lower_found_xy.second << "," << lower_found_z << "," << truthPV_x -> at(0) << "," << truthPV_y -> at(0) << "," << truthPV_z -> at(0) << "," << lower_gx.size() << "," << hitfromG4P_X -> size() << endl;
            upper_gx.clear(); upper_gy.clear(); upper_gz.clear();
            lower_gx.clear(); lower_gy.clear(); lower_gz.clear();
            upper_eta_value.clear();  upper_phi_value.clear();  upper_layer.clear();
            lower_eta_value.clear();  lower_phi_value.clear();  lower_layer.clear();
        }
    }

}