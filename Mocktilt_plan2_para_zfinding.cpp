#include <iostream>
#include <limits>
#include <fstream>
#include <string>
#include <vector>
#include <sstream>
#include <numeric>
#include <math.h>
#include <algorithm>
#include <omp.h>
using namespace std;

/*********************************************************************
 *                      GLOBAL VARIABLES;
 * ******************************************************************/
float zmin = -21.94;
float zmax = 20.22;        // From sPHENIX paper, the stave's length is aroung 27.12 cm;
float scanstep = 0.0025;  // unit: cm; 
int bins = (zmax - zmin)/scanstep + 1;
float threshold_eta_xl = 0.0003; float threshold_phi_xl = 0.0003;
float threshold_eta_l = 0.0008; float threshold_phi_l = 0.0008;
float threshold_eta_s = 0.008; float threshold_phi_s = 0.008;
float ls_cut = 630;

/*********************************************************************
 *              SCAN WITH ALL THE GUESS OF Z COORDINATES
 * ******************************************************************/
float zScan (const vector<float> &gx, const vector<float> &gy, const vector<float> &gz, const float &etaThreshold, const float &phiThreshold) {
    vector<float> cnt(bins);

    vector<float> phi_values, relative_r_values;
    for (int j = 0; j < gx.size(); j++) {
        float relative_x = gx[j];
        float relative_y = gy[j];
        float phi = atan2(relative_y, relative_x);
        float relative_r = sqrt(relative_x * relative_x + relative_y * relative_y);
        
        phi_values.push_back(phi);
        relative_r_values.push_back(relative_r);
    }
    
    #pragma omp parallel for
    for (int s = 0; s < bins; s++) {
        int counter = 0;
        float x1, y1, z1, rg1, rv1, theta1, eta1, phi1;
        float x2, y2, z2, rg2, rv2, theta2, eta2, phi2;
        int layer1, layer2;
        vector<float> eta_0, eta_1, eta_2;
        vector<float> phi_0, phi_1, phi_2;
        #pragma omp parallel for
        for (int i = 0; i < gx.size(); i++) {
            // float x = gx[i];
            // float y = gy[i];
            float z = gz[i] - (zmin + s*scanstep);
            // float phi = atan2(y, x);
            // float rv = sqrt(x*x + y*y); float rg = sqrt(gx[i]*gx[i] + gy[i]*gy[i]);
            float theta = atan2(relative_r_values[i], z);
            float eta;
            if (z >= 0) {
                eta = -log(tan(theta/2));
            }
            else {
                eta = log(tan((M_PI - theta)/2));
            }
            if (relative_r_values[i] >= 2.37 && relative_r_values[i] <= 2.80) {
			    eta_0.push_back(eta);
                phi_0.push_back(phi_values[i]);
		    }
		    else if (relative_r_values[i] >= 3.14 && relative_r_values[i] <= 3.59) {
			    eta_1.push_back(eta);
                phi_1.push_back(phi_values[i]);
		    }
		    else if (relative_r_values[i] >= 3.91 && relative_r_values[i] <= 4.3925) {
			    eta_2.push_back(eta);
                phi_2.push_back(phi_values[i]);
		    }
        }
        #pragma omp parallel for reduction(+:counter)
        for (int i = 0; i < eta_0.size(); i++) {
            for (int j = 0; j < eta_1.size(); j++) {
                float d_eta = abs(eta_0[i] - eta_1[j]);
                float d_phi = abs(phi_0[i] - phi_1[j]);
                bool within_deta_range = d_eta < etaThreshold;
                bool within_dphi_range = d_phi < phiThreshold;
                if (within_deta_range && within_dphi_range) {
                    counter++;
                }
            }
        }
        #pragma omp parallel for reduction(+:counter)
        for (int i = 0; i < eta_0.size(); i++) {
            for (int j = 0; j < eta_2.size(); j++) {
                float d_eta = abs(eta_0[i] - eta_2[j]);
                float d_phi = abs(phi_0[i] - phi_2[j]);
                bool within_deta_range = d_eta < etaThreshold;
                bool within_dphi_range = d_phi < phiThreshold;
                if (within_deta_range && within_dphi_range) {
                    counter++;
                }
            }
        }
        #pragma omp parallel for reduction(+:counter)
        for (int i = 0; i < eta_1.size(); i++) {
            for (int j = 0; j < eta_2.size(); j++) {
                float d_eta = abs(eta_1[i] - eta_2[j]);
                float d_phi = abs(phi_1[i] - phi_2[j]);
                bool within_deta_range = d_eta < etaThreshold;
                bool within_dphi_range = d_phi < phiThreshold;
                if (within_deta_range && within_dphi_range) {
                    counter++;
                }
            }
        }
        eta_0.clear(); std::vector<float>().swap(eta_0);
        eta_1.clear(); std::vector<float>().swap(eta_1);
        eta_2.clear(); std::vector<float>().swap(eta_2);
        phi_0.clear(); std::vector<float>().swap(phi_0);
        phi_1.clear(); std::vector<float>().swap(phi_1);
        phi_2.clear(); std::vector<float>().swap(phi_2);
        cnt[s] = counter;
    }
    
    auto most = max_element(cnt.begin(), cnt.end());
    size_t max_i = most - cnt.begin();
    float ctz = ((zmin + max_i*scanstep)*cnt[max_i] + (zmin + (max_i - 1)*scanstep)*cnt[max_i - 1] + (zmin + (max_i + 1)*scanstep)*cnt[max_i + 1]) 
                / (cnt[max_i] + cnt[max_i - 1] + cnt[max_i + 1]);

    cnt.clear(); phi_values.clear(); relative_r_values.clear();
    std::vector<float>().swap(cnt);
    std::vector<float>().swap(phi_values);
    std::vector<float>().swap(relative_r_values);

    return ctz;
}

/*********************************************************************
 *                      MAIN PROGRAMME;
 * ******************************************************************/
int main (int argc, char *argv[]) {
    omp_set_num_threads(10);
    ifstream myfile("plan2_tilted_upper_everyThing.txt");
    // ifstream myfile("plan2_tilted_upper.txt");
	if (!myfile .is_open()) {
		cout << "Unable to open sorted file." << endl;
		system("read -n 1 -s -p \"Press any key to continue...\" echo");
		exit(1);
	}
    ofstream outFile("plan2_para_foundZ.txt", std::ios_base::app);
	if (!outFile .is_open()) {
		cout << "Unable to open modified file." << endl;
		system("read -n 1 -s -p \"Press any key to continue...\" echo");
		exit(1);
	}

    int event_num;
    int indicator = 990;    int cntr = 0;
    float global_x, global_y, global_z, vertex_x, vertex_y, vertex_z;
    vector<float> gx, gy, gz, vx, vy, vz;
    vector<int> evt;
    string line_in, fragment;
    int starter = atoi(argv[1]);
    while(getline(myfile, line_in)) {
        stringstream str(line_in);
        getline(str, fragment, ',');
        event_num = stoi(fragment);
        getline(str, fragment, ',');
        global_x = stod(fragment);
        getline(str, fragment, ',');
        global_y = stod(fragment);
        getline(str, fragment, ',');
        global_z = stod(fragment);
        
        if (indicator != event_num) {
            float re;
            if (cntr >= starter) {
                if (gx.size() > ls_cut) {
                    re = zScan(gx, gy, gz, threshold_eta_l, threshold_phi_l);
                }else {
                    re = zScan(gx, gy, gz, threshold_eta_s, threshold_phi_s);
                }
                outFile << indicator << "," << re << endl;
            }
            cntr++;
            indicator = event_num;
            gx.clear();
            gy.clear();
            gz.clear();
            evt.clear();
            std::vector<float>().swap(gx);
            std::vector<float>().swap(gy);
            std::vector<float>().swap(gz);
            std::vector<int>().swap(evt);
        }
        gx.push_back(global_x);
        gy.push_back(global_y);
        gz.push_back(global_z);
        evt.push_back(event_num);
    }
    return 0;
}