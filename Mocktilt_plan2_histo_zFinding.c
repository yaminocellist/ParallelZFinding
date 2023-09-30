#include "TH1.h"
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
//////////////////////////////////////////////////////////////////////

/*************************************
 * @brief Remove NaN values in vetor
 * 
 * @param vec 
 * @return ** void 
 ************************************/

void remove_nan(std::vector<float>& vec) {
    vec.erase(
        std::remove_if(
            vec.begin(),
            vec.end(),
            [](float value) { return std::isnan(value); }
        ),
        vec.end()
    );
}

/*********************************************************************
 *                  Another Scanning method for zFinding:
 * ******************************************************************/
float zScan2 (const int &evt, const vector<float> &gx, const vector<float> &gy, const vector<float> &gz, const float &etaThreshold, const float &phiThreshold) {
    // /*
    vector<float> phi_values, relative_r_values;
    vector<float> eta_0, eta_1, eta_2;
    vector<float> phi_0, phi_1, phi_2;
    TH1D *h = new TH1D("", "", bins, zmin - scanstep/2, zmax + scanstep/2);
    for (int j = 0; j < gx.size(); j++) {
        float relative_x = gx[j];
        float relative_y = gy[j];
        float phi = atan2(relative_y, relative_x);
        float global_r = sqrt(relative_x * relative_x + relative_y * relative_y);
        float relative_r = global_r;
        
        if (relative_r >= 2.37 && relative_r <= 2.80) {
		    eta_0.push_back(1.);            phi_0.push_back(1.);
            eta_1.push_back(std::nan(""));  phi_1.push_back(std::nan(""));
            eta_2.push_back(std::nan(""));  phi_2.push_back(std::nan(""));
		}
		else if (relative_r >= 3.14 && relative_r <= 3.59) {
		    eta_0.push_back(std::nan(""));  phi_0.push_back(std::nan(""));
            eta_1.push_back(1.);            phi_1.push_back(1.);
            eta_2.push_back(std::nan(""));  phi_2.push_back(std::nan(""));
		}
		else if (relative_r >= 3.91 && relative_r <= 4.3925) {
		    eta_0.push_back(std::nan(""));  phi_0.push_back(std::nan(""));
            eta_1.push_back(std::nan(""));  phi_1.push_back(std::nan(""));
            eta_2.push_back(1.);            phi_2.push_back(1.);
		}
        else {
            eta_0.push_back(std::nan(""));  phi_0.push_back(std::nan(""));
            eta_1.push_back(std::nan(""));  phi_1.push_back(std::nan(""));
            eta_2.push_back(std::nan(""));  phi_2.push_back(std::nan(""));
        }

        phi_values.push_back(phi);
        relative_r_values.push_back(relative_r);
    }
    // std::cout << eta_0.size() << "," << eta_1.size() << "," << eta_2.size() << "," << phi_0.size() << "," << phi_1.size() << "," << phi_2.size() << std::endl;

    for (int s = 0; s < bins; s++) {
        std::vector<float> eta_0Dummy = eta_0;
        std::vector<float> eta_1Dummy = eta_1;
        std::vector<float> eta_2Dummy = eta_2;
        std::vector<float> phi_0Dummy = phi_0;
        std::vector<float> phi_1Dummy = phi_1;
        std::vector<float> phi_2Dummy = phi_2;
        for (int i = 0; i < gx.size(); i++) {
            // double x = gx[i];
            // double y = gy[i];
            float z = gz[i] - (zmin + s*scanstep);
            // double phi = atan2(y, x);
            // double rv = sqrt(x*x + y*y); double rg = sqrt(gx[i]*gx[i] + gy[i]*gy[i]);
            float theta = atan2(relative_r_values[i], z);
            float eta;
            if (z >= 0) {
                eta = -log(tan(theta/2));
            }
            else {
                eta = log(tan((M_PI - theta)/2));
            }
            eta_0Dummy[i] *= eta;   eta_1Dummy[i] *= eta;   eta_2Dummy[i] *= eta;
            phi_0Dummy[i] *= phi_values[i]; phi_1Dummy[i] *= phi_values[i]; phi_2Dummy[i] *= phi_values[i];
        }

        remove_nan(eta_0Dummy); remove_nan(eta_1Dummy); remove_nan(eta_2Dummy);
        remove_nan(phi_0Dummy); remove_nan(phi_1Dummy); remove_nan(phi_2Dummy);

        // eta_0Dummy.erase(std::remove_if(eta_0Dummy.begin(), eta_0Dummy.end(), [](float value) { return std::isnan(value); } ), eta_0Dummy.end());
        // eta_1Dummy.erase(std::remove_if(eta_1Dummy.begin(), eta_1Dummy.end(), [](float value) { return std::isnan(value); } ), eta_1Dummy.end());
        // eta_2Dummy.erase(std::remove_if(eta_2Dummy.begin(), eta_2Dummy.end(), [](float value) { return std::isnan(value); } ), eta_2Dummy.end());
        // phi_0Dummy.erase(std::remove_if(phi_0Dummy.begin(), phi_0Dummy.end(), [](float value) { return std::isnan(value); } ), phi_0Dummy.end());
        // phi_1Dummy.erase(std::remove_if(phi_1Dummy.begin(), phi_1Dummy.end(), [](float value) { return std::isnan(value); } ), phi_1Dummy.end());
        // phi_2Dummy.erase(std::remove_if(phi_2Dummy.begin(), phi_2Dummy.end(), [](float value) { return std::isnan(value); } ), phi_2Dummy.end());

        // for(float val : eta_0Dummy) {
        //     std::cout << val << " ";
        // }
        // std::cout << std::endl;

        for (int i = 0; i < eta_0Dummy.size(); i++) {
            for (int j = 0; j < eta_1Dummy.size(); j++) {
                float d_eta = abs(eta_0Dummy[i] - eta_1Dummy[j]);
                float d_phi = abs(phi_0Dummy[i] - phi_1Dummy[j]);
                bool within_deta_range, within_dphi_range;
                within_deta_range = d_eta < etaThreshold;
                within_dphi_range = d_phi < phiThreshold;
                if (within_deta_range && within_dphi_range ) {
                    h -> Fill(zmin + s*scanstep);
                }
            }
        }
        for (int i = 0; i < eta_0Dummy.size(); i++) {
            for (int j = 0; j < eta_2Dummy.size(); j++) {
                float d_eta = abs(eta_0Dummy[i] - eta_2Dummy[j]);
                float d_phi = abs(phi_0Dummy[i] - phi_2Dummy[j]);
                bool within_deta_range, within_dphi_range;
                within_deta_range = d_eta < etaThreshold;
                within_dphi_range = d_phi < phiThreshold;
                if (within_deta_range && within_dphi_range ) {
                    h -> Fill(zmin + s*scanstep);
                }
            }
        }
        for (int i = 0; i < eta_1Dummy.size(); i++) {
            for (int j = 0; j < eta_2Dummy.size(); j++) {
                float d_eta = abs(eta_1Dummy[i] - eta_2Dummy[j]);
                float d_phi = abs(phi_1Dummy[i] - phi_2Dummy[j]);
                bool within_deta_range, within_dphi_range;
                within_deta_range = d_eta < etaThreshold;
                within_dphi_range = d_phi < phiThreshold;
                if (within_deta_range && within_dphi_range ) {
                    h -> Fill(zmin + s*scanstep);
                }
            }
        }
        eta_0Dummy.clear(); std::vector<float>().swap(eta_0Dummy);
        eta_1Dummy.clear(); std::vector<float>().swap(eta_1Dummy);
        eta_2Dummy.clear(); std::vector<float>().swap(eta_2Dummy);
        phi_0Dummy.clear(); std::vector<float>().swap(phi_0Dummy);
        phi_1Dummy.clear(); std::vector<float>().swap(phi_1Dummy);
        phi_2Dummy.clear(); std::vector<float>().swap(phi_2Dummy);
    }

    int binmax = h -> GetMaximumBin(); 
    float xlabel0 = h -> GetXaxis() -> GetBinCenter(binmax);
    int entry0 = h -> GetBinContent(binmax);
    float xlabel_1 = h -> GetXaxis() -> GetBinCenter(binmax - 1); double entry_1 = h -> GetBinContent(binmax - 1);
    float xlabel1 = h -> GetXaxis() -> GetBinCenter(binmax + 1); double entry1 = h -> GetBinContent(binmax + 1);
    float ctz = (xlabel0*entry0 + xlabel_1*entry_1 + xlabel1*entry1) / (entry0 + entry1 + entry_1);

    TCanvas *can2 = new TCanvas("c2","c2",0,50,1800,550);
    h -> Draw();
    h -> SetFillColor(kYellow - 7);
    h -> SetLineWidth(1);
    h -> SetFillStyle(1001);
    h -> GetXaxis() -> SetTitle("Z coordinate [cm]");
    h -> GetXaxis() -> SetTitleSize(.05);
    h -> GetXaxis() -> SetLabelSize(.04);
    h -> GetXaxis() -> CenterTitle(true);
    h -> GetXaxis() -> SetNdivisions(31, 5, 0);
    h -> GetXaxis() -> SetTitleOffset(.8);
    h -> GetYaxis() -> SetTitle("# of Counts");
    h -> GetYaxis() -> SetTitleSize(.05);
    h -> GetYaxis() -> SetLabelSize(.04);
    h -> GetYaxis() -> SetTitleOffset(.62);
    h -> GetYaxis() -> CenterTitle(true);
    h -> SetTitle("Displacement at the farthest end is 0.4 mm");
    h -> SetMinimum(0);
    int max_entry = h -> GetBinContent(h -> FindBin(ctz));
    TLine *l1 = new TLine(ctz, 0, ctz, max_entry);
    l1 -> SetLineColor(kRed);
    l1 -> SetLineStyle(2);
    l1 -> SetLineWidth(1);
    l1 -> Draw("same");
    TLegend *lg = new TLegend(0.12, 0.8, 0.72, 0.9);
    if (gx.size() > ls_cut) {
        lg -> AddEntry(h, Form("%2.4fcm, Scan_step=%1.0f#mum, threshold_eta |#Delta#eta|<%1.5f, |#Delta#phi|<%1.5f", ctz, scanstep*1e4, threshold_eta_l, threshold_phi_l), "f");
    }else {
        lg -> AddEntry(h, Form("%2.4fcm, Scan_step=%1.0f#mum, threshold_eta |#Delta#eta|<%1.5f, |#Delta#phi|<%1.5f", ctz, scanstep*1e4, threshold_eta_s, threshold_phi_s), "f");
    }
    
    gStyle -> SetLegendTextSize(.03);
    lg -> Draw("same");
    gPad -> SetGrid(1,1); gPad -> Update();
    can2 -> SaveAs(Form("histoZScan/foundz_%d_04.png", evt));
    delete can2;
    delete h;   h = 0;
    phi_values.clear(); relative_r_values.clear();
    std::vector<float>().swap(phi_values);
    std::vector<float>().swap(relative_r_values);
    eta_0.clear(); std::vector<float>().swap(eta_0);
    eta_1.clear(); std::vector<float>().swap(eta_1);
    eta_2.clear(); std::vector<float>().swap(eta_2);
    phi_0.clear(); std::vector<float>().swap(phi_0);
    phi_1.clear(); std::vector<float>().swap(phi_1);
    phi_2.clear(); std::vector<float>().swap(phi_2);
    return ctz;
}

/*********************************************************************
 *                  Scan for found_z vertex:
 * ******************************************************************/
float zScan (const int &evt, const vector<float> &gx, const vector<float> &gy, const vector<float> &gz, const float &etaThreshold, const float &phiThreshold) {
    vector<float> phi_values, global_r_values, relative_r_values;
    TH1D *h = new TH1D("", "", bins, zmin - scanstep/2, zmax + scanstep/2);
    for (int j = 0; j < gx.size(); j++) {
        float relative_x = gx[j];
        float relative_y = gy[j];
        float phi = atan2(relative_y, relative_x);
        float global_r = sqrt(relative_x * relative_x + relative_y * relative_y);
        float relative_r = global_r;
        
        phi_values.push_back(phi);
        global_r_values.push_back(global_r);
        relative_r_values.push_back(relative_r);
    }

    for (int s = 0; s < bins; s++) {
        vector<float> eta_0, eta_1, eta_2;
        vector<float> phi_0, phi_1, phi_2;
        
        for (int i = 0; i < gx.size(); i++) {
            // double x = gx[i];
            // double y = gy[i];
            float z = gz[i] - (zmin + s*scanstep);
            // double phi = atan2(y, x);
            // double rv = sqrt(x*x + y*y); double rg = sqrt(gx[i]*gx[i] + gy[i]*gy[i]);
            float theta = atan2(relative_r_values[i], z);
            float eta;
            if (z >= 0) {
                eta = -log(tan(theta/2));
            }
            else {
                eta = log(tan((M_PI - theta)/2));
            }
            if (global_r_values[i] >= 2.37 && global_r_values[i] <= 2.80) {
			    eta_0.push_back(eta);
                phi_0.push_back(phi_values[i]);
		    }
		    else if (global_r_values[i] >= 3.14 && global_r_values[i] <= 3.59) {
			    eta_1.push_back(eta);
                phi_1.push_back(phi_values[i]);
		    }
		    else if (global_r_values[i] >= 3.91 && global_r_values[i] <= 4.3925) {
			    eta_2.push_back(eta);
                phi_2.push_back(phi_values[i]);
		    }
        }

        for (int i = 0; i < eta_0.size(); i++) {
            for (int j = 0; j < eta_1.size(); j++) {
                float d_eta = abs(eta_0[i] - eta_1[j]);
                float d_phi = abs(phi_0[i] - phi_1[j]);
                bool within_deta_range, within_dphi_range;
                within_deta_range = d_eta < etaThreshold;
                within_dphi_range = d_phi < phiThreshold;
                if (within_deta_range && within_dphi_range) {
                    h -> Fill(zmin + s*scanstep);
                }
            }
        }
        for (int i = 0; i < eta_0.size(); i++) {
            for (int j = 0; j < eta_2.size(); j++) {
                float d_eta = abs(eta_0[i] - eta_2[j]);
                float d_phi = abs(phi_0[i] - phi_2[j]);
                bool within_deta_range, within_dphi_range;
                within_deta_range = d_eta < etaThreshold;
                within_dphi_range = d_phi < phiThreshold;
                if (within_deta_range && within_dphi_range) {
                    h -> Fill(zmin + s*scanstep);
                }
            }
        }
        for (int i = 0; i < eta_1.size(); i++) {
            for (int j = 0; j < eta_2.size(); j++) {
                float d_eta = abs(eta_1[i] - eta_2[j]);
                float d_phi = abs(phi_1[i] - phi_2[j]);
                bool within_deta_range, within_dphi_range;
                within_deta_range = d_eta < etaThreshold;
                within_dphi_range = d_phi < phiThreshold;
                if (within_deta_range && within_dphi_range) {
                    h -> Fill(zmin + s*scanstep);
                }
            }
        }
        eta_0.clear(); std::vector<float>().swap(eta_0);
        eta_1.clear(); std::vector<float>().swap(eta_1);
        eta_2.clear(); std::vector<float>().swap(eta_2);
        phi_0.clear(); std::vector<float>().swap(phi_0);
        phi_1.clear(); std::vector<float>().swap(phi_1);
        phi_2.clear(); std::vector<float>().swap(phi_2);
    }

    int binmax = h -> GetMaximumBin(); 
    float xlabel0 = h -> GetXaxis() -> GetBinCenter(binmax);
    int entry0 = h -> GetBinContent(binmax);
    float xlabel_1 = h -> GetXaxis() -> GetBinCenter(binmax - 1); double entry_1 = h -> GetBinContent(binmax - 1);
    float xlabel1 = h -> GetXaxis() -> GetBinCenter(binmax + 1); double entry1 = h -> GetBinContent(binmax + 1);
    float ctz = (xlabel0*entry0 + xlabel_1*entry_1 + xlabel1*entry1) / (entry0 + entry1 + entry_1);

    TCanvas *can2 = new TCanvas("c2","c2",0,50,1800,550);
    h -> Draw();
    h -> SetFillColor(kYellow - 7);
    h -> SetLineWidth(1);
    h -> SetFillStyle(1001);
    h -> GetXaxis() -> SetTitle("Z coordinate [cm]");
    h -> GetXaxis() -> SetTitleSize(.05);
    h -> GetXaxis() -> SetLabelSize(.04);
    h -> GetXaxis() -> CenterTitle(true);
    h -> GetXaxis() -> SetNdivisions(31, 5, 0);
    h -> GetXaxis() -> SetTitleOffset(.8);
    h -> GetYaxis() -> SetTitle("# of Counts");
    h -> GetYaxis() -> SetTitleSize(.05);
    h -> GetYaxis() -> SetLabelSize(.04);
    h -> GetYaxis() -> SetTitleOffset(.62);
    h -> GetYaxis() -> CenterTitle(true);
    h -> SetTitle("Displacement at the farthest end is 0.4 mm");
    h -> SetMinimum(0);
    int max_entry = h -> GetBinContent(h -> FindBin(ctz));
    TLine *l1 = new TLine(ctz, 0, ctz, max_entry);
    l1 -> SetLineColor(kRed);
    l1 -> SetLineStyle(2);
    l1 -> SetLineWidth(1);
    l1 -> Draw("same");
    TLegend *lg = new TLegend(0.12, 0.8, 0.72, 0.9);
    if (gx.size() > ls_cut) {
        lg -> AddEntry(h, Form("%2.4fcm, Scan_step=%1.0f#mum, threshold_eta |#Delta#eta|<%1.5f, |#Delta#phi|<%1.5f", ctz, scanstep*1e4, threshold_eta_l, threshold_phi_l), "f");
    }else {
        lg -> AddEntry(h, Form("%2.4fcm, Scan_step=%1.0f#mum, threshold_eta |#Delta#eta|<%1.5f, |#Delta#phi|<%1.5f", ctz, scanstep*1e4, threshold_eta_s, threshold_phi_s), "f");
    }
    
    gStyle -> SetLegendTextSize(.03);
    lg -> Draw("same");
    gPad -> SetGrid(1,1); gPad -> Update();
    can2 -> SaveAs(Form("histoZScan/foundz_%d_04.png", evt));
    delete can2;
    delete h;   h = 0;
    phi_values.clear(); global_r_values.clear(); relative_r_values.clear();
    std::vector<float>().swap(phi_values);
    std::vector<float>().swap(relative_r_values);
    std::vector<float>().swap(global_r_values);
    return ctz;
}

void Mocktilt_plan2_histo_zFinding (Int_t target) {
    ifstream myfile("plan2_tilted_upper_everyThing_0.4.txt");
	if (!myfile .is_open()) {
		cout << "Unable to open sorted file." << endl;
		system("read -n 1 -s -p \"Press any key to continue...\" echo");
		exit(1);
	}
    ofstream outFile("plan2_histo_foundZ_0.4.txt", std::ios_base::app);
	if (!outFile .is_open()) {
		cout << "Unable to open modified file." << endl;
		system("read -n 1 -s -p \"Press any key to continue...\" echo");
		exit(1);
	}

    int event_num;
    int indicator = 990;    int cntr = 0;
    float global_x, global_y, global_z;
    vector<float> gx, gy, gz;
    vector<int> evt;
    string line_in, fragment;
    while(getline(myfile, line_in)) {
        stringstream str(line_in);
        getline(str, fragment, ',');
        event_num = stoi(fragment);

        // getline(str, fragment, ','); getline(str, fragment, ','); getline(str, fragment, ',');
        // getline(str, fragment, ','); getline(str, fragment, ','); getline(str, fragment, ',');
        // getline(str, fragment, ','); getline(str, fragment, ','); getline(str, fragment, ',');
        // getline(str, fragment, ','); getline(str, fragment, ','); getline(str, fragment, ',');
        // getline(str, fragment, ','); getline(str, fragment, ','); getline(str, fragment, ',');
        // getline(str, fragment, ','); getline(str, fragment, ','); getline(str, fragment, ',');

        getline(str, fragment, ',');
        global_x = stof(fragment);
        getline(str, fragment, ',');
        global_y = stof(fragment);
        getline(str, fragment, ',');
        global_z = stof(fragment);
        
        // getline(str, fragment, ','); getline(str, fragment, ','); getline(str, fragment, ',');
        // getline(str, fragment, ','); getline(str, fragment, ','); getline(str, fragment, ',');
        
        // getline(str, fragment, ',');
        // vertex_x = stod(fragment);
        // getline(str, fragment, ',');
        // vertex_y = stod(fragment);
        // getline(str, fragment, ',');
        // vertex_z = stod(fragment);

        if (indicator != event_num) {
            if (cntr >= target) {
                float re;
                if (gx.size() > ls_cut) {
                    re = zScan2(evt[0], gx, gy, gz, threshold_eta_l, threshold_phi_l);
                }else {
                    re = zScan2(evt[0], gx, gy, gz, threshold_eta_s, threshold_phi_s);
                }                
                outFile << indicator << "," << re << "," << gx.size() << endl;
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
}