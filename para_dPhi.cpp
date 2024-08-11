#include "Fit/BinData.h"
#include "HFitInterface.h"
#include "Fit/Fitter.h"
#include "Math/WrappedMultiTF1.h"
#include "Math/WrappedParamFunction.h"
#include "TROOT.h"
#include "TVirtualFitter.h"
#include <TStopwatch.h>
#include "headers/histogramPlotting.h"
#include "TCanvas.h"
#include "TApplication.h"
#include <thread>
#include <mutex>
#include <chrono>

TH1 *h[100];
TH1D *h_ZonOne[20];
TH1D *h_CenonOne[14];
std::mutex m_mutex;
int batch_number = 50;   // 28147*8 = 225176
int N = 1000;
double range_min = -M_PI;
double range_max = M_PI;
double bin_width = (range_max - range_min) / N;
TH1D *h_dPhi_nomix = new TH1D("", "", N, range_min, range_max);
TH1D *h_Background_dPhi = new TH1D("", "", N, range_min, range_max);

void dPhiAccumulator (int id, std::vector<int> index, std::vector<double> MBD_cen, TBranch *branch11, TBranch* branch16, std::vector<float>* ClusPhi, std::vector<int>* ClusLayer) {
    double phi, dPhi, phi_0, phi_1;
    std::vector<double> Phi0, Phi1;
    
    int local_index;
    std::unique_lock<std::mutex> lock(m_mutex);
    for (int i = 0; i < batch_number; i++) {
        local_index = i + batch_number*id;
        branch16->GetEntry(index[local_index]);   branch11->GetEntry(index[local_index]);
        // std::cout << MBD_true_z[local_index] << std::endl;
        // std::cout << id << ", " << local_index << std::endl;
        // std::cout << id << ", " << local_index << ", " << index[local_index] << ", " << MBD_true_z[local_index] << ", " << MBD_cen[local_index] << ", " << ClusPhi->size() << ", " << ClusLayer->size() << std::endl;
        Phi0.reserve(ClusPhi->size());
        Phi1.reserve(ClusPhi->size());
        for (int j = 0; j < ClusPhi->size(); j++) {
            phi = ClusPhi->at(j);
            if (ClusLayer->at(j) == 3 || ClusLayer->at(j) == 4) {
                Phi0.push_back(phi);
            }
            else {
                Phi1.push_back(phi);
            }
        }
        for (int k = 0; k < Phi0.size(); k++) {
            for (int l = 0; l < Phi1.size(); l++) {
                dPhi = Phi0[k] - Phi1[l];
                if (dPhi > M_PI)    dPhi -= 2 * M_PI;
                if (dPhi < -M_PI)   dPhi += 2 * M_PI;
                h_dPhi_nomix->Fill(dPhi);
            }
        }
        Phi0.clear();   Phi1.clear();
    }
    lock.unlock();
}

void dPhi_in_bins_of_Centrality_ver1 (int id, int target, std::vector<int> index, std::vector<double> MBD_true_z, std::vector<double> MBD_cen, TBranch *branch11, TBranch* branch16, std::vector<float>* ClusPhi, std::vector<int>* ClusLayer) {
    h_CenonOne[id] = new TH1D(Form("dPhi of %f to %f", static_cast<double>(id)*0.05, static_cast<double>(id + 1)*0.05), Form("dPhi of %f to %f;dPhi;# of counts", static_cast<double>(id)*0.05, static_cast<double>(id + 1)*0.05), N, range_min, range_max);
    double phi, dPhi, phi_0, phi_1;
    std::vector<double> Phi0, Phi1;
    std::unique_lock<std::mutex> lock(m_mutex);
    double cen_lower_range  = static_cast<double>(id)*0.05;
    double cen_higher_range = static_cast<double>(id+1)*0.05;
    for (int i = 0; i < target; i++) {
        if (MBD_cen[i] >= cen_lower_range && MBD_cen[i] <= cen_higher_range) {
            branch16->GetEntry(index[i]);   // ClusPhi;
            branch11->GetEntry(index[i]);   // ClusLayer;
            Phi0.reserve(ClusPhi->size());  Phi1.reserve(ClusPhi->size());
            for (int j = 0; j < ClusPhi->size(); j++) {
                double phi = ClusPhi->at(j);
                if (ClusLayer->at(j) == 3 || ClusLayer->at(j) == 4) {
                    Phi0.push_back(phi);
                }
                else {
                    Phi1.push_back(phi);
                }
            }

            for (int k = 0; k < Phi0.size(); k++) {
                for (int l = 0; l < Phi1.size(); l++) {
                    dPhi = Phi0[k] - Phi1[l];
                    if (dPhi > M_PI)    dPhi = dPhi - 2*M_PI;
                    if (dPhi < -M_PI)   dPhi = dPhi + 2*M_PI;
                    h_CenonOne[id] -> Fill(dPhi);
                }
            }

            Phi0.clear();   Phi1.clear();
        }
    }
    lock.unlock();
}

void dPhi_in_bins_of_Z_vtx_ver1 (int id, int target, std::vector<int> index, std::vector<double> MBD_true_z, std::vector<double> MBD_cen, TBranch *branch11, TBranch* branch16, std::vector<float>* ClusPhi, std::vector<int>* ClusLayer) {
    h_ZonOne[id] = new TH1D(Form("dPhi of %f to %f", static_cast<double>(id)*0.05, static_cast<double>(id + 1)*0.05), Form("dPhi of %f to %f;dPhi;# of counts", static_cast<double>(id)*0.05, static_cast<double>(id + 1)*0.05), N, range_min, range_max);
    double phi, dPhi, phi_0, phi_1;
    std::vector<double> Phi0, Phi1;
    std::unique_lock<std::mutex> lock(m_mutex);
    double z_lower_range  = -6 - static_cast<double>(id);
    double z_higher_range = -5 - static_cast<double>(id);
    h_ZonOne[id] = new TH1D(Form("dPhi of %1.0f to %1.0f", z_lower_range, z_higher_range), Form("dPhi of %1.0f to %1.0f;dPhi;# of counts", z_lower_range, z_higher_range), N, range_min, range_max);
    for (int i = 0; i < target; i++) {
        if (MBD_true_z[i] >= z_lower_range && MBD_true_z[i] <= z_higher_range) {
            branch16->GetEntry(index[i]);   // ClusPhi;
            branch11->GetEntry(index[i]);   // ClusLayer;
            Phi0.reserve(ClusPhi->size());  Phi1.reserve(ClusPhi->size());
            for (int j = 0; j < ClusPhi->size(); j++) {
                double phi = ClusPhi->at(j);
                if (ClusLayer->at(j) == 3 || ClusLayer->at(j) == 4) {
                    Phi0.push_back(phi);
                }
                else {
                    Phi1.push_back(phi);
                }
            }

            for (int k = 0; k < Phi0.size(); k++) {
                for (int l = 0; l < Phi1.size(); l++) {
                    dPhi = Phi0[k] - Phi1[l];
                    if (dPhi > M_PI)    dPhi = dPhi - 2*M_PI;
                    if (dPhi < -M_PI)   dPhi = dPhi + 2*M_PI;
                    h_ZonOne[id] -> Fill(dPhi);
                }
            }

            Phi0.clear();   Phi1.clear();
        }
    }
    lock.unlock();
}

void dPhi_mixing (int id, int batch_size, const std::vector<std::vector<double>>& event_Phi0, const std::vector<std::vector<double>>& event_Phi1) {
    int local_index;
    double dPhi;
    std::unique_lock<std::mutex> lock(m_mutex);
    for (int i = 0; i < batch_size; i++) {
        local_index = i + batch_size*id;
        for (int j = local_index; j < event_Phi1.size(); j++) {
            // Process pairs from event_Phi0[i] and event_Phi1[j]
            for (double phi0 : event_Phi0[local_index]) {
                for (double phi1 : event_Phi1[j]) {
                    dPhi = phi0 - phi1;
                    if (dPhi > M_PI) dPhi -= 2 * M_PI;
                    if (dPhi < -M_PI) dPhi += 2 * M_PI;
                    h_Background_dPhi->Fill(dPhi);
                }
            }
            // Process pairs from event_Phi0[j] and event_Phi1[i] to ensure symmetry
            for (double phi0 : event_Phi0[j]) {
                for (double phi1 : event_Phi1[local_index]) {
                    dPhi = phi0 - phi1;
                    if (dPhi > M_PI) dPhi -= 2 * M_PI;
                    if (dPhi < -M_PI) dPhi += 2 * M_PI;
                    h_Background_dPhi->Fill(dPhi);
                }
            }
        }
    }
    lock.unlock();
}

void testFitter_ver2 (int id) {
    auto func = (TF1*)gROOT->GetFunction("gaus");
    std::unique_lock<std::mutex> lock(m_mutex); // lock the mutex;
    for (int i = 0; i < 5; i++) {
        func -> SetParameters(1000, 0, 1);
        double zeros[] = {0, 0, 0};
        func -> SetParErrors(zeros);

        h[i + id*5] -> Fit(func, (i==0?"N":"NQ"), "", -3, 3);
        auto x = func -> GetParameter(1);
        auto y = func -> GetParameter(2);
        printf("Now: %d, %d, #%d of histograms, %.18le, %.18le\n", id, i, i + id*5, x, y);
    }
    lock.unlock();
}

int main(int argc, char* argv[]) {
    // Start the stopwatch:
    TStopwatch timer;   timer.Start();
    // Get and print the current PC time
    auto now = std::chrono::system_clock::now();
    auto now_c = std::chrono::system_clock::to_time_t(now);
    std::cout << "Run started at: " << std::put_time(std::localtime(&now_c), "%Y-%m-%d %H:%M:%S") << std::endl;
    TApplication theApp("App", &argc, argv);
    std::vector<std::string> method;
    // <METHOD> <TARGET> <CEN_LOW> <CEN_HIGH> <Z_LOW> <Z_HIGH>
    std::string func = "NIL";
    int target       = 100;
    double cen_low   = 0.0;
    double cen_high  = 0.7;
    double z_low     = 25.;
    double z_high    = 5.; 

    if (argc > 1) {
        if (argc > 2) target   = std::stoi(argv[2]);
        if (argc > 3) cen_low  = std::stod(argv[3]);
        if (argc > 4) cen_high = std::stod(argv[4]);
        if (argc > 5) z_low    = std::stod(argv[5]);
        if (argc > 6) z_high   = std::stod(argv[6]);
        for (int i = 1; i < argc; i++) {
            std::string substring = argv[i];
            method.push_back(substring);
        }
    }
    // for (const auto& option : method) {
    //     std::cout << "Option: " << option << std::endl;
    // }

    TFile *file = TFile::Open("../../External/Data_CombinedNtuple_Run20869_HotDead_BCO_ADC_Survey.root");
    if (!file || file->IsZombie()) {
        std::cerr << "Error opening file!" << std::endl;
        exit(1);
    }
    TTree *EventTree = (TTree*)file->Get("EventTree");
    if (!EventTree) {
        std::cerr << "Error getting TTree!" << std::endl;
        exit(1);
    }

    std::string filePath = "../zFindingResults/new_foundZ_DCAfit_0_16_-001_001_step2e-1_z_25_5.txt";
    std::ifstream myfile(filePath);
    if (!myfile.is_open()){
		std::cout << "Unable to open linelabel" << std::endl;
		system("read -n 1 -s -p \"Press any key to continue...\" echo");
		exit(1);
 	}

    std::string line, value;
    std::vector<int> index, event, NHits;
    std::vector<double> foundZ, MBD_true_z, MBD_cen;

    getline(myfile, line);
    while (getline(myfile, line)) {
        std::stringstream data(line);
        getline(data, value, ',');  int i     = std::stoi(value);
        getline(data, value, ',');  int e     = std::stoi(value);
        getline(data, value, ',');  int N     = std::stoi(value);
        getline(data, value, ',');  double f  = std::stod(value);
        getline(data, value, ',');  double Mz = std::stod(value);
        getline(data, value, ',');  double Mc = std::stod(value);

        if (method[0] == "mix") {
            if (Mz >= -z_low && Mz <= -z_high && Mc >= cen_low && Mc <= cen_high) {
                index.push_back(i);     event.push_back(e);         NHits.push_back(N);
                foundZ.push_back(f);    MBD_true_z.push_back(Mz);   MBD_cen.push_back(Mc);
            }
        }
        else {
            if (Mz >= -25. && Mz <= -5. && Mc <= 0.7) {
                index.push_back(i);     event.push_back(e);         NHits.push_back(N);
                foundZ.push_back(f);    MBD_true_z.push_back(Mz);   MBD_cen.push_back(Mc);
            }
        }
        
    }
    std::cout << "Selected size: " << event.size() << std::endl;
    Long64_t nEntries = EventTree -> GetEntries();
    int event25, NClus;
    float MBD_centrality, MBD_z_vtx;
    std::vector<float> *ClusX     = nullptr;
    std::vector<float> *ClusY     = nullptr;
    std::vector<float> *ClusZ     = nullptr;
    std::vector<int>   *ClusLayer = nullptr;
    std::vector<float> *ClusR     = nullptr;    // FYI only;
    std::vector<float> *ClusPhi   = nullptr;    // FYI only;
    std::vector<float> *ClusEta   = nullptr;    // FYI only;

    TBranch *branch25 = (TBranch*)EventTree->GetListOfBranches()->At(25);    // event25;
    TBranch *branch10 = (TBranch*)EventTree->GetListOfBranches()->At(10);    // NClus;
    TBranch *branch11 = (TBranch*)EventTree->GetListOfBranches()->At(11);    // ClusLayer;
    TBranch *branch12 = (TBranch*)EventTree->GetListOfBranches()->At(12);    // ClusX;
    TBranch *branch13 = (TBranch*)EventTree->GetListOfBranches()->At(13);    // ClusY;
    TBranch *branch14 = (TBranch*)EventTree->GetListOfBranches()->At(14);    // ClusZ;
    TBranch *branch30 = (TBranch*)EventTree->GetListOfBranches()->At(30);    // MBD_centrality;
    TBranch *branch31 = (TBranch*)EventTree->GetListOfBranches()->At(31);    // MBD_z_vtx;
    TBranch *branch15 = (TBranch*)EventTree->GetListOfBranches()->At(15);    // ClusR;
    TBranch *branch16 = (TBranch*)EventTree->GetListOfBranches()->At(16);    // ClusPhi;
    TBranch *branch17 = (TBranch*)EventTree->GetListOfBranches()->At(17);    // ClusEta;
    
    branch25->SetAddress(&event25);
    branch10->SetAddress(&NClus);
    branch11->SetAddress(&ClusLayer);
    branch12->SetAddress(&ClusX);
    branch13->SetAddress(&ClusY);
    branch14->SetAddress(&ClusZ);
    branch15->SetAddress(&ClusR);
    branch16->SetAddress(&ClusPhi);
    branch17->SetAddress(&ClusEta);
    branch30->SetAddress(&MBD_centrality);
    branch31->SetAddress(&MBD_z_vtx);

    target = target > event.size() ? event.size() : target;

    std::cout << target << "," << cen_low << "," << cen_high << "," << z_low << "," << z_high << std::endl;
    
    ROOT::EnableThreadSafety();
    // Use 'sysctl -n hw.logicalcpu' to determine max number of threads:
    if (method[0] == "nomix") {
        std::thread thsafe[8];
        std::cout<<"multi-thready safe:"<<std::endl;
        for(int i = 0; i < 8; ++i)     thsafe[i]= std::thread(dPhiAccumulator,i,std::cref(index),std::cref(MBD_cen),branch11,branch16,ClusPhi,ClusLayer);
        for(int i = 0; i < 8; ++i)     thsafe[i].join();
        angularPlot1D(h_dPhi_nomix, method, "dPhi of unmixed");

        TFile* file = TFile::Open("dPhi.root", "UPDATE");
        if (!file || file -> IsZombie()) {
            std::cerr << "You messed up!" << std::endl;
            exit(1);
        }
        h_dPhi_nomix -> Write("nomix", TObject::kOverwrite);
    }
    else if (method[0] == "perCen") {
        std::thread thsafe[14];
        std::cout << "safe dPhi of different centralities" << std::endl;
        for (int i = 0; i < 14; i++)
            thsafe[i] = std::thread(dPhi_in_bins_of_Centrality_ver1,i,target,std::cref(index),std::cref(MBD_true_z),std::cref(MBD_cen),branch11,branch16,ClusPhi,ClusLayer);

        for (int i = 0; i < 14; i++)
            thsafe[i].join();

        std::vector<TH1D*> h(h_CenonOne, h_CenonOne + 14);
        ArrayPlot1D_Rescale(h, method, "dPhi_per_centralities_rescale");
    }
    else if (method[0] == "perZ") {
        std::thread thsafe[20];
        std::cout << "safe dPhi of different Z vertices" << std::endl;
        for (int i = 0; i < 20; i++)
            thsafe[i] = std::thread(dPhi_in_bins_of_Z_vtx_ver1,i,target,std::cref(index),std::cref(MBD_true_z),std::cref(MBD_cen),branch11,branch16,ClusPhi,ClusLayer);

        for (int i = 0; i < 20; i++)
            thsafe[i].join();

        std::vector<TH1D*> h(h_ZonOne, h_ZonOne + 20);
        ArrayPlot1D_Rescale_ver2(h, method, "dPhi_per_Z_vtx_rescale");
    }
    else if (method[0] == "mix") {
        double phi, dPhi, phi_0, phi_1; int clus_layer;
        std::vector<double> Phi0, Phi1;
        std::vector<std::vector <double>> event_Phi0, event_Phi1;
        // Signal dPhi is unmixed events' dPhi:
        TH1D *h_Signal_dPhi = new TH1D("dPhi of unmixed", Form("dPhi of unmixed %d events;dPhi value;# of counts", target), N, range_min, range_max);
        // Fill Signals:
        for (int i = 0; i < target; i++) {
            branch16->GetEntry(index[i]);   branch11->GetEntry(index[i]);
            event_Phi0.push_back(std::vector <double>());   event_Phi1.push_back(std::vector <double>());
            for (int j = 0; j < ClusPhi->size(); j++) {
                phi        = ClusPhi->at(j);
                clus_layer = ClusLayer->at(j);
                if (clus_layer == 3 || clus_layer == 4) {
                    Phi0.push_back(phi);
                    event_Phi0[i].push_back(phi);
                }
                else {
                    Phi1.push_back(phi);
                    event_Phi1[i].push_back(phi);
                }
            }
            for (int k = 0; k < Phi0.size(); k++) {
                for (int l = 0; l < Phi1.size(); l++) {
                    dPhi = Phi0[k] - Phi1[l];
                    if (dPhi > M_PI)    dPhi = dPhi - 2*M_PI;
                    if (dPhi < -M_PI)   dPhi = dPhi + 2*M_PI;
                    h_Signal_dPhi -> Fill(dPhi);
                }
            }
            Phi0.clear();   Phi1.clear();
        }

        std::thread thsafe[8];
        std::cout << "Mixing events: \n" << std::endl;
        for (int i = 0; i < 8; i++)
            thsafe[i] = std::thread(dPhi_mixing,i,target/8,event_Phi0,event_Phi1);
        for (int i = 0; i < 8; i++)
            thsafe[i].join();


        TCanvas *c1 = new TCanvas("c1", "dPhi Histogram", 1920, 1056);
        double phi_range_low = -2.4, phi_range_high = -1.8;
    int bin_range_low = h_Background_dPhi->FindBin(phi_range_low), bin_range_high = h_Background_dPhi->FindBin(phi_range_high);
    double max_unmixed = -1, max_mixed = -1, current_binContent;
    for (int bin = bin_range_low; bin <= bin_range_high; bin++) {
        current_binContent = h_Signal_dPhi->GetBinContent(bin);
        if (max_unmixed < current_binContent)  max_unmixed = current_binContent;
        current_binContent = h_Background_dPhi->GetBinContent(bin);
        if (max_mixed < current_binContent)  max_mixed = current_binContent;
    }
    h_Signal_dPhi -> Add(h_Signal_dPhi, max_mixed/max_unmixed - 1);
    h_Signal_dPhi -> Draw("SAME");

    max_mixed   = h_Background_dPhi->GetBinContent(h_Background_dPhi->GetMaximumBin());
    max_unmixed = h_Signal_dPhi->GetBinContent(h_Signal_dPhi->GetMaximumBin());
    std::cout << h_Signal_dPhi->GetBinContent(h_Signal_dPhi->GetMaximumBin()) << ", " << h_Background_dPhi->GetBinContent(h_Background_dPhi->GetMaximumBin())  << std::endl;
    if (max_unmixed > max_mixed) {
        h_Signal_dPhi -> GetYaxis() -> SetRangeUser(1e7, max_unmixed*1.2);
        h_Background_dPhi -> GetYaxis() -> SetRangeUser(1e7, max_unmixed*1.2);
    }   
    else {
        h_Signal_dPhi -> GetYaxis() -> SetRangeUser(1e7, max_mixed*1.2);
        h_Background_dPhi -> GetYaxis() -> SetRangeUser(1e7, max_mixed*1.2);
    }                      
    h_Background_dPhi -> Draw("SAME");
    h_Background_dPhi -> SetLineColor(2);

    double pi = TMath::Pi();
    int bin_min  = 1;  // The first bin
    int bin_max  = h_Signal_dPhi->GetNbinsX();  // The last bin
    // Calculate bin positions for each label
    int bin_pi   = bin_max;
    int bin_0    = bin_min + (bin_max - bin_min)/2;
    int binPi_2  = bin_min + 3*(bin_max - bin_min)/4;
    int bin_pi_2 = bin_min + (bin_max - bin_min)/4;
    // Set the labels at the calculated positions
    h_Background_dPhi->GetXaxis()->SetBinLabel(bin_0, "0");
    h_Background_dPhi->GetXaxis()->SetBinLabel(bin_pi_2, "#frac{-#pi}{2}");
    h_Background_dPhi->GetXaxis()->SetBinLabel(binPi_2, "#frac{#pi}{2}");
    h_Background_dPhi->GetXaxis()->SetBinLabel(bin_pi, "#pi");
    h_Background_dPhi->GetXaxis()->SetBinLabel(bin_min, "-#pi");
    // Ensure the custom labels are displayed by setting the number of divisions
    h_Background_dPhi->GetXaxis()->SetNdivisions(9, 0, 0, kFALSE);
    h_Background_dPhi->GetXaxis()->SetLabelSize(0.04);
    // Update histogram to refresh the axis
    // h[0]->Draw("HIST");
    h_Background_dPhi->GetXaxis()->LabelsOption("h"); // Draw the labels vertically

    h_Signal_dPhi->GetXaxis()->SetBinLabel(bin_0, "0");
    h_Signal_dPhi->GetXaxis()->SetBinLabel(bin_pi_2, "#frac{-#pi}{2}");
    h_Signal_dPhi->GetXaxis()->SetBinLabel(binPi_2, "#frac{#pi}{2}");
    h_Signal_dPhi->GetXaxis()->SetBinLabel(bin_pi, "#pi");
    h_Signal_dPhi->GetXaxis()->SetBinLabel(bin_min, "-#pi");
    // Ensure the custom labels are displayed by setting the number of divisions
    h_Signal_dPhi->GetXaxis()->SetNdivisions(9, 0, 0, kFALSE);
    h_Signal_dPhi->GetXaxis()->SetLabelSize(0.04);
    // Update histogram to refresh the axis
    // h[0]->Draw("HIST");
    h_Signal_dPhi->GetXaxis()->LabelsOption("h"); // Draw the labels vertically

    h_Signal_dPhi -> GetXaxis() -> CenterTitle(true);
    h_Signal_dPhi -> GetYaxis() -> CenterTitle(true);
    h_Background_dPhi -> GetXaxis() -> CenterTitle(true);
    h_Background_dPhi -> GetYaxis() -> CenterTitle(true);

    // TLegend *lg = new TLegend(0.12, 0.8, 0.33, 0.9);
    // lg -> AddEntry(h_Background_dPhi, Form("z_vtx between %2.2f and %2.2f", z_lower_range, z_upper_range), "l");
    // gStyle -> SetLegendTextSize(.023);
    // lg->Draw("same");
    
        c1 -> SaveAs("../../External/zFindingPlots/dPhi_mixed.png");
        backgroundCancelling_dPhi(h_Background_dPhi, h_Signal_dPhi, method, target);
    }
    // TCanvas *c1 = new TCanvas("c1", "dPhi Histogram", 1920, 1056);
    // h_dPhi_nomix -> Draw();
    // c1->Update();    c1->Modified();

    // Stop the stopwatch and print the runtime:
    timer.Stop();   timer.Print();

    theApp.Run();
    return 0;
}
