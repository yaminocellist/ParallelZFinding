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
    TApplication theApp("App", &argc, argv);
    std::string method;
    if (argc == 2)   method = argv[1];
    std::vector<std::string> sub_options;
    if (!method.empty()) {
        sub_options = splitString(method, '_');
        if (sub_options.size() == 2) {
            sub_options.push_back("");
        }
        else if (sub_options.size() == 1) {
            sub_options.push_back("");
            sub_options.push_back("");
        }
    }
    else {
        sub_options.push_back("");
    }

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

        if (Mz >= -25. && Mz <= -5. && Mc <= 0.7) {
            index.push_back(i);     event.push_back(e);         NHits.push_back(N);
            foundZ.push_back(f);    MBD_true_z.push_back(Mz);   MBD_cen.push_back(Mc);
        }
    }

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

    // if (sub_options[0] == "foundz")      INTTZAnalysisLite(filePath, sub_options);
    // if (sub_options[0] == "dPhi")        INTTdPhiAnalysis(tree, filePath, sub_options);
    // if (sub_options[0] == "dPhiMix")     INTTMixingEvent(tree, filePath, sub_options);
    // if (sub_options[0] == "dEta")        INTTdEtaAnalysis(tree, filePath, sub_options);

    int target = std::stoi(sub_options[1]) > event.size() ? event.size() : std::stoi(sub_options[1]);
    ROOT::EnableThreadSafety();
    // Use 'sysctl -n hw.logicalcpu' to determine max number of threads:
    if (method == "nomix") {
        std::thread thsafe[8];
        std::cout<<"multi-thready safe:"<<std::endl;
        for(int i = 0; i < 8; ++i)     thsafe[i]= std::thread(dPhiAccumulator,i,std::cref(index),std::cref(MBD_cen),branch11,branch16,ClusPhi,ClusLayer);
        for(int i = 0; i < 8; ++i)     thsafe[i].join();
        angularPlot1D(h_dPhi_nomix, sub_options, "dPhi of unmixed");

        TFile* file = TFile::Open("dPhi.root", "UPDATE");
        if (!file || file -> IsZombie()) {
            std::cerr << "You messed up!" << std::endl;
            exit(1);
        }
        h_dPhi_nomix -> Write("nomix", TObject::kOverwrite);
    }
    else if (sub_options[0] == "perCen") {
        std::thread thsafe[14];
        std::cout << "safe dPhi of different centralities" << std::endl;
        for (int i = 0; i < 14; i++)
            thsafe[i] = std::thread(dPhi_in_bins_of_Centrality_ver1,i,target,std::cref(index),std::cref(MBD_true_z),std::cref(MBD_cen),branch11,branch16,ClusPhi,ClusLayer);

        for (int i = 0; i < 14; i++)
            thsafe[i].join();

        std::vector<TH1D*> h(h_CenonOne, h_CenonOne + 14);
        ArrayPlot1D_Rescale(h, sub_options, "dPhi_per_centralities_rescale");
    }
    else if (sub_options[0] == "perZ") {
        std::thread thsafe[20];
        std::cout << "safe dPhi of different centralities" << std::endl;
        for (int i = 0; i < 20; i++)
            thsafe[i] = std::thread(dPhi_in_bins_of_Z_vtx_ver1,i,target,std::cref(index),std::cref(MBD_true_z),std::cref(MBD_cen),branch11,branch16,ClusPhi,ClusLayer);

        for (int i = 0; i < 20; i++)
            thsafe[i].join();

        std::vector<TH1D*> h(h_ZonOne, h_ZonOne + 20);
        ArrayPlot1D_Rescale_ver2(h, sub_options, "dPhi_per_Z_vtx_rescale");
    }
    // TCanvas *c1 = new TCanvas("c1", "dPhi Histogram", 1920, 1056);
    // h_dPhi_nomix -> Draw();
    // c1->Update();    c1->Modified();

    // Stop the stopwatch and print the runtime:
    timer.Stop();   timer.Print();

    theApp.Run();
    return 0;
}
