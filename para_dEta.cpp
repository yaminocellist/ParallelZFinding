#include "Math/WrappedMultiTF1.h"
#include "Math/WrappedParamFunction.h"
#include "TROOT.h"
#include "TVirtualFitter.h"
#include "TCanvas.h"
#include "TApplication.h"
#include <TStopwatch.h>
#include <thread>
#include <mutex>
#include "headers/histogramPlotting.h"
#include "headers/CommonParaFuncs.h"

std::mutex m_mutex;
int N = 1000;
double range_min = -2*M_PI;
double range_max = 2*M_PI;

double bin_width = (range_max - range_min) / N;
TH1D *h_Eta_diff = new TH1D("", "", N, range_min*1e9, range_max*1e9);
TH1D *h_dEta_nomix = new TH1D("dEta of no-mixing events", ";dEta;# of counts", N, range_min, range_max);
TH2D *h_Eta_Z = new TH2D("Eta v.s. Z", "Eta distributions of all events;Eta Value;MBD Z VTX [cm]", N, range_min, range_max, 20, -25., -5.);
TH2D *h_dEta_Z = new TH2D("dEta v.s. Z", ";dEta Value;MBD Z VTX [cm]", N, range_min, range_max, 20, -25., -5.);
TH2D *h_dEta_cen = new TH2D("dEta v.s. centrality", ";dEta Value;Centrality", N, range_min, range_max, 14, .0, .7);
TH1D *h_dEta_Background = new TH1D("", "", N, range_min, range_max);

TH1D* h_dEta_cenOnOne[14];
TH1D* h_dEta_ZOnOne[20];

void Eta_diff_gathering (int id, int chunck_size, std::vector<int> index, TBranch *branch11, TBranch* branch14,TBranch* branch15, TBranch* branch17, std::vector<int>* ClusLayer, std::vector<float>* ClusZ, std::vector<float>* ClusR, std::vector<float>* ClusEta) {
    double z_vtx, cen, dZ, R, theta, eta, dEta;
    int local_index;
    std::unique_lock<std::mutex> lock(m_mutex);
    for (int i = 0; i < chunck_size; i++) {
        local_index = i + chunck_size*id;
        branch11->GetEntry(index[local_index]);   // ClusLayer;
        branch14->GetEntry(index[local_index]);   // ClusZ;
        branch15->GetEntry(index[local_index]);   // ClusR;
        branch17->GetEntry(index[local_index]);   // ClusEta;
        for (int j = 0; j < ClusZ->size(); j++) {
            dZ       = ClusZ->at(j);
            R        = ClusR->at(j);
            theta    = std::atan2(R, dZ);
            if (dZ >= 0)    eta = -std::log(std::tan(theta/2));
            if (dZ <  0)    eta = std::log(std::tan((M_PI - theta)/2));
            h_Eta_diff -> Fill(eta - ClusEta->at(j));
        }
    }
    lock.unlock();
}

void dEta_per_centrality (
    const int &id, 
    const int &target, 
    const std::vector<int>    &index,
    const std::vector<double> &MBD_true_z,
    const std::vector<double> &MBD_cen,
    TBranch* branch11, 
    TBranch* branch14, 
    TBranch* branch15,
    const std::vector<int>*   const ClusLayer, 
    const std::vector<float>* const ClusZ, 
    const std::vector<float>* const ClusR
) {
    double cen_lower_range = static_cast<double>(id)*0.05;
    double cen_upper_range = static_cast<double>(id + 1)*0.05;
    h_dEta_cenOnOne[id] = new TH1D(Form("%2.2f < Centrality < %2.2f", cen_lower_range, cen_upper_range), ";dEta value;# of counts", N, range_min, range_max);
    double z_vtx, cen, dZ, R, theta, eta, dEta;
    std::vector<double> Eta0, Eta1;
    std::unique_lock<std::mutex> lock(m_mutex);
    for (int i = 0; i < target; i++) {
        cen = MBD_cen[i];
        z_vtx = MBD_true_z[i];
        if (cen >= cen_lower_range && cen <= cen_upper_range) {
            branch11->GetEntry(index[i]);   // ClusLayer;
            branch14->GetEntry(index[i]);   // ClusZ;
            branch15->GetEntry(index[i]);   // ClusR;
            for (int j = 0; j < ClusZ->size(); j++) {
                dZ       = ClusZ->at(j) - z_vtx;
                R        = ClusR->at(j);
                theta    = std::atan2(R, dZ);
                if (dZ >= 0)    eta = -std::log(std::tan(theta/2));
                if (dZ <  0)    eta = std::log(std::tan((M_PI - theta)/2));
                if (ClusLayer->at(j) == 3 || ClusLayer->at(j) == 4) {
                    Eta0.push_back(eta);
                }
                else {
                    Eta1.push_back(eta);
                }
            }
            for (double eta0: Eta0) {
                for (double eta1: Eta1) {
                    dEta = eta0 - eta1;
                    h_dEta_cenOnOne[id] -> Fill(dEta);
                }
            }
            Eta0.clear();   Eta1.clear();
        }
    }
    lock.unlock();
}

void dEta_per_centrality2 (
    const int &id, 
    const int &target, 
    const std::vector<int>    &index,
    const std::vector<double> &MBD_true_z,
    const std::vector<double> &MBD_cen,
    TBranch* branch11, 
    TBranch* branch14, 
    TBranch* branch15,
    TBranch* branch16,
    const std::vector<int>*   const ClusLayer, 
    const std::vector<float>* const ClusZ, 
    const std::vector<float>* const ClusR,
    const std::vector<float>* const ClusPhi
) {
    double cen_lower_range = static_cast<double>(id)*0.05;
    double cen_upper_range = static_cast<double>(id + 1)*0.05;
    h_dEta_cenOnOne[id] = new TH1D(Form("%2.2f < Centrality < %2.2f", cen_lower_range, cen_upper_range), ";dEta value;# of counts", N, range_min, range_max);
    double z_vtx, cen, dZ, R, theta, eta, dEta, phi, dPhi;
    std::vector<EtaWithPhi> Eta0, Eta1;
    std::unique_lock<std::mutex> lock(m_mutex);
    for (int i = 0; i < target; i++) {
        cen = MBD_cen[i];
        z_vtx = MBD_true_z[i];
        if (cen >= cen_lower_range && cen <= cen_upper_range) {
            branch11->GetEntry(index[i]);   // ClusLayer;
            branch14->GetEntry(index[i]);   // ClusZ;
            branch15->GetEntry(index[i]);   // ClusR;
            branch16->GetEntry(index[i]);   // ClusPhi;
            for (int j = 0; j < ClusZ->size(); j++) {
                dZ       = ClusZ->at(j) - z_vtx;
                R        = ClusR->at(j);
                theta    = std::atan2(R, dZ);
                phi      = ClusPhi->at(j);
                if (dZ >= 0)    eta = -std::log(std::tan(theta/2));
                if (dZ <  0)    eta = std::log(std::tan((M_PI - theta)/2));
                if (ClusLayer->at(j) == 3 || ClusLayer->at(j) == 4) {
                    Eta0.emplace_back(eta, phi);
                }
                else {
                    Eta1.emplace_back(eta, phi);
                }
            }
            for (const EtaWithPhi& eta0: Eta0) {
                for (const EtaWithPhi& eta1: Eta1) {
                    dEta = eta0.eta_value - eta1.eta_value;
                    dPhi = eta0.phi_value - eta1.phi_value;
                    if (std::abs(dPhi) < dPhi_cut) {
                        h_dEta_cenOnOne[id] -> Fill(dEta);
                    }
                }
            }
            Eta0.clear();   Eta1.clear();
        }
    }
    lock.unlock();
}

void dEta_per_Z (
    const int &id, 
    const int &target, 
    const std::vector<int>    &index,
    const std::vector<double> &MBD_true_z,
    TBranch* branch11, 
    TBranch* branch14, 
    TBranch* branch15,
    const std::vector<int>*   const ClusLayer, 
    const std::vector<float>* const ClusZ, 
    const std::vector<float>* const ClusR
) {
    double z_lower_range = -25. + static_cast<double>(id);
    double z_upper_range = -24. + static_cast<double>(id);
    h_dEta_ZOnOne[id] = new TH1D(Form("%2.2f < MBD Z Vtx < %2.2f", z_lower_range, z_upper_range), ";dEta value;# of counts", N, range_min, range_max);
    double z_vtx, cen, dZ, R, theta, eta, dEta;
    std::vector<double> Eta0, Eta1;
    std::unique_lock<std::mutex> lock(m_mutex);
    for (int i = 0; i < target; i++) {
        z_vtx = MBD_true_z[i];
        if (z_vtx >= z_lower_range && z_vtx <= z_upper_range) {
            branch11->GetEntry(index[i]);   // ClusLayer;
            branch14->GetEntry(index[i]);   // ClusZ;
            branch15->GetEntry(index[i]);   // ClusR;
            int local_size = ClusZ->size();
            // std::cout << local_size << std::endl;
            // Eta0.reserve(local_size);   Eta1.reserve(local_size);
            for (int j = 0; j < local_size; j++) {
                dZ       = ClusZ->at(j) - z_vtx;
                R        = ClusR->at(j);
                theta    = std::atan2(R, dZ);
                if (dZ >= 0)    eta = -std::log(std::tan(theta/2));
                if (dZ <  0)    eta = std::log(std::tan((M_PI - theta)/2));
                if (ClusLayer->at(j) == 3 || ClusLayer->at(j) == 4) {
                    Eta0.push_back(eta);
                }
                else {
                    Eta1.push_back(eta);
                }
            }
            for (double eta0: Eta0) {
                for (double eta1: Eta1) {
                    dEta = eta0 - eta1;
                    h_dEta_ZOnOne[id] -> Fill(dEta);
                }
            }
            Eta0.clear();   Eta1.clear();
        }
    }
    lock.unlock();
}

void dEta_per_Z2 (
    const int &id, 
    const int &target, 
    const std::vector<int>    &index,
    const std::vector<double> &MBD_true_z,
    TBranch* branch11, 
    TBranch* branch14, 
    TBranch* branch15,
    TBranch* branch16,
    const std::vector<int>*   const ClusLayer, 
    const std::vector<float>* const ClusZ, 
    const std::vector<float>* const ClusR,
    const std::vector<float>* const ClusPhi
) {
    double z_lower_range = -25. + static_cast<double>(id);
    double z_upper_range = -24. + static_cast<double>(id);
    h_dEta_ZOnOne[id] = new TH1D(Form("%2.2f < MBD Z Vtx < %2.2f", z_lower_range, z_upper_range), ";dEta value;# of counts", N, range_min, range_max);
    double z_vtx, cen, dZ, R, theta, eta, dEta, phi, dPhi;
    std::vector<EtaWithPhi> Eta0, Eta1;
    std::unique_lock<std::mutex> lock(m_mutex);
    for (int i = 0; i < target; i++) {
        z_vtx = MBD_true_z[i];
        if (z_vtx >= z_lower_range && z_vtx <= z_upper_range) {
            branch11->GetEntry(index[i]);   // ClusLayer;
            branch14->GetEntry(index[i]);   // ClusZ;
            branch15->GetEntry(index[i]);   // ClusR;
            branch16->GetEntry(index[i]);   // ClusPhi;
            int local_size = ClusZ->size();
            // std::cout << local_size << std::endl;
            // Eta0.reserve(local_size);   Eta1.reserve(local_size);
            for (int j = 0; j < local_size; j++) {
                dZ       = ClusZ->at(j) - z_vtx;
                R        = ClusR->at(j);
                theta    = std::atan2(R, dZ);
                phi      = ClusPhi->at(j);
                if (dZ >= 0)    eta = -std::log(std::tan(theta/2));
                if (dZ <  0)    eta = std::log(std::tan((M_PI - theta)/2));
                if (ClusLayer->at(j) == 3 || ClusLayer->at(j) == 4) {
                    Eta0.emplace_back(eta, phi);
                }
                else {
                    Eta1.emplace_back(eta, phi);
                }
            }
            for (const EtaWithPhi& eta0: Eta0) {
                for (const EtaWithPhi& eta1: Eta1) {
                    dEta = eta0.eta_value - eta1.eta_value;
                    dPhi = eta0.phi_value - eta1.phi_value;
                    if (std::abs(dPhi) < dPhi_cut) {
                        h_dEta_ZOnOne[id] -> Fill(dEta);
                    }
                }
            }
            Eta0.clear();   Eta1.clear();
        }
    }
    lock.unlock();
}

void dEtaNoMix (
    const int &id,
    const int &chunck_size,
    const std::vector<int> &index,
    const std::vector<double> &MBD_true_z,
    const std::vector<double> &MBD_cen,
    TBranch* branch11,
    TBranch* branch14,
    TBranch* branch15,
    const std::vector<int>*   const ClusLayer, 
    const std::vector<float>* const ClusZ,
    const std::vector<float>* const ClusR
) {
    double z_vtx, cen, dZ, R, theta, eta, dEta;
    std::vector<double> Eta0, Eta1;
    int local_index;
    // std::cout << chunck_size << std::endl;
    std::unique_lock<std::mutex> lock(m_mutex);
    for (int i = 0; i < chunck_size; i++) {
        local_index = i + chunck_size*id;
        // std::cout << local_index << std::endl;
        branch11->GetEntry(index[local_index]);   // ClusLayer;
        branch14->GetEntry(index[local_index]);   // ClusZ;
        branch15->GetEntry(index[local_index]);   // ClusR;
        z_vtx = MBD_true_z[local_index];
        cen = MBD_cen[local_index];
        for (int j = 0; j < ClusZ->size(); j++) {
            dZ       = ClusZ->at(j) - z_vtx;
            R        = ClusR->at(j);
            theta    = std::atan2(R, dZ);
            if (dZ >= 0)    eta = -std::log(std::tan(theta/2));
            if (dZ <  0)    eta = std::log(std::tan((M_PI - theta)/2));
            // h_Eta_Z -> Fill(eta, z_vtx);
            if (ClusLayer->at(j) == 3 || ClusLayer->at(j) == 4) {
                Eta0.push_back(eta);
            }
            else {
                Eta1.push_back(eta);
            }
        }
        for (int k = 0; k < Eta0.size(); k++) {
            for (int l = 0; l < Eta1.size(); l++) {
                dEta = Eta0[k] - Eta1[l];
                // h_dEta_nomix -> Fill(dEta);
                h_dEta_Z     -> Fill(dEta, z_vtx);
                h_dEta_cen   -> Fill(dEta, cen);
            }
        }
        Eta0.clear();   Eta1.clear();
    }
    lock.unlock();
}

void dEtaNoMix2 (
    const int &id,
    const int &chunck_size,
    const std::vector<int> &index,
    const std::vector<double> &MBD_true_z,
    const std::vector<double> &MBD_cen,
    TBranch* branch11,
    TBranch* branch14,
    TBranch* branch15,
    TBranch* branch16,
    const std::vector<int>*   const ClusLayer, 
    const std::vector<float>* const ClusZ,
    const std::vector<float>* const ClusR,
    const std::vector<float>* const ClusPhi
) {
    double z_vtx, cen, dZ, R, theta, eta, dEta, phi, dPhi;
    std::vector<EtaWithPhi> Eta0, Eta1;
    int local_index;
    std::unique_lock<std::mutex> lock(m_mutex);
    for (int i = 0; i < chunck_size; i++) {
        local_index = i + chunck_size*id;
        branch11->GetEntry(index[local_index]);   // ClusLayer;
        branch14->GetEntry(index[local_index]);   // ClusZ;
        branch15->GetEntry(index[local_index]);   // ClusR;
        branch16->GetEntry(index[local_index]);   // ClusPhi;
        z_vtx = MBD_true_z[local_index];
        cen = MBD_cen[local_index];
        for (int j = 0; j < ClusZ->size(); j++) {
            dZ       = ClusZ->at(j) - z_vtx;
            R        = ClusR->at(j);
            theta    = std::atan2(R, dZ);
            phi      = ClusPhi->at(j);
            if (dZ >= 0)    eta = -std::log(std::tan(theta/2));
            if (dZ <  0)    eta = std::log(std::tan((M_PI - theta)/2));
            // h_Eta_Z -> Fill(eta, z_vtx);
            if (ClusLayer->at(j) == 3 || ClusLayer->at(j) == 4) {
                Eta0.emplace_back(eta, phi);
            }
            else {
                Eta1.emplace_back(eta, phi);
            }
        }
        for (int k = 0; k < Eta0.size(); k++) {
            for (int l = 0; l < Eta1.size(); l++) {
                dEta = Eta0[k].eta_value - Eta1[l].eta_value;
                dPhi = Eta0[k].phi_value - Eta1[l].phi_value;
                if (std::abs(dPhi) < dPhi_cut) {
                    // h_dEta_nomix -> Fill(dEta);
                    h_dEta_Z     -> Fill(dEta, z_vtx);
                    // h_dEta_cen   -> Fill(dEta, cen);
                }
            }
        }
        Eta0.clear();   Eta1.clear();
    }
    lock.unlock();
}

void dEtaMixing (
    const int &starter,
    const int &ender,
    const std::vector<std::vector<double>> &event_Eta0,
    const std::vector<std::vector<double>> &event_Eta1
) {
    double dEta;
    std::unique_lock<std::mutex> lock(m_mutex);
    for (int i = starter; i <= ender; i++) {
        for (int j = i; j < event_Eta1.size(); j++) {
            for (double eta0: event_Eta0[i]) {
                for (double eta1: event_Eta1[j]) {
                    dEta = eta0 - eta1;
                    h_dEta_Background -> Fill(dEta);
                }
            }
            for (double eta1: event_Eta1[i]) {
                for (double eta0: event_Eta0[j]) {
                    dEta = eta0 - eta1;
                    h_dEta_Background -> Fill(dEta);
                }
            }
        }
    }
    lock.unlock();
}

void dEtaMixing2 (
    const int &starter,
    const int &ender,
    const std::vector<std::vector<EtaWithPhi>> &event_Eta0,
    const std::vector<std::vector<EtaWithPhi>> &event_Eta1
) {
    double dEta, dPhi;
    std::unique_lock<std::mutex> lock(m_mutex);
    for (int i = starter; i <= ender; i++) {
        for (int j = i; j < event_Eta1.size(); j++) {
            for (const EtaWithPhi& eta0: event_Eta0[i]) {
                for (const EtaWithPhi& eta1: event_Eta1[j]) {
                    dEta = eta0.eta_value - eta1.eta_value;
                    dPhi = eta0.phi_value - eta1.phi_value;
                    if (std::abs(dPhi) < dPhi_cut)  h_dEta_Background -> Fill(dEta);
                }
            }
            for (const EtaWithPhi& eta1: event_Eta1[i]) {
                for (const EtaWithPhi& eta0: event_Eta0[j]) {
                    dEta = eta0.eta_value - eta1.eta_value;
                    dPhi = eta0.phi_value - eta1.phi_value;
                    if (std::abs(dPhi) < dPhi_cut)  h_dEta_Background -> Fill(dEta);
                }
            }
        }
    }
    lock.unlock();
}

int main(int argc, char* argv[]) {
    // Start the stopwatch:
    TStopwatch timer;   timer.Start();
    current_PC_time();
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
    std::vector<float> *ClusR     = nullptr;    
    std::vector<float> *ClusPhi   = nullptr;    
    std::vector<float> *ClusEta   = nullptr;

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
    std::cout << "Num of events: " << target << ", centrality greater than " << cen_low << ", smaller than " << cen_high << ", z vtx greater than " << -z_low << ", smaller than " << -z_high << std::endl;
    ROOT::EnableThreadSafety();
    if (method[0] == "nomix") {
        std::thread thsafe[8];
        std::cout<<"dEta without mixing:"<<std::endl;
        // for(int i = 0; i < 8; ++i)     thsafe[i]= std::thread(dEtaNoMix,i,target/8,std::cref(index),std::cref(MBD_true_z),std::cref(MBD_cen),branch11,branch14,branch15,ClusLayer,ClusZ,ClusR);
        for(int i = 0; i < 8; ++i)     thsafe[i]= std::thread(dEtaNoMix2,i,target/8,std::cref(index),std::cref(MBD_true_z),std::cref(MBD_cen),branch11,branch14,branch15,branch16,ClusLayer,ClusZ,ClusR,ClusPhi);
        for(int i = 0; i < 8; ++i)     thsafe[i].join();

        // save_histo_toRootFile(h_dEta_nomix, method, "dEta");
        // h_dEta_nomix->SetTitle(Form("dEta of %d no mixing events with dPhi cut < %0.2f", target, dPhi_cut));
        // angularPlot1D(h_dEta_nomix, method, "dEta of unmixed with dPhi cut");
        // angularPlot2D(h_Eta_Z, method, "Eta v.s. Z distribution");
        h_dEta_Z->SetTitle(Form("dEta v.s. MBD Z of %d no mixing events with dPhi cut < %0.2f", target, dPhi_cut));
        angularPlot2D(h_dEta_Z, method, "dEta v.s. Z distribution with dPhi cut");
        // h_dEta_cen->SetTitle(Form("dEta v.s. MBD centrality of %d no mixing events < %0.2f", target, dPhi_cut));
        // angularPlot2D(h_dEta_cen, method, "dEta v.s. Centrality distribution with dPhi cut");
    }
    else if (method[0] == "san") {
        std::thread thsafe[8];
        std::cout<<"Sanity Check:"<<std::endl;
        for(int i = 0; i < 8; ++i)     thsafe[i]= std::thread(Eta_diff_gathering,i,target/8,std::cref(index),branch11,branch14,branch15,branch17,ClusLayer,ClusZ,ClusR,ClusEta);
        for(int i = 0; i < 8; ++i)     thsafe[i].join();

        h_Eta_diff->SetTitle(Form("Eta sanity check of %d no mixing events", target));
        angularPlot1D(h_Eta_diff, method, "Eta sanity check");
        histogramOverflowCheck(h_Eta_diff);
    }
    else if (method[0] == "perCen") {
        std::thread thsafe[14];
        std::cout << "dEta of different centralities:" << std::endl;
        std::cout << target << std::endl;
        // for(int i = 0; i < 14; i++)     thsafe[i] = std::thread(dEta_per_centrality,i,target,std::cref(index),std::cref(MBD_true_z),std::cref(MBD_cen),branch11,branch14,branch15,ClusLayer,ClusZ,ClusR);
        for(int i = 0; i < 14; i++)     thsafe[i] = std::thread(dEta_per_centrality2,i,target,std::cref(index),std::cref(MBD_true_z),std::cref(MBD_cen),branch11,branch14,branch15,branch16,ClusLayer,ClusZ,ClusR,ClusPhi);
        for(int i = 0; i < 14; i++)     thsafe[i].join();
        
        for (int i = 0; i < 14; i++)
            h_dEta_cenOnOne[i]->SetTitle(Form("dEta of %d events in different centrality ranges with dPhi cut < %0.2f", target, dPhi_cut));
        std::vector<TH1D*> h(h_dEta_cenOnOne, h_dEta_cenOnOne + 14);
        // ArrayPlot_dEta_1D_Logy(h, method, "dEta per centralities with dPhi cut");
        ArrayPlot1D_Rescale_dEta(h, method, "dEta per centralities Rescaled with dPhi cut");
    }
    else if (method[0] == "perZ") {
        std::thread thsafe[20];
        std::cout << "dEta of different centralities:" << std::endl;
        // for(int i = 0; i < 20; i++)     thsafe[i] = std::thread(dEta_per_Z,i,target,std::cref(index),std::cref(MBD_true_z),branch11,branch14,branch15,ClusLayer,ClusZ,ClusR);
        for(int i = 0; i < 20; i++)     thsafe[i] = std::thread(dEta_per_Z2,i,target,std::cref(index),std::cref(MBD_true_z),branch11,branch14,branch15,branch16,ClusLayer,ClusZ,ClusR,ClusPhi);
        for(int i = 0; i < 20; i++)     thsafe[i].join();
        
        for (int i = 0; i < 20; i++)
            h_dEta_ZOnOne[i]->SetTitle(Form("dEta of %d events in different MBD Z Vtx ranges, with dPhi cut < %0.2f", target, dPhi_cut));
        TH1D* dummyHistogram = new TH1D("dummy", "Dummy for empty Histogram", N, range_min, range_max);
        std::vector<TH1D*> h;
        for (int i = 0; i < 20; i++)
            h.push_back(h_dEta_ZOnOne[i]);
        // ArrayPlot_dEta_1D_Logy(h, method, "dEta per Z");
        ArrayPlot1D_Rescale_dEta(h, method, "dEta per Z Rescaled with dPhi cut");
        // TCanvas *c1 = new TCanvas("c1dEta", "dPhi Histogram", 1920, 1056);
        // c1->Update();    c1->Modified();
    }
    else if (method[0] == "mix") {
        double z_vtx, cen, dZ, R, theta, eta, dEta, phi, dPhi;
        std::vector<EtaWithPhi> Eta0, Eta1;
        std::vector<std::vector<EtaWithPhi>> event_Eta0, event_Eta1;
        // std::vector<double> Eta0, Eta1;
        // std::vector<std::vector<double>> event_Eta0, event_Eta1;
        TH1D *h_dEta_Signal = new TH1D("no mixing events", ";dEta value;# of counts", N, range_min, range_max);
        for (int i = 0; i < target; i++) {
            branch11->GetEntry(index[i]);   // ClusLayer;
            branch14->GetEntry(index[i]);   // ClusZ;
            branch15->GetEntry(index[i]);   // ClusR;
            branch16->GetEntry(index[i]);   // ClusPhi;
            z_vtx = MBD_true_z[i];
            // event_Eta0.push_back(std::vector <double>());   event_Eta1.push_back(std::vector <double>());
            event_Eta0.emplace_back(std::vector <EtaWithPhi>());   event_Eta1.emplace_back(std::vector <EtaWithPhi>());
            for (int j = 0; j < ClusZ->size(); j++) {
                dZ    = ClusZ->at(j) - z_vtx;
                R     = ClusR->at(j);
                theta = std::atan2(R, dZ);
                phi   = ClusPhi->at(j);
                if (dZ >= 0)    eta = -std::log(std::tan(theta/2));
                if (dZ <  0)    eta = std::log(std::tan((M_PI - theta)/2));
                if (ClusLayer->at(j) == 3 || ClusLayer->at(j) == 4) {
                    Eta0.emplace_back(eta, phi);
                    event_Eta0[i].emplace_back(eta, phi);
                }
                else {
                    Eta1.emplace_back(eta, phi);
                    event_Eta1[i].emplace_back(eta, phi);
                }
            }
            for (const EtaWithPhi& eta0: Eta0) {
                for (const EtaWithPhi& eta1: Eta1) {
                    dEta = eta0.eta_value - eta1.eta_value;
                    dPhi = eta0.phi_value - eta1.phi_value;
                    if (std::abs(dPhi) < dPhi_cut)  h_dEta_Signal -> Fill(dEta);
                }
            }
            Eta0.clear();   Eta1.clear();
        }

        std::vector<int> boundaries = readCsvToVector("/Users/yaminocellist/codeGarage/boundaries.csv");

        std::thread thsafe[8];
        // for(int i = 0; i < 8; i++)    thsafe[i] = std::thread(dEtaMixing,boundaries[2*i], boundaries[2*i+1],event_Eta0,event_Eta1);
        for(int i = 0; i < 8; i++)    thsafe[i] = std::thread(dEtaMixing2,boundaries[2*i], boundaries[2*i+1],event_Eta0,event_Eta1);
        for(int i = 0; i < 8; i++)    thsafe[i].join();

        TCanvas *c1 = new TCanvas("c1", "dPhi Histogram", 1920, 1056);
        c1 -> Divide(1, 2);
        c1 -> cd(1);

        // double eta_range_low = -3.4, eta_range_high = -3.1;
        // // double eta_range_low = -2., eta_range_high = -1.9375;
        // int bin_range_low = h_dEta_Background->FindBin(eta_range_low), bin_range_high = h_dEta_Background->FindBin(eta_range_high);
        // double max_unmixed = -1, max_mixed = -1, current_binContent;
        // for (int bin = bin_range_low; bin <= bin_range_high; bin++) {
        //     current_binContent = h_dEta_Signal->GetBinContent(bin);
        //     if (max_unmixed < current_binContent)  max_unmixed = current_binContent;
        //     current_binContent = h_dEta_Background->GetBinContent(bin);
        //     if (max_mixed < current_binContent)  max_mixed = current_binContent;
        // }
        // double N = max_mixed/max_unmixed;

        double N1 = h_dEta_Background -> Integral(h_dEta_Background->FindFixBin(-2*M_PI), h_dEta_Background->FindFixBin(-0.2), "") + 
                    h_dEta_Background -> Integral(h_dEta_Background->FindFixBin(0.1), h_dEta_Background->FindFixBin(2*M_PI), "");
        double N2 = h_dEta_Signal -> Integral(h_dEta_Signal->FindFixBin(-2*M_PI), h_dEta_Signal->FindFixBin(-0.2), "") + 
                    h_dEta_Signal -> Integral(h_dEta_Signal->FindFixBin(0.1), h_dEta_Signal->FindFixBin(2*M_PI), "");
        double N  = N1/N2;

        h_dEta_Signal -> Add(h_dEta_Signal, N-1);
        h_dEta_Signal -> SetLineColor(2);   h_dEta_Signal -> SetLineWidth(3);
        h_dEta_Signal -> Draw("SAME");
        h_dEta_Signal -> SetTitle(Form("dEta of %d mixing events, %2.2f < centrality < %2.2f, %2.2fcm < MBD Z VTX < %2.2fcm", target, cen_low, cen_high, -z_low, -z_high));
        h_dEta_Signal -> GetXaxis() -> SetTitleSize(.05);
        h_dEta_Signal -> GetYaxis() -> SetTitleSize(.05);
        h_dEta_Signal->GetXaxis()->CenterTitle(true);   h_dEta_Signal->GetYaxis()->CenterTitle(true);
        h_dEta_Background->Draw("same");    h_dEta_Background -> SetLineWidth(3);
        h_dEta_Background->GetXaxis()->SetRangeUser(-M_PI*5/4, M_PI*5/4);
        h_dEta_Signal->GetXaxis()->SetRangeUser(-M_PI*5/4, M_PI*5/4);

        TLegend *lg = new TLegend(0.12, 0.8, 0.36, 0.9);
        lg -> AddEntry(h_dEta_Signal, "dEta of non-mixing events with dPhi cut", "l");
        lg -> AddEntry(h_dEta_Background, "dEta of mixing events with dPhi cut", "l");
        gStyle -> SetLegendTextSize(.041);
        lg->Draw("same");

        c1 -> cd(2);
        TH1D* h_dEta_Normalized = (TH1D*) h_dEta_Signal->Clone("Normalized hBackground");
        h_dEta_Normalized -> SetTitle("Background Subtracted");
        h_dEta_Normalized -> Add(h_dEta_Background, -1);
        h_dEta_Normalized -> Sumw2();
        h_dEta_Normalized -> Draw("HIST SAME");
        h_dEta_Normalized -> Draw("e1psame");
        h_dEta_Normalized -> GetXaxis() -> SetTitleSize(.05);
        h_dEta_Normalized -> SetFillColor(kYellow - 7);
        h_dEta_Normalized -> SetFillStyle(1001);
        c1->Update();    c1->Modified();
        c1 -> SaveAs(Form("../../External/zFindingPlots/dEta_mixed_%d_%2.2f_%2.2f_%2.2f_%2.2f_with_dPhi_cut.png", target, cen_low, cen_high, -z_low, -z_high));
    }

    // Stop the stopwatch and print the runtime:
    timer.Stop();   timer.Print();

    theApp.Run();
    return 0;
}