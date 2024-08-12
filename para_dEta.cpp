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
TH1D *h_dEta_nomix = new TH1D("", "", N, range_min, range_max);
TH1D *h_Background_dEta = new TH1D("", "", N, range_min, range_max);

void dEtaNoMix (int id, int chunck_size, std::vector<int> index, std::vector<double> MBD_true_z, TBranch *branch11, TBranch* branch14,TBranch* branch15, std::vector<int>* ClusLayer, std::vector<float>* ClusZ, std::vector<float>* ClusR) {
    double z_vtx, dZ, R, theta, eta, dEta;
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
        for (int k = 0; k < Eta0.size(); k++) {
            for (int l = 0; l < Eta1.size(); l++) {
                dEta = Eta0[k] - Eta1[l];
                h_dEta_nomix -> Fill(dEta);
                // h_dPhi_Z       -> Fill(dPhi, fz);
                // h_dPhi_cen     -> Fill(dPhi, cen);
            }
        }
        Eta0.clear();   Eta1.clear();
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
    if (method[0] == "nomix") {
        std::thread thsafe[8];
        std::cout<<"multi-thready safe:"<<std::endl;
        for(int i = 0; i < 8; ++i)     thsafe[i]= std::thread(dEtaNoMix,i,target/8,std::cref(index),std::cref(MBD_true_z),branch11,branch14,branch15,ClusLayer,ClusZ,ClusR);
        for(int i = 0; i < 8; ++i)     thsafe[i].join();

        save_histo_toRootFile(h_dEta_nomix, method, "dEta");
        angularPlot1D(h_dEta_nomix, method, "dEta of unmixed");   
    }

    // Stop the stopwatch and print the runtime:
    timer.Stop();   timer.Print();

    theApp.Run();
    return 0;
}