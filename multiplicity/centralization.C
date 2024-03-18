#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <vector>
#include <stdlib.h>
#include "TH1.h"
#include "TGraph.h"
#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TLatex.h"

using namespace std;

/*********************************************************************************
 *  Beginning of the main algorithm
 *  Dealing with data ony whose VZERO value is up to 5% of the whole distribution
 * *******************************************************************************/

//*** c:\root\bin\thisroot.bat
//*** 1V1MIT\Week3

void centralization () {
    ifstream myfile("../week5/V40_50.txt");
    ifstream phiRaw("../week3/phiRaw.txt");
    ifstream etaRaw("../week3/etaRaw.txt");
    
    if (!myfile.is_open()){
		  cout << "Unable to open linelabel" << endl;
		  system("read -n 1 -s -p \"Press any key to continue...\" echo");
		  exit(1);
 	}
 	  if (!etaRaw.is_open()){
		cout << "Unable to open etaRaw" << endl;
		system("read -n 1 -s -p \"Press any key to continue...\" echo");
		exit(1);
 	}
    if (!phiRaw.is_open()){
		cout << "Unable to open phiRaw" << endl;
		system("read -n 1 -s -p \"Press any key to continue...\" echo");
		exit(1);
 	}
  /*******************************
	*       Sets Variable bins;
	* ******************************/
  const double pi = 3.142;
	/* Double_t binEdges[481];
	double dx_1 = 0.01*pi;
	binEdges[0] = -2*pi;
	for (int i = 1; i <= 190; i++) {
    	binEdges[i] = binEdges[0] + i*dx_1;
  }
	double dx_2 = 0.002*pi;
	for (int j = 1; j <= 100; j++) {
    	binEdges[j+190] = binEdges[190] + j*dx_2;
  }
	for (int l = 1; l <= 190; l++) {
    	binEdges[l + 290] = binEdges[290] + l*dx_1;
  } */

  Double_t binEdges[481];
	double dx_1 = 0.01;
	binEdges[0] = -2;
	for (int i = 1; i <= 190; i++) {
    	binEdges[i] = binEdges[0] + i*dx_1;
  }
	double dx_2 = 0.002;
	for (int j = 1; j <= 100; j++) {
    	binEdges[j+190] = binEdges[190] + j*dx_2;
  }
	for (int l = 1; l <= 190; l++) {
    	binEdges[l + 290] = binEdges[290] + l*dx_1;
  }
    
    TH1D *HSignal = new TH1D("HSignal & HBackground","HSignal & HBackground; Values; Entries",480, binEdges);
    HSignal -> SetMarkerStyle(20);
    HSignal -> SetLineColor(kRed);
    TH1D *HBackground = new TH1D("HBackground","HBackground; Values; Entries",480, binEdges);
    HBackground -> SetMarkerStyle(21);
    HBackground -> SetLineColor(kBlue);

    TH1D *HEta = new TH1D("HEta","Delta Eta; Values; Entries",400,-2,2);
    HEta -> SetMarkerStyle(20);
    HEta -> SetLineColor(kRed);
    TH1D *HPhi = new TH1D("HPhi","Delta Phi; Values; Entries",400,-6.284,6.284);
    HPhi -> SetMarkerStyle(21);
    HPhi -> SetLineColor(kBlue);

    //* 'line_0' indicates on which line stores "hits VZERO" of layer_0;
    //* 'line_1' indicates on which line stores "hits VZERO" of layer_1;
    int line_0, line_1;
    //* 'hits_0' indicates how many lines of data in layer_0;
    //* ‘hits_1’ indicates how many lines of data in layer_1;
    double hits_0, hits_1;
    /*
     * 'eta_0' is a vector stores all eta values in layer_0;
     * 'eta_1' is a vector stores all eta values in layer_1;
     * 'phi_0' is a vector stores all phi values in layer_0;
     * 'phi_1' is a vector stores all phi values in layer_1;
     */
    
    /*
    * 'linestream' temporally stores 'line_0' or 'line_1' for 'myfile';
    * 'dummy' temporally stores any useless line;
    * 'dummy1' temporally stores HOW MANY 'hits' in layer_0;
    * 'dummy2' temporally stores HOW MANY 'hits' in layer_1;
    * 'dummy_eta_0' stores data of eta in layer_0;
    * 'dummy_eta_1' stores data of eta in layer_1;
    * 'dummy_phi_0' stores data of phi in layer_0;
    * 'dummy_phi_1' stores data of phi in layer_1;
    */
    string linestream, dummy, dummy1, dummy2, dummy_eta_0, dummy_eta_1, dummy_phi_0, dummy_phi_1;
    double vzero;
    int indicator = 1;
    double events = 0;

    line_1 = 0;
    hits_1 = 0;

    double phiS = 0.1;
    double phiB = 0.2;

    double Eta0, Phi0, Eta1, Phi1, delphi, deleta, absdphi, absdeta;
    while (indicator <= 1999) {
        getline(myfile, linestream);
        line_0 = stoi(linestream);
        vector<double> eta_0, eta_1, phi_0, phi_1;
    //* throw off useless lines;
        for (int i = 1; i < line_0 - line_1 - hits_1; i++) {
            getline(etaRaw, dummy);
            getline(phiRaw, dummy);
        }
        getline(myfile, linestream);
        line_1 = stoi(linestream);
        //cout << "here comes data between the " << line_0 << "th line and the " << line_1 << "th line." <<endl;
    //* here begins next useful parcel; also double check VZERO value;
        getline(etaRaw, dummy1);
        hits_0 = stod(dummy1);
        getline(phiRaw, dummy2);
        vzero = stod(dummy2);   
        /* if (vzero < 15164) {
        exit(0);
        } */

        for (int ii = 1; ii <= hits_0; ii++) {
            getline(etaRaw, dummy_eta_0);
            eta_0.push_back(stod(dummy_eta_0));
            getline(phiRaw, dummy_phi_0);
            phi_0.push_back(stod(dummy_phi_0));
        }
        getline(etaRaw, dummy1);
        hits_1 = stod(dummy1);
        getline(phiRaw, dummy2);
        /* if (stod(dummy2) == vzero) {
            vzero = stod(dummy2);
        }
        else {
          cout << "You've made mistake." << endl;
          exit(0);
        } */
        /* if (vzero < 15164) {
        exit(0);
        } */
        for (int j = 1; j <= hits_1; j++) {
            getline(etaRaw, dummy_eta_1);
            eta_1.push_back(stod(dummy_eta_1));
            getline(phiRaw, dummy_phi_1);
            phi_1.push_back(stod(dummy_phi_1));
        }

        for (int k = 0; k < hits_0; k++) {
          
        Eta0 = eta_0[k];
        Phi0 = phi_0[k];
          
          for (int l = 0; l < hits_1; l++) {
            Eta1 = eta_1[l];
            
            Phi1 = phi_1[l];
            
            deleta = eta_0[k] - eta_1[l];
            delphi = phi_0[k] - phi_1[l];
            absdphi = abs(delphi);
            absdeta = abs(deleta);
            //HEta -> Fill(deleta);
            //HPhi -> Fill(delphi); 
            if (absdphi <= phiS){
              HSignal -> Fill(deleta);
            }
            else  {
              HBackground -> Fill(deleta);
            }
            /* if (absdeta < 0.05*pi){
              HSignal -> Fill(delphi);
            }
            else if (absdeta < 0.2*pi) {
              HBackground -> Fill(delphi);
            } */
          }
        } 
    



    /* for (auto it = eta_1.begin(); it != eta_1.end(); it++){
        cout << *it << endl;
    } */
        indicator++;
    }
    
    // double N1 = HBackground -> GetEntries();
    // double N2 = HSignal -> GetEntries();

    double N1 = HBackground -> Integral(HBackground->FindFixBin(-2), HBackground->FindFixBin(-phiS), "") + 
                HBackground -> Integral(HBackground->FindFixBin(phiS), HBackground->FindFixBin(2), "");
    double N2 = HSignal -> Integral(HSignal->FindFixBin(-2), HSignal->FindFixBin(-phiS), "") + 
                HSignal -> Integral(HSignal->FindFixBin(phiS), HSignal->FindFixBin(2), "");
    double N  = N2/N1;

    // double N1 = HBackground -> Integral(HBackground->FindFixBin(-2*pi), HBackground->FindFixBin(-0.05*pi), "") + 
    //             HBackground -> Integral(HBackground->FindFixBin(0.05*pi), HBackground->FindFixBin(2*pi), "");
    // double N2 = HSignal -> Integral(HSignal->FindFixBin(-2*pi), HSignal->FindFixBin(-0.05*pi), "") + 
    //             HSignal -> Integral(HSignal->FindFixBin(0.05*pi), HSignal->FindFixBin(2*pi), "");
    


    TCanvas *c1 = new TCanvas();
    c1 -> Divide(1, 2);
    c1 -> cd(1);
    
    TH1D* hb = (TH1D*) HBackground -> Clone("Signal & Normalized Background");
    hb -> Add(HBackground, N - 1);
    hb -> Sumw2();
    hb -> SetFillColor(kRed - 7);
    hb -> SetMarkerSize(0.55);
    hb -> SetFillStyle(1001);
    
    // hb -> GetXaxis() -> SetRange(hb -> FindFixBin(-0.098),hb -> FindFixBin(0.098));
    hb -> GetXaxis() -> SetRange(hb -> FindFixBin(-phiS),hb -> FindFixBin(phiS));
    hb -> SetTitle(Form("Centralized Histogram for %5.0d Events", indicator - 1));
    hb -> SetLabelSize(0.055, "xy");

    HSignal -> Sumw2();
    HSignal -> SetFillColor(kCyan - 7);
    HSignal -> SetMarkerSize(0.55);
    HSignal -> SetFillStyle(3144);
    HSignal -> Draw("HIST SAME");
    HSignal -> Draw("e1psame");
    hb -> Draw("HIST SAME");
    hb -> Draw("e1psame");
    // HSignal -> GetXaxis() -> SetRange(HSignal -> FindFixBin(-0.098),HSignal -> FindFixBin(0.098));
    HSignal -> GetXaxis() -> SetRange(HSignal -> FindFixBin(-phiS+dx_2),HSignal -> FindFixBin(phiS));
    HSignal -> SetTitle(Form("Centralized Histogram for %5.0d Events", indicator - 1));
    HSignal -> SetLabelSize(0.055, "xy");

    TLegend *legend = new TLegend(0.29, 0.15, 0.71, 0.28);
    legend -> AddEntry(HSignal, Form("HSignal for |#Delta#phi| < %1.2f", phiS), "f");
    legend -> AddEntry(hb, Form("Normalized HBackground for %1.2f < |#Delta#phi| < %1.2f", phiS, phiB), "f");
    //legend -> SetHeader("Red Region:");
    gStyle -> SetLegendTextSize(0.05);
    legend -> Draw("same");
    gStyle->SetEndErrorSize(2);
    gStyle->SetErrorX(0.5);

    gPad ->SetLogy(1);
    c1 -> cd(2);

    TH1D* hDiff = (TH1D*) HSignal->Clone("hDiff");
    hDiff -> Add(HBackground, -N);
    double h = hDiff -> Integral(hDiff->FindFixBin(-phiS), hDiff->FindFixBin(phiS), "");
    double Ratio = h/(indicator - 1);
    hDiff -> Sumw2();
    hDiff -> Draw("HIST SAME");
    //hDiff -> SetLineColor(kBlue);
    hDiff -> SetFillColor(kYellow - 7);
    hDiff -> SetFillStyle(1001);
    hDiff -> SetMarkerSize(0.7);
    hDiff -> SetTitle(Form("Average Multiplicity Density is  %5.2f", Ratio/2));
    //hDiff -> SetMarkerColor(kRed);
    hDiff -> Draw("e1psame");
    hDiff -> SetLabelSize(0.05, "xy");
    hDiff -> GetXaxis() -> SetRange(hDiff -> FindFixBin(-phiS+dx_2),hDiff -> FindFixBin(phiS));
    TLegend *legend2 = new TLegend(0.11, 0.8, 0.469, 0.9);
    legend2 -> AddEntry(hDiff, Form("Normalized in range of %1.2f < |#Delta#phi| < %1.2f", phiS, phiB), "f");
    //legend -> SetHeader("Red Region:");
    gStyle -> SetLegendTextSize(0.055);
    legend2 -> Draw("same");

    gStyle->SetEndErrorSize(2);
    gStyle->SetErrorX(0.5);
    //gPad ->SetLogy(1);

    


    TLine *l = new TLine(-2,0,2,0);
	  l ->Draw("same"); 
    l -> SetLineColor(kGreen);
    cout << "Ratio is " << Ratio << endl;
    cout << "Well done." << endl;
}