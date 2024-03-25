#include "./headers/Analysis.h"
void MBD_Analysis() {}

void load(std::string method = "1")
{
    TFile *f = TFile::Open("./rootData/metrics.root");
    if (!f || f->IsZombie()) {
        std::cerr << "Error opening file" << std::endl;
        return;
    }
    TTree *MetricTree;
    f->GetObject("MetricTree", MetricTree);
    if (!MetricTree) {
        std::cerr << "Tree not found" << std::endl;
        return;
    }
    
    TCanvas* globalCanvas = new TCanvas("globalCanvas", "Canvas with 2 Square Pads", 1200, 600); // Width, Height
    
    TH1F* h1 = nullptr;
    if (globalCanvas == nullptr) {
        globalCanvas = new TCanvas("globalCanvas", "Analysis Results", 800, 600);
    } else {
        globalCanvas->Clear();
    }
    if (h1 != nullptr) {
        delete h1;
        h1 = nullptr;
    }

    Int_t Event, Npart, NpartFromSource;
    Float_t foundZ_vtx, trueXfromSource, trueYfromSource, trueZfromSource;
    Float_t centralityBimp, centralityImpactparam, centralityMbdquantity, centralityMbd;
    MetricTree -> SetBranchAddress("Event", &Event);
    MetricTree -> SetBranchAddress("foundZ_vtx", &foundZ_vtx);
    MetricTree -> SetBranchAddress("Npart", &Npart);
    MetricTree -> SetBranchAddress("trueXfromSource", &trueXfromSource);
    MetricTree -> SetBranchAddress("trueYfromSource", &trueYfromSource);
    MetricTree -> SetBranchAddress("trueZfromSource", &trueZfromSource);
    MetricTree -> SetBranchAddress("NpartFromSource", &NpartFromSource);
    MetricTree -> SetBranchAddress("centralityBimp", &centralityBimp);
    MetricTree -> SetBranchAddress("centralityImpactparam", &centralityImpactparam);
    MetricTree -> SetBranchAddress("centralityMbdquantity", &centralityMbdquantity);
    MetricTree -> SetBranchAddress("centralityMbd", &centralityMbd);

    TGraph *g1 = new TGraph();
    for (Long64_t i = 0; i < MetricTree->GetEntries(); ++i) {
        MetricTree->GetEntry(i);
        // std::cout << Event << std::endl;
        if (centralityMbd >= MBD_lower && centralityMbd <= MBD_upper && trueZfromSource >= -25. && trueZfromSource <= -15.) {
            g1->SetPoint(g1->GetN(), trueZfromSource*10, (foundZ_vtx - trueZfromSource)*10);
        }
        
    }
    if (method == "project") {
        // double ymin = g1->GetHistogram()->GetMinimum();
        // double ymax = g1->GetHistogram()->GetMaximum();
        int nbins = 1000;
        h1 = new TH1F("ProjectionY", "Histogram of vtx_z_resolution vs vtx_z 1D projection;Found_z - True_z [mm];Entries", nbins, -60, +60);
        int npoints = g1->GetN();
        double x, y;
        for (int i = 0; i < npoints; i++) {
            g1->GetPoint(i, x, y);
            h1->Fill(y);
        }
        h1 -> Draw();
        h1 -> GetXaxis() -> CenterTitle(true);
        h1 -> GetYaxis() -> CenterTitle(true);
        h1 -> GetXaxis() -> SetRangeUser(-60, 60); // Setting x range;
        h1 -> GetYaxis() -> SetRangeUser(1, 100.);
        h1 -> SetFillColor(kYellow - 7);
        h1 -> SetLineWidth(1);
        h1 -> SetFillStyle(1001);
        static TLine *l1 = new TLine();
        l1 -> SetLineColor(kRed);
        l1 -> SetLineStyle(0);
        l1 -> SetLineWidth(2);
        l1 -> SetX1(0);   l1 -> SetY1(0);
        l1 -> SetX2(0);   l1 -> SetY2(h1->GetMaximum());
        l1 -> Draw("same");
        gPad -> SetLogy();
        TLegend *legend = new TLegend(0.15,0.8,0.35,0.9);
        // legend->SetHeader("Legend","C"); // option "C" allows to center the header
        legend->AddEntry(g1,Form("Centrality %2.0f-%2.0f", MBD_lower, MBD_upper),"f");
        legend->SetTextSize(0.04);
        legend->Draw("same");
        globalCanvas->Update();
    }
    else if (method == "slice") {
        // double ymin = g1->GetHistogram()->GetMinimum();
        // double ymax = g1->GetHistogram()->GetMaximum();
        TGraph *g2 = new TGraph();
        for (Long64_t i = 0; i < MetricTree->GetEntries(); ++i) {
            MetricTree->GetEntry(i);
            // std::cout << Event << std::endl;
            if (centralityMbd >= MBD_lower && centralityMbd <= MBD_upper && trueZfromSource >= -22.5 && trueZfromSource <= -21.5) {
                g2->SetPoint(g2->GetN(), trueZfromSource*10, (foundZ_vtx - trueZfromSource)*10);
            }   
        }
        int nbins = 200;
        h1 = new TH1F("Partial_ProjectionY", "Projection around True_z = -22cm;Found_z - True_z [mm];Entries", nbins, -60, +60);
        int npoints = g2->GetN();
        double x, y;
        for (int i = 0; i < npoints; i++) {
            g2->GetPoint(i, x, y);
            h1->Fill(y);
        }
        h1 -> Draw();
        h1 -> GetXaxis() -> CenterTitle(true);
        h1 -> GetYaxis() -> CenterTitle(true);
        h1 -> GetXaxis() -> SetRangeUser(-60, 60); // Setting x range;
        h1 -> GetYaxis() -> SetRangeUser(1, 30.);
        h1 -> SetFillColor(kYellow - 7);
        h1 -> SetLineWidth(1);
        h1 -> SetFillStyle(1001);
        static TLine *l1 = new TLine();
        l1 -> SetLineColor(kRed);
        l1 -> SetLineStyle(0);
        l1 -> SetLineWidth(2);
        l1 -> SetX1(0);   l1 -> SetY1(0);
        l1 -> SetX2(0);   l1 -> SetY2(h1->GetMaximum());
        l1 -> Draw("same");
        // gPad -> SetLogy();
        TLegend *legend = new TLegend(0.15,0.8,0.35,0.9);
        // legend->SetHeader("Legend","C"); // option "C" allows to center the header
        legend->AddEntry(g1,Form("Centrality %2.0f-%2.0f", MBD_lower, MBD_upper),"f");
        legend->SetTextSize(0.04);
        legend->Draw("same");
        globalCanvas->Update();
    }
    else if (method == "allXY") {
        globalCanvas->Clear();
        globalCanvas->Divide(2, 1); // Divide the canvas into 2 pads: 2 columns, 1 row
        TGraph *g2 = new TGraph();
        TGraph *g3 = new TGraph();
        for (Long64_t i = 0; i < MetricTree->GetEntries(); ++i) {
            MetricTree->GetEntry(i);
            // std::cout << Event << std::endl;
            if (foundZ_vtx - trueZfromSource > 0) {
                g2->SetPoint(g2->GetN(), trueXfromSource*10, trueYfromSource*10);
            } else {
                g3->SetPoint(g2->GetN(), trueXfromSource*10, trueYfromSource*10);
            }
        }

        globalCanvas->cd(1);
        gPad->SetMargin(0.1, 0.1, 0.1, 0.1);
        gPad->SetGrid(1, 1);
        if (g2->GetN() > 0) { // Check if there are points to draw
            g2->Draw("AP"); 
            g2->GetXaxis()->SetLimits(-0.5, 0.5); // X axis limits
            g2->SetMinimum(-0.5); // Y axis minimum
            g2->SetMaximum(0.5);  // Y axis maximum
            g2 -> SetMarkerStyle(29);
            g2 -> SetMarkerSize(0.9);
            g2 -> SetMarkerColor(kBlue - 7);
            g2 -> SetLineWidth(3);
            g2 -> SetLineColor(kWhite);
            gStyle -> SetTitleW(0.7);  //per cent of the pad width
            gStyle -> SetTitleH(0.08); //per cent of the pad height
            g2 -> SetTitle("For resolution > 0");
            g2 -> GetXaxis() -> SetTitle("True X Vertex [mm]");
            g2 -> GetXaxis() -> SetTitleSize(0.05);
            g2 -> GetXaxis() -> SetLabelSize(0.025);
            g2 -> GetXaxis() -> CenterTitle(true);
            g2 -> GetYaxis() -> SetTitle("True Y Vertex [mm]");
            g2 -> GetYaxis() -> SetTitleSize(0.05);
            g2 -> GetYaxis() -> SetLabelSize(0.025);
            g2 -> GetYaxis() -> CenterTitle(true);
            g2 -> GetYaxis() -> SetTitleOffset(0.8); 
            g2 -> GetXaxis() -> SetTitleOffset(0.8); 
        } else {
            std::cerr << "No points to display." << std::endl;
        }
        globalCanvas->cd(2); // Change to the second pad
        gPad->SetMargin(0.1, 0.1, 0.1, 0.1);
        gPad->SetGrid(1, 1);
        if (g3->GetN() > 0) {
            g3->Draw("AP"); // Draw the graph with axes and points
            g3->GetXaxis()->SetLimits(-0.5, 0.5); // X axis limits
            g3->SetMinimum(-0.5); // Y axis minimum 
            g3->SetMaximum(0.5);  // Y axis maximum
            g3->GetXaxis()->SetLimits(-0.5, 0.5); // X axis limits
            g3->SetMinimum(-0.5); // Y axis minimum
            g3->SetMaximum(0.5);  // Y axis maximum
            g3 -> SetMarkerStyle(29);
            g3 -> SetMarkerSize(0.9);
            g3 -> SetMarkerColor(kRed - 7);
            g3 -> SetLineWidth(3);
            g3 -> SetLineColor(kWhite);
            gStyle -> SetTitleW(0.7);  //per cent of the pad width
            gStyle -> SetTitleH(0.08); //per cent of the pad height
            g3 -> SetTitle("For resolution < 0");
            g3 -> GetXaxis() -> SetTitle("True X Vertex [mm]");
            g3 -> GetXaxis() -> SetTitleSize(0.05);
            g3 -> GetXaxis() -> SetLabelSize(0.025);
            g3 -> GetXaxis() -> CenterTitle(true);
            g3 -> GetYaxis() -> SetTitle("True Y Vertex [mm]");
            g3 -> GetYaxis() -> SetTitleSize(0.05);
            g3 -> GetYaxis() -> SetLabelSize(0.025);
            g3 -> GetYaxis() -> CenterTitle(true);
            g3 -> GetYaxis() -> SetTitleOffset(0.8); 
            g3 -> GetXaxis() -> SetTitleOffset(0.8); 
        } else {
            std::cerr << "No points in the second graph to display." << std::endl;
        }

        globalCanvas->Update();
    }
    else if (method == "sliceXY") {
        globalCanvas->Clear();
        globalCanvas->Divide(2, 1); // Divide the canvas into 2 pads: 2 columns, 1 row
        TGraph *g2 = new TGraph();
        TGraph *g3 = new TGraph();
        for (Long64_t i = 0; i < MetricTree->GetEntries(); ++i) {
            MetricTree->GetEntry(i);
            if (trueZfromSource >= -22.5 && trueZfromSource <= -21.5) {
                if (foundZ_vtx - trueZfromSource > 0) {
                    g2->SetPoint(g2->GetN(), trueXfromSource*10, trueYfromSource*10);
                } else {
                    g3->SetPoint(g2->GetN(), trueXfromSource*10, trueYfromSource*10);
                }
            }
            
        }

        globalCanvas->cd(1);
        gPad->SetMargin(0.1, 0.1, 0.1, 0.1);
        gPad->SetGrid(1, 1);
        if (g2->GetN() > 0) { // Check if there are points to draw
            g2->Draw("AP"); 
            g2->GetXaxis()->SetLimits(-0.5, 0.5); // X axis limits
            g2->SetMinimum(-0.5); // Y axis minimum
            g2->SetMaximum(0.5);  // Y axis maximum
            g2 -> SetMarkerStyle(29);
            g2 -> SetMarkerSize(0.9);
            g2 -> SetMarkerColor(kBlue - 7);
            g2 -> SetLineWidth(3);
            g2 -> SetLineColor(kWhite);
            gStyle -> SetTitleW(0.7);  //per cent of the pad width
            gStyle -> SetTitleH(0.08); //per cent of the pad height
            g2 -> SetTitle("For resolution > 0");
            g2 -> GetXaxis() -> SetTitle("True X Vertex [mm]");
            g2 -> GetXaxis() -> SetTitleSize(0.05);
            g2 -> GetXaxis() -> SetLabelSize(0.025);
            g2 -> GetXaxis() -> CenterTitle(true);
            g2 -> GetYaxis() -> SetTitle("True Y Vertex [mm]");
            g2 -> GetYaxis() -> SetTitleSize(0.05);
            g2 -> GetYaxis() -> SetLabelSize(0.025);
            g2 -> GetYaxis() -> CenterTitle(true);
            g2 -> GetYaxis() -> SetTitleOffset(0.8); 
            g2 -> GetXaxis() -> SetTitleOffset(0.8); 
        } else {
            std::cerr << "No points to display." << std::endl;
        }
        globalCanvas->cd(2); // Change to the second pad
        gPad->SetMargin(0.1, 0.1, 0.1, 0.1);
        gPad->SetGrid(1, 1);
        if (g3->GetN() > 0) {
            g3->Draw("AP"); // Draw the graph with axes and points
            g3->GetXaxis()->SetLimits(-0.5, 0.5); // X axis limits
            g3->SetMinimum(-0.5); // Y axis minimum 
            g3->SetMaximum(0.5);  // Y axis maximum
            g3->GetXaxis()->SetLimits(-0.5, 0.5); // X axis limits
            g3->SetMinimum(-0.5); // Y axis minimum
            g3->SetMaximum(0.5);  // Y axis maximum
            g3 -> SetMarkerStyle(29);
            g3 -> SetMarkerSize(0.9);
            g3 -> SetMarkerColor(kRed - 7);
            g3 -> SetLineWidth(3);
            g3 -> SetLineColor(kWhite);
            gStyle -> SetTitleW(0.7);  //per cent of the pad width
            gStyle -> SetTitleH(0.08); //per cent of the pad height
            g3 -> SetTitle("For resolution < 0");
            g3 -> GetXaxis() -> SetTitle("True X Vertex [mm]");
            g3 -> GetXaxis() -> SetTitleSize(0.05);
            g3 -> GetXaxis() -> SetLabelSize(0.025);
            g3 -> GetXaxis() -> CenterTitle(true);
            g3 -> GetYaxis() -> SetTitle("True Y Vertex [mm]");
            g3 -> GetYaxis() -> SetTitleSize(0.05);
            g3 -> GetYaxis() -> SetLabelSize(0.025);
            g3 -> GetYaxis() -> CenterTitle(true);
            g3 -> GetYaxis() -> SetTitleOffset(0.8); 
            g3 -> GetXaxis() -> SetTitleOffset(0.8); 
        } else {
            std::cerr << "No points in the second graph to display." << std::endl;
        }

        globalCanvas->Update();
    }
    else {
        globalCanvas->Clear();
        globalCanvas->cd();

        if (g1->GetN() > 0) { // Check if there are points to draw
            g1->Draw("AP"); // Specify the drawing options as needed
            globalCanvas->Update(); // Make sure to update the canvas to reflect the drawn graph
        } else {
            std::cerr << "No points in the graph to display." << std::endl;
        }
        g1 -> SetMarkerStyle(29);
        g1 -> SetMarkerSize(0.5);
        g1 -> SetMarkerColor(kBlue - 7);
        g1 -> SetLineWidth(3);
        g1 -> SetLineColor(kWhite);
        gStyle -> SetTitleW(0.7);  //per cent of the pad width
        gStyle -> SetTitleH(0.08); //per cent of the pad height
        g1 -> GetXaxis() -> SetTitle("True Z Vertex [mm]");
        g1 -> GetXaxis() -> SetTitleSize(0.06);
        g1 -> GetXaxis() -> SetLabelSize(0.04);
        g1 -> GetXaxis() -> CenterTitle(true);
        g1 -> GetYaxis() -> SetTitle("Found_z - True_z [mm]");
        g1 -> GetYaxis() -> SetTitleSize(0.06);
        g1 -> GetYaxis() -> SetLabelSize(0.025);
        g1 -> GetYaxis() -> CenterTitle(true);
        g1 -> SetMinimum(-100); // Setting y range;
        g1 -> SetMaximum(100);  // Setting y range;
        g1 -> GetYaxis() -> SetTitleOffset(0.8); 
        g1 -> GetXaxis() -> SetTitleOffset(0.8); 
        // g1 -> GetXaxis() -> SetLimits(0, 6000); // Setting x range;
        g1 -> Draw("AP SAME");
        gPad->SetGrid(5, 2); gPad->Update();

        TLegend *legend = new TLegend(0.35,0.8,0.65,0.9);
        legend->AddEntry(g1,Form("Centrality %2.0f-%2.0f", MBD_lower, MBD_upper),"f");
        legend->SetTextSize(0.04);
        legend->Draw("same");
        globalCanvas->Update();
    }
    // delete g1;
}