#include <TFile.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TLine.h>
#include <iostream>

void dPhi_Fitter() {
    // Open the .root file
    TFile *inputFile = new TFile("../External/zFindingPlots/hDiff_output.root", "READ");
    if (!inputFile || inputFile->IsZombie()) {
        std::cerr << "Error: Cannot open .root file!" << std::endl;
        return;
    }

    // Retrieve the histogram
    TH1D *hDiff = (TH1D*)inputFile->Get("Background Subtracted Signal");
    if (!hDiff) {
        std::cerr << "Error: Cannot find the target histogram!" << std::endl;
        inputFile->Close();
        return;
    }

    // Create a canvas for plotting
    TCanvas *can = new TCanvas("can", "Background Subtracted Signal", 1920, 1056);
    
    // Set some drawing options (optional)
    hDiff->SetLineColor(kBlue);
    hDiff->SetLineWidth(2);

    // hDiff->SetMarkerStyle(0);  // Disable markers
    // hDiff->SetMarkerSize(0);   // Ensure no marker size

    hDiff->Draw("HIST");

    // Add additional processing or fitting here if needed
    // Example: fitting a Gaussian function to the histogram
    // hDiff->Fit("gaus");
    
    // Optional: Draw a line at y=0 for reference
    TLine *line = new TLine(hDiff->GetXaxis()->GetXmin(), 0, hDiff->GetXaxis()->GetXmax(), 0);
    line->SetLineColor(kRed);
    line->Draw("same");

    // Optional: Add a legend
    TLegend *legend = new TLegend(0.7, 0.8, 0.9, 0.9);
    legend->AddEntry(hDiff, "Background Subtracted Signal", "l");
    legend->Draw();

    // Show the plot
    can->Update();

    // Close the input file
    // inputFile->Close();
    
}
