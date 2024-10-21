#include <TFile.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TLine.h>
#include <TStopwatch.h>

#include <iostream>

double ChiTwo (const TH1D * const hDiff, const double &p1, const double &p2, const double &low_range) {
    int nbins = hDiff->GetNbinsX();
    double x_value, y_value;
    double ChiSum = 0.;
    for (int i = 1; i < nbins; i++) {
        x_value = hDiff->GetBinCenter(i);
        y_value = hDiff->GetBinContent(i);
        // if (x_value >= (hDiff->GetXaxis()->GetXmin())/5.)   break;
        if (x_value >= low_range)   break;

        double calculated = p1*p1*cos(2*x_value)+p2 - y_value;
        // double bin_error = hDiff->GetBinError(i);
        double bin_error = 1.;
        ChiSum += (calculated*calculated)/(bin_error*bin_error);
        // ChiSum += ((p1*p1*cos(2*x_value)+p2 - y_value)*(p1*p1*cos(2*x_value)+p2 - y_value))/((hDiff->GetBinError(i))*(hDiff->GetBinError(i)));
    }
    return ChiSum;
}

// Actually this is a Least-Square test:
double ChiTwo_ver2 (const TH1D * const hDiff, const double &p1, const double &p2, const double &low_range, const double &high_range) {
    int nbins = hDiff->GetNbinsX();
    double x_value, y_value;
    double ChiSum = 0.;
    for (int i = 1; i < nbins; i++) {
        x_value = hDiff->GetBinCenter(i);
        y_value = hDiff->GetBinContent(i);
        // if (x_value <= (hDiff->GetXaxis()->GetXmin())/5. || x_value >= (hDiff->GetXaxis()->GetXmax())/5.) {
        if (x_value <= low_range || x_value >= high_range) {
            double calculated = p1*p1*cos(2*x_value)+p2 - y_value;
            ChiSum += calculated*calculated;
        }
    }
    return ChiSum;
}

// "Standard" Chi-Squared test:
double ChiTwo_ver3 (const TH1D * const hDiff, const double &p1, const double &p2, const double &low_range, const double &high_range) {
    int nbins = hDiff->GetNbinsX();
    double x_value, y_value;
    double ChiSum = 0.;
    for (int i = 1; i < nbins; i++) {
        x_value = hDiff->GetBinCenter(i);
        y_value = hDiff->GetBinContent(i);
        // if (x_value <= (hDiff->GetXaxis()->GetXmin())/5. || x_value >= (hDiff->GetXaxis()->GetXmax())/5.) {
        if (x_value <= low_range || x_value >= high_range) {
            double calculated = p1*p1*cos(2*x_value)+p2 - y_value;
            ChiSum += (calculated*calculated)/y_value;
        }
    }
    return ChiSum;
}

void fit_the_hist (TH1D *hDiff, TCanvas *can) {
    double left_subrange_min = hDiff->GetXaxis()->GetXmin();
    // double left_subrange_max = (hDiff->GetXaxis()->GetXmin())/5.;
    double left_subrange_max = -0.2;
    TF1 *bottom_noise = new TF1("bottom_noise", "[0]*[0]*cos(2*x) + [1]", left_subrange_min, left_subrange_max);

    double right_subrange_max = hDiff->GetXaxis()->GetXmax();
    // double right_subrange_min = (hDiff->GetXaxis()->GetXmax())/5.;
    double right_subrange_min = +0.2;
    TF1 *bottom_noise_2 = new TF1("bottom_noise", "[0]*[0]*cos(2*x) + [1]", right_subrange_min, right_subrange_max);

    double rst;
    double min_rst = std::numeric_limits<double>::max();
    double par1 = 800, par2 = -40e3;
    for (double p_1 = 800.0; p_1 < 1500.0; p_1++) {
        for (double p_2 = -400000.0; p_2 < 400000.0; p_2 += 200.0) {
            rst = ChiTwo_ver2(hDiff, p_1, p_2, left_subrange_max, right_subrange_min);
            if (min_rst > rst) {
                min_rst = rst;
                par1 = p_1;
                par2 = p_2;
            }
        }
    }
    std::cout << min_rst << ", " << par1 << ", " << par2 << std::endl;

    bottom_noise -> SetLineColor(kTeal - 7);
    bottom_noise -> SetLineWidth(10);
    bottom_noise -> SetParameter(0, par1);
    bottom_noise -> SetParameter(1, par2);
    bottom_noise -> Draw("same");
    bottom_noise_2 -> SetLineColor(kTeal - 7);
    bottom_noise_2 -> SetLineWidth(10);
    bottom_noise_2 -> SetParameter(0, par1);
    bottom_noise_2 -> SetParameter(1, par2);
    bottom_noise_2 -> Draw("same");

    can -> cd(2);
    int nBins = hDiff->GetNbinsX();
    TH1D *hSubtracted = new TH1D("hSubtracted", "Subtracted Histogram", nBins, hDiff->GetXaxis()->GetXmin(), hDiff->GetXaxis()->GetXmax());
    TF1 *bottomNoise = new TF1("bottomNoise", "[0]*[0]*cos(2*x) + [1]", hDiff->GetXaxis()->GetXmin(), hDiff->GetXaxis()->GetXmax());
    bottomNoise->SetParameters(par1, par2);
    // Loop over each bin
    for (int i = 1; i <= nBins; ++i) {
        double binCenter = hDiff->GetBinCenter(i);
        double fitValue = bottomNoise->Eval(binCenter);
        double newBinContent = hDiff->GetBinContent(i) - fitValue;
        hSubtracted->SetBinContent(i, newBinContent);
    }

    for (int i = 1; i <= nBins; ++i) {
        std::cout << hSubtracted->GetBinCenter(i) << ",  " << hSubtracted->GetBinContent(i) << std::endl;
    }

    hSubtracted->Draw();
    hSubtracted->SetFillStyle(3003);
    hSubtracted->SetFillColor(kBlue - 7);
    hSubtracted->GetXaxis()->SetTitle("dPhi value");
    hSubtracted->GetXaxis()->CenterTitle(true);
    hSubtracted->GetXaxis()->SetTitleSize(0.05);
    hSubtracted->GetXaxis()->SetTitleOffset(.8);
    hSubtracted->GetYaxis()->SetTitle("# of count");
    hSubtracted->GetYaxis()->CenterTitle(true);
    hSubtracted->GetYaxis()->SetTitleSize(0.05);
    hSubtracted->GetYaxis()->SetTitleOffset(.8);

    TLine *line_horizontal = new TLine(hSubtracted->GetXaxis()->GetXmin(), 0, hSubtracted->GetXaxis()->GetXmax(), 0);
    line_horizontal->SetLineColor(kRed);
    line_horizontal->Draw("same");

    // Show the plot
    can->Update();
}

void fit_the_hist_ver2 (TH1D *hDiff, TCanvas *can) {
    for (int i = 1; i<hDiff->GetNbinsX(); i++)  hDiff->SetBinError(i, 1.);
    double par1 = 800, par2 = -95200;
    double left_subrange_min = hDiff->GetXaxis()->GetXmin();
    double left_subrange_max = (hDiff->GetXaxis()->GetXmin())/5.;
    TF1 *bottom_noise = new TF1("bottom_noise", "[0]*[0]*cos(2*x) + [1]", left_subrange_min, left_subrange_max);

    double right_subrange_max = hDiff->GetXaxis()->GetXmax();
    double right_subrange_min = (hDiff->GetXaxis()->GetXmax())/5.;
    TF1 *bottom_noise_2 = new TF1("bottom_noise", "[0]*[0]*cos(2*x) + [1]", right_subrange_min, right_subrange_max);
    bottom_noise -> SetLineColor(7);
    bottom_noise -> SetLineWidth(10);
    bottom_noise -> SetParameter(0, par1);
    bottom_noise -> SetParameter(1, par2);
    bottom_noise -> Draw("same");
    bottom_noise_2 -> SetLineColor(8);
    bottom_noise_2 -> SetLineWidth(10);
    bottom_noise_2 -> SetParameter(0, par1);
    bottom_noise_2 -> SetParameter(1, par2);
    hDiff->Fit(bottom_noise, "RMS");
    hDiff->Fit(bottom_noise_2, "RMS+");
    bottom_noise_2 -> Draw("same");
    // Show the plot
    can->Update();
}

void dPhi_Fitter() {
    TStopwatch timer;   timer.Start();

    TFile *inputFile = new TFile("../External/zFindingPlots/hDiff_output.root", "READ");
    if (!inputFile || inputFile->IsZombie()) {
        std::cerr << "Error: Cannot open .root file!" << std::endl;
        return;
    }

    TH1D *hDiff = (TH1D*)inputFile->Get("Background Subtracted Signal");
    if (!hDiff) {
        std::cerr << "Error: Cannot find the target histogram!" << std::endl;
        inputFile->Close();
        return;
    }

    TCanvas *can = new TCanvas("can", "Background Subtracted Signal", 1920, 1056);
    can -> Divide(1, 2);
    can -> cd(1);
    hDiff->SetLineColor(kBlue);
    hDiff->SetLineWidth(2);
    // hDiff->SetMarkerStyle(0);  // Disable markers
    // hDiff->SetMarkerSize(10);   // Ensure no marker size
    hDiff->Draw("hist");

    // Add additional processing or fitting here if needed
    // Example: fitting a Gaussian function to the histogram
    // hDiff->Fit("gaus");
    
    // Cosmetics:
    TLine *line_horizontal = new TLine(hDiff->GetXaxis()->GetXmin(), 0, hDiff->GetXaxis()->GetXmax(), 0);
    line_horizontal->SetLineColor(kRed);
    line_horizontal->Draw("same");
    // TLine *line_left = new TLine((hDiff->GetXaxis()->GetXmin())/5., hDiff->GetMinimum(), (hDiff->GetXaxis()->GetXmin())/5., hDiff->GetMaximum());
    // line_left->SetLineColor(kRed);
    // line_left->Draw("same");
    TLegend *legend = new TLegend(0.7, 0.8, 0.9, 0.9);
    legend->AddEntry(hDiff, "Background Subtracted Signal", "l");
    legend->Draw();

    fit_the_hist(hDiff, can);
    // fit_the_hist_ver2(hDiff, can);

    // Stop the stopwatch and print the runtime:
    timer.Stop();   timer.Print();
}
