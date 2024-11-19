#include <TFile.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TLine.h>
#include <TStopwatch.h>

#include <iostream>
#include "./headers/globalDefinitions.h"
#include "./headers/histogramPlotting.h"

double rootFit_ChiTwo (const TH1D * const hDiff, const double &p1, const double &p2, const double &low_range, const double &high_range) {
    int nbins = hDiff->GetNbinsX();
    double x_value, y_value;
    double ChiSum = 0.;
    for (int i = 1; i < nbins; i++) {
        x_value = hDiff->GetBinCenter(i);
        y_value = hDiff->GetBinContent(i);
        if (x_value <= low_range || x_value >= high_range) {
            double calculated = p1*p1*cos(2*x_value)+p2 - y_value;
            double bin_error = hDiff->GetBinError(i);
            ChiSum += (calculated*calculated)/(bin_error*bin_error);
        }
    }
    return ChiSum;
}

double LeastSquare (const TH1D * const hDiff, const double &p1, const double &p2, const double &low_range, const double &high_range) {
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
double ChiSquared (const TH1D * const hDiff, const double &p1, const double &p2, const double &low_range, const double &high_range) {
    int nbins = hDiff->GetNbinsX();
    double x_value, y_value;
    double ChiSum = 0.;
    for (int i = 1; i < nbins; i++) {
        x_value = hDiff->GetBinCenter(i);
        y_value = hDiff->GetBinContent(i);
        // if (x_value <= (hDiff->GetXaxis()->GetXmin())/5. || x_value >= (hDiff->GetXaxis()->GetXmax())/5.) {
        if (x_value <= low_range || x_value >= high_range) {
            double expected = p1*p1*cos(2*x_value)+p2;
            double calculated = expected - y_value;
            ChiSum += (calculated*calculated)/expected;
        }
    }
    return ChiSum;
}

double TotalLeastSquare (const TH1D * const hDiff, const double &p1, const double &p2, const double &low_range, const double &high_range) {
    int nbins = hDiff->GetNbinsX();
    double x_value, y_value;
    double ChiSum = 0.;
    for (int i = 1; i < nbins; i++) {
        x_value = hDiff->GetBinCenter(i);
        y_value = hDiff->GetBinContent(i);
        // if (x_value <= (hDiff->GetXaxis()->GetXmin())/5. || x_value >= (hDiff->GetXaxis()->GetXmax())/5.) {
        if (x_value <= low_range || x_value >= high_range) {
            double numerator   = p1*p1*cos(2*x_value)+p2 - y_value;
            double dydx        = -2*p1*p1*sin(2*x_value);
            double denominator = 1 + dydx*dydx;
            ChiSum += (numerator*numerator)/denominator;
        }
    }
    return ChiSum;
}

void fit_the_hist (
    // const std::string &rootFilePath,
    const std::vector<std::string> &rootFilePaths,
    const TString &option,
    const int &target
) {
    std::string rootFilePath = rootFilePaths[target];
    printBlue(rootFilePath);
    TFile *inputFile = new TFile(rootFilePath.c_str(), "READ");
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

    std::string hDiff_title = hDiff->GetTitle();
    std::vector<std::string> splitted = splitString(hDiff_title, ' ');
    int numberOfEvents = std::stoi(splitted[3]);
    double left_subrange_min = hDiff->GetXaxis()->GetXmin();
    // double left_subrange_max = (hDiff->GetXaxis()->GetXmin())/5.;
    double left_subrange_max = -abs_fit_range;
    TF1 *bottom_noise = new TF1("bottom_noise", "[0]*[0]*cos(2*x) + [1]", left_subrange_min, left_subrange_max);

    double right_subrange_max = hDiff->GetXaxis()->GetXmax();
    // double right_subrange_min = (hDiff->GetXaxis()->GetXmax())/5.;
    double right_subrange_min = abs_fit_range;
    TF1 *bottom_noise_2 = new TF1("bottom_noise", "[0]*[0]*cos(2*x) + [1]", right_subrange_min, right_subrange_max);

    TCanvas *can = new TCanvas("can", "Background Subtracted Signal", 3024, 1964);
    can -> Divide(1, 2);
    can -> cd(1);
    hDiff->SetLineColor(kBlue);
    hDiff->SetLineWidth(2);
    // hDiff->SetMarkerStyle(0);  // Disable markers
    // hDiff->SetMarkerSize(10);   // Ensure no marker size
    hDiff->Draw("hist");
    
    // Cosmetics:
    TLine *line_horizontal = new TLine(hDiff->GetXaxis()->GetXmin(), 0, hDiff->GetXaxis()->GetXmax(), 0);
    line_horizontal->SetLineColor(kRed);
    line_horizontal->Draw("same");
    // TLine *line_left = new TLine((hDiff->GetXaxis()->GetXmin())/5., hDiff->GetMinimum(), (hDiff->GetXaxis()->GetXmin())/5., hDiff->GetMaximum());
    // line_left->SetLineColor(kRed);
    // line_left->Draw("same");
    TLegend *legend = new TLegend(0.13, 0.75, 0.29, 0.9);
    legend->AddEntry(hDiff, Form("Fit with %s", option.Data()), "l");
    legend->AddEntry(bottom_noise, Form("Fit range ~ [%.2f, %.2f]",left_subrange_max, right_subrange_min), "l");
    legend->Draw();

    double rst;
    double min_rst = std::numeric_limits<double>::max();
    double par1 = 100., par2 = -30e4;
    double start_p1 = par1, start_p2 = par2;

    std::function<double(const TH1D *const, const double &, const double &, const double &, const double &)> myFitFunction;
    if (option == "TLS")    myFitFunction = TotalLeastSquare;
    if (option == "LS")     myFitFunction = LeastSquare;
    if (option == "root")   myFitFunction = rootFit_ChiTwo;
    if (option == "Chi2")   myFitFunction = ChiSquared;

    for (double p_1 = start_p1; p_1 < 1800.0; p_1++) {
            for (double p_2 = start_p2; p_2 < 80.e4; p_2 += 100.0) {
            rst = myFitFunction(hDiff, p_1, p_2, left_subrange_max, right_subrange_min);
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

    std::vector<std::string> cutInfos = splitString(rootFilePath, '_');
    std::string cenRange1 = cutInfos[3];
    std::string cenRange2 = cutInfos[4].substr(0, cutInfos[4].length()-5);
    TH1D *hFiltered = new TH1D("hFiltered", Form("Filtered Histogram of %s < Centrality < %s;dPhi value;# of count", cenRange1.c_str(), cenRange2.c_str()), nBins, hDiff->GetXaxis()->GetXmin(), hDiff->GetXaxis()->GetXmax());
    TF1 *bottomNoise = new TF1("bottomNoise", "[0]*[0]*cos(2*x) + [1]", hDiff->GetXaxis()->GetXmin(), hDiff->GetXaxis()->GetXmax());
    bottomNoise->SetParameters(par1, par2);
    // Loop over each bin
    for (int i = 1; i <= nBins; ++i) {
        double binCenter = hDiff->GetBinCenter(i);
        double fitValue = bottomNoise->Eval(binCenter);
        double newBinContent = hDiff->GetBinContent(i) - fitValue;
        hFiltered->SetBinContent(i, newBinContent);
    }

    // for (int i = 1; i <= nBins; ++i) {
    //     std::cout << hFiltered->GetBinCenter(i) << ",  " << hFiltered->GetBinContent(i) << std::endl;
    // }

    hFiltered->Draw();
    hFiltered->SetFillStyle(3003);
    hFiltered->SetFillColor(kBlue - 7);
    hFiltered->GetXaxis()->CenterTitle(true);
    hFiltered->GetXaxis()->SetTitleSize(0.05);
    hFiltered->GetXaxis()->SetTitleOffset(.8);
    hFiltered->GetYaxis()->CenterTitle(true);
    hFiltered->GetYaxis()->SetTitleSize(0.05);
    hFiltered->GetYaxis()->SetTitleOffset(.8);

    std::cout << left_subrange_min << ", " << left_subrange_max << ", " << right_subrange_min << ", " << right_subrange_max << std::endl;
    double Signal = hFiltered->Integral(hFiltered->FindFixBin(left_subrange_max),hFiltered->FindFixBin(right_subrange_min));
    double Background = hFiltered->Integral(hFiltered->FindFixBin(left_subrange_min), hFiltered->FindFixBin(left_subrange_max))
                      + hFiltered->Integral(hFiltered->FindFixBin(right_subrange_min), hFiltered->FindFixBin(right_subrange_max));
    double Multiplicity = Signal/numberOfEvents;
    double SNR = Signal/Background;
    std::cout << Signal << ", " << Background << ", " << SNR << ", " << Multiplicity << ", " << numberOfEvents << std::endl;

    TLine *line_horizon = new TLine(hFiltered->GetXaxis()->GetXmin(), 0, hFiltered->GetXaxis()->GetXmax(), 0);
    line_horizon->SetLineColor(kRed);
    line_horizon->Draw("same");
    TLegend *legend2 = new TLegend(0.13, 0.75, 0.29, 0.9);
    legend2->AddEntry(hFiltered, Form("Multiplicity Density: %.2f", Multiplicity), "f");
    legend2->AddEntry(hFiltered, Form("Signal/Background Ratio: %.2f",SNR), "l");
    legend2->Draw("same");

    // Show the plot
    can->Update();
}

void fit_hists_all (
    const std::vector<std::string> &rootFilePaths,
    const TString &option
) {
    int i = 0;
    for (const std::string &rootFilePath : rootFilePaths) {
        std::vector<std::string> cutInfos = splitString(rootFilePath, '_');
        std::string cenRange1 = cutInfos[8];
        std::string cenRange2 = cutInfos[9].substr(0, cutInfos[9].length()-5);
        TFile *inputFile = new TFile(rootFilePath.c_str(), "READ");
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

        TCanvas *can = new TCanvas(Form("can%d", i), Form("Background Subtracted Signal%d", i), 3024, 1964);
        can -> Divide(1, 2);
        can -> cd(1);
        hDiff->SetLineColor(kBlue);
        hDiff->SetLineWidth(2);
        // hDiff->SetMarkerStyle(0);  // Disable markers
        // hDiff->SetMarkerSize(10);   // Ensure no marker size
        hDiff->Draw("hist");

        double left_subrange_min = hDiff->GetXaxis()->GetXmin();
        // double left_subrange_max = (hDiff->GetXaxis()->GetXmin())/5.;
        double left_subrange_max = -abs_fit_range;
        TF1 *bottom_noise = new TF1("bottom_noise", "[0]*[0]*cos(2*x) + [1]", left_subrange_min, left_subrange_max);

        double right_subrange_max = hDiff->GetXaxis()->GetXmax();
        // double right_subrange_min = (hDiff->GetXaxis()->GetXmax())/5.;
        double right_subrange_min = +abs_fit_range;
        TF1 *bottom_noise_2 = new TF1("bottom_noise", "[0]*[0]*cos(2*x) + [1]", right_subrange_min, right_subrange_max);
        
        // Cosmetics:
        TLine *line_horizontal = new TLine(hDiff->GetXaxis()->GetXmin(), 0, hDiff->GetXaxis()->GetXmax(), 0);
        line_horizontal->SetLineColor(kRed);
        line_horizontal->Draw("same");
        TLegend *legend = new TLegend(0.13, 0.75, 0.29, 0.9);
        legend->AddEntry(hDiff, Form("Fit with %s", option.Data()), "l");
        legend->AddEntry(bottom_noise, Form("Fit range ~ [%.2f, %.2f]",left_subrange_max, right_subrange_min), "l");
        legend->Draw();

        double rst;
        double min_rst = std::numeric_limits<double>::max();
        double par1 = 100., par2 = -40e4;
        double start_p1 = par1, start_p2 = par2;

        std::function<double(const TH1D *const, const double &, const double &, const double &, const double &)> myFitFunction;
        if (option == "TLS")    myFitFunction = TotalLeastSquare;
        if (option == "LS")     myFitFunction = LeastSquare;
        if (option == "root")   myFitFunction = rootFit_ChiTwo;
        if (option == "Chi2")   myFitFunction = ChiSquared;
        for (double p_1 = start_p1; p_1 < 1800.0; p_1++) {
            for (double p_2 = start_p2; p_2 < 40.e4; p_2 += 100.0) {
                rst = myFitFunction(hDiff, p_1, p_2, left_subrange_max, right_subrange_min);
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
        // TH1D *hFiltered = new TH1D("hFiltered", "Subtracted Histogram", nBins, hDiff->GetXaxis()->GetXmin(), hDiff->GetXaxis()->GetXmax());
        TH1D *hFiltered = new TH1D("hFiltered", Form("Filtered Histogram of %s < Centrality < %s;dPhi value;# of count", cenRange1.c_str(), cenRange2.c_str()), nBins, hDiff->GetXaxis()->GetXmin(), hDiff->GetXaxis()->GetXmax());
        TF1 *bottomNoise = new TF1("bottomNoise", "[0]*[0]*cos(2*x) + [1]", hDiff->GetXaxis()->GetXmin(), hDiff->GetXaxis()->GetXmax());
        bottomNoise->SetParameters(par1, par2);
        // Loop over each bin
        for (int i = 1; i <= nBins; ++i) {
            double binCenter = hDiff->GetBinCenter(i);
            double fitValue = bottomNoise->Eval(binCenter);
            double newBinContent = hDiff->GetBinContent(i) - fitValue;
            hFiltered->SetBinContent(i, newBinContent);
        }

        hFiltered->Draw();
        hFiltered->SetFillStyle(3003);
        hFiltered->SetFillColor(kBlue - 7);
        hFiltered->GetXaxis()->SetTitle("dPhi value");
        hFiltered->GetXaxis()->CenterTitle(true);
        hFiltered->GetXaxis()->SetTitleSize(0.05);
        hFiltered->GetXaxis()->SetTitleOffset(.8);
        hFiltered->GetYaxis()->SetTitle("# of count");
        hFiltered->GetYaxis()->CenterTitle(true);
        hFiltered->GetYaxis()->SetTitleSize(0.05);
        hFiltered->GetYaxis()->SetTitleOffset(.8);

        std::string hDiff_title = hDiff->GetTitle();
        std::vector<std::string> splitted = splitString(hDiff_title, ' ');
        int numberOfEvents = std::stoi(splitted[3]);
        std::cout << left_subrange_min << ", " << left_subrange_max << ", " << right_subrange_min << ", " << right_subrange_max << std::endl;
        double Signal = hFiltered->Integral(hFiltered->FindFixBin(left_subrange_max),hFiltered->FindFixBin(right_subrange_min));
        double Background = hFiltered->Integral(hFiltered->FindFixBin(left_subrange_min), hFiltered->FindFixBin(left_subrange_max))
                        + hFiltered->Integral(hFiltered->FindFixBin(right_subrange_min), hFiltered->FindFixBin(right_subrange_max));
        double Multiplicity = Signal/numberOfEvents;
        double SNR = Signal/Background;
        std::cout << Signal << ", " << Background << ", " << SNR << ", " << Multiplicity << ", " << numberOfEvents << std::endl;

        TLine *line_horizon = new TLine(hFiltered->GetXaxis()->GetXmin(), 0, hFiltered->GetXaxis()->GetXmax(), 0);
        line_horizon->SetLineColor(kRed);
        line_horizon->Draw("same");
        TLegend *legend2 = new TLegend(0.13, 0.75, 0.29, 0.9);
        legend2->AddEntry(hFiltered, Form("Multiplicity Density: %.2f", Multiplicity), "f");
        legend2->AddEntry(hFiltered, Form("Signal/Background Ratio: %.2f",SNR), "l");
        legend2->Draw("same");

        can->Update();
        std::string fileName = rootFilePath.substr(0, rootFilePath.length()-5);
        can -> SaveAs(Form("%s_%s.png", fileName.c_str(), option.Data()));
        i++;
    }
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

void dPhi_Fitter(std::string opt = "") {
    TStopwatch timer;   timer.Start();

    std::string dirPath = "../External/zFindingPlots";  // The directory containing the files
    std::string filePrefix = "hDiff_with_dEta_cut";  // Prefix of the file you are looking for
    std::vector<std::string> rootFilePaths = findRootFiles(dirPath, filePrefix);
    if (rootFilePaths.empty()) {
        std::cerr << "Error: Cannot find any .root files with the prefix hDiff_.." << filePrefix << "!" << std::endl;
        exit(1);
    }

    printRed("Found .root files:");
    printSeparation();
    for (const std::string &rootFilePath : rootFilePaths) {
        std::cout << rootFilePath << std::endl;
    }
    printSeparation();
    std::cout << opt << std::endl;
    // // TFile *inputFile = new TFile("../External/zFindingPlots/hDiff_output.root", "READ");
    // TFile *inputFile = new TFile(rootFilePath.c_str(), "READ");
    // if (!inputFile || inputFile->IsZombie()) {
    //     std::cerr << "Error: Cannot open .root file!" << std::endl;
    //     return;
    // }

    // TH1D *hDiff = (TH1D*)inputFile->Get("Background Subtracted Signal");
    // if (!hDiff) {
    //     std::cerr << "Error: Cannot find the target histogram!" << std::endl;
    //     inputFile->Close();
    //     return;
    // }

    // TCanvas *can = new TCanvas("can", "Background Subtracted Signal", 1920, 1056);
    // can -> Divide(1, 2);
    // can -> cd(1);
    // hDiff->SetLineColor(kBlue);
    // hDiff->SetLineWidth(2);
    // // hDiff->SetMarkerStyle(0);  // Disable markers
    // // hDiff->SetMarkerSize(10);   // Ensure no marker size
    // hDiff->Draw("hist");

    // // Add additional processing or fitting here if needed
    // // Example: fitting a Gaussian function to the histogram
    // // hDiff->Fit("gaus");
    
    // // Cosmetics:
    // TLine *line_horizontal = new TLine(hDiff->GetXaxis()->GetXmin(), 0, hDiff->GetXaxis()->GetXmax(), 0);
    // line_horizontal->SetLineColor(kRed);
    // line_horizontal->Draw("same");
    // // TLine *line_left = new TLine((hDiff->GetXaxis()->GetXmin())/5., hDiff->GetMinimum(), (hDiff->GetXaxis()->GetXmin())/5., hDiff->GetMaximum());
    // // line_left->SetLineColor(kRed);
    // // line_left->Draw("same");
    // TLegend *legend = new TLegend(0.7, 0.8, 0.9, 0.9);
    // legend->AddEntry(hDiff, "Background Subtracted Signal", "l");
    // legend->Draw();

    std::vector<std::string> selections = splitString(opt, '_');
    if (selections[0] == "single")    fit_the_hist(rootFilePaths, selections[1], std::stoi(selections[2]));
    if (selections[0] == "all")       fit_hists_all(rootFilePaths, selections[1]);
    if (selections[0] == "A")         plotMultiCen();
    // fit_the_hist_ver2(hDiff, can);

    // Stop the stopwatch and print the runtime:
    timer.Stop();   timer.Print();
}
