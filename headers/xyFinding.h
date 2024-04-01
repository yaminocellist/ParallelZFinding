#include "globalDefinitions.h"

double gMedian(const TH1D * h1) {
   int n = h1->GetXaxis()->GetNbins();
   std::vector<double>  x(n);
   h1->GetXaxis()->GetCenter( &x[0] );
   const double * y = h1->GetArray();
   // exclude underflow/overflows from bin content array y
   return TMath::Median(n, &x[0], &y[1]);
}

double angularDistance (const int & evt, const vector<myTrackletMember> &t0, const vector<myTrackletMember> &t1) {
    const int n1 = 100;         // number of bins in [0, 0.1]
    const int n2 = 10;          // number of bins in [0.1, 1.1], [1.1, 2.1], ...
    const int n3 = 8 * n2 + n1; // total number of bins
    float unequal_bins[n3 + 1]; // array of bin edges
    // Set bin edges
    for (int i = 0; i < n1; ++i) {
        unequal_bins[i] = 0.1 * static_cast<double>(i) / n1;  // divide [0, 0.1] into n1 bins
    }
    for (int i = 0; i < 8; ++i) {
        for (int j = 0; j < n2; ++j) {
            unequal_bins[n1 + i * n2 + j] = 0.1 + static_cast<double>(i) + static_cast<double>(j) / n2;  // divide [0.1 + i, 1.1 + i] into n2 bins
        }
    }
    unequal_bins[n3] = 8.1;  // set the last bin edge
    TH1F *h = new TH1F(Form("histo_%d", evt), Form("Angular Distance Histogram of Event %d", evt), n3, unequal_bins);
    // TH1F *h = new TH1F("", "", 200, 0., 8.);
    // std::cout << t0.size() << "," << t1.size() << std::endl;
    for (int i = 0; i < t0.size(); i++) {
        for (int j = 0; j < t1.size(); j++) {
            double angular_distance = std::sqrt((t0[i].eta - t1[j].eta)*(t0[i].eta - t1[j].eta) + (t0[i].phi - t1[j].phi)*(t0[i].phi - t1[j].phi));
            h -> Fill(angular_distance);
        }
    }
    // std::cout << "Hi" << std::endl;
    int bin1 = h->FindBin(0.01);  // find the bin number corresponding to 0.0
    int bin2 = h->FindBin(0.02);  // find the bin number corresponding to 0.02
    double minContent = h->GetBinContent(bin1);
    int minBin = bin1;  
    for (int i = bin1 + 1; i <= bin2; ++i) {
        double currentContent = h -> GetBinContent(i);
        if (currentContent < minContent) {
            minContent = currentContent;
            minBin = i;
        }
    }
    double rst = (h -> GetBinCenter(minBin)) - 0.1/(n1*2);

    TCanvas *can2 = new TCanvas("c2","c2",0,50,1800,550);
    h -> Draw();
    h -> SetFillColor(kYellow - 7);
    h -> SetLineWidth(1);
    h -> SetFillStyle(1001);
    h -> GetXaxis() -> SetTitle("Angular Distance");
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
    h -> GetXaxis() -> SetRangeUser(0, 0.2); // Setting x range;
    h -> GetYaxis()->SetRangeUser(0, 1000);
    gPad -> SetGrid(1,1); gPad -> Update();
    // gPad -> SetLogx();
    can2 -> SaveAs(Form("../External/xyFindingPlots/angular_distance_%d.png", evt));
    delete can2;
    h -> Reset("ICESM");    delete h;   h = 0;
    
    // std::cout << rst << std::endl;
    return rst;
}

std::pair<double, double> findXY (const int & evt, const vector<myTrackletMember> &t0, const vector<myTrackletMember> &t1, const double &angular_distance_cut_squared, const double &found_z) {
    TH2D *h_intersection = new TH2D("", "", 1800, -4.5, 4.5, 1800, -4.5, 4.5);
    h_intersection -> Sumw2(kFALSE);
    for (int i = 0; i < t0.size(); i++) {
        for (int j = 0; j < t1.size(); j++) {
            double a_d_squared = (t0[i].eta - t1[j].eta)*(t0[i].eta - t1[j].eta) + 
                                 (t0[i].phi - t1[j].phi)*(t0[i].phi - t1[j].phi);
            if (a_d_squared <= angular_distance_cut_squared) {
                double t = (found_z - t0[i].z) / (t1[j].z - t0[i].z);
                double x_at_found_z = t0[i].x + t * (t1[j].x - t0[i].x);
                double y_at_found_z = t0[i].y + t * (t1[j].y - t0[i].y);
                h_intersection -> Fill(x_at_found_z, y_at_found_z);
            }       
        }
    }

    TH1D *h_proj_x = h_intersection -> ProjectionX("h_proj_x");
    TH1D *h_proj_y = h_intersection -> ProjectionY("h_proj_y");
    double median_x = gMedian(h_proj_x);
    double median_y = gMedian(h_proj_y);
    TCanvas *canvas = new TCanvas("c", "c", 0,50,1000,1000);
    gPad->SetMargin(0.1, 0.1, 0.1, 0.1);
    gPad->SetGrid(1, 1);

    h_intersection -> Draw("colz");
    h_intersection -> GetYaxis() -> SetRangeUser(-3., 3.);
    h_intersection -> GetXaxis() -> SetRangeUser(-3., 3.);
    h_intersection -> GetXaxis() -> SetTitle("X (cm)");
    h_intersection -> GetYaxis() -> SetTitle("Y (cm)");
    h_intersection -> GetYaxis() -> SetTitleOffset(1);
    h_intersection -> GetXaxis() -> SetTitleSize(0.05);
    h_intersection -> GetYaxis() -> SetTitleSize(0.05); 
    h_intersection -> SetTitle(Form("#splitline{Data from event #%d}{      Angular Separation cut at %1.3f}", evt, sqrt(angular_distance_cut_squared)));
    gStyle -> SetTitleSize(0.032, "t");
    h_intersection -> GetXaxis() -> CenterTitle(true);
    h_intersection -> GetYaxis() -> CenterTitle(true);
    Color_t colorS = kTeal - 7;
    TMarker *markerS = new TMarker(0, 0, 21); // 8 is the default marker style for TH2D
    markerS->SetMarkerColor(colorS);
    h_intersection -> SetMarkerColor(colorS);
    TLegend *lg = new TLegend(0.12, 0.85, 0.7, 0.9);
    lg -> AddEntry(markerS, Form("Found x and y vertex at: %1.5f, %1.5f (cm)", median_x, median_y), "P");
    gStyle -> SetLegendTextSize(.02);
    lg -> Draw("same");
    canvas -> SaveAs(Form("../External/xyFindingPlots/event_%d.png", evt));
    delete canvas;
    h_intersection -> Reset("ICESM");
    h_proj_y -> Reset("ICESM");
    h_proj_x -> Reset("ICESM");
    delete h_intersection;  delete h_proj_x;    delete h_proj_y;
    h_intersection = 0; h_proj_x = 0;   h_proj_y = 0;
    std::cout << median_x << "," << median_y << std::endl;
    return std::make_pair(median_x, median_y);
}