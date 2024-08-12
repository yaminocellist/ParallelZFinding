#include "Math/WrappedMultiTF1.h"
#include "Math/WrappedParamFunction.h"
#include "TROOT.h"
#include "TVirtualFitter.h"
#include <TStopwatch.h>
#include <thread>
#include <mutex>
#include <chrono>
#include "headers/histogramPlotting.h"
#include "TCanvas.h"
#include "TApplication.h"

std::mutex m_mutex;
int N = 1000;
double range_min = -2*M_PI;
double range_max = 2*M_PI;

double bin_width = (range_max - range_min) / N;
TH1D *h_dEta_nomix = new TH1D("", "", N, range_min, range_max);
TH1D *h_Background_dEta = new TH1D("", "", N, range_min, range_max);
