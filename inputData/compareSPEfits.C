// Per-channel comparison of log-normal fit parameters: SPE-only pulse injection
// (fitSPEonly.C output) vs the offline data fits (allFits_v4.root).
// One plot = one channel, one parameter, two normalised histograms.
//
// Selection, the same on both sides:
//   data: status_ln == 0 && 0.5 <= redChi2_ln <= 1.5 (already harvested before fitting)
//   sim : passHarvest && status_ln == 0 && 0.5 <= redChi2_ln <= 1.5
//
// Output in outDir: <param>_per_channel.pdf (one page per channel) for U0, tau,
// sigma, peakHeight, fitRiseNs_ln, fitFallNs_ln; compareSPEfits.root with all
// histograms; compareSPEfits_means.csv with the per-channel means.
//
//   root -l -b -q 'compareSPEfits.C+("SPEonly_2kEvent_fits.root","allFits_v4.root","SPEfitCompare")'
#include "TCanvas.h"
#include "TFile.h"
#include "TH1D.h"
#include "TLegend.h"
#include "TMath.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TTree.h"

#include <algorithm>
#include <cstdio>
#include <iostream>
#include <string>
#include <vector>

namespace {
const int kNChan = 96;
const std::vector<std::string> kParams = {"U0", "tau", "sigma", "peakHeight", "fitRiseNs_ln", "fitFallNs_ln"};
const std::vector<std::string> kTitles = {"U_{0} (mV)", "#tau (ns)", "#sigma", "peak height (mV)",
                                          "fitted 20-80% rise (ns)", "fitted 80-20% fall (ns)"};

// values[param][chan]
using Values = std::vector<std::vector<std::vector<float>>>;

void readTree(const char* fileName, bool isSim, Values& v) {
    TFile* f = TFile::Open(fileName, "READ");
    TTree* t = (TTree*)f->Get("speFits");
    Int_t chan, statusLN, pass = 1;
    Double_t redChi2LN, x[6];
    t->SetBranchStatus("*", 0);
    for (const char* b : {"chan", "status_ln", "redChi2_ln"}) t->SetBranchStatus(b, 1);
    t->SetBranchAddress("chan", &chan);
    t->SetBranchAddress("status_ln", &statusLN);
    t->SetBranchAddress("redChi2_ln", &redChi2LN);
    if (isSim) {
        t->SetBranchStatus("passHarvest", 1);
        t->SetBranchAddress("passHarvest", &pass);
    }
    for (size_t p = 0; p < kParams.size(); ++p) {
        t->SetBranchStatus(kParams[p].c_str(), 1);
        t->SetBranchAddress(kParams[p].c_str(), &x[p]);
    }
    const Long64_t n = t->GetEntries();
    for (Long64_t i = 0; i < n; ++i) {
        t->GetEntry(i);
        if (!pass || statusLN != 0 || redChi2LN < 0.5 || redChi2LN > 1.5) continue;
        for (size_t p = 0; p < kParams.size(); ++p) v[p][chan].push_back(x[p]);
    }
    f->Close();
}

double quantile(std::vector<float> a, double q) {
    if (a.empty()) return 0;
    size_t k = std::min(a.size() - 1, size_t(q * a.size()));
    std::nth_element(a.begin(), a.begin() + k, a.end());
    return a[k];
}

double mean(const std::vector<float>& a) {
    double s = 0;
    for (float x : a) s += x;
    return a.empty() ? 0 : s / a.size();
}
}  // namespace

void compareSPEfits(const char* simFile = "SPEonly_2kEvent_fits.root",
                    const char* dataFile = "allFits_v4.root",
                    const char* outDir = "SPEfitCompare") {
    Values sim(kParams.size(), std::vector<std::vector<float>>(kNChan));
    Values data(kParams.size(), std::vector<std::vector<float>>(kNChan));
    readTree(simFile, true, sim);
    readTree(dataFile, false, data);

    gSystem->mkdir(outDir, kTRUE);
    gStyle->SetOptStat(0);
    TFile out(Form("%s/compareSPEfits.root", outDir), "RECREATE");
    FILE* csv = fopen(Form("%s/compareSPEfits_means.csv", outDir), "w");
    fprintf(csv, "param,chan,nSim,meanSim,nData,meanData,ratioSimOverData\n");

    TCanvas c("c", "c", 800, 600);
    for (size_t p = 0; p < kParams.size(); ++p) {
        TString pdf = Form("%s/%s_per_channel.pdf", outDir, kParams[p].c_str());
        c.Print(pdf + "[");
        for (int ch = 0; ch < kNChan; ++ch) {
            const std::vector<float>& s = sim[p][ch];
            const std::vector<float>& d = data[p][ch];
            // Range from the 0.5% and 99.5% quantiles of both samples, padded 10%
            double lo = std::min(quantile(s, 0.005), quantile(d, 0.005));
            double hi = std::max(quantile(s, 0.995), quantile(d, 0.995));
            const double pad = 0.1 * (hi - lo > 0 ? hi - lo : 1.0);
            lo -= pad;
            hi += pad;

            TH1D* hd = new TH1D(Form("data_%s_ch%02d", kParams[p].c_str(), ch),
                                Form("Channel %d;%s;normalised", ch, kTitles[p].c_str()), 80, lo, hi);
            TH1D* hs = new TH1D(Form("sim_%s_ch%02d", kParams[p].c_str(), ch),
                                Form("Channel %d;%s;normalised", ch, kTitles[p].c_str()), 80, lo, hi);
            for (float x : d) hd->Fill(x);
            for (float x : s) hs->Fill(x);
            out.cd();
            hd->Write();
            hs->Write();
            if (hd->Integral() > 0) hd->Scale(1.0 / hd->Integral());
            if (hs->Integral() > 0) hs->Scale(1.0 / hs->Integral());

            hd->SetLineColor(kGray + 2);
            hd->SetFillColorAlpha(kGray, 0.5);
            hs->SetLineColor(kRed + 1);
            hs->SetLineWidth(2);
            hd->SetMaximum(1.15 * std::max(hd->GetMaximum(), hs->GetMaximum()));
            hd->Draw("HIST");
            hs->Draw("HIST SAME");

            const double ms = mean(s), md = mean(d);
            TLegend leg(0.50, 0.74, 0.89, 0.89);
            leg.SetBorderSize(0);
            leg.SetFillStyle(0);
            leg.AddEntry(hs, Form("pulse injection (N=%zu, mean %.3g)", s.size(), ms), "l");
            leg.AddEntry(hd, Form("offline fitting data (N=%zu, mean %.3g)", d.size(), md), "f");
            leg.Draw();
            c.Print(pdf, Form("Title:ch%02d", ch));
            fprintf(csv, "%s,%d,%zu,%.5g,%zu,%.5g,%.4f\n", kParams[p].c_str(), ch, s.size(), ms, d.size(), md,
                    md != 0 ? ms / md : 0.0);
            delete hd;
            delete hs;
        }
        c.Print(pdf + "]");
        std::cout << "Wrote " << pdf << std::endl;
    }
    fclose(csv);
    out.Close();
}
