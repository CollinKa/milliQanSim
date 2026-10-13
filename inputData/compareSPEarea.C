// Pulse-area check of the SPE-only injection (waveinject_slab_SPEonly.C output):
// the area of each injected waveform vs the area it was assigned.
//
// Waveform area = sum of all 1024 bins x 2.5 ns (pVs), for every pulse, with no
// pulse finding and no selection. Compared per channel with
//   - injectedArea: the area drawn for that pulse from the SPEmeans.txt Gaussian
//   - the input Gaussian itself (SPEmeans.txt mean, sigma), not truncated, scaled
//     to the number of pulses. The CSV also gives the mean and RMS of the
//     Gaussian truncated to the sampling range [0, 30000], which the injected
//     areas should match.
//
// Output in outDir:
//   SPEarea_per_channel.pdf        one page per channel: waveform area and the input Gaussian
//   SPEarea_diff_per_channel.pdf   one page per channel: waveform area - injected area (pVs)
//   compareSPEarea_summary.csv     per-channel means and RMS; channel 8 also split into 8A/8B
//   area_chNN.png, diff_chNN.png   previews for channels 0, 8 and 72
//
//   root -l -b -q 'compareSPEarea.C+("SPEonly_2kEvent_waveinjected.root","SPEmeans.txt","SPEareaCompare")'
#include "TCanvas.h"
#include "TF1.h"
#include "TFile.h"
#include "TH1D.h"
#include "TLegend.h"
#include "TMath.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TTree.h"

#include <cmath>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

namespace {
struct Stat {
    double n = 0, s = 0, s2 = 0;
    void add(double x) { n++; s += x; s2 += x * x; }
    double mean() const { return n ? s / n : 0; }
    double rms() const { return n > 1 ? std::sqrt(std::max(0.0, s2 / n - mean() * mean())) : 0; }
};
}  // namespace

void compareSPEarea(const char* inputFile = "SPEonly_2kEvent_waveinjected.root",
                    const char* speMeansFile = "SPEmeans.txt",
                    const char* outDir = "SPEareaCompare") {
    const int nDigitizers = 6, nChannelsPerDigitizer = 16, nBins = 1024, nChan = 96;
    const double binWidth = 2.5, areaMax = 30000;  // same sampling range as the injection

    std::vector<double> gMean(nChan), gSigma(nChan);
    {
        std::ifstream in(speMeansFile);
        std::string line;
        std::getline(in, line);  // header
        while (std::getline(in, line)) {
            std::istringstream iss(line);
            int chan;
            double mean, meanErr, sigma;
            iss >> chan >> mean >> meanErr >> sigma;
            gMean[chan] = mean * 1e3;
            gSigma[chan] = sigma * 1e3;
        }
    }

    TFile* f = TFile::Open(inputFile, "READ");
    TTree* t = (TTree*)f->Get("Events");
    static Float_t waveform[nDigitizers][nChannelsPerDigitizer][nBins];
    Float_t injectedArea[nDigitizers][nChannelsPerDigitizer];
    t->SetBranchAddress("waveform", waveform);
    t->SetBranchAddress("injectedArea", injectedArea);

    // area[ch], injected[ch] per pulse
    std::vector<std::vector<float>> area(nChan), injected(nChan);
    const Long64_t nEvents = t->GetEntries();
    for (Long64_t i = 0; i < nEvents; ++i) {
        t->GetEntry(i);
        for (int ch = 0; ch < nChan; ++ch) {
            const int d = ch / nChannelsPerDigitizer, c = ch % nChannelsPerDigitizer;
            double sum = 0;
            for (int b = 0; b < nBins; ++b) sum += waveform[d][c][b];
            area[ch].push_back(sum * binWidth);
            injected[ch].push_back(injectedArea[d][c]);
        }
    }
    f->Close();

    gSystem->mkdir(outDir, kTRUE);
    gStyle->SetOptStat(0);
    FILE* csv = fopen(Form("%s/compareSPEarea_summary.csv", outDir), "w");
    fprintf(csv, "chan,n,gausMean,gausSigma,truncGausMean,truncGausRMS,injectedMean,injectedRMS,waveAreaMean,waveAreaRMS,"
                 "ratioOfMeans,diffMean,diffRMS,injMeanOverInput,injRMSOverInput,waveMeanOverInput,waveRMSOverInput\n");

    TCanvas c("c", "c", 800, 600);
    TString pdfArea = Form("%s/SPEarea_per_channel.pdf", outDir);
    TString pdfRatio = Form("%s/SPEarea_diff_per_channel.pdf", outDir);
    c.Print(pdfArea + "[");
    for (int ch = 0; ch < nChan; ++ch) {
        // Sampling range of the injection
        const double sLo = 0.0, sHi = areaMax;
        const double hi = std::min(areaMax, gMean[ch] + 5 * gSigma[ch]);
        const double lo = -0.05 * hi;
        TH1D hWave(Form("wave_ch%02d", ch), Form("Channel %d;pulse area (pVs);pulses", ch), 100, lo, hi);
        Stat sWave, sInj, sDiff;
        for (size_t k = 0; k < area[ch].size(); ++k) {
            hWave.Fill(area[ch][k]);
            sWave.add(area[ch][k]);
            sInj.add(injected[ch][k]);
            sDiff.add(area[ch][k] - injected[ch][k]);
        }

        // Input Gaussian (SPEmeans.txt mean, sigma), not truncated, scaled to the
        // number of pulses: what the injected areas are meant to follow.
        TF1 g(Form("g_ch%02d", ch), "[0]*TMath::Gaus(x,[1],[2],1)", lo, hi);
        g.SetParameters(area[ch].size() * hWave.GetBinWidth(1), gMean[ch], gSigma[ch]);
        g.SetNpx(500);
        // Mean and RMS of the Gaussian truncated to the sampling range: what the
        // injected areas should average to.
        TF1 g1("g1", "TMath::Gaus(x,[0],[1])", sLo, sHi);
        TF1 gx("gx", "x*TMath::Gaus(x,[0],[1])", sLo, sHi);
        TF1 gxx("gxx", "x*x*TMath::Gaus(x,[0],[1])", sLo, sHi);
        for (TF1* fn : {&g1, &gx, &gxx}) fn->SetParameters(gMean[ch], gSigma[ch]);
        const double i0 = g1.Integral(sLo, sHi);
        const double tMean = gx.Integral(sLo, sHi) / i0;
        const double tRMS = std::sqrt(std::max(0.0, gxx.Integral(sLo, sHi) / i0 - tMean * tMean));

        hWave.SetLineColor(kRed + 1);
        hWave.SetLineWidth(2);
        g.SetLineColor(kBlack);
        g.SetLineWidth(2);
        hWave.SetMaximum(1.2 * std::max(hWave.GetMaximum(), g.GetMaximum(lo, hi)));
        hWave.Draw("HIST");
        g.Draw("SAME");
        TLegend leg(0.45, 0.70, 0.89, 0.89);
        leg.SetBorderSize(0);
        leg.SetFillStyle(0);
        leg.AddEntry(&hWave, Form("waveform area (mean %.0f)", sWave.mean()), "l");
        leg.AddEntry(&g, Form("input Gaussian %.0f #pm %.0f", gMean[ch], gSigma[ch]), "l");
        leg.AddEntry((TObject*)nullptr, Form("sampled on [0, %.0f]", areaMax), "");
        leg.Draw();
        c.Print(pdfArea, Form("Title:ch%02d", ch));
        if (ch == 0 || ch == 8 || ch == 72) c.SaveAs(Form("%s/area_ch%02d.png", outDir, ch));

        fprintf(csv, "%d,%.0f,%.1f,%.1f,%.1f,%.1f,%.1f,%.1f,%.1f,%.1f,%.5f,%.1f,%.1f,%.5f,%.5f,%.5f,%.5f\n", ch, sWave.n,
                gMean[ch], gSigma[ch], tMean, tRMS, sInj.mean(), sInj.rms(), sWave.mean(), sWave.rms(),
                sWave.mean() / sInj.mean(), sDiff.mean(), sDiff.rms(), sInj.mean() / gMean[ch], sInj.rms() / gSigma[ch],
                sWave.mean() / gMean[ch], sWave.rms() / gSigma[ch]);
        if (ch == 8) {
            // Even events use 8A, odd 8B (waveinject_slab_SPEonly.C)
            for (int type = 0; type < 2; ++type) {
                Stat sw, si, sd;
                for (size_t k = type; k < area[ch].size(); k += 2) {
                    sw.add(area[ch][k]);
                    si.add(injected[ch][k]);
                    sd.add(area[ch][k] - injected[ch][k]);
                }
                fprintf(csv, "%s,%.0f,%.1f,%.1f,%.1f,%.1f,%.1f,%.1f,%.1f,%.1f,%.5f,%.1f,%.1f,%.5f,%.5f,%.5f,%.5f\n",
                        type ? "8B" : "8A", sw.n, gMean[ch], gSigma[ch], tMean, tRMS, si.mean(), si.rms(), sw.mean(),
                        sw.rms(), sw.mean() / si.mean(), sd.mean(), sd.rms(), si.mean() / gMean[ch],
                        si.rms() / gSigma[ch], sw.mean() / gMean[ch], sw.rms() / gSigma[ch]);
            }
        }
    }
    c.Print(pdfArea + "]");

    c.Print(pdfRatio + "[");
    for (int ch = 0; ch < nChan; ++ch) {
        TH1D hD(Form("diff_ch%02d", ch), Form("Channel %d;waveform area - injected area (pVs);pulses", ch), 100, -2000, 2000);
        TH1D hA(Form("diffA_ch%02d", ch), "", 100, -2000, 2000);  // channel 8 only: 8A
        TH1D hB(Form("diffB_ch%02d", ch), "", 100, -2000, 2000);  // channel 8 only: 8B
        double sumW = 0, sumI = 0;
        for (size_t k = 0; k < area[ch].size(); ++k) {
            const double dA = area[ch][k] - injected[ch][k];
            hD.Fill(dA);
            if (ch == 8) (k % 2 == 0 ? hA : hB).Fill(dA);
            sumW += area[ch][k];
            sumI += injected[ch][k];
        }
        hD.SetLineColor(kRed + 1);
        hD.SetLineWidth(2);
        hD.Draw("HIST");
        TLegend leg(0.50, 0.72, 0.89, 0.89);
        leg.SetBorderSize(0);
        leg.SetFillStyle(0);
        leg.AddEntry(&hD, Form("all (mean %.0f, rms %.0f pVs)", hD.GetMean(), hD.GetRMS()), "l");
        leg.AddEntry((TObject*)nullptr, Form("mean waveform / mean injected = %.4f", sumI > 0 ? sumW / sumI : 0.0), "");
        if (ch == 8) {
            hA.SetLineColor(kAzure + 2);
            hB.SetLineColor(kGreen + 2);
            hA.SetLineWidth(2);
            hB.SetLineWidth(2);
            hA.Draw("HIST SAME");
            hB.Draw("HIST SAME");
            leg.AddEntry(&hA, Form("8A (mean %.0f pVs)", hA.GetMean()), "l");
            leg.AddEntry(&hB, Form("8B (mean %.0f pVs)", hB.GetMean()), "l");
        }
        leg.Draw();
        c.Print(pdfRatio, Form("Title:ch%02d", ch));
        if (ch == 0 || ch == 8 || ch == 72) c.SaveAs(Form("%s/diff_ch%02d.png", outDir, ch));
    }
    c.Print(pdfRatio + "]");
    fclose(csv);
    std::cout << "Wrote " << pdfArea << ", " << pdfRatio << " and " << outDir << "/compareSPEarea_summary.csv"
              << std::endl;
}
