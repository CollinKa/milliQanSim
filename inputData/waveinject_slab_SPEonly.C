// SPE-only slab pulse injection, for validating the log-normal SPE templates.
// No Geant4 input: every event has exactly one PE in each of the 96 channels
// (summing is ignored, one channel = one PE). Channel 8 alternates between its
// two PMT shapes: even events use 8A (fast tau), odd events 8B (slow tau).
// No QE rescaling, since only SPEs are studied.
//
// Pulse building follows the per-hit path of waveinject_slab.C: sample the SPE
// area from SPEmeans.txt, scale the template by area / templateArea in bins
// 500-660, shift by the hit time, zero outside 500-660, add noise, clamp.
// makeLogNormalTemplate and its constants are copied from waveinject_slab.C;
// keep the two in sync.
//
// Output tree "Events": waveform[6][16][1024], eventWeight, and
// injectedArea[6][16] (pVs, the sampled area per channel, for the harvest-window
// selection when comparing with the v4 data fits).
//
//   root -l -b -q 'waveinject_slab_SPEonly.C+("SPEonly_2kEvent_waveinjected.root", 2000)'
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TF1.h"
#include "TRandom3.h"
#include "TString.h"

#include <cmath>
#include <cstring>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

// Log-normal SPE template on the digitizer binning (1024 bins x 2.5 ns):
// f(t) = U0 * exp(-0.5 * (ln((t-t0)/tau) / sigma)^2) + ped for t > t0, else ped.
const double kTemplateU0 = 50.0;    // mV
const double kTemplateT0 = 1285.0;  // ns
const double kTemplatePed = 0.0;    // mV

TH1F* makeLogNormalTemplate(const char* name, double tau, double sigma) {
    TH1F* h = new TH1F(name, name, 1024, 0, 2560);
    h->SetDirectory(nullptr);
    for (int bin = 1; bin <= h->GetNbinsX(); ++bin) {
        double dt = h->GetBinCenter(bin) - kTemplateT0;
        double value = kTemplatePed;
        if (dt > 0) value += kTemplateU0 * std::exp(-0.5 * std::pow(std::log(dt / tau) / sigma, 2));
        h->SetBinContent(bin, value);
    }
    return h;
}

void waveinject_slab_SPEonly(
    TString outputFile = "SPEonly_2kEvent_waveinjected.root",
    Long64_t nEvents = 2000,
    double hitTime = 0.0,  // ns, same for every PE; early so the pulse stays inside bins 500-660
    TString lognormalFile = "SPElognormal.txt",
    TString speMeansFile = "SPEmeans.txt") {

    std::cout << "Outputting file " << outputFile << std::endl;
    std::cout << "Events " << nEvents << ", hit time " << hitTime << " ns" << std::endl;
    std::cout << "Using log-normal SPE shapes " << lognormalFile << std::endl;
    std::cout << "Using SPE means file " << speMeansFile << std::endl;

    const int nDigitizers = 6;
    const int nChannelsPerDigitizer = 16;
    const int nSPEChan = 96;
    const int nBins = 1024;
    double binWidth = 2.5;
    double rms_noise = 1;
    const double nVsToPVs = 1e3;  // SPEmeans.txt is in nVs; sampling uses pVs

    TFile* outfile = new TFile(outputFile, "RECREATE");
    Float_t waveform[nDigitizers][nChannelsPerDigitizer][nBins] = {{{0}}};
    Float_t injectedArea[nDigitizers][nChannelsPerDigitizer] = {{0}};
    Double_t eventWeight = 1;

    TTree* injectedTree = new TTree("Events", "SPE-only injected digitizer waveforms");
    injectedTree->Branch("waveform", waveform, Form("waveform[%d][%d][%d]/F", nDigitizers, nChannelsPerDigitizer, nBins));
    injectedTree->Branch("eventWeight", &eventWeight, "eventWeight/D");
    injectedTree->Branch("injectedArea", injectedArea, Form("injectedArea[%d][%d]/F", nDigitizers, nChannelsPerDigitizer));

    // Per-channel SPE area Gaussians from SPEmeans.txt (mean/sigma in nVs -> pVs).
    // Range [0, 30000] as in waveinject_slab.C: TF1::GetRandom only samples
    // inside the range, and 30000 pVs is >= mean + 3 sigma for every channel.
    std::vector<TF1*> fit(nSPEChan, nullptr);
    {
        std::ifstream in(speMeansFile.Data());
        if (!in.is_open()) {
            std::cerr << "Failed to open SPE means file: " << speMeansFile << std::endl;
            return;
        }
        std::string line;
        std::getline(in, line);  // header
        while (std::getline(in, line)) {
            std::istringstream iss(line);
            int chan = -1;
            double mean = 0, meanErr = 0, sigma = 0;
            iss >> chan >> mean >> meanErr >> sigma;
            double m = mean * nVsToPVs, s = sigma * nVsToPVs;
            fit[chan] = new TF1(Form("fit_chan%d", chan), "gaus(0)", 0, 30000);
            fit[chan]->SetParameters(1.0, m, s);
        }
    }

    // Per-channel log-normal SPE templates from SPElognormal.txt; 8A/8B are
    // channel 8's two PMT shapes, row 8 its mean (used for the area A).
    std::vector<TH1F*> speTemplate(nSPEChan, nullptr);
    TH1F* speTemplate8A = nullptr;
    TH1F* speTemplate8B = nullptr;
    {
        std::ifstream in(lognormalFile.Data());
        if (!in.is_open()) {
            std::cerr << "Failed to open log-normal SPE file: " << lognormalFile << std::endl;
            return;
        }
        std::string line;
        std::getline(in, line);  // header
        while (std::getline(in, line)) {
            std::istringstream iss(line);
            std::string label;
            double tau = 0, tauErr = 0, sigma = 0;
            iss >> label >> tau >> tauErr >> sigma;
            TH1F* h = makeLogNormalTemplate(Form("spe_template_%s", label.c_str()), tau, sigma);
            if (label == "8A") speTemplate8A = h;
            else if (label == "8B") speTemplate8B = h;
            else speTemplate[std::stoi(label)] = h;
        }
    }

    // Template area A per channel in pVs (mV x ns), over the whole record.
    // Channel 8 uses the mean-tau template (row 8) for 8A and 8B as well.
    std::vector<double> templateArea(nSPEChan, 0.0);
    for (int ch = 0; ch < nSPEChan; ++ch) {
        templateArea[ch] = speTemplate[ch]->Integral(1, speTemplate[ch]->GetNbinsX()) * binWidth;
    }

    // Max values: 96 entries of 1250, as in waveinject_slab.C (b679b4d)
    std::vector<double> maxValues(nSPEChan, 1250);
    maxValues[23] = maxValues[37] = maxValues[38] = 1255;

    TRandom3 randGen(2004);
    const int integer_shift = static_cast<int>(hitTime / binWidth);
    const double fractional_shift = fmod(hitTime, binWidth) / binWidth;

    for (Long64_t i = 0; i < nEvents; i++) {
        if (i % (nEvents / 10) == 0) std::cout << "Processing Event " << i << "..." << std::endl;
        memset(waveform, 0, sizeof(waveform));
        memset(injectedArea, 0, sizeof(injectedArea));

        for (int ch = 0; ch < nSPEChan; ++ch) {
            int digitizer = ch / nChannelsPerDigitizer;
            int channel = ch % nChannelsPerDigitizer;

            TH1F* pulse_shape = speTemplate[ch];
            if (ch == 8) pulse_shape = (i % 2 == 0) ? speTemplate8A : speTemplate8B;

            TH1F* new_waveform = (TH1F*)pulse_shape->Clone();
            double event_area = fit[ch]->GetRandom();
            injectedArea[digitizer][channel] = event_area;

            // Scale only the bins within [500, 660]
            double scaleFactor = event_area / templateArea[ch];
            for (int bin = 500; bin <= 660; ++bin) {
                new_waveform->SetBinContent(bin, new_waveform->GetBinContent(bin) * scaleFactor);
            }

            for (int bin = 0; bin < nBins; ++bin) {
                int shifted_bin = bin - integer_shift;
                if (shifted_bin >= nBins || shifted_bin < 1) continue;
                double value = (1.0 - fractional_shift) * new_waveform->GetBinContent(shifted_bin) +
                               fractional_shift * new_waveform->GetBinContent(shifted_bin - 1);
                waveform[digitizer][channel][bin] += value;
            }
            for (int bin = 0; bin < integer_shift; ++bin) waveform[digitizer][channel][bin] = 0.0;
            delete new_waveform;

            for (int bin = 0; bin < nBins; bin++) {
                if (bin < 500 || bin > 660) waveform[digitizer][channel][bin] = 0;
                waveform[digitizer][channel][bin] += randGen.Gaus(0, rms_noise);
                if (waveform[digitizer][channel][bin] > maxValues[ch]) waveform[digitizer][channel][bin] = maxValues[ch];
                if (waveform[digitizer][channel][bin] < -50) waveform[digitizer][channel][bin] = -50;
            }
        }
        injectedTree->Fill();
    }

    outfile->cd();
    injectedTree->Write();
    outfile->Close();
    for (int ch = 0; ch < nSPEChan; ++ch) {
        delete fit[ch];
        delete speTemplate[ch];
    }
    delete speTemplate8A;
    delete speTemplate8B;
}
