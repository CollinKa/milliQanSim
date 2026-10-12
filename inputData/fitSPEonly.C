// Fit the SPE-only injected pulses (waveinject_slab_SPEonly.C output) exactly
// as the v4 data was fitted, so the two can be compared channel by channel.
//
// Per channel per event:
//  1. Offline pulse finding, ported from OfflineFactory::findPulses /
//     processChannel (milliqanOffline SlabPulseInjectionSim, pulseFindingSlab.json:
//     highThresh 10 mV, lowThresh 1 mV, nConsecSamples 3, nConsecSamplesEnd 3),
//     on the dynamic-pedestal-subtracted wave. Gives the offline time (fit seed),
//     area, height, riseSamples, fallSamples, duration, npulses, pickupFlagTight.
//  2. The data selection: extract_ids.C (npulses == 1, riseSamples >= 3,
//     fallSamples >= 3, pickupFlagTight == 0) and the make_spe_harvest.C area
//     window 1000*(mean +/- 1 sigma) from SPEmeans.txt, both edges exclusive.
//     prePulseMean and the timeFit_module_calibrated window are not applied: the
//     sim baseline is flat 1 mV noise and every pulse sits at ~1285 ns.
//  3. The v4 fit: fitPulse from runSpeFit_core.h (verbatim runSpeFit.cc),
//     log-normal model, --pre 200 --post 300, seeded with the offline time.
//     isSlab = false because the injected pulses are already positive.
//
// Every channel with a pulse is fitted and written; passHarvest marks the ones
// that pass step 2, so the selection can be applied (or not) when plotting.
//
//   root -l -b -q 'fitSPEonly.C+("SPEonly_2kEvent_waveinjected.root","SPEonly_2kEvent_fits.root")'
#include "runSpeFit_core.h"

#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"

#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

namespace {

// pulseFindingSlab.json
constexpr double kOffHighThresh = 10.0;  // mV
constexpr double kOffLowThresh = 1.0;    // mV
constexpr int kOffNConsec = 3;
constexpr int kOffNConsecEnd = 3;
constexpr double kOffSidebandEnd = 50.0;  // ns, findPulses starts at sideband_range[1] + 1

// OfflineFactory::findPulses, line for line.
std::vector<std::pair<float, float>> offlineFindPulses(TH1D *wave) {
    std::vector<std::pair<float, float>> bounds;
    int istart = wave->FindBin(kOffSidebandEnd + 1);
    bool inpulse = false;
    int nover = 0, nunder = 0;
    int i_begin = istart;
    int i_stop_searching = wave->GetNbinsX() - kOffNConsec;
    int i_stop_final_pulse = wave->GetNbinsX();
    for (int i = istart; i < i_stop_searching || (inpulse && i < i_stop_final_pulse); i++) {
        float v = wave->GetBinContent(i);
        if (!inpulse) {
            if (v < kOffLowThresh) {
                nover = 0;
                i_begin = i;
            } else if (v >= kOffHighThresh) {
                nover++;
            } else {
                i_begin = i;
            }
            if (nover >= kOffNConsec) {
                inpulse = true;
                nunder = 0;
            }
        } else {
            if (v < kOffHighThresh) nunder++;
            else nunder = 0;
            if (nunder >= kOffNConsecEnd || i == (i_stop_final_pulse - 1)) {
                bounds.push_back({(float)wave->GetBinLowEdge(i_begin), (float)wave->GetBinLowEdge(i + 1) - 0.01f});
                inpulse = false;
                nover = 0;
                nunder = 0;
                i_begin = i;
            }
        }
    }
    return bounds;
}

struct OfflinePulse {
    float time = -1, area = 0, height = 0, duration = 0;
    int riseSamples = 0, fallSamples = 0, npulses = 0;
    bool pickupFlagTight = false;
};

// OfflineFactory::processChannel for the first pulse (ipulse 0).
OfflinePulse offlinePulse(TH1D *wave) {
    OfflinePulse p;
    std::vector<std::pair<float, float>> pulseBounds = offlineFindPulses(wave);
    p.npulses = pulseBounds.size();
    if (p.npulses == 0) return p;
    const std::pair<float, float> b = pulseBounds[0];
    wave->SetAxisRange(b.first, b.second);
    float height = wave->GetMaximum();
    int above20 = -2, above80 = -1;
    for (int iStart = wave->FindBin(b.first); iStart <= wave->FindBin(b.second); iStart++) {
        if (wave->GetBinContent(iStart) > height * 0.2 && above20 < 0) above20 = iStart;
        if (wave->GetBinContent(iStart) > height * 0.8 && above80 < 0) above80 = iStart;
        if (above80 >= 0) break;
    }
    int riseSamples = above80 - above20;
    above20 = -2;
    above80 = -1;
    for (int iFall = wave->FindBin(b.second); iFall >= wave->FindBin(b.first); iFall--) {
        if (wave->GetBinContent(iFall) > height * 0.2 && above20 < 0) above20 = iFall;
        if (wave->GetBinContent(iFall) > height * 0.8 && above80 < 0) above80 = iFall;
        if (above80 >= 0) break;
    }
    int fallSamples = above20 - above80;
    float area = wave->Integral("width");
    float duration = b.second - b.first;
    bool qual = true;
    if (fallSamples < 2) qual = false;
    if (qual && riseSamples < 2) qual = false;
    if (qual && (height > 17. && height < 100.) && (0.001 * area < 0.04 * (height - 17.))) qual = false;
    if (qual && height < 25. && duration < 4.6 * (height - 12.)) qual = false;
    if (qual && height >= 25. && duration < 60. && fallSamples < 6) qual = false;
    if (qual && 0.001 * area < 0.4 && duration < 150. * (0.001 * area)) qual = false;
    if (qual && 0.001 * area >= 0.4 && duration < 60.) qual = false;
    bool qual_tight = false;
    if (!qual && riseSamples >= fallSamples && b.first > 70.0 && b.first < 2400.0 && p.npulses < 3 && duration < 100.0) qual_tight = true;
    wave->GetXaxis()->SetRange(0, 0);

    p.time = b.first;
    p.area = area;
    p.height = height;
    p.duration = duration;
    p.riseSamples = riseSamples;
    p.fallSamples = fallSamples;
    p.pickupFlagTight = qual_tight;
    return p;
}

}  // namespace

void fitSPEonly(const char *inputFile = "SPEonly_2kEvent_waveinjected.root",
                const char *outputFile = "SPEonly_2kEvent_fits.root",
                const char *speMeansFile = "SPEmeans.txt") {
    const int nDigitizers = 6, nChannelsPerDigitizer = 16, nBins = 1024, nSPEChan = 96;
    const double preWindow = 200.0, postWindow = 300.0;  // v4: --pre 200 --post 300

    // Harvest windows, make_spe_harvest.C with nSigma = 1
    std::vector<double> winLow(nSPEChan), winHigh(nSPEChan);
    {
        std::ifstream in(speMeansFile);
        std::string line;
        std::getline(in, line);  // header
        while (std::getline(in, line)) {
            std::istringstream iss(line);
            int chan;
            double mean, meanErr, sigma;
            iss >> chan >> mean >> meanErr >> sigma;
            winLow[chan] = 1000.0 * (mean - sigma);
            winHigh[chan] = 1000.0 * (mean + sigma);
        }
    }

    TFile *inFile = TFile::Open(inputFile, "READ");
    TTree *inTree = (TTree *)inFile->Get("Events");
    static Float_t waveform[nDigitizers][nChannelsPerDigitizer][nBins];
    Float_t injectedArea[nDigitizers][nChannelsPerDigitizer];
    inTree->SetBranchAddress("waveform", waveform);
    inTree->SetBranchAddress("injectedArea", injectedArea);

    TH1::AddDirectory(kFALSE);
    TFile *outFile = TFile::Open(outputFile, "RECREATE");
    TTree fits("speFits", "Fits of SPE-only injected pulses, v4 fit");
    Int_t tEvent, tChan, tCh8Type, tNpulses, tRise, tFall, tPickupTight, tPass, tStatusLN, tPedSamples;
    Double_t tInjArea, tOffArea, tOffHeight, tSeedT0, tDuration;
    Double_t tU0, tTau, tSigma, tPeak, tT0LN, tPedLN, tFitRiseLN, tFitFallLN, tRedChi2LN, tPedRMS;
    fits.Branch("event", &tEvent, "event/I");
    fits.Branch("chan", &tChan, "chan/I");
    fits.Branch("ch8Type", &tCh8Type, "ch8Type/I");  // channel 8: 0 = 8A (even event), 1 = 8B; -1 otherwise
    fits.Branch("injectedArea", &tInjArea, "injectedArea/D");
    fits.Branch("area", &tOffArea, "area/D");        // offline pulse-finder area
    fits.Branch("height", &tOffHeight, "height/D");  // offline pulse-finder height
    fits.Branch("npulses", &tNpulses, "npulses/I");
    fits.Branch("riseSamples", &tRise, "riseSamples/I");
    fits.Branch("fallSamples", &tFall, "fallSamples/I");
    fits.Branch("duration", &tDuration, "duration/D");
    fits.Branch("pickupFlagTight", &tPickupTight, "pickupFlagTight/I");
    fits.Branch("passHarvest", &tPass, "passHarvest/I");
    fits.Branch("seedT0", &tSeedT0, "seedT0/D");
    fits.Branch("status_ln", &tStatusLN, "status_ln/I");
    fits.Branch("redChi2_ln", &tRedChi2LN, "redChi2_ln/D");
    fits.Branch("U0", &tU0, "U0/D");
    fits.Branch("tau", &tTau, "tau/D");
    fits.Branch("sigma", &tSigma, "sigma/D");
    fits.Branch("peakHeight", &tPeak, "peakHeight/D");
    fits.Branch("t0_ln", &tT0LN, "t0_ln/D");
    fits.Branch("ped_ln", &tPedLN, "ped_ln/D");
    fits.Branch("fitRiseNs_ln", &tFitRiseLN, "fitRiseNs_ln/D");
    fits.Branch("fitFallNs_ln", &tFitFallLN, "fitFallNs_ln/D");
    fits.Branch("pedestalRMS", &tPedRMS, "pedestalRMS/D");
    fits.Branch("pedestalSamples", &tPedSamples, "pedestalSamples/I");

    const Long64_t nEvents = inTree->GetEntries();
    Long64_t nNoPulse = 0, nFitted = 0, nPass = 0, nPassGood = 0;
    for (Long64_t i = 0; i < nEvents; ++i) {
        inTree->GetEntry(i);
        for (int ch = 0; ch < nSPEChan; ++ch) {
            const int d = ch / nChannelsPerDigitizer, c = ch % nChannelsPerDigitizer;
            TH1D wave(Form("w_%lld_%d", i, ch), "", nBins, 0, nBins * 2.5);
            for (int b = 0; b < nBins; ++b) wave.SetBinContent(b + 1, waveform[d][c][b]);

            // Offline view: dynamic pedestal subtracted, then pulse finding
            TH1D offWave(wave);
            applyDynamicPedestal(&offWave);
            OfflinePulse p = offlinePulse(&offWave);
            if (p.npulses == 0) { ++nNoPulse; continue; }

            FitResult r;
            const bool ok = fitPulse(&wave, true, p.time, preWindow, postWindow, false, 1, r,
                                     nullptr, Form("e%lld_c%d", i, ch), false);
            delete r.func;
            if (!ok) continue;
            ++nFitted;

            tEvent = i;
            tChan = ch;
            tCh8Type = (ch == 8) ? int(i % 2) : -1;
            tInjArea = injectedArea[d][c];
            tOffArea = p.area;
            tOffHeight = p.height;
            tNpulses = p.npulses;
            tRise = p.riseSamples;
            tFall = p.fallSamples;
            tDuration = p.duration;
            tPickupTight = p.pickupFlagTight;
            tPass = (p.npulses == 1 && p.riseSamples >= 3 && p.fallSamples >= 3 && !p.pickupFlagTight &&
                     p.area > winLow[ch] && p.area < winHigh[ch]);
            tSeedT0 = p.time;
            tStatusLN = r.status;
            tRedChi2LN = r.ndf > 0 ? r.chi2 / r.ndf : -1.0;
            tU0 = r.U0;
            tTau = r.tauLN;
            tSigma = r.sigma;
            tPeak = r.peakHeight;
            tT0LN = r.t0;
            tPedLN = r.ped;
            tFitRiseLN = r.fitRiseNs;
            tFitFallLN = r.fitFallNs;
            tPedRMS = r.pedestalRMS;
            tPedSamples = r.pedestalSamples;
            fits.Fill();
            if (tPass) {
                ++nPass;
                if (tStatusLN == 0 && tRedChi2LN >= 0.5 && tRedChi2LN <= 1.5) ++nPassGood;
            }
        }
        if (i % (nEvents / 10) == 0) std::cout << "Event " << i << " / " << nEvents << std::endl;
    }

    outFile->cd();
    fits.Write();
    outFile->Close();
    inFile->Close();
    std::cout << "channels with no offline pulse: " << nNoPulse << std::endl;
    std::cout << "fitted: " << nFitted << ", pass harvest: " << nPass
              << ", pass harvest && status_ln==0 && 0.5<=redChi2_ln<=1.5: " << nPassGood << std::endl;
    std::cout << "Wrote " << outputFile << std::endl;
}
