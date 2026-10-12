// Verbatim copy of lines 16-528 of milliqanOffline Run3Detector/src/runSpeFit.cc,
// branch SlabPulseInjectionSim at 2fb617aa (the fit that produced allFits_v4.root),
// minus the MilliDAQ include. Do not edit: refresh it from the source instead.
#pragma once
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <limits>
#include <fstream>
#include <sstream>
#include <iostream>
#include <string>
#include <vector>

#include "TCanvas.h"
#include "TF1.h"
#include "TFile.h"
#include "TFitResult.h"
#include "TFitResultPtr.h"
#include "TH1D.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TROOT.h"
#include "TString.h"
#include "TStyle.h"
#include "TTree.h"


namespace {

// Pulse finding defaults, matching configuration/pulseFinding/pulseFindingSlab.json
constexpr double kHighThreshold = 10.0;  // mV, starts a pulse
constexpr int kNConsecSamples = 3;       // samples above threshold to confirm
constexpr double kSidebandStart = 0.0;   // ns, noise estimation window
constexpr double kSidebandEnd = 50.0;    // ns

// Dynamic pedestal, matching OfflineFactory::prepareWave exactly.
constexpr int kPedestalTotalSamples = 400;        // first 1000 ns, avoids trigger
constexpr int kPedestalConsecutiveSamples = 16;   // block length
constexpr double kPedestalGranularity = 0.25;     // mV, MPV binning
constexpr double kRmsVariationMax = 4.0;          // mV, block must be this quiet
constexpr double kPedestalVariationMax = 80.0;    // mV, block mean must be this small

// Seed time constants, from an earlier by-hand fit of a slab SPE pulse.
constexpr double kSeedTauRise = 1.14;  // ns
constexpr double kSeedTauFall = 1.53;  // ns

char *getCmdOption(char **begin, char **end, const std::string &option) {
  char **itr = std::find(begin, end, option);
  if (itr != end && ++itr != end) return *itr;
  return 0;
}

bool cmdOptionExists(char **begin, char **end, const std::string &option) {
  return std::find(begin, end, option) != end;
}

// The model is unchanged:
//
//   f(t) = A*(exp(-(t-t0)/tau_fall) - exp(-(t-t0)/tau_rise)) + ped,  t > t0
//
// but it is fitted in terms of the pulse area Q and the gap dTau, where
// tau_fall = tau_rise + dTau and A = Q/dTau. Two reasons:
//
//  - As tau_fall approaches tau_rise the bracket collapses toward zero and A
//    diverges to compensate, so only the product A*dTau is determined. That
//    product is the area, so fitting Q directly removes the degeneracy. Fitting
//    A instead leaves it 45% uncertain even when the curve lies on the data.
//  - Fitting dTau > 0 guarantees tau_fall > tau_rise, instead of letting the
//    fit swap the two and invert the pulse.
//
// A and tau_fall are recovered afterwards, with errors propagated through the
// covariance matrix.
double speShape(double *x, double *p) {
  const double t = x[0];
  const double Q = p[0];
  const double t0 = p[1];
  const double tauRise = p[2];
  const double dTau = p[3];
  const double ped = p[4];
  if (t <= t0) return ped;
  const double tauFall = tauRise + dTau;
  if (dTau <= 0.0 || tauRise <= 0.0) return ped;
  return (Q / dTau) *
             (std::exp(-(t - t0) / tauFall) - std::exp(-(t - t0) / tauRise)) +
         ped;
}

// Log-normal pulse shape:
//
//   f(t) = U0 * exp(-0.5 * (ln((t-t0)/tau) / sigma)^2) + ped,   t > t0
//
// The peak is U0 and it occurs at t = t0 + tau, so amplitude, peak time and
// asymmetry are carried by three nearly independent parameters. Unlike the
// difference of exponentials, whose fall/rise 20-80% width ratio cannot go
// below 5.54, this shape covers any ratio: sigma 0.5 gives 3.4, matching the
// measured riseSamples/fallSamples of about 0.29.
double logNormalShape(double *x, double *p) {
  const double t = x[0];
  const double U0 = p[0];
  const double t0 = p[1];
  const double tau = p[2];
  const double sigma = p[3];
  const double ped = p[4];
  const double dt = t - t0;
  if (dt <= 1e-9 || tau <= 0.0 || sigma <= 0.0) return ped;
  const double z = std::log(dt / tau) / sigma;
  return U0 * std::exp(-0.5 * z * z) + ped;
}

// 20%-80% rise and 80%-20% fall widths of a fitted curve, in ns, so the fit can
// be compared directly against the pulse finder's riseSamples and fallSamples
// (which use the same levels) without assuming a functional form.
void measureFittedWidths(TF1 &fit, double lo, double hi, double ped,
                         double &riseNs, double &fallNs) {
  riseNs = fallNs = -1.0;
  const int nStep = 4000;
  const double step = (hi - lo) / nStep;
  double peak = 0.0, tPeak = lo;
  for (int i = 0; i <= nStep; ++i) {
    const double t = lo + i * step;
    const double v = fit.Eval(t) - ped;
    if (v > peak) {
      peak = v;
      tPeak = t;
    }
  }
  if (peak <= 0.0) return;
  auto cross = [&](double from, double to, double level) {
    const int n = 4000;
    const double h = (to - from) / n;
    for (int i = 0; i <= n; ++i) {
      const double t = from + i * h;
      if ((fit.Eval(t) - ped) >= level) return t;
    }
    return std::numeric_limits<double>::quiet_NaN();
  };
  const double t20r = cross(lo, tPeak, 0.2 * peak);
  const double t80r = cross(lo, tPeak, 0.8 * peak);
  const double t80f = cross(hi, tPeak, 0.8 * peak);
  const double t20f = cross(hi, tPeak, 0.2 * peak);
  // check if the cross values are valid.t80r == t80r is false if t80r is NaN.
  // This is safe gard that ensure the amplitude is positive.
  if (t80r == t80r && t20r == t20r) riseNs = t80r - t20r;
  if (t80f == t80f && t20f == t20f) fallNs = t80f - t20f;
}

// Value of the bracket at its maximum. Converts between the fitted A and the
// peak height of the pulse, which are not the same thing.
double shapePeakFactor(double tauRise, double tauFall) {
  if (tauFall <= tauRise) return 1.0;
  const double dPeak =
      (tauRise * tauFall / (tauFall - tauRise)) * std::log(tauFall / tauRise);
  return std::exp(-dPeak / tauFall) - std::exp(-dPeak / tauRise);
}

// Dynamic pedestal, a direct port of OfflineFactory::prepareWave. The waveform
// is scanned in blocks of 16 samples over the first 400 samples (1000 ns, before
// the trigger region). Blocks that are quiet enough (RMS < 4 mV) and not already
// sitting on signal (|mean| < 80 mV) contribute their mean, rounded to 0.25 mV,
// to a histogram. The pedestal is that histogram's most probable value, not the
// average: an average would be dragged up by any block containing signal.
// Returns the value subtracted.
double applyDynamicPedestal(TH1D *wave) {
  const int nBins =
      std::min(kPedestalTotalSamples, wave->GetNbinsX());
  TH1D histTemp("temp", "temp",
                1 + int(kPedestalVariationMax / kPedestalGranularity + 1e-3) * 2,
                -kPedestalVariationMax - kPedestalGranularity / 2,
                kPedestalVariationMax + kPedestalGranularity / 2);
  for (int ibin = 1; ibin + kPedestalConsecutiveSamples - 1 <= nBins;
       ibin += kPedestalConsecutiveSamples) {
    double sum = 0.0;
    for (int jbin = ibin; jbin < ibin + kPedestalConsecutiveSamples; ++jbin) {
      sum += wave->GetBinContent(jbin);
    }
    const double mean = sum / kPedestalConsecutiveSamples;
    double var = 0.0;
    for (int jbin = ibin; jbin < ibin + kPedestalConsecutiveSamples; ++jbin) {
      var += std::pow(wave->GetBinContent(jbin) - mean, 2.0);
    }
    const double rms = std::fabs(std::sqrt(var / kPedestalConsecutiveSamples));
    if (std::fabs(mean) < kPedestalVariationMax && rms < kRmsVariationMax) {
      histTemp.Fill(std::round(mean / kPedestalGranularity) *
                    kPedestalGranularity);
    }
  }

  double pedestal = 0.0;
  if (histTemp.GetEntries()) {
    pedestal = histTemp.GetBinCenter(histTemp.GetMaximumBin());
  }
  for (int bin = 1; bin <= wave->GetNbinsX(); ++bin) {
    wave->SetBinContent(bin, wave->GetBinContent(bin) - pedestal);
  }
  return pedestal;
}

// Mean and RMS of the samples in [lo, hi], returning the sample count. Run
// after the pedestal subtraction, so the mean should sit near zero.
int measureRange(TH1D *wave, double lo, double hi, double &mean, double &rms) {
  const int firstBin = wave->FindBin(lo);
  const int lastBin = wave->FindBin(hi);
  double sum = 0.0;
  double sum2 = 0.0;
  int n = 0;
  for (int bin = firstBin; bin <= lastBin; ++bin) {
    const double v = wave->GetBinContent(bin);
    sum += v;
    sum2 += v * v;
    ++n;
  }
  if (n == 0) {
    mean = 0.0;
    rms = 1.0;
    return 0;
  }
  mean = sum / n;
  const double variance = sum2 / n - mean * mean;
  rms = variance > 0.0 ? std::sqrt(variance) : 1.0;   // notice that although it is called rms, it is actually the standard deviation. I just follow the convention from the OfflineFactory.
  if (rms <= 0.0) rms = 1.0;
  return n;
}

// The 0-50 ns sideband, matching OfflineFactory::measureSideband. Recorded for
// comparison; it is no longer the fit uncertainty, see fitPulse.
void measureSideband(TH1D *wave, double &mean, double &rms) {
  measureRange(wave, kSidebandStart, kSidebandEnd, mean, rms);
}

// The pedestal region runs from the start of the record to the trigger window,
// which the stage-1 selection guarantees contains no pulse, stopping a guard
// short of the pulse in case one starts early.
constexpr double kPreTriggerEndNs = 1000.0;
constexpr double kPedGuardNs = 50.0;
constexpr int kMinPedSamples = 64;

// First sample where the waveform stays above threshold for kNConsecSamples.
// This approximates OfflineFactory::findPulses well enough to locate the
// ipulse==0 pulse, which is the one the selection kept.
int findFirstPulseBin(TH1D *wave, double baseline, double threshold) {
  const int nBins = wave->GetNbinsX();
  for (int bin = 1; bin <= nBins - kNConsecSamples; ++bin) {
    bool allAbove = true;
    for (int k = 0; k < kNConsecSamples; ++k) {
      if (wave->GetBinContent(bin + k) - baseline < threshold) {
        allAbove = false;
        break;
      }
    }
    if (allAbove) return bin;
  }
  return -1;
}

// Rough time constants straight from the data: the rise from 10% to 90% of the
// peak, and the decay from the peak down to 1/e of it. Both are crude, but they
// put the fit in the right decade, which fixed seeds do not.
void estimateTaus(TH1D *wave, double baseline, double fitLow, double fitHigh,
                  double &tauRise, double &tauFall) {
  const int lowBin = wave->FindBin(fitLow);
  const int highBin = wave->FindBin(fitHigh);
  int peakBin = lowBin;
  double peak = 0.0;
  for (int bin = lowBin; bin <= highBin; ++bin) {
    const double v = wave->GetBinContent(bin) - baseline;
    if (v > peak) {
      peak = v;
      peakBin = bin;
    }
  }
  if (peak <= 0.0) return;

  int bin10 = peakBin;
  int bin90 = peakBin;
  for (int bin = peakBin; bin >= lowBin; --bin) {
    const double v = wave->GetBinContent(bin) - baseline;
    if (v > 0.9 * peak) bin90 = bin;
    if (v > 0.1 * peak) bin10 = bin;
    else break;
  }
  const double riseTime =
      wave->GetBinLowEdge(bin90) - wave->GetBinLowEdge(bin10);
  // 10-90% of an exponential rise spans about 2.2 time constants.
  if (riseTime > 0.0) tauRise = std::max(0.1, riseTime / 2.2);

  int binE = peakBin;
  for (int bin = peakBin; bin <= highBin; ++bin) {
    binE = bin;
    if (wave->GetBinContent(bin) - baseline < peak / std::exp(1.0)) break;
  }
  const double fallTime =
      wave->GetBinLowEdge(binE) - wave->GetBinLowEdge(peakBin);
  if (fallTime > 0.0) tauFall = std::max(tauRise * 1.5, fallTime);
  if (tauFall <= tauRise) tauFall = tauRise * 3.0;
}

}  // namespace

// One pulse's worth of fit output.
struct FitResult {
  double Q = 0, errQ = 0, A = 0, errA = 0;
  double t0 = 0, errT0 = 0;
  double tauRise = 0, errTauRise = 0, tauFall = 0, errTauFall = 0;
  double dTau = 0, errDTau = 0, ped = 0, errPed = 0;
  double U0 = 0, errU0 = 0, tauLN = 0, errTauLN = 0, sigma = 0, errSigma = 0;
  double chi2 = 0, peakHeight = 0, baselineRMS = 0, dynPedestal = 0;
  // Noise over the whole pedestal region, which is what the fit uses.
  // baselineRMS above is the 0-50 ns sideband, kept so the two can be compared.
  double pedestalMean = 0, pedestalRMS = 0;
  int pedestalSamples = 0;
  double fitRiseNs = -1, fitFallNs = -1;
  TF1 *func = nullptr;  // owned by the caller, for overlay drawing
  int ndf = 0, status = -1, model = 0;
  bool ok = false;
};

// Pedestal-correct, locate the pulse, fit, and fill result. When canvas is not
// null the waveform and fitted curve are drawn into it.
// model 0 = difference of exponentials, 1 = log-normal.
bool fitPulse(TH1D *wave, bool haveT0, double userT0, double preWindow,
              double postWindow, bool isSlab, int model, FitResult &r,
              TCanvas *canvas, const char *tag, bool verbose) {
  r.model = model;
  if (isSlab) wave->Scale(-1);
  r.dynPedestal = applyDynamicPedestal(wave);

  double baseline = 0.0, baselineRMS = 1.0;
  measureSideband(wave, baseline, baselineRMS);
  r.baselineRMS = baselineRMS;

  double pulseStart = 0.0;
  if (haveT0) {
    pulseStart = userT0;
  } else {
    const int pulseBin = findFirstPulseBin(wave, baseline, kHighThreshold);
    if (pulseBin < 0) return false;
    pulseStart = wave->GetBinLowEdge(pulseBin);
  }

  // Clamp to the record, so pulses near either end still fit over whatever is
  // available instead of running off the histogram.
  const double axisLow = wave->GetXaxis()->GetXmin();
  const double axisHigh = wave->GetXaxis()->GetXmax();
  const double fitLow = std::max(axisLow, pulseStart - preWindow);
  const double fitHigh = std::min(axisHigh, pulseStart + postWindow);
  if (fitHigh - fitLow < 20.0) return false;

  // The per-bin uncertainty is measured over the pedestal region, not the
  // 0-50 ns sideband.
  //
  // The sideband is twenty samples spanning 50 ns, while the fit spans 500 ns.
  // Noise slower than about 20 MHz cannot develop within 50 ns, so the
  // sideband is structurally blind to it, and that is exactly where the
  // difference between pulses this model fits and pulses it does not turns out
  // to lie: measured over 100 of each, the failing ones carry 1.30x the noise
  // power at 1-5 MHz and 1.24x at 5-20 MHz, but only 1.02x above 20 MHz. An
  // uncertainty measured over 50 ns therefore describes the wrong thing, which
  // is why redChi2 showed no dependence on baselineRMS at all.
  //
  // This is not a tuned choice. The region is every sample before the trigger
  // window, which the stage-1 selection guarantees is pulse free, and the rule
  // is to use all of it rather than some fraction picked for its effect on
  // chi2. The sideband stays in baselineRMS so the two can be compared.
  double pedMean = baseline, pedRMS = baselineRMS;
  const double pedEnd = std::min(kPreTriggerEndNs, pulseStart - kPedGuardNs);
  r.pedestalSamples = measureRange(wave, axisLow, pedEnd, pedMean, pedRMS);
  if (r.pedestalSamples < kMinPedSamples) {
    // Too little pedestal to measure: fall back rather than fit against a
    // meaningless uncertainty.
    pedMean = baseline;
    pedRMS = baselineRMS;
  }
  r.pedestalMean = pedMean;
  r.pedestalRMS = pedRMS;
  // From here these are what the fit sees: the bin errors below, the ped
  // bound, the parameter seeds and peakHeight.
  baseline = pedMean;
  baselineRMS = pedRMS;

  for (int bin = 1; bin <= wave->GetNbinsX(); ++bin) {
    wave->SetBinError(bin, baselineRMS);
  }

  double peakHeight = 0.0;
  const int lowBin = wave->FindBin(fitLow);
  const int highBin = wave->FindBin(fitHigh);
  for (int bin = lowBin; bin <= highBin; ++bin) {
    peakHeight = std::max(peakHeight, wave->GetBinContent(bin) - baseline);
  }
  r.peakHeight = peakHeight;
  if (peakHeight <= 0.0) return false;

  TF1 *fitp =
      new TF1(Form("%s_%s", model == 1 ? "speFitLN" : "speFit", tag),
              model == 1 ? logNormalShape : speShape, fitLow, fitHigh, 5);
  TF1 &fit = *fitp;
  r.func = fitp;

  if (model == 1) {
    // Peak time seeds tau directly, since the log-normal peaks at t0 + tau.
    double tPeak = pulseStart;
    double best = 0.0;
    for (int bin = lowBin; bin <= highBin; ++bin) {
      const double v = wave->GetBinContent(bin) - baseline;
      if (v > best) {
        best = v;
        tPeak = wave->GetBinCenter(bin);
      }
    }
    const double seedTau = std::max(1.0, tPeak - pulseStart);
    fit.SetParNames("U0", "t0", "tau", "sigma", "ped");
    // sigma 0.5 reproduces the measured riseSamples/fallSamples of ~0.29.
    fit.SetParameters(peakHeight, pulseStart, seedTau, 0.5, baseline);
    fit.SetParLimits(0, 0.0, 100.0 * std::max(1.0, peakHeight));
    fit.SetParLimits(1, fitLow, fitHigh);
    fit.SetParLimits(2, 0.5, 500.0);
    fit.SetParLimits(3, 0.05, 3.0);
  } else {
    double seedTauRise = kSeedTauRise, seedTauFall = kSeedTauFall;
    estimateTaus(wave, baseline, fitLow, fitHigh, seedTauRise, seedTauFall);
    const double seedFactor = shapePeakFactor(seedTauRise, seedTauFall);
    const double seedA = seedFactor > 0.0 ? peakHeight / seedFactor : peakHeight;
    const double seedDTau = std::max(0.1, seedTauFall - seedTauRise);
    fit.SetParNames("Q", "t0", "tau_rise", "dTau", "ped");
    fit.SetParameters(seedA * seedDTau, pulseStart, seedTauRise, seedDTau,
                      baseline);
    fit.SetParLimits(0, 0.0, 1000.0 * std::max(1.0, peakHeight) * seedDTau);
    fit.SetParLimits(1, fitLow, fitHigh);
    fit.SetParLimits(2, 0.05, 200.0);
    fit.SetParLimits(3, 0.01, 500.0);
  }
  fit.SetParLimits(4, baseline - 5.0 * baselineRMS,
                   baseline + 5.0 * baselineRMS);
  fit.SetLineColor(kRed);
  fit.SetLineStyle(2);

  TFitResultPtr res = wave->Fit(&fit, verbose ? "RS" : "RQS");
  r.status = int(res);

  r.t0 = fit.GetParameter(1);  r.errT0 = fit.GetParError(1);
  r.ped = fit.GetParameter(4); r.errPed = fit.GetParError(4);
  r.chi2 = fit.GetChisquare();
  r.ndf = fit.GetNDF();
  measureFittedWidths(fit, fitLow, fitHigh, r.ped, r.fitRiseNs, r.fitFallNs);

  if (model == 1) {
    r.U0 = fit.GetParameter(0);    r.errU0 = fit.GetParError(0);
    r.tauLN = fit.GetParameter(2); r.errTauLN = fit.GetParError(2);
    r.sigma = fit.GetParameter(3); r.errSigma = fit.GetParError(3);
  } else {
    r.Q = fit.GetParameter(0);       r.errQ = fit.GetParError(0);
    r.tauRise = fit.GetParameter(2); r.errTauRise = fit.GetParError(2);
    r.dTau = fit.GetParameter(3);    r.errDTau = fit.GetParError(3);

    // Recover the parameters of the A form. Q and dTau are strongly
    // correlated, so the covariance term is not optional. This is error
    // propagation, not a correlation measurement: the shape parameters are
    // compared across channels by plotting them, not through Minuit.
    double covQD = 0.0, covRD = 0.0;
    if (res.Get()) {
      covQD = res->CovMatrix(0, 3);
      covRD = res->CovMatrix(2, 3);
    }
    r.tauFall = r.tauRise + r.dTau;
    r.errTauFall = std::sqrt(std::max(
        0.0, r.errTauRise * r.errTauRise + r.errDTau * r.errDTau + 2.0 * covRD));
    if (r.dTau > 0.0) {
      r.A = r.Q / r.dTau;
      const double dAdQ = 1.0 / r.dTau;
      const double dAdD = -r.Q / (r.dTau * r.dTau);
      r.errA = std::sqrt(std::max(0.0, dAdQ * dAdQ * r.errQ * r.errQ +
                                           dAdD * dAdD * r.errDTau * r.errDTau +
                                           2.0 * dAdQ * dAdD * covQD));
    }
  }
  r.ok = true;

  if (canvas) {
    canvas->cd();
    canvas->SetGrid();
    TH1D *draw = (TH1D *)wave->Clone(Form("h_%s", tag));
    draw->GetXaxis()->SetRangeUser(fitLow - 20.0, fitHigh + 20.0);
    draw->SetTitle(Form("%s;Time (ns);Amplitude (mV)", tag));
    draw->SetLineColor(kBlack);
    draw->SetDirectory(nullptr);
    draw->Draw("HIST");
    fit.DrawCopy("SAME");
    TLatex latex;
    latex.SetNDC();
    latex.SetTextSize(0.032);
    if (model == 1) {
      latex.DrawLatex(0.58, 0.86, Form("U_{0} = %.2f #pm %.2f mV", r.U0, r.errU0));
      latex.DrawLatex(0.58, 0.82,
                      Form("#tau = %.2f #pm %.2f ns", r.tauLN, r.errTauLN));
      latex.DrawLatex(0.58, 0.78,
                      Form("#sigma = %.3f #pm %.3f", r.sigma, r.errSigma));
      latex.DrawLatex(0.58, 0.74,
                      Form("rise/fall = %.3f",
                           r.fitFallNs > 0 ? r.fitRiseNs / r.fitFallNs : 0.0));
    } else {
      latex.DrawLatex(0.58, 0.86, Form("A = %.2f #pm %.2f", r.A, r.errA));
      latex.DrawLatex(0.58, 0.82, Form("Q = %.1f #pm %.1f pVs", r.Q, r.errQ));
      latex.DrawLatex(0.58, 0.78,
                      Form("#tau_{rise} = %.2f #pm %.2f ns", r.tauRise, r.errTauRise));
      latex.DrawLatex(0.58, 0.74,
                      Form("#tau_{fall} = %.2f #pm %.2f ns", r.tauFall, r.errTauFall));
    }
    latex.DrawLatex(0.58, 0.70, Form("t_{0} = %.2f #pm %.2f ns", r.t0, r.errT0));
    latex.DrawLatex(0.58, 0.66, Form("ped = %.3f #pm %.3f", r.ped, r.errPed));
    latex.DrawLatex(0.58, 0.62,
                    Form("#chi^{2}/ndf = %.2f (status %d)",
                         r.ndf > 0 ? r.chi2 / r.ndf : 0.0, r.status));
  }
  return true;
}
