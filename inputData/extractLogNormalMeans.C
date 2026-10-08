// Per-channel mean log-normal shape (tau, sigma) for the SPE template used by
// waveinject_slab.C. Reads the SPE fit tree (speFits) and writes a text file in
// the same spirit as SPEmeans.txt.
//
// Selection: status_ln == 0 && minRedChi2 <= redChi2_ln <= maxRedChi2.
// Channel 8 has two tau populations (two PMTs on one summing amp), so it gets
// three rows: "8" over all its pulses, "8A" for tau < ch8TauSplit (fast) and
// "8B" for tau >= ch8TauSplit (slow).
//
//   root -l -b -q 'extractLogNormalMeans.C("allFits_v4.root")'
//   root -l -b -q 'extractLogNormalMeans.C("allFits_v4.root","SPElognormal.txt",0.5,1.5,50)'
#include "TFile.h"
#include "TTree.h"

#include <cmath>
#include <cstdio>
#include <iostream>
#include <string>
#include <vector>

namespace {
struct Acc {
    double n = 0, sTau = 0, sTau2 = 0, sSig = 0, sSig2 = 0;
    void add(double tau, double sigma) {
        n += 1;
        sTau += tau;
        sTau2 += tau * tau;
        sSig += sigma;
        sSig2 += sigma * sigma;
    }
    // mean and error on the mean
    void get(double s, double s2, double& mean, double& err) const {
        mean = s / n;
        double var = s2 / n - mean * mean;
        err = (n > 1 && var > 0) ? std::sqrt(var / (n - 1)) : 0;
    }
};

void writeRow(FILE* out, const std::string& label, const Acc& a) {
    if (a.n == 0) {
        std::cerr << "WARNING: no selected fits for chan " << label << ", row skipped" << std::endl;
        return;
    }
    double tau, tauErr, sigma, sigmaErr;
    a.get(a.sTau, a.sTau2, tau, tauErr);
    a.get(a.sSig, a.sSig2, sigma, sigmaErr);
    std::fprintf(out, "%s %.4f %.4f %.5f %.5f %.0f\n", label.c_str(), tau, tauErr, sigma, sigmaErr, a.n);
}
}  // namespace

int extractLogNormalMeans(const char* fitFile,
                          const char* outFile = "SPElognormal.txt",
                          double minRedChi2 = 0.5,
                          double maxRedChi2 = 1.5,
                          double ch8TauSplit = 50.0) {
    const int nChan = 96;

    TFile* f = TFile::Open(fitFile, "READ");
    if (!f || f->IsZombie()) {
        std::cerr << "ERROR: could not open " << fitFile << std::endl;
        return 1;
    }
    auto* t = dynamic_cast<TTree*>(f->Get("speFits"));
    if (!t) {
        std::cerr << "ERROR: no speFits tree in " << fitFile << std::endl;
        return 2;
    }

    // tau and sigma are always the log-normal's. In a --model both file the
    // generic status/redChi2 belong to the difference of exponentials, so use
    // the log-normal's own pair whenever the file has it.
    Int_t chan = 0, status = 0;
    Double_t redChi2 = 0, tau = 0, sigma = 0;
    const bool haveLN = t->GetBranch("status_ln") != nullptr;
    t->SetBranchStatus("*", 0);
    for (const char* b : {"chan", "tau", "sigma"}) t->SetBranchStatus(b, 1);
    t->SetBranchAddress("chan", &chan);
    t->SetBranchAddress("tau", &tau);
    t->SetBranchAddress("sigma", &sigma);
    const char* statusBr = haveLN ? "status_ln" : "status";
    const char* chi2Br = haveLN ? "redChi2_ln" : "redChi2";
    t->SetBranchStatus(statusBr, 1);
    t->SetBranchStatus(chi2Br, 1);
    t->SetBranchAddress(statusBr, &status);
    t->SetBranchAddress(chi2Br, &redChi2);
    std::cout << "Selecting " << statusBr << " == 0 && " << minRedChi2 << " <= " << chi2Br
              << " <= " << maxRedChi2 << std::endl;

    std::vector<Acc> acc(nChan);
    Acc ch8A, ch8B;
    const Long64_t n = t->GetEntries();
    for (Long64_t i = 0; i < n; ++i) {
        t->GetEntry(i);
        if (status != 0 || redChi2 < minRedChi2 || redChi2 > maxRedChi2) continue;
        if (chan < 0 || chan >= nChan) continue;
        acc[chan].add(tau, sigma);
        if (chan == 8) (tau < ch8TauSplit ? ch8A : ch8B).add(tau, sigma);
    }

    FILE* out = std::fopen(outFile, "w");
    if (!out) {
        std::cerr << "ERROR: could not write " << outFile << std::endl;
        return 3;
    }
    std::fprintf(out, "chan tauMean tauMeanErr sigmaMean sigmaMeanErr nFits\n");
    for (int c = 0; c < nChan; ++c) {
        writeRow(out, std::to_string(c), acc[c]);
        if (c == 8) {
            writeRow(out, "8A", ch8A);
            writeRow(out, "8B", ch8B);
        }
    }
    std::fclose(out);
    std::cout << "Wrote " << outFile << " (8A: tau < " << ch8TauSplit << " ns, 8B: tau >= "
              << ch8TauSplit << " ns)" << std::endl;
    f->Close();
    return 0;
}
