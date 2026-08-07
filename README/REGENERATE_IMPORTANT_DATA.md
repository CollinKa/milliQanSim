# Regenerate important slab SPE / wave-injection data

Large `.root` outputs (merged SPE pulses, 1M-event waveinjection cases, cosmic
waveinjected samples, plot dumps) should **not** be pushed to GitHub. Keep the
small inputs + scripts in git (or locally), and regenerate the heavy products
with the steps below.

Related setup notes for this Mac / Geant4 11 branch:
[`README_mac_geant4_11.md`](../README_mac_geant4_11.md).

Flattening Geant4 sim → lightweight tree:
[`flatlightwithphotonslab_latest.py`](./flatlightwithphotonslab_latest.py).

---

## 0. Environment

```bash
cd /Users/haoliangzheng/Desktop/CERN/milliQanSim_Geant4P11/milliQanSim
source ./setup_mac_g4p11.sh
# or: conda activate milliqan-sim-g4p11
```

Confirm ROOT is available (`root -l -b -q`). Wave-injection scripts load
`../build/libMilliQanCore.dylib`, so build the sim first if that library is
missing:

```bash
bash beamSetupSlab.sh
```

---

## 1. Small inputs you need (do not regenerate from nothing)

| File | Role |
|------|------|
| `inputData/modified_waveform.root` | SPE waveform template histogram (read-only) |
| `inputData/SPEmeans.txt` | Per-channel SPE area Gaussian mean/sigma **[nVs]** |
| `inputData/SPE_hists.root` | Offline SPE area/height hists to compare against |
| `inputData/offline_pulse_area_chan3.txt` | Optional Chan3 offline area overlay (pVs) |

`SPEmeans.txt` format:

```text
chan SPEmean SPEmeanErr SPEsigma SPEsigmaErr
0 1.726 0.01 0.863 0.019
...
```

If you update offline SPE fits, copy a fresh `SPEmeans.txt` here before
re-running wave injection / per-channel studies.

---

## 2. Offline SPE fits → `SPEmeans.txt` / fit ROOT summaries

Scripts:

- `inputData/OfflineSlabSPEFit.C`
- `inputData/OfflineSlabSPEFit_SPEplotter.C`

These expect offline selected-pulse ROOT files under a local directory
(default in the script: `~/Downloads/SlabSPEStudy_merged_pulses_spe/`).

```bash
cd inputData
# smoke: one run
root -l -b -q 'OfflineSlabSPEFit.C("smoke")'
root -l -b -q 'OfflineSlabSPEFit_SPEplotter.C("smoke")'

# full: hadd all runs then fit (produces large hadd ROOT — keep local)
root -l -b -q 'OfflineSlabSPEFit.C("all")'
root -l -b -q 'OfflineSlabSPEFit_SPEplotter.C("all")'
```

Typical local outputs (do **not** push the large hadd):

- `OfflineSlabSPEFit_smoke_Run1605.root` + `*_fitSummary.csv`
- `OfflineSlabSPEFit_all.root` + `*_fitSummary.csv`
- `SlabSPEStudy_merged_pulses_spe_hadd.root` (large)

Update `SPEmeans.txt` from the fit summary before using it in wave injection.

---

## 3. SPE area template study (single Gaussian)

Script: `inputData/SlabSPEareaStudy.C`

Uses `modified_waveform.root`, samples SPE area from a Gaussian, scales the
template, writes area/height, fits area-vs-height, and can overlay offline Chan3.

```bash
cd inputData
root -l -b -q 'SlabSPEareaStudy.C+'

# refit / redraw only (no resampling)
root -l -b -q -e 'gROOT->ProcessLine(".L SlabSPEareaStudy.C+"); fitSlabSPEarea();'
root -l -b -q -e 'gROOT->ProcessLine(".L SlabSPEareaStudy.C+"); compareAreaDistributions();'
```

Outputs:

- `SlabSPEarea.root`
- `SlabSPEarea_areaVsHeight.pdf` / `.png`
- `SlabSPEarea_offline_vs_sim.pdf` / `.png`

---

## 4. Per-channel SPE area / height study

Script: `inputData/SlabSPEareaStudy_perChannel.C`

Reads `SPEmeans.txt`, samples per channel, compares to `SPE_hists.root`.

```bash
cd inputData
# default: offline-style pulse window, all chans 0..95
root -l -b -q 'SlabSPEareaStudy_perChannel.C+'

# full-waveform integral
root -l -b -q 'SlabSPEareaStudy_perChannel.C+(20000,0,95,"full")'

# smoke: one channel
root -l -b -q 'SlabSPEareaStudy_perChannel.C+(5000,3,3,"full")'

# redraw comparisons only
root -l -b -q -e 'gROOT->ProcessLine(".L SlabSPEareaStudy_perChannel.C+"); comparePerChannelToSPEHists();'
root -l -b -q -e 'gROOT->ProcessLine(".L SlabSPEareaStudy_perChannel.C+"); comparePerChannelHeight();'
root -l -b -q -e 'gROOT->ProcessLine(".L SlabSPEareaStudy_perChannel.C+"); comparePerChannelHeightVsArea();'
```

Outputs (local / optional in git if small):

- `SlabSPEarea_perChannel.root` / `SlabSPEarea_perChannel_fullIntegral.root`
- CSVs and `*_plots/` directories

---

## 5. Geant4 slab sim → flat tree → wave injection

### 5a. Build / run slab beam muons

```bash
bash beamSetupSlab.sh
# smoke output: build/beamMuonSlab_MilliQan.root
```

For larger samples, edit `runMac/mcp_novis.mac` (or your production mac) and run:

```bash
cd build
./MilliQan ../runMac/mcp_novis.mac
```

### 5b. Flatten to lightweight tree

Edit input/output paths at the bottom of
`README/flatlightwithphotonslab_latest.py`, then:

```bash
python README/flatlightwithphotonslab_latest.py
```

On macOS this loads `build/libMilliQanCore.dylib` (see comments in the script).

### 5c. Inject SPE waveforms (production path)

Script: `inputData/waveinject_slab.C`

Uses `modified_waveform.root` + per-channel Gaussians from `SPEmeans.txt`.

```bash
cd inputData
root -l -b -q 'waveinject_slab.C+("../build/beamMuonSlab_1kEvent.root","beamMuonSlab_1kEvent_waveinjected.root","modified_waveform.root","SPEmeans.txt")'
```

### 5d. Extract pulse area / height (offline-style)

Script: `inputData/extractWaveformAreaHeight.C`

```bash
cd inputData
root -l -b -q 'extractWaveformAreaHeight.C+("beamMuonSlab_1kEvent_waveinjected.root","beamMuonSlab_1kEvent_waveinjected_areaHeight.root")'
# or cosmic example defaults:
root -l -b -q 'extractWaveformAreaHeight.C+("cosmicSlabCorr.root","cosmicSlabCorr_areaHeight.root")'
```

---

## 6. Large pulse-injection study cases (local only)

Script: `inputData/waveinject_slab_InjectionStudy.C`

Produces ~400 MB files under `inputData/pulseInjectionStudy/`
(`1MCd109center_waveinjectedCase*.root`). **Do not push these.**

Point `inputFile` / `outputFile` in the script (or pass arguments) at your
1M-event Geant4 file, then:

```bash
cd inputData
mkdir -p pulseInjectionStudy
root -l -b -q 'waveinject_slab_InjectionStudy.C+'
```

Same rule for other multi-hundred-MB products:

- `inputData/SlabSPEStudy_merged_pulses_hadd.root`
- `inputData/SlabSPEStudy_merged_pulses_spe_hadd.root`
- `rawCosmicSim/MilliQan.root`

Store them on disk / CERNBox / EOS, not in git.

---

## 7. Suggested regenerate order

1. Update offline SPE fits → refresh `SPEmeans.txt` (and optionally `SPE_hists.root`).
2. `SlabSPEareaStudy.C` / `SlabSPEareaStudy_perChannel.C` to validate template + per-channel SPE.
3. Run / obtain Geant4 slab `MilliQan.root`.
4. Flatten with `flatlightwithphotonslab_latest.py` if you need the flat tree.
5. `waveinject_slab.C` (+ `extractWaveformAreaHeight.C`) for analysis samples.
6. Keep large ROOT products out of git; document paths in notes or a local config.

---

## 8. What belongs in git vs local

**Useful to keep / push (small):**

- This README and `README_mac_geant4_11.md`
- `README/flatlightwithphotonslab_latest.py`
- Analysis scripts under `inputData/*.C` (optional)
- `modified_waveform.root`, `SPEmeans.txt`, `SPE_hists.root`, small CSVs/txt

**Keep local / regenerate (large):**

- Any `*_hadd.root`, `pulseInjectionStudy/*.root`, multi-MB cosmic/beam waveinjected files
- Plot galleries under `inputData/*_plots/`
- ACLiC build leftovers: `*_C.so`, `*_C.d`, `*_rdict.pcm`
