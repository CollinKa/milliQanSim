# Per-channel log-normal SPE template for slab pulse injection

Branch `SlabPulseInjectLogNormalSPE`. `waveinject_slab.C` no longer uses the
single measured `average_waveform` from `modified_waveform.root`. Each summing
channel gets its own SPE shape, built from the per-pulse log-normal fits of the
slab SPE study.

## Files

All files are in `inputData/`.

| File | What it is |
|---|---|
| `extractLogNormalMeans.C` | Reads the SPE fit tree (`speFits`) and writes the mean `tau` and `sigma` per channel to `SPElognormal.txt`. |
| `SPElognormal.txt` | The per-channel shapes: `chan tauMean tauMeanErr sigmaMean sigmaMeanErr nFits`. Rows `0`–`95`, plus `8A` and `8B`. |
| `waveinject_slab.C` | Pulse injection. Builds one template per channel from `SPElognormal.txt`. |
| `SPEmeans.txt` | Per-channel SPE **area** Gaussians (nVs). Unchanged, still used to sample the pulse area. |

The fit tree is `allFits_v4.root`: 2.40M pulses, produced on lxplus by
`speFit.exe --model both`. The local copy is
`~/Documents/Working on lxplux/spe_results/allFits_v4.root`, and the EOS copy is
`/eos/user/h/hazheng/SlabSPEFits/allFits_v4.root`. How the fits were made is
documented in `~/Documents/Working on lxplux/condorTest/SPE_STUDY.md`.

## How to run

From `inputData/`, with the env from `setup_mac_g4p11.sh` sourced:

```bash
# 1. per-channel means -> SPElognormal.txt (about 10 s)
root -l -b -q 'extractLogNormalMeans.C+("/Users/haoliangzheng/Documents/Working on lxplux/spe_results/allFits_v4.root")'

# 2. pulse injection (about 1 min for 1k beam-muon events)
root -l -b -q 'waveinject_slab.C("../build/beamMuonSlab_1kEvent.root","beamMuonSlab_1kEvent_waveinjected.root")'
```

`extractLogNormalMeans.C` takes these arguments:
`(fitFile, outFile="SPElognormal.txt", minRedChi2=0.5, maxRedChi2=1.5, ch8TauSplit=50)`.

`waveinject_slab.C` takes
`(inputFile, outputFile, lognormalFile="SPElognormal.txt", speMeansFile="SPEmeans.txt")`.
The third argument used to be the waveform template ROOT file. Run it
**interpreted**, without the `+`, because ACLiC cannot link `libMilliQanCore.dylib`
through `R__LOAD_LIBRARY`.

## Decisions

1. **Shape.** The template is
   `f(t) = U0 * exp(-0.5 * (ln((t-t0)/tau)/sigma)^2) + ped`, with
   **U0 = 50 mV, ped = 0** and **t0 = 1285 ns**. Only `tau` and `sigma` change
   per channel. The histogram uses the digitizer binning: 1024 bins of 2.5 ns,
   0–2560 ns, with the function evaluated at each bin centre.

2. **Why t0 = 1285 ns.** That is where the old `average_waveform` starts its
   rise. Before it the template is a slow foot under 4% of peak. At 1285 ns it
   crosses 5% of its peak (32.9 at bin 529, 1320 ns) and climbs steeply
   (bin 515 is 5.1%, bin 520 is 35%). Putting t0 there keeps the hit-time shift
   (`initial_hit_time / binWidth`) aligning pulses exactly as before. The
   log-normal peak now lands at t0 + tau: about 1323 ns for the fastest channel
   and 1381 ns for channel 92.

3. **Per-channel mean.** The mean of `tau` and of `sigma` over fits with
   `status_ln == 0 && 0.5 <= redChi2_ln <= 1.5`, which is 1.77M of 2.40M
   pulses. The log-normal's own status and chi2 are used, not the generic
   `status`/`redChi2`, which belong to the difference-of-exponentials model in
   a `--model both` file.

4. **Channel 8 has two PMTs with different tau.** Every channel is a summing amp
   of two PMTs, and only channel 8 shows two tau peaks: 7568 fits at 37.9 ns
   (sigma 0.704) and 9298 at 63.7 ns (sigma 0.577). The split at 50 ns sits in
   the empty valley between them (49–55 ns). The text file has three rows:
   - `8A`: tau < 50, fast
   - `8B`: tau >= 50, slow
   - `8`: both together, tau 52.2, sigma 0.634

5. **Which channel-8 PMT is which.** Channel 8 is fed by sim PMTs 34 and 35, the
   front pair of data slab 4. **PMT 34 gets 8A (fast) and PMT 35 gets 8B
   (slow).** The data cannot tell which physical PMT is the fast one, so this
   assignment is a choice.

6. **Channel 8 above 5000 hits uses row `8`.** In that branch all hits of a PMT
   are summed into one pulse, and both PMTs are equally likely to see the
   photons, so the mean shape is used. Pulses this large saturate anyway
   (above 1250 mV, and longer than the 2560 ns record).

7. **Where the template is chosen.** Right after the channel check
   (`if (remappedPMT < 0 || ... || !fit[remappedPMT]) continue;`), a local
   `TH1F* pulse_shape` is set, so both existing
   `TH1F* new_waveform = (TH1F*)pulse_shape->Clone();` lines work unchanged.
   This one spot covers both the >5000-hit and the per-hit paths.
   `pmtHitsMap` is keyed by the **sim** PMT number, which is one physical PMT.
   That is why channel 8 can be chosen per PMT there, without changing the
   grouping loop.

8. **Channel 92 is treated as high gain.** Its `SPEmeans.txt` area
   (14.4 nVs, about 4× other channels) is kept as is, and it gets its own
   template (tau 95.3, sigma 0.451).

9. **Scaling to data: `scaleFactor = area / A`.** This replaces
   `area / 828.03 * (SPE_height / 59.0878)`. `area` is sampled from the
   channel's `SPEmeans.txt` Gaussian, in pVs: one sample per photoelectron, or
   their sum above 5000 hits. `A` is that channel's template area in pVs over
   the whole record, `Integral(1, nBins) * 2.5 ns` on the U0 = 50 mV, ped = 0,
   t0 = 1285 ns template. It is computed once per channel before the event loop
   (`templateArea[ch]`). Units: mV per bin × ns per bin = mV·ns = pVs. It
   agrees with the closed form `U0 · τ · σ · √(2π) · e^(σ²/2)`. Only bins
   500–660 are scaled and kept, and they hold over 99% of the area, so the
   injected pulse is under 1% short of the sampled area.
   - **Channel 8 takes `A` from the mean-tau template (row 8) for 8A and 8B
     too.** The 8A template has about 15% less area than row 8 and the 8B
     template about 8% more, so 8A and 8B pulses come out at the same height
     with different areas. That matches data, where U0 is 48.0 mV for the fast
     peak and 49.6 mV for the slow one.
   - `SPE_height` and the old constants are no longer used; the old lines are
     kept as comments.

10. **Max pulse height list has 96 entries.** `waveinject_slab.C` starts from
   commit `b679b4d` ("change mas pulse height to associated with 96 channels",
   2026-08-27). The `084ca95` version listed only 92 `maxValues`, so channels
   92–95 read past the end of the list.

11. **SPE area sampling range is `[0, 30000]` pVs** (was `[0, 5000]`).
    `TF1::GetRandom()` only samples inside the TF1 range. The 5000 came from
    the bars, whose SPE area is much smaller. The slab SPE areas are about
    5 nVs, and up to 14.4 nVs for channel 92, so 5000 cut off most of the
    distribution for 55 of 96 channels. For example, channel 92's drawn mean was
    3317 pVs instead of 14371, and channel 8's 3029 instead of 4389.
    - 30000 is at least mean + 3σ for every channel. The closest is channel 92
      at +2.98σ, which removes 0.14% of its draws and lowers its mean by 0.1%.
    - The lower edge at 0 means no negative area is ever drawn. Channels whose
      σ is about half their mean therefore come out slightly high
      (channel 0: 1726 → 1772).
    - Keep the upper edge above the largest mean + 3σ if `SPEmeans.txt` is
      redone.

## Checks done

- Sim PMTs 34 and 35 map to `remappedPMT` 8. Neighbouring PMTs map to 73
  (32/33) and 9 (36/37).
- Every template peaks at 50.0 mV at t0 + tau (8A at t0+38.8 ns, 8B at +63.8,
  92 at +96.2). Under 1% of each template's area lies outside bins 500–660.
- The full-record template area equals the closed form
  `U0 · τ · σ · √(2π) · e^(σ²/2)`: 5677 pVs for channel 0 and 5066 pVs for
  channel 8 (row 8).
- `beamMuonSlab_1kEvent.root` (1000 events) runs through. 96 + 2 templates are
  built, and 207 events have a channel-8 pulse above 5 mV.

## Validation test plan: SPE-only injection vs the v4 data fits

Status: steps 1–3 done (see "Progress" below).

**Goal.** Check that the injected SPE pulses look like the real SPE pulses the
templates were built from. Inject SPEs with no Geant4 input, fit them exactly as
the data was fitted, and compare the fitted quantities with `allFits_v4.root`.

**Everything runs locally.**
- No MilliDAQ conversion and no ntuplizer: the fit reads the injected
  `waveform[6][16][1024]` array directly.
- `runSpeFit.cc` is only copied down from lxplus as a reference. Nothing is
  changed or run on lxplus.

### Steps

1. **`waveinject_slab_SPEonly.C`** (new, `inputData/`).
   - 2000 events, with exactly one PE in every one of the 96 channels per event.
     Summing is ignored: one channel is one PE.
   - Channel 8 alternates between its two PMT shapes: even events use 8A
     (fast tau) and odd events 8B (slow tau), so each gets exactly 1000.
   - No QE line (`if (randGen.Uniform() > calibration) continue;` is dropped),
     since only SPEs matter here and nPE is not studied.
   - The pulse area is sampled from `SPEmeans.txt` and scaled with
     `scaleFactor = area / templateArea[ch]`, as in `waveinject_slab.C`.
   - Same output tree as `waveinject_slab.C` (`waveform`, `eventWeight`), plus
     a branch with the injected area per channel, used for the harvest-window
     selection in step 3.
2. **Local fit macro.** Reads the SPE-only output and fits every channel of
   every event with the fit from `runSpeFit.cc`, using the v4 settings:
   - the log-normal model
   - `--pre 200 --post 300`
   - per-bin error equal to the `pedestalRMS` of the pre-trigger samples
   - the same parameter bounds
   - the seed taken from the known injection time

   Writes a tree with the v4 branch names (`U0`, `tau`, `sigma`, `peakHeight`,
   `fitRiseNs_ln`, `fitFallNs_ln`, `status_ln`, `redChi2_ln`, ...).
3. **Comparison macro.** Per channel, sim vs v4 data:
   - **Selection, the same on both sides.** For the sim, the harvest window
     `1000·(mean ± sigma)` from `SPEmeans.txt` applied to the injected area;
     the data already passed it before fitting. Then
     `status_ln == 0 && 0.5 <= redChi2_ln <= 1.5`.
   - **Compared:** height (`U0`, `peakHeight`), rise and fall time
     (`fitRiseNs_ln`, `fitFallNs_ln`), `tau` and `sigma`.
   - **Not compared: pulse area.** The data area in this study is the offline
     pulse-finder `area`, which the v4 fit tree does not store. An area computed
     on the sim waveform would use a different definition.

### Expected limits

- **Noise.** The sim noise is 1 mV white (`rms_noise`), while data noise has
  per-channel RMS and spectral lines. The chi2 cut therefore keeps different
  fractions on the two sides, so compare means and shapes, not pass rates.
- **Shape spread.** Each channel has one fixed shape, so the sim `tau` and
  `sigma` are narrow peaks broadened only by noise. The test checks the mean
  shape, not the data's spread (tau RMS 5–14 ns).
- **Height vs area.** The sim height is exactly proportional to the sampled
  area, so the sim height–area correlation is tighter than in data.
- **Expected height gap.** The heights should come out about 16% low on
  average (see "Not done yet"). This test measures that gap per channel with the data's own
  fit.
- **Hit time.** Inject at a fixed early hit time (0–20 ns) so the whole pulse
  stays inside bins 500–660 (see the late-hit item below).

### Progress

- **Step 1 is written:** `waveinject_slab_SPEonly.C`, not committed yet.
  - Run it compiled:
    `root -l -b -q 'waveinject_slab_SPEonly.C+("SPEonly_2kEvent_waveinjected.root", 2000)'`.
    It needs no sim library, so ACLiC works. 2000 events take about 30 s.
  - Arguments: `(outputFile, nEvents=2000, hitTime=0 ns, lognormalFile, speMeansFile)`.
  - The SPE area Gaussian is sampled on `[0, 30000]` pVs, the same as
    `waveinject_slab.C` (decision 11).
  - Check on 2000 events:

    | Channel | injected area (pVs) | waveform integral (pVs) | mean max (mV) |
    |---|---|---|---|
    | 0 | 1792 | 1782 | 17.2 |
    | 8 | 4689 | 4479 | 47.3 |
    | 92 | 14488 | 14404 | 122.4 |

  - Channel 8 split exactly 1000/1000. The 8A events peak at 1328 ns and the
    8B events at 1349 ns.
  - Channel 8's waveform integral is 4% below its injected area by design: 8A
    and 8B are scaled with row 8's area.
- **Step 2 is done: `fitSPEonly.C` (+ `runSpeFit_core.h`), run locally.**
  - `runSpeFit_core.h` is lines 16–528 of `runSpeFit.cc` copied verbatim
    (minus the MilliDAQ include). The source is the local `milliqanOffline`
    checkout, branch `SlabPulseInjectionSim` at `2fb617aa`. The lxplus
    `~/condorTest/speFit.exe` that produced v4 was built 19 minutes after that
    commit, and its strings match it.
  - Before fitting, every pulse goes through a port of the offline pulse
    finder: `OfflineFactory::findPulses` and `processChannel`, with
    `pulseFindingSlab.json` (high 10 mV, low 1 mV, 3/3 consecutive samples).
    That gives the offline `time` (the fit seed), `area`, `height`,
    `riseSamples`, `fallSamples`, `duration` and `pickupFlagTight`.
  - `passHarvest` applies the data selection:
    - from `extract_ids.C`: `npulses == 1`, `riseSamples >= 3`,
      `fallSamples >= 3`, `pickupFlagTight == 0`
    - from `make_spe_harvest.C`: `1000·(mean ± 1σ)` on the offline area,
      edges exclusive

    `prePulseMean` and the calibrated-time window are not applied; the sim
    pulses always pass them.
  - The fit is `fitPulse`, log-normal, `--pre 200 --post 300`, with
    `isSlab = false`, since the injected pulses are already positive.
  - Run: `root -l -b -q 'fitSPEonly.C+("SPEonly_2kEvent_waveinjected.root","SPEonly_2kEvent_fits.root")'`.
    2000 events take about 2.5 min.
  - Result: 192000 channel-events in total. 11233 had no offline pulse (below
    the 10 mV threshold), 180767 were fitted, and 125662 pass the harvest
    selection, of which 125602 also pass the χ² cut. That is 1046–1996 per
    channel. Mean sim χ²_red is 1.008.
  - The fit recovers each template's τ and σ (channel 0: τ 79.43 vs 79.40,
    σ 0.504 vs 0.503).
- **Step 3 is done: `compareSPEfits.C`.**
  - Run: `root -l -b -q 'compareSPEfits.C+("SPEonly_2kEvent_fits.root","<path>/allFits_v4.root","SPEfitCompare")'`.
    It takes about 15 s.
  - Output in `SPEfitCompare/`:
    - one PDF per parameter (`U0`, `tau`, `sigma`, `peakHeight`,
      `fitRiseNs_ln`, `fitFallNs_ln`), one page per channel, each page
      showing the normalised "pulse injection" and "offline fitting data"
      histograms
    - `compareSPEfits.root`, with all histograms
    - `compareSPEfits_means.csv`, with the per-channel means
  - Mean of sim/data per-channel mean ratios:

    | Parameter | mean | rms | extremes |
    |---|---|---|---|
    | `tau` | 1.000 | 0.004 | 0.965 (ch8) |
    | `sigma` | 1.000 | 0.002 | 1.014 (ch8) |
    | `U0` | 0.913 | 0.058 | 0.698 (ch93) … 1.019 (ch17) |
    | `peakHeight` | 0.924 | 0.047 | 0.759 (ch93) … 1.038 (ch17) |
    | `fitRiseNs_ln` | 1.033 | 0.044 | up to 1.201 (ch72) |
    | `fitFallNs_ln` | 1.030 | 0.035 | up to 1.171 (ch72) |

  - **τ and σ match the data means by construction,** but the sim
    distributions are narrow peaks, as expected for one fixed shape per
    channel.
  - **The rise and fall times are on average 3% longer in the sim, and up to
    20% longer on some channels.** A template built from the mean τ and mean σ
    need not have the mean rise and fall time. This is not yet checked.
  - `fitFallNs_ln` is negative in both samples. That is how
    `measureFittedWidths` defines it (`t80f - t20f`), not an error.
  - **U0 is about 9% low on average, worst on small-area channels.** See the
    height-gap item in "Not done yet".

- **Area check: `compareSPEarea.C`.** This compares the SPE-only injection with
  its own input, so v4 is not needed. The waveform area is the sum of all 1024
  bins × 2.5 ns, for every injected pulse, with no selection.
  - Run: `root -l -b -q 'compareSPEarea.C+("SPEonly_2kEvent_waveinjected.root","SPEmeans.txt","SPEareaCompare")'`.
    It takes about 10 s.
  - Output in `SPEareaCompare/`:
    - `SPEarea_per_channel.pdf`: waveform area vs the input Gaussian
      (`SPEmeans.txt` mean ± σ, not truncated)
    - `SPEarea_diff_per_channel.pdf`: waveform − injected, per pulse
    - `compareSPEarea_summary.csv`
    - PNG previews for channels 0, 8 and 72

    The pulse-by-pulse ratio was dropped: for tiny injected areas the 80 pVs
    of noise makes it meaningless.
  - **The injected areas follow the assigned Gaussian.** Their mean matches the
    truncated-Gaussian mean to −0.02% on average, at most 2.4% (channel 74),
    which is the statistics of 2000 draws.
  - **The waveform area keeps 0.992 of the injected area on average**
    (RMS 0.005 over the 95 channels other than 8). The loss is exactly the part
    of the template beyond bin 660, which the injection zeroes. Channel 72 has
    the longest τ (122 ns) and keeps 0.970, matching the 0.969 of its template
    that lies inside bins 500–660.
  - **The spread of waveform − injected is 82 pVs on average,** the
    expected noise: sqrt(1024) × 1 mV × 2.5 ns = 80 pVs.
  - **Channel 8 is off by design:** 8A keeps 0.84 of its injected area
    (−738 pVs on average) and 8B 1.07 (+319 pVs), because both are scaled with
    row 8's area A. The two together keep 0.955.

- **SPE area sampling range: the group chose `[0, 30000]`.** That is what
  `waveinject_slab.C` and `waveinject_slab_SPEonly.C` use. The alternative
  `[mean − σ, mean + σ]` was tried with 2000 events and dropped. The test code
  for it was not kept.

  Value / input `SPEmeans.txt` value, over the 95 channels other than 8:

  | range | injected mean | injected RMS | waveform mean | waveform RMS |
  |---|---|---|---|---|
  | `[0, 30000]` (chosen) | 1.022 (0.994–1.094) | 0.951 (0.854–0.993) | 1.014 | 0.945 |
  | mean ± 1σ (dropped) | 1.000 (0.986–1.014) | 0.538 (0.523–0.555) | 0.992 | 0.537 |

  - **`[0, 30000]`:** no negative areas. The mean is biased up by 2.2% on
    average and up to 9.4% (channel 67), because the negative tail is cut, and
    the RMS is 5% low on average.
  - **mean ± 1σ:** the mean is exact and there are no negative areas, but the
    RMS is only 0.54σ (as expected for a Gaussian truncated at ±1σ, 0.5396),
    and the shape is flat-topped with hard edges. That would make summed
    multi-PE areas about half as wide as they should be.

### Recreate the validation plots

The output files are not committed. To regenerate them, run from
`inputData/` with the env from `setup_mac_g4p11.sh`. The whole chain takes
about 3 min.

```bash
# 1. SPE-only injection, 2000 events (~30 s). About 740 MB, keep it out of git.
root -l -b -q 'waveinject_slab_SPEonly.C+("SPEonly_2kEvent_waveinjected.root", 2000)'

# 2. Offline pulse finding + harvest selection + v4 log-normal fit (~2.5 min, ~15 MB)
root -l -b -q 'fitSPEonly.C+("SPEonly_2kEvent_waveinjected.root","SPEonly_2kEvent_fits.root")'

# 3. Per-channel comparison with the data fits (~15 s) -> SPEfitCompare/
root -l -b -q 'compareSPEfits.C+("SPEonly_2kEvent_fits.root","/Users/haoliangzheng/Documents/Working on lxplux/spe_results/allFits_v4.root","SPEfitCompare")'
```

- Step 1 reproduces the same events. The noise uses `TRandom3 randGen(2004)`,
  and the area draws (`TF1::GetRandom`) use ROOT's `gRandom` with its fixed
  default seed.
- `allFits_v4.root` is not in this repo. A local copy is in
  `~/Documents/Working on lxplux/spe_results/`, and the EOS copy is
  `/eos/user/h/hazheng/SlabSPEFits/allFits_v4.root`.
- For the area check, run
  `root -l -b -q 'compareSPEarea.C+("SPEonly_2kEvent_waveinjected.root","SPEmeans.txt","SPEareaCompare")'`
  on the step 1 output (about 10 s).
- Step 3 writes `SPEfitCompare/<param>_per_channel.pdf` for `U0`, `tau`,
  `sigma`, `peakHeight`, `fitRiseNs_ln` and `fitFallNs_ln`, one page per
  channel. It also writes `compareSPEfits.root` and `compareSPEfits_means.csv`.
- `SPEheight_vs_dataU0.csv` (the height-gap numbers in "Not done yet") came
  from a one-off check macro that is not in the repo. For each channel it
  used:
  - the injected height `50 × <drawn area> / A`, where `<drawn area>` is the
    mean of 100k `TF1::GetRandom()` draws on `[0, 5000]` and on `[0, 30000]`
  - the data U0 from `per_channel_v4/ln_params_per_channel.csv`

## Not done yet

- [ ] **Injected SPE heights are about 16% below data.** For a mean-area SPE,
  the injected height is `50 * <drawn area> / A`. The drawn area is the mean of
  the `SPEmeans.txt` Gaussian on `[0, 30000]`. Against the mean fitted U0 in data
  (`per_channel_v4/ln_params_per_channel.csv`), the ratio is 0.838 on average,
  with an RMS of 0.075 over 96 channels. Per-channel numbers are in
  `SPEheight_vs_dataU0.csv` (it also has the old `[0, 5000]` column, which
  averages 0.687).
  - Examples: channel 0 is 15.6 vs 23.3 mV (0.67), channel 8 is 45.5 vs 49.1
    (0.93), channel 92 is 120.8 vs 122.2 (0.99).
  - Lowest ratios: channels 93 (0.61), 72 (0.63), 80 (0.64), 0 (0.67).
  - **The gap is larger for small-area channels:** the ratio is 0.77 for
    SPEmean < 2500 pVs (17 channels), 0.83 for 2500–4000 (46) and 0.89 for
    ≥ 4000 (33). The correlation is 0.56.

  **Cause, confirmed:** the offline `area` is `Integral("width")` between the
  pulse bounds only. The bounds start at the last sample below the **10 mV**
  high threshold and end after 3 samples below it, so the whole tail under
  10 mV is left out. On the sim pulses, offline area / injected area is 0.66
  for channel 0, 0.75 for channel 9, 0.82 for channel 8 and 0.96 for
  channel 92. `SPEmeans.txt` therefore holds offline (partial) areas, and
  scaling a full template to them gives pulses that are too small. Small
  pulses lose the most. Fix: rescale `SPEmeans.txt` to full-pulse areas, or
  scale by `area / (offline-style area of the template)`.
- [ ] **Check the >5000-hit saturation claim** (height above 1250 mV, pulse
  longer than 2560 ns) on real high-hit events.
- [ ] **Noise and zeroing on summed channels.**
  - In the >5000-hit branch, the line
    `for (...) waveform[digitizer][channel][bin] = 0.0;` wipes the partner
    PMT's pulse if that PMT was processed first.
  - The per-hit branch adds noise once per PMT, so a summed channel gets it
    twice. That is accepted for now, since two PMTs feed the summing amp.
  - The per-hit branch also zeroes bins outside 500–660 after the partner PMT
    was added.
- [ ] **Only bins 500–660 are scaled.** In the >5000-hit branch the template
  tail beyond bin 660 (t > 1650 ns) stays at its unscaled 50 mV-normalised
  value. That is under 1% of the area, but it is not scaled with the pulse.
- [ ] **The calibration array `cali`** (QE / 10) is still the bar-detector-style
  list noted as a TODO in the file header.
- [ ] **The channel 92 high-gain assumption is not yet confirmed.** Check the
  raw stage-1 area spectrum for a peak near 3 nVs. If there is one, the
  calibration is wrong and channel 92 needs re-harvesting.
- [ ] **Late hits lose their tail in the per-hit path.** After the time shift,
  bins outside 500–660 are zeroed. A hit later than about 100 ns therefore
  loses its tail, and one later than about 360 ns is erased entirely. This
  affects the Geant4 injection as well as the SPE-only test.
- [x] **Validation test** (see "Validation test plan" above). Steps 1–3 are
  done. The open follow-ups are the height gap and the rise/fall excess.
- [ ] **Commit the remaining files.** The `waveinject_slab.C` changes are
  committed and pushed (`5e9916f`, `8b35e96`), but three files are still
  untracked: `extractLogNormalMeans.C`, `SPElognormal.txt` and this file.
  `waveinject_slab.C` reads `SPElognormal.txt` by default, so a fresh clone of
  the branch cannot run it until `SPElognormal.txt` is committed.

## Unnecessary checks to remove (do at the end)

Guards against input that never happens with our files. Each item is named by
its code, because line numbers shift as the file is edited.

### `waveinject_slab.C`: log-normal loading and channel 8 (added on this branch)

- [ ] `if (chan < 0 || chan >= nSPEChan) { delete h; continue; }`: the file
  only has 0–95, 8A and 8B.
- [ ] `line.empty() ||` in `if (line.empty() || line.find("chan") != std::string::npos) continue;`:
  there are no blank lines. Keep the `line.find("chan")` part, which skips
  the header.
- [ ] `if (!(iss >> label >> tau >> tauErr >> sigma)) continue;`: every row has
  the same format, so plain `iss >> label >> tau >> tauErr >> sigma;` works.
- [ ] The `nMissing` / 8A / 8B "Log-normal SPE file is incomplete" abort block:
  `extractLogNormalMeans.C` always writes all 98 rows.
- [ ] `if (!in.is_open()) { ... return; }` for `lognormalFile`: optional. If the
  path is wrong, nothing is loaded and the `templateArea` loop crashes on a
  null pointer. Keep this one (or the `nMissing` block) for a clear message,
  or drop both.
- [ ] `if (PMT_number == 34) ... else if (PMT_number == 35) ...`: channel 8 only
  comes from PMTs 34 and 35, so this can be
  `pulse_shape = (PMT_number == 34) ? speTemplate8A : speTemplate8B;`.

Keep these, they are not checks:
- `h->SetDirectory(nullptr)`: without it all 98 templates get written to the
  output file.
- `kTemplatePed = 0.0`: it records the ped = 0 decision.

### `waveinject_slab.C`: `SPEmeans.txt` loading and the hit loop (from `084ca95`)

- [ ] `if (line.empty()) continue;`: no blank lines.
- [ ] `if (!(iss >> chan >> mean >> meanErr >> sigma >> sigmaErr)) continue;`:
  fixed format.
- [ ] `if (chan < 0 || chan >= nSPEChan) continue;`: only 0–95.
- [ ] `if (!(mean > 0.0 && sigma > 0.0)) continue;`: the same kind as the
  removed tau/sigma check.
- [ ] The `nRead == 0` "No valid SPE means loaded" abort: only triggers on an
  empty file.
- [ ] `if (remappedPMT < 0 || remappedPMT >= nSPEChan || !fit[remappedPMT]) continue;`
  in the hit loop. The mapping always gives 0–95 (`pmtPosition` asserts), and
  all 96 `fit[]` exist. It also comes after `cali[remappedPMT]` has already
  read with that index, so it guards nothing.

### `extractLogNormalMeans.C`

- [ ] `if (chan < 0 || chan >= nChan) continue;`
- [ ] The `a.n == 0` "no selected fits ... row skipped" warning in `writeRow`.
- [ ] The `haveLN` fallback to generic `status`/`redChi2`: `allFits_v4.root`
  always has `status_ln` and `redChi2_ln`, so use those branches directly.
- [ ] The `(n > 1 && var > 0)` guard in the error-on-the-mean calculation.
