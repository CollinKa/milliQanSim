# milliQanSim slab beam muons on Apple Silicon with Geant4 11.2

This branch uses a separate native ARM64 Conda environment and a separate
Geant4 installation, so it does not replace the Geant4 10.7 setup used by the
older `MacArmVer` checkout.

Verified on this Mac:

- Conda environment: `milliqan-sim-g4p11`
- Geant4: 11.2.2, ARM64, multithreaded
- Geant4 Qt and Qt3D support: enabled
- ROOT: 6.38.04
- CMake: 3.31.8
- Geant4 install: `/Users/haoliangzheng/software/geant4-11.2.2-qt-install`

## Activate this branch's environment

From the repository root:

```bash
cd /Users/haoliangzheng/Desktop/CERN/milliQanSim_Geant4P11/milliQanSim
source ./setup_mac_g4p11.sh
```

The script activates `milliqan-sim-g4p11`, sources Geant4 11.2.2, and sets
`Geant4_DIR` to the matching CMake package. Check the selected version with:

```bash
geant4-config --version
geant4-config --has-feature qt
```

The expected answers are `11.2.2` and `yes`.

## Build and run one slab beam-muon event

```bash
bash beamSetupSlab.sh
```

This selects `particlesMu.ini`, installs the slab source variants into the
active source paths, builds `build/MilliQan`, and runs one event using
`runMac/mcp_novis.mac`.

The smoke-test output is:

```text
build/beamMuonSlab_MilliQan.root
```

To rerun without reconfiguring:

```bash
cd /Users/haoliangzheng/Desktop/CERN/milliQanSim_Geant4P11/milliQanSim/build
./MilliQan ../runMac/mcp_novis.mac
```

## Open the Qt slab visualization

After `beamSetupSlab.sh` has selected and built the slab configuration:

```bash
cd /Users/haoliangzheng/Desktop/CERN/milliQanSim_Geant4P11/milliQanSim/build
./MilliQan ../runMac/visQT.mac -GUI
```

The program opens an `OGLSQt` window and draws the detector. In the Qt session
input, generate and display one beam-muon event with:

```text
/run/beamOn 1
/vis/viewer/flush
```

The initial geometry construction prints many overlap checks, so allow it to
finish before concluding that the window is unresponsive.

To check the event count:

```bash
python -c "import ROOT; f=ROOT.TFile.Open('beamMuonSlab_MilliQan.root'); print(f.Get('Events').GetEntries())"
```

## Switch back to the old Geant4 10.7 checkout

Changing Git branches does not change the active libraries. Open a fresh
terminal, or explicitly activate the old environment and source the old
Geant4 installation before working in the old checkout:

```bash
source /opt/miniconda3/etc/profile.d/conda.sh
conda activate milliqan-sim
source /Users/haoliangzheng/software/geant4-10.7.4-qt-install/bin/geant4.sh
cd /Users/haoliangzheng/Desktop/CERN/milliQanSim
```

Use `milliqan-sim-g4p11` only with this Geant4 11 checkout, and use
`milliqan-sim` only with the old Geant4 10.7 checkout. Keeping separate build
directories prevents CMake from retaining the wrong Geant4 paths.

## Known non-fatal warnings

The one-event run reports the existing 0.25 mm `airGapPanel_physic` geometry
overlap and several Geant4 shutdown cache warnings. The run still exits with
status 0 and writes a readable one-entry ROOT tree. Review the geometry overlap
before treating this configuration as production-validated physics.
