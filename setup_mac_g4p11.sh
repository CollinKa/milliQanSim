#!/usr/bin/env bash
# Activate the conda env and Geant4 11.2.2 used to build this checkout on macOS.
# Source it from the repository root before building or running:
#   source ./setup_mac_g4p11.sh
# Override the defaults if your env or Geant4 install live elsewhere:
#   MQ_CONDA_ENV=<env> MQ_GEANT4_INSTALL=<prefix> source ./setup_mac_g4p11.sh

if [[ "${BASH_SOURCE[0]}" == "${0}" ]]; then
    echo "Please source this script: source ./setup_mac_g4p11.sh"
    exit 1
fi

MQ_CONDA_ENV=${MQ_CONDA_ENV:-milliqan-sim-g4p11}
MQ_GEANT4_INSTALL=${MQ_GEANT4_INSTALL:-$HOME/software/geant4-11.2.2-qt-install}

if [[ -n "${CONDA_EXE:-}" ]]; then
    conda_base=$(dirname "$(dirname "$CONDA_EXE")")
else
    conda_base=/opt/miniconda3
fi
source "$conda_base/etc/profile.d/conda.sh"
conda activate "$MQ_CONDA_ENV" || return 1

if [[ ! -f "$MQ_GEANT4_INSTALL/bin/geant4.sh" ]]; then
    echo "No Geant4 install at $MQ_GEANT4_INSTALL; set MQ_GEANT4_INSTALL or run build_geant4_mac.sh"
    return 1
fi
source "$MQ_GEANT4_INSTALL/bin/geant4.sh"

export Geant4_DIR=$MQ_GEANT4_INSTALL/lib/cmake/Geant4

echo "Activated $MQ_CONDA_ENV with Geant4 $(geant4-config --version) from $MQ_GEANT4_INSTALL"
