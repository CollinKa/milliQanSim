#!/usr/bin/env bash
# Build Geant4 11.2.2 (Qt + Qt3D, multithreaded, datasets) against a conda env on macOS arm64.
# The env must already exist, e.g.
#   conda create -n milliqan-sim-g4p11 -c conda-forge python=3.11 "cmake<4" ninja git root boost qt-main qt
# Usage:
#   bash build_geant4_mac.sh [conda-env] [install-prefix] [source-dir]
# Defaults match setup_mac_g4p11.sh. The source and build trees go in [source-dir]
# (default ~/software/src), not /tmp, which macOS clears.
set -euo pipefail

ENV_NAME=${1:-milliqan-sim-g4p11}
PREFIX=${2:-$HOME/software/geant4-11.2.2-qt-install}
SRC_ROOT=${3:-$HOME/software/src}
G4_VER=11.2.2

# conda's compiler activate.d scripts reference unset variables, so relax -u here
set +u
if [[ -n "${CONDA_EXE:-}" ]]; then
    conda_base=$(dirname "$(dirname "$CONDA_EXE")")
else
    conda_base=/opt/miniconda3
fi
source "$conda_base/etc/profile.d/conda.sh"
conda activate "$ENV_NAME"
set -u

cmake_major=$(cmake --version | awk 'NR == 1 { split($3, v, "."); print v[1] }')
if (( cmake_major >= 4 )); then
    echo "CMake $cmake_major.x found; Geant4 $G4_VER needs cmake<4 (conda install -n $ENV_NAME -c conda-forge 'cmake<4')"
    exit 1
fi

mkdir -p "$SRC_ROOT"
cd "$SRC_ROOT"
if [[ ! -d geant4-v$G4_VER ]]; then
    curl -sSL -o geant4-v$G4_VER.tar.gz \
        https://gitlab.cern.ch/geant4/geant4/-/archive/v$G4_VER/geant4-v$G4_VER.tar.gz
    tar xzf geant4-v$G4_VER.tar.gz
fi

# Always configure from a clean build tree: cached flags survive a reconfigure
build_dir=geant4-v$G4_VER-build
rm -rf "$build_dir"
mkdir "$build_dir"
cd "$build_dir"

# CMAKE_PREFIX_PATH lets CMake find Qt5 in the env; OpenGL/X11 is off because
# the XQuartz GLX viewer does not work on macOS, so the viewer is OGLSQt.
args=(
    -G Ninja
    -DCMAKE_INSTALL_PREFIX="$PREFIX"
    -DCMAKE_PREFIX_PATH="$CONDA_PREFIX"
    -DCMAKE_OSX_ARCHITECTURES=arm64
    -DGEANT4_INSTALL_DATA=ON
    -DGEANT4_BUILD_MULTITHREADED=ON
    -DGEANT4_USE_QT=ON
    -DGEANT4_USE_QT3D=ON
    -DGEANT4_USE_OPENGL_X11=OFF
    -DGEANT4_USE_SYSTEM_ZLIB=ON
)

start=$(date +%s)
echo "Configuring (log: $PWD/configure.log)"
cmake "${args[@]}" ../geant4-v$G4_VER > configure.log 2>&1
echo "Building (log: $PWD/build.log)"
ninja -j"$(sysctl -n hw.ncpu)" > build.log 2>&1
echo "Installing to $PREFIX"
ninja install > install.log 2>&1

# The env's paths are baked into the install, so record which env it belongs to
{
    echo "Geant4 $G4_VER built $(date '+%F %T') in conda env $ENV_NAME ($CONDA_PREFIX)"
    echo "cmake ${args[*]} <src>"
} > "$PREFIX/BUILD_INFO.txt"

echo "Done in $(( ($(date +%s) - start) / 60 )) min"
echo "Geant4 $("$PREFIX/bin/geant4-config" --version), Qt: $("$PREFIX/bin/geant4-config" --has-feature qt)"
