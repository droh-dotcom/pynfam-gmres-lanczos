# PyNFAM GMRES and Lanczos Solvers

GMRES and Lanczos implementations for charge-changing quasiparticle random-phase approximation (QRPA) strength-function calculations, built on [PyNFAM](https://github.com/mld1812/pynfam).

This repository accompanies the manuscript *Lanczos Method for Quasiparticle Random-Phase Approximation Strength Functions in Atomic Nuclei* by Dong Min Roh, Chao Yang, Matthew L. Dai, and Jonathan Engel.

## Contents

- `Makefile`, `pnfam_setup.f90`, and `pnfam_solver.f90`: replacement files for PyNFAM's `exes/pnfam/` directory.
- `run_pynfam_GMRES.py`: example input script for GMRES calculations.
- `run_pynfam_LANCZOS.py`: example input script for Lanczos calculations.

These files extend an existing PyNFAM installation; the remaining PyNFAM source files are required.

## Requirements

Follow the [PyNFAM installation instructions](https://github.com/mld1812/pynfam#readme) to install its Python dependencies and build HFBTHO. The supplied Makefile defaults to GNU Fortran (`gfortran`) and requires GNU Make, BLAS, LAPACK, and the HFBTHO library and module files.

Before building pnFAM, adjust `L_LINEARALGEBRA`, `L_HFBTHO`, and `I_HFBTHO` in the supplied Makefile for your system. The default HFBTHO paths refer to `exes/hfbtho_blocking_ba9f1a8/src/hfbtho/`. Place the HFBTHO executable required by PyNFAM in `exes/` as well.

## Installation

Download this repository into a directory named `pynfam-gmres-lanczos` and clone PyNFAM alongside it:

```bash
git clone https://github.com/mld1812/pynfam.git
```

The directory layout should be:

```text
parent-directory/
├── pynfam-gmres-lanczos/
│   ├── Makefile
│   ├── pnfam_setup.f90
│   ├── pnfam_solver.f90
│   ├── run_pynfam_GMRES.py
│   └── run_pynfam_LANCZOS.py
└── pynfam/
```

From the parent directory, install the Python package and copy the replacement files and example scripts:

```bash
cd pynfam
python3 -m pip install .
cp ../pynfam-gmres-lanczos/Makefile exes/pnfam/
cp ../pynfam-gmres-lanczos/pnfam_setup.f90 exes/pnfam/
cp ../pynfam-gmres-lanczos/pnfam_solver.f90 exes/pnfam/
cp ../pynfam-gmres-lanczos/run_pynfam_GMRES.py .
cp ../pynfam-gmres-lanczos/run_pynfam_LANCZOS.py .
```

The copy commands replace the corresponding upstream files. Complete the HFBTHO setup described above before proceeding.

## Build and run

The solver is selected **at compile time** with `SOLVER=GMRES` or `SOLVER=LANCZOS`. Rebuild and copy the executable whenever you switch solvers; the Python script alone does not select the solver.

### GMRES

Starting from the PyNFAM root directory:

```bash
cd exes/pnfam
make clean
make pnfam_main.x SOLVER=GMRES
cp pnfam_main.x ../
cd ../..
python3 run_pynfam_GMRES.py
```

### Lanczos

Starting from the PyNFAM root directory:

```bash
cd exes/pnfam
make clean
make pnfam_main.x SOLVER=LANCZOS
cp pnfam_main.x ../
cd ../..
python3 run_pynfam_LANCZOS.py
```

The explicit `pnfam_main.x` target builds the executable only. Plain `make SOLVER=GMRES` or `make SOLVER=LANCZOS` also requests the library and documentation targets.

## Inputs and results

Run the scripts from the PyNFAM root so that `./exes` and `./tests` resolve correctly. Edit `pynfam_inputs` and `override_settings` in each script to change the nucleus, model space, interaction, energy grid, or iteration settings. In the Lanczos example, `max_iter` sets the number of Lanczos iterations, and the contour requests a single evaluation at zero frequency.

The supplied scripts write results under:

```text
tests/TEST_GMRES/
tests/TEST_LANCZOS/
```

Lanczos additionally writes `Theta.txt`, `RF0.txt`, and `RF1.txt` in the individual pnFAM working directories. These contain projected eigenvalues and transition amplitudes for reconstructing the strength function as described in the manuscript. Use these quantities for the Lanczos reconstruction; the conventional strength entries printed by pnFAM in this mode are not the reconstructed profile.

The supplied scripts use a small two-shell example with three protons and three neutrons. They must be configured with the appropriate inputs to reproduce the manuscript's production calculations.

## Acknowledgment

Please cite the accompanying manuscript and acknowledge the original [PyNFAM project](https://github.com/mld1812/pynfam) when using this work. Retain the original attribution notices in the source files.
