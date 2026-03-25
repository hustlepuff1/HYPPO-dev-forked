# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

HYPPO is an implicit coupled transonic/supersonic/hypersonic flow solver built on **foam-extend 4.0**. It solves the compressible Navier-Stokes equations (density, momentum, total energy) using a block-coupled implicit approach with Rusanov inviscid fluxes, Diamond gradient scheme, and centered viscous fluxes.

## Build Commands

```bash
# Build everything (BCs library + solver)
cd src && ./Allwmake

# Build only boundary conditions library
cd src && wmake libso BCs

# Build only solver
cd src && wmake solver
```

Output: solver binary to `$FOAM_USER_APPBIN/hyppo`, BC library to `$FOAM_USER_LIBBIN/libhyperCoupledFoam`.

## Dependencies

- **foam-extend 4.1** — uses fvBlockMatrix (block-coupled system), not compatible with standard OpenFOAM
- **OpenFOAM 7** — needed for `blockMesh` (foam-extend's blockMesh can't handle `#calc` entries in blockMeshDict)
- **Mutation++** — chemical kinetics / thermodynamics library, expected at `$HOME/PATO/src/thirdParty/mutation++/`
- **Eigen** — bundled within Mutation++ at `mutation++/thirdparty/eigen`
- **PATO** — provides Mutation++ data files at `$HOME/PATO/data/ThermoTransportChemistry/mutation++/`

## Environment Setup

The user has bash aliases configured in `~/.bashrc`:
- `loadHyppo` — loads foam-extend 4.1 + Mutation++ environment
- `of7` — loads OpenFOAM 7 (needed for `blockMesh`)

## Running Tutorials

Meshing and solving require **different OpenFOAM environments**:

```bash
# Step 1: Generate mesh with OpenFOAM 7 (supports #calc in blockMeshDict)
of7
cd tutorials/halfSphere/<CaseName>
blockMesh

# Step 2: In a NEW terminal, run solver with foam-extend
loadHyppo
cd tutorials/halfSphere/<CaseName>
hyppo

# Parallel run
decomposePar -force
mpirun -np <N> hyppo -parallel
reconstructPar

# Clean case
./Allclean
```

Note: The `Allrun` scripts in some tutorials are outdated and may not work.

## Architecture

### Solver (`src/solver/`)

The main time loop is in `hyppo.C`. It orchestrates field creation, equation assembly, and the block-coupled solve. Key include files:

- **Equations** (`rhoUE/`): `rhoEqn.H`, `rhoUEqn.H`, `rhoEEqn.H` — conservation of mass, momentum, and energy with implicit Jacobian terms assembled into a `fvBlockMatrix<vector6>`.
- **Numerics** (`numerics/`): `Diamond.C/.H` (gradient computation for non-orthogonal meshes), `rusanovFlux.H` (Rusanov inviscid flux with dissipation), `CenteredFluxViscosity.H` (viscous flux).
- **Gas models** — selected via `constant/transportProperties` `model` keyword:
  - `"none"` — perfect gas, constant properties
  - `"sutherland"` — ideal gas with Sutherland viscosity law
  - `"readTable"` — real gas via Mutation++ lookup tables (air5 or air11)
  - `"finiteRate"` — finite-rate chemistry, 5 species (N2, N, O2, O, NO)
- **Chemistry** (`finiteRate/`): species transport equations with production rates and Jacobians from Mutation++.
- **Turbulence** (`turbulence/`): Baldwin-Lomax algebraic model, optionally activated after a user-specified time.
- **Simulation control** (`simControl/`): adaptive CFL ramping — starts at ~0.3, increases 20% every 10 iterations up to user max; auto-restarts with reduced timestep on negative internal energy.
- **Table lookup** (`readTable/`): `ReadTable` class for interpolating equilibrium thermodynamic properties.

### Boundary Conditions (`src/BCs/`)

Compiled into `libhyperCoupledFoam`:

- `maxwellSlipU` — Maxwell slip velocity for rarefied flows (thermal creep, curvature effects)
- `smoluchowskiJumpT` — Smoluchowski temperature jump for rarefied flows
- `fixedRho` — fixed density BC
- `mixedFixedValueSlip` — generic mixed slip base class

### Tutorial Case Structure

```
CaseName/
├── 0/                    # Initial fields (U, p, T)
├── constant/
│   ├── transportProperties   # Gas model selection
│   ├── constantProperties    # Cp, gamma, mu, Pr
│   ├── turbulenceProperties  # Turbulence model config
│   └── fluidComposition      # Species (for chemistry models)
└── system/
    ├── controlDict           # Time control, output
    ├── fvSchemes             # Numerical schemes
    ├── fvSolution            # Solver settings (GMRES + Cholesky)
    └── blockMeshDict         # Mesh generation
```

## Key Implementation Details

- The solver uses foam-extend's `fvBlockMatrix<vector6>` for block-coupled implicit solving of (ρ, ρU, ρE) — this is fundamentally different from segregated OpenFOAM solvers.
- GMRES with Cholesky preconditioner is the default linear solver.
- Temporal discretization: implicit Euler (1st or 2nd order via `backward` scheme).
- Jacobians for the implicit coupling are computed via finite differences in the equation files.
- Platform-specific compiler flags exist in `Make/options` for Linux (g++) and Darwin (c++).

## Bugs Fixed (2026-03-25)

### 1. Unconditional Mutation++ initialization (`createChemistryFields.H`)
`Mutation::Mixture` was created for all gas models, even `"none"` (perfect gas). Crashed immediately when Mutation++ data files weren't needed. Fixed by wrapping initialization in `if (model == "finiteRate" || model == "readTable")` using a `Mutation::Mixture*` pointer, with local references (`Mutation::Mixture& mix = *mixPtr`) inside conditional blocks.

### 2. Missing `return false` in `testEnergy()` (`simControl.C:172-183`)
The function checked if internal energy was negative but had no `return false` for the valid-energy path. Undefined behavior caused the solver to always think energy was wrong, shrinking the timestep to zero on every run. This was the critical bug — it broke all gas models. Fixed by adding `return false;` at the end of the function.

### 3. Mass fraction check ran for all models (`createChemistryFields.H`)
The `Y_O + Y_NO + Y_N + Y_N2 + Y_O2 != 1` check ran unconditionally. Wrapped in `if (model == "finiteRate")`.
