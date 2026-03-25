# HYPPO-dev (Forked)

This is a forked version of **HYPPO**, an implicit coupled supersonic/hypersonic flow solver built on [foam-extend 4.1](https://sourceforge.net/projects/foam-extend/) and [Mutation++](https://github.com/mutationpp/Mutationpp).

## Credits

- **Original repository**: [https://gitlab.com/PATO/hyppo-dev](https://gitlab.com/PATO/hyppo-dev)
- **Original authors**: Jeremy Chevalier, Dylan Gomesse
- **Professor / Project lead**: Jean Lachaud — [https://jeanlachaud.com/](https://jeanlachaud.com/)

Full credit goes to the original authors and the PATO team for the development of HYPPO.

## What is HYPPO?

HYPPO solves the compressible Navier-Stokes equations for supersonic and hypersonic flows using a block-coupled implicit approach. It computes the density, momentum, and total energy using an implicit Rusanov scheme with GMRES solver and Cholesky preconditioner.

### Gas Models

Selected via `constant/transportProperties`:

| Model | Keyword | Description |
|-------|---------|-------------|
| Constant properties | `none` | Perfect ideal gas with constant viscosity, Cp, Pr |
| Sutherland | `sutherland` | Ideal gas with temperature-dependent viscosity |
| Equilibrium chemistry | `readTable` | Real gas via Mutation++ lookup tables (air5 or air11) |
| Finite-rate chemistry | `finiteRate` | 5-species air (N2, N, O2, O, NO) via Mutation++ |

### Numerical Methods

- **Inviscid fluxes**: Rusanov scheme
- **Viscous fluxes**: Centered scheme
- **Gradient computation**: Diamond scheme (handles non-orthogonal meshes)
- **Time integration**: Implicit Euler (1st/2nd order)
- **Linear solver**: GMRES with Cholesky preconditioner
- **CFL control**: Adaptive — constant for 100 iterations, then increases ~25% every 10 iterations up to user-specified max

### Turbulence

Optional Baldwin-Lomax algebraic model, activated after a user-specified time.

## Dependencies

- **foam-extend 4.1** — block-coupled matrix system (`fvBlockMatrix`)
- **OpenFOAM 7** — needed for mesh generation (`blockMesh`) since `blockMeshDict` files use `#calc` entries unsupported by foam-extend
- **Mutation++** — chemical kinetics and thermodynamic properties (via PATO)
- **Eigen** — linear algebra (bundled with Mutation++)

## Building

```bash
cd src
./Allwmake
```

This builds the boundary conditions library (`libhyperCoupledFoam`) and the solver binary (`hyppo`).

## Running a Tutorial

Meshing and solving require different OpenFOAM environments:

```bash
# Terminal 1: generate mesh with OpenFOAM 7
source /path/to/openfoam7/etc/bashrc
cd tutorials/halfSphere/<CaseName>
blockMesh

# Terminal 2: run solver with foam-extend 4.1
source /path/to/foam-extend-4.1/etc/bashrc
export MPP_DATA_DIRECTORY=/path/to/mutation++/data
export LD_LIBRARY_PATH=/path/to/mutation++/install/lib:$LD_LIBRARY_PATH
cd tutorials/halfSphere/<CaseName>
hyppo
```

### Available Tutorials

| Case | Mach | Altitude | Model | Approx. Runtime |
|------|------|----------|-------|-----------------|
| PerfectGasMach10 | 10 | 20 km | Ideal gas (`none`) | ~15s |
| ChemicalEquilibriumMach10 | 10 | 20 km | Equilibrium air5 (`readTable`) | ~1 min |
| ChemicalEquilibriumMach50 | 56 | 50 km | Equilibrium air11 (`readTable`) | ~6 min |
| FiniteRateMach5 | 5 | — | Sutherland (`sutherland`) | ~1 min |

## Bugs Fixed in This Fork

The original code had several bugs that prevented the tutorials from running. The following fixes were applied:

### 1. Unconditional Mutation++ initialization (`createChemistryFields.H`)

The `Mutation::Mixture` object was created for **all** gas models, including `"none"` (perfect gas) and `"sutherland"`, which do not require chemistry. This caused an immediate crash when Mutation++ tried to load data files that weren't needed.

**Fix**: Wrapped the Mutation++ initialization in a conditional so it only runs for `"finiteRate"` and `"readTable"` models. Changed from a stack-allocated object to a pointer (`Mutation::Mixture* mixPtr`) with local references inside conditional blocks.

### 2. Missing `return false` in `testEnergy()` (`simControl.C`)

This was the critical bug. The `testEnergy()` function checked whether internal energy was negative, but **did not return `false`** when energy was valid. Due to undefined behavior (falling off the end of a non-void function), the solver always interpreted this as "energy is negative", causing it to continuously reduce the timestep until it crashed with `"Lower Courant number needed"`.

This bug affected **every gas model** — even with a perfectly valid simulation, the solver would spiral the timestep down to zero.

**Fix**: Added `return false;` at the end of the function.

### 3. Unconditional mass fraction check (`createChemistryFields.H`)

The check `Y_O + Y_NO + Y_N + Y_N2 + Y_O2 != 1` ran for all models, but is only meaningful for finite-rate chemistry.

**Fix**: Wrapped in `if (model == "finiteRate")`.
