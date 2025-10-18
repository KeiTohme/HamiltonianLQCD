# Hamiltonian LQCD Documentation

## Overview

This is a comprehensive implementation of Hamiltonian Lattice Quantum Chromodynamics (LQCD) with a focus on real-time gluon dynamics. The program is written primarily in Fortran 90/95 and supports various lattice geometries, fermion formulations, and parallelization options.

## Features

### 1. Lattice Types
- **Square Lattice**: Hypercubic lattice with periodic boundary conditions
- **Hexagonal Lattice**: 2D hexagonal lattice structure

### 2. Fermion Formulations
- **Wilson Fermions**: Standard Wilson discretization with improved locality
- **Staggered Fermions**: Kogut-Susskind formulation with reduced fermion doubling

### 3. Gauge Group
- SU(N) gauge theory with arbitrary N
- Pre-implemented generators for SU(2) and SU(3)
- Generic implementation for higher SU(N) groups

### 4. Hamiltonian Formulation
- Electric field energy (kinetic term)
- Magnetic field energy (plaquette term)
- Real-time evolution using symplectic integrators

### 5. Parallelization
- **OpenMP**: Shared-memory parallelization (optional)
- **MPI**: Distributed-memory parallelization (optional)
- Both can be enabled simultaneously for hybrid parallelization

### 6. Computational Methods
- **Standard**: Classical simulation on CPU
- **Tensor Network**: Framework for tensor network methods
- **Quantum Computer**: Framework for quantum computing implementation

## Building the Code

### Prerequisites
- Fortran 90/95 compiler (gfortran, ifort, etc.)
- Optional: MPI library for distributed parallelization

### Basic Build
```bash
make
```

### Build with OpenMP
```bash
make openmp
# or
make USE_OPENMP=yes
```

### Build with MPI
```bash
make mpi
# or
make USE_MPI=yes
```

### Build with both OpenMP and MPI
```bash
make hybrid
```

### Clean Build
```bash
make clean
make rebuild
```

## Running the Code

### Basic Usage
```bash
./hamiltonian_lqcd input_file.dat [output_file.dat]
```

### Examples
```bash
# Small test case
./hamiltonian_lqcd examples/input_small_test.dat test_output.dat

# Square lattice simulation
./hamiltonian_lqcd examples/input_square_lattice.dat square_results.dat

# Hexagonal lattice simulation
./hamiltonian_lqcd examples/input_hexagonal_lattice.dat hex_results.dat
```

### Running with MPI
```bash
mpirun -np 4 ./hamiltonian_lqcd input.dat output.dat
```

## Input File Format

The input file uses a simple key-value format with comments starting with `#`.

### Parameters

#### Physical Parameters
- `gauge_coupling` (or `g`): Gauge coupling constant
- `fermion_mass` (or `m`): Fermion mass
- `num_colors` (or `Nc`): Number of colors (typically 3 for QCD)

#### Lattice Parameters
- `spacetime_dim` (or `D_euclid`): Space-time dimension
- `lattice_size` (or `Nsize`): Comma-separated lattice sizes for each dimension
- `lattice_type`: Either `square` or `hexagonal`

#### Fermion Parameters
- `fermion_type`: Either `wilson` or `staggered`

#### Time Evolution Parameters
- `time_step` (or `dt`): Time step for evolution
- `num_time_steps`: Number of time steps to simulate

#### Parallelization Options
- `use_openmp`: Enable OpenMP (true/false)
- `use_mpi`: Enable MPI (true/false)
- `num_threads`: Number of OpenMP threads

#### Computational Method
- `comp_method`: Method to use (`standard`, `tensor_network`, or `quantum`)

## Output

The program generates a data file with the following columns:
1. Step number
2. Time
3. Total energy (Hamiltonian)
4. Average plaquette value

This data can be used to analyze the time evolution of the gauge field configuration.

## Theory Background

### Hamiltonian Formulation

The Hamiltonian for lattice gauge theory is:

```
H = (g²/2) Σ_x Σ_i E_i^a(x)² - (1/g²) Σ_plaq Re[Tr(U_plaq)]
```

where:
- E_i^a(x) is the electric field (conjugate momentum to gauge links)
- U_plaq is the gauge link product around plaquettes
- g is the gauge coupling constant

### Time Evolution

Real-time evolution is performed using Hamilton's equations:
- dU/dt = g² * E * U
- dE/dt = -∂S/∂U

We use a symplectic (leapfrog) integrator to preserve the Hamiltonian structure.

### Wilson Fermions

Wilson fermions add a term to remove fermion doubling:
```
D_Wilson = m + Σ_μ [(1-γ_μ) U_μ(x) ψ(x+μ) + (1+γ_μ) U_μ†(x-μ) ψ(x-μ)]/2 - (r/2) Σ_μ ∇_μ² ψ
```

### Staggered Fermions

Staggered fermions use phase factors to reduce fermion doubling:
```
D_staggered = m χ(x) + Σ_μ η_μ(x) [U_μ(x) χ(x+μ) - U_μ†(x-μ) χ(x-μ)]/2
```

where η_μ(x) = (-1)^(x_1 + ... + x_(μ-1))

## Group Theory and Representations

### SU(N) Generators

The program implements SU(N) Lie algebra generators:
- For SU(2): Pauli matrices σ_i / 2
- For SU(3): Gell-Mann matrices λ_a / 2
- For general SU(N): Orthogonal traceless Hermitian matrices

### Representation Composition

The gauge field transforms in the adjoint representation of SU(N), while fermions transform in the fundamental representation. The program correctly implements:
- Fundamental representation (for fermions)
- Adjoint representation (for gauge fields)
- Tensor product decomposition for composite systems

## Advanced Features

### Tensor Network Methods

The framework is designed to support tensor network methods for quantum state representation. Future implementations could include:
- Matrix Product States (MPS)
- Projected Entangled Pair States (PEPS)
- Tree Tensor Networks (TTN)

### Quantum Computer Implementation

The code structure allows for future quantum computing implementations:
- Gate-based quantum circuits for time evolution
- Variational quantum algorithms
- Quantum phase estimation

## Performance Considerations

### OpenMP Parallelization
- Site-level parallelization for local operations
- Reduction operations for global observables
- Thread-safe gauge field updates

### MPI Parallelization
- Domain decomposition of lattice
- Halo exchange for boundary communications
- Collective operations for global measurements

## References

1. Kogut, J., & Susskind, L. (1975). Hamiltonian formulation of Wilson's lattice gauge theories. Physical Review D, 11(2), 395.

2. Wilson, K. G. (1974). Confinement of quarks. Physical review D, 10(8), 2445.

3. Rothe, H. J. (2012). Lattice gauge theories: an introduction (Vol. 43). World Scientific Publishing Company.

4. Gattringer, C., & Lang, C. B. (2010). Quantum chromodynamics on the lattice: an introductory presentation (Vol. 788). Springer Science & Business Media.

## Contact and Support

For questions, issues, or contributions, please refer to the repository:
https://github.com/KeiTohme/HamiltonianLQCD
