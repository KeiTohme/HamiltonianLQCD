# Lattice QCD Hamiltonian Formalism Program

This is a Fortran 90/95 implementation of lattice QCD in the Hamiltonian formalism, designed to study real-time gluon dynamics.

## Features

- **Lattice Geometries**: Square and hexagonal lattices
- **Fermion Types**: Wilson and Staggered fermions
- **Gauge Groups**: SU(N) with arbitrary N (optimized for SU(2) and SU(3))
- **Parallelization**: OpenMP and MPI support
- **Time Evolution**: Leapfrog integrator (with hooks for RK4 and symplectic methods)
- **Interfaces**: Prepared for tensor network and quantum computer backends

## Directory Structure

```
lattice_qcd/
├── src/                    # Source code
│   ├── main.f90           # Main program
│   ├── parameter_module.f90    # Parameter management
│   ├── lattice_module.f90      # Lattice geometry
│   ├── gauge_field_module.f90  # Gauge field operations
│   ├── fermion_module.f90      # Fermion implementations
│   ├── hamiltonian_module.f90  # Hamiltonian construction
│   ├── time_evolution_module.f90 # Time evolution algorithms
│   └── parallel_module.f90     # Parallelization support
├── input/                  # Input parameter files
│   └── parameters.inp     # Sample parameter file
├── output/                # Output directory
├── build/                 # Build directory (created by make)
├── Makefile              # Build configuration
└── README.md             # This file
```

## Building the Program

### Requirements

- Fortran 90/95 compiler (gfortran, ifort, etc.)
- OpenMP support (optional)
- MPI implementation (optional, e.g., OpenMPI, MPICH)

### Compilation

```bash
# Serial version with OpenMP
make

# MPI+OpenMP version
make mpi

# Debug version
make debug

# Clean build
make clean
```

## Running the Program

```bash
# Serial/OpenMP version
./lattice_qcd_hamiltonian input/parameters.inp

# MPI version
mpirun -np 4 ./lattice_qcd_hamiltonian_mpi input/parameters.inp
```

## Input Parameters

The program reads parameters from an input file. See `input/parameters.inp` for an example.

### Key Parameters:

- `g`: Gauge coupling constant
- `m_fermion`: Fermion mass
- `Nc`: Number of colors (e.g., 3 for QCD)
- `D_euclid`: Spacetime dimension (2, 3, or 4)
- `Nsize`: Lattice size in each dimension
- `lattice_type`: "square" or "hexagonal" (hexagonal only for 2D)
- `fermion_type`: "wilson" or "staggered"
- `dt`: Time step for evolution
- `n_steps`: Number of time evolution steps

### Parallelization Options:

- `use_mpi`: Enable MPI parallelization
- `use_openmp`: Enable OpenMP parallelization
- `n_threads`: Number of OpenMP threads

## Physics Background

### Hamiltonian Formalism

The program implements the Hamiltonian formulation of lattice QCD:

```
H = H_gauge + H_fermion
```

where:
- `H_gauge`: Pure gauge Hamiltonian (electric + magnetic energy)
- `H_fermion`: Fermion kinetic term and gauge-fermion interaction

### Gauge Field Dynamics

The gauge fields are represented by:
- Link variables: U_μ(x) ∈ SU(N)
- Electric field: E_μ(x) (conjugate momentum to A_μ)

### Fermion Implementations

1. **Wilson Fermions**: 
   - Removes fermion doubling
   - Breaks chiral symmetry explicitly
   - 4-component spinors

2. **Staggered Fermions**:
   - Reduces fermion doubling
   - Preserves remnant chiral symmetry
   - 1-component per site

## Output

The program outputs:
- Energy evolution: `output/energy.dat`
- Gauge configurations: `output/gauge_config_*.dat`
- Fermion configurations: `output/fermion_config_*.dat`

## Extending the Code

### Adding Tensor Network Support

The code has hooks for tensor network methods in:
- `hamiltonian_module.f90`: `setup_tensor_network_hamiltonian`
- Implement actual TN backend integration as needed

### Adding Quantum Computer Support

Quantum computer interface prepared in:
- `hamiltonian_module.f90`: `setup_quantum_computer_hamiltonian`
- Implement quantum circuit generation as needed

### Adding New Observables

Add measurement routines in:
- `time_evolution_module.f90`: `measure_observables`

## Mathematical Details

### SU(N) Generators

The code implements general SU(N) generators:
- SU(2): Pauli matrices σᵢ/2
- SU(3): Gell-Mann matrices λₐ/2
- SU(N): Generalized construction

### Group Integration

Link updates use exponential map:
```
U_μ(x) → exp(ig dt E_μ(x)) U_μ(x)
```

with reunitarization to maintain SU(N) constraint.

## Performance Considerations

- OpenMP parallelization over lattice sites
- MPI domain decomposition for large lattices
- Optimized matrix operations for small N

## Known Limitations

- Hexagonal lattice only implemented for 2D
- Tensor network and quantum computer interfaces are placeholders
- Some advanced integrators (RK4, Symplectic-4) not fully implemented

## References

1. Kogut, J. B., & Susskind, L. (1975). "Hamiltonian formulation of Wilson's lattice gauge theories"
2. Creutz, M. (1983). "Quarks, gluons and lattices"
3. Rothe, H. J. (2012). "Lattice gauge theories: an introduction"

## License

This code is provided for educational and research purposes.