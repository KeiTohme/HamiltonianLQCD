# Developer Guide for Hamiltonian LQCD

## Code Architecture

### Module Structure

The code is organized in a modular fashion with clear separation of concerns:

```
src/
├── main.f90                      # Main program driver
├── modules/
│   └── parameters_mod.f90        # Parameter management and I/O
├── lattice/
│   └── lattice_mod.f90          # Lattice structure and neighbors
├── gauge/
│   ├── su_n_mod.f90             # SU(N) gauge group operations
│   └── hamiltonian_mod.f90      # Hamiltonian and observables
├── fermions/
│   ├── wilson_fermion_mod.f90   # Wilson fermion formulation
│   └── staggered_fermion_mod.f90 # Staggered fermion formulation
├── time_evolution/
│   └── time_evolution_mod.f90   # Real-time evolution algorithms
└── parallel/
    ├── openmp_mod.f90           # OpenMP parallelization
    └── mpi_mod.f90              # MPI parallelization
```

### Module Dependencies

```
main.f90
  ├─> parameters_mod
  ├─> lattice_mod
  │     └─> parameters_mod
  ├─> su_n_mod
  │     ├─> parameters_mod
  │     └─> lattice_mod
  ├─> wilson_fermion_mod
  │     ├─> parameters_mod
  │     ├─> lattice_mod
  │     └─> su_n_mod
  ├─> staggered_fermion_mod
  │     ├─> parameters_mod
  │     ├─> lattice_mod
  │     └─> su_n_mod
  ├─> hamiltonian_mod
  │     ├─> parameters_mod
  │     ├─> lattice_mod
  │     └─> su_n_mod
  ├─> time_evolution_mod
  │     ├─> parameters_mod
  │     ├─> lattice_mod
  │     ├─> su_n_mod
  │     └─> hamiltonian_mod
  └─> parallel modules (openmp_mod, mpi_mod)
        └─> parameters_mod
```

## Key Modules

### parameters_mod.f90

**Purpose**: Manage simulation parameters and configuration

**Key Components**:
- `read_parameters(filename)`: Read input file
- `print_parameters()`: Display current parameters
- `cleanup_parameters()`: Free memory

**Global Variables**:
- Physical: `gauge_coupling`, `fermion_mass`, `num_colors`
- Lattice: `spacetime_dim`, `lattice_size`, `total_sites`, `lattice_type`
- Evolution: `time_step`, `num_time_steps`
- Parallel: `use_openmp`, `use_mpi`, `num_threads`

### lattice_mod.f90

**Purpose**: Define lattice structure and neighbor relations

**Key Types**:
```fortran
type :: lattice_site
    integer :: site_index
    integer, allocatable :: coords(:)
    integer, allocatable :: neighbors(:,:)
end type lattice_site
```

**Key Functions**:
- `initialize_lattice()`: Set up lattice structure
- `initialize_square_lattice()`: Hypercubic lattice
- `initialize_hexagonal_lattice()`: Hexagonal lattice (2D)
- `index_to_coords(idx, coords)`: Convert linear index to coordinates
- `coords_to_index(coords)`: Convert coordinates to linear index

### su_n_mod.f90

**Purpose**: SU(N) gauge group operations

**Key Types**:
```fortran
type :: su_n_matrix
    complex(dp), allocatable :: elements(:,:)
end type su_n_matrix
```

**Global Arrays**:
- `gauge_links(site, direction)`: Gauge link variables
- `gauge_momenta(site, direction)`: Conjugate momenta (E fields)
- `su_n_generators(a)`: Lie algebra generators

**Key Functions**:
- `initialize_gauge_fields()`: Initialize gauge configuration
- `initialize_su_n_generators()`: Set up generators for SU(2), SU(3), or general SU(N)
- `matrix_multiply(A, B)`: Matrix multiplication
- `matrix_dagger(A)`: Hermitian conjugate
- `trace(A)`: Matrix trace

### wilson_fermion_mod.f90

**Purpose**: Wilson fermion discretization

**Key Arrays**:
- `psi(site, spin, color)`: Fermion field
- `psi_bar(site, spin, color)`: Conjugate fermion field

**Key Functions**:
- `apply_wilson_dirac_operator(psi_in, psi_out)`: Apply D_Wilson
- `apply_projector_minus/plus(mu, spinor)`: Gamma matrix projectors
- `apply_gauge_link(site, mu, spinor)`: Gauge covariant derivative

### staggered_fermion_mod.f90

**Purpose**: Staggered fermion discretization

**Key Arrays**:
- `chi(site, color)`: Staggered fermion field (one component per site)
- `chi_bar(site, color)`: Conjugate field

**Key Functions**:
- `apply_staggered_dirac_operator(chi_in, chi_out)`: Apply D_staggered
- `staggered_phase(site, mu)`: Compute eta_mu(x) phase factors

### hamiltonian_mod.f90

**Purpose**: Hamiltonian formulation and observables

**Key Functions**:
- `compute_hamiltonian()`: Total Hamiltonian
- `compute_electric_energy()`: Kinetic term from E fields
- `compute_magnetic_energy()`: Potential term from plaquettes
- `compute_plaquette(site, mu, nu)`: Single plaquette value
- `compute_wilson_loops(R, T, value)`: Wilson loop observables

### time_evolution_mod.f90

**Purpose**: Real-time evolution using Hamiltonian dynamics

**Key Functions**:
- `evolve_gauge_fields()`: Main evolution loop
- `symplectic_step()`: Leapfrog integrator step
- `update_momenta(dt)`: Update E fields
- `update_links(dt)`: Update U fields
- `compute_gauge_force(site, mu, force)`: Force from gauge action
- `reunitarize(U)`: Restore SU(N) unitarity

## Adding New Features

### Adding a New Observable

1. **Define the observable in `hamiltonian_mod.f90`**:

```fortran
function compute_new_observable() result(obs_value)
    implicit none
    real(dp) :: obs_value
    integer :: site, mu, nu
    
    obs_value = 0.0_dp
    
    ! Your computation here
    do site = 1, total_sites
        ! ...
    end do
    
end function compute_new_observable
```

2. **Add measurement in `time_evolution_mod.f90`**:

```fortran
! In measure_observables subroutine
real(dp) :: new_obs
new_obs = compute_new_observable()
! Store or output new_obs
```

3. **Update output if needed**:

```fortran
! In save_results subroutine
write(unit_num, '(I6,4E16.8)') i-1, (i-1)*time_step, &
    energy_history(i), plaquette_history(i), new_obs_history(i)
```

### Adding a New Fermion Type

1. **Create new module** `src/fermions/new_fermion_mod.f90`:

```fortran
module new_fermion_mod
    use parameters_mod
    use lattice_mod
    use su_n_mod
    implicit none
    
    ! Define fermion fields
    complex(dp), allocatable :: fermion_field(:,:,:)
    
contains

    subroutine initialize_new_fermions()
        ! Allocate and initialize fields
    end subroutine initialize_new_fermions
    
    subroutine apply_new_dirac_operator(psi_in, psi_out)
        ! Implement Dirac operator
    end subroutine apply_new_dirac_operator
    
    subroutine cleanup_new_fermions()
        ! Free memory
    end subroutine cleanup_new_fermions

end module new_fermion_mod
```

2. **Add to main program** `src/main.f90`:

```fortran
use new_fermion_mod

! In initialization section
select case (trim(fermion_type))
    case ('wilson')
        call initialize_wilson_fermions()
    case ('staggered')
        call initialize_staggered_fermions()
    case ('new_type')
        call initialize_new_fermions()
end select
```

3. **Update Makefile**:

```makefile
FERMION_SOURCES = $(FERMION_DIR)/wilson_fermion_mod.f90 \
                  $(FERMION_DIR)/staggered_fermion_mod.f90 \
                  $(FERMION_DIR)/new_fermion_mod.f90
```

### Adding a New Lattice Type

1. **Implement in `lattice_mod.f90`**:

```fortran
subroutine initialize_new_lattice_type()
    implicit none
    integer :: site_idx
    
    allocate(lattice(total_sites))
    
    do site_idx = 1, total_sites
        ! Set up site structure
        ! Define neighbors according to new geometry
    end do
    
end subroutine initialize_new_lattice_type
```

2. **Add case in `initialize_lattice()`**:

```fortran
select case (trim(lattice_type))
    case ('square')
        call initialize_square_lattice()
    case ('hexagonal')
        call initialize_hexagonal_lattice()
    case ('new_type')
        call initialize_new_lattice_type()
end select
```

### Implementing Improved Actions

To implement improved gauge actions (Symanzik, Iwasaki, etc.):

1. **Add new plaquette functions in `hamiltonian_mod.f90`**:

```fortran
function compute_rectangle(site, mu, nu) result(rect_value)
    ! Compute 1x2 rectangle
end function compute_rectangle

subroutine compute_improved_magnetic_energy()
    ! Use weighted sum of plaquettes and rectangles
    magnetic_energy = -c0 * plaquette_sum - c1 * rectangle_sum
end subroutine compute_improved_magnetic_energy
```

2. **Update force computation**:

```fortran
subroutine compute_improved_gauge_force(site, mu, force)
    ! Include staples from rectangles
end subroutine compute_improved_gauge_force
```

## Parallelization Guidelines

### OpenMP Parallelization

**Where to parallelize**:
- Loops over lattice sites for local operations
- Gauge force computation
- Observable measurements

**Example**:
```fortran
!$OMP PARALLEL DO PRIVATE(site, mu) REDUCTION(+:total_sum)
do site = 1, total_sites
    do mu = 1, spacetime_dim
        ! Local computation
        total_sum = total_sum + local_value
    end do
end do
!$OMP END PARALLEL DO
```

**Considerations**:
- Avoid race conditions in gauge link updates
- Use proper reduction clauses for global sums
- Consider thread-safe random number generation

### MPI Parallelization

**Domain Decomposition Strategy**:
1. Partition lattice into sub-domains
2. Assign each sub-domain to an MPI rank
3. Implement halo exchange for boundary links

**Implementation outline**:
```fortran
subroutine decompose_lattice()
    ! Divide lattice among MPI ranks
    ! Assign local sites to each rank
end subroutine decompose_lattice

subroutine exchange_boundaries()
    ! Send/receive boundary link data
    ! Use MPI_Isend/MPI_Irecv for non-blocking communication
end subroutine exchange_boundaries

subroutine gather_observables()
    ! Use MPI_Reduce for global observables
    call MPI_REDUCE(local_energy, global_energy, 1, &
                   MPI_DOUBLE_PRECISION, MPI_SUM, 0, &
                   MPI_COMM_WORLD, ierr)
end subroutine gather_observables
```

## Testing and Validation

### Unit Tests

Test individual components:

1. **Lattice structure**: Verify neighbor relations
2. **SU(N) operations**: Check unitarity and algebra
3. **Fermion operators**: Test Hermiticity properties
4. **Integrators**: Verify energy conservation

### Integration Tests

1. **Small lattices**: Run quick tests on 4⁴ or smaller
2. **Energy conservation**: Check ΔE/E < tolerance
3. **Gauge invariance**: Verify observables are gauge invariant

### Benchmark Tests

1. **Scaling tests**: Performance vs lattice size
2. **Parallel efficiency**: Speedup with OpenMP/MPI
3. **Comparison**: Validate against known results

## Coding Standards

### Fortran Style

```fortran
! Use explicit typing
implicit none

! Clear variable names
integer :: lattice_site_index
real(dp) :: gauge_coupling_constant

! Consistent indentation (4 spaces)
do i = 1, n
    if (condition) then
        ! code
    end if
end do

! Comment complex algorithms
! This implements the symplectic integrator
call update_momenta(dt/2)
call update_links(dt)
call update_momenta(dt/2)
```

### Memory Management

```fortran
! Always deallocate allocated arrays
if (allocated(array)) deallocate(array)

! Check allocation status
allocate(array(n), stat=ierr)
if (ierr /= 0) then
    print *, 'Allocation failed'
    stop
end if
```

### Error Handling

```fortran
! Check file operations
open(unit=10, file='input.dat', status='old', iostat=ierr)
if (ierr /= 0) then
    print *, 'Error opening file'
    return
end if

! Validate parameters
if (time_step <= 0.0_dp) then
    print *, 'Error: time_step must be positive'
    stop
end if
```

## Performance Optimization

### Compiler Flags

```makefile
# Optimization
FCFLAGS = -O3 -march=native -funroll-loops

# Debugging
FCFLAGS = -O0 -g -fcheck=all -fbacktrace

# Profiling
FCFLAGS = -O2 -pg
```

### Hot Spots

Profile the code to identify bottlenecks:

```bash
# Compile with profiling
make FCFLAGS="-O2 -pg"

# Run simulation
./hamiltonian_lqcd input.dat

# Analyze profile
gprof hamiltonian_lqcd gmon.out > profile.txt
```

Common hot spots:
1. Gauge force computation
2. Matrix operations
3. Reunitarization

### Memory Optimization

- Use contiguous memory layouts
- Minimize allocations inside loops
- Consider structure-of-arrays vs array-of-structures

## Future Extensions

### Tensor Network Methods

1. Implement Matrix Product State (MPS) representation
2. Add TEBD or TDVP time evolution
3. Optimize for 1D+1 or 2D+1 systems

### Quantum Computing

1. Implement gate decomposition for gauge links
2. Add quantum circuit construction
3. Interface with quantum simulators (Qiskit, Cirq)

### Advanced Features

- [ ] Adaptive time stepping
- [ ] Multiple gauge actions (Symanzik, Iwasaki)
- [ ] Fermion determinant computation
- [ ] Topology (instantons, monopoles)
- [ ] Multi-grid methods
- [ ] GPU acceleration

## Contributing

### Workflow

1. Fork the repository
2. Create a feature branch
3. Make changes with clear commits
4. Add tests for new features
5. Update documentation
6. Submit pull request

### Code Review Checklist

- [ ] Code compiles without warnings
- [ ] Tests pass
- [ ] Documentation updated
- [ ] Performance acceptable
- [ ] Memory leaks checked
- [ ] Style consistent

## References

### Books
- "Lattice Gauge Theories: An Introduction" - H.J. Rothe
- "Quantum Chromodynamics on the Lattice" - Gattringer & Lang

### Papers
- Kogut & Susskind (1975) - Hamiltonian formulation
- Wilson (1974) - Lattice gauge theory
- Nielsen & Ninomiya (1981) - Fermion doubling theorem

### Numerical Methods
- Hairer et al. - "Geometric Numerical Integration"
- Press et al. - "Numerical Recipes"

## Support

- GitHub Issues: Report bugs and request features
- Discussions: Ask questions and share ideas
- Documentation: Check docs/ directory for more info
