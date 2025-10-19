# Hamiltonian LQCD Usage Guide

## Quick Start

### 1. Build the Program

Basic build without parallelization:
```bash
make
```

Build with OpenMP support:
```bash
make openmp
```

Build with MPI support (requires MPI library):
```bash
make mpi
```

Build with both OpenMP and MPI:
```bash
make hybrid
```

### 2. Run a Test Simulation

```bash
./hamiltonian_lqcd examples/input_small_test.dat test_output.dat
```

This runs a small 2D test case that completes quickly.

## Input File Configuration

### Example Input File

```
# Physical parameters
gauge_coupling = 1.0        # Gauge coupling constant g
fermion_mass = 0.1          # Fermion mass m
num_colors = 3              # Number of colors Nc (SU(3) for QCD)

# Lattice parameters
spacetime_dim = 4           # Space-time dimension D_euclid
lattice_size = 8,8,8,8      # Lattice size Nsize for each dimension
lattice_type = square       # Lattice type: 'square' or 'hexagonal'

# Fermion type
fermion_type = wilson       # Fermion type: 'wilson' or 'staggered'

# Time evolution parameters
time_step = 0.01            # Time step dt
num_time_steps = 100        # Number of time steps

# Parallelization options
use_openmp = true           # Use OpenMP parallelization
use_mpi = false             # Use MPI parallelization
num_threads = 4             # Number of OpenMP threads

# Computational method
comp_method = standard      # Method: 'standard', 'tensor_network', or 'quantum'
```

### Parameter Details

#### Physical Parameters

**gauge_coupling** (or **g**)
- Type: Real number
- Default: 1.0
- Description: The gauge coupling constant that determines the strength of the gauge interaction
- Typical range: 0.1 - 2.0

**fermion_mass** (or **m**)
- Type: Real number
- Default: 0.1
- Description: The mass of the fermion field in lattice units
- Note: Small masses require larger lattices to avoid finite-size effects

**num_colors** (or **Nc**)
- Type: Integer
- Default: 3
- Description: Number of colors in SU(Nc) gauge group
- Common values: 2 (SU(2)), 3 (SU(3) for QCD)

#### Lattice Parameters

**spacetime_dim** (or **D_euclid**)
- Type: Integer
- Default: 4
- Description: Number of space-time dimensions
- Note: 4D for realistic QCD, 2D or 3D for testing

**lattice_size** (or **Nsize**)
- Type: Comma-separated integers
- Default: 8,8,8,8
- Description: Number of sites in each dimension
- Format: N1,N2,N3,...,ND where D = spacetime_dim
- Example: "16,16,16,32" for 16³×32 lattice

**lattice_type**
- Type: String
- Default: 'square'
- Options:
  - 'square': Hypercubic lattice with periodic boundaries
  - 'hexagonal': Hexagonal lattice (2D only)

#### Fermion Parameters

**fermion_type**
- Type: String
- Default: 'wilson'
- Options:
  - 'wilson': Wilson fermions with improved locality
  - 'staggered': Staggered (Kogut-Susskind) fermions

#### Time Evolution Parameters

**time_step** (or **dt**)
- Type: Real number
- Default: 0.01
- Description: Time step for evolution in lattice units
- Note: Smaller values give more accurate but slower evolution

**num_time_steps**
- Type: Integer
- Default: 100
- Description: Total number of time steps to simulate
- Total time = time_step × num_time_steps

#### Parallelization Parameters

**use_openmp**
- Type: Boolean (true/false)
- Default: false
- Description: Enable OpenMP shared-memory parallelization
- Note: Requires compilation with OpenMP support

**use_mpi**
- Type: Boolean (true/false)
- Default: false
- Description: Enable MPI distributed-memory parallelization
- Note: Requires MPI library and compilation with MPI support

**num_threads**
- Type: Integer
- Default: 1
- Description: Number of OpenMP threads to use
- Note: Only relevant if use_openmp = true

#### Computational Method

**comp_method**
- Type: String
- Default: 'standard'
- Options:
  - 'standard': Classical CPU-based simulation
  - 'tensor_network': Framework for tensor network methods
  - 'quantum': Framework for quantum computing implementation

## Running Simulations

### Standard Run

```bash
./hamiltonian_lqcd input.dat output.dat
```

If output filename is not specified, results are saved to `time_evolution.dat`.

### Running with MPI

```bash
mpirun -np 4 ./hamiltonian_lqcd input.dat output.dat
```

Make sure to compile with MPI support first:
```bash
make mpi
```

### Setting OpenMP Threads

If using OpenMP, you can control the number of threads via environment variable:

```bash
export OMP_NUM_THREADS=8
./hamiltonian_lqcd input.dat output.dat
```

Or set it in the input file:
```
use_openmp = true
num_threads = 8
```

## Output Interpretation

### Console Output

The program prints:
1. **Banner and initialization**: Program start, parameter summary
2. **Lattice setup**: Lattice type and size confirmation
3. **Gauge field initialization**: SU(N) group information
4. **Initial Hamiltonian**: Electric and magnetic energy components
5. **Time evolution progress**: Step number, time, energy, and average plaquette
6. **Results**: File save confirmation and final status

### Output Data File

The output file contains columns:
1. **Step**: Time step number (0 to num_time_steps)
2. **Time**: Physical time in lattice units
3. **Energy**: Total Hamiltonian energy
4. **Avg_Plaquette**: Average plaquette value (normalized)

Example:
```
# Time evolution results
# Step, Time, Energy, Avg_Plaquette
     0  0.00000000E+00 -0.32000000E+02  0.10000000E+01
     1  0.05000000E-01 -0.32000000E+02  0.10000000E+01
     2  0.10000000E+00 -0.32000000E+02  0.99987453E+00
     ...
```

### Analyzing Results

You can plot the results using tools like Python, gnuplot, or MATLAB:

**Python example:**
```python
import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt('time_evolution.dat')
time = data[:, 1]
energy = data[:, 2]
plaquette = data[:, 3]

plt.figure(figsize=(12, 4))

plt.subplot(131)
plt.plot(time, energy)
plt.xlabel('Time')
plt.ylabel('Energy')
plt.title('Energy Evolution')

plt.subplot(132)
plt.plot(time, plaquette)
plt.xlabel('Time')
plt.ylabel('Average Plaquette')
plt.title('Plaquette Evolution')

plt.subplot(133)
plt.plot(time, energy - energy[0])
plt.xlabel('Time')
plt.ylabel('Energy Change')
plt.title('Energy Conservation')

plt.tight_layout()
plt.savefig('results.png')
plt.show()
```

**Gnuplot example:**
```gnuplot
set terminal png size 1200,400
set output 'results.png'
set multiplot layout 1,3

set title "Energy Evolution"
plot 'time_evolution.dat' u 2:3 w l title 'Energy'

set title "Plaquette Evolution"
plot 'time_evolution.dat' u 2:4 w l title 'Avg Plaquette'

set title "Energy Conservation"
plot 'time_evolution.dat' u 2:(column(3)-column(3)[1]) w l title 'ΔE'

unset multiplot
```

## Physical Observables

### Energy Components

- **Electric Energy**: Kinetic energy from gauge field momenta (E field)
- **Magnetic Energy**: Potential energy from plaquette action (B field)
- **Total Energy**: Should be approximately conserved in real-time evolution

### Plaquette

The average plaquette measures the average of gauge link products around elementary squares:
- Value of 1.0: Perfect alignment (weak field limit)
- Lower values: Stronger field fluctuations
- Useful for monitoring thermalization and field configuration

### Wilson Loops

Wilson loops of size R×T measure the potential between static quarks:
- Related to confinement in QCD
- Can be computed for analysis of confinement properties

## Performance Tips

### Lattice Size

- Start with small lattices (4⁴ or 8⁴) for testing
- Computational cost scales as V × Nsteps where V is the lattice volume
- Memory usage scales as V × Nc² × Ndim

### Time Step

- Smaller time steps give better energy conservation
- Typical values: 0.01 - 0.05 in lattice units
- Check energy conservation to validate time step choice

### Parallelization

- **Small lattices (< 16⁴)**: OpenMP is sufficient
- **Medium lattices (16⁴ - 32⁴)**: OpenMP + MPI hybrid
- **Large lattices (> 32⁴)**: MPI with domain decomposition

Optimal thread count typically equals number of physical cores.

## Troubleshooting

### Build Issues

**Problem**: Cannot find gfortran
```
Solution: Install gfortran
- Ubuntu/Debian: sudo apt-get install gfortran
- macOS: brew install gcc
- Windows: Install MinGW or WSL
```

**Problem**: OpenMP not working
```
Solution: Compile with OpenMP flag
make openmp
```

**Problem**: MPI not found
```
Solution: Install MPI library
- Ubuntu/Debian: sudo apt-get install libopenmpi-dev
- macOS: brew install open-mpi
Then compile with: make mpi
```

### Runtime Issues

**Problem**: Program crashes or gives incorrect results
```
Check:
1. Lattice size is not too large for available memory
2. Time step is not too large (try reducing dt)
3. Input file format is correct (no syntax errors)
```

**Problem**: Poor performance
```
Solutions:
1. Enable OpenMP: make openmp
2. Reduce lattice size or number of time steps
3. Use optimized compiler flags (-O3)
4. For large simulations, use MPI parallelization
```

**Problem**: Energy not conserved
```
Solutions:
1. Reduce time step (dt)
2. Check numerical stability
3. Increase lattice spacing (reduce coupling)
```

## Advanced Usage

### Custom Initial Conditions

Modify the `initialize_gauge_fields()` function in `su_n_mod.f90` to set custom initial gauge configurations.

### Additional Observables

Extend the `measure_observables()` function in `time_evolution_mod.f90` to compute additional quantities.

### Different Actions

The code currently uses the Wilson gauge action. Other actions (Symanzik, Iwasaki, etc.) can be implemented by modifying `hamiltonian_mod.f90`.

## Example Workflows

### 1. Thermalization Study

```bash
# Start with random configuration
# Evolve for sufficient time
./hamiltonian_lqcd examples/input_square_lattice.dat therm_output.dat

# Plot energy vs time to check thermalization
# Expect energy to reach equilibrium after some time
```

### 2. Confinement Study

```bash
# Run multiple simulations with different couplings
# Measure Wilson loops
# Analyze potential V(R) ~ σR + const
```

### 3. Scaling Study

```bash
# Run on different lattice sizes
# Compare physical observables
# Check for finite-size effects
```

## References and Further Reading

1. **Lattice QCD introductions**:
   - Rothe, H. J. (2012). "Lattice Gauge Theories: An Introduction"
   - Gattringer & Lang (2010). "Quantum Chromodynamics on the Lattice"

2. **Hamiltonian formulation**:
   - Kogut & Susskind (1975). "Hamiltonian formulation of Wilson's lattice gauge theories"

3. **Numerical methods**:
   - Hairer et al. (2006). "Geometric Numerical Integration"

4. **Advanced topics**:
   - Tensor networks for gauge theories
   - Quantum computing for lattice field theory

## Support

For questions, issues, or contributions:
- GitHub Issues: https://github.com/KeiTohme/HamiltonianLQCD/issues
- Documentation: See docs/ directory in the repository
