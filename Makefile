# ==============================================================================
# Makefile for Hamiltonian LQCD
# ==============================================================================

# Compiler settings
FC = gfortran
FCFLAGS = -O2 -Wall -Wextra -fcheck=bounds
LDFLAGS = 

# Optional flags
OPENMP_FLAG = -fopenmp
MPI_FC = mpif90
MPI_FLAG = -DUSE_MPI

# Target executable
TARGET = hamiltonian_lqcd

# Source directories
SRC_DIR = src
MOD_DIR = src/modules
LATTICE_DIR = src/lattice
GAUGE_DIR = src/gauge
FERMION_DIR = src/fermions
TIME_DIR = src/time_evolution
PARALLEL_DIR = src/parallel

# Object directory
OBJ_DIR = obj

# Module files
MOD_SOURCES = $(MOD_DIR)/parameters_mod.f90

# Lattice files
LATTICE_SOURCES = $(LATTICE_DIR)/lattice_mod.f90

# Gauge files
GAUGE_SOURCES = $(GAUGE_DIR)/su_n_mod.f90 \
                $(GAUGE_DIR)/hamiltonian_mod.f90

# Fermion files
FERMION_SOURCES = $(FERMION_DIR)/wilson_fermion_mod.f90 \
                  $(FERMION_DIR)/staggered_fermion_mod.f90

# Time evolution files
TIME_SOURCES = $(TIME_DIR)/time_evolution_mod.f90

# Parallel files
PARALLEL_SOURCES = $(PARALLEL_DIR)/openmp_mod.f90 \
                   $(PARALLEL_DIR)/mpi_mod.f90

# Main program
MAIN_SOURCE = $(SRC_DIR)/main.f90

# All sources in compilation order
SOURCES = $(MOD_SOURCES) \
          $(LATTICE_SOURCES) \
          $(GAUGE_SOURCES) \
          $(FERMION_SOURCES) \
          $(TIME_SOURCES) \
          $(PARALLEL_SOURCES) \
          $(MAIN_SOURCE)

# Object files
OBJECTS = $(patsubst $(SRC_DIR)/%.f90,$(OBJ_DIR)/%.o,$(SOURCES))

# Build options
USE_OPENMP ?= yes
USE_MPI ?= no

# Add OpenMP flags if enabled
ifeq ($(USE_OPENMP),yes)
    FCFLAGS += $(OPENMP_FLAG) -D_OPENMP
endif

# Use MPI compiler if enabled
ifeq ($(USE_MPI),yes)
    FC = $(MPI_FC)
    FCFLAGS += $(MPI_FLAG)
endif

# Default target
all: directories $(TARGET)

# Create necessary directories
directories:
	@mkdir -p $(OBJ_DIR)
	@mkdir -p $(OBJ_DIR)/modules
	@mkdir -p $(OBJ_DIR)/lattice
	@mkdir -p $(OBJ_DIR)/gauge
	@mkdir -p $(OBJ_DIR)/fermions
	@mkdir -p $(OBJ_DIR)/time_evolution
	@mkdir -p $(OBJ_DIR)/parallel

# Link executable
$(TARGET): $(OBJECTS)
	$(FC) $(FCFLAGS) -o $@ $^ $(LDFLAGS)
	@echo ""
	@echo "Build successful! Executable: $(TARGET)"
	@echo ""

# Compile object files
$(OBJ_DIR)/modules/%.o: $(MOD_DIR)/%.f90
	$(FC) $(FCFLAGS) -c $< -o $@

$(OBJ_DIR)/lattice/%.o: $(LATTICE_DIR)/%.f90
	$(FC) $(FCFLAGS) -c $< -o $@

$(OBJ_DIR)/gauge/%.o: $(GAUGE_DIR)/%.f90
	$(FC) $(FCFLAGS) -c $< -o $@

$(OBJ_DIR)/fermions/%.o: $(FERMION_DIR)/%.f90
	$(FC) $(FCFLAGS) -c $< -o $@

$(OBJ_DIR)/time_evolution/%.o: $(TIME_DIR)/%.f90
	$(FC) $(FCFLAGS) -c $< -o $@

$(OBJ_DIR)/parallel/%.o: $(PARALLEL_DIR)/%.f90
	$(FC) $(FCFLAGS) -c $< -o $@

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.f90
	$(FC) $(FCFLAGS) -c $< -o $@

# Build with OpenMP
openmp:
	$(MAKE) USE_OPENMP=yes

# Build with MPI
mpi:
	$(MAKE) USE_MPI=yes

# Build with both OpenMP and MPI
hybrid:
	$(MAKE) USE_OPENMP=yes USE_MPI=yes

# Clean build files
clean:
	rm -rf $(OBJ_DIR)
	rm -f *.mod
	rm -f $(TARGET)
	@echo "Clean complete"

# Clean and rebuild
rebuild: clean all

# Help message
help:
	@echo "Hamiltonian LQCD Makefile"
	@echo ""
	@echo "Targets:"
	@echo "  all     - Build with default settings"
	@echo "  openmp  - Build with OpenMP support"
	@echo "  mpi     - Build with MPI support"
	@echo "  hybrid  - Build with both OpenMP and MPI"
	@echo "  clean   - Remove build files"
	@echo "  rebuild - Clean and rebuild"
	@echo "  help    - Show this help message"
	@echo ""
	@echo "Options:"
	@echo "  USE_OPENMP=yes/no  (default: yes)"
	@echo "  USE_MPI=yes/no     (default: no)"
	@echo ""
	@echo "Examples:"
	@echo "  make"
	@echo "  make openmp"
	@echo "  make USE_OPENMP=no"
	@echo "  make mpi"
	@echo "  make hybrid"

.PHONY: all directories openmp mpi hybrid clean rebuild help
