# =============================================================================
# Makefile for Lattice QCD Hamiltonian Formalism
# =============================================================================

# Compiler options
FC = gfortran
FFLAGS = -O3 -Wall -fcheck=bounds -fbacktrace
LDFLAGS = 

# Directories
SRCDIR = src
MODDIR = src/modules
OBJDIR = build
BINDIR = .

# Executable name
EXEC = lattice_qcd

# Optional features
USE_OPENMP = yes
USE_MPI = no

# OpenMP flags
ifeq ($(USE_OPENMP),yes)
    FFLAGS += -fopenmp
    LDFLAGS += -fopenmp
endif

# MPI flags
ifeq ($(USE_MPI),yes)
    FC = mpif90
    FFLAGS += -DUSE_MPI
endif

# Source files
MODULES = \
    $(MODDIR)/mod_parameters.f90 \
    $(MODDIR)/mod_su_algebra.f90 \
    $(MODDIR)/mod_lattice.f90 \
    $(MODDIR)/mod_gauge_field.f90 \
    $(MODDIR)/mod_wilson_fermion.f90 \
    $(MODDIR)/mod_staggered_fermion.f90 \
    $(MODDIR)/mod_hamiltonian.f90 \
    $(MODDIR)/mod_time_evolution.f90 \
    $(MODDIR)/mod_parallel.f90 \
    $(MODDIR)/mod_tensor_network.f90 \
    $(MODDIR)/mod_quantum_interface.f90

MAIN = $(SRCDIR)/main.f90

# Object files
MOD_OBJS = $(patsubst $(MODDIR)/%.f90,$(OBJDIR)/%.o,$(MODULES))
MAIN_OBJ = $(OBJDIR)/main.o

ALL_OBJS = $(MOD_OBJS) $(MAIN_OBJ)

# =============================================================================
# Targets
# =============================================================================

.PHONY: all clean dirs help

all: dirs $(BINDIR)/$(EXEC)

# Create necessary directories
dirs:
	@mkdir -p $(OBJDIR)
	@mkdir -p output

# Link executable
$(BINDIR)/$(EXEC): $(ALL_OBJS)
	@echo "Linking $@"
	$(FC) $(LDFLAGS) -o $@ $(ALL_OBJS)
	@echo "Build successful!"

# Compile modules (with dependencies)
$(OBJDIR)/mod_parameters.o: $(MODDIR)/mod_parameters.f90
	@echo "Compiling $<"
	$(FC) $(FFLAGS) -J$(OBJDIR) -c $< -o $@

$(OBJDIR)/mod_su_algebra.o: $(MODDIR)/mod_su_algebra.f90 $(OBJDIR)/mod_parameters.o
	@echo "Compiling $<"
	$(FC) $(FFLAGS) -J$(OBJDIR) -c $< -o $@

$(OBJDIR)/mod_lattice.o: $(MODDIR)/mod_lattice.f90 $(OBJDIR)/mod_parameters.o
	@echo "Compiling $<"
	$(FC) $(FFLAGS) -J$(OBJDIR) -c $< -o $@

$(OBJDIR)/mod_gauge_field.o: $(MODDIR)/mod_gauge_field.f90 \
    $(OBJDIR)/mod_parameters.o $(OBJDIR)/mod_lattice.o $(OBJDIR)/mod_su_algebra.o
	@echo "Compiling $<"
	$(FC) $(FFLAGS) -J$(OBJDIR) -c $< -o $@

$(OBJDIR)/mod_wilson_fermion.o: $(MODDIR)/mod_wilson_fermion.f90 \
    $(OBJDIR)/mod_parameters.o $(OBJDIR)/mod_lattice.o \
    $(OBJDIR)/mod_gauge_field.o $(OBJDIR)/mod_su_algebra.o
	@echo "Compiling $<"
	$(FC) $(FFLAGS) -J$(OBJDIR) -c $< -o $@

$(OBJDIR)/mod_staggered_fermion.o: $(MODDIR)/mod_staggered_fermion.f90 \
    $(OBJDIR)/mod_parameters.o $(OBJDIR)/mod_lattice.o \
    $(OBJDIR)/mod_gauge_field.o $(OBJDIR)/mod_su_algebra.o
	@echo "Compiling $<"
	$(FC) $(FFLAGS) -J$(OBJDIR) -c $< -o $@

$(OBJDIR)/mod_hamiltonian.o: $(MODDIR)/mod_hamiltonian.f90 \
    $(OBJDIR)/mod_parameters.o $(OBJDIR)/mod_lattice.o \
    $(OBJDIR)/mod_gauge_field.o $(OBJDIR)/mod_su_algebra.o \
    $(OBJDIR)/mod_wilson_fermion.o $(OBJDIR)/mod_staggered_fermion.o
	@echo "Compiling $<"
	$(FC) $(FFLAGS) -J$(OBJDIR) -c $< -o $@

$(OBJDIR)/mod_time_evolution.o: $(MODDIR)/mod_time_evolution.f90 \
    $(OBJDIR)/mod_parameters.o $(OBJDIR)/mod_lattice.o \
    $(OBJDIR)/mod_gauge_field.o $(OBJDIR)/mod_hamiltonian.o \
    $(OBJDIR)/mod_su_algebra.o
	@echo "Compiling $<"
	$(FC) $(FFLAGS) -J$(OBJDIR) -c $< -o $@

$(OBJDIR)/mod_parallel.o: $(MODDIR)/mod_parallel.f90 $(OBJDIR)/mod_parameters.o
	@echo "Compiling $<"
	$(FC) $(FFLAGS) -J$(OBJDIR) -c $< -o $@

$(OBJDIR)/mod_tensor_network.o: $(MODDIR)/mod_tensor_network.f90 \
    $(OBJDIR)/mod_parameters.o $(OBJDIR)/mod_lattice.o $(OBJDIR)/mod_gauge_field.o
	@echo "Compiling $<"
	$(FC) $(FFLAGS) -J$(OBJDIR) -c $< -o $@

$(OBJDIR)/mod_quantum_interface.o: $(MODDIR)/mod_quantum_interface.f90 \
    $(OBJDIR)/mod_parameters.o $(OBJDIR)/mod_lattice.o \
    $(OBJDIR)/mod_gauge_field.o $(OBJDIR)/mod_su_algebra.o
	@echo "Compiling $<"
	$(FC) $(FFLAGS) -J$(OBJDIR) -c $< -o $@

# Compile main program
$(OBJDIR)/main.o: $(SRCDIR)/main.f90 $(MOD_OBJS)
	@echo "Compiling $<"
	$(FC) $(FFLAGS) -J$(OBJDIR) -c $< -o $@

# Clean build files
clean:
	@echo "Cleaning build files"
	rm -rf $(OBJDIR)/*.o $(OBJDIR)/*.mod
	rm -f $(BINDIR)/$(EXEC)
	@echo "Clean complete"

# Help message
help:
	@echo "Lattice QCD Hamiltonian Formalism - Build System"
	@echo ""
	@echo "Available targets:"
	@echo "  make          - Build the program (default)"
	@echo "  make clean    - Remove build files"
	@echo "  make help     - Show this help message"
	@echo ""
	@echo "Compilation options:"
	@echo "  USE_OPENMP=yes/no  - Enable/disable OpenMP (default: yes)"
	@echo "  USE_MPI=yes/no     - Enable/disable MPI (default: no)"
	@echo ""
	@echo "Examples:"
	@echo "  make                          - Build with OpenMP"
	@echo "  make USE_OPENMP=no            - Build without OpenMP"
	@echo "  make USE_MPI=yes              - Build with MPI"
	@echo "  make USE_OPENMP=yes USE_MPI=yes - Build with both"
	@echo ""
	@echo "Running the program:"
	@echo "  ./lattice_qcd [input_file]"
	@echo "  Default input file: input/parameters.inp"
