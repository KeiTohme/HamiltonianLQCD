# Makefile for Lattice QCD Hamiltonian Program

# Compiler settings
FC = gfortran
MPIFC = mpifort

# Compiler flags
FFLAGS = -O3 -Wall -fcheck=bounds
DEBUGFLAGS = -g -O0 -fbacktrace -fcheck=all
OMPFLAGS = -fopenmp
MPIFLAGS = -DUSE_MPI

# Directories
SRCDIR = src
BUILDDIR = build
BINDIR = .

# Source files
SOURCES = parameter_module.f90 \
          lattice_module.f90 \
          gauge_field_module.f90 \
          fermion_module.f90 \
          hamiltonian_module.f90 \
          parallel_module.f90 \
          tensor_network_module.f90 \
          quantum_computer_module.f90 \
          time_evolution_module.f90 \
          main.f90

# Object files
OBJECTS = $(SOURCES:%.f90=$(BUILDDIR)/%.o)

# Module dependencies
$(BUILDDIR)/main.o: $(BUILDDIR)/parameter_module.o \
                    $(BUILDDIR)/lattice_module.o \
                    $(BUILDDIR)/gauge_field_module.o \
                    $(BUILDDIR)/fermion_module.o \
                    $(BUILDDIR)/hamiltonian_module.o \
                    $(BUILDDIR)/parallel_module.o \
                    $(BUILDDIR)/time_evolution_module.o

$(BUILDDIR)/time_evolution_module.o: $(BUILDDIR)/parameter_module.o \
                                     $(BUILDDIR)/lattice_module.o \
                                     $(BUILDDIR)/gauge_field_module.o \
                                     $(BUILDDIR)/fermion_module.o \
                                     $(BUILDDIR)/hamiltonian_module.o

$(BUILDDIR)/hamiltonian_module.o: $(BUILDDIR)/parameter_module.o \
                                  $(BUILDDIR)/lattice_module.o \
                                  $(BUILDDIR)/gauge_field_module.o \
                                  $(BUILDDIR)/fermion_module.o

$(BUILDDIR)/tensor_network_module.o: $(BUILDDIR)/parameter_module.o \
                                    $(BUILDDIR)/lattice_module.o \
                                    $(BUILDDIR)/gauge_field_module.o \
                                    $(BUILDDIR)/fermion_module.o

$(BUILDDIR)/quantum_computer_module.o: $(BUILDDIR)/parameter_module.o \
                                      $(BUILDDIR)/lattice_module.o \
                                      $(BUILDDIR)/gauge_field_module.o \
                                      $(BUILDDIR)/fermion_module.o

$(BUILDDIR)/fermion_module.o: $(BUILDDIR)/parameter_module.o \
                              $(BUILDDIR)/lattice_module.o \
                              $(BUILDDIR)/gauge_field_module.o

$(BUILDDIR)/gauge_field_module.o: $(BUILDDIR)/parameter_module.o \
                                  $(BUILDDIR)/lattice_module.o

$(BUILDDIR)/lattice_module.o: $(BUILDDIR)/parameter_module.o

$(BUILDDIR)/parallel_module.o: $(BUILDDIR)/parameter_module.o \
                               $(BUILDDIR)/lattice_module.o \
                               $(BUILDDIR)/gauge_field_module.o \
                               $(BUILDDIR)/fermion_module.o

# Targets
.PHONY: all serial parallel mpi clean debug

# Default target - serial with OpenMP
all: serial

# Serial version (with optional OpenMP)
serial: FC := $(FC)
serial: FFLAGS += $(OMPFLAGS)
serial: $(BINDIR)/lattice_qcd_hamiltonian

# Parallel version with MPI and OpenMP
mpi: FC := $(MPIFC)
mpi: FFLAGS += $(OMPFLAGS) $(MPIFLAGS)
mpi: $(BINDIR)/lattice_qcd_hamiltonian_mpi

# Debug version
debug: FFLAGS = $(DEBUGFLAGS) $(OMPFLAGS)
debug: $(BINDIR)/lattice_qcd_hamiltonian_debug

# Build executable
$(BINDIR)/lattice_qcd_hamiltonian: $(OBJECTS)
	@mkdir -p $(BINDIR)
	$(FC) $(FFLAGS) -o $@ $^

$(BINDIR)/lattice_qcd_hamiltonian_mpi: $(OBJECTS)
	@mkdir -p $(BINDIR)
	$(FC) $(FFLAGS) -o $@ $^

$(BINDIR)/lattice_qcd_hamiltonian_debug: $(OBJECTS)
	@mkdir -p $(BINDIR)
	$(FC) $(FFLAGS) -o $@ $^

# Build objects
$(BUILDDIR)/%.o: $(SRCDIR)/%.f90
	@mkdir -p $(BUILDDIR)
	$(FC) $(FFLAGS) -c -o $@ $< -J$(BUILDDIR)

# Clean
clean:
	rm -rf $(BUILDDIR)
	rm -f $(BINDIR)/lattice_qcd_hamiltonian*
	rm -rf output/*

# Create necessary directories
dirs:
	@mkdir -p $(BUILDDIR)
	@mkdir -p output

# Run test
test: serial
	@mkdir -p output
	./lattice_qcd_hamiltonian input/parameters.inp

# Help
help:
	@echo "Available targets:"
	@echo "  make         - Build serial version with OpenMP"
	@echo "  make mpi     - Build MPI+OpenMP version"
	@echo "  make debug   - Build debug version"
	@echo "  make test    - Run test with sample parameters"
	@echo "  make clean   - Clean build files"
	@echo "  make help    - Show this help"