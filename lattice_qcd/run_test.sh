#!/bin/bash

# Test script for Lattice QCD Hamiltonian program

echo "=== Lattice QCD Hamiltonian Test Script ==="
echo

# Create output directory if it doesn't exist
mkdir -p output

# Clean previous build
echo "Cleaning previous build..."
make clean

# Build the program
echo "Building program..."
make

if [ $? -ne 0 ]; then
    echo "Build failed!"
    exit 1
fi

echo "Build successful!"
echo

# Run test with sample parameters
echo "Running test simulation..."
./lattice_qcd_hamiltonian input/parameters.inp

if [ $? -eq 0 ]; then
    echo
    echo "Test completed successfully!"
    echo "Check output directory for results:"
    ls -la output/
else
    echo "Test failed!"
    exit 1
fi