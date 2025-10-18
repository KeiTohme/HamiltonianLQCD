#!/usr/bin/env python3
"""
Plot results from Hamiltonian LQCD simulation

Usage:
    python plot_results.py <data_file> [output_image]
    
Example:
    python plot_results.py time_evolution.dat results.png
"""

import sys
import numpy as np
import matplotlib.pyplot as plt

def plot_results(data_file, output_file='results.png'):
    """
    Plot time evolution results from Hamiltonian LQCD simulation
    
    Parameters:
    -----------
    data_file : str
        Path to the data file
    output_file : str
        Path to save the plot
    """
    
    # Load data
    try:
        data = np.loadtxt(data_file)
    except Exception as e:
        print(f"Error loading data file: {e}")
        return
    
    # Extract columns
    step = data[:, 0]
    time = data[:, 1]
    energy = data[:, 2]
    plaquette = data[:, 3]
    
    # Create figure with subplots
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    fig.suptitle('Hamiltonian LQCD Simulation Results', fontsize=16, fontweight='bold')
    
    # Plot 1: Energy vs Time
    ax1 = axes[0, 0]
    ax1.plot(time, energy, 'b-', linewidth=2, label='Total Energy')
    ax1.set_xlabel('Time (lattice units)', fontsize=12)
    ax1.set_ylabel('Energy', fontsize=12)
    ax1.set_title('Energy Evolution', fontsize=14, fontweight='bold')
    ax1.grid(True, alpha=0.3)
    ax1.legend(fontsize=10)
    
    # Plot 2: Average Plaquette vs Time
    ax2 = axes[0, 1]
    ax2.plot(time, plaquette, 'r-', linewidth=2, label='Avg Plaquette')
    ax2.set_xlabel('Time (lattice units)', fontsize=12)
    ax2.set_ylabel('Average Plaquette', fontsize=12)
    ax2.set_title('Plaquette Evolution', fontsize=14, fontweight='bold')
    ax2.grid(True, alpha=0.3)
    ax2.legend(fontsize=10)
    ax2.set_ylim([min(plaquette) * 0.99, max(plaquette) * 1.01])
    
    # Plot 3: Energy Conservation (relative change)
    ax3 = axes[1, 0]
    energy_change = (energy - energy[0]) / abs(energy[0]) * 100  # Percentage change
    ax3.plot(time, energy_change, 'g-', linewidth=2, label='Relative Energy Change')
    ax3.set_xlabel('Time (lattice units)', fontsize=12)
    ax3.set_ylabel('ΔE / E₀ (%)', fontsize=12)
    ax3.set_title('Energy Conservation', fontsize=14, fontweight='bold')
    ax3.grid(True, alpha=0.3)
    ax3.axhline(y=0, color='k', linestyle='--', alpha=0.5)
    ax3.legend(fontsize=10)
    
    # Plot 4: Phase space (Energy vs Plaquette)
    ax4 = axes[1, 1]
    scatter = ax4.scatter(plaquette, energy, c=time, cmap='viridis', 
                         s=50, alpha=0.6, edgecolors='black', linewidth=0.5)
    ax4.set_xlabel('Average Plaquette', fontsize=12)
    ax4.set_ylabel('Energy', fontsize=12)
    ax4.set_title('Phase Space (colored by time)', fontsize=14, fontweight='bold')
    ax4.grid(True, alpha=0.3)
    cbar = plt.colorbar(scatter, ax=ax4)
    cbar.set_label('Time', fontsize=10)
    
    # Add statistics as text
    stats_text = f"Statistics:\n"
    stats_text += f"Time range: {time[0]:.4f} - {time[-1]:.4f}\n"
    stats_text += f"Energy range: {energy.min():.4f} - {energy.max():.4f}\n"
    stats_text += f"Energy std dev: {energy.std():.4e}\n"
    stats_text += f"Plaquette range: {plaquette.min():.6f} - {plaquette.max():.6f}\n"
    stats_text += f"Max energy change: {abs(energy_change).max():.4f}%"
    
    fig.text(0.02, 0.02, stats_text, fontsize=9, family='monospace',
             verticalalignment='bottom', bbox=dict(boxstyle='round', 
             facecolor='wheat', alpha=0.5))
    
    # Adjust layout and save
    plt.tight_layout(rect=[0, 0.08, 1, 0.96])
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"Plot saved to: {output_file}")
    
    # Print statistics
    print("\n" + "="*60)
    print("SIMULATION STATISTICS")
    print("="*60)
    print(stats_text)
    print("="*60)
    
    # Show plot
    try:
        plt.show()
    except:
        pass  # May fail in non-interactive environment

def main():
    """Main function"""
    
    if len(sys.argv) < 2:
        print("Usage: python plot_results.py <data_file> [output_image]")
        print("\nExample:")
        print("  python plot_results.py time_evolution.dat results.png")
        sys.exit(1)
    
    data_file = sys.argv[1]
    output_file = sys.argv[2] if len(sys.argv) > 2 else 'results.png'
    
    plot_results(data_file, output_file)

if __name__ == '__main__':
    main()
