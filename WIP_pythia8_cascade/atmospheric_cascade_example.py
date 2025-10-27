#!/usr/bin/env python3
"""
Atmospheric cascade simulation using PythiaCascadePython.

This script demonstrates atmospheric shower simulation similar to 
the C++ main483.cc example, but using the pure Python implementation.
"""

import numpy as np
import matplotlib.pyplot as plt
from WIP_pythia8_cascade.pythia_cascade_python import PythiaCascadePython
from chromo.kinematics import FixedTarget
from chromo.constants import GeV


class AtmosphereConfig:
    """Configuration for atmospheric model."""
    
    def __init__(self, legend, do_exponential, do_heavy_ion, zenith_angle, e_kin_min):
        self.legend = legend
        self.do_exponential = do_exponential
        self.do_heavy_ion = do_heavy_ion
        self.zenith_angle = zenith_angle
        self.e_kin_min = e_kin_min
    
    def initial_height(self):
        """Get initial height above ground."""
        return MEDIUM_HEIGHT if self.do_exponential else H
    
    def get_depth(self, h):
        """Get medium depth at specified height."""
        return (1.0 / np.cos(self.zenith_angle) * 
                (RHO_AIR * H * np.exp(-h / H) if self.do_exponential 
                 else RHO_AIR * (self.initial_height() - h)))


# Constants (converted from main483.cc)
KM2MM = 1e6  # mm per km
kg_m3_to_GeV_c2mm1mb1 = 5.60958865e26 * 0.001 * 1e-31

# Medium parameters
M_AIR = 0.9315
MEDIUM_DENSITY = 1.225 * kg_m3_to_GeV_c2mm1mb1 / M_AIR
RHO_AIR = 1.225e-4  # g cm^-2 mm^-1
MEDIUM_HEIGHT = 100 * KM2MM
H = 10.4 * KM2MM
MP = 0.938272  # Proton mass in GeV


def simulate_atmospheric_cascade():
    """
    Simulate atmospheric cascades with different configurations.
    
    This is a simplified version of the main483.cc atmospheric shower simulation.
    """
    print("Atmospheric Cascade Simulation")
    print("=" * 50)
    
    # Configuration cases
    configurations = [
        AtmosphereConfig("Uniform p/n atmosphere", False, False, 0.0, 0.2),
        AtmosphereConfig("Uniform nitrogen", False, True, 0.0, 0.2),
        AtmosphereConfig("Exponential nitrogen", True, True, 0.0, 0.2),
        AtmosphereConfig("Exponential nitrogen at 45°", True, True, np.pi/4, 0.2),
    ]
    
    # Primary proton energy
    p_primary = 1e6 * GeV  # 1 EeV
    n_events = 5  # Reduced for demonstration
    ecm_min = 10.5  # Minimum CM energy (GeV)
    
    print(f"Primary energy: {p_primary/GeV:.0e} GeV")
    print(f"Number of events per case: {n_events}")
    
    results = {}
    
    for i_case, config in enumerate(configurations):
        print(f"\nCase {i_case + 1}: {config.legend}")
        print("-" * 40)
        
        # Set up cascade generator
        # Use nitrogen-14 as target for heavy ion cases, proton for p/n gas
        if config.do_heavy_ion:
            target = (7, 14)  # N-14
        else:
            target = (1, 1)   # proton
            
        try:
            evt_kin = FixedTarget(p_primary, "p", target)
            cascade = PythiaCascadePython(
                evt_kin,
                seed=42 + i_case,
                max_energy=2e6,  # 2 EeV max
                rapid_decays=True,
            )
            
            # Get cross section
            cs = cascade._cross_section()
            print(f"Cross section: {cs.total:.2f} mb")
            
            case_results = {
                'interactions': [],
                'hadrons': [],
                'muons': [],
                'neutrinos': [],
            }
            
            # Generate events
            for i_event in range(n_events):
                print(f"  Event {i_event + 1}/{n_events}", end="")
                
                try:
                    # Generate collision
                    success = cascade._generate()
                    if not success:
                        print(" - failed")
                        continue
                    
                    event_data = cascade._current_event_data
                    n_particles = len(event_data['pid'])
                    print(f" - {n_particles} particles")
                    
                    # Analyze particles
                    analyze_event(event_data, config, case_results)
                    
                except Exception as e:
                    print(f" - error: {e}")
                    continue
            
            results[config.legend] = case_results
            
        except Exception as e:
            print(f"Failed to initialize case: {e}")
            continue
    
    # Simple analysis and plotting
    plot_results(results, configurations)
    
    return results


def analyze_event(event_data, config, case_results):
    """Analyze a single event for atmospheric shower properties."""
    
    # Count different particle types
    hadron_count = 0
    muon_count = 0
    neutrino_count = 0
    interaction_count = 0
    
    for i, (pid, status, energy, mass) in enumerate(zip(
        event_data['pid'], event_data['status'], 
        event_data['en'], event_data['m']
    )):
        # Skip if below energy threshold
        if energy - mass < config.e_kin_min:
            continue
            
        # Count interactions (target particles)
        if status in [-181, -182]:
            interaction_count += 1
            
        # Count final state particles
        elif status > 0:
            abs_pid = abs(pid)
            
            # Hadrons
            if is_hadron(abs_pid):
                hadron_count += 1
                
            # Muons
            elif abs_pid == 13:
                muon_count += 1
                
            # Neutrinos
            elif abs_pid in [12, 14, 16]:
                neutrino_count += 1
    
    case_results['interactions'].append(interaction_count)
    case_results['hadrons'].append(hadron_count)
    case_results['muons'].append(muon_count)
    case_results['neutrinos'].append(neutrino_count)


def is_hadron(abs_pid):
    """Check if particle is a hadron."""
    # Simplified hadron identification
    return (abs_pid > 100 and abs_pid not in [11, 12, 13, 14, 15, 16, 22, 130, 310])


def plot_results(results, configurations):
    """Create simple plots of the results."""
    try:
        fig, axes = plt.subplots(2, 2, figsize=(12, 10))
        fig.suptitle("Atmospheric Cascade Results", fontsize=16)
        
        colors = ['red', 'blue', 'black', 'green']
        
        # Plot histograms for each particle type
        particle_types = [
            ('interactions', 'Number of Interactions'),
            ('hadrons', 'Number of Hadrons'),
            ('muons', 'Number of Muons'),  
            ('neutrinos', 'Number of Neutrinos'),
        ]
        
        for idx, (particle_type, title) in enumerate(particle_types):
            ax = axes[idx // 2, idx % 2]
            
            for i, (config_name, case_data) in enumerate(results.items()):
                if particle_type in case_data and case_data[particle_type]:
                    counts = case_data[particle_type]
                    ax.hist(counts, bins=max(1, len(set(counts))), 
                           alpha=0.7, label=config_name, color=colors[i % len(colors)])
            
            ax.set_title(title)
            ax.set_xlabel('Count')
            ax.set_ylabel('Frequency')
            ax.legend()
            ax.grid(True, alpha=0.3)
        
        plt.tight_layout()
        plt.savefig('atmospheric_cascade_results.png', dpi=150, bbox_inches='tight')
        print(f"\nResults plot saved as 'atmospheric_cascade_results.png'")
        
    except ImportError:
        print("\nMatplotlib not available, skipping plots")
    except Exception as e:
        print(f"\nPlotting failed: {e}")


def demonstrate_cross_sections():
    """Demonstrate cross section calculations."""
    print("\nCross Section Demonstration")
    print("=" * 30)
    
    energies = [10, 100, 1000, 10000]  # GeV
    targets = [
        ((1, 1), "proton"),
        ((6, 12), "carbon-12"),
        ((7, 14), "nitrogen-14"),
        ((8, 16), "oxygen-16"),
    ]
    
    print(f"{'Energy (GeV)':<12} {'Target':<12} {'σ_total (mb)':<15}")
    print("-" * 45)
    
    for energy in energies:
        for target, name in targets:
            try:
                evt_kin = FixedTarget(energy * GeV, "p", target)
                cascade = PythiaCascadePython(evt_kin, seed=42)
                cs = cascade._cross_section()
                
                print(f"{energy:<12.0f} {name:<12} {cs.total:<15.2f}")
                
            except Exception as e:
                print(f"{energy:<12.0f} {name:<12} {'Error':<15}")


def main():
    """Main function."""
    print("PythiaCascadePython Atmospheric Simulation")
    print("=" * 60)
    
    # Demonstrate cross sections first
    demonstrate_cross_sections()
    
    # Run atmospheric cascade simulation
    results = simulate_atmospheric_cascade()
    
    # Print summary
    print("\nSummary:")
    print("=" * 20)
    for config_name, data in results.items():
        if data['hadrons']:
            avg_hadrons = np.mean(data['hadrons'])
            avg_muons = np.mean(data['muons'])
            avg_interactions = np.mean(data['interactions'])
            
            print(f"{config_name}:")
            print(f"  Avg interactions: {avg_interactions:.1f}")
            print(f"  Avg hadrons: {avg_hadrons:.1f}")
            print(f"  Avg muons: {avg_muons:.1f}")


if __name__ == "__main__":
    main()
