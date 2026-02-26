#!/usr/bin/env python3
"""
Master runner script for the Tidal-Induced Stellar Warp study.

Runs:
  1. Baseline simulation (isolated disk, no tidal field)
  2. Fiducial tidal simulation
  3. Parameter study (varying satellite mass and orbital inclination)
  4. Analysis and visualization
"""

import os
import sys
import json
import copy
import time
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from simulation import DEFAULT_CONFIG, run_simulation
from analysis import warp_evolution_summary, warp_profile, fourier_bending_modes
from visualize import generate_all_figures

OUTPUT_BASE = os.path.join(os.path.dirname(__file__), "..", "output")
FIG_DIR = os.path.join(os.path.dirname(__file__), "..", "figures")


def run_fiducial(config):
    """Run baseline and fiducial tidal simulations."""
    config = config.copy()
    config["output_dir"] = OUTPUT_BASE

    print("\n" + "=" * 70)
    print("PHASE 1: BASELINE (Isolated Disk)")
    print("=" * 70)
    run_simulation(config, include_tidal=False, label="baseline")

    print("\n" + "=" * 70)
    print("PHASE 2: FIDUCIAL TIDAL SIMULATION")
    print("=" * 70)
    run_simulation(config, include_tidal=True, label="tidal")


def run_parameter_study(base_config):
    """Run parameter study varying satellite mass and inclination."""

    # --- Vary satellite mass ---
    M_sat_values = [0.5, 1.0, 2.0, 4.0, 8.0]
    print("\n" + "=" * 70)
    print("PHASE 3a: PARAMETER STUDY - Satellite Mass")
    print("=" * 70)

    for M_sat in M_sat_values:
        config = base_config.copy()
        config["output_dir"] = OUTPUT_BASE
        config["M_sat"] = M_sat
        label = f"param_Msat_{M_sat:.1f}"
        print(f"\n--- M_sat = {M_sat} x 1e10 Msun ---")
        run_simulation(config, include_tidal=True, label=label)

    # --- Vary orbital inclination ---
    incl_values = [15.0, 30.0, 45.0, 60.0, 75.0, 90.0]
    print("\n" + "=" * 70)
    print("PHASE 3b: PARAMETER STUDY - Orbital Inclination")
    print("=" * 70)

    for incl in incl_values:
        config = base_config.copy()
        config["output_dir"] = OUTPUT_BASE
        config["incl_sat"] = incl
        label = f"param_incl_{incl:.0f}"
        print(f"\n--- inclination = {incl} deg ---")
        run_simulation(config, include_tidal=True, label=label)


def analyze_results():
    """Run analysis and produce RESULTS.json."""
    print("\n" + "=" * 70)
    print("PHASE 4: ANALYSIS")
    print("=" * 70)

    results = {"config": DEFAULT_CONFIG.copy(), "metrics": {}}

    # Analyze baseline
    baseline_dir = os.path.join(OUTPUT_BASE, "baseline")
    if os.path.exists(baseline_dir):
        bl_summary = warp_evolution_summary(baseline_dir)
        results["metrics"]["baseline"] = {
            "max_warp_amplitude": round(
                float(np.max(bl_summary["max_warp_amplitude"])), 4
            ),
            "final_z_rms": round(float(bl_summary["z_rms"][-1]), 4),
        }
        print(
            f"  Baseline: max warp = {results['metrics']['baseline']['max_warp_amplitude']:.4f} kpc"
        )

    # Analyze tidal
    tidal_dir = os.path.join(OUTPUT_BASE, "tidal")
    if os.path.exists(tidal_dir):
        td_summary = warp_evolution_summary(tidal_dir)
        results["metrics"]["tidal"] = {
            "max_warp_amplitude": round(
                float(np.max(td_summary["max_warp_amplitude"])), 4
            ),
            "final_warp_amplitude": round(
                float(td_summary["max_warp_amplitude"][-1]), 4
            ),
            "final_z_rms": round(float(td_summary["z_rms"][-1]), 4),
            "mean_z_outer_peak": round(float(np.max(td_summary["mean_z_outer"])), 4),
        }
        print(
            f"  Tidal: max warp = {results['metrics']['tidal']['max_warp_amplitude']:.4f} kpc"
        )

        # Warp amplification factor
        if "baseline" in results["metrics"]:
            bl_max = results["metrics"]["baseline"]["max_warp_amplitude"]
            td_max = results["metrics"]["tidal"]["max_warp_amplitude"]
            if bl_max > 0:
                results["metrics"]["warp_amplification"] = round(td_max / bl_max, 4)
            else:
                results["metrics"]["warp_amplification"] = float("inf")

    # Parameter study results
    M_sat_values = [0.5, 1.0, 2.0, 4.0, 8.0]
    param_Msat_results = {}
    for M_sat in M_sat_values:
        d = os.path.join(OUTPUT_BASE, f"param_Msat_{M_sat:.1f}")
        if os.path.exists(d):
            s = warp_evolution_summary(d)
            param_Msat_results[str(M_sat)] = {
                "max_warp_amplitude": round(float(np.max(s["max_warp_amplitude"])), 4),
            }
    results["metrics"]["param_study_Msat"] = param_Msat_results

    incl_values = [15.0, 30.0, 45.0, 60.0, 75.0, 90.0]
    param_incl_results = {}
    for incl in incl_values:
        d = os.path.join(OUTPUT_BASE, f"param_incl_{incl:.0f}")
        if os.path.exists(d):
            s = warp_evolution_summary(d)
            param_incl_results[str(incl)] = {
                "max_warp_amplitude": round(float(np.max(s["max_warp_amplitude"])), 4),
            }
    results["metrics"]["param_study_inclination"] = param_incl_results

    # Save RESULTS.json
    results_path = os.path.join(OUTPUT_BASE, "..", "RESULTS.json")
    with open(results_path, "w") as f:
        json.dump(results, f, indent=2, default=str)
    print(f"\n  Results saved to {results_path}")

    return results


def main():
    t_start = time.time()

    config = DEFAULT_CONFIG.copy()

    # Use moderate resolution for full run
    config["N_particles"] = 20000
    config["dt"] = 2.0
    config["T_total"] = 2000.0
    config["snap_interval"] = 50.0

    # Phase 1 & 2: Baseline + Fiducial
    run_fiducial(config)

    # Phase 3: Parameter study (also at moderate resolution for speed)
    run_parameter_study(config)

    # Phase 4: Analysis
    results = analyze_results()

    # Phase 5: Visualization
    print("\n" + "=" * 70)
    print("PHASE 5: VISUALIZATION")
    print("=" * 70)
    generate_all_figures(OUTPUT_BASE, FIG_DIR, config)

    # Also generate parameter study plots
    from visualize import plot_parameter_study

    M_sat_values = [0.5, 1.0, 2.0, 4.0, 8.0]
    M_sat_dirs = [
        os.path.join(OUTPUT_BASE, f"param_Msat_{v:.1f}") for v in M_sat_values
    ]
    plot_parameter_study(M_sat_dirs, M_sat_values, "M_sat [1e10 Msun]", FIG_DIR)

    incl_values = [15.0, 30.0, 45.0, 60.0, 75.0, 90.0]
    incl_dirs = [os.path.join(OUTPUT_BASE, f"param_incl_{v:.0f}") for v in incl_values]
    plot_parameter_study(incl_dirs, incl_values, "Inclination [deg]", FIG_DIR)

    total_time = time.time() - t_start
    print(f"\n{'=' * 70}")
    print(f"ALL DONE. Total wall time: {total_time:.1f} s ({total_time / 60:.1f} min)")
    print(f"{'=' * 70}")


if __name__ == "__main__":
    main()
