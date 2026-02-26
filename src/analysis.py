#!/usr/bin/env python3
"""
Analysis module for Tidal-Induced Stellar Warp Simulation
=========================================================
Computes warp profiles, bending modes, Fourier decomposition,
and diagnostic quantities from simulation snapshots.
"""

import numpy as np
import os
import json
import glob


def load_snapshot(snap_path):
    """Load a single snapshot."""
    data = np.load(snap_path)
    return {
        "pos": data["pos"],
        "vel": data["vel"],
        "t": float(data["t"]),
        "mass_pp": float(data["mass_pp"]),
    }


def load_all_snapshots(sim_dir):
    """Load all snapshots from a simulation directory."""
    files = sorted(glob.glob(os.path.join(sim_dir, "snap_*.npz")))
    snapshots = []
    for f in files:
        snapshots.append(load_snapshot(f))
    return snapshots


def warp_profile(pos, R_bins=None, n_phi_bins=12):
    """
    Compute the warp profile: mean z as a function of R and azimuthal angle.

    Returns:
        R_mid : radial bin centers
        z_mean : mean z height in each radial bin
        z_std : std of z in each radial bin
        warp_amplitude : amplitude of m=1 bending mode at each R
        warp_phase : phase of m=1 bending mode at each R
    """
    x, y, z = pos[:, 0], pos[:, 1], pos[:, 2]
    R = np.sqrt(x**2 + y**2)
    phi = np.arctan2(y, x)

    if R_bins is None:
        R_bins = np.linspace(0.5, 15.0, 30)

    R_mid = 0.5 * (R_bins[:-1] + R_bins[1:])
    n_R = len(R_mid)

    z_mean = np.zeros(n_R)
    z_std = np.zeros(n_R)
    warp_amplitude = np.zeros(n_R)
    warp_phase = np.zeros(n_R)
    n_in_bin = np.zeros(n_R, dtype=int)

    for i in range(n_R):
        mask = (R >= R_bins[i]) & (R < R_bins[i + 1])
        n_in_bin[i] = np.sum(mask)

        if n_in_bin[i] < 10:
            continue

        z_ring = z[mask]
        phi_ring = phi[mask]

        z_mean[i] = np.mean(z_ring)
        z_std[i] = np.std(z_ring)

        # Fourier decomposition for m=1 (warp) mode
        # z(phi) = A * cos(phi - phi_0)
        # = a1*cos(phi) + b1*sin(phi)
        a1 = 2.0 * np.mean(z_ring * np.cos(phi_ring))
        b1 = 2.0 * np.mean(z_ring * np.sin(phi_ring))

        warp_amplitude[i] = np.sqrt(a1**2 + b1**2)
        warp_phase[i] = np.arctan2(b1, a1)

    return {
        "R_mid": R_mid,
        "z_mean": z_mean,
        "z_std": z_std,
        "warp_amplitude": warp_amplitude,
        "warp_phase": warp_phase,
        "n_in_bin": n_in_bin,
    }


def fourier_bending_modes(pos, R_bins=None, m_max=4):
    """
    Compute Fourier bending modes m=0,1,2,...m_max at each radius.

    The m=1 mode is the classic warp; m=0 is vertical displacement;
    m=2 is a saddle-like mode (corrugation).
    """
    x, y, z = pos[:, 0], pos[:, 1], pos[:, 2]
    R = np.sqrt(x**2 + y**2)
    phi = np.arctan2(y, x)

    if R_bins is None:
        R_bins = np.linspace(0.5, 15.0, 30)

    R_mid = 0.5 * (R_bins[:-1] + R_bins[1:])
    n_R = len(R_mid)

    amplitudes = np.zeros((m_max + 1, n_R))
    phases = np.zeros((m_max + 1, n_R))

    for i in range(n_R):
        mask = (R >= R_bins[i]) & (R < R_bins[i + 1])
        if np.sum(mask) < 10:
            continue

        z_ring = z[mask]
        phi_ring = phi[mask]

        for m in range(m_max + 1):
            if m == 0:
                amplitudes[m, i] = np.abs(np.mean(z_ring))
            else:
                am = 2.0 * np.mean(z_ring * np.cos(m * phi_ring))
                bm = 2.0 * np.mean(z_ring * np.sin(m * phi_ring))
                amplitudes[m, i] = np.sqrt(am**2 + bm**2)
                phases[m, i] = np.arctan2(bm, am)

    return {
        "R_mid": R_mid,
        "amplitudes": amplitudes,
        "phases": phases,
    }


def compute_angular_momentum_profile(pos, vel, mass_pp, R_bins=None):
    """Compute the tilt of angular momentum vector in radial bins."""
    x, y, z = pos[:, 0], pos[:, 1], pos[:, 2]
    R = np.sqrt(x**2 + y**2)

    if R_bins is None:
        R_bins = np.linspace(0.5, 15.0, 30)

    R_mid = 0.5 * (R_bins[:-1] + R_bins[1:])
    n_R = len(R_mid)

    L_tilt_x = np.zeros(n_R)  # tilt angles
    L_tilt_y = np.zeros(n_R)

    for i in range(n_R):
        mask = (R >= R_bins[i]) & (R < R_bins[i + 1])
        if np.sum(mask) < 10:
            continue

        p = pos[mask]
        v = vel[mask]

        # Angular momentum: L = m * (r x v)
        L = mass_pp * np.cross(p, v)
        L_total = np.sum(L, axis=0)
        L_mag = np.linalg.norm(L_total)

        if L_mag > 0:
            # Tilt angle from z-axis
            L_tilt_x[i] = np.degrees(np.arctan2(L_total[0], L_total[2]))
            L_tilt_y[i] = np.degrees(np.arctan2(L_total[1], L_total[2]))

    return {
        "R_mid": R_mid,
        "L_tilt_x": L_tilt_x,
        "L_tilt_y": L_tilt_y,
    }


def compute_energy(pos, vel, mass_pp, config):
    """Compute kinetic, potential (halo), and total energy."""
    KE = 0.5 * mass_pp * np.sum(vel**2)

    # Halo potential energy
    r = np.sqrt(np.sum(pos**2, axis=1))
    from simulation import nfw_potential

    phi_halo = nfw_potential(r, config["M_halo_vir"], config["c_halo"], config["R_vir"])
    PE_halo = mass_pp * np.sum(phi_halo)

    return {
        "KE": float(KE),
        "PE_halo": float(PE_halo),
        "E_total": float(KE + PE_halo),
    }


def warp_evolution_summary(sim_dir, R_bins=None):
    """
    Compute warp properties across all snapshots.
    Returns time series of warp amplitude, phase, and vertical thickness.
    """
    snapshots = load_all_snapshots(sim_dir)

    times = []
    max_warp_amplitudes = []
    mean_z_outer = []
    z_rms_all = []

    warp_profiles_list = []

    for snap in snapshots:
        t = snap["t"]
        pos = snap["pos"]

        wp = warp_profile(pos, R_bins=R_bins)

        times.append(t)
        max_warp_amplitudes.append(np.max(wp["warp_amplitude"]))

        # Mean z of outer disk (R > 8 kpc)
        outer_mask = wp["R_mid"] > 8.0
        if np.any(outer_mask):
            mean_z_outer.append(np.mean(np.abs(wp["z_mean"][outer_mask])))
        else:
            mean_z_outer.append(0.0)

        z_rms_all.append(np.std(pos[:, 2]))
        warp_profiles_list.append(wp)

    return {
        "times": np.array(times),
        "max_warp_amplitude": np.array(max_warp_amplitudes),
        "mean_z_outer": np.array(mean_z_outer),
        "z_rms": np.array(z_rms_all),
        "profiles": warp_profiles_list,
    }


if __name__ == "__main__":
    import sys

    if len(sys.argv) < 2:
        print("Usage: python analysis.py <sim_dir>")
        sys.exit(1)

    sim_dir = sys.argv[1]
    print(f"Analyzing simulation in {sim_dir}")

    summary = warp_evolution_summary(sim_dir)

    print(f"\nTime range: {summary['times'][0]:.0f} - {summary['times'][-1]:.0f} Myr")
    print(f"Max warp amplitude: {np.max(summary['max_warp_amplitude']):.4f} kpc")
    print(f"Final z_rms: {summary['z_rms'][-1]:.4f} kpc")
    print(f"Max outer disk <|z|>: {np.max(summary['mean_z_outer']):.4f} kpc")
