#!/usr/bin/env python3
"""
Tilted-Ring Model for Tidal-Induced Stellar Warp
=================================================

The disk is decomposed into N concentric rings, each described by:
  - R_i: radius
  - M_i: mass
  - L_i: angular momentum vector (direction encodes tilt)

The rings interact via:
  1. Self-gravity coupling: neighboring rings exert gravitational torques
     that resist differential tilting (bending stiffness)
  2. Dark matter halo torque: an oblate/spherical halo provides a restoring
     torque (differential precession)
  3. External tidal torque: from an orbiting satellite galaxy

This follows the formalism of:
  - Sparke & Casertano (1988) for the bending wave dynamics
  - López-Corredoira et al. (2002) for the tidal torque calculation
  - Shen & Sellwood (2006) for the ring coupling

Physics:
  Each ring's angular momentum precesses according to:
    dL_i/dt = T_self_i + T_halo_i + T_tidal_i

  where the torques are computed from the gravitational interactions.

Units (Galactic):
  Length: kpc
  Mass: 1e10 Msun
  Time: Myr
  G = 4.498e-6 kpc^3 / (1e10 Msun * Myr^2)

Author: OpenResearch (automated)
Date: 2026-02-26
"""

import numpy as np
import json
import os
import time as walltime

# ============================================================
# Physical constants
# ============================================================
G_GRAV = 4.498502151575286e-6  # kpc^3 / (1e10 Msun * Myr^2)

# ============================================================
# Default configuration
# ============================================================
DEFAULT_CONFIG = {
    # --- Disk parameters ---
    "N_rings": 60,  # number of concentric rings
    "M_disk": 5.0,  # total disk mass [1e10 Msun]
    "R_disk": 3.5,  # disk scale length [kpc]
    "z_disk": 0.3,  # disk scale height [kpc]
    "R_min": 0.5,  # inner ring radius [kpc]
    "R_max": 18.0,  # outer ring radius [kpc]
    # --- Halo parameters (NFW) ---
    "M_halo_vir": 100.0,  # virial mass [1e10 Msun]
    "c_halo": 12.0,  # concentration
    "R_vir": 200.0,  # virial radius [kpc]
    "q_halo": 0.95,  # halo flattening (1=spherical, <1=oblate)
    # --- Satellite parameters ---
    "M_sat": 10.0,  # satellite mass [1e10 Msun]
    "R_sat_orbit": 25.0,  # semi-major axis [kpc]
    "sat_eccentricity": 0.3,  # eccentricity
    "incl_sat": 45.0,  # orbital inclination [deg]
    "Omega_sat_period": 800.0,  # orbital period [Myr]
    "sat_softening": 2.0,  # softening length [kpc]
    # --- Numerical parameters ---
    "dt": 1.0,  # time step [Myr]
    "T_total": 3000.0,  # total time [Myr]
    "snap_interval": 25.0,  # snapshot interval [Myr]
    "seed": 42,
    # --- Output ---
    "output_dir": "../output",
}


# ============================================================
# Satellite orbit
# ============================================================


def satellite_position(t, config):
    """Satellite on an inclined elliptical orbit (Kepler)."""
    period = config["Omega_sat_period"]
    omega_mean = 2.0 * np.pi / period
    a = config["R_sat_orbit"]
    ecc = config.get("sat_eccentricity", 0.3)
    incl = np.radians(config["incl_sat"])

    M = omega_mean * t
    E = M
    for _ in range(15):
        E = M + ecc * np.sin(E)

    cos_f = (np.cos(E) - ecc) / (1 - ecc * np.cos(E))
    sin_f = (np.sqrt(1 - ecc**2) * np.sin(E)) / (1 - ecc * np.cos(E))
    f = np.arctan2(sin_f, cos_f)

    r = a * (1 - ecc**2) / (1 + ecc * np.cos(f))

    x_orb = r * np.cos(f)
    y_orb = r * np.sin(f)

    x_sat = x_orb
    y_sat = y_orb * np.cos(incl)
    z_sat = y_orb * np.sin(incl)

    return np.array([x_sat, y_sat, z_sat])


# ============================================================
# Ring model setup
# ============================================================


def setup_rings(config):
    """
    Initialize the ring model.

    Returns:
      R : array of ring radii (N_rings,)
      M : array of ring masses (N_rings,)
      L_hat : array of angular momentum unit vectors (N_rings, 3)
             Initially all pointing in z-direction (flat disk)
      Omega : array of angular velocities (N_rings,) [rad/Myr]
      Sigma : surface density at each ring radius
    """
    N = config["N_rings"]
    R_min = config["R_min"]
    R_max = config["R_max"]
    Rd = config["R_disk"]
    M_disk = config["M_disk"]

    # Logarithmic spacing gives better resolution in inner disk
    R = np.linspace(R_min, R_max, N)
    dR = R[1] - R[0]

    # Surface density: exponential disk
    Sigma = (M_disk / (2.0 * np.pi * Rd**2)) * np.exp(-R / Rd)

    # Ring mass: M_i = 2*pi*R_i * Sigma_i * dR
    M = 2.0 * np.pi * R * Sigma * dR

    # Angular momentum direction: initially all in z-direction (flat disk)
    L_hat = np.zeros((N, 3))
    L_hat[:, 2] = 1.0  # all rings in z-direction

    # Circular velocity from NFW halo + disk
    M_vir = config["M_halo_vir"]
    c = config["c_halo"]
    R_vir = config["R_vir"]
    Rs = R_vir / c

    A = M_vir / (np.log(1.0 + c) - c / (1.0 + c))
    x = R / Rs
    M_enc_halo = A * (np.log(1.0 + x) - x / (1.0 + x))
    M_enc_disk = M_disk * (1.0 - (1.0 + R / Rd) * np.exp(-R / Rd))

    v_circ = np.sqrt(G_GRAV * (M_enc_halo + M_enc_disk) / R)

    # Angular velocity
    Omega = v_circ / R  # rad/Myr

    # Angular momentum magnitude per unit mass
    L_mag = v_circ * R  # specific angular momentum kpc^2/Myr

    return {
        "R": R,
        "dR": dR,
        "M": M,
        "L_hat": L_hat,
        "Omega": Omega,
        "Sigma": Sigma,
        "v_circ": v_circ,
        "L_mag": L_mag,
    }


# ============================================================
# Torque calculations
# ============================================================


def ring_coupling_torque(rings, config):
    """
    Compute the self-gravity torque between rings.

    When neighboring rings are mutually tilted, they exert a gravitational
    torque on each other that tends to align them (bending stiffness).

    The coupling torque on ring i from ring j is proportional to:
      T_ij ~ G * M_i * M_j / |R_i - R_j| * sin(angle_between_L_i_and_L_j)

    This implements the Sparke & Casertano (1988) coupling:
      T_i = sum_j C_ij * (L_hat_j - L_hat_i)

    where C_ij depends on the mutual gravitational interaction.

    For a thin disk, the coupling coefficient between rings i and j is:
      C_ij = -pi * G * Sigma_i * Sigma_j * R_i * R_j * dR^2 * K(R_i, R_j)

    where K involves elliptic integrals, but for |i-j| = 1 the dominant term is:
      C ~ -G * M_i * M_j / (2 * |R_i - R_j|)
    """
    N = len(rings["R"])
    R = rings["R"]
    M = rings["M"]
    L_hat = rings["L_hat"]
    Sigma = rings["Sigma"]
    dR = rings["dR"]
    h = config["z_disk"]

    torque = np.zeros((N, 3))

    # Compute coupling for all ring pairs (dominated by nearest neighbors)
    for i in range(N):
        for j in range(N):
            if i == j:
                continue

            # Coupling coefficient (softened)
            d_ij = np.abs(R[i] - R[j])
            r_eff = np.sqrt(d_ij**2 + h**2)  # softened by disk thickness

            # Gravitational coupling strength
            # Using the thin-disk approximation for ring-ring interaction
            # C_ij = G * M_i * M_j / r_eff^3 * R_i * R_j (torque per unit misalignment)
            C_ij = G_GRAV * M[i] * M[j] * min(R[i], R[j]) / (r_eff**2 * max(R[i], R[j]))

            # Torque: proportional to cross product of L vectors (drives alignment)
            # This produces a precession of ring i toward ring j's orientation
            delta_L = L_hat[j] - L_hat[i]

            # The torque is perpendicular to L_hat_i and in the direction of delta_L
            # T_i = C_ij * (L_hat_i x (L_hat_j x L_hat_i)) / |L_i|
            # Simplified: T_i = C_ij * (L_hat_j - (L_hat_j . L_hat_i) * L_hat_i) / M_i*R_i*Omega_i
            cross = np.cross(L_hat[i], np.cross(L_hat[j], L_hat[i]))

            torque[i] += C_ij * cross

    return torque


def ring_coupling_torque_fast(rings, config):
    """
    Optimized ring coupling torque using vectorized operations.
    Considers coupling between all ring pairs within n_couple neighbors.
    """
    N = len(rings["R"])
    R = rings["R"]
    M = rings["M"]
    L_hat = rings["L_hat"]
    h = config["z_disk"]

    torque = np.zeros((N, 3))

    n_couple = min(8, N // 2)

    for di in range(1, n_couple + 1):
        # Vectorized: all pairs (i, i+di) at once
        i_idx = np.arange(N - di)
        j_idx = i_idx + di

        d_ij = R[j_idx] - R[i_idx]
        r_eff = np.sqrt(d_ij**2 + h**2)

        C_ij = G_GRAV * M[i_idx] * M[j_idx] * R[i_idx] / (r_eff**2 * R[j_idx])

        # Torque on ring i from ring j: C * (L_i x (L_j x L_i))
        # L_j x L_i
        LjxLi = np.cross(L_hat[j_idx], L_hat[i_idx])
        cross_ij = np.cross(L_hat[i_idx], LjxLi)
        torque[i_idx] += C_ij[:, np.newaxis] * cross_ij

        # Torque on ring j from ring i
        C_ji = C_ij  # same magnitude
        LixLj = np.cross(L_hat[i_idx], L_hat[j_idx])
        cross_ji = np.cross(L_hat[j_idx], LixLj)
        torque[j_idx] += C_ji[:, np.newaxis] * cross_ji

    return torque


def halo_torque(rings, config):
    """
    Torque from a flattened dark matter halo.

    An oblate halo (q < 1) exerts a torque that causes differential
    precession of tilted rings. This is a key ingredient for warp dynamics:
    it provides a restoring torque that opposes the warp but also
    causes the line of nodes to precess.

    For a flattened NFW halo with axis ratio q:
      T_halo,i = -nu_p^2(R_i) * L_i * sin(theta_i) * cos(theta_i) * phi_hat

    where nu_p is the precession frequency and theta_i is the tilt angle.

    In the small-angle limit:
      dL_hat/dt |_halo = -nu_p^2(R) * (L_hat_z_component - 1) correction

    For a slightly oblate halo, the precession frequency is:
      nu_p^2 ~ (1 - q^2) * G*M_enc(R) / R^3
    """
    N = len(rings["R"])
    R = rings["R"]
    L_hat = rings["L_hat"]
    q = config["q_halo"]
    M_vir = config["M_halo_vir"]
    c = config["c_halo"]
    R_vir = config["R_vir"]
    Rs = R_vir / c

    A = M_vir / (np.log(1.0 + c) - c / (1.0 + c))
    x = R / Rs
    M_enc = A * (np.log(1.0 + x) - x / (1.0 + x))

    # Flattening parameter
    epsilon = 1 - q**2  # ~0.1 for q=0.95

    # Halo precession frequency squared
    nu_p_sq = epsilon * G_GRAV * M_enc / R**3

    # Vectorized halo torque: T = nu_p^2 * M * L_mag * (z_hat x L_hat)
    z_hat = np.array([0, 0, 1.0])

    # z_hat x L_hat for all rings
    cross_zL = np.cross(z_hat, L_hat)  # (N, 3)

    prefactor = nu_p_sq * rings["M"] * rings["L_mag"]  # (N,)
    torque = prefactor[:, np.newaxis] * cross_zL

    return torque


def tidal_torque(rings, t, config):
    """
    Tidal torque from the satellite galaxy on each ring.

    The tidal potential of a point mass at r_sat produces a torque on a
    ring of radius R:

      T_tidal = (3/2) * G * M_sat * M_ring * R^2 / d^5 *
                (d_vec x L_hat) * (d_vec . L_hat)

    where d is the distance from the galactic center to the satellite.

    This is the quadrupole (l=2) tidal torque which is the dominant term
    for warp generation.
    """
    r_sat = satellite_position(t, config)
    M_sat = config["M_sat"]
    eps = config["sat_softening"]

    N = len(rings["R"])
    R = rings["R"]
    M = rings["M"]
    L_hat = rings["L_hat"]

    d = np.sqrt(np.sum(r_sat**2) + eps**2)
    d_hat = r_sat / d

    # Vectorized tidal torque
    dot_dL = L_hat @ d_hat  # (N,)
    cross_dL = np.cross(d_hat, L_hat)  # (N, 3)

    T_mag = 1.5 * G_GRAV * M_sat * M * R**2 / d**3  # (N,)

    torque = (T_mag * dot_dL)[:, np.newaxis] * cross_dL

    return torque


# ============================================================
# Time integration
# ============================================================


def compute_dLhat_dt(rings, t, config, include_tidal=True):
    """
    Compute the time derivative of L_hat for each ring.

    dL_hat_i/dt = T_total_i / (M_i * R_i^2 * Omega_i)

    where the factor M_i * R_i^2 * Omega_i = |L_i| is the angular
    momentum magnitude.
    """
    N = len(rings["R"])

    # Compute all torques
    T_self = ring_coupling_torque_fast(rings, config)
    T_halo = halo_torque(rings, config)
    T_total = T_self + T_halo

    if include_tidal:
        T_tidal = tidal_torque(rings, t, config)
        T_total += T_tidal

    # Convert torque to angular velocity of L_hat
    # dL_hat/dt = T / |L| where |L| = M * v_circ * R = M * L_mag
    L_magnitude = rings["M"] * rings["L_mag"]  # (N,)
    L_magnitude = np.maximum(L_magnitude, 1e-30)
    dL_hat = T_total / L_magnitude[:, np.newaxis]

    return dL_hat, T_self, T_halo, T_total


def normalize_L_hat(L_hat):
    """Ensure all L_hat vectors are unit vectors."""
    norms = np.linalg.norm(L_hat, axis=1, keepdims=True)
    norms = np.maximum(norms, 1e-10)
    return L_hat / norms


def run_ring_simulation(config, include_tidal=True, label="tidal"):
    """
    Run the tilted-ring warp simulation.
    """
    print(f"=" * 60)
    print(f"Ring Model Simulation: {label}")
    print(
        f"N_rings={config['N_rings']}, T={config['T_total']} Myr, dt={config['dt']} Myr"
    )
    print(f"Tidal perturbation: {'ON' if include_tidal else 'OFF'}")
    print(f"M_sat={config['M_sat']} x 1e10 Msun, R_orbit={config['R_sat_orbit']} kpc")
    print(f"Inclination={config['incl_sat']} deg, q_halo={config['q_halo']}")
    print(f"=" * 60)

    # Initialize rings
    rings = setup_rings(config)

    dt = config["dt"]
    T = config["T_total"]
    snap_dt = config["snap_interval"]
    N_steps = int(T / dt)

    out_dir = os.path.join(config["output_dir"], label)
    os.makedirs(out_dir, exist_ok=True)

    # Save config
    with open(os.path.join(out_dir, "config.json"), "w") as f:
        json.dump(config, f, indent=2)

    # Storage
    snap_times = []
    snap_count = 0

    # Save initial state
    np.savez_compressed(
        os.path.join(out_dir, f"ring_{snap_count:04d}.npz"),
        R=rings["R"],
        M=rings["M"],
        L_hat=rings["L_hat"].copy(),
        Omega=rings["Omega"],
        Sigma=rings["Sigma"],
        v_circ=rings["v_circ"],
        L_mag=rings["L_mag"],
        t=0.0,
    )
    snap_times.append(0.0)
    snap_count += 1

    t_wall_start = walltime.time()

    for step in range(1, N_steps + 1):
        t = step * dt

        # RK4 integration for accuracy
        L_hat_0 = rings["L_hat"].copy()

        # k1
        dL1, _, _, _ = compute_dLhat_dt(rings, t - dt, config, include_tidal)

        # k2
        rings["L_hat"] = normalize_L_hat(L_hat_0 + 0.5 * dt * dL1)
        dL2, _, _, _ = compute_dLhat_dt(rings, t - 0.5 * dt, config, include_tidal)

        # k3
        rings["L_hat"] = normalize_L_hat(L_hat_0 + 0.5 * dt * dL2)
        dL3, _, _, _ = compute_dLhat_dt(rings, t - 0.5 * dt, config, include_tidal)

        # k4
        rings["L_hat"] = normalize_L_hat(L_hat_0 + dt * dL3)
        dL4, T_self, T_halo, T_total = compute_dLhat_dt(rings, t, config, include_tidal)

        # Update
        rings["L_hat"] = normalize_L_hat(
            L_hat_0 + (dt / 6.0) * (dL1 + 2 * dL2 + 2 * dL3 + dL4)
        )

        # Save snapshot
        if step % int(snap_dt / dt) == 0:
            # Compute warp diagnostics
            tilt_angle = np.arccos(np.clip(rings["L_hat"][:, 2], -1, 1))
            max_tilt = np.degrees(np.max(tilt_angle))

            # Outer disk tilt (R > 10 kpc)
            outer = rings["R"] > 10
            outer_tilt = np.degrees(np.mean(tilt_angle[outer])) if np.any(outer) else 0

            np.savez_compressed(
                os.path.join(out_dir, f"ring_{snap_count:04d}.npz"),
                R=rings["R"],
                M=rings["M"],
                L_hat=rings["L_hat"].copy(),
                Omega=rings["Omega"],
                Sigma=rings["Sigma"],
                v_circ=rings["v_circ"],
                L_mag=rings["L_mag"],
                t=t,
            )
            snap_times.append(t)
            snap_count += 1

            elapsed = walltime.time() - t_wall_start
            eta = elapsed / step * (N_steps - step)

            # Satellite position
            r_sat = satellite_position(t, config)
            d_sat = np.linalg.norm(r_sat)

            print(
                f"  t={t:6.0f} Myr  |  max_tilt={max_tilt:6.3f} deg  |  "
                f"outer_tilt={outer_tilt:6.3f} deg  |  d_sat={d_sat:5.1f} kpc  |  "
                f"wall={elapsed:.1f}s  ETA={eta:.1f}s"
            )

    np.save(os.path.join(out_dir, "snap_times.npy"), np.array(snap_times))

    total_wall = walltime.time() - t_wall_start
    print(f"\nSimulation '{label}' complete. Wall time: {total_wall:.1f}s")
    print(f"Saved {snap_count} snapshots to {out_dir}/")

    return out_dir


# ============================================================
# Generate N-body particles from ring model (for visualization)
# ============================================================


def rings_to_particles(ring_data, N_particles=50000, seed=42):
    """
    Generate particle positions from the ring model state.
    Each ring contributes particles proportional to its mass,
    distributed as a thin annulus tilted according to L_hat.
    """
    rng = np.random.default_rng(seed)
    R = ring_data["R"]
    M = ring_data["M"]
    L_hat = ring_data["L_hat"]
    dR = R[1] - R[0]

    N_rings = len(R)
    total_mass = np.sum(M)

    all_pos = []

    for i in range(N_rings):
        # Number of particles in this ring
        n_i = max(1, int(N_particles * M[i] / total_mass))

        # Generate particles in the ring's plane
        phi = rng.uniform(0, 2 * np.pi, n_i)
        r_i = R[i] + rng.normal(0, dR * 0.3, n_i)  # small radial spread
        z_i = rng.normal(0, 0.3, n_i)  # vertical spread (scale height)

        # Positions in ring's local frame (before tilt)
        x_local = r_i * np.cos(phi)
        y_local = r_i * np.sin(phi)
        z_local = z_i

        # Rotate to align with L_hat[i]
        # The ring's normal is L_hat[i]; we need to rotate from z-axis to L_hat[i]
        L = L_hat[i]
        z_axis = np.array([0, 0, 1.0])

        if np.allclose(L, z_axis, atol=1e-8):
            # No rotation needed
            pos_ring = np.column_stack([x_local, y_local, z_local])
        elif np.allclose(L, -z_axis, atol=1e-8):
            # 180 degree rotation
            pos_ring = np.column_stack([x_local, y_local, -z_local])
        else:
            # Rotation axis: z x L
            rot_axis = np.cross(z_axis, L)
            rot_axis /= np.linalg.norm(rot_axis)

            # Rotation angle
            cos_angle = np.dot(z_axis, L)
            angle = np.arccos(np.clip(cos_angle, -1, 1))

            # Rodrigues' rotation formula
            K = np.array(
                [
                    [0, -rot_axis[2], rot_axis[1]],
                    [rot_axis[2], 0, -rot_axis[0]],
                    [-rot_axis[1], rot_axis[0], 0],
                ]
            )

            R_mat = np.eye(3) + np.sin(angle) * K + (1 - cos_angle) * K @ K

            local_pos = np.column_stack([x_local, y_local, z_local])
            pos_ring = (R_mat @ local_pos.T).T

        all_pos.append(pos_ring)

    return np.vstack(all_pos)


# ============================================================
# Entry point
# ============================================================

if __name__ == "__main__":
    import sys

    config = DEFAULT_CONFIG.copy()

    if "--test" in sys.argv:
        config["T_total"] = 500.0
        config["snap_interval"] = 50.0
        print("*** TEST MODE ***")

    # Run baseline
    if "--baseline" in sys.argv or "--all" in sys.argv:
        run_ring_simulation(config, include_tidal=False, label="baseline")

    # Run tidal
    if "--tidal" in sys.argv or "--all" in sys.argv:
        run_ring_simulation(config, include_tidal=True, label="tidal")

    if len(sys.argv) == 1:
        run_ring_simulation(config, include_tidal=False, label="baseline")
        run_ring_simulation(config, include_tidal=True, label="tidal")
