#!/usr/bin/env python3
"""
Tidal-Induced Stellar Warp Simulation
======================================
N-body simulation of a stellar disk embedded in a dark matter halo,
subject to tidal perturbation from an orbiting satellite galaxy.

Physics:
  - Stellar disk: self-gravitating N-body particles
  - Dark matter halo: rigid NFW potential
  - Satellite: point-mass on a prescribed orbit, exerting tidal forces
  - Integration: Leapfrog (kick-drift-kick) with softened gravity

Units (Galactic):
  - Length: kpc
  - Mass: 1e10 Msun
  - Time: Myr
  - G = 4.498e-6 kpc^3 / (1e10 Msun * Myr^2)   [from astropy]

Author: OpenResearch (automated)
Date: 2026-02-26
"""

import numpy as np
from scipy.ndimage import map_coordinates
import json
import os
import time as walltime

# ============================================================
# Physical constants in our unit system
# ============================================================
G_GRAV = 4.498502151575286e-6  # kpc^3 / (1e10 Msun * Myr^2)

# ============================================================
# Configuration dataclass (as dict for JSON serialization)
# ============================================================
DEFAULT_CONFIG = {
    # --- Disk parameters ---
    "N_particles": 50000,  # number of stellar particles
    "M_disk": 5.0,  # disk mass [1e10 Msun]
    "R_disk": 3.5,  # disk scale length [kpc]
    "z_disk": 0.3,  # disk scale height [kpc]
    "R_trunc": 15.0,  # disk truncation radius [kpc]
    # --- Halo parameters (NFW) ---
    "M_halo_vir": 100.0,  # virial mass [1e10 Msun]
    "c_halo": 12.0,  # concentration parameter
    "R_vir": 200.0,  # virial radius [kpc]
    # --- Satellite parameters ---
    "M_sat": 10.0,  # satellite mass [1e10 Msun] (massive perturber)
    "R_sat_orbit": 25.0,  # semi-major axis [kpc] (pericenter ~17 kpc with e=0.3)
    "sat_eccentricity": 0.3,  # orbital eccentricity
    "incl_sat": 45.0,  # orbital inclination wrt disk plane [deg]
    "Omega_sat_period": 800.0,  # orbital period [Myr]
    "sat_softening": 2.0,  # satellite softening length [kpc]
    # --- Numerical parameters ---
    "softening": 0.3,  # gravitational softening [kpc]
    "pm_Ng": 128,  # PM grid cells per dimension
    "dt": 1.0,  # time step [Myr]
    "T_total": 2000.0,  # total integration time [Myr]
    "snap_interval": 50.0,  # snapshot interval [Myr]
    "tree_theta": 0.7,  # Barnes-Hut opening angle
    "seed": 42,
    # --- Output ---
    "output_dir": "../output",
}


# ============================================================
# Potential models
# ============================================================


def nfw_potential(r, M_vir, c, R_vir):
    """NFW halo potential."""
    Rs = R_vir / c
    rho0_factor = M_vir / (4.0 * np.pi * Rs**3 * (np.log(1.0 + c) - c / (1.0 + c)))
    x = r / Rs
    phi = -4.0 * np.pi * G_GRAV * rho0_factor * Rs**2 * np.log(1.0 + x) / x
    return phi


def nfw_acceleration(pos, M_vir, c, R_vir):
    """NFW halo acceleration (vectorized)."""
    Rs = R_vir / c
    r = np.sqrt(np.sum(pos**2, axis=1))
    r = np.maximum(r, 1e-10)  # avoid division by zero

    A = M_vir / (np.log(1.0 + c) - c / (1.0 + c))
    x = r / Rs

    # M_enclosed(r) for NFW
    M_enc = A * (np.log(1.0 + x) - x / (1.0 + x))

    # a = -G*M_enc(r)/r^2 * r_hat = -G*M_enc(r)/r^3 * r_vec
    acc_mag = -G_GRAV * M_enc / r**3
    acc = acc_mag[:, np.newaxis] * pos
    return acc


def satellite_position(t, config):
    """
    Compute satellite position on an inclined elliptical orbit.

    Uses a Keplerian orbit with eccentricity to allow close pericentric passages
    that are far more effective at inducing warps than circular orbits.
    """
    period = config["Omega_sat_period"]
    omega_mean = 2.0 * np.pi / period
    R0 = config["R_sat_orbit"]  # semi-major axis
    ecc = config.get("sat_eccentricity", 0.5)
    incl = np.radians(config["incl_sat"])

    # Solve Kepler's equation M = E - e*sin(E) iteratively
    M = omega_mean * t  # mean anomaly
    E = M  # initial guess for eccentric anomaly
    for _ in range(10):
        E = M + ecc * np.sin(E)

    # True anomaly
    cos_f = (np.cos(E) - ecc) / (1 - ecc * np.cos(E))
    sin_f = (np.sqrt(1 - ecc**2) * np.sin(E)) / (1 - ecc * np.cos(E))
    f = np.arctan2(sin_f, cos_f)

    # Orbital radius
    r = R0 * (1 - ecc**2) / (1 + ecc * np.cos(f))

    # Position in orbital plane
    x_orb = r * np.cos(f)
    y_orb = r * np.sin(f)

    # Rotate by inclination around x-axis
    x_sat = x_orb
    y_sat = y_orb * np.cos(incl)
    z_sat = y_orb * np.sin(incl)

    return np.array([x_sat, y_sat, z_sat])


def satellite_acceleration(pos, t, config):
    """Tidal acceleration from satellite (direct + indirect term).

    The indirect term (-a_sat at disk center) ensures we work in the
    non-inertial frame centered on the host galaxy.
    """
    M_sat = config["M_sat"]
    eps = config["sat_softening"]

    r_sat = satellite_position(t, config)  # (3,)

    # Direct term: acceleration of each particle toward satellite
    dr = r_sat[np.newaxis, :] - pos  # (N, 3)
    dist2 = np.sum(dr**2, axis=1) + eps**2
    dist3 = dist2**1.5
    a_direct = G_GRAV * M_sat * dr / dist3[:, np.newaxis]

    # Indirect term: acceleration of galaxy center toward satellite
    dist_sat = np.sqrt(np.sum(r_sat**2) + eps**2)
    a_indirect = -G_GRAV * M_sat * r_sat / dist_sat**3

    return a_direct + a_indirect[np.newaxis, :]


# ============================================================
# Self-gravity via Particle-Mesh (PM) method with FFT
# ============================================================


class PMGrid:
    """
    Particle-Mesh gravity solver using FFT-based Poisson solver.

    This is O(N + Ng^3 log Ng) instead of O(N^2), enabling much
    larger particle counts. Uses Cloud-In-Cell (CIC) interpolation
    for mass assignment and force interpolation.
    """

    def __init__(self, Ng=64, box_size=40.0, softening=0.1):
        """
        Parameters
        ----------
        Ng : int
            Number of grid cells per dimension
        box_size : float
            Full box size in kpc (centered on origin, so [-L/2, L/2])
        softening : float
            Gravitational softening length in kpc
        """
        self.Ng = Ng
        self.L = box_size
        self.dx = box_size / Ng
        self.softening = softening

        # Precompute Green's function in Fourier space
        self._precompute_greens_function()

    def _precompute_greens_function(self):
        """Precompute the Fourier-space Green's function for the Poisson equation."""
        Ng = self.Ng
        dx = self.dx

        # Wave numbers
        kx = np.fft.fftfreq(Ng, d=dx) * 2 * np.pi
        ky = np.fft.fftfreq(Ng, d=dx) * 2 * np.pi
        kz = np.fft.fftfreq(Ng, d=dx) * 2 * np.pi

        KX, KY, KZ = np.meshgrid(kx, ky, kz, indexing="ij")

        # |k|^2 with softening
        k2 = KX**2 + KY**2 + KZ**2
        eps_k = (2 * np.pi * self.softening / self.L) ** 2  # softening in k-space

        # Green's function: -1/k^2 (with softening)
        # Avoid division by zero at k=0
        self.green_k = np.zeros_like(k2)
        nonzero = k2 > 0
        self.green_k[nonzero] = -1.0 / (k2[nonzero] + eps_k)

        # Store wavenumbers for gradient computation
        self.KX = KX
        self.KY = KY
        self.KZ = KZ

    def deposit_mass(self, pos, mass_per_particle, N_particles):
        """
        Deposit particle masses onto grid using Cloud-In-Cell (CIC).
        Returns density field rho(x,y,z).
        """
        Ng = self.Ng
        L = self.L
        dx = self.dx

        # Map positions to grid coordinates [0, Ng)
        # Positions are in [-L/2, L/2], map to [0, L]
        gx = (pos[:, 0] + L / 2) / dx
        gy = (pos[:, 1] + L / 2) / dx
        gz = (pos[:, 2] + L / 2) / dx

        # CIC: each particle contributes to 8 neighboring cells
        ix = np.floor(gx).astype(int)
        iy = np.floor(gy).astype(int)
        iz = np.floor(gz).astype(int)

        # Fractional positions within cell
        fx = gx - ix
        fy = gy - iy
        fz = gz - iz

        # Periodic boundary wrapping
        ix = ix % Ng
        iy = iy % Ng
        iz = iz % Ng
        ix1 = (ix + 1) % Ng
        iy1 = (iy + 1) % Ng
        iz1 = (iz + 1) % Ng

        rho = np.zeros((Ng, Ng, Ng), dtype=np.float64)

        # CIC weights for 8 corners
        total_mass = mass_per_particle * N_particles
        w = mass_per_particle / dx**3  # mass per unit volume

        for (dix, wx), (diy, wy), (diz, wz) in [
            ((ix, 1 - fx), (iy, 1 - fy), (iz, 1 - fz)),
            ((ix1, fx), (iy, 1 - fy), (iz, 1 - fz)),
            ((ix, 1 - fx), (iy1, fy), (iz, 1 - fz)),
            ((ix, 1 - fx), (iy, 1 - fy), (iz1, fz)),
            ((ix1, fx), (iy1, fy), (iz, 1 - fz)),
            ((ix1, fx), (iy, 1 - fy), (iz1, fz)),
            ((ix, 1 - fx), (iy1, fy), (iz1, fz)),
            ((ix1, fx), (iy1, fy), (iz1, fz)),
        ]:
            weight = w * wx * wy * wz
            np.add.at(rho, (dix, diy, diz), weight)

        return rho

    def solve_potential(self, rho):
        """Solve Poisson equation: nabla^2 phi = 4*pi*G*rho."""
        rho_k = np.fft.fftn(rho)
        phi_k = 4.0 * np.pi * G_GRAV * self.green_k * rho_k
        phi = np.real(np.fft.ifftn(phi_k))
        return phi

    def compute_acceleration_field(self, phi):
        """Compute acceleration field a = -grad(phi) using FFT differentiation."""
        phi_k = np.fft.fftn(phi)

        # a_x = -d(phi)/dx = -i*kx*phi_k in Fourier space
        ax = np.real(np.fft.ifftn(-1j * self.KX * phi_k))
        ay = np.real(np.fft.ifftn(-1j * self.KY * phi_k))
        az = np.real(np.fft.ifftn(-1j * self.KZ * phi_k))

        return ax, ay, az

    def interpolate_acceleration(self, pos, ax_field, ay_field, az_field):
        """
        Interpolate acceleration field to particle positions using CIC.
        """
        Ng = self.Ng
        L = self.L
        dx = self.dx
        N = len(pos)

        # Map positions to grid coordinates
        gx = (pos[:, 0] + L / 2) / dx
        gy = (pos[:, 1] + L / 2) / dx
        gz = (pos[:, 2] + L / 2) / dx

        # CIC interpolation: same weights as deposit, but read instead of write
        ix = np.floor(gx).astype(int)
        iy = np.floor(gy).astype(int)
        iz = np.floor(gz).astype(int)

        fx = gx - ix
        fy = gy - iy
        fz = gz - iz

        ix = ix % Ng
        iy = iy % Ng
        iz = iz % Ng
        ix1 = (ix + 1) % Ng
        iy1 = (iy + 1) % Ng
        iz1 = (iz + 1) % Ng

        acc = np.zeros((N, 3), dtype=np.float64)

        for field_idx, field in enumerate([ax_field, ay_field, az_field]):
            val = np.zeros(N)
            for (dix, wx), (diy, wy), (diz, wz) in [
                ((ix, 1 - fx), (iy, 1 - fy), (iz, 1 - fz)),
                ((ix1, fx), (iy, 1 - fy), (iz, 1 - fz)),
                ((ix, 1 - fx), (iy1, fy), (iz, 1 - fz)),
                ((ix, 1 - fx), (iy, 1 - fy), (iz1, fz)),
                ((ix1, fx), (iy1, fy), (iz, 1 - fz)),
                ((ix1, fx), (iy, 1 - fy), (iz1, fz)),
                ((ix, 1 - fx), (iy1, fy), (iz1, fz)),
                ((ix1, fx), (iy1, fy), (iz1, fz)),
            ]:
                val += field[dix, diy, diz] * wx * wy * wz
            acc[:, field_idx] = val

        return acc

    def compute_self_gravity(self, pos, mass_per_particle, N_particles):
        """
        Full PM self-gravity computation pipeline.

        Returns acceleration array (N, 3).
        """
        # 1. Deposit mass onto grid
        rho = self.deposit_mass(pos, mass_per_particle, N_particles)

        # 2. Solve Poisson equation
        phi = self.solve_potential(rho)

        # 3. Compute acceleration field
        ax, ay, az = self.compute_acceleration_field(phi)

        # 4. Interpolate to particle positions
        acc = self.interpolate_acceleration(pos, ax, ay, az)

        return acc


# Global PM grid instance (lazily initialized)
_pm_grid = None


def get_pm_grid(config):
    """Get or create the PM grid solver."""
    global _pm_grid
    if _pm_grid is None:
        # Box must be large enough to contain the disk
        box_size = 2.2 * config["R_trunc"]  # some padding beyond truncation
        Ng = config.get("pm_Ng", 128)  # default 128^3 grid
        _pm_grid = PMGrid(Ng=Ng, box_size=box_size, softening=config["softening"])
    return _pm_grid


def self_gravity_pm(pos, mass_per_particle, config):
    """Compute self-gravity using Particle-Mesh method."""
    pm = get_pm_grid(config)

    # Clip particles to box (particles that escape the box are ignored)
    L2 = pm.L / 2 - pm.dx
    mask = (
        (np.abs(pos[:, 0]) < L2) & (np.abs(pos[:, 1]) < L2) & (np.abs(pos[:, 2]) < L2)
    )

    if np.sum(mask) < 100:
        return np.zeros_like(pos)

    # Compute on in-box particles
    acc_full = np.zeros_like(pos)
    pos_in = pos[mask]

    acc_in = pm.compute_self_gravity(pos_in, mass_per_particle, len(pos_in))
    acc_full[mask] = acc_in

    return acc_full


# ============================================================
# Initial conditions: exponential disk + isothermal z-profile
# ============================================================


def _vertical_frequency_MN(R_cyl, M_disk, a, b, M_vir, c_halo, R_vir):
    """
    Compute vertical oscillation frequency nu_z at radius R in the midplane
    for combined Miyamoto-Nagai disk + NFW halo.

    nu_z^2 = d^2(Phi)/dz^2 |_{z=0}
    """
    # Miyamoto-Nagai: d^2 phi/dz^2 at z=0
    # phi_MN = -G*M / sqrt(R^2 + (a+sqrt(z^2+b^2))^2)
    # At z=0: (a+b), denom = (R^2 + (a+b)^2)^(1/2)
    ab = a + b
    denom2 = R_cyl**2 + ab**2
    denom = np.sqrt(denom2)

    # d^2 phi / dz^2 at z=0 for MN
    # = G*M * (a+b)/b * [1/(R^2+(a+b)^2)^(3/2) - 3*(a+b)^2 / ((R^2+(a+b)^2)^(5/2) * b)]
    # Simplified:
    # nu_z^2 = G*M * (a+b) / (b * denom^3) * (1 - 3*(a+b)^2/denom^2 * ... )
    # Let me use numerical differentiation instead for robustness
    dz = 0.001  # kpc
    pos_0 = np.column_stack([R_cyl, np.zeros_like(R_cyl), np.zeros_like(R_cyl)])
    pos_p = np.column_stack([R_cyl, np.zeros_like(R_cyl), np.full_like(R_cyl, dz)])
    pos_m = np.column_stack([R_cyl, np.zeros_like(R_cyl), np.full_like(R_cyl, -dz)])

    # MN acceleration z-component at z=+dz and z=-dz
    az_p = miyamoto_nagai_acceleration(pos_p, M_disk, a, b)[:, 2]
    az_m = miyamoto_nagai_acceleration(pos_m, M_disk, a, b)[:, 2]

    # NFW acceleration z-component
    az_nfw_p = nfw_acceleration(pos_p, M_vir, c_halo, R_vir)[:, 2]
    az_nfw_m = nfw_acceleration(pos_m, M_vir, c_halo, R_vir)[:, 2]

    # Total d(a_z)/dz ≈ (a_z(+dz) - a_z(-dz)) / (2*dz)
    # nu_z^2 = -d(a_z)/dz (since a_z = -d(phi)/dz, so d(a_z)/dz = -d^2(phi)/dz^2 = -nu_z^2)
    daz_dz = ((az_p + az_nfw_p) - (az_m + az_nfw_m)) / (2 * dz)
    nu_z_sq = -daz_dz
    nu_z_sq = np.maximum(nu_z_sq, 1e-20)  # ensure positive

    return np.sqrt(nu_z_sq)


def generate_disk_ic(config):
    """
    Generate initial conditions for a stellar disk in equilibrium with
    the combined Miyamoto-Nagai disk + NFW halo potential.

    Surface density: Sigma(R) = Sigma_0 * exp(-R/Rd)
    Vertical distribution: consistent with the total potential
    Velocities: circular velocity + correct velocity dispersions from Jeans equations
    """
    rng = np.random.default_rng(config["seed"])
    N = config["N_particles"]
    Rd = config["R_disk"]
    z0 = config["z_disk"]
    M_disk = config["M_disk"]
    R_trunc = config["R_trunc"]

    # --- Sample radial positions from exponential distribution ---
    R = np.zeros(N)
    count = 0
    while count < N:
        R_cand = rng.exponential(Rd, size=N * 2)
        mask = R_cand < R_trunc
        R_cand = R_cand[mask]
        n_take = min(len(R_cand), N - count)
        R[count : count + n_take] = R_cand[:n_take]
        count += n_take

    # Azimuthal angles
    phi = rng.uniform(0, 2 * np.pi, N)

    # Convert to Cartesian
    x = R * np.cos(phi)
    y = R * np.sin(phi)

    # --- Compute vertical frequency and set z-distribution consistently ---
    M_vir = config["M_halo_vir"]
    c = config["c_halo"]
    R_vir = config["R_vir"]
    a_MN = config["R_disk"]
    b_MN = config["z_disk"]

    nu_z = _vertical_frequency_MN(R, M_disk, a_MN, b_MN, M_vir, c, R_vir)

    # Vertical positions: Gaussian distribution with scale height z0 (= b of MN potential)
    # For a near-harmonic potential, z ~ N(0, z0) is a good approximation
    z = rng.normal(0, z0, N)

    pos = np.column_stack([x, y, z])

    # --- Compute circular velocities from total potential ---
    # v_circ^2 = R * |d(Phi)/dR| at z=0
    # Use the actual acceleration from both potentials
    pos_midplane = np.column_stack([R, np.zeros(N), np.zeros(N)])

    # Radial direction acceleration at midplane
    a_MN_mid = miyamoto_nagai_acceleration(pos_midplane, M_disk, a_MN, b_MN)
    a_NFW_mid = nfw_acceleration(pos_midplane, M_vir, c, R_vir)
    a_total_R = a_MN_mid + a_NFW_mid

    # a_R = a_x for particles on x-axis (pos_midplane has y=0, z=0)
    # v_circ^2 = -R * a_R
    a_radial = a_total_R[:, 0]  # x-component = radial for particles on x-axis
    v_circ_sq = -R * a_radial
    v_circ_sq = np.maximum(v_circ_sq, 0.0)
    v_circ = np.sqrt(v_circ_sq)

    # Vertical velocity dispersion: sigma_vz = nu_z * z0
    sigma_vz = nu_z * z0
    sigma_vz = np.maximum(sigma_vz, 1e-5)

    # Radial velocity dispersion: Toomre Q ~ 1.5
    # sigma_R ~ sigma_vz * kappa / nu_z ~ sigma_vz for flat rotation curve
    sigma_R = sigma_vz * 1.4  # slightly warmer radially than vertically
    sigma_R = np.maximum(sigma_R, 1e-5)

    # Asymmetric drift correction
    v_mean = np.sqrt(np.maximum(v_circ_sq - sigma_R**2, 0.1 * v_circ_sq))

    # Velocity components (cylindrical -> Cartesian)
    vR = rng.normal(0, sigma_R)
    vphi = v_mean + rng.normal(0, sigma_R * 0.7)  # azimuthal dispersion ~ 0.7 * sigma_R

    vx = -vphi * np.sin(phi) + vR * np.cos(phi)
    vy = vphi * np.cos(phi) + vR * np.sin(phi)
    vz = rng.normal(0, sigma_vz)

    vel = np.column_stack([vx, vy, vz])

    # Particle mass
    mass_per_particle = M_disk / N

    return pos, vel, mass_per_particle


# ============================================================
# Leapfrog integrator (KDK scheme)
# ============================================================


def miyamoto_nagai_acceleration(pos, M_disk, a, b):
    """
    Miyamoto-Nagai disk potential acceleration.

    phi = -G*M / sqrt(R^2 + (a + sqrt(z^2 + b^2))^2)

    This provides a smooth, analytic disk potential that doesn't suffer
    from numerical heating issues of particle-based self-gravity.

    Parameters
    ----------
    pos : (N, 3) array
    M_disk : disk mass [1e10 Msun]
    a : radial scale length [kpc]
    b : vertical scale height [kpc]
    """
    x, y, z = pos[:, 0], pos[:, 1], pos[:, 2]
    R2 = x**2 + y**2
    zb = np.sqrt(z**2 + b**2)
    azb = a + zb
    denom = (R2 + azb**2) ** 1.5

    ax = -G_GRAV * M_disk * x / denom
    ay = -G_GRAV * M_disk * y / denom
    az = -G_GRAV * M_disk * z * azb / (zb * denom)

    return np.column_stack([ax, ay, az])


def compute_total_acceleration(pos, mass_per_particle, t, config, include_tidal=True):
    """
    Compute total acceleration on all particles.

    Uses analytic potentials for halo and disk (avoiding numerical heating),
    plus the external tidal field from the satellite.
    """
    # 1. NFW halo
    acc = nfw_acceleration(pos, config["M_halo_vir"], config["c_halo"], config["R_vir"])

    # 2. Analytic disk potential (Miyamoto-Nagai)
    # a ~ R_disk, b ~ z_disk are the scale parameters
    acc += miyamoto_nagai_acceleration(
        pos, config["M_disk"], config["R_disk"], config["z_disk"]
    )

    # 3. Tidal force from satellite
    if include_tidal:
        acc += satellite_acceleration(pos, t, config)

    return acc


def run_simulation(config, include_tidal=True, label="tidal"):
    """
    Main simulation loop.

    Parameters
    ----------
    config : dict
        Simulation configuration
    include_tidal : bool
        Whether to include the satellite tidal field
    label : str
        Label for output files
    """
    # Reset PM grid for fresh simulation
    global _pm_grid
    _pm_grid = None

    print(f"=" * 60)
    print(f"Starting simulation: {label}")
    print(
        f"N={config['N_particles']}, T={config['T_total']} Myr, dt={config['dt']} Myr"
    )
    print(f"Tidal perturbation: {'ON' if include_tidal else 'OFF'}")
    print(f"=" * 60)

    # Generate initial conditions
    pos, vel, mass_pp = generate_disk_ic(config)

    dt = config["dt"]
    T = config["T_total"]
    snap_dt = config["snap_interval"]
    N_steps = int(T / dt)

    out_dir = os.path.join(config["output_dir"], label)
    os.makedirs(out_dir, exist_ok=True)

    # Save config
    with open(os.path.join(out_dir, "config.json"), "w") as f:
        json.dump(config, f, indent=2)

    # Initial acceleration (for leapfrog)
    t = 0.0
    acc = compute_total_acceleration(pos, mass_pp, t, config, include_tidal)

    # Storage for snapshots
    snap_times = []
    snap_count = 0

    # Save IC
    np.savez_compressed(
        os.path.join(out_dir, f"snap_{snap_count:04d}.npz"),
        pos=pos,
        vel=vel,
        t=t,
        mass_pp=mass_pp,
    )
    snap_times.append(t)
    snap_count += 1

    t_wall_start = walltime.time()

    for step in range(1, N_steps + 1):
        # Kick (half step)
        vel += 0.5 * dt * acc

        # Drift (full step)
        pos += dt * vel

        # Update time
        t = step * dt

        # Compute new acceleration
        acc = compute_total_acceleration(pos, mass_pp, t, config, include_tidal)

        # Kick (half step)
        vel += 0.5 * dt * acc

        # Save snapshot
        if step % int(snap_dt / dt) == 0:
            np.savez_compressed(
                os.path.join(out_dir, f"snap_{snap_count:04d}.npz"),
                pos=pos.copy(),
                vel=vel.copy(),
                t=t,
                mass_pp=mass_pp,
            )
            snap_times.append(t)
            snap_count += 1

            elapsed = walltime.time() - t_wall_start
            eta = elapsed / step * (N_steps - step)
            print(
                f"  t = {t:.0f} Myr  |  snap {snap_count}  |  "
                f"wall: {elapsed:.1f}s  |  ETA: {eta:.1f}s  |  "
                f"z_rms = {np.std(pos[:, 2]):.4f} kpc"
            )

    # Save snapshot index
    np.save(os.path.join(out_dir, "snap_times.npy"), np.array(snap_times))

    total_wall = walltime.time() - t_wall_start
    print(f"\nSimulation '{label}' complete. Wall time: {total_wall:.1f}s")
    print(f"Saved {snap_count} snapshots to {out_dir}/")

    return out_dir


# ============================================================
# Entry point
# ============================================================
if __name__ == "__main__":
    import sys

    config = DEFAULT_CONFIG.copy()

    # Allow quick test with fewer particles
    if "--test" in sys.argv:
        config["N_particles"] = 5000
        config["T_total"] = 500.0
        config["dt"] = 2.0
        config["snap_interval"] = 50.0
        print("*** TEST MODE: reduced resolution ***")

    if "--fast" in sys.argv:
        config["N_particles"] = 20000
        config["T_total"] = 2000.0
        config["dt"] = 2.0
        config["snap_interval"] = 50.0
        print("*** FAST MODE ***")

    # Run baseline (isolated disk)
    if "--baseline" in sys.argv or "--all" in sys.argv:
        run_simulation(config, include_tidal=False, label="baseline")

    # Run tidal simulation
    if "--tidal" in sys.argv or "--all" in sys.argv:
        run_simulation(config, include_tidal=True, label="tidal")

    # Default: run both
    if len(sys.argv) == 1:
        run_simulation(config, include_tidal=False, label="baseline")
        run_simulation(config, include_tidal=True, label="tidal")
