# Tidal-Induced Stellar Warp: Ring Model Simulation Study

## Objective
Investigate the formation mechanism of stellar disk warps induced by tidal interactions 
with a satellite galaxy, using a tilted-ring dynamical model.

## Method
We employ a tilted-ring model (Sparke & Casertano 1988 formalism) where the stellar disk 
is decomposed into 60 concentric rings (R = 0.5--18 kpc). Each ring evolves under three torques:

1. **Self-gravity coupling**: Gravitational torques between misaligned rings resist 
   differential tilting and propagate bending waves inward (ring-ring interaction, 
   softened by disk thickness h = 0.3 kpc, coupling up to 8 neighbors).
2. **Halo torque**: A slightly oblate NFW halo (q = 0.95, M_vir = 10^12 Msun, c = 12) 
   causes differential precession of tilted rings.
3. **Tidal torque**: The quadrupole tidal field from an orbiting satellite 
   (M_sat = 10^11 Msun, semi-major axis a = 25 kpc, e = 0.3, inclination i = 45 deg, 
   period P = 800 Myr) drives warp excitation.

Integration: 4th-order Runge-Kutta, dt = 1 Myr, T_total = 3000 Myr.

## Key Results

| Metric | Value |
|--------|-------|
| Peak warp tilt (fiducial) | 1.70 deg |
| Final outer disk tilt (R > 10 kpc) | 1.34 deg |
| Baseline max tilt | 0.00 deg |
| Warp growth pattern | Stepwise, synchronized with pericentric passages |

### 1. Warp Formation Mechanism (Fig. 1)

The simulation reveals a clear **impulsive tidal excitation** mechanism:

- **Panel (a)**: The warp tilt profile grows monotonically with radius at all times, 
  confirming that the outer disk is most susceptible to tidal torquing (T_tidal ~ R^2).
  The inner disk (R < 4 kpc) remains nearly flat due to the disk's bending stiffness 
  (self-gravity coupling dominates over tidal torque there).

- **Panel (b)**: The line of nodes (phase of the tilt vector) shows systematic precession 
  over time, rotating from ~ -22 deg to ~ -2 deg across 3 Gyr. This precession is 
  driven by the combination of the oblate halo torque and the orbital motion of the 
  satellite. Importantly, the line of nodes is nearly constant with radius at any given 
  time, indicating a coherent (non-wound-up) warp -- consistent with the "modified tilt" 
  mode predicted by Sparke & Casertano (1988).

- **Panel (c)**: The warp amplitude grows in a characteristic **staircase pattern**, with 
  each step corresponding to a pericentric passage of the satellite (period ~ 800 Myr, 
  marked by vertical dashed lines). Between passages, the amplitude plateaus. This 
  demonstrates that warp growth is impulsive rather than continuous, with each close 
  encounter delivering a discrete "kick" of tidal torque.

### 2. Torque Decomposition (Fig. 7)

The torque analysis at t = 3000 Myr reveals a clear hierarchy:

- **Tidal torque** dominates at all radii R > 2 kpc (by 1--3 orders of magnitude over 
  self-gravity), scaling as ~ R^2 and reaching ~ 10^-8 at R = 18 kpc.
- **Self-gravity coupling** is significant only in the inner disk (R < 5 kpc), where it 
  acts as a bending stiffness that resists differential tilting and keeps the inner disk 
  aligned. The characteristic bump at R ~ 3--5 kpc reflects the peak of the disk surface 
  density (exponential scale length R_d = 3.5 kpc).
- **Halo torque** is subdominant by ~3 orders of magnitude, explaining why the q_halo 
  parameter study showed negligible sensitivity. This is expected: the halo torque scales 
  as (1 - q^2) ~ 0.1, while the tidal torque from a 10^11 Msun satellite at 20 kpc 
  pericenter is intrinsically much stronger.

### 3. Parameter Dependencies

**Satellite mass (Fig. 4)**: Warp amplitude scales nearly linearly with M_sat over the 
range 1--50 x 10^10 Msun (0.17 deg to 8.55 deg peak tilt). This is expected from the 
tidal torque formula T ~ G * M_sat * M_ring * R^2 / d^3. The time evolution panel shows 
that more massive perturbers produce faster warp growth with the same staircase pattern.

**Orbital inclination (Fig. 5)**: Peak warp tilt peaks sharply at i ~ 45 deg (1.70 deg) 
and drops to near zero at i = 90 deg (0.18 deg). This is the theoretically predicted 
sin(2i) dependence of the z-component of the tidal torque: the tidal field must have 
a component both along and perpendicular to the disk plane to generate a net tilting 
torque. At i = 90 deg, the satellite orbit is polar and the time-averaged torque nearly 
cancels; at i = 0 deg, the orbit is coplanar and no out-of-plane torque exists. The 
optimal angle is sin(2i) maximized at i = 45 deg. The time evolution panel also reveals 
interesting behavior: i = 90 deg shows oscillatory warp (the torque reverses sign each 
half-orbit), while intermediate inclinations produce steady monotonic growth.

**Halo flattening**: No significant sensitivity was found for q = 0.8--1.0 in our 
parameter regime. This is because the halo torque is 3 orders of magnitude weaker than 
the tidal torque. Halo flattening would become relevant only for much weaker perturbers 
or after the satellite is removed (the halo torque would then control the warp precession 
and long-term evolution/winding).

## Physical Interpretation

The warp formation in our simulation follows a three-stage process:

1. **Tidal excitation**: During each pericentric passage (d_peri = 17.5 kpc), the 
   satellite delivers a strong quadrupole tidal torque that preferentially tilts the 
   outer disk rings (T ~ R^2 / d^3).

2. **Bending wave propagation**: Self-gravity coupling between adjacent rings 
   transmits the perturbation inward as a bending wave. However, the wave speed is 
   limited by the disk surface density, so the inner disk responds more slowly than 
   the outer disk, creating the characteristic radially increasing tilt profile.

3. **Coherent precession**: The entire warp precesses as a nearly coherent structure 
   (constant line of nodes with radius), maintained by the combined action of 
   self-gravity coupling and the halo potential. This is consistent with the "modified 
   tilt mode" of Sparke & Casertano (1988).

## Relevance to Observations

- **Milky Way warp**: The observed MW warp has an amplitude of ~1--2 deg at R ~ 15 kpc, 
  consistent with our fiducial model (1.34 deg outer tilt). The LMC (M ~ 1--2 x 10^11 Msun, 
  d ~ 50 kpc) is a plausible driver (Tsuchiya 2002; Bennett+ 2021).

- **S-shaped warp morphology**: Our model produces the characteristic radially increasing 
  tilt observed in most warped galaxies (Garcia-Ruiz+ 2002).

- **Warp precession**: The predicted precession of the line of nodes (~20 deg over 3 Gyr) 
  is consistent with recent measurements of the MW warp precession (Huang+ 2024).

## Conclusions

The tilted-ring simulation confirms that **tidal torquing by a satellite galaxy is an 
efficient and physically well-motivated mechanism for generating stellar disk warps**. 
The key findings are:

1. Warp growth is impulsive, driven by discrete pericentric passages
2. The outer disk is preferentially warped (T ~ R^2), while the inner disk is stabilized 
   by self-gravity coupling
3. Warp amplitude scales linearly with satellite mass and peaks at orbital inclination ~ 45 deg
4. The warp precesses coherently, maintained by disk self-gravity
5. Halo flattening is subdominant in the presence of a strong tidal perturber
