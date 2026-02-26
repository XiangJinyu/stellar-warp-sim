# Tidal-Induced Stellar Warp Simulation

A tilted-ring dynamical model for studying the formation mechanism of stellar disk warps driven by tidal interactions with satellite galaxies.

## Physics

Galactic disk warps — where the outer disk bends out of the midplane — are ubiquitous in spiral galaxies. This project investigates how tidal torques from an orbiting satellite excite and sustain these warps, using the Sparke & Casertano (1988) tilted-ring formalism.

The disk is decomposed into 60 concentric rings (R = 0.5–18 kpc). Each ring's angular momentum vector evolves under three torques:

| Torque | Source | Role |
|--------|--------|------|
| **Self-gravity coupling** | Ring–ring gravitational interaction | Bending stiffness; propagates warp inward |
| **Halo torque** | Oblate NFW dark matter halo | Differential precession of tilted rings |
| **Tidal torque** | Orbiting satellite (quadrupole field) | Drives warp excitation |

## Key Results

- **Impulsive excitation**: Warp grows in a staircase pattern, with each step synchronized to a pericentric passage of the satellite
- **Radial structure**: Outer disk is preferentially warped (T ∝ R²); inner disk stabilized by self-gravity
- **Mass scaling**: Warp amplitude scales linearly with satellite mass (0.17°–8.55° for M_sat = 10¹⁰–5×10¹¹ M☉)
- **Inclination dependence**: Peak warp at orbital inclination ~45° (sin 2i dependence), near-zero at 90°
- **Coherent precession**: Line of nodes precesses at ~7 deg/Gyr, consistent with Gaia measurements (Poggio+ 2020)
- **Torque hierarchy**: Tidal ≫ self-gravity ≫ halo (at R > 2 kpc)

Fiducial model: peak tilt = 1.70°, outer disk tilt = 1.34° — at the lower end of the observed Milky Way stellar warp (2–5°; Chen+ 2019, Skowron+ 2019). See [REPORT.md](REPORT.md) for full discussion of caveats including missing halo dynamical response.

## Project Structure

```
stellar_warp_sim/
├── src/
│   ├── ring_model.py          # Core simulation: tilted-ring model, RK4 integrator
│   ├── run_ring_study.py      # Master runner: baseline + tidal + parameter study + analysis + plots
│   ├── simulation.py          # Earlier N-body attempt (deprecated, kept for reference)
│   ├── analysis.py            # N-body analysis utilities
│   └── visualize.py           # N-body visualization utilities
├── output/                    # Simulation snapshots (~2500 .npz files across 20 runs)
│   ├── baseline/              # Isolated disk (no tidal field)
│   ├── tidal/                 # Fiducial tidal simulation
│   ├── param_Msat_*/          # Satellite mass parameter study (6 values)
│   ├── param_incl_*/          # Orbital inclination study (7 values)
│   └── param_qhalo_*/         # Halo flattening study (5 values)
├── figures/                   # 7 publication-quality figures
├── paper/                     # Full MNRAS-format LaTeX manuscript + compiled PDF
├── REPORT.md                  # Full English report with embedded figures
├── REPORT_CN.md               # Full Chinese report (中文报告)
├── RESULTS.json               # Quantitative results for all runs
└── README.md
```

## Figures

| Figure | Description |
|--------|-------------|
| [fig1](figures/fig1_warp_evolution.png) | Warp tilt profile, line of nodes, and amplitude evolution |
| [fig2](figures/fig2_3d_warp.png) | 3D visualization of warped disk at four epochs |
| [fig3](figures/fig3_disk_views.png) | Face-on and edge-on projected density maps |
| [fig4](figures/fig4_param_Msat.png) | Parameter study: satellite mass |
| [fig5](figures/fig5_param_inclination.png) | Parameter study: orbital inclination |
| [fig6](figures/fig6_param_qhalo.png) | Parameter study: halo flattening |
| [fig7](figures/fig7_torque_decomposition.png) | Torque decomposition (tidal vs self-gravity vs halo) |

## Quick Start

```bash
# Run the full study (baseline + tidal + parameter study + analysis + figures)
cd src
python run_ring_study.py

# Or run individual simulations
python ring_model.py --all          # baseline + fiducial tidal
python ring_model.py --test --tidal # quick test (500 Myr)
```

### Requirements

- Python 3.8+
- NumPy, SciPy, Matplotlib

No GPU or special hardware required. The full study (20 simulations × 3 Gyr) completes in ~2 minutes on a laptop.

## Units

| Quantity | Unit |
|----------|------|
| Length | kpc |
| Mass | 10¹⁰ M☉ |
| Time | Myr |
| G | 4.498 × 10⁻⁶ kpc³ / (10¹⁰ M☉ · Myr²) |

## Fiducial Model Parameters

| Parameter | Value | Description |
|-----------|-------|-------------|
| M_disk | 5 × 10¹⁰ M☉ | Disk mass |
| R_disk | 3.5 kpc | Exponential scale length |
| M_halo | 10¹² M☉ | NFW halo virial mass |
| q_halo | 0.95 | Halo axis ratio |
| M_sat | 10¹¹ M☉ | Satellite mass |
| a_sat | 25 kpc | Satellite semi-major axis |
| e_sat | 0.3 | Orbital eccentricity |
| i_sat | 45° | Orbital inclination |
| P_sat | 800 Myr | Orbital period |

## References

- Sparke, L. S. & Casertano, S. 1988, MNRAS, 234, 873
- Chen, X., et al. 2019, Nature Astronomy, 3, 320
- Skowron, D. M., et al. 2019, Science, 365, 478
- Poggio, E., et al. 2020, Nature Astronomy, 4, 590
- Sellwood, J. A. & Debattista, V. P. 2021, MNRAS, 510, 2532
- Vasiliev, E. 2023, Galaxies, 11, 59

See [REPORT.md](REPORT.md) for the complete reference list and detailed analysis.

## Paper

A full manuscript in MNRAS format is included in the `paper/` directory:

- **[Compiled PDF](paper/tidal_stellar_warp_paper.pdf)** -- ready-to-read version
- `paper/main.tex` -- main LaTeX source (uses `\input` for each section)
- `paper/sec_introduction.tex` -- Introduction: observational context, 4 warp mechanisms, scope
- `paper/sec_method.tex` -- Method: tilted-ring formalism, NFW halo, 3 torque derivations
- `paper/sec_results.tex` -- Results: fiducial model, torque decomposition, 3 parameter studies
- `paper/sec_discussion.tex` -- Discussion: MW comparison, precession, caveats
- `paper/sec_conclusions.tex` -- Conclusions: 7 key findings
- `paper/references.bib` -- 25 bibliography entries

To compile locally:
```bash
cd paper
pdflatex main && bibtex main && pdflatex main && pdflatex main
```

## License

This project is released for academic and educational use.
