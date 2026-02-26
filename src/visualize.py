#!/usr/bin/env python3
"""
Visualization module for Tidal-Induced Stellar Warp Simulation
===============================================================
Generates publication-quality figures.
"""

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import os
import sys

# Add src to path
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from analysis import (
    load_all_snapshots,
    load_snapshot,
    warp_profile,
    fourier_bending_modes,
    compute_angular_momentum_profile,
    warp_evolution_summary,
)
from simulation import satellite_position, DEFAULT_CONFIG


def setup_style():
    """Setup publication-quality plot style."""
    plt.rcParams.update(
        {
            "font.size": 12,
            "axes.labelsize": 14,
            "axes.titlesize": 14,
            "legend.fontsize": 10,
            "xtick.labelsize": 11,
            "ytick.labelsize": 11,
            "figure.dpi": 150,
            "savefig.dpi": 150,
            "savefig.bbox": "tight",
        }
    )


# ============================================================
# Figure 1: Edge-on and face-on views at selected times
# ============================================================


def plot_disk_views(sim_dir, config, fig_dir, times_to_plot=None):
    """Plot face-on and edge-on views of the disk at selected times."""
    setup_style()
    snapshots = load_all_snapshots(sim_dir)

    if times_to_plot is None:
        times_to_plot = [0, 500, 1000, 1500, 2000]

    # Find closest snapshots
    all_times = np.array([s["t"] for s in snapshots])

    fig, axes = plt.subplots(2, len(times_to_plot), figsize=(4 * len(times_to_plot), 8))

    for j, t_target in enumerate(times_to_plot):
        idx = np.argmin(np.abs(all_times - t_target))
        snap = snapshots[idx]
        pos = snap["pos"]
        t = snap["t"]

        x, y, z = pos[:, 0], pos[:, 1], pos[:, 2]

        # Face-on (x-y)
        ax = axes[0, j]
        ax.hist2d(
            x, y, bins=200, range=[[-18, 18], [-18, 18]], norm=LogNorm(), cmap="inferno"
        )

        # Plot satellite position
        r_sat = satellite_position(t, config)
        ax.plot(
            r_sat[0], r_sat[1], "c*", ms=12, mew=1.5, markeredgecolor="white", zorder=10
        )

        ax.set_xlim(-18, 18)
        ax.set_ylim(-18, 18)
        ax.set_aspect("equal")
        ax.set_title(f"t = {t:.0f} Myr")
        if j == 0:
            ax.set_ylabel("y [kpc]")
        ax.set_xlabel("x [kpc]")

        # Edge-on (x-z)
        ax = axes[1, j]
        ax.hist2d(
            x,
            z,
            bins=[200, 100],
            range=[[-18, 18], [-5, 5]],
            norm=LogNorm(),
            cmap="inferno",
        )
        ax.set_xlim(-18, 18)
        ax.set_ylim(-5, 5)
        ax.set_aspect("equal")
        if j == 0:
            ax.set_ylabel("z [kpc]")
        ax.set_xlabel("x [kpc]")

    fig.suptitle("Disk Evolution: Face-on (top) and Edge-on (bottom)", y=1.02)
    plt.tight_layout()
    plt.savefig(os.path.join(fig_dir, "fig1_disk_views.png"))
    plt.close()
    print("  Saved fig1_disk_views.png")


# ============================================================
# Figure 2: Warp profile evolution
# ============================================================


def plot_warp_profiles(sim_dir, fig_dir, label="tidal"):
    """Plot warp amplitude and phase as functions of R at different times."""
    setup_style()
    summary = warp_evolution_summary(sim_dir)

    fig, axes = plt.subplots(1, 3, figsize=(16, 5))

    # Select time snapshots to plot
    times = summary["times"]
    indices = np.linspace(0, len(times) - 1, 8, dtype=int)

    cmap = plt.cm.viridis
    colors = [cmap(i / len(indices)) for i in range(len(indices))]

    for k, idx in enumerate(indices):
        wp = summary["profiles"][idx]
        t = times[idx]
        c = colors[k]

        # Warp amplitude (m=1)
        axes[0].plot(wp["R_mid"], wp["warp_amplitude"], color=c, label=f"t={t:.0f}")

        # Mean z displacement
        axes[1].plot(wp["R_mid"], wp["z_mean"], color=c, label=f"t={t:.0f}")

        # Warp phase (line of nodes)
        # Only plot where amplitude is significant
        sig_mask = wp["warp_amplitude"] > 0.01
        if np.any(sig_mask):
            axes[2].plot(
                wp["R_mid"][sig_mask],
                np.degrees(wp["warp_phase"][sig_mask]),
                "o-",
                color=c,
                ms=3,
                label=f"t={t:.0f}",
            )

    axes[0].set_xlabel("R [kpc]")
    axes[0].set_ylabel("Warp amplitude (m=1) [kpc]")
    axes[0].set_title("Warp Amplitude")
    axes[0].legend(fontsize=8, ncol=2)
    axes[0].set_xlim(0, 15)

    axes[1].set_xlabel("R [kpc]")
    axes[1].set_ylabel("<z> [kpc]")
    axes[1].set_title("Mean Vertical Displacement")
    axes[1].axhline(0, color="gray", ls="--", lw=0.5)
    axes[1].set_xlim(0, 15)

    axes[2].set_xlabel("R [kpc]")
    axes[2].set_ylabel("Warp phase [deg]")
    axes[2].set_title("Line of Nodes")
    axes[2].set_xlim(0, 15)
    axes[2].legend(fontsize=8, ncol=2)

    plt.suptitle(f"Warp Profile Evolution ({label})", y=1.02)
    plt.tight_layout()
    plt.savefig(os.path.join(fig_dir, f"fig2_warp_profiles_{label}.png"))
    plt.close()
    print(f"  Saved fig2_warp_profiles_{label}.png")


# ============================================================
# Figure 3: Warp amplitude time evolution
# ============================================================


def plot_warp_time_evolution(tidal_dir, baseline_dir, fig_dir):
    """Compare warp amplitude evolution between tidal and baseline runs."""
    setup_style()

    fig, axes = plt.subplots(1, 3, figsize=(16, 5))

    for sim_dir, label, color in [
        (tidal_dir, "Tidal", "C0"),
        (baseline_dir, "Isolated", "C1"),
    ]:
        if not os.path.exists(sim_dir):
            continue
        summary = warp_evolution_summary(sim_dir)

        axes[0].plot(
            summary["times"],
            summary["max_warp_amplitude"],
            color=color,
            label=label,
            lw=2,
        )
        axes[1].plot(
            summary["times"], summary["mean_z_outer"], color=color, label=label, lw=2
        )
        axes[2].plot(summary["times"], summary["z_rms"], color=color, label=label, lw=2)

    axes[0].set_xlabel("Time [Myr]")
    axes[0].set_ylabel("Max warp amplitude [kpc]")
    axes[0].set_title("Peak m=1 Bending Amplitude")
    axes[0].legend()

    axes[1].set_xlabel("Time [Myr]")
    axes[1].set_ylabel("<|z|> (R>8 kpc) [kpc]")
    axes[1].set_title("Mean Outer Disk Displacement")
    axes[1].legend()

    axes[2].set_xlabel("Time [Myr]")
    axes[2].set_ylabel("z_rms [kpc]")
    axes[2].set_title("Vertical Disk Thickness")
    axes[2].legend()

    plt.tight_layout()
    plt.savefig(os.path.join(fig_dir, "fig3_warp_time_evolution.png"))
    plt.close()
    print("  Saved fig3_warp_time_evolution.png")


# ============================================================
# Figure 4: Fourier bending modes
# ============================================================


def plot_bending_modes(sim_dir, fig_dir, label="tidal"):
    """Plot m=0,1,2,3 bending mode amplitudes at final snapshot."""
    setup_style()

    snapshots = load_all_snapshots(sim_dir)

    # Plot at 3 times: early, mid, final
    times_frac = [0.0, 0.5, 1.0]

    fig, axes = plt.subplots(1, 3, figsize=(16, 5))

    for ax_idx, frac in enumerate(times_frac):
        snap_idx = int(frac * (len(snapshots) - 1))
        snap = snapshots[snap_idx]
        pos = snap["pos"]
        t = snap["t"]

        modes = fourier_bending_modes(pos)

        for m in range(5):
            style = ["-", "-", "--", ":", "-."][m]
            axes[ax_idx].plot(
                modes["R_mid"], modes["amplitudes"][m], ls=style, lw=2, label=f"m={m}"
            )

        axes[ax_idx].set_xlabel("R [kpc]")
        axes[ax_idx].set_ylabel("Amplitude [kpc]")
        axes[ax_idx].set_title(f"t = {t:.0f} Myr")
        axes[ax_idx].legend()
        axes[ax_idx].set_xlim(0, 15)

    plt.suptitle(f"Fourier Bending Modes ({label})", y=1.02)
    plt.tight_layout()
    plt.savefig(os.path.join(fig_dir, f"fig4_bending_modes_{label}.png"))
    plt.close()
    print(f"  Saved fig4_bending_modes_{label}.png")


# ============================================================
# Figure 5: Angular momentum tilt profile
# ============================================================


def plot_angular_momentum_tilt(sim_dir, fig_dir, label="tidal"):
    """Plot tilt of angular momentum vector as function of R."""
    setup_style()

    snapshots = load_all_snapshots(sim_dir)

    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    times = np.array([s["t"] for s in snapshots])
    indices = np.linspace(0, len(times) - 1, 6, dtype=int)
    cmap = plt.cm.plasma
    colors = [cmap(i / len(indices)) for i in range(len(indices))]

    for k, idx in enumerate(indices):
        snap = snapshots[idx]
        Lp = compute_angular_momentum_profile(snap["pos"], snap["vel"], snap["mass_pp"])
        t = snap["t"]
        c = colors[k]

        axes[0].plot(Lp["R_mid"], Lp["L_tilt_x"], color=c, label=f"t={t:.0f}")
        axes[1].plot(Lp["R_mid"], Lp["L_tilt_y"], color=c, label=f"t={t:.0f}")

    axes[0].set_xlabel("R [kpc]")
    axes[0].set_ylabel("L_x tilt [deg]")
    axes[0].set_title("Angular Momentum Tilt (x-component)")
    axes[0].legend(fontsize=8)
    axes[0].axhline(0, color="gray", ls="--", lw=0.5)
    axes[0].set_xlim(0, 15)

    axes[1].set_xlabel("R [kpc]")
    axes[1].set_ylabel("L_y tilt [deg]")
    axes[1].set_title("Angular Momentum Tilt (y-component)")
    axes[1].legend(fontsize=8)
    axes[1].axhline(0, color="gray", ls="--", lw=0.5)
    axes[1].set_xlim(0, 15)

    plt.suptitle(f"Disk Tilt Profile ({label})", y=1.02)
    plt.tight_layout()
    plt.savefig(os.path.join(fig_dir, f"fig5_angular_momentum_{label}.png"))
    plt.close()
    print(f"  Saved fig5_angular_momentum_{label}.png")


# ============================================================
# Figure 6: 3D visualization of warped disk
# ============================================================


def plot_3d_warp(sim_dir, config, fig_dir):
    """3D scatter plot showing the warped disk morphology."""
    setup_style()

    snapshots = load_all_snapshots(sim_dir)
    snap = snapshots[-1]  # final snapshot
    pos = snap["pos"]
    t = snap["t"]

    x, y, z = pos[:, 0], pos[:, 1], pos[:, 2]
    R = np.sqrt(x**2 + y**2)

    # Subsample for clarity
    rng = np.random.default_rng(42)
    idx = rng.choice(len(pos), min(10000, len(pos)), replace=False)

    fig = plt.figure(figsize=(12, 8))
    ax = fig.add_subplot(111, projection="3d")

    scatter = ax.scatter(
        x[idx],
        y[idx],
        z[idx],
        c=z[idx],
        cmap="coolwarm",
        s=1,
        alpha=0.5,
        vmin=-2,
        vmax=2,
    )

    # Plot satellite
    r_sat = satellite_position(t, config)
    ax.scatter(
        *r_sat,
        c="cyan",
        s=200,
        marker="*",
        edgecolors="black",
        linewidths=1,
        zorder=10,
        label="Satellite",
    )

    ax.set_xlabel("x [kpc]")
    ax.set_ylabel("y [kpc]")
    ax.set_zlabel("z [kpc]")
    ax.set_title(f"Warped Stellar Disk at t = {t:.0f} Myr")
    ax.set_zlim(-5, 5)
    plt.colorbar(scatter, label="z [kpc]", shrink=0.6)
    ax.legend()

    plt.savefig(os.path.join(fig_dir, "fig6_3d_warp.png"))
    plt.close()
    print("  Saved fig6_3d_warp.png")


# ============================================================
# Figure 7: Parameter study summary
# ============================================================


def plot_parameter_study(param_dirs, param_values, param_name, fig_dir):
    """
    Plot warp amplitude vs a varied parameter.
    param_dirs: list of simulation directories
    param_values: corresponding parameter values
    param_name: name of varied parameter
    """
    setup_style()

    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    max_amps = []
    final_amps = []

    for sim_dir in param_dirs:
        if not os.path.exists(sim_dir):
            max_amps.append(np.nan)
            final_amps.append(np.nan)
            continue

        summary = warp_evolution_summary(sim_dir)
        max_amps.append(np.max(summary["max_warp_amplitude"]))
        final_amps.append(summary["max_warp_amplitude"][-1])

    axes[0].plot(param_values, max_amps, "o-", lw=2, ms=8)
    axes[0].set_xlabel(param_name)
    axes[0].set_ylabel("Peak warp amplitude [kpc]")
    axes[0].set_title("Maximum Warp Amplitude")

    axes[1].plot(param_values, final_amps, "s-", lw=2, ms=8, color="C1")
    axes[1].set_xlabel(param_name)
    axes[1].set_ylabel("Final warp amplitude [kpc]")
    axes[1].set_title("Final Warp Amplitude")

    plt.suptitle(f"Parameter Study: {param_name}", y=1.02)
    plt.tight_layout()
    plt.savefig(os.path.join(fig_dir, f"fig7_param_study_{param_name}.png"))
    plt.close()
    print(f"  Saved fig7_param_study_{param_name}.png")


# ============================================================
# Main: generate all figures
# ============================================================


def generate_all_figures(output_base, fig_dir, config=None):
    """Generate all analysis figures."""
    if config is None:
        config = DEFAULT_CONFIG.copy()

    os.makedirs(fig_dir, exist_ok=True)

    tidal_dir = os.path.join(output_base, "tidal")
    baseline_dir = os.path.join(output_base, "baseline")

    print("\nGenerating figures...")

    if os.path.exists(tidal_dir):
        plot_disk_views(tidal_dir, config, fig_dir)
        plot_warp_profiles(tidal_dir, fig_dir, "tidal")
        plot_bending_modes(tidal_dir, fig_dir, "tidal")
        plot_angular_momentum_tilt(tidal_dir, fig_dir, "tidal")
        plot_3d_warp(tidal_dir, config, fig_dir)

    if os.path.exists(baseline_dir):
        plot_warp_profiles(baseline_dir, fig_dir, "baseline")
        plot_bending_modes(baseline_dir, fig_dir, "baseline")

    if os.path.exists(tidal_dir) and os.path.exists(baseline_dir):
        plot_warp_time_evolution(tidal_dir, baseline_dir, fig_dir)

    print("\nAll figures generated!")


if __name__ == "__main__":
    output_base = os.path.join(os.path.dirname(__file__), "..", "output")
    fig_dir = os.path.join(os.path.dirname(__file__), "..", "figures")
    generate_all_figures(output_base, fig_dir)
