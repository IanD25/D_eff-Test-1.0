#!/usr/bin/env python3
"""
DS Phase 3D-1: 3D Ising Universality Test
=====================================================================
Kill test: Does the SV degeneracy swap occur in the 3D Ising model
(a different universality class from 2D)?

3D Ising: ν ≈ 0.630, γ ≈ 1.237, η ≈ 0.036
6 neighbors → 6×6 FIM → 6 singular values
T_c ≈ 4.5115 (Ferrenberg & Landau)

If the degeneracy swap transfers from 2D to 3D, the phenomenon generalizes
across universality classes. If not, the generalization hypothesis is dead.

Adapted from ising_fisher_phase_transition.py (Phase 2 v2, Manhattan distance).
"""

import os
import sys
import json
import time
import datetime
from collections import deque, Counter

import numpy as np
from scipy.optimize import curve_fit
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

OUTPUT_DIR = "phase3d1_results"
RAW_DIR = os.path.join(OUTPUT_DIR, "raw_data")
os.makedirs(RAW_DIR, exist_ok=True)

T_C_3D = 4.5115   # 3D Ising critical temperature (Ferrenberg & Landau)
T_C_2D = 2.0 / np.log(1 + np.sqrt(2))  # For reference comparisons

np.random.seed(42)

# Hardcoded Phase 2 (2D Ising, N=128, v2 Manhattan) SV profiles for comparison
PHASE2_SV_PROFILES = {
    1.50: [1.000, 0.225, 0.225, 0.040],
    1.02: [1.000, 1.000, 0.344, 0.090],
    1.00: [1.000, 1.000, 0.542, 0.079],
    0.90: [1.000, 0.479, 0.479, 0.223],
}


# ============================================================
# MODULE 1: 3D Ising Monte Carlo (Wolff cluster algorithm)
# ============================================================

def initialize_lattice_3d(L, mode='random'):
    """Returns L×L×L array of spins in {+1, -1}."""
    if mode == 'cold':
        return np.ones((L, L, L), dtype=np.int8)
    elif mode == 'hot' or mode == 'random':
        return np.where(np.random.random((L, L, L)) < 0.5, 1, -1).astype(np.int8)
    else:
        raise ValueError(f"Unknown mode: {mode}")


def wolff_step_3d(lattice, T, J=1.0):
    """
    Execute one Wolff cluster flip on 3D periodic lattice.
    6 neighbors: ±x, ±y, ±z with periodic boundary conditions.
    """
    L = lattice.shape[0]
    p_add = 1.0 - np.exp(-2.0 * J / T)

    # Choose random seed
    si = np.random.randint(L)
    sj = np.random.randint(L)
    sk = np.random.randint(L)
    seed_spin = lattice[si, sj, sk]

    # BFS cluster growth
    visited = np.zeros((L, L, L), dtype=np.bool_)
    visited[si, sj, sk] = True
    queue = deque()
    queue.append((si, sj, sk))
    cluster = [(si, sj, sk)]

    while queue:
        ci, cj, ck = queue.popleft()
        # Check 6 neighbors with periodic BC
        for di, dj, dk in [(-1,0,0), (1,0,0), (0,-1,0), (0,1,0), (0,0,-1), (0,0,1)]:
            ni = (ci + di) % L
            nj = (cj + dj) % L
            nk = (ck + dk) % L
            if not visited[ni, nj, nk] and lattice[ni, nj, nk] == seed_spin:
                if np.random.random() < p_add:
                    visited[ni, nj, nk] = True
                    queue.append((ni, nj, nk))
                    cluster.append((ni, nj, nk))

    # Flip cluster
    for ci, cj, ck in cluster:
        lattice[ci, cj, ck] = -lattice[ci, cj, ck]

    return lattice


def wolff_sweep_3d(lattice, T, n_steps):
    """Execute n_steps Wolff cluster flips."""
    for _ in range(n_steps):
        wolff_step_3d(lattice, T)
    return lattice


def equilibrate_3d(lattice, T, n_sweeps=200):
    """Equilibrate 3D lattice at temperature T using Wolff steps."""
    for _ in range(n_sweeps):
        wolff_step_3d(lattice, T)
    return lattice


# ============================================================
# MODULE 2: Correlation Function Measurement (3D)
# ============================================================

def measure_magnetization_3d(lattice):
    """Mean absolute magnetization: |sum(s)| / L^3"""
    return np.abs(np.mean(lattice))


def measure_energy_3d(lattice, J=1.0):
    """
    Energy per spin: -J * sum(s_i * s_j) / L^3 over NN pairs.
    3D: 3 bonds per site (right, down, forward).
    """
    L = lattice.shape[0]
    L3 = L * L * L
    lat = lattice.astype(np.float64)
    e = -J * np.sum(lat * np.roll(lat, -1, axis=0))   # x-direction
    e += -J * np.sum(lat * np.roll(lat, -1, axis=1))   # y-direction
    e += -J * np.sum(lat * np.roll(lat, -1, axis=2))   # z-direction
    return e / L3


def _build_distance_grid_3d(L):
    """
    Pre-compute 3D Manhattan distance grid for periodic lattice.
    Returns array of shape (L, L, L) where dist[dx, dy, dz] = Manhattan distance.
    """
    ix = np.arange(L)
    dx = np.minimum(ix, L - ix)
    # Build 3D distance grid via broadcasting
    dist = dx[:, None, None] + dx[None, :, None] + dx[None, None, :]
    return dist


def measure_correlations_3d(lattice, max_r=None):
    """
    Measure connected G(r) using FFT-based autocorrelation on 3D periodic lattice.
    Returns G[r] for r = 0, 1, ..., max_r using Manhattan distance.
    """
    L = lattice.shape[0]
    if max_r is None:
        max_r = L // 2

    lat_float = lattice.astype(np.float64)
    m = np.mean(lat_float)
    m_sq = m * m

    # FFT-based autocorrelation (3D)
    F = np.fft.fftn(lat_float)
    power = F * np.conj(F)
    G_full = np.real(np.fft.ifftn(power)) / (L * L * L)

    # Radial average using Manhattan distance
    dist_grid = _build_distance_grid_3d(L)

    G_r = np.zeros(max_r + 1)
    counts = np.zeros(max_r + 1)

    for r in range(max_r + 1):
        mask = (dist_grid == r)
        if np.any(mask):
            G_r[r] = np.mean(G_full[mask]) - m_sq
            counts[r] = np.sum(mask)

    return G_r


def accumulate_correlations_3d(lattice, T, n_configs=500, n_sweeps_between=5):
    """
    Measure G(r) averaged over n_configs independent configurations.
    Returns averaged G(r) array.
    """
    L = lattice.shape[0]
    max_r = L // 2
    G_sum = np.zeros(max_r + 1)

    for i in range(n_configs):
        for _ in range(n_sweeps_between):
            wolff_step_3d(lattice, T)
        G_r = measure_correlations_3d(lattice, max_r=max_r)
        G_sum += G_r[:max_r + 1]

    return G_sum / n_configs


def estimate_correlation_length(G_r):
    """
    Fit G(r) ~ A * exp(-r/xi) for r in [1, len(G_r)//2].
    Returns (xi, R^2).
    """
    max_fit = len(G_r) // 2
    r_vals = np.arange(1, max_fit + 1)
    g_vals = G_r[1:max_fit + 1]

    # Filter positive values for log fit
    pos = g_vals > 0
    if np.sum(pos) < 3:
        return float('inf'), 0.0

    r_pos = r_vals[pos]
    log_g = np.log(g_vals[pos])

    try:
        def exp_decay(r, log_A, inv_xi):
            return log_A - r * inv_xi

        popt, _ = curve_fit(exp_decay, r_pos, log_g, p0=[0.0, 0.1], maxfev=5000)
        inv_xi = popt[1]
        if inv_xi <= 0:
            return float('inf'), 0.0
        xi = 1.0 / inv_xi

        # R^2
        predicted = exp_decay(r_pos, *popt)
        ss_res = np.sum((log_g - predicted) ** 2)
        ss_tot = np.sum((log_g - np.mean(log_g)) ** 2)
        r_sq = 1.0 - ss_res / ss_tot if ss_tot > 0 else 0.0

        return xi, r_sq
    except Exception:
        return float('inf'), 0.0


# ============================================================
# MODULE 3: Thermal Fisher Kernel and FIM Construction (3D)
# ============================================================

def build_thermal_kernel_3d(G_r, L, v0):
    """
    Build thermal probability distribution p_{v0}(u; T) using |G(r)|.
    v0 = (i0, j0, k0) tuple on 3D periodic lattice.
    Returns flat array of shape (L^3,).
    """
    i0, j0, k0 = v0
    max_r = len(G_r) - 1

    # Build distance from v0 to all sites (3D Manhattan, periodic)
    ix = np.arange(L)
    di = np.minimum(np.abs(ix - i0), L - np.abs(ix - i0))
    dj = np.minimum(np.abs(ix - j0), L - np.abs(ix - j0))
    dk = np.minimum(np.abs(ix - k0), L - np.abs(ix - k0))

    # 3D distance grid via broadcasting
    dist = di[:, None, None] + dj[None, :, None] + dk[None, None, :]

    # Build weights from |G(r)|
    weights = np.zeros((L, L, L), dtype=np.float64)
    for r in range(max_r + 1):
        mask = (dist == r)
        weights[mask] = np.abs(G_r[r])

    # For distances beyond max_r, use last available value
    beyond = dist > max_r
    if np.any(beyond):
        weights[beyond] = np.abs(G_r[max_r])

    # Flatten and normalize
    p = weights.flatten()
    total = np.sum(p)
    if total < 1e-30:
        p = np.ones(L ** 3) / (L ** 3)
    else:
        p = p / total

    return p


def compute_FIM_thermal_3d(G_r, L, v0):
    """
    Compute 6×6 Fisher Information Matrix at vertex v0 using thermal kernel.
    v0 = (i0, j0, k0) tuple. Returns 6×6 FIM matrix.
    """
    i0, j0, k0 = v0

    # p_{v0}
    p_v0 = build_thermal_kernel_3d(G_r, L, v0)
    log_p_v0 = np.log(p_v0 + 1e-30)

    # 6 neighbors (periodic BC): ±x, ±y, ±z
    neighbors = [
        ((i0 - 1) % L, j0, k0),           # -x
        ((i0 + 1) % L, j0, k0),           # +x
        (i0, (j0 - 1) % L, k0),           # -y
        (i0, (j0 + 1) % L, k0),           # +y
        (i0, j0, (k0 - 1) % L),           # -z
        (i0, j0, (k0 + 1) % L),           # +z
    ]
    k = len(neighbors)

    # Score vectors
    score_vectors = np.zeros((k, L ** 3))
    for j, w in enumerate(neighbors):
        p_w = build_thermal_kernel_3d(G_r, L, w)
        log_p_w = np.log(p_w + 1e-30)
        score_vectors[j, :] = log_p_w - log_p_v0

    # FIM: F_ij = sum_u s_i(u) * s_j(u) * p_{v0}(u)
    weighted_scores = score_vectors * np.sqrt(p_v0)[np.newaxis, :]
    FIM = weighted_scores @ weighted_scores.T

    return FIM


def gap_based_rank(sv_norm):
    """Gap-based rank from normalized singular values."""
    if len(sv_norm) <= 1:
        return 1
    ratios = sv_norm[1:] / np.maximum(sv_norm[:-1], 1e-15)
    return int(np.argmin(ratios)) + 1


def participation_ratio(svs):
    """PR = (sum sv_i)^2 / sum(sv_i^2)"""
    s = np.sum(svs)
    s2 = np.sum(svs ** 2)
    if s2 < 1e-30:
        return 0.0
    return s * s / s2


def disorder_index(sv_norm, rank):
    """eta = sv[rank] / sv[rank-1] — ratio of (d+1)-th to d-th SV."""
    if rank >= len(sv_norm) or rank < 1:
        return 0.0
    return sv_norm[rank] / sv_norm[rank - 1] if sv_norm[rank - 1] > 1e-15 else 0.0


def sv_profile(FIM):
    """Returns normalized SV profile. Shape: (6,) for 3D."""
    svs = np.linalg.svd(FIM, compute_uv=False)
    if svs[0] < 1e-30:
        return svs
    return svs / svs[0]


def fisher_diagnostics_3d(G_r, L, n_samples=30):
    """
    Sample n_samples random vertices on the L^3 lattice.
    For each: compute FIM, rank, PR, eta, sv_profile.
    Returns diagnostics dict.
    """
    ranks = []
    prs = []
    etas = []
    sv_profiles = []
    gap_ratios = []

    # Sample random interior vertices
    margin = max(2, L // 8)
    sampled = 0
    attempts = 0
    while sampled < n_samples and attempts < n_samples * 5:
        i = np.random.randint(margin, L - margin)
        j = np.random.randint(margin, L - margin)
        k = np.random.randint(margin, L - margin)
        attempts += 1

        FIM = compute_FIM_thermal_3d(G_r, L, (i, j, k))
        svs = np.linalg.svd(FIM, compute_uv=False)

        if svs[0] < 1e-30:
            continue

        sv_n = svs / svs[0]
        r = gap_based_rank(sv_n)
        pr = participation_ratio(svs)
        eta = disorder_index(sv_n, r)
        gr = svs[0] / svs[1] if len(svs) > 1 and svs[1] > 1e-30 else float('inf')

        ranks.append(r)
        prs.append(pr)
        etas.append(eta)
        sv_profiles.append(sv_n)
        gap_ratios.append(gr)
        sampled += 1

    if not ranks:
        return None

    sv_arr = np.array(sv_profiles)
    gr_arr = np.array([g for g in gap_ratios if g < 1e10])

    return {
        'rank_mean': float(np.mean(ranks)),
        'rank_std': float(np.std(ranks)),
        'rank_distribution': dict(Counter(ranks)),
        'pr_mean': float(np.mean(prs)),
        'pr_std': float(np.std(prs)),
        'eta_mean': float(np.mean(etas)),
        'eta_std': float(np.std(etas)),
        'sv_profile_mean': sv_arr.mean(axis=0).tolist(),
        'sv_profile_std': sv_arr.std(axis=0).tolist(),
        'gap_ratio_mean': float(np.mean(gr_arr)) if len(gr_arr) > 0 else float('inf'),
    }


# ============================================================
# MODULE 4: Temperature Sweep and Diagnostics
# ============================================================

def measure_macroscopic_3d(lattice, T, n_configs=500, n_sweeps_between=5):
    """
    Measure macroscopic observables over n_configs decorrelated configurations.
    3D: L^3 sites, 3 bonds per site for energy.
    """
    L = lattice.shape[0]
    L3 = L * L * L

    ms = []
    m2s = []
    m4s = []
    Es = []
    E2s = []

    for _ in range(n_configs):
        for _ in range(n_sweeps_between):
            wolff_step_3d(lattice, T)
        m_val = np.mean(lattice.astype(np.float64))
        ms.append(np.abs(m_val))
        m2s.append(m_val ** 2)
        m4s.append(m_val ** 4)

        e_val = measure_energy_3d(lattice)
        Es.append(e_val)
        E2s.append(e_val ** 2)

    m_mean = np.mean(ms)
    m2_mean = np.mean(m2s)
    m4_mean = np.mean(m4s)
    E_mean = np.mean(Es)
    E2_mean = np.mean(E2s)

    chi = L3 * (m2_mean - m_mean ** 2)
    C = L3 * (E2_mean - E_mean ** 2) / (T * T)
    U_binder = 1.0 - m4_mean / (3.0 * m2_mean ** 2) if m2_mean > 1e-15 else 0.0

    return {
        'magnetization': float(m_mean),
        'susceptibility': float(chi),
        'energy': float(E_mean),
        'specific_heat': float(C),
        'binder': float(U_binder),
    }


def temperature_sweep_3d(L, temperatures, n_eq_sweeps=200, n_configs=500):
    """
    Full temperature sweep for lattice size L (3D).
    Returns list of result dicts, one per temperature.
    """
    results = []

    for t_idx, T in enumerate(temperatures):
        t0 = time.time()
        T_over_Tc = T / T_C_3D

        # Initialize: hot start above Tc, cold start below
        mode = 'hot' if T > T_C_3D else 'cold'
        lattice = initialize_lattice_3d(L, mode=mode)

        # Equilibrate
        equilibrate_3d(lattice, T, n_sweeps=n_eq_sweeps)

        # Near Tc: use more configs for better statistics
        nc = n_configs
        ns = 30
        if 0.98 <= T_over_Tc <= 1.05:
            nc = min(n_configs * 2, 1000)
            ns = 50

        # Accumulate correlations
        lat_corr = lattice.copy()
        G_r = accumulate_correlations_3d(lat_corr, T, n_configs=nc, n_sweeps_between=5)

        # Correlation length
        xi, xi_r2 = estimate_correlation_length(G_r)

        # Fisher diagnostics
        fisher = fisher_diagnostics_3d(G_r, L, n_samples=ns)

        # Macroscopic observables (fresh pass from equilibrated state)
        lat_macro = lattice.copy()
        macro = measure_macroscopic_3d(lat_macro, T, n_configs=nc, n_sweeps_between=5)

        dt = time.time() - t0

        result = {
            'T': float(T),
            'T_over_Tc': float(T_over_Tc),
            'L': L,
            'G_r': G_r.tolist(),
            'xi': float(xi) if np.isfinite(xi) else -1.0,
            'xi_fit_R2': float(xi_r2),
            'time_sec': float(dt),
        }
        if fisher is not None:
            result.update({f'fisher_{k}': v for k, v in fisher.items()})
        result.update({f'macro_{k}': v for k, v in macro.items()})

        rank_str = f"{fisher['rank_mean']:.2f}" if fisher else "N/A"
        eta_str = f"{fisher['eta_mean']:.3f}" if fisher else "N/A"
        pr_str = f"{fisher['pr_mean']:.3f}" if fisher else "N/A"
        chi_str = f"{macro['susceptibility']:.1f}"

        print(f"  [{t_idx+1}/{len(temperatures)}] T={T:.4f} (T/Tc={T_over_Tc:.3f}): "
              f"rank={rank_str} eta={eta_str} PR={pr_str} chi={chi_str}  [{dt:.1f}s]")

        results.append(result)

    return results


# ============================================================
# MODULE 5: Plots and Results Report
# ============================================================

def plot_diagnostic_panel(results, L):
    """Figure 1: 4-panel diagnostic plot for 3D Ising."""
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    fig.suptitle(f"Fisher Diagnostic Panel — 3D Ising L={L}", fontsize=14, fontweight='bold')

    T_ratios = [r['T_over_Tc'] for r in results]

    # (a) D_eff (rank) vs T/Tc
    ax = axes[0, 0]
    ranks = [r.get('fisher_rank_mean', np.nan) for r in results]
    rank_stds = [r.get('fisher_rank_std', 0) for r in results]
    ax.errorbar(T_ratios, ranks, yerr=rank_stds, fmt='o-', color='tab:blue', capsize=3)
    ax.axhline(4, color='green', ls='--', alpha=0.5, label='D_eff=4 (predicted 3D)')
    ax.axhline(3, color='gray', ls='--', alpha=0.5, label='D_eff=3 (2D value)')
    ax.axvline(1.0, color='red', ls='--', alpha=0.5, label='T_c')
    ax.set_xlabel('T / T_c')
    ax.set_ylabel('D_eff (gap-based rank)')
    ax.set_title('(a) Effective Dimension')
    ax.legend(fontsize=8)

    # (b) eta vs T/Tc
    ax = axes[0, 1]
    etas = [r.get('fisher_eta_mean', np.nan) for r in results]
    eta_stds = [r.get('fisher_eta_std', 0) for r in results]
    ax.errorbar(T_ratios, etas, yerr=eta_stds, fmt='s-', color='tab:orange', capsize=3)
    ax.axhline(0.23, color='blue', ls='--', alpha=0.4, label='Torus η=0.23')
    ax.axhline(0.68, color='green', ls='--', alpha=0.4, label='RGG η=0.68')
    ax.axhline(0.93, color='red', ls='--', alpha=0.4, label='ER η=0.93')
    ax.axvline(1.0, color='red', ls='--', alpha=0.5)
    ax.set_xlabel('T / T_c')
    ax.set_ylabel('Disorder Index η')
    ax.set_title('(b) Disorder Index')
    ax.legend(fontsize=8)

    # (c) Gap ratio vs T/Tc
    ax = axes[1, 0]
    grs = [r.get('fisher_gap_ratio_mean', np.nan) for r in results]
    ax.plot(T_ratios, grs, 'D-', color='tab:green')
    ax.axhline(1.0, color='gray', ls='--', alpha=0.5, label='No gap')
    ax.axvline(1.0, color='red', ls='--', alpha=0.5)
    ax.set_xlabel('T / T_c')
    ax.set_ylabel('Gap Ratio (sv₁/sv₂)')
    ax.set_title('(c) Gap Ratio')
    ax.set_yscale('log')
    ax.legend(fontsize=8)

    # (d) Susceptibility vs T/Tc
    ax = axes[1, 1]
    chis = [r.get('macro_susceptibility', np.nan) for r in results]
    ax.plot(T_ratios, chis, '^-', color='tab:red')
    ax.axvline(1.0, color='red', ls='--', alpha=0.5, label='T_c')
    ax.set_xlabel('T / T_c')
    ax.set_ylabel('Susceptibility χ')
    ax.set_title('(d) Susceptibility')
    ax.legend(fontsize=8)

    plt.tight_layout()
    fname = os.path.join(OUTPUT_DIR, f"diagnostic_panel_L{L}.png")
    plt.savefig(fname, dpi=150)
    plt.close()
    print(f"  Saved {fname}")


def plot_sv_evolution(results, L):
    """Figure 2: SV profile evolution at representative temperatures (6 bars)."""
    target_ratios = [1.50, 1.10, 1.02, 1.00, 0.90]
    fig, axes = plt.subplots(2, 3, figsize=(18, 10))
    fig.suptitle(f"SV Profile Evolution — 3D Ising L={L} (6 SVs)", fontsize=14, fontweight='bold')
    axes_flat = axes.flatten()

    for idx, t_target in enumerate(target_ratios):
        ax = axes_flat[idx]
        best = min(results, key=lambda r: abs(r['T_over_Tc'] - t_target))
        sv_mean = np.array(best.get('fisher_sv_profile_mean', [0]*6))
        sv_std = np.array(best.get('fisher_sv_profile_std', [0]*6))
        rank = best.get('fisher_rank_mean', 0)
        eta = best.get('fisher_eta_mean', 0)

        # Color by temperature regime
        if t_target >= 1.10:
            colors = ['tab:blue'] * len(sv_mean)
        elif t_target >= 1.00:
            colors = ['tab:purple'] * len(sv_mean)
        else:
            colors = ['tab:red'] * len(sv_mean)

        ax.bar(range(1, len(sv_mean) + 1), sv_mean, yerr=sv_std, capsize=4,
               color=colors, alpha=0.7, edgecolor='black')
        ax.set_xlabel('SV Index')
        ax.set_ylabel('Normalized SV')
        ax.set_title(f"T/Tc={best['T_over_Tc']:.3f}  rank={rank:.1f}  η={eta:.2f}")
        ax.set_ylim(0, 1.15)
        ax.set_xticks(range(1, len(sv_mean) + 1))

    # Last panel: mean SV profile at T_c with extra detail
    ax = axes_flat[5]
    tc_result = min(results, key=lambda r: abs(r['T_over_Tc'] - 1.0))
    sv_mean = np.array(tc_result.get('fisher_sv_profile_mean', [0]*6))
    sv_std = np.array(tc_result.get('fisher_sv_profile_std', [0]*6))
    ax.bar(range(1, len(sv_mean) + 1), sv_mean, yerr=sv_std, capsize=4,
           color='tab:purple', alpha=0.7, edgecolor='black')
    ax.set_xlabel('SV Index')
    ax.set_ylabel('Normalized SV')
    ax.set_title(f"T=T_c (detail): mean ± std")
    ax.set_ylim(0, 1.15)
    ax.set_xticks(range(1, len(sv_mean) + 1))

    # Add SV values as text annotations on this panel
    for i, (m, s) in enumerate(zip(sv_mean, sv_std)):
        ax.text(i + 1, m + s + 0.03, f"{m:.3f}", ha='center', va='bottom', fontsize=8)

    plt.tight_layout()
    fname = os.path.join(OUTPUT_DIR, f"sv_profile_evolution_L{L}.png")
    plt.savefig(fname, dpi=150)
    plt.close()
    print(f"  Saved {fname}")


def plot_sv_comparison_2d_vs_3d(results_3d, L):
    """
    Figure 3: Side-by-side 2D vs 3D SV profile comparison at matched T/T_c.
    Left: 2D Ising (hardcoded from Phase 2 N=128 v2).
    Right: 3D Ising from this experiment.
    """
    comparison_ratios = [1.50, 1.02, 1.00, 0.90]
    fig, axes = plt.subplots(len(comparison_ratios), 2, figsize=(14, 16))
    fig.suptitle("SV Profile Comparison: 2D Ising vs 3D Ising", fontsize=14, fontweight='bold')

    for row, t_target in enumerate(comparison_ratios):
        # Left: 2D
        ax_2d = axes[row, 0]
        sv_2d = PHASE2_SV_PROFILES[t_target]
        ax_2d.bar(range(1, 5), sv_2d, color='tab:blue', alpha=0.7, edgecolor='black')
        ax_2d.set_ylabel('Normalized SV')
        ax_2d.set_title(f"2D Ising T/Tc={t_target:.2f}")
        ax_2d.set_ylim(0, 1.15)
        ax_2d.set_xticks(range(1, 5))
        ax_2d.set_xlabel('SV Index')
        # Annotate degeneracy
        if abs(sv_2d[1] - sv_2d[2]) < 0.05:
            ax_2d.annotate('SV2≈SV3', xy=(2.5, max(sv_2d[1], sv_2d[2]) + 0.05),
                          ha='center', fontsize=8, color='blue')
        if abs(sv_2d[0] - sv_2d[1]) < 0.05:
            ax_2d.annotate('SV1≈SV2', xy=(1.5, sv_2d[0] + 0.05),
                          ha='center', fontsize=8, color='red')

        # Right: 3D
        ax_3d = axes[row, 1]
        best = min(results_3d, key=lambda r: abs(r['T_over_Tc'] - t_target))
        sv_3d = np.array(best.get('fisher_sv_profile_mean', [0]*6))
        sv_3d_std = np.array(best.get('fisher_sv_profile_std', [0]*6))
        ax_3d.bar(range(1, 7), sv_3d, yerr=sv_3d_std, capsize=3,
                 color='tab:orange', alpha=0.7, edgecolor='black')
        ax_3d.set_ylabel('Normalized SV')
        ax_3d.set_title(f"3D Ising T/Tc={best['T_over_Tc']:.3f}")
        ax_3d.set_ylim(0, 1.15)
        ax_3d.set_xticks(range(1, 7))
        ax_3d.set_xlabel('SV Index')
        # Annotate degeneracies in 3D
        # Check for triple degeneracy SV2=SV3=SV4
        if len(sv_3d) >= 4:
            if abs(sv_3d[1] - sv_3d[2]) < 0.05 and abs(sv_3d[2] - sv_3d[3]) < 0.05:
                ax_3d.annotate('SV2≈SV3≈SV4', xy=(3, max(sv_3d[1:4]) + 0.05),
                              ha='center', fontsize=8, color='blue')
            # Check if SV1≈SV2 (swap)
            if abs(sv_3d[0] - sv_3d[1]) < 0.05:
                ax_3d.annotate('SV1≈SV2', xy=(1.5, sv_3d[0] + 0.05),
                              ha='center', fontsize=8, color='red')

    plt.tight_layout()
    fname = os.path.join(OUTPUT_DIR, "sv_comparison_2d_vs_3d.png")
    plt.savefig(fname, dpi=150)
    plt.close()
    print(f"  Saved {fname}")


def plot_leading_indicator(results, L):
    """Figure 4: Leading indicator comparison — eta vs chi."""
    T_ratios = [r['T_over_Tc'] for r in results]
    etas = [r.get('fisher_eta_mean', np.nan) for r in results]
    chis = [r.get('macro_susceptibility', np.nan) for r in results]

    # eta onset: exceeds mean(eta_high_T) + 2*std(eta_high_T)
    high_T_etas = [e for e, t in zip(etas, T_ratios) if t >= 1.3 and not np.isnan(e)]
    if high_T_etas:
        eta_threshold = np.mean(high_T_etas) + 2 * np.std(high_T_etas)
    else:
        eta_threshold = 0.35

    # chi onset: exceeds 10% of peak
    valid_chis = [c for c in chis if not np.isnan(c)]
    chi_threshold = 0.1 * max(valid_chis) if valid_chis else 0

    T_fisher = None
    T_chi = None

    # Scan from high T downward
    for t, e in sorted(zip(T_ratios, etas), reverse=True):
        if not np.isnan(e) and e > eta_threshold and T_fisher is None:
            T_fisher = t

    for t, c in sorted(zip(T_ratios, chis), reverse=True):
        if not np.isnan(c) and c > chi_threshold and T_chi is None:
            T_chi = c  # BUG FIX: should be t not c
            T_chi = t

    fig, ax1 = plt.subplots(figsize=(10, 6))
    ax2 = ax1.twinx()

    ax1.plot(T_ratios, etas, 's-', color='tab:blue', label='Fisher η', markersize=6)
    ax2.plot(T_ratios, chis, '^-', color='tab:red', label='Susceptibility χ', markersize=6)

    ax1.axhline(eta_threshold, color='tab:blue', ls=':', alpha=0.4)
    ax2.axhline(chi_threshold, color='tab:red', ls=':', alpha=0.4)

    if T_fisher is not None:
        ax1.axvline(T_fisher, color='tab:blue', ls='--', alpha=0.6, label=f'η onset: T/Tc={T_fisher:.3f}')
    if T_chi is not None:
        ax1.axvline(T_chi, color='tab:red', ls='--', alpha=0.6, label=f'χ onset: T/Tc={T_chi:.3f}')

    ax1.axvline(1.0, color='black', ls='--', alpha=0.3)

    ax1.set_xlabel('T / T_c')
    ax1.set_ylabel('Disorder Index η', color='tab:blue')
    ax2.set_ylabel('Susceptibility χ', color='tab:red')

    delta = (T_fisher - T_chi) if T_fisher is not None and T_chi is not None else 0
    if delta > 0.01:
        verdict = f"Fisher LEADING by ΔT/Tc = {delta:.3f}"
    elif delta < -0.01:
        verdict = f"Fisher LAGGING by ΔT/Tc = {-delta:.3f}"
    else:
        verdict = "Simultaneous onset"

    ax1.set_title(f"Leading Indicator Test — 3D Ising L={L}: {verdict}")

    lines1, labels1 = ax1.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    ax1.legend(lines1 + lines2, labels1 + labels2, loc='upper left', fontsize=8)

    plt.tight_layout()
    fname = os.path.join(OUTPUT_DIR, f"leading_indicator_L{L}.png")
    plt.savefig(fname, dpi=150)
    plt.close()
    print(f"  Saved {fname}")

    return T_fisher, T_chi, verdict


def plot_pr_rmax(results, L):
    """Figure 5: PR vs r_max at 3 temperatures."""
    target_ratios = [1.50, 1.02, 1.00]
    colors = ['tab:blue', 'tab:orange', 'tab:red']
    labels = ['T/Tc=1.50', 'T/Tc=1.02', 'T/Tc=1.00']

    fig, ax = plt.subplots(figsize=(10, 6))

    for t_target, color, label in zip(target_ratios, colors, labels):
        best = min(results, key=lambda r: abs(r['T_over_Tc'] - t_target))
        G_r_full = np.array(best['G_r'])

        # Compute PR at different r_max truncations
        r_max_values = list(range(3, len(G_r_full), 2))
        pr_values = []

        for r_max in r_max_values:
            G_r_trunc = G_r_full[:r_max + 1]
            prs = []
            for _ in range(10):
                margin = max(2, L // 8)
                i = np.random.randint(margin, L - margin)
                j = np.random.randint(margin, L - margin)
                k = np.random.randint(margin, L - margin)
                FIM = compute_FIM_thermal_3d(G_r_trunc, L, (i, j, k))
                svs = np.linalg.svd(FIM, compute_uv=False)
                if svs[0] > 1e-30:
                    prs.append(participation_ratio(svs))
            pr_values.append(np.mean(prs) if prs else 0)

        ax.plot(r_max_values, pr_values, 'o-', color=color, label=label, markersize=4)

    ax.axhline(3.0, color='gray', ls=':', alpha=0.3, label='d=3')
    ax.axhline(4.0, color='gray', ls='--', alpha=0.3, label='d+1=4')
    ax.set_xlabel('r_max (correlation truncation range)')
    ax.set_ylabel('Participation Ratio')
    ax.set_title(f'PR vs Correlation Range — 3D Ising L={L}')
    ax.legend()
    plt.tight_layout()
    fname = os.path.join(OUTPUT_DIR, f"pr_rmax_thermal_L{L}.png")
    plt.savefig(fname, dpi=150)
    plt.close()
    print(f"  Saved {fname}")


def plot_finite_size_scaling(all_results):
    """Figure 6: Finite-size scaling — all L overlaid."""
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))
    fig.suptitle("Finite-Size Scaling — 3D Ising", fontsize=14, fontweight='bold')

    colors = {16: 'tab:blue', 32: 'tab:orange', 48: 'tab:red'}

    for L_val, results in sorted(all_results.items()):
        T_ratios = [r['T_over_Tc'] for r in results]
        ranks = [r.get('fisher_rank_mean', np.nan) for r in results]
        etas = [r.get('fisher_eta_mean', np.nan) for r in results]

        ax1.plot(T_ratios, ranks, 'o-', color=colors.get(L_val, 'black'),
                 label=f'L={L_val}', markersize=5)
        ax2.plot(T_ratios, etas, 's-', color=colors.get(L_val, 'black'),
                 label=f'L={L_val}', markersize=5)

    ax1.axvline(1.0, color='black', ls='--', alpha=0.3)
    ax1.axhline(4, color='green', ls='--', alpha=0.3, label='D=4 predicted')
    ax1.axhline(3, color='gray', ls=':', alpha=0.3, label='D=3 (2D value)')
    ax1.set_xlabel('T / T_c')
    ax1.set_ylabel('D_eff (rank mean)')
    ax1.set_title('Effective Dimension')
    ax1.legend()

    ax2.axvline(1.0, color='black', ls='--', alpha=0.3)
    ax2.axhline(0.23, color='blue', ls='--', alpha=0.3)
    ax2.axhline(0.68, color='green', ls='--', alpha=0.3)
    ax2.set_xlabel('T / T_c')
    ax2.set_ylabel('Disorder Index η')
    ax2.set_title('Disorder Index')
    ax2.legend()

    plt.tight_layout()
    fname = os.path.join(OUTPUT_DIR, "finite_size_scaling.png")
    plt.savefig(fname, dpi=150)
    plt.close()
    print(f"  Saved {fname}")


def identify_degeneracy_pattern(sv_profile_mean):
    """
    Identify degeneracy pattern in 6-component SV profile.
    Returns a string describing the pattern.
    """
    sv = np.array(sv_profile_mean)
    if len(sv) < 6:
        return "INCOMPLETE"

    tol = 0.05  # tolerance for degeneracy detection

    patterns = []

    # Check SV1≈SV2
    if abs(sv[0] - sv[1]) < tol:
        patterns.append("SV1≈SV2")
    # Check SV2≈SV3
    if abs(sv[1] - sv[2]) < tol:
        patterns.append("SV2≈SV3")
    # Check SV3≈SV4
    if abs(sv[2] - sv[3]) < tol:
        patterns.append("SV3≈SV4")
    # Check SV4≈SV5
    if abs(sv[3] - sv[4]) < tol:
        patterns.append("SV4≈SV5")
    # Check SV5≈SV6
    if abs(sv[4] - sv[5]) < tol:
        patterns.append("SV5≈SV6")

    # Check triple degeneracies
    if abs(sv[1] - sv[2]) < tol and abs(sv[2] - sv[3]) < tol:
        patterns.append("TRIPLE:SV2≈SV3≈SV4")
    if abs(sv[0] - sv[1]) < tol and abs(sv[1] - sv[2]) < tol:
        patterns.append("TRIPLE:SV1≈SV2≈SV3")

    if not patterns:
        return "NO DEGENERACY"

    return "; ".join(patterns)


def write_results_report(all_results, verdicts, total_time):
    """Write the full results markdown report for Phase 3D-1."""
    lines = []
    lines.append("# Phase 3D-1 Results: 3D Ising Universality Test\n")
    lines.append(f"Generated: {datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")

    # Section 1: Parameters
    lines.append("## 1. System Parameters\n")
    lines.append(f"- T_c (3D) = {T_C_3D:.4f}")
    lines.append(f"- T_c (2D) = {T_C_2D:.6f} (reference)")
    lines.append(f"- Lattice sizes: {sorted(all_results.keys())}")
    lines.append(f"- Temperature points: 12 (T/Tc from 1.50 to 0.90)")
    lines.append(f"- FIM dimensionality: 6×6 (6 neighbors)")
    lines.append(f"- Distance metric: 3D Manhattan")
    lines.append(f"- Total runtime: {total_time:.1f}s\n")

    # Section 2: Correlation Function
    lines.append("## 2. Correlation Function Behavior\n")
    for L_val in sorted(all_results.keys()):
        lines.append(f"### L = {L_val}\n")
        lines.append("| T/Tc | ξ | R² (exp fit) | Notes |")
        lines.append("|------|---|-------------|-------|")
        for r in all_results[L_val]:
            xi_str = f"{r['xi']:.2f}" if r['xi'] > 0 else "∞"
            r2_str = f"{r['xi_fit_R2']:.3f}"
            note = ""
            if r['xi_fit_R2'] < 0.85:
                note = "Power law onset"
            if r['xi'] < 0:
                xi_str = "divergent"
                note = "Power law / divergent"
            lines.append(f"| {r['T_over_Tc']:.3f} | {xi_str} | {r2_str} | {note} |")
        lines.append("")

    # Section 3: Fisher Diagnostics Table
    lines.append("## 3. Fisher Diagnostics vs Temperature\n")
    for L_val in sorted(all_results.keys()):
        lines.append(f"### L = {L_val}\n")
        lines.append("| T/Tc | Rank (mean±std) | η (mean±std) | PR (mean) | Gap Ratio | χ | C |")
        lines.append("|------|----------------|-------------|-----------|-----------|---|---|")
        for r in all_results[L_val]:
            rank_m = r.get('fisher_rank_mean', 0)
            rank_s = r.get('fisher_rank_std', 0)
            eta_m = r.get('fisher_eta_mean', 0)
            eta_s = r.get('fisher_eta_std', 0)
            pr_m = r.get('fisher_pr_mean', 0)
            gr = r.get('fisher_gap_ratio_mean', 0)
            chi = r.get('macro_susceptibility', 0)
            C_val = r.get('macro_specific_heat', 0)
            lines.append(f"| {r['T_over_Tc']:.3f} | {rank_m:.2f}±{rank_s:.2f} | "
                        f"{eta_m:.3f}±{eta_s:.3f} | {pr_m:.3f} | {gr:.2f} | "
                        f"{chi:.1f} | {C_val:.2f} |")
        lines.append("")

    # Section 4: SV Profile Evolution — THE PRIMARY RESULT
    lines.append("## 4. SV Profile Evolution — THE PRIMARY RESULT\n")
    primary_L = 32 if 32 in all_results else min(all_results.keys())
    lines.append(f"SV profiles (mean normalized) at key temperatures (L={primary_L} primary):\n")

    if primary_L in all_results:
        for target in [1.50, 1.30, 1.10, 1.02, 1.00, 0.90]:
            best = min(all_results[primary_L], key=lambda r: abs(r['T_over_Tc'] - target))
            sv = best.get('fisher_sv_profile_mean', [])
            sv_str = ", ".join([f"{v:.3f}" for v in sv])
            pattern = identify_degeneracy_pattern(sv)
            lines.append(f"- T/Tc={best['T_over_Tc']:.3f}: [{sv_str}]  → {pattern}")
        lines.append("")

        # SWAP ASSESSMENT
        lines.append("### Degeneracy Swap Assessment\n")
        # Compare high-T and near-Tc patterns
        high_T = min(all_results[primary_L], key=lambda r: abs(r['T_over_Tc'] - 1.50))
        near_Tc = min(all_results[primary_L], key=lambda r: abs(r['T_over_Tc'] - 1.02))
        at_Tc = min(all_results[primary_L], key=lambda r: abs(r['T_over_Tc'] - 1.00))

        pat_high = identify_degeneracy_pattern(high_T.get('fisher_sv_profile_mean', []))
        pat_near = identify_degeneracy_pattern(near_Tc.get('fisher_sv_profile_mean', []))
        pat_at = identify_degeneracy_pattern(at_Tc.get('fisher_sv_profile_mean', []))

        lines.append(f"- High T (T/Tc=1.50): {pat_high}")
        lines.append(f"- Near Tc (T/Tc≈1.02): {pat_near}")
        lines.append(f"- At Tc (T/Tc=1.00): {pat_at}")
        lines.append("")

        swap_observed = (pat_high != pat_near) or (pat_high != pat_at)
        if swap_observed:
            lines.append("**RESULT: DEGENERACY SWAP OBSERVED** ✓")
        else:
            lines.append("**RESULT: NO DEGENERACY SWAP OBSERVED** ✗")
        lines.append("")

    # Section 5: 2D vs 3D Comparison
    lines.append("## 5. 2D vs 3D Comparison\n")
    lines.append("| Observable | 2D Ising (Phase 2) | 3D Ising (Phase 3D-1) |")
    lines.append("|---|---|---|")
    if primary_L in all_results:
        tc_3d = min(all_results[primary_L], key=lambda r: abs(r['T_over_Tc'] - 1.0))
        below_3d = min(all_results[primary_L], key=lambda r: abs(r['T_over_Tc'] - 0.9))
        rank_3d = tc_3d.get('fisher_rank_mean', 0)
        eta_tc_3d = tc_3d.get('fisher_eta_mean', 0)
        eta_below_3d = below_3d.get('fisher_eta_mean', 0)
        chi_peak_t = max(all_results[primary_L], key=lambda r: r.get('macro_susceptibility', 0))
        chi_peak_ratio = chi_peak_t['T_over_Tc']

        lines.append(f"| Rank | 3 | {rank_3d:.1f} |")
        lines.append(f"| η at T_c (primary) | 0.147 | {eta_tc_3d:.3f} |")
        lines.append(f"| η below T_c | 0.465 | {eta_below_3d:.3f} |")
        lines.append(f"| Swap observed? | YES | {'YES' if swap_observed else 'NO'} |")

        # Swap range
        if swap_observed:
            # Find swap range by looking at where pattern changes
            lines.append(f"| Swap range (T/T_c) | 1.10–1.04 | SEE PROFILES |")
        else:
            lines.append(f"| Swap range (T/T_c) | 1.10–1.04 | N/A |")

        lines.append(f"| χ peak (T/T_c) | 1.005 | {chi_peak_ratio:.3f} |")

        if primary_L in verdicts:
            T_f, T_c_onset, v = verdicts[primary_L]
            lines.append(f"| Leading indicator | YES (ΔT/Tc=0.13) | {v} |")
    lines.append("")

    # Section 6: Leading Indicator Assessment
    lines.append("## 6. Leading Indicator Assessment\n")
    for L_val, (T_f, T_c_onset, verdict) in sorted(verdicts.items()):
        t_f_str = f"{T_f:.3f}" if T_f is not None else "N/A"
        t_c_str = f"{T_c_onset:.3f}" if T_c_onset is not None else "N/A"
        lines.append(f"- L={L_val}: T_Fisher={t_f_str}, T_chi={t_c_str} → **{verdict}**")
    lines.append("")

    # Section 7: Falsifiability Assessment
    lines.append("## 7. Falsifiability Assessment\n")
    lines.append("| ID | Prediction | Confidence | Result |")
    lines.append("|---|---|---|---|")

    if primary_L in all_results:
        # P3D-1: Rank = 4
        ranks_all = [r.get('fisher_rank_mean', 0) for r in all_results[primary_L]]
        rank_mode = Counter([round(r) for r in ranks_all]).most_common(1)[0][0]
        p1 = "PASS" if rank_mode == 4 else f"FAIL (rank={rank_mode})"
        lines.append(f"| P3D-1 | Rank = 4 at all T | Tier 3 ‡_p | {p1} |")

        # P3D-2: Swap observed (KILL TEST)
        p2 = "PASS" if swap_observed else "**FAIL — KILL**"
        lines.append(f"| P3D-2 | SV degeneracy swap occurs | Tier 3 ‡_p | {p2} |")

        # P3D-3: η minimum at Tc
        tc_eta = tc_3d.get('fisher_eta_mean', 1.0)
        other_etas = [r.get('fisher_eta_mean', 1.0) for r in all_results[primary_L]
                     if abs(r['T_over_Tc'] - 1.0) > 0.01]
        eta_min_at_tc = tc_eta <= min(other_etas) if other_etas else False
        p3 = "PASS" if eta_min_at_tc else f"FAIL (η_Tc={tc_eta:.3f}, min_other={min(other_etas):.3f})"
        lines.append(f"| P3D-3 | η minimum at T_c | Tier 3 ‡_p | {p3} |")

        # P3D-4: Fisher leads susceptibility
        if primary_L in verdicts:
            _, _, v = verdicts[primary_L]
            p4 = "PASS" if "LEADING" in v else f"FAIL ({v})"
            lines.append(f"| P3D-4 | Fisher leads susceptibility | Tier 3 ‡_p | {p4} |")
    lines.append("")

    # Section 8: Runtime Log
    lines.append("## 8. Runtime Log\n")
    for L_val in sorted(all_results.keys()):
        times = [r['time_sec'] for r in all_results[L_val]]
        lines.append(f"- L={L_val}: {sum(times):.1f}s total ({np.mean(times):.1f}s per T point)")
    lines.append(f"- **Grand total: {total_time:.1f}s**")

    report_path = os.path.join(OUTPUT_DIR, "PHASE3D1_3D_ISING_RESULTS.md")
    with open(report_path, 'w') as f:
        f.write('\n'.join(lines))
    print(f"  Saved {report_path}")


# ============================================================
# MAIN
# ============================================================

def main():
    print("=" * 70)
    print("DS Phase 3D-1: 3D Ising Universality Test")
    print("=" * 70)
    print(f"T_c (3D) = {T_C_3D:.4f}")
    print(f"Output: {OUTPUT_DIR}/")
    print()

    temperatures = T_C_3D * np.array([
        1.50, 1.30, 1.20, 1.15, 1.10, 1.07,
        1.04, 1.02, 1.01, 1.005, 1.000, 0.900,
    ])

    # Lattice sizes: L=16 (fast), L=32 (primary), L=48 (stretch)
    lattice_sizes = [16, 32]
    all_results = {}
    verdicts = {}
    t_total_start = time.time()

    for L_val in lattice_sizes:
        print(f"\n{'=' * 70}")
        print(f"  3D Lattice L = {L_val}  ({L_val**3} sites)")
        print(f"{'=' * 70}")

        results = temperature_sweep_3d(L_val, temperatures,
                                        n_eq_sweeps=200, n_configs=500)
        all_results[L_val] = results

        # Save raw data
        raw_path = os.path.join(RAW_DIR, f"results_L{L_val}.json")
        with open(raw_path, 'w') as f:
            json.dump(results, f, indent=2, default=str)
        print(f"  Raw data saved to {raw_path}")

    elapsed = time.time() - t_total_start
    print(f"\nL=16 + L=32 completed in {elapsed:.1f}s")

    # Decide whether to run L=48
    if elapsed < 900:  # 15 min budget for L=16+L=32
        print(f"\nWithin budget ({elapsed:.0f}s < 900s). Running L=48 stretch goal...")
        L_val = 48
        print(f"\n{'=' * 70}")
        print(f"  3D Lattice L = {L_val}  ({L_val**3} sites)")
        print(f"{'=' * 70}")

        # Reduce n_configs if tight on time
        nc = 300 if elapsed > 600 else 500
        results = temperature_sweep_3d(L_val, temperatures,
                                        n_eq_sweeps=200, n_configs=nc)
        all_results[L_val] = results

        raw_path = os.path.join(RAW_DIR, f"results_L{L_val}.json")
        with open(raw_path, 'w') as f:
            json.dump(results, f, indent=2, default=str)
        print(f"  Raw data saved to {raw_path}")
    else:
        print(f"\nOver budget ({elapsed:.0f}s >= 900s). Skipping L=48.")

    total_time = time.time() - t_total_start

    # ---- PLOTS ----
    print(f"\n{'=' * 70}")
    print("  GENERATING PLOTS")
    print(f"{'=' * 70}")

    for L_val in sorted(all_results.keys()):
        plot_diagnostic_panel(all_results[L_val], L_val)
        plot_sv_evolution(all_results[L_val], L_val)
        T_f, T_c_onset, verdict = plot_leading_indicator(all_results[L_val], L_val)
        verdicts[L_val] = (T_f, T_c_onset, verdict)

    # SV comparison with 2D (using L=32 primary)
    primary_L = 32 if 32 in all_results else min(all_results.keys())
    plot_sv_comparison_2d_vs_3d(all_results[primary_L], primary_L)

    # PR vs r_max for L=32 (skip L=16 — too small for meaningful r_max range)
    if 32 in all_results:
        plot_pr_rmax(all_results[32], 32)

    # Finite-size scaling
    if len(all_results) >= 2:
        plot_finite_size_scaling(all_results)

    # Results report
    write_results_report(all_results, verdicts, total_time)

    print(f"\n{'=' * 70}")
    print(f"  TOTAL TIME: {total_time:.1f}s")
    print(f"{'=' * 70}")

    # Print headline result
    print(f"\n{'*' * 70}")
    primary_results = all_results.get(primary_L, [])
    if primary_results:
        high_T = min(primary_results, key=lambda r: abs(r['T_over_Tc'] - 1.50))
        near_Tc = min(primary_results, key=lambda r: abs(r['T_over_Tc'] - 1.02))
        at_Tc = min(primary_results, key=lambda r: abs(r['T_over_Tc'] - 1.00))

        pat_high = identify_degeneracy_pattern(high_T.get('fisher_sv_profile_mean', []))
        pat_near = identify_degeneracy_pattern(near_Tc.get('fisher_sv_profile_mean', []))
        pat_at = identify_degeneracy_pattern(at_Tc.get('fisher_sv_profile_mean', []))

        swap = (pat_high != pat_near) or (pat_high != pat_at)

        print(f"  HEADLINE: {'DEGENERACY SWAP OBSERVED' if swap else 'NO SWAP — GENERALIZATION DEAD'}")
        print(f"  High T pattern: {pat_high}")
        print(f"  Near Tc pattern: {pat_near}")
        print(f"  At Tc pattern: {pat_at}")
    print(f"{'*' * 70}")


if __name__ == '__main__':
    main()
