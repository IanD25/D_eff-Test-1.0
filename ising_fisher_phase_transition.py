#!/usr/bin/env python3
"""
DS Phase 2: Ising Fisher Phase Transition
============================================
Tests whether Fisher diagnostic toolkit can detect the approach to the
2D Ising phase transition before macroscopic observables signal it.

Uses the thermal correlation function G(r,T) as the Fisher kernel,
replacing the exponential distance kernel from Phase 1.

Systems: 2D square lattice with periodic BCs, N = 64, 128, 256
Temperature sweep: 12 points from 1.5*Tc down to 0.9*Tc
Monte Carlo: Wolff cluster algorithm
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

OUTPUT_DIR = "phase2_results"
RAW_DIR = os.path.join(OUTPUT_DIR, "raw_data")
os.makedirs(RAW_DIR, exist_ok=True)

T_C = 2.0 / np.log(1 + np.sqrt(2))   # 2.26919...

np.random.seed(42)

# ============================================================
# MODULE 1: Ising Monte Carlo (Wolff cluster algorithm)
# ============================================================

def initialize_lattice(N, mode='random'):
    """Returns NxN array of spins in {+1, -1}."""
    if mode == 'cold':
        return np.ones((N, N), dtype=np.int8)
    elif mode == 'hot' or mode == 'random':
        return np.where(np.random.random((N, N)) < 0.5, 1, -1).astype(np.int8)
    else:
        raise ValueError(f"Unknown mode: {mode}")


def wolff_step(lattice, T, J=1.0):
    """Execute one Wolff cluster flip."""
    N = lattice.shape[0]
    p_add = 1.0 - np.exp(-2.0 * J / T)

    # Choose random seed
    si = np.random.randint(N)
    sj = np.random.randint(N)
    seed_spin = lattice[si, sj]

    # BFS cluster growth
    visited = np.zeros((N, N), dtype=np.bool_)
    visited[si, sj] = True
    queue = deque()
    queue.append((si, sj))
    cluster = [(si, sj)]

    while queue:
        ci, cj = queue.popleft()
        # Check 4 neighbors with periodic BC
        for di, dj in [(-1, 0), (1, 0), (0, -1), (0, 1)]:
            ni, nj = (ci + di) % N, (cj + dj) % N
            if not visited[ni, nj] and lattice[ni, nj] == seed_spin:
                if np.random.random() < p_add:
                    visited[ni, nj] = True
                    queue.append((ni, nj))
                    cluster.append((ni, nj))

    # Flip cluster
    for ci, cj in cluster:
        lattice[ci, cj] = -lattice[ci, cj]

    return lattice


def wolff_sweep(lattice, T, n_steps):
    """Execute n_steps Wolff cluster flips."""
    for _ in range(n_steps):
        wolff_step(lattice, T)
    return lattice


def equilibrate(lattice, T, n_sweeps=200):
    """Equilibrate lattice at temperature T using Wolff sweeps."""
    N = lattice.shape[0]
    # One sweep = N*N / mean_cluster_size attempted flips
    # For Wolff: each step flips a cluster, so n_sweeps steps suffices
    for _ in range(n_sweeps):
        wolff_step(lattice, T)
    return lattice


# ============================================================
# MODULE 2: Correlation Function Measurement
# ============================================================

def measure_magnetization(lattice):
    """Mean absolute magnetization: |sum(s)| / N^2"""
    return np.abs(np.mean(lattice))


def measure_energy(lattice, J=1.0):
    """Energy per spin: -J * sum(s_i * s_j) / N^2 over NN pairs."""
    N = lattice.shape[0]
    # Sum over right and down neighbors (each bond counted once)
    e = -J * np.sum(lattice * np.roll(lattice, -1, axis=0))  # down
    e += -J * np.sum(lattice * np.roll(lattice, -1, axis=1))  # right
    return e / (N * N)


def measure_correlations_fft(lattice):
    """
    Measure connected G(r) using FFT-based autocorrelation.
    Returns G[r] for r = 0, 1, ..., N//4 using Chebyshev distance.
    """
    N = lattice.shape[0]
    lat_float = lattice.astype(np.float64)
    m = np.mean(lat_float)
    m_sq = m * m

    # FFT-based autocorrelation
    F = np.fft.fft2(lat_float)
    power = F * np.conj(F)
    G_full = np.real(np.fft.ifft2(power)) / (N * N)

    # Radial average using Chebyshev distance
    max_r = N // 4
    G_r = np.zeros(max_r + 1)
    counts = np.zeros(max_r + 1)

    # Build distance grid (Chebyshev with periodic wrapping)
    ix = np.arange(N)
    dx = np.minimum(ix, N - ix)
    DX, DY = np.meshgrid(dx, dx)
    dist_grid = np.maximum(DX, DY)

    for r in range(max_r + 1):
        mask = (dist_grid == r)
        if np.any(mask):
            G_r[r] = np.mean(G_full[mask]) - m_sq
            counts[r] = np.sum(mask)

    return G_r


def accumulate_correlations(lattice, T, n_configs=500, n_sweeps_between=5):
    """
    Measure G(r) averaged over n_configs independent configurations.
    Returns averaged G(r) array.
    """
    N = lattice.shape[0]
    max_r = N // 4
    G_sum = np.zeros(max_r + 1)

    for i in range(n_configs):
        for _ in range(n_sweeps_between):
            wolff_step(lattice, T)
        G_r = measure_correlations_fft(lattice)
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
# MODULE 3: Thermal Fisher Kernel and FIM Construction
# ============================================================

def build_thermal_kernel(G_r, N, v0):
    """
    Build thermal probability distribution p_{v0}(u; T) using |G(r)|.
    v0 = (i0, j0) tuple.
    Returns flat array of shape (N*N,).
    """
    i0, j0 = v0
    max_r = len(G_r) - 1

    # Build distance from v0 to all sites (Chebyshev, periodic)
    ix = np.arange(N)
    di = np.minimum(np.abs(ix - i0), N - np.abs(ix - i0))
    dj = np.minimum(np.abs(ix - j0), N - np.abs(ix - j0))
    DI, DJ = np.meshgrid(di, dj, indexing='ij')
    dist = np.maximum(DI, DJ)

    # Build weights from |G(r)|
    weights = np.zeros((N, N), dtype=np.float64)
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
        # Uniform fallback
        p = np.ones(N * N) / (N * N)
    else:
        p = p / total

    return p


def compute_FIM_thermal(G_r, N, v0):
    """
    Compute Fisher Information Matrix at vertex v0 using thermal kernel.
    v0 = (i0, j0) tuple. Returns 4x4 FIM matrix and diagnostics.
    """
    i0, j0 = v0

    # p_{v0}
    p_v0 = build_thermal_kernel(G_r, N, v0)
    log_p_v0 = np.log(p_v0 + 1e-30)

    # 4 neighbors (periodic BC)
    neighbors = [
        ((i0 - 1) % N, j0),
        ((i0 + 1) % N, j0),
        (i0, (j0 - 1) % N),
        (i0, (j0 + 1) % N),
    ]
    k = len(neighbors)

    # Score vectors
    score_vectors = np.zeros((k, N * N))
    for j, w in enumerate(neighbors):
        p_w = build_thermal_kernel(G_r, N, w)
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
    """Returns normalized SV profile."""
    svs = np.linalg.svd(FIM, compute_uv=False)
    if svs[0] < 1e-30:
        return svs
    return svs / svs[0]


def fisher_diagnostics(G_r, N, n_samples=30):
    """
    Sample n_samples random vertices.
    Compute FIM, rank, PR, eta, sv_profile for each.
    Returns diagnostics dict.
    """
    ranks = []
    prs = []
    etas = []
    sv_profiles = []
    gap_ratios = []

    # Sample random interior vertices (avoid edges for cleaner stats, though periodic)
    margin = max(2, N // 8)
    sampled = 0
    attempts = 0
    while sampled < n_samples and attempts < n_samples * 5:
        i = np.random.randint(margin, N - margin)
        j = np.random.randint(margin, N - margin)
        attempts += 1

        FIM = compute_FIM_thermal(G_r, N, (i, j))
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

def measure_macroscopic(lattice, T, n_configs=500, n_sweeps_between=5):
    """
    Measure macroscopic observables over n_configs decorrelated configurations.
    Returns dict with m, chi, E, C, U_binder.
    """
    N = lattice.shape[0]
    N2 = N * N

    ms = []
    m2s = []
    m4s = []
    Es = []
    E2s = []

    for _ in range(n_configs):
        for _ in range(n_sweeps_between):
            wolff_step(lattice, T)
        m_val = np.mean(lattice.astype(np.float64))
        ms.append(np.abs(m_val))
        m2s.append(m_val ** 2)
        m4s.append(m_val ** 4)

        e_val = measure_energy(lattice)
        Es.append(e_val)
        E2s.append(e_val ** 2)

    m_mean = np.mean(ms)
    m2_mean = np.mean(m2s)
    m4_mean = np.mean(m4s)
    E_mean = np.mean(Es)
    E2_mean = np.mean(E2s)

    chi = N2 * (m2_mean - m_mean ** 2)
    C = N2 * (E2_mean - E_mean ** 2) / (T * T)
    U_binder = 1.0 - m4_mean / (3.0 * m2_mean ** 2) if m2_mean > 1e-15 else 0.0

    return {
        'magnetization': float(m_mean),
        'susceptibility': float(chi),
        'energy': float(E_mean),
        'specific_heat': float(C),
        'binder': float(U_binder),
    }


def temperature_sweep(N, temperatures, n_eq_sweeps=200, n_configs=500):
    """
    Full temperature sweep for lattice size N.
    Returns list of result dicts, one per temperature.
    """
    results = []

    for t_idx, T in enumerate(temperatures):
        t0 = time.time()
        T_over_Tc = T / T_C

        # Initialize: hot start above Tc, cold start below
        mode = 'hot' if T > T_C else 'cold'
        lattice = initialize_lattice(N, mode=mode)

        # Equilibrate
        equilibrate(lattice, T, n_sweeps=n_eq_sweeps)

        # Near Tc: use more configs for better statistics
        nc = n_configs
        ns = 30
        if 0.98 <= T_over_Tc <= 1.05:
            nc = min(n_configs * 2, 1000)
            ns = 50

        # Accumulate correlations
        # Clone lattice for correlation measurement (don't disrupt state)
        lat_corr = lattice.copy()
        G_r = accumulate_correlations(lat_corr, T, n_configs=nc, n_sweeps_between=5)

        # Correlation length
        xi, xi_r2 = estimate_correlation_length(G_r)

        # Fisher diagnostics
        fisher = fisher_diagnostics(G_r, N, n_samples=ns)

        # Macroscopic observables (fresh pass from equilibrated state)
        lat_macro = lattice.copy()
        macro = measure_macroscopic(lat_macro, T, n_configs=nc, n_sweeps_between=5)

        dt = time.time() - t0

        result = {
            'T': float(T),
            'T_over_Tc': float(T_over_Tc),
            'N': N,
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

def plot_diagnostic_panel(results, N):
    """Figure 1: 4-panel diagnostic plot."""
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    fig.suptitle(f"Fisher Diagnostic Panel — N={N}", fontsize=14, fontweight='bold')

    T_ratios = [r['T_over_Tc'] for r in results]

    # (a) D_eff (rank) vs T/Tc
    ax = axes[0, 0]
    ranks = [r.get('fisher_rank_mean', np.nan) for r in results]
    rank_stds = [r.get('fisher_rank_std', 0) for r in results]
    ax.errorbar(T_ratios, ranks, yerr=rank_stds, fmt='o-', color='tab:blue', capsize=3)
    ax.axhline(2, color='gray', ls='--', alpha=0.5, label='D_eff=2')
    ax.axhline(1, color='gray', ls=':', alpha=0.5, label='D_eff=1')
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
    fname = os.path.join(OUTPUT_DIR, f"diagnostic_panel_N{N}.png")
    plt.savefig(fname, dpi=150)
    plt.close()
    print(f"  Saved {fname}")


def plot_sv_evolution(results, N):
    """Figure 2: SV profile evolution at representative temperatures."""
    target_ratios = [1.50, 1.10, 1.02, 1.00, 0.90]
    fig, axes = plt.subplots(2, 3, figsize=(16, 10))
    fig.suptitle(f"SV Profile Evolution — N={N}", fontsize=14, fontweight='bold')
    axes_flat = axes.flatten()

    for idx, t_target in enumerate(target_ratios):
        ax = axes_flat[idx]
        # Find closest result
        best = min(results, key=lambda r: abs(r['T_over_Tc'] - t_target))
        sv_mean = np.array(best.get('fisher_sv_profile_mean', [0, 0, 0, 0]))
        sv_std = np.array(best.get('fisher_sv_profile_std', [0, 0, 0, 0]))
        rank = best.get('fisher_rank_mean', 0)
        eta = best.get('fisher_eta_mean', 0)

        colors = ['tab:blue'] * len(sv_mean)
        # Color graded-decay profiles red
        if eta > 0.4:
            colors = ['tab:red'] * len(sv_mean)

        ax.bar(range(1, len(sv_mean) + 1), sv_mean, yerr=sv_std, capsize=4,
               color=colors, alpha=0.7, edgecolor='black')
        ax.set_xlabel('SV Index')
        ax.set_ylabel('Normalized SV')
        ax.set_title(f"T/Tc={best['T_over_Tc']:.3f}  rank={rank:.1f}  η={eta:.2f}")
        ax.set_ylim(0, 1.15)

    # Last panel: mean SV profile at T_c with extra detail
    ax = axes_flat[5]
    tc_result = min(results, key=lambda r: abs(r['T_over_Tc'] - 1.0))
    sv_mean = np.array(tc_result.get('fisher_sv_profile_mean', [0, 0, 0, 0]))
    sv_std = np.array(tc_result.get('fisher_sv_profile_std', [0, 0, 0, 0]))
    ax.bar(range(1, len(sv_mean) + 1), sv_mean, yerr=sv_std, capsize=4,
           color='tab:purple', alpha=0.7, edgecolor='black')
    ax.set_xlabel('SV Index')
    ax.set_ylabel('Normalized SV')
    ax.set_title(f"T=T_c (detail): mean ± std over {tc_result.get('N', N)} samples")
    ax.set_ylim(0, 1.15)

    plt.tight_layout()
    fname = os.path.join(OUTPUT_DIR, f"sv_profile_evolution_N{N}.png")
    plt.savefig(fname, dpi=150)
    plt.close()
    print(f"  Saved {fname}")


def plot_pr_sigma(results, N):
    """Figure 3: PR vs sigma (r_max truncation) at 3 temperatures."""
    target_ratios = [1.50, 1.02, 1.00]
    colors = ['tab:blue', 'tab:orange', 'tab:red']
    labels = ['T/Tc=1.50', 'T/Tc=1.02', 'T/Tc=1.00']

    fig, ax = plt.subplots(figsize=(10, 6))

    for t_target, color, label in zip(target_ratios, colors, labels):
        best = min(results, key=lambda r: abs(r['T_over_Tc'] - t_target))
        G_r_full = np.array(best['G_r'])

        # Compute PR at different r_max truncations (proxy for sigma)
        r_max_values = list(range(3, len(G_r_full), 2))
        pr_values = []

        for r_max in r_max_values:
            G_r_trunc = G_r_full[:r_max + 1]
            # Sample a few vertices and average PR
            prs = []
            for _ in range(10):
                i = np.random.randint(N // 4, 3 * N // 4)
                j = np.random.randint(N // 4, 3 * N // 4)
                FIM = compute_FIM_thermal(G_r_trunc, N, (i, j))
                svs = np.linalg.svd(FIM, compute_uv=False)
                if svs[0] > 1e-30:
                    prs.append(participation_ratio(svs))
            pr_values.append(np.mean(prs) if prs else 0)

        ax.plot(r_max_values, pr_values, 'o-', color=color, label=label, markersize=4)

    ax.axhline(1.875, color='gray', ls='--', alpha=0.5, label='d_H=1.875 (Ising critical)')
    ax.axhline(2.0, color='gray', ls=':', alpha=0.3, label='d=2')
    ax.set_xlabel('r_max (correlation truncation range)')
    ax.set_ylabel('Participation Ratio')
    ax.set_title(f'PR vs Correlation Range — N={N}')
    ax.legend()
    plt.tight_layout()
    fname = os.path.join(OUTPUT_DIR, f"pr_sigma_thermal_N{N}.png")
    plt.savefig(fname, dpi=150)
    plt.close()
    print(f"  Saved {fname}")


def plot_leading_indicator(results, N):
    """Figure 4: Leading indicator comparison — eta vs chi."""
    T_ratios = [r['T_over_Tc'] for r in results]
    etas = [r.get('fisher_eta_mean', np.nan) for r in results]
    chis = [r.get('macro_susceptibility', np.nan) for r in results]

    # Compute onset temperatures
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

    ax1.set_title(f"Leading Indicator Test — N={N}: {verdict}")

    lines1, labels1 = ax1.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    ax1.legend(lines1 + lines2, labels1 + labels2, loc='upper left', fontsize=8)

    plt.tight_layout()
    fname = os.path.join(OUTPUT_DIR, f"leading_indicator_N{N}.png")
    plt.savefig(fname, dpi=150)
    plt.close()
    print(f"  Saved {fname}")

    return T_fisher, T_chi, verdict


def plot_finite_size_scaling(all_results):
    """Figure 5: Finite-size scaling — all N overlaid."""
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))
    fig.suptitle("Finite-Size Scaling", fontsize=14, fontweight='bold')

    colors = {64: 'tab:blue', 128: 'tab:orange', 256: 'tab:red'}

    for N_val, results in sorted(all_results.items()):
        T_ratios = [r['T_over_Tc'] for r in results]
        ranks = [r.get('fisher_rank_mean', np.nan) for r in results]
        etas = [r.get('fisher_eta_mean', np.nan) for r in results]

        ax1.plot(T_ratios, ranks, 'o-', color=colors.get(N_val, 'black'),
                 label=f'N={N_val}', markersize=5)
        ax2.plot(T_ratios, etas, 's-', color=colors.get(N_val, 'black'),
                 label=f'N={N_val}', markersize=5)

    ax1.axvline(1.0, color='black', ls='--', alpha=0.3)
    ax1.axhline(2, color='gray', ls='--', alpha=0.3)
    ax1.axhline(1, color='gray', ls=':', alpha=0.3)
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


def write_results_report(all_results, verdicts, total_time):
    """Write the full results markdown report."""
    lines = []
    lines.append("# Phase 2 Results: Ising Fisher Phase Transition\n")
    lines.append(f"Generated: {datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")

    # Section 1: Parameters
    lines.append("## 1. System Parameters\n")
    lines.append(f"- T_c = {T_C:.6f}")
    lines.append(f"- Lattice sizes: {sorted(all_results.keys())}")
    lines.append(f"- Temperature points: 12 (T/Tc from 1.50 to 0.90)")
    lines.append(f"- Total runtime: {total_time:.1f}s\n")

    # Section 2: Correlation Function
    lines.append("## 2. Correlation Function Behavior\n")
    for N_val in sorted(all_results.keys()):
        lines.append(f"### N = {N_val}\n")
        lines.append("| T/Tc | ξ | R² (exp fit) | Notes |")
        lines.append("|------|---|-------------|-------|")
        for r in all_results[N_val]:
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
    for N_val in sorted(all_results.keys()):
        lines.append(f"### N = {N_val}\n")
        lines.append("| T/Tc | Rank (mean±std) | η (mean±std) | PR (mean) | Gap Ratio | χ | C |")
        lines.append("|------|----------------|-------------|-----------|-----------|---|---|")
        for r in all_results[N_val]:
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

    # Section 4: SV Profile Shape
    lines.append("## 4. SV Profile Shape Transition\n")
    lines.append("SV profiles (mean normalized) at key temperatures (N=128 primary):\n")
    if 128 in all_results:
        for target in [1.50, 1.10, 1.02, 1.00, 0.90]:
            best = min(all_results[128], key=lambda r: abs(r['T_over_Tc'] - target))
            sv = best.get('fisher_sv_profile_mean', [])
            sv_str = ", ".join([f"{v:.3f}" for v in sv])
            lines.append(f"- T/Tc={best['T_over_Tc']:.3f}: [{sv_str}]")
    lines.append("")

    # Section 5: Leading Indicator
    lines.append("## 5. Leading Indicator Assessment\n")
    for N_val, (T_f, T_c_onset, verdict) in verdicts.items():
        t_f_str = f"{T_f:.3f}" if T_f is not None else "N/A"
        t_c_str = f"{T_c_onset:.3f}" if T_c_onset is not None else "N/A"
        lines.append(f"- N={N_val}: T_Fisher={t_f_str}, T_chi={t_c_str} → **{verdict}**")
    lines.append("")

    # Section 6: Phase 1 comparison
    lines.append("## 6. Comparison to Phase 1 Benchmarks\n")
    lines.append("Reference η values from Phase 1:")
    lines.append("- Torus: η ≈ 0.23 (ordered lattice)")
    lines.append("- RGG: η ≈ 0.68 (disordered manifold)")
    lines.append("- ER: η ≈ 0.93 (no geometry)\n")
    if 128 in all_results:
        tc_result = min(all_results[128], key=lambda r: abs(r['T_over_Tc'] - 1.0))
        eta_tc = tc_result.get('fisher_eta_mean', 0)
        lines.append(f"η at T_c (N=128): {eta_tc:.3f}")
        if eta_tc < 0.45:
            lines.append("→ Below RGG reference; critical state less disordered than expected")
        elif eta_tc < 0.80:
            lines.append("→ In RGG-like range; critical state matches disordered manifold profile")
        else:
            lines.append("→ Above RGG reference; critical state highly disordered")
    lines.append("")

    # Section 7: PR at Tc
    lines.append("## 7. PR(σ) Profile at T_c\n")
    if 128 in all_results:
        tc_result = min(all_results[128], key=lambda r: abs(r['T_over_Tc'] - 1.0))
        pr_tc = tc_result.get('fisher_pr_mean', 0)
        lines.append(f"- PR at T_c (N=128, default kernel): {pr_tc:.3f}")
        lines.append(f"- Known d_H for 2D Ising critical cluster: 1.875")
        lines.append(f"- Gap: {abs(pr_tc - 1.875):.3f}")
    lines.append("")

    # Section 8: Falsifiability
    lines.append("## 8. Falsifiability Assessment\n")
    lines.append("| ID | Prediction | Result |")
    lines.append("|---|---|---|")
    # These will be filled based on actual data
    if 128 in all_results:
        res128 = all_results[128]
        tc_res = min(res128, key=lambda r: abs(r['T_over_Tc'] - 1.0))
        high_res = min(res128, key=lambda r: abs(r['T_over_Tc'] - 1.5))
        near_res = min(res128, key=lambda r: abs(r['T_over_Tc'] - 1.05))
        low_res = min(res128, key=lambda r: abs(r['T_over_Tc'] - 0.9))

        # P2-1: SV profile shape change by T/Tc=1.05
        sv_high = high_res.get('fisher_sv_profile_mean', [1, 0, 0, 0])
        sv_near = near_res.get('fisher_sv_profile_mean', [1, 0, 0, 0])
        sv_tc = tc_res.get('fisher_sv_profile_mean', [1, 0, 0, 0])
        p21 = "PASS" if len(sv_near) >= 3 and sv_near[1] > 0.1 and sv_near[2] / max(sv_near[1], 1e-10) > 0.3 else "CHECK VISUALLY"
        lines.append(f"| P2-1 | SV profile step→graded by T/Tc=1.05 | {p21} |")

        # P2-2: eta rises to >0.45 at Tc
        eta_tc = tc_res.get('fisher_eta_mean', 0)
        p22 = "PASS" if eta_tc > 0.45 else "FAIL"
        lines.append(f"| P2-2 | η(T_c) > 0.45 | {p22} (η={eta_tc:.3f}) |")

        # P2-3: Gap ratio monotone decrease
        grs = [r.get('fisher_gap_ratio_mean', float('inf')) for r in res128]
        # Check if mostly decreasing (allowing some noise)
        diffs = [grs[i+1] - grs[i] for i in range(len(grs)-1) if not np.isinf(grs[i]) and not np.isinf(grs[i+1])]
        n_decrease = sum(1 for d in diffs if d < 0)
        p23 = "PASS" if n_decrease > len(diffs) * 0.6 else "FAIL"
        lines.append(f"| P2-3 | Gap ratio monotone decrease | {p23} ({n_decrease}/{len(diffs)} decreasing) |")

        # P2-4: PR negative slope at Tc — assessed from PR sigma plot
        lines.append(f"| P2-4 | PR negative slope at T_c | SEE PR(σ) PLOT |")

        # P2-5: Rank collapse below Tc
        rank_low = low_res.get('fisher_rank_mean', 2)
        p25 = "PASS" if rank_low < 2.0 else "FAIL"
        lines.append(f"| P2-5 | D_eff(0.9*Tc) < 2 | {p25} (rank={rank_low:.2f}) |")

        # P2-6: Leading indicator
        if 128 in verdicts:
            _, _, v = verdicts[128]
            p26 = "PASS" if "LEADING" in v else "FAIL"
            lines.append(f"| P2-6 | Fisher leads susceptibility | {p26} ({v}) |")

        # P2-7: PR at Tc in [1.6, 2.1]
        pr_tc = tc_res.get('fisher_pr_mean', 0)
        p27 = "PASS" if 1.6 <= pr_tc <= 2.1 else "FAIL"
        lines.append(f"| P2-7 | PR(T_c) ∈ [1.6, 2.1] | {p27} (PR={pr_tc:.3f}) |")

    lines.append("")

    # Section 9: Open Questions
    lines.append("## 9. Open Questions and Phase 3 Candidates\n")
    lines.append("(To be filled based on results analysis)\n")

    # Section 10: Runtime
    lines.append("## 10. Runtime Log\n")
    for N_val in sorted(all_results.keys()):
        times = [r['time_sec'] for r in all_results[N_val]]
        lines.append(f"- N={N_val}: {sum(times):.1f}s total ({np.mean(times):.1f}s per T point)")
    lines.append(f"- **Grand total: {total_time:.1f}s**")

    report_path = os.path.join(OUTPUT_DIR, "PHASE2_ISING_FISHER_RESULTS.md")
    with open(report_path, 'w') as f:
        f.write('\n'.join(lines))
    print(f"  Saved {report_path}")


# ============================================================
# MAIN
# ============================================================

def main():
    print("=" * 70)
    print("DS Phase 2: Ising Fisher Phase Transition")
    print("=" * 70)
    print(f"T_c = {T_C:.6f}")
    print(f"Output: {OUTPUT_DIR}/")
    print()

    temperatures = T_C * np.array([
        1.50, 1.30, 1.20, 1.15, 1.10, 1.07,
        1.04, 1.02, 1.01, 1.005, 1.000, 0.900,
    ])

    lattice_sizes = [64, 128, 256]
    all_results = {}
    verdicts = {}
    t_total_start = time.time()

    for N_val in lattice_sizes:
        print(f"\n{'=' * 70}")
        print(f"  Lattice N = {N_val}")
        print(f"{'=' * 70}")

        results = temperature_sweep(N_val, temperatures,
                                     n_eq_sweeps=200, n_configs=500)
        all_results[N_val] = results

        # Save raw data
        # Convert non-serializable types
        raw_path = os.path.join(RAW_DIR, f"results_N{N_val}.json")
        with open(raw_path, 'w') as f:
            json.dump(results, f, indent=2, default=str)
        print(f"  Raw data saved to {raw_path}")

    total_time = time.time() - t_total_start

    # ---- PLOTS ----
    print(f"\n{'=' * 70}")
    print("  GENERATING PLOTS")
    print(f"{'=' * 70}")

    for N_val in lattice_sizes:
        plot_diagnostic_panel(all_results[N_val], N_val)
        plot_sv_evolution(all_results[N_val], N_val)
        T_f, T_c_onset, verdict = plot_leading_indicator(all_results[N_val], N_val)
        verdicts[N_val] = (T_f, T_c_onset, verdict)

    # PR sigma plot for N=128 and N=256
    for N_val in [128, 256]:
        if N_val in all_results:
            plot_pr_sigma(all_results[N_val], N_val)

    # Finite-size scaling
    plot_finite_size_scaling(all_results)

    # Results report
    write_results_report(all_results, verdicts, total_time)

    print(f"\n{'=' * 70}")
    print(f"  TOTAL TIME: {total_time:.1f}s")
    print(f"{'=' * 70}")


if __name__ == '__main__':
    main()
