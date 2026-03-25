#!/usr/bin/env python3
"""
Phase 4B-FSS: q=4 Potts Finite-Size Scaling — N=256 and N=512
================================================================
STANDALONE script — no relative imports, designed for remote/background execution.

Tests whether SV₂/SV₁ at T_c increases with system size N, confirming that
the low value at N=128 (0.521) is due to logarithmic corrections to scaling
rather than genuine first-order behavior.

q=4 Potts on 2D square lattice: marginal transition at T_c = 1/ln(3) ≈ 0.9102
with multiplicative logarithmic corrections (Salas & Sokal 1997).

If SV₂/SV₁ increases with N → log corrections confirmed, q=4 is continuous
If SV₂/SV₁ stays flat or decreases → genuine anomaly, needs investigation

Reference values (Phase 2, N=128):
  q=2 (Ising):  SV₂/SV₁ = 1.000 (continuous)
  q=4:          SV₂/SV₁ = 0.521 (N=128, this study)
  q=5:          SV₂/SV₁ = 0.691 (weak first-order)
  q=10:         SV₂/SV₁ = 0.374 (strong first-order)

Expected runtimes (Apple M4 / RTX 2000):
  N=256: ~30-40 min
  N=512: ~2-3 hours

Author: Ian Darling / Claude Code
Date: 2026-03-25
"""

import os
import sys
import time
import csv
import json
import datetime
import numpy as np
from collections import deque, Counter

try:
    from scipy.optimize import curve_fit
    HAS_SCIPY = True
except ImportError:
    HAS_SCIPY = False

try:
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    HAS_MATPLOTLIB = True
except ImportError:
    HAS_MATPLOTLIB = False

np.random.seed(42)

# ─── Configuration ────────────────────────────────────────────────────
Q = 4
T_C_Q4 = 1.0 / np.log(1.0 + np.sqrt(4.0))  # = 1/ln(3) ≈ 0.9102

# System sizes to test (N=128 reference included for direct comparison)
LATTICE_SIZES = [128, 256, 512]

# Temperature sweep: focused near T_c for scaling study
# Fewer points than N=128 full sweep, concentrated where it matters
T_OVER_TC_VALUES = np.array([
    0.70, 0.80, 0.90,          # Below T_c (3 points)
    0.94, 0.96, 0.98, 0.99,    # Approaching T_c (4 points)
    1.00, 1.01, 1.02,          # At/near T_c (3 points)
    1.04, 1.06, 1.10,          # Just above T_c (3 points)
    1.20, 1.35, 1.50,          # Well above T_c (3 points)
])
T_VALUES = T_OVER_TC_VALUES * T_C_Q4

# MC parameters — scale with N for comparable statistics
MC_PARAMS = {
    128:  {'n_eq': 300,  'n_configs': 500,  'n_fim': 30},
    256:  {'n_eq': 400,  'n_configs': 500,  'n_fim': 30},
    512:  {'n_eq': 500,  'n_configs': 400,  'n_fim': 25},
}

# Phase 2 reference values
REFERENCE = {
    'q=2 (Ising, N=128)': 1.000,
    'q=4 (N=128)':        0.521,
    'q=5 (N=128)':        0.691,
    'q=10 (N=128)':       0.374,
}

RESULTS_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                           "results_fss")


# =====================================================================
# MODULE 1: Potts Monte Carlo (Wolff cluster algorithm)
# =====================================================================

def initialize_potts_lattice(N, q=4, mode='random'):
    """Returns N×N array of spins in {0, 1, ..., q-1}."""
    if mode == 'cold':
        return np.zeros((N, N), dtype=np.int8)
    elif mode in ('hot', 'random'):
        return np.random.randint(0, q, size=(N, N)).astype(np.int8)
    else:
        raise ValueError(f"Unknown mode: {mode}")


def wolff_step_potts(lattice, T, q=4, J=1.0):
    """
    One Wolff cluster flip for q-state Potts model.
    P_add = 1 - exp(-J/T). Cluster: same-state neighbors bonded with P_add.
    Flip: reassign entire cluster to a uniformly random OTHER state.
    """
    N = lattice.shape[0]
    p_add = 1.0 - np.exp(-J / T)

    si = np.random.randint(N)
    sj = np.random.randint(N)
    seed_state = lattice[si, sj]

    visited = np.zeros((N, N), dtype=np.bool_)
    visited[si, sj] = True
    queue = deque()
    queue.append((si, sj))
    cluster = [(si, sj)]

    while queue:
        ci, cj = queue.popleft()
        for di, dj in [(-1, 0), (1, 0), (0, -1), (0, 1)]:
            ni, nj = (ci + di) % N, (cj + dj) % N
            if not visited[ni, nj] and lattice[ni, nj] == seed_state:
                if np.random.random() < p_add:
                    visited[ni, nj] = True
                    queue.append((ni, nj))
                    cluster.append((ni, nj))

    new_state = np.random.randint(0, q - 1)
    if new_state >= seed_state:
        new_state += 1
    new_state = np.int8(new_state)

    for ci, cj in cluster:
        lattice[ci, cj] = new_state

    return lattice


def equilibrate_potts(lattice, T, q=4, n_sweeps=300):
    """Equilibrate Potts lattice at temperature T."""
    for _ in range(n_sweeps):
        wolff_step_potts(lattice, T, q)
    return lattice


# =====================================================================
# MODULE 2: Potts Correlation Function Measurement
# =====================================================================

def measure_potts_correlations(lattice, q=4, max_r=None):
    """
    Measure Potts correlation function G(r) = <δ(s₀, s_r)> - 1/q.
    FFT approach: autocorrelation of indicator(lattice==s), summed over q states.
    Radially averaged using Manhattan distance.
    """
    N = lattice.shape[0]
    if max_r is None:
        max_r = N // 4

    N2 = N * N
    autocorr_total = np.zeros((N, N), dtype=np.float64)

    for s in range(q):
        indicator = (lattice == s).astype(np.float64)
        F = np.fft.fft2(indicator)
        power = F * np.conj(F)
        autocorr_s = np.real(np.fft.ifft2(power)) / N2
        autocorr_total += autocorr_s

    ix = np.arange(N)
    dx = np.minimum(ix, N - ix)
    DX, DY = np.meshgrid(dx, dx)
    dist_grid = DX + DY

    G_r = np.zeros(max_r + 1)
    for r in range(max_r + 1):
        mask = (dist_grid == r)
        if np.any(mask):
            G_r[r] = np.mean(autocorr_total[mask]) - 1.0 / q

    return G_r


def accumulate_potts_correlations(lattice, T, q=4, n_configs=500,
                                   n_sweeps_between=5):
    """Measure G(r) averaged over n_configs configurations."""
    N = lattice.shape[0]
    max_r = N // 4
    G_sum = np.zeros(max_r + 1)

    for i in range(n_configs):
        for _ in range(n_sweeps_between):
            wolff_step_potts(lattice, T, q)
        G_r = measure_potts_correlations(lattice, q=q, max_r=max_r)
        G_sum += G_r[:max_r + 1]

    return G_sum / n_configs


def estimate_correlation_length(G_r):
    """
    Fit G(r) ~ A * exp(-r/xi). Returns (xi, R²).
    Falls back to simple estimate if scipy unavailable.
    """
    max_fit = len(G_r) // 2
    r_vals = np.arange(1, max_fit + 1)
    g_vals = G_r[1:max_fit + 1]

    pos = g_vals > 0
    if np.sum(pos) < 3:
        return float('inf'), 0.0

    r_pos = r_vals[pos]
    log_g = np.log(g_vals[pos])

    if HAS_SCIPY:
        try:
            def exp_decay(r, log_A, inv_xi):
                return log_A - r * inv_xi
            popt, _ = curve_fit(exp_decay, r_pos, log_g, p0=[0.0, 0.1], maxfev=5000)
            inv_xi = popt[1]
            if inv_xi <= 0:
                return float('inf'), 0.0
            xi = 1.0 / inv_xi
            predicted = exp_decay(r_pos, *popt)
            ss_res = np.sum((log_g - predicted) ** 2)
            ss_tot = np.sum((log_g - np.mean(log_g)) ** 2)
            r_sq = 1.0 - ss_res / ss_tot if ss_tot > 0 else 0.0
            return xi, r_sq
        except Exception:
            pass

    # Fallback: simple slope estimate from first and last positive points
    if len(r_pos) >= 2:
        slope = (log_g[-1] - log_g[0]) / (r_pos[-1] - r_pos[0])
        if slope < 0:
            return -1.0 / slope, 0.5
    return float('inf'), 0.0


# =====================================================================
# MODULE 3: Thermal Fisher Kernel and FIM (verbatim Phase 2)
# =====================================================================

def build_thermal_kernel(G_r, N, v0):
    """
    Build thermal probability distribution p_{v0}(u; T) using |G(r)|.
    v0 = (i0, j0) tuple. Returns flat array of shape (N*N,).
    Manhattan distance, periodic BC.
    """
    i0, j0 = v0
    max_r = len(G_r) - 1

    ix = np.arange(N)
    di = np.minimum(np.abs(ix - i0), N - np.abs(ix - i0))
    dj = np.minimum(np.abs(ix - j0), N - np.abs(ix - j0))
    DI, DJ = np.meshgrid(di, dj, indexing='ij')
    dist = DI + DJ

    weights = np.zeros((N, N), dtype=np.float64)
    for r in range(max_r + 1):
        mask = (dist == r)
        weights[mask] = np.abs(G_r[r])

    beyond = dist > max_r
    if np.any(beyond):
        weights[beyond] = np.abs(G_r[max_r])

    p = weights.flatten()
    total = np.sum(p)
    if total < 1e-30:
        p = np.ones(N * N) / (N * N)
    else:
        p = p / total

    return p


def compute_FIM_thermal(G_r, N, v0):
    """
    Compute 4×4 Fisher Information Matrix at vertex v0.
    v0 = (i0, j0) tuple. Manhattan distance, periodic BC.
    """
    i0, j0 = v0

    p_v0 = build_thermal_kernel(G_r, N, v0)
    log_p_v0 = np.log(p_v0 + 1e-30)

    neighbors = [
        ((i0 - 1) % N, j0),
        ((i0 + 1) % N, j0),
        (i0, (j0 - 1) % N),
        (i0, (j0 + 1) % N),
    ]
    k = len(neighbors)

    score_vectors = np.zeros((k, N * N))
    for j, w in enumerate(neighbors):
        p_w = build_thermal_kernel(G_r, N, w)
        log_p_w = np.log(p_w + 1e-30)
        score_vectors[j, :] = log_p_w - log_p_v0

    weighted_scores = score_vectors * np.sqrt(p_v0)[np.newaxis, :]
    FIM = weighted_scores @ weighted_scores.T

    return FIM


# =====================================================================
# MODULE 4: FIM Diagnostics
# =====================================================================

def gap_based_rank(sv_norm):
    """Gap-based rank from normalized singular values."""
    if len(sv_norm) <= 1:
        return 1
    ratios = sv_norm[1:] / np.maximum(sv_norm[:-1], 1e-15)
    return int(np.argmin(ratios)) + 1


def participation_ratio(svs):
    """PR = (sum sv_i)² / sum(sv_i²)"""
    s = np.sum(svs)
    s2 = np.sum(svs ** 2)
    if s2 < 1e-30:
        return 0.0
    return s * s / s2


def disorder_index(sv_norm, rank):
    """eta = sv[rank] / sv[rank-1]"""
    if rank >= len(sv_norm) or rank < 1:
        return 0.0
    return sv_norm[rank] / sv_norm[rank - 1] if sv_norm[rank - 1] > 1e-15 else 0.0


def fisher_diagnostics(G_r, N, n_samples=30):
    """
    Sample n_samples random vertices. Compute FIM, rank, PR, eta, sv_profile.
    Returns diagnostics dict.
    """
    ranks = []
    prs = []
    etas = []
    sv_profiles_list = []
    gap_ratios = []

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
        sv_profiles_list.append(sv_n)
        gap_ratios.append(gr)
        sampled += 1

    if not ranks:
        return None

    sv_arr = np.array(sv_profiles_list)

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
    }


# =====================================================================
# MODULE 5: Macroscopic Observables
# =====================================================================

def measure_potts_energy(lattice, J=1.0):
    """Energy per spin for Potts: -J * (same-state neighbor pairs) / N²."""
    N = lattice.shape[0]
    same_right = np.sum(lattice == np.roll(lattice, -1, axis=1))
    same_down = np.sum(lattice == np.roll(lattice, -1, axis=0))
    return -J * (same_right + same_down) / (N * N)


def measure_potts_macroscopic(lattice, T, q=4, n_configs=500,
                               n_sweeps_between=5):
    """
    Measure macroscopic Potts observables.
    Order parameter: m = (q * max_fraction - 1) / (q - 1).
    """
    N = lattice.shape[0]
    N2 = N * N

    ms, m2s, Es, E2s = [], [], [], []

    for _ in range(n_configs):
        for _ in range(n_sweeps_between):
            wolff_step_potts(lattice, T, q)

        counts = np.bincount(lattice.flatten(), minlength=q)
        max_frac = np.max(counts) / N2
        m_val = (q * max_frac - 1.0) / (q - 1.0)
        ms.append(m_val)
        m2s.append(m_val ** 2)

        e_val = measure_potts_energy(lattice)
        Es.append(e_val)
        E2s.append(e_val ** 2)

    m_mean = np.mean(ms)
    E_mean = np.mean(Es)
    chi = N2 * (np.mean(m2s) - m_mean ** 2)

    return {
        'magnetization': float(m_mean),
        'susceptibility': float(chi),
        'energy': float(E_mean),
    }


# =====================================================================
# MODULE 6: Single Temperature Runner
# =====================================================================

def run_temperature(N, T_val, t_idx, n_temps, params):
    """Run full Potts q=4 analysis at one temperature for lattice size N."""
    t0 = time.time()
    T_over_Tc = T_val / T_C_Q4

    mode = 'hot' if T_val > T_C_Q4 else 'cold'
    lattice = initialize_potts_lattice(N, q=Q, mode=mode)
    equilibrate_potts(lattice, T_val, q=Q, n_sweeps=params['n_eq'])

    nc = params['n_configs']
    ns = params['n_fim']
    if 0.95 <= T_over_Tc <= 1.05:
        nc = min(nc * 2, 1000)
        ns = min(ns + 20, 50)

    # Correlations
    lat_corr = lattice.copy()
    G_r = accumulate_potts_correlations(lat_corr, T_val, q=Q,
                                         n_configs=nc, n_sweeps_between=5)
    xi, xi_r2 = estimate_correlation_length(G_r)

    # FIM diagnostics
    fisher = fisher_diagnostics(G_r, N, n_samples=ns)

    # Macroscopic observables
    lat_macro = lattice.copy()
    macro = measure_potts_macroscopic(lat_macro, T_val, q=Q,
                                       n_configs=nc, n_sweeps_between=5)

    dt = time.time() - t0

    sv_raw = fisher['sv_profile_mean'] if fisher else [0] * 4
    sv_profile = sv_raw.tolist() if hasattr(sv_raw, 'tolist') else list(sv_raw)
    sv2sv1_val = sv_profile[1] if len(sv_profile) > 1 else 0.0

    result = {
        'N': N,
        'T': T_val,
        'T_over_Tc': T_over_Tc,
        'sv2sv1_mean': sv2sv1_val,
        'sv2sv1_std': fisher['sv_profile_std'][1] if fisher and len(fisher['sv_profile_std']) > 1 else 0.0,
        'rank_mean': fisher['rank_mean'] if fisher else 0.0,
        'eta_mean': fisher['eta_mean'] if fisher else 0.0,
        'pr_mean': fisher['pr_mean'] if fisher else 0.0,
        'xi': xi if np.isfinite(xi) else -1.0,
        'chi': macro['susceptibility'],
        'energy': macro['energy'],
        'sv_profile': sv_profile,
        'time_sec': dt,
    }

    print(f"  [{t_idx+1}/{n_temps}] N={N} T={T_val:.4f} (T/Tc={T_over_Tc:.3f}): "
          f"SV₂/SV₁={result['sv2sv1_mean']:.4f} rank={result['rank_mean']:.2f} "
          f"η={result['eta_mean']:.4f} PR={result['pr_mean']:.3f} "
          f"ξ={xi:.1f} χ={macro['susceptibility']:.1f} [{dt:.1f}s]",
          flush=True)

    return result


# =====================================================================
# MODULE 7: Finite-Size Scaling Analysis
# =====================================================================

def analyze_scaling(all_results):
    """
    Compare SV₂/SV₁ at T_c across system sizes.
    Returns analysis text.
    """
    lines = []
    lines.append("=" * 70)
    lines.append("PHASE 4B-FSS: FINITE-SIZE SCALING ANALYSIS")
    lines.append("=" * 70)
    lines.append(f"q = {Q}, T_c = {T_C_Q4:.6f}")
    lines.append(f"Lattice sizes: {LATTICE_SIZES}")
    lines.append(f"Timestamp: {datetime.datetime.now().isoformat()}")
    lines.append("")

    # Extract SV₂/SV₁ at T_c for each N
    sv_at_tc = {}
    for N_val in LATTICE_SIZES:
        n_results = [r for r in all_results if r['N'] == N_val]
        if not n_results:
            continue
        T_ratios = np.array([r['T_over_Tc'] for r in n_results])
        tc_idx = np.argmin(np.abs(T_ratios - 1.0))
        sv_at_tc[N_val] = {
            'sv2sv1': n_results[tc_idx]['sv2sv1_mean'],
            'sv2sv1_std': n_results[tc_idx]['sv2sv1_std'],
            'T_over_Tc': n_results[tc_idx]['T_over_Tc'],
            'rank': n_results[tc_idx]['rank_mean'],
            'pr': n_results[tc_idx]['pr_mean'],
            'xi': n_results[tc_idx]['xi'],
        }

    lines.append("SV₂/SV₁ at T_c by system size:")
    lines.append("-" * 50)
    lines.append(f"{'N':>6} {'SV₂/SV₁':>10} {'±':>8} {'T/Tc':>7} {'rank':>6} {'PR':>7} {'ξ':>8}")
    for N_val in sorted(sv_at_tc.keys()):
        d = sv_at_tc[N_val]
        lines.append(f"{N_val:6d} {d['sv2sv1']:10.4f} {d['sv2sv1_std']:8.4f} "
                     f"{d['T_over_Tc']:7.3f} {d['rank']:6.2f} {d['pr']:7.3f} "
                     f"{d['xi']:8.1f}")
    lines.append("")

    # Scaling trend
    sizes = sorted(sv_at_tc.keys())
    if len(sizes) >= 2:
        vals = [sv_at_tc[n]['sv2sv1'] for n in sizes]
        increasing = all(vals[i] <= vals[i+1] for i in range(len(vals)-1))
        decreasing = all(vals[i] >= vals[i+1] for i in range(len(vals)-1))

        lines.append("SCALING TREND:")
        for i in range(1, len(sizes)):
            delta = vals[i] - vals[i-1]
            pct = 100 * delta / vals[i-1] if vals[i-1] > 0 else 0
            lines.append(f"  N={sizes[i-1]}→{sizes[i]}: "
                        f"SV₂/SV₁ {vals[i-1]:.4f} → {vals[i]:.4f} "
                        f"(Δ={delta:+.4f}, {pct:+.1f}%)")
        lines.append("")

        if increasing:
            lines.append("VERDICT: SV₂/SV₁ INCREASES with N → logarithmic corrections confirmed")
            lines.append("  q=4 transition is continuous (as theory predicts)")
            lines.append("  N=128 value of 0.521 was suppressed by finite-size log corrections")
        elif decreasing:
            lines.append("VERDICT: SV₂/SV₁ DECREASES with N → ANOMALOUS")
            lines.append("  This would contradict theory — needs investigation")
        else:
            lines.append("VERDICT: NON-MONOTONIC trend — inconclusive")
            lines.append("  May need larger system sizes or more statistics")

    lines.append("")

    # Compare with Phase 2 reference spectrum
    lines.append("REFERENCE SPECTRUM (Phase 2, N=128):")
    lines.append(f"  q=2 (Ising):  1.000  (continuous)")
    for N_val in sorted(sv_at_tc.keys()):
        d = sv_at_tc[N_val]
        lines.append(f"  q=4 (N={N_val}):  {d['sv2sv1']:.4f}  ← THIS STUDY")
    lines.append(f"  q=5:          0.691  (weak first-order)")
    lines.append(f"  q=10:         0.374  (strong first-order)")
    lines.append("")

    # Key question: does q=4 rise above q=5 at large N?
    if sizes:
        largest_N = max(sizes)
        largest_sv = sv_at_tc[largest_N]['sv2sv1']
        if largest_sv > 0.691:
            lines.append(f"✓ At N={largest_N}, q=4 SV₂/SV₁ ({largest_sv:.4f}) > q=5 (0.691)")
            lines.append("  → q-ordering restored, transition is continuous")
        elif largest_sv > 0.60:
            lines.append(f"~ At N={largest_N}, q=4 SV₂/SV₁ ({largest_sv:.4f}) approaching q=5 (0.691)")
            lines.append("  → Trend suggests continuous, may need N=1024 to confirm")
        else:
            lines.append(f"✗ At N={largest_N}, q=4 SV₂/SV₁ ({largest_sv:.4f}) still < q=5 (0.691)")
            lines.append("  → Log corrections are severe, or genuine anomaly")

    lines.append("")
    lines.append("=" * 70)

    return "\n".join(lines)


def make_scaling_plot(all_results):
    """Generate finite-size scaling plots."""
    if not HAS_MATPLOTLIB:
        print("  [matplotlib not available, skipping plots]", flush=True)
        return

    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    fig.suptitle(f"Phase 4B-FSS: q=4 Potts Finite-Size Scaling (T_c={T_C_Q4:.4f})",
                fontsize=14, fontweight='bold')

    colors = {128: '#2196F3', 256: '#4CAF50', 512: '#FF9800', 1024: '#F44336'}

    # A: SV₂/SV₁ vs T/T_c for each N
    ax = axes[0, 0]
    for N_val in LATTICE_SIZES:
        n_results = sorted([r for r in all_results if r['N'] == N_val],
                          key=lambda r: r['T'])
        if not n_results:
            continue
        Ts = [r['T_over_Tc'] for r in n_results]
        svs = [r['sv2sv1_mean'] for r in n_results]
        errs = [r['sv2sv1_std'] for r in n_results]
        c = colors.get(N_val, '#795548')
        ax.errorbar(Ts, svs, yerr=errs, fmt='o-', color=c, capsize=3,
                   label=f'N={N_val}', markersize=4)
    ax.axvline(1.0, color='red', linestyle='--', alpha=0.7, label='T_c')
    ax.axhline(0.691, color='orange', linestyle=':', alpha=0.5, label='q=5 ref')
    ax.axhline(1.000, color='green', linestyle=':', alpha=0.5, label='q=2 ref')
    ax.set_xlabel('T / T_c'); ax.set_ylabel('SV₂/SV₁')
    ax.set_title('FIM Spectral Ratio vs Temperature')
    ax.legend(fontsize=8); ax.grid(True, alpha=0.3)

    # B: SV₂/SV₁ at T_c vs N (the key scaling plot)
    ax = axes[0, 1]
    Ns_plot = []
    sv_tc_plot = []
    sv_tc_err = []
    for N_val in LATTICE_SIZES:
        n_results = [r for r in all_results if r['N'] == N_val]
        if not n_results:
            continue
        T_ratios = np.array([r['T_over_Tc'] for r in n_results])
        tc_idx = np.argmin(np.abs(T_ratios - 1.0))
        Ns_plot.append(N_val)
        sv_tc_plot.append(n_results[tc_idx]['sv2sv1_mean'])
        sv_tc_err.append(n_results[tc_idx]['sv2sv1_std'])

    if Ns_plot:
        ax.errorbar(Ns_plot, sv_tc_plot, yerr=sv_tc_err, fmt='s-',
                    color='#2196F3', capsize=5, markersize=8, linewidth=2)
        ax.axhline(0.691, color='orange', linestyle=':', alpha=0.5, label='q=5 (N=128)')
        ax.axhline(1.000, color='green', linestyle=':', alpha=0.5, label='q=2 (N=128)')
        ax.set_xlabel('N (lattice size)'); ax.set_ylabel('SV₂/SV₁ at T_c')
        ax.set_title('Finite-Size Scaling of SV₂/SV₁')
        ax.set_xscale('log', base=2)
        ax.set_xticks(Ns_plot)
        ax.set_xticklabels([str(n) for n in Ns_plot])
        ax.legend(fontsize=8); ax.grid(True, alpha=0.3)

    # C: Correlation length at T_c vs N
    ax = axes[1, 0]
    for N_val in LATTICE_SIZES:
        n_results = sorted([r for r in all_results if r['N'] == N_val],
                          key=lambda r: r['T'])
        if not n_results:
            continue
        Ts = [r['T_over_Tc'] for r in n_results]
        xis = [r['xi'] for r in n_results]
        valid = [(t, x) for t, x in zip(Ts, xis) if 0 < x < N_val * 10]
        if valid:
            c = colors.get(N_val, '#795548')
            ax.plot([t for t, x in valid], [x for t, x in valid], 'o-',
                   color=c, label=f'N={N_val}', markersize=4)
    ax.axvline(1.0, color='red', linestyle='--', alpha=0.7)
    ax.set_xlabel('T / T_c'); ax.set_ylabel('ξ')
    ax.set_title('Correlation Length')
    ax.set_yscale('log')
    ax.legend(fontsize=8); ax.grid(True, alpha=0.3)

    # D: Susceptibility
    ax = axes[1, 1]
    for N_val in LATTICE_SIZES:
        n_results = sorted([r for r in all_results if r['N'] == N_val],
                          key=lambda r: r['T'])
        if not n_results:
            continue
        Ts = [r['T_over_Tc'] for r in n_results]
        chis = [r['chi'] for r in n_results]
        c = colors.get(N_val, '#795548')
        ax.plot(Ts, chis, 'o-', color=c, label=f'N={N_val}', markersize=4)
    ax.axvline(1.0, color='red', linestyle='--', alpha=0.7)
    ax.set_xlabel('T / T_c'); ax.set_ylabel('χ')
    ax.set_title('Susceptibility')
    ax.legend(fontsize=8); ax.grid(True, alpha=0.3)

    plt.tight_layout()
    path = os.path.join(RESULTS_DIR, "phase4b_fss_scaling.png")
    plt.savefig(path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  Plot saved: {path}", flush=True)


# =====================================================================
# MAIN
# =====================================================================

def main():
    os.makedirs(RESULTS_DIR, exist_ok=True)

    # Parse command-line: optionally specify which sizes to run
    # Usage: python potts4_finite_size_scaling.py [128] [256] [512]
    if len(sys.argv) > 1:
        sizes = [int(x) for x in sys.argv[1:] if x.isdigit()]
        if sizes:
            global LATTICE_SIZES
            LATTICE_SIZES = sizes

    print("=" * 70, flush=True)
    print("PHASE 4B-FSS: q=4 POTTS FINITE-SIZE SCALING STUDY", flush=True)
    print("=" * 70, flush=True)
    print(f"q = {Q}, T_c = {T_C_Q4:.6f}", flush=True)
    print(f"Lattice sizes: {LATTICE_SIZES}", flush=True)
    print(f"T sweep: {len(T_VALUES)} values, T/T_c from {T_OVER_TC_VALUES[0]:.2f} to {T_OVER_TC_VALUES[-1]:.2f}", flush=True)
    print(f"Start: {datetime.datetime.now().isoformat()}", flush=True)
    print(flush=True)

    print("Reference values (Phase 2, N=128):", flush=True)
    for label, sv in REFERENCE.items():
        print(f"  {label}: SV₂/SV₁ = {sv:.3f}", flush=True)
    print(flush=True)

    all_results = []
    t_global_start = time.time()

    for N_val in LATTICE_SIZES:
        params = MC_PARAMS.get(N_val, {'n_eq': 500, 'n_configs': 400, 'n_fim': 25})
        n_temps = len(T_VALUES)

        print(f"\n{'─' * 60}", flush=True)
        print(f"  N = {N_val} ({N_val}×{N_val} = {N_val**2:,} spins)", flush=True)
        print(f"  MC: {params['n_eq']} eq + {params['n_configs']} configs, "
              f"{params['n_fim']} FIM samples", flush=True)
        print(f"{'─' * 60}", flush=True)

        t_n_start = time.time()

        for i, T_val in enumerate(T_VALUES):
            np.random.seed(42 + N_val + i)  # Reproducible per (N, T)
            res = run_temperature(N_val, T_val, i, n_temps, params)
            all_results.append(res)

            # Save intermediate results after each temperature
            # (safety: if script crashes, partial results are preserved)
            interim_path = os.path.join(RESULTS_DIR,
                                         f"phase4b_fss_N{N_val}_interim.csv")
            n_results = [r for r in all_results if r['N'] == N_val]
            with open(interim_path, 'w', newline='') as f:
                writer = csv.DictWriter(f, fieldnames=[
                    'N', 'T', 'T_over_Tc', 'sv2sv1_mean', 'sv2sv1_std',
                    'rank_mean', 'eta_mean', 'pr_mean', 'xi', 'chi',
                    'energy', 'time_sec'
                ])
                writer.writeheader()
                for r in n_results:
                    row = {k: f"{v:.6f}" if isinstance(v, float) else v
                           for k, v in r.items() if k != 'sv_profile'}
                    writer.writerow(row)

        t_n_total = time.time() - t_n_start
        print(f"\n  N={N_val} complete: {t_n_total:.1f}s ({t_n_total/60:.1f} min)", flush=True)

    t_global = time.time() - t_global_start
    print(f"\n{'=' * 70}", flush=True)
    print(f"ALL SIZES COMPLETE: {t_global:.1f}s ({t_global/60:.1f} min)", flush=True)
    print(f"{'=' * 70}", flush=True)

    # Write combined CSV
    csv_path = os.path.join(RESULTS_DIR, "phase4b_fss_all_results.csv")
    with open(csv_path, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=[
            'N', 'T', 'T_over_Tc', 'sv2sv1_mean', 'sv2sv1_std',
            'rank_mean', 'eta_mean', 'pr_mean', 'xi', 'chi',
            'energy', 'time_sec'
        ])
        writer.writeheader()
        for r in all_results:
            row = {k: f"{v:.6f}" if isinstance(v, float) else v
                   for k, v in r.items() if k != 'sv_profile'}
            writer.writerow(row)
    print(f"Results CSV: {csv_path}", flush=True)

    # Scaling analysis
    analysis = analyze_scaling(all_results)
    analysis_path = os.path.join(RESULTS_DIR, "phase4b_fss_analysis.txt")
    with open(analysis_path, 'w') as f:
        f.write(analysis)
    print(f"Analysis: {analysis_path}", flush=True)

    # Plots
    make_scaling_plot(all_results)

    # Save JSON summary for programmatic access
    summary = {}
    for N_val in LATTICE_SIZES:
        n_results = [r for r in all_results if r['N'] == N_val]
        if not n_results:
            continue
        T_ratios = np.array([r['T_over_Tc'] for r in n_results])
        tc_idx = np.argmin(np.abs(T_ratios - 1.0))
        summary[str(N_val)] = {
            'sv2sv1_at_tc': n_results[tc_idx]['sv2sv1_mean'],
            'sv2sv1_std': n_results[tc_idx]['sv2sv1_std'],
            'T_over_Tc': n_results[tc_idx]['T_over_Tc'],
            'rank': n_results[tc_idx]['rank_mean'],
            'pr': n_results[tc_idx]['pr_mean'],
            'total_time_sec': sum(r['time_sec'] for r in n_results),
        }

    json_path = os.path.join(RESULTS_DIR, "phase4b_fss_summary.json")
    with open(json_path, 'w') as f:
        json.dump({
            'q': Q,
            'T_c': T_C_Q4,
            'lattice_sizes': LATTICE_SIZES,
            'reference': REFERENCE,
            'results': summary,
            'total_runtime_sec': t_global,
            'timestamp': datetime.datetime.now().isoformat(),
        }, f, indent=2)
    print(f"JSON summary: {json_path}", flush=True)

    # Print analysis to stdout
    print(flush=True)
    print(analysis, flush=True)


if __name__ == "__main__":
    main()
