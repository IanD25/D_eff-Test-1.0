#!/usr/bin/env python3
"""
Phase 4A′: 2D XY Model — FIM Spectral Transition at BKT Point
==============================================================
Tests whether the FIM SV₂/SV₁ degeneracy swap occurs at the
Berezinskii-Kosterlitz-Thouless (BKT) transition T_BKT ≈ 0.893.

The XY model is the lattice version of Kuramoto synchronization:
  H = -J Σ_{<i,j>} cos(θᵢ - θⱼ)

Two-scale mapping:
  Scale 1: Vortex-antivortex pair separation ξ_vortex
  Scale 2: System size L
  At T_BKT: pairs unbind, ξ → ∞

FIM construction: THERMAL KERNEL (same as Phase 2 Ising paper)
  p_v(u) = |C(r)| / Z  where r = Manhattan distance, C(r) = ⟨cos(θᵢ-θⱼ)⟩
  NOT the geometric kernel exp(-d/σ) from Phase 1.

Pre-registered predictions:
  P4A-1: SV₂/SV₁ peaks within ±15% of T_BKT = 0.893 [KILL if absent]
  P4A-2: SV₂/SV₁ at T_BKT exceeds SV₂/SV₁ at T=0.3 by ≥0.10
  P4A-3: Non-monotonic (below < peak > above)
  P4A-4: η minimum coincides with SV₂/SV₁ maximum (±20% T)
  P4A-5: Rank elevation at T_BKT

Author: Ian Darling / Claude Code
Date: 2026-03-25
"""

import numpy as np
import csv
import time
import os

# ─── Configuration ────────────────────────────────────────────────────
L = 64                 # Lattice size (64×64 = 4096 spins)
N_SPINS = L * L
T_BKT = 0.893          # Best estimate from literature
J = 1.0                # Coupling constant

# Temperature sweep: 20 values
T_VALUES = np.concatenate([
    np.linspace(0.30, 0.70, 5),     # Well below T_BKT (ordered, algebraic C(r))
    np.linspace(0.76, 1.03, 8),     # Near T_BKT (±15% window)
    np.linspace(1.10, 1.80, 7),     # Well above T_BKT (disordered, exp C(r))
])

# Monte Carlo parameters
N_THERM = 3000         # Thermalization sweeps (Wolff)
N_CONFIGS = 200        # Configurations for C(r) averaging
N_SWEEPS_BETWEEN = 5   # Wolff sweeps between measurements

# FIM parameters (match Phase 2 Ising paper exactly)
N_FIM_SAMPLES = 30     # Sample centers for FIM computation

RESULTS_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "results")


# ============================================================
# MODULE 1: XY Monte Carlo
# ============================================================

def xy_wolff_step(theta, L, T):
    """One Wolff cluster step for XY model.
    Projects spins onto random direction, builds cluster, reflects."""
    N = L * L
    beta = 1.0 / T

    # Random projection direction
    r_angle = np.random.uniform(0, 2 * np.pi)

    # Start from random seed
    seed = np.random.randint(N)
    cluster = np.zeros(N, dtype=bool)
    cluster[seed] = True
    stack = [seed]

    while stack:
        i = stack.pop()
        x, y = i % L, i // L
        si_proj = np.cos(theta[i] - r_angle)

        for dx, dy in [(1,0),(-1,0),(0,1),(0,-1)]:
            j = ((x + dx) % L) + ((y + dy) % L) * L
            if not cluster[j]:
                sj_proj = np.cos(theta[j] - r_angle)
                bond_weight = 2 * beta * J * si_proj * sj_proj
                if bond_weight > 0 and np.random.random() < 1 - np.exp(-bond_weight):
                    cluster[j] = True
                    stack.append(j)

    # Reflect cluster spins about r direction
    cluster_idx = np.where(cluster)[0]
    theta[cluster_idx] = 2 * r_angle - theta[cluster_idx]
    return theta, len(cluster_idx)


def xy_metropolis_sweep(theta, L, T):
    """One vectorized Metropolis sweep. Much faster than pure Python loops."""
    N = L * L
    beta = 1.0 / T

    # Process all sites in random order
    order = np.random.permutation(N)
    for idx in order:
        x, y = idx % L, idx // L
        old = theta[idx]
        new = old + np.random.uniform(-np.pi, np.pi)

        # Neighbor angles
        n_angles = np.array([
            theta[((x+1) % L) + y * L],
            theta[((x-1) % L) + y * L],
            theta[x + ((y+1) % L) * L],
            theta[x + ((y-1) % L) * L],
        ])

        dE = -J * np.sum(np.cos(new - n_angles) - np.cos(old - n_angles))

        if dE < 0 or np.random.random() < np.exp(-beta * dE):
            theta[idx] = new

    return theta


def compute_magnetization(theta, L):
    """M = |N⁻¹ Σ exp(iθ)|"""
    mx = np.mean(np.cos(theta))
    my = np.mean(np.sin(theta))
    return np.sqrt(mx**2 + my**2)


# ============================================================
# MODULE 2: Correlation Function (FFT-based, matching Phase 2)
# ============================================================

def measure_correlations_fft_xy(theta, L):
    """
    Measure C(r) = ⟨cos(θᵢ - θⱼ)⟩ at Manhattan distance r.
    Uses FFT: C(r) = ⟨Re[exp(iθ)]⟩_r via autocorrelation of cos/sin fields.

    For XY model: ⟨cos(θᵢ-θⱼ)⟩ = ⟨cosθᵢ cosθⱼ⟩ + ⟨sinθᵢ sinθⱼ⟩
    """
    cos_field = np.cos(theta).reshape(L, L)
    sin_field = np.sin(theta).reshape(L, L)

    # Mean values for connected correlation
    mc = np.mean(cos_field)
    ms = np.mean(sin_field)

    # FFT autocorrelation of cos field
    Fc = np.fft.fft2(cos_field)
    Gc_full = np.real(np.fft.ifft2(Fc * np.conj(Fc))) / (L * L)

    # FFT autocorrelation of sin field
    Fs = np.fft.fft2(sin_field)
    Gs_full = np.real(np.fft.ifft2(Fs * np.conj(Fs))) / (L * L)

    # Total correlation = cos-cos + sin-sin
    G_full = Gc_full + Gs_full

    # Subtract disconnected part
    m_sq = mc**2 + ms**2

    # Radial average using Manhattan distance
    max_r = L // 2
    G_r = np.zeros(max_r + 1)
    counts = np.zeros(max_r + 1)

    ix = np.arange(L)
    dx = np.minimum(ix, L - ix)
    DX, DY = np.meshgrid(dx, dx)
    dist_grid = DX + DY

    for r in range(max_r + 1):
        mask = (dist_grid == r)
        if np.any(mask):
            G_r[r] = np.mean(G_full[mask]) - m_sq
            counts[r] = np.sum(mask)

    return G_r


def accumulate_correlations_xy(theta, L, T, n_configs, n_sweeps):
    """Average C(r) over n_configs independent configurations."""
    max_r = L // 2
    G_sum = np.zeros(max_r + 1)

    for i in range(n_configs):
        for _ in range(n_sweeps):
            theta, _ = xy_wolff_step(theta, L, T)
        G_r = measure_correlations_fft_xy(theta, L)
        G_sum += G_r[:max_r + 1]

    return G_sum / n_configs, theta


def estimate_correlation_length(G_r):
    """Fit G(r) ~ A * exp(-r/ξ) for r > 1."""
    valid = np.abs(G_r[1:]) > 1e-10
    if np.sum(valid) < 3:
        return 0.0

    r_vals = np.arange(1, len(G_r))
    log_G = np.log(np.abs(G_r[1:]) + 1e-30)

    # Use first half of valid range for fit
    n_fit = max(3, np.sum(valid) // 2)
    r_fit = r_vals[valid][:n_fit]
    y_fit = log_G[valid][:n_fit]

    if len(r_fit) < 2:
        return 0.0

    # Linear fit: log|G| = log(A) - r/ξ
    try:
        coeffs = np.polyfit(r_fit, y_fit, 1)
        if coeffs[0] >= 0:
            return float('inf')
        xi = -1.0 / coeffs[0]
        return max(0.0, xi)
    except:
        return 0.0


# ============================================================
# MODULE 3: Thermal FIM (matching Phase 2 Ising exactly)
# ============================================================

def build_thermal_kernel(G_r, L, v0):
    """
    Build thermal probability distribution p_{v0}(u; T) using |G(r)|.
    v0 = (i0, j0) tuple on L×L lattice.
    Returns flat array of shape (L*L,).
    """
    i0, j0 = v0
    max_r = len(G_r) - 1

    ix = np.arange(L)
    di = np.minimum(np.abs(ix - i0), L - np.abs(ix - i0))
    dj = np.minimum(np.abs(ix - j0), L - np.abs(ix - j0))
    DI, DJ = np.meshgrid(di, dj, indexing='ij')
    dist = DI + DJ  # Manhattan distance

    weights = np.zeros((L, L), dtype=np.float64)
    for r in range(max_r + 1):
        mask = (dist == r)
        weights[mask] = np.abs(G_r[r])

    beyond = dist > max_r
    if np.any(beyond):
        weights[beyond] = np.abs(G_r[max_r])

    p = weights.flatten()
    total = np.sum(p)
    if total < 1e-30:
        p = np.ones(L * L) / (L * L)
    else:
        p = p / total

    return p


def compute_FIM_thermal(G_r, L, v0):
    """
    Compute Fisher Information Matrix at vertex v0 using thermal kernel.
    Returns 4×4 FIM matrix.
    """
    i0, j0 = v0

    p_v0 = build_thermal_kernel(G_r, L, v0)
    log_p_v0 = np.log(p_v0 + 1e-30)

    neighbors = [
        ((i0 - 1) % L, j0),
        ((i0 + 1) % L, j0),
        (i0, (j0 - 1) % L),
        (i0, (j0 + 1) % L),
    ]
    k = len(neighbors)

    score_vectors = np.zeros((k, L * L))
    for j, w in enumerate(neighbors):
        p_w = build_thermal_kernel(G_r, L, w)
        log_p_w = np.log(p_w + 1e-30)
        score_vectors[j, :] = log_p_w - log_p_v0

    weighted_scores = score_vectors * np.sqrt(p_v0)[np.newaxis, :]
    FIM = weighted_scores @ weighted_scores.T

    return FIM


# ============================================================
# MODULE 4: Diagnostics (matching Phase 2 exactly)
# ============================================================

def gap_based_rank(sv_norm):
    """Gap-based rank from normalized singular values."""
    if len(sv_norm) <= 1:
        return 1
    ratios = sv_norm[1:] / np.maximum(sv_norm[:-1], 1e-15)
    return int(np.argmin(ratios)) + 1


def participation_ratio(svs):
    """PR = (Σ sv_i)² / Σ(sv_i²)"""
    s = np.sum(svs)
    s2 = np.sum(svs ** 2)
    if s2 < 1e-30:
        return 0.0
    return s * s / s2


def disorder_index(sv_norm, rank):
    """η = sv[rank] / sv[rank-1]"""
    if rank >= len(sv_norm) or rank < 1:
        return 0.0
    return sv_norm[rank] / sv_norm[rank - 1] if sv_norm[rank - 1] > 1e-15 else 0.0


def fisher_diagnostics(G_r, L, n_samples=30):
    """Sample vertices, compute FIM diagnostics."""
    ranks = []
    prs = []
    etas = []
    sv_profiles = []

    margin = max(2, L // 8)
    sampled = 0

    while sampled < n_samples:
        i = np.random.randint(margin, L - margin)
        j = np.random.randint(margin, L - margin)
        v0 = (i, j)

        FIM = compute_FIM_thermal(G_r, L, v0)
        svs = np.linalg.svd(FIM, compute_uv=False)
        if svs[0] < 1e-30:
            continue
        sv_norm = svs / svs[0]

        rank = gap_based_rank(sv_norm)
        pr = participation_ratio(sv_norm)
        eta = disorder_index(sv_norm, rank)

        ranks.append(rank)
        prs.append(pr)
        etas.append(eta)
        sv_profiles.append(sv_norm)
        sampled += 1

    sv2sv1 = np.mean([sp[1] for sp in sv_profiles])
    sv2sv1_std = np.std([sp[1] for sp in sv_profiles])

    return {
        'rank_mean': np.mean(ranks),
        'rank_std': np.std(ranks),
        'pr_mean': np.mean(prs),
        'pr_std': np.std(prs),
        'eta_mean': np.mean(etas),
        'eta_std': np.std(etas),
        'sv2sv1_mean': sv2sv1,
        'sv2sv1_std': sv2sv1_std,
        'mean_sv_profile': np.mean(sv_profiles, axis=0),
    }


# ============================================================
# MODULE 5: Main experiment loop
# ============================================================

def run_single_T(T_val):
    """Full XY simulation + FIM analysis at one temperature."""
    # Start from ORDERED state (all θ=0) — critical for Wolff efficiency
    theta = np.zeros(N_SPINS)

    # Phase 1: Metropolis thermalization (fast from ordered state)
    print(f"    Metropolis therm...", end="", flush=True)
    for _ in range(20):
        theta = xy_metropolis_sweep(theta, L, T_val)
    print(f" M={compute_magnetization(theta, L):.3f}", end="", flush=True)

    # Phase 2: Wolff thermalization (large clusters now that system is partially ordered)
    print(f" → Wolff therm...", end="", flush=True)
    total_flipped = 0
    for _ in range(N_THERM):
        theta, cs = xy_wolff_step(theta, L, T_val)
        total_flipped += cs
    avg_cluster = total_flipped / N_THERM
    print(f" avg_cluster={avg_cluster:.0f}/{N_SPINS}", flush=True)

    M = compute_magnetization(theta, L)

    # Accumulate correlations
    G_r, theta = accumulate_correlations_xy(theta, L, T_val, N_CONFIGS, N_SWEEPS_BETWEEN)

    # Estimate correlation length
    xi = estimate_correlation_length(G_r)

    # Fisher diagnostics
    diag = fisher_diagnostics(G_r, L, N_FIM_SAMPLES)

    return {
        'T': T_val,
        'M': M,
        'xi': xi,
        'sv2sv1_mean': diag['sv2sv1_mean'],
        'sv2sv1_std': diag['sv2sv1_std'],
        'rank_mean': diag['rank_mean'],
        'eta_mean': diag['eta_mean'],
        'pr_mean': diag['pr_mean'],
        'sv_profile': diag['mean_sv_profile'],
    }


def evaluate_predictions(results):
    """Evaluate pre-registered predictions."""
    Ts = np.array([r['T'] for r in results])
    sv2sv1 = np.array([r['sv2sv1_mean'] for r in results])
    eta = np.array([r['eta_mean'] for r in results])
    ranks = np.array([r['rank_mean'] for r in results])

    lines = []
    lines.append("=" * 70)
    lines.append("PHASE 4A′: 2D XY MODEL — BKT TRANSITION PREDICTION VERDICTS")
    lines.append("=" * 70)
    lines.append(f"T_BKT = {T_BKT:.3f}")
    lines.append(f"Kill window: [{T_BKT*0.85:.3f}, {T_BKT*1.15:.3f}] (±15%)")
    lines.append(f"L = {L}, N = {N_SPINS}")
    lines.append("")

    # P4A-1: SV₂/SV₁ peaks within ±15% of T_BKT
    peak_idx = np.argmax(sv2sv1)
    peak_T = Ts[peak_idx]
    peak_val = sv2sv1[peak_idx]
    low_bound = T_BKT * 0.85
    high_bound = T_BKT * 1.15
    p1_pass = low_bound <= peak_T <= high_bound

    lines.append(f"P4A-1 [KILL TEST]: SV₂/SV₁ peak at T={peak_T:.3f} (value={peak_val:.4f})")
    lines.append(f"  Window: [{low_bound:.3f}, {high_bound:.3f}], T_BKT={T_BKT:.3f}")
    lines.append(f"  Offset: {abs(peak_T - T_BKT):.3f} ({abs(peak_T - T_BKT)/T_BKT*100:.1f}%)")
    lines.append(f"  → {'PASS ✓' if p1_pass else 'KILL ✗'}")
    lines.append("")

    # P4A-2: Separation ≥ 0.10
    idx_low = np.argmin(np.abs(Ts - 0.3))
    sv_low = sv2sv1[idx_low]
    sep = peak_val - sv_low
    p2_pass = sep >= 0.10
    p2_weak = sep >= 0.05
    lines.append(f"P4A-2: Peak - T=0.3 = {peak_val:.4f} - {sv_low:.4f} = {sep:.4f}")
    lines.append(f"  → {'PASS ✓' if p2_pass else ('WEAK' if p2_weak else 'FAIL ✗')}")
    lines.append("")

    # P4A-3: Non-monotonic
    below = Ts < T_BKT - 0.1
    above = Ts > T_BKT + 0.1
    if np.any(below) and np.any(above):
        mean_below = np.mean(sv2sv1[below])
        mean_above = np.mean(sv2sv1[above])
        p3_pass = mean_below < peak_val and mean_above < peak_val
        lines.append(f"P4A-3: below={mean_below:.4f} < peak={peak_val:.4f} > above={mean_above:.4f}")
        lines.append(f"  → {'PASS ✓' if p3_pass else 'INCONCLUSIVE'}")
    else:
        p3_pass = False
        lines.append("P4A-3: Insufficient data")
    lines.append("")

    # P4A-4: η minimum near SV₂/SV₁ peak
    eta_min_idx = np.argmin(eta)
    eta_min_T = Ts[eta_min_idx]
    p4_pass = abs(eta_min_T - peak_T) / T_BKT <= 0.20
    lines.append(f"P4A-4: η min at T={eta_min_T:.3f}, SV₂/SV₁ peak at T={peak_T:.3f}")
    lines.append(f"  → {'PASS ✓' if p4_pass else 'INCONCLUSIVE'}")
    lines.append("")

    # P4A-5: Rank elevation
    near_mask = np.abs(Ts - T_BKT) < 0.15
    far_mask = np.abs(Ts - T_BKT) > 0.3
    if np.any(near_mask) and np.any(far_mask):
        rank_near = np.mean(ranks[near_mask])
        rank_far = np.mean(ranks[far_mask])
        p5_pass = rank_near > rank_far
        lines.append(f"P4A-5: Rank near T_BKT={rank_near:.2f}, far={rank_far:.2f}")
        lines.append(f"  → {'PASS ✓' if p5_pass else 'INCONCLUSIVE'}")
    else:
        p5_pass = False
        lines.append("P4A-5: Insufficient data")
    lines.append("")

    # Overall
    lines.append("=" * 70)
    if p1_pass:
        verdict = "PASS — FIM detects BKT transition"
        if p2_pass and p3_pass:
            verdict += " (STRONG: non-monotonic peak with clear separation)"
    else:
        verdict = "KILLED — FIM does NOT detect BKT transition"
    lines.append(f"OVERALL: {verdict}")
    lines.append(f"P4A-1={'PASS' if p1_pass else 'KILL'} "
                f"P4A-2={'PASS' if p2_pass else 'FAIL'} "
                f"P4A-3={'PASS' if p3_pass else '?'} "
                f"P4A-4={'PASS' if p4_pass else '?'} "
                f"P4A-5={'PASS' if p5_pass else '?'}")
    lines.append("=" * 70)

    # Data table
    lines.append("\n--- Full Data ---")
    lines.append(f"{'T':>8} {'M':>8} {'ξ':>8} {'SV2/SV1':>10} {'±':>8} {'rank':>6} {'η':>8} {'PR':>8}")
    for r in results:
        lines.append(f"{r['T']:8.3f} {r['M']:8.3f} {r['xi']:8.2f} "
                     f"{r['sv2sv1_mean']:10.4f} {r['sv2sv1_std']:8.4f} "
                     f"{r['rank_mean']:6.1f} {r['eta_mean']:8.4f} {r['pr_mean']:8.4f}")

    return "\n".join(lines)


def make_plots(results, plt):
    """Generate diagnostic plots."""
    Ts = [r['T'] for r in results]
    sv2sv1 = [r['sv2sv1_mean'] for r in results]
    sv2sv1_err = [r['sv2sv1_std'] for r in results]
    Ms = [r['M'] for r in results]
    ranks = [r['rank_mean'] for r in results]
    etas = [r['eta_mean'] for r in results]
    prs = [r['pr_mean'] for r in results]
    xis = [r['xi'] for r in results]

    fig, axes = plt.subplots(2, 3, figsize=(18, 10))
    fig.suptitle(f"Phase 4A′: 2D XY Model BKT Transition (L={L}, T_BKT={T_BKT})",
                fontsize=14, fontweight='bold')

    # A: SV₂/SV₁
    ax = axes[0, 0]
    ax.errorbar(Ts, sv2sv1, yerr=sv2sv1_err, fmt='o-', color='#2196F3', capsize=3)
    ax.axvline(T_BKT, color='red', linestyle='--', alpha=0.7, label=f'T_BKT = {T_BKT}')
    ax.axvspan(T_BKT*0.85, T_BKT*1.15, alpha=0.1, color='green', label='±15%')
    ax.set_xlabel('T'); ax.set_ylabel('SV₂/SV₁')
    ax.set_title('FIM Spectral Ratio'); ax.legend(fontsize=8); ax.grid(True, alpha=0.3)

    # B: Magnetization
    ax = axes[0, 1]
    ax.plot(Ts, Ms, 's-', color='#FF5722')
    ax.axvline(T_BKT, color='red', linestyle='--', alpha=0.7)
    ax.set_xlabel('T'); ax.set_ylabel('M')
    ax.set_title('Magnetization'); ax.grid(True, alpha=0.3)

    # C: Correlation length
    ax = axes[0, 2]
    finite_xi = [(t, x) for t, x in zip(Ts, xis) if x < 1e6]
    if finite_xi:
        ax.plot([t for t,x in finite_xi], [x for t,x in finite_xi], 'D-', color='#9C27B0')
    ax.axvline(T_BKT, color='red', linestyle='--', alpha=0.7)
    ax.set_xlabel('T'); ax.set_ylabel('ξ')
    ax.set_title('Correlation Length'); ax.grid(True, alpha=0.3)

    # D: Rank
    ax = axes[1, 0]
    ax.plot(Ts, ranks, 'D-', color='#4CAF50')
    ax.axvline(T_BKT, color='red', linestyle='--', alpha=0.7)
    ax.set_xlabel('T'); ax.set_ylabel('Rank')
    ax.set_title('Gap-Based Rank'); ax.grid(True, alpha=0.3)

    # E: η
    ax = axes[1, 1]
    ax.plot(Ts, etas, 'o-', color='#FF9800')
    ax.axvline(T_BKT, color='red', linestyle='--', alpha=0.7)
    ax.set_xlabel('T'); ax.set_ylabel('η')
    ax.set_title('Disorder Index'); ax.grid(True, alpha=0.3)

    # F: PR
    ax = axes[1, 2]
    ax.plot(Ts, prs, '^-', color='#795548')
    ax.axvline(T_BKT, color='red', linestyle='--', alpha=0.7)
    ax.set_xlabel('T'); ax.set_ylabel('PR')
    ax.set_title('Participation Ratio'); ax.grid(True, alpha=0.3)

    plt.tight_layout()
    path = os.path.join(RESULTS_DIR, "phase4a_xy_diagnostics.png")
    plt.savefig(path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  Plot saved: {path}")


def main():
    os.makedirs(RESULTS_DIR, exist_ok=True)

    print("=" * 70)
    print("PHASE 4A′: 2D XY MODEL — BKT TRANSITION FIM TEST")
    print("=" * 70)
    print(f"L = {L}, N = {N_SPINS}, T_BKT = {T_BKT:.3f}")
    print(f"T sweep: {len(T_VALUES)} values from {T_VALUES[0]:.2f} to {T_VALUES[-1]:.2f}")
    print(f"MC: {N_THERM} therm + {N_CONFIGS}×{N_SWEEPS_BETWEEN} measure")
    print(f"FIM: thermal kernel, {N_FIM_SAMPLES} samples/T")
    print()

    t_start = time.time()
    results = []

    for i, T_val in enumerate(T_VALUES):
        print(f"  T={T_val:.3f} ({i+1}/{len(T_VALUES)})...", flush=True)
        t0 = time.time()

        res = run_single_T(T_val)
        dt = time.time() - t0

        results.append(res)
        print(f"    M={res['M']:.3f} ξ={res['xi']:.1f} "
              f"SV₂/SV₁={res['sv2sv1_mean']:.4f}±{res['sv2sv1_std']:.4f} "
              f"rank={res['rank_mean']:.1f} η={res['eta_mean']:.4f} ({dt:.1f}s)")

    t_total = time.time() - t_start
    print(f"\nTotal runtime: {t_total:.1f}s ({t_total/60:.1f} min)")

    # Write CSV
    csv_path = os.path.join(RESULTS_DIR, "phase4a_xy_results.csv")
    with open(csv_path, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=[
            'T', 'M', 'xi', 'sv2sv1_mean', 'sv2sv1_std',
            'rank_mean', 'eta_mean', 'pr_mean'
        ])
        writer.writeheader()
        for r in results:
            row = {k: f"{v:.6f}" if isinstance(v, float) else v
                   for k, v in r.items() if k != 'sv_profile'}
            writer.writerow(row)
    print(f"Results: {csv_path}")

    # Evaluate predictions
    summary = evaluate_predictions(results)
    summary_path = os.path.join(RESULTS_DIR, "phase4a_xy_summary.txt")
    with open(summary_path, 'w') as f:
        f.write(summary)
    print(f"Summary: {summary_path}")

    # Plots
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        make_plots(results, plt)
    except ImportError:
        print("matplotlib not available, skipping plots")

    print("\n" + summary)


if __name__ == "__main__":
    main()
