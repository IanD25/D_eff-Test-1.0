#!/usr/bin/env python3
"""
Phase 4C: 3D Anderson Metal-Insulator Transition — FIM Test
=============================================================
Tests whether the FIM SV₂/SV₁ detects the Anderson localization transition
in a 3D tight-binding model with random on-site disorder.

The Anderson model:
  H = Σᵢ εᵢ|i⟩⟨i| − t Σ_{⟨i,j⟩} |i⟩⟨j|
  εᵢ ~ Uniform[-W/2, W/2], t = 1, 3D cubic lattice with periodic BC.

Critical disorder: W_c ≈ 16.5 ± 0.5 (well-established numerically).
  W < W_c: metallic (extended wavefunctions, IPR ~ 1/N)
  W > W_c: insulating (localized wavefunctions, IPR ~ O(1))

The two-scale commensurability:
  Scale 1: Localization length ξ_loc (diverges at W_c from above)
  Scale 2: System size L (fixed)
  At W_c: ξ_loc ~ L → commensurability → FIM should detect transition

Correlation kernel: C(r) = ⟨|ψ(i)|·|ψ(j)|⟩ averaged over eigenstates near E=0,
replaces thermal G(r,T) from Phase 2.

Pre-registered predictions:
  P4C-1: SV₂/SV₁ peaks in W ∈ [14.5, 18.5] (W_c ± 12%)     [KILL if absent]
  P4C-2: SV₂/SV₁(W=2) > SV₂/SV₁(W=24)                       [metallic > insulating]
  P4C-3: η minimum aligns with SV₂/SV₁ maximum
  P4C-4: IPR transition midpoint ↔ SV₂/SV₁ peak (±20% W)
  P4C-5: Rank elevation at W_c vs W=2 and W=24

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
from collections import Counter

import scipy.sparse as sp
from scipy.sparse.linalg import eigsh

try:
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    HAS_MATPLOTLIB = True
except ImportError:
    HAS_MATPLOTLIB = False

np.random.seed(42)

# ─── Configuration ────────────────────────────────────────────────────
L = 16                    # Lattice dimension: L×L×L = 4096 sites
N_SITES = L ** 3          # Total sites
T_HOP = 1.0              # Hopping amplitude (energy scale)
W_C = 16.5                # Critical disorder (literature value)

# Disorder sweep: 18 values, concentrated near W_c
W_VALUES = np.array([
    2.0, 4.0, 6.0, 8.0, 10.0,          # Metallic regime (5 points)
    12.0, 13.5, 15.0, 15.5,             # Approaching W_c (4 points)
    16.0, 16.5, 17.0, 17.5,             # At/near W_c (4 points)
    18.5, 20.0, 21.5, 23.0, 24.0,       # Insulating regime (5 points)
])

# Eigenstates: k=50 near E=0 (band center, transition is sharpest)
K_EIG = 50

# Disorder realizations per W value
N_REALIZATIONS = 10

# FIM sampling
N_FIM_SAMPLES = 30        # Random sites for FIM computation
N_CORR_SOURCES = 200      # Random source sites for C(r) measurement

RESULTS_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "results")


# =====================================================================
# MODULE 1: Anderson Hamiltonian
# =====================================================================

def site_index(x, y, z, L):
    """Map 3D coordinates to flat index (periodic BC)."""
    return (x % L) * L * L + (y % L) * L + (z % L)


def build_anderson_hamiltonian(L, W, t=1.0, seed=None):
    """
    Build sparse Anderson Hamiltonian on L×L×L cubic lattice.
    H = Σᵢ εᵢ|i⟩⟨i| − t Σ_{⟨i,j⟩} |i⟩⟨j|
    εᵢ ~ Uniform[-W/2, W/2], periodic boundary conditions.

    Returns scipy.sparse CSR matrix (N×N where N=L³).
    """
    if seed is not None:
        np.random.seed(seed)

    N = L ** 3
    H = sp.lil_matrix((N, N), dtype=np.float64)

    # On-site disorder
    disorder = np.random.uniform(-W / 2.0, W / 2.0, size=N)
    for i in range(N):
        H[i, i] = disorder[i]

    # Nearest-neighbor hopping (6 neighbors in 3D cubic, periodic BC)
    for x in range(L):
        for y in range(L):
            for z in range(L):
                i = site_index(x, y, z, L)
                for dx, dy, dz in [(1,0,0), (-1,0,0),
                                    (0,1,0), (0,-1,0),
                                    (0,0,1), (0,0,-1)]:
                    j = site_index(x+dx, y+dy, z+dz, L)
                    H[i, j] = -t

    return H.tocsr()


# =====================================================================
# MODULE 2: Eigenstate Analysis
# =====================================================================

def compute_eigenstates(H, k=50):
    """
    Compute k eigenstates nearest to E=0 using shift-invert.
    Returns (eigenvalues, eigenvectors) where eigvecs[:,i] is the i-th state.
    """
    try:
        eigenvalues, eigenvectors = eigsh(H, k=k, sigma=0.0, which='LM')
        # Sort by |E|
        order = np.argsort(np.abs(eigenvalues))
        return eigenvalues[order], eigenvectors[:, order]
    except Exception as e:
        print(f"  [eigsh failed: {e}, trying smaller k]", flush=True)
        try:
            eigenvalues, eigenvectors = eigsh(H, k=k//2, sigma=0.0, which='LM')
            order = np.argsort(np.abs(eigenvalues))
            return eigenvalues[order], eigenvectors[:, order]
        except Exception as e2:
            print(f"  [eigsh failed again: {e2}]", flush=True)
            return np.array([]), np.array([])


def compute_ipr(eigenvectors):
    """
    Inverse Participation Ratio: IPR = Σᵢ |ψ(i)|⁴
    Returns array of IPR values, one per eigenstate.
    Extended state: IPR ~ 1/N. Localized state: IPR ~ O(1).
    """
    probs = np.abs(eigenvectors) ** 2  # |ψ(i)|² for each state
    ipr = np.sum(probs ** 2, axis=0)   # Σᵢ |ψ(i)|⁴
    return ipr


# =====================================================================
# MODULE 3: Wavefunction Correlation Function C(r)
# =====================================================================

def measure_wavefunction_correlations(eigenvectors, L, n_sources=200, max_r=None):
    """
    Measure C(r) = ⟨|ψ(i)|·|ψ(j)|⟩ as a function of Manhattan distance r.

    Averaged over eigenstates and random source sites.
    This replaces the thermal G(r,T) from the Potts/Ising tests.

    For each eigenstate:
      amplitudes[i] = |ψ(i)|
      For each source site i, compute |ψ(i)|·|ψ(j)| for all j,
      bin by Manhattan distance r = |Δx|+|Δy|+|Δz|.
    """
    N = L ** 3
    if max_r is None:
        max_r = (3 * L) // 4  # Manhattan distance can go up to 3*L/2

    n_states = eigenvectors.shape[1]

    # Precompute 3D coordinates for all sites
    coords = np.zeros((N, 3), dtype=np.int32)
    for i in range(N):
        coords[i, 0] = i // (L * L)
        coords[i, 1] = (i // L) % L
        coords[i, 2] = i % L

    # Accumulate C(r)
    C_sum = np.zeros(max_r + 1, dtype=np.float64)
    C_count = np.zeros(max_r + 1, dtype=np.int64)

    for s in range(n_states):
        amplitudes = np.abs(eigenvectors[:, s])

        # Random source sites
        sources = np.random.choice(N, size=min(n_sources, N), replace=False)

        for src in sources:
            sx, sy, sz = coords[src]
            amp_src = amplitudes[src]

            # Manhattan distances from source to all sites (periodic BC)
            dx = np.minimum(np.abs(coords[:, 0] - sx), L - np.abs(coords[:, 0] - sx))
            dy = np.minimum(np.abs(coords[:, 1] - sy), L - np.abs(coords[:, 1] - sy))
            dz = np.minimum(np.abs(coords[:, 2] - sz), L - np.abs(coords[:, 2] - sz))
            dists = dx + dy + dz

            products = amp_src * amplitudes

            for r in range(max_r + 1):
                mask = (dists == r)
                if np.any(mask):
                    C_sum[r] += np.sum(products[mask])
                    C_count[r] += np.sum(mask)

    # Average
    C_r = np.zeros(max_r + 1)
    for r in range(max_r + 1):
        if C_count[r] > 0:
            C_r[r] = C_sum[r] / C_count[r]

    return C_r


def measure_wavefunction_correlations_fast(eigenvectors, L, max_r=None):
    """
    Fast C(r) using FFT autocorrelation of |ψ| on 3D lattice.
    Analogous to the FFT correlation measurement in Phase 2 Ising/Potts.

    For each eigenstate, compute autocorrelation of |ψ(r)| field,
    then radially average using Manhattan distance.
    """
    N = L ** 3
    if max_r is None:
        max_r = (3 * L) // 4

    n_states = eigenvectors.shape[1]

    # Precompute Manhattan distance grid
    ix = np.arange(L)
    dx = np.minimum(ix, L - ix)
    DX, DY, DZ = np.meshgrid(dx, dx, dx, indexing='ij')
    dist_grid = DX + DY + DZ

    C_r_total = np.zeros(max_r + 1)

    for s in range(n_states):
        # Reshape to 3D
        amp = np.abs(eigenvectors[:, s]).reshape((L, L, L))

        # FFT autocorrelation: C(Δr) = IFFT(|FFT(amp)|²) / N
        F = np.fft.fftn(amp)
        power = F * np.conj(F)
        autocorr = np.real(np.fft.ifftn(power)) / N

        # Radial average using Manhattan distance
        for r in range(max_r + 1):
            mask = (dist_grid == r)
            if np.any(mask):
                C_r_total[r] += np.mean(autocorr[mask])

    # Average over eigenstates
    C_r = C_r_total / max(n_states, 1)

    # Subtract uniform background: C_connected(r) = C(r) - <|ψ|>²
    # (analogous to subtracting 1/q in Potts)
    if len(C_r) > 0 and C_r[0] > 0:
        # Background ~ (1/N)² * N = 1/N for extended state
        bg = np.mean(C_r[max_r//2:]) if max_r > 2 else 0.0
        C_r_connected = C_r - bg
        # Don't let it go negative at small r
        C_r_connected = np.maximum(C_r_connected, 0.0)
        return C_r_connected

    return C_r


# =====================================================================
# MODULE 4: Thermal-style FIM on 3D Lattice
# =====================================================================

def build_correlation_kernel_3d(C_r, L, v0):
    """
    Build probability distribution from C(r) centered at v0 = (x0, y0, z0).
    3D Manhattan distance, periodic BC.
    Returns flat array of shape (L³,).
    """
    x0, y0, z0 = v0
    max_r = len(C_r) - 1
    N = L ** 3

    ix = np.arange(L)
    dx = np.minimum(np.abs(ix - x0), L - np.abs(ix - x0))
    dy = np.minimum(np.abs(ix - y0), L - np.abs(ix - y0))
    dz = np.minimum(np.abs(ix - z0), L - np.abs(ix - z0))
    DX, DY, DZ = np.meshgrid(dx, dy, dz, indexing='ij')
    dist = DX + DY + DZ

    weights = np.zeros((L, L, L), dtype=np.float64)
    for r in range(max_r + 1):
        mask = (dist == r)
        weights[mask] = np.abs(C_r[r])

    beyond = dist > max_r
    if np.any(beyond):
        weights[beyond] = np.abs(C_r[max_r])

    p = weights.flatten()
    total = np.sum(p)
    if total < 1e-30:
        p = np.ones(N) / N
    else:
        p = p / total

    return p


def compute_FIM_3d(C_r, L, v0):
    """
    Compute 6×6 Fisher Information Matrix at vertex v0 on 3D cubic lattice.
    6 neighbors (±x, ±y, ±z). Score vectors from log-probability differences.
    """
    x0, y0, z0 = v0
    N = L ** 3

    p_v0 = build_correlation_kernel_3d(C_r, L, v0)
    log_p_v0 = np.log(p_v0 + 1e-30)

    neighbors = [
        ((x0 - 1) % L, y0, z0),
        ((x0 + 1) % L, y0, z0),
        (x0, (y0 - 1) % L, z0),
        (x0, (y0 + 1) % L, z0),
        (x0, y0, (z0 - 1) % L),
        (x0, y0, (z0 + 1) % L),
    ]
    k = len(neighbors)

    score_vectors = np.zeros((k, N))
    for j, w in enumerate(neighbors):
        p_w = build_correlation_kernel_3d(C_r, L, w)
        log_p_w = np.log(p_w + 1e-30)
        score_vectors[j, :] = log_p_w - log_p_v0

    weighted_scores = score_vectors * np.sqrt(p_v0)[np.newaxis, :]
    FIM = weighted_scores @ weighted_scores.T

    return FIM


# =====================================================================
# MODULE 5: FIM Diagnostics
# =====================================================================

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
    """eta = sv[rank] / sv[rank-1]"""
    if rank >= len(sv_norm) or rank < 1:
        return 0.0
    return sv_norm[rank] / sv_norm[rank - 1] if sv_norm[rank - 1] > 1e-15 else 0.0


def fisher_diagnostics_3d(C_r, L, n_samples=30):
    """
    Sample n_samples random vertices on 3D lattice.
    Compute FIM (6×6), SVD, rank, PR, eta, sv_profile.
    """
    ranks = []
    prs = []
    etas = []
    sv_profiles_list = []

    margin = max(1, L // 8)
    sampled = 0
    attempts = 0
    while sampled < n_samples and attempts < n_samples * 5:
        x = np.random.randint(margin, L - margin)
        y = np.random.randint(margin, L - margin)
        z = np.random.randint(margin, L - margin)
        attempts += 1

        FIM = compute_FIM_3d(C_r, L, (x, y, z))
        svs = np.linalg.svd(FIM, compute_uv=False)

        if svs[0] < 1e-30:
            continue

        sv_n = svs / svs[0]
        r = gap_based_rank(sv_n)
        pr = participation_ratio(svs)
        eta = disorder_index(sv_n, r)

        ranks.append(r)
        prs.append(pr)
        etas.append(eta)
        sv_profiles_list.append(sv_n)
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
# MODULE 6: Single Disorder Strength Runner
# =====================================================================

def run_disorder_strength(W, w_idx, n_W, n_realizations=10):
    """Run full Anderson analysis at one disorder strength W."""
    t0 = time.time()

    all_ipr = []
    all_sv_profiles = []
    all_ranks = []
    all_prs = []
    all_etas = []

    for real in range(n_realizations):
        seed = 42 + w_idx * 1000 + real

        # Build Hamiltonian
        H = build_anderson_hamiltonian(L, W, t=T_HOP, seed=seed)

        # Diagonalize near E=0
        evals, evecs = compute_eigenstates(H, k=K_EIG)
        if len(evals) == 0:
            continue

        # IPR
        ipr = compute_ipr(evecs)
        all_ipr.extend(ipr.tolist())

        # Wavefunction correlation C(r) — use fast FFT version
        C_r = measure_wavefunction_correlations_fast(evecs, L)

        # FIM diagnostics
        fisher = fisher_diagnostics_3d(C_r, L, n_samples=N_FIM_SAMPLES)
        if fisher is not None:
            sv_raw = fisher['sv_profile_mean']
            all_sv_profiles.append(sv_raw)
            all_ranks.append(fisher['rank_mean'])
            all_prs.append(fisher['pr_mean'])
            all_etas.append(fisher['eta_mean'])

    dt = time.time() - t0

    # Aggregate across realizations
    if not all_sv_profiles:
        return None

    sv_arr = np.array(all_sv_profiles)
    ipr_arr = np.array(all_ipr)

    # SV₂/SV₁ from mean SV profile
    sv_mean = sv_arr.mean(axis=0)
    sv_std = sv_arr.std(axis=0)
    sv2sv1 = sv_mean[1] if len(sv_mean) > 1 else 0.0
    sv2sv1_std = sv_std[1] if len(sv_std) > 1 else 0.0

    result = {
        'W': W,
        'W_over_Wc': W / W_C,
        'sv2sv1_mean': float(sv2sv1),
        'sv2sv1_std': float(sv2sv1_std),
        'rank_mean': float(np.mean(all_ranks)),
        'eta_mean': float(np.mean(all_etas)),
        'pr_mean': float(np.mean(all_prs)),
        'ipr_mean': float(np.mean(ipr_arr)),
        'ipr_std': float(np.std(ipr_arr)),
        'ipr_median': float(np.median(ipr_arr)),
        'n_realizations': len(all_sv_profiles),
        'n_eigenstates': len(all_ipr),
        'sv_profile': sv_mean.tolist(),
        'time_sec': dt,
    }

    print(f"  [{w_idx+1}/{n_W}] W={W:5.1f} (W/Wc={W/W_C:.3f}): "
          f"SV₂/SV₁={sv2sv1:.4f}±{sv2sv1_std:.4f} "
          f"rank={np.mean(all_ranks):.2f} η={np.mean(all_etas):.4f} "
          f"PR={np.mean(all_prs):.3f} IPR={np.mean(ipr_arr):.4f} "
          f"[{dt:.1f}s, {len(all_sv_profiles)} real]", flush=True)

    return result


# =====================================================================
# MODULE 7: Prediction Evaluation
# =====================================================================

def evaluate_predictions(results):
    """Evaluate pre-registered predictions."""
    W_vals = np.array([r['W'] for r in results])
    sv2sv1 = np.array([r['sv2sv1_mean'] for r in results])
    ranks = np.array([r['rank_mean'] for r in results])
    etas = np.array([r['eta_mean'] for r in results])
    iprs = np.array([r['ipr_mean'] for r in results])

    lines = []
    lines.append("=" * 70)
    lines.append("PHASE 4C: ANDERSON METAL-INSULATOR TRANSITION — VERDICTS")
    lines.append("=" * 70)
    lines.append(f"L = {L}, N = {N_SITES}, W_c = {W_C}")
    lines.append(f"Eigenstates: {K_EIG} near E=0, {N_REALIZATIONS} realizations/W")
    lines.append("")

    # P4C-1 [KILL]: SV₂/SV₁ peaks in W ∈ [14.5, 18.5]
    peak_idx = np.argmax(sv2sv1)
    peak_W = W_vals[peak_idx]
    peak_sv = sv2sv1[peak_idx]
    p1_pass = 14.5 <= peak_W <= 18.5
    lines.append(f"P4C-1 [KILL TEST]: SV₂/SV₁ peak at W={peak_W:.1f} (value={peak_sv:.4f})")
    lines.append(f"  Window: [14.5, 18.5], W_c={W_C}")
    lines.append(f"  Offset: {abs(peak_W - W_C):.1f} ({100*abs(peak_W - W_C)/W_C:.1f}%)")
    lines.append(f"  → {'PASS ✓' if p1_pass else 'KILL ✗'}")
    lines.append("")

    # P4C-2: SV₂/SV₁(W=2) > SV₂/SV₁(W=24)
    sv_low = sv2sv1[np.argmin(np.abs(W_vals - 2.0))]
    sv_high = sv2sv1[np.argmin(np.abs(W_vals - 24.0))]
    p2_pass = sv_low > sv_high
    lines.append(f"P4C-2: SV₂/SV₁(W≈2)={sv_low:.4f} vs SV₂/SV₁(W≈24)={sv_high:.4f}")
    lines.append(f"  → {'PASS ✓' if p2_pass else 'FAIL ✗'}")
    lines.append("")

    # P4C-3: η minimum aligns with SV₂/SV₁ maximum
    eta_min_idx = np.argmin(etas)
    eta_min_W = W_vals[eta_min_idx]
    p3_align = abs(eta_min_W - peak_W) / W_C < 0.20
    lines.append(f"P4C-3: η minimum at W={eta_min_W:.1f}, SV₂/SV₁ peak at W={peak_W:.1f}")
    lines.append(f"  Separation: {abs(eta_min_W - peak_W):.1f} ({100*abs(eta_min_W - peak_W)/W_C:.1f}% of W_c)")
    lines.append(f"  → {'PASS ✓' if p3_align else 'INCONCLUSIVE'}")
    lines.append("")

    # P4C-4: IPR transition midpoint ↔ SV₂/SV₁ peak (±20%)
    # IPR midpoint: halfway between min (metallic) and max (insulating)
    ipr_mid = (np.min(iprs) + np.max(iprs)) / 2.0
    ipr_mid_idx = np.argmin(np.abs(iprs - ipr_mid))
    ipr_mid_W = W_vals[ipr_mid_idx]
    p4_align = abs(ipr_mid_W - peak_W) / W_C < 0.20
    lines.append(f"P4C-4: IPR midpoint at W≈{ipr_mid_W:.1f}, SV₂/SV₁ peak at W={peak_W:.1f}")
    lines.append(f"  Separation: {abs(ipr_mid_W - peak_W):.1f} ({100*abs(ipr_mid_W - peak_W)/W_C:.1f}% of W_c)")
    lines.append(f"  → {'PASS ✓' if p4_align else 'INCONCLUSIVE'}")
    lines.append("")

    # P4C-5: Rank elevation at W_c vs endpoints
    rank_at_peak = ranks[peak_idx]
    rank_low = ranks[np.argmin(np.abs(W_vals - 2.0))]
    rank_high = ranks[np.argmin(np.abs(W_vals - 24.0))]
    p5_pass = rank_at_peak > rank_low and rank_at_peak > rank_high
    lines.append(f"P4C-5: Rank at peak W={peak_W:.1f}: {rank_at_peak:.2f} "
                 f"vs W≈2: {rank_low:.2f}, W≈24: {rank_high:.2f}")
    lines.append(f"  → {'PASS ✓' if p5_pass else 'INCONCLUSIVE'}")
    lines.append("")

    # Overall classification
    lines.append("=" * 70)
    if p1_pass:
        lines.append("OVERALL: PASS — FIM detects Anderson metal-insulator transition")
    else:
        lines.append("OVERALL: KILL — FIM does NOT detect Anderson transition")
    verdict_str = (f"P4C-1={'PASS' if p1_pass else 'KILL'} "
                   f"P4C-2={'PASS' if p2_pass else 'FAIL'} "
                   f"P4C-3={'PASS' if p3_align else '?'} "
                   f"P4C-4={'PASS' if p4_align else '?'} "
                   f"P4C-5={'PASS' if p5_pass else '?'}")
    lines.append(verdict_str)
    lines.append("=" * 70)

    # Full data table
    lines.append("\n--- Full Data ---")
    lines.append(f"{'W':>6} {'W/Wc':>6} {'SV2/SV1':>10} {'±':>8} {'rank':>6} "
                 f"{'η':>8} {'PR':>7} {'IPR':>9} {'IPR_std':>9} {'n_real':>6}")
    for r in results:
        lines.append(f"{r['W']:6.1f} {r['W_over_Wc']:6.3f} "
                     f"{r['sv2sv1_mean']:10.4f} {r['sv2sv1_std']:8.4f} "
                     f"{r['rank_mean']:6.2f} {r['eta_mean']:8.4f} "
                     f"{r['pr_mean']:7.3f} {r['ipr_mean']:9.6f} "
                     f"{r['ipr_std']:9.6f} {r['n_realizations']:6d}")

    return "\n".join(lines)


# =====================================================================
# MODULE 8: Plots
# =====================================================================

def make_plots(results):
    """Generate diagnostic plots."""
    if not HAS_MATPLOTLIB:
        print("  [matplotlib not available, skipping plots]", flush=True)
        return

    W_vals = [r['W'] for r in results]
    sv2sv1 = [r['sv2sv1_mean'] for r in results]
    sv2sv1_err = [r['sv2sv1_std'] for r in results]
    ranks = [r['rank_mean'] for r in results]
    etas = [r['eta_mean'] for r in results]
    prs = [r['pr_mean'] for r in results]
    iprs = [r['ipr_mean'] for r in results]
    ipr_errs = [r['ipr_std'] for r in results]

    fig, axes = plt.subplots(2, 3, figsize=(18, 10))
    fig.suptitle(f"Phase 4C: 3D Anderson Transition (L={L}, W_c={W_C})",
                fontsize=14, fontweight='bold')

    # A: SV₂/SV₁ vs W
    ax = axes[0, 0]
    ax.errorbar(W_vals, sv2sv1, yerr=sv2sv1_err, fmt='o-', color='#2196F3',
               capsize=3, label='SV₂/SV₁')
    ax.axvline(W_C, color='red', linestyle='--', alpha=0.7, label=f'W_c={W_C}')
    ax.axvspan(14.5, 18.5, alpha=0.1, color='red', label='Kill window')
    ax.set_xlabel('W (disorder)'); ax.set_ylabel('SV₂/SV₁')
    ax.set_title('FIM Spectral Ratio'); ax.legend(fontsize=7); ax.grid(True, alpha=0.3)

    # B: IPR vs W
    ax = axes[0, 1]
    ax.errorbar(W_vals, iprs, yerr=ipr_errs, fmt='s-', color='#4CAF50',
               capsize=3, label='IPR')
    ax.axvline(W_C, color='red', linestyle='--', alpha=0.7)
    ax.set_xlabel('W (disorder)'); ax.set_ylabel('IPR')
    ax.set_title('Inverse Participation Ratio'); ax.legend(fontsize=7)
    ax.set_yscale('log'); ax.grid(True, alpha=0.3)

    # C: SV₂/SV₁ and IPR together (dual axis)
    ax = axes[0, 2]
    ax.plot(W_vals, sv2sv1, 'o-', color='#2196F3', label='SV₂/SV₁')
    ax.set_xlabel('W (disorder)'); ax.set_ylabel('SV₂/SV₁', color='#2196F3')
    ax.axvline(W_C, color='red', linestyle='--', alpha=0.7)
    ax2 = ax.twinx()
    ax2.plot(W_vals, iprs, 's-', color='#4CAF50', label='IPR')
    ax2.set_ylabel('IPR', color='#4CAF50')
    ax.set_title('SV₂/SV₁ vs IPR'); ax.grid(True, alpha=0.3)

    # D: Rank
    ax = axes[1, 0]
    ax.plot(W_vals, ranks, 'D-', color='#FF9800')
    ax.axvline(W_C, color='red', linestyle='--', alpha=0.7)
    ax.set_xlabel('W'); ax.set_ylabel('Rank')
    ax.set_title('Gap-Based Rank'); ax.grid(True, alpha=0.3)

    # E: η (disorder index)
    ax = axes[1, 1]
    ax.plot(W_vals, etas, 'o-', color='#9C27B0')
    ax.axvline(W_C, color='red', linestyle='--', alpha=0.7)
    ax.set_xlabel('W'); ax.set_ylabel('η')
    ax.set_title('Disorder Index'); ax.grid(True, alpha=0.3)

    # F: PR
    ax = axes[1, 2]
    ax.plot(W_vals, prs, '^-', color='#795548')
    ax.axvline(W_C, color='red', linestyle='--', alpha=0.7)
    ax.set_xlabel('W'); ax.set_ylabel('PR')
    ax.set_title('Participation Ratio'); ax.grid(True, alpha=0.3)

    plt.tight_layout()
    path = os.path.join(RESULTS_DIR, "phase4c_anderson_diagnostics.png")
    plt.savefig(path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  Plot saved: {path}", flush=True)

    # SV profile at peak
    peak_idx = np.argmax([r['sv2sv1_mean'] for r in results])
    sv_prof = results[peak_idx]['sv_profile']
    fig2, ax2 = plt.subplots(figsize=(8, 5))
    ax2.bar(range(1, len(sv_prof)+1), sv_prof, color='#2196F3', alpha=0.8)
    ax2.axhline(0.01, color='red', linestyle='--', alpha=0.5, label='Threshold 0.01')
    ax2.set_xlabel('SV Index'); ax2.set_ylabel('Normalized SV')
    ax2.set_title(f'Anderson SV Profile at W={results[peak_idx]["W"]:.1f} '
                  f'(SV₂/SV₁={sv_prof[1]:.4f})')
    ax2.legend()
    path2 = os.path.join(RESULTS_DIR, "phase4c_sv_profile_at_peak.png")
    plt.savefig(path2, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  Plot saved: {path2}", flush=True)


# =====================================================================
# MAIN
# =====================================================================

def main():
    os.makedirs(RESULTS_DIR, exist_ok=True)

    print("=" * 70, flush=True)
    print("PHASE 4C: 3D ANDERSON METAL-INSULATOR TRANSITION — FIM TEST", flush=True)
    print("=" * 70, flush=True)
    print(f"L = {L}, N = {N_SITES} sites, W_c = {W_C}", flush=True)
    print(f"W sweep: {len(W_VALUES)} values from {W_VALUES[0]:.1f} to {W_VALUES[-1]:.1f}", flush=True)
    print(f"Eigenstates: {K_EIG} near E=0, {N_REALIZATIONS} realizations/W", flush=True)
    print(f"FIM: {N_FIM_SAMPLES} samples/realization", flush=True)
    print(f"Start: {datetime.datetime.now().isoformat()}", flush=True)
    print(flush=True)

    print("Two-scale commensurability:", flush=True)
    print(f"  Scale 1: Localization length ξ_loc (diverges at W_c)", flush=True)
    print(f"  Scale 2: System size L = {L}", flush=True)
    print(f"  At W_c: ξ_loc ~ L → commensurability → FIM peak expected", flush=True)
    print(flush=True)

    t_start = time.time()
    results = []

    for i, W in enumerate(W_VALUES):
        res = run_disorder_strength(W, i, len(W_VALUES), N_REALIZATIONS)
        if res is not None:
            results.append(res)

        # Interim save
        interim_path = os.path.join(RESULTS_DIR, "phase4c_interim.csv")
        with open(interim_path, 'w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=[
                'W', 'W_over_Wc', 'sv2sv1_mean', 'sv2sv1_std',
                'rank_mean', 'eta_mean', 'pr_mean',
                'ipr_mean', 'ipr_std', 'n_realizations', 'time_sec'
            ])
            writer.writeheader()
            for r in results:
                row = {k: f"{v:.6f}" if isinstance(v, float) else v
                       for k, v in r.items()
                       if k not in ('sv_profile', 'ipr_median', 'n_eigenstates',
                                     'W_over_Wc')}
                row['W_over_Wc'] = f"{r['W_over_Wc']:.6f}"
                writer.writerow(row)

    t_total = time.time() - t_start
    print(f"\nTotal runtime: {t_total:.1f}s ({t_total/60:.1f} min)", flush=True)

    # Write final CSV
    csv_path = os.path.join(RESULTS_DIR, "phase4c_anderson_results.csv")
    with open(csv_path, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=[
            'W', 'W_over_Wc', 'sv2sv1_mean', 'sv2sv1_std',
            'rank_mean', 'eta_mean', 'pr_mean',
            'ipr_mean', 'ipr_std', 'ipr_median',
            'n_realizations', 'n_eigenstates', 'time_sec'
        ])
        writer.writeheader()
        for r in results:
            row = {k: f"{v:.6f}" if isinstance(v, float) else v
                   for k, v in r.items() if k != 'sv_profile'}
            writer.writerow(row)
    print(f"Results: {csv_path}", flush=True)

    # Evaluate predictions
    summary = evaluate_predictions(results)
    summary_path = os.path.join(RESULTS_DIR, "phase4c_anderson_summary.txt")
    with open(summary_path, 'w') as f:
        f.write(summary)
    print(f"Summary: {summary_path}", flush=True)

    # Plots
    make_plots(results)

    # JSON summary
    peak_idx = np.argmax([r['sv2sv1_mean'] for r in results])
    json_path = os.path.join(RESULTS_DIR, "phase4c_anderson_summary.json")
    with open(json_path, 'w') as f:
        json.dump({
            'L': L,
            'N_sites': N_SITES,
            'W_c': W_C,
            'K_eig': K_EIG,
            'n_realizations': N_REALIZATIONS,
            'peak_W': float(results[peak_idx]['W']),
            'peak_sv2sv1': float(results[peak_idx]['sv2sv1_mean']),
            'peak_ipr': float(results[peak_idx]['ipr_mean']),
            'total_runtime_sec': t_total,
            'timestamp': datetime.datetime.now().isoformat(),
        }, f, indent=2)

    print(flush=True)
    print(summary, flush=True)


if __name__ == "__main__":
    main()
