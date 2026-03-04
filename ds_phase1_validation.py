#!/usr/bin/env python3
"""
DS Framework Phase 1 Validation
================================
Tests that effective dimensionality D_eff = rank(Fisher Information Matrix)
converges with geometric and spectral dimension estimates on benchmark
torus lattices where the true answer is known.

Test systems:
  - 2D flat torus (200x200): true dimension = 2.0
  - 3D flat torus (25x25x25): true dimension = 3.0

Three measurement routes:
  1. Growth dimension (BFS ball volumes)
  2. Spectral dimension (Laplacian eigenvalue Weyl law)
  3. Fisher information rank (DS framework core claim)

Phase 1 PASSES if all three routes converge within specified tolerances.
"""

import os
import time
import datetime
from collections import deque

import numpy as np
import networkx as nx
from scipy.sparse.linalg import eigsh
from scipy.sparse import diags
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# ============================================================================
# Configuration
# ============================================================================

OUTPUT_DIR = "phase1_results"
np.random.seed(42)

SYSTEMS = {
    '2d': {
        'name': '2D Torus (200x200)',
        'prefix': 'torus2d',
        'n': 200,
        'dim': 2,
        'true_d': 2.0,
        'expected_nodes': 40000,
        'expected_degree': 4,
        'gate_lo': 1.85,
        'gate_hi': 2.15,
        'agreement_gate': 0.10,
        'growth_samples': 30,
        'fisher_samples': 20,
        'n_eigenvalues': 500,
    },
    '3d': {
        'name': '3D Torus (25x25x25)',
        'prefix': 'torus3d',
        'n': 25,
        'dim': 3,
        'true_d': 3.0,
        'expected_nodes': 15625,
        'expected_degree': 6,
        'gate_lo': 2.70,
        'gate_hi': 3.30,
        'agreement_gate': 0.15,
        'growth_samples': 15,
        'fisher_samples': 20,
        'n_eigenvalues': 500,
    },
}

# ============================================================================
# Torus Construction
# ============================================================================

def make_torus_2d(n):
    """Build 2D flat torus manually. Every vertex has degree 4."""
    G = nx.Graph()
    G.add_nodes_from(range(n * n))
    for i in range(n):
        for j in range(n):
            v = i * n + j
            G.add_edge(v, i * n + (j + 1) % n)
            G.add_edge(v, ((i + 1) % n) * n + j)
    return G


def make_torus_3d(n):
    """Build 3D flat torus manually. Every vertex has degree 6."""
    G = nx.Graph()
    n2 = n * n
    G.add_nodes_from(range(n * n * n))
    for i in range(n):
        for j in range(n):
            for k in range(n):
                v = i * n2 + j * n + k
                G.add_edge(v, i * n2 + j * n + (k + 1) % n)
                G.add_edge(v, i * n2 + ((j + 1) % n) * n + k)
                G.add_edge(v, ((i + 1) % n) * n2 + j * n + k)
    return G


def validate_torus(G, cfg):
    """Validate torus topology: node count and degree uniformity."""
    n_nodes = G.number_of_nodes()
    assert n_nodes == cfg['expected_nodes'], \
        f"Expected {cfg['expected_nodes']} nodes, got {n_nodes}"
    degrees = set(dict(G.degree()).values())
    assert degrees == {cfg['expected_degree']}, \
        f"Expected all degree {cfg['expected_degree']}, got {degrees}"
    print(f"  Validated: {n_nodes} nodes, all degree {cfg['expected_degree']}")


# ============================================================================
# Route 1: Growth Dimension
# ============================================================================

def bfs_distances(G, source):
    """BFS returning dict of distances. O(V+E)."""
    dist = {source: 0}
    queue = deque([source])
    while queue:
        v = queue.popleft()
        for w in G.neighbors(v):
            if w not in dist:
                dist[w] = dist[v] + 1
                queue.append(w)
    return dist


def compute_ball_volumes(G, n_samples=30):
    """Compute average ball volumes V(r) for r = 0..maxR."""
    nodes = list(G.nodes())
    samples = np.random.choice(nodes, size=min(n_samples, len(nodes)), replace=False)

    d0 = bfs_distances(G, samples[0])
    max_dist = max(d0.values())
    maxR = max_dist // 2

    all_volumes = []
    for i, v in enumerate(samples):
        dists = bfs_distances(G, v)
        dist_vals = np.array(list(dists.values()))
        volumes = np.array([np.sum(dist_vals <= r) for r in range(maxR + 1)])
        all_volumes.append(volumes)
        if (i + 1) % 10 == 0:
            print(f"    BFS {i+1}/{len(samples)}")

    mean_volumes = np.mean(all_volumes, axis=0)
    return mean_volumes, maxR


def estimate_growth_dimension(mean_volumes, maxR):
    """Estimate growth dimension via log-log fit and point-by-point delta."""
    rMin = max(3, maxR // 4)
    rMax = (maxR * 2) // 3

    # Method 1: Log-log linear regression
    radii = np.arange(rMin, rMax + 1, dtype=float)
    log_r = np.log(radii)
    log_v = np.log(mean_volumes[rMin:rMax + 1].astype(float))
    coeffs = np.polyfit(log_r, log_v, 1)
    slope = coeffs[0]

    predicted = np.polyval(coeffs, log_r)
    ss_res = np.sum((log_v - predicted) ** 2)
    ss_tot = np.sum((log_v - np.mean(log_v)) ** 2)
    r_squared = 1.0 - ss_res / ss_tot

    # Method 2: Point-by-point delta
    deltas = []
    for r in range(1, maxR):
        if mean_volumes[r] > 0 and mean_volumes[r + 1] > 0:
            d = (np.log(float(mean_volumes[r + 1])) - np.log(float(mean_volumes[r]))) / \
                (np.log(float(r + 1)) - np.log(float(r)))
            deltas.append((r, d))

    delta_in_range = [d for r, d in deltas if rMin <= r <= rMax]
    mean_delta = np.mean(delta_in_range) if delta_in_range else float('nan')

    return {
        'fit_dimension': slope,
        'delta_dimension': mean_delta,
        'r_squared': r_squared,
        'rMin': rMin,
        'rMax': rMax,
        'maxR': maxR,
        'agreement': abs(slope - mean_delta),
        'deltas': deltas,
        'mean_volumes': mean_volumes,
    }


# ============================================================================
# Route 2: Spectral Dimension
# ============================================================================

def compute_spectral_dimension(G, n_eigenvalues=500):
    """Estimate dimension from Weyl law: N(lambda) ~ lambda^{d/2}."""
    A = nx.adjacency_matrix(G).astype(float)
    n = G.number_of_nodes()
    degree = list(dict(G.degree()).values())[0]

    D = diags([degree] * n, 0, format='csr')
    L = D - A

    k = min(n_eigenvalues + 1, n - 2)

    # Use shift-invert mode for numerical stability
    try:
        eigenvalues = eigsh(L, k=k, sigma=0.01, return_eigenvectors=False)
    except Exception:
        print("    shift-invert failed, trying which='SM'...")
        eigenvalues = eigsh(L, k=k, which='SM', return_eigenvectors=False)

    eigenvalues = np.sort(np.real(eigenvalues))
    eigenvalues = eigenvalues[eigenvalues > 1e-10]

    n_eigs = len(eigenvalues)
    indices = np.arange(1, n_eigs + 1, dtype=float)

    fit_start = max(5, n_eigs // 10)
    fit_end = (n_eigs * 4) // 5

    log_lam = np.log(eigenvalues[fit_start:fit_end])
    log_N = np.log(indices[fit_start:fit_end])

    coeffs = np.polyfit(log_lam, log_N, 1)
    slope = coeffs[0]
    spectral_dim = 2.0 * slope

    predicted = np.polyval(coeffs, log_lam)
    ss_res = np.sum((log_N - predicted) ** 2)
    ss_tot = np.sum((log_N - np.mean(log_N)) ** 2)
    r_squared = 1.0 - ss_res / ss_tot

    return {
        'spectral_dimension': spectral_dim,
        'weyl_slope': slope,
        'r_squared': r_squared,
        'n_eigenvalues_used': fit_end - fit_start,
        'fit_range': (fit_start, fit_end),
        'eigenvalues': eigenvalues,
    }


# ============================================================================
# Route 3: Fisher Information Rank
# ============================================================================

def compute_fisher_dimension(G, sigma=3.0, n_samples=20):
    """
    Estimate effective dimension via Fisher Information Rank.
    D_eff = rank(FIM) in neighbor directions.
    """
    nodes = list(G.nodes())
    n = len(nodes)

    samples = [nodes[i] for i in np.random.choice(n, size=n_samples, replace=False)]

    all_ranks = []
    all_threshold_ranks = []
    all_singular_values = []
    all_participation_ratios = []

    for s_idx, v0 in enumerate(samples):
        if (s_idx + 1) % 5 == 0:
            print(f"    Fisher sample {s_idx+1}/{n_samples}")

        # BFS distances from v0
        dists_v0 = bfs_distances(G, v0)
        dist_arr_v0 = np.array([dists_v0[u] for u in nodes], dtype=float)
        log_unnorm_v0 = -dist_arr_v0 / sigma
        log_Z_v0 = np.logaddexp.reduce(log_unnorm_v0)
        log_p_v0 = log_unnorm_v0 - log_Z_v0
        p_v0 = np.exp(log_p_v0)

        neighbors = list(G.neighbors(v0))
        k = len(neighbors)

        # Score vectors for each neighbor direction
        score_vectors = np.zeros((k, n))
        for j, w in enumerate(neighbors):
            dists_w = bfs_distances(G, w)
            dist_arr_w = np.array([dists_w[u] for u in nodes], dtype=float)
            log_unnorm_w = -dist_arr_w / sigma
            log_Z_w = np.logaddexp.reduce(log_unnorm_w)
            log_p_w = log_unnorm_w - log_Z_w
            score_vectors[j, :] = log_p_w - log_p_v0

        # Fisher Information Matrix in neighbor basis
        weighted_scores = score_vectors * np.sqrt(p_v0)[np.newaxis, :]
        F = weighted_scores @ weighted_scores.T

        # Effective rank via singular values
        sv = np.linalg.svd(F, compute_uv=False)
        sv_norm = sv / sv[0] if sv[0] > 0 else sv

        # Gap-based rank detection: find largest relative drop in SVs
        # ratio[i] = sv_norm[i+1] / sv_norm[i]; rank = argmin(ratio) + 1
        if len(sv_norm) > 1 and sv_norm[0] > 0:
            ratios = sv_norm[1:] / np.maximum(sv_norm[:-1], 1e-15)
            gap_rank = int(np.argmin(ratios) + 1)
        else:
            gap_rank = 1

        # Also keep fixed-threshold rank for diagnostics
        threshold = 0.01
        threshold_rank = int(np.sum(sv_norm > threshold))

        eff_rank = gap_rank

        # Participation ratio
        if np.sum(sv ** 2) > 0:
            participation_ratio = (np.sum(sv)) ** 2 / np.sum(sv ** 2)
        else:
            participation_ratio = 0.0

        all_ranks.append(eff_rank)
        all_threshold_ranks.append(threshold_rank)
        all_singular_values.append(sv_norm)
        all_participation_ratios.append(participation_ratio)

    mean_rank = np.mean(all_ranks)
    median_rank = np.median(all_ranks)
    mean_participation = np.mean(all_participation_ratios)

    # Average singular value profile
    max_len = max(len(sv) for sv in all_singular_values)
    padded = np.zeros((len(all_singular_values), max_len))
    for i, sv in enumerate(all_singular_values):
        padded[i, :len(sv)] = sv
    mean_sv = np.mean(padded, axis=0)

    return {
        'fisher_dimension_mean': mean_rank,
        'fisher_dimension_median': median_rank,
        'fisher_participation_ratio': mean_participation,
        'mean_singular_values': mean_sv,
        'individual_ranks': all_ranks,
        'sigma': sigma,
        'threshold': threshold,
    }


def fisher_sigma_sweep(G, sigmas=[1.5, 2.0, 3.0, 5.0, 8.0], n_samples=10):
    """Run Fisher dimension at multiple sigma values for sensitivity check."""
    results = {}
    for sigma in sigmas:
        print(f"    sigma = {sigma}")
        r = compute_fisher_dimension(G, sigma=sigma, n_samples=n_samples)
        results[sigma] = {
            'mean_rank': r['fisher_dimension_mean'],
            'participation_ratio': r['fisher_participation_ratio'],
            'sv_profile': r['mean_singular_values'],
        }
    return results


# ============================================================================
# Plotting
# ============================================================================

def plot_growth_delta(deltas, true_d, rMin, rMax, prefix, output_dir):
    """Plot delta(r) vs r with reference line at true dimension."""
    fig, ax = plt.subplots(figsize=(10, 6))
    rs = [r for r, d in deltas]
    ds = [d for r, d in deltas]
    ax.plot(rs, ds, 'b.-', alpha=0.7, markersize=4)
    ax.axhline(y=true_d, color='r', linestyle='--', linewidth=2, label=f'True d = {true_d}')
    ax.axvspan(rMin, rMax, alpha=0.15, color='green', label=f'Fitting range [{rMin}, {rMax}]')
    ax.set_xlabel('Radius r', fontsize=12)
    ax.set_ylabel('Local dimension Δ(r)', fontsize=12)
    ax.set_title(f'Growth Dimension Δ(r) — {prefix}', fontsize=14)
    ax.legend(fontsize=11)
    ax.grid(True, alpha=0.3)
    ax.set_ylim(max(0, true_d - 2), true_d + 2)
    path = os.path.join(output_dir, f'{prefix}_growth_delta.png')
    fig.savefig(path, dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f"    Saved {path}")


def plot_growth_loglog(mean_volumes, rMin, rMax, slope, true_d, prefix, output_dir):
    """Plot log V(r) vs log r with fitted line."""
    fig, ax = plt.subplots(figsize=(10, 6))
    valid = np.arange(1, len(mean_volumes))
    log_r = np.log(valid.astype(float))
    log_v = np.log(mean_volumes[1:].astype(float))
    ax.plot(log_r, log_v, 'b.', alpha=0.5, markersize=3)

    # Fitted line in range
    fit_r = np.arange(rMin, rMax + 1, dtype=float)
    fit_log_r = np.log(fit_r)
    fit_log_v = np.log(mean_volumes[rMin:rMax + 1].astype(float))
    coeffs = np.polyfit(fit_log_r, fit_log_v, 1)
    ax.plot(fit_log_r, np.polyval(coeffs, fit_log_r), 'r-', linewidth=2,
            label=f'Fit: slope = {slope:.4f} (true = {true_d})')

    ax.set_xlabel('log(r)', fontsize=12)
    ax.set_ylabel('log(V(r))', fontsize=12)
    ax.set_title(f'Growth Dimension log-log — {prefix}', fontsize=14)
    ax.legend(fontsize=11)
    ax.grid(True, alpha=0.3)
    path = os.path.join(output_dir, f'{prefix}_growth_loglog.png')
    fig.savefig(path, dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f"    Saved {path}")


def plot_spectral_weyl(eigenvalues, fit_start, fit_end, weyl_slope, spectral_dim,
                       true_d, prefix, output_dir):
    """Plot log N(lambda) vs log(lambda) with fitted line."""
    fig, ax = plt.subplots(figsize=(10, 6))
    n_eigs = len(eigenvalues)
    indices = np.arange(1, n_eigs + 1, dtype=float)
    log_lam = np.log(eigenvalues)
    log_N = np.log(indices)
    ax.plot(log_lam, log_N, 'b.', alpha=0.4, markersize=3)

    # Fitted line
    fit_log_lam = log_lam[fit_start:fit_end]
    fit_log_N = log_N[fit_start:fit_end]
    coeffs = np.polyfit(fit_log_lam, fit_log_N, 1)
    ax.plot(fit_log_lam, np.polyval(coeffs, fit_log_lam), 'r-', linewidth=2,
            label=f'Weyl slope = {weyl_slope:.4f} → d = {spectral_dim:.4f} (true = {true_d})')

    ax.set_xlabel('log(λ)', fontsize=12)
    ax.set_ylabel('log(N(λ))', fontsize=12)
    ax.set_title(f'Spectral Dimension (Weyl Law) — {prefix}', fontsize=14)
    ax.legend(fontsize=11)
    ax.grid(True, alpha=0.3)
    path = os.path.join(output_dir, f'{prefix}_spectral_weyl.png')
    fig.savefig(path, dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f"    Saved {path}")


def plot_fisher_sv(mean_sv, true_d, prefix, output_dir):
    """Bar chart of normalized singular values of FIM."""
    fig, ax = plt.subplots(figsize=(8, 6))
    x = np.arange(1, len(mean_sv) + 1)
    colors = ['#2196F3' if i < true_d else '#BDBDBD' for i in range(len(mean_sv))]
    ax.bar(x, mean_sv, color=colors, edgecolor='black', linewidth=0.5)
    ax.axhline(y=0.01, color='r', linestyle='--', alpha=0.7, label='Threshold (0.01)')
    ax.set_xlabel('Singular Value Index', fontsize=12)
    ax.set_ylabel('Normalized Singular Value', fontsize=12)
    ax.set_title(f'Fisher Information Singular Values — {prefix}', fontsize=14)
    ax.set_xticks(x)
    ax.legend(fontsize=11)
    ax.grid(True, alpha=0.3, axis='y')
    path = os.path.join(output_dir, f'{prefix}_fisher_sv.png')
    fig.savefig(path, dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f"    Saved {path}")


def plot_convergence(estimates, true_d, gate_lo, gate_hi, prefix, output_dir):
    """Bar chart comparing all three route estimates."""
    fig, ax = plt.subplots(figsize=(8, 6))
    labels = ['Route 1\nGrowth', 'Route 2\nSpectral', 'Route 3\nFisher']
    vals = [estimates['growth'], estimates['spectral'], estimates['fisher']]
    colors = ['#4CAF50' if gate_lo <= v <= gate_hi else '#F44336' for v in vals]

    bars = ax.bar(labels, vals, color=colors, edgecolor='black', linewidth=0.8, width=0.5)

    ax.axhline(y=true_d, color='black', linestyle='-', linewidth=2, label=f'True d = {true_d}')
    ax.axhspan(gate_lo, gate_hi, alpha=0.12, color='green', label=f'Gate [{gate_lo}, {gate_hi}]')

    for bar, val in zip(bars, vals):
        ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.02,
                f'{val:.3f}', ha='center', va='bottom', fontsize=12, fontweight='bold')

    ax.set_ylabel('Dimension Estimate', fontsize=12)
    ax.set_title(f'Convergence Summary — {prefix}', fontsize=14)
    ax.legend(fontsize=11, loc='upper right')
    ax.grid(True, alpha=0.3, axis='y')
    ax.set_ylim(min(gate_lo - 0.3, min(vals) - 0.3), max(gate_hi + 0.3, max(vals) + 0.3))
    path = os.path.join(output_dir, f'{prefix}_convergence.png')
    fig.savefig(path, dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f"    Saved {path}")


# ============================================================================
# Main Execution
# ============================================================================

def run_system(key, cfg):
    """Run all three routes on one test system."""
    print(f"\n{'='*70}")
    print(f"  {cfg['name']}")
    print(f"{'='*70}")

    results = {'config': cfg}

    # --- Step 0: Build torus ---
    print(f"\n[Step 0] Building {cfg['name']}...")
    t0 = time.time()
    if cfg['dim'] == 2:
        G = make_torus_2d(cfg['n'])
    else:
        G = make_torus_3d(cfg['n'])
    validate_torus(G, cfg)
    print(f"  Construction time: {time.time()-t0:.1f}s")

    # --- Step 1: Route 1 — Growth Dimension ---
    print(f"\n[Route 1] Growth Dimension...")
    t0 = time.time()
    mean_volumes, maxR = compute_ball_volumes(G, n_samples=cfg['growth_samples'])
    growth = estimate_growth_dimension(mean_volumes, maxR)
    growth_time = time.time() - t0

    growth_est = growth['fit_dimension']
    gate_pass = cfg['gate_lo'] <= growth_est <= cfg['gate_hi']
    print(f"  Fit dimension:   {growth_est:.4f}")
    print(f"  Delta dimension: {growth['delta_dimension']:.4f}")
    print(f"  R-squared:       {growth['r_squared']:.6f}")
    print(f"  Fit range:       [{growth['rMin']}, {growth['rMax']}] (maxR={maxR})")
    print(f"  Gate [{cfg['gate_lo']}, {cfg['gate_hi']}]: {'PASS' if gate_pass else 'FAIL'}")
    print(f"  Time: {growth_time:.1f}s")
    results['growth'] = growth
    results['growth_pass'] = gate_pass

    plot_growth_delta(growth['deltas'], cfg['true_d'], growth['rMin'], growth['rMax'],
                      cfg['prefix'], OUTPUT_DIR)
    plot_growth_loglog(mean_volumes, growth['rMin'], growth['rMax'],
                       growth_est, cfg['true_d'], cfg['prefix'], OUTPUT_DIR)

    # --- Step 2: Route 2 — Spectral Dimension ---
    print(f"\n[Route 2] Spectral Dimension...")
    t0 = time.time()
    spectral = compute_spectral_dimension(G, n_eigenvalues=cfg['n_eigenvalues'])
    spectral_time = time.time() - t0

    spectral_est = spectral['spectral_dimension']
    gate_pass = cfg['gate_lo'] <= spectral_est <= cfg['gate_hi']
    print(f"  Spectral dimension: {spectral_est:.4f}")
    print(f"  Weyl slope:         {spectral['weyl_slope']:.4f}")
    print(f"  R-squared:          {spectral['r_squared']:.6f}")
    print(f"  Eigenvalues used:   {spectral['n_eigenvalues_used']}")
    print(f"  Gate [{cfg['gate_lo']}, {cfg['gate_hi']}]: {'PASS' if gate_pass else 'FAIL'}")
    print(f"  Time: {spectral_time:.1f}s")
    results['spectral'] = spectral
    results['spectral_pass'] = gate_pass

    plot_spectral_weyl(spectral['eigenvalues'], spectral['fit_range'][0],
                       spectral['fit_range'][1], spectral['weyl_slope'],
                       spectral_est, cfg['true_d'], cfg['prefix'], OUTPUT_DIR)

    # --- Step 3: Route 3 — Fisher Information Rank ---
    print(f"\n[Route 3] Fisher Information Rank...")
    t0 = time.time()
    fisher = compute_fisher_dimension(G, sigma=3.0, n_samples=cfg['fisher_samples'])
    fisher_time = time.time() - t0

    fisher_est = fisher['fisher_dimension_mean']
    gate_pass = cfg['gate_lo'] <= fisher_est <= cfg['gate_hi']
    print(f"  Fisher dimension (mean):   {fisher_est:.4f}")
    print(f"  Fisher dimension (median): {fisher['fisher_dimension_median']:.4f}")
    print(f"  Participation ratio:       {fisher['fisher_participation_ratio']:.4f}")
    print(f"  Singular values:           {np.array2string(fisher['mean_singular_values'], precision=4)}")
    print(f"  Individual ranks:          {fisher['individual_ranks']}")
    print(f"  sigma = {fisher['sigma']}, threshold = {fisher['threshold']}")
    print(f"  Gate [{cfg['gate_lo']}, {cfg['gate_hi']}]: {'PASS' if gate_pass else 'FAIL'}")
    print(f"  Time: {fisher_time:.1f}s")
    results['fisher'] = fisher
    results['fisher_pass'] = gate_pass

    plot_fisher_sv(fisher['mean_singular_values'], cfg['true_d'],
                   cfg['prefix'], OUTPUT_DIR)

    # --- Sigma sensitivity sweep ---
    print(f"\n  [Fisher sigma sweep]")
    sweep = fisher_sigma_sweep(G, sigmas=[1.5, 2.0, 3.0, 5.0, 8.0], n_samples=10)
    for sigma, sr in sweep.items():
        print(f"    sigma={sigma}: rank={sr['mean_rank']:.2f}, "
              f"participation={sr['participation_ratio']:.3f}")
    results['fisher_sweep'] = sweep

    # --- Convergence ---
    estimates = {
        'growth': growth['fit_dimension'],
        'spectral': spectral['spectral_dimension'],
        'fisher': fisher['fisher_dimension_mean'],
    }
    vals = list(estimates.values())
    max_diff = max(vals) - min(vals)
    agreement_pass = max_diff < cfg['agreement_gate']

    print(f"\n[Convergence]")
    print(f"  Route 1 (Growth):   {estimates['growth']:.4f}")
    print(f"  Route 2 (Spectral): {estimates['spectral']:.4f}")
    print(f"  Route 3 (Fisher):   {estimates['fisher']:.4f}")
    print(f"  Max pairwise diff:  {max_diff:.4f}")
    print(f"  Agreement gate (<{cfg['agreement_gate']}): {'PASS' if agreement_pass else 'FAIL'}")

    results['estimates'] = estimates
    results['max_diff'] = max_diff
    results['agreement_pass'] = agreement_pass

    plot_convergence(estimates, cfg['true_d'], cfg['gate_lo'], cfg['gate_hi'],
                     cfg['prefix'], OUTPUT_DIR)

    return results


def generate_report(all_results, total_time):
    """Generate PHASE1_RESULTS.md report."""
    lines = []
    lines.append("# DS Framework Phase 1 Validation Results")
    lines.append(f"## Date: {datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    lines.append(f"## Runtime: {total_time:.1f} seconds")
    lines.append("")

    overall_pass = True
    any_pass = False

    for key in ['2d', '3d']:
        r = all_results[key]
        cfg = r['config']
        est = r['estimates']
        lines.append(f"### {cfg['name']}")
        lines.append("")
        lines.append("| Route | Method | Estimate | Gate [{:.2f}, {:.2f}] | Status |".format(
            cfg['gate_lo'], cfg['gate_hi']))
        lines.append("|-------|--------|----------|-------------------|--------|")

        growth_pass = r['growth_pass']
        spectral_pass = r['spectral_pass']
        fisher_pass = r['fisher_pass']

        lines.append(f"| 1 | Growth Dimension | {est['growth']:.4f} "
                     f"(R²={r['growth']['r_squared']:.4f}) | "
                     f"[{cfg['gate_lo']}, {cfg['gate_hi']}] | "
                     f"{'PASS' if growth_pass else 'FAIL'} |")
        lines.append(f"| 2 | Spectral Dimension | {est['spectral']:.4f} "
                     f"(R²={r['spectral']['r_squared']:.4f}) | "
                     f"[{cfg['gate_lo']}, {cfg['gate_hi']}] | "
                     f"{'PASS' if spectral_pass else 'FAIL'} |")
        lines.append(f"| 3 | Fisher Info Rank | {est['fisher']:.4f} "
                     f"(PR={r['fisher']['fisher_participation_ratio']:.3f}) | "
                     f"[{cfg['gate_lo']}, {cfg['gate_hi']}] | "
                     f"{'PASS' if fisher_pass else 'FAIL'} |")
        lines.append("")
        lines.append(f"Cross-route agreement: {r['max_diff']:.4f} "
                     f"(gate: < {cfg['agreement_gate']}) -> "
                     f"{'PASS' if r['agreement_pass'] else 'FAIL'}")
        lines.append("")

        all_pass = growth_pass and spectral_pass and fisher_pass and r['agreement_pass']
        if all_pass:
            any_pass = True
        else:
            overall_pass = False

    # Overall verdict
    if overall_pass and any_pass:
        verdict = "PASS"
    elif any_pass:
        verdict = "PARTIAL"
    else:
        verdict = "FAIL"

    lines.append(f"### Phase 1 Overall: **{verdict}**")
    lines.append("")
    lines.append("- **PASS**: All individual gates pass AND both convergence gates pass")
    lines.append("- **PARTIAL**: Some routes pass, others marginal")
    lines.append("- **FAIL**: Core route(s) fail")
    lines.append("")

    # Diagnostics
    lines.append("### Diagnostics")
    lines.append("")
    for key in ['2d', '3d']:
        r = all_results[key]
        cfg = r['config']
        lines.append(f"**{cfg['name']} — Fisher sigma sweep:**")
        lines.append("")
        lines.append("| sigma | Mean Rank | Participation Ratio |")
        lines.append("|-------|-----------|---------------------|")
        for sigma, sr in r['fisher_sweep'].items():
            lines.append(f"| {sigma} | {sr['mean_rank']:.2f} | {sr['participation_ratio']:.3f} |")
        lines.append("")
        lines.append(f"**{cfg['name']} — Fisher SV profile:** "
                     f"{np.array2string(r['fisher']['mean_singular_values'], precision=4)}")
        lines.append("")

    report = "\n".join(lines)
    path = os.path.join(OUTPUT_DIR, "PHASE1_RESULTS.md")
    with open(path, 'w') as f:
        f.write(report)
    print(f"\nReport saved to {path}")
    return report


def main():
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    t_start = time.time()

    all_results = {}
    for key in ['2d', '3d']:
        all_results[key] = run_system(key, SYSTEMS[key])

    total_time = time.time() - t_start
    report = generate_report(all_results, total_time)

    print("\n" + "=" * 70)
    print("  PHASE 1 COMPLETE")
    print("=" * 70)
    print(report)


if __name__ == '__main__':
    main()
