#!/usr/bin/env python3
"""
DS Phase 1 Extension: Sierpinski Gasket + 4D Torus
====================================================
Extends Phase 1 validation with:
  - Test 2A: Sierpinski gasket (L=6, L=7) — fractal with non-integer dimension
  - Test 3A: 4D flat torus (n=12) — integer dimension d=4

The Sierpinski test is the CRITICAL test: Fisher rank can only return integers,
but the gasket has Hausdorff dimension log(3)/log(2) ≈ 1.585. The question is
what the FIM singular value profile looks like on a fractal, and whether the
participation ratio yields a meaningful non-integer dimension estimate.

Known dimensions of Sierpinski gasket:
  - Hausdorff: d_H = log(3)/log(2) ≈ 1.5849
  - Spectral:  d_S = 2*log(3)/log(5) ≈ 1.3652
  - Walk:      d_W = log(5)/log(2) ≈ 2.3219
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

# Reference dimensions for Sierpinski gasket
D_HAUSDORFF = np.log(3) / np.log(2)       # ≈ 1.5849
D_SPECTRAL  = 2 * np.log(3) / np.log(5)   # ≈ 1.3652
D_WALK      = np.log(5) / np.log(2)        # ≈ 2.3219


# ============================================================================
# Shared Utility Functions (from ds_phase1_validation.py, self-contained)
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


def compute_ball_volumes(G, n_samples=30, sample_nodes=None):
    """Compute average ball volumes V(r) for r = 0..maxR."""
    nodes = list(G.nodes())
    if sample_nodes is not None:
        samples = sample_nodes[:n_samples]
    else:
        samples = [nodes[i] for i in np.random.choice(len(nodes),
                   size=min(n_samples, len(nodes)), replace=False)]

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


def estimate_growth_dimension(mean_volumes, maxR, rMin_frac=0.25, rMax_frac=0.67):
    """Estimate growth dimension via log-log fit and point-by-point delta."""
    rMin = max(3, int(maxR * rMin_frac))
    rMax = int(maxR * rMax_frac)

    radii = np.arange(rMin, rMax + 1, dtype=float)
    log_r = np.log(radii)
    log_v = np.log(mean_volumes[rMin:rMax + 1].astype(float))
    coeffs = np.polyfit(log_r, log_v, 1)
    slope = coeffs[0]

    predicted = np.polyval(coeffs, log_r)
    ss_res = np.sum((log_v - predicted) ** 2)
    ss_tot = np.sum((log_v - np.mean(log_v)) ** 2)
    r_squared = 1.0 - ss_res / ss_tot

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
        'deltas': deltas,
        'mean_volumes': mean_volumes,
    }


def compute_spectral_dimension_eigsh(G, n_eigenvalues=300):
    """Estimate spectral dimension from Laplacian eigenvalues using sparse solver."""
    A = nx.adjacency_matrix(G).astype(float)
    n = G.number_of_nodes()
    degrees = np.array([G.degree(v) for v in G.nodes()])
    D = diags(degrees.astype(float), 0, format='csr')
    L = D - A

    k = min(n_eigenvalues + 1, n - 2)

    try:
        eigenvalues = eigsh(L, k=k, sigma=0.01, return_eigenvectors=False)
    except Exception:
        print("    shift-invert failed, trying which='SM'...")
        eigenvalues = eigsh(L, k=k, which='SM', return_eigenvectors=False,
                           maxiter=5000, tol=1e-6)

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


def torus_eigenvalues_exact(side_n, dim):
    """Compute exact Laplacian eigenvalues for a d-dimensional torus."""
    freqs = 2.0 * np.pi * np.arange(side_n) / side_n
    cosines = 1.0 - np.cos(freqs)

    if dim == 2:
        eigs = 2.0 * (cosines[:, None] + cosines[None, :]).ravel()
    elif dim == 3:
        eigs = 2.0 * (cosines[:, None, None] + cosines[None, :, None]
                       + cosines[None, None, :]).ravel()
    elif dim == 4:
        eigs = 2.0 * (cosines[:, None, None, None] + cosines[None, :, None, None]
                       + cosines[None, None, :, None] + cosines[None, None, None, :]).ravel()
    else:
        raise ValueError(f"Unsupported dimension {dim}")

    eigs = np.sort(eigs)
    eigs = eigs[eigs > 1e-10]
    return eigs


def compute_spectral_dimension_analytical(side_n, dim, n_eigenvalues=500):
    """Spectral dimension from analytical torus eigenvalues."""
    all_eigs = torus_eigenvalues_exact(side_n, dim)
    eigenvalues = all_eigs[:n_eigenvalues]

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


def compute_fisher_single_vertex(G, v0, nodes, sigma=3.0):
    """Compute Fisher Information for a single vertex. Returns sv_norm, gap_rank, participation_ratio."""
    n = len(nodes)
    dists_v0 = bfs_distances(G, v0)
    dist_arr_v0 = np.array([dists_v0[u] for u in nodes], dtype=float)
    log_unnorm_v0 = -dist_arr_v0 / sigma
    log_Z_v0 = np.logaddexp.reduce(log_unnorm_v0)
    log_p_v0 = log_unnorm_v0 - log_Z_v0
    p_v0 = np.exp(log_p_v0)

    neighbors = list(G.neighbors(v0))
    k = len(neighbors)

    score_vectors = np.zeros((k, n))
    for j, w in enumerate(neighbors):
        dists_w = bfs_distances(G, w)
        dist_arr_w = np.array([dists_w[u] for u in nodes], dtype=float)
        log_unnorm_w = -dist_arr_w / sigma
        log_Z_w = np.logaddexp.reduce(log_unnorm_w)
        log_p_w = log_unnorm_w - log_Z_w
        score_vectors[j, :] = log_p_w - log_p_v0

    weighted_scores = score_vectors * np.sqrt(p_v0)[np.newaxis, :]
    F = weighted_scores @ weighted_scores.T

    sv = np.linalg.svd(F, compute_uv=False)
    sv_norm = sv / sv[0] if sv[0] > 0 else sv

    # Gap-based rank
    if len(sv_norm) > 1 and sv_norm[0] > 0:
        ratios = sv_norm[1:] / np.maximum(sv_norm[:-1], 1e-15)
        gap_rank = int(np.argmin(ratios) + 1)
    else:
        gap_rank = 1

    # Participation ratio
    if np.sum(sv ** 2) > 0:
        participation_ratio = (np.sum(sv)) ** 2 / np.sum(sv ** 2)
    else:
        participation_ratio = 0.0

    # All gap ratios
    gap_ratios = list(ratios) if len(sv_norm) > 1 else []

    return {
        'sv_norm': sv_norm,
        'gap_rank': gap_rank,
        'participation_ratio': participation_ratio,
        'gap_ratios': gap_ratios,
    }


def compute_fisher_sampled(G, sigma=3.0, n_samples=20, sample_nodes=None):
    """Compute Fisher dimension over sampled vertices."""
    nodes = list(G.nodes())
    n = len(nodes)

    if sample_nodes is not None:
        samples = sample_nodes[:n_samples]
    else:
        samples = [nodes[i] for i in np.random.choice(n, size=n_samples, replace=False)]

    all_results = []
    for s_idx, v0 in enumerate(samples):
        if (s_idx + 1) % 5 == 0:
            print(f"    Fisher sample {s_idx+1}/{len(samples)}")
        r = compute_fisher_single_vertex(G, v0, nodes, sigma)
        all_results.append(r)

    ranks = [r['gap_rank'] for r in all_results]
    prs = [r['participation_ratio'] for r in all_results]
    svs = [r['sv_norm'] for r in all_results]

    max_len = max(len(sv) for sv in svs)
    padded = np.zeros((len(svs), max_len))
    for i, sv in enumerate(svs):
        padded[i, :len(sv)] = sv
    mean_sv = np.mean(padded, axis=0)
    std_sv = np.std(padded, axis=0)

    return {
        'fisher_dimension_mean': np.mean(ranks),
        'fisher_dimension_median': np.median(ranks),
        'fisher_participation_ratio': np.mean(prs),
        'fisher_participation_std': np.std(prs),
        'mean_singular_values': mean_sv,
        'std_singular_values': std_sv,
        'individual_ranks': ranks,
        'individual_prs': prs,
        'individual_svs': svs,
        'sigma': sigma,
    }


def compute_fisher_all_vertices(G, sigma=3.0, vertex_filter=None):
    """Compute Fisher at every vertex (or filtered subset). Returns per-vertex results."""
    nodes = list(G.nodes())
    if vertex_filter is not None:
        targets = vertex_filter
    else:
        targets = nodes

    all_results = {}
    for i, v in enumerate(targets):
        if (i + 1) % 100 == 0:
            print(f"    Fisher all-vertex: {i+1}/{len(targets)}")
        r = compute_fisher_single_vertex(G, v, nodes, sigma)
        all_results[v] = r

    prs = [r['participation_ratio'] for r in all_results.values()]
    ranks = [r['gap_rank'] for r in all_results.values()]

    return {
        'per_vertex': all_results,
        'all_prs': prs,
        'all_ranks': ranks,
        'mean_pr': np.mean(prs),
        'std_pr': np.std(prs),
        'mean_rank': np.mean(ranks),
    }


# ============================================================================
# Sierpinski Gasket Construction
# ============================================================================

def make_sierpinski_gasket(level):
    """Build Sierpinski gasket graph at given level, recursively.

    Returns (G, corners) where corners = (top, bottom_left, bottom_right).
    """
    if level == 0:
        G = nx.Graph()
        G.add_edges_from([(0, 1), (1, 2), (0, 2)])
        return G, (0, 1, 2)

    # Get previous level
    G_prev, (top_p, bl_p, br_p) = make_sierpinski_gasket(level - 1)
    m = G_prev.number_of_nodes()

    # Create three copies
    GA = G_prev.copy()
    GB = nx.relabel_nodes(G_prev, {v: v + m for v in G_prev.nodes()})
    GC = nx.relabel_nodes(G_prev, {v: v + 2 * m for v in G_prev.nodes()})

    # Corners for each copy
    # A: top of triangle -> corners same as original
    a_top, a_bl, a_br = top_p, bl_p, br_p
    # B: bottom-left -> offset by m
    b_top, b_bl, b_br = top_p + m, bl_p + m, br_p + m
    # C: bottom-right -> offset by 2m
    c_top, c_bl, c_br = top_p + 2 * m, bl_p + 2 * m, br_p + 2 * m

    # Merge into one graph
    G = nx.compose(nx.compose(GA, GB), GC)

    # Merge corner pairs:
    # A's bottom-left = B's top
    G = nx.contracted_nodes(G, a_bl, b_top, self_loops=False)
    # A's bottom-right = C's top
    G = nx.contracted_nodes(G, a_br, c_top, self_loops=False)
    # B's bottom-right = C's bottom-left
    G = nx.contracted_nodes(G, b_br, c_bl, self_loops=False)

    # New corners of SG(L)
    new_top = a_top
    new_bl = b_bl
    new_br = c_br

    # Relabel to contiguous integers
    mapping = {old: new for new, old in enumerate(sorted(G.nodes()))}
    G = nx.relabel_nodes(G, mapping)
    new_top = mapping[new_top]
    new_bl = mapping[new_bl]
    new_br = mapping[new_br]

    return G, (new_top, new_bl, new_br)


def validate_sierpinski(G, corners, level):
    """Validate Sierpinski gasket topology."""
    n_nodes = G.number_of_nodes()
    expected = (3 ** (level + 1) + 3) // 2
    assert n_nodes == expected, f"Expected {expected} vertices, got {n_nodes}"

    degrees = dict(G.degree())
    deg2 = [v for v, d in degrees.items() if d == 2]
    deg4 = [v for v, d in degrees.items() if d == 4]
    other = [v for v, d in degrees.items() if d not in (2, 4)]
    assert len(deg2) == 3, f"Expected 3 degree-2 vertices (corners), got {len(deg2)}"
    assert len(other) == 0, f"Found vertices with unexpected degrees: {set(degrees[v] for v in other)}"
    assert set(deg2) == set(corners), f"Degree-2 vertices should be corners"

    assert nx.is_connected(G), "Graph should be connected"

    diameter = nx.diameter(G)
    expected_diam = 2 ** level
    assert abs(diameter - expected_diam) <= 1, \
        f"Expected diameter ~{expected_diam}, got {diameter}"

    print(f"  Validated: {n_nodes} vertices, {len(deg2)} corners (deg 2), "
          f"{len(deg4)} interior (deg 4), diameter {diameter}")
    return deg4


# ============================================================================
# 4D Torus Construction
# ============================================================================

def make_torus_4d(n):
    """Build 4D flat torus. Every vertex has degree 8."""
    G = nx.Graph()
    n2 = n * n
    n3 = n * n * n
    total = n ** 4
    G.add_nodes_from(range(total))

    for i1 in range(n):
        for i2 in range(n):
            for i3 in range(n):
                for i4 in range(n):
                    v = i1 * n3 + i2 * n2 + i3 * n + i4
                    # 4 neighbor pairs (±1 in each dimension)
                    G.add_edge(v, i1 * n3 + i2 * n2 + i3 * n + (i4 + 1) % n)
                    G.add_edge(v, i1 * n3 + i2 * n2 + ((i3 + 1) % n) * n + i4)
                    G.add_edge(v, i1 * n3 + ((i2 + 1) % n) * n2 + i3 * n + i4)
                    G.add_edge(v, ((i1 + 1) % n) * n3 + i2 * n2 + i3 * n + i4)
    return G


def validate_torus_4d(G, n):
    """Validate 4D torus topology."""
    n_nodes = G.number_of_nodes()
    assert n_nodes == n ** 4, f"Expected {n**4} nodes, got {n_nodes}"
    degrees = set(dict(G.degree()).values())
    assert degrees == {8}, f"Expected all degree 8, got {degrees}"
    print(f"  Validated: {n_nodes} nodes, all degree 8")


# ============================================================================
# Plotting Functions
# ============================================================================

def plot_growth_delta(deltas, ref_dims, rMin, rMax, prefix, output_dir, title=None):
    """Plot delta(r) vs r with reference lines."""
    fig, ax = plt.subplots(figsize=(10, 6))
    rs = [r for r, d in deltas]
    ds = [d for r, d in deltas]
    ax.plot(rs, ds, 'b.-', alpha=0.7, markersize=4)
    colors = ['r', 'g', 'purple']
    for i, (label, val) in enumerate(ref_dims.items()):
        ax.axhline(y=val, color=colors[i % len(colors)], linestyle='--',
                   linewidth=2, label=f'{label} = {val:.4f}')
    ax.axvspan(rMin, rMax, alpha=0.15, color='green', label=f'Fitting range [{rMin}, {rMax}]')
    ax.set_xlabel('Radius r', fontsize=12)
    ax.set_ylabel('Local dimension Δ(r)', fontsize=12)
    ax.set_title(title or f'Growth Dimension Δ(r) — {prefix}', fontsize=14)
    ax.legend(fontsize=10)
    ax.grid(True, alpha=0.3)
    path = os.path.join(output_dir, f'{prefix}_growth_delta.png')
    fig.savefig(path, dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f"    Saved {path}")


def plot_growth_loglog(mean_volumes, rMin, rMax, slope, ref_dims, prefix, output_dir, title=None):
    """Plot log V(r) vs log r with fitted line."""
    fig, ax = plt.subplots(figsize=(10, 6))
    valid = np.arange(1, len(mean_volumes))
    log_r = np.log(valid.astype(float))
    log_v = np.log(mean_volumes[1:].astype(float))
    ax.plot(log_r, log_v, 'b.', alpha=0.5, markersize=3)

    fit_r = np.arange(rMin, rMax + 1, dtype=float)
    fit_log_r = np.log(fit_r)
    fit_log_v = np.log(mean_volumes[rMin:rMax + 1].astype(float))
    coeffs = np.polyfit(fit_log_r, fit_log_v, 1)
    first_ref = list(ref_dims.values())[0]
    ax.plot(fit_log_r, np.polyval(coeffs, fit_log_r), 'r-', linewidth=2,
            label=f'Fit: slope = {slope:.4f} (ref = {first_ref:.4f})')

    ax.set_xlabel('log(r)', fontsize=12)
    ax.set_ylabel('log(V(r))', fontsize=12)
    ax.set_title(title or f'Growth Dimension log-log — {prefix}', fontsize=14)
    ax.legend(fontsize=11)
    ax.grid(True, alpha=0.3)
    path = os.path.join(output_dir, f'{prefix}_growth_loglog.png')
    fig.savefig(path, dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f"    Saved {path}")


def plot_spectral_weyl(eigenvalues, fit_start, fit_end, weyl_slope, spectral_dim,
                       ref_dims, prefix, output_dir, title=None):
    """Plot Weyl law fit."""
    fig, ax = plt.subplots(figsize=(10, 6))
    n_eigs = len(eigenvalues)
    indices = np.arange(1, n_eigs + 1, dtype=float)
    ax.plot(np.log(eigenvalues), np.log(indices), 'b.', alpha=0.4, markersize=3)

    log_lam = np.log(eigenvalues[fit_start:fit_end])
    log_N = np.log(indices[fit_start:fit_end])
    coeffs = np.polyfit(log_lam, log_N, 1)
    first_ref_name = list(ref_dims.keys())[0]
    first_ref_val = list(ref_dims.values())[0]
    ax.plot(log_lam, np.polyval(coeffs, log_lam), 'r-', linewidth=2,
            label=f'Weyl slope={weyl_slope:.4f} → d={spectral_dim:.4f} (ref {first_ref_name}={first_ref_val:.4f})')

    ax.set_xlabel('log(λ)', fontsize=12)
    ax.set_ylabel('log(N(λ))', fontsize=12)
    ax.set_title(title or f'Spectral Dimension (Weyl Law) — {prefix}', fontsize=14)
    ax.legend(fontsize=10)
    ax.grid(True, alpha=0.3)
    path = os.path.join(output_dir, f'{prefix}_spectral_weyl.png')
    fig.savefig(path, dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f"    Saved {path}")


def plot_fisher_sv(mean_sv, std_sv, ref_dims, prefix, output_dir, title=None):
    """Bar chart of Fisher singular values with error bars."""
    fig, ax = plt.subplots(figsize=(8, 6))
    x = np.arange(1, len(mean_sv) + 1)
    ax.bar(x, mean_sv, color='#2196F3', edgecolor='black', linewidth=0.5,
           yerr=std_sv if std_sv is not None else None, capsize=4, ecolor='red')
    ax.axhline(y=0.01, color='gray', linestyle=':', alpha=0.5, label='Threshold 0.01')
    ax.set_xlabel('Singular Value Index', fontsize=12)
    ax.set_ylabel('Normalized Singular Value', fontsize=12)
    ax.set_title(title or f'Fisher Information Singular Values — {prefix}', fontsize=14)
    ax.set_xticks(x)
    ax.legend(fontsize=10)
    ax.grid(True, alpha=0.3, axis='y')
    path = os.path.join(output_dir, f'{prefix}_fisher_sv.png')
    fig.savefig(path, dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f"    Saved {path}")


def plot_fisher_sv_sigma_sweep(sweep_results, prefix, output_dir):
    """Panel plot of Fisher SV profiles at multiple sigma values."""
    sigmas = sorted(sweep_results.keys())
    n_panels = len(sigmas)
    cols = 3
    rows = (n_panels + cols - 1) // cols
    fig, axes = plt.subplots(rows, cols, figsize=(5 * cols, 4 * rows))
    axes = np.atleast_2d(axes)

    for idx, sigma in enumerate(sigmas):
        row, col = idx // cols, idx % cols
        ax = axes[row, col]
        sv = sweep_results[sigma]['mean_sv']
        std = sweep_results[sigma].get('std_sv', None)
        x = np.arange(1, len(sv) + 1)
        ax.bar(x, sv, color='#2196F3', edgecolor='black', linewidth=0.5,
               yerr=std, capsize=3, ecolor='red')
        ax.set_title(f'σ = {sigma}', fontsize=11)
        ax.set_ylim(0, 1.15)
        ax.set_xticks(x)
        ax.grid(True, alpha=0.3, axis='y')

    # Turn off unused axes
    for idx in range(n_panels, rows * cols):
        row, col = idx // cols, idx % cols
        axes[row, col].set_visible(False)

    fig.suptitle(f'Fisher SV Profiles — Sigma Sweep — {prefix}', fontsize=14)
    fig.tight_layout()
    path = os.path.join(output_dir, f'{prefix}_fisher_sv_sigma_sweep.png')
    fig.savefig(path, dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f"    Saved {path}")


def plot_participation_histogram(prs, ref_dims, prefix, output_dir, title=None):
    """Histogram of participation ratios with reference dimension lines."""
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.hist(prs, bins=30, color='#2196F3', edgecolor='black', alpha=0.8, density=True)
    colors = ['r', 'g', 'purple']
    for i, (label, val) in enumerate(ref_dims.items()):
        ax.axvline(x=val, color=colors[i % len(colors)], linestyle='--',
                   linewidth=2, label=f'{label} = {val:.4f}')
    mean_pr = np.mean(prs)
    ax.axvline(x=mean_pr, color='black', linestyle='-', linewidth=2,
               label=f'Mean PR = {mean_pr:.4f}')
    ax.set_xlabel('Participation Ratio', fontsize=12)
    ax.set_ylabel('Density', fontsize=12)
    ax.set_title(title or f'Fisher Participation Ratio Distribution — {prefix}', fontsize=14)
    ax.legend(fontsize=10)
    ax.grid(True, alpha=0.3)
    path = os.path.join(output_dir, f'{prefix}_fisher_participation_hist.png')
    fig.savefig(path, dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f"    Saved {path}")


def plot_convergence_fractal(estimates, ref_dims, prefix, output_dir, title=None):
    """Bar chart comparing route estimates with multiple reference lines."""
    fig, ax = plt.subplots(figsize=(9, 6))
    labels = list(estimates.keys())
    vals = list(estimates.values())
    ax.bar(labels, vals, color='#4CAF50', edgecolor='black', linewidth=0.8, width=0.5)

    colors = ['black', 'red', 'purple']
    for i, (label, val) in enumerate(ref_dims.items()):
        ax.axhline(y=val, color=colors[i % len(colors)], linestyle='--',
                   linewidth=2, label=f'{label} = {val:.4f}')

    for i, (label, val) in enumerate(zip(labels, vals)):
        ax.text(i, val + 0.02, f'{val:.3f}', ha='center', va='bottom',
                fontsize=11, fontweight='bold')

    ax.set_ylabel('Dimension Estimate', fontsize=12)
    ax.set_title(title or f'Convergence Summary — {prefix}', fontsize=14)
    ax.legend(fontsize=10, loc='upper right')
    ax.grid(True, alpha=0.3, axis='y')
    path = os.path.join(output_dir, f'{prefix}_convergence.png')
    fig.savefig(path, dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f"    Saved {path}")


# ============================================================================
# Sierpinski Gasket Runner
# ============================================================================

def run_sierpinski(level):
    """Run all three routes on Sierpinski gasket at given level."""
    prefix = f'sierpinski_L{level}'
    print(f"\n{'='*70}")
    print(f"  Sierpinski Gasket Level {level}")
    print(f"{'='*70}")

    results = {'level': level, 'prefix': prefix}

    # --- Build ---
    print(f"\n[Step 0] Building Sierpinski gasket level {level}...")
    t0 = time.time()
    G, corners = make_sierpinski_gasket(level)
    interior = validate_sierpinski(G, corners, level)
    print(f"  Construction time: {time.time()-t0:.1f}s")
    results['n_nodes'] = G.number_of_nodes()
    results['n_interior'] = len(interior)

    # --- Route 1: Growth Dimension ---
    print(f"\n[Route 1] Growth Dimension...")
    t0 = time.time()
    sample_interior = [interior[i] for i in
                       np.random.choice(len(interior), size=min(25, len(interior)), replace=False)]
    mean_volumes, maxR = compute_ball_volumes(G, n_samples=25, sample_nodes=sample_interior)
    growth = estimate_growth_dimension(mean_volumes, maxR, rMin_frac=0.20, rMax_frac=0.60)
    growth_time = time.time() - t0

    print(f"  Growth dimension:  {growth['fit_dimension']:.4f}")
    print(f"  Delta dimension:   {growth['delta_dimension']:.4f}")
    print(f"  R²:                {growth['r_squared']:.6f}")
    print(f"  Fit range:         [{growth['rMin']}, {growth['rMax']}] (maxR={maxR})")
    gate_lo, gate_hi = 1.40, 1.75
    gate_pass = gate_lo <= growth['fit_dimension'] <= gate_hi
    print(f"  Gate [{gate_lo}, {gate_hi}]: {'PASS' if gate_pass else 'FAIL'}")
    print(f"  Time: {growth_time:.1f}s")
    results['growth'] = growth
    results['growth_pass'] = gate_pass

    ref_dims = {'d_H': D_HAUSDORFF, 'd_S': D_SPECTRAL}
    plot_growth_delta(growth['deltas'], ref_dims, growth['rMin'], growth['rMax'],
                      prefix, OUTPUT_DIR)
    plot_growth_loglog(mean_volumes, growth['rMin'], growth['rMax'],
                       growth['fit_dimension'], ref_dims, prefix, OUTPUT_DIR)

    # --- Route 2: Spectral Dimension ---
    print(f"\n[Route 2] Spectral Dimension...")
    t0 = time.time()
    n_eigs = min(300, G.number_of_nodes() - 2)
    spectral = compute_spectral_dimension_eigsh(G, n_eigenvalues=n_eigs)
    spectral_time = time.time() - t0

    print(f"  Spectral dimension: {spectral['spectral_dimension']:.4f}")
    print(f"  Weyl slope:         {spectral['weyl_slope']:.4f}")
    print(f"  R²:                 {spectral['r_squared']:.6f}")
    gate_lo_s, gate_hi_s = 1.20, 1.55
    gate_pass_s = gate_lo_s <= spectral['spectral_dimension'] <= gate_hi_s
    print(f"  Gate [{gate_lo_s}, {gate_hi_s}]: {'PASS' if gate_pass_s else 'FAIL'}")
    print(f"  Time: {spectral_time:.1f}s")
    results['spectral'] = spectral
    results['spectral_pass'] = gate_pass_s

    plot_spectral_weyl(spectral['eigenvalues'], spectral['fit_range'][0],
                       spectral['fit_range'][1], spectral['weyl_slope'],
                       spectral['spectral_dimension'], ref_dims, prefix, OUTPUT_DIR)

    # --- Route 3: Fisher Information ---
    print(f"\n[Route 3] Fisher Information (sampled)...")
    t0 = time.time()
    fisher = compute_fisher_sampled(G, sigma=3.0, n_samples=min(25, len(interior)),
                                     sample_nodes=sample_interior)
    fisher_time = time.time() - t0

    print(f"  Fisher rank (mean):       {fisher['fisher_dimension_mean']:.4f}")
    print(f"  Fisher rank (median):     {fisher['fisher_dimension_median']:.4f}")
    print(f"  Participation ratio:      {fisher['fisher_participation_ratio']:.4f} "
          f"± {fisher['fisher_participation_std']:.4f}")
    print(f"  Mean SVs:                 {np.array2string(fisher['mean_singular_values'], precision=4)}")
    print(f"  SV std:                   {np.array2string(fisher['std_singular_values'], precision=4)}")
    print(f"  Individual ranks:         {fisher['individual_ranks']}")
    print(f"  Time: {fisher_time:.1f}s")
    results['fisher'] = fisher

    plot_fisher_sv(fisher['mean_singular_values'], fisher['std_singular_values'],
                   ref_dims, prefix, OUTPUT_DIR)

    # --- Fisher sigma sweep ---
    print(f"\n  [Fisher sigma sweep]")
    sweep = {}
    for sigma in [1.5, 2.0, 3.0, 5.0, 8.0, 12.0]:
        print(f"    sigma = {sigma}")
        r = compute_fisher_sampled(G, sigma=sigma, n_samples=min(15, len(interior)),
                                    sample_nodes=sample_interior)
        sweep[sigma] = {
            'mean_rank': r['fisher_dimension_mean'],
            'mean_pr': r['fisher_participation_ratio'],
            'std_pr': r['fisher_participation_std'],
            'mean_sv': r['mean_singular_values'],
            'std_sv': r['std_singular_values'],
        }
        print(f"      rank={sweep[sigma]['mean_rank']:.2f}, "
              f"PR={sweep[sigma]['mean_pr']:.4f} ± {sweep[sigma]['std_pr']:.4f}")
    results['fisher_sweep'] = sweep
    plot_fisher_sv_sigma_sweep(sweep, prefix, OUTPUT_DIR)

    # --- Fisher per-vertex map (for L=6 or small graphs) ---
    if G.number_of_nodes() <= 5000:
        print(f"\n  [Fisher per-vertex map] ({len(interior)} interior vertices)")
        t0 = time.time()
        per_vertex = compute_fisher_all_vertices(G, sigma=3.0, vertex_filter=interior)
        print(f"  Per-vertex Fisher: mean PR = {per_vertex['mean_pr']:.4f} "
              f"± {per_vertex['std_pr']:.4f}")
        print(f"  Per-vertex ranks: {np.unique(per_vertex['all_ranks'], return_counts=True)}")
        print(f"  Time: {time.time()-t0:.1f}s")
        results['per_vertex'] = per_vertex

        plot_participation_histogram(per_vertex['all_prs'], ref_dims, prefix, OUTPUT_DIR)

    # --- Convergence ---
    estimates = {
        'Route 1\nGrowth': growth['fit_dimension'],
        'Route 2\nSpectral': spectral['spectral_dimension'],
        'Route 3\nFisher PR': fisher['fisher_participation_ratio'],
    }
    results['estimates'] = estimates

    plot_convergence_fractal(estimates, ref_dims, prefix, OUTPUT_DIR)

    return results


# ============================================================================
# 4D Torus Runner
# ============================================================================

def run_torus_4d(n=12):
    """Run all three routes on 4D torus."""
    prefix = 'torus4d'
    print(f"\n{'='*70}")
    print(f"  4D Torus ({n}x{n}x{n}x{n})")
    print(f"{'='*70}")

    results = {'n': n, 'prefix': prefix}

    # --- Build ---
    print(f"\n[Step 0] Building 4D torus ({n}^4 = {n**4} vertices)...")
    t0 = time.time()
    G = make_torus_4d(n)
    validate_torus_4d(G, n)
    print(f"  Construction time: {time.time()-t0:.1f}s")

    # --- Route 1: Growth Dimension ---
    print(f"\n[Route 1] Growth Dimension...")
    t0 = time.time()
    mean_volumes, maxR = compute_ball_volumes(G, n_samples=15)
    growth = estimate_growth_dimension(mean_volumes, maxR)
    growth_time = time.time() - t0

    gate_lo, gate_hi = 3.70, 4.30
    gate_pass = gate_lo <= growth['fit_dimension'] <= gate_hi
    print(f"  Growth dimension:  {growth['fit_dimension']:.4f}")
    print(f"  R²:                {growth['r_squared']:.6f}")
    print(f"  Fit range:         [{growth['rMin']}, {growth['rMax']}]")
    print(f"  Gate [{gate_lo}, {gate_hi}]: {'PASS' if gate_pass else 'FAIL'}")
    print(f"  Time: {growth_time:.1f}s")
    results['growth'] = growth
    results['growth_pass'] = gate_pass

    ref_dims = {'True d': 4.0}
    plot_growth_delta(growth['deltas'], ref_dims, growth['rMin'], growth['rMax'],
                      prefix, OUTPUT_DIR)
    plot_growth_loglog(mean_volumes, growth['rMin'], growth['rMax'],
                       growth['fit_dimension'], ref_dims, prefix, OUTPUT_DIR)

    # --- Route 2: Spectral Dimension (analytical) ---
    print(f"\n[Route 2] Spectral Dimension (analytical)...")
    t0 = time.time()
    spectral = compute_spectral_dimension_analytical(n, 4, n_eigenvalues=500)
    spectral_time = time.time() - t0

    gate_lo_s, gate_hi_s = 3.70, 4.30
    gate_pass_s = gate_lo_s <= spectral['spectral_dimension'] <= gate_hi_s
    print(f"  Spectral dimension: {spectral['spectral_dimension']:.4f}")
    print(f"  R²:                 {spectral['r_squared']:.6f}")
    print(f"  Gate [{gate_lo_s}, {gate_hi_s}]: {'PASS' if gate_pass_s else 'FAIL'}")
    print(f"  Time: {spectral_time:.1f}s")
    results['spectral'] = spectral
    results['spectral_pass'] = gate_pass_s

    plot_spectral_weyl(spectral['eigenvalues'], spectral['fit_range'][0],
                       spectral['fit_range'][1], spectral['weyl_slope'],
                       spectral['spectral_dimension'], ref_dims, prefix, OUTPUT_DIR)

    # --- Route 3: Fisher Information ---
    print(f"\n[Route 3] Fisher Information Rank...")
    t0 = time.time()
    fisher = compute_fisher_sampled(G, sigma=3.0, n_samples=20)
    fisher_time = time.time() - t0

    gate_pass_f = fisher['fisher_dimension_mean'] == 4.0
    print(f"  Fisher rank (mean):   {fisher['fisher_dimension_mean']:.4f}")
    print(f"  Participation ratio:  {fisher['fisher_participation_ratio']:.4f}")
    print(f"  Mean SVs:             {np.array2string(fisher['mean_singular_values'], precision=4)}")
    print(f"  Individual ranks:     {fisher['individual_ranks']}")
    print(f"  Gate (rank=4): {'PASS' if gate_pass_f else 'FAIL'}")
    print(f"  Time: {fisher_time:.1f}s")
    results['fisher'] = fisher
    results['fisher_pass'] = gate_pass_f

    plot_fisher_sv(fisher['mean_singular_values'], fisher['std_singular_values'],
                   ref_dims, prefix, OUTPUT_DIR)

    # --- Fisher sigma sweep ---
    print(f"\n  [Fisher sigma sweep]")
    sweep = {}
    for sigma in [1.5, 2.0, 3.0, 5.0, 8.0]:
        print(f"    sigma = {sigma}")
        r = compute_fisher_sampled(G, sigma=sigma, n_samples=10)
        sweep[sigma] = {
            'mean_rank': r['fisher_dimension_mean'],
            'mean_pr': r['fisher_participation_ratio'],
            'mean_sv': r['mean_singular_values'],
        }
        print(f"      rank={sweep[sigma]['mean_rank']:.2f}, PR={sweep[sigma]['mean_pr']:.4f}")
    results['fisher_sweep'] = sweep

    # --- Convergence ---
    estimates = {
        'growth': growth['fit_dimension'],
        'spectral': spectral['spectral_dimension'],
        'fisher': fisher['fisher_dimension_mean'],
    }
    vals = list(estimates.values())
    max_diff = max(vals) - min(vals)
    agreement_gate = 0.20
    agreement_pass = max_diff < agreement_gate

    print(f"\n[Convergence]")
    print(f"  Route 1 (Growth):   {estimates['growth']:.4f}")
    print(f"  Route 2 (Spectral): {estimates['spectral']:.4f}")
    print(f"  Route 3 (Fisher):   {estimates['fisher']:.4f}")
    print(f"  Max pairwise diff:  {max_diff:.4f}")
    print(f"  Agreement gate (<{agreement_gate}): {'PASS' if agreement_pass else 'FAIL'}")
    results['estimates'] = estimates
    results['max_diff'] = max_diff
    results['agreement_pass'] = agreement_pass

    convergence_estimates = {
        'Route 1\nGrowth': estimates['growth'],
        'Route 2\nSpectral': estimates['spectral'],
        'Route 3\nFisher': estimates['fisher'],
    }
    gate_band = {'gate_lo': gate_lo, 'gate_hi': gate_hi}

    # Convergence plot
    fig, ax = plt.subplots(figsize=(8, 6))
    labels = list(convergence_estimates.keys())
    vals = list(convergence_estimates.values())
    colors_bar = ['#4CAF50' if gate_lo <= v <= gate_hi else '#F44336' for v in vals]
    bars = ax.bar(labels, vals, color=colors_bar, edgecolor='black', linewidth=0.8, width=0.5)
    ax.axhline(y=4.0, color='black', linestyle='-', linewidth=2, label='True d = 4.0')
    ax.axhspan(gate_lo, gate_hi, alpha=0.12, color='green',
               label=f'Gate [{gate_lo}, {gate_hi}]')
    for bar, val in zip(bars, vals):
        ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.02,
                f'{val:.3f}', ha='center', va='bottom', fontsize=12, fontweight='bold')
    ax.set_ylabel('Dimension Estimate', fontsize=12)
    ax.set_title('Convergence Summary — 4D Torus', fontsize=14)
    ax.legend(fontsize=11)
    ax.grid(True, alpha=0.3, axis='y')
    path = os.path.join(OUTPUT_DIR, 'torus4d_convergence.png')
    fig.savefig(path, dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f"    Saved {path}")

    return results


# ============================================================================
# Report Generation
# ============================================================================

def generate_extension_report(sierpinski_results, torus4d_results, total_time):
    """Generate PHASE1_EXTENSION_RESULTS.md."""
    lines = []
    lines.append("# DS Phase 1 Extension Results")
    lines.append(f"## Date: {datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    lines.append(f"## Runtime: {total_time:.1f} seconds")
    lines.append("")

    # Sierpinski
    lines.append("## Test 2A: Sierpinski Gasket")
    lines.append("")
    lines.append(f"Reference dimensions: d_H = {D_HAUSDORFF:.4f}, d_S = {D_SPECTRAL:.4f}")
    lines.append("")

    for sr in sierpinski_results:
        L = sr['level']
        lines.append(f"### Level {L} ({sr['n_nodes']} vertices, {sr['n_interior']} interior)")
        lines.append("")
        lines.append("| Route | Method | Estimate | Gate | Status |")
        lines.append("|-------|--------|----------|------|--------|")

        g = sr['growth']
        lines.append(f"| 1 | Growth Dimension | {g['fit_dimension']:.4f} (R²={g['r_squared']:.4f}) "
                     f"| [1.40, 1.75] | {'PASS' if sr['growth_pass'] else 'FAIL'} |")

        s = sr['spectral']
        lines.append(f"| 2 | Spectral Dimension | {s['spectral_dimension']:.4f} "
                     f"(R²={s['r_squared']:.4f}) | [1.20, 1.55] | {'PASS' if sr['spectral_pass'] else 'FAIL'} |")

        f = sr['fisher']
        lines.append(f"| 3 | Fisher PR (mean) | {f['fisher_participation_ratio']:.4f} "
                     f"± {f['fisher_participation_std']:.4f} | — | — |")
        lines.append(f"| 3 | Fisher Rank (mode) | {np.median(f['individual_ranks']):.0f} | — | — |")
        lines.append("")

        lines.append(f"Fisher SV profile: {np.array2string(f['mean_singular_values'], precision=4)}")
        lines.append(f"Fisher SV std:     {np.array2string(f['std_singular_values'], precision=4)}")
        lines.append("")

        if 'per_vertex' in sr:
            pv = sr['per_vertex']
            lines.append(f"Per-vertex Fisher (all {sr['n_interior']} interior): "
                        f"mean PR = {pv['mean_pr']:.4f} ± {pv['std_pr']:.4f}")
            ranks, counts = np.unique(pv['all_ranks'], return_counts=True)
            lines.append(f"Per-vertex rank distribution: {dict(zip(ranks.tolist(), counts.tolist()))}")
            lines.append("")

        # Sigma sweep
        lines.append("**Sigma sensitivity:**")
        lines.append("")
        lines.append("| sigma | Mean Rank | Mean PR | PR Std |")
        lines.append("|-------|-----------|---------|--------|")
        for sigma, sv in sr['fisher_sweep'].items():
            lines.append(f"| {sigma} | {sv['mean_rank']:.2f} | {sv['mean_pr']:.4f} | {sv['std_pr']:.4f} |")
        lines.append("")

    # 4D Torus
    lines.append("## Test 3A: 4D Torus")
    lines.append("")
    r4 = torus4d_results
    lines.append(f"### 4D Torus (n={r4['n']}, {r4['n']**4} vertices)")
    lines.append("")
    lines.append("| Route | Method | Estimate | Gate | Status |")
    lines.append("|-------|--------|----------|------|--------|")

    g4 = r4['growth']
    lines.append(f"| 1 | Growth Dimension | {g4['fit_dimension']:.4f} (R²={g4['r_squared']:.4f}) "
                 f"| [3.70, 4.30] | {'PASS' if r4['growth_pass'] else 'FAIL'} |")
    s4 = r4['spectral']
    lines.append(f"| 2 | Spectral Dimension | {s4['spectral_dimension']:.4f} "
                 f"(R²={s4['r_squared']:.4f}) | [3.70, 4.30] | {'PASS' if r4['spectral_pass'] else 'FAIL'} |")
    f4 = r4['fisher']
    lines.append(f"| 3 | Fisher Rank | {f4['fisher_dimension_mean']:.4f} | rank=4 | "
                 f"{'PASS' if r4['fisher_pass'] else 'FAIL'} |")
    lines.append("")
    lines.append(f"Cross-route agreement: {r4['max_diff']:.4f} "
                 f"(gate < 0.20) → {'PASS' if r4['agreement_pass'] else 'FAIL'}")
    lines.append("")
    lines.append(f"Fisher SV profile: {np.array2string(f4['mean_singular_values'], precision=4)}")
    lines.append("")

    # Sigma sweep
    lines.append("**4D Torus sigma sweep:**")
    lines.append("")
    lines.append("| sigma | Mean Rank | PR |")
    lines.append("|-------|-----------|-----|")
    for sigma, sv in r4['fisher_sweep'].items():
        lines.append(f"| {sigma} | {sv['mean_rank']:.2f} | {sv['mean_pr']:.4f} |")
    lines.append("")

    report = "\n".join(lines)
    path = os.path.join(OUTPUT_DIR, "PHASE1_EXTENSION_RESULTS.md")
    with open(path, 'w') as f:
        f.write(report)
    print(f"\nReport saved to {path}")
    return report


# ============================================================================
# Main
# ============================================================================

def main():
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    t_start = time.time()

    # --- Sierpinski Gasket ---
    sierpinski_results = []
    for level in [6, 7]:
        r = run_sierpinski(level)
        sierpinski_results.append(r)

    # --- 4D Torus ---
    torus4d_results = run_torus_4d(n=12)

    total_time = time.time() - t_start
    report = generate_extension_report(sierpinski_results, torus4d_results, total_time)

    print("\n" + "=" * 70)
    print("  PHASE 1 EXTENSION COMPLETE")
    print("=" * 70)
    print(report)


if __name__ == '__main__':
    main()
