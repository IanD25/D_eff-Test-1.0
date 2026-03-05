#!/usr/bin/env python3
"""
DS Phase 1 Extension v2: Extended Sigma Sweep + Sierpinski Carpet
==================================================================
Test A: Push sigma sweep on gasket L=7/L=8 to sigma=50 — does PR → d_H?
Test B: Sierpinski carpet (d_H ≈ 1.893) — does PR-vs-sigma generalize?

The MONEY QUESTION: if gasket PR → 1.585 and carpet PR → 1.893 at large sigma,
then Fisher PR recovers Hausdorff dimension on self-similar fractals.
"""

import os
import sys
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

OUTPUT_DIR = "phase1_results"
np.random.seed(42)

# Reference dimensions
GASKET_D_H = np.log(3) / np.log(2)        # ≈ 1.5849
GASKET_D_S = 2 * np.log(3) / np.log(5)    # ≈ 1.3652
CARPET_D_H = np.log(8) / np.log(3)        # ≈ 1.8928
CARPET_D_S = 1.805                          # numerical estimate
CARPET_D_W = 2.097

SIGMA_SWEEP = [1.5, 2.0, 3.0, 5.0, 8.0, 12.0, 16.0, 20.0, 25.0, 30.0, 40.0, 50.0]

# ============================================================================
# Shared Utilities (self-contained)
# ============================================================================

def bfs_distances(G, source):
    dist = {source: 0}
    queue = deque([source])
    while queue:
        v = queue.popleft()
        for w in G.neighbors(v):
            if w not in dist:
                dist[w] = dist[v] + 1
                queue.append(w)
    return dist


def compute_fisher_single_vertex(G, v0, nodes, sigma=3.0):
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

    # Check for degenerate FIM
    frob_norm = np.linalg.norm(F, 'fro')
    if frob_norm < 1e-12:
        return None  # FIM degenerate

    sv = np.linalg.svd(F, compute_uv=False)
    sv_norm = sv / sv[0] if sv[0] > 0 else sv

    if len(sv_norm) > 1 and sv_norm[0] > 0:
        ratios = sv_norm[1:] / np.maximum(sv_norm[:-1], 1e-15)
        gap_rank = int(np.argmin(ratios) + 1)
    else:
        gap_rank = 1

    if np.sum(sv ** 2) > 0:
        pr = (np.sum(sv)) ** 2 / np.sum(sv ** 2)
    else:
        pr = 0.0

    return {'sv_norm': sv_norm, 'gap_rank': gap_rank, 'participation_ratio': pr}


def fisher_sigma_sweep(G, sigmas, sample_nodes, n_samples=20):
    """Run Fisher PR at multiple sigma values. Returns dict of sigma -> stats."""
    nodes = list(G.nodes())
    samples = sample_nodes[:n_samples]
    results = {}

    for sigma in sigmas:
        prs, ranks, svs = [], [], []
        n_degenerate = 0
        for v in samples:
            r = compute_fisher_single_vertex(G, v, nodes, sigma)
            if r is None:
                n_degenerate += 1
            else:
                prs.append(r['participation_ratio'])
                ranks.append(r['gap_rank'])
                svs.append(r['sv_norm'])

        if len(prs) < len(samples) // 2:
            print(f"    sigma={sigma}: FIM DEGENERATE ({n_degenerate}/{len(samples)} samples)")
            results[sigma] = {'degenerate': True, 'n_degenerate': n_degenerate}
        else:
            max_len = max(len(sv) for sv in svs) if svs else 0
            if max_len > 0:
                padded = np.zeros((len(svs), max_len))
                for i, sv in enumerate(svs):
                    padded[i, :len(sv)] = sv
                mean_sv = np.mean(padded, axis=0)
            else:
                mean_sv = np.array([])

            results[sigma] = {
                'degenerate': False,
                'mean_pr': np.mean(prs),
                'std_pr': np.std(prs),
                'mean_rank': np.mean(ranks),
                'mean_sv': mean_sv,
                'n_degenerate': n_degenerate,
                'n_valid': len(prs),
            }
            deg_str = f" ({n_degenerate} degen)" if n_degenerate > 0 else ""
            print(f"    sigma={sigma}: PR={results[sigma]['mean_pr']:.4f} "
                  f"± {results[sigma]['std_pr']:.4f}, rank={results[sigma]['mean_rank']:.2f}{deg_str}")

    return results


def compute_fisher_all_vertices(G, sigma, vertex_filter):
    """Compute Fisher at every vertex in filter. Returns per-vertex results."""
    nodes = list(G.nodes())
    prs, ranks = [], []
    for i, v in enumerate(vertex_filter):
        if (i + 1) % 200 == 0:
            print(f"    Per-vertex: {i+1}/{len(vertex_filter)}")
        r = compute_fisher_single_vertex(G, v, nodes, sigma)
        if r is not None:
            prs.append(r['participation_ratio'])
            ranks.append(r['gap_rank'])
    return {
        'prs': prs, 'ranks': ranks,
        'mean_pr': np.mean(prs), 'std_pr': np.std(prs),
    }


# ============================================================================
# Growth & Spectral Dimension (for carpet)
# ============================================================================

def compute_ball_volumes(G, n_samples=30, sample_nodes=None):
    nodes = list(G.nodes())
    if sample_nodes is not None:
        samples = sample_nodes[:n_samples]
    else:
        samples = [nodes[i] for i in np.random.choice(len(nodes),
                   size=min(n_samples, len(nodes)), replace=False)]
    d0 = bfs_distances(G, samples[0])
    maxR = max(d0.values()) // 2
    all_volumes = []
    for i, v in enumerate(samples):
        dists = bfs_distances(G, v)
        dist_vals = np.array(list(dists.values()))
        volumes = np.array([np.sum(dist_vals <= r) for r in range(maxR + 1)])
        all_volumes.append(volumes)
        if (i + 1) % 10 == 0:
            print(f"    BFS {i+1}/{len(samples)}")
    return np.mean(all_volumes, axis=0), maxR


def estimate_growth_dimension(mean_volumes, maxR, rMin_frac=0.20, rMax_frac=0.60):
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
    return {
        'fit_dimension': slope, 'r_squared': r_squared,
        'rMin': rMin, 'rMax': rMax, 'maxR': maxR, 'deltas': deltas,
        'mean_volumes': mean_volumes,
    }


def compute_spectral_dimension_eigsh(G, n_eigenvalues=300):
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
        'spectral_dimension': spectral_dim, 'weyl_slope': slope,
        'r_squared': r_squared, 'n_eigenvalues_used': fit_end - fit_start,
        'fit_range': (fit_start, fit_end), 'eigenvalues': eigenvalues,
    }


# ============================================================================
# Sierpinski Gasket Construction (reused from extension v1)
# ============================================================================

def make_sierpinski_gasket(level):
    if level == 0:
        G = nx.Graph()
        G.add_edges_from([(0, 1), (1, 2), (0, 2)])
        return G, (0, 1, 2)
    G_prev, (top_p, bl_p, br_p) = make_sierpinski_gasket(level - 1)
    m = G_prev.number_of_nodes()
    GA = G_prev.copy()
    GB = nx.relabel_nodes(G_prev, {v: v + m for v in G_prev.nodes()})
    GC = nx.relabel_nodes(G_prev, {v: v + 2 * m for v in G_prev.nodes()})
    a_top, a_bl, a_br = top_p, bl_p, br_p
    b_top, b_bl, b_br = top_p + m, bl_p + m, br_p + m
    c_top, c_bl, c_br = top_p + 2 * m, bl_p + 2 * m, br_p + 2 * m
    G = nx.compose(nx.compose(GA, GB), GC)
    G = nx.contracted_nodes(G, a_bl, b_top, self_loops=False)
    G = nx.contracted_nodes(G, a_br, c_top, self_loops=False)
    G = nx.contracted_nodes(G, b_br, c_bl, self_loops=False)
    new_top, new_bl, new_br = a_top, b_bl, c_br
    mapping = {old: new for new, old in enumerate(sorted(G.nodes()))}
    G = nx.relabel_nodes(G, mapping)
    return G, (mapping[new_top], mapping[new_bl], mapping[new_br])


# ============================================================================
# Sierpinski Carpet Construction
# ============================================================================

def is_removed(x, y, L):
    """Check if position (x,y) is in a removed square at any level."""
    for k in range(L):
        dx = (x // (3**k)) % 3
        dy = (y // (3**k)) % 3
        if dx == 1 and dy == 1:
            return True
    return False


def make_sierpinski_carpet(L):
    """Build the Sierpinski carpet graph at level L."""
    side = 3 ** L
    surviving = set()
    for x in range(side):
        for y in range(side):
            if not is_removed(x, y, L):
                surviving.add((x, y))

    G = nx.Graph()
    for (x, y) in surviving:
        G.add_node((x, y))
        if (x + 1, y) in surviving:
            G.add_edge((x, y), (x + 1, y))
        if (x, y + 1) in surviving:
            G.add_edge((x, y), (x, y + 1))
    return G


def validate_carpet(G, L):
    """Validate Sierpinski carpet topology."""
    n_nodes = G.number_of_nodes()
    expected = 8 ** L
    assert n_nodes == expected, f"Expected {expected} vertices, got {n_nodes}"

    degrees = dict(G.degree())
    deg_counts = {}
    for v, d in degrees.items():
        deg_counts[d] = deg_counts.get(d, 0) + 1
    assert all(d in {1, 2, 3, 4} for d in deg_counts), \
        f"Unexpected degrees: {set(deg_counts.keys())}"
    assert nx.is_connected(G), "Graph should be connected"
    assert min(degrees.values()) >= 1, "Found isolated vertices"

    deg4 = [v for v, d in degrees.items() if d == 4]
    print(f"  Validated: {n_nodes} vertices, connected")
    print(f"  Degree distribution: {dict(sorted(deg_counts.items()))}")
    print(f"  Interior (deg 4): {len(deg4)} vertices")
    return deg4, deg_counts


# ============================================================================
# Plotting
# ============================================================================

def plot_pr_vs_sigma(sweep_data, title, filename, ref_lines):
    """Plot PR vs sigma with reference lines. sweep_data = dict of label -> sweep_results."""
    fig, ax = plt.subplots(figsize=(12, 7))
    colors = ['#2196F3', '#FF5722', '#4CAF50', '#9C27B0']
    markers = ['o', 's', '^', 'D']

    for i, (label, sweep) in enumerate(sweep_data.items()):
        sigmas, prs, errs = [], [], []
        for sigma in sorted(sweep.keys()):
            r = sweep[sigma]
            if not r.get('degenerate', False):
                sigmas.append(sigma)
                prs.append(r['mean_pr'])
                errs.append(r['std_pr'])
        ax.errorbar(sigmas, prs, yerr=errs, marker=markers[i % len(markers)],
                    color=colors[i % len(colors)], linewidth=2, markersize=6,
                    capsize=3, label=label)

    line_styles = ['-', '--', '-.', ':']
    line_colors = ['red', 'red', 'blue', 'blue']
    for i, (rlabel, rval) in enumerate(ref_lines.items()):
        ax.axhline(y=rval, color=line_colors[i % len(line_colors)],
                   linestyle=line_styles[i % len(line_styles)],
                   linewidth=1.5, alpha=0.8, label=f'{rlabel} = {rval:.4f}')

    ax.set_xscale('log')
    ax.set_xlabel('σ (log scale)', fontsize=13)
    ax.set_ylabel('Participation Ratio', fontsize=13)
    ax.set_title(title, fontsize=14)
    ax.legend(fontsize=10, loc='best')
    ax.grid(True, alpha=0.3)
    path = os.path.join(OUTPUT_DIR, filename)
    fig.savefig(path, dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f"    Saved {path}")


def plot_sv_panel(sweep, prefix, output_dir, title=None):
    """Panel plot of SV profiles across sigma values."""
    sigmas = sorted([s for s in sweep if not sweep[s].get('degenerate', False)])
    n = len(sigmas)
    cols = min(4, n)
    rows = (n + cols - 1) // cols
    fig, axes = plt.subplots(rows, cols, figsize=(4 * cols, 3.5 * rows))
    if rows == 1 and cols == 1:
        axes = np.array([[axes]])
    elif rows == 1:
        axes = axes[np.newaxis, :]
    elif cols == 1:
        axes = axes[:, np.newaxis]

    for idx, sigma in enumerate(sigmas):
        row, col = idx // cols, idx % cols
        ax = axes[row, col]
        sv = sweep[sigma]['mean_sv']
        x = np.arange(1, len(sv) + 1)
        ax.bar(x, sv, color='#2196F3', edgecolor='black', linewidth=0.5)
        ax.set_title(f'σ = {sigma}', fontsize=10)
        ax.set_ylim(0, 1.15)
        ax.set_xticks(x)
        ax.grid(True, alpha=0.3, axis='y')

    for idx in range(n, rows * cols):
        row, col = idx // cols, idx % cols
        axes[row, col].set_visible(False)

    fig.suptitle(title or f'SV Profiles — Sigma Sweep — {prefix}', fontsize=13)
    fig.tight_layout()
    path = os.path.join(output_dir, f'{prefix}_sv_profiles_full_sweep.png')
    fig.savefig(path, dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f"    Saved {path}")


def plot_participation_histogram(prs, ref_dims, prefix, output_dir, title=None):
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.hist(prs, bins=40, color='#2196F3', edgecolor='black', alpha=0.8, density=True)
    colors = ['r', 'g', 'purple']
    for i, (label, val) in enumerate(ref_dims.items()):
        ax.axvline(x=val, color=colors[i % len(colors)], linestyle='--',
                   linewidth=2, label=f'{label} = {val:.4f}')
    ax.axvline(x=np.mean(prs), color='black', linestyle='-', linewidth=2,
               label=f'Mean = {np.mean(prs):.4f}')
    ax.set_xlabel('Participation Ratio', fontsize=12)
    ax.set_ylabel('Density', fontsize=12)
    ax.set_title(title or f'PR Distribution — {prefix}', fontsize=14)
    ax.legend(fontsize=10)
    ax.grid(True, alpha=0.3)
    path = os.path.join(output_dir, f'{prefix}_fisher_participation_hist.png')
    fig.savefig(path, dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f"    Saved {path}")


def plot_growth_delta(deltas, ref_dims, rMin, rMax, prefix, output_dir):
    fig, ax = plt.subplots(figsize=(10, 6))
    rs = [r for r, d in deltas]
    ds = [d for r, d in deltas]
    ax.plot(rs, ds, 'b.-', alpha=0.7, markersize=4)
    colors = ['r', 'g', 'purple']
    for i, (label, val) in enumerate(ref_dims.items()):
        ax.axhline(y=val, color=colors[i % len(colors)], linestyle='--',
                   linewidth=2, label=f'{label} = {val:.4f}')
    ax.axvspan(rMin, rMax, alpha=0.15, color='green', label=f'Fit range [{rMin}, {rMax}]')
    ax.set_xlabel('Radius r', fontsize=12)
    ax.set_ylabel('Δ(r)', fontsize=12)
    ax.set_title(f'Growth Dimension — {prefix}', fontsize=14)
    ax.legend(fontsize=10)
    ax.grid(True, alpha=0.3)
    path = os.path.join(output_dir, f'{prefix}_growth_delta.png')
    fig.savefig(path, dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f"    Saved {path}")


def plot_growth_loglog(mean_volumes, rMin, rMax, slope, ref_dim, prefix, output_dir):
    fig, ax = plt.subplots(figsize=(10, 6))
    valid = np.arange(1, len(mean_volumes))
    ax.plot(np.log(valid.astype(float)), np.log(mean_volumes[1:].astype(float)), 'b.', alpha=0.5, markersize=3)
    fit_r = np.arange(rMin, rMax + 1, dtype=float)
    fit_log_r = np.log(fit_r)
    fit_log_v = np.log(mean_volumes[rMin:rMax + 1].astype(float))
    coeffs = np.polyfit(fit_log_r, fit_log_v, 1)
    ax.plot(fit_log_r, np.polyval(coeffs, fit_log_r), 'r-', linewidth=2,
            label=f'Fit: slope = {slope:.4f} (ref = {ref_dim:.4f})')
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
                       ref_dim, prefix, output_dir):
    fig, ax = plt.subplots(figsize=(10, 6))
    n_eigs = len(eigenvalues)
    indices = np.arange(1, n_eigs + 1, dtype=float)
    ax.plot(np.log(eigenvalues), np.log(indices), 'b.', alpha=0.4, markersize=3)
    log_lam = np.log(eigenvalues[fit_start:fit_end])
    log_N = np.log(indices[fit_start:fit_end])
    coeffs = np.polyfit(log_lam, log_N, 1)
    ax.plot(log_lam, np.polyval(coeffs, log_lam), 'r-', linewidth=2,
            label=f'slope={weyl_slope:.4f} → d={spectral_dim:.4f} (ref={ref_dim:.4f})')
    ax.set_xlabel('log(λ)', fontsize=12)
    ax.set_ylabel('log(N(λ))', fontsize=12)
    ax.set_title(f'Spectral Dimension (Weyl) — {prefix}', fontsize=14)
    ax.legend(fontsize=10)
    ax.grid(True, alpha=0.3)
    path = os.path.join(output_dir, f'{prefix}_spectral_weyl.png')
    fig.savefig(path, dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f"    Saved {path}")


def plot_fisher_sv(mean_sv, std_sv, prefix, output_dir, title=None):
    fig, ax = plt.subplots(figsize=(8, 6))
    x = np.arange(1, len(mean_sv) + 1)
    ax.bar(x, mean_sv, color='#2196F3', edgecolor='black', linewidth=0.5,
           yerr=std_sv, capsize=4, ecolor='red')
    ax.set_xlabel('SV Index', fontsize=12)
    ax.set_ylabel('Normalized SV', fontsize=12)
    ax.set_title(title or f'Fisher SV Profile — {prefix}', fontsize=14)
    ax.set_xticks(x)
    ax.grid(True, alpha=0.3, axis='y')
    path = os.path.join(output_dir, f'{prefix}_fisher_sv.png')
    fig.savefig(path, dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f"    Saved {path}")


def plot_convergence_fractal(estimates, ref_dims, prefix, output_dir):
    fig, ax = plt.subplots(figsize=(9, 6))
    labels = list(estimates.keys())
    vals = list(estimates.values())
    ax.bar(labels, vals, color='#4CAF50', edgecolor='black', linewidth=0.8, width=0.5)
    colors = ['black', 'red', 'purple']
    for i, (label, val) in enumerate(ref_dims.items()):
        ax.axhline(y=val, color=colors[i % len(colors)], linestyle='--',
                   linewidth=2, label=f'{label} = {val:.4f}')
    for i, val in enumerate(vals):
        ax.text(i, val + 0.02, f'{val:.3f}', ha='center', va='bottom', fontsize=11, fontweight='bold')
    ax.set_ylabel('Dimension Estimate', fontsize=12)
    ax.set_title(f'Convergence — {prefix}', fontsize=14)
    ax.legend(fontsize=10)
    ax.grid(True, alpha=0.3, axis='y')
    path = os.path.join(output_dir, f'{prefix}_convergence.png')
    fig.savefig(path, dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f"    Saved {path}")


# ============================================================================
# Test A: Extended Gasket Sigma Sweep
# ============================================================================

def run_gasket_sigma_sweep():
    """Extended sigma sweep on gasket L=7 and L=8."""
    print(f"\n{'='*70}")
    print("  TEST A: Extended Gasket Sigma Sweep")
    print(f"{'='*70}")

    results = {}

    for level in [7, 8]:
        print(f"\n--- Sierpinski Gasket L={level} ---")
        t0 = time.time()
        G, corners = make_sierpinski_gasket(level)
        n = G.number_of_nodes()
        expected = (3 ** (level + 1) + 3) // 2
        assert n == expected, f"Expected {expected}, got {n}"
        degrees = dict(G.degree())
        interior = [v for v, d in degrees.items() if d == 4]
        print(f"  Built: {n} vertices, {len(interior)} interior, {time.time()-t0:.1f}s")

        # Sample interior vertices
        sample_interior = [interior[i] for i in
                           np.random.choice(len(interior), size=min(20, len(interior)), replace=False)]

        # Extended sigma sweep
        print(f"\n  [Sigma sweep — {len(SIGMA_SWEEP)} values]")
        sweep = fisher_sigma_sweep(G, SIGMA_SWEEP, sample_interior, n_samples=20)
        results[f'L{level}'] = {'sweep': sweep, 'level': level, 'n_nodes': n}

        # SV panel plot
        plot_sv_panel(sweep, f'sierpinski_L{level}',  OUTPUT_DIR,
                      f'SV Profiles — Gasket L={level}')

        # Per-vertex Fisher at sigma=3.0 for L=8
        if level == 8:
            print(f"\n  [Per-vertex Fisher L=8 at sigma=3.0]")
            t0 = time.time()
            pv = compute_fisher_all_vertices(G, sigma=3.0, vertex_filter=interior)
            print(f"  Per-vertex: mean PR = {pv['mean_pr']:.4f} ± {pv['std_pr']:.4f}, "
                  f"time={time.time()-t0:.1f}s")
            ranks_u, ranks_c = np.unique(pv['ranks'], return_counts=True)
            print(f"  Rank distribution: {dict(zip(ranks_u.tolist(), ranks_c.tolist()))}")
            results['L8_pervertex'] = pv

            plot_participation_histogram(
                pv['prs'], {'d_H': GASKET_D_H, 'd_S': GASKET_D_S},
                'sierpinski_L8', OUTPUT_DIR,
                'PR Distribution — Gasket L=8 (all interior)')

    # PR-vs-sigma plot (gasket only)
    sweep_data = {}
    for key in ['L7', 'L8']:
        sweep_data[f'Gasket {key}'] = results[key]['sweep']

    plot_pr_vs_sigma(
        sweep_data,
        'PR vs σ — Sierpinski Gasket (Extended Sweep)',
        'sierpinski_pr_vs_sigma_extended.png',
        {'Gasket d_H': GASKET_D_H, 'Gasket d_S': GASKET_D_S}
    )

    return results


# ============================================================================
# Test B: Sierpinski Carpet
# ============================================================================

def run_sierpinski_carpet():
    """Build and test Sierpinski carpet at L=3 and L=4."""
    print(f"\n{'='*70}")
    print("  TEST B: Sierpinski Carpet")
    print(f"{'='*70}")

    results = {}

    for level in [3, 4]:
        prefix = f'carpet_L{level}'
        print(f"\n--- Sierpinski Carpet L={level} ---")
        t0 = time.time()
        G = make_sierpinski_carpet(level)
        deg4, deg_counts = validate_carpet(G, level)
        print(f"  Construction time: {time.time()-t0:.1f}s")
        results[f'L{level}'] = {'level': level, 'n_nodes': G.number_of_nodes(),
                                 'deg_counts': deg_counts, 'n_deg4': len(deg4)}

        # Sample interior (deg 4) vertices
        sample_deg4 = [deg4[i] for i in
                       np.random.choice(len(deg4), size=min(30, len(deg4)), replace=False)]

        if level == 4:
            # Full routes only on L=4

            # Route 1: Growth
            print(f"\n  [Route 1] Growth Dimension...")
            t0 = time.time()
            mean_vol, maxR = compute_ball_volumes(G, n_samples=30, sample_nodes=sample_deg4)
            growth = estimate_growth_dimension(mean_vol, maxR)
            gate_lo, gate_hi = 1.70, 2.10
            gpass = gate_lo <= growth['fit_dimension'] <= gate_hi
            print(f"  Growth dim: {growth['fit_dimension']:.4f} (R²={growth['r_squared']:.4f})")
            print(f"  Gate [{gate_lo}, {gate_hi}]: {'PASS' if gpass else 'FAIL'}")
            print(f"  Time: {time.time()-t0:.1f}s")
            results['L4']['growth'] = growth
            results['L4']['growth_pass'] = gpass

            ref_dims = {'d_H': CARPET_D_H, 'd_S': CARPET_D_S}
            plot_growth_delta(growth['deltas'], ref_dims, growth['rMin'],
                              growth['rMax'], prefix, OUTPUT_DIR)
            plot_growth_loglog(mean_vol, growth['rMin'], growth['rMax'],
                               growth['fit_dimension'], CARPET_D_H, prefix, OUTPUT_DIR)

            # Route 2: Spectral
            print(f"\n  [Route 2] Spectral Dimension...")
            t0 = time.time()
            spectral = compute_spectral_dimension_eigsh(G, n_eigenvalues=300)
            gate_lo_s, gate_hi_s = 1.60, 2.00
            spass = gate_lo_s <= spectral['spectral_dimension'] <= gate_hi_s
            print(f"  Spectral dim: {spectral['spectral_dimension']:.4f} "
                  f"(R²={spectral['r_squared']:.4f})")
            print(f"  Gate [{gate_lo_s}, {gate_hi_s}]: {'PASS' if spass else 'FAIL'}")
            print(f"  Time: {time.time()-t0:.1f}s")
            results['L4']['spectral'] = spectral
            results['L4']['spectral_pass'] = spass

            plot_spectral_weyl(spectral['eigenvalues'], spectral['fit_range'][0],
                               spectral['fit_range'][1], spectral['weyl_slope'],
                               spectral['spectral_dimension'], CARPET_D_S, prefix, OUTPUT_DIR)

        # Fisher sigma sweep (both levels)
        print(f"\n  [Fisher sigma sweep L={level}]")
        sweep = fisher_sigma_sweep(G, SIGMA_SWEEP, sample_deg4, n_samples=20)
        results[f'L{level}']['sweep'] = sweep

        if level == 4:
            # Fisher SV at sigma=3
            r3 = sweep.get(3.0, {})
            if not r3.get('degenerate', True):
                plot_fisher_sv(r3['mean_sv'], None, prefix, OUTPUT_DIR,
                               f'Fisher SV — Carpet L=4 (σ=3)')

            plot_sv_panel(sweep, prefix, OUTPUT_DIR, f'SV Profiles — Carpet L={level}')

            # Per-vertex Fisher
            print(f"\n  [Per-vertex Fisher L=4 at sigma=3.0]")
            t0 = time.time()
            pv = compute_fisher_all_vertices(G, sigma=3.0, vertex_filter=deg4)
            print(f"  Per-vertex (deg4): mean PR = {pv['mean_pr']:.4f} ± {pv['std_pr']:.4f}, "
                  f"time={time.time()-t0:.1f}s")
            ranks_u, ranks_c = np.unique(pv['ranks'], return_counts=True)
            print(f"  Rank distribution: {dict(zip(ranks_u.tolist(), ranks_c.tolist()))}")
            results['L4_pervertex'] = pv

            plot_participation_histogram(
                pv['prs'], {'d_H': CARPET_D_H, 'd_S': CARPET_D_S},
                prefix, OUTPUT_DIR,
                'PR Distribution — Carpet L=4 (deg-4 interior)')

            # PR-vs-sigma for carpet alone
            plot_pr_vs_sigma(
                {f'Carpet L={level}': sweep},
                f'PR vs σ — Carpet L={level}',
                f'{prefix}_pr_vs_sigma.png',
                {'Carpet d_H': CARPET_D_H, 'Carpet d_S': CARPET_D_S}
            )

            # Convergence plot
            pr3 = sweep[3.0]['mean_pr'] if not sweep[3.0].get('degenerate', True) else 0
            plot_convergence_fractal(
                {'Growth': growth['fit_dimension'],
                 'Spectral': spectral['spectral_dimension'],
                 'Fisher PR\n(σ=3)': pr3},
                {'d_H': CARPET_D_H, 'd_S': CARPET_D_S},
                prefix, OUTPUT_DIR)

    return results


# ============================================================================
# The Money Plot
# ============================================================================

def make_money_plot(gasket_sweep, carpet_sweep):
    """Gasket vs Carpet PR-vs-sigma overlay — the most important plot."""
    print(f"\n  [Money Plot: Gasket vs Carpet PR-vs-sigma]")

    sweep_data = {}
    # Use L=8 gasket if available, else L=7
    if 'L8' in gasket_sweep:
        sweep_data['Gasket (L=8)'] = gasket_sweep['L8']['sweep']
    elif 'L7' in gasket_sweep:
        sweep_data['Gasket (L=7)'] = gasket_sweep['L7']['sweep']
    if 'L4' in carpet_sweep:
        sweep_data['Carpet (L=4)'] = carpet_sweep['L4']['sweep']

    plot_pr_vs_sigma(
        sweep_data,
        'Fisher PR vs σ — Gasket vs Carpet',
        'fractal_pr_vs_sigma_comparison.png',
        {
            'Gasket d_H': GASKET_D_H,
            'Gasket d_S': GASKET_D_S,
            'Carpet d_H': CARPET_D_H,
            'Carpet d_S': CARPET_D_S,
        }
    )


# ============================================================================
# Report
# ============================================================================

def generate_report(gasket_results, carpet_results, total_time):
    lines = []
    lines.append("# DS Phase 1 Extension v2 Results")
    lines.append(f"## Date: {datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    lines.append(f"## Runtime: {total_time:.1f} seconds")
    lines.append("")

    # Gasket sigma sweep
    lines.append("## Test A: Extended Gasket Sigma Sweep")
    lines.append("")
    lines.append("| sigma | L=7 PR | L=7 std | L=8 PR | L=8 std |")
    lines.append("|-------|--------|---------|--------|---------|")
    for sigma in SIGMA_SWEEP:
        r7 = gasket_results.get('L7', {}).get('sweep', {}).get(sigma, {})
        r8 = gasket_results.get('L8', {}).get('sweep', {}).get(sigma, {})
        pr7 = f"{r7['mean_pr']:.4f}" if not r7.get('degenerate', True) else "DEGEN"
        std7 = f"{r7['std_pr']:.4f}" if not r7.get('degenerate', True) else "—"
        pr8 = f"{r8['mean_pr']:.4f}" if not r8.get('degenerate', True) else "DEGEN"
        std8 = f"{r8['std_pr']:.4f}" if not r8.get('degenerate', True) else "—"
        lines.append(f"| {sigma} | {pr7} | {std7} | {pr8} | {std8} |")
    lines.append("")
    lines.append(f"Reference: d_H = {GASKET_D_H:.4f}, d_S = {GASKET_D_S:.4f}")
    lines.append("")

    if 'L8_pervertex' in gasket_results:
        pv = gasket_results['L8_pervertex']
        lines.append(f"L=8 per-vertex Fisher (sigma=3.0): "
                     f"mean PR = {pv['mean_pr']:.4f} ± {pv['std_pr']:.4f}")
        prs = pv['prs']
        lines.append(f"  min={min(prs):.4f}, Q1={np.percentile(prs,25):.4f}, "
                     f"median={np.median(prs):.4f}, Q3={np.percentile(prs,75):.4f}, "
                     f"max={max(prs):.4f}")
        lines.append("")

    # Carpet
    lines.append("## Test B: Sierpinski Carpet")
    lines.append("")
    if 'L4' in carpet_results and 'growth' in carpet_results['L4']:
        g = carpet_results['L4']['growth']
        s = carpet_results['L4']['spectral']
        lines.append(f"### Carpet L=4 ({carpet_results['L4']['n_nodes']} vertices)")
        lines.append("")
        lines.append("| Route | Estimate | Gate | Status |")
        lines.append("|-------|----------|------|--------|")
        lines.append(f"| Growth | {g['fit_dimension']:.4f} (R²={g['r_squared']:.4f}) "
                     f"| [1.70, 2.10] | {'PASS' if carpet_results['L4']['growth_pass'] else 'FAIL'} |")
        lines.append(f"| Spectral | {s['spectral_dimension']:.4f} (R²={s['r_squared']:.4f}) "
                     f"| [1.60, 2.00] | {'PASS' if carpet_results['L4']['spectral_pass'] else 'FAIL'} |")
        lines.append("")

    # Carpet sigma sweep
    for lvl in ['L3', 'L4']:
        if lvl in carpet_results and 'sweep' in carpet_results[lvl]:
            lines.append(f"**Carpet {lvl} sigma sweep:**")
            lines.append("")
            lines.append("| sigma | PR | std |")
            lines.append("|-------|-----|-----|")
            for sigma in SIGMA_SWEEP:
                r = carpet_results[lvl]['sweep'].get(sigma, {})
                if not r.get('degenerate', True):
                    lines.append(f"| {sigma} | {r['mean_pr']:.4f} | {r['std_pr']:.4f} |")
                else:
                    lines.append(f"| {sigma} | DEGEN | — |")
            lines.append("")

    if 'L4_pervertex' in carpet_results:
        pv = carpet_results['L4_pervertex']
        lines.append(f"Carpet L=4 per-vertex (deg4, sigma=3.0): "
                     f"mean PR = {pv['mean_pr']:.4f} ± {pv['std_pr']:.4f}")
        lines.append("")

    lines.append(f"Reference: Carpet d_H = {CARPET_D_H:.4f}, d_S = {CARPET_D_S:.4f}")
    lines.append("")

    report = "\n".join(lines)
    path = os.path.join(OUTPUT_DIR, "PHASE1_EXTENSION_V2_RESULTS.md")
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

    gasket_results = run_gasket_sigma_sweep()
    carpet_results = run_sierpinski_carpet()
    make_money_plot(gasket_results, carpet_results)

    total_time = time.time() - t_start
    report = generate_report(gasket_results, carpet_results, total_time)

    print("\n" + "=" * 70)
    print("  PHASE 1 EXTENSION v2 COMPLETE")
    print("=" * 70)
    print(report)


if __name__ == '__main__':
    main()
