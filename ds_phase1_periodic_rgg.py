#!/usr/bin/env python3
"""
DS Phase 1: Periodic RGG + Boundary Analysis — Resolving the d+1 Puzzle
========================================================================
Test A: RGG on periodic domain (flat torus) — eliminates density gradient.
        If Fisher rank drops from d+1 to d, the bounded-domain density
        gradient caused the extra informative direction.
Test B: Interior-vs-boundary stratification on bounded RGG —
        spatial confirmation of the density-gradient hypothesis.
"""

import os
import time
import datetime
from collections import deque, Counter

import numpy as np
import networkx as nx
from scipy.spatial import cKDTree, KDTree
from scipy.sparse.linalg import eigsh
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

OUTPUT_DIR = "phase1_results"
np.random.seed(42)

# ============================================================================
# BFS — optimized using adjacency list
# ============================================================================

def bfs_distances_fast(adj, source, n):
    dist = np.full(n, -1, dtype=np.int32)
    dist[source] = 0
    queue = deque([source])
    while queue:
        v = queue.popleft()
        d = dist[v] + 1
        for w in adj[v]:
            if dist[w] == -1:
                dist[w] = d
                queue.append(w)
    return dist


def build_adjacency_list(G):
    nodes = list(G.nodes())
    node_to_idx = {v: i for i, v in enumerate(nodes)}
    n = len(nodes)
    adj = [[] for _ in range(n)]
    for u, v in G.edges():
        i, j = node_to_idx[u], node_to_idx[v]
        adj[i].append(j)
        adj[j].append(i)
    return adj, nodes, node_to_idx, n


# ============================================================================
# Graph construction
# ============================================================================

def make_periodic_rgg(n_points, d, r):
    """Build RGG with periodic BCs using cKDTree boxsize."""
    points = np.random.uniform(0, 1, size=(n_points, d))
    tree = cKDTree(points, boxsize=1.0)
    pairs = tree.query_pairs(r)
    G = nx.Graph()
    G.add_nodes_from(range(n_points))
    G.add_edges_from(pairs)
    return G, points


def make_bounded_rgg(n_points, d, r):
    """Build standard bounded RGG in [0,1]^d (no periodic BCs)."""
    points = np.random.uniform(0, 1, size=(n_points, d))
    tree = KDTree(points)
    pairs = tree.query_pairs(r)
    G = nx.Graph()
    G.add_nodes_from(range(n_points))
    G.add_edges_from(pairs)
    return G, points


def validate_graph(G, name):
    """Validate and report graph statistics."""
    n = G.number_of_nodes()
    m = G.number_of_edges()
    connected = nx.is_connected(G)

    if connected:
        G_lcc = G
        lcc_frac = 1.0
    else:
        components = sorted(nx.connected_components(G), key=len, reverse=True)
        lcc = components[0]
        lcc_frac = len(lcc) / n
        G_lcc = G.subgraph(lcc).copy()
        mapping = {v: i for i, v in enumerate(G_lcc.nodes())}
        G_lcc = nx.relabel_nodes(G_lcc, mapping)

    degrees_lcc = [d for _, d in G_lcc.degree()]
    clustering = nx.average_clustering(G_lcc)

    adj_tmp, _, _, n_tmp = build_adjacency_list(G_lcc)
    samples = np.random.choice(n_tmp, min(30, n_tmp), replace=False)
    max_dist = 0
    for s in samples:
        dists = bfs_distances_fast(adj_tmp, s, n_tmp)
        md = np.max(dists[dists >= 0])
        if md > max_dist:
            max_dist = md
    diameter = max_dist

    stats = {
        'name': name,
        'n_original': n,
        'n_lcc': G_lcc.number_of_nodes(),
        'm_lcc': G_lcc.number_of_edges(),
        'connected': connected,
        'lcc_frac': lcc_frac,
        'degree_min': min(degrees_lcc),
        'degree_max': max(degrees_lcc),
        'degree_mean': np.mean(degrees_lcc),
        'degree_median': np.median(degrees_lcc),
        'degree_std': np.std(degrees_lcc),
        'clustering': clustering,
        'diameter': diameter,
    }

    print(f"\n  {name}:")
    print(f"    Vertices: {stats['n_lcc']} (LCC {stats['lcc_frac']*100:.1f}%)")
    print(f"    Edges: {stats['m_lcc']}")
    print(f"    Connected: {stats['connected']}")
    print(f"    Degree: min={stats['degree_min']}, max={stats['degree_max']}, "
          f"mean={stats['degree_mean']:.1f}, median={stats['degree_median']:.0f}, "
          f"std={stats['degree_std']:.1f}")
    print(f"    Clustering: {stats['clustering']:.4f}")
    print(f"    Diameter: {stats['diameter']} (est.)")

    return G_lcc, stats


# ============================================================================
# Route 1: Growth dimension
# ============================================================================

def growth_dimension(adj, n, n_samples=30, label=""):
    samples = np.random.choice(n, min(n_samples, n), replace=False)
    all_vols = {}
    max_dists = []

    for s in samples:
        dists = bfs_distances_fast(adj, s, n)
        dists_valid = dists[dists >= 0]
        maxR = int(np.max(dists_valid))
        max_dists.append(maxR)
        for r in range(1, maxR + 1):
            vol = int(np.sum(dists_valid <= r))
            if r not in all_vols:
                all_vols[r] = []
            all_vols[r].append(vol)

    avg_maxR = int(np.mean(max_dists))
    rs = sorted(all_vols.keys())
    mean_vols = [np.mean(all_vols[r]) for r in rs]

    rMin = max(2, avg_maxR // 5)
    rMax = avg_maxR * 2 // 3
    mask = np.array([(r >= rMin and r <= rMax) for r in rs])
    if np.sum(mask) < 3:
        rMin, rMax = 2, max(rs) - 1
        mask = np.array([(r >= rMin and r <= rMax) for r in rs])

    log_r = np.log(np.array(rs, dtype=float))
    log_v = np.log(np.array(mean_vols, dtype=float))
    log_r_fit = log_r[mask]
    log_v_fit = log_v[mask]

    if len(log_r_fit) < 2:
        return {'dim': np.nan, 'R2': 0, 'maxR': avg_maxR, 'rs': rs, 'vols': mean_vols}

    coeffs = np.polyfit(log_r_fit, log_v_fit, 1)
    dim = coeffs[0]
    predicted = np.polyval(coeffs, log_r_fit)
    ss_res = np.sum((log_v_fit - predicted) ** 2)
    ss_tot = np.sum((log_v_fit - np.mean(log_v_fit)) ** 2)
    R2 = 1 - ss_res / ss_tot if ss_tot > 0 else 0

    print(f"    Growth dim {label}: {dim:.3f} (R²={R2:.4f}, [{rMin},{rMax}], maxR={avg_maxR})")
    return {'dim': dim, 'R2': R2, 'maxR': avg_maxR, 'rMin': rMin, 'rMax': rMax,
            'rs': rs, 'vols': mean_vols, 'coeffs': coeffs}


# ============================================================================
# Route 2: Spectral dimension
# ============================================================================

def spectral_dimension(G, n_eigs=300, label=""):
    n = G.number_of_nodes()
    L = nx.laplacian_matrix(G).astype(float)
    k_request = min(n_eigs, n - 2)
    try:
        eigenvalues = eigsh(L, k=k_request, which='SM', return_eigenvectors=False,
                           sigma=0.01, mode='normal', maxiter=5000)
        eigenvalues = np.sort(eigenvalues)
    except Exception:
        try:
            eigenvalues = eigsh(L, k=k_request, which='SM', return_eigenvectors=False, maxiter=5000)
            eigenvalues = np.sort(eigenvalues)
        except Exception as e2:
            print(f"    Spectral: failed ({e2})")
            return {'dim': np.nan, 'R2': 0}

    pos = eigenvalues[eigenvalues > 1e-8]
    if len(pos) < 10:
        return {'dim': np.nan, 'R2': 0}

    ranks = np.arange(1, len(pos) + 1)
    i_start = len(pos) // 4
    i_end = 3 * len(pos) // 4
    if i_end - i_start < 5:
        i_start, i_end = 0, len(pos)

    log_lam = np.log(pos[i_start:i_end])
    log_rank = np.log(ranks[i_start:i_end].astype(float))
    coeffs = np.polyfit(log_lam, log_rank, 1)
    d_spec = 2.0 * coeffs[0]
    predicted = np.polyval(coeffs, log_lam)
    ss_res = np.sum((log_rank - predicted) ** 2)
    ss_tot = np.sum((log_rank - np.mean(log_rank)) ** 2)
    R2 = 1 - ss_res / ss_tot if ss_tot > 0 else 0

    print(f"    Spectral dim {label}: {d_spec:.3f} (R²={R2:.4f})")
    return {'dim': d_spec, 'R2': R2, 'eigenvalues': pos}


# ============================================================================
# Route 3: Fisher Information — per-vertex rank and PR
# ============================================================================

def compute_fisher_per_vertex(adj, v0_idx, n, sigma):
    dist_v0 = bfs_distances_fast(adj, v0_idx, n)
    log_unnorm_v0 = -dist_v0.astype(float) / sigma
    log_unnorm_v0[dist_v0 < 0] = -1000.0
    log_Z_v0 = np.logaddexp.reduce(log_unnorm_v0)
    log_p_v0 = log_unnorm_v0 - log_Z_v0
    p_v0 = np.exp(log_p_v0)

    neighbors = adj[v0_idx]
    k = len(neighbors)
    if k < 2:
        return None

    score_vectors = np.zeros((k, n))
    for j, w_idx in enumerate(neighbors):
        dist_w = bfs_distances_fast(adj, w_idx, n)
        log_unnorm_w = -dist_w.astype(float) / sigma
        log_unnorm_w[dist_w < 0] = -1000.0
        log_Z_w = np.logaddexp.reduce(log_unnorm_w)
        log_p_w = log_unnorm_w - log_Z_w
        score_vectors[j, :] = log_p_w - log_p_v0

    weighted_scores = score_vectors * np.sqrt(p_v0)[np.newaxis, :]
    F = weighted_scores @ weighted_scores.T

    frob_norm = np.linalg.norm(F, 'fro')
    if frob_norm < 1e-12:
        return None

    sv = np.linalg.svd(F, compute_uv=False)
    sv_norm = sv / sv[0] if sv[0] > 0 else sv

    if len(sv_norm) > 1:
        ratios = sv_norm[1:] / sv_norm[:-1]
        ratios[ratios == 0] = 1e-15
        gap_rank = int(np.argmin(ratios)) + 1
    else:
        gap_rank = 1

    if np.sum(sv ** 2) > 0:
        pr = (np.sum(sv)) ** 2 / np.sum(sv ** 2)
    else:
        pr = 0.0

    return {'sv': sv, 'sv_norm': sv_norm, 'rank': gap_rank, 'pr': pr, 'degree': k}


def fisher_analysis(adj, n, sigma, sample_indices, label=""):
    results = []
    n_degenerate = 0
    for v_idx in sample_indices:
        r = compute_fisher_per_vertex(adj, v_idx, n, sigma)
        if r is None:
            n_degenerate += 1
        else:
            r['vertex_idx'] = v_idx
            results.append(r)

    if not results:
        print(f"    Fisher {label}: ALL DEGENERATE")
        return None

    ranks = [r['rank'] for r in results]
    prs = [r['pr'] for r in results]
    degrees = [r['degree'] for r in results]

    print(f"    Fisher {label} (σ={sigma}): rank mean={np.mean(ranks):.2f} "
          f"median={np.median(ranks):.0f} std={np.std(ranks):.2f}, "
          f"PR mean={np.mean(prs):.2f} std={np.std(prs):.2f}, "
          f"deg mean={np.mean(degrees):.1f}")

    return {
        'results': results, 'ranks': ranks, 'prs': prs,
        'degrees': degrees, 'n_degenerate': n_degenerate, 'sigma': sigma,
    }


def fisher_sigma_sweep(adj, n, sigmas, sample_indices, label=""):
    sweep = {}
    for sigma in sigmas:
        fa = fisher_analysis(adj, n, sigma, sample_indices, label=f"{label} σ={sigma}")
        if fa is not None:
            sweep[sigma] = fa
    return sweep


# ============================================================================
# Plotting helpers
# ============================================================================

def get_modal_degree(fisher_data):
    degrees = [r['degree'] for r in fisher_data['results']]
    return Counter(degrees).most_common(1)[0][0]


def get_sv_profiles_at_degree(fisher_data, deg):
    """Get SV profiles at given degree, or closest available."""
    profiles = [r['sv_norm'] for r in fisher_data['results'] if r['degree'] == deg]
    if not profiles:
        all_deg = sorted(set(r['degree'] for r in fisher_data['results']))
        closest = min(all_deg, key=lambda d: abs(d - deg))
        profiles = [r['sv_norm'] for r in fisher_data['results'] if r['degree'] == closest]
        deg = closest
    return profiles, deg


def plot_degree_comparison(deg_bounded, deg_periodic, d, filename):
    """Overlay degree histograms: bounded vs periodic."""
    fig, ax = plt.subplots(figsize=(9, 5))
    all_degs = deg_bounded + deg_periodic
    bins = range(min(all_degs), max(all_degs) + 2)
    ax.hist(deg_bounded, bins=bins, align='left', alpha=0.5, edgecolor='black',
            label=f'Bounded (mean={np.mean(deg_bounded):.1f}, std={np.std(deg_bounded):.1f})', color='coral')
    ax.hist(deg_periodic, bins=bins, align='left', alpha=0.5, edgecolor='black',
            label=f'Periodic (mean={np.mean(deg_periodic):.1f}, std={np.std(deg_periodic):.1f})', color='steelblue')
    ax.set_xlabel('Degree')
    ax.set_ylabel('Count')
    ax.set_title(f'Degree Distribution: Bounded vs Periodic RGG d={d}')
    ax.legend()
    plt.tight_layout()
    plt.savefig(os.path.join(OUTPUT_DIR, filename), dpi=150)
    plt.close()


def plot_fisher_rank_hist(ranks, true_d, name, filename):
    fig, ax = plt.subplots(figsize=(8, 5))
    bins = range(min(ranks), max(ranks) + 2)
    ax.hist(ranks, bins=bins, align='left', edgecolor='black', alpha=0.7)
    ax.axvline(true_d, color='red', linestyle='--', linewidth=2, label=f'True d={true_d}')
    ax.set_xlabel('Fisher Gap-Based Rank')
    ax.set_ylabel('Count')
    ax.set_title(f'Fisher Rank Distribution: {name}')
    ax.legend()
    plt.tight_layout()
    plt.savefig(os.path.join(OUTPUT_DIR, filename), dpi=150)
    plt.close()


def plot_fisher_pr_hist(prs, true_d, name, filename):
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.hist(prs, bins=20, edgecolor='black', alpha=0.7)
    ax.axvline(true_d, color='red', linestyle='--', linewidth=2, label=f'True d={true_d}')
    ax.set_xlabel('Fisher Participation Ratio')
    ax.set_ylabel('Count')
    ax.set_title(f'Fisher PR Distribution: {name}')
    ax.legend()
    plt.tight_layout()
    plt.savefig(os.path.join(OUTPUT_DIR, filename), dpi=150)
    plt.close()


def plot_sv_profile(fisher_data, name, filename, target_deg=None):
    """SV profile bar chart at modal (or target) degree."""
    if target_deg is None:
        target_deg = get_modal_degree(fisher_data)
    profiles, actual_deg = get_sv_profiles_at_degree(fisher_data, target_deg)
    if not profiles:
        return actual_deg

    max_len = max(len(sv) for sv in profiles)
    padded = np.zeros((len(profiles), max_len))
    for i, sv in enumerate(profiles):
        padded[i, :len(sv)] = sv
    mean_sv = np.mean(padded, axis=0)
    std_sv = np.std(padded, axis=0)

    fig, ax = plt.subplots(figsize=(10, 5))
    x = np.arange(1, max_len + 1)
    ax.bar(x, mean_sv, yerr=std_sv, capsize=3, edgecolor='black', alpha=0.7)
    ax.set_xlabel('Singular Value Index')
    ax.set_ylabel('Normalized SV')
    ax.set_title(f'Mean SV Profile (degree={actual_deg}, n={len(profiles)}): {name}')
    plt.tight_layout()
    plt.savefig(os.path.join(OUTPUT_DIR, filename), dpi=150)
    plt.close()
    return actual_deg


def plot_pr_vs_degree(fisher_data, true_d, name, filename):
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.scatter(fisher_data['degrees'], fisher_data['prs'], alpha=0.5, s=20)
    ax.axhline(true_d, color='red', linestyle='--', linewidth=2, label=f'True d={true_d}')
    corr = np.corrcoef(fisher_data['degrees'], fisher_data['prs'])[0, 1]
    ax.set_xlabel('Vertex Degree')
    ax.set_ylabel('Fisher PR')
    ax.set_title(f'PR vs Degree (r={corr:.3f}): {name}')
    ax.legend()
    plt.tight_layout()
    plt.savefig(os.path.join(OUTPUT_DIR, filename), dpi=150)
    plt.close()
    return corr


def plot_growth(growth, true_d, name, filename):
    fig, ax = plt.subplots(figsize=(8, 5))
    rs = np.array(growth['rs'], dtype=float)
    vols = np.array(growth['vols'], dtype=float)
    ax.plot(np.log(rs), np.log(vols), 'b.-', label='Data')
    if 'coeffs' in growth:
        r_fit = np.linspace(np.log(growth.get('rMin', 2)), np.log(growth.get('rMax', max(rs))), 50)
        ax.plot(r_fit, np.polyval(growth['coeffs'], r_fit), 'r--', linewidth=2,
                label=f"Fit: slope={growth['dim']:.3f} (R²={growth['R2']:.4f})")
    ax.set_xlabel('log(r)')
    ax.set_ylabel('log(V(r))')
    ax.set_title(f'Growth Dimension: {name} (true d={true_d})')
    ax.legend()
    plt.tight_layout()
    plt.savefig(os.path.join(OUTPUT_DIR, filename), dpi=150)
    plt.close()


# ============================================================================
# Three-panel comparison: bounded vs periodic vs ER
# ============================================================================

def plot_three_panel_ranks(bounded_ranks, periodic_ranks, er_ranks, true_d, filename):
    """Three-panel rank comparison."""
    fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(16, 5))

    all_ranks = bounded_ranks + periodic_ranks + (er_ranks if er_ranks else [])
    bins = range(min(all_ranks), max(all_ranks) + 2)

    ax1.hist(bounded_ranks, bins=bins, align='left', edgecolor='black', alpha=0.7, color='coral')
    ax1.axvline(true_d, color='red', linestyle='--', linewidth=2)
    ax1.set_title(f'Bounded RGG d={true_d}\nrank mode={Counter(bounded_ranks).most_common(1)[0][0]}')
    ax1.set_xlabel('Fisher Rank')
    ax1.set_ylabel('Count')

    ax2.hist(periodic_ranks, bins=bins, align='left', edgecolor='black', alpha=0.7, color='steelblue')
    ax2.axvline(true_d, color='red', linestyle='--', linewidth=2)
    ax2.set_title(f'Periodic RGG d={true_d}\nrank mode={Counter(periodic_ranks).most_common(1)[0][0]}')
    ax2.set_xlabel('Fisher Rank')

    if er_ranks:
        ax3.hist(er_ranks, bins=bins, align='left', edgecolor='black', alpha=0.7, color='gray')
        ax3.axvline(true_d, color='red', linestyle='--', linewidth=2)
        ax3.set_title(f'ER (no geometry)\nrank mode={Counter(er_ranks).most_common(1)[0][0]}')
    else:
        ax3.text(0.5, 0.5, 'ER data\nnot available', ha='center', va='center', transform=ax3.transAxes)
        ax3.set_title('ER (no geometry)')
    ax3.set_xlabel('Fisher Rank')

    plt.suptitle(f'Fisher Rank: Bounded vs Periodic vs ER (d={true_d})', fontsize=14, y=1.02)
    plt.tight_layout()
    plt.savefig(os.path.join(OUTPUT_DIR, filename), dpi=150, bbox_inches='tight')
    plt.close()


# ============================================================================
# Boundary analysis plots
# ============================================================================

def plot_boundary_analysis_bars(interior_data, boundary_data, d, filename):
    """Grouped bar chart: interior vs boundary Fisher rank and PR."""
    fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(15, 5))

    # Rank comparison
    categories = ['Interior', 'Boundary']
    rank_means = [np.mean(interior_data['ranks']), np.mean(boundary_data['ranks'])]
    rank_stds = [np.std(interior_data['ranks']), np.std(boundary_data['ranks'])]
    bars = ax1.bar(categories, rank_means, yerr=rank_stds, capsize=5, color=['steelblue', 'coral'],
                   edgecolor='black', alpha=0.8)
    ax1.axhline(d, color='red', linestyle='--', linewidth=2, label=f'd={d}')
    ax1.set_ylabel('Fisher Rank (mean)')
    ax1.set_title(f'Fisher Rank: Interior vs Boundary (d={d})')
    ax1.legend()
    for bar, val in zip(bars, rank_means):
        ax1.text(bar.get_x() + bar.get_width()/2., bar.get_height() + 0.05,
                f'{val:.2f}', ha='center', va='bottom', fontweight='bold')

    # PR comparison
    pr_means = [np.mean(interior_data['prs']), np.mean(boundary_data['prs'])]
    pr_stds = [np.std(interior_data['prs']), np.std(boundary_data['prs'])]
    bars = ax2.bar(categories, pr_means, yerr=pr_stds, capsize=5, color=['steelblue', 'coral'],
                   edgecolor='black', alpha=0.8)
    ax2.axhline(d, color='red', linestyle='--', linewidth=2, label=f'd={d}')
    ax2.set_ylabel('Fisher PR (mean)')
    ax2.set_title(f'Fisher PR: Interior vs Boundary (d={d})')
    ax2.legend()
    for bar, val in zip(bars, pr_means):
        ax2.text(bar.get_x() + bar.get_width()/2., bar.get_height() + 0.05,
                f'{val:.2f}', ha='center', va='bottom', fontweight='bold')

    # Degree comparison
    deg_means = [np.mean(interior_data['degrees']), np.mean(boundary_data['degrees'])]
    deg_stds = [np.std(interior_data['degrees']), np.std(boundary_data['degrees'])]
    bars = ax3.bar(categories, deg_means, yerr=deg_stds, capsize=5, color=['steelblue', 'coral'],
                   edgecolor='black', alpha=0.8)
    ax3.set_ylabel('Degree (mean)')
    ax3.set_title(f'Degree: Interior vs Boundary (d={d})')
    for bar, val in zip(bars, deg_means):
        ax3.text(bar.get_x() + bar.get_width()/2., bar.get_height() + 0.2,
                f'{val:.1f}', ha='center', va='bottom', fontweight='bold')

    plt.tight_layout()
    plt.savefig(os.path.join(OUTPUT_DIR, filename), dpi=150)
    plt.close()


def plot_spatial_scatter_by_rank(points, fisher_results, d_true, filename):
    """2D scatter plot of vertex positions colored by Fisher rank."""
    fig, ax = plt.subplots(figsize=(8, 8))

    positions = []
    ranks = []
    for r in fisher_results:
        v_idx = r['vertex_idx']
        positions.append(points[v_idx])
        ranks.append(r['rank'])

    positions = np.array(positions)
    ranks = np.array(ranks)

    unique_ranks = sorted(set(ranks))
    colors_map = {d_true: 'blue', d_true + 1: 'red'}
    for rank_val in unique_ranks:
        mask = ranks == rank_val
        color = colors_map.get(rank_val, 'green')
        ax.scatter(positions[mask, 0], positions[mask, 1], c=color, s=40,
                  label=f'Rank {rank_val} (n={np.sum(mask)})', alpha=0.7, edgecolors='black', linewidths=0.5)

    # Draw boundary margin
    margin = 0.15
    rect = plt.Rectangle((margin, margin), 1 - 2*margin, 1 - 2*margin,
                         fill=False, edgecolor='gray', linestyle='--', linewidth=2,
                         label=f'Interior (margin={margin})')
    ax.add_patch(rect)

    ax.set_xlim(-0.02, 1.02)
    ax.set_ylim(-0.02, 1.02)
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_title(f'Spatial Distribution of Fisher Rank (d={d_true})')
    ax.set_aspect('equal')
    ax.legend(loc='upper right')
    plt.tight_layout()
    plt.savefig(os.path.join(OUTPUT_DIR, filename), dpi=150)
    plt.close()


# ============================================================================
# Test A: Periodic RGG
# ============================================================================

def run_periodic_rgg(d, n_points, r, sigmas=[2.0, 3.0, 5.0, 8.0]):
    label = f"Periodic RGG d={d}"
    print(f"\n{'='*70}")
    print(f"  TEST A: {label} (N={n_points}, r={r})")
    print(f"{'='*70}")

    print("  Building periodic RGG...")
    G, points = make_periodic_rgg(n_points, d, r)
    G_lcc, stats = validate_graph(G, label)

    adj, nodes, node_to_idx, n = build_adjacency_list(G_lcc)
    degrees = [len(adj[i]) for i in range(n)]

    # Fisher: sample 50 random vertices (no restriction needed — no boundary depletion)
    fisher_samples = np.random.choice(n, min(50, n), replace=False).tolist()
    print(f"  Fisher samples: {len(fisher_samples)} random vertices")

    # Route 1: Growth
    print("\n  Route 1: Growth dimension")
    growth = growth_dimension(adj, n, n_samples=30, label=label)
    plot_growth(growth, d, label, f"periodic_rgg_{d}d_growth.png")

    # Route 2: Spectral
    print("\n  Route 2: Spectral dimension")
    spectral = spectral_dimension(G_lcc, n_eigs=300, label=label)

    # Route 3: Fisher (main at σ=2.0)
    print(f"\n  Route 3: Fisher (σ={sigmas[0]}, {len(fisher_samples)} samples)")
    fisher_main = fisher_analysis(adj, n, sigma=sigmas[0], sample_indices=fisher_samples, label=label)

    if fisher_main:
        plot_fisher_rank_hist(fisher_main['ranks'], d, label, f"periodic_rgg_{d}d_fisher_rank_hist.png")
        plot_fisher_pr_hist(fisher_main['prs'], d, label, f"periodic_rgg_{d}d_fisher_pr_hist.png")
        modal_deg = plot_sv_profile(fisher_main, label, f"periodic_rgg_{d}d_sv_profile.png")
        pr_deg_corr = plot_pr_vs_degree(fisher_main, d, label, f"periodic_rgg_{d}d_pr_vs_degree.png")
    else:
        modal_deg = None
        pr_deg_corr = None

    # Sigma sweep
    print(f"\n  Fisher sigma sweep: {sigmas}")
    sweep = fisher_sigma_sweep(adj, n, sigmas, fisher_samples, label=label)

    return {
        'stats': stats, 'growth': growth, 'spectral': spectral,
        'fisher_main': fisher_main, 'fisher_sweep': sweep,
        'modal_degree': modal_deg, 'pr_deg_corr': pr_deg_corr,
        'degrees': degrees, 'points': points,
    }


# ============================================================================
# Test B: Boundary analysis on bounded RGG
# ============================================================================

def classify_vertex(position, margin=0.15):
    return 'interior' if all(margin <= c <= 1 - margin for c in position) else 'boundary'


def run_boundary_analysis(d, n_points, r, sigma=2.0, n_per_group=50, margin=0.15):
    label = f"Bounded RGG d={d} boundary analysis"
    print(f"\n{'='*70}")
    print(f"  TEST B: {label}")
    print(f"{'='*70}")

    print("  Building bounded RGG...")
    G, points = make_bounded_rgg(n_points, d, r)

    # Take LCC
    if not nx.is_connected(G):
        lcc = max(nx.connected_components(G), key=len)
        G = G.subgraph(lcc).copy()
        lcc_nodes = sorted(lcc)
        # Remap points
        points = points[lcc_nodes]
        mapping = {v: i for i, v in enumerate(lcc_nodes)}
        G = nx.relabel_nodes(G, mapping)

    adj, nodes, node_to_idx, n = build_adjacency_list(G)

    # Classify all vertices
    categories = {}
    for i in range(n):
        cat = classify_vertex(points[i], margin=margin)
        if cat not in categories:
            categories[cat] = []
        categories[cat].append(i)

    n_interior = len(categories.get('interior', []))
    n_boundary = len(categories.get('boundary', []))
    print(f"  Interior: {n_interior} ({100*n_interior/n:.1f}%), Boundary: {n_boundary} ({100*n_boundary/n:.1f}%)")

    # Sample from each group
    interior_indices = categories.get('interior', [])
    boundary_indices = categories.get('boundary', [])
    np.random.shuffle(interior_indices)
    np.random.shuffle(boundary_indices)
    interior_sample = interior_indices[:n_per_group]
    boundary_sample = boundary_indices[:n_per_group]

    print(f"  Sampling {len(interior_sample)} interior + {len(boundary_sample)} boundary vertices")

    # Run Fisher on each group
    print(f"\n  Fisher on interior vertices (σ={sigma}):")
    fisher_interior = fisher_analysis(adj, n, sigma=sigma, sample_indices=interior_sample,
                                      label=f"Interior d={d}")

    print(f"\n  Fisher on boundary vertices (σ={sigma}):")
    fisher_boundary = fisher_analysis(adj, n, sigma=sigma, sample_indices=boundary_sample,
                                      label=f"Boundary d={d}")

    # Also get full-sample results
    all_sample = interior_sample + boundary_sample
    fisher_all = fisher_analysis(adj, n, sigma=sigma, sample_indices=all_sample,
                                 label=f"Full d={d}")

    # Plots
    if fisher_interior and fisher_boundary:
        plot_boundary_analysis_bars(fisher_interior, fisher_boundary, d,
                                   f"bounded_rgg_{d}d_boundary_analysis_bars.png")

    if d == 2 and fisher_all:
        plot_spatial_scatter_by_rank(points, fisher_all['results'], d,
                                    f"bounded_rgg_{d}d_boundary_analysis_spatial.png")

    return {
        'fisher_interior': fisher_interior,
        'fisher_boundary': fisher_boundary,
        'fisher_all': fisher_all,
        'n_interior': n_interior,
        'n_boundary': n_boundary,
        'n_total': n,
        'points': points,
        'margin': margin,
    }


# ============================================================================
# Results report
# ============================================================================

def write_results(periodic_2d, periodic_3d, boundary_2d, boundary_3d,
                  bounded_2d_degrees=None, bounded_3d_degrees=None):
    lines = []
    lines.append("# Phase 1: Periodic RGG + Boundary Analysis — Results")
    lines.append(f"\nGenerated: {datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    lines.append("")

    # Test A: Periodic RGG
    lines.append("## Test A: Periodic Random Geometric Graphs")
    lines.append("")

    for tag, pdata, true_d in [("d=2", periodic_2d, 2), ("d=3", periodic_3d, 3)]:
        lines.append(f"### Periodic RGG {tag}")
        lines.append("")
        s = pdata['stats']
        lines.append(f"- **Vertices (LCC):** {s['n_lcc']} ({s['lcc_frac']*100:.1f}%)")
        lines.append(f"- **Edges:** {s['m_lcc']}")
        lines.append(f"- **Degree:** min={s['degree_min']}, max={s['degree_max']}, "
                     f"mean={s['degree_mean']:.1f}, median={s['degree_median']:.0f}, std={s['degree_std']:.1f}")
        lines.append(f"- **Clustering:** {s['clustering']:.4f}")
        lines.append(f"- **Diameter:** {s['diameter']} (est.)")
        lines.append("")

        g = pdata['growth']
        lines.append(f"**Growth Dimension:** {g['dim']:.3f} (R²={g['R2']:.4f})")
        gate = [true_d - 0.3, true_d + 0.3]
        gpass = gate[0] <= g['dim'] <= gate[1]
        lines.append(f"  Gate [{gate[0]:.1f}, {gate[1]:.1f}]: **{'PASS' if gpass else 'FAIL'}**")
        lines.append("")

        sp = pdata['spectral']
        lines.append(f"**Spectral Dimension:** {sp['dim']:.3f} (R²={sp['R2']:.4f})")
        lines.append("")

        if pdata['fisher_main']:
            fm = pdata['fisher_main']
            lines.append(f"**Fisher Rank Distribution:**")
            lines.append(f"  Mean={np.mean(fm['ranks']):.2f}, Median={np.median(fm['ranks']):.0f}, "
                        f"Std={np.std(fm['ranks']):.2f}")
            rank_counter = Counter(fm['ranks'])
            lines.append(f"  Distribution: {dict(sorted(rank_counter.items()))}")
            lines.append("")

            lines.append(f"**Fisher PR Distribution:**")
            lines.append(f"  Mean={np.mean(fm['prs']):.3f}, Median={np.median(fm['prs']):.3f}, "
                        f"Std={np.std(fm['prs']):.3f}")
            lines.append("")

            if pdata['pr_deg_corr'] is not None:
                lines.append(f"**PR-Degree Correlation:** r = {pdata['pr_deg_corr']:.3f}")
                lines.append("")

            if pdata['fisher_sweep']:
                lines.append("**Sigma Sweep:**")
                lines.append("")
                lines.append("| σ | PR mean | PR std | Rank mean | Rank mode |")
                lines.append("|---|---------|--------|-----------|-----------|")
                for sigma, data in sorted(pdata['fisher_sweep'].items()):
                    mode = Counter(data['ranks']).most_common(1)[0][0]
                    lines.append(f"| {sigma:.1f} | {np.mean(data['prs']):.3f} | {np.std(data['prs']):.3f} "
                               f"| {np.mean(data['ranks']):.2f} | {mode} |")
                lines.append("")

        lines.append("---")
        lines.append("")

    # Comparison with bounded RGG (from previous run)
    lines.append("### Comparison: Periodic vs Bounded RGG")
    lines.append("")
    lines.append("| Metric | Bounded d=2 | Periodic d=2 | Bounded d=3 | Periodic d=3 |")
    lines.append("|--------|-------------|--------------|-------------|--------------|")

    if periodic_2d['fisher_main'] and periodic_3d['fisher_main']:
        # Bounded data from previous run (hardcoded from results)
        lines.append(f"| Growth dim | 1.573 | {periodic_2d['growth']['dim']:.3f} | 2.309 | {periodic_3d['growth']['dim']:.3f} |")
        lines.append(f"| Fisher rank (mean) | 2.90 | {np.mean(periodic_2d['fisher_main']['ranks']):.2f} "
                    f"| 3.60 | {np.mean(periodic_3d['fisher_main']['ranks']):.2f} |")
        lines.append(f"| Fisher rank (mode) | 3 | {Counter(periodic_2d['fisher_main']['ranks']).most_common(1)[0][0]} "
                    f"| 4 | {Counter(periodic_3d['fisher_main']['ranks']).most_common(1)[0][0]} |")
        lines.append(f"| Fisher PR (σ=2) | 3.210 | {np.mean(periodic_2d['fisher_main']['prs']):.3f} "
                    f"| 3.982 | {np.mean(periodic_3d['fisher_main']['prs']):.3f} |")
        lines.append(f"| PR-deg corr | 0.325 | {periodic_2d['pr_deg_corr']:.3f} "
                    f"| 0.464 | {periodic_3d['pr_deg_corr']:.3f} |")
        lines.append(f"| Degree std | 3.9 | {periodic_2d['stats']['degree_std']:.1f} "
                    f"| 4.4 | {periodic_3d['stats']['degree_std']:.1f} |")
    lines.append("")

    # Test B: Boundary analysis
    lines.append("## Test B: Interior vs Boundary Analysis (Bounded RGG)")
    lines.append("")

    for tag, bdata, true_d in [("d=2", boundary_2d, 2), ("d=3", boundary_3d, 3)]:
        lines.append(f"### Bounded RGG {tag}")
        lines.append("")
        lines.append(f"- Interior: {bdata['n_interior']} vertices ({100*bdata['n_interior']/bdata['n_total']:.1f}%), "
                    f"margin={bdata['margin']}")
        lines.append(f"- Boundary: {bdata['n_boundary']} vertices ({100*bdata['n_boundary']/bdata['n_total']:.1f}%)")
        lines.append("")

        fi = bdata['fisher_interior']
        fb = bdata['fisher_boundary']
        fa = bdata['fisher_all']

        if fi and fb:
            lines.append(f"| Metric | Interior | Boundary | Full Sample |")
            lines.append(f"|--------|----------|----------|-------------|")
            lines.append(f"| N vertices | {len(fi['ranks'])} | {len(fb['ranks'])} | {len(fa['ranks'])} |")
            lines.append(f"| Mean degree | {np.mean(fi['degrees']):.1f} | {np.mean(fb['degrees']):.1f} | {np.mean(fa['degrees']):.1f} |")
            lines.append(f"| Fisher rank (mean) | {np.mean(fi['ranks']):.2f} | {np.mean(fb['ranks']):.2f} | {np.mean(fa['ranks']):.2f} |")
            fi_mode = Counter(fi['ranks']).most_common(1)[0][0]
            fb_mode = Counter(fb['ranks']).most_common(1)[0][0]
            lines.append(f"| Fisher rank (mode) | {fi_mode} | {fb_mode} | {Counter(fa['ranks']).most_common(1)[0][0]} |")
            lines.append(f"| Fisher PR (mean) | {np.mean(fi['prs']):.3f} | {np.mean(fb['prs']):.3f} | {np.mean(fa['prs']):.3f} |")
            lines.append(f"| Fisher PR (std) | {np.std(fi['prs']):.3f} | {np.std(fb['prs']):.3f} | {np.std(fa['prs']):.3f} |")
            lines.append("")

        lines.append("---")
        lines.append("")

    path = os.path.join(OUTPUT_DIR, "PHASE1_PERIODIC_RGG_RESULTS.md")
    with open(path, 'w') as f:
        f.write('\n'.join(lines))
    print(f"\n  Results written to {path}")
    return path


# ============================================================================
# MAIN
# ============================================================================

if __name__ == "__main__":
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    t0 = time.time()

    print("=" * 70)
    print("DS Phase 1: Periodic RGG + Boundary Analysis")
    print("=" * 70)

    # ---- Test A1: Periodic RGG d=2 ----
    periodic_2d = run_periodic_rgg(d=2, n_points=5000, r=0.03, sigmas=[2.0, 3.0, 5.0, 8.0])

    # ---- Test A2: Periodic RGG d=3 ----
    periodic_3d = run_periodic_rgg(d=3, n_points=5000, r=0.09, sigmas=[2.0, 3.0, 5.0, 8.0])

    # ---- Test B1: Boundary analysis d=2 ----
    boundary_2d = run_boundary_analysis(d=2, n_points=5000, r=0.03, sigma=2.0,
                                        n_per_group=50, margin=0.15)

    # ---- Test B2: Boundary analysis d=3 ----
    boundary_3d = run_boundary_analysis(d=3, n_points=5000, r=0.09, sigma=2.0,
                                        n_per_group=50, margin=0.15)

    # ---- Degree comparison plots ----
    print(f"\n{'='*70}")
    print("  COMPARISON PLOTS")
    print(f"{'='*70}")

    # We need bounded RGG degree distributions for comparison
    # Re-build bounded RGGs to get their degree distributions
    print("  Building bounded RGGs for degree comparison...")
    np.random.seed(42)
    G_bounded_2d, _ = make_bounded_rgg(5000, 2, 0.03)
    bounded_2d_degrees = [d for _, d in G_bounded_2d.degree()]
    G_bounded_3d, _ = make_bounded_rgg(5000, 3, 0.09)
    bounded_3d_degrees = [d for _, d in G_bounded_3d.degree()]

    plot_degree_comparison(bounded_2d_degrees, periodic_2d['degrees'], 2,
                          "periodic_rgg_2d_degree_comparison.png")
    print("  Saved periodic_rgg_2d_degree_comparison.png")

    plot_degree_comparison(bounded_3d_degrees, periodic_3d['degrees'], 3,
                          "periodic_rgg_3d_degree_comparison.png")
    print("  Saved periodic_rgg_3d_degree_comparison.png")

    # Three-panel rank comparison (bounded vs periodic vs ER)
    # Get bounded RGG ranks from previous run results (hardcoded)
    # Bounded d=2: rank distribution {1: 1, 2: 4, 3: 44, 4: 1} => 50 samples
    bounded_2d_ranks = [1] + [2]*4 + [3]*44 + [4]*1
    bounded_3d_ranks = [1]*1 + [2]*2 + [3]*13 + [4]*34
    er_ranks = [1]*30  # ER n=1000 all had rank 1

    if periodic_2d['fisher_main']:
        plot_three_panel_ranks(bounded_2d_ranks, periodic_2d['fisher_main']['ranks'],
                              er_ranks, 2, "rgg_periodic_vs_bounded_vs_er_ranks_2d.png")
        print("  Saved rgg_periodic_vs_bounded_vs_er_ranks_2d.png")

    if periodic_3d['fisher_main']:
        plot_three_panel_ranks(bounded_3d_ranks, periodic_3d['fisher_main']['ranks'],
                              er_ranks, 3, "rgg_periodic_vs_bounded_vs_er_ranks_3d.png")
        print("  Saved rgg_periodic_vs_bounded_vs_er_ranks_3d.png")

    # ---- Results report ----
    write_results(periodic_2d, periodic_3d, boundary_2d, boundary_3d)

    elapsed = time.time() - t0
    print(f"\n{'='*70}")
    print(f"  TOTAL TIME: {elapsed:.1f}s")
    print(f"{'='*70}")
