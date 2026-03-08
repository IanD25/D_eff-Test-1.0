#!/usr/bin/env python3
"""
DS Phase 1: Block-Spin Coarse-Graining Test (P7 / DPI Validation)
==================================================================
Tests whether D_eff is monotonically non-increasing under coarse-graining.

Test 1: 2D torus 128×128 with 2×2 block-spin (6 levels)
Test 2: 2D torus 81×81 with 3×3 block-spin (3 levels)
Test 3: 3D torus 32^3 with 2×2×2 block-spin (4 levels)
"""

import os
import time
import datetime
from collections import deque
from itertools import product

import numpy as np
import networkx as nx
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

OUTPUT_DIR = "phase1_results"
np.random.seed(42)

EPSILON = 0.05  # P7 monotonicity tolerance

# ============================================================================
# Graph Construction
# ============================================================================

def make_torus_2d(n):
    G = nx.Graph()
    for i in range(n):
        for j in range(n):
            v = i * n + j
            G.add_node(v)
            # Right neighbor
            G.add_edge(v, i * n + (j + 1) % n)
            # Down neighbor
            G.add_edge(v, ((i + 1) % n) * n + j)
    return G


def make_torus_3d(n):
    G = nx.Graph()
    for i in range(n):
        for j in range(n):
            for k in range(n):
                v = i * n * n + j * n + k
                G.add_node(v)
                G.add_edge(v, i * n * n + j * n + (k + 1) % n)
                G.add_edge(v, i * n * n + ((j + 1) % n) * n + k)
                G.add_edge(v, ((i + 1) % n) * n * n + j * n + k)
    return G


# ============================================================================
# Block-Spin Coarse-Graining
# ============================================================================

def block_spin_2d(G_fine, n_fine):
    """Coarse-grain a 2D torus by 2×2 block-spin."""
    assert n_fine % 2 == 0
    n_coarse = n_fine // 2

    fine_to_coarse = {}
    for I in range(n_coarse):
        for J in range(n_coarse):
            coarse_v = I * n_coarse + J
            for di in range(2):
                for dj in range(2):
                    fine_i = (2 * I + di) % n_fine
                    fine_j = (2 * J + dj) % n_fine
                    fine_v = fine_i * n_fine + fine_j
                    fine_to_coarse[fine_v] = coarse_v

    G_coarse = nx.Graph()
    G_coarse.add_nodes_from(range(n_coarse * n_coarse))
    for u, v in G_fine.edges():
        cu, cv = fine_to_coarse[u], fine_to_coarse[v]
        if cu != cv:
            G_coarse.add_edge(cu, cv)

    return G_coarse, n_coarse


def block_spin_2d_3x3(G_fine, n_fine):
    """Coarse-grain a 2D torus by 3×3 block-spin."""
    assert n_fine % 3 == 0
    n_coarse = n_fine // 3

    fine_to_coarse = {}
    for I in range(n_coarse):
        for J in range(n_coarse):
            coarse_v = I * n_coarse + J
            for di in range(3):
                for dj in range(3):
                    fine_i = (3 * I + di) % n_fine
                    fine_j = (3 * J + dj) % n_fine
                    fine_v = fine_i * n_fine + fine_j
                    fine_to_coarse[fine_v] = coarse_v

    G_coarse = nx.Graph()
    G_coarse.add_nodes_from(range(n_coarse * n_coarse))
    for u, v in G_fine.edges():
        cu, cv = fine_to_coarse[u], fine_to_coarse[v]
        if cu != cv:
            G_coarse.add_edge(cu, cv)

    return G_coarse, n_coarse


def block_spin_3d(G_fine, n_fine):
    """Coarse-grain a 3D torus by 2×2×2 block-spin."""
    assert n_fine % 2 == 0
    n_coarse = n_fine // 2

    fine_to_coarse = {}
    for I in range(n_coarse):
        for J in range(n_coarse):
            for K in range(n_coarse):
                coarse_v = I * n_coarse * n_coarse + J * n_coarse + K
                for di in range(2):
                    for dj in range(2):
                        for dk in range(2):
                            fi = (2 * I + di) % n_fine
                            fj = (2 * J + dj) % n_fine
                            fk = (2 * K + dk) % n_fine
                            fine_v = fi * n_fine * n_fine + fj * n_fine + fk
                            fine_to_coarse[fine_v] = coarse_v

    G_coarse = nx.Graph()
    G_coarse.add_nodes_from(range(n_coarse ** 3))
    for u, v in G_fine.edges():
        cu, cv = fine_to_coarse[u], fine_to_coarse[v]
        if cu != cv:
            G_coarse.add_edge(cu, cv)

    return G_coarse, n_coarse


# ============================================================================
# BFS and Growth Dimension
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


def compute_growth_dimension(G, n_samples=20):
    nodes = list(G.nodes())
    n = len(nodes)
    if n < 10:
        return {'fit_dimension': float('nan'), 'r_squared': 0, 'reliable': False}

    samples = [nodes[i] for i in
               np.random.choice(n, size=min(n_samples, n), replace=False)]

    d0 = bfs_distances(G, samples[0])
    maxR = max(d0.values()) // 2
    if maxR < 3:
        return {'fit_dimension': float('nan'), 'r_squared': 0, 'reliable': False}

    all_volumes = []
    for v in samples:
        dists = bfs_distances(G, v)
        dist_vals = np.array(list(dists.values()))
        volumes = np.array([np.sum(dist_vals <= r) for r in range(maxR + 1)])
        all_volumes.append(volumes)

    mean_vol = np.mean(all_volumes, axis=0)

    rMin = max(2, int(maxR * 0.20))
    rMax = max(rMin + 2, int(maxR * 0.60))
    rMax = min(rMax, maxR)

    if rMax - rMin < 3:
        return {'fit_dimension': float('nan'), 'r_squared': 0, 'reliable': False,
                'note': f'fit range too small: [{rMin},{rMax}]'}

    radii = np.arange(rMin, rMax + 1, dtype=float)
    log_r = np.log(radii)
    log_v = np.log(np.maximum(mean_vol[rMin:rMax + 1], 1).astype(float))
    coeffs = np.polyfit(log_r, log_v, 1)
    slope = coeffs[0]
    predicted = np.polyval(coeffs, log_r)
    ss_res = np.sum((log_v - predicted) ** 2)
    ss_tot = np.sum((log_v - np.mean(log_v)) ** 2)
    r_squared = 1.0 - ss_res / ss_tot if ss_tot > 0 else 0

    reliable = (rMax - rMin >= 4) and (r_squared > 0.95)
    return {'fit_dimension': slope, 'r_squared': r_squared, 'reliable': reliable,
            'rMin': rMin, 'rMax': rMax, 'maxR': maxR}


# ============================================================================
# Spectral Dimension (Analytical for Tori)
# ============================================================================

def torus_eigenvalues_exact(side_n, dim):
    """Analytical eigenvalues for a d-dimensional torus with side length n."""
    if dim == 2:
        eigs = []
        for k1 in range(side_n):
            for k2 in range(side_n):
                lam = 2 * (1 - np.cos(2 * np.pi * k1 / side_n)) + \
                      2 * (1 - np.cos(2 * np.pi * k2 / side_n))
                eigs.append(lam)
    elif dim == 3:
        eigs = []
        for k1 in range(side_n):
            for k2 in range(side_n):
                for k3 in range(side_n):
                    lam = 2 * (1 - np.cos(2 * np.pi * k1 / side_n)) + \
                          2 * (1 - np.cos(2 * np.pi * k2 / side_n)) + \
                          2 * (1 - np.cos(2 * np.pi * k3 / side_n))
                    eigs.append(lam)
    else:
        raise ValueError(f"Unsupported dim={dim}")

    eigs = np.sort(np.array(eigs))
    return eigs


def compute_spectral_dimension(side_n, dim):
    """Spectral dimension from analytical eigenvalues."""
    if side_n < 3:
        return {'spectral_dimension': float('nan'), 'reliable': False}

    eigenvalues = torus_eigenvalues_exact(side_n, dim)
    eigenvalues = eigenvalues[eigenvalues > 1e-10]
    n_eigs = len(eigenvalues)
    if n_eigs < 10:
        return {'spectral_dimension': float('nan'), 'reliable': False}

    indices = np.arange(1, n_eigs + 1, dtype=float)
    fit_start = max(5, n_eigs // 10)
    fit_end = (n_eigs * 4) // 5
    if fit_end - fit_start < 5:
        return {'spectral_dimension': float('nan'), 'reliable': False}

    log_lam = np.log(eigenvalues[fit_start:fit_end])
    log_N = np.log(indices[fit_start:fit_end])
    coeffs = np.polyfit(log_lam, log_N, 1)
    slope = coeffs[0]
    spectral_dim = 2.0 * slope
    predicted = np.polyval(coeffs, log_lam)
    ss_res = np.sum((log_N - predicted) ** 2)
    ss_tot = np.sum((log_N - np.mean(log_N)) ** 2)
    r_squared = 1.0 - ss_res / ss_tot if ss_tot > 0 else 0

    reliable = r_squared > 0.95
    return {'spectral_dimension': spectral_dim, 'weyl_slope': slope,
            'r_squared': r_squared, 'reliable': reliable,
            'eigenvalues': eigenvalues, 'fit_range': (fit_start, fit_end)}


# ============================================================================
# Fisher Information
# ============================================================================

def compute_fisher_vertex(G, v0, nodes, sigma):
    n = len(nodes)
    dists_v0 = bfs_distances(G, v0)
    dist_arr = np.array([dists_v0[u] for u in nodes], dtype=float)
    log_unnorm = -dist_arr / sigma
    log_Z = np.logaddexp.reduce(log_unnorm)
    log_p = log_unnorm - log_Z
    p = np.exp(log_p)

    neighbors = list(G.neighbors(v0))
    k = len(neighbors)
    scores = np.zeros((k, n))
    for j, w in enumerate(neighbors):
        dw = bfs_distances(G, w)
        dw_arr = np.array([dw[u] for u in nodes], dtype=float)
        lw = -dw_arr / sigma
        lw_Z = np.logaddexp.reduce(lw)
        lw_p = lw - lw_Z
        scores[j, :] = lw_p - log_p

    weighted = scores * np.sqrt(p)[np.newaxis, :]
    F = weighted @ weighted.T

    if np.linalg.norm(F, 'fro') < 1e-12:
        return None

    sv = np.linalg.svd(F, compute_uv=False)
    sv_norm = sv / sv[0] if sv[0] > 0 else sv

    # Gap-based rank
    if len(sv_norm) > 1:
        ratios = sv_norm[1:] / np.maximum(sv_norm[:-1], 1e-15)
        gap_rank = int(np.argmin(ratios) + 1)
    else:
        gap_rank = 1

    pr = (np.sum(sv)) ** 2 / np.sum(sv ** 2) if np.sum(sv ** 2) > 0 else 0

    return {'sv_norm': sv_norm, 'gap_rank': gap_rank, 'pr': pr}


def compute_fisher_stats(G, sigma, n_samples=20):
    """Compute Fisher rank and PR statistics."""
    nodes = list(G.nodes())
    n = len(nodes)
    if n < 5:
        return {'mean_rank': float('nan'), 'mean_pr': float('nan'),
                'std_pr': 0, 'reliable': False, 'mean_sv': np.array([])}

    n_s = min(n_samples, n)
    samples = [nodes[i] for i in np.random.choice(n, size=n_s, replace=False)]

    ranks, prs, svs = [], [], []
    for v in samples:
        r = compute_fisher_vertex(G, v, nodes, sigma)
        if r is not None:
            ranks.append(r['gap_rank'])
            prs.append(r['pr'])
            svs.append(r['sv_norm'])

    if len(prs) == 0:
        return {'mean_rank': float('nan'), 'mean_pr': float('nan'),
                'std_pr': 0, 'reliable': False, 'mean_sv': np.array([])}

    max_len = max(len(sv) for sv in svs)
    padded = np.zeros((len(svs), max_len))
    for i, sv in enumerate(svs):
        padded[i, :len(sv)] = sv
    mean_sv = np.mean(padded, axis=0)

    return {
        'mean_rank': np.mean(ranks),
        'mean_pr': np.mean(prs),
        'std_pr': np.std(prs),
        'mean_sv': mean_sv,
        'reliable': len(prs) >= n_s // 2,
        'n_samples': len(prs),
    }


# ============================================================================
# Full Measurement at One Level
# ============================================================================

def measure_level(G, side_n, dim, level, sigma_list=[2.0, 3.0, 5.0, 8.0]):
    """Run all measurements on a graph at one coarse-graining level."""
    n = G.number_of_nodes()
    print(f"  Level {level}: {side_n}{'²' if dim == 2 else '³'} = {n} vertices")

    # Growth dimension
    growth = compute_growth_dimension(G, n_samples=20)
    g_str = f"{growth['fit_dimension']:.4f}" if not np.isnan(growth['fit_dimension']) else "N/A"
    g_rel = "✓" if growth.get('reliable', False) else "✗"
    print(f"    Growth: {g_str} (R²={growth.get('r_squared', 0):.4f}) [{g_rel}]")

    # Spectral dimension
    spectral = compute_spectral_dimension(side_n, dim)
    s_str = f"{spectral['spectral_dimension']:.4f}" if not np.isnan(spectral['spectral_dimension']) else "N/A"
    s_rel = "✓" if spectral.get('reliable', False) else "✗"
    print(f"    Spectral: {s_str} [{s_rel}]")

    # Fisher at multiple sigma
    fisher_by_sigma = {}
    for sigma in sigma_list:
        fs = compute_fisher_stats(G, sigma, n_samples=min(20, n))
        fisher_by_sigma[sigma] = fs
        if not np.isnan(fs['mean_pr']):
            print(f"    Fisher σ={sigma}: rank={fs['mean_rank']:.2f}, "
                  f"PR={fs['mean_pr']:.4f} ± {fs['std_pr']:.4f}")
        else:
            print(f"    Fisher σ={sigma}: N/A (too small)")

    return {
        'level': level, 'side_n': side_n, 'dim': dim, 'n_nodes': n,
        'growth': growth, 'spectral': spectral,
        'fisher_by_sigma': fisher_by_sigma,
    }


# ============================================================================
# P7 Monotonicity Check
# ============================================================================

def check_monotonicity(levels_data, sigma=3.0):
    """Check P7 monotonicity across coarse-graining levels."""
    results = []
    prev_pr = None
    all_pass = True
    worst_violation = 0

    for data in levels_data:
        fs = data['fisher_by_sigma'].get(sigma, {})
        pr = fs.get('mean_pr', float('nan'))
        rank = fs.get('mean_rank', float('nan'))
        growth = data['growth']['fit_dimension']
        spectral = data['spectral']['spectral_dimension']

        if prev_pr is not None and not np.isnan(pr) and not np.isnan(prev_pr):
            delta_pr = pr - prev_pr
            if delta_pr > EPSILON:
                status = "FAIL"
                all_pass = False
            elif delta_pr > 0:
                status = "MARGINAL"
            else:
                status = "PASS"
            worst_violation = max(worst_violation, delta_pr)
        else:
            delta_pr = None
            status = "—"

        results.append({
            'level': data['level'], 'side_n': data['side_n'],
            'n_nodes': data['n_nodes'],
            'growth': growth, 'spectral': spectral,
            'rank': rank, 'pr': pr,
            'delta_pr': delta_pr, 'status': status,
            'reliable': fs.get('reliable', False),
        })
        if not np.isnan(pr):
            prev_pr = pr

    return results, all_pass, worst_violation


# ============================================================================
# Plotting
# ============================================================================

def plot_deff_vs_level(levels_data, ref_dim, title, filename, sigma=3.0):
    """D_eff vs coarse-graining level — the KEY plot."""
    fig, ax = plt.subplots(figsize=(10, 7))

    levels = [d['level'] for d in levels_data]
    growth = [d['growth']['fit_dimension'] for d in levels_data]
    spectral = [d['spectral']['spectral_dimension'] for d in levels_data]
    fisher_pr = [d['fisher_by_sigma'].get(sigma, {}).get('mean_pr', float('nan'))
                 for d in levels_data]

    # Filter NaN for plotting
    valid_g = [(l, g) for l, g in zip(levels, growth) if not np.isnan(g)]
    valid_s = [(l, s) for l, s in zip(levels, spectral) if not np.isnan(s)]
    valid_f = [(l, f) for l, f in zip(levels, fisher_pr) if not np.isnan(f)]

    if valid_g:
        ax.plot(*zip(*valid_g), 'o-', color='#2196F3', linewidth=2, markersize=8,
                label='Growth Dim')
    if valid_s:
        ax.plot(*zip(*valid_s), 's-', color='#FF5722', linewidth=2, markersize=8,
                label='Spectral Dim')
    if valid_f:
        ax.plot(*zip(*valid_f), '^-', color='#4CAF50', linewidth=2, markersize=8,
                label=f'Fisher PR (σ={sigma})')

    ax.axhline(y=ref_dim, color='black', linestyle='--', linewidth=2, alpha=0.5,
               label=f'd = {ref_dim}')

    # Annotate with lattice sizes
    for d in levels_data:
        side = d['side_n']
        dim = d['dim']
        label = f"{side}{'²' if dim == 2 else '³'}"
        ax.annotate(label, (d['level'], ref_dim - 0.15), ha='center', fontsize=9, color='gray')

    ax.set_xlabel('Coarse-Graining Level', fontsize=13)
    ax.set_ylabel('Dimension Estimate', fontsize=13)
    ax.set_title(title, fontsize=14)
    ax.legend(fontsize=10, loc='best')
    ax.grid(True, alpha=0.3)
    ax.set_xticks(levels)
    path = os.path.join(OUTPUT_DIR, filename)
    fig.savefig(path, dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f"  Saved {path}")


def plot_sv_profiles(levels_data, title, filename, sigma=3.0):
    """Panel plot of Fisher SV profiles at each coarse-graining level."""
    n = len(levels_data)
    fig, axes = plt.subplots(1, n, figsize=(3.5 * n, 4))
    if n == 1:
        axes = [axes]

    for idx, data in enumerate(levels_data):
        ax = axes[idx]
        fs = data['fisher_by_sigma'].get(sigma, {})
        sv = fs.get('mean_sv', np.array([]))
        side = data['side_n']
        dim = data['dim']

        if len(sv) > 0:
            x = np.arange(1, len(sv) + 1)
            ax.bar(x, sv, color='#2196F3', edgecolor='black', linewidth=0.5)
            ax.set_ylim(0, 1.15)
            ax.set_xticks(x)
        else:
            ax.text(0.5, 0.5, 'N/A', ha='center', va='center', transform=ax.transAxes)

        ax.set_title(f'L{data["level"]}: {side}{"²" if dim == 2 else "³"}', fontsize=10)
        ax.grid(True, alpha=0.3, axis='y')

    fig.suptitle(title, fontsize=13)
    fig.tight_layout()
    path = os.path.join(OUTPUT_DIR, filename)
    fig.savefig(path, dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f"  Saved {path}")


def plot_pr_vs_sigma_by_level(levels_data, title, filename):
    """PR-vs-sigma at each coarse-graining level."""
    fig, ax = plt.subplots(figsize=(10, 7))
    colors = plt.cm.viridis(np.linspace(0, 0.9, len(levels_data)))

    for idx, data in enumerate(levels_data):
        sigmas, prs = [], []
        for sigma in sorted(data['fisher_by_sigma'].keys()):
            fs = data['fisher_by_sigma'][sigma]
            if not np.isnan(fs.get('mean_pr', float('nan'))):
                sigmas.append(sigma)
                prs.append(fs['mean_pr'])
        if sigmas:
            side = data['side_n']
            dim = data['dim']
            label = f'L{data["level"]}: {side}{"²" if dim == 2 else "³"}'
            ax.plot(sigmas, prs, 'o-', color=colors[idx], linewidth=2,
                    markersize=6, label=label)

    ax.set_xlabel('σ', fontsize=13)
    ax.set_ylabel('Participation Ratio', fontsize=13)
    ax.set_title(title, fontsize=14)
    ax.legend(fontsize=9, loc='best')
    ax.grid(True, alpha=0.3)
    path = os.path.join(OUTPUT_DIR, filename)
    fig.savefig(path, dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f"  Saved {path}")


# ============================================================================
# Test 1: 2D Torus, 2×2 Block-Spin
# ============================================================================

def run_test_2d_2x2():
    print(f"\n{'='*70}")
    print("  TEST 1: 2D Torus, 2×2 Block-Spin (128×128 → 4×4)")
    print(f"{'='*70}")

    # Sanity check: block-spin 128→64 vs direct construction
    print("\n  [Sanity Check: block-spin 128→64 vs direct construction]")
    G128 = make_torus_2d(128)
    G64_bs, n_bs = block_spin_2d(G128, 128)
    G64_direct = make_torus_2d(64)

    n_bs_nodes = G64_bs.number_of_nodes()
    n_direct_nodes = G64_direct.number_of_nodes()
    assert n_bs_nodes == n_direct_nodes == 4096, \
        f"Node count mismatch: bs={n_bs_nodes}, direct={n_direct_nodes}"

    deg_bs = sorted(set(dict(G64_bs.degree()).values()))
    deg_direct = sorted(set(dict(G64_direct.degree()).values()))
    assert deg_bs == deg_direct == [4], \
        f"Degree mismatch: bs={deg_bs}, direct={deg_direct}"

    # Check diameter
    d_bs = bfs_distances(G64_bs, 0)
    d_direct = bfs_distances(G64_direct, 0)
    diam_bs = max(d_bs.values())
    diam_direct = max(d_direct.values())
    assert diam_bs == diam_direct, f"Diameter mismatch: bs={diam_bs}, direct={diam_direct}"

    print(f"  ✓ Block-spin 128→64 matches direct construction")
    print(f"    Nodes: {n_bs_nodes}, Degrees: {deg_bs}, Diameter: {diam_bs}")

    # Build hierarchy using direct construction (simpler, identical result)
    K = 7  # 2^7 = 128
    levels_data = []
    for level in range(K - 1):  # levels 0 through 5
        side = 2 ** (K - level)
        if side < 4:
            break
        G = make_torus_2d(side)
        data = measure_level(G, side, dim=2, level=level)
        levels_data.append(data)

    # Monotonicity check
    mono_results, all_pass, worst = check_monotonicity(levels_data, sigma=3.0)
    print(f"\n  P7 Monotonicity (σ=3.0): {'PASS' if all_pass else 'FAIL'}")
    print(f"  Worst ΔPR: {worst:+.4f} (tolerance: {EPSILON})")

    # Plots
    plot_deff_vs_level(levels_data, 2.0,
                       'D_eff vs Coarse-Graining Level — 2D Torus (2×2 blocking)',
                       'coarsegrain_2d_deff_vs_level.png')
    plot_sv_profiles(levels_data,
                     'Fisher SV Profiles — 2D Torus Coarse-Graining Hierarchy',
                     'coarsegrain_2d_sv_profiles.png')
    plot_pr_vs_sigma_by_level(levels_data,
                              'PR vs σ at Each Coarse-Graining Level — 2D Torus',
                              'coarsegrain_2d_pr_vs_sigma_by_level.png')

    return levels_data, mono_results, all_pass


# ============================================================================
# Test 2: 2D Torus, 3×3 Block-Spin
# ============================================================================

def run_test_2d_3x3():
    print(f"\n{'='*70}")
    print("  TEST 2: 2D Torus, 3×3 Block-Spin (81×81 → 9×9)")
    print(f"{'='*70}")

    # Sanity check: block-spin 81→27 vs direct construction
    print("\n  [Sanity Check: block-spin 81→27 vs direct construction]")
    G81 = make_torus_2d(81)
    G27_bs, n_bs = block_spin_2d_3x3(G81, 81)
    G27_direct = make_torus_2d(27)

    n_bs_nodes = G27_bs.number_of_nodes()
    n_direct_nodes = G27_direct.number_of_nodes()
    assert n_bs_nodes == n_direct_nodes == 729, \
        f"Node count mismatch: bs={n_bs_nodes}, direct={n_direct_nodes}"

    deg_bs = sorted(set(dict(G27_bs.degree()).values()))
    deg_direct = sorted(set(dict(G27_direct.degree()).values()))
    assert deg_bs == deg_direct == [4], \
        f"Degree mismatch: bs={deg_bs}, direct={deg_direct}"

    print(f"  ✓ Block-spin 81→27 matches direct construction")
    print(f"    Nodes: {n_bs_nodes}, Degrees: {deg_bs}")

    # Build hierarchy
    sides = [81, 27, 9]  # 3^4, 3^3, 3^2
    levels_data = []
    for level, side in enumerate(sides):
        G = make_torus_2d(side)
        data = measure_level(G, side, dim=2, level=level)
        levels_data.append(data)

    # Monotonicity check
    mono_results, all_pass, worst = check_monotonicity(levels_data, sigma=3.0)
    print(f"\n  P7 Monotonicity (σ=3.0): {'PASS' if all_pass else 'FAIL'}")
    print(f"  Worst ΔPR: {worst:+.4f} (tolerance: {EPSILON})")

    # Plots
    plot_deff_vs_level(levels_data, 2.0,
                       'D_eff vs Coarse-Graining Level — 2D Torus (3×3 blocking)',
                       'coarsegrain_2d_3x3_comparison.png')

    return levels_data, mono_results, all_pass


# ============================================================================
# Test 3: 3D Torus, 2×2×2 Block-Spin
# ============================================================================

def run_test_3d_2x2x2():
    print(f"\n{'='*70}")
    print("  TEST 3: 3D Torus, 2×2×2 Block-Spin (32³ → 4³)")
    print(f"{'='*70}")

    # Sanity check: block-spin 32→16 vs direct construction
    print("\n  [Sanity Check: block-spin 32→16 vs direct construction]")
    G32 = make_torus_3d(32)
    G16_bs, n_bs = block_spin_3d(G32, 32)
    G16_direct = make_torus_3d(16)

    n_bs_nodes = G16_bs.number_of_nodes()
    n_direct_nodes = G16_direct.number_of_nodes()
    assert n_bs_nodes == n_direct_nodes == 4096, \
        f"Node count mismatch: bs={n_bs_nodes}, direct={n_direct_nodes}"

    deg_bs = sorted(set(dict(G16_bs.degree()).values()))
    deg_direct = sorted(set(dict(G16_direct.degree()).values()))
    assert deg_bs == deg_direct == [6], \
        f"Degree mismatch: bs={deg_bs}, direct={deg_direct}"

    print(f"  ✓ Block-spin 32³→16³ matches direct construction")
    print(f"    Nodes: {n_bs_nodes}, Degrees: {deg_bs}")

    # Build hierarchy
    levels_data = []
    for level, side in enumerate([32, 16, 8, 4]):
        G = make_torus_3d(side)
        data = measure_level(G, side, dim=3, level=level)
        levels_data.append(data)

    # Monotonicity check
    mono_results, all_pass, worst = check_monotonicity(levels_data, sigma=3.0)
    print(f"\n  P7 Monotonicity (σ=3.0): {'PASS' if all_pass else 'FAIL'}")
    print(f"  Worst ΔPR: {worst:+.4f} (tolerance: {EPSILON})")

    # Plots
    plot_deff_vs_level(levels_data, 3.0,
                       'D_eff vs Coarse-Graining Level — 3D Torus (2×2×2 blocking)',
                       'coarsegrain_3d_deff_vs_level.png')
    plot_sv_profiles(levels_data,
                     'Fisher SV Profiles — 3D Torus Coarse-Graining Hierarchy',
                     'coarsegrain_3d_sv_profiles.png')
    plot_pr_vs_sigma_by_level(levels_data,
                              'PR vs σ at Each Coarse-Graining Level — 3D Torus',
                              'coarsegrain_3d_pr_vs_sigma_by_level.png')

    return levels_data, mono_results, all_pass


# ============================================================================
# Report Generation
# ============================================================================

def format_mono_table(mono_results, dim_label):
    lines = []
    lines.append(f"### {dim_label}")
    lines.append("")
    lines.append("| Level | Size | Growth | Spectral | Fisher Rank | Fisher PR (σ=3) | ΔPR | Status |")
    lines.append("|-------|------|--------|----------|-------------|-----------------|-----|--------|")
    for r in mono_results:
        side = r['side_n']
        dim = 2 if '2D' in dim_label or '3×3' in dim_label else 3
        size = f"{side}{'²' if dim == 2 else '³'}"
        g = f"{r['growth']:.3f}" if not np.isnan(r['growth']) else "N/A"
        s = f"{r['spectral']:.3f}" if not np.isnan(r['spectral']) else "N/A"
        rank = f"{r['rank']:.2f}" if not np.isnan(r['rank']) else "N/A"
        pr = f"{r['pr']:.4f}" if not np.isnan(r['pr']) else "N/A"
        delta = f"{r['delta_pr']:+.4f}" if r['delta_pr'] is not None else "—"
        lines.append(f"| {r['level']} | {size} | {g} | {s} | {rank} | {pr} | {delta} | {r['status']} |")
    lines.append("")
    return "\n".join(lines)


def generate_report(test1, test2, test3, total_time):
    lines = []
    lines.append("# DS Phase 1: Coarse-Graining Test (P7 / DPI) Results")
    lines.append(f"## Date: {datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    lines.append(f"## Runtime: {total_time:.1f} seconds")
    lines.append("")

    # Test 1
    data1, mono1, pass1 = test1
    lines.append(format_mono_table(mono1, "2D Torus, 2×2 Blocking (128² → 4²)"))
    lines.append(f"**P7 Verdict: {'PASS' if pass1 else 'FAIL'}**")
    lines.append("")

    # Test 2
    data2, mono2, pass2 = test2
    lines.append(format_mono_table(mono2, "2D Torus, 3×3 Blocking (81² → 9²)"))
    lines.append(f"**P7 Verdict: {'PASS' if pass2 else 'FAIL'}**")
    lines.append("")

    # Test 3
    data3, mono3, pass3 = test3
    lines.append(format_mono_table(mono3, "3D Torus, 2×2×2 Blocking (32³ → 4³)"))
    lines.append(f"**P7 Verdict: {'PASS' if pass3 else 'FAIL'}**")
    lines.append("")

    # Overall
    lines.append("## Overall P7 Assessment")
    lines.append("")
    all_pass = pass1 and pass2 and pass3
    if all_pass:
        lines.append("**ALL TESTS PASS.** D_eff is monotonically non-increasing under coarse-graining")
        lines.append("across all tested configurations (2D/2×2, 2D/3×3, 3D/2×2×2).")
    else:
        fails = []
        if not pass1: fails.append("2D/2×2")
        if not pass2: fails.append("2D/3×3")
        if not pass3: fails.append("3D/2×2×2")
        lines.append(f"**FAILURES in: {', '.join(fails)}**")
    lines.append("")

    # Sanity checks
    lines.append("## Block-Spin Sanity Checks")
    lines.append("- 2D 128→64 (2×2): ✓ Matches direct torus construction")
    lines.append("- 2D 81→27 (3×3): ✓ Matches direct torus construction")
    lines.append("- 3D 32→16 (2×2×2): ✓ Matches direct torus construction")
    lines.append("")

    report = "\n".join(lines)
    path = os.path.join(OUTPUT_DIR, "PHASE1_COARSEGRAIN_RESULTS.md")
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

    test1 = run_test_2d_2x2()
    test2 = run_test_2d_3x3()
    test3 = run_test_3d_2x2x2()

    total_time = time.time() - t_start
    report = generate_report(test1, test2, test3, total_time)

    print("\n" + "=" * 70)
    print(f"  COARSE-GRAINING TEST COMPLETE — {total_time:.1f}s")
    print("=" * 70)
    print(report)


if __name__ == '__main__':
    main()
