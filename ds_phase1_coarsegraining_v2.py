#!/usr/bin/env python3
"""
DS Phase 1: Coarse-Graining with Scale-Aware Sigma (v2)
========================================================
Re-runs the block-spin coarse-graining hierarchy with σ scaled proportionally
to lattice side length, so σ/side = constant across all levels.

This resolves the PR divergence seen in v1 (where fixed σ on shrinking lattice
caused PR to increase, appearing to violate P7/DPI).

2D hierarchy: 128² → 64² → 32² → 16² → 8² → 4²
3D hierarchy: 32³ → 16³ → 8³ → 4³

Three sigma schedules:
  B: σ/side = 12.5%  (primary — widest reliable range)
  C: σ/side = 6.25%  (robustness check)
  Fixed: σ = 3.0     (original — for comparison)
"""

import os
import time
import datetime
from collections import deque

import numpy as np
import networkx as nx
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

OUTPUT_DIR = "phase1_results"
np.random.seed(42)

SIGMA_MIN_RELIABLE = 2.0  # below this, Fisher results are unreliable

# ============================================================================
# Graph Construction
# ============================================================================

def make_torus_2d(n):
    G = nx.Graph()
    for i in range(n):
        for j in range(n):
            v = i * n + j
            G.add_node(v)
            G.add_edge(v, i * n + (j + 1) % n)
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


def block_spin_3d(G_fine, n_fine):
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
# BFS
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


# ============================================================================
# Growth Dimension
# ============================================================================

def compute_growth_dimension(G, n_samples=20):
    nodes = list(G.nodes())
    n = len(nodes)
    if n < 10:
        return {'fit_dimension': float('nan'), 'r_squared': 0, 'reliable': False}

    samples = [nodes[i] for i in np.random.choice(n, size=min(n_samples, n), replace=False)]
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
        return {'fit_dimension': float('nan'), 'r_squared': 0, 'reliable': False}

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
    return {'fit_dimension': slope, 'r_squared': r_squared, 'reliable': reliable}


# ============================================================================
# Spectral Dimension (Analytical)
# ============================================================================

def torus_eigenvalues_exact(side_n, dim):
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
    return np.sort(np.array(eigs))


def compute_spectral_dimension(side_n, dim):
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
    spectral_dim = 2.0 * coeffs[0]
    predicted = np.polyval(coeffs, log_lam)
    ss_res = np.sum((log_N - predicted) ** 2)
    ss_tot = np.sum((log_N - np.mean(log_N)) ** 2)
    r_squared = 1.0 - ss_res / ss_tot if ss_tot > 0 else 0

    reliable = r_squared > 0.95
    return {'spectral_dimension': spectral_dim, 'r_squared': r_squared, 'reliable': reliable}


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

    if len(sv_norm) > 1:
        ratios = sv_norm[1:] / np.maximum(sv_norm[:-1], 1e-15)
        gap_rank = int(np.argmin(ratios) + 1)
    else:
        gap_rank = 1

    pr = (np.sum(sv)) ** 2 / np.sum(sv ** 2) if np.sum(sv ** 2) > 0 else 0

    return {'sv_norm': sv_norm, 'gap_rank': gap_rank, 'pr': pr}


def compute_fisher_stats(G, sigma, n_samples=20):
    nodes = list(G.nodes())
    n = len(nodes)
    if n < 5:
        return {'mean_rank': float('nan'), 'mean_pr': float('nan'),
                'std_pr': 0, 'reliable': False}

    n_s = min(n_samples, n)
    samples = [nodes[i] for i in np.random.choice(n, size=n_s, replace=False)]

    ranks, prs = [], []
    for v in samples:
        r = compute_fisher_vertex(G, v, nodes, sigma)
        if r is not None:
            ranks.append(r['gap_rank'])
            prs.append(r['pr'])

    if len(prs) == 0:
        return {'mean_rank': float('nan'), 'mean_pr': float('nan'),
                'std_pr': 0, 'reliable': False}

    return {
        'mean_rank': np.mean(ranks),
        'mean_pr': np.mean(prs),
        'std_pr': np.std(prs),
        'reliable': len(prs) >= n_s // 2,
        'n_samples': len(prs),
        'ranks': ranks,
        'prs': prs,
    }


# ============================================================================
# Run hierarchy
# ============================================================================

def run_hierarchy(dim, start_side, schedules, n_fisher_samples=20):
    """
    Run full coarse-graining hierarchy.

    schedules: dict of {name: ratio} where σ = ratio * side at each level.
               Also include {'fixed': sigma_value} for constant-σ comparison.
    """
    print(f"\n{'='*70}")
    print(f"  {dim}D Coarse-Graining Hierarchy (start={start_side})")
    print(f"{'='*70}")

    # Build hierarchy
    if dim == 2:
        G = make_torus_2d(start_side)
    else:
        G = make_torus_3d(start_side)

    hierarchy = [(G, start_side)]
    side = start_side
    while side > 2:
        if side % 2 != 0:
            break
        if dim == 2:
            G_coarse, side_coarse = block_spin_2d(hierarchy[-1][0], side)
        else:
            G_coarse, side_coarse = block_spin_3d(hierarchy[-1][0], side)
        hierarchy.append((G_coarse, side_coarse))
        side = side_coarse

    n_levels = len(hierarchy)
    print(f"  Levels: {n_levels} ({' → '.join([str(s) for _, s in hierarchy])})")

    # Results storage
    results = {
        'dim': dim,
        'levels': [],
        'sides': [],
        'growth_dims': [],
        'spectral_dims': [],
    }

    # Add storage for each schedule
    for sname in schedules:
        results[f'fisher_{sname}'] = {
            'sigmas': [], 'ranks': [], 'prs': [], 'reliable': []
        }

    # Measure at each level
    for level, (G_level, side_n) in enumerate(hierarchy):
        n_verts = G_level.number_of_nodes()
        print(f"\n  Level {level}: side={side_n}, vertices={n_verts}")
        results['levels'].append(level)
        results['sides'].append(side_n)

        # Growth dimension
        growth = compute_growth_dimension(G_level, n_samples=n_fisher_samples)
        gd = growth['fit_dimension']
        results['growth_dims'].append(gd)
        if not np.isnan(gd):
            print(f"    Growth dim: {gd:.4f} (R²={growth.get('r_squared',0):.4f})")
        else:
            print(f"    Growth dim: N/A (too small)")

        # Spectral dimension
        spectral = compute_spectral_dimension(side_n, dim)
        sd = spectral.get('spectral_dimension', float('nan'))
        results['spectral_dims'].append(sd)
        if not np.isnan(sd):
            print(f"    Spectral dim: {sd:.4f} (R²={spectral.get('r_squared',0):.4f})")
        else:
            print(f"    Spectral dim: N/A")

        # Fisher under each schedule
        for sname, sval in schedules.items():
            if sname == 'fixed':
                sigma = sval  # constant σ
            else:
                sigma = sval * side_n  # ratio × side

            is_reliable = sigma >= SIGMA_MIN_RELIABLE
            fisher = compute_fisher_stats(G_level, sigma, n_samples=n_fisher_samples)

            results[f'fisher_{sname}']['sigmas'].append(sigma)
            results[f'fisher_{sname}']['ranks'].append(fisher['mean_rank'])
            results[f'fisher_{sname}']['prs'].append(fisher['mean_pr'])
            results[f'fisher_{sname}']['reliable'].append(is_reliable and fisher.get('reliable', False))

            rel_tag = "✓" if is_reliable else "✗ (σ<2)"
            if not np.isnan(fisher['mean_pr']):
                print(f"    Fisher [{sname}] σ={sigma:.3f}: rank={fisher['mean_rank']:.2f} "
                      f"PR={fisher['mean_pr']:.3f} (±{fisher['std_pr']:.3f}) [{rel_tag}]")
            else:
                print(f"    Fisher [{sname}] σ={sigma:.3f}: DEGENERATE [{rel_tag}]")

    return results


# ============================================================================
# Plotting
# ============================================================================

def plot_combined_chart(results, filename):
    """
    The key deliverable: D_eff vs coarse-graining level with both
    constant-σ and scaled-σ Fisher PR lines.
    """
    dim = results['dim']
    levels = results['levels']
    sides = results['sides']
    true_d = float(dim)

    fig, ax = plt.subplots(figsize=(12, 7))

    # Growth dimension (blue circles)
    gd = np.array(results['growth_dims'])
    mask_gd = ~np.isnan(gd)
    if np.any(mask_gd):
        ax.plot(np.array(levels)[mask_gd], gd[mask_gd], 'bo-', markersize=8,
                linewidth=2, label='Growth Dim', zorder=5)

    # Spectral dimension (red squares)
    sd = np.array(results['spectral_dims'])
    mask_sd = ~np.isnan(sd)
    if np.any(mask_sd):
        ax.plot(np.array(levels)[mask_sd], sd[mask_sd], 'rs-', markersize=8,
                linewidth=2, label='Spectral Dim', zorder=5)

    # Fisher PR — fixed σ (green dashed)
    if 'fisher_fixed' in results:
        fdata = results['fisher_fixed']
        pr_fixed = np.array(fdata['prs'])
        mask_f = ~np.isnan(pr_fixed)
        if np.any(mask_f):
            ax.plot(np.array(levels)[mask_f], pr_fixed[mask_f], 'g^--', markersize=8,
                    linewidth=2, label='Fisher PR (σ=3.0, fixed)', zorder=5, alpha=0.7)

    # Fisher PR — Schedule B (green solid, primary)
    if 'fisher_B' in results:
        fdata = results['fisher_B']
        pr_B = np.array(fdata['prs'])
        reliable_B = np.array(fdata['reliable'])
        mask_B = ~np.isnan(pr_B)
        if np.any(mask_B):
            # Plot reliable points solid, unreliable as open
            rel_mask = mask_B & reliable_B
            unrel_mask = mask_B & ~reliable_B
            if np.any(rel_mask):
                ax.plot(np.array(levels)[rel_mask], pr_B[rel_mask], 'gv-', markersize=10,
                        linewidth=2.5, label='Fisher PR (σ/side=12.5%, scaled)',
                        zorder=6, markerfacecolor='green')
            if np.any(unrel_mask):
                ax.plot(np.array(levels)[unrel_mask], pr_B[unrel_mask], 'gv', markersize=10,
                        linewidth=2.5, zorder=6, markerfacecolor='none', markeredgewidth=2,
                        label='Fisher PR (scaled, σ<2 unreliable)')

    # Fisher PR — Schedule C (orange, robustness check)
    if 'fisher_C' in results:
        fdata = results['fisher_C']
        pr_C = np.array(fdata['prs'])
        reliable_C = np.array(fdata['reliable'])
        mask_C = ~np.isnan(pr_C)
        if np.any(mask_C):
            rel_mask = mask_C & reliable_C
            if np.any(rel_mask):
                ax.plot(np.array(levels)[rel_mask], pr_C[rel_mask], 'D-', color='orange',
                        markersize=8, linewidth=2, label='Fisher PR (σ/side=6.25%, scaled)',
                        zorder=5, alpha=0.8)

    # Fisher Rank — Schedule B (as a separate indicator)
    if 'fisher_B' in results:
        fdata = results['fisher_B']
        ranks_B = np.array(fdata['ranks'])
        reliable_B = np.array(fdata['reliable'])
        rel_mask = ~np.isnan(ranks_B) & reliable_B
        if np.any(rel_mask):
            ax.plot(np.array(levels)[rel_mask], ranks_B[rel_mask], 'k+', markersize=12,
                    markeredgewidth=2, label='Fisher Rank (scaled)', zorder=7)

    # True dimension line
    ax.axhline(true_d, color='gray', linestyle='--', linewidth=1.5,
               label=f'd = {dim}', zorder=1)

    # Unreliable region shading
    if 'fisher_B' in results:
        fdata = results['fisher_B']
        for i, (lev, sigma, rel) in enumerate(zip(levels, fdata['sigmas'], fdata['reliable'])):
            if not rel:
                ax.axvspan(lev - 0.4, lev + 0.4, alpha=0.1, color='red', zorder=0)

    # Formatting
    ax.set_xlabel('Coarse-Graining Level', fontsize=13)
    ax.set_ylabel('Estimated Dimension', fontsize=13)
    ax.set_title(f'{dim}D Torus: D_eff vs Coarse-Graining Level\n'
                 f'(Constant σ vs Scale-Aware σ)', fontsize=14)
    ax.set_xticks(levels)
    ax.set_xticklabels([f'L{l}\n{s}{"²" if dim == 2 else "³"}' for l, s in zip(levels, sides)])
    ax.legend(loc='upper left', fontsize=10)
    ax.grid(True, alpha=0.3)

    # Add sigma annotations for Schedule B
    if 'fisher_B' in results:
        fdata = results['fisher_B']
        for i, (lev, sigma) in enumerate(zip(levels, fdata['sigmas'])):
            ax.annotate(f'σ={sigma:.1f}', xy=(lev, 0.5), fontsize=7,
                       ha='center', color='green', alpha=0.6)

    plt.tight_layout()
    plt.savefig(os.path.join(OUTPUT_DIR, filename), dpi=150)
    plt.close()
    print(f"  Saved {filename}")


# ============================================================================
# Results report
# ============================================================================

def write_results(results_2d, results_3d):
    lines = []
    lines.append("# Phase 1: Coarse-Graining with Scale-Aware Sigma — Results")
    lines.append(f"\nGenerated: {datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    lines.append("")

    for tag, res in [("2D", results_2d), ("3D", results_3d)]:
        dim = res['dim']
        lines.append(f"## {tag} Torus Hierarchy")
        lines.append("")

        # Main table
        lines.append("### All Routes — Comparison Table")
        lines.append("")
        header = "| Level | Side | Growth | Spectral |"
        sep = "|-------|------|--------|----------|"
        for sname in ['fixed', 'B', 'C']:
            key = f'fisher_{sname}'
            if key in res:
                header += f" σ ({sname}) | Rank ({sname}) | PR ({sname}) |"
                sep += "---------|------------|----------|"
        lines.append(header)
        lines.append(sep)

        for i, (level, side) in enumerate(zip(res['levels'], res['sides'])):
            gd = res['growth_dims'][i]
            sd = res['spectral_dims'][i]
            gd_s = f"{gd:.3f}" if not np.isnan(gd) else "N/A"
            sd_s = f"{sd:.3f}" if not np.isnan(sd) else "N/A"
            row = f"| {level} | {side} | {gd_s} | {sd_s} |"

            for sname in ['fixed', 'B', 'C']:
                key = f'fisher_{sname}'
                if key in res:
                    fdata = res[key]
                    sigma = fdata['sigmas'][i]
                    rank = fdata['ranks'][i]
                    pr = fdata['prs'][i]
                    rel = fdata['reliable'][i]
                    rel_mark = "" if rel else "*"
                    rank_s = f"{rank:.1f}{rel_mark}" if not np.isnan(rank) else "N/A"
                    pr_s = f"{pr:.3f}{rel_mark}" if not np.isnan(pr) else "N/A"
                    row += f" {sigma:.3f} | {rank_s} | {pr_s} |"

            lines.append(row)

        lines.append("")
        lines.append("*Asterisk (*) indicates σ < 2.0 — results may be unreliable.*")
        lines.append("")

        # P7 monotonicity check for scaled sigma
        key_B = 'fisher_B'
        if key_B in res:
            fdata = res[key_B]
            prs = fdata['prs']
            reliable = fdata['reliable']

            # Check monotonicity on reliable levels only
            rel_prs = [(i, pr) for i, (pr, rel) in enumerate(zip(prs, reliable))
                       if rel and not np.isnan(pr)]

            if len(rel_prs) >= 2:
                monotone = all(rel_prs[j+1][1] <= rel_prs[j][1] + 0.05
                              for j in range(len(rel_prs) - 1))
                pr_first = rel_prs[0][1]
                pr_last = rel_prs[-1][1]
                delta = pr_last - pr_first

                lines.append(f"### P7 Monotonicity (Schedule B, reliable levels)")
                lines.append(f"- PR at first reliable level: {pr_first:.3f}")
                lines.append(f"- PR at last reliable level: {pr_last:.3f}")
                lines.append(f"- Change: {delta:+.3f}")
                lines.append(f"- Monotonically non-increasing (ε=0.05): **{'PASS' if monotone else 'FAIL'}**")
                lines.append("")

        lines.append("---")
        lines.append("")

    path = os.path.join(OUTPUT_DIR, "PHASE1_COARSEGRAIN_V2_RESULTS.md")
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
    print("DS Phase 1: Coarse-Graining with Scale-Aware Sigma (v2)")
    print("=" * 70)

    # Sigma schedules
    schedules_2d = {
        'fixed': 3.0,    # constant σ = 3.0
        'B': 0.125,      # σ/side = 12.5%
        'C': 0.0625,     # σ/side = 6.25%
    }

    schedules_3d = {
        'fixed': 3.0,
        'B': 0.125,
        'C': 0.0625,
    }

    # 2D hierarchy
    results_2d = run_hierarchy(dim=2, start_side=128, schedules=schedules_2d, n_fisher_samples=20)

    # 3D hierarchy
    results_3d = run_hierarchy(dim=3, start_side=32, schedules=schedules_3d, n_fisher_samples=20)

    # Plots
    print(f"\n{'='*70}")
    print("  GENERATING CHARTS")
    print(f"{'='*70}")

    plot_combined_chart(results_2d, "coarsegrain_2d_deff_vs_level_scaled_sigma.png")
    plot_combined_chart(results_3d, "coarsegrain_3d_deff_vs_level_scaled_sigma.png")

    # Results report
    write_results(results_2d, results_3d)

    elapsed = time.time() - t0
    print(f"\n{'='*70}")
    print(f"  TOTAL TIME: {elapsed:.1f}s")
    print(f"{'='*70}")
