#!/usr/bin/env python3
"""
DS Phase 1 Extension v4: Gasket L=10 + Carpet L=6 with sigma to 500
=====================================================================
The decisive test. Can the carpet PR cross d_H = 1.893 at sigma ~ 270-500?

Gasket L=10: 88,575 vertices, diameter ~1024. sigma=500 is ~49% of diameter.
Carpet L=6:  262,144 vertices, diameter ~1458. sigma=500 is ~34% of diameter.

Also re-runs L=9 gasket and L=5 carpet at extended sigma for comparison.
"""

import os
import sys
import time
import datetime
from collections import deque

import numpy as np
import networkx as nx
from scipy.optimize import curve_fit
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

OUTPUT_DIR = "phase1_results"
np.random.seed(42)

# Reference dimensions
GASKET_D_H = np.log(3) / np.log(2)        # ≈ 1.5849
GASKET_D_S = 2 * np.log(3) / np.log(5)    # ≈ 1.3652
CARPET_D_H = np.log(8) / np.log(3)        # ≈ 1.8928
CARPET_D_S = 1.805

# Extended sigma sweep — push to 500
SIGMA_SWEEP = [1.5, 2.0, 3.0, 5.0, 8.0, 12.0, 16.0, 20.0, 25.0, 30.0,
               40.0, 50.0, 65.0, 80.0, 100.0, 130.0, 170.0, 220.0, 280.0, 360.0, 450.0]

# ============================================================================
# BFS — optimized for large graphs using adjacency list
# ============================================================================

def bfs_distances_fast(adj, source, n):
    """BFS using pre-built adjacency list. Returns distance array."""
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
    """Convert networkx graph to adjacency list with integer node IDs."""
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
# Fisher Information — optimized for large graphs
# ============================================================================

def compute_fisher_single_vertex_fast(adj, v0_idx, n, sigma):
    """Compute Fisher PR at vertex v0_idx using pre-built adjacency list."""
    dist_v0 = bfs_distances_fast(adj, v0_idx, n)

    log_unnorm_v0 = -dist_v0.astype(float) / sigma
    # Handle unreachable nodes (shouldn't happen on connected graph)
    log_unnorm_v0[dist_v0 < 0] = -1000.0
    log_Z_v0 = np.logaddexp.reduce(log_unnorm_v0)
    log_p_v0 = log_unnorm_v0 - log_Z_v0
    p_v0 = np.exp(log_p_v0)

    neighbors = adj[v0_idx]
    k = len(neighbors)
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

    if np.sum(sv ** 2) > 0:
        pr = (np.sum(sv)) ** 2 / np.sum(sv ** 2)
    else:
        pr = 0.0

    return {'sv_norm': sv_norm, 'participation_ratio': pr}


def fisher_sigma_sweep_fast(adj, n, sigmas, sample_indices, n_samples=20):
    """Run Fisher PR sweep using fast adjacency list BFS."""
    samples = sample_indices[:n_samples]
    results = {}

    for sigma in sigmas:
        prs, svs = [], []
        n_degenerate = 0
        for v_idx in samples:
            r = compute_fisher_single_vertex_fast(adj, v_idx, n, sigma)
            if r is None:
                n_degenerate += 1
            else:
                prs.append(r['participation_ratio'])
                svs.append(r['sv_norm'])

        if len(prs) < len(samples) // 2:
            print(f"    sigma={sigma:7.1f}: FIM DEGENERATE ({n_degenerate}/{len(samples)} samples)")
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
                'mean_sv': mean_sv,
                'n_degenerate': n_degenerate,
                'n_valid': len(prs),
            }
            deg_str = f" ({n_degenerate} degen)" if n_degenerate > 0 else ""
            print(f"    sigma={sigma:7.1f}: PR={results[sigma]['mean_pr']:.4f} "
                  f"± {results[sigma]['std_pr']:.4f}{deg_str}")

    return results


# ============================================================================
# Sierpinski Gasket Construction
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
    for k in range(L):
        dx = (x // (3**k)) % 3
        dy = (y // (3**k)) % 3
        if dx == 1 and dy == 1:
            return True
    return False


def make_sierpinski_carpet(L):
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


# ============================================================================
# Growth Dimension
# ============================================================================

def compute_growth_dimension(adj, n, sample_indices, n_samples=30):
    """Compute growth dimension using fast BFS."""
    samples = sample_indices[:n_samples]

    # Get maxR from first sample
    d0 = bfs_distances_fast(adj, samples[0], n)
    maxR = int(np.max(d0[d0 >= 0])) // 2

    all_volumes = []
    for i, v_idx in enumerate(samples):
        dist = bfs_distances_fast(adj, v_idx, n)
        volumes = np.array([np.sum(dist[dist >= 0] <= r) for r in range(maxR + 1)])
        all_volumes.append(volumes)
        if (i + 1) % 10 == 0:
            print(f"    Growth BFS {i+1}/{len(samples)}")

    mean_vol = np.mean(all_volumes, axis=0)

    rMin = max(3, int(maxR * 0.20))
    rMax = int(maxR * 0.60)
    radii = np.arange(rMin, rMax + 1, dtype=float)
    log_r = np.log(radii)
    log_v = np.log(mean_vol[rMin:rMax + 1].astype(float))
    coeffs = np.polyfit(log_r, log_v, 1)
    slope = coeffs[0]
    predicted = np.polyval(coeffs, log_r)
    ss_res = np.sum((log_v - predicted) ** 2)
    ss_tot = np.sum((log_v - np.mean(log_v)) ** 2)
    r_squared = 1.0 - ss_res / ss_tot

    deltas = []
    for r in range(1, maxR):
        if mean_vol[r] > 0 and mean_vol[r + 1] > 0:
            d = (np.log(float(mean_vol[r + 1])) - np.log(float(mean_vol[r]))) / \
                (np.log(float(r + 1)) - np.log(float(r)))
            deltas.append((r, d))

    return {
        'fit_dimension': slope, 'r_squared': r_squared,
        'rMin': rMin, 'rMax': rMax, 'maxR': maxR,
        'deltas': deltas, 'mean_volumes': mean_vol,
    }


# ============================================================================
# Extrapolation
# ============================================================================

def power_law_model(sigma, d_inf, A, alpha):
    return d_inf + A * sigma**(-alpha)


def fit_and_extrapolate(sigmas, prs, label, d_H, d_S):
    """Fit PR(sigma) = d_inf + A*sigma^(-alpha) and extrapolate."""
    print(f"\n  --- Extrapolation: {label} ---")
    print(f"  d_H = {d_H:.4f}, d_S = {d_S:.4f}")

    results = {}
    for start_idx, fit_label in [(0, 'all'), (4, 'sigma>=8'), (8, 'sigma>=30')]:
        if start_idx >= len(sigmas):
            continue
        s = np.array(sigmas[start_idx:])
        p = np.array(prs[start_idx:])
        if len(s) < 4:
            continue
        try:
            popt, pcov = curve_fit(power_law_model, s, p,
                                   p0=[d_H, 2.0, 0.5], maxfev=10000,
                                   bounds=([0.0, 0, 0.01], [5.0, 50.0, 5.0]))
            perr = np.sqrt(np.diag(pcov))
            d_inf, A, alpha = popt
            d_inf_err = perr[0]

            residuals = p - power_law_model(s, *popt)
            ss_res = np.sum(residuals**2)
            ss_tot = np.sum((p - np.mean(p))**2)
            r2 = 1 - ss_res/ss_tot

            pred = {sig: power_law_model(sig, *popt) for sig in [500, 1000, 5000]}

            print(f"    Fit ({fit_label}): d_inf = {d_inf:.4f} ± {d_inf_err:.4f}, "
                  f"A = {A:.4f}, alpha = {alpha:.4f}, R² = {r2:.6f}")
            print(f"      σ=500→{pred[500]:.4f}, σ=1000→{pred[1000]:.4f}, σ=5000→{pred[5000]:.4f}")

            if d_inf < d_H:
                sigma_cross = (A / (d_H - d_inf))**(1/alpha)
                print(f"      PR crosses d_H at σ ≈ {sigma_cross:.1f}")
            else:
                excess = d_inf - d_H
                print(f"      d_inf ABOVE d_H by {excess:.4f}")

            results[fit_label] = {
                'd_inf': d_inf, 'd_inf_err': d_inf_err,
                'A': A, 'alpha': alpha, 'r2': r2, 'predictions': pred
            }
        except Exception as e:
            print(f"    Fit ({fit_label}): FAILED - {e}")
            results[fit_label] = None

    return results


# ============================================================================
# Plotting
# ============================================================================

def plot_pr_vs_sigma(sweep_data, title, filename, ref_lines, extrapolation=None):
    """Plot PR vs sigma with reference lines and optional extrapolation curves."""
    fig, ax = plt.subplots(figsize=(14, 8))
    colors = ['#2196F3', '#FF5722', '#4CAF50', '#9C27B0', '#FF9800', '#795548']
    markers = ['o', 's', '^', 'D', 'v', 'p']

    for i, (label, sweep) in enumerate(sweep_data.items()):
        sigmas, prs, errs = [], [], []
        for sigma in sorted(sweep.keys()):
            r = sweep[sigma]
            if not r.get('degenerate', False):
                sigmas.append(sigma)
                prs.append(r['mean_pr'])
                errs.append(r['std_pr'])
        ax.errorbar(sigmas, prs, yerr=errs, marker=markers[i % len(markers)],
                    color=colors[i % len(colors)], linewidth=2, markersize=5,
                    capsize=3, label=label)

    # Extrapolation curves
    if extrapolation:
        for key, fit in extrapolation.items():
            if fit and 'sigma>=8' in fit and fit['sigma>=8']:
                f = fit['sigma>=8']
                sigma_ext = np.logspace(0, np.log10(5000), 300)
                pr_ext = power_law_model(sigma_ext, f['d_inf'], f['A'], f['alpha'])
                color = '#2196F3' if 'gasket' in key.lower() else '#FF5722'
                ax.plot(sigma_ext, pr_ext, ':', color=color, linewidth=1, alpha=0.5,
                        label=f'{key} extrap → d∞={f["d_inf"]:.3f}')

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
    ax.legend(fontsize=8, loc='best')
    ax.grid(True, alpha=0.3)
    path = os.path.join(OUTPUT_DIR, filename)
    fig.savefig(path, dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f"    Saved {path}")


def plot_sv_panel(sweep, prefix, output_dir, title=None):
    sigmas = sorted([s for s in sweep if not sweep[s].get('degenerate', False)])
    n = len(sigmas)
    if n == 0:
        return
    cols = min(5, n)
    rows = (n + cols - 1) // cols
    fig, axes = plt.subplots(rows, cols, figsize=(3.5 * cols, 3 * rows))
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
        ax.set_title(f'σ={sigma}', fontsize=9)
        ax.set_ylim(0, 1.15)
        ax.set_xticks(x)
        ax.grid(True, alpha=0.3, axis='y')

    for idx in range(n, rows * cols):
        row, col = idx // cols, idx % cols
        axes[row, col].set_visible(False)

    fig.suptitle(title or f'SV Profiles — {prefix}', fontsize=12)
    fig.tight_layout()
    path = os.path.join(output_dir, f'{prefix}_sv_profiles_full_sweep.png')
    fig.savefig(path, dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f"    Saved {path}")


def plot_growth_delta(deltas, ref_dims, rMin, rMax, prefix, output_dir):
    fig, ax = plt.subplots(figsize=(10, 6))
    rs = [r for r, d in deltas]
    ds = [d for r, d in deltas]
    ax.plot(rs, ds, 'b.-', alpha=0.7, markersize=3)
    for i, (label, val) in enumerate(ref_dims.items()):
        ax.axhline(y=val, color=['r', 'g', 'purple'][i % 3], linestyle='--',
                   linewidth=2, label=f'{label} = {val:.4f}')
    ax.axvspan(rMin, rMax, alpha=0.15, color='green', label=f'Fit [{rMin},{rMax}]')
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
    ax.plot(np.log(valid.astype(float)), np.log(mean_volumes[1:].astype(float)),
            'b.', alpha=0.5, markersize=2)
    fit_r = np.arange(rMin, rMax + 1, dtype=float)
    fit_log_r = np.log(fit_r)
    fit_log_v = np.log(mean_volumes[rMin:rMax + 1].astype(float))
    coeffs = np.polyfit(fit_log_r, fit_log_v, 1)
    ax.plot(fit_log_r, np.polyval(coeffs, fit_log_r), 'r-', linewidth=2,
            label=f'slope = {slope:.4f} (ref = {ref_dim:.4f})')
    ax.set_xlabel('log(r)', fontsize=12)
    ax.set_ylabel('log(V(r))', fontsize=12)
    ax.set_title(f'Growth Dimension log-log — {prefix}', fontsize=14)
    ax.legend(fontsize=11)
    ax.grid(True, alpha=0.3)
    path = os.path.join(output_dir, f'{prefix}_growth_loglog.png')
    fig.savefig(path, dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f"    Saved {path}")


# ============================================================================
# Main
# ============================================================================

def main():
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    t_start = time.time()

    all_results = {}

    # ==================================================================
    # GASKET L=9 and L=10
    # ==================================================================
    print(f"\n{'='*70}")
    print("  GASKET: L=9 and L=10 — Sigma Sweep to 450")
    print(f"{'='*70}")

    for level in [9, 10]:
        print(f"\n--- Sierpinski Gasket L={level} ---")
        t0 = time.time()
        G, corners = make_sierpinski_gasket(level)
        n_nodes = G.number_of_nodes()
        expected = (3 ** (level + 1) + 3) // 2
        assert n_nodes == expected, f"Expected {expected}, got {n_nodes}"
        build_time = time.time() - t0
        print(f"  Built: {n_nodes} vertices, {build_time:.1f}s")

        # Convert to fast adjacency list
        t0 = time.time()
        adj, nodes, node_to_idx, n = build_adjacency_list(G)
        degrees = [len(adj[i]) for i in range(n)]
        interior_idx = [i for i in range(n) if degrees[i] == 4]
        print(f"  Adjacency list: {n} vertices, {len(interior_idx)} interior, {time.time()-t0:.1f}s")

        # Diameter
        d0 = bfs_distances_fast(adj, 0, n)
        far = int(np.argmax(d0))
        d1 = bfs_distances_fast(adj, far, n)
        diameter = int(np.max(d1))
        print(f"  Diameter: {diameter}")

        # Free the networkx graph to save memory
        del G

        # Sample interior vertices
        sample_idx = [interior_idx[i] for i in
                      np.random.choice(len(interior_idx), size=min(20, len(interior_idx)), replace=False)]

        # Sigma sweep
        print(f"\n  [Sigma sweep — {len(SIGMA_SWEEP)} values]")
        t0 = time.time()
        sweep = fisher_sigma_sweep_fast(adj, n, SIGMA_SWEEP, sample_idx, n_samples=20)
        print(f"  Sweep time: {time.time()-t0:.1f}s")

        all_results[f'gasket_L{level}'] = {
            'sweep': sweep, 'level': level, 'n_nodes': n_nodes,
            'diameter': diameter, 'n_interior': len(interior_idx)
        }

        plot_sv_panel(sweep, f'sierpinski_L{level}', OUTPUT_DIR,
                      f'SV Profiles — Gasket L={level} (σ to 450)')

        # Free memory
        del adj

    # Gasket comparison plot
    gasket_sweep_data = {}
    for key in ['gasket_L9', 'gasket_L10']:
        gasket_sweep_data[key.replace('gasket_', 'Gasket ')] = all_results[key]['sweep']

    plot_pr_vs_sigma(
        gasket_sweep_data,
        'PR vs σ — Gasket L=9 vs L=10 (Extended to σ=450)',
        'sierpinski_pr_vs_sigma_v4.png',
        {'Gasket d_H': GASKET_D_H, 'Gasket d_S': GASKET_D_S}
    )

    # ==================================================================
    # CARPET L=5 and L=6
    # ==================================================================
    print(f"\n{'='*70}")
    print("  CARPET: L=5 and L=6 — Sigma Sweep to 450")
    print(f"{'='*70}")

    for level in [5, 6]:
        prefix = f'carpet_L{level}'
        print(f"\n--- Sierpinski Carpet L={level} ---")
        t0 = time.time()
        G = make_sierpinski_carpet(level)
        n_nodes = G.number_of_nodes()
        expected = 8 ** level
        assert n_nodes == expected, f"Expected {expected}, got {n_nodes}"
        build_time = time.time() - t0
        print(f"  Built: {n_nodes} vertices, {build_time:.1f}s")

        # Convert to fast adjacency list
        t0 = time.time()
        adj, nodes, node_to_idx, n = build_adjacency_list(G)
        degrees_arr = [len(adj[i]) for i in range(n)]
        deg_counts = {}
        for d in degrees_arr:
            deg_counts[d] = deg_counts.get(d, 0) + 1
        deg4_idx = [i for i in range(n) if degrees_arr[i] == 4]
        print(f"  Adjacency list: {n} verts, {len(deg4_idx)} interior (deg 4), {time.time()-t0:.1f}s")
        print(f"  Degree dist: {dict(sorted(deg_counts.items()))}")

        # Diameter
        d0 = bfs_distances_fast(adj, 0, n)
        far = int(np.argmax(d0))
        d1 = bfs_distances_fast(adj, far, n)
        diameter = int(np.max(d1))
        print(f"  Diameter: {diameter}")

        all_results[f'carpet_L{level}'] = {
            'level': level, 'n_nodes': n_nodes, 'diameter': diameter,
            'deg_counts': deg_counts, 'n_deg4': len(deg4_idx)
        }

        sample_deg4 = [deg4_idx[i] for i in
                       np.random.choice(len(deg4_idx), size=min(30, len(deg4_idx)), replace=False)]

        # Growth dimension for L=6
        if level == 6:
            print(f"\n  [Growth Dimension]")
            t0 = time.time()
            growth = compute_growth_dimension(adj, n, sample_deg4, n_samples=30)
            gate_lo, gate_hi = 1.70, 2.10
            gpass = gate_lo <= growth['fit_dimension'] <= gate_hi
            print(f"  Growth dim: {growth['fit_dimension']:.4f} (R²={growth['r_squared']:.4f})")
            print(f"  Gate [{gate_lo}, {gate_hi}]: {'PASS' if gpass else 'FAIL'}")
            print(f"  Time: {time.time()-t0:.1f}s")
            all_results['carpet_L6']['growth'] = growth
            all_results['carpet_L6']['growth_pass'] = gpass

            ref_dims = {'d_H': CARPET_D_H, 'd_S': CARPET_D_S}
            plot_growth_delta(growth['deltas'], ref_dims, growth['rMin'],
                              growth['rMax'], prefix, OUTPUT_DIR)
            plot_growth_loglog(growth['mean_volumes'], growth['rMin'], growth['rMax'],
                               growth['fit_dimension'], CARPET_D_H, prefix, OUTPUT_DIR)

        # Free networkx graph
        del G

        # Sigma sweep
        print(f"\n  [Fisher sigma sweep L={level}]")
        t0 = time.time()
        sweep = fisher_sigma_sweep_fast(adj, n, SIGMA_SWEEP, sample_deg4, n_samples=20)
        print(f"  Sweep time: {time.time()-t0:.1f}s")
        all_results[f'carpet_L{level}']['sweep'] = sweep

        plot_sv_panel(sweep, prefix, OUTPUT_DIR, f'SV Profiles — Carpet L={level} (σ to 450)')

        del adj

    # Carpet comparison plot
    carpet_sweep_data = {}
    for key in ['carpet_L5', 'carpet_L6']:
        carpet_sweep_data[key.replace('carpet_', 'Carpet ')] = all_results[key]['sweep']

    plot_pr_vs_sigma(
        carpet_sweep_data,
        'PR vs σ — Carpet L=5 vs L=6 (Extended to σ=450)',
        'carpet_pr_vs_sigma_v4.png',
        {'Carpet d_H': CARPET_D_H, 'Carpet d_S': CARPET_D_S}
    )

    # ==================================================================
    # EXTRAPOLATION
    # ==================================================================
    print(f"\n{'='*70}")
    print("  EXTRAPOLATION ANALYSIS")
    print(f"{'='*70}")

    all_extrap = {}

    for sys_key, d_H, d_S in [('gasket_L10', GASKET_D_H, GASKET_D_S),
                                ('gasket_L9', GASKET_D_H, GASKET_D_S),
                                ('carpet_L6', CARPET_D_H, CARPET_D_S),
                                ('carpet_L5', CARPET_D_H, CARPET_D_S)]:
        if sys_key in all_results and 'sweep' in all_results[sys_key]:
            sw = all_results[sys_key]['sweep']
            sigmas, prs = [], []
            for sigma in sorted(sw.keys()):
                if not sw[sigma].get('degenerate', False):
                    sigmas.append(sigma)
                    prs.append(sw[sigma]['mean_pr'])
            fit = fit_and_extrapolate(sigmas, prs, sys_key, d_H, d_S)
            all_extrap[sys_key] = {'sigmas': sigmas, 'prs': prs, 'fits': fit}

    # ==================================================================
    # L vs L AGREEMENT TABLES
    # ==================================================================
    print(f"\n{'='*70}")
    print("  L=9 vs L=10 AGREEMENT (Gasket)")
    print(f"{'='*70}")
    sw9 = all_results.get('gasket_L9', {}).get('sweep', {})
    sw10 = all_results.get('gasket_L10', {}).get('sweep', {})
    print(f"\n  {'sigma':>7}  {'L=9 PR':>8}  {'L=10 PR':>8}  {'Δ':>8}  {'Note':>12}")
    for sigma in SIGMA_SWEEP:
        r9, r10 = sw9.get(sigma, {}), sw10.get(sigma, {})
        if not r9.get('degenerate', True) and not r10.get('degenerate', True):
            p9, p10 = r9['mean_pr'], r10['mean_pr']
            delta = p10 - p9
            note = "✓ CLOSE" if abs(delta) < 0.02 else "DIVERGING" if abs(delta) > 0.05 else "~moderate"
            print(f"  {sigma:7.1f}  {p9:8.4f}  {p10:8.4f}  {delta:+8.4f}  {note:>12}")

    print(f"\n{'='*70}")
    print("  L=5 vs L=6 AGREEMENT (Carpet)")
    print(f"{'='*70}")
    sw5 = all_results.get('carpet_L5', {}).get('sweep', {})
    sw6 = all_results.get('carpet_L6', {}).get('sweep', {})
    print(f"\n  {'sigma':>7}  {'L=5 PR':>8}  {'L=6 PR':>8}  {'Δ':>8}  {'Note':>12}")
    for sigma in SIGMA_SWEEP:
        r5, r6 = sw5.get(sigma, {}), sw6.get(sigma, {})
        if not r5.get('degenerate', True) and not r6.get('degenerate', True):
            p5, p6 = r5['mean_pr'], r6['mean_pr']
            delta = p6 - p5
            note = "✓ CLOSE" if abs(delta) < 0.02 else "DIVERGING" if abs(delta) > 0.05 else "~moderate"
            print(f"  {sigma:7.1f}  {p5:8.4f}  {p6:8.4f}  {delta:+8.4f}  {note:>12}")

    # ==================================================================
    # MONEY PLOT v4
    # ==================================================================
    print(f"\n  [Money Plot v4: Gasket L=10 vs Carpet L=6]")
    money_sweep = {}
    money_sweep['Gasket L=10 (89k)'] = all_results['gasket_L10']['sweep']
    money_sweep['Carpet L=6 (262k)'] = all_results['carpet_L6']['sweep']

    plot_pr_vs_sigma(
        money_sweep,
        'Fisher PR vs σ — Gasket L=10 vs Carpet L=6 (σ to 450)',
        'fractal_pr_vs_sigma_v4.png',
        {
            'Gasket d_H': GASKET_D_H,
            'Gasket d_S': GASKET_D_S,
            'Carpet d_H': CARPET_D_H,
            'Carpet d_S': CARPET_D_S,
        },
        extrapolation=all_extrap
    )

    # Comprehensive multi-level comparison plot
    all_sweep = {}
    for key in ['gasket_L9', 'gasket_L10']:
        all_sweep[key.replace('_', ' ').title()] = all_results[key]['sweep']
    for key in ['carpet_L5', 'carpet_L6']:
        all_sweep[key.replace('_', ' ').title()] = all_results[key]['sweep']

    plot_pr_vs_sigma(
        all_sweep,
        'Fisher PR vs σ — All Systems (v4 Comprehensive)',
        'fractal_pr_vs_sigma_v4_all.png',
        {
            'Gasket d_H': GASKET_D_H,
            'Gasket d_S': GASKET_D_S,
            'Carpet d_H': CARPET_D_H,
            'Carpet d_S': CARPET_D_S,
        }
    )

    # ==================================================================
    # REPORT
    # ==================================================================
    total_time = time.time() - t_start

    lines = []
    lines.append("# DS Phase 1 Extension v4 Results")
    lines.append(f"## Date: {datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    lines.append(f"## Runtime: {total_time:.1f} seconds")
    lines.append("")

    # System specs
    lines.append("## Systems")
    lines.append("")
    lines.append("| System | Vertices | Diameter | Interior (deg 4) |")
    lines.append("|--------|----------|----------|-------------------|")
    for key in ['gasket_L9', 'gasket_L10', 'carpet_L5', 'carpet_L6']:
        r = all_results[key]
        interior = r.get('n_interior', r.get('n_deg4', '?'))
        lines.append(f"| {key} | {r['n_nodes']:,} | {r['diameter']} | {interior:,} |")
    lines.append("")

    # Gasket table
    lines.append("## Gasket PR-vs-Sigma (L=9, L=10)")
    lines.append("")
    lines.append("| sigma | L=9 PR | L=9 std | L=10 PR | L=10 std | Δ(L10−L9) |")
    lines.append("|-------|--------|---------|---------|----------|-----------|")
    for sigma in SIGMA_SWEEP:
        r9 = sw9.get(sigma, {})
        r10 = sw10.get(sigma, {})
        p9 = f"{r9['mean_pr']:.4f}" if not r9.get('degenerate', True) else "DEGEN"
        s9 = f"{r9['std_pr']:.4f}" if not r9.get('degenerate', True) else "—"
        p10 = f"{r10['mean_pr']:.4f}" if not r10.get('degenerate', True) else "DEGEN"
        s10 = f"{r10['std_pr']:.4f}" if not r10.get('degenerate', True) else "—"
        if not r9.get('degenerate', True) and not r10.get('degenerate', True):
            delta = f"{r10['mean_pr'] - r9['mean_pr']:+.4f}"
        else:
            delta = "—"
        lines.append(f"| {sigma} | {p9} | {s9} | {p10} | {s10} | {delta} |")
    lines.append("")
    lines.append(f"Reference: d_H = {GASKET_D_H:.4f}, d_S = {GASKET_D_S:.4f}")
    lines.append("")

    # Carpet table
    lines.append("## Carpet PR-vs-Sigma (L=5, L=6)")
    lines.append("")

    if 'growth' in all_results.get('carpet_L6', {}):
        g = all_results['carpet_L6']['growth']
        lines.append(f"### Carpet L=6: {all_results['carpet_L6']['n_nodes']:,} vertices, diameter {all_results['carpet_L6']['diameter']}")
        lines.append(f"Growth dimension: {g['fit_dimension']:.4f} (R²={g['r_squared']:.4f}), Gate: {'PASS' if all_results['carpet_L6']['growth_pass'] else 'FAIL'}")
        lines.append("")

    lines.append("| sigma | L=5 PR | L=5 std | L=6 PR | L=6 std | Δ(L6−L5) |")
    lines.append("|-------|--------|---------|--------|---------|----------|")
    for sigma in SIGMA_SWEEP:
        r5 = sw5.get(sigma, {})
        r6 = sw6.get(sigma, {})
        p5 = f"{r5['mean_pr']:.4f}" if not r5.get('degenerate', True) else "DEGEN"
        s5 = f"{r5['std_pr']:.4f}" if not r5.get('degenerate', True) else "—"
        p6 = f"{r6['mean_pr']:.4f}" if not r6.get('degenerate', True) else "DEGEN"
        s6 = f"{r6['std_pr']:.4f}" if not r6.get('degenerate', True) else "—"
        if not r5.get('degenerate', True) and not r6.get('degenerate', True):
            delta = f"{r6['mean_pr'] - r5['mean_pr']:+.4f}"
        else:
            delta = "—"
        lines.append(f"| {sigma} | {p5} | {s5} | {p6} | {s6} | {delta} |")
    lines.append("")
    lines.append(f"Reference: d_H = {CARPET_D_H:.4f}, d_S = {CARPET_D_S:.4f}")
    lines.append("")

    # Extrapolation table
    lines.append("## Extrapolation: PR(σ) = d_∞ + A·σ^(-α)")
    lines.append("")
    lines.append("| System | Fit Range | d_∞ | ± err | α | R² | Crosses d_H at σ ≈ |")
    lines.append("|--------|-----------|-----|-------|---|-----|---------------------|")
    for sys_key in ['gasket_L10', 'gasket_L9', 'carpet_L6', 'carpet_L5']:
        if sys_key in all_extrap:
            for fit_range, fit_data in all_extrap[sys_key]['fits'].items():
                if fit_data:
                    d_H_ref = GASKET_D_H if 'gasket' in sys_key else CARPET_D_H
                    if fit_data['d_inf'] < d_H_ref:
                        cross = f"{(fit_data['A'] / (d_H_ref - fit_data['d_inf']))**(1/fit_data['alpha']):.0f}"
                    else:
                        cross = "Never"
                    lines.append(f"| {sys_key} | {fit_range} | {fit_data['d_inf']:.4f} | "
                                f"{fit_data['d_inf_err']:.4f} | {fit_data['alpha']:.4f} | "
                                f"{fit_data['r2']:.6f} | {cross} |")
    lines.append("")

    report = "\n".join(lines)
    path = os.path.join(OUTPUT_DIR, "PHASE1_EXTENSION_V4_RESULTS.md")
    with open(path, 'w') as f:
        f.write(report)
    print(f"\nReport saved to {path}")

    print("\n" + "=" * 70)
    print(f"  PHASE 1 EXTENSION v4 COMPLETE — {total_time:.1f}s")
    print("=" * 70)
    print(report)


if __name__ == '__main__':
    main()
