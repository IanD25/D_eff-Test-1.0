#!/usr/bin/env python3
"""
DS Phase 1 Extension v3: Gasket L=9 + Carpet L=5 with sigma to 100
====================================================================
Push to larger systems to resolve whether Fisher PR converges to d_H.
Also includes curve-fit extrapolation analysis.

Gasket L=9: 29,526 vertices, diameter ~512. sigma=100 is ~20% of diameter.
Carpet L=5: 32,768 vertices, diameter ~240. sigma=100 is ~42% of diameter.
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
CARPET_D_S = 1.805                          # numerical estimate

# Extended sigma sweep — push to 100
SIGMA_SWEEP = [1.5, 2.0, 3.0, 5.0, 8.0, 12.0, 16.0, 20.0, 25.0, 30.0, 40.0, 50.0, 65.0, 80.0, 100.0]

# ============================================================================
# Shared Utilities
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

    frob_norm = np.linalg.norm(F, 'fro')
    if frob_norm < 1e-12:
        return None

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
            print(f"    sigma={sigma:6.1f}: PR={results[sigma]['mean_pr']:.4f} "
                  f"± {results[sigma]['std_pr']:.4f}, rank={results[sigma]['mean_rank']:.2f}{deg_str}")

    return results


# ============================================================================
# Growth & Spectral Dimension
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
        try:
            eigenvalues = eigsh(L, k=k, which='SM', return_eigenvectors=False,
                               maxiter=5000, tol=1e-6)
        except Exception as e:
            print(f"    eigsh failed entirely: {e}")
            return None
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
# Curve Fitting / Extrapolation
# ============================================================================

def power_law_model(sigma, d_inf, A, alpha):
    return d_inf + A * sigma**(-alpha)


def fit_and_extrapolate(sigmas, prs, label, d_H, d_S):
    """Fit PR(sigma) = d_inf + A*sigma^(-alpha) and extrapolate."""
    print(f"\n  --- Extrapolation: {label} ---")
    print(f"  d_H = {d_H:.4f}, d_S = {d_S:.4f}")

    results = {}
    for start_idx, fit_label in [(0, 'all'), (2, 'sigma>=3'), (4, 'sigma>=8')]:
        s = np.array(sigmas[start_idx:])
        p = np.array(prs[start_idx:])
        try:
            popt, pcov = curve_fit(power_law_model, s, p,
                                   p0=[d_H, 2.0, 0.5], maxfev=10000,
                                   bounds=([0.1, 0, 0.01], [5.0, 30.0, 5.0]))
            perr = np.sqrt(np.diag(pcov))
            d_inf, A, alpha = popt
            d_inf_err = perr[0]

            residuals = p - power_law_model(s, *popt)
            ss_res = np.sum(residuals**2)
            ss_tot = np.sum((p - np.mean(p))**2)
            r2 = 1 - ss_res/ss_tot

            pred = {sig: power_law_model(sig, *popt) for sig in [200, 500, 1000, 5000]}

            print(f"    Fit ({fit_label}): d_inf = {d_inf:.4f} ± {d_inf_err:.4f}, "
                  f"A = {A:.4f}, alpha = {alpha:.4f}, R² = {r2:.6f}")
            print(f"      Predictions: σ=200→{pred[200]:.4f}, σ=500→{pred[500]:.4f}, "
                  f"σ=1000→{pred[1000]:.4f}, σ=5000→{pred[5000]:.4f}")

            if d_inf < d_H:
                sigma_cross = (A / (d_H - d_inf))**(1/alpha)
                print(f"      PR crosses d_H at σ ≈ {sigma_cross:.1f}")
            else:
                excess = d_inf - d_H
                print(f"      d_inf ABOVE d_H by {excess:.4f} — PR never reaches d_H in this model")

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

def plot_pr_vs_sigma(sweep_data, title, filename, ref_lines):
    fig, ax = plt.subplots(figsize=(14, 8))
    colors = ['#2196F3', '#FF5722', '#4CAF50', '#9C27B0', '#FF9800']
    markers = ['o', 's', '^', 'D', 'v']

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
    ax.legend(fontsize=9, loc='best')
    ax.grid(True, alpha=0.3)
    path = os.path.join(OUTPUT_DIR, filename)
    fig.savefig(path, dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f"    Saved {path}")


def plot_extrapolation(sigmas_dict, prs_dict, fits_dict, ref_dims, filename, title):
    """Plot actual data + extrapolation curves."""
    fig, ax = plt.subplots(figsize=(14, 8))
    colors = {'gasket': '#2196F3', 'carpet': '#FF5722'}

    for key in sigmas_dict:
        s = sigmas_dict[key]
        p = prs_dict[key]
        ax.scatter(s, p, color=colors.get(key, 'gray'), s=50, zorder=5,
                   label=f'{key} (measured)')

        # Plot best extrapolation curve
        fit = fits_dict.get(key, {})
        best_fit = fit.get('sigma>=3') or fit.get('all')
        if best_fit:
            sigma_ext = np.logspace(np.log10(1.5), np.log10(5000), 200)
            pr_ext = power_law_model(sigma_ext, best_fit['d_inf'],
                                     best_fit['A'], best_fit['alpha'])
            ax.plot(sigma_ext, pr_ext, '--', color=colors.get(key, 'gray'),
                    linewidth=1.5, alpha=0.7,
                    label=f'{key} fit → d∞={best_fit["d_inf"]:.3f}±{best_fit["d_inf_err"]:.3f}')

    # Reference lines
    for label, val in ref_dims.items():
        style = '-' if 'd_H' in label else '--'
        color = 'red' if 'Gasket' in label else 'blue'
        ax.axhline(y=val, color=color, linestyle=style, linewidth=1.5, alpha=0.7,
                   label=f'{label} = {val:.4f}')

    ax.set_xscale('log')
    ax.set_xlabel('σ (log scale)', fontsize=13)
    ax.set_ylabel('Participation Ratio', fontsize=13)
    ax.set_title(title, fontsize=14)
    ax.legend(fontsize=9, loc='upper right')
    ax.grid(True, alpha=0.3)
    ax.set_xlim(1, 6000)
    path = os.path.join(OUTPUT_DIR, filename)
    fig.savefig(path, dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f"    Saved {path}")


def plot_sv_panel(sweep, prefix, output_dir, title=None):
    sigmas = sorted([s for s in sweep if not sweep[s].get('degenerate', False)])
    n = len(sigmas)
    cols = min(5, n)
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

    fig.suptitle(title or f'SV Profiles — {prefix}', fontsize=13)
    fig.tight_layout()
    path = os.path.join(output_dir, f'{prefix}_sv_profiles_full_sweep.png')
    fig.savefig(path, dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f"    Saved {path}")


def plot_growth_delta(deltas, ref_dims, rMin, rMax, prefix, output_dir):
    fig, ax = plt.subplots(figsize=(10, 6))
    rs = [r for r, d in deltas]
    ds = [d for r, d in deltas]
    ax.plot(rs, ds, 'b.-', alpha=0.7, markersize=4)
    colors_ref = ['r', 'g', 'purple']
    for i, (label, val) in enumerate(ref_dims.items()):
        ax.axhline(y=val, color=colors_ref[i % len(colors_ref)], linestyle='--',
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
    ax.plot(np.log(valid.astype(float)), np.log(mean_volumes[1:].astype(float)),
            'b.', alpha=0.5, markersize=3)
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


def plot_convergence_fractal(estimates, ref_dims, prefix, output_dir):
    fig, ax = plt.subplots(figsize=(9, 6))
    labels = list(estimates.keys())
    vals = list(estimates.values())
    ax.bar(labels, vals, color='#4CAF50', edgecolor='black', linewidth=0.8, width=0.5)
    colors_ref = ['black', 'red', 'purple']
    for i, (label, val) in enumerate(ref_dims.items()):
        ax.axhline(y=val, color=colors_ref[i % len(colors_ref)], linestyle='--',
                   linewidth=2, label=f'{label} = {val:.4f}')
    for i, val in enumerate(vals):
        ax.text(i, val + 0.02, f'{val:.3f}', ha='center', va='bottom',
                fontsize=11, fontweight='bold')
    ax.set_ylabel('Dimension Estimate', fontsize=12)
    ax.set_title(f'Convergence — {prefix}', fontsize=14)
    ax.legend(fontsize=10)
    ax.grid(True, alpha=0.3, axis='y')
    path = os.path.join(output_dir, f'{prefix}_convergence.png')
    fig.savefig(path, dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f"    Saved {path}")


# ============================================================================
# Main: Gasket L=8 + L=9
# ============================================================================

def run_gasket():
    print(f"\n{'='*70}")
    print("  GASKET: L=8 and L=9 — Sigma Sweep to 100")
    print(f"{'='*70}")

    results = {}

    for level in [8, 9]:
        print(f"\n--- Sierpinski Gasket L={level} ---")
        t0 = time.time()
        G, corners = make_sierpinski_gasket(level)
        n = G.number_of_nodes()
        expected = (3 ** (level + 1) + 3) // 2
        assert n == expected, f"Expected {expected}, got {n}"
        degrees = dict(G.degree())
        interior = [v for v, d in degrees.items() if d == 4]
        build_time = time.time() - t0
        print(f"  Built: {n} vertices, {len(interior)} interior, {build_time:.1f}s")

        # Compute diameter estimate via BFS from a corner
        d_corner = bfs_distances(G, corners[0])
        diameter = max(d_corner.values())
        print(f"  Diameter: {diameter}")

        sample_interior = [interior[i] for i in
                           np.random.choice(len(interior), size=min(20, len(interior)), replace=False)]

        print(f"\n  [Sigma sweep — {len(SIGMA_SWEEP)} values]")
        t0 = time.time()
        sweep = fisher_sigma_sweep(G, SIGMA_SWEEP, sample_interior, n_samples=20)
        print(f"  Sweep time: {time.time()-t0:.1f}s")
        results[f'L{level}'] = {'sweep': sweep, 'level': level, 'n_nodes': n,
                                'diameter': diameter}

        plot_sv_panel(sweep, f'sierpinski_L{level}', OUTPUT_DIR,
                      f'SV Profiles — Gasket L={level} (σ to 100)')

    # PR-vs-sigma plot
    sweep_data = {}
    for key in ['L8', 'L9']:
        sweep_data[f'Gasket {key}'] = results[key]['sweep']

    plot_pr_vs_sigma(
        sweep_data,
        'PR vs σ — Sierpinski Gasket L=8 vs L=9 (Extended to σ=100)',
        'sierpinski_pr_vs_sigma_v3.png',
        {'Gasket d_H': GASKET_D_H, 'Gasket d_S': GASKET_D_S}
    )

    return results


# ============================================================================
# Main: Carpet L=4 + L=5
# ============================================================================

def run_carpet():
    print(f"\n{'='*70}")
    print("  CARPET: L=4 and L=5 — Sigma Sweep to 100")
    print(f"{'='*70}")

    results = {}

    for level in [4, 5]:
        prefix = f'carpet_L{level}'
        print(f"\n--- Sierpinski Carpet L={level} ---")
        t0 = time.time()
        G = make_sierpinski_carpet(level)
        n = G.number_of_nodes()
        expected = 8 ** level
        assert n == expected, f"Expected {expected}, got {n}"
        assert nx.is_connected(G), "Graph should be connected"
        degrees = dict(G.degree())
        deg_counts = {}
        for v, d in degrees.items():
            deg_counts[d] = deg_counts.get(d, 0) + 1
        deg4 = [v for v, d in degrees.items() if d == 4]
        build_time = time.time() - t0
        print(f"  Built: {n} vertices, {len(deg4)} interior (deg 4), {build_time:.1f}s")
        print(f"  Degree distribution: {dict(sorted(deg_counts.items()))}")

        # Diameter estimate
        nodes_list = list(G.nodes())
        d0 = bfs_distances(G, nodes_list[0])
        far_node = max(d0, key=d0.get)
        d1 = bfs_distances(G, far_node)
        diameter = max(d1.values())
        print(f"  Diameter: {diameter}")

        results[f'L{level}'] = {
            'level': level, 'n_nodes': n, 'deg_counts': deg_counts,
            'n_deg4': len(deg4), 'diameter': diameter
        }

        sample_deg4 = [deg4[i] for i in
                       np.random.choice(len(deg4), size=min(30, len(deg4)), replace=False)]

        if level == 5:
            # Growth dimension
            print(f"\n  [Route 1] Growth Dimension...")
            t0 = time.time()
            mean_vol, maxR = compute_ball_volumes(G, n_samples=30, sample_nodes=sample_deg4)
            growth = estimate_growth_dimension(mean_vol, maxR)
            gate_lo, gate_hi = 1.70, 2.10
            gpass = gate_lo <= growth['fit_dimension'] <= gate_hi
            print(f"  Growth dim: {growth['fit_dimension']:.4f} (R²={growth['r_squared']:.4f})")
            print(f"  Gate [{gate_lo}, {gate_hi}]: {'PASS' if gpass else 'FAIL'}")
            print(f"  Time: {time.time()-t0:.1f}s")
            results['L5']['growth'] = growth
            results['L5']['growth_pass'] = gpass

            ref_dims = {'d_H': CARPET_D_H, 'd_S': CARPET_D_S}
            plot_growth_delta(growth['deltas'], ref_dims, growth['rMin'],
                              growth['rMax'], prefix, OUTPUT_DIR)
            plot_growth_loglog(mean_vol, growth['rMin'], growth['rMax'],
                               growth['fit_dimension'], CARPET_D_H, prefix, OUTPUT_DIR)

            # Spectral dimension
            print(f"\n  [Route 2] Spectral Dimension...")
            t0 = time.time()
            spectral = compute_spectral_dimension_eigsh(G, n_eigenvalues=500)
            if spectral is not None:
                gate_lo_s, gate_hi_s = 1.60, 2.00
                spass = gate_lo_s <= spectral['spectral_dimension'] <= gate_hi_s
                print(f"  Spectral dim: {spectral['spectral_dimension']:.4f} "
                      f"(R²={spectral['r_squared']:.4f})")
                print(f"  Gate [{gate_lo_s}, {gate_hi_s}]: {'PASS' if spass else 'FAIL'}")
                print(f"  Time: {time.time()-t0:.1f}s")
                results['L5']['spectral'] = spectral
                results['L5']['spectral_pass'] = spass

                plot_spectral_weyl(spectral['eigenvalues'], spectral['fit_range'][0],
                                   spectral['fit_range'][1], spectral['weyl_slope'],
                                   spectral['spectral_dimension'], CARPET_D_S, prefix, OUTPUT_DIR)
            else:
                print(f"  Spectral: SKIPPED (eigsh failed)")
                results['L5']['spectral'] = None
                results['L5']['spectral_pass'] = False

        # Fisher sigma sweep (both levels)
        print(f"\n  [Fisher sigma sweep L={level}]")
        t0 = time.time()
        sweep = fisher_sigma_sweep(G, SIGMA_SWEEP, sample_deg4, n_samples=20)
        print(f"  Sweep time: {time.time()-t0:.1f}s")
        results[f'L{level}']['sweep'] = sweep

        plot_sv_panel(sweep, prefix, OUTPUT_DIR, f'SV Profiles — Carpet L={level} (σ to 100)')

    # PR-vs-sigma plot
    sweep_data = {}
    for key in ['L4', 'L5']:
        sweep_data[f'Carpet {key}'] = results[key]['sweep']

    plot_pr_vs_sigma(
        sweep_data,
        'PR vs σ — Sierpinski Carpet L=4 vs L=5 (Extended to σ=100)',
        'carpet_pr_vs_sigma_v3.png',
        {'Carpet d_H': CARPET_D_H, 'Carpet d_S': CARPET_D_S}
    )

    # L=5 convergence plot
    if 'growth' in results.get('L5', {}):
        g = results['L5']['growth']
        s = results['L5'].get('spectral', {})
        s_dim = s['spectral_dimension'] if s else 0
        sweep5 = results['L5']['sweep']
        pr3 = sweep5[3.0]['mean_pr'] if not sweep5.get(3.0, {}).get('degenerate', True) else 0
        plot_convergence_fractal(
            {'Growth': g['fit_dimension'],
             'Spectral': s_dim,
             'Fisher PR\n(σ=3)': pr3},
            {'d_H': CARPET_D_H, 'd_S': CARPET_D_S},
            'carpet_L5', OUTPUT_DIR)

    return results


# ============================================================================
# Main
# ============================================================================

def main():
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    t_start = time.time()

    # Run both systems
    gasket_results = run_gasket()
    carpet_results = run_carpet()

    # ----------------------------------------------------------------
    # Extrapolation analysis
    # ----------------------------------------------------------------
    print(f"\n{'='*70}")
    print("  EXTRAPOLATION ANALYSIS")
    print(f"{'='*70}")

    all_extrap = {}

    # Gasket L=9 (most reliable)
    if 'L9' in gasket_results:
        sw9 = gasket_results['L9']['sweep']
        sigmas_g, prs_g = [], []
        for sigma in sorted(sw9.keys()):
            if not sw9[sigma].get('degenerate', False):
                sigmas_g.append(sigma)
                prs_g.append(sw9[sigma]['mean_pr'])
        fit_g = fit_and_extrapolate(sigmas_g, prs_g, 'Gasket L=9', GASKET_D_H, GASKET_D_S)
        all_extrap['gasket'] = fit_g

    # Also fit L=8 for comparison
    if 'L8' in gasket_results:
        sw8 = gasket_results['L8']['sweep']
        sigmas_g8, prs_g8 = [], []
        for sigma in sorted(sw8.keys()):
            if not sw8[sigma].get('degenerate', False):
                sigmas_g8.append(sigma)
                prs_g8.append(sw8[sigma]['mean_pr'])
        fit_g8 = fit_and_extrapolate(sigmas_g8, prs_g8, 'Gasket L=8', GASKET_D_H, GASKET_D_S)
        all_extrap['gasket_L8'] = fit_g8

    # Carpet L=5 (most reliable)
    if 'L5' in carpet_results:
        sw5 = carpet_results['L5']['sweep']
        sigmas_c, prs_c = [], []
        for sigma in sorted(sw5.keys()):
            if not sw5[sigma].get('degenerate', False):
                sigmas_c.append(sigma)
                prs_c.append(sw5[sigma]['mean_pr'])
        fit_c = fit_and_extrapolate(sigmas_c, prs_c, 'Carpet L=5', CARPET_D_H, CARPET_D_S)
        all_extrap['carpet'] = fit_c

    # Carpet L=4 for comparison
    if 'L4' in carpet_results:
        sw4 = carpet_results['L4']['sweep']
        sigmas_c4, prs_c4 = [], []
        for sigma in sorted(sw4.keys()):
            if not sw4[sigma].get('degenerate', False):
                sigmas_c4.append(sigma)
                prs_c4.append(sw4[sigma]['mean_pr'])
        fit_c4 = fit_and_extrapolate(sigmas_c4, prs_c4, 'Carpet L=4', CARPET_D_H, CARPET_D_S)
        all_extrap['carpet_L4'] = fit_c4

    # ----------------------------------------------------------------
    # Money plot: Gasket L=9 vs Carpet L=5 overlay
    # ----------------------------------------------------------------
    print(f"\n  [Money Plot v3: Gasket L=9 vs Carpet L=5]")
    money_sweep = {}
    if 'L9' in gasket_results:
        money_sweep['Gasket (L=9, 29k verts)'] = gasket_results['L9']['sweep']
    if 'L5' in carpet_results:
        money_sweep['Carpet (L=5, 33k verts)'] = carpet_results['L5']['sweep']

    plot_pr_vs_sigma(
        money_sweep,
        'Fisher PR vs σ — Gasket L=9 vs Carpet L=5 (σ to 100)',
        'fractal_pr_vs_sigma_v3.png',
        {
            'Gasket d_H': GASKET_D_H,
            'Gasket d_S': GASKET_D_S,
            'Carpet d_H': CARPET_D_H,
            'Carpet d_S': CARPET_D_S,
        }
    )

    # Extrapolation plot
    if 'L9' in gasket_results and 'L5' in carpet_results:
        plot_extrapolation(
            {'gasket': sigmas_g, 'carpet': sigmas_c},
            {'gasket': prs_g, 'carpet': prs_c},
            all_extrap,
            {
                'Gasket d_H': GASKET_D_H,
                'Gasket d_S': GASKET_D_S,
                'Carpet d_H': CARPET_D_H,
                'Carpet d_S': CARPET_D_S,
            },
            'fractal_pr_extrapolation.png',
            'PR Extrapolation — Power Law Fit to σ=5000'
        )

    # ----------------------------------------------------------------
    # L=8 vs L=9 agreement analysis
    # ----------------------------------------------------------------
    print(f"\n{'='*70}")
    print("  L=8 vs L=9 AGREEMENT (Gasket)")
    print(f"{'='*70}")
    if 'L8' in gasket_results and 'L9' in gasket_results:
        sw8 = gasket_results['L8']['sweep']
        sw9 = gasket_results['L9']['sweep']
        print(f"\n  {'sigma':>6}  {'L8 PR':>8}  {'L9 PR':>8}  {'Δ(L9-L8)':>10}  {'Agreement':>10}")
        for sigma in SIGMA_SWEEP:
            r8 = sw8.get(sigma, {})
            r9 = sw9.get(sigma, {})
            if not r8.get('degenerate', True) and not r9.get('degenerate', True):
                p8, p9 = r8['mean_pr'], r9['mean_pr']
                delta = p9 - p8
                agree = "✓ CLOSE" if abs(delta) < 0.02 else "DIVERGING" if abs(delta) > 0.05 else "~moderate"
                print(f"  {sigma:6.1f}  {p8:8.4f}  {p9:8.4f}  {delta:+10.4f}  {agree:>10}")

    # L=4 vs L=5 agreement (Carpet)
    print(f"\n{'='*70}")
    print("  L=4 vs L=5 AGREEMENT (Carpet)")
    print(f"{'='*70}")
    if 'L4' in carpet_results and 'L5' in carpet_results:
        sw4 = carpet_results['L4']['sweep']
        sw5 = carpet_results['L5']['sweep']
        print(f"\n  {'sigma':>6}  {'L4 PR':>8}  {'L5 PR':>8}  {'Δ(L5-L4)':>10}  {'Agreement':>10}")
        for sigma in SIGMA_SWEEP:
            r4 = sw4.get(sigma, {})
            r5 = sw5.get(sigma, {})
            if not r4.get('degenerate', True) and not r5.get('degenerate', True):
                p4, p5 = r4['mean_pr'], r5['mean_pr']
                delta = p5 - p4
                agree = "✓ CLOSE" if abs(delta) < 0.02 else "DIVERGING" if abs(delta) > 0.05 else "~moderate"
                print(f"  {sigma:6.1f}  {p4:8.4f}  {p5:8.4f}  {delta:+10.4f}  {agree:>10}")

    # ----------------------------------------------------------------
    # Report
    # ----------------------------------------------------------------
    total_time = time.time() - t_start

    lines = []
    lines.append("# DS Phase 1 Extension v3 Results")
    lines.append(f"## Date: {datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    lines.append(f"## Runtime: {total_time:.1f} seconds")
    lines.append("")

    # Gasket results
    lines.append("## Gasket Sigma Sweep (L=8 and L=9, σ to 100)")
    lines.append("")
    lines.append("| sigma | L=8 PR | L=8 std | L=9 PR | L=9 std | Δ(L9−L8) |")
    lines.append("|-------|--------|---------|--------|---------|----------|")
    sw8 = gasket_results.get('L8', {}).get('sweep', {})
    sw9 = gasket_results.get('L9', {}).get('sweep', {})
    for sigma in SIGMA_SWEEP:
        r8 = sw8.get(sigma, {})
        r9 = sw9.get(sigma, {})
        p8 = f"{r8['mean_pr']:.4f}" if not r8.get('degenerate', True) else "DEGEN"
        s8 = f"{r8['std_pr']:.4f}" if not r8.get('degenerate', True) else "—"
        p9 = f"{r9['mean_pr']:.4f}" if not r9.get('degenerate', True) else "DEGEN"
        s9 = f"{r9['std_pr']:.4f}" if not r9.get('degenerate', True) else "—"
        if not r8.get('degenerate', True) and not r9.get('degenerate', True):
            delta = f"{r9['mean_pr'] - r8['mean_pr']:+.4f}"
        else:
            delta = "—"
        lines.append(f"| {sigma} | {p8} | {s8} | {p9} | {s9} | {delta} |")
    lines.append("")
    lines.append(f"Reference: d_H = {GASKET_D_H:.4f}, d_S = {GASKET_D_S:.4f}")
    lines.append(f"L=8 diameter: {gasket_results.get('L8', {}).get('diameter', '?')}")
    lines.append(f"L=9 diameter: {gasket_results.get('L9', {}).get('diameter', '?')}")
    lines.append("")

    # Carpet results
    lines.append("## Carpet Results (L=4 and L=5, σ to 100)")
    lines.append("")

    if 'L5' in carpet_results:
        lines.append(f"### Carpet L=5 ({carpet_results['L5']['n_nodes']} vertices)")
        lines.append(f"Diameter: {carpet_results['L5']['diameter']}")
        lines.append(f"Degree distribution: {dict(sorted(carpet_results['L5']['deg_counts'].items()))}")
        lines.append("")

        if 'growth' in carpet_results['L5']:
            g = carpet_results['L5']['growth']
            lines.append(f"Growth dimension: {g['fit_dimension']:.4f} (R²={g['r_squared']:.4f})")
            lines.append(f"  Gate [1.70, 2.10]: {'PASS' if carpet_results['L5']['growth_pass'] else 'FAIL'}")

        if carpet_results['L5'].get('spectral'):
            s = carpet_results['L5']['spectral']
            lines.append(f"Spectral dimension: {s['spectral_dimension']:.4f} (R²={s['r_squared']:.4f})")
            lines.append(f"  Gate [1.60, 2.00]: {'PASS' if carpet_results['L5']['spectral_pass'] else 'FAIL'}")
        lines.append("")

    lines.append("| sigma | L=4 PR | L=4 std | L=5 PR | L=5 std | Δ(L5−L4) |")
    lines.append("|-------|--------|---------|--------|---------|----------|")
    sw4 = carpet_results.get('L4', {}).get('sweep', {})
    sw5 = carpet_results.get('L5', {}).get('sweep', {})
    for sigma in SIGMA_SWEEP:
        r4 = sw4.get(sigma, {})
        r5 = sw5.get(sigma, {})
        p4 = f"{r4['mean_pr']:.4f}" if not r4.get('degenerate', True) else "DEGEN"
        s4 = f"{r4['std_pr']:.4f}" if not r4.get('degenerate', True) else "—"
        p5 = f"{r5['mean_pr']:.4f}" if not r5.get('degenerate', True) else "DEGEN"
        s5 = f"{r5['std_pr']:.4f}" if not r5.get('degenerate', True) else "—"
        if not r4.get('degenerate', True) and not r5.get('degenerate', True):
            delta = f"{r5['mean_pr'] - r4['mean_pr']:+.4f}"
        else:
            delta = "—"
        lines.append(f"| {sigma} | {p4} | {s4} | {p5} | {s5} | {delta} |")
    lines.append("")
    lines.append(f"Reference: d_H = {CARPET_D_H:.4f}, d_S = {CARPET_D_S:.4f}")
    lines.append("")

    # Extrapolation summary
    lines.append("## Extrapolation: PR(σ) = d_∞ + A·σ^(-α)")
    lines.append("")
    lines.append("| System | Fit Range | d_∞ | ± err | A | α | R² |")
    lines.append("|--------|-----------|-----|-------|---|---|-----|")
    for sys_key, sys_label in [('gasket', 'Gasket L=9'), ('gasket_L8', 'Gasket L=8'),
                                ('carpet', 'Carpet L=5'), ('carpet_L4', 'Carpet L=4')]:
        if sys_key in all_extrap:
            for fit_range, fit_data in all_extrap[sys_key].items():
                if fit_data:
                    lines.append(f"| {sys_label} | {fit_range} | {fit_data['d_inf']:.4f} | "
                                f"{fit_data['d_inf_err']:.4f} | {fit_data['A']:.4f} | "
                                f"{fit_data['alpha']:.4f} | {fit_data['r2']:.6f} |")
    lines.append("")

    report = "\n".join(lines)
    path = os.path.join(OUTPUT_DIR, "PHASE1_EXTENSION_V3_RESULTS.md")
    with open(path, 'w') as f:
        f.write(report)
    print(f"\nReport saved to {path}")

    print("\n" + "=" * 70)
    print("  PHASE 1 EXTENSION v3 COMPLETE")
    print(f"  Total time: {total_time:.1f}s")
    print("=" * 70)
    print(report)


if __name__ == '__main__':
    main()
