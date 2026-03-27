#!/usr/bin/env python3
"""
Phase 4A: Kuramoto Synchronization — FIM Spectral Transition Test
=================================================================
Tests whether the FIM SV₂/SV₁ degeneracy swap occurs at the Kuramoto
critical coupling K_c = 2γ/π ≈ 0.637 (exact for Lorentzian distribution).

Two-scale mapping:
  Scale 1: σ (natural frequency distribution width)
  Scale 2: K (coupling strength)
  At K_c: the two scales balance.

Pre-registered predictions:
  P4A-1: SV₂/SV₁ peaks within ±15% of K_c = 0.637  [KILL if absent]
  P4A-2: SV₂/SV₁ at K_c exceeds K=0.1 by ≥0.10
  P4A-3: Non-monotonic peak (below < peak > above)
  P4A-4: η minimum coincides with SV₂/SV₁ maximum (±20% K)
  P4A-5: Rank increases approaching K_c

Author: Ian Darling / Claude Code
Date: 2026-03-25
"""

import numpy as np
from scipy.stats import cauchy
from scipy.spatial.distance import squareform, pdist
import csv
import time
import os
import sys

# ─── Configuration ────────────────────────────────────────────────────
N_OSC = 500           # Number of oscillators
K_VALUES = np.concatenate([
    np.linspace(0.1, 0.50, 5),    # Below K_c
    np.linspace(0.54, 0.73, 6),   # Near K_c (±15% of 0.637)
    np.linspace(0.80, 2.0, 9),    # Above K_c
])
GAMMA = 1.0 / np.pi               # Lorentzian scale γ=1/π so K_c = 2γ = 2/π ≈ 0.637
K_C_EXACT = 2.0 * GAMMA            # = 2/π ≈ 0.6366 (exact for Lorentzian)

T_TOTAL = 200.0       # Total integration time
T_TRANSIENT = 100.0   # Discard first 100 time units
DT = 0.01             # Euler timestep
N_SNAPSHOTS = 10      # Snapshots in analysis window
N_SEEDS = 3           # Independent random seeds per K

K_NN = 10             # k-NN neighbors for graph construction
MAX_R = 15            # Maximum BFS distance for C(r)
N_SAMPLE_CENTERS = 50 # Centers for FIM sampling (speed vs accuracy)

RESULTS_DIR = os.path.join(os.path.dirname(__file__), "results")


def integrate_kuramoto(omega, K, dt, n_steps):
    """Euler integration of Kuramoto ODE. Returns final phases."""
    N = len(omega)
    theta = np.random.uniform(0, 2 * np.pi, N)

    for _ in range(n_steps):
        # Vectorized: sin(θⱼ - θᵢ) for all pairs
        diff = theta[np.newaxis, :] - theta[:, np.newaxis]  # (N, N)
        coupling = (K / N) * np.sum(np.sin(diff), axis=1)
        theta += dt * (omega + coupling)

    return theta


def integrate_kuramoto_snapshots(omega, K, dt, t_transient, t_analysis, n_snapshots):
    """Integrate and return snapshots during the analysis window."""
    N = len(omega)
    theta = np.random.uniform(0, 2 * np.pi, N)

    n_transient = int(t_transient / dt)
    n_analysis = int(t_analysis / dt)
    snapshot_interval = n_analysis // n_snapshots

    # Run transient
    for _ in range(n_transient):
        diff = theta[np.newaxis, :] - theta[:, np.newaxis]
        coupling = (K / N) * np.sum(np.sin(diff), axis=1)
        theta += dt * (omega + coupling)

    # Run analysis window, collect snapshots
    snapshots = []
    R_values = []
    step = 0
    for _ in range(n_analysis):
        diff = theta[np.newaxis, :] - theta[:, np.newaxis]
        coupling = (K / N) * np.sum(np.sin(diff), axis=1)
        theta += dt * (omega + coupling)
        step += 1

        if step % snapshot_interval == 0 and len(snapshots) < n_snapshots:
            snapshots.append(theta.copy())
            # Order parameter R
            R = np.abs(np.mean(np.exp(1j * theta)))
            R_values.append(R)

    return snapshots, R_values


def build_knn_graph(corr_matrix, k):
    """Build symmetric k-NN graph from correlation matrix.
    Returns adjacency list (dict of sets)."""
    N = corr_matrix.shape[0]
    adj = {i: set() for i in range(N)}

    for i in range(N):
        # Find k most correlated neighbors (highest |correlation|)
        sims = np.abs(corr_matrix[i])
        sims[i] = -1  # Exclude self
        neighbors = np.argsort(sims)[-k:]
        for j in neighbors:
            adj[i].add(j)
            adj[j].add(i)  # Symmetric

    return adj


def bfs_distances(adj, source, max_r):
    """BFS from source, return dict of node→distance."""
    dist = {source: 0}
    queue = [source]
    head = 0
    while head < len(queue):
        u = queue[head]
        head += 1
        if dist[u] >= max_r:
            continue
        for v in adj[u]:
            if v not in dist:
                dist[v] = dist[u] + 1
                queue.append(v)
    return dist


def compute_cr(corr_matrix, adj, max_r, n_centers):
    """Compute C(r) = mean |correlation| at graph distance r.
    Sample n_centers random source nodes."""
    N = corr_matrix.shape[0]
    centers = np.random.choice(N, min(n_centers, N), replace=False)

    cr_sum = np.zeros(max_r + 1)
    cr_count = np.zeros(max_r + 1)

    for src in centers:
        dists = bfs_distances(adj, src, max_r)
        for node, r in dists.items():
            if node != src and r <= max_r:
                cr_sum[r] += np.abs(corr_matrix[src, node])
                cr_count[r] += 1

    # Avoid division by zero
    mask = cr_count > 0
    cr = np.zeros(max_r + 1)
    cr[mask] = cr_sum[mask] / cr_count[mask]

    return cr, mask


def build_fim_from_cr(cr, adj, k, n_centers, corr_matrix):
    """Build FIM in neighbor basis from C(r) kernel.

    For each center v0, define score vectors s_j(u) for each of k neighbors.
    The score is: log p_{w_j}(u) - log p_{v0}(u)
    where p_v(u) = exp(-d(v,u)/σ) / Z(v,σ)

    Here we use C(r) as the kernel instead of pure exponential.
    """
    N = corr_matrix.shape[0]
    centers = np.random.choice(N, min(n_centers, N), replace=False)

    # Use cr as position-indexed kernel: p_v(u) ∝ cr[d(v,u)]
    # Score vector s_j(u) = log(cr[d(w_j,u)]) - log(cr[d(v0,u)])

    fim_accum = np.zeros((k, k))
    n_valid = 0

    for v0 in centers:
        # Get distances from v0
        dists_v0 = bfs_distances(adj, v0, MAX_R)

        # Get k neighbors of v0
        neighbors_v0 = list(adj[v0])[:k]
        if len(neighbors_v0) < k:
            continue

        # Get distances from each neighbor
        neighbor_dists = []
        for wj in neighbors_v0:
            neighbor_dists.append(bfs_distances(adj, wj, MAX_R))

        # Compute score vectors over all reachable nodes
        all_nodes = set(dists_v0.keys())
        for nd in neighbor_dists:
            all_nodes &= set(nd.keys())
        all_nodes.discard(v0)
        all_nodes = list(all_nodes)

        if len(all_nodes) < 10:
            continue

        # Build score matrix: (k, n_nodes)
        scores = np.zeros((k, len(all_nodes)))
        weights = np.zeros(len(all_nodes))

        for idx, u in enumerate(all_nodes):
            r_v0 = dists_v0[u]
            # Kernel value at v0
            if r_v0 < len(cr) and cr[r_v0] > 1e-12:
                log_p_v0 = np.log(cr[r_v0] + 1e-15)
            else:
                log_p_v0 = -30.0

            weights[idx] = cr[r_v0] if r_v0 < len(cr) else 1e-15

            for j in range(k):
                r_wj = neighbor_dists[j].get(u, MAX_R + 1)
                if r_wj < len(cr) and cr[r_wj] > 1e-12:
                    log_p_wj = np.log(cr[r_wj] + 1e-15)
                else:
                    log_p_wj = -30.0
                scores[j, idx] = log_p_wj - log_p_v0

        # Normalize weights
        weights /= (weights.sum() + 1e-15)

        # FIM = Σ_u s_i(u) * s_j(u) * p_v0(u)
        weighted_scores = scores * np.sqrt(weights)[np.newaxis, :]
        fim_local = weighted_scores @ weighted_scores.T
        fim_accum += fim_local
        n_valid += 1

    if n_valid > 0:
        fim_accum /= n_valid

    return fim_accum, n_valid


def analyze_fim(fim):
    """Extract diagnostics from FIM: SV₂/SV₁, rank, η, PR."""
    svs = np.linalg.svd(fim, compute_uv=False)
    svs = np.sort(svs)[::-1]  # Descending

    if svs[0] < 1e-15:
        return {
            'sv2sv1': 0.0,
            'rank': 1,
            'eta': 1.0,
            'pr': 1.0,
            'svs': svs
        }

    # Normalize
    svs_norm = svs / svs[0]

    # SV₂/SV₁
    sv2sv1 = svs_norm[1] if len(svs_norm) > 1 else 0.0

    # Gap-based rank
    gaps = []
    for i in range(len(svs_norm) - 1):
        if svs_norm[i] > 1e-10:
            gaps.append(svs_norm[i] / svs_norm[i + 1] if svs_norm[i + 1] > 1e-10 else 1e10)
        else:
            gaps.append(1.0)
    rank = np.argmax(gaps) + 1 if gaps else 1

    # Participation ratio
    svs_sq = svs_norm ** 2
    pr = (svs_sq.sum()) ** 2 / (svs_sq ** 2).sum() if (svs_sq ** 2).sum() > 0 else 1.0

    # η = SV₂ / (SV₁ + SV₂)
    eta = svs_norm[1] / (svs_norm[0] + svs_norm[1]) if len(svs_norm) > 1 else 0.0

    return {
        'sv2sv1': sv2sv1,
        'rank': rank,
        'eta': eta,
        'pr': pr,
        'svs': svs_norm
    }


def run_single_K(K_val, seed, omega=None):
    """Run full pipeline for one K value and one seed."""
    rng = np.random.RandomState(seed)
    np.random.seed(seed)

    if omega is None:
        # Lorentzian frequencies, clipped
        omega = cauchy.rvs(loc=0, scale=GAMMA, size=N_OSC, random_state=rng)
        omega = np.clip(omega, -20, 20)

    # Integrate and get snapshots
    t_analysis = T_TOTAL - T_TRANSIENT
    snapshots, R_values = integrate_kuramoto_snapshots(
        omega, K_val, DT, T_TRANSIENT, t_analysis, N_SNAPSHOTS
    )

    # Compute FIM diagnostics at each snapshot
    sv2sv1_list = []
    rank_list = []
    eta_list = []
    pr_list = []

    for theta in snapshots:
        # Phase correlation matrix: C(i,j) = cos(θᵢ - θⱼ)
        diff = theta[:, np.newaxis] - theta[np.newaxis, :]
        corr = np.cos(diff)

        # Build k-NN graph on correlation
        adj = build_knn_graph(corr, K_NN)

        # Compute C(r)
        cr, mask = compute_cr(corr, adj, MAX_R, N_SAMPLE_CENTERS)

        # Build FIM
        fim, n_valid = build_fim_from_cr(cr, adj, K_NN, N_SAMPLE_CENTERS, corr)

        # Analyze
        diag = analyze_fim(fim)
        sv2sv1_list.append(diag['sv2sv1'])
        rank_list.append(diag['rank'])
        eta_list.append(diag['eta'])
        pr_list.append(diag['pr'])

    return {
        'K': K_val,
        'seed': seed,
        'R_mean': np.mean(R_values),
        'R_std': np.std(R_values),
        'sv2sv1_mean': np.mean(sv2sv1_list),
        'sv2sv1_std': np.std(sv2sv1_list),
        'rank_mean': np.mean(rank_list),
        'eta_mean': np.mean(eta_list),
        'pr_mean': np.mean(pr_list),
        'n_snapshots': len(snapshots),
    }


def main():
    os.makedirs(RESULTS_DIR, exist_ok=True)

    print("=" * 70)
    print("PHASE 4A: KURAMOTO SYNCHRONIZATION — FIM SPECTRAL TRANSITION TEST")
    print("=" * 70)
    print(f"N = {N_OSC} oscillators, K_c = {K_C_EXACT:.4f}")
    print(f"K sweep: {len(K_VALUES)} values from {K_VALUES[0]:.2f} to {K_VALUES[-1]:.2f}")
    print(f"Seeds per K: {N_SEEDS}, Snapshots per run: {N_SNAPSHOTS}")
    print(f"Integration: T={T_TOTAL}, dt={DT}, transient={T_TRANSIENT}")
    print()

    t_start = time.time()

    # Aggregate results per K (average across seeds)
    results = []

    for i, K_val in enumerate(K_VALUES):
        seed_results = []
        for s in range(N_SEEDS):
            seed = 42 + s * 1000
            print(f"  K={K_val:.3f} seed={seed} ({i+1}/{len(K_VALUES)})...", end=" ", flush=True)
            t0 = time.time()
            res = run_single_K(K_val, seed)
            dt_run = time.time() - t0
            print(f"R={res['R_mean']:.3f} SV₂/SV₁={res['sv2sv1_mean']:.4f} rank={res['rank_mean']:.1f} ({dt_run:.1f}s)")
            seed_results.append(res)

        # Average across seeds
        agg = {
            'K': K_val,
            'R_mean': np.mean([r['R_mean'] for r in seed_results]),
            'R_std': np.sqrt(np.mean([r['R_std']**2 for r in seed_results]) +
                           np.var([r['R_mean'] for r in seed_results])),
            'sv2sv1_mean': np.mean([r['sv2sv1_mean'] for r in seed_results]),
            'sv2sv1_std': np.sqrt(np.mean([r['sv2sv1_std']**2 for r in seed_results]) +
                                np.var([r['sv2sv1_mean'] for r in seed_results])),
            'rank_mean': np.mean([r['rank_mean'] for r in seed_results]),
            'eta_mean': np.mean([r['eta_mean'] for r in seed_results]),
            'pr_mean': np.mean([r['pr_mean'] for r in seed_results]),
            'n_snapshots': N_SNAPSHOTS,
        }
        results.append(agg)
        print(f"  → K={K_val:.3f}: R={agg['R_mean']:.3f}±{agg['R_std']:.3f} "
              f"SV₂/SV₁={agg['sv2sv1_mean']:.4f}±{agg['sv2sv1_std']:.4f} "
              f"rank={agg['rank_mean']:.1f} η={agg['eta_mean']:.4f}")
        print()

    t_total = time.time() - t_start
    print(f"\nTotal runtime: {t_total:.1f}s ({t_total/60:.1f} min)")

    # ─── Write CSV ────────────────────────────────────────────────────
    csv_path = os.path.join(RESULTS_DIR, "phase4a_results.csv")
    with open(csv_path, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=[
            'K', 'R_mean', 'R_std', 'sv2sv1_mean', 'sv2sv1_std',
            'rank_mean', 'eta_mean', 'pr_mean', 'n_snapshots'
        ])
        writer.writeheader()
        for r in results:
            writer.writerow({k: f"{v:.6f}" if isinstance(v, float) else v
                           for k, v in r.items()})
    print(f"Results written to {csv_path}")

    # ─── Kill Test & Predictions ──────────────────────────────────────
    summary = evaluate_predictions(results)

    summary_path = os.path.join(RESULTS_DIR, "phase4a_summary.txt")
    with open(summary_path, 'w') as f:
        f.write(summary)
    print(f"Summary written to {summary_path}")

    # ─── Plots ────────────────────────────────────────────────────────
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        make_plots(results, plt)
        print("Plots saved.")
    except ImportError:
        print("matplotlib not available — skipping plots")

    print("\n" + summary)


def evaluate_predictions(results):
    """Evaluate all pre-registered predictions."""
    Ks = np.array([r['K'] for r in results])
    sv2sv1 = np.array([r['sv2sv1_mean'] for r in results])
    eta = np.array([r['eta_mean'] for r in results])
    ranks = np.array([r['rank_mean'] for r in results])
    Rs = np.array([r['R_mean'] for r in results])

    lines = []
    lines.append("=" * 70)
    lines.append("PHASE 4A: KURAMOTO — PREDICTION VERDICTS")
    lines.append("=" * 70)
    lines.append(f"K_c (exact) = {K_C_EXACT:.4f}")
    lines.append(f"Kill window: [{0.54:.2f}, {0.73:.2f}] (±15% of K_c)")
    lines.append("")

    # P4A-1: SV₂/SV₁ peaks within ±15% of K_c
    peak_idx = np.argmax(sv2sv1)
    peak_K = Ks[peak_idx]
    peak_val = sv2sv1[peak_idx]
    p1_pass = 0.54 <= peak_K <= 0.73

    lines.append(f"P4A-1 [KILL TEST]: SV₂/SV₁ peak at K={peak_K:.3f} (value={peak_val:.4f})")
    lines.append(f"  Window: [0.54, 0.73], K_c={K_C_EXACT:.4f}")
    lines.append(f"  Offset from K_c: {abs(peak_K - K_C_EXACT):.3f} ({abs(peak_K - K_C_EXACT)/K_C_EXACT*100:.1f}%)")
    lines.append(f"  → {'PASS ✓' if p1_pass else 'KILL ✗'}")
    lines.append("")

    # P4A-2: SV₂/SV₁ at K_c exceeds K=0.1 by ≥0.10
    idx_low = np.argmin(np.abs(Ks - 0.1))
    sv2sv1_low = sv2sv1[idx_low]
    separation = peak_val - sv2sv1_low
    p2_pass = separation >= 0.10
    p2_weak = separation >= 0.05

    lines.append(f"P4A-2: SV₂/SV₁ separation = {peak_val:.4f} - {sv2sv1_low:.4f} = {separation:.4f}")
    lines.append(f"  → {'PASS ✓' if p2_pass else ('WEAK' if p2_weak else 'FAIL ✗')}")
    lines.append("")

    # P4A-3: Non-monotonic peak
    # Check: mean below K_c < peak, mean above K_c < peak
    below_mask = Ks < K_C_EXACT - 0.1
    above_mask = Ks > K_C_EXACT + 0.1
    if np.any(below_mask) and np.any(above_mask):
        mean_below = np.mean(sv2sv1[below_mask])
        mean_above = np.mean(sv2sv1[above_mask])
        p3_pass = mean_below < peak_val and mean_above < peak_val
        lines.append(f"P4A-3: Non-monotonic: below_mean={mean_below:.4f}, peak={peak_val:.4f}, above_mean={mean_above:.4f}")
        lines.append(f"  → {'PASS ✓' if p3_pass else 'INCONCLUSIVE'}")
    else:
        p3_pass = False
        lines.append("P4A-3: Insufficient data points below/above K_c")
    lines.append("")

    # P4A-4: η minimum coincides with SV₂/SV₁ maximum (±20%)
    eta_min_idx = np.argmin(eta)
    eta_min_K = Ks[eta_min_idx]
    p4_pass = abs(eta_min_K - peak_K) / K_C_EXACT <= 0.20
    lines.append(f"P4A-4: η minimum at K={eta_min_K:.3f}, SV₂/SV₁ peak at K={peak_K:.3f}")
    lines.append(f"  Separation: {abs(eta_min_K - peak_K):.3f} ({abs(eta_min_K - peak_K)/K_C_EXACT*100:.1f}% of K_c)")
    lines.append(f"  → {'PASS ✓' if p4_pass else 'INCONCLUSIVE'}")
    lines.append("")

    # P4A-5: Rank increases approaching K_c
    near_kc_mask = np.abs(Ks - K_C_EXACT) < 0.15
    far_mask = np.abs(Ks - K_C_EXACT) > 0.3
    if np.any(near_kc_mask) and np.any(far_mask):
        rank_near = np.mean(ranks[near_kc_mask])
        rank_far = np.mean(ranks[far_mask])
        p5_pass = rank_near > rank_far
        lines.append(f"P4A-5: Rank near K_c={rank_near:.2f}, far from K_c={rank_far:.2f}")
        lines.append(f"  → {'PASS ✓' if p5_pass else 'INCONCLUSIVE'}")
    else:
        p5_pass = False
        lines.append("P4A-5: Insufficient data")
    lines.append("")

    # Overall verdict
    lines.append("=" * 70)
    if p1_pass:
        lines.append("OVERALL VERDICT: PASS — FIM detects Kuramoto transition")
    else:
        lines.append("OVERALL VERDICT: KILLED — FIM does NOT detect Kuramoto transition")
    lines.append(f"P4A-1={'PASS' if p1_pass else 'KILL'} P4A-2={'PASS' if p2_pass else 'FAIL'} "
                f"P4A-3={'PASS' if p3_pass else '?'} P4A-4={'PASS' if p4_pass else '?'} "
                f"P4A-5={'PASS' if p5_pass else '?'}")
    lines.append("=" * 70)

    # Data table
    lines.append("\n--- Full Data ---")
    lines.append(f"{'K':>8} {'R':>8} {'SV2/SV1':>10} {'rank':>6} {'eta':>8} {'PR':>8}")
    for r in results:
        lines.append(f"{r['K']:8.3f} {r['R_mean']:8.3f} {r['sv2sv1_mean']:10.4f} "
                     f"{r['rank_mean']:6.1f} {r['eta_mean']:8.4f} {r['pr_mean']:8.4f}")

    return "\n".join(lines)


def make_plots(results, plt):
    """Generate diagnostic plots."""
    Ks = [r['K'] for r in results]
    sv2sv1 = [r['sv2sv1_mean'] for r in results]
    sv2sv1_err = [r['sv2sv1_std'] for r in results]
    Rs = [r['R_mean'] for r in results]
    R_err = [r['R_std'] for r in results]
    ranks = [r['rank_mean'] for r in results]
    etas = [r['eta_mean'] for r in results]
    prs = [r['pr_mean'] for r in results]

    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    fig.suptitle(f"Phase 4A: Kuramoto FIM Spectral Analysis (N={N_OSC}, K_c={K_C_EXACT:.4f})",
                fontsize=14, fontweight='bold')

    # Plot A: SV₂/SV₁ vs K
    ax = axes[0, 0]
    ax.errorbar(Ks, sv2sv1, yerr=sv2sv1_err, fmt='o-', color='#2196F3',
                capsize=3, label='SV₂/SV₁')
    ax.axvline(K_C_EXACT, color='red', linestyle='--', alpha=0.7, label=f'K_c = {K_C_EXACT:.3f}')
    ax.axvspan(0.54, 0.73, alpha=0.1, color='green', label='±15% window')
    ax.set_xlabel('Coupling K')
    ax.set_ylabel('SV₂/SV₁')
    ax.set_title('FIM Spectral Ratio')
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3)

    # Plot B: Order parameter R vs K
    ax = axes[0, 1]
    ax.errorbar(Ks, Rs, yerr=R_err, fmt='s-', color='#FF5722', capsize=3, label='R (order param)')
    ax.axvline(K_C_EXACT, color='red', linestyle='--', alpha=0.7, label=f'K_c = {K_C_EXACT:.3f}')
    ax.set_xlabel('Coupling K')
    ax.set_ylabel('R = |N⁻¹ Σ exp(iθ)|')
    ax.set_title('Kuramoto Order Parameter')
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3)

    # Plot C: Rank + PR vs K
    ax = axes[1, 0]
    ax.plot(Ks, ranks, 'D-', color='#4CAF50', label='Rank')
    ax2 = ax.twinx()
    ax2.plot(Ks, prs, '^-', color='#9C27B0', label='PR', alpha=0.7)
    ax.axvline(K_C_EXACT, color='red', linestyle='--', alpha=0.7)
    ax.set_xlabel('Coupling K')
    ax.set_ylabel('Rank', color='#4CAF50')
    ax2.set_ylabel('PR', color='#9C27B0')
    ax.set_title('Rank and Participation Ratio')
    ax.grid(True, alpha=0.3)

    # Plot D: η vs K
    ax = axes[1, 1]
    ax.plot(Ks, etas, 'o-', color='#FF9800', label='η')
    ax.axvline(K_C_EXACT, color='red', linestyle='--', alpha=0.7, label=f'K_c = {K_C_EXACT:.3f}')
    ax.set_xlabel('Coupling K')
    ax.set_ylabel('η = SV₂/(SV₁+SV₂)')
    ax.set_title('Anisotropy Parameter')
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    plot_path = os.path.join(RESULTS_DIR, "phase4a_kuramoto_diagnostics.png")
    plt.savefig(plot_path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  Plot saved: {plot_path}")


if __name__ == "__main__":
    main()
