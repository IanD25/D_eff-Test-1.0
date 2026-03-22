"""
DS Phase 3F: Percolation Threshold FDS Test
Ontology Kill-Shot Experiment

Version: 1.0
Date: March 2026
Author: Ian Darling
Executor: Claude Code

KILL CONDITION:  SV2/SV1 < 0.70 at p_c across both D=2 and D=3 → ontology abandoned
CONFIRM CONDITION: SV2/SV1 >= 0.85 at p_c → State 2 confirmed, ontology survives

FIM construction is VERBATIM from Phase 1 Geometric regime:
  kernel: p(u|v) = exp(-d(u,v)/sigma) / Z
  score:  s_j(u) = log p(u|w_j) - log p(u|v)
  FIM:    F = (S * sqrt(p)) @ (S * sqrt(p)).T
  SVD:    sigma_values = svd(F); d_eff = argmin(sv_norm[i+1]/sv_norm[i]) + 1
"""

import json
import os
import sys
import time
from collections import defaultdict
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np

# ============================================================================
# CONSTANTS
# ============================================================================

SIGMA_DEFAULT = 3.0
N_REALIZATIONS = 5
N_CENTERS = 20
MIN_GIANT_FACTOR = 10  # giant must be >= MIN_GIANT_FACTOR * N_CENTERS

# 2D sweep parameters
SWEEP_2D = {
    "L_values": [32, 64, 128, 256],
    "p_values": np.linspace(0.35, 0.65, 13).tolist(),
    "p_c": 0.5,
}

# 3D sweep parameters
SWEEP_3D = {
    "L_values": [16, 32, 64],
    "p_values": np.linspace(0.18, 0.32, 11).tolist(),
    "p_c": 0.24881,
}

# Known ground truth
GROUND_TRUTH = {
    2: {"p_c": 0.5, "d_f": 91 / 48, "nu": 4 / 3},
    3: {"p_c": 0.24881, "d_f": 2.523, "nu": 0.876},
}

# Kill/confirm thresholds
KILL_THRESHOLD = 0.70
CONFIRM_THRESHOLD = 0.85

RESULTS_DIR = Path(__file__).parent.parent / "results"
RAW_DIR = RESULTS_DIR / "raw_data"


# ============================================================================
# MODULE 1: LATTICE AND PERCOLATION
# ============================================================================


def build_percolation_lattice(L, dim, p, seed=None):
    """
    Build a bond percolation lattice graph with periodic boundary conditions.

    Parameters:
        L    : int, linear size
        dim  : int, 2 or 3
        p    : float, bond occupation probability
        seed : int or None

    Returns:
        G       : nx.Graph (occupied bonds only)
        giant   : set of node indices in largest connected component
        giant_fraction : len(giant) / L^dim
    """
    rng = np.random.default_rng(seed)
    N = L ** dim
    G = nx.Graph()
    G.add_nodes_from(range(N))

    if dim == 2:
        for i in range(L):
            for j in range(L):
                node = i * L + j
                # Horizontal bond: (i, j) -- (i, (j+1) mod L)
                right = i * L + (j + 1) % L
                if rng.random() < p:
                    G.add_edge(node, right)
                # Vertical bond: (i, j) -- ((i+1) mod L, j)
                down = ((i + 1) % L) * L + j
                if rng.random() < p:
                    G.add_edge(node, down)
    elif dim == 3:
        for i in range(L):
            for j in range(L):
                for k in range(L):
                    node = i * L * L + j * L + k
                    # +x
                    nx_node = i * L * L + j * L + (k + 1) % L
                    if rng.random() < p:
                        G.add_edge(node, nx_node)
                    # +y
                    ny_node = i * L * L + ((j + 1) % L) * L + k
                    if rng.random() < p:
                        G.add_edge(node, ny_node)
                    # +z
                    nz_node = ((i + 1) % L) * L * L + j * L + k
                    if rng.random() < p:
                        G.add_edge(node, nz_node)
    else:
        raise ValueError(f"dim must be 2 or 3, got {dim}")

    # Find giant component
    if G.number_of_edges() == 0:
        return G, set(), 0.0

    components = sorted(nx.connected_components(G), key=len, reverse=True)
    giant = components[0]
    giant_fraction = len(giant) / N

    return G, giant, giant_fraction


# ============================================================================
# MODULE 2: GRAPH DISTANCES (BFS on percolating cluster)
# ============================================================================


def bfs_distances(G, source):
    """BFS shortest path distances from source. Returns dict {node: distance}."""
    return dict(nx.single_source_shortest_path_length(G, source))


def select_sample_centers(G_sub, giant_nodes, n_samples, rng):
    """
    Select n_samples center nodes from the giant component.
    Prefer nodes with degree >= 2 (avoid boundary leaves).
    Falls back to all nodes if fewer than n_samples have degree >= 2.
    """
    # Filter to degree >= 2
    high_deg = [n for n in giant_nodes if G_sub.degree(n) >= 2]
    if len(high_deg) >= n_samples:
        indices = rng.choice(len(high_deg), size=n_samples, replace=False)
        return [high_deg[i] for i in indices]

    # Fallback: use all available nodes
    giant_list = list(giant_nodes)
    if len(giant_list) <= n_samples:
        return giant_list
    indices = rng.choice(len(giant_list), size=n_samples, replace=False)
    return [giant_list[i] for i in indices]


# ============================================================================
# MODULE 3: FIM CONSTRUCTION (VERBATIM FROM PHASE 1)
# ============================================================================


def compute_fim_at_center(G_sub, center, sigma, nodes_list):
    """
    Compute FIM for a single center vertex on the giant component subgraph.

    This is the Phase 1 geometric regime pipeline:
      1. BFS distances from center to all nodes
      2. Exponential kernel → log-space probability distribution
      3. For each neighbor: score vector = log p_w - log p_v
      4. FIM = (S * sqrt(p)) @ (S * sqrt(p)).T
      5. SVD → diagnostics

    Parameters:
        G_sub      : nx.Graph, giant component subgraph
        center     : node index
        sigma      : float, kernel scale parameter
        nodes_list : list, ordered node indices in giant component

    Returns:
        dict with keys: sv1, sv2, sv3, sv4, sv2_sv1_ratio,
                        rank, pr, eta, n_neighbors, sv_profile
        or None if center has degree < 2
    """
    neighbors = list(G_sub.neighbors(center))
    k = len(neighbors)
    if k < 2:
        return None

    n = len(nodes_list)
    node_to_idx = {node: i for i, node in enumerate(nodes_list)}

    # BFS distances from center
    dists_center = bfs_distances(G_sub, center)
    dist_arr_center = np.full(n, 1e6)  # unreachable = large distance
    for node, d in dists_center.items():
        if node in node_to_idx:
            dist_arr_center[node_to_idx[node]] = d

    # Log-space probability distribution at center
    log_unnorm_v0 = -dist_arr_center / sigma
    log_Z_v0 = np.logaddexp.reduce(log_unnorm_v0)
    log_p_v0 = log_unnorm_v0 - log_Z_v0
    p_v0 = np.exp(log_p_v0)

    # Score vectors for each neighbor direction
    score_vectors = np.zeros((k, n))
    for j, w in enumerate(neighbors):
        dists_w = bfs_distances(G_sub, w)
        dist_arr_w = np.full(n, 1e6)
        for node, d in dists_w.items():
            if node in node_to_idx:
                dist_arr_w[node_to_idx[node]] = d

        log_unnorm_w = -dist_arr_w / sigma
        log_Z_w = np.logaddexp.reduce(log_unnorm_w)
        log_p_w = log_unnorm_w - log_Z_w
        score_vectors[j, :] = log_p_w - log_p_v0

    # Fisher Information Matrix
    weighted_scores = score_vectors * np.sqrt(p_v0)[np.newaxis, :]
    F = weighted_scores @ weighted_scores.T

    # SVD
    sv = np.linalg.svd(F, compute_uv=False)
    if sv[0] <= 0:
        return None

    sv_norm = sv / sv[0]

    # Gap-based rank
    if len(sv_norm) > 1:
        ratios = sv_norm[1:] / np.maximum(sv_norm[:-1], 1e-15)
        gap_rank = int(np.argmin(ratios) + 1)
    else:
        gap_rank = 1

    # Participation ratio
    pr = float(np.sum(sv) ** 2 / np.sum(sv ** 2)) if np.sum(sv ** 2) > 0 else 0.0

    # Disorder index eta
    if gap_rank < len(sv_norm):
        eta = float(sv_norm[gap_rank] / sv_norm[gap_rank - 1])
    else:
        eta = 0.0

    # Extract top 4 SVs (pad with 0 if fewer)
    sv_top = np.zeros(4)
    for i in range(min(4, len(sv_norm))):
        sv_top[i] = sv_norm[i]

    sv2_sv1 = float(sv_top[1]) if sv_top[0] > 0 else 0.0

    return {
        "sv1": float(sv_top[0]),
        "sv2": float(sv_top[1]),
        "sv3": float(sv_top[2]),
        "sv4": float(sv_top[3]),
        "sv2_sv1_ratio": sv2_sv1,
        "rank": gap_rank,
        "pr": pr,
        "eta": eta,
        "n_neighbors": k,
        "sv_profile": sv_norm.tolist(),
    }


# ============================================================================
# MODULE 4: SWEEP ACROSS p VALUES
# ============================================================================


def run_sweep(dim, L_values, p_values, n_realizations, n_centers, sigma,
              results_dir=None):
    """
    Run the full percolation FDS sweep for one dimension.

    For each (L, p): build n_realizations percolation lattices.
    For each realization: sample n_centers from giant component, compute FIM.
    Save intermediate results after each (L, p) pair.

    Returns:
        all_results: list of dicts (one per data point)
    """
    all_results = []
    total_pairs = len(L_values) * len(p_values)
    pair_idx = 0

    for L in L_values:
        for p in p_values:
            pair_idx += 1
            t0 = time.time()

            for r in range(n_realizations):
                seed = L * 10000 + int(p * 10000) + r
                G, giant, giant_frac = build_percolation_lattice(L, dim, p, seed=seed)

                # Check minimum giant component size
                if len(giant) < MIN_GIANT_FACTOR * n_centers:
                    all_results.append({
                        "dim": dim, "L": L, "p": round(p, 6),
                        "realization": r, "center": -1,
                        "giant_size": len(giant),
                        "giant_fraction": giant_frac,
                        "sv1": np.nan, "sv2": np.nan,
                        "sv3": np.nan, "sv4": np.nan,
                        "sv2_sv1_ratio": np.nan,
                        "rank": np.nan, "pr": np.nan, "eta": np.nan,
                        "below_threshold": True,
                    })
                    continue

                # Build giant component subgraph
                G_sub = G.subgraph(giant).copy()
                nodes_list = sorted(giant)

                # Select sample centers
                rng = np.random.default_rng(seed + 9999)
                centers = select_sample_centers(G_sub, giant, n_centers, rng)

                for c_idx, center in enumerate(centers):
                    result = compute_fim_at_center(G_sub, center, sigma, nodes_list)

                    if result is None:
                        continue

                    all_results.append({
                        "dim": dim, "L": L, "p": round(p, 6),
                        "realization": r, "center": int(center),
                        "giant_size": len(giant),
                        "giant_fraction": giant_frac,
                        "sv1": result["sv1"],
                        "sv2": result["sv2"],
                        "sv3": result["sv3"],
                        "sv4": result["sv4"],
                        "sv2_sv1_ratio": result["sv2_sv1_ratio"],
                        "rank": result["rank"],
                        "pr": result["pr"],
                        "eta": result["eta"],
                        "below_threshold": False,
                    })

            elapsed = time.time() - t0
            # Count valid results for this (L, p) pair
            valid = [r for r in all_results
                     if r["L"] == L and r["p"] == round(p, 6)
                     and not r.get("below_threshold", False)]
            n_valid = len(valid)
            mean_sv2 = np.nanmean([r["sv2_sv1_ratio"] for r in valid]) if valid else np.nan

            print(f"  [{pair_idx}/{total_pairs}] D={dim} L={L} p={p:.4f} "
                  f"| {n_valid} samples | mean SV2/SV1={mean_sv2:.4f} "
                  f"| {elapsed:.1f}s")

            # Save intermediate results
            if results_dir:
                raw_file = results_dir / f"results_{dim}D_L{L}.json"
                L_results = [r for r in all_results if r["L"] == L]
                with open(raw_file, "w") as f:
                    json.dump(L_results, f, indent=2, default=str)

    return all_results


# ============================================================================
# MODULE 5: PLOTS AND RESULTS REPORT
# ============================================================================


def aggregate_by_L_p(results, dim):
    """Aggregate results by (L, p): compute mean and std of SV2/SV1, PR, rank."""
    filtered = [r for r in results if r["dim"] == dim and not r.get("below_threshold")]
    grouped = defaultdict(list)
    for r in filtered:
        grouped[(r["L"], r["p"])].append(r)

    agg = {}
    for (L, p), recs in grouped.items():
        sv2_vals = [r["sv2_sv1_ratio"] for r in recs if not np.isnan(r["sv2_sv1_ratio"])]
        pr_vals = [r["pr"] for r in recs if not np.isnan(r["pr"])]
        rank_vals = [r["rank"] for r in recs if not np.isnan(r["rank"])]
        eta_vals = [r["eta"] for r in recs if not np.isnan(r["eta"])]
        giant_vals = [r["giant_fraction"] for r in recs]

        agg[(L, p)] = {
            "sv2_mean": np.mean(sv2_vals) if sv2_vals else np.nan,
            "sv2_std": np.std(sv2_vals) if sv2_vals else np.nan,
            "pr_mean": np.mean(pr_vals) if pr_vals else np.nan,
            "pr_std": np.std(pr_vals) if pr_vals else np.nan,
            "rank_mean": np.mean(rank_vals) if rank_vals else np.nan,
            "eta_mean": np.mean(eta_vals) if eta_vals else np.nan,
            "giant_mean": np.mean(giant_vals) if giant_vals else np.nan,
            "n_samples": len(sv2_vals),
        }
    return agg


def plot_sv2_vs_p(results, dim, sweep_config, out_dir):
    """Figure 1/2: SV2/SV1 vs p for each L, with kill/confirm lines."""
    agg = aggregate_by_L_p(results, dim)
    p_c = sweep_config["p_c"]

    fig, ax = plt.subplots(figsize=(10, 6))

    for L in sweep_config["L_values"]:
        ps = sorted(set(p for (l, p) in agg if l == L))
        means = [agg[(L, p)]["sv2_mean"] for p in ps if (L, p) in agg]
        stds = [agg[(L, p)]["sv2_std"] for p in ps if (L, p) in agg]
        ps_valid = [p for p in ps if (L, p) in agg]

        ax.errorbar(ps_valid, means, yerr=stds, marker="o", label=f"L={L}",
                     capsize=3, linewidth=1.5, markersize=4)

    ax.axvline(p_c, color="gray", linestyle="--", alpha=0.7, label=f"p_c={p_c}")
    ax.axhline(KILL_THRESHOLD, color="red", linestyle="--", alpha=0.7,
               label=f"Kill threshold ({KILL_THRESHOLD})")
    ax.axhline(CONFIRM_THRESHOLD, color="green", linestyle="--", alpha=0.7,
               label=f"Confirm threshold ({CONFIRM_THRESHOLD})")

    ax.set_xlabel("Bond probability p", fontsize=12)
    ax.set_ylabel("SV₂/SV₁", fontsize=12)
    ax.set_title(f"{dim}D Percolation: SV₂/SV₁ vs Bond Probability p", fontsize=14)
    ax.legend(fontsize=10)
    ax.grid(True, alpha=0.3)
    ax.set_ylim(0, 1.1)

    fig.tight_layout()
    fig.savefig(out_dir / f"sv2_sv1_vs_p_{dim}D.png", dpi=150)
    plt.close(fig)
    print(f"  Saved: sv2_sv1_vs_p_{dim}D.png")


def plot_sv_profile_at_threshold(results, dim, L_target, sweep_config, out_dir):
    """Figure 3/4: SV profile bar chart at below/at/above threshold."""
    p_c = sweep_config["p_c"]
    ps = sorted(sweep_config["p_values"])

    # Find below, at, above threshold p values
    p_below = max(p for p in ps if p < p_c - 0.05) if any(p < p_c - 0.05 for p in ps) else ps[0]
    p_at = min(ps, key=lambda x: abs(x - p_c))
    p_above = min(p for p in ps if p > p_c + 0.05) if any(p > p_c + 0.05 for p in ps) else ps[-1]

    fig, axes = plt.subplots(1, 3, figsize=(14, 5))
    labels_p = [
        (p_below, f"Below (p={p_below:.2f})"),
        (p_at, f"At threshold (p={p_at:.2f})"),
        (p_above, f"Above (p={p_above:.2f})"),
    ]

    for ax_idx, (p_val, label) in enumerate(labels_p):
        ax = axes[ax_idx]
        recs = [r for r in results
                if r["dim"] == dim and r["L"] == L_target
                and abs(r["p"] - p_val) < 0.001
                and not r.get("below_threshold")]

        if not recs:
            ax.text(0.5, 0.5, "No data", ha="center", va="center", transform=ax.transAxes)
            ax.set_title(label)
            continue

        # Average SV profile across all samples
        sv_means = [
            np.nanmean([r["sv1"] for r in recs]),
            np.nanmean([r["sv2"] for r in recs]),
            np.nanmean([r["sv3"] for r in recs]),
            np.nanmean([r["sv4"] for r in recs]),
        ]

        colors = ["#2196F3", "#FF9800", "#4CAF50", "#9C27B0"]
        bars = ax.bar(["SV₁", "SV₂", "SV₃", "SV₄"], sv_means, color=colors, alpha=0.8)
        ax.set_ylim(0, 1.15)
        ax.set_title(label, fontsize=11)
        ax.set_ylabel("Normalized SV" if ax_idx == 0 else "")

        # Annotate values
        for bar, val in zip(bars, sv_means):
            if not np.isnan(val):
                ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.02,
                        f"{val:.3f}", ha="center", va="bottom", fontsize=9)

    fig.suptitle(f"{dim}D Percolation SV Profile at L={L_target}", fontsize=14)
    fig.tight_layout()
    fig.savefig(out_dir / f"sv_profile_{dim}D_at_threshold.png", dpi=150)
    plt.close(fig)
    print(f"  Saved: sv_profile_{dim}D_at_threshold.png")


def plot_rank_and_pr(results, dim, sweep_config, out_dir):
    """Figure 5/6: Gap-based rank and PR vs p for each L."""
    agg = aggregate_by_L_p(results, dim)
    p_c = sweep_config["p_c"]

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))

    for L in sweep_config["L_values"]:
        ps = sorted(set(p for (l, p) in agg if l == L))
        ranks = [agg[(L, p)]["rank_mean"] for p in ps if (L, p) in agg]
        prs = [agg[(L, p)]["pr_mean"] for p in ps if (L, p) in agg]
        ps_valid = [p for p in ps if (L, p) in agg]

        ax1.plot(ps_valid, ranks, marker="o", label=f"L={L}", linewidth=1.5, markersize=4)
        ax2.plot(ps_valid, prs, marker="s", label=f"L={L}", linewidth=1.5, markersize=4)

    for ax in [ax1, ax2]:
        ax.axvline(p_c, color="gray", linestyle="--", alpha=0.7)
        ax.grid(True, alpha=0.3)
        ax.legend(fontsize=9)

    ax1.set_xlabel("Bond probability p")
    ax1.set_ylabel("Gap-based rank (d_eff)")
    ax1.set_title(f"{dim}D: Gap-Based Rank vs p")

    ax2.set_xlabel("Bond probability p")
    ax2.set_ylabel("Participation Ratio (PR)")
    ax2.set_title(f"{dim}D: Participation Ratio vs p")

    # Add d_f reference line
    d_f = GROUND_TRUTH[dim]["d_f"]
    ax2.axhline(d_f, color="red", linestyle=":", alpha=0.7, label=f"d_f={d_f:.3f}")
    ax2.legend(fontsize=9)

    fig.tight_layout()
    fig.savefig(out_dir / f"rank_and_pr_{dim}D.png", dpi=150)
    plt.close(fig)
    print(f"  Saved: rank_and_pr_{dim}D.png")


def plot_giant_fraction(results, sweep_2d, sweep_3d, out_dir):
    """Figure 7: Giant component fraction vs p for both dimensions."""
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))

    for dim, sweep_config, ax in [(2, sweep_2d, ax1), (3, sweep_3d, ax2)]:
        agg = aggregate_by_L_p(results, dim)
        # Also include below-threshold data for giant fraction
        all_filtered = [r for r in results if r["dim"] == dim]
        giant_agg = defaultdict(list)
        for r in all_filtered:
            giant_agg[(r["L"], r["p"])].append(r["giant_fraction"])

        for L in sweep_config["L_values"]:
            ps = sorted(set(p for (l, p) in giant_agg if l == L))
            means = [np.mean(giant_agg[(L, p)]) for p in ps if (L, p) in giant_agg]
            ps_valid = [p for p in ps if (L, p) in giant_agg]
            ax.plot(ps_valid, means, marker="o", label=f"L={L}",
                    linewidth=1.5, markersize=4)

        ax.axvline(sweep_config["p_c"], color="gray", linestyle="--", alpha=0.7,
                   label=f"p_c={sweep_config['p_c']}")
        ax.set_xlabel("Bond probability p")
        ax.set_ylabel("Giant component fraction")
        ax.set_title(f"{dim}D: Giant Fraction vs p")
        ax.legend(fontsize=9)
        ax.grid(True, alpha=0.3)

    fig.tight_layout()
    fig.savefig(out_dir / "giant_fraction_vs_p.png", dpi=150)
    plt.close(fig)
    print("  Saved: giant_fraction_vs_p.png")


def plot_finite_size_sv2(results, sweep_2d, sweep_3d, out_dir):
    """Figure 8: SV2/SV1 at p_c vs L for both dimensions."""
    fig, ax = plt.subplots(figsize=(8, 6))

    for dim, sweep_config, marker, color in [
        (2, sweep_2d, "o", "#2196F3"),
        (3, sweep_3d, "s", "#FF5722"),
    ]:
        agg = aggregate_by_L_p(results, dim)
        p_c = sweep_config["p_c"]

        # Find closest p to p_c
        ps = sorted(set(p for (_, p) in agg))
        p_nearest = min(ps, key=lambda x: abs(x - p_c)) if ps else None
        if p_nearest is None:
            continue

        Ls = sorted(sweep_config["L_values"])
        sv2_at_pc = []
        sv2_err = []
        Ls_valid = []

        for L in Ls:
            if (L, p_nearest) in agg:
                sv2_at_pc.append(agg[(L, p_nearest)]["sv2_mean"])
                sv2_err.append(agg[(L, p_nearest)]["sv2_std"])
                Ls_valid.append(L)

        if Ls_valid:
            ax.errorbar(Ls_valid, sv2_at_pc, yerr=sv2_err, marker=marker,
                        color=color, label=f"{dim}D (p≈{p_nearest:.4f})",
                        capsize=3, linewidth=2, markersize=6)

    ax.axhline(KILL_THRESHOLD, color="red", linestyle="--", alpha=0.7,
               label=f"Kill ({KILL_THRESHOLD})")
    ax.axhline(CONFIRM_THRESHOLD, color="green", linestyle="--", alpha=0.7,
               label=f"Confirm ({CONFIRM_THRESHOLD})")

    ax.set_xlabel("Lattice size L", fontsize=12)
    ax.set_ylabel("SV₂/SV₁ at p_c", fontsize=12)
    ax.set_title("Finite-Size Scaling: SV₂/SV₁ at Threshold", fontsize=14)
    ax.legend(fontsize=10)
    ax.grid(True, alpha=0.3)
    ax.set_xscale("log", base=2)

    fig.tight_layout()
    fig.savefig(out_dir / "finite_size_sv2_at_threshold.png", dpi=150)
    plt.close(fig)
    print("  Saved: finite_size_sv2_at_threshold.png")


def generate_report(results, out_dir):
    """Generate PHASE3F_PERCOLATION_RESULTS.md"""
    lines = []
    lines.append("# Phase 3F: Percolation FDS Results\n")
    lines.append(f"Generated: {time.strftime('%Y-%m-%d %H:%M:%S')}\n")
    lines.append(f"Sigma: {SIGMA_DEFAULT} | Realizations: {N_REALIZATIONS} "
                 f"| Centers/realization: {N_CENTERS}\n\n")

    # ---- Gate Verification ----
    lines.append("## Gate Verification\n\n")

    for dim in [2, 3]:
        p_c = GROUND_TRUTH[dim]["p_c"]
        agg = aggregate_by_L_p(results, dim)

        # Gate A: Giant component forms near p_c
        ps = sorted(set(p for (_, p) in agg))
        if ps:
            p_nearest = min(ps, key=lambda x: abs(x - p_c))
            Ls = sorted(set(l for (l, _) in agg))
            giant_at_pc = []
            for L in Ls:
                if (L, p_nearest) in agg:
                    giant_at_pc.append(agg[(L, p_nearest)]["giant_mean"])

            gate_a = all(g > 0.3 for g in giant_at_pc) if giant_at_pc else False
            lines.append(f"**Gate A ({dim}D):** Giant fraction at p≈{p_nearest:.4f}: "
                         f"{[f'{g:.3f}' for g in giant_at_pc]} → "
                         f"{'PASS' if gate_a else 'FAIL'}\n\n")

        # Gate B: PR at threshold vs known d_f
        d_f = GROUND_TRUTH[dim]["d_f"]
        pr_at_pc = []
        for L in Ls:
            if (L, p_nearest) in agg:
                pr_at_pc.append((L, agg[(L, p_nearest)]["pr_mean"]))

        if pr_at_pc:
            lines.append(f"**Gate B ({dim}D):** PR at threshold (d_f={d_f:.3f}):\n")
            for L, pr in pr_at_pc:
                pct_off = abs(pr - d_f) / d_f * 100 if d_f > 0 else 0
                lines.append(f"  - L={L}: PR={pr:.3f} ({pct_off:.1f}% from d_f)\n")
            lines.append("\n")

    # ---- Primary Result: Kill Condition ----
    lines.append("## Primary Result: Kill Condition Assessment\n\n")

    kill_triggered = True  # Assume kill unless proven otherwise
    for dim in [2, 3]:
        p_c = GROUND_TRUTH[dim]["p_c"]
        agg = aggregate_by_L_p(results, dim)
        ps = sorted(set(p for (_, p) in agg))
        if not ps:
            lines.append(f"**{dim}D:** No valid data\n\n")
            continue

        p_nearest = min(ps, key=lambda x: abs(x - p_c))
        Ls = sorted(set(l for (l, _) in agg))

        lines.append(f"### {dim}D (p_c = {p_c}, measured at p = {p_nearest:.4f})\n\n")
        lines.append("| L | SV₂/SV₁ | PR | Rank | η | Verdict |\n")
        lines.append("|---|---------|-----|------|---|--------|\n")

        for L in Ls:
            if (L, p_nearest) not in agg:
                continue
            a = agg[(L, p_nearest)]
            sv2 = a["sv2_mean"]
            pr = a["pr_mean"]
            rank = a["rank_mean"]
            eta = a["eta_mean"]

            if sv2 >= CONFIRM_THRESHOLD:
                verdict = "✅ CONFIRM"
                kill_triggered = False
            elif sv2 < KILL_THRESHOLD:
                verdict = "❌ KILL"
            else:
                verdict = "⚠️ INCONCLUSIVE"
                kill_triggered = False

            lines.append(f"| {L} | {sv2:.4f} | {pr:.3f} | {rank:.1f} | "
                         f"{eta:.4f} | {verdict} |\n")

        lines.append("\n")

    # ---- Overall Verdict ----
    lines.append("## Overall Verdict\n\n")

    # Check all (dim, L) pairs at threshold
    all_sv2 = []
    for dim in [2, 3]:
        agg = aggregate_by_L_p(results, dim)
        ps = sorted(set(p for (_, p) in agg))
        if not ps:
            continue
        p_c = GROUND_TRUTH[dim]["p_c"]
        p_nearest = min(ps, key=lambda x: abs(x - p_c))
        for L in sorted(set(l for (l, _) in agg)):
            if (L, p_nearest) in agg:
                all_sv2.append(agg[(L, p_nearest)]["sv2_mean"])

    if not all_sv2:
        lines.append("**NO DATA — cannot assess.**\n\n")
    elif all(sv2 < KILL_THRESHOLD for sv2 in all_sv2):
        lines.append("**🔴 KILL SHOT TRIGGERED:** SV₂/SV₁ < 0.70 across all "
                     "lattice sizes in both dimensions.\n\n")
        lines.append("The two-dimensional ontology is **ABANDONED**. The SV swap "
                     "is a thermal artifact, not a universal signature of "
                     "dimensional coupling transitions.\n\n")
    elif all(sv2 >= CONFIRM_THRESHOLD for sv2 in all_sv2):
        lines.append("**🟢 ONTOLOGY SURVIVES:** SV₂/SV₁ ≥ 0.85 at threshold "
                     "across all lattice sizes in both dimensions.\n\n")
        lines.append("State 2 isotropic profile confirmed at percolation threshold. "
                     "The SV swap is NOT a thermal artifact — it appears in a "
                     "non-thermal geometric phase transition.\n\n")
    else:
        max_sv2 = max(all_sv2)
        min_sv2 = min(all_sv2)
        lines.append(f"**🟡 INCONCLUSIVE:** SV₂/SV₁ ranges from {min_sv2:.4f} to "
                     f"{max_sv2:.4f}. Between kill and confirm thresholds.\n\n")

    # ---- Prediction Outcomes ----
    lines.append("## Pre-Registered Prediction Outcomes\n\n")

    # P3F-1: SV2/SV1 >= 0.85 at p_c
    if all_sv2:
        mean_all = np.mean(all_sv2)
        lines.append(f"**P3F-1 (PRIMARY):** Mean SV₂/SV₁ at p_c = {mean_all:.4f} → "
                     f"{'PASS' if mean_all >= CONFIRM_THRESHOLD else 'FAIL' if mean_all < KILL_THRESHOLD else 'PARTIAL'}\n\n")

    # P3F-3: Peak at p_c
    for dim in [2, 3]:
        agg = aggregate_by_L_p(results, dim)
        Ls = sorted(set(l for (l, _) in agg))
        if not Ls:
            continue
        L_largest = Ls[-1]
        ps_for_L = sorted(set(p for (l, p) in agg if l == L_largest))
        if not ps_for_L:
            continue
        sv2_by_p = {p: agg[(L_largest, p)]["sv2_mean"]
                    for p in ps_for_L if (L_largest, p) in agg}
        if sv2_by_p:
            p_max = max(sv2_by_p, key=sv2_by_p.get)
            p_c = GROUND_TRUTH[dim]["p_c"]
            p_step = abs(ps_for_L[1] - ps_for_L[0]) if len(ps_for_L) > 1 else 0.05
            peak_near_pc = abs(p_max - p_c) <= p_step * 1.5
            lines.append(f"**P3F-3 ({dim}D):** SV₂/SV₁ peak at p={p_max:.4f} "
                         f"(p_c={p_c:.4f}) → "
                         f"{'PASS' if peak_near_pc else 'FAIL'}\n\n")

    # P3F-4: Finite-size scaling
    for dim in [2, 3]:
        agg = aggregate_by_L_p(results, dim)
        p_c = GROUND_TRUTH[dim]["p_c"]
        ps = sorted(set(p for (_, p) in agg))
        if not ps:
            continue
        p_nearest = min(ps, key=lambda x: abs(x - p_c))
        Ls = sorted(set(l for (l, _) in agg))
        sv2_by_L = [(L, agg[(L, p_nearest)]["sv2_mean"])
                    for L in Ls if (L, p_nearest) in agg]
        if len(sv2_by_L) >= 2:
            monotone = all(sv2_by_L[i][1] <= sv2_by_L[i + 1][1]
                           for i in range(len(sv2_by_L) - 1))
            lines.append(f"**P3F-4 ({dim}D):** SV₂/SV₁ at p_c vs L: "
                         f"{[(L, f'{v:.4f}') for L, v in sv2_by_L]} → "
                         f"{'PASS (monotonically increasing)' if monotone else 'FAIL'}\n\n")

    # ---- Next Steps ----
    lines.append("## Next Step Recommendation\n\n")
    if all_sv2 and all(sv2 < KILL_THRESHOLD for sv2 in all_sv2):
        lines.append("- FDS thermal regime is domain-specific to correlation-based transitions.\n")
        lines.append("- Paper 2 stands. Broader ontological claim does not.\n")
        lines.append("- Next: submit papers as-is.\n")
    elif all_sv2 and all(sv2 >= CONFIRM_THRESHOLD for sv2 in all_sv2):
        lines.append("- Next test: Erdős-Rényi giant component transition.\n")
        lines.append("  - Random graph analog: no spatial structure at all.\n")
        lines.append("  - If swap appears: claim is about graph topology generically.\n")
        lines.append("  - If swap disappears: claim requires spatial embedding.\n")
    else:
        lines.append("- Run L=512 in 2D to resolve inconclusive range.\n")
        lines.append("- Check sigma sensitivity: run sigma=2.0 and sigma=5.0.\n")

    report = "".join(lines)
    report_path = out_dir / "PHASE3F_PERCOLATION_RESULTS.md"
    with open(report_path, "w") as f:
        f.write(report)
    print(f"  Saved: PHASE3F_PERCOLATION_RESULTS.md")

    return report


# ============================================================================
# MAIN
# ============================================================================


def main():
    """Run the full Phase 3F percolation experiment."""
    print("=" * 70)
    print("Phase 3F: Percolation Threshold FDS Test")
    print("Ontology Kill-Shot Experiment")
    print("=" * 70)
    print()

    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    RAW_DIR.mkdir(parents=True, exist_ok=True)

    all_results = []

    # ---- 2D Sweep ----
    print("━" * 50)
    print("2D Bond Percolation (square lattice, periodic BC)")
    print(f"  L values: {SWEEP_2D['L_values']}")
    print(f"  p range: {SWEEP_2D['p_values'][0]:.2f} – {SWEEP_2D['p_values'][-1]:.2f} "
          f"({len(SWEEP_2D['p_values'])} steps)")
    print(f"  p_c = {SWEEP_2D['p_c']} (exact)")
    print(f"  σ = {SIGMA_DEFAULT}, {N_REALIZATIONS} realizations, "
          f"{N_CENTERS} centers each")
    print("━" * 50)
    t0 = time.time()

    results_2d = run_sweep(
        dim=2,
        L_values=SWEEP_2D["L_values"],
        p_values=SWEEP_2D["p_values"],
        n_realizations=N_REALIZATIONS,
        n_centers=N_CENTERS,
        sigma=SIGMA_DEFAULT,
        results_dir=RAW_DIR,
    )
    all_results.extend(results_2d)
    t_2d = time.time() - t0
    print(f"\n2D complete in {t_2d:.1f}s ({t_2d / 60:.1f} min)\n")

    # ---- 3D Sweep ----
    # Start with L=16 and L=32 only; L=64 if time allows
    L_3d = SWEEP_3D["L_values"][:2]  # [16, 32]
    remaining_budget = 20 * 60 - t_2d  # 20 min total budget minus 2D time

    print("━" * 50)
    print("3D Bond Percolation (cubic lattice, periodic BC)")
    print(f"  L values: {L_3d} (L=64 if time allows)")
    print(f"  p range: {SWEEP_3D['p_values'][0]:.2f} – {SWEEP_3D['p_values'][-1]:.2f} "
          f"({len(SWEEP_3D['p_values'])} steps)")
    print(f"  p_c = {SWEEP_3D['p_c']}")
    print(f"  Remaining budget: {remaining_budget:.0f}s ({remaining_budget / 60:.1f} min)")
    print("━" * 50)
    t1 = time.time()

    results_3d = run_sweep(
        dim=3,
        L_values=L_3d,
        p_values=SWEEP_3D["p_values"],
        n_realizations=N_REALIZATIONS,
        n_centers=N_CENTERS,
        sigma=SIGMA_DEFAULT,
        results_dir=RAW_DIR,
    )
    all_results.extend(results_3d)
    t_3d_partial = time.time() - t1
    print(f"\n3D (L≤32) complete in {t_3d_partial:.1f}s ({t_3d_partial / 60:.1f} min)\n")

    # Try L=64 if we have time
    elapsed_total = time.time() - t0
    if elapsed_total < 15 * 60:  # If we're under 15 min total
        print("  Time available — running 3D L=64...")
        t2 = time.time()
        results_3d_64 = run_sweep(
            dim=3,
            L_values=[64],
            p_values=SWEEP_3D["p_values"],
            n_realizations=N_REALIZATIONS,
            n_centers=N_CENTERS,
            sigma=SIGMA_DEFAULT,
            results_dir=RAW_DIR,
        )
        all_results.extend(results_3d_64)
        t_3d_64 = time.time() - t2
        print(f"\n3D L=64 complete in {t_3d_64:.1f}s ({t_3d_64 / 60:.1f} min)\n")
    else:
        print(f"  Skipping 3D L=64 (elapsed {elapsed_total / 60:.1f} min > 15 min)\n")

    # ---- Generate Plots ----
    print("━" * 50)
    print("Generating figures...")
    print("━" * 50)

    plot_sv2_vs_p(all_results, 2, SWEEP_2D, RESULTS_DIR)
    plot_sv2_vs_p(all_results, 3, SWEEP_3D, RESULTS_DIR)
    plot_sv_profile_at_threshold(all_results, 2, 128, SWEEP_2D, RESULTS_DIR)
    plot_sv_profile_at_threshold(all_results, 3, 32, SWEEP_3D, RESULTS_DIR)
    plot_rank_and_pr(all_results, 2, SWEEP_2D, RESULTS_DIR)
    plot_rank_and_pr(all_results, 3, SWEEP_3D, RESULTS_DIR)
    plot_giant_fraction(all_results, SWEEP_2D, SWEEP_3D, RESULTS_DIR)
    plot_finite_size_sv2(all_results, SWEEP_2D, SWEEP_3D, RESULTS_DIR)

    # ---- Generate Report ----
    print("━" * 50)
    print("Generating results report...")
    print("━" * 50)

    report = generate_report(all_results, RESULTS_DIR)

    # ---- Save all raw data ----
    with open(RAW_DIR / "all_results.json", "w") as f:
        json.dump(all_results, f, indent=2, default=str)
    print("  Saved: raw_data/all_results.json")

    total_time = time.time() - t0
    print()
    print("=" * 70)
    print(f"Phase 3F COMPLETE — Total time: {total_time:.0f}s ({total_time / 60:.1f} min)")
    print(f"Results in: {RESULTS_DIR}")
    print("=" * 70)

    # Print the verdict to stdout
    print()
    for line in report.split("\n"):
        if "KILL SHOT" in line or "ONTOLOGY SURVIVES" in line or "INCONCLUSIVE" in line:
            print(line)
            break


if __name__ == "__main__":
    main()
