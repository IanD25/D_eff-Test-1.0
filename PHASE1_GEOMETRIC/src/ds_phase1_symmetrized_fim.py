#!/usr/bin/env python3
"""
DS Phase 1: Symmetrized Fisher Information Matrix
===================================================
Tests whether symmetrizing score vectors (pairing each neighbor with its
most anti-parallel counterpart) resolves the d+1 puzzle on random geometric graphs.

Systems tested:
  1. 2D torus (200×200) — control
  2. 3D torus (50×50×50) — control
  3. Periodic RGG d=2 (N=5000) — KEY TEST
  4. Periodic RGG d=3 (N=5000) — KEY TEST
  5. ER (n=1000) — negative control
  6. Sierpinski gasket L=7 — fractal check
"""

import os
import time
import datetime
from collections import deque, Counter

import numpy as np
import networkx as nx
from scipy.spatial import cKDTree
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

OUTPUT_DIR = "phase1_results"
np.random.seed(42)

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
# Gap-based rank
# ============================================================================

def gap_based_rank(sv_norm):
    if len(sv_norm) <= 1:
        return 1
    ratios = sv_norm[1:] / np.maximum(sv_norm[:-1], 1e-15)
    return int(np.argmin(ratios)) + 1


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


def make_periodic_rgg(n_points, d, r):
    points = np.random.uniform(0, 1, size=(n_points, d))
    tree = cKDTree(points, boxsize=1.0)
    pairs = tree.query_pairs(r)
    G = nx.Graph()
    G.add_nodes_from(range(n_points))
    G.add_edges_from(pairs)
    return G, points


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
# Symmetrized Fisher Information
# ============================================================================

def compute_symmetrized_fisher(adj, v0_idx, n, sigma):
    """
    Compute both standard and symmetrized FIM at vertex v0_idx.
    Uses adjacency list for fast BFS.
    """
    neighbors = adj[v0_idx]
    k = len(neighbors)
    if k < 2:
        return None

    # BFS from v0
    dist_v0 = bfs_distances_fast(adj, v0_idx, n)
    log_unnorm_v0 = -dist_v0.astype(float) / sigma
    log_unnorm_v0[dist_v0 < 0] = -1000.0
    log_Z_v0 = np.logaddexp.reduce(log_unnorm_v0)
    log_p_v0 = log_unnorm_v0 - log_Z_v0
    p_v0 = np.exp(log_p_v0)

    # Score vectors for all neighbors
    score_vectors = np.zeros((k, n))
    for j, w_idx in enumerate(neighbors):
        dist_w = bfs_distances_fast(adj, w_idx, n)
        log_unnorm_w = -dist_w.astype(float) / sigma
        log_unnorm_w[dist_w < 0] = -1000.0
        log_Z_w = np.logaddexp.reduce(log_unnorm_w)
        log_p_w = log_unnorm_w - log_Z_w
        score_vectors[j, :] = log_p_w - log_p_v0

    # Standard FIM
    weighted_scores = score_vectors * np.sqrt(p_v0)[np.newaxis, :]
    F_standard = weighted_scores @ weighted_scores.T

    if np.linalg.norm(F_standard, 'fro') < 1e-12:
        return None

    sv_standard = np.linalg.svd(F_standard, compute_uv=False)
    sv_standard_norm = sv_standard / sv_standard[0] if sv_standard[0] > 0 else sv_standard
    rank_standard = gap_based_rank(sv_standard_norm)

    # PR standard
    pr_standard = (np.sum(sv_standard))**2 / np.sum(sv_standard**2) if np.sum(sv_standard**2) > 0 else 0

    # Cosine similarity matrix (from Gram matrix F_standard)
    diag = np.sqrt(np.maximum(np.diag(F_standard), 1e-30))
    cosine_matrix = F_standard / np.outer(diag, diag)

    # Find anti-parallel partner for each neighbor
    cosine_for_pairing = cosine_matrix.copy()
    np.fill_diagonal(cosine_for_pairing, np.inf)
    opposites = np.argmin(cosine_for_pairing, axis=1)

    pairing_cosines = np.array([cosine_matrix[j, opposites[j]] for j in range(k)])
    mean_pairing_quality = np.mean(np.abs(pairing_cosines))

    # Symmetrized score vectors (all-k version)
    sym_scores = np.zeros((k, n))
    for j in range(k):
        sym_scores[j, :] = 0.5 * (score_vectors[j, :] - score_vectors[opposites[j], :])

    # Symmetrized FIM
    weighted_sym = sym_scores * np.sqrt(p_v0)[np.newaxis, :]
    F_sym = weighted_sym @ weighted_sym.T

    if np.linalg.norm(F_sym, 'fro') < 1e-12:
        return {
            'standard_svs': sv_standard_norm, 'standard_rank': rank_standard,
            'standard_pr': pr_standard,
            'symmetrized_svs': np.zeros(1), 'symmetrized_rank': 0,
            'symmetrized_pr': 0,
            'disorder_index': 0, 'pairing_quality': mean_pairing_quality,
            'pairing_cosines': pairing_cosines, 'degree': k,
        }

    sv_sym = np.linalg.svd(F_sym, compute_uv=False)
    sv_sym_norm = sv_sym / sv_sym[0] if sv_sym[0] > 1e-12 else sv_sym
    rank_sym = gap_based_rank(sv_sym_norm)

    pr_sym = (np.sum(sv_sym))**2 / np.sum(sv_sym**2) if np.sum(sv_sym**2) > 0 else 0

    # Disorder index
    if rank_sym < len(sv_standard_norm):
        disorder_index = sv_standard_norm[rank_sym] / sv_standard_norm[rank_sym - 1] if sv_standard_norm[rank_sym - 1] > 0 else 0
    else:
        disorder_index = 0.0

    return {
        'standard_svs': sv_standard_norm,
        'standard_rank': rank_standard,
        'standard_pr': pr_standard,
        'symmetrized_svs': sv_sym_norm,
        'symmetrized_rank': rank_sym,
        'symmetrized_pr': pr_sym,
        'disorder_index': disorder_index,
        'pairing_quality': mean_pairing_quality,
        'pairing_cosines': pairing_cosines,
        'degree': k,
    }


# ============================================================================
# Run analysis on a system
# ============================================================================

def run_system(G, name, sigma, sample_indices, true_d=None):
    """Run symmetrized Fisher analysis on a system."""
    print(f"\n  {name}: {G.number_of_nodes()} vertices, {G.number_of_edges()} edges")
    adj, nodes, node_to_idx, n = build_adjacency_list(G)

    results = []
    for v_idx in sample_indices:
        r = compute_symmetrized_fisher(adj, v_idx, n, sigma)
        if r is not None:
            r['vertex_idx'] = v_idx
            results.append(r)

    if not results:
        print(f"    ALL DEGENERATE")
        return None

    std_ranks = [r['standard_rank'] for r in results]
    sym_ranks = [r['symmetrized_rank'] for r in results]
    std_prs = [r['standard_pr'] for r in results]
    sym_prs = [r['symmetrized_pr'] for r in results]
    disorders = [r['disorder_index'] for r in results]
    pairings = [r['pairing_quality'] for r in results]

    print(f"    Standard  rank: mean={np.mean(std_ranks):.2f} median={np.median(std_ranks):.0f} "
          f"mode={Counter(std_ranks).most_common(1)[0][0]} "
          f"PR={np.mean(std_prs):.3f}")
    print(f"    Symmetrized rank: mean={np.mean(sym_ranks):.2f} median={np.median(sym_ranks):.0f} "
          f"mode={Counter(sym_ranks).most_common(1)[0][0]} "
          f"PR={np.mean(sym_prs):.3f}")
    print(f"    Disorder η: mean={np.mean(disorders):.3f} std={np.std(disorders):.3f}")
    print(f"    Pairing quality: mean={np.mean(pairings):.3f}")

    return {
        'name': name,
        'true_d': true_d,
        'sigma': sigma,
        'results': results,
        'std_ranks': std_ranks,
        'sym_ranks': sym_ranks,
        'std_prs': std_prs,
        'sym_prs': sym_prs,
        'disorders': disorders,
        'pairings': pairings,
    }


# ============================================================================
# Plotting
# ============================================================================

def plot_rank_comparison(data, filename):
    """Overlay standard vs symmetrized rank histograms."""
    fig, ax = plt.subplots(figsize=(9, 5))
    all_ranks = data['std_ranks'] + data['sym_ranks']
    bins = range(min(all_ranks), max(all_ranks) + 2)
    ax.hist(data['std_ranks'], bins=bins, align='left', alpha=0.5, edgecolor='black',
            label=f"Standard (mean={np.mean(data['std_ranks']):.2f})", color='coral')
    ax.hist(data['sym_ranks'], bins=bins, align='left', alpha=0.5, edgecolor='black',
            label=f"Symmetrized (mean={np.mean(data['sym_ranks']):.2f})", color='steelblue')
    if data['true_d'] is not None:
        ax.axvline(data['true_d'], color='red', linestyle='--', linewidth=2,
                  label=f"True d={data['true_d']}")
    ax.set_xlabel('Fisher Gap-Based Rank')
    ax.set_ylabel('Count')
    ax.set_title(f"Standard vs Symmetrized Rank: {data['name']}")
    ax.legend()
    plt.tight_layout()
    plt.savefig(os.path.join(OUTPUT_DIR, filename), dpi=150)
    plt.close()


def plot_sv_comparison(data, filename):
    """Side-by-side SV profiles: standard vs symmetrized."""
    results = data['results']
    degrees = [r['degree'] for r in results]
    modal_deg = Counter(degrees).most_common(1)[0][0]

    std_svs = [r['standard_svs'] for r in results if r['degree'] == modal_deg]
    sym_svs = [r['symmetrized_svs'] for r in results if r['degree'] == modal_deg]

    if not std_svs:
        closest = min(set(degrees), key=lambda d: abs(d - modal_deg))
        std_svs = [r['standard_svs'] for r in results if r['degree'] == closest]
        sym_svs = [r['symmetrized_svs'] for r in results if r['degree'] == closest]
        modal_deg = closest

    if not std_svs:
        return

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))

    # Standard
    max_len = max(len(sv) for sv in std_svs)
    padded = np.zeros((len(std_svs), max_len))
    for i, sv in enumerate(std_svs):
        padded[i, :len(sv)] = sv
    mean_sv = np.mean(padded, axis=0)
    std_sv = np.std(padded, axis=0)
    x = np.arange(1, max_len + 1)
    ax1.bar(x, mean_sv, yerr=std_sv, capsize=3, edgecolor='black', alpha=0.7, color='coral')
    ax1.set_xlabel('SV Index')
    ax1.set_ylabel('Normalized SV')
    ax1.set_title(f'Standard FIM (deg={modal_deg}, n={len(std_svs)})')

    # Symmetrized
    max_len2 = max(len(sv) for sv in sym_svs)
    padded2 = np.zeros((len(sym_svs), max_len2))
    for i, sv in enumerate(sym_svs):
        padded2[i, :len(sv)] = sv
    mean_sv2 = np.mean(padded2, axis=0)
    std_sv2 = np.std(padded2, axis=0)
    x2 = np.arange(1, max_len2 + 1)
    ax2.bar(x2, mean_sv2, yerr=std_sv2, capsize=3, edgecolor='black', alpha=0.7, color='steelblue')
    ax2.set_xlabel('SV Index')
    ax2.set_ylabel('Normalized SV')
    ax2.set_title(f'Symmetrized FIM (deg={modal_deg}, n={len(sym_svs)})')

    plt.suptitle(f'SV Profile Comparison: {data["name"]}', fontsize=13, y=1.02)
    plt.tight_layout()
    plt.savefig(os.path.join(OUTPUT_DIR, filename), dpi=150, bbox_inches='tight')
    plt.close()


def plot_disorder_histogram(data, filename):
    """Histogram of disorder index η."""
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.hist(data['disorders'], bins=20, edgecolor='black', alpha=0.7)
    ax.axvline(np.mean(data['disorders']), color='red', linestyle='--',
              label=f"Mean η={np.mean(data['disorders']):.3f}")
    ax.set_xlabel('Disorder Index η')
    ax.set_ylabel('Count')
    ax.set_title(f"Disorder Index Distribution: {data['name']}")
    ax.legend()
    plt.tight_layout()
    plt.savefig(os.path.join(OUTPUT_DIR, filename), dpi=150)
    plt.close()


def plot_pairing_quality_comparison(all_data, filename):
    """Bar chart of pairing quality across all systems."""
    fig, ax = plt.subplots(figsize=(10, 5))
    names = [d['name'] for d in all_data if d is not None]
    qualities = [np.mean(d['pairings']) for d in all_data if d is not None]
    colors = ['steelblue'] * len(names)

    bars = ax.bar(range(len(names)), qualities, color=colors, edgecolor='black', alpha=0.8)
    ax.set_xticks(range(len(names)))
    ax.set_xticklabels(names, rotation=30, ha='right', fontsize=9)
    ax.set_ylabel('Mean Pairing Quality |cos|')
    ax.set_title('Anti-Parallel Pairing Quality Across Systems')
    ax.set_ylim(0, 1.05)

    for bar, val in zip(bars, qualities):
        ax.text(bar.get_x() + bar.get_width()/2., bar.get_height() + 0.02,
                f'{val:.3f}', ha='center', va='bottom', fontsize=10, fontweight='bold')

    plt.tight_layout()
    plt.savefig(os.path.join(OUTPUT_DIR, filename), dpi=150)
    plt.close()


# ============================================================================
# Results report
# ============================================================================

def write_results(all_data):
    lines = []
    lines.append("# Phase 1: Symmetrized Fisher Information Matrix — Results")
    lines.append(f"\nGenerated: {datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    lines.append("")

    # Summary table
    lines.append("## Summary Table")
    lines.append("")
    lines.append("| System | Std Rank (mode) | Sym Rank (mode) | Std PR | Sym PR | Disorder η | Pairing Quality |")
    lines.append("|--------|----------------|-----------------|--------|--------|------------|-----------------|")

    for data in all_data:
        if data is None:
            continue
        std_mode = Counter(data['std_ranks']).most_common(1)[0][0]
        sym_mode = Counter(data['sym_ranks']).most_common(1)[0][0]
        lines.append(f"| {data['name']} | {std_mode} (mean={np.mean(data['std_ranks']):.2f}) | "
                    f"{sym_mode} (mean={np.mean(data['sym_ranks']):.2f}) | "
                    f"{np.mean(data['std_prs']):.3f} | {np.mean(data['sym_prs']):.3f} | "
                    f"{np.mean(data['disorders']):.3f} | {np.mean(data['pairings']):.3f} |")

    lines.append("")

    # Detailed results per system
    for data in all_data:
        if data is None:
            continue
        lines.append(f"## {data['name']}")
        lines.append("")
        lines.append(f"- σ = {data['sigma']}")
        lines.append(f"- True d = {data['true_d']}")
        lines.append(f"- N samples = {len(data['results'])}")
        lines.append("")

        lines.append(f"**Standard rank distribution:** {dict(sorted(Counter(data['std_ranks']).items()))}")
        lines.append(f"**Symmetrized rank distribution:** {dict(sorted(Counter(data['sym_ranks']).items()))}")
        lines.append("")

        lines.append(f"**Standard PR:** mean={np.mean(data['std_prs']):.3f}, std={np.std(data['std_prs']):.3f}")
        lines.append(f"**Symmetrized PR:** mean={np.mean(data['sym_prs']):.3f}, std={np.std(data['sym_prs']):.3f}")
        lines.append("")

        lines.append(f"**Disorder index η:** mean={np.mean(data['disorders']):.3f}, std={np.std(data['disorders']):.3f}")
        lines.append(f"**Pairing quality:** mean={np.mean(data['pairings']):.3f}, std={np.std(data['pairings']):.3f}")
        lines.append("")
        lines.append("---")
        lines.append("")

    path = os.path.join(OUTPUT_DIR, "PHASE1_SYMMETRIZED_FIM_RESULTS.md")
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
    print("DS Phase 1: Symmetrized Fisher Information Matrix")
    print("=" * 70)

    all_data = []
    sigma_main = 3.0

    # --- System 1: 2D Torus (200×200) ---
    print(f"\n{'='*70}")
    print("  System 1: 2D Torus (200×200)")
    print(f"{'='*70}")
    G = make_torus_2d(200)
    adj, nodes, node_to_idx, n = build_adjacency_list(G)
    samples = np.random.choice(n, 20, replace=False).tolist()
    data_torus2d = run_system(G, "2D Torus 200²", sigma_main, samples, true_d=2)
    all_data.append(data_torus2d)
    del G  # free memory

    # --- System 2: 3D Torus (50×50×50) ---
    print(f"\n{'='*70}")
    print("  System 2: 3D Torus (50×50×50)")
    print(f"{'='*70}")
    G = make_torus_3d(50)
    adj, nodes, node_to_idx, n = build_adjacency_list(G)
    samples = np.random.choice(n, 20, replace=False).tolist()
    data_torus3d = run_system(G, "3D Torus 50³", sigma_main, samples, true_d=3)
    all_data.append(data_torus3d)
    del G

    # --- System 3: Periodic RGG d=2 (N=5000) ---
    print(f"\n{'='*70}")
    print("  System 3: Periodic RGG d=2 (N=5000)")
    print(f"{'='*70}")
    G, _ = make_periodic_rgg(5000, 2, 0.03)
    if not nx.is_connected(G):
        lcc = max(nx.connected_components(G), key=len)
        G = G.subgraph(lcc).copy()
        mapping = {v: i for i, v in enumerate(sorted(G.nodes()))}
        G = nx.relabel_nodes(G, mapping)
    adj, nodes, node_to_idx, n = build_adjacency_list(G)
    degrees = [len(adj[i]) for i in range(n)]
    high_deg = [i for i in range(n) if degrees[i] >= 8]
    np.random.shuffle(high_deg)
    samples = high_deg[:50]
    data_rgg2d = run_system(G, "Periodic RGG d=2", sigma_main, samples, true_d=2)
    all_data.append(data_rgg2d)

    # Also run at σ=5.0 for robustness
    print("  (σ=5.0 robustness check)")
    data_rgg2d_s5 = run_system(G, "Periodic RGG d=2 (σ=5)", 5.0, samples, true_d=2)
    del G

    # --- System 4: Periodic RGG d=3 (N=5000) ---
    print(f"\n{'='*70}")
    print("  System 4: Periodic RGG d=3 (N=5000)")
    print(f"{'='*70}")
    G, _ = make_periodic_rgg(5000, 3, 0.09)
    if not nx.is_connected(G):
        lcc = max(nx.connected_components(G), key=len)
        G = G.subgraph(lcc).copy()
        mapping = {v: i for i, v in enumerate(sorted(G.nodes()))}
        G = nx.relabel_nodes(G, mapping)
    adj, nodes, node_to_idx, n = build_adjacency_list(G)
    degrees = [len(adj[i]) for i in range(n)]
    high_deg = [i for i in range(n) if degrees[i] >= 8]
    np.random.shuffle(high_deg)
    samples = high_deg[:50]
    data_rgg3d = run_system(G, "Periodic RGG d=3", sigma_main, samples, true_d=3)
    all_data.append(data_rgg3d)

    # σ=5.0 robustness
    print("  (σ=5.0 robustness check)")
    data_rgg3d_s5 = run_system(G, "Periodic RGG d=3 (σ=5)", 5.0, samples, true_d=3)
    del G

    # --- System 5: ER (n=1000) ---
    print(f"\n{'='*70}")
    print("  System 5: ER (n=1000, p=0.015)")
    print(f"{'='*70}")
    G = nx.erdos_renyi_graph(1000, 0.015, seed=42)
    if not nx.is_connected(G):
        lcc = max(nx.connected_components(G), key=len)
        G = G.subgraph(lcc).copy()
        mapping = {v: i for i, v in enumerate(sorted(G.nodes()))}
        G = nx.relabel_nodes(G, mapping)
    samples = np.random.choice(G.number_of_nodes(), 20, replace=False).tolist()
    data_er = run_system(G, "ER n=1000", sigma_main, samples, true_d=None)
    all_data.append(data_er)
    del G

    # --- System 6: Sierpinski Gasket L=7 ---
    print(f"\n{'='*70}")
    print("  System 6: Sierpinski Gasket L=7")
    print(f"{'='*70}")
    G, _ = make_sierpinski_gasket(7)
    adj, nodes, node_to_idx, n = build_adjacency_list(G)
    degrees = [len(adj[i]) for i in range(n)]
    # Interior vertices have degree 4
    interior = [i for i in range(n) if degrees[i] == 4]
    np.random.shuffle(interior)
    samples = interior[:30]
    GASKET_D_H = np.log(3) / np.log(2)
    data_gasket = run_system(G, "Gasket L=7", sigma_main, samples, true_d=GASKET_D_H)
    all_data.append(data_gasket)
    del G

    # ============================================================
    # PLOTS
    # ============================================================
    print(f"\n{'='*70}")
    print("  GENERATING PLOTS")
    print(f"{'='*70}")

    # 1. Rank comparison: RGG d=2
    if data_rgg2d:
        plot_rank_comparison(data_rgg2d, "symmetrized_rgg2d_rank_comparison.png")
        print("  Saved symmetrized_rgg2d_rank_comparison.png")

    # 2. Rank comparison: RGG d=3
    if data_rgg3d:
        plot_rank_comparison(data_rgg3d, "symmetrized_rgg3d_rank_comparison.png")
        print("  Saved symmetrized_rgg3d_rank_comparison.png")

    # 3. SV comparison: RGG d=2
    if data_rgg2d:
        plot_sv_comparison(data_rgg2d, "symmetrized_rgg2d_sv_comparison.png")
        print("  Saved symmetrized_rgg2d_sv_comparison.png")

    # 4. Disorder histogram: RGG d=2
    if data_rgg2d:
        plot_disorder_histogram(data_rgg2d, "symmetrized_rgg2d_disorder_hist.png")
        print("  Saved symmetrized_rgg2d_disorder_hist.png")

    # 5. Pairing quality comparison
    plot_pairing_quality_comparison(all_data, "symmetrized_pairing_quality_comparison.png")
    print("  Saved symmetrized_pairing_quality_comparison.png")

    # Results report
    write_results(all_data)

    elapsed = time.time() - t0
    print(f"\n{'='*70}")
    print(f"  TOTAL TIME: {elapsed:.1f}s")
    print(f"{'='*70}")
