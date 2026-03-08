#!/usr/bin/env python3
"""
Phase 3B: Financial Correlation Network — Standalone Fisher Diagnostics
=======================================================================
Same scientific computation as the LEAN algorithm, but uses yfinance for
data access. This avoids LEAN deployment issues while producing identical
Fisher diagnostic results.

Computes weekly Fisher diagnostics on a k-NN graph built from rolling
correlation matrices of ~100 S&P 500 stocks, 2005-2024.

Output: fisher_financial_results.json (same format as LEAN ObjectStore output)
"""

import numpy as np
import json
import os
import sys
import time
from collections import deque
from datetime import datetime, timedelta

# ============================================================================
# Data acquisition via yfinance
# ============================================================================

def download_data(tickers, start='2004-06-01', end='2024-12-31'):
    """Download daily adjusted close prices for all tickers."""
    import yfinance as yf

    print(f"Downloading data for {len(tickers)} tickers...")
    print(f"  Period: {start} to {end}")

    # Download in batches to avoid rate limits
    all_data = {}
    batch_size = 20
    for i in range(0, len(tickers), batch_size):
        batch = tickers[i:i + batch_size]
        print(f"  Batch {i // batch_size + 1}: {batch[:5]}...")
        try:
            data = yf.download(batch, start=start, end=end,
                             auto_adjust=True, progress=False)
            if len(batch) == 1:
                # Single ticker returns Series not DataFrame
                all_data[batch[0]] = data['Close']
            else:
                for ticker in batch:
                    if ticker in data['Close'].columns:
                        all_data[ticker] = data['Close'][ticker]
        except Exception as e:
            print(f"  Warning: batch download failed: {e}")
            # Try individual downloads
            for ticker in batch:
                try:
                    t = yf.Ticker(ticker)
                    hist = t.history(start=start, end=end)
                    if not hist.empty:
                        all_data[ticker] = hist['Close']
                except Exception:
                    print(f"    Skipped {ticker}")

    print(f"  Downloaded {len(all_data)} tickers successfully")
    return all_data


# ============================================================================
# Graph construction and C(r) computation
# ============================================================================

def build_knn_graph(corr_matrix, k):
    """Build k-NN graph from correlation matrix. Symmetrize."""
    n = corr_matrix.shape[0]
    adjacency = np.zeros((n, n), dtype=bool)

    for i in range(n):
        row = corr_matrix[i].copy()
        row[i] = -np.inf
        top_k = np.argsort(row)[-k:]
        for j in top_k:
            adjacency[i, j] = True
            adjacency[j, i] = True

    neighbors = {}
    for i in range(n):
        neighbors[i] = list(np.where(adjacency[i])[0])

    return adjacency, neighbors


def bfs_distances(nbr_list, source, n, max_dist=15):
    """BFS from source using adjacency list. Pure Python."""
    dist = np.full(n, -1, dtype=int)
    dist[source] = 0
    queue = deque([source])
    while queue:
        node = queue.popleft()
        if dist[node] >= max_dist:
            continue
        for neighbor in nbr_list[node]:
            if dist[neighbor] == -1:
                dist[neighbor] = dist[node] + 1
                queue.append(neighbor)
    return dist


def largest_component_size(adjacency, n):
    """Find size of largest connected component via BFS."""
    visited = set()
    max_size = 0
    for start in range(n):
        if start in visited:
            continue
        queue = deque([start])
        component = set()
        while queue:
            node = queue.popleft()
            if node in component:
                continue
            component.add(node)
            visited.add(node)
            for neighbor in range(n):
                if adjacency[node, neighbor] and neighbor not in component:
                    queue.append(neighbor)
        max_size = max(max_size, len(component))
    return max_size


def compute_Cr(corr_matrix, adjacency, n, n_samples=50):
    """Compute radially-averaged correlation C(r) on the k-NN graph."""
    nbr_list = [list(np.where(adjacency[i])[0]) for i in range(n)]
    sample_idx = np.random.choice(n, min(n_samples, n), replace=False)
    max_r = 15

    C_r = np.zeros(max_r + 1)
    counts = np.zeros(max_r + 1)

    for src in sample_idx:
        dist = bfs_distances(nbr_list, int(src), n, max_r)
        for dst in range(n):
            if src == dst:
                continue
            d = dist[dst]
            if 0 < d <= max_r:
                C_r[d] += corr_matrix[src, dst]
                counts[d] += 1

    mask = counts > 0
    C_r[mask] /= counts[mask]
    C_r[~mask] = 0.0

    # Find effective max_r
    effective_max = max_r
    for r in range(max_r, 0, -1):
        if counts[r] > 0:
            effective_max = r
            break

    return C_r[:effective_max + 1]


# ============================================================================
# Fisher diagnostics (same as Phase 2)
# ============================================================================

def build_kernel(C_r, distances, n):
    """Build probability distribution from C(r) kernel and distances."""
    weights = np.zeros(n)
    for u in range(n):
        d = int(distances[u])
        if 0 <= d < len(C_r):
            weights[u] = abs(C_r[d])

    total = np.sum(weights)
    if total < 1e-15:
        return None
    return weights / total


def fisher_diagnostics(C_r, adjacency, neighbors, n, n_samples=30):
    """Compute Fisher diagnostics on the financial graph."""
    nbr_list = [list(np.where(adjacency[i])[0]) for i in range(n)]
    sample_idx = np.random.choice(n, min(n_samples, n), replace=False)

    ranks = []
    etas = []
    prs = []
    gap_ratios = []
    sv2_sv1s = []
    sv_profiles = []

    for v0 in sample_idx:
        v0 = int(v0)
        nbrs = neighbors.get(v0, [])
        if len(nbrs) < 2:
            continue

        dist_v0 = bfs_distances(nbr_list, v0, n, 15)
        p_v0 = build_kernel(C_r, dist_v0, n)
        if p_v0 is None:
            continue

        k = len(nbrs)
        score_vectors = np.zeros((k, n))

        valid = True
        for j, wj in enumerate(nbrs):
            dist_wj = bfs_distances(nbr_list, int(wj), n, 15)
            p_wj = build_kernel(C_r, dist_wj, n)
            if p_wj is None:
                valid = False
                break
            eps = 1e-12
            score_vectors[j] = np.log(p_wj + eps) - np.log(p_v0 + eps)

        if not valid:
            continue

        weighted_scores = score_vectors * np.sqrt(p_v0)[np.newaxis, :]
        FIM = weighted_scores @ weighted_scores.T

        try:
            sv = np.linalg.svd(FIM, compute_uv=False)
        except np.linalg.LinAlgError:
            continue

        sv = np.sort(sv)[::-1]
        if sv[0] < 1e-15:
            continue

        sv_norm = sv / sv[0]

        # Gap-based rank
        rank = 1
        best_gap = 0
        for i in range(len(sv) - 1):
            if sv[i + 1] < 1e-10:
                rank = i + 1
                break
            gap = sv[i] / sv[i + 1]
            if gap > best_gap:
                best_gap = gap
                rank = i + 1

        pr = (np.sum(sv)) ** 2 / (np.sum(sv ** 2) + 1e-15)

        if rank < len(sv):
            eta = sv[rank] / (sv[rank - 1] + 1e-15)
        else:
            eta = 0.0

        sv2_sv1 = float(sv_norm[1]) if len(sv_norm) > 1 else 0.0
        gap_ratio = sv[0] / (sv[1] + 1e-15) if len(sv) > 1 else float('inf')

        ranks.append(rank)
        etas.append(eta)
        prs.append(pr)
        gap_ratios.append(gap_ratio)
        sv2_sv1s.append(sv2_sv1)
        sv_profiles.append(sv_norm[:min(k, 8)].tolist())

    if len(ranks) == 0:
        return None

    return {
        'rank': float(np.mean(ranks)),
        'rank_std': float(np.std(ranks)),
        'eta': float(np.mean(etas)),
        'eta_std': float(np.std(etas)),
        'pr': float(np.mean(prs)),
        'sv2_sv1': float(np.mean(sv2_sv1s)),
        'sv2_sv1_std': float(np.std(sv2_sv1s)),
        'gap_ratio': float(np.mean(gap_ratios)),
        'sv_profile': np.mean(sv_profiles, axis=0).tolist() if sv_profiles else [],
        'n_valid': len(ranks)
    }


# ============================================================================
# Single-window computation
# ============================================================================

def compute_single_window(returns_df, window_end_idx, window_size,
                           k_nn=10, n_fisher_samples=30, n_cr_samples=50):
    """
    For a single rolling window ending at window_end_idx:
    compute correlation matrix -> k-NN graph -> C(r) -> Fisher diagnostics.
    """
    start_idx = max(0, window_end_idx - window_size)
    window_returns = returns_df.iloc[start_idx:window_end_idx]

    # Drop assets with too much missing data
    valid_cols = window_returns.columns[window_returns.notna().sum() >= int(window_size * 0.7)]
    window_returns = window_returns[valid_cols].dropna(axis=0, how='all')

    if window_returns.shape[1] < 40 or len(window_returns) < int(window_size * 0.5):
        return None

    # Fill remaining NaN
    window_returns = window_returns.ffill().bfill().dropna(axis=1)
    if window_returns.shape[1] < 40:
        return None

    n = window_returns.shape[1]

    # Correlation matrix
    corr_matrix = window_returns.corr().values
    np.fill_diagonal(corr_matrix, 1.0)
    corr_matrix = np.nan_to_num(corr_matrix, nan=0.0)

    # k-NN graph
    adjacency, neighbors = build_knn_graph(corr_matrix, k_nn)

    # Check connectivity
    lcs = largest_component_size(adjacency, n)
    if lcs < n * 0.8:
        return None

    # C(r)
    C_r = compute_Cr(corr_matrix, adjacency, n, n_cr_samples)
    if C_r is None or len(C_r) < 3:
        return None

    # Fisher diagnostics
    diagnostics = fisher_diagnostics(C_r, adjacency, neighbors,
                                      n, n_fisher_samples)
    if diagnostics is not None:
        diagnostics['n_assets'] = n
        diagnostics['graph_size'] = lcs
        diagnostics['cr_length'] = len(C_r)
        diagnostics['cr_values'] = C_r.tolist()

    return diagnostics


# ============================================================================
# Main computation loop
# ============================================================================

def main():
    import pandas as pd

    # --- Configuration ---
    tickers = [
        # Tech
        "AAPL", "MSFT", "GOOG", "AMZN", "INTC", "CSCO", "ORCL", "IBM",
        "ADBE", "TXN", "QCOM", "AVGO",
        # Finance
        "JPM", "BAC", "WFC", "GS", "MS", "C", "USB", "PNC", "BK",
        "AXP", "MET", "PRU", "ALL", "TRV",
        # Health
        "JNJ", "UNH", "PFE", "MRK", "LLY", "TMO", "ABT", "MDT", "BMY",
        "AMGN", "GILD", "CVS", "CI",
        # Consumer
        "PG", "KO", "PEP", "WMT", "COST", "HD", "MCD", "NKE", "SBUX", "TGT",
        "CL", "GIS", "SYY",
        # Industrial
        "GE", "CAT", "HON", "MMM", "UPS", "RTX", "LMT", "DE", "EMR", "ITW",
        "GD", "NSC", "CSX", "WM", "ETN",
        # Energy
        "XOM", "CVX", "COP", "SLB", "EOG", "PSX", "VLO", "MPC", "OXY", "HAL",
        # Utilities
        "NEE", "DUK", "SO", "D", "AEP", "EXC", "SRE",
        # Telecom/Media
        "VZ", "T", "CMCSA", "DIS", "NFLX",
    ]
    ref_tickers = ["SPY", "^VIX"]

    windows = [30, 60, 90, 120]
    k_nn = 10
    n_fisher_samples = 30
    n_cr_samples = 50
    step_days = 5  # weekly (every 5 trading days)

    # --- Download data ---
    all_tickers = tickers + ref_tickers
    raw_data = download_data(all_tickers, start='2004-06-01', end='2024-12-31')

    # Build price DataFrame
    price_df = pd.DataFrame(raw_data)
    price_df = price_df.sort_index()
    print(f"\nPrice data: {price_df.shape[0]} days x {price_df.shape[1]} tickers")
    print(f"Date range: {price_df.index[0]} to {price_df.index[-1]}")

    # Extract reference data
    spy_prices = price_df['SPY'].values if 'SPY' in price_df.columns else np.full(len(price_df), np.nan)
    vix_values = price_df['^VIX'].values if '^VIX' in price_df.columns else np.full(len(price_df), np.nan)

    # Filter to analysis tickers only (exclude SPY and VIX)
    analysis_cols = [t for t in tickers if t in price_df.columns]
    analysis_df = price_df[analysis_cols]
    print(f"Analysis tickers: {len(analysis_cols)}")

    # Compute log returns
    returns_df = np.log(analysis_df / analysis_df.shift(1))
    returns_df = returns_df.iloc[1:]  # drop first NaN row

    # --- Main computation loop ---
    max_window = max(windows)
    start_idx = max_window + 30  # buffer after warmup
    total_steps = (len(returns_df) - start_idx) // step_days

    print(f"\nComputation: {total_steps} weekly steps")
    print(f"  Starting from index {start_idx} ({returns_df.index[start_idx].strftime('%Y-%m-%d')})")

    dates = []
    spy_list = []
    vix_list = []
    results = {w: [] for w in windows}

    t0 = time.time()
    step = 0

    for idx in range(start_idx, len(returns_df), step_days):
        date = returns_df.index[idx]
        dates.append(date.strftime('%Y-%m-%d'))

        # Record SPY and VIX
        # Align indices: returns_df is offset by 1 from price_df
        price_idx = idx + 1  # offset for the shift
        if price_idx < len(spy_prices):
            spy_list.append(float(spy_prices[price_idx]) if not np.isnan(spy_prices[price_idx]) else float('nan'))
            vix_list.append(float(vix_values[price_idx]) if not np.isnan(vix_values[price_idx]) else float('nan'))
        else:
            spy_list.append(float('nan'))
            vix_list.append(float('nan'))

        for window in windows:
            try:
                result = compute_single_window(
                    returns_df, idx, window,
                    k_nn=k_nn,
                    n_fisher_samples=n_fisher_samples,
                    n_cr_samples=n_cr_samples
                )
                results[window].append(result)
            except Exception as e:
                results[window].append(None)

        step += 1
        if step % 50 == 0:
            elapsed = time.time() - t0
            r90 = results[90][-1]
            sv_str = f"SV2/SV1={r90['sv2_sv1']:.3f}" if r90 else "None"
            print(f"  [{step}/{total_steps}] {date.strftime('%Y-%m-%d')} "
                  f"({elapsed:.0f}s) {sv_str}")

    elapsed = time.time() - t0
    print(f"\nCompleted {step} steps in {elapsed:.0f}s ({elapsed/60:.1f} min)")

    # --- Save results ---
    output = {
        'dates': dates,
        'vix': vix_list,
        'spy': spy_list,
        'windows': windows,
        'k_nn': k_nn,
        'n_fisher_samples': n_fisher_samples,
        'n_cr_samples': n_cr_samples,
        'tickers': tickers,
        'n_analysis_tickers': len(analysis_cols),
        'data_source': 'yfinance',
    }

    for w in windows:
        cleaned = []
        for r in results[w]:
            if r is None:
                cleaned.append({})
            else:
                cleaned.append(r)
        output[f'results_w{w}'] = cleaned

    outdir = '/tmp/D_eff-Test-1.0/phase3b_results/raw_data'
    os.makedirs(outdir, exist_ok=True)
    outpath = os.path.join(outdir, 'fisher_financial_results.json')

    with open(outpath, 'w') as f:
        json.dump(output, f)

    print(f"\nResults saved to: {outpath}")
    print(f"JSON size: {os.path.getsize(outpath) / 1024:.1f} KB")

    # Summary
    for w in windows:
        valid = sum(1 for r in results[w] if r and 'sv2_sv1' in r)
        print(f"  Window {w}d: {valid}/{len(results[w])} valid results")

    r90_valid = [r for r in results[90] if r and 'sv2_sv1' in r]
    if r90_valid:
        sv_vals = [r['sv2_sv1'] for r in r90_valid]
        print(f"\nSV2/SV1 (90d) summary:")
        print(f"  Mean: {np.mean(sv_vals):.4f}")
        print(f"  Std:  {np.std(sv_vals):.4f}")
        print(f"  Min:  {np.min(sv_vals):.4f}")
        print(f"  Max:  {np.max(sv_vals):.4f}")


if __name__ == '__main__':
    main()
