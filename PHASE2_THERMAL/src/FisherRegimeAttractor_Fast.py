#!/usr/bin/env python3
"""
DS Phase 3B-3: Fisher Regime Attractor — Fast Demo Version
Optimized for quick execution with synthetic data.
"""

import numpy as np
import pandas as pd
from scipy.sparse.csgraph import shortest_path
from scipy import sparse
from scipy.cluster.vq import kmeans2
import json
from datetime import datetime, timedelta
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')

# Configuration - OPTIMIZED FOR SPEED
COMPUTE_FREQUENCY = 20  # every 20 trading days (monthly instead of weekly)
WINDOWS = [90]  # Only 90d window for speed
LOOKBACK_ZSCORE = 252
MIN_CLUSTER_SIZE = 8  # Reduced from 12
MAX_CLUSTERS = 15  # Reduced from 20
N_FISHER_SAMPLES = 10  # Reduced from 20
UNIVERSE_SIZE = 50  # Reduced from 200 for demo

class FisherRegimeAttractor:
    """Fast demo version of Fisher Regime Attractor."""
    
    def __init__(self, start_date='2006-01-01', end_date='2024-12-31'):
        self.start_date = pd.to_datetime(start_date)
        self.end_date = pd.to_datetime(end_date)
        self.fra_history = []
        self.cluster_history = []
        self.sv2sv1_histories = {w: [] for w in WINDOWS}
        self.cluster_heat_histories = {}
        self.results_dir = Path("phase3b3_results")
        self.results_dir.mkdir(exist_ok=True)
        print(f"FRA initialized: {start_date} to {end_date}")
    
    def create_synthetic_data(self):
        """Create synthetic price data with crisis patterns."""
        print("Creating synthetic data with crisis patterns...")
        
        dates = pd.date_range(start=self.start_date, end=self.end_date, freq='B')
        n_tickers = UNIVERSE_SIZE
        tickers = [f'TICK{i:02d}' for i in range(n_tickers)]
        
        np.random.seed(42)
        price_data = pd.DataFrame(index=dates, columns=tickers)
        
        # Create base factors with crisis events
        t = np.arange(len(dates))
        
        # Factor 1: Market factor with crisis dips
        factor1 = np.sin(t / 500) * 0.02
        # Add crisis events (2008, 2020)
        crisis_2008 = np.exp(-((t - 500)**2) / 50000) * -0.15
        crisis_2020 = np.exp(-((t - 3700)**2) / 50000) * -0.12
        factor1 = factor1 + crisis_2008 + crisis_2020
        
        # Factor 2: Sector rotation
        factor2 = np.sin(t / 300) * 0.015
        
        # Factor 3: Volatility regime
        factor3 = np.sin(t / 200) * 0.01
        
        # Generate prices
        for i, ticker in enumerate(tickers):
            weights = np.random.rand(3)
            weights = weights / weights.sum()
            
            factor_returns = (weights[0] * factor1 + 
                            weights[1] * factor2 + 
                            weights[2] * factor3)
            
            idio_returns = np.random.randn(len(dates)) * 0.005
            total_returns = factor_returns + idio_returns
            
            prices = 100 * np.exp(np.cumsum(total_returns))
            price_data[ticker] = prices
        
        print(f"Created {len(dates)} dates, {n_tickers} tickers")
        return price_data
    
    def compute_regime(self, price_data):
        """Main computation loop."""
        print(f"Computing regime...")
        
        all_dates = price_data.index
        weekly_dates = all_dates[::COMPUTE_FREQUENCY]
        
        for i, date in enumerate(weekly_dates):
            if i % 5 == 0:
                print(f"  {date.date()} ({i+1}/{len(weekly_dates)})")
            
            lookback_start = date - timedelta(days=120)
            mask = (price_data.index >= lookback_start) & (price_data.index <= date)
            recent_data = price_data[mask]
            
            valid_tickers = recent_data.columns[recent_data.notna().sum() >= 60].tolist()
            if len(valid_tickers) < 30:
                continue
            
            for window in WINDOWS:
                try:
                    window_start = date - timedelta(days=window + 10)
                    window_mask = (price_data.index >= window_start) & (price_data.index <= date)
                    window_data = price_data.loc[window_mask, valid_tickers]
                    
                    if len(window_data) < window * 0.7:
                        continue
                    
                    returns = np.log(window_data / window_data.shift(1)).dropna()
                    if len(returns) < window * 0.5:
                        continue
                    
                    corr_matrix = returns.corr().values
                    n = corr_matrix.shape[0]
                    np.fill_diagonal(corr_matrix, 1.0)
                    corr_matrix = np.nan_to_num(corr_matrix, nan=0.0)
                    
                    # Market Fisher
                    sv2sv1 = self._compute_market_fisher(corr_matrix, n)
                    fra = self._compute_fra(sv2sv1, window)
                    
                    # Clustering
                    clusters = self._dynamic_clusters(corr_matrix, n)
                    cluster_heats = self._compute_cluster_heats(corr_matrix, clusters, n)
                    asset_scores = self._compute_asset_temperatures(
                        corr_matrix, clusters, cluster_heats, n, valid_tickers, fra
                    )
                    
                    # Store
                    unique_clusters = np.unique(clusters)
                    self.cluster_history.append({
                        'date': str(date.date()),
                        'fra': fra,
                        'n_clusters': len(unique_clusters),
                        'cluster_heats': cluster_heats,
                        'top_10_assets': asset_scores[:10],
                        'cluster_sizes': {int(c): int(np.sum(clusters == c)) 
                                        for c in unique_clusters},
                    })
                    
                    self.fra_history.append({
                        'date': date,
                        'fra': {window: fra}
                    })
                
                except Exception as e:
                    pass
        
        print(f"Complete: {len(self.fra_history)} points")
    
    def _compute_market_fisher(self, corr_matrix, n):
        """Market-level Fisher SV2/SV1."""
        if n < 15:
            return 0.0
        
        k = min(8, n - 1)
        adjacency, neighbors = self._build_knn(corr_matrix, k)
        C_r = self._compute_Cr(corr_matrix, adjacency, n, min(30, n))
        if C_r is None or len(C_r) < 2:
            return 0.0
        
        n_samples = min(N_FISHER_SAMPLES, n)
        sample_idx = np.random.choice(n, n_samples, replace=False)
        adj_sparse = sparse.csr_matrix(adjacency.astype(float))
        
        sv2_sv1s = []
        for v0 in sample_idx:
            result = self._compute_single_fim(C_r, adj_sparse, neighbors, v0, n)
            if result:
                sv2_sv1s.append(result['sv2_sv1'])
        
        return float(np.mean(sv2_sv1s)) if sv2_sv1s else 0.0
    
    def _compute_fra(self, sv2sv1, window):
        """Convert SV2/SV1 to FRA [-1, +1]."""
        hist = self.sv2sv1_histories[window]
        hist.append(sv2sv1)
        
        if len(hist) < 26:
            return 0.0
        
        recent = np.array(hist[-LOOKBACK_ZSCORE:])
        z = (sv2sv1 - np.mean(recent)) / (np.std(recent) + 1e-8)
        return float(np.tanh(-z / 2.0))
    
    def _build_knn(self, corr_matrix, k):
        """Build k-NN graph."""
        n = corr_matrix.shape[0]
        adjacency = np.zeros((n, n), dtype=bool)
        for i in range(n):
            row = corr_matrix[i].copy()
            row[i] = -np.inf
            top_k = np.argsort(row)[-k:]
            for j in top_k:
                adjacency[i, j] = True
                adjacency[j, i] = True
        neighbors = {i: list(np.where(adjacency[i])[0]) for i in range(n)}
        return adjacency, neighbors
    
    def _compute_Cr(self, corr_matrix, adjacency, n, n_samples):
        """Compute C(r)."""
        adj_sparse = sparse.csr_matrix(adjacency.astype(float))
        sample_idx = np.random.choice(n, min(n_samples, n), replace=False)
        dist_matrix = shortest_path(adj_sparse, method='D', indices=sample_idx, directed=False)
        
        max_r = int(np.nanmax(dist_matrix[np.isfinite(dist_matrix)]))
        if max_r < 2:
            return None
        max_r = min(max_r, 8)
        
        C_r = np.zeros(max_r + 1)
        counts = np.zeros(max_r + 1)
        for idx_i, src in enumerate(sample_idx):
            for dst in range(n):
                if src == dst:
                    continue
                d = dist_matrix[idx_i, dst]
                if np.isfinite(d) and int(d) <= max_r:
                    C_r[int(d)] += corr_matrix[src, dst]
                    counts[int(d)] += 1
        mask = counts > 0
        C_r[mask] /= counts[mask]
        return C_r
    
    def _compute_single_fim(self, C_r, adj_sparse, neighbors, v0, n):
        """Compute FIM at vertex."""
        nbrs = neighbors.get(v0, [])
        if len(nbrs) < 2:
            return None
        
        sources = [v0] + list(nbrs)
        dists = shortest_path(adj_sparse, method='D', indices=sources, directed=False)
        
        p_v0 = self._build_kernel(C_r, dists[0], n)
        if p_v0 is None:
            return None
        
        k = len(nbrs)
        score_vectors = np.zeros((k, n))
        eps = 1e-12
        
        for j in range(k):
            p_wj = self._build_kernel(C_r, dists[j + 1], n)
            if p_wj is None:
                return None
            score_vectors[j] = np.log(p_wj + eps) - np.log(p_v0 + eps)
        
        weighted = score_vectors * np.sqrt(p_v0)[np.newaxis, :]
        FIM = weighted @ weighted.T
        
        sv = np.linalg.svd(FIM, compute_uv=False)
        sv = np.sort(sv)[::-1]
        sv_norm = sv / (sv[0] + 1e-15)
        
        sv2_sv1 = float(sv_norm[1]) if len(sv_norm) > 1 else 0.0
        return {'sv2_sv1': sv2_sv1}
    
    def _build_kernel(self, C_r, distances, n):
        """Build kernel from C(r)."""
        weights = np.zeros(n)
        for u in range(n):
            d = int(distances[u])
            if d < len(C_r) and np.isfinite(distances[u]):
                weights[u] = abs(C_r[d])
        total = np.sum(weights)
        return weights / total if total > 1e-15 else None
    
    def _dynamic_clusters(self, corr_matrix, n):
        """Spectral clustering with adaptive k."""
        if n < MIN_CLUSTER_SIZE * 2:
            return np.zeros(n, dtype=int)
        
        dist = np.sqrt(2.0 * np.clip(1.0 - corr_matrix, 0, 2))
        np.fill_diagonal(dist, 0)
        
        median_d = np.median(dist[dist > 0])
        if median_d == 0:
            return np.zeros(n, dtype=int)
        
        affinity = np.exp(-dist**2 / (2 * median_d**2))
        np.fill_diagonal(affinity, 0)
        
        D = np.diag(np.sum(affinity, axis=1))
        L = D - affinity
        D_inv_sqrt = np.diag(1.0 / (np.sqrt(np.diag(D)) + 1e-10))
        L_norm = D_inv_sqrt @ L @ D_inv_sqrt
        
        try:
            eigenvalues, eigenvectors = np.linalg.eigh(L_norm)
        except:
            return np.zeros(n, dtype=int)
        
        max_k = min(MAX_CLUSTERS, n // MIN_CLUSTER_SIZE)
        if max_k < 2:
            return np.zeros(n, dtype=int)
        
        gaps = np.diff(eigenvalues[1:max_k+1])
        k = np.argmax(gaps) + 2
        k = max(k, 2)
        
        features = eigenvectors[:, 1:k+1]
        norms = np.linalg.norm(features, axis=1, keepdims=True) + 1e-10
        features = features / norms
        
        try:
            _, labels = kmeans2(features, k, minit='++')
        except:
            return np.zeros(n, dtype=int)
        
        return self._merge_small_clusters(labels, corr_matrix, MIN_CLUSTER_SIZE)
    
    def _merge_small_clusters(self, labels, corr_matrix, min_size):
        """Merge small clusters."""
        unique, counts = np.unique(labels, return_counts=True)
        
        for c, count in zip(unique, counts):
            if count < min_size:
                members = np.where(labels == c)[0]
                best_target = -1
                best_corr = -1
                for other_c in unique:
                    if other_c == c:
                        continue
                    other_members = np.where(labels == other_c)[0]
                    if len(other_members) < min_size:
                        continue
                    mean_corr = np.mean(corr_matrix[np.ix_(members, other_members)])
                    if mean_corr > best_corr:
                        best_corr = mean_corr
                        best_target = other_c
                
                if best_target >= 0:
                    labels[members] = best_target
        
        unique_new = np.unique(labels)
        mapping = {old: new for new, old in enumerate(unique_new)}
        return np.array([mapping[l] for l in labels])
    
    def _compute_cluster_heats(self, corr_matrix, labels, n):
        """Per-cluster Fisher temperature."""
        heats = {}
        unique_clusters = np.unique(labels)
        
        for c in unique_clusters:
            members = np.where(labels == c)[0]
            if len(members) < MIN_CLUSTER_SIZE:
                heats[int(c)] = 0.0
                continue
            
            sub_corr = corr_matrix[np.ix_(members, members)]
            n_c = len(members)
            k_nn = max(3, n_c // 3)
            
            sv2sv1 = self._compute_subgraph_fisher(sub_corr, n_c, k_nn)
            
            key = f'cluster_{c}'
            if key not in self.cluster_heat_histories:
                self.cluster_heat_histories[key] = []
            hist = self.cluster_heat_histories[key]
            hist.append(sv2sv1)
            
            if len(hist) < 13:
                heats[int(c)] = 0.0
            else:
                recent = np.array(hist[-26:])
                z = (sv2sv1 - np.mean(recent)) / (np.std(recent) + 1e-8)
                heats[int(c)] = float(z)
        
        return heats
    
    def _compute_subgraph_fisher(self, corr_matrix, n, k):
        """Fisher for subgraph."""
        if n < 5:
            return 0.0
        
        adjacency, neighbors = self._build_knn(corr_matrix, k)
        C_r = self._compute_Cr(corr_matrix, adjacency, n, min(15, n))
        if C_r is None or len(C_r) < 2:
            return 0.0
        
        n_samples = min(8, n)
        sample_idx = np.random.choice(n, n_samples, replace=False)
        adj_sparse = sparse.csr_matrix(adjacency.astype(float))
        
        sv2_sv1s = []
        for v0 in sample_idx:
            result = self._compute_single_fim(C_r, adj_sparse, neighbors, v0, n)
            if result:
                sv2_sv1s.append(result['sv2_sv1'])
        
        return float(np.mean(sv2_sv1s)) if sv2_sv1s else 0.0
    
    def _compute_asset_temperatures(self, corr_matrix, labels, cluster_heats, 
                                     n, ticker_names, fra):
        """Per-asset Fisher temperature."""
        asset_scores = []
        unique_clusters = np.unique(labels)
        
        for c in unique_clusters:
            members = np.where(labels == c)[0]
            if len(members) < MIN_CLUSTER_SIZE:
                continue
            
            sub_corr = corr_matrix[np.ix_(members, members)]
            n_c = len(members)
            k_nn = max(3, n_c // 3)
            
            asset_sv2sv1 = self._compute_per_asset_fisher(sub_corr, n_c, k_nn)
            
            if not asset_sv2sv1:
                continue
            
            mean_sv = np.mean(list(asset_sv2sv1.values()))
            cluster_heat = cluster_heats.get(int(c), 0.0)
            
            for local_idx, sv in asset_sv2sv1.items():
                global_idx = members[local_idx]
                if global_idx < len(ticker_names):
                    name = str(ticker_names[global_idx])
                    micro = sv - mean_sv
                    
                    composite = (0.4 * fra + 0.4 * cluster_heat + 0.2 * micro)
                    
                    asset_scores.append({
                        'name': name,
                        'composite': float(composite),
                        'macro_fra': float(fra),
                        'meso_heat': float(cluster_heat),
                        'micro_temp': float(micro),
                        'cluster': int(c),
                        'sv2_sv1': float(sv)
                    })
        
        asset_scores.sort(key=lambda x: abs(x['composite']), reverse=True)
        return asset_scores
    
    def _compute_per_asset_fisher(self, corr_matrix, n, k):
        """Fisher for each asset."""
        if n < 5:
            return {}
        
        adjacency, neighbors = self._build_knn(corr_matrix, k)
        C_r = self._compute_Cr(corr_matrix, adjacency, n, min(15, n))
        if C_r is None or len(C_r) < 2:
            return {}
        
        adj_sparse = sparse.csr_matrix(adjacency.astype(float))
        
        asset_sv2sv1 = {}
        n_samples = min(10, n)
        sample_idx = np.random.choice(n, n_samples, replace=False)
        
        for v0 in sample_idx:
            result = self._compute_single_fim(C_r, adj_sparse, neighbors, v0, n)
            if result:
                asset_sv2sv1[v0] = result['sv2_sv1']
        
        return asset_sv2sv1
    
    def save_results(self):
        """Save results."""
        print("Saving results...")
        
        fra_data = []
        for result in self.fra_history:
            fra_data.append({
                'date': str(result['date'].date()),
                'fra_90d': result['fra'].get(90, 0.0),
            })
        
        with open(self.results_dir / 'fra_history.json', 'w') as f:
            json.dump(fra_data, f, indent=2)
        
        with open(self.results_dir / 'cluster_history.json', 'w') as f:
            json.dump(self.cluster_history, f, indent=2, default=str)
        
        print(f"Results saved to {self.results_dir}/")
        return fra_data, self.cluster_history
    
    def generate_report(self):
        """Generate summary report."""
        report = []
        report.append("# DS Phase 3B-3: Fisher Regime Attractor Results")
        report.append(f"**Generated:** {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
        report.append(f"**Period:** {self.start_date.date()} to {self.end_date.date()}")
        report.append(f"**Data points:** {len(self.fra_history)} computations")
        report.append("")
        
        if self.fra_history:
            fra_90d = [r['fra'].get(90, 0) for r in self.fra_history if 90 in r['fra']]
            if fra_90d:
                report.append("## FRA Statistics (90-day window)")
                report.append(f"- Mean: {np.mean(fra_90d):.3f}")
                report.append(f"- Std: {np.std(fra_90d):.3f}")
                report.append(f"- Min: {np.min(fra_90d):.3f}")
                report.append(f"- Max: {np.max(fra_90d):.3f}")
                report.append(f"- Negative regime (FRA < -0.3): {100 * np.mean(np.array(fra_90d) < -0.3):.1f}%")
                report.append("")
        
        if self.cluster_history:
            n_clusters = [h['n_clusters'] for h in self.cluster_history]
            report.append("## Cluster Statistics")
            report.append(f"- Mean clusters: {np.mean(n_clusters):.1f}")
            report.append(f"- Std: {np.std(n_clusters):.1f}")
            report.append(f"- Range: {np.min(n_clusters)} to {np.max(n_clusters)}")
            report.append("")
        
        report.append("## Files Generated")
        report.append("- `fra_history.json`: FRA values")
        report.append("- `cluster_history.json`: Cluster data")
        report.append("- `analysis_report.md`: This report")
        
        report_path = self.results_dir / "analysis_report.md"
        with open(report_path, 'w') as f:
            f.write('\n'.join(report))
        
        print(f"Report saved to {report_path}")
        return report_path

def main():
    print("=" * 60)
    print("DS Phase 3B-3: Fisher Regime Attractor (Fast Demo)")
    print("=" * 60)
    
    fra_system = FisherRegimeAttractor()
    price_data = fra_system.create_synthetic_data()
    fra_system.compute_regime(price_data)
    fra_system.save_results()
    fra_system.generate_report()
    
    print("\n" + "=" * 60)
    print("ANALYSIS COMPLETE")
    print("=" * 60)

if __name__ == "__main__":
    main()