#!/usr/bin/env python3
"""
DS Phase 3B-3: Fisher Regime Attractor — Multi-Scale Gradient System
Standalone Python implementation for analysis of historical data.

This implements the Fisher Regime Attractor (FRA) system:
1. FRA(t) ∈ [-1, +1] - market position on criticality spectrum
2. Cluster heat vector - gradient scores for dynamic communities  
3. Asset temperature ranking - every asset scored by local Fisher state

All outputs are continuous gradients, no binary signals.
"""

import numpy as np
import pandas as pd
from scipy.sparse.csgraph import shortest_path
from scipy import sparse
from scipy.cluster.vq import kmeans2
from scipy.spatial.distance import pdist, squareform
import json
from datetime import datetime, timedelta
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')

# Configuration
COMPUTE_FREQUENCY = 5  # every 5 trading days (weekly)
WINDOWS = [30, 60, 90]  # multi-scale windows
LOOKBACK_ZSCORE = 252  # 1-year rolling z-score baseline
MIN_CLUSTER_SIZE = 12  # minimum for valid Fisher
MAX_CLUSTERS = 20  # cap on spectral clustering k
N_FISHER_SAMPLES = 20  # per cluster (reduced for speed)
UNIVERSE_SIZE = 200  # top 200 most liquid names

class FisherRegimeAttractor:
    """
    Main class implementing the Fisher Regime Attractor system.
    """
    
    def __init__(self, data_path=None, start_date='2006-01-01', end_date='2024-12-31'):
        """
        Initialize the FRA system.
        
        Args:
            data_path: Path to price data CSV (columns: date, ticker, close)
            start_date: Start date for analysis
            end_date: End date for analysis
        """
        self.data_path = data_path
        self.start_date = pd.to_datetime(start_date)
        self.end_date = pd.to_datetime(end_date)
        
        # Storage
        self.fra_history = []  # List of dicts with date, FRA values
        self.cluster_history = []  # List of dicts with cluster info
        self.asset_scores_history = []  # List of dicts with asset scores
        
        # Rolling histories for z-scoring
        self.sv2sv1_histories = {w: [] for w in WINDOWS}
        self.cluster_heat_histories = {}  # cluster_id -> list
        
        # Results directory
        self.results_dir = Path("phase3b3_results")
        self.results_dir.mkdir(exist_ok=True)
        
        print(f"Fisher Regime Attractor initialized: {start_date} to {end_date}")
        print(f"Windows: {WINDOWS}, Min cluster size: {MIN_CLUSTER_SIZE}")
    
    def load_data(self):
        """
        Load price data from CSV.
        Expected format: CSV with columns [date, ticker, close]
        """
        if not self.data_path:
            print("No data path provided. Using synthetic data for demonstration.")
            return self._create_synthetic_data()
        
        print(f"Loading data from {self.data_path}...")
        df = pd.read_csv(self.data_path)
        df['date'] = pd.to_datetime(df['date'])
        
        # Filter by date range
        mask = (df['date'] >= self.start_date) & (df['date'] <= self.end_date)
        df = df[mask].copy()
        
        # Pivot to wide format (dates × tickers)
        price_data = df.pivot(index='date', columns='ticker', values='close')
        
        # Filter to top N by average dollar volume (simplified: by non-NaN count)
        non_nan_counts = price_data.notna().sum()
        top_tickers = non_nan_counts.nlargest(UNIVERSE_SIZE).index.tolist()
        price_data = price_data[top_tickers]
        
        print(f"Loaded {len(price_data)} dates, {len(top_tickers)} tickers")
        return price_data
    
    def _create_synthetic_data(self):
        """
        Create synthetic price data for demonstration.
        In practice, you would load real price data.
        """
        print("Creating synthetic data for demonstration...")
        
        # Generate dates
        dates = pd.date_range(start=self.start_date, end=self.end_date, freq='B')
        
        # Create synthetic tickers
        n_tickers = min(UNIVERSE_SIZE, 50)  # Smaller for demo
        tickers = [f'TICKER_{i:03d}' for i in range(n_tickers)]
        
        # Base returns with some correlation structure
        np.random.seed(42)
        base_returns = np.random.randn(len(dates), 3) * 0.01
        
        # Each ticker is a weighted combination of base factors + noise
        weights = np.random.rand(n_tickers, 3)
        weights = weights / weights.sum(axis=1, keepdims=True)
        
        price_data = pd.DataFrame(index=dates, columns=tickers)
        
        for i, ticker in enumerate(tickers):
            # Factor returns
            factor_returns = base_returns @ weights[i]
            # Add idiosyncratic noise
            idio_returns = np.random.randn(len(dates)) * 0.005
            total_returns = factor_returns + idio_returns
            
            # Convert to prices (start at 100)
            prices = 100 * np.exp(np.cumsum(total_returns))
            price_data[ticker] = prices
        
        print(f"Created synthetic data: {len(dates)} dates, {n_tickers} tickers")
        return price_data
    
    def compute_regime(self, price_data):
        """
        Main computation loop: weekly FRA + clusters + asset scores.
        """
        print(f"Computing regime from {price_data.index[0].date()} to {price_data.index[-1].date()}")
        
        # Get weekly dates (every 5 trading days)
        all_dates = price_data.index
        weekly_dates = all_dates[::COMPUTE_FREQUENCY]
        
        for i, date in enumerate(weekly_dates):
            if i % 20 == 0:
                print(f"  Processing {date.date()} ({i+1}/{len(weekly_dates)})")
            
            # Get active tickers with sufficient data
            lookback_end = date
            lookback_start = date - timedelta(days=120)  # Enough for 90d window
            
            mask = (price_data.index >= lookback_start) & (price_data.index <= lookback_end)
            recent_data = price_data[mask]
            
            # Filter tickers with enough non-NaN values
            valid_tickers = recent_data.columns[recent_data.notna().sum() >= 60].tolist()
            
            if len(valid_tickers) < 50:
                continue
            
            # Store results for this date
            date_results = {
                'date': date,
                'fra': {},
                'clusters': {},
                'cluster_heats': {},
                'asset_scores': []
            }
            
            # Compute for each window
            fra_values = {}
            for window in WINDOWS:
                try:
                    # Get window data
                    window_end = date
                    window_start = date - timedelta(days=window + 10)  # Buffer
                    window_mask = (price_data.index >= window_start) & (price_data.index <= window_end)
                    window_data = price_data.loc[window_mask, valid_tickers]
                    
                    if len(window_data) < window * 0.7:
                        continue
                    
                    # Compute returns
                    returns = np.log(window_data / window_data.shift(1)).dropna()
                    if len(returns) < window * 0.5:
                        continue
                    
                    # Correlation matrix
                    corr_matrix = returns.corr().values
                    n = corr_matrix.shape[0]
                    np.fill_diagonal(corr_matrix, 1.0)
                    corr_matrix = np.nan_to_num(corr_matrix, nan=0.0)
                    
                    # Market-level Fisher
                    sv2sv1 = self._compute_market_fisher(corr_matrix, n)
                    
                    # Compute FRA
                    fra = self._compute_fra(sv2sv1, window)
                    fra_values[window] = fra
                    
                    # For 90d window, also compute clusters and asset scores
                    if window == 90:
                        # Dynamic clustering
                        clusters = self._dynamic_clusters(corr_matrix, n)
                        
                        # Cluster heat vector
                        cluster_heats = self._compute_cluster_heats(corr_matrix, clusters, n)
                        
                        # Asset temperature ranking
                        asset_scores = self._compute_asset_temperatures(
                            corr_matrix, clusters, cluster_heats, n, 
                            valid_tickers, fra
                        )
                        
                        date_results['clusters'] = clusters
                        date_results['cluster_heats'] = cluster_heats
                        date_results['asset_scores'] = asset_scores
                        
                        # Store cluster info
                        unique_clusters = np.unique(clusters)
                        cluster_sizes = {int(c): int(np.sum(clusters == c)) 
                                        for c in unique_clusters}
                        
                        self.cluster_history.append({
                            'date': str(date.date()),
                            'fra': fra,
                            'n_clusters': len(unique_clusters),
                            'cluster_heats': cluster_heats,
                            'top_10_assets': asset_scores[:10] if asset_scores else [],
                            'cluster_sizes': cluster_sizes,
                        })
                
                except Exception as e:
                    print(f"  Error at {date.date()} window={window}: {str(e)}")
                    continue
            
            # Store FRA values
            date_results['fra'] = fra_values
            
            # Compute derivatives if we have history
            if len(self.fra_history) > 0:
                prev_fra = self.fra_history[-1]['fra'].get(90, 0)
                if 90 in fra_values:
                    velocity = fra_values[90] - prev_fra
                    date_results['fra_velocity'] = velocity
            
            # Compute momentum if we have both 30d and 90d
            if 30 in fra_values and 90 in fra_values:
                momentum = fra_values[30] - fra_values[90]
                date_results['fra_momentum'] = momentum
            
            self.fra_history.append(date_results)
        
        print(f"Computation complete: {len(self.fra_history)} weekly points")
    
    def _compute_market_fisher(self, corr_matrix, n):
        """
        Full-market Fisher SV2/SV1.
        Build k-NN graph (k=10), compute C(r), run Fisher on 20 samples.
        """
        if n < 20:
            return 0.0
        
        # k-NN graph
        k = min(10, n - 1)
        adjacency, neighbors = self._build_knn(corr_matrix, k)
        
        # C(r)
        C_r = self._compute_Cr(corr_matrix, adjacency, n, min(50, n))
        if C_r is None or len(C_r) < 2:
            return 0.0
        
        # Fisher diagnostics (sampled)
        n_samples = min(N_FISHER_SAMPLES, n)
        sample_idx = np.random.choice(n, n_samples, replace=False)
        
        adj_sparse = sparse.csr_matrix(adjacency.astype(float))
        
        sv2_sv1s = []
        for v0 in sample_idx:
            result = self._compute_single_fim(C_r, adj_sparse, neighbors, v0, n)
            if result is not None:
                sv2_sv1s.append(result['sv2_sv1'])
        
        if len(sv2_sv1s) == 0:
            return 0.0
        
        return float(np.mean(sv2_sv1s))
    
    def _compute_fra(self, sv2sv1, window):
        """
        Convert SV2/SV1 to FRA on [-1, +1].
        Uses rolling z-score, sign-flipped, tanh-compressed.
        """
        # Append to rolling history
        hist = self.sv2sv1_histories[window]
        hist.append(sv2sv1)
        
        if len(hist) < 52:  # need 1 year minimum
            return 0.0
        
        # Use trailing LOOKBACK_ZSCORE points
        recent = np.array(hist[-LOOKBACK_ZSCORE:])
        z = (sv2sv1 - np.mean(recent)) / (np.std(recent) + 1e-8)
        fra = float(np.tanh(-z / 2.0))  # sign-flip + compress
        return fra
    
    def _build_knn(self, corr_matrix, k):
        """Build k-NN graph from correlation matrix."""
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
        """Compute C(r) - correlation as function of graph distance."""
        adj_sparse = sparse.csr_matrix(adjacency.astype(float))
        sample_idx = np.random.choice(n, min(n_samples, n), replace=False)
        dist_matrix = shortest_path(adj_sparse, method='D', indices=sample_idx, directed=False)
        
        max_r = int(np.nanmax(dist_matrix[np.isfinite(dist_matrix)]))
        if max_r < 2:
            return None
        max_r = min(max_r, 10)
        
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
        """Compute FIM at a single vertex."""
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
        if total < 1e-15:
            return None
        return weights / total
    
    def _dynamic_clusters(self, corr_matrix, n):
        """
        Spectral clustering with adaptive k.
        Returns array of cluster assignments, length n.
        """
        if n < MIN_CLUSTER_SIZE * 2:
            return np.zeros(n, dtype=int)
        
        # Distance matrix
        dist = np.sqrt(2.0 * np.clip(1.0 - corr_matrix, 0, 2))
        np.fill_diagonal(dist, 0)
        
        # Affinity
        median_d = np.median(dist[dist > 0])
        if median_d == 0:
            return np.zeros(n, dtype=int)
        
        affinity = np.exp(-dist**2 / (2 * median_d**2))
        np.fill_diagonal(affinity, 0)
        
        # Laplacian
        D = np.diag(np.sum(affinity, axis=1))
        L = D - affinity
        D_inv_sqrt = np.diag(1.0 / (np.sqrt(np.diag(D)) + 1e-10))
        L_norm = D_inv_sqrt @ L @ D_inv_sqrt
        
        # Eigendecomposition (smallest eigenvalues)
        try:
            eigenvalues, eigenvectors = np.linalg.eigh(L_norm)
        except:
            return np.zeros(n, dtype=int)
        
        # Eigengap: find k
        max_k = min(MAX_CLUSTERS, n // MIN_CLUSTER_SIZE)
        if max_k < 2:
            return np.zeros(n, dtype=int)
        
        gaps = np.diff(eigenvalues[1:max_k+1])
        k = np.argmax(gaps) + 2  # +2 because we skipped eigenvalue 0
        k = max(k, 2)
        
        # K-means on top-k eigenvectors
        features = eigenvectors[:, 1:k+1]
        # Normalize rows
        norms = np.linalg.norm(features, axis=1, keepdims=True) + 1e-10
        features = features / norms
        
        try:
            _, labels = kmeans2(features, k, minit='++')
        except:
            return np.zeros(n, dtype=int)
        
        # Merge small clusters
        labels = self._merge_small_clusters(labels, corr_matrix, MIN_CLUSTER_SIZE)
        
        return labels
    
    def _merge_small_clusters(self, labels, corr_matrix, min_size):
        """Merge clusters below min_size into nearest neighbor."""
        unique, counts = np.unique(labels, return_counts=True)
        
        for c, count in zip(unique, counts):
            if count < min_size:
                # Find nearest cluster by mean correlation
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
        
        # Relabel consecutively
        unique_new = np.unique(labels)
        mapping = {old: new for new, old in enumerate(unique_new)}
        return np.array([mapping[l] for l in labels])
    
    def _compute_cluster_heats(self, corr_matrix, labels, n):
        """
        Per-cluster Fisher SV2/SV1, z-scored against cluster history.
        Returns dict: cluster_id → heat (continuous, signed).
        """
        heats = {}
        unique_clusters = np.unique(labels)
        
        for c in unique_clusters:
            members = np.where(labels == c)[0]
            if len(members) < MIN_CLUSTER_SIZE:
                heats[int(c)] = 0.0
                continue
            
            # Sub-correlation matrix
            sub_corr = corr_matrix[np.ix_(members, members)]
            n_c = len(members)
            k_nn = max(4, n_c // 3)
            
            # Fisher on sub-graph
            sv2sv1 = self._compute_subgraph_fisher(sub_corr, n_c, k_nn)
            
            # Z-score against this cluster's history
            key = f'cluster_{c}'
            if key not in self.cluster_heat_histories:
                self.cluster_heat_histories[key] = []
            hist = self.cluster_heat_histories[key]
            hist.append(sv2sv1)
            
            if len(hist) < 26:  # need 6 months
                heats[int(c)] = 0.0
            else:
                recent = np.array(hist[-52:])
                z = (sv2sv1 - np.mean(recent)) / (np.std(recent) + 1e-8)
                heats[int(c)] = float(z)
        
        return heats
    
    def _compute_subgraph_fisher(self, corr_matrix, n, k):
        """Compute Fisher SV2/SV1 for a subgraph."""
        if n < 7:
            return 0.0
        
        adjacency, neighbors = self._build_knn(corr_matrix, k)
        C_r = self._compute_Cr(corr_matrix, adjacency, n, min(20, n))
        if C_r is None or len(C_r) < 2:
            return 0.0
        
        n_samples = min(N_FISHER_SAMPLES, n)
        sample_idx = np.random.choice(n, n_samples, replace=False)
        adj_sparse = sparse.csr_matrix(adjacency.astype(float))
        
        sv2_sv1s = []
        for v0 in sample_idx:
            result = self._compute_single_fim(C_r, adj_sparse, neighbors, v0, n)
            if result is not None:
                sv2_sv1s.append(result['sv2_sv1'])
        
        if len(sv2_sv1s) == 0:
            return 0.0
        
        return float(np.mean(sv2_sv1s))
    
    def _compute_asset_temperatures(self, corr_matrix, labels, cluster_heats, 
                                     n, ticker_names, fra):
        """
        Per-asset Fisher temperature: deviation from cluster mean.
        Returns list of (asset_name, composite_score, components) dicts.
        """
        asset_scores = []
        unique_clusters = np.unique(labels)
        
        for c in unique_clusters:
            members = np.where(labels == c)[0]
            if len(members) < MIN_CLUSTER_SIZE:
                continue
            
            sub_corr = corr_matrix[np.ix_(members, members)]
            n_c = len(members)
            k_nn = max(4, n_c // 3)
            
            # Per-asset SV2/SV1
            asset_sv2sv1 = self._compute_per_asset_fisher(sub_corr, n_c, k_nn)
            
            if not asset_sv2sv1:
                continue
            
            mean_sv = np.mean(list(asset_sv2sv1.values()))
            cluster_heat = cluster_heats.get(int(c), 0.0)
            
            for local_idx, sv in asset_sv2sv1.items():
                global_idx = members[local_idx]
                if global_idx < len(ticker_names):
                    name = str(ticker_names[global_idx])
                    micro = sv - mean_sv  # deviation from cluster mean
                    
                    # Composite score
                    composite = (0.4 * fra 
                               + 0.4 * cluster_heat 
                               + 0.2 * micro)
                    
                    asset_scores.append({
                        'name': name,
                        'composite': float(composite),
                        'macro_fra': float(fra),
                        'meso_heat': float(cluster_heat),
                        'micro_temp': float(micro),
                        'cluster': int(c),
                        'sv2_sv1': float(sv)
                    })
        
        # Sort by |composite| descending
        asset_scores.sort(key=lambda x: abs(x['composite']), reverse=True)
        return asset_scores
    
    def _compute_per_asset_fisher(self, corr_matrix, n, k):
        """
        Compute Fisher SV2/SV1 for each asset in a subgraph.
        Returns dict: local_index -> sv2_sv1.
        """
        if n < 7:
            return {}
        
        adjacency, neighbors = self._build_knn(corr_matrix, k)
        C_r = self._compute_Cr(corr_matrix, adjacency, n, min(20, n))
        if C_r is None or len(C_r) < 2:
            return {}
        
        adj_sparse = sparse.csr_matrix(adjacency.astype(float))
        
        asset_sv2sv1 = {}
        # Sample up to 20 assets per cluster for speed
        n_samples = min(20, n)
        sample_idx = np.random.choice(n, n_samples, replace=False)
        
        for v0 in sample_idx:
            result = self._compute_single_fim(C_r, adj_sparse, neighbors, v0, n)
            if result is not None:
                asset_sv2sv1[v0] = result['sv2_sv1']
        
        return asset_sv2sv1
    
    def save_results(self):
        """Save all results to JSON files."""
        print("Saving results...")
        
        # FRA history
        fra_data = []
        for result in self.fra_history:
            fra_data.append({
                'date': str(result['date'].date()),
                'fra_30d': result['fra'].get(30, 0.0),
                'fra_60d': result['fra'].get(60, 0.0),
                'fra_90d': result['fra'].get(90, 0.0),
                'fra_velocity': result.get('fra_velocity', 0.0),
                'fra_momentum': result.get('fra_momentum', 0.0),
            })
        
        with open(self.results_dir / 'fra_history.json', 'w') as f:
            json.dump(fra_data, f, indent=2)
        
        # Cluster history
        with open(self.results_dir / 'cluster_history.json', 'w') as f:
            json.dump(self.cluster_history, f, indent=2, default=str)
        
        print(f"Results saved to {self.results_dir}/")
        return fra_data, self.cluster_history
    
    def generate_report(self):
        """Generate a summary report of the analysis."""
        report = []
        report.append("# DS Phase 3B-3: Fisher Regime Attractor Results")
        report.append(f"**Generated:** {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
        report.append(f"**Period:** {self.start_date.date()} to {self.end_date.date()}")
        report.append(f"**Data points:** {len(self.fra_history)} weekly computations")
        report.append("")
        
        # FRA statistics
        if self.fra_history:
            fra_90d = [r['fra'].get(90, 0) for r in self.fra_history if 90 in r['fra']]
            if fra_90d:
                report.append("## FRA Statistics (90-day window)")
                report.append(f"- Mean: {np.mean(fra_90d):.3f}")
                report.append(f"- Std: {np.std(fra_90d):.3f}")
                report.append(f"- Min: {np.min(fra_90d):.3f}")
                report.append(f"- Max: {np.max(fra_90d):.3f}")
                report.append(f"- Time in negative regime (FRA < -0.3): {100 * np.mean(np.array(fra_90d) < -0.3):.1f}%")
                report.append(f"- Time in positive regime (FRA > 0.3): {100 * np.mean(np.array(fra_90d) > 0.3):.1f}%")
                report.append("")
        
        # Cluster statistics
        if self.cluster_history:
            n_clusters = [h['n_clusters'] for h in self.cluster_history]
            report.append("## Cluster Statistics")
            report.append(f"- Mean clusters per week: {np.mean(n_clusters):.1f}")
            report.append(f"- Std clusters: {np.std(n_clusters):.1f}")
            report.append(f"- Min clusters: {np.min(n_clusters)}")
            report.append(f"- Max clusters: {np.max(n_clusters)}")
            report.append("")
        
        # Top assets across history
        if self.cluster_history:
            all_top_assets = []
            for h in self.cluster_history:
                all_top_assets.extend(h.get('top_10_assets', []))
            
            if all_top_assets:
                # Count frequency in top 10
                asset_counts = {}
                for asset in all_top_assets:
                    name = asset.get('name', 'Unknown')
                    asset_counts[name] = asset_counts.get(name, 0) + 1
                
                top_10_frequent = sorted(asset_counts.items(), key=lambda x: x[1], reverse=True)[:10]
                report.append("## Most Frequent Top Assets")
                report.append("| Asset | Weeks in Top 10 |")
                report.append("|-------|-----------------|")
                for name, count in top_10_frequent:
                    report.append(f"| {name} | {count} |")
                report.append("")
        
        report.append("## Files Generated")
        report.append("- `fra_history.json`: FRA values and derivatives")
        report.append("- `cluster_history.json`: Cluster assignments and heats")
        report.append("- `analysis_report.md`: This report")
        
        report_path = self.results_dir / "analysis_report.md"
        with open(report_path, 'w') as f:
            f.write('\n'.join(report))
        
        print(f"Report saved to {report_path}")
        return report_path

def main():
    """
    Main function to run the Fisher Regime Attractor analysis.
    """
    print("=" * 60)
    print("DS Phase 3B-3: Fisher Regime Attractor")
    print("Multi-Scale Gradient System")
    print("=" * 60)
    
    # Initialize the system
    fra_system = FisherRegimeAttractor(
        data_path=None,  # Set to your price data CSV path
        start_date='2006-01-01',
        end_date='2024-12-31'
    )
    
    # Load data
    price_data = fra_system.load_data()
    
    # Compute regime
    fra_system.compute_regime(price_data)
    
    # Save results
    fra_system.save_results()
    
    # Generate report
    fra_system.generate_report()
    
    print("\n" + "=" * 60)
    print("ANALYSIS COMPLETE")
    print("=" * 60)
    print("\nNext steps:")
    print("1. Use the generated JSON files for visualization")
    print("2. Run the post-processing script to create figures")
    print("3. Test the pre-registered predictions")
    print("\nTo use with real data:")
    print("1. Export price data as CSV with columns: date, ticker, close")
    print("2. Update data_path in main() function")
    print("3. Run: python FisherRegimeAttractor.py")

if __name__ == "__main__":
    main()