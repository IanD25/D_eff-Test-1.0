"""
DS Phase 3B-3: Fisher Regime Attractor — QuantConnect LEAN Algorithm
Multi-Scale Gradient System for Market Regime Analysis

This is the QuantConnect-compatible version that inherits from QCAlgorithm.
"""

from AlgorithmImports import *
import numpy as np
from scipy.sparse.csgraph import shortest_path
from scipy import sparse
from scipy.cluster.vq import kmeans2
import json

# Configuration
WINDOWS = [90]  # Use 90d only for speed (can add 30d, 60d later)
LOOKBACK_ZSCORE = 252
MIN_CLUSTER_SIZE = 12
MAX_CLUSTERS = 20
N_FISHER_SAMPLES = 20

class FisherRegimeAttractor(QCAlgorithm):
    """
    Fisher Regime Attractor algorithm for QuantConnect LEAN.
    Computes market regime, cluster dynamics, and asset temperatures.
    """
    
    def initialize(self):
        """Initialize the algorithm."""
        self.set_start_date(2006, 1, 1)
        self.set_end_date(2024, 12, 31)
        self.set_cash(100000)
        
        # Add SPY for reference
        self.spy = self.add_equity("SPY", Resolution.DAILY).symbol
        
        # Add VIX
        self.vix = self.add_data(CBOE, "VIX", Resolution.DAILY).symbol
        
        # Universe: top 200 by dollar volume
        self.add_universe(self.coarse_selection)
        self.universe_settings.resolution = Resolution.DAILY
        
        # Schedule weekly computation (every Monday)
        self.schedule.on(
            self.date_rules.every(DayOfWeek.MONDAY),
            self.time_rules.after_market_open("SPY", 30),
            self.compute_regime
        )
        
        # Storage for results
        self.fra_history = []
        self.cluster_history = []
        self.sv2sv1_history = []
        
        # LEAN Charts
        fra_chart = Chart("Fisher Regime Attractor")
        fra_chart.add_series(Series("FRA (90d)", SeriesType.LINE, 0))
        fra_chart.add_series(Series("FRA Velocity", SeriesType.LINE, 1))
        self.add_chart(fra_chart)
        
        cluster_chart = Chart("Cluster Dynamics")
        cluster_chart.add_series(Series("N Clusters", SeriesType.LINE, 0))
        cluster_chart.add_series(Series("Max |Heat|", SeriesType.LINE, 1))
        self.add_chart(cluster_chart)
        
        self.set_warm_up(timedelta(days=400))
        
        self.log("FisherRegimeAttractor initialized")
    
    def coarse_selection(self, coarse):
        """Select top 200 by dollar volume."""
        if self.time.month % 3 != 1 or self.time.day > 7:
            return Universe.UNCHANGED
        
        sorted_by_volume = sorted(
            [x for x in coarse if x.has_fundamental_data and x.price > 5],
            key=lambda x: x.dollar_volume, reverse=True
        )
        return [x.symbol for x in sorted_by_volume[:200]]
    
    def compute_regime(self):
        """Weekly computation: FRA + clusters + asset scores."""
        if self.is_warming_up:
            return
        
        date = self.time
        
        # Get active securities
        active = [s for s in self.active_securities.keys() 
                 if s != self.spy and s != self.vix
                 and self.securities[s].price > 0]
        
        if len(active) < 50:
            self.log(f"{date.date()}: Insufficient active securities ({len(active)})")
            return
        
        try:
            # Compute for 90d window
            window = 90
            
            # Get price history
            history = self.history(active, window, Resolution.DAILY)
            if history.empty:
                return
            
            closes = history['close'].unstack(level=0)
            closes = closes.dropna(axis=1, thresh=int(window * 0.8))
            if closes.shape[1] < 50:
                return
            
            # Compute returns
            returns = np.log(closes / closes.shift(1)).dropna()
            if len(returns) < window * 0.7:
                return
            
            # Correlation matrix
            corr_matrix = returns.corr().values
            n = corr_matrix.shape[0]
            np.fill_diagonal(corr_matrix, 1.0)
            corr_matrix = np.nan_to_num(corr_matrix, nan=0.0)
            
            ticker_names = list(returns.columns)
            
            # --- Market-level Fisher ---
            sv2sv1 = self._compute_market_fisher(corr_matrix, n)
            
            # --- FRA Computation ---
            fra = self._compute_fra(sv2sv1, window)
            
            # --- Dynamic Clustering ---
            clusters = self._dynamic_clusters(corr_matrix, n)
            
            # --- Cluster Heat Vector ---
            cluster_heats = self._compute_cluster_heats(corr_matrix, clusters, n)
            
            # --- Asset Temperature Ranking ---
            asset_scores = self._compute_asset_temperatures(
                corr_matrix, clusters, cluster_heats, n, ticker_names, fra
            )
            
            # --- Store and Plot ---
            unique_clusters = np.unique(clusters)
            n_clusters = len(unique_clusters)
            max_heat = max(abs(h) for h in cluster_heats.values()) if cluster_heats else 0
            
            # Plot FRA
            self.plot("Fisher Regime Attractor", "FRA (90d)", fra)
            
            # Plot velocity
            if len(self.fra_history) > 0:
                velocity = fra - self.fra_history[-1]
                self.plot("Fisher Regime Attractor", "FRA Velocity", velocity)
            
            # Plot clusters
            self.plot("Cluster Dynamics", "N Clusters", n_clusters)
            self.plot("Cluster Dynamics", "Max |Heat|", max_heat)
            
            # Store results
            self.fra_history.append(fra)
            
            self.cluster_history.append({
                'date': str(date.date()),
                'fra': fra,
                'n_clusters': n_clusters,
                'cluster_heats': cluster_heats,
                'top_10_assets': asset_scores[:10] if asset_scores else [],
                'cluster_sizes': {int(c): int(np.sum(clusters == c)) 
                                for c in unique_clusters},
            })
            
            # Log summary
            top_3 = asset_scores[:3] if asset_scores else []
            top_str = ', '.join([f"{a['name']}={a['composite']:.3f}" for a in top_3])
            self.log(f"{date.date()}: FRA={fra:.3f}, k={n_clusters}, "
                    f"maxH={max_heat:.2f}, top=[{top_str}]")
        
        except Exception as e:
            self.log(f"Error at {date.date()}: {str(e)}")
    
    def _compute_market_fisher(self, corr_matrix, n):
        """Market-level Fisher SV2/SV1."""
        if n < 20:
            return 0.0
        
        k = min(10, n - 1)
        adjacency, neighbors = self._build_knn(corr_matrix, k)
        C_r = self._compute_Cr(corr_matrix, adjacency, n, min(50, n))
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
        
        return float(np.mean(sv2_sv1s)) if sv2_sv1s else 0.0
    
    def _compute_fra(self, sv2sv1, window):
        """Convert SV2/SV1 to FRA on [-1, +1]."""
        self.sv2sv1_history.append(sv2sv1)
        
        if len(self.sv2sv1_history) < 52:
            return 0.0
        
        recent = np.array(self.sv2sv1_history[-LOOKBACK_ZSCORE:])
        z = (sv2sv1 - np.mean(recent)) / (np.std(recent) + 1e-8)
        fra = float(np.tanh(-z / 2.0))
        return fra
    
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
            k_nn = max(4, n_c // 3)
            
            sv2sv1 = self._compute_subgraph_fisher(sub_corr, n_c, k_nn)
            heats[int(c)] = float(sv2sv1)
        
        return heats
    
    def _compute_subgraph_fisher(self, corr_matrix, n, k):
        """Fisher for subgraph."""
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
            k_nn = max(4, n_c // 3)
            
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
        if n < 7:
            return {}
        
        adjacency, neighbors = self._build_knn(corr_matrix, k)
        C_r = self._compute_Cr(corr_matrix, adjacency, n, min(20, n))
        if C_r is None or len(C_r) < 2:
            return {}
        
        adj_sparse = sparse.csr_matrix(adjacency.astype(float))
        
        asset_sv2sv1 = {}
        n_samples = min(20, n)
        sample_idx = np.random.choice(n, n_samples, replace=False)
        
        for v0 in sample_idx:
            result = self._compute_single_fim(C_r, adj_sparse, neighbors, v0, n)
            if result is not None:
                asset_sv2sv1[v0] = result['sv2_sv1']
        
        return asset_sv2sv1
    
    def on_end_of_algorithm(self):
        """Save results to ObjectStore."""
        output = {
            'fra_history': self.fra_history,
            'cluster_history': self.cluster_history,
        }
        self.object_store.save("fisher_regime_attractor_results", 
                                json.dumps(output, default=str))
        self.log(f"Saved {len(self.fra_history)} weeks of results to ObjectStore")
