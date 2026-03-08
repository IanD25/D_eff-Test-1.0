# region imports
from AlgorithmImports import *
import numpy as np
from collections import deque
import json
import math
# endregion

# ============================================================================
# Phase 3B: Financial Correlation Network — Fisher Diagnostic Algorithm
# ============================================================================
# Uses the DS Framework's Fisher Information Matrix diagnostics on
# financial correlation networks. Tests whether SV2/SV1 rises before
# known market crises (2008, 2020).
#
# Architecture: QCAlgorithm subclass that computes weekly Fisher diagnostics
# on a k-NN graph built from rolling correlation matrices of ~100 S&P 500
# stocks. Results saved to ObjectStore for post-processing.
# ============================================================================

class FisherDiagnosticAlgorithm(QCAlgorithm):

    def initialize(self):
        """
        LEAN required initialization.
        1. Set backtest period: 2005-01-01 to 2024-12-31
        2. Add equity universe: ~100 liquid S&P 500 stocks
        3. Schedule weekly Fisher computation
        4. Initialize storage structures
        """
        self.set_start_date(2005, 1, 1)
        self.set_end_date(2024, 12, 31)
        self.set_cash(100000)  # nominal - algorithm does not trade

        # --- Universe: ~100 highly liquid S&P 500 constituents ---
        # Fixed list chosen for sector diversity and long data availability
        self.tickers = [
            # Tech (12)
            "AAPL", "MSFT", "GOOG", "AMZN", "INTC", "CSCO", "ORCL", "IBM",
            "ADBE", "TXN", "QCOM", "AVGO",
            # Finance (14)
            "JPM", "BAC", "WFC", "GS", "MS", "C", "USB", "PNC", "BK",
            "AXP", "MET", "PRU", "ALL", "TRV",
            # Health (13)
            "JNJ", "UNH", "PFE", "MRK", "LLY", "TMO", "ABT", "MDT", "BMY",
            "AMGN", "GILD", "CVS", "CI",
            # Consumer (13)
            "PG", "KO", "PEP", "WMT", "COST", "HD", "MCD", "NKE", "SBUX", "TGT",
            "CL", "GIS", "SYY",
            # Industrial (15)
            "GE", "CAT", "HON", "MMM", "UPS", "RTX", "LMT", "DE", "EMR", "ITW",
            "GD", "NSC", "CSX", "WM", "ETN",
            # Energy (10)
            "XOM", "CVX", "COP", "SLB", "EOG", "PSX", "VLO", "MPC", "OXY", "HAL",
            # Utilities (7)
            "NEE", "DUK", "SO", "D", "AEP", "EXC", "SRE",
            # Telecom/Media (5)
            "VZ", "T", "CMCSA", "DIS", "NFLX",
            # Index reference
            "SPY",
        ]

        self.symbols = []
        for ticker in self.tickers:
            try:
                equity = self.add_equity(ticker, Resolution.DAILY)
                equity.set_data_normalization_mode(DataNormalizationMode.ADJUSTED)
                self.symbols.append(equity.symbol)
            except Exception as e:
                self.log(f"Could not add {ticker}: {str(e)}")

        # --- Try to add VIX ---
        self.vix_symbol = None
        try:
            vix_data = self.add_data(CBOE, "VIX", Resolution.DAILY)
            self.vix_symbol = vix_data.symbol
        except Exception as e:
            self.log(f"VIX via CBOE not available: {str(e)}. Will skip VIX.")

        # --- Schedule weekly Fisher computation ---
        # Run every Monday 30 min after market open
        self.schedule.on(
            self.date_rules.every(DayOfWeek.MONDAY),
            self.time_rules.after_market_open("SPY", 30),
            self.compute_fisher_diagnostics
        )

        # --- Parameters ---
        self.windows = [30, 60, 90, 120]
        self.k_nn = 10  # reduced from 12 to ensure deeper graph structure with ~100 nodes
        self.n_fisher_samples = 30
        self.n_cr_samples = 50  # reduced for speed

        # --- Results storage ---
        self.results = {w: [] for w in self.windows}
        self.dates = []
        self.vix_values = []
        self.spy_values = []
        self.computation_count = 0

        # --- LEAN Charts ---
        sv_chart = Chart("Fisher SV2/SV1")
        sv_chart.add_series(Series("SV2/SV1 (90d)", SeriesType.LINE, 0))
        sv_chart.add_series(Series("VIX (scaled)", SeriesType.LINE, 1))
        self.add_chart(sv_chart)

        eta_chart = Chart("Fisher Eta")
        eta_chart.add_series(Series("Eta (90d)", SeriesType.LINE, 0))
        self.add_chart(eta_chart)

        rank_chart = Chart("Fisher Rank")
        rank_chart.add_series(Series("Rank (90d)", SeriesType.LINE, 0))
        self.add_chart(rank_chart)

        multi_chart = Chart("Multi-Scale SV2/SV1")
        for w in self.windows:
            multi_chart.add_series(Series(f"SV2/SV1 ({w}d)", SeriesType.LINE, 0))
        self.add_chart(multi_chart)

        self.set_warm_up(timedelta(days=150))

    # ========================================================================
    # Core computation: scheduled weekly
    # ========================================================================

    def compute_fisher_diagnostics(self):
        """
        Scheduled weekly: compute Fisher diagnostics at all four windows.
        """
        if self.is_warming_up:
            return

        current_date = self.time

        # Get list of symbols with active data (exclude SPY from analysis universe)
        active_symbols = []
        for s in self.symbols:
            if s.value == "SPY":
                continue
            try:
                sec = self.securities[s]
                if sec.has_data and sec.price > 0:
                    active_symbols.append(s)
            except Exception:
                continue

        if len(active_symbols) < 40:
            return  # Not enough assets

        # Record SPY and VIX
        try:
            spy_price = float(self.securities["SPY"].price)
        except Exception:
            spy_price = float('nan')

        vix_price = float('nan')
        if self.vix_symbol is not None:
            try:
                if self.securities[self.vix_symbol].has_data:
                    vix_price = float(self.securities[self.vix_symbol].price)
            except Exception:
                pass

        self.dates.append(str(current_date.date()))
        self.spy_values.append(spy_price)
        self.vix_values.append(vix_price)

        for window in self.windows:
            try:
                result = self._compute_single_window(active_symbols, window)
                self.results[window].append(result)

                # Plot 90-day primary window metrics
                if window == 90 and result is not None:
                    self.plot("Fisher SV2/SV1", "SV2/SV1 (90d)", result['sv2_sv1'])
                    if not math.isnan(vix_price):
                        self.plot("Fisher SV2/SV1", "VIX (scaled)", vix_price / 100.0)
                    self.plot("Fisher Eta", "Eta (90d)", result['eta'])
                    self.plot("Fisher Rank", "Rank (90d)", result['rank'])

                # Plot all windows on multi-scale chart
                if result is not None:
                    self.plot("Multi-Scale SV2/SV1", f"SV2/SV1 ({window}d)",
                              result['sv2_sv1'])

            except Exception as e:
                self.log(f"Error at {current_date} window={window}: {str(e)}")
                self.results[window].append(None)

        # Progress logging (every 10 computations)
        self.computation_count += 1
        if self.computation_count % 10 == 0:
            r90 = self.results[90][-1]
            if r90 is not None:
                self.log(f"[{self.computation_count}] {current_date.date()}: "
                         f"rank={r90['rank']:.1f}, eta={r90['eta']:.3f}, "
                         f"SV2/SV1={r90['sv2_sv1']:.3f}, "
                         f"n_assets={len(active_symbols)}, VIX={vix_price:.1f}")
            else:
                self.log(f"[{self.computation_count}] {current_date.date()}: "
                         f"90d window returned None, n_assets={len(active_symbols)}")

    # ========================================================================
    # Single-window computation
    # ========================================================================

    def _compute_single_window(self, symbols, window):
        """
        For a single rolling window: correlation matrix -> k-NN graph ->
        C(r) -> Fisher diagnostics. Returns dict or None.
        """
        # Step 1: Get historical returns
        history = self.history(symbols, window, Resolution.DAILY)
        if history.empty:
            return None

        # Pivot to (dates x assets) close price matrix
        try:
            closes = history['close'].unstack(level=0)
        except Exception:
            return None

        # Drop assets with too much missing data
        closes = closes.dropna(axis=1, thresh=int(window * 0.7))
        if closes.shape[1] < 40:
            return None

        # Forward-fill remaining gaps, then drop any remaining NaN rows
        closes = closes.ffill().dropna()

        # Compute log returns
        returns = np.log(closes / closes.shift(1)).dropna()
        if len(returns) < int(window * 0.5):
            return None

        # Step 2: Correlation matrix
        corr_matrix = returns.corr().values  # N x N numpy array
        n = corr_matrix.shape[0]

        # Clean up
        np.fill_diagonal(corr_matrix, 1.0)
        corr_matrix = np.nan_to_num(corr_matrix, nan=0.0)

        # Step 3: k-NN graph
        adjacency, neighbors = self._build_knn_graph(corr_matrix, self.k_nn)

        # Check graph connectivity
        largest_component = self._largest_component_size(adjacency, n)
        if largest_component < n * 0.8:
            # Graph too fragmented
            return None

        # Step 4: C(r) via sampled BFS
        C_r = self._compute_Cr(corr_matrix, adjacency, n, self.n_cr_samples)
        if C_r is None or len(C_r) < 3:
            return None

        # Step 5: Fisher diagnostics
        diagnostics = self._fisher_diagnostics(C_r, adjacency, neighbors,
                                                n, self.n_fisher_samples)
        if diagnostics is not None:
            diagnostics['n_assets'] = n
            diagnostics['graph_size'] = largest_component
            diagnostics['cr_length'] = len(C_r)
            diagnostics['cr_values'] = C_r.tolist()

        return diagnostics

    # ========================================================================
    # Graph construction
    # ========================================================================

    def _build_knn_graph(self, corr_matrix, k):
        """
        Build k-NN graph from correlation matrix. Connect each asset to
        its k most-correlated peers. Symmetrize.
        """
        n = corr_matrix.shape[0]
        adjacency = np.zeros((n, n), dtype=bool)

        for i in range(n):
            row = corr_matrix[i].copy()
            row[i] = -np.inf  # exclude self
            top_k = np.argsort(row)[-k:]
            for j in top_k:
                adjacency[i, j] = True
                adjacency[j, i] = True  # symmetrize

        # Build neighbors dict
        neighbors = {}
        for i in range(n):
            neighbors[i] = list(np.where(adjacency[i])[0])

        return adjacency, neighbors

    def _largest_component_size(self, adjacency, n):
        """Find size of largest connected component via BFS."""
        visited = set()
        max_size = 0
        for start in range(n):
            if start in visited:
                continue
            # BFS
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

    # ========================================================================
    # C(r) computation via BFS
    # ========================================================================

    def _bfs_distances(self, adjacency, source, n, max_dist=15):
        """BFS from source, return distance array. Pure Python."""
        dist = np.full(n, -1, dtype=int)
        dist[source] = 0
        queue = deque([source])
        while queue:
            node = queue.popleft()
            if dist[node] >= max_dist:
                continue
            for neighbor in range(n):
                if adjacency[node, neighbor] and dist[neighbor] == -1:
                    dist[neighbor] = dist[node] + 1
                    queue.append(neighbor)
        return dist

    def _bfs_distances_fast(self, neighbors_list, source, n, max_dist=15):
        """BFS from source using adjacency list (faster). Pure Python."""
        dist = np.full(n, -1, dtype=int)
        dist[source] = 0
        queue = deque([source])
        while queue:
            node = queue.popleft()
            if dist[node] >= max_dist:
                continue
            for neighbor in neighbors_list[node]:
                if dist[neighbor] == -1:
                    dist[neighbor] = dist[node] + 1
                    queue.append(neighbor)
        return dist

    def _compute_Cr(self, corr_matrix, adjacency, n, n_samples):
        """
        Compute radially-averaged correlation C(r) on the k-NN graph.
        """
        # Build adjacency list for faster BFS
        nbr_list = []
        for i in range(n):
            nbr_list.append(list(np.where(adjacency[i])[0]))

        # Sample vertices
        sample_idx = np.random.choice(n, min(n_samples, n), replace=False)
        max_r = 15

        C_r = np.zeros(max_r + 1)
        counts = np.zeros(max_r + 1)

        for src in sample_idx:
            dist = self._bfs_distances_fast(nbr_list, int(src), n, max_r)
            for dst in range(n):
                if src == dst:
                    continue
                d = dist[dst]
                if 0 < d <= max_r:
                    C_r[d] += corr_matrix[src, dst]
                    counts[d] += 1

        # Average
        mask = counts > 0
        C_r[mask] /= counts[mask]
        C_r[~mask] = 0.0

        # Find effective max_r (last distance with data)
        effective_max = max_r
        for r in range(max_r, 0, -1):
            if counts[r] > 0:
                effective_max = r
                break

        return C_r[:effective_max + 1]

    # ========================================================================
    # Fisher diagnostics (adapted from Phase 2)
    # ========================================================================

    def _fisher_diagnostics(self, C_r, adjacency, neighbors, n, n_samples):
        """
        Compute Fisher diagnostics on the financial graph using C(r) as kernel.
        """
        # Build adjacency list for BFS
        nbr_list = []
        for i in range(n):
            nbr_list.append(list(np.where(adjacency[i])[0]))

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

            # BFS distances from v0
            dist_v0 = self._bfs_distances_fast(nbr_list, v0, n, 15)

            # Build kernel for v0
            p_v0 = self._build_kernel(C_r, dist_v0, n)
            if p_v0 is None:
                continue

            # Build kernels for each neighbor, compute score vectors
            k = len(nbrs)
            score_vectors = np.zeros((k, n))

            valid = True
            for j, wj in enumerate(nbrs):
                dist_wj = self._bfs_distances_fast(nbr_list, int(wj), n, 15)
                p_wj = self._build_kernel(C_r, dist_wj, n)
                if p_wj is None:
                    valid = False
                    break
                # Score vector: log p_wj - log p_v0
                eps = 1e-12
                score_vectors[j] = np.log(p_wj + eps) - np.log(p_v0 + eps)

            if not valid:
                continue

            # FIM: weighted inner product
            weighted_scores = score_vectors * np.sqrt(p_v0)[np.newaxis, :]
            FIM = weighted_scores @ weighted_scores.T  # k x k

            # SVD
            try:
                sv = np.linalg.svd(FIM, compute_uv=False)
            except np.linalg.LinAlgError:
                continue

            sv = np.sort(sv)[::-1]  # descending
            if sv[0] < 1e-15:
                continue

            sv_norm = sv / sv[0]  # normalize

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

            # Participation ratio
            pr = (np.sum(sv)) ** 2 / (np.sum(sv ** 2) + 1e-15)

            # Disorder index (eta)
            if rank < len(sv):
                eta = sv[rank] / (sv[rank - 1] + 1e-15)
            else:
                eta = 0.0

            # SV2/SV1
            sv2_sv1 = float(sv_norm[1]) if len(sv_norm) > 1 else 0.0

            # Gap ratio
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

    def _build_kernel(self, C_r, distances, n):
        """Build probability distribution from C(r) kernel and distances."""
        weights = np.zeros(n)
        for u in range(n):
            d = int(distances[u])
            if 0 <= d < len(C_r):
                weights[u] = abs(C_r[d])
            # else: weight stays 0 (unreachable vertices)

        total = np.sum(weights)
        if total < 1e-15:
            return None
        return weights / total

    # ========================================================================
    # End of algorithm: save results
    # ========================================================================

    def on_end_of_algorithm(self):
        """Save all results to ObjectStore for post-processing."""
        output = {
            'dates': self.dates,
            'vix': self.vix_values,
            'spy': self.spy_values,
            'windows': self.windows,
            'k_nn': self.k_nn,
            'n_fisher_samples': self.n_fisher_samples,
            'n_cr_samples': self.n_cr_samples,
            'tickers': self.tickers,
        }

        for w in self.windows:
            # Convert None entries to empty dicts for JSON
            cleaned = []
            for r in self.results[w]:
                if r is None:
                    cleaned.append({})
                else:
                    cleaned.append(r)
            output[f'results_w{w}'] = cleaned

        # Save to ObjectStore
        json_str = json.dumps(output)
        self.object_store.save("fisher_financial_results", json_str)
        self.log(f"Saved {len(self.dates)} data points to ObjectStore")
        if self.dates:
            self.log(f"Date range: {self.dates[0]} to {self.dates[-1]}")

        # Summary statistics
        r90 = output.get('results_w90', [])
        valid_90 = [r for r in r90 if r and 'sv2_sv1' in r]
        if valid_90:
            sv_vals = [r['sv2_sv1'] for r in valid_90]
            self.log(f"SV2/SV1 (90d) — mean={np.mean(sv_vals):.3f}, "
                     f"std={np.std(sv_vals):.3f}, "
                     f"max={np.max(sv_vals):.3f}, "
                     f"n_valid={len(valid_90)}/{len(r90)}")
