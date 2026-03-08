# DS Phase 3B-2: Sector-Level Fisher Decomposition
# LEAN Algorithm for QuantConnect
# Version: 1.0
# Target: Identify sector-level stress before market-wide propagation

from AlgorithmImports import *
import numpy as np
from scipy.sparse.csgraph import shortest_path
from scipy import sparse
import json

class FisherSectorAlgorithm(QCAlgorithm):
    
    def initialize(self):
        self.set_start_date(2005, 1, 1)
        self.set_end_date(2024, 12, 31)
        self.set_cash(100000)
        
        # --- Sector definitions (GICS-based) ---
        self.sector_map = {
            'Technology': ["AAPL","MSFT","GOOG","AMZN","INTC","CSCO","ORCL","IBM",
                          "ADBE","TXN","QCOM","AVGO"],
            'Financials': ["JPM","BAC","WFC","GS","MS","C","USB","PNC","BK",
                          "AXP","MET","PRU","ALL","TRV"],
            'HealthCare': ["JNJ","UNH","PFE","MRK","LLY","TMO","ABT","MDT","BMY",
                          "AMGN","GILD","CVS","CI"],
            'Consumer':   ["PG","KO","PEP","WMT","COST","HD","MCD","NKE","SBUX",
                          "TGT","CL","GIS","SYY"],
            'Industrials':["GE","CAT","HON","MMM","UPS","RTX","LMT","DE","EMR",
                          "ITW","GD","NSC","CSX","WM","ETN"],
            'Energy':     ["XOM","CVX","COP","SLB","EOG","PSX","VLO","MPC","OXY","HAL"],
            'Utilities':  ["NEE","DUK","SO","D","AEP","EXC","SRE"],
        }
        
        # Flatten all tickers + add SPY
        all_tickers = set()
        for tickers in self.sector_map.values():
            all_tickers.update(tickers)
        all_tickers.add("SPY")
        
        self.symbols = {}
        for ticker in all_tickers:
            eq = self.add_equity(ticker, Resolution.DAILY)
            eq.set_data_normalization_mode(DataNormalizationMode.ADJUSTED)
            self.symbols[ticker] = eq.symbol
        
        # Add VIX
        self.vix = self.add_data(CBOE, "VIX", Resolution.DAILY).symbol
        
        # Schedule weekly computation (every Monday)
        self.schedule.on(
            self.date_rules.every(DayOfWeek.MONDAY),
            self.time_rules.after_market_open("SPY", 30),
            self.compute_all_fisher
        )
        
        # Storage for results
        self.results = {
            'dates': [],
            'market': [],      # market-level SV2/SV1
            'market_rank': [],
            'market_pr': [],
            'vix': [],
            'spy': [],
        }
        
        # Initialize sector storage
        for sector in self.sector_map:
            self.results[f'sector_{sector}'] = []        # sector SV2/SV1 (90d)
            self.results[f'sector_{sector}_30d'] = []    # sector SV2/SV1 (30d)
            self.results[f'spread_{sector}'] = []        # 30d-90d spread
            self.results[f'rank_{sector}'] = []          # sector rank
            self.results[f'pr_{sector}'] = []            # sector participation ratio
            self.results[f'assets_{sector}'] = []        # per-asset Fisher temps
            self.results[f'n_assets_{sector}'] = []      # number of active assets
        
        # LEAN charts
        sector_chart = Chart("Sector SV2/SV1")
        for sector in self.sector_map:
            sector_chart.add_series(Series(sector, SeriesType.LINE, 0))
        sector_chart.add_series(Series("Market", SeriesType.LINE, 0))
        self.add_chart(sector_chart)
        
        spread_chart = Chart("Sector Spread (30d-90d)")
        for sector in self.sector_map:
            spread_chart.add_series(Series(sector, SeriesType.LINE, 0))
        self.add_chart(spread_chart)
        
        vix_chart = Chart("VIX")
        vix_chart.add_series(Series("VIX", SeriesType.LINE, 0))
        self.add_chart(vix_chart)
        
        # Warmup period (150 days for 90-day window + buffer)
        self.set_warm_up(timedelta(days=150))
        
        self.log("FisherSectorAlgorithm initialized: 7 sectors, 2005-2024")
    
    def compute_all_fisher(self):
        """Weekly: compute market-level and all sector-level Fisher diagnostics."""
        if self.is_warming_up:
            return
        
        date = self.time
        self.results['dates'].append(str(date.date()))
        
        # SPY and VIX
        spy_price = self.securities["SPY"].price if self.securities["SPY"].has_data else float('nan')
        self.results['spy'].append(float(spy_price))
        
        vix_val = self.securities[self.vix].price if self.securities[self.vix].has_data else float('nan')
        self.results['vix'].append(float(vix_val))
        if not np.isnan(vix_val):
            self.plot("VIX", "VIX", vix_val)
        
        # --- Market level (90d, full universe, k=10) ---
        all_active = [t for t in self.symbols.keys() if t != "SPY" 
                      and self.securities[self.symbols[t]].has_data 
                      and self.securities[self.symbols[t]].price > 0]
        
        mkt_result = self._compute_fisher(all_active, window=90, k_nn=10)
        mkt_sv = mkt_result['sv2_sv1'] if mkt_result else float('nan')
        mkt_rank = mkt_result['rank'] if mkt_result else float('nan')
        mkt_pr = mkt_result['pr'] if mkt_result else float('nan')
        
        self.results['market'].append(mkt_sv)
        self.results['market_rank'].append(mkt_rank)
        self.results['market_pr'].append(mkt_pr)
        
        if not np.isnan(mkt_sv):
            self.plot("Sector SV2/SV1", "Market", mkt_sv)
        
        # --- Sector level ---
        for sector, tickers in self.sector_map.items():
            active = [t for t in tickers 
                     if t in self.symbols and self.securities[self.symbols[t]].has_data 
                     and self.securities[self.symbols[t]].price > 0]
            
            n_active = len(active)
            self.results[f'n_assets_{sector}'].append(n_active)
            
            # Determine k based on sector size
            k = 5 if n_active >= 10 else 4
            
            # 90d window
            r90 = self._compute_fisher(active, window=90, k_nn=k) if n_active >= 7 else None
            sv90 = r90['sv2_sv1'] if r90 else float('nan')
            rank90 = r90['rank'] if r90 else float('nan')
            pr90 = r90['pr'] if r90 else float('nan')
            
            self.results[f'sector_{sector}'].append(sv90)
            self.results[f'rank_{sector}'].append(rank90)
            self.results[f'pr_{sector}'].append(pr90)
            
            # 30d window for spread
            r30 = self._compute_fisher(active, window=30, k_nn=k) if n_active >= 7 else None
            sv30 = r30['sv2_sv1'] if r30 else float('nan')
            self.results[f'sector_{sector}_30d'].append(sv30)
            
            # Spread = 30d - 90d
            spread = (sv30 - sv90) if (not np.isnan(sv30) and not np.isnan(sv90)) else float('nan')
            self.results[f'spread_{sector}'].append(spread)
            
            # Per-asset Fisher temperatures (90d)
            asset_temps = {}
            if r90 and 'asset_sv2_sv1' in r90:
                asset_temps = r90['asset_sv2_sv1']
            self.results[f'assets_{sector}'].append(asset_temps)
            
            # Plot
            if not np.isnan(sv90):
                self.plot("Sector SV2/SV1", sector, sv90)
            if not np.isnan(spread):
                self.plot("Sector Spread (30d-90d)", sector, spread)
        
        # Log summary
        fin_sv = self.results['sector_Financials'][-1]
        tech_sv = self.results['sector_Technology'][-1]
        energy_sv = self.results['sector_Energy'][-1]
        
        self.log(f"{date.date()}: Mkt={mkt_sv:.3f}, Fin={fin_sv:.3f}, Tech={tech_sv:.3f}, Energy={energy_sv:.3f}, VIX={vix_val:.1f}")
    
    def _compute_fisher(self, tickers, window, k_nn):
        """
        Compute Fisher diagnostics for a set of tickers.
        Returns dict with sv2_sv1, rank, eta, pr, gap_ratio, 
        and asset_sv2_sv1 (per-asset Fisher temperatures).
        
        Adapted from Phase 3B — same pipeline, parameterized k_nn.
        """
        symbols = [self.symbols[t] for t in tickers if t in self.symbols]
        if len(symbols) < 7:
            return None
        
        try:
            history = self.history(symbols, window, Resolution.DAILY)
            if history.empty:
                return None
            
            closes = history['close'].unstack(level=0)
            closes = closes.dropna(axis=1, thresh=int(window * 0.8))
            if closes.shape[1] < 7:
                return None
            
            returns = np.log(closes / closes.shift(1)).dropna()
            if len(returns) < window * 0.7:
                return None
            
            corr_matrix = returns.corr().values
            n = corr_matrix.shape[0]
            np.fill_diagonal(corr_matrix, 1.0)
            corr_matrix = np.nan_to_num(corr_matrix, nan=0.0)
            
            # Column names for asset mapping
            col_names = list(returns.columns)
            
            # k-NN graph
            adjacency, neighbors = self._build_knn(corr_matrix, min(k_nn, n - 1))
            
            # C(r)
            C_r = self._compute_Cr(corr_matrix, adjacency, n, min(50, n))
            if C_r is None or len(C_r) < 2:
                return None
            
            # Fisher diagnostics (sampled)
            n_samples = min(20, n)
            sample_idx = np.random.choice(n, n_samples, replace=False)
            
            adj_sparse = sparse.csr_matrix(adjacency.astype(float))
            
            ranks, etas, prs, gap_ratios, sv2_sv1s = [], [], [], [], []
            asset_sv2_sv1 = {}  # per-asset Fisher temperature
            
            for v0 in sample_idx:
                result = self._compute_single_fim(C_r, adj_sparse, neighbors, v0, n)
                if result is None:
                    continue
                ranks.append(result['rank'])
                etas.append(result['eta'])
                prs.append(result['pr'])
                gap_ratios.append(result['gap_ratio'])
                sv2_sv1s.append(result['sv2_sv1'])
                
                # Map back to ticker name
                if v0 < len(col_names):
                    ticker_name = str(col_names[v0])
                    asset_sv2_sv1[ticker_name] = result['sv2_sv1']
            
            if len(ranks) == 0:
                return None
            
            return {
                'sv2_sv1': float(np.mean(sv2_sv1s)),
                'rank': float(np.mean(ranks)),
                'eta': float(np.mean(etas)),
                'pr': float(np.mean(prs)),
                'gap_ratio': float(np.mean(gap_ratios)),
                'asset_sv2_sv1': asset_sv2_sv1,
                'n_valid': len(ranks)
            }
        except Exception as e:
            self.log(f"Fisher computation error: {str(e)}")
            return None
    
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
        
        pr = (np.sum(sv))**2 / (np.sum(sv**2) + 1e-15)
        eta = sv[rank] / (sv[rank-1] + 1e-15) if rank < len(sv) else 0.0
        sv2_sv1 = float(sv_norm[1]) if len(sv_norm) > 1 else 0.0
        gap_ratio = float(sv[0] / (sv[1] + 1e-15)) if len(sv) > 1 else float('inf')
        
        return {'rank': rank, 'eta': eta, 'pr': pr, 'sv2_sv1': sv2_sv1, 'gap_ratio': gap_ratio}
    
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
    
    def on_end_of_algorithm(self):
        """Save all results to ObjectStore."""
        self.object_store.save("fisher_sector_results", json.dumps(self.results))
        self.log(f"Saved {len(self.results['dates'])} dates to ObjectStore")
        
        # Summary statistics
        n_dates = len(self.results['dates'])
        mkt_valid = np.sum(~np.isnan(self.results['market']))
        self.log(f"Market-level valid: {mkt_valid}/{n_dates} ({100*mkt_valid/n_dates:.1f}%)")
        
        for sector in self.sector_map:
            sector_valid = np.sum(~np.isnan(self.results[f'sector_{sector}']))
            self.log(f"{sector} valid: {sector_valid}/{n_dates} ({100*sector_valid/n_dates:.1f}%)")
