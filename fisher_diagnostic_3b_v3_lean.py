from AlgorithmImports import *
import numpy as np
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import shortest_path
import json, csv, io

class FisherDiagnostic3B(QCAlgorithm):

    def initialize(self):
        self.set_start_date(2005, 1, 1)
        self.set_end_date(2024, 12, 31)
        self.set_cash(100000)

        self.tickers = [
            "AAPL","MSFT","GOOG","AMZN","JPM","BAC","WFC","GS","MS","C",
            "JNJ","UNH","PFE","MRK","LLY","ABT","TMO","MDT","BMY","AMGN",
            "PG","KO","PEP","WMT","COST","HD","MCD","NKE","SBUX","TGT",
            "XOM","CVX","COP","SLB","EOG","VLO","OXY","HAL",
            "GE","CAT","HON","MMM","UPS","LMT","DE","EMR","ITW","GD",
            "NEE","DUK","SO","D","AEP","EXC","SRE",
            "V","MA","AXP","BLK","MET","PRU","ALL","TRV",
            "INTC","CSCO","ORCL","IBM","ADBE","TXN","QCOM",
            "GILD","CVS","CI","CL","GIS","SYY",
            "VZ","T","CMCSA","DIS",
            "RTX","NSC","CSX","WM","ETN",
            "USB","PNC","BK",
            "LOW","TJX","CRM","AVGO","NFLX","PYPL",
        ]

        self.symbols = {}
        for t in self.tickers:
            try:
                self.symbols[t] = self.add_equity(t, Resolution.DAILY).symbol
            except:
                pass

        self.spy = self.add_equity("SPY", Resolution.DAILY).symbol

        # Try adding VIX — may need custom data depending on LEAN setup
        try:
            self.vix = self.add_data(CBOE, "VIX", Resolution.DAILY).symbol
            self.has_vix = True
        except:
            self.has_vix = False
            self.log("VIX not available as custom data — will record NaN")

        self.windows = [30, 60, 90]
        self.k_nn = 10
        self.n_fisher_samples = 20
        self.n_cr_samples = 50
        self.step = 5  # compute every 5 trading days
        self.day_count = 0

        self.results = []
        self.gate_log = []

        self.set_warm_up(timedelta(days=150))

    def on_data(self, data):
        if self.is_warming_up:
            return

        self.day_count += 1
        if self.day_count % self.step != 0:
            return

        date_str = str(self.time.date())

        # Get SPY and VIX
        spy_price = float(self.securities[self.spy].price) if self.spy in self.securities else None
        vix_price = None
        if self.has_vix:
            try:
                vix_price = float(self.securities[self.vix].price)
            except:
                pass

        row = {"date": date_str, "spy": spy_price, "vix": vix_price}

        for w in self.windows:
            try:
                result = self._compute_window(w)
                if result:
                    row[f"w{w}_sv2sv1"] = round(result["sv2sv1"], 4)
                    row[f"w{w}_rank"] = round(result["rank"], 2)
                    row[f"w{w}_eta"] = round(result["eta"], 4)
                    row[f"w{w}_pr"] = round(result["pr"], 3)
                    row[f"w{w}_n_valid"] = result["n_valid"]
                    self.gate_log.append(f"{date_str} w{w}: OK n={result['n_valid']}")
                else:
                    for k in ["sv2sv1", "rank", "eta", "pr", "n_valid"]:
                        row[f"w{w}_{k}"] = None
                    self.gate_log.append(f"{date_str} w{w}: FAIL (no result)")
            except Exception as e:
                for k in ["sv2sv1", "rank", "eta", "pr", "n_valid"]:
                    row[f"w{w}_{k}"] = None
                self.gate_log.append(f"{date_str} w{w}: ERROR {str(e)[:80]}")

        self.results.append(row)

        if len(self.results) % 50 == 0:
            self.log(f"Processed {len(self.results)} dates through {date_str}")

    def _compute_window(self, window):
        """Full pipeline: history -> correlation -> graph -> C(r) -> FIM -> diagnostics."""
        syms = list(self.symbols.values())
        history = self.history(syms, window, Resolution.DAILY)
        if history.empty:
            return None

        closes = history["close"].unstack(level=0)
        closes = closes.dropna(axis=1, thresh=int(window * 0.8))
        if closes.shape[1] < 50:
            return None

        returns = np.log(closes / closes.shift(1)).dropna()
        if len(returns) < window * 0.7:
            return None

        corr = returns.corr().values
        n = corr.shape[0]
        np.fill_diagonal(corr, 1.0)
        corr = np.nan_to_num(corr, nan=0.0)

        # k-NN graph
        adj = self._build_knn(corr, min(self.k_nn, n - 1))

        # C(r)
        C_r = self._compute_cr(corr, adj, n)
        if C_r is None or len(C_r) < 3:
            return None

        # Fisher diagnostics
        return self._fisher(C_r, adj, n)

    def _build_knn(self, corr, k):
        """k-NN graph from correlation. Returns bool adjacency matrix."""
        n = corr.shape[0]
        adj = np.zeros((n, n), dtype=bool)
        for i in range(n):
            row = corr[i].copy()
            row[i] = -np.inf
            top_k = np.argsort(row)[-k:]
            for j in top_k:
                adj[i, j] = True
                adj[j, i] = True
        return adj

    def _compute_cr(self, corr, adj, n):
        """Radially-averaged correlation C(r) on graph."""
        adj_sp = csr_matrix(adj.astype(float))
        samples = np.random.choice(n, min(self.n_cr_samples, n), replace=False)
        dist = shortest_path(adj_sp, method="D", indices=samples, directed=False)

        max_r = int(np.nanmax(dist[np.isfinite(dist)]))
        if max_r < 3:
            return None
        max_r = min(max_r, 15)

        C_r = np.zeros(max_r + 1)
        counts = np.zeros(max_r + 1)

        for idx_i, src in enumerate(samples):
            for dst in range(n):
                if src == dst:
                    continue
                d = dist[idx_i, dst]
                if np.isfinite(d) and int(d) <= max_r:
                    C_r[int(d)] += corr[src, dst]
                    counts[int(d)] += 1

        mask = counts > 0
        C_r[mask] /= counts[mask]
        return C_r

    def _fisher(self, C_r, adj, n):
        """Compute FIM diagnostics from C(r) kernel on graph."""
        adj_sp = csr_matrix(adj.astype(float))
        neighbors = {i: list(np.where(adj[i])[0]) for i in range(n)}
        samples = np.random.choice(n, min(self.n_fisher_samples, n), replace=False)

        ranks, etas, prs, sv2sv1s = [], [], [], []

        for v0 in samples:
            nbrs = neighbors.get(v0, [])
            if len(nbrs) < 2:
                continue

            sources = [v0] + list(nbrs)
            try:
                dists = shortest_path(adj_sp, method="D", indices=sources, directed=False)
            except:
                continue

            # Build kernel for v0
            p_v0 = self._kernel(C_r, dists[0], n)
            if p_v0 is None:
                continue

            # Score vectors
            k = len(nbrs)
            scores = np.zeros((k, n))
            valid = True
            for j in range(k):
                p_wj = self._kernel(C_r, dists[j + 1], n)
                if p_wj is None:
                    valid = False
                    break
                with np.errstate(divide="ignore", invalid="ignore"):
                    log_ratio = np.log(p_wj) - np.log(p_v0)
                    log_ratio = np.nan_to_num(log_ratio, nan=0.0, posinf=0.0, neginf=0.0)
                scores[j] = log_ratio

            if not valid:
                continue

            # FIM
            weighted = scores * np.sqrt(p_v0)[np.newaxis, :]
            F = weighted @ weighted.T
            sv = np.linalg.svd(F, compute_uv=False)
            if sv[0] < 1e-15:
                continue
            sv = sv / sv[0]

            # Diagnostics
            gaps = sv[:-1] / (sv[1:] + 1e-15)
            rank = int(np.argmax(gaps) + 1)
            eta = float(sv[min(rank, len(sv) - 1)] / (sv[rank - 1] + 1e-15))
            pr = float((np.sum(sv)) ** 2 / (np.sum(sv ** 2) + 1e-15))
            s2s1 = float(sv[1] / (sv[0] + 1e-15)) if len(sv) > 1 else 0.0

            ranks.append(rank)
            etas.append(eta)
            prs.append(pr)
            sv2sv1s.append(s2s1)

        if len(ranks) == 0:
            return None

        return {
            "sv2sv1": float(np.mean(sv2sv1s)),
            "rank": float(np.mean(ranks)),
            "eta": float(np.mean(etas)),
            "pr": float(np.mean(prs)),
            "n_valid": len(ranks),
        }

    def _kernel(self, C_r, distances, n):
        """Build probability kernel from C(r) and BFS distances."""
        weights = np.zeros(n)
        for u in range(n):
            d = distances[u]
            if np.isfinite(d) and int(d) < len(C_r):
                weights[u] = abs(C_r[int(d)])
        total = np.sum(weights)
        if total < 1e-15:
            return None
        return weights / total

    def on_end_of_algorithm(self):
        """Save CSV results and gate log."""
        if not self.results:
            self.log("No results to save")
            return

        # Build CSV
        fields = ["date", "spy", "vix"]
        for w in self.windows:
            fields += [f"w{w}_sv2sv1", f"w{w}_rank", f"w{w}_eta", f"w{w}_pr", f"w{w}_n_valid"]

        buf = io.StringIO()
        writer = csv.DictWriter(buf, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        for row in self.results:
            writer.writerow(row)

        csv_str = buf.getvalue()
        self.object_store.save("phase3b_results_csv", csv_str)
        self.log(f"Saved {len(self.results)} rows to ObjectStore:phase3b_results_csv")

        # Gate log
        gate_str = "\n".join(self.gate_log)
        self.object_store.save("phase3b_gate_log", gate_str)
        self.log(f"Saved {len(self.gate_log)} gate entries")

        # Summary stats
        self.log("=" * 72)
        self.log("PHASE 3B v3: FISHER DIAGNOSTIC — END OF RUN SUMMARY")
        self.log("=" * 72)
        self.log(f"Period: {self.results[0]['date']} to {self.results[-1]['date']}")
        self.log(f"Total computation dates: {len(self.results)}")
        self.log(f"Universe: {len(self.symbols)} assets resolved")

        for w in self.windows:
            sv_vals = [r.get(f"w{w}_sv2sv1") for r in self.results if r.get(f"w{w}_sv2sv1") is not None]
            rk_vals = [r.get(f"w{w}_rank") for r in self.results if r.get(f"w{w}_rank") is not None]
            eta_vals = [r.get(f"w{w}_eta") for r in self.results if r.get(f"w{w}_eta") is not None]
            pr_vals = [r.get(f"w{w}_pr") for r in self.results if r.get(f"w{w}_pr") is not None]
            nv_vals = [r.get(f"w{w}_n_valid") for r in self.results if r.get(f"w{w}_n_valid") is not None]
            null_ct = sum(1 for r in self.results if r.get(f"w{w}_sv2sv1") is None)

            if sv_vals:
                self.log(f"--- Window {w}d ({len(sv_vals)} valid / {null_ct} null) ---")
                self.log(f"  SV2/SV1: mean={np.mean(sv_vals):.4f}  std={np.std(sv_vals):.4f}  "
                         f"min={np.min(sv_vals):.4f}  max={np.max(sv_vals):.4f}")
                self.log(f"  Rank:    mean={np.mean(rk_vals):.2f}  std={np.std(rk_vals):.2f}  "
                         f"min={np.min(rk_vals):.2f}  max={np.max(rk_vals):.2f}")
                self.log(f"  Eta:     mean={np.mean(eta_vals):.4f}  std={np.std(eta_vals):.4f}  "
                         f"min={np.min(eta_vals):.4f}  max={np.max(eta_vals):.4f}")
                self.log(f"  PR:      mean={np.mean(pr_vals):.3f}  std={np.std(pr_vals):.3f}  "
                         f"min={np.min(pr_vals):.3f}  max={np.max(pr_vals):.3f}")
                self.log(f"  N_valid: mean={np.mean(nv_vals):.1f}  min={np.min(nv_vals)}  max={np.max(nv_vals)}")

                # 2sigma threshold crossings for SV2/SV1
                arr = np.array(sv_vals)
                mean_sv = np.mean(arr)
                std_sv = np.std(arr)
                threshold = mean_sv + 2 * std_sv
                dates_valid = [r["date"] for r in self.results if r.get(f"w{w}_sv2sv1") is not None]
                crossings = [(dates_valid[i], sv_vals[i]) for i in range(len(sv_vals)) if sv_vals[i] > threshold]
                self.log(f"  2sigma threshold: {threshold:.4f}  crossings: {len(crossings)}/{len(sv_vals)} "
                         f"({100*len(crossings)/len(sv_vals):.1f}%)")
                if crossings:
                    self.log(f"  First crossing: {crossings[0][0]} (val={crossings[0][1]:.4f})")
                    self.log(f"  Last crossing:  {crossings[-1][0]} (val={crossings[-1][1]:.4f})")

        # VIX correlation
        sv90 = []
        vix_vals = []
        for r in self.results:
            s = r.get("w90_sv2sv1")
            v = r.get("vix")
            if s is not None and v is not None:
                sv90.append(s)
                vix_vals.append(v)
        if len(sv90) > 10:
            corr_sv_vix = float(np.corrcoef(sv90, vix_vals)[0, 1])
            self.log(f"--- Correlations ---")
            self.log(f"  SV2/SV1(90d) vs VIX: r = {corr_sv_vix:.3f}")

        # Gate summary
        ok_count = sum(1 for g in self.gate_log if "OK" in g)
        fail_count = sum(1 for g in self.gate_log if "FAIL" in g)
        err_count = sum(1 for g in self.gate_log if "ERROR" in g)
        self.log(f"--- Gate Summary ---")
        self.log(f"  OK: {ok_count}  FAIL: {fail_count}  ERROR: {err_count}")

        self.log("=" * 72)
        self.log("Upload phase3b_results_csv from ObjectStore to chat for analysis.")
        self.log("=" * 72)
