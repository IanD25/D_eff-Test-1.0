from AlgorithmImports import *
import numpy as np
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import shortest_path
from collections import deque
import csv, io

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

        try:
            self.vix = self.add_data(CBOE, "VIX", Resolution.DAILY).symbol
            self.has_vix = True
        except:
            self.has_vix = False
            self.log("VIX not available as custom data — will record NaN")

        self.windows = [30, 60, 90]
        self.k_nn = 10
        self.n_fisher_samples = 15   # SPEED: reduced from 20 (~25% faster Fisher)
        self.n_cr_samples = 30       # SPEED: reduced from 50 (~40% faster C(r))
        self.step = 5
        self.day_count = 0

        # Rolling threshold tracking (~252 trading days = 50 steps at step=5)
        self.rolling_n = 50
        self.rolling_min_periods = 20   # need ≥20 periods before threshold is valid
        self.sv_history = {w: deque(maxlen=self.rolling_n) for w in self.windows}
        self.above_thresh  = {w: False for w in self.windows}
        self.rolling_crossings = {w: [] for w in self.windows}  # (dir, date, sv, thresh)

        # Crisis reference dates for proximity logging
        self.crisis_dates = {
            "Bear Stearns": "2008-03-14",
            "Lehman":       "2008-09-15",
            "COVID onset":  "2020-02-24",
            "COVID bottom": "2020-03-23",
        }

        self.results  = []
        self.gate_log = []

        self.set_warm_up(timedelta(days=150))

    # ------------------------------------------------------------------
    def on_data(self, data):
        if self.is_warming_up:
            return

        self.day_count += 1
        if self.day_count % self.step != 0:
            return

        date_str = str(self.time.date())

        spy_price = float(self.securities[self.spy].price) if self.spy in self.securities else None
        vix_price = None
        if self.has_vix:
            try:
                vix_price = float(self.securities[self.vix].price)
            except:
                pass

        row = {"date": date_str, "spy": spy_price, "vix": vix_price}

        # SPEED: fetch history once for the largest window, slice for smaller ones
        max_w = max(self.windows)
        syms  = list(self.symbols.values())
        try:
            full_hist = self.history(syms, max_w, Resolution.DAILY)
        except Exception as e:
            for w in self.windows:
                for k in ["sv2sv1","rank","eta","pr","n_valid"]:
                    row[f"w{w}_{k}"] = None
                self.gate_log.append(f"{date_str} w{w}: FAIL (history error {str(e)[:40]})")
            self.results.append(row)
            return

        if full_hist.empty:
            for w in self.windows:
                for k in ["sv2sv1","rank","eta","pr","n_valid"]:
                    row[f"w{w}_{k}"] = None
                self.gate_log.append(f"{date_str} w{w}: FAIL (empty history)")
            self.results.append(row)
            return

        for w in self.windows:
            try:
                result = self._compute_from_history(full_hist, w)
                if result:
                    sv = round(result["sv2sv1"], 4)
                    row[f"w{w}_sv2sv1"] = sv
                    row[f"w{w}_rank"]   = round(result["rank"], 2)
                    row[f"w{w}_eta"]    = round(result["eta"],  4)
                    row[f"w{w}_pr"]     = round(result["pr"],   3)
                    row[f"w{w}_n_valid"] = result["n_valid"]
                    self.gate_log.append(f"{date_str} w{w}: OK n={result['n_valid']}")
                    self._check_rolling_threshold(w, sv, date_str)
                else:
                    for k in ["sv2sv1","rank","eta","pr","n_valid"]:
                        row[f"w{w}_{k}"] = None
                    self.gate_log.append(f"{date_str} w{w}: FAIL (no result)")
            except Exception as e:
                for k in ["sv2sv1","rank","eta","pr","n_valid"]:
                    row[f"w{w}_{k}"] = None
                self.gate_log.append(f"{date_str} w{w}: ERROR {str(e)[:80]}")

        self.results.append(row)

        if len(self.results) % 50 == 0:
            self.log(f"Processed {len(self.results)} dates through {date_str}")

    # ------------------------------------------------------------------
    def _compute_from_history(self, full_hist, window):
        """Slice pre-fetched 90d history to `window` days and run Fisher pipeline."""
        closes = full_hist["close"].unstack(level=0)

        # Slice to requested window (take last `window` rows)
        if len(closes) > window:
            closes = closes.iloc[-window:]

        closes = closes.dropna(axis=1, thresh=int(window * 0.8))
        if closes.shape[1] < 50:
            return None

        returns = np.log(closes / closes.shift(1)).dropna()
        if len(returns) < window * 0.7:
            return None

        corr = returns.corr().values
        n    = corr.shape[0]
        np.fill_diagonal(corr, 1.0)
        corr = np.nan_to_num(corr, nan=0.0)

        adj = self._build_knn(corr, min(self.k_nn, n - 1))
        C_r = self._compute_cr(corr, adj, n)
        if C_r is None or len(C_r) < 3:
            return None

        return self._fisher(C_r, adj, n)

    # ------------------------------------------------------------------
    def _check_rolling_threshold(self, window, sv, date_str):
        """
        Maintain rolling 252-day (≈50 periods) mean+2σ threshold per window.
        Log every UP crossing with crisis proximity info. Append sv AFTER check.
        """
        hist = self.sv_history[window]

        if len(hist) >= self.rolling_min_periods:
            arr    = np.array(hist)
            r_mean = float(np.mean(arr))
            r_std  = float(np.std(arr))
            thresh = r_mean + 2.0 * r_std

            was_above = self.above_thresh[window]
            is_above  = sv > thresh

            if is_above and not was_above:
                # ── Crossed UP ──────────────────────────────────────────
                msg = (f"[ROLL-UP w{window}] {date_str}  "
                       f"SV={sv:.4f} > thresh={thresh:.4f}  "
                       f"(roll_mean={r_mean:.4f} roll_std={r_std:.4f} n={len(arr)})")
                self.log(msg)
                # Check proximity to every reference crisis
                from datetime import date as dt
                cross_dt = dt.fromisoformat(date_str)
                for crisis_name, crisis_str in self.crisis_dates.items():
                    crisis_dt = dt.fromisoformat(crisis_str)
                    delta = (crisis_dt - cross_dt).days
                    if 0 <= delta <= 60:
                        marker = "*** P3B-1 CANDIDATE (<= 30d) ***" if delta <= 30 else f"({delta}d before)"
                        self.log(f"  --> {delta}d before {crisis_name}  {marker}")
                self.rolling_crossings[window].append(("UP", date_str, sv, thresh))

            elif not is_above and was_above:
                # ── Crossed DOWN ────────────────────────────────────────
                self.log(f"[ROLL-DN w{window}] {date_str}  "
                         f"SV={sv:.4f} < thresh={thresh:.4f}")
                self.rolling_crossings[window].append(("DOWN", date_str, sv, thresh))

            self.above_thresh[window] = is_above

        # Always append AFTER the check (threshold based on history BEFORE this point)
        hist.append(sv)

    # ------------------------------------------------------------------
    def _build_knn(self, corr, k):
        n   = corr.shape[0]
        adj = np.zeros((n, n), dtype=bool)
        for i in range(n):
            row    = corr[i].copy()
            row[i] = -np.inf
            top_k  = np.argsort(row)[-k:]
            for j in top_k:
                adj[i, j] = True
                adj[j, i] = True
        return adj

    def _compute_cr(self, corr, adj, n):
        adj_sp  = csr_matrix(adj.astype(float))
        samples = np.random.choice(n, min(self.n_cr_samples, n), replace=False)
        dist    = shortest_path(adj_sp, method="D", indices=samples, directed=False)

        max_r = int(np.nanmax(dist[np.isfinite(dist)]))
        if max_r < 3:
            return None
        max_r = min(max_r, 15)

        C_r    = np.zeros(max_r + 1)
        counts = np.zeros(max_r + 1)

        for idx_i, src in enumerate(samples):
            for dst in range(n):
                if src == dst:
                    continue
                d = dist[idx_i, dst]
                if np.isfinite(d) and int(d) <= max_r:
                    C_r[int(d)]    += corr[src, dst]
                    counts[int(d)] += 1

        mask = counts > 0
        C_r[mask] /= counts[mask]
        return C_r

    def _fisher(self, C_r, adj, n):
        adj_sp    = csr_matrix(adj.astype(float))
        neighbors = {i: list(np.where(adj[i])[0]) for i in range(n)}
        samples   = np.random.choice(n, min(self.n_fisher_samples, n), replace=False)

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

            p_v0 = self._kernel(C_r, dists[0], n)
            if p_v0 is None:
                continue

            k      = len(nbrs)
            scores = np.zeros((k, n))
            valid  = True
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

            weighted = scores * np.sqrt(p_v0)[np.newaxis, :]
            F  = weighted @ weighted.T
            sv = np.linalg.svd(F, compute_uv=False)
            if sv[0] < 1e-15:
                continue
            sv = sv / sv[0]

            gaps = sv[:-1] / (sv[1:] + 1e-15)
            rank = int(np.argmax(gaps) + 1)
            eta  = float(sv[min(rank, len(sv) - 1)] / (sv[rank - 1] + 1e-15))
            pr   = float((np.sum(sv)) ** 2 / (np.sum(sv ** 2) + 1e-15))
            s2s1 = float(sv[1] / (sv[0] + 1e-15)) if len(sv) > 1 else 0.0

            ranks.append(rank)
            etas.append(eta)
            prs.append(pr)
            sv2sv1s.append(s2s1)

        if len(ranks) == 0:
            return None

        return {
            "sv2sv1": float(np.mean(sv2sv1s)),
            "rank":   float(np.mean(ranks)),
            "eta":    float(np.mean(etas)),
            "pr":     float(np.mean(prs)),
            "n_valid": len(ranks),
        }

    def _kernel(self, C_r, distances, n):
        weights = np.zeros(n)
        for u in range(n):
            d = distances[u]
            if np.isfinite(d) and int(d) < len(C_r):
                weights[u] = abs(C_r[int(d)])
        total = np.sum(weights)
        if total < 1e-15:
            return None
        return weights / total

    # ------------------------------------------------------------------
    def on_end_of_algorithm(self):
        if not self.results:
            self.log("No results to save")
            return

        # Save CSV
        fields = ["date", "spy", "vix"]
        for w in self.windows:
            fields += [f"w{w}_sv2sv1", f"w{w}_rank", f"w{w}_eta", f"w{w}_pr", f"w{w}_n_valid"]

        buf    = io.StringIO()
        writer = csv.DictWriter(buf, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        for row in self.results:
            writer.writerow(row)
        self.object_store.save("phase3b_results_csv", buf.getvalue())
        self.log(f"Saved {len(self.results)} rows to ObjectStore:phase3b_results_csv")

        gate_str = "\n".join(self.gate_log)
        self.object_store.save("phase3b_gate_log", gate_str)
        self.log(f"Saved {len(self.gate_log)} gate entries")

        # ── Summary stats ─────────────────────────────────────────────
        self.log("=" * 72)
        self.log("PHASE 3B v3: FISHER DIAGNOSTIC — END OF RUN SUMMARY")
        self.log("=" * 72)
        self.log(f"Period: {self.results[0]['date']} to {self.results[-1]['date']}")
        self.log(f"Total computation dates: {len(self.results)}")
        self.log(f"Universe: {len(self.symbols)} assets resolved")

        for w in self.windows:
            sv_vals  = [r.get(f"w{w}_sv2sv1") for r in self.results if r.get(f"w{w}_sv2sv1") is not None]
            rk_vals  = [r.get(f"w{w}_rank")   for r in self.results if r.get(f"w{w}_rank")   is not None]
            eta_vals = [r.get(f"w{w}_eta")    for r in self.results if r.get(f"w{w}_eta")    is not None]
            pr_vals  = [r.get(f"w{w}_pr")     for r in self.results if r.get(f"w{w}_pr")     is not None]
            nv_vals  = [r.get(f"w{w}_n_valid") for r in self.results if r.get(f"w{w}_n_valid") is not None]
            null_ct  = sum(1 for r in self.results if r.get(f"w{w}_sv2sv1") is None)

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
                self.log(f"  N_valid: mean={np.mean(nv_vals):.1f}  "
                         f"min={np.min(nv_vals)}  max={np.max(nv_vals)}")

        # ── VIX correlation ───────────────────────────────────────────
        sv90 = [r.get("w90_sv2sv1") for r in self.results if r.get("w90_sv2sv1") is not None]
        vix_vals = [r.get("vix") for r in self.results if r.get("w90_sv2sv1") is not None and r.get("vix") is not None]
        # Pair them correctly
        paired_sv, paired_vix = [], []
        for r in self.results:
            if r.get("w90_sv2sv1") is not None and r.get("vix") is not None:
                paired_sv.append(r["w90_sv2sv1"])
                paired_vix.append(r["vix"])
        if len(paired_sv) > 10:
            self.log(f"--- Correlations ---")
            self.log(f"  SV2/SV1(90d) vs VIX: r = {float(np.corrcoef(paired_sv, paired_vix)[0,1]):.3f}")

        # ── Gate summary ──────────────────────────────────────────────
        ok_ct   = sum(1 for g in self.gate_log if "OK"    in g)
        fail_ct = sum(1 for g in self.gate_log if "FAIL"  in g)
        err_ct  = sum(1 for g in self.gate_log if "ERROR" in g)
        self.log(f"--- Gate Summary ---")
        self.log(f"  OK: {ok_ct}  FAIL: {fail_ct}  ERROR: {err_ct}")

        # ── Rolling threshold crossing summary (P3B-1 kill test) ─────
        self.log("=" * 72)
        self.log("ROLLING THRESHOLD CROSSING SUMMARY  (P3B-1 KILL TEST)")
        self.log(f"Rolling window: {self.rolling_n} periods  "
                 f"(~{self.rolling_n * self.step} trading days ≈ 252d)")
        self.log(f"Crisis reference dates: {self.crisis_dates}")

        from datetime import date as dt

        p3b1_pass = False  # will set True if any UP crossing is within 30d of a crisis

        for w in self.windows:
            ups = [(d, v, t) for (direction, d, v, t) in self.rolling_crossings[w]
                   if direction == "UP"]
            dns = [(d, v, t) for (direction, d, v, t) in self.rolling_crossings[w]
                   if direction == "DOWN"]
            self.log(f"--- w{w}: {len(ups)} UP crossings, {len(dns)} DOWN crossings ---")
            for cross_date_str, sv_val, thresh_val in ups:
                cross_dt = dt.fromisoformat(cross_date_str)
                # Find nearest subsequent DOWN crossing (episode end)
                end_str = "ongoing"
                for (dd, dv, dt2) in dns:
                    if dd > cross_date_str:
                        end_str = dd
                        break
                self.log(f"  UP  {cross_date_str} → {end_str}  "
                         f"SV={sv_val:.4f} thresh={thresh_val:.4f}")
                # Crisis proximity
                for crisis_name, crisis_str in self.crisis_dates.items():
                    crisis_dt_obj = dt.fromisoformat(crisis_str)
                    delta = (crisis_dt_obj - cross_dt).days
                    if 0 <= delta <= 90:
                        verdict = ""
                        if delta <= 30:
                            verdict = "  *** P3B-1 PASS ***"
                            p3b1_pass = True
                        elif delta <= 60:
                            verdict = "  (within 60d)"
                        self.log(f"    {delta:3d}d before {crisis_name}{verdict}")

        self.log("--- P3B-1 KILL TEST VERDICT ---")
        if p3b1_pass:
            self.log("  RESULT: PASS — at least one rolling UP crossing within 30d of a crisis")
        else:
            self.log("  RESULT: FAIL — no rolling UP crossing found within 30d of any crisis")
            self.log("  (check individual crossing dates above against 60-90d windows)")

        self.log("=" * 72)
