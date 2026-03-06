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

        from datetime import date as dt, timedelta

        gate_str = "\n".join(self.gate_log)
        self.object_store.save("phase3b_gate_log", gate_str)

        # ── Summary stats ─────────────────────────────────────────────
        self.log("=" * 72)
        self.log("PHASE 3B v3: FISHER DIAGNOSTIC — END OF RUN SUMMARY")
        self.log("=" * 72)
        self.log(f"Period: {self.results[0]['date']} to {self.results[-1]['date']}")
        self.log(f"Total computation dates: {len(self.results)}")
        self.log(f"Universe: {len(self.symbols)} assets resolved")

        for w in self.windows:
            sv_vals  = [r[f"w{w}_sv2sv1"] for r in self.results if r.get(f"w{w}_sv2sv1") is not None]
            rk_vals  = [r[f"w{w}_rank"]   for r in self.results if r.get(f"w{w}_rank")   is not None]
            eta_vals = [r[f"w{w}_eta"]    for r in self.results if r.get(f"w{w}_eta")    is not None]
            pr_vals  = [r[f"w{w}_pr"]     for r in self.results if r.get(f"w{w}_pr")     is not None]
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

        # ── VIX correlation ───────────────────────────────────────────
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
        self.log(f"--- Gate Summary: OK={ok_ct}  FAIL={fail_ct}  ERROR={err_ct} ---")

        # ══════════════════════════════════════════════════════════════
        # P3B-1  ROLLING THRESHOLD CROSSING SUMMARY
        # ══════════════════════════════════════════════════════════════
        self.log("=" * 72)
        self.log("P3B-1: ROLLING THRESHOLD CROSSINGS")
        self.log(f"Rolling window: ~{self.rolling_n * self.step} trading days  "
                 f"Crises: {self.crisis_dates}")

        p3b1_pass = False

        for w in self.windows:
            ups = [(d, v, t) for (direction, d, v, t) in self.rolling_crossings[w] if direction == "UP"]
            dns = [(d, v, t) for (direction, d, v, t) in self.rolling_crossings[w] if direction == "DOWN"]
            self.log(f"--- w{w}: {len(ups)} UP  {len(dns)} DOWN ---")
            for cross_date_str, sv_val, thresh_val in ups:
                cross_dt = dt.fromisoformat(cross_date_str)
                end_str = next((dd for (dd, dv, dt2) in dns if dd > cross_date_str), "ongoing")
                self.log(f"  UP {cross_date_str}->{end_str} SV={sv_val:.4f} thr={thresh_val:.4f}")
                for crisis_name, crisis_str in self.crisis_dates.items():
                    delta = (dt.fromisoformat(crisis_str) - cross_dt).days
                    if 0 <= delta <= 90:
                        tag = " *** P3B-1 PASS ***" if delta <= 30 else f" ({delta}d)"
                        if delta <= 30:
                            p3b1_pass = True
                        self.log(f"    {delta:3d}d before {crisis_name}{tag}")

        self.log(f"P3B-1 VERDICT: {'PASS' if p3b1_pass else 'FAIL'}")

        # ══════════════════════════════════════════════════════════════
        # P3B-4  YEARLY SV2/SV1 PEAKS
        # ══════════════════════════════════════════════════════════════
        self.log("=" * 72)
        self.log("P3B-4: YEARLY SV2/SV1 PEAKS")
        self.log("year  w30_peak  w30_date    w60_peak  w60_date    w90_peak  w90_date")
        for year in range(2005, 2025):
            ys = str(year)
            yr = [r for r in self.results if r["date"][:4] == ys]
            parts = []
            for w in self.windows:
                vals = [(r["date"], r[f"w{w}_sv2sv1"]) for r in yr if r.get(f"w{w}_sv2sv1") is not None]
                if vals:
                    pd_, pv = max(vals, key=lambda x: x[1])
                    parts.append(f"{pv:.4f}  {pd_}")
                else:
                    parts.append("N/A       N/A       ")
            self.log(f"{ys}  {'  '.join(parts)}")

        # Direct P3B-4 verdict: 2008 peak vs 2018 peak per window
        for w in self.windows:
            sv08 = [r[f"w{w}_sv2sv1"] for r in self.results if r["date"][:4] == "2008" and r.get(f"w{w}_sv2sv1") is not None]
            sv18 = [r[f"w{w}_sv2sv1"] for r in self.results if r["date"][:4] == "2018" and r.get(f"w{w}_sv2sv1") is not None]
            if sv08 and sv18:
                p08, p18 = max(sv08), max(sv18)
                v = "PASS" if p08 > p18 else "FAIL"
                self.log(f"P3B-4 w{w}: 2008_peak={p08:.4f}  2018_peak={p18:.4f}  -> {v}")

        # ══════════════════════════════════════════════════════════════
        # P3B-5  SHORT-WINDOW LEADS LONG-WINDOW
        # ══════════════════════════════════════════════════════════════
        self.log("=" * 72)
        self.log("P3B-5: FIRST ROLLING-THRESHOLD CROSSING IN 90d PRE-CRISIS")
        self.log("(positive lead = shorter window fires earlier)")

        key_crises = [
            ("Bear_Stearns", "2008-03-14"),
            ("Lehman",       "2008-09-15"),
            ("COVID_bottom", "2020-03-23"),
        ]

        for crisis_name, crisis_str in key_crises:
            crisis_dt_obj = dt.fromisoformat(crisis_str)
            lookback = (crisis_dt_obj - timedelta(days=90)).isoformat()[:10]
            first = {}
            for w in self.windows:
                for direction, d, sv, thresh in self.rolling_crossings[w]:
                    if direction == "UP" and lookback <= d <= crisis_str:
                        first[w] = d
                        break
            self.log(f"--- {crisis_name} ({crisis_str}) ---")
            for w in self.windows:
                fc = first.get(w, "none")
                if fc != "none":
                    lead = (crisis_dt_obj - dt.fromisoformat(fc)).days
                    self.log(f"  w{w}: {fc}  ({lead}d before)")
                else:
                    self.log(f"  w{w}: no crossing in 90d window")
            fc30 = first.get(30)
            fc90 = first.get(90)
            if fc30 and fc90:
                lead = (dt.fromisoformat(fc90) - dt.fromisoformat(fc30)).days
                v = "PASS" if lead >= 5 else f"FAIL (lead={lead}d)"
                self.log(f"  P3B-5 w30 leads w90 by {lead}d -> {v}")
            elif fc30 and not fc90:
                self.log(f"  P3B-5 w30 fired, w90 did not -> PASS")
            else:
                self.log(f"  P3B-5 insufficient data")

        # ══════════════════════════════════════════════════════════════
        # P3B-2 / P3B-3  CRISIS WINDOW DATA TABLES + DIRECT VERDICTS
        # Logs CSV rows for 90d pre-crisis to 30d post-crisis
        # ══════════════════════════════════════════════════════════════
        self.log("=" * 72)
        self.log("P3B-2/P3B-3: CRISIS WINDOW DATA")
        self.log("date,w30_sv,w30_rk,w30_eta,w60_sv,w60_rk,w60_eta,w90_sv,w90_rk,w90_eta")

        crisis_windows = [
            ("Bear_Stearns", "2008-03-14", 90, 30),
            ("Lehman",       "2008-09-15", 90, 30),
            ("COVID_bottom", "2020-03-23", 90, 30),
        ]

        for crisis_name, crisis_str, pre_days, post_days in crisis_windows:
            crisis_dt_obj = dt.fromisoformat(crisis_str)
            w_start = (crisis_dt_obj - timedelta(days=pre_days)).isoformat()[:10]
            w_end   = (crisis_dt_obj + timedelta(days=post_days)).isoformat()[:10]
            self.log(f"# {crisis_name} {w_start} to {w_end}")
            for row in self.results:
                d = row["date"]
                if w_start <= d <= w_end:
                    self.log(
                        f"{d},"
                        f"{row.get('w30_sv2sv1','')},{row.get('w30_rank','')},{row.get('w30_eta','')},"
                        f"{row.get('w60_sv2sv1','')},{row.get('w60_rank','')},{row.get('w60_eta','')},"
                        f"{row.get('w90_sv2sv1','')},{row.get('w90_rank','')},{row.get('w90_eta','')}"
                    )

            # ── P3B-2: Rank trend in 60d pre-crisis ──────────────────
            pre60_start = (crisis_dt_obj - timedelta(days=60)).isoformat()[:10]
            pre_rows = [r for r in self.results if pre60_start <= r["date"] <= crisis_str]
            for w in [30, 90]:
                rk = [r[f"w{w}_rank"] for r in pre_rows if r.get(f"w{w}_rank") is not None]
                if len(rk) >= 4:
                    n = len(rk)
                    early = float(np.mean(rk[:n // 2]))
                    late  = float(np.mean(rk[n // 2:]))
                    delta_rk = late - early
                    v = "PASS" if delta_rk > 0 else "FAIL"
                    self.log(f"P3B-2 {crisis_name} w{w}: rank early={early:.2f} late={late:.2f} "
                             f"delta={delta_rk:+.2f} -> {v}")

            # ── P3B-3: η extremum within ±30d of crisis ──────────────
            eta30_start = (crisis_dt_obj - timedelta(days=30)).isoformat()[:10]
            eta30_end   = (crisis_dt_obj + timedelta(days=30)).isoformat()[:10]
            surr_start  = (crisis_dt_obj - timedelta(days=90)).isoformat()[:10]
            surr_end    = (crisis_dt_obj + timedelta(days=90)).isoformat()[:10]
            for w in [30, 90]:
                eta_win  = [r[f"w{w}_eta"] for r in self.results
                            if eta30_start <= r["date"] <= eta30_end and r.get(f"w{w}_eta") is not None]
                eta_surr = [r[f"w{w}_eta"] for r in self.results
                            if surr_start  <= r["date"] <= surr_end  and r.get(f"w{w}_eta") is not None]
                if eta_win and len(eta_surr) >= 4:
                    mu, sigma = float(np.mean(eta_surr)), float(np.std(eta_surr))
                    emin, emax = min(eta_win), max(eta_win)
                    extremum = (emin < mu - sigma) or (emax > mu + sigma)
                    v = "PASS" if extremum else "FAIL"
                    self.log(f"P3B-3 {crisis_name} w{w}: eta_win=[{emin:.4f},{emax:.4f}] "
                             f"surr={mu:.4f}+/-{sigma:.4f} -> {v}")

        # ══════════════════════════════════════════════════════════════
        # FINAL VERDICT SUMMARY
        # ══════════════════════════════════════════════════════════════
        self.log("=" * 72)
        self.log("PREDICTION VERDICT SUMMARY")
        self.log("(search log for P3B-1/P3B-2/P3B-3/P3B-4/P3B-5 PASS/FAIL)")
        self.log("=" * 72)
