from AlgorithmImports import *
import numpy as np
from collections import deque
import json
import math

# ══════════════════════════════════════════════════════════════════════
# GICS Sector Map — ~250 common S&P 500 tickers
# ══════════════════════════════════════════════════════════════════════
_GICS = {
    'Technology': [
        'AAPL','MSFT','GOOG','GOOGL','AMZN','META','NVDA','AVGO','ADBE','CRM',
        'CSCO','ACN','TXN','INTC','ORCL','IBM','QCOM','AMD','AMAT','INTU',
        'MU','ADI','LRCX','KLAC','SNPS','CDNS','NOW','PANW','FTNT','MRVL',
        'NXPI','ON','MPWR','GEN','CTSH','IT','ANSS','KEYS','ZBRA','TDY',
    ],
    'Financials': [
        'JPM','BAC','WFC','GS','MS','C','USB','PNC','BK','AXP','MET','PRU',
        'ALL','TRV','AIG','CME','ICE','SPGI','MSCI','BLK','SCHW','CB','AON',
        'MMC','AFL','CINF','HBAN','KEY','CFG','RF','FITB','MTB','NTRS','STT',
        'BRO','RJF','SIVB','ZION','CMA','SBNY','FRC',
    ],
    'HealthCare': [
        'JNJ','UNH','PFE','MRK','LLY','TMO','ABT','MDT','BMY','AMGN','GILD',
        'CVS','CI','ELV','HUM','ISRG','SYK','BDX','BSX','DXCM','REGN','VRTX',
        'ZTS','IDXX','IQV','A','BAX','EW','HOLX','ALGN','TECH','MTD','WAT',
    ],
    'ConsumerDisc': [
        'HD','MCD','NKE','SBUX','TGT','LOW','TJX','ROST','DHI','LEN','GM','F',
        'TSLA','BKNG','EBAY','MAR','HLT','YUM','DPZ','BBY','POOL','GPC',
        'APTV','BWA','CCL','RCL','MGM','WYNN','DRI','CMG',
    ],
    'ConsumerStap': [
        'PG','KO','PEP','WMT','COST','CL','GIS','SYY','MO','PM','MNST','EL',
        'STZ','KMB','CHD','K','SJM','HRL','ADM','KR','HSY','MKC','CPB','CAG',
        'CLX','BF.B','TAP',
    ],
    'Industrials': [
        'GE','CAT','HON','MMM','UPS','RTX','LMT','DE','EMR','ITW','GD','NSC',
        'CSX','WM','ETN','ROK','FDX','BA','GWW','SWK','PH','ODFL','CMI','IR',
        'DOV','AME','FAST','XYL','TT','CARR','OTIS','GNRC','PCAR','J','LHX',
        'HWM','TXT','LDOS','NOC','HII',
    ],
    'Energy': [
        'XOM','CVX','COP','SLB','EOG','PSX','VLO','MPC','OXY','HAL','DVN',
        'FANG','APA','BKR','HES','CTRA','MRO','TRGP','WMB','KMI','OKE',
    ],
    'Utilities': [
        'NEE','DUK','SO','D','AEP','EXC','SRE','ED','ES','WEC','XEL','PPL',
        'DTE','ETR','ATO','CMS','CNP','FE','EVRG','NI','PNW','AES',
    ],
    'RealEstate': [
        'PLD','AMT','EQIX','CCI','PSA','O','WELL','SPG','DLR','AVB','EQR',
        'VTR','ARE','BXP','SLG','REG','FRT','KIM','UDR','ESS','MAA','CPT',
    ],
    'CommServices': [
        'VZ','T','CMCSA','DIS','NFLX','TMUS','CHTR','EA','TTWO','ATVI','WBD',
        'PARA','FOX','FOXA','LYV','IPG','OMC','MTCH',
    ],
    'Materials': [
        'LIN','APD','SHW','ECL','NEM','FCX','NUE','VMC','MLM','DD','DOW','CE',
        'PPG','ALB','EMN','IFF','FMC','CF','MOS','PKG','IP','SEE','AVY','BLL',
    ],
}
GICS_MAP = {}
for _sect, _tickers in _GICS.items():
    for _t in _tickers:
        GICS_MAP[_t] = _sect

# ══════════════════════════════════════════════════════════════════════
# Algorithm
# ══════════════════════════════════════════════════════════════════════
class FisherRegimeAttractor(QCAlgorithm):

    # --- tunables ---
    WINDOWS         = [30, 90]
    UNIVERSE_SIZE   = 200
    K_NN_MARKET     = 10
    N_FISHER_SAMPLES= 20
    N_CR_SAMPLES    = 30
    MIN_CLUSTER_SIZE= 12
    MAX_CLUSTERS    = 15
    ZSCORE_LOOKBACK = 52      # weeks (1 year)

    # ──────────────────────────────────────────────────────────────
    # Lifecycle
    # ──────────────────────────────────────────────────────────────
    def initialize(self):
        self.set_start_date(2006, 1, 1)
        self.set_end_date(2024, 12, 31)
        self.set_cash(100000)

        # Dynamic universe
        self._universe_set = False
        self.add_universe(self.coarse_filter)
        self.universe_settings.resolution = Resolution.DAILY

        # Reference instruments
        spy = self.add_equity("SPY", Resolution.DAILY)
        spy.set_data_normalization_mode(DataNormalizationMode.ADJUSTED)
        self.spy = spy.symbol

        self.vix_symbol = None
        try:
            vix = self.add_data(CBOE, "VIX", Resolution.DAILY)
            self.vix_symbol = vix.symbol
        except:
            self.log("[INIT] VIX via CBOE not available")

        # Weekly schedule
        self.schedule.on(
            self.date_rules.every(DayOfWeek.MONDAY),
            self.time_rules.after_market_open("SPY", 30),
            self.compute_regime)

        # ── Charts ──
        c1 = Chart("Fisher Regime Attractor")
        c1.add_series(Series("FRA 90d",   SeriesType.LINE, 0))
        c1.add_series(Series("FRA 30d",   SeriesType.LINE, 0))
        c1.add_series(Series("Velocity",  SeriesType.LINE, 1))
        c1.add_series(Series("Momentum",  SeriesType.LINE, 1))
        self.add_chart(c1)

        c2 = Chart("Market SV2/SV1")
        c2.add_series(Series("SV2/SV1 90d", SeriesType.LINE, 0))
        c2.add_series(Series("SV2/SV1 30d", SeriesType.LINE, 0))
        self.add_chart(c2)

        c3 = Chart("Cluster Diagnostics")
        c3.add_series(Series("N Clusters", SeriesType.LINE, 0))
        c3.add_series(Series("Max Heat",   SeriesType.LINE, 1))
        c3.add_series(Series("Mean Heat",  SeriesType.LINE, 1))
        self.add_chart(c3)

        c4 = Chart("VIX")
        c4.add_series(Series("VIX", SeriesType.LINE, 0))
        self.add_chart(c4)

        # ── State ──
        self.sv_history   = {w: [] for w in self.WINDOWS}
        self.fra_hist_90  = []
        self.fra_hist_30  = []
        self.prev_fra_90  = None
        self.prev_regime  = None
        self.cluster_heat_buffers = {}   # cid → deque
        self.prev_labels  = None
        self.prev_symbols = None
        self.step = 0

        # end-of-run accumulators
        self.run_dates    = []
        self.run_fra90    = []
        self.run_fra30    = []
        self.run_vel      = []
        self.run_mom      = []
        self.run_regimes  = []
        self.run_vix      = []
        self.run_nclusters= []
        self.run_maxheat  = []
        self.run_nassets  = []
        self.asset_top10_count = {}      # ticker → count of top-10 appearances
        self.asset_composite_sum = {}    # ticker → cumulative |composite|

        self.set_warm_up(timedelta(days=400))

    # ──────────────────────────────────────────────────────────────
    # Universe selection (quarterly rebalance by dollar volume)
    # ──────────────────────────────────────────────────────────────
    def coarse_filter(self, coarse):
        if self._universe_set and (self.time.month % 3 != 1 or self.time.day > 7):
            return Universe.UNCHANGED
        self._universe_set = True
        filtered = sorted(
            [x for x in coarse if x.has_fundamental_data and x.price > 5.0],
            key=lambda x: x.dollar_volume, reverse=True)
        selected = [x.symbol for x in filtered[:self.UNIVERSE_SIZE]]
        self.log(f"[UNIVERSE] {self.time.date()}: selected {len(selected)} assets "
                 f"(top dollar volume)")
        return selected

    # ──────────────────────────────────────────────────────────────
    # Main weekly computation
    # ──────────────────────────────────────────────────────────────
    def compute_regime(self):
        if self.is_warming_up:
            return

        active = [s for s in self.active_securities.keys()
                  if s != self.spy
                  and (self.vix_symbol is None or s != self.vix_symbol)
                  and self.securities[s].price > 0]
        if len(active) < 50:
            return

        self.step += 1
        date = self.time

        # VIX
        vix_val = float('nan')
        if self.vix_symbol and self.securities[self.vix_symbol].has_data:
            try:
                vix_val = float(self.securities[self.vix_symbol].price)
            except:
                pass

        fra_map     = {}
        sv_map      = {}
        clusters_90 = None
        heats_90    = {}
        scores_90   = []
        names_90    = []
        n_used      = 0

        for window in self.WINDOWS:
            try:
                history = self.history(active, window, Resolution.DAILY)
                if history.empty:
                    continue
                closes = history['close'].unstack(level=0)
                closes = closes.dropna(axis=1, thresh=int(window * 0.8))
                if closes.shape[1] < 50:
                    continue

                returns = np.log(closes / closes.shift(1)).dropna()
                if len(returns) < window * 0.6:
                    continue

                corr = returns.corr().values
                asset_names = [str(c).split()[0] for c in returns.columns]
                n = len(asset_names)
                np.fill_diagonal(corr, 1.0)
                corr = np.nan_to_num(corr, nan=0.0)

                # Market-level Fisher SV2/SV1
                mkt_sv = self._market_fisher(corr, n)
                if mkt_sv is None:
                    continue
                sv_map[window] = mkt_sv

                # FRA
                fra = self._compute_fra(mkt_sv, window)
                fra_map[window] = fra

                # Plot raw SV2/SV1
                self.plot("Market SV2/SV1", f"SV2/SV1 {window}d", mkt_sv)

                # 90d window: clustering + heat + asset temps
                if window == 90:
                    n_used = n
                    names_90 = asset_names
                    try:
                        clusters_90 = self._spectral_clusters(corr, n)
                        heats_90 = self._cluster_heats(corr, clusters_90, n,
                                                        asset_names)
                        scores_90 = self._asset_temperatures(
                            corr, clusters_90, heats_90, n,
                            asset_names, fra)
                    except Exception as e:
                        self.log(f"[WARN] Clustering/heat error {date.date()}: "
                                 f"{str(e)[:150]}")
            except Exception as e:
                self.log(f"[ERROR] w={window} {date.date()}: {str(e)[:200]}")

        # ── Derived signals ──
        fra_90 = fra_map.get(90, float('nan'))
        fra_30 = fra_map.get(30, float('nan'))

        if not np.isnan(fra_90):
            self.fra_hist_90.append(fra_90)
            self.plot("Fisher Regime Attractor", "FRA 90d", fra_90)

        if not np.isnan(fra_30):
            self.fra_hist_30.append(fra_30)
            self.plot("Fisher Regime Attractor", "FRA 30d", fra_30)

        vel = 0.0
        if self.prev_fra_90 is not None and not np.isnan(fra_90):
            vel = fra_90 - self.prev_fra_90
            self.plot("Fisher Regime Attractor", "Velocity", vel)

        mom = float('nan')
        if not np.isnan(fra_30) and not np.isnan(fra_90):
            mom = fra_30 - fra_90
            self.plot("Fisher Regime Attractor", "Momentum", mom)

        n_cl = len(np.unique(clusters_90)) if clusters_90 is not None else 0
        max_heat = max((abs(h) for h in heats_90.values()), default=0.0)
        mean_heat = float(np.mean(list(heats_90.values()))) if heats_90 else 0.0

        if n_cl > 0:
            self.plot("Cluster Diagnostics", "N Clusters", n_cl)
            self.plot("Cluster Diagnostics", "Max Heat", max_heat)
            self.plot("Cluster Diagnostics", "Mean Heat", mean_heat)
        if not np.isnan(vix_val):
            self.plot("VIX", "VIX", vix_val)

        # ── Regime classification ──
        regime = self._regime_name(fra_90)

        # ── Threshold alerts ──
        if self.prev_fra_90 is not None and not np.isnan(fra_90):
            for thr in [-0.76, -0.50, -0.30, 0.30, 0.50, 0.76]:
                prev = self.prev_fra_90
                if (prev > thr >= fra_90) or (prev < thr <= fra_90):
                    d = "v" if fra_90 < prev else "^"
                    old_r = self._regime_name(prev)
                    self.log(f"  *** FRA CROSSED {thr:+.2f} {d} "
                             f"[{date.date()}] — {old_r} -> {regime} ***")

        # ── Diagnostic log ──
        hottest_gics = "N/A"
        if heats_90 and clusters_90 is not None:
            hot_c = max(heats_90, key=lambda c: abs(heats_90[c]))
            members = np.where(clusters_90 == hot_c)[0]
            gics_cnt = {}
            for m in members:
                if m < len(names_90):
                    sec = GICS_MAP.get(names_90[m], 'Other')
                    gics_cnt[sec] = gics_cnt.get(sec, 0) + 1
            if gics_cnt:
                top_sec = max(gics_cnt, key=gics_cnt.get)
                hottest_gics = f"{top_sec}({gics_cnt[top_sec]})"

        top3_str = ""
        if scores_90:
            top3 = scores_90[:3]
            top3_str = ", ".join(
                f"{a['name']}={a['composite']:+.3f}" for a in top3)

        mom_str = f"{mom:+.4f}" if not np.isnan(mom) else "  N/A "
        self.log(
            f"[W{self.step:>4}] {date.date()} | "
            f"FRA={fra_90:+.4f} [{regime:>12}] | "
            f"v={vel:+.4f} m={mom_str} | "
            f"K={n_cl} Hmax={max_heat:.2f} [{hottest_gics}] | "
            f"N={n_used} VIX={vix_val:.1f} | "
            f"Top: {top3_str}")

        # ── Accumulate for end-of-run ──
        self.run_dates.append(str(date.date()))
        self.run_fra90.append(fra_90)
        self.run_fra30.append(fra_30)
        self.run_vel.append(vel)
        self.run_mom.append(mom if not np.isnan(mom) else 0.0)
        self.run_regimes.append(regime)
        self.run_vix.append(vix_val)
        self.run_nclusters.append(n_cl)
        self.run_maxheat.append(max_heat)
        self.run_nassets.append(n_used)

        for a in scores_90[:10]:
            t = a['name']
            self.asset_top10_count[t] = self.asset_top10_count.get(t, 0) + 1
            self.asset_composite_sum[t] = (
                self.asset_composite_sum.get(t, 0.0) + abs(a['composite']))

        if not np.isnan(fra_90):
            self.prev_fra_90 = fra_90
        self.prev_regime = regime

    # ──────────────────────────────────────────────────────────────
    # FRA Construction
    # ──────────────────────────────────────────────────────────────
    def _compute_fra(self, sv2sv1, window):
        hist = self.sv_history[window]
        hist.append(sv2sv1)
        if len(hist) < self.ZSCORE_LOOKBACK:
            return 0.0
        recent = np.array(hist[-self.ZSCORE_LOOKBACK:])
        m = np.mean(recent)
        s = np.std(recent)
        if s < 1e-10:
            return 0.0
        z = (sv2sv1 - m) / s
        return float(np.tanh(-z / 2.0))

    # ──────────────────────────────────────────────────────────────
    # Market-level Fisher SV2/SV1
    # ──────────────────────────────────────────────────────────────
    def _market_fisher(self, corr, n):
        k = min(self.K_NN_MARKET, n - 1)
        adj, nbrs = self._build_knn(corr, n, k)
        Cr = self._compute_Cr(corr, adj, nbrs, n,
                               min(self.N_CR_SAMPLES, n))
        if Cr is None or len(Cr) < 2:
            return None
        n_samp = min(self.N_FISHER_SAMPLES, n)
        indices = np.random.choice(n, n_samp, replace=False)
        sv_vals = []
        for v0 in indices:
            r = self._single_fim(Cr, adj, nbrs, v0, n)
            if r is not None:
                sv_vals.append(r)
        return float(np.mean(sv_vals)) if sv_vals else None

    # ──────────────────────────────────────────────────────────────
    # Spectral Clustering
    # ──────────────────────────────────────────────────────────────
    def _spectral_clusters(self, corr, n):
        # Distance
        dist = np.sqrt(np.clip(2.0 * (1.0 - corr), 0, 4))
        np.fill_diagonal(dist, 0)

        # Affinity (Gaussian kernel)
        med_d = np.median(dist[dist > 0])
        if med_d < 1e-10:
            med_d = 1.0
        A = np.exp(-dist**2 / (2.0 * med_d**2))
        np.fill_diagonal(A, 0)

        # Normalized Laplacian
        deg = A.sum(axis=1)
        deg_inv_sqrt = np.where(deg > 1e-10, 1.0 / np.sqrt(deg), 0.0)
        D_is = np.diag(deg_inv_sqrt)
        L_sym = np.eye(n) - D_is @ A @ D_is

        # Eigendecomposition
        eigvals, eigvecs = np.linalg.eigh(L_sym)

        # Eigengap → k
        max_k = min(self.MAX_CLUSTERS, n // self.MIN_CLUSTER_SIZE)
        if max_k < 2:
            return np.zeros(n, dtype=int)
        # Use eigenvalues 1..max_k (skip 0th which ≈ 0)
        gap_range = min(max_k, len(eigvals) - 2)
        if gap_range < 2:
            return np.zeros(n, dtype=int)
        gaps = np.diff(eigvals[1:gap_range + 2])
        k = int(np.argmax(gaps)) + 2
        k = max(2, min(k, max_k))

        # K-means on first k eigenvectors (skip eigvec 0)
        features = eigvecs[:, 1:k + 1].copy()
        norms = np.linalg.norm(features, axis=1, keepdims=True) + 1e-10
        features /= norms
        labels = self._kmeans(features, k)

        # Merge small clusters
        labels = self._merge_small(labels, corr)
        return labels

    def _kmeans(self, X, k, n_iter=30, n_init=3):
        n, d = X.shape
        best_labels = np.zeros(n, dtype=int)
        best_inertia = float('inf')
        for _ in range(n_init):
            # k-means++ init
            centers = np.empty((k, d))
            centers[0] = X[np.random.randint(n)]
            for c in range(1, k):
                dists = np.array([np.sum((X - centers[j])**2, axis=1)
                                  for j in range(c)])
                min_dists = np.min(dists, axis=0)
                total = min_dists.sum()
                if total < 1e-15:
                    centers[c] = X[np.random.randint(n)]
                else:
                    probs = min_dists / total
                    centers[c] = X[np.random.choice(n, p=probs)]
            for _ in range(n_iter):
                d2 = np.array([np.sum((X - centers[j])**2, axis=1)
                               for j in range(k)])
                labels = np.argmin(d2, axis=0)
                new_c = np.empty_like(centers)
                for j in range(k):
                    mask = labels == j
                    new_c[j] = X[mask].mean(axis=0) if mask.any() else X[np.random.randint(n)]
                if np.allclose(centers, new_c, atol=1e-8):
                    break
                centers = new_c
            inertia = sum(np.sum((X[labels == j] - centers[j])**2)
                          for j in range(k))
            if inertia < best_inertia:
                best_inertia = inertia
                best_labels = labels.copy()
        return best_labels

    def _merge_small(self, labels, corr):
        for _ in range(5):    # iterate in case of cascading merges
            uniq, counts = np.unique(labels, return_counts=True)
            small = uniq[counts < self.MIN_CLUSTER_SIZE]
            if len(small) == 0:
                break
            large = uniq[counts >= self.MIN_CLUSTER_SIZE]
            if len(large) == 0:
                break
            for s in small:
                s_members = np.where(labels == s)[0]
                best_c, best_r = -1, -2.0
                for lg in large:
                    lg_members = np.where(labels == lg)[0]
                    r = float(np.mean(corr[np.ix_(s_members, lg_members)]))
                    if r > best_r:
                        best_r = r
                        best_c = lg
                if best_c >= 0:
                    labels[s_members] = best_c
        # Relabel to contiguous
        for new_id, old_id in enumerate(np.unique(labels)):
            labels[labels == old_id] = new_id
        return labels

    # ──────────────────────────────────────────────────────────────
    # Cluster Heat
    # ──────────────────────────────────────────────────────────────
    def _cluster_heats(self, corr, labels, n, asset_names):
        heats = {}
        for c in np.unique(labels):
            members = np.where(labels == c)[0]
            nc = len(members)
            if nc < 7:
                heats[int(c)] = 0.0
                continue
            sub_corr = corr[np.ix_(members, members)]
            k_nn = max(4, nc // 3)
            sv = self._subgraph_fisher(sub_corr, nc, k_nn)
            if sv is None:
                heats[int(c)] = 0.0
                continue
            # z-score against this cluster's history
            cid = int(c)
            if cid not in self.cluster_heat_buffers:
                self.cluster_heat_buffers[cid] = deque(maxlen=104)
            buf = self.cluster_heat_buffers[cid]
            buf.append(sv)
            if len(buf) < 26:
                heats[cid] = 0.0
            else:
                arr = np.array(buf)[-self.ZSCORE_LOOKBACK:]
                m, s = np.mean(arr), np.std(arr)
                heats[cid] = float((sv - m) / (s + 1e-10))
        return heats

    # ──────────────────────────────────────────────────────────────
    # Asset Temperatures
    # ──────────────────────────────────────────────────────────────
    def _asset_temperatures(self, corr, labels, heats, n,
                             asset_names, fra):
        scores = []
        for c in np.unique(labels):
            members = np.where(labels == c)[0]
            nc = len(members)
            if nc < 7:
                continue
            sub_corr = corr[np.ix_(members, members)]
            k_nn = max(4, nc // 3)
            adj, nbrs = self._build_knn(sub_corr, nc, k_nn)
            Cr = self._compute_Cr(sub_corr, adj, nbrs, nc,
                                   min(nc, self.N_CR_SAMPLES))
            if Cr is None or len(Cr) < 2:
                continue
            # Per-asset SV2/SV1
            asset_svs = {}
            for i in range(nc):
                r = self._single_fim(Cr, adj, nbrs, i, nc)
                if r is not None:
                    asset_svs[i] = r
            if not asset_svs:
                continue
            mean_sv = float(np.mean(list(asset_svs.values())))
            std_sv  = max(float(np.std(list(asset_svs.values()))), 1e-10)
            c_heat  = heats.get(int(c), 0.0)

            for local_i, sv in asset_svs.items():
                global_i = members[local_i]
                if global_i >= len(asset_names):
                    continue
                micro = (sv - mean_sv) / std_sv
                composite = 0.4 * fra + 0.4 * c_heat + 0.2 * micro
                name = asset_names[global_i]
                scores.append({
                    'name':      name,
                    'composite': float(composite),
                    'fra':       float(fra),
                    'heat':      float(c_heat),
                    'micro':     float(micro),
                    'sv2sv1':    float(sv),
                    'cluster':   int(c),
                    'gics':      GICS_MAP.get(name, 'Other'),
                })
        scores.sort(key=lambda x: abs(x['composite']), reverse=True)
        return scores

    # ──────────────────────────────────────────────────────────────
    # Subgraph Fisher (for cluster heat — averaged SV2/SV1)
    # ──────────────────────────────────────────────────────────────
    def _subgraph_fisher(self, sub_corr, nc, k_nn):
        adj, nbrs = self._build_knn(sub_corr, nc, k_nn)
        Cr = self._compute_Cr(sub_corr, adj, nbrs, nc,
                               min(nc, self.N_CR_SAMPLES))
        if Cr is None or len(Cr) < 2:
            return None
        n_samp = min(self.N_FISHER_SAMPLES, nc)
        indices = np.random.choice(nc, n_samp, replace=False)
        svs = []
        for v0 in indices:
            r = self._single_fim(Cr, adj, nbrs, v0, nc)
            if r is not None:
                svs.append(r)
        return float(np.mean(svs)) if svs else None

    # ──────────────────────────────────────────────────────────────
    # Fisher Pipeline (verbatim Phase 3B)
    # ──────────────────────────────────────────────────────────────
    def _build_knn(self, corr, n, k):
        adj = np.zeros((n, n), dtype=bool)
        for i in range(n):
            row = corr[i].copy()
            row[i] = -np.inf
            top_k = np.argsort(row)[-k:]
            for j in top_k:
                adj[i, j] = True
                adj[j, i] = True
        nbrs = {i: list(np.where(adj[i])[0]) for i in range(n)}
        return adj, nbrs

    def _bfs(self, adj, src, n, max_d=15):
        dist = np.full(n, -1, dtype=int)
        dist[src] = 0
        queue = deque([src])
        while queue:
            u = queue.popleft()
            if dist[u] >= max_d:
                continue
            for v in range(n):
                if adj[u, v] and dist[v] < 0:
                    dist[v] = dist[u] + 1
                    queue.append(v)
        return dist

    def _compute_Cr(self, corr, adj, nbrs, n, n_samples):
        samples = np.random.choice(n, min(n_samples, n), replace=False)
        max_r = 0
        all_dists = {}
        for src in samples:
            d = self._bfs(adj, src, n)
            all_dists[src] = d
            valid = d[d > 0]
            if len(valid) > 0:
                mr = int(valid.max())
                max_r = max(max_r, mr)
        if max_r < 2:
            return None
        max_r = min(max_r, 10)
        Cr = np.zeros(max_r + 1)
        counts = np.zeros(max_r + 1)
        for src in samples:
            d = all_dists[src]
            for dst in range(n):
                if src == dst:
                    continue
                dd = d[dst]
                if 0 < dd <= max_r:
                    Cr[dd] += corr[src, dst]
                    counts[dd] += 1
        mask = counts > 0
        Cr[mask] /= counts[mask]
        return Cr

    def _single_fim(self, Cr, adj, nbrs, v0, n):
        nb = nbrs.get(v0, [])
        if len(nb) < 2:
            return None
        d_v0 = self._bfs(adj, v0, n)
        p_v0 = self._kernel(Cr, d_v0, n)
        if p_v0 is None:
            return None
        k = len(nb)
        score = np.zeros((k, n))
        eps = 1e-12
        for j, wj in enumerate(nb):
            d_wj = self._bfs(adj, wj, n)
            p_wj = self._kernel(Cr, d_wj, n)
            if p_wj is None:
                return None
            score[j] = np.log(p_wj + eps) - np.log(p_v0 + eps)
        weighted = score * np.sqrt(p_v0)[np.newaxis, :]
        FIM = weighted @ weighted.T
        sv = np.linalg.svd(FIM, compute_uv=False)
        sv = np.sort(sv)[::-1]
        if sv[0] < 1e-15:
            return None
        return float(sv[1] / sv[0]) if len(sv) > 1 else 0.0

    def _kernel(self, Cr, dist, n):
        w = np.zeros(n)
        for u in range(n):
            d = dist[u]
            if 0 <= d < len(Cr):
                w[u] = abs(Cr[d])
        total = w.sum()
        if total < 1e-15:
            return None
        return w / total

    # ──────────────────────────────────────────────────────────────
    # Regime name
    # ──────────────────────────────────────────────────────────────
    def _regime_name(self, fra):
        if np.isnan(fra):
            return "N/A"
        if fra < -0.76:  return "CRISIS"
        if fra < -0.50:  return "High Stress"
        if fra < -0.30:  return "Stress"
        if fra < -0.10:  return "Mild Stress"
        if fra <  0.10:  return "Equilibrium"
        if fra <  0.30:  return "Mild Disp"
        if fra <  0.50:  return "Dispersion"
        if fra <  0.76:  return "High Disp"
        return "EXTREME DISP"

    # ──────────────────────────────────────────────────────────────
    # End-of-algorithm: KPI summary + ObjectStore export
    # ──────────────────────────────────────────────────────────────
    def on_end_of_algorithm(self):
        N = len(self.run_dates)
        if N == 0:
            self.log("[END] No data collected.")
            return

        sep = "=" * 80
        self.log(sep)
        self.log("PHASE 3B-3: FISHER REGIME ATTRACTOR — END OF RUN SUMMARY")
        self.log(sep)
        self.log(f"Period: {self.run_dates[0]} to {self.run_dates[-1]}")
        self.log(f"Weekly computation dates: {N}")
        self.log(f"Windows: {self.WINDOWS}")
        self.log(f"Universe target: {self.UNIVERSE_SIZE}")
        self.log(f"Min cluster size: {self.MIN_CLUSTER_SIZE}")
        self.log(f"Fisher samples: {self.N_FISHER_SAMPLES}")

        # ── FRA Statistics ──
        fra90 = np.array([x for x in self.run_fra90 if not np.isnan(x)])
        fra30 = np.array([x for x in self.run_fra30 if not np.isnan(x)])

        self.log(f"\n--- FRA Statistics ---")
        if len(fra90) > 0:
            min_i = int(np.argmin(fra90))
            max_i = int(np.argmax(fra90))
            valid_dates = [d for d, f in zip(self.run_dates, self.run_fra90)
                           if not np.isnan(f)]
            self.log(f"  FRA 90d: mean={np.mean(fra90):+.4f}  "
                     f"std={np.std(fra90):.4f}  "
                     f"min={np.min(fra90):+.4f} ({valid_dates[min_i]})  "
                     f"max={np.max(fra90):+.4f} ({valid_dates[max_i]})")
            self.log(f"    5th%={np.percentile(fra90, 5):+.4f}  "
                     f"25th%={np.percentile(fra90, 25):+.4f}  "
                     f"50th%={np.percentile(fra90, 50):+.4f}  "
                     f"75th%={np.percentile(fra90, 75):+.4f}  "
                     f"95th%={np.percentile(fra90, 95):+.4f}")
        if len(fra30) > 0:
            self.log(f"  FRA 30d: mean={np.mean(fra30):+.4f}  "
                     f"std={np.std(fra30):.4f}  "
                     f"min={np.min(fra30):+.4f}  "
                     f"max={np.max(fra30):+.4f}")

        vel = np.array(self.run_vel)
        mom = np.array(self.run_mom)
        self.log(f"  Velocity: mean={np.mean(vel):+.4f}  "
                 f"std={np.std(vel):.4f}")
        self.log(f"  Momentum: mean={np.mean(mom):+.4f}  "
                 f"std={np.std(mom):.4f}")

        # ── Regime Distribution ──
        self.log(f"\n--- Regime Distribution ---")
        regime_counts = {}
        for r in self.run_regimes:
            regime_counts[r] = regime_counts.get(r, 0) + 1
        for r_name in ["CRISIS", "High Stress", "Stress", "Mild Stress",
                        "Equilibrium", "Mild Disp", "Dispersion",
                        "High Disp", "EXTREME DISP", "N/A"]:
            cnt = regime_counts.get(r_name, 0)
            if cnt > 0:
                self.log(f"  {r_name:>14}: {cnt:>4} weeks ({100*cnt/N:.1f}%)")

        # ── Crisis Episodes (FRA < -0.50) ──
        self.log(f"\n--- Crisis Episodes (FRA < -0.50) ---")
        in_crisis = False
        c_start = None
        c_min = 0.0
        episode = 0
        for i in range(N):
            f = self.run_fra90[i]
            if np.isnan(f):
                continue
            if f < -0.50 and not in_crisis:
                in_crisis = True
                c_start = self.run_dates[i]
                c_min = f
            elif f < -0.50 and in_crisis:
                c_min = min(c_min, f)
            elif f >= -0.50 and in_crisis:
                in_crisis = False
                episode += 1
                dur = i - next(j for j in range(N)
                               if self.run_dates[j] == c_start)
                self.log(f"  {episode}. {c_start} to {self.run_dates[i]} "
                         f"({dur} wks, min FRA={c_min:+.3f})")
        if in_crisis:
            episode += 1
            self.log(f"  {episode}. {c_start} to {self.run_dates[-1]} "
                     f"(ongoing, min FRA={c_min:+.3f})")
        if episode == 0:
            self.log("  (none)")

        # ── Cluster Statistics ──
        ncl = [x for x in self.run_nclusters if x > 0]
        self.log(f"\n--- Cluster Statistics ---")
        if ncl:
            self.log(f"  Mean K: {np.mean(ncl):.1f}  "
                     f"Std: {np.std(ncl):.1f}  "
                     f"Min: {min(ncl)}  Max: {max(ncl)}")

        # ── VIX Correlation ──
        self.log(f"\n--- Market Indicator Correlations ---")
        vix_arr = np.array(self.run_vix)
        pairs = [(fra90, [v for v, f in zip(self.run_vix, self.run_fra90)
                          if not np.isnan(f) and not np.isnan(v)],
                  "FRA_90d vs VIX"),
                 (fra30, [v for v, f in zip(self.run_vix, self.run_fra30)
                          if not np.isnan(f) and not np.isnan(v)],
                  "FRA_30d vs VIX")]
        for arr1, vix_sub, label in pairs:
            if len(arr1) > 50 and len(vix_sub) >= len(arr1):
                r = np.corrcoef(arr1, np.array(vix_sub[:len(arr1)]))[0, 1]
                self.log(f"  {label}: r = {r:+.3f}")

        # ── Most Active Assets ──
        self.log(f"\n--- Most Active Assets (top-10 by total |composite|) ---")
        sorted_assets = sorted(self.asset_composite_sum.items(),
                               key=lambda x: x[1], reverse=True)[:10]
        for i, (ticker, total) in enumerate(sorted_assets):
            cnt = self.asset_top10_count.get(ticker, 0)
            avg = total / max(cnt, 1)
            gics = GICS_MAP.get(ticker, 'Other')
            self.log(f"  {i+1:>2}. {ticker:>6} [{gics:>14}]: "
                     f"top-10 in {cnt} weeks, "
                     f"avg|composite|={avg:.3f}")

        # ── Prediction Quick-Check ──
        self.log(f"\n--- Prediction Pre-Assessment ---")
        # P3B3-1: FRA < -0.3 within 60d before crisis
        crises_check = [("2008 GFC", "2008-09-15"),
                        ("2020 COVID", "2020-03-23")]
        for name, peak_str in crises_check:
            found = False
            for i in range(N):
                d = self.run_dates[i]
                f = self.run_fra90[i]
                if np.isnan(f):
                    continue
                if d <= peak_str and d >= peak_str[:4] + "-01-01":
                    if f < -0.3:
                        self.log(f"  P3B3-1 ({name}): FRA={f:+.3f} on {d} "
                                 f"(before {peak_str}) -> CANDIDATE PASS")
                        found = True
                        break
            if not found:
                self.log(f"  P3B3-1 ({name}): No FRA < -0.3 found "
                         f"before {peak_str}")

        self.log(sep)

        # ── ObjectStore Export ──
        output = {
            'dates':      self.run_dates,
            'fra_90d':    self.run_fra90,
            'fra_30d':    self.run_fra30,
            'velocity':   self.run_vel,
            'momentum':   self.run_mom,
            'regimes':    self.run_regimes,
            'vix':        self.run_vix,
            'n_clusters': self.run_nclusters,
            'max_heat':   self.run_maxheat,
            'n_assets':   self.run_nassets,
            'params': {
                'windows':    self.WINDOWS,
                'universe':   self.UNIVERSE_SIZE,
                'k_nn':       self.K_NN_MARKET,
                'n_samples':  self.N_FISHER_SAMPLES,
                'min_cluster':self.MIN_CLUSTER_SIZE,
                'max_clusters':self.MAX_CLUSTERS,
            }
        }
        try:
            self.object_store.save("fisher_regime_attractor_results",
                                    json.dumps(output, default=str))
            self.log(f"[END] Saved {N} results to ObjectStore")
        except Exception as e:
            self.log(f"[END] ObjectStore save failed: {str(e)[:100]}")
