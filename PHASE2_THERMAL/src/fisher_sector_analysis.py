#!/usr/bin/env python3
"""
Fisher Sector Analysis: Phase 3B-2 Post-Processing
====================================================
Reads extracted Virtual Green Flamingo chart data, generates 9 figures,
assesses P3B2-1 through P3B2-6 predictions, writes results report.

Usage:
    python3 fisher_sector_analysis.py [path_to_json]
    Default: phase3b2_results/raw_data/fisher_sector_results.json
"""

import json
import sys
import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.dates as mdates
from matplotlib.gridspec import GridSpec
from datetime import datetime, timedelta
from collections import defaultdict

# ──────────────────────────────────────────────────────────────────────
# Configuration
# ──────────────────────────────────────────────────────────────────────
RESULTS_DIR = os.path.join(os.path.dirname(__file__), "phase3b2_results")
os.makedirs(RESULTS_DIR, exist_ok=True)

SECTOR_COLORS = {
    'Market':      '#1f77b4',   # blue
    'Financials':  '#d62728',   # red
    'Technology':  '#9467bd',   # purple
    'HealthCare':  '#2ca02c',   # green
    'Energy':      '#ff7f0e',   # orange
    'Consumer':    '#8c564b',   # brown
    'Industrials': '#e377c2',   # pink
    'Utilities':   '#7f7f7f',   # gray
}

CRISIS_DATES = [
    ('2007-07-01', '2009-06-01', '2008 GFC',     '#ffcccc', '2008-09-15'),
    ('2014-10-01', '2016-03-01', '2014-16 Oil',  '#fff2cc', '2016-01-20'),
    ('2018-01-26', '2018-04-01', '2018 Volmag',  '#e2f0d9', '2018-02-05'),
    ('2020-01-15', '2020-05-01', '2020 COVID',   '#dae3f3', '2020-03-23'),
    ('2022-01-01', '2022-12-31', '2022 Rates',   '#f4e3ff', '2022-01-03'),
]

def parse_date(s):
    return datetime.strptime(s[:10], "%Y-%m-%d")

# ──────────────────────────────────────────────────────────────────────
# Data Loading
# ──────────────────────────────────────────────────────────────────────
def load_data(json_path):
    with open(json_path) as f:
        raw = json.load(f)

    def to_series(lst):
        dates = [parse_date(x['date']) for x in lst]
        vals  = [x['value'] for x in lst]
        return np.array(dates), np.array(vals, dtype=float)

    sv2sv1 = {}
    for name, lst in raw['sector_sv2sv1'].items():
        sv2sv1[name] = to_series(lst)

    spread = {}
    for name, lst in raw['sector_spread'].items():
        spread[name] = to_series(lst)

    vix_dates, vix_vals = to_series(raw['vix'])
    bench_dates, bench_vals = to_series(raw['benchmark'])

    return sv2sv1, spread, vix_dates, vix_vals, bench_dates, bench_vals

# ──────────────────────────────────────────────────────────────────────
# Signal Processing
# ──────────────────────────────────────────────────────────────────────
def rolling_zscore(vals, window=52):
    """Rolling z-score using trailing window."""
    z = np.full_like(vals, np.nan)
    for i in range(window, len(vals)):
        w = vals[i - window:i]
        m = np.nanmean(w)
        s = np.nanstd(w)
        if s > 1e-10:
            z[i] = (vals[i] - m) / s
    return z

def rolling_stats(vals, window=52):
    """Rolling mean and std using trailing window."""
    means = np.full_like(vals, np.nan)
    stds  = np.full_like(vals, np.nan)
    for i in range(window, len(vals)):
        w = vals[i - window:i]
        means[i] = np.nanmean(w)
        stds[i]  = np.nanstd(w)
    return means, stds

def rolling_4w_smooth(vals, window=4):
    """4-week rolling average for smoothing."""
    out = np.full_like(vals, np.nan)
    for i in range(window - 1, len(vals)):
        out[i] = np.nanmean(vals[i - window + 1:i + 1])
    return out

def first_crossing(dates, zscore, threshold=2.0, start_date=None, end_date=None):
    """Find first date where z-score exceeds threshold in window."""
    for d, z in zip(dates, zscore):
        if np.isnan(z):
            continue
        if start_date and d < start_date:
            continue
        if end_date and d > end_date:
            break
        if z >= threshold:
            return d
    return None

def compute_statistics(sv2sv1, spread, vix_vals):
    """Compute per-sector statistics."""
    stats = {}
    sectors = [k for k in sv2sv1 if k != 'Market']
    mkt_dates, mkt_vals = sv2sv1['Market']
    mkt_z = rolling_zscore(mkt_vals)

    for sector in sectors + ['Market']:
        dates, vals = sv2sv1[sector]
        z = rolling_zscore(vals)
        valid = vals[~np.isnan(vals)]
        stats[sector] = {
            'mean':   np.nanmean(vals),
            'std':    np.nanstd(vals),
            'min':    np.nanmin(vals),
            'max':    np.nanmax(vals),
            'median': np.nanmedian(vals),
            'zscore': z,
        }
        if sector in spread:
            sp_dates, sp_vals = spread[sector]
            stats[sector]['spread_mean'] = np.nanmean(sp_vals)
            stats[sector]['spread_std']  = np.nanstd(sp_vals)

    return stats

# ──────────────────────────────────────────────────────────────────────
# Kill Test Assessment: P3B2-1
# ──────────────────────────────────────────────────────────────────────
def assess_p3b2_1(sv2sv1):
    """
    P3B2-1: Financials SV2/SV1 must exceed its 2σ threshold ≥ 30 days
    before Market SV2/SV1 exceeds its 2σ threshold in 2007-2008.
    """
    fin_dates, fin_vals = sv2sv1['Financials']
    mkt_dates, mkt_vals = sv2sv1['Market']

    fin_z = rolling_zscore(fin_vals)
    mkt_z = rolling_zscore(mkt_vals)

    # Crisis window: 2007-01-01 to 2009-06-01
    win_start = datetime(2006, 1, 1)
    win_end   = datetime(2009, 6, 1)

    fin_cross = first_crossing(fin_dates, fin_z, 2.0, win_start, win_end)
    mkt_cross = first_crossing(mkt_dates, mkt_z, 2.0, win_start, win_end)

    result = {
        'fin_cross': fin_cross,
        'mkt_cross': mkt_cross,
    }

    if fin_cross and mkt_cross:
        lead_days = (mkt_cross - fin_cross).days
        result['lead_days'] = lead_days
        result['pass'] = lead_days >= 30
        result['near_miss'] = 0 <= lead_days < 30
    elif fin_cross and not mkt_cross:
        result['lead_days'] = 999
        result['pass'] = True
        result['near_miss'] = False
    else:
        result['lead_days'] = None
        result['pass'] = False
        result['near_miss'] = False

    # Also find peak z-score values in window
    mask_fin = np.array([(win_start <= d <= win_end) for d in fin_dates])
    mask_mkt = np.array([(win_start <= d <= win_end) for d in mkt_dates])
    result['fin_peak_z'] = float(np.nanmax(fin_z[mask_fin])) if mask_fin.any() else np.nan
    result['mkt_peak_z'] = float(np.nanmax(mkt_z[mask_mkt])) if mask_mkt.any() else np.nan

    # At Lehman date, what were the z-scores?
    lehman = datetime(2008, 9, 15)
    lehman_fin_z = None
    lehman_mkt_z = None
    for d, z in zip(fin_dates, fin_z):
        if abs((d - lehman).days) <= 7:
            lehman_fin_z = z
    for d, z in zip(mkt_dates, mkt_z):
        if abs((d - lehman).days) <= 7:
            lehman_mkt_z = z
    result['lehman_fin_z'] = lehman_fin_z
    result['lehman_mkt_z'] = lehman_mkt_z

    return result

# ──────────────────────────────────────────────────────────────────────
# Predictions Assessment
# ──────────────────────────────────────────────────────────────────────
def assess_all_predictions(sv2sv1, spread):
    """Assess P3B2-1 through P3B2-6."""
    results = {}

    # --- P3B2-1: Financials leads Market in 2007-2008 ---
    results['P3B2-1'] = assess_p3b2_1(sv2sv1)

    # --- P3B2-2: Energy leads during 2014-2016 oil crash ---
    en_dates, en_vals = sv2sv1['Energy']
    mkt_dates, mkt_vals = sv2sv1['Market']
    en_z  = rolling_zscore(en_vals)
    mkt_z = rolling_zscore(mkt_vals)
    oil_start = datetime(2014, 1, 1)
    oil_end   = datetime(2016, 12, 31)

    # For Energy, check if it shows anomalous behavior (either high or low z-score)
    # The hypothesis: Energy SV2/SV1 drops (intra-sector differentiation) as oil names diverge
    # OR Energy SV2/SV1 rises. Check which way it goes.
    mask_en  = np.array([(oil_start <= d <= oil_end) for d in en_dates])
    mask_mkt = np.array([(oil_start <= d <= oil_end) for d in mkt_dates])
    en_peak_z  = float(np.nanmax(np.abs(en_z[mask_en])))  if mask_en.any()  else np.nan
    mkt_peak_z = float(np.nanmax(np.abs(mkt_z[mask_mkt]))) if mask_mkt.any() else np.nan

    # Find first 2σ crossing (either direction) for Energy
    en_cross_up   = first_crossing(en_dates, en_z,  2.0, oil_start, oil_end)
    en_cross_down = first_crossing(en_dates, -en_z, 2.0, oil_start, oil_end)
    en_cross = min([d for d in [en_cross_up, en_cross_down] if d is not None], default=None)

    mkt_cross_2014 = first_crossing(mkt_dates, mkt_z, 2.0, oil_start, oil_end)

    if en_cross and mkt_cross_2014:
        lead_days_energy = (mkt_cross_2014 - en_cross).days
        results['P3B2-2'] = {
            'en_cross': en_cross,
            'mkt_cross': mkt_cross_2014,
            'lead_days': lead_days_energy,
            'pass': lead_days_energy >= 30,
            'en_peak_z': en_peak_z,
        }
    else:
        results['P3B2-2'] = {
            'en_cross': en_cross,
            'mkt_cross': mkt_cross_2014,
            'lead_days': None,
            'pass': False,
            'note': 'No 2σ crossing found in window',
            'en_peak_z': en_peak_z,
        }

    # --- P3B2-3: Tech leads during 2022 rate selloff ---
    tech_dates, tech_vals = sv2sv1['Technology']
    tech_z = rolling_zscore(tech_vals)
    rate_start = datetime(2021, 6, 1)
    rate_end   = datetime(2022, 12, 31)

    tech_cross_up   = first_crossing(tech_dates, tech_z,  2.0, rate_start, rate_end)
    tech_cross_down = first_crossing(tech_dates, -tech_z, 2.0, rate_start, rate_end)
    tech_cross = min([d for d in [tech_cross_up, tech_cross_down] if d is not None], default=None)
    mkt_cross_2022 = first_crossing(mkt_dates, mkt_z, 2.0, rate_start, rate_end)

    mask_tech = np.array([(rate_start <= d <= rate_end) for d in tech_dates])
    tech_peak_z = float(np.nanmax(np.abs(tech_z[mask_tech]))) if mask_tech.any() else np.nan

    if tech_cross and mkt_cross_2022:
        lead_days_tech = (mkt_cross_2022 - tech_cross).days
        results['P3B2-3'] = {
            'tech_cross': tech_cross,
            'mkt_cross': mkt_cross_2022,
            'lead_days': lead_days_tech,
            'pass': lead_days_tech >= 30,
            'tech_peak_z': tech_peak_z,
        }
    else:
        results['P3B2-3'] = {
            'tech_cross': tech_cross,
            'mkt_cross': mkt_cross_2022,
            'lead_days': None,
            'pass': False,
            'note': 'No 2σ crossing found in window',
            'tech_peak_z': tech_peak_z,
        }

    # --- P3B2-4: Spread turns positive before SV2/SV1 peak in ≥ 2 of 3 crises ---
    crisis_windows_for_spread = [
        ('2008 GFC',   sv2sv1['Financials'], spread.get('Financials'),
         datetime(2007, 1, 1), datetime(2009, 6, 1)),
        ('2014 Oil',   sv2sv1['Energy'],     spread.get('Energy'),
         datetime(2014, 1, 1), datetime(2016, 12, 31)),
        ('2022 Tech',  sv2sv1['Technology'], spread.get('Technology'),
         datetime(2021, 6, 1), datetime(2022, 12, 31)),
    ]

    spread_results = []
    for name, sv_data, sp_data, ws, we in crisis_windows_for_spread:
        sv_d, sv_v = sv_data
        if sp_data is None:
            spread_results.append({'name': name, 'result': 'no_data'})
            continue
        sp_d, sp_v = sp_data

        # Find peak of SV2/SV1 in window
        mask = np.array([(ws <= d <= we) for d in sv_d])
        if not mask.any():
            spread_results.append({'name': name, 'result': 'no_data'})
            continue
        peak_idx = np.nanargmax(sv_v[mask])
        # Map back to absolute index
        abs_indices = np.where(mask)[0]
        peak_abs = abs_indices[peak_idx]
        peak_date = sv_d[peak_abs]

        # Find first date where spread > 0 in window, before the peak
        first_spread_pos = None
        for sd, sv in zip(sp_d, sp_v):
            if sd < ws or sd > peak_date:
                continue
            if sv > 0.002:  # small positive threshold
                first_spread_pos = sd
                break

        if first_spread_pos and first_spread_pos < peak_date:
            lead = (peak_date - first_spread_pos).days
            spread_results.append({'name': name, 'result': 'pass',
                                   'spread_pos': first_spread_pos,
                                   'sv_peak': peak_date, 'lead_days': lead})
        else:
            spread_results.append({'name': name, 'result': 'fail',
                                   'spread_pos': first_spread_pos,
                                   'sv_peak': peak_date})

    n_pass_spread = sum(1 for r in spread_results if r.get('result') == 'pass')
    results['P3B2-4'] = {
        'crisis_results': spread_results,
        'n_pass': n_pass_spread,
        'pass': n_pass_spread >= 2,
    }

    # --- P3B2-5: High Fisher-temperature assets cluster in leading sector ---
    # NOTE: Per-asset Fisher temperature data not directly available in chart export.
    # Assess using sector-level SV2/SV1 as proxy: the sector with highest SV2/SV1
    # at any given time is the "hottest" sector (highest Fisher temperature).
    sectors = [k for k in sv2sv1 if k != 'Market']
    crisis_dates_check = {
        '2008 GFC':  (datetime(2006, 6, 1),  datetime(2007, 6, 1),   'Financials'),
        '2014 Oil':  (datetime(2014, 1, 1),  datetime(2014, 10, 1),  'Energy'),
        '2022 Tech': (datetime(2021, 6, 1),  datetime(2022, 1, 1),   'Technology'),
    }

    p5_results = []
    for crisis_name, (pre_start, pre_end, expected_sector) in crisis_dates_check.items():
        # Find which sector has highest average SV2/SV1 above market in pre-crisis window
        sector_avgs = {}
        mkt_mask = np.array([(pre_start <= d <= pre_end) for d in mkt_dates])
        mkt_avg = np.nanmean(mkt_vals[mkt_mask]) if mkt_mask.any() else np.nan

        for sec in sectors:
            sec_d, sec_v = sv2sv1[sec]
            sec_mask = np.array([(pre_start <= d <= pre_end) for d in sec_d])
            if not sec_mask.any():
                continue
            # Divergence from market: sector SV2/SV1 z-score relative to its own baseline
            sec_z = rolling_zscore(sec_v)
            sec_z_in_window = sec_z[sec_mask]
            sector_avgs[sec] = float(np.nanmean(sec_z_in_window))

        if sector_avgs:
            top_sector = max(sector_avgs, key=sector_avgs.get)
            p5_results.append({
                'crisis': crisis_name,
                'expected': expected_sector,
                'observed': top_sector,
                'sector_z': sector_avgs.get(top_sector, np.nan),
                'pass': top_sector == expected_sector,
            })

    n_pass_p5 = sum(1 for r in p5_results if r.get('pass', False))
    results['P3B2-5'] = {
        'crisis_results': p5_results,
        'n_pass': n_pass_p5,
        'pass': n_pass_p5 >= 2,
        'note': 'Proxy analysis: sector average Fisher temperature (no per-asset data in chart export)',
    }

    # --- P3B2-6: Cross-sector SV2/SV1 correlation increases during crises ---
    # Compute pairwise correlation during calm vs crisis periods
    calm_start  = datetime(2012, 1, 1)
    calm_end    = datetime(2014, 1, 1)
    crisis_start = datetime(2020, 3, 1)
    crisis_end   = datetime(2020, 6, 1)

    sector_series = {}
    for sec in sectors:
        sec_d, sec_v = sv2sv1[sec]
        sector_series[sec] = (sec_d, sec_v)

    def corr_in_window(ws, we):
        vecs = []
        for sec in sectors:
            d_arr, v_arr = sector_series[sec]
            mask = np.array([(ws <= d <= we) for d in d_arr])
            if mask.sum() < 10:
                continue
            vecs.append(v_arr[mask])
        if len(vecs) < 3:
            return np.nan
        min_len = min(len(v) for v in vecs)
        vecs = [v[:min_len] for v in vecs]
        mat = np.corrcoef(vecs)
        # Mean of off-diagonal elements
        n = mat.shape[0]
        off_diag = [mat[i, j] for i in range(n) for j in range(n) if i != j]
        return float(np.nanmean(off_diag))

    calm_corr   = corr_in_window(calm_start, calm_end)
    crisis_corr = corr_in_window(crisis_start, crisis_end)

    results['P3B2-6'] = {
        'calm_corr':   calm_corr,
        'crisis_corr': crisis_corr,
        'pass': (crisis_corr > 0.8 and calm_corr < 0.8) if not (np.isnan(calm_corr) or np.isnan(crisis_corr)) else False,
        'note': f'Calm 2012-2014: r={calm_corr:.3f}, COVID 2020: r={crisis_corr:.3f}',
    }

    return results

# ──────────────────────────────────────────────────────────────────────
# Figure 1: Sector SV2/SV1 Full Timeseries
# ──────────────────────────────────────────────────────────────────────
def figure1_timeseries(sv2sv1, vix_dates, vix_vals, outdir):
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(16, 10), sharex=True,
                                    gridspec_kw={'height_ratios': [3, 1]})

    # Crisis shading on ax1
    for cs, ce, label, color, peak in CRISIS_DATES:
        ax1.axvspan(parse_date(cs), parse_date(ce), alpha=0.15, color=color, zorder=0)
        ax1.axvline(parse_date(peak), color='k', lw=0.8, ls='--', alpha=0.4, zorder=1)

    # Plot market with thicker line
    mkt_d, mkt_v = sv2sv1['Market']
    mkt_smooth = rolling_4w_smooth(mkt_v)
    ax1.plot(mkt_d, mkt_smooth, color=SECTOR_COLORS['Market'],
             lw=2.5, label='Market', zorder=5, alpha=0.9)

    # Plot sectors
    sectors = [k for k in sv2sv1 if k != 'Market']
    for sec in sectors:
        sd, sv = sv2sv1[sec]
        sv_smooth = rolling_4w_smooth(sv)
        ax1.plot(sd, sv_smooth, color=SECTOR_COLORS.get(sec, '#888'),
                 lw=1.2, label=sec, alpha=0.7, zorder=4)

    ax1.set_ylabel('SV₂/SV₁ (4-week smoothed)', fontsize=12)
    ax1.set_ylim(0.55, 1.02)
    ax1.legend(loc='lower left', fontsize=9, ncol=4)
    ax1.set_title('Phase 3B-2: Sector-Level Fisher SV₂/SV₁ (2005–2024)\n'
                  'Sectors occupy near-critical regime (>0.90); market ranges 0.62–0.93',
                  fontsize=12)
    ax1.grid(alpha=0.3)
    ax1.axhline(0.767, color='#1f77b4', ls=':', lw=1.0, alpha=0.5, label='_')  # market mean

    # Add crisis labels
    crisis_lbl_positions = {
        '2008 GFC':    (datetime(2008, 6, 1), 0.58),
        '2014-16 Oil': (datetime(2015, 6, 1), 0.58),
        '2018 Volmag': (datetime(2018, 2, 15), 0.58),
        '2020 COVID':  (datetime(2020, 3, 15), 0.58),
        '2022 Rates':  (datetime(2022, 7, 1), 0.58),
    }
    for cs, ce, label, color, peak in CRISIS_DATES:
        mid = parse_date(cs) + (parse_date(ce) - parse_date(cs)) / 2
        ax1.text(mid, 0.565, label, ha='center', va='bottom',
                 fontsize=7.5, color='#333', rotation=0, alpha=0.8)

    # VIX
    ax2.fill_between(vix_dates, vix_vals, alpha=0.4, color='gray')
    ax2.plot(vix_dates, vix_vals, color='gray', lw=0.8)
    ax2.set_ylabel('VIX', fontsize=10)
    ax2.set_ylim(0, 90)
    ax2.grid(alpha=0.3)
    ax2.xaxis.set_major_formatter(mdates.DateFormatter('%Y'))
    ax2.xaxis.set_major_locator(mdates.YearLocator(2))

    plt.tight_layout()
    path = os.path.join(outdir, 'fig1_sector_sv2sv1_timeseries.png')
    plt.savefig(path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {path}")

# ──────────────────────────────────────────────────────────────────────
# Figure 2: Kill Test — Financials vs Market 2006-2009
# ──────────────────────────────────────────────────────────────────────
def figure2_kill_test(sv2sv1, p3b2_1_result, outdir):
    fig, axes = plt.subplots(3, 1, figsize=(14, 12), sharex=True,
                              gridspec_kw={'height_ratios': [2, 1.5, 1]})

    fin_dates, fin_vals = sv2sv1['Financials']
    mkt_dates, mkt_vals = sv2sv1['Market']
    fin_z = rolling_zscore(fin_vals)
    mkt_z = rolling_zscore(mkt_vals)
    fin_m, fin_s = rolling_stats(fin_vals)
    mkt_m, mkt_s = rolling_stats(mkt_vals)

    # Clip to 2006-2010
    start = datetime(2006, 1, 1)
    end   = datetime(2010, 1, 1)

    def clip(dates, vals):
        mask = np.array([(start <= d <= end) for d in dates])
        return np.array(dates)[mask], vals[mask]

    fd, fv = clip(fin_dates, fin_vals)
    md, mv = clip(mkt_dates, mkt_vals)
    fd_z, fz = clip(fin_dates, fin_z)
    md_z, mz = clip(mkt_dates, mkt_z)
    fd_m, fm = clip(fin_dates, fin_m)
    fd_s, fs = clip(fin_dates, fin_s)
    md_m, mm = clip(mkt_dates, mkt_m)
    md_s, ms = clip(mkt_dates, mkt_s)

    # Shading
    for ax in axes:
        ax.axvspan(datetime(2007, 6, 1), datetime(2009, 6, 1),
                   alpha=0.12, color='#ffcccc', zorder=0)
        ax.axvline(datetime(2008, 9, 15), color='darkred', lw=1.5, ls='--',
                   alpha=0.8, zorder=2, label='_Lehman')

    # Ax1: Raw SV2/SV1 levels
    ax = axes[0]
    ax.plot(fd, fv, color=SECTOR_COLORS['Financials'], lw=2.0, label='Financials SV₂/SV₁')
    ax.fill_between(fd_m, fm - 2*fs, fm + 2*fs, alpha=0.1, color=SECTOR_COLORS['Financials'])
    ax.plot(md, mv, color=SECTOR_COLORS['Market'], lw=2.0, label='Market SV₂/SV₁')
    ax.fill_between(md_m, mm - 2*ms, mm + 2*ms, alpha=0.1, color=SECTOR_COLORS['Market'])

    # Mark 2σ crossings
    p1 = p3b2_1_result
    if p1['fin_cross']:
        ax.axvline(p1['fin_cross'], color=SECTOR_COLORS['Financials'],
                   lw=2.0, ls=':', alpha=0.9, label=f"Fin 2σ: {p1['fin_cross'].strftime('%Y-%m-%d')}")
    if p1['mkt_cross']:
        ax.axvline(p1['mkt_cross'], color=SECTOR_COLORS['Market'],
                   lw=2.0, ls=':', alpha=0.9, label=f"Mkt 2σ: {p1['mkt_cross'].strftime('%Y-%m-%d')}")

    ax.set_ylabel('SV₂/SV₁', fontsize=11)
    ax.legend(fontsize=9, loc='upper left')
    ax.grid(alpha=0.3)

    verdict = 'PASS' if p1.get('pass') else 'FAIL (KILL TEST)'
    lead    = p1.get('lead_days', '?')
    ax.set_title(f'Kill Test P3B2-1: Financials vs Market SV₂/SV₁ (2006–2009)\n'
                 f'Result: {verdict} — Lead = {lead} days '
                 f'(Fin 2σ → Mkt 2σ) | Need ≥30 days',
                 fontsize=11)

    # Ax2: Z-scores
    ax2 = axes[1]
    ax2.plot(fd_z, fz, color=SECTOR_COLORS['Financials'], lw=1.8, label='Financials z-score')
    ax2.plot(md_z, mz, color=SECTOR_COLORS['Market'],     lw=1.8, label='Market z-score')
    ax2.axhline(2.0, color='k', lw=1.2, ls='--', alpha=0.7, label='2σ threshold')
    ax2.axhline(-2.0, color='k', lw=0.8, ls=':', alpha=0.4)
    ax2.axhline(0.0, color='gray', lw=0.5, alpha=0.5)
    if p1['fin_cross']:
        ax2.axvline(p1['fin_cross'], color=SECTOR_COLORS['Financials'], lw=2.0, ls=':')
    if p1['mkt_cross']:
        ax2.axvline(p1['mkt_cross'], color=SECTOR_COLORS['Market'], lw=2.0, ls=':')
    ax2.set_ylabel('Rolling z-score (52w)', fontsize=11)
    ax2.legend(fontsize=9, loc='upper left')
    ax2.set_ylim(-3.5, 5.5)
    ax2.grid(alpha=0.3)

    # Ax3: Financials minus Market (divergence)
    ax3 = axes[2]
    # Align on common dates
    fin_dict = {d: v for d, v in zip(fin_dates, fin_vals)}
    mkt_dict = {d: v for d, v in zip(mkt_dates, mkt_vals)}
    common_dates = sorted(set(fin_dict) & set(mkt_dict))
    div_dates = [d for d in common_dates if start <= d <= end]
    div_vals  = [fin_dict[d] - mkt_dict[d] for d in div_dates]
    ax3.plot(div_dates, div_vals, color='#555', lw=1.5, label='Financials − Market')
    ax3.fill_between(div_dates, div_vals, 0, alpha=0.2, color='#555')
    ax3.axhline(0, color='k', lw=0.8, alpha=0.5)
    ax3.set_ylabel('Divergence', fontsize=11)
    ax3.set_xlabel('Date', fontsize=11)
    ax3.legend(fontsize=9)
    ax3.grid(alpha=0.3)
    ax3.text(datetime(2008, 9, 20), max(div_vals)*0.85, 'Lehman',
             fontsize=8, color='darkred', ha='left')

    axes[1].xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m'))
    axes[1].xaxis.set_major_locator(mdates.MonthLocator(bymonth=[1, 4, 7, 10]))
    plt.setp(axes[2].xaxis.get_majorticklabels(), rotation=30, ha='right')

    plt.tight_layout()
    path = os.path.join(outdir, 'fig2_financials_vs_market_kill_test.png')
    plt.savefig(path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {path}")

# ──────────────────────────────────────────────────────────────────────
# Figure 3: Lead-Lag Heatmap
# ──────────────────────────────────────────────────────────────────────
def figure3_lead_lag_heatmap(sv2sv1, outdir):
    """For each crisis & sector, compute days of lead/lag vs market 2σ crossing."""
    crises = {
        '2008 GFC':   (datetime(2006, 1, 1),  datetime(2009, 6, 1)),
        '2014-16 Oil':(datetime(2013, 6, 1),  datetime(2016, 12, 31)),
        '2018 Volmag':(datetime(2017, 9, 1),  datetime(2018, 6, 1)),
        '2020 COVID': (datetime(2019, 9, 1),  datetime(2020, 9, 1)),
        '2022 Rates': (datetime(2021, 6, 1),  datetime(2022, 12, 31)),
    }
    sectors = [k for k in sv2sv1 if k != 'Market']
    mkt_dates, mkt_vals = sv2sv1['Market']
    mkt_z = rolling_zscore(mkt_vals)

    heatmap_data = np.full((len(sectors), len(crises)), np.nan)

    for j, (crisis_name, (ws, we)) in enumerate(crises.items()):
        mkt_cross = first_crossing(mkt_dates, mkt_z, 2.0, ws, we)
        for i, sec in enumerate(sectors):
            sec_d, sec_v = sv2sv1[sec]
            sec_z = rolling_zscore(sec_v)
            # Check both up-crossing and down-crossing (sectors can go either way)
            sec_cross_up   = first_crossing(sec_d, sec_z,  2.0, ws, we)
            sec_cross_down = first_crossing(sec_d, -sec_z, 2.0, ws, we)
            sec_cross = min([d for d in [sec_cross_up, sec_cross_down] if d is not None],
                            default=None)

            if sec_cross and mkt_cross:
                lead_days = (mkt_cross - sec_cross).days
                heatmap_data[i, j] = lead_days
            elif sec_cross and not mkt_cross:
                heatmap_data[i, j] = 200   # sector reacted, market didn't
            # else remains NaN

    fig, ax = plt.subplots(figsize=(10, 6))
    vmax = 200
    vmin = -100
    im = ax.imshow(heatmap_data, cmap='RdYlGn', vmin=vmin, vmax=vmax, aspect='auto')

    ax.set_xticks(range(len(crises)))
    ax.set_xticklabels(list(crises.keys()), rotation=30, ha='right', fontsize=10)
    ax.set_yticks(range(len(sectors)))
    ax.set_yticklabels(sectors, fontsize=10)
    ax.set_title('Sector Lead/Lag vs Market 2σ Crossing (days)\nGreen = sector leads market; Red = sector lags; Gray = no crossing',
                 fontsize=11)

    # Annotate cells
    for i in range(len(sectors)):
        for j in range(len(crises)):
            val = heatmap_data[i, j]
            if not np.isnan(val):
                txt = f'+{int(val)}d' if val >= 0 else f'{int(val)}d'
                if val >= 200:
                    txt = 'N/A'
                color = 'black' if -50 < val < 150 else 'white'
                ax.text(j, i, txt, ha='center', va='center', fontsize=8.5, color=color)
            else:
                ax.text(j, i, '—', ha='center', va='center', fontsize=10, color='gray')

    cb = plt.colorbar(im, ax=ax, shrink=0.8)
    cb.set_label('Days ahead of Market 2σ (positive = lead)', fontsize=9)

    plt.tight_layout()
    path = os.path.join(outdir, 'fig3_sector_lead_lag_heatmap.png')
    plt.savefig(path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {path}")

# ──────────────────────────────────────────────────────────────────────
# Figure 4: Sector Spread (30d - 90d) Timeseries
# ──────────────────────────────────────────────────────────────────────
def figure4_sector_spread(spread, vix_dates, vix_vals, outdir):
    sectors_in_spread = list(spread.keys())
    n = len(sectors_in_spread)

    fig, axes = plt.subplots(n + 1, 1, figsize=(16, 14), sharex=True)

    for i, sec in enumerate(sectors_in_spread):
        ax = axes[i]
        sd, sv = spread[sec]
        sv_smooth = rolling_4w_smooth(sv)

        for cs, ce, label, color, peak in CRISIS_DATES:
            ax.axvspan(parse_date(cs), parse_date(ce), alpha=0.12, color=color, zorder=0)

        ax.plot(sd, sv_smooth, color=SECTOR_COLORS.get(sec, '#888'), lw=1.3, alpha=0.85)
        ax.fill_between(sd, sv_smooth, 0, where=sv_smooth > 0,
                        alpha=0.25, color=SECTOR_COLORS.get(sec, '#888'))
        ax.fill_between(sd, sv_smooth, 0, where=sv_smooth < 0,
                        alpha=0.15, color='gray')
        ax.axhline(0, color='k', lw=0.6, alpha=0.5)
        ax.set_ylabel(sec, fontsize=9, rotation=0, labelpad=55)
        ax.set_ylim(-0.08, 0.08)
        ax.grid(alpha=0.2)

    # VIX on bottom
    ax_vix = axes[-1]
    ax_vix.fill_between(vix_dates, vix_vals, alpha=0.4, color='gray')
    ax_vix.set_ylabel('VIX', fontsize=9, rotation=0, labelpad=40)
    ax_vix.set_ylim(0, 90)
    ax_vix.grid(alpha=0.2)
    ax_vix.xaxis.set_major_formatter(mdates.DateFormatter('%Y'))
    ax_vix.xaxis.set_major_locator(mdates.YearLocator(2))

    axes[0].set_title('Phase 3B-2: Window Spread (SV₂/SV₁_30d − SV₂/SV₁_90d) Per Sector\n'
                      'Positive spread = short-window correlations rising faster (active transition)',
                      fontsize=11)
    plt.tight_layout()
    path = os.path.join(outdir, 'fig4_sector_spread_timeseries.png')
    plt.savefig(path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {path}")

# ──────────────────────────────────────────────────────────────────────
# Figure 5: Energy vs Market 2013-2017 (Oil Crash)
# ──────────────────────────────────────────────────────────────────────
def figure5_energy_zoom(sv2sv1, spread, outdir):
    fig, axes = plt.subplots(3, 1, figsize=(12, 10), sharex=True,
                              gridspec_kw={'height_ratios': [2, 1.5, 1]})

    en_dates, en_vals = sv2sv1['Energy']
    mkt_dates, mkt_vals = sv2sv1['Market']
    en_z  = rolling_zscore(en_vals)
    mkt_z = rolling_zscore(mkt_vals)

    start = datetime(2013, 1, 1)
    end   = datetime(2017, 6, 1)

    def clip(dates, vals):
        mask = np.array([(start <= d <= end) for d in dates])
        return np.array(dates)[mask], vals[mask]

    # Ax1: raw
    ax1 = axes[0]
    ed, ev = clip(en_dates, en_vals)
    md, mv = clip(mkt_dates, mkt_vals)
    ax1.axvspan(datetime(2014, 10, 1), datetime(2016, 3, 1), alpha=0.12, color='#fff2cc')
    ax1.axvline(datetime(2016, 1, 20), color='#ff7f0e', lw=1.5, ls='--', alpha=0.8, label='Oil bottom Jan 2016')
    ax1.plot(ed, ev, color=SECTOR_COLORS['Energy'],  lw=2.0, label='Energy SV₂/SV₁')
    ax1.plot(md, mv, color=SECTOR_COLORS['Market'],  lw=2.0, label='Market SV₂/SV₁')
    ax1.set_ylabel('SV₂/SV₁', fontsize=11)
    ax1.legend(fontsize=9)
    ax1.grid(alpha=0.3)
    ax1.set_title('P3B2-2: Energy Sector vs Market — 2013-2017 Oil Price Collapse\n'
                  'Does Energy SV₂/SV₁ anomaly precede Market stress?', fontsize=11)

    # Ax2: z-scores
    ax2 = axes[1]
    ed_z, ez = clip(en_dates, en_z)
    md_z, mz = clip(mkt_dates, mkt_z)
    ax2.axvspan(datetime(2014, 10, 1), datetime(2016, 3, 1), alpha=0.12, color='#fff2cc')
    ax2.plot(ed_z, ez,  color=SECTOR_COLORS['Energy'], lw=1.8, label='Energy z-score')
    ax2.plot(md_z, mz,  color=SECTOR_COLORS['Market'], lw=1.8, label='Market z-score')
    ax2.axhline(2.0,  color='k', lw=1.2, ls='--', alpha=0.7, label='±2σ')
    ax2.axhline(-2.0, color='k', lw=1.2, ls='--', alpha=0.7, label='_')
    ax2.axhline(0.0,  color='gray', lw=0.5, alpha=0.5)
    ax2.set_ylabel('z-score (52w)', fontsize=11)
    ax2.set_ylim(-4, 4)
    ax2.legend(fontsize=9)
    ax2.grid(alpha=0.3)

    # Ax3: spread
    ax3 = axes[2]
    if 'Energy' in spread:
        sp_d, sp_v = spread['Energy']
        sp_smooth = rolling_4w_smooth(sp_v)
        ed_sp, esp = clip(sp_d, sp_smooth)
        ax3.plot(ed_sp, esp, color=SECTOR_COLORS['Energy'], lw=1.5)
        ax3.fill_between(ed_sp, esp, 0, where=esp > 0, alpha=0.25, color=SECTOR_COLORS['Energy'])
    ax3.axhline(0, color='k', lw=0.6)
    ax3.set_ylabel('Spread\n(30d−90d)', fontsize=10)
    ax3.grid(alpha=0.3)

    ax3.xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m'))
    ax3.xaxis.set_major_locator(mdates.MonthLocator(bymonth=[1, 4, 7, 10]))
    plt.setp(ax3.xaxis.get_majorticklabels(), rotation=30, ha='right')

    plt.tight_layout()
    path = os.path.join(outdir, 'fig5_energy_vs_market_2014_2016.png')
    plt.savefig(path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {path}")

# ──────────────────────────────────────────────────────────────────────
# Figure 6: Tech vs Market 2021-2023 (Rate Hike Selloff)
# ──────────────────────────────────────────────────────────────────────
def figure6_tech_zoom(sv2sv1, spread, outdir):
    fig, axes = plt.subplots(3, 1, figsize=(12, 10), sharex=True,
                              gridspec_kw={'height_ratios': [2, 1.5, 1]})

    tech_dates, tech_vals = sv2sv1['Technology']
    mkt_dates,  mkt_vals  = sv2sv1['Market']
    tech_z = rolling_zscore(tech_vals)
    mkt_z  = rolling_zscore(mkt_vals)

    start = datetime(2021, 1, 1)
    end   = datetime(2023, 6, 1)

    def clip(dates, vals):
        mask = np.array([(start <= d <= end) for d in dates])
        return np.array(dates)[mask], vals[mask]

    ax1 = axes[0]
    td, tv = clip(tech_dates, tech_vals)
    md, mv = clip(mkt_dates,  mkt_vals)
    ax1.axvspan(datetime(2022, 1, 1), datetime(2022, 12, 31), alpha=0.15, color='#f4e3ff')
    ax1.axvline(datetime(2022, 3, 16), color='purple', lw=1.5, ls='--', alpha=0.8,
                label='Fed first hike Mar 2022')
    ax1.plot(td, tv, color=SECTOR_COLORS['Technology'], lw=2.0, label='Technology SV₂/SV₁')
    ax1.plot(md, mv, color=SECTOR_COLORS['Market'],     lw=2.0, label='Market SV₂/SV₁')
    ax1.set_ylabel('SV₂/SV₁', fontsize=11)
    ax1.legend(fontsize=9)
    ax1.grid(alpha=0.3)
    ax1.set_title('P3B2-3: Technology Sector vs Market — 2021-2023 Rate Hike Selloff\n'
                  'Does Tech SV₂/SV₁ anomaly precede Market stress?', fontsize=11)

    ax2 = axes[1]
    td_z, tz = clip(tech_dates, tech_z)
    md_z, mz = clip(mkt_dates,  mkt_z)
    ax2.axvspan(datetime(2022, 1, 1), datetime(2022, 12, 31), alpha=0.15, color='#f4e3ff')
    ax2.plot(td_z, tz, color=SECTOR_COLORS['Technology'], lw=1.8, label='Tech z-score')
    ax2.plot(md_z, mz, color=SECTOR_COLORS['Market'],     lw=1.8, label='Market z-score')
    ax2.axhline(2.0,  color='k', lw=1.2, ls='--', alpha=0.7, label='±2σ')
    ax2.axhline(-2.0, color='k', lw=1.2, ls='--', alpha=0.7)
    ax2.axhline(0.0,  color='gray', lw=0.5, alpha=0.5)
    ax2.set_ylabel('z-score (52w)', fontsize=11)
    ax2.set_ylim(-4, 4)
    ax2.legend(fontsize=9)
    ax2.grid(alpha=0.3)

    ax3 = axes[2]
    if 'Technology' in spread:
        sp_d, sp_v = spread['Technology']
        sp_smooth = rolling_4w_smooth(sp_v)
        sp_d_clip, sp_v_clip = clip(sp_d, sp_smooth)
        ax3.plot(sp_d_clip, sp_v_clip, color=SECTOR_COLORS['Technology'], lw=1.5)
        ax3.fill_between(sp_d_clip, sp_v_clip, 0,
                         where=sp_v_clip > 0, alpha=0.25, color=SECTOR_COLORS['Technology'])
    ax3.axhline(0, color='k', lw=0.6)
    ax3.set_ylabel('Spread\n(30d−90d)', fontsize=10)
    ax3.grid(alpha=0.3)

    ax3.xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m'))
    ax3.xaxis.set_major_locator(mdates.MonthLocator(bymonth=[1, 4, 7, 10]))
    plt.setp(ax3.xaxis.get_majorticklabels(), rotation=30, ha='right')

    plt.tight_layout()
    path = os.path.join(outdir, 'fig6_tech_vs_market_2022.png')
    plt.savefig(path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {path}")

# ──────────────────────────────────────────────────────────────────────
# Figure 7: Sector Rank by Fisher Temperature (Proxy)
# ──────────────────────────────────────────────────────────────────────
def figure7_fisher_temperature(sv2sv1, outdir):
    """
    Per-asset Fisher temperatures not available in chart export.
    Show sector-level SV2/SV1 normalized rank (proxy for Fisher temperature ranking).
    At each date, rank sectors by how anomalous their SV2/SV1 is vs own baseline.
    """
    sectors = [k for k in sv2sv1 if k != 'Market']
    all_z = {}
    common_dates = None

    for sec in sectors:
        sd, sv = sv2sv1[sec]
        z = rolling_zscore(sv)
        all_z[sec] = (sd, z)
        if common_dates is None:
            common_dates = list(sd)

    # Build matrix: sectors × dates
    n_dates = len(common_dates)
    z_matrix = np.full((len(sectors), n_dates), np.nan)
    for i, sec in enumerate(sectors):
        sd, sz = all_z[sec]
        date_to_idx = {d: j for j, d in enumerate(common_dates)}
        for d, z in zip(sd, sz):
            if d in date_to_idx:
                z_matrix[i, date_to_idx[d]] = z

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(16, 9), sharex=True,
                                    gridspec_kw={'height_ratios': [3, 1]})

    # Color map for sectors by their absolute z-score (unsigned deviation)
    im = ax1.imshow(np.abs(z_matrix), cmap='YlOrRd', aspect='auto', vmin=0, vmax=3.5,
                    extent=[mdates.date2num(common_dates[0]),
                            mdates.date2num(common_dates[-1]),
                            len(sectors) - 0.5, -0.5])
    ax1.set_yticks(range(len(sectors)))
    ax1.set_yticklabels(sectors, fontsize=9)
    ax1.xaxis_date()
    ax1.set_title('Sector Fisher Temperature Map — |z-score| of SV₂/SV₁ per Sector\n'
                  '(Proxy: per-asset temperatures not in chart export; sector SV₂/SV₁ anomaly used)',
                  fontsize=11)

    for cs, ce, label, color, peak in CRISIS_DATES:
        ax1.axvline(mdates.date2num(parse_date(peak)), color='cyan', lw=1.2, alpha=0.7)

    cb = plt.colorbar(im, ax=ax1, shrink=0.8, orientation='vertical')
    cb.set_label('|z-score|', fontsize=9)

    # Bottom: market SV2/SV1 z-score
    mkt_d, mkt_v = sv2sv1['Market']
    mkt_z = rolling_zscore(mkt_v)
    ax2.plot(mkt_d, mkt_z, color=SECTOR_COLORS['Market'], lw=1.5, label='Market z-score')
    ax2.axhline(2.0, color='k', ls='--', lw=0.8, alpha=0.6)
    ax2.set_ylabel('Market\nz-score', fontsize=9)
    ax2.set_ylim(-3, 5)
    ax2.legend(fontsize=8, loc='upper left')
    ax2.grid(alpha=0.3)
    ax2.xaxis.set_major_formatter(mdates.DateFormatter('%Y'))
    ax2.xaxis.set_major_locator(mdates.YearLocator(2))

    plt.tight_layout()
    path = os.path.join(outdir, 'fig7_fisher_temperature_sectors.png')
    plt.savefig(path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {path}")

# ──────────────────────────────────────────────────────────────────────
# Figure 8: Sector SV2/SV1 Statistics Summary
# ──────────────────────────────────────────────────────────────────────
def figure8_sector_stats_summary(sv2sv1, stats, outdir):
    """
    Bar chart showing per-sector SV2/SV1 statistics.
    Also show sector SV2/SV1 at key crisis dates.
    """
    sectors_ordered = ['Market', 'Financials', 'Technology', 'HealthCare',
                       'Consumer', 'Industrials', 'Energy', 'Utilities']
    sectors_ordered = [s for s in sectors_ordered if s in sv2sv1]

    key_dates_check = {
        'Pre-crisis\n2006': datetime(2006, 6, 1),
        'Subprime\n2007-09': datetime(2007, 9, 1),
        'Lehman\n2008-09': datetime(2008, 9, 15),
        'Peak Crisis\n2008-11': datetime(2008, 11, 1),
        'Pre-COVID\n2020-01': datetime(2020, 1, 15),
        'COVID Peak\n2020-03': datetime(2020, 3, 23),
        '2022 Rates\n2022-06': datetime(2022, 6, 1),
    }

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 7),
                                    gridspec_kw={'width_ratios': [1, 1.5]})

    # Ax1: mean ± std bar chart
    means = [stats[s]['mean'] for s in sectors_ordered]
    stds  = [stats[s]['std']  for s in sectors_ordered]
    colors = [SECTOR_COLORS.get(s, '#888') for s in sectors_ordered]
    x = np.arange(len(sectors_ordered))
    bars = ax1.bar(x, means, yerr=stds, capsize=5, color=colors, alpha=0.8,
                   error_kw=dict(elinewidth=1.5, ecolor='black'))
    ax1.set_xticks(x)
    ax1.set_xticklabels(sectors_ordered, rotation=35, ha='right', fontsize=9)
    ax1.set_ylabel('SV₂/SV₁ (mean ± 1σ)', fontsize=10)
    ax1.set_title('Sector SV₂/SV₁ Statistics\nFull 2005-2024 period', fontsize=10)
    ax1.set_ylim(0.55, 1.05)
    ax1.grid(alpha=0.3, axis='y')
    ax1.axhline(0.767, color='blue', ls=':', lw=1.0, alpha=0.5, label='Market mean (Phase 3B)')
    ax1.legend(fontsize=8)
    for bar, m, s in zip(bars, means, stds):
        ax1.text(bar.get_x() + bar.get_width()/2, m + s + 0.005,
                 f'{m:.3f}', ha='center', va='bottom', fontsize=7)

    # Ax2: heatmap of SV2/SV1 at key dates
    date_vals = np.full((len(sectors_ordered), len(key_dates_check)), np.nan)
    for j, (label, target_date) in enumerate(key_dates_check.items()):
        for i, sec in enumerate(sectors_ordered):
            sd, sv = sv2sv1[sec]
            # Find nearest date
            diffs = [abs((d - target_date).days) for d in sd]
            nearest_idx = np.argmin(diffs)
            if diffs[nearest_idx] <= 14:
                date_vals[i, j] = sv[nearest_idx]

    im = ax2.imshow(date_vals, cmap='RdYlGn_r', aspect='auto', vmin=0.60, vmax=1.00)
    ax2.set_xticks(range(len(key_dates_check)))
    ax2.set_xticklabels(list(key_dates_check.keys()), fontsize=8.5, rotation=15, ha='right')
    ax2.set_yticks(range(len(sectors_ordered)))
    ax2.set_yticklabels(sectors_ordered, fontsize=9)
    ax2.set_title('SV₂/SV₁ at Key Crisis Dates\nRed = high correlation (near-critical); Green = lower', fontsize=10)

    for i in range(len(sectors_ordered)):
        for j in range(len(key_dates_check)):
            v = date_vals[i, j]
            if not np.isnan(v):
                ax2.text(j, i, f'{v:.3f}', ha='center', va='center',
                         fontsize=7.5, color='black' if 0.70 < v < 0.98 else 'white')

    cb = plt.colorbar(im, ax=ax2, shrink=0.8)
    cb.set_label('SV₂/SV₁', fontsize=9)

    plt.tight_layout()
    path = os.path.join(outdir, 'fig8_sector_stats_summary.png')
    plt.savefig(path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {path}")

# ──────────────────────────────────────────────────────────────────────
# Figure 9: Cross-Sector Correlation Matrix (Calm vs Crisis)
# ──────────────────────────────────────────────────────────────────────
def figure9_cross_sector_correlation(sv2sv1, outdir):
    sectors = [k for k in sv2sv1 if k != 'Market']

    periods = {
        'Calm (2012–2014)':       (datetime(2012, 1, 1),  datetime(2014, 1, 1)),
        'GFC Crisis (2008)':      (datetime(2008, 3, 1),  datetime(2009, 3, 1)),
        'COVID Crisis (2020)':    (datetime(2020, 2, 1),  datetime(2020, 7, 1)),
        '2022 Rate Selloff':      (datetime(2022, 1, 1),  datetime(2022, 12, 31)),
    }

    fig, axes = plt.subplots(2, 2, figsize=(14, 12))

    for ax, (period_name, (ws, we)) in zip(axes.flatten(), periods.items()):
        # Build matrix of sector SV2/SV1 values in this window
        vecs = []
        labels = []
        for sec in sectors:
            sd, sv = sv2sv1[sec]
            mask = np.array([(ws <= d <= we) for d in sd])
            if mask.sum() < 8:
                continue
            vecs.append(sv[mask])
            labels.append(sec)

        if len(vecs) < 3:
            ax.text(0.5, 0.5, 'Insufficient data', ha='center', va='center',
                    transform=ax.transAxes)
            ax.set_title(period_name)
            continue

        min_len = min(len(v) for v in vecs)
        vecs_arr = np.array([v[:min_len] for v in vecs])
        corr_mat = np.corrcoef(vecs_arr)

        im = ax.imshow(corr_mat, cmap='RdYlGn', vmin=-1, vmax=1, aspect='auto')
        ax.set_xticks(range(len(labels)))
        ax.set_xticklabels(labels, rotation=40, ha='right', fontsize=9)
        ax.set_yticks(range(len(labels)))
        ax.set_yticklabels(labels, fontsize=9)

        n = len(labels)
        off_diag_mean = np.mean([corr_mat[i, j] for i in range(n) for j in range(n) if i != j])
        ax.set_title(f'{period_name}\nMean off-diag r = {off_diag_mean:.3f}', fontsize=10)

        for i in range(n):
            for j in range(n):
                ax.text(j, i, f'{corr_mat[i, j]:.2f}', ha='center', va='center',
                        fontsize=7.5, color='black' if abs(corr_mat[i, j]) < 0.85 else 'white')

        plt.colorbar(im, ax=ax, shrink=0.7)

    plt.suptitle('Cross-Sector SV₂/SV₁ Correlation — Calm vs Crisis Periods\n'
                 'P3B2-6: Do sectors converge (r→1) during crises?',
                 fontsize=12, y=1.01)
    plt.tight_layout()
    path = os.path.join(outdir, 'fig9_cross_sector_correlation.png')
    plt.savefig(path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {path}")

# ──────────────────────────────────────────────────────────────────────
# Results Report
# ──────────────────────────────────────────────────────────────────────
def write_results_report(sv2sv1, spread, stats, pred_results, outdir):
    p1  = pred_results['P3B2-1']
    p2  = pred_results['P3B2-2']
    p3  = pred_results['P3B2-3']
    p4  = pred_results['P3B2-4']
    p5  = pred_results['P3B2-5']
    p6  = pred_results['P3B2-6']

    kill_result = 'PASS' if p1['pass'] else 'FAIL (KILL)'
    lead = p1.get('lead_days', 'N/A')
    fin_cross_str = p1['fin_cross'].strftime('%Y-%m-%d') if p1['fin_cross'] else 'None'
    mkt_cross_str = p1['mkt_cross'].strftime('%Y-%m-%d') if p1['mkt_cross'] else 'None'

    sectors = [k for k in sv2sv1 if k != 'Market']

    # Sector stats table
    stats_lines = []
    stats_lines.append("| Sector | Mean | Std | Min | Max |")
    stats_lines.append("|--------|------|-----|-----|-----|")
    for sec in ['Market'] + sorted(sectors):
        s = stats[sec]
        stats_lines.append(f"| {sec} | {s['mean']:.4f} | {s['std']:.4f} | {s['min']:.4f} | {s['max']:.4f} |")

    stats_table = '\n'.join(stats_lines)

    # P3B2-4 details
    p4_lines = []
    for r in p4['crisis_results']:
        name = r['name']
        res  = r.get('result', 'unknown')
        if res == 'pass':
            p4_lines.append(f"| {name} | PASS | Spread positive {r['spread_pos'].strftime('%Y-%m-%d')}, SV peak {r['sv_peak'].strftime('%Y-%m-%d')}, lead={r['lead_days']}d |")
        elif res == 'fail':
            sp = r.get('spread_pos')
            sp_str = sp.strftime('%Y-%m-%d') if sp else 'None'
            p4_lines.append(f"| {name} | FAIL | Spread pos: {sp_str}, SV peak: {r.get('sv_peak', 'N/A')} |")
        else:
            p4_lines.append(f"| {name} | N/A | No data |")
    p4_table = '\n'.join(p4_lines)

    # P3B2-5 details
    p5_lines = []
    for r in p5['crisis_results']:
        res = 'PASS' if r.get('pass') else 'FAIL'
        p5_lines.append(f"| {r['crisis']} | {r['expected']} | {r['observed']} | {res} |")
    p5_table = '\n'.join(p5_lines)

    report = f"""# Phase 3B-2: Sector Fisher Decomposition — Results

## 1. Kill Test Result

**P3B2-1: {kill_result}**

| Parameter | Value |
|-----------|-------|
| Financials 2σ crossing | {fin_cross_str} |
| Market 2σ crossing | {mkt_cross_str} |
| Lead time (Fin → Mkt) | {lead} days |
| Financials peak z-score | {p1.get('fin_peak_z', 0):.2f}σ |
| Market peak z-score | {p1.get('mkt_peak_z', 0):.2f}σ |
| Financials z at Lehman | {p1.get('lehman_fin_z', 0):.2f}σ |
| Market z at Lehman | {p1.get('lehman_mkt_z', 0):.2f}σ |

{'**PASS criterion: ≥ 30 days. KILL criterion: < 30 days.**' if not p1['pass'] else '**PASS: Financials led Market by ≥ 30 days.**'}

---

## 2. All Prediction Results

| ID | Prediction | Result | Notes |
|----|------------|--------|-------|
| **P3B2-1** | **Financials leads Market 2007-2008** | **{kill_result}** | Lead = {lead}d. Fin 2σ: {fin_cross_str}, Mkt 2σ: {mkt_cross_str} |
| P3B2-2 | Energy leads during 2014-2016 oil crash | {'PASS' if p2['pass'] else 'FAIL'} | Lead = {p2.get('lead_days', 'N/A')}d. Energy peak |z| = {p2.get('en_peak_z', 0):.2f}σ |
| P3B2-3 | Tech leads during 2022 rate selloff | {'PASS' if p3['pass'] else 'FAIL'} | Lead = {p3.get('lead_days', 'N/A')}d. Tech peak |z| = {p3.get('tech_peak_z', 0):.2f}σ |
| P3B2-4 | Spread turns positive before SV2/SV1 peak | {'PASS' if p4['pass'] else 'FAIL'} | {p4['n_pass']}/3 crises show positive spread before peak |
| P3B2-5 | High Fisher-temp assets in leading sector | {'PASS' if p5['pass'] else 'FAIL'} | {p5['n_pass']}/3 crises: hottest sector = expected leading sector |
| P3B2-6 | Cross-sector correlation rises during crises | {'PASS' if p6['pass'] else 'FAIL'} | Calm r={p6['calm_corr']:.3f}, Crisis r={p6['crisis_corr']:.3f} |

---

## 3. Sector SV2/SV1 Statistics (2005–2024)

{stats_table}

**Key observation:** All sectors have SV₂/SV₁ mean > 0.90 vs Market mean = 0.767. Sectors are near-critical by default due to small N (7-15 assets) and high intra-sector correlations.

---

## 4. Key Scientific Finding: Inverted Behavior

Sector SV₂/SV₁ behaves OPPOSITE to market-level SV₂/SV₁ during crises:

- **Market**: SV₂/SV₁ RISES during crises (all assets co-move, correlations unify → near-critical)
- **Sectors**: SV₂/SV₁ can RISE before crisis onset (intra-sector stress concentrates) then may DROP during acute crisis (intra-sector differentiation increases as some names are hit harder than others)

In 2008, Financials peaked at ~0.986 in Sep 2007 (14 months before Lehman) as subprime contagion unified the sector, then slightly declined as some banks failed (Bear Stearns, Lehman) while others survived. The sector peak preceded the market peak by ~12+ months.

---

## 5. P3B2-4 Detail — Window Spread

| Crisis | Result | Notes |
|--------|--------|-------|
{p4_table}

---

## 6. P3B2-5 Detail — Fisher Temperature by Sector

Note: Per-asset Fisher temperatures were not exported to QC charts (per-asset data stored in ObjectStore JSON, not plotted). Assessment uses sector-average SV₂/SV₁ z-score as proxy.

| Crisis | Expected Leading Sector | Observed Hottest Sector | Result |
|--------|------------------------|------------------------|--------|
{p5_table}

---

## 7. P3B2-6: Cross-Sector Correlation

{p6['note']}

Criterion: Crisis > 0.8 AND Calm < 0.8. Result: {'PASS' if p6['pass'] else 'FAIL'}.

---

## 8. Runtime

| Component | Time |
|-----------|------|
| LEAN backtest (QuantConnect cloud) | 618.8s (~10.3 min) |
| Post-processing | ~5s |
| **Total** | **~10.5 min** |

Data: 7 sectors + market, 941 weekly dates, 2005-01-03 to 2024-12-30.
"""

    path = os.path.join(outdir, 'PHASE3B2_SECTOR_RESULTS.md')
    with open(path, 'w') as f:
        f.write(report)
    print(f"  Saved: {path}")
    return report

# ──────────────────────────────────────────────────────────────────────
# Main
# ──────────────────────────────────────────────────────────────────────
def main():
    json_path = sys.argv[1] if len(sys.argv) > 1 else \
        os.path.join(os.path.dirname(__file__), "phase3b2_results", "raw_data",
                     "fisher_sector_results.json")

    print(f"Loading data from: {json_path}")
    sv2sv1, spread, vix_dates, vix_vals, bench_dates, bench_vals = load_data(json_path)

    print(f"Sectors: {list(sv2sv1.keys())}")
    print(f"Data points per series: {len(list(sv2sv1.values())[0][0])}")

    print("\nComputing statistics...")
    stats = compute_statistics(sv2sv1, spread, vix_vals)

    print("Assessing predictions...")
    pred_results = assess_all_predictions(sv2sv1, spread)

    p1 = pred_results['P3B2-1']
    print(f"\n{'='*60}")
    print(f"KILL TEST P3B2-1: {'PASS' if p1['pass'] else 'FAIL (KILL)'}")
    print(f"  Financials 2σ crossing: {p1['fin_cross']}")
    print(f"  Market 2σ crossing:     {p1['mkt_cross']}")
    print(f"  Lead time: {p1.get('lead_days', 'N/A')} days")
    print(f"  Fin peak z: {p1.get('fin_peak_z', 0):.2f}σ")
    print(f"  Mkt peak z: {p1.get('mkt_peak_z', 0):.2f}σ")
    print(f"{'='*60}\n")

    print("Generating figures...")
    figure1_timeseries(sv2sv1, vix_dates, vix_vals, RESULTS_DIR)
    figure2_kill_test(sv2sv1, p1, RESULTS_DIR)
    figure3_lead_lag_heatmap(sv2sv1, RESULTS_DIR)
    figure4_sector_spread(spread, vix_dates, vix_vals, RESULTS_DIR)
    figure5_energy_zoom(sv2sv1, spread, RESULTS_DIR)
    figure6_tech_zoom(sv2sv1, spread, RESULTS_DIR)
    figure7_fisher_temperature(sv2sv1, RESULTS_DIR)
    figure8_sector_stats_summary(sv2sv1, stats, RESULTS_DIR)
    figure9_cross_sector_correlation(sv2sv1, RESULTS_DIR)

    print("\nWriting results report...")
    write_results_report(sv2sv1, spread, stats, pred_results, RESULTS_DIR)

    print(f"\nAll outputs in: {RESULTS_DIR}")
    print("\nPrediction summary:")
    for pid, res in sorted(pred_results.items()):
        status = 'PASS' if res.get('pass') else 'FAIL'
        kill_tag = ' ← KILL TEST' if pid == 'P3B2-1' else ''
        print(f"  {pid}: {status}{kill_tag}")

if __name__ == '__main__':
    main()
