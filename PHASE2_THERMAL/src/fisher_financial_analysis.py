#!/usr/bin/env python3
"""
Phase 3B: Financial Correlation Network — Post-Processing Script
================================================================
Reads the JSON output from the LEAN FisherDiagnosticAlgorithm and generates:
  - 9 matplotlib figures
  - PHASE3B_FINANCIAL_RESULTS.md handback report

Usage:
    python3 fisher_financial_analysis.py <path_to_json>

The JSON file is saved by the LEAN algorithm to ObjectStore. Copy it to a
local path before running this script.
"""

import json
import sys
import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from datetime import datetime, timedelta

# ============================================================================
# Helper functions
# ============================================================================

def parse_dates(date_strings):
    """Convert list of date strings to datetime objects."""
    dates = []
    for d in date_strings:
        try:
            dates.append(datetime.strptime(d.strip(), "%Y-%m-%d"))
        except ValueError:
            try:
                dates.append(datetime.strptime(d.strip()[:10], "%Y-%m-%d"))
            except ValueError:
                dates.append(None)
    return dates


def extract_metric(results_list, key, default=float('nan')):
    """Extract a metric from results list, handling None/empty entries."""
    values = []
    for r in results_list:
        if r and key in r:
            values.append(float(r[key]))
        else:
            values.append(default)
    return np.array(values)


def rolling_zscore(values, window=252):
    """Compute rolling z-score: (x - rolling_mean) / rolling_std."""
    n = len(values)
    zscores = np.full(n, np.nan)
    for i in range(window, n):
        segment = values[i - window:i]
        valid = segment[~np.isnan(segment)]
        if len(valid) > 10:
            mu = np.mean(valid)
            sigma = np.std(valid)
            if sigma > 1e-10:
                zscores[i] = (values[i] - mu) / sigma
    return zscores


def find_crisis_dates():
    """Return key crisis reference dates."""
    return {
        '2008_lehman': datetime(2008, 9, 15),
        '2020_covid': datetime(2020, 3, 11),
        '2018_volmageddon': datetime(2018, 2, 5),
        '2011_eu_debt': datetime(2011, 8, 5),
        '2015_china': datetime(2015, 8, 24),
        '2022_rate_hike': datetime(2022, 6, 13),
    }


# ============================================================================
# Figure generators
# ============================================================================

def figure1_sv2_sv1_vs_vix(dates, sv2_sv1, vix, outdir):
    """Figure 1: HEADLINE — SV2/SV1 (90d) vs VIX with crisis markers."""
    fig, ax1 = plt.subplots(figsize=(14, 6))

    # SV2/SV1
    color1 = 'tab:blue'
    ax1.set_xlabel('Date')
    ax1.set_ylabel('SV2/SV1 (90d)', color=color1)
    ax1.plot(dates, sv2_sv1, color=color1, alpha=0.7, linewidth=0.8, label='SV2/SV1')
    ax1.tick_params(axis='y', labelcolor=color1)

    # VIX on secondary axis
    ax2 = ax1.twinx()
    color2 = 'tab:red'
    ax2.set_ylabel('VIX', color=color2)
    vix_clean = np.where(np.isnan(vix), 0, vix)
    ax2.plot(dates, vix_clean, color=color2, alpha=0.5, linewidth=0.5, label='VIX')
    ax2.tick_params(axis='y', labelcolor=color2)

    # Crisis markers
    crises = find_crisis_dates()
    for name, date in crises.items():
        if dates[0] <= date <= dates[-1]:
            ax1.axvline(x=date, color='gray', linestyle='--', alpha=0.5, linewidth=0.8)
            ax1.text(date, ax1.get_ylim()[1] * 0.95, name.split('_')[0],
                     rotation=90, fontsize=7, va='top', ha='right', alpha=0.6)

    # Rolling z-score threshold
    zscores = rolling_zscore(sv2_sv1, window=52)  # ~1 year of weekly data
    threshold_mask = zscores > 2.0
    if np.any(threshold_mask):
        ax1.scatter(np.array(dates)[threshold_mask],
                   sv2_sv1[threshold_mask],
                   color='orange', s=15, zorder=5, alpha=0.7,
                   label='SV2/SV1 > 2σ above trailing mean')

    ax1.set_title('Phase 3B: SV2/SV1 (90d) vs VIX — Financial Fisher Diagnostic')
    ax1.legend(loc='upper left', fontsize=8)
    ax1.xaxis.set_major_formatter(mdates.DateFormatter('%Y'))
    ax1.xaxis.set_major_locator(mdates.YearLocator(2))
    fig.tight_layout()
    fig.savefig(os.path.join(outdir, 'sv2_sv1_vs_vix.png'), dpi=150)
    plt.close(fig)
    print("  Figure 1: sv2_sv1_vs_vix.png")


def figure2_diagnostic_panel(dates, data_90, outdir):
    """Figure 2: 4-panel diagnostic time series."""
    sv2_sv1 = extract_metric(data_90, 'sv2_sv1')
    eta = extract_metric(data_90, 'eta')
    gap_ratio = extract_metric(data_90, 'gap_ratio')
    rank = extract_metric(data_90, 'rank')

    fig, axes = plt.subplots(4, 1, figsize=(14, 10), sharex=True)

    axes[0].plot(dates, sv2_sv1, 'b-', linewidth=0.8)
    axes[0].set_ylabel('SV2/SV1')
    axes[0].set_title('Fisher Diagnostics (90-day window)')

    axes[1].plot(dates, eta, 'g-', linewidth=0.8)
    axes[1].set_ylabel('η (disorder index)')

    axes[2].plot(dates, gap_ratio, 'r-', linewidth=0.8)
    axes[2].set_ylabel('Gap Ratio (SV1/SV2)')
    axes[2].set_yscale('log')

    axes[3].plot(dates, rank, 'k-', linewidth=0.8)
    axes[3].set_ylabel('Rank')
    axes[3].set_xlabel('Date')

    # Crisis markers on all panels
    crises = find_crisis_dates()
    for ax in axes:
        for name, date in crises.items():
            if dates[0] <= date <= dates[-1]:
                ax.axvline(x=date, color='gray', linestyle='--', alpha=0.4, linewidth=0.5)

    for ax in axes:
        ax.xaxis.set_major_formatter(mdates.DateFormatter('%Y'))
        ax.xaxis.set_major_locator(mdates.YearLocator(2))

    fig.tight_layout()
    fig.savefig(os.path.join(outdir, 'diagnostic_panel_timeseries.png'), dpi=150)
    plt.close(fig)
    print("  Figure 2: diagnostic_panel_timeseries.png")


def figure3_multiscale(dates, all_results, windows, outdir):
    """Figure 3: SV2/SV1 at all 4 windows overlaid."""
    fig, ax = plt.subplots(figsize=(14, 6))

    colors = ['tab:blue', 'tab:green', 'tab:orange', 'tab:red']
    for i, w in enumerate(windows):
        sv2_sv1 = extract_metric(all_results[w], 'sv2_sv1')
        ax.plot(dates, sv2_sv1, color=colors[i % len(colors)],
                linewidth=0.7, alpha=0.8, label=f'{w}d window')

    crises = find_crisis_dates()
    for name, date in crises.items():
        if dates[0] <= date <= dates[-1]:
            ax.axvline(x=date, color='gray', linestyle='--', alpha=0.4, linewidth=0.5)

    ax.set_title('Multi-Scale SV2/SV1 — All Windows')
    ax.set_xlabel('Date')
    ax.set_ylabel('SV2/SV1')
    ax.legend(fontsize=9)
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%Y'))
    ax.xaxis.set_major_locator(mdates.YearLocator(2))

    fig.tight_layout()
    fig.savefig(os.path.join(outdir, 'multiscale_sv2_sv1.png'), dpi=150)
    plt.close(fig)
    print("  Figure 3: multiscale_sv2_sv1.png")


def figure4_pr_window_profile(dates, all_results, windows, outdir):
    """Figure 4: PR vs window length at key dates."""
    crises = find_crisis_dates()
    key_dates_map = {
        'Pre-2008 (2008-06)': datetime(2008, 6, 1),
        '2008 Crisis (2008-10)': datetime(2008, 10, 1),
        'Post-2008 (2009-06)': datetime(2009, 6, 1),
        'Pre-COVID (2020-01)': datetime(2020, 1, 15),
        'COVID Crisis (2020-03)': datetime(2020, 3, 15),
        'Calm (2017-06)': datetime(2017, 6, 1),
    }

    fig, ax = plt.subplots(figsize=(10, 6))
    colors = ['blue', 'red', 'green', 'cyan', 'magenta', 'gray']

    for ci, (label, target_date) in enumerate(key_dates_map.items()):
        # Find closest date
        best_idx = None
        best_dist = float('inf')
        for i, d in enumerate(dates):
            dist = abs((d - target_date).days)
            if dist < best_dist:
                best_dist = dist
                best_idx = i

        if best_idx is None or best_dist > 30:
            continue

        prs = []
        for w in windows:
            r = all_results[w][best_idx] if best_idx < len(all_results[w]) else None
            if r and 'pr' in r:
                prs.append(r['pr'])
            else:
                prs.append(float('nan'))

        ax.plot(windows, prs, 'o-', color=colors[ci % len(colors)],
                label=label, markersize=6, linewidth=1.5)

    ax.set_xlabel('Rolling Window (days)')
    ax.set_ylabel('Participation Ratio')
    ax.set_title('Participation Ratio vs Window Length at Key Dates')
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3)

    fig.tight_layout()
    fig.savefig(os.path.join(outdir, 'pr_window_profile.png'), dpi=150)
    plt.close(fig)
    print("  Figure 4: pr_window_profile.png")


def figure_crisis_zoom(dates, data_90, spy_values, vix_values,
                        crisis_date, crisis_name, filename, outdir,
                        window_days=180):
    """Generic crisis zoom figure."""
    start = crisis_date - timedelta(days=window_days)
    end = crisis_date + timedelta(days=window_days)

    # Filter to window
    mask = [(start <= d <= end) for d in dates]
    idx = [i for i, m in enumerate(mask) if m]

    if len(idx) < 5:
        print(f"  Skipping {filename} — insufficient data near {crisis_name}")
        return

    zoom_dates = [dates[i] for i in idx]
    zoom_sv2 = extract_metric([data_90[i] for i in idx], 'sv2_sv1')
    zoom_spy = np.array([spy_values[i] for i in idx])
    zoom_vix = np.array([vix_values[i] for i in idx])

    fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(12, 8), sharex=True)

    # SPY
    ax1.plot(zoom_dates, zoom_spy, 'k-', linewidth=1)
    ax1.set_ylabel('SPY Price')
    ax1.set_title(f'Crisis Zoom: {crisis_name}')
    ax1.axvline(x=crisis_date, color='red', linestyle='--', alpha=0.7)

    # SV2/SV1
    ax2.plot(zoom_dates, zoom_sv2, 'b-', linewidth=1.2)
    ax2.set_ylabel('SV2/SV1 (90d)')
    ax2.axvline(x=crisis_date, color='red', linestyle='--', alpha=0.7)

    # Rolling z-score for this window
    full_sv2 = extract_metric(data_90, 'sv2_sv1')
    full_zscores = rolling_zscore(full_sv2, window=52)
    zoom_z = np.array([full_zscores[i] if i < len(full_zscores) else np.nan for i in idx])
    ax2_twin = ax2.twinx()
    ax2_twin.plot(zoom_dates, zoom_z, 'orange', linewidth=0.8, alpha=0.7)
    ax2_twin.axhline(y=2.0, color='orange', linestyle=':', alpha=0.5)
    ax2_twin.set_ylabel('Z-score', color='orange')

    # VIX
    ax3.plot(zoom_dates, zoom_vix, 'r-', linewidth=1)
    ax3.set_ylabel('VIX')
    ax3.set_xlabel('Date')
    ax3.axvline(x=crisis_date, color='red', linestyle='--', alpha=0.7)

    for ax in [ax1, ax2, ax3]:
        ax.xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m'))

    fig.tight_layout()
    fig.savefig(os.path.join(outdir, filename), dpi=150)
    plt.close(fig)
    print(f"  Figure: {filename}")


def figure8_sv_profile_snapshots(dates, data_90, outdir):
    """Figure 8: SV bar charts at 6 key dates."""
    key_dates_map = {
        'Pre-2008\n(2008-06)': datetime(2008, 6, 1),
        '2008 Crisis\n(2008-10)': datetime(2008, 10, 1),
        'Recovery\n(2009-06)': datetime(2009, 6, 1),
        'Calm\n(2017-06)': datetime(2017, 6, 1),
        'Pre-COVID\n(2020-01)': datetime(2020, 1, 15),
        'COVID\n(2020-03)': datetime(2020, 3, 15),
    }

    fig, axes = plt.subplots(2, 3, figsize=(14, 8))
    axes_flat = axes.flatten()

    for ai, (label, target_date) in enumerate(key_dates_map.items()):
        ax = axes_flat[ai]

        # Find closest date
        best_idx = None
        best_dist = float('inf')
        for i, d in enumerate(dates):
            dist = abs((d - target_date).days)
            if dist < best_dist:
                best_dist = dist
                best_idx = i

        if best_idx is None or best_dist > 30:
            ax.set_title(f'{label}\nNo data', fontsize=9)
            continue

        r = data_90[best_idx]
        if r and 'sv_profile' in r and len(r['sv_profile']) > 0:
            profile = r['sv_profile']
            x = np.arange(len(profile))
            ax.bar(x, profile, color='steelblue', alpha=0.8)
            ax.set_title(f'{label}\nSV2/SV1={r["sv2_sv1"]:.3f}', fontsize=9)
            ax.set_xlabel('SV index')
            ax.set_ylabel('Normalized SV')
            ax.set_ylim(0, 1.1)
        else:
            ax.set_title(f'{label}\nNo SV data', fontsize=9)

    fig.suptitle('SV Profile Snapshots at Key Dates (90d window)', fontsize=12)
    fig.tight_layout()
    fig.savefig(os.path.join(outdir, 'sv_profile_snapshots.png'), dpi=150)
    plt.close(fig)
    print("  Figure 8: sv_profile_snapshots.png")


def figure9_correlation_function(dates, data_90, outdir):
    """Figure 9: C(r) at key dates."""
    key_dates_map = {
        'Pre-2008 (2008-06)': datetime(2008, 6, 1),
        '2008 Crisis (2008-10)': datetime(2008, 10, 1),
        'Calm (2017-06)': datetime(2017, 6, 1),
        'COVID (2020-03)': datetime(2020, 3, 15),
    }

    fig, ax = plt.subplots(figsize=(10, 6))
    colors = ['blue', 'red', 'gray', 'magenta']

    for ci, (label, target_date) in enumerate(key_dates_map.items()):
        best_idx = None
        best_dist = float('inf')
        for i, d in enumerate(dates):
            dist = abs((d - target_date).days)
            if dist < best_dist:
                best_dist = dist
                best_idx = i

        if best_idx is None or best_dist > 30:
            continue

        r = data_90[best_idx]
        if r and 'cr_values' in r:
            cr = r['cr_values']
            ax.plot(range(len(cr)), cr, 'o-', color=colors[ci % len(colors)],
                    label=label, markersize=4, linewidth=1.2)

    ax.set_xlabel('Graph Distance r')
    ax.set_ylabel('C(r) — Average Correlation')
    ax.set_title('Correlation Function C(r) at Key Dates')
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3)
    ax.axhline(y=0, color='black', linestyle='-', alpha=0.3)

    fig.tight_layout()
    fig.savefig(os.path.join(outdir, 'correlation_function_key_dates.png'), dpi=150)
    plt.close(fig)
    print("  Figure 9: correlation_function_key_dates.png")


# ============================================================================
# Prediction assessment
# ============================================================================

def assess_predictions(dates, all_results, windows, spy_values, vix_values):
    """Assess all pre-registered predictions."""
    crises = find_crisis_dates()
    data_90 = all_results[90]
    sv2_sv1_90 = extract_metric(data_90, 'sv2_sv1')
    rank_90 = extract_metric(data_90, 'rank')
    eta_90 = extract_metric(data_90, 'eta')

    # Rolling z-score (52-week trailing)
    zscores = rolling_zscore(sv2_sv1_90, window=52)

    results = {}

    # --- P3B-1: SV2/SV1 rises before crashes (KILL TEST) ---
    p3b1_pass = False
    for crisis_name, crisis_date in [('2008_lehman', crises['2008_lehman']),
                                       ('2020_covid', crises['2020_covid'])]:
        # Check 30 days before crisis
        for i, d in enumerate(dates):
            days_before = (crisis_date - d).days
            if 0 < days_before <= 30:
                if not np.isnan(zscores[i]) and zscores[i] > 2.0:
                    p3b1_pass = True
                    results[f'P3B-1_{crisis_name}'] = {
                        'pass': True,
                        'date': str(d.date()),
                        'zscore': float(zscores[i]),
                        'sv2_sv1': float(sv2_sv1_90[i]),
                        'days_before': days_before
                    }
                    break

    # Also check 60-day window for softer signal
    for crisis_name, crisis_date in [('2008_lehman', crises['2008_lehman']),
                                       ('2020_covid', crises['2020_covid'])]:
        for i, d in enumerate(dates):
            days_before = (crisis_date - d).days
            if 0 < days_before <= 60:
                if not np.isnan(zscores[i]) and zscores[i] > 1.5:
                    key = f'P3B-1_{crisis_name}_60d'
                    if key not in results:
                        results[key] = {
                            'pass': True,
                            'date': str(d.date()),
                            'zscore': float(zscores[i]),
                            'sv2_sv1': float(sv2_sv1_90[i]),
                            'days_before': days_before
                        }
                    break

    results['P3B-1_verdict'] = 'PASS' if p3b1_pass else 'FAIL (KILL)'

    # --- P3B-2: Rank transitions before crashes ---
    p3b2_pass = False
    for crisis_name, crisis_date in [('2008_lehman', crises['2008_lehman']),
                                       ('2020_covid', crises['2020_covid'])]:
        for i, d in enumerate(dates):
            days_before = (crisis_date - d).days
            if 0 < days_before <= 60:
                if not np.isnan(rank_90[i]) and rank_90[i] > 2:
                    p3b2_pass = True
                    results[f'P3B-2_{crisis_name}'] = {
                        'pass': True,
                        'date': str(d.date()),
                        'rank': float(rank_90[i]),
                        'days_before': days_before
                    }
                    break

    results['P3B-2_verdict'] = 'PASS' if p3b2_pass else 'FAIL'

    # --- P3B-3: eta feature near crashes ---
    p3b3_pass = False
    for crisis_name, crisis_date in [('2008_lehman', crises['2008_lehman']),
                                       ('2020_covid', crises['2020_covid'])]:
        for i, d in enumerate(dates):
            days_diff = abs((d - crisis_date).days)
            if days_diff <= 30 and not np.isnan(eta_90[i]):
                # Check if local extremum (compare to neighbors)
                if i > 0 and i < len(eta_90) - 1:
                    is_max = (eta_90[i] > eta_90[i - 1] and eta_90[i] > eta_90[i + 1])
                    is_min = (eta_90[i] < eta_90[i - 1] and eta_90[i] < eta_90[i + 1])
                    if is_max or is_min:
                        p3b3_pass = True
                        results[f'P3B-3_{crisis_name}'] = {
                            'pass': True,
                            'date': str(d.date()),
                            'eta': float(eta_90[i]),
                            'type': 'maximum' if is_max else 'minimum'
                        }
                        break

    results['P3B-3_verdict'] = 'PASS' if p3b3_pass else 'FAIL'

    # --- P3B-4: 2008 peak > 2018 peak ---
    p3b4_pass = False
    # Find max SV2/SV1 within 60 days of each crisis
    crisis_peaks = {}
    for crisis_name, crisis_date in [('2008', crises['2008_lehman']),
                                       ('2018', crises['2018_volmageddon'])]:
        peak = -np.inf
        for i, d in enumerate(dates):
            if abs((d - crisis_date).days) <= 60:
                if not np.isnan(sv2_sv1_90[i]):
                    peak = max(peak, sv2_sv1_90[i])
        crisis_peaks[crisis_name] = peak

    if crisis_peaks.get('2008', -np.inf) > crisis_peaks.get('2018', -np.inf):
        p3b4_pass = True
    results['P3B-4'] = {
        'pass': p3b4_pass,
        '2008_peak': crisis_peaks.get('2008', float('nan')),
        '2018_peak': crisis_peaks.get('2018', float('nan')),
    }
    results['P3B-4_verdict'] = 'PASS' if p3b4_pass else 'FAIL'

    # --- P3B-5: Short-window leads long-window ---
    p3b5_pass = False
    sv2_30 = extract_metric(all_results[30], 'sv2_sv1')
    sv2_90_arr = sv2_sv1_90
    zscores_30 = rolling_zscore(sv2_30, window=52)
    zscores_90 = rolling_zscore(sv2_90_arr, window=52)

    for crisis_name, crisis_date in [('2008_lehman', crises['2008_lehman']),
                                       ('2020_covid', crises['2020_covid'])]:
        # Find first date where 30d z-score > 1.5 within 90 days before crisis
        date_30 = None
        date_90 = None
        for i, d in enumerate(dates):
            days_before = (crisis_date - d).days
            if 0 < days_before <= 90:
                if date_30 is None and not np.isnan(zscores_30[i]) and zscores_30[i] > 1.5:
                    date_30 = d
                if date_90 is None and not np.isnan(zscores_90[i]) and zscores_90[i] > 1.5:
                    date_90 = d

        if date_30 is not None and date_90 is not None:
            lead_days = (date_90 - date_30).days
            if lead_days >= 5:
                p3b5_pass = True
                results[f'P3B-5_{crisis_name}'] = {
                    'pass': True,
                    'lead_days': lead_days,
                    'date_30d': str(date_30.date()),
                    'date_90d': str(date_90.date()),
                }

    results['P3B-5_verdict'] = 'PASS' if p3b5_pass else 'FAIL'

    return results


# ============================================================================
# Report generator
# ============================================================================

def write_results_report(dates, all_results, windows, spy_values, vix_values,
                          predictions, outdir, json_path):
    """Write PHASE3B_FINANCIAL_RESULTS.md."""
    data_90 = all_results[90]
    sv2_sv1_90 = extract_metric(data_90, 'sv2_sv1')
    rank_90 = extract_metric(data_90, 'rank')
    eta_90 = extract_metric(data_90, 'eta')

    valid_mask = ~np.isnan(sv2_sv1_90)
    n_valid = np.sum(valid_mask)
    n_total = len(dates)

    report = f"""# Phase 3B Results: Financial Correlation Network Fisher Diagnostics

Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}

## 1. Data Summary

- Date range: {dates[0].strftime('%Y-%m-%d')} to {dates[-1].strftime('%Y-%m-%d')}
- Total computation dates: {n_total}
- Valid 90d results: {n_valid} ({100*n_valid/max(n_total,1):.1f}%)
- Windows: {windows}
- k_nn: 10
- Source: {json_path}

## 2. SV2/SV1 (90d) Summary Statistics

| Metric | Value |
|--------|-------|
| Mean | {np.nanmean(sv2_sv1_90):.4f} |
| Std | {np.nanstd(sv2_sv1_90):.4f} |
| Min | {np.nanmin(sv2_sv1_90):.4f} |
| Max | {np.nanmax(sv2_sv1_90):.4f} |
| Median | {np.nanmedian(sv2_sv1_90):.4f} |

## 3. Rank Summary (90d)

| Metric | Value |
|--------|-------|
| Mean | {np.nanmean(rank_90):.2f} |
| Std | {np.nanstd(rank_90):.2f} |
| Min | {np.nanmin(rank_90):.1f} |
| Max | {np.nanmax(rank_90):.1f} |

## 4. Prediction Assessment

| ID | Prediction | Result | Notes |
|---|---|---|---|
| **P3B-1** | **SV2/SV1 rises before crashes** | **{predictions['P3B-1_verdict']}** | Kill test |
| P3B-2 | Rank transitions before crashes | {predictions['P3B-2_verdict']} | |
| P3B-3 | eta feature near crashes | {predictions['P3B-3_verdict']} | |
| P3B-4 | 2008 peak > 2018 peak | {predictions['P3B-4_verdict']} | |
| P3B-5 | Short-window leads long-window | {predictions['P3B-5_verdict']} | |

### P3B-1 Detail (Kill Test)
"""

    for key in sorted(predictions.keys()):
        if key.startswith('P3B-1_') and key != 'P3B-1_verdict':
            val = predictions[key]
            if isinstance(val, dict):
                report += f"\n**{key}:** "
                for k2, v2 in val.items():
                    report += f"{k2}={v2}, "
                report += "\n"

    report += f"""
### P3B-4 Detail
- 2008 peak SV2/SV1: {predictions['P3B-4'].get('2008_peak', 'N/A')}
- 2018 peak SV2/SV1: {predictions['P3B-4'].get('2018_peak', 'N/A')}

## 5. Crisis-Proximate Behavior

### Near 2008 Lehman (Sep 15, 2008)
"""
    # Find data near 2008
    lehman = datetime(2008, 9, 15)
    for i, d in enumerate(dates):
        days = (d - lehman).days
        if -60 <= days <= 30:
            r = data_90[i]
            if r and 'sv2_sv1' in r:
                report += f"| {d.strftime('%Y-%m-%d')} | SV2/SV1={r['sv2_sv1']:.3f} | rank={r['rank']:.1f} | eta={r['eta']:.3f} |\n"

    report += f"""
### Near COVID (Mar 11, 2020)
"""
    covid = datetime(2020, 3, 11)
    for i, d in enumerate(dates):
        days = (d - covid).days
        if -60 <= days <= 30:
            r = data_90[i]
            if r and 'sv2_sv1' in r:
                report += f"| {d.strftime('%Y-%m-%d')} | SV2/SV1={r['sv2_sv1']:.3f} | rank={r['rank']:.1f} | eta={r['eta']:.3f} |\n"

    report += """
## 6. Figures

| # | File | Description |
|---|------|-------------|
| 1 | sv2_sv1_vs_vix.png | HEADLINE: SV2/SV1 vs VIX with crisis markers |
| 2 | diagnostic_panel_timeseries.png | 4-panel diagnostic time series |
| 3 | multiscale_sv2_sv1.png | Multi-scale SV2/SV1 |
| 4 | pr_window_profile.png | PR vs window at key dates |
| 5 | crisis_zoom_2008.png | 2008 crisis zoom |
| 6 | crisis_zoom_2020.png | COVID crisis zoom |
| 7 | crisis_zoom_2018.png | 2018 Volmageddon zoom |
| 8 | sv_profile_snapshots.png | SV bar charts at key dates |
| 9 | correlation_function_key_dates.png | C(r) at key dates |

## 7. Comparison to Physics Results

| System | Context | SV2/SV1 at transition | Swap? |
|--------|---------|----------------------|-------|
| 2D Ising | T_c (continuous) | 1.000 | YES |
| 3D Ising | T_c (continuous) | 0.956 | YES |
| Potts q=5 | T_c (1st-order, weak) | 0.691 | NO |
| Potts q=10 | T_c (1st-order, strong) | 0.374 | NO |
| Financial | Crisis (if applicable) | See above | TBD |

"""

    report_path = os.path.join(outdir, 'PHASE3B_FINANCIAL_RESULTS.md')
    with open(report_path, 'w') as f:
        f.write(report)
    print(f"  Report: PHASE3B_FINANCIAL_RESULTS.md")


# ============================================================================
# Main
# ============================================================================

def main():
    if len(sys.argv) < 2:
        print("Usage: python3 fisher_financial_analysis.py <path_to_json>")
        print("  The JSON file is the output from the LEAN FisherDiagnosticAlgorithm.")
        sys.exit(1)

    json_path = sys.argv[1]
    print(f"Loading results from: {json_path}")

    with open(json_path, 'r') as f:
        data = json.load(f)

    # Parse data
    dates = parse_dates(data['dates'])
    vix = np.array(data['vix'])
    spy = np.array(data['spy'])
    windows = data['windows']

    all_results = {}
    for w in windows:
        key = f'results_w{w}'
        if key in data:
            all_results[w] = data[key]
        else:
            all_results[w] = []

    print(f"  Dates: {len(dates)} ({dates[0].strftime('%Y-%m-%d')} to {dates[-1].strftime('%Y-%m-%d')})")
    for w in windows:
        valid = sum(1 for r in all_results[w] if r and 'sv2_sv1' in r)
        print(f"  Window {w}d: {valid}/{len(all_results[w])} valid results")

    # Create output directory
    outdir = os.path.join(os.path.dirname(json_path), '..', 'phase3b_results')
    outdir = os.path.abspath(outdir)
    os.makedirs(outdir, exist_ok=True)
    os.makedirs(os.path.join(outdir, 'raw_data'), exist_ok=True)

    # Copy raw JSON
    import shutil
    shutil.copy2(json_path, os.path.join(outdir, 'raw_data', 'fisher_financial_results.json'))

    # Extract 90d data
    data_90 = all_results[90]
    sv2_sv1_90 = extract_metric(data_90, 'sv2_sv1')

    print("\nGenerating figures...")

    # Figure 1: SV2/SV1 vs VIX
    figure1_sv2_sv1_vs_vix(dates, sv2_sv1_90, vix, outdir)

    # Figure 2: Diagnostic panel
    figure2_diagnostic_panel(dates, data_90, outdir)

    # Figure 3: Multi-scale
    figure3_multiscale(dates, all_results, windows, outdir)

    # Figure 4: PR vs window
    figure4_pr_window_profile(dates, all_results, windows, outdir)

    # Figure 5: Crisis zoom 2008
    crises = find_crisis_dates()
    figure_crisis_zoom(dates, data_90, spy, vix,
                        crises['2008_lehman'], '2008 Lehman Brothers',
                        'crisis_zoom_2008.png', outdir)

    # Figure 6: Crisis zoom 2020
    figure_crisis_zoom(dates, data_90, spy, vix,
                        crises['2020_covid'], '2020 COVID-19',
                        'crisis_zoom_2020.png', outdir)

    # Figure 7: Crisis zoom 2018
    figure_crisis_zoom(dates, data_90, spy, vix,
                        crises['2018_volmageddon'], '2018 Volmageddon',
                        'crisis_zoom_2018.png', outdir)

    # Figure 8: SV profile snapshots
    figure8_sv_profile_snapshots(dates, data_90, outdir)

    # Figure 9: Correlation function
    figure9_correlation_function(dates, data_90, outdir)

    # Assess predictions
    print("\nAssessing predictions...")
    predictions = assess_predictions(dates, all_results, windows, spy, vix)
    for key in sorted(predictions.keys()):
        if key.endswith('_verdict'):
            print(f"  {key}: {predictions[key]}")

    # Write report
    print("\nWriting report...")
    write_results_report(dates, all_results, windows, spy, vix,
                          predictions, outdir, json_path)

    print(f"\nAll outputs saved to: {outdir}")
    print("Done.")


if __name__ == '__main__':
    main()
