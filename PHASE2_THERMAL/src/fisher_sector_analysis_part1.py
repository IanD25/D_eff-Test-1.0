#!/usr/bin/env python3
"""
DS Phase 3B-2: Sector-Level Fisher Analysis
Post-processing script for LEAN backtest results

Generates 9 figures and tests 6 predictions.
Kill test: P3B2-1 - Financials must lead Market in 2007-2008.
"""

import json
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from datetime import datetime
from scipy import stats
from pathlib import Path

# Configuration
RESULTS_DIR = Path("phase3b2_results")
RESULTS_DIR.mkdir(exist_ok=True)

# Crisis periods for analysis
CRISIS_PERIODS = {
    '2008_Financial': ('2007-01-01', '2009-12-31'),
    '2014_Oil': ('2014-01-01', '2016-12-31'),
    '2018_Volmageddon': ('2017-10-01', '2018-06-30'),
    '2020_COVID': ('2019-10-01', '2020-12-31'),
    '2022_RateHike': ('2021-10-01', '2023-03-31'),
}

# Key crisis dates for markers
CRISIS_DATES = {
    '2008-09-15': 'Lehman',
    '2015-01-01': 'Oil Crash',
    '2018-02-05': 'Volmageddon',
    '2020-03-16': 'COVID',
    '2022-06-01': 'Rate Hike',
}

SECTORS = ['Technology', 'Financials', 'HealthCare', 'Consumer', 
           'Industrials', 'Energy', 'Utilities']

def load_results(json_path):
    """Load results from LEAN ObjectStore JSON."""
    with open(json_path, 'r') as f:
        data = json.load(f)
    
    # Convert to DataFrame
    df = pd.DataFrame()
    df['date'] = pd.to_datetime(data['dates'])
    df['market'] = data['market']
    df['market_rank'] = data['market_rank']
    df['market_pr'] = data['market_pr']
    df['vix'] = data['vix']
    df['spy'] = data['spy']
    
    # Add sector data
    for sector in SECTORS:
        df[f'sector_{sector}'] = data[f'sector_{sector}']
        df[f'sector_{sector}_30d'] = data[f'sector_{sector}_30d']
        df[f'spread_{sector}'] = data[f'spread_{sector}']
        df[f'rank_{sector}'] = data[f'rank_{sector}']
        df[f'pr_{sector}'] = data[f'pr_{sector}']
        df[f'n_assets_{sector}'] = data[f'n_assets_{sector}']
    
    # Add asset temperature data (as dict column)
    df['asset_temps'] = [{} for _ in range(len(df))]
    for i, date in enumerate(data['dates']):
        for sector in SECTORS:
            temps = data[f'assets_{sector}'][i]
            if temps:
                for ticker, temp in temps.items():
                    df.at[i, 'asset_temps'][ticker] = temp
    
    return df, data

def compute_z_scores(df, series_name, window=252):
    """Compute rolling z-scores for a series."""
    series = df[series_name].copy()
    rolling_mean = series.rolling(window=window, min_periods=60).mean()
    rolling_std = series.rolling(window=window, min_periods=60).std()
    z_scores = (series - rolling_mean) / rolling_std
    return z_scores

def find_crossing_dates(z_series, threshold=2):
    """Find dates where z-score crosses threshold."""
    above = z_series > threshold
    crossings = above.astype(int).diff()
    up_crossings = df[df.index.isin(crossings[crossings == 1].index)]
    down_crossings = df[df.index.isin(crossings[crossings == -1].index)]
    return up_crossings, down_crossings

def figure1_sector_timeseries(df):
    """Figure 1: All 7 sector SV2/SV1 + market, overlaid."""
    fig, ax = plt.subplots(figsize=(14, 8))
    
    # Plot market
    ax.plot(df['date'], df['market'], label='Market', color='black', linewidth=2.5)
    
    # Plot sectors with distinct colors
    colors = plt.cm.Set3(np.linspace(0, 1, len(SECTORS)))
    for sector, color in zip(SECTORS, colors):
        ax.plot(df['date'], df[f'sector_{sector}'], label=sector, 
                color=color, alpha=0.8, linewidth=1.5)
    
    # Add crisis markers
    for crisis_date, label in CRISIS_DATES.items():
        date = pd.to_datetime(crisis_date)
        if date >= df['date'].min() and date <= df['date'].max():
            ax.axvline(date, color='red', alpha=0.3, linestyle='--')
            ax.text(date, ax.get_ylim()[1]*0.95, label, rotation=90, 
                   verticalalignment='top', fontsize=9)
    
    ax.set_xlabel('Date')
    ax.set_ylabel('SV2/SV1 (90-day window)')
    ax.set_title('Sector-Level Fisher SV2/SV1 (2005-2024)')
    ax.legend(loc='upper left', ncol=2)
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(RESULTS_DIR / 'sector_sv2_sv1_timeseries.png', dpi=300)
    plt.close()
    
    print("✓ Figure 1: sector_sv2_sv1_timeseries.png")

def figure2_financials_vs_market_2008(df):
    """Figure 2: KILL TEST - Financials vs Market 2007-2009."""
    mask = (df['date'] >= '2006-01-01') & (df['date'] <= '2009-12-31')
    df_sub = df[mask].copy()
    
    # Compute z-scores
    df_sub['market_z'] = compute_z_scores(df_sub, 'market')
    df_sub['financials_z'] = compute_z_scores(df_sub, 'sector_Financials')
    
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(14, 10), sharex=True)
    
    # Panel A: Raw SV2/SV1
    ax1.plot(df_sub['date'], df_sub['market'], label='Market', color='black', linewidth=2)
    ax1.plot(df_sub['date'], df_sub['sector_Financials'], label='Financials', 
             color='blue', linewidth=2)
    ax1.axvline(pd.to_datetime('2008-09-15'), color='red', alpha=0.5, linestyle='--', label='Lehman')
    ax1.set_ylabel('SV2/SV1')
    ax1.set_title('Financials vs Market SV2/SV1 (2006-2009)')
    ax1.legend()
    ax1.grid(True, alpha=0.3)
    
    # Panel B: Z-scores
    ax2.plot(df_sub['date'], df_sub['market_z'], label='Market Z', color='black', linewidth=2)
    ax2.plot(df_sub['date'], df_sub['financials_z'], label='Financials Z', 
             color='blue', linewidth=2)
    ax2.axhline(2, color='red', alpha=0.5, linestyle='--', label='2σ threshold')
    ax2.axhline(-2, color='red', alpha=0.5, linestyle='--')
    ax2.axvline(pd.to_datetime('2008-09-15'), color='red', alpha=0.5, linestyle='--')
    ax2.set_ylabel('Z-score')
    ax2.set_xlabel('Date')
    ax2.legend()
    ax2.grid(True, alpha=0.3)
    
    # Find crossing dates
    fin_up, _ = find_crossing_dates(df_sub['financials_z'], threshold=2)
    mkt_up, _ = find_crossing_dates(df_sub['market_z'], threshold=2)
    
    if not fin_up.empty and not mkt_up.empty:
        fin_first = fin_up['date'].min()
        mkt_first = mkt_up['date'].min()
        lead_days = (mkt_first - fin_first).days
        
        ax2.text(0.02, 0.95, f'Financials leads by {lead_days} days', 
                transform=ax2.transAxes, fontsize=10, 
                bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))
        
        # Mark crossing points
        ax2.scatter(fin_first, 2, color='blue', s=100, zorder=5, 
                   label='Fin > 2σ')
        ax2.scatter(mkt_first, 2, color='black', s=100, zorder=5, 
                   label='Mkt > 2σ')
    
    plt.tight_layout()
    plt.savefig(RESULTS_DIR / 'financials_vs_market_2007_2008.png', dpi=300)
    plt.close()
    
    print("✓ Figure 2: financials_vs_market_2007_2008.png")
    
    # Return lead/lag analysis for kill test
    if not fin_up.empty and not mkt_up.empty:
        return lead_days
    return None

def figure3_sector_lead_lag_heatmap(df):
    """Figure 3: Heatmap of sector lead/lag for each crisis."""
    # Compute z-scores for all sectors
    z_data = {}
    for sector in SECTORS:
        z_data[sector] = compute_z_scores(df, f'sector_{sector}')
    
    # For each crisis period, find first crossing > 2σ
    lead_lag_matrix = pd.DataFrame(index=SECTORS, 
                                   columns=list(CRISIS_PERIODS.keys()))
    
    for crisis, (start, end) in CRISIS_PERIODS.items():
        mask = (df['date'] >= start) & (df['date'] <= end)
        if not mask.any():
            continue
            
        for sector in SECTORS:
            z_series = z_data[sector][mask]
            above = z_series > 2
            if above.any():
                first_crossing = df.loc[mask & above, 'date'].min()
                # Days from start of crisis period
                start_date = pd.to_datetime(start)
                lead_lag = (first_crossing - start_date).days
                lead_lag_matrix.loc[sector, crisis] = lead_lag
            else:
                lead_lag_matrix.loc[sector, crisis] = np.nan
    
    # Create heatmap
    fig, ax = plt.subplots(figsize=(12, 8))
    im = ax.imshow(lead_lag_matrix.astype(float).values, cmap='RdYlGn_r', 
                   aspect='auto')
    
    # Add text
    for i in range(len(SECTORS)):
        for j in range(len(CRISIS_PERIODS)):
            value = lead_lag_matrix.iloc[i, j]
            if not np.isnan(value):
                ax.text(j, i, f'{int(value)}', ha='center', va='center', 
                       color='black', fontsize=9)
    
    ax.set_xticks(range(len(CRISIS_PERIODS)))
    ax.set_xticklabels([c.replace('_', '\n') for c in CRISIS_PERIODS.keys()], 
                       rotation=45, ha='right')
    ax.set_yticks(range(len(SECTORS)))
    ax.set_yticklabels(SECTORS)
    ax.set_title('Sector Lead/Lag (Days to exceed 2σ after crisis start)')
    
    plt.colorbar(im, ax=ax, label='Days')
    plt.tight_layout()
    plt.savefig(RESULTS_DIR / 'sector_lead_lag_heatmap.png', dpi=300)
    plt.close()
    
    print("✓ Figure 3: sector_lead_lag_heatmap.png")
    return lead_lag_matrix