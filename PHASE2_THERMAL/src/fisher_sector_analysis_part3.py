def figure8_sector_rank_comparison(df):
    """Figure 8: Sector rank vs market rank over time."""
    fig, axes = plt.subplots(4, 2, figsize=(16, 12))
    axes = axes.flatten()
    
    for idx, sector in enumerate(SECTORS):
        ax = axes[idx]
        
        # Plot sector rank
        ax.plot(df['date'], df[f'rank_{sector}'], label=f'{sector} Rank', 
                color='blue', linewidth=1.5, alpha=0.7)
        
        # Plot market rank
        ax.plot(df['date'], df['market_rank'], label='Market Rank', 
                color='red', linewidth=1.5, alpha=0.7)
        
        # Add crisis markers
        for crisis_date in CRISIS_DATES.keys():
            date = pd.to_datetime(crisis_date)
            if date >= df['date'].min() and date <= df['date'].max():
                ax.axvline(date, color='gray', alpha=0.2, linestyle='--')
        
        ax.set_title(f'{sector} vs Market Rank')
        ax.set_ylabel('Rank')
        ax.grid(True, alpha=0.3)
        
        if idx == 0:
            ax.legend(loc='upper right')
    
    # Remove empty subplot
    axes[-1].axis('off')
    
    plt.suptitle('Sector Fisher Rank vs Market Rank (Lower = More Unified)', fontsize=14)
    plt.tight_layout()
    plt.savefig(RESULTS_DIR / 'sector_rank_comparison.png', dpi=300)
    plt.close()
    
    print("✓ Figure 8: sector_rank_comparison.png")

def figure9_cross_sector_correlation(df):
    """Figure 9: Correlation matrix of sector SV2/SV1 time series."""
    # Extract sector SV2/SV1 series
    sector_series = {}
    for sector in SECTORS:
        series = df[f'sector_{sector}'].copy()
        # Fill NaN with forward/backward fill
        series = series.fillna(method='ffill').fillna(method='bfill')
        sector_series[sector] = series
    
    # Create correlation matrix
    sector_df = pd.DataFrame(sector_series)
    corr_matrix = sector_df.corr()
    
    # Compute crisis vs calm correlations
    crisis_corrs = {}
    calm_corrs = {}
    
    for crisis, (start, end) in CRISIS_PERIODS.items():
        mask = (df['date'] >= start) & (df['date'] <= end)
        if mask.sum() > 10:  # Need enough data points
            crisis_df = sector_df[mask]
            crisis_corr = crisis_df.corr().values.mean()
            crisis_corrs[crisis] = crisis_corr
    
    # Calm periods (exclude crisis periods)
    calm_mask = pd.Series(True, index=df.index)
    for start, end in CRISIS_PERIODS.values():
        calm_mask &= ~((df['date'] >= start) & (df['date'] <= end))
    
    if calm_mask.sum() > 10:
        calm_df = sector_df[calm_mask]
        calm_corr = calm_df.corr().values.mean()
    
    # Plot correlation matrix
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))
    
    # Panel A: Full period correlation
    im1 = ax1.imshow(corr_matrix.values, cmap='RdBu_r', vmin=-1, vmax=1)
    ax1.set_xticks(range(len(SECTORS)))
    ax1.set_yticks(range(len(SECTORS)))
    ax1.set_xticklabels([s[:4] for s in SECTORS], rotation=45, ha='right')
    ax1.set_yticklabels([s[:4] for s in SECTORS])
    ax1.set_title('Cross-Sector SV2/SV1 Correlation (2005-2024)')
    
    # Add correlation values
    for i in range(len(SECTORS)):
        for j in range(len(SECTORS)):
            ax1.text(j, i, f'{corr_matrix.iloc[i, j]:.2f}', 
                    ha='center', va='center', color='black', fontsize=8)
    
    plt.colorbar(im1, ax=ax1, label='Correlation')
    
    # Panel B: Crisis vs calm comparison
    if crisis_corrs and 'calm_corr' in locals():
        crises = list(crisis_corrs.keys())
        crisis_values = [crisis_corrs[c] for c in crises]
        
        x_pos = np.arange(len(crises) + 1)
        values = crisis_values + [calm_corr]
        labels = [c.replace('_', '\n') for c in crises] + ['Calm\nPeriod']
        
        bars = ax2.bar(x_pos, values, color=['red']*len(crises) + ['blue'])
        ax2.axhline(0.5, color='gray', linestyle='--', alpha=0.5, label='0.5 threshold')
        ax2.axhline(0.8, color='gray', linestyle='--', alpha=0.5, label='0.8 threshold')
        
        ax2.set_xticks(x_pos)
        ax2.set_xticklabels(labels, rotation=45, ha='right')
        ax2.set_ylabel('Mean Pairwise Correlation')
        ax2.set_title('Crisis vs Calm: Sector Correlation Convergence')
        ax2.legend()
        ax2.grid(True, alpha=0.3, axis='y')
    
    plt.tight_layout()
    plt.savefig(RESULTS_DIR / 'cross_sector_correlation_of_sv2sv1.png', dpi=300)
    plt.close()
    
    print("✓ Figure 9: cross_sector_correlation_of_sv2sv1.png")
    
    return corr_matrix, crisis_corrs, calm_corr if 'calm_corr' in locals() else None

def test_predictions(df, lead_lag_matrix, corr_matrix, crisis_corrs, calm_corr):
    """Test the 6 pre-registered predictions."""
    results = {}
    
    # P3B2-1: Financials leads Market in 2007-2008 (KILL TEST)
    mask_2008 = (df['date'] >= '2007-01-01') & (df['date'] <= '2009-12-31')
    df_2008 = df[mask_2008].copy()
    
    df_2008['market_z'] = compute_z_scores(df_2008, 'market')
    df_2008['financials_z'] = compute_z_scores(df_2008, 'sector_Financials')
    
    fin_up, _ = find_crossing_dates(df_2008['financials_z'], threshold=2)
    mkt_up, _ = find_crossing_dates(df_2008['market_z'], threshold=2)
    
    if not fin_up.empty and not mkt_up.empty:
        fin_first = fin_up['date'].min()
        mkt_first = mkt_up['date'].min()
        lead_days = (mkt_first - fin_first).days
        
        results['P3B2-1'] = {
            'lead_days': lead_days,
            'pass': lead_days >= 30,
            'kill': lead_days < 30
        }
    else:
        results['P3B2-1'] = {
            'lead_days': None,
            'pass': False,
            'kill': True  # No crossing = fail
        }
    
    # P3B2-2: Energy leads during 2014-2016 oil crash
    mask_2014 = (df['date'] >= '2014-01-01') & (df['date'] <= '2016-12-31')
    df_2014 = df[mask_2014].copy()
    
    df_2014['market_z'] = compute_z_scores(df_2014, 'market')
    df_2014['energy_z'] = compute_z_scores(df_2014, 'sector_Energy')
    
    energy_up, _ = find_crossing_dates(df_2014['energy_z'], threshold=2)
    mkt_up, _ = find_crossing_dates(df_2014['market_z'], threshold=2)
    
    if not energy_up.empty and not mkt_up.empty:
        energy_first = energy_up['date'].min()
        mkt_first = mkt_up['date'].min()
        lead_days = (mkt_first - energy_first).days
        
        results['P3B2-2'] = {
            'lead_days': lead_days,
            'pass': lead_days >= 30
        }
    else:
        results['P3B2-2'] = {
            'lead_days': None,
            'pass': False
        }
    
    # P3B2-3: Tech leads during 2022 rate selloff
    mask_2022 = (df['date'] >= '2021-10-01') & (df['date'] <= '2023-03-31')
    df_2022 = df[mask_2022].copy()
    
    df_2022['market_z'] = compute_z_scores(df_2022, 'market')
    df_2022['tech_z'] = compute_z_scores(df_2022, 'sector_Technology')
    
    tech_up, _ = find_crossing_dates(df_2022['tech_z'], threshold=2)
    mkt_up, _ = find_crossing_dates(df_2022['market_z'], threshold=2)
    
    if not tech_up.empty and not mkt_up.empty:
        tech_first = tech_up['date'].min()
        mkt_first = mkt_up['date'].min()
        lead_days = (mkt_first - tech_first).days
        
        results['P3B2-3'] = {
            'lead_days': lead_days,
            'pass': lead_days >= 30
        }
    else:
        results['P3B2-3'] = {
            'lead_days': None,
            'pass': False
        }
    
    # P3B2-4: Spread > 0 precedes SV2/SV1 peak
    crises_to_test = ['2008_Financial', '2014_Oil', '2020_COVID']
    spread_lead_success = 0
    
    for crisis in crises_to_test:
        start, end = CRISIS_PERIODS[crisis]
        mask = (df['date'] >= start) & (df['date'] <= end)
        if not mask.any():
            continue
        
        # Find leading sector for this crisis
        if crisis in lead_lag_matrix.columns:
            sector_lead = lead_lag_matrix[crisis].dropna()
            if not sector_lead.empty:
                leading_sector = sector_lead.idxmin()
                
                # Check if spread > 0 precedes peak
                crisis_df = df[mask].copy()
                sector_series = crisis_df[f'sector_{leading_sector}']
                spread_series = crisis_df[f'spread_{leading_sector}']
                
                # Find peak
                peak_idx = sector_series.idxmax()
                peak_date = crisis_df.loc[peak_idx, 'date']
                
                # Check 2 weeks before peak
                two_weeks_before = peak_date - pd.Timedelta(days=14)
                mask_before = crisis_df['date'] <= two_weeks_before
                
                if mask_before.any():
                    spread_before = spread_series[mask_before]
                    if (spread_before > 0).any():
                        spread_lead_success += 1
    
    results['P3B2-4'] = {
        'successes': spread_lead_success,
        'total_crises': len(crises_to_test),
        'pass': spread_lead_success >= 2
    }
    
    # P3B2-5: High Fisher-temperature assets cluster in leading sector
    # (Simplified test - would need more detailed asset-level analysis)
    results['P3B2-5'] = {
        'pass': 'NEEDS_ASSET_LEVEL_ANALYSIS',
        'note': 'Requires detailed asset temperature tracking per sector'
    }
    
    # P3B2-6: Cross-sector correlation increases during crises
    if crisis_corrs and calm_corr is not None:
        crisis_avg = np.mean(list(crisis_corrs.values()))
        results['P3B2-6'] = {
            'crisis_avg': crisis_avg,
            'calm_avg': calm_corr,
            'pass': crisis_avg > 0.8 and calm_corr < 0.5
        }
    else:
        results['P3B2-6'] = {
            'pass': False,
            'note': 'Insufficient data for correlation analysis'
        }
    
    return results

def generate_results_report(df, lead_lag_matrix, corr_matrix, prediction_results):
    """Generate PHASE3B2_SECTOR_RESULTS.md report."""
    report = []
    report.append("# DS Phase 3B-2: Sector-Level Fisher Decomposition Results")
    report.append(f"**Generated:** {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    report.append(f"**Period:** {df['date'].min().date()} to {df['date'].max().date()}")
    report.append(f"**Total dates:** {len(df)}")
    report.append("")
    
    # Validation Gates
    report.append("## 1. Validation Gates")
    report.append("")
    
    # Gate A: Sector subgraph viability
    valid_sectors = []
    for sector in SECTORS:
        valid_pct = (1 - df[f'sector_{sector}'].isna().mean()) * 100
        valid_sectors.append((sector, valid_pct))
    
    report.append("### Gate A: Sector Subgraph Viability")
    report.append("| Sector | Valid % | Status |")
    report.append("|--------|---------|--------|")
    for sector, pct in valid_sectors:
        status = "✓" if pct >= 90 else "⚠" if pct >= 70 else "✗"
        report.append(f"| {sector} | {pct:.1f}% | {status} |")
    
    viable_count = sum(1 for _, pct in valid_sectors if pct >= 90)
    report.append(f"\n**Result:** {viable_count}/7 sectors viable (≥90% valid)")
    report.append(f"**Gate A:** {'PASS' if viable_count >= 5 else 'FAIL'}")
    report.append("")
    
    # Gate B: Sector SV2/SV1 variance
    report.append("### Gate B: Sector SV2/SV1 Variance")
    report.append("| Sector | Std Dev | Corr with Market | Status |")
    report.append("|--------|---------|------------------|--------|")
    for sector in SECTORS:
        std = df[f'sector_{sector}'].std()
        corr = df[f'sector_{sector}'].corr(df['market'])
        std_ok = std > 0.05
        corr_ok = abs(corr) < 0.95
        status = "✓" if std_ok and corr_ok else "⚠" if std_ok else "✗"
        report.append(f"| {sector} | {std:.3f} | {corr:.3f} | {status} |")
    
    report.append(f"\n**Gate B:** {'PASS' if all(std > 0.05 for std in [df[f'sector_{s}'].std() for s in SECTORS]) else 'FAIL'}")
    report.append("")
    
    # Prediction Results
    report.append("## 2. Pre-Registered Predictions")
    report.append("")
    
    for pred_id, result in prediction_results.items():
        report.append(f"### {pred_id}")
        
        if pred_id == 'P3B2-1':
            report.append("**KILL TEST: Financials leads Market in 2007-2008**")
            if result['lead_days'] is not None:
                report.append(f"- Financials leads by {result['lead_days']} days")
                report.append(f"- **Result:** {'PASS' if result['pass'] else 'FAIL'}")
                report.append(f"- **Kill condition:** {'TRIGGERED' if result['kill'] else 'NOT TRIGGERED'}")
            else:
                report.append("- No 2σ crossing detected")
                report.append("- **Result: FAIL**")
                report.append("- **Kill condition: TRIGGERED**")
        
        elif pred_id == 'P3B2-2':
            report.append("**Energy leads during 2014-2016 oil crash**")
            if result['lead_days'] is not None:
                report.append(f"- Energy leads by {result['lead_days']} days")
                report.append(f"- **Result:** {'PASS' if result['pass'] else 'FAIL'}")
            else:
                report.append("- No 2σ crossing detected")
                report.append("- **Result: FAIL**")
        
        elif pred_id == 'P3B2-3':
            report.append("**Tech leads during 2022 rate selloff**")
            if result['lead_days'] is not None:
                report.append(f"- Technology leads by {result['lead_days']} days")
                report.append(f"- **Result:** {'PASS' if result['pass'] else 'FAIL'}")
            else:
                report.append("- No 2σ crossing detected")
                report.append("- **Result: FAIL**")
        
        elif pred_id == 'P3B2-4':
            report.append("**Sector spread > 0 precedes SV2/SV1 peak**")
            report.append(f"- Successes: {result['successes']}/{result['total_crises']} crises")
            report.append(f"- **Result:** {'PASS' if result['pass'] else 'FAIL'}")
        
        elif pred_id == 'P3B2-5':
            report.append("**High Fisher-temperature assets cluster in leading sector**")
            report.append(f"- **Result:** {result['pass']}")
            if 'note' in result:
                report.append(f"- Note: {result['note']}")
        
        elif pred_id == 'P3B2-6':
            report.append("**Cross-sector correlation increases during crises**")
            if 'crisis_avg' in result:
                report.append(f"- Crisis average: {result['crisis_avg']:.3f}")
                report.append(f"- Calm average: {result['calm_avg']:.3f}")
                report.append(f"- **Result:** {'PASS' if result['pass'] else 'FAIL'}")
            else:
                report.append(f"- **Result:** {result['pass']}")
                if 'note' in result:
                    report.append(f"- Note: {result['note']}")
        
        report.append("")
    
    # Sector Statistics
    report.append("## 3. Sector Statistics")
    report.append("")
    report.append("| Sector | Mean SV2/SV1 | Std Dev | Min | Max | Valid % |")
    report.append("|--------|--------------|---------|-----|-----|---------|")
    for sector in SECTORS:
        series = df[f'sector_{sector}']
        valid_pct = (1 - series.isna().mean()) * 100
        report.append(f"| {sector} | {series.mean():.3f} | {series.std():.3f} | "
                     f"{series.min():.3f} | {series.max():.3f} | {valid_pct:.1f}% |")
    
    # Lead/Lag Summary
    report.append("\n## 4. Crisis Lead/Lag Summary")
    report.append("")
    report.append("Days for sector to exceed 2σ after crisis start:")
    report.append("")
    
    for crisis in CRISIS_PERIODS:
        if crisis in lead_lag_matrix.columns:
            report.append(f"### {crisis.replace('_', ' ')}")
            sorted_sectors = lead_lag_matrix[crisis].dropna().sort_values()
            for sector, days in sorted_sectors.items():
                report.append(f"- {sector}: {int(days)} days")
            report.append("")
    
    # Figures Generated
    report.append("## 5. Figures Generated")
    report.append("")
    figures = [
        "1. sector_sv2_sv1_timeseries.png",
        "2. financials_vs_market_2007_2008.png",
        "3. sector_lead_lag_heatmap.png", 
        "4. sector_spread_timeseries.png",
        "5. energy_vs_market_2014_2016.png",
        "6. tech_vs_market_2022.png",
        "7. top_fisher_temperature_assets.png",
        "8. sector_rank_comparison.png",
        "9. cross_sector_correlation_of_sv2sv1.png"
    ]
    for fig in figures:
        report.append(f"- {fig}")
    
    # Overall Assessment
    report.append("\n## 6. Overall Assessment")
    report.append("")
    
    kill_test = prediction_results.get('P3B2-1', {})
    if kill_test.get('kill', False):
        report.append("**❌ KILL CONDITION TRIGGERED**")
        report.append("Financials did not lead Market in 2007-2008.")
        report.append("Sector decomposition hypothesis fails for the clearest case.")
    else:
        report.append("**✅ KILL CONDITION NOT TRIGGERED**")
        report.append("Financials led Market in 2007-2008.")
        report.append("Sector decomposition hypothesis survives initial test.")
    
    # Write report
    report_path = RESULTS_DIR / "PHASE3B2_SECTOR_RESULTS.md"
    with open(report_path, 'w') as f:
        f.write('\n'.join(report))
    
    print(f"✓ Results report: {report_path}")
    return report_path

def main():
    """Main analysis pipeline."""
    print("DS Phase 3B-2: Sector-Level Fisher Analysis")
    print("=" * 50)
    
    # Load results (simulate - in practice would load from ObjectStore)
    print("Loading results...")
    # For now, create a dummy structure - in practice you'd load from:
    # json_path = "fisher_sector_results.json"
    # df, raw_data = load_results(json_path)
    
    print("Note: Analysis script ready. To use:")
    print("1. Run FisherSectorAlgorithm.py in QuantConnect")
    print("2. Download results from ObjectStore as JSON")
    print("3. Update json_path in main() function")
    print("4. Run: python fisher_sector_analysis.py")
    
    # The analysis functions are all defined and ready to use
    print("\nAnalysis functions available:")
    print("- figure1_sector_timeseries(df)")
    print("- figure2_financials_vs_market_2008(df)")
    print("- figure3_sector_lead_lag_heatmap(df)")
    print("- figure4_sector_spread_timeseries(df)")
    print("- figure5_energy_vs_market_2014_2016(df)")
    print("- figure6_tech_vs_market_2022(df)")
    print("- figure7_top_fisher_temperature_assets(df, raw_data)")
    print("- figure8_sector_rank_comparison(df)")
    print("- figure9_cross_sector_correlation(df)")
    print("- test_predictions(df, lead_lag_matrix, corr_matrix, crisis_corrs, calm_corr)")
    print("- generate_results_report(df, lead_lag_matrix, corr_matrix, prediction_results)")

if __name__ == "__main__":
    main()