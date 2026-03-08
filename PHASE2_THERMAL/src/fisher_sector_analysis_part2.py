def figure4_sector_spread_timeseries(df):
    """Figure 4: 30d-90d spread per sector."""
    fig, axes = plt.subplots(4, 2, figsize=(16, 12))
    axes = axes.flatten()
    
    for idx, sector in enumerate(SECTORS):
        ax = axes[idx]
        spread_series = df[f'spread_{sector}'].copy()
        
        # 4-week rolling average for smoothing
        spread_smooth = spread_series.rolling(window=4, min_periods=1).mean()
        
        ax.plot(df['date'], spread_smooth, color='purple', linewidth=1.5)
        ax.axhline(0, color='gray', alpha=0.5, linestyle='--')
        ax.fill_between(df['date'], 0, spread_smooth, 
                       where=spread_smooth > 0, color='red', alpha=0.3, label='Positive')
        ax.fill_between(df['date'], 0, spread_smooth, 
                       where=spread_smooth < 0, color='green', alpha=0.3, label='Negative')
        
        # Add crisis markers
        for crisis_date in CRISIS_DATES.keys():
            date = pd.to_datetime(crisis_date)
            if date >= df['date'].min() and date <= df['date'].max():
                ax.axvline(date, color='red', alpha=0.2, linestyle='--')
        
        ax.set_title(f'{sector} Spread (30d-90d)')
        ax.set_ylabel('Spread')
        ax.grid(True, alpha=0.3)
        
        if idx == 0:
            ax.legend()
    
    # Remove empty subplot
    axes[-1].axis('off')
    
    plt.suptitle('Sector Window Spread: SV2/SV1(30d) - SV2/SV1(90d)', fontsize=14)
    plt.tight_layout()
    plt.savefig(RESULTS_DIR / 'sector_spread_timeseries.png', dpi=300)
    plt.close()
    
    print("✓ Figure 4: sector_spread_timeseries.png")

def figure5_energy_vs_market_2014_2016(df):
    """Figure 5: Energy sector during oil price crash."""
    mask = (df['date'] >= '2014-01-01') & (df['date'] <= '2016-12-31')
    df_sub = df[mask].copy()
    
    # Compute z-scores
    df_sub['market_z'] = compute_z_scores(df_sub, 'market')
    df_sub['energy_z'] = compute_z_scores(df_sub, 'sector_Energy')
    
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(14, 8), sharex=True)
    
    # Panel A: Raw values
    ax1.plot(df_sub['date'], df_sub['market'], label='Market', color='black', linewidth=2)
    ax1.plot(df_sub['date'], df_sub['sector_Energy'], label='Energy', 
             color='orange', linewidth=2)
    ax1.axvline(pd.to_datetime('2015-01-01'), color='red', alpha=0.5, 
                linestyle='--', label='Oil Crash')
    ax1.set_ylabel('SV2/SV1')
    ax1.set_title('Energy vs Market SV2/SV1 (2014-2016 Oil Price Crash)')
    ax1.legend()
    ax1.grid(True, alpha=0.3)
    
    # Panel B: Z-scores
    ax2.plot(df_sub['date'], df_sub['market_z'], label='Market Z', color='black', linewidth=2)
    ax2.plot(df_sub['date'], df_sub['energy_z'], label='Energy Z', 
             color='orange', linewidth=2)
    ax2.axhline(2, color='red', alpha=0.5, linestyle='--', label='2σ threshold')
    ax2.axhline(-2, color='red', alpha=0.5, linestyle='--')
    ax2.axvline(pd.to_datetime('2015-01-01'), color='red', alpha=0.5, linestyle='--')
    ax2.set_ylabel('Z-score')
    ax2.set_xlabel('Date')
    ax2.legend()
    ax2.grid(True, alpha=0.3)
    
    # Find crossing dates
    energy_up, _ = find_crossing_dates(df_sub['energy_z'], threshold=2)
    mkt_up, _ = find_crossing_dates(df_sub['market_z'], threshold=2)
    
    if not energy_up.empty and not mkt_up.empty:
        energy_first = energy_up['date'].min()
        mkt_first = mkt_up['date'].min()
        lead_days = (mkt_first - energy_first).days
        
        ax2.text(0.02, 0.95, f'Energy leads by {lead_days} days', 
                transform=ax2.transAxes, fontsize=10, 
                bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))
    
    plt.tight_layout()
    plt.savefig(RESULTS_DIR / 'energy_vs_market_2014_2016.png', dpi=300)
    plt.close()
    
    print("✓ Figure 5: energy_vs_market_2014_2016.png")

def figure6_tech_vs_market_2022(df):
    """Figure 6: Tech sector during 2022 rate hike selloff."""
    mask = (df['date'] >= '2021-10-01') & (df['date'] <= '2023-03-31')
    df_sub = df[mask].copy()
    
    # Compute z-scores
    df_sub['market_z'] = compute_z_scores(df_sub, 'market')
    df_sub['tech_z'] = compute_z_scores(df_sub, 'sector_Technology')
    
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(14, 8), sharex=True)
    
    # Panel A: Raw values
    ax1.plot(df_sub['date'], df_sub['market'], label='Market', color='black', linewidth=2)
    ax1.plot(df_sub['date'], df_sub['sector_Technology'], label='Technology', 
             color='blue', linewidth=2)
    ax1.axvline(pd.to_datetime('2022-06-01'), color='red', alpha=0.5, 
                linestyle='--', label='Rate Hike Peak')
    ax1.set_ylabel('SV2/SV1')
    ax1.set_title('Technology vs Market SV2/SV1 (2022 Rate Hike Selloff)')
    ax1.legend()
    ax1.grid(True, alpha=0.3)
    
    # Panel B: Z-scores
    ax2.plot(df_sub['date'], df_sub['market_z'], label='Market Z', color='black', linewidth=2)
    ax2.plot(df_sub['date'], df_sub['tech_z'], label='Technology Z', 
             color='blue', linewidth=2)
    ax2.axhline(2, color='red', alpha=0.5, linestyle='--', label='2σ threshold')
    ax2.axhline(-2, color='red', alpha=0.5, linestyle='--')
    ax2.axvline(pd.to_datetime('2022-06-01'), color='red', alpha=0.5, linestyle='--')
    ax2.set_ylabel('Z-score')
    ax2.set_xlabel('Date')
    ax2.legend()
    ax2.grid(True, alpha=0.3)
    
    # Find crossing dates
    tech_up, _ = find_crossing_dates(df_sub['tech_z'], threshold=2)
    mkt_up, _ = find_crossing_dates(df_sub['market_z'], threshold=2)
    
    if not tech_up.empty and not mkt_up.empty:
        tech_first = tech_up['date'].min()
        mkt_first = mkt_up['date'].min()
        lead_days = (mkt_first - tech_first).days
        
        ax2.text(0.02, 0.95, f'Technology leads by {lead_days} days', 
                transform=ax2.transAxes, fontsize=10, 
                bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))
    
    plt.tight_layout()
    plt.savefig(RESULTS_DIR / 'tech_vs_market_2022.png', dpi=300)
    plt.close()
    
    print("✓ Figure 6: tech_vs_market_2022.png")

def figure7_top_fisher_temperature_assets(df, raw_data):
    """Figure 7: Top-10 highest Fisher-temperature assets per week."""
    # Extract top assets each week
    top_assets_by_week = []
    
    for i, date in enumerate(df['date']):
        asset_temps = df.at[i, 'asset_temps']
        if asset_temps:
            # Sort by temperature (descending)
            sorted_temps = sorted(asset_temps.items(), key=lambda x: x[1], reverse=True)
            top_10 = sorted_temps[:10]
            top_assets_by_week.append({
                'date': date,
                'top_assets': top_10,
                'top_sector': None
            })
    
    # Create time series of top assets
    all_assets = set()
    for week in top_assets_by_week:
        for ticker, _ in week['top_assets']:
            all_assets.add(ticker)
    
    # Map tickers to sectors
    ticker_to_sector = {}
    for sector, tickers in raw_data.get('sector_map', {}).items():
        for ticker in tickers:
            ticker_to_sector[ticker] = sector
    
    # Create presence matrix
    presence_matrix = pd.DataFrame(index=df['date'], columns=list(all_assets))
    for week in top_assets_by_week:
        date = week['date']
        for ticker, temp in week['top_assets']:
            presence_matrix.at[date, ticker] = temp
    
    # Plot top 20 most frequent assets
    freq_counts = presence_matrix.count().sort_values(ascending=False)
    top_20_assets = freq_counts.head(20).index
    
    fig, ax = plt.subplots(figsize=(16, 10))
    
    # Create heatmap
    heatmap_data = presence_matrix[top_20_assets].T
    
    # Sort by sector
    heatmap_data['sector'] = heatmap_data.index.map(ticker_to_sector)
    heatmap_data = heatmap_data.sort_values('sector', ascending=False)
    heatmap_data = heatmap_data.drop('sector', axis=1)
    
    im = ax.imshow(heatmap_data.values, cmap='YlOrRd', aspect='auto', 
                   interpolation='nearest')
    
    # Add sector color coding
    sector_colors = {'Technology': 'blue', 'Financials': 'red', 
                     'HealthCare': 'green', 'Consumer': 'orange',
                     'Industrials': 'purple', 'Energy': 'brown',
                     'Utilities': 'gray'}
    
    for i, ticker in enumerate(heatmap_data.index):
        sector = ticker_to_sector.get(ticker, 'Unknown')
        color = sector_colors.get(sector, 'black')
        ax.text(-0.5, i, sector[:3], ha='right', va='center', 
               color=color, fontsize=8, fontweight='bold')
    
    ax.set_yticks(range(len(heatmap_data)))
    ax.set_yticklabels(heatmap_data.index)
    ax.set_xlabel('Date')
    ax.set_title('Top Fisher-Temperature Assets by Week (Color = Temperature)')
    
    # Add crisis markers
    xticks = ax.get_xticks()
    date_labels = [df['date'].iloc[int(x)].strftime('%Y-%m') 
                   if x < len(df) else '' for x in xticks]
    ax.set_xticklabels(date_labels, rotation=45, ha='right')
    
    plt.colorbar(im, ax=ax, label='Fisher Temperature')
    plt.tight_layout()
    plt.savefig(RESULTS_DIR / 'top_fisher_temperature_assets.png', dpi=300)
    plt.close()
    
    print("✓ Figure 7: top_fisher_temperature_assets.png")