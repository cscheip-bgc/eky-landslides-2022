"""
Precipitation analysis and plotting module.

This module processes precipitation gauge data from eastern Kentucky
and generates rolling sum precipitation plots.

Data source: https://www.weather.gov/wrh/Climate?wfo=jkl
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import sys

# Project config
sys.path.append(str(Path(__file__).resolve().parent.parent))
import config as cfg


def load_precipitation_data(data_path=None):
    """
    Load precipitation data from CSV file.

    Parameters
    ----------
    data_path : str or Path, optional
        Path to the precipitation data CSV file (defaults to
        ``config.precip_gauge_csv``)

    Returns
    -------
    pd.DataFrame
        DataFrame containing precipitation data from multiple gauges
    """
    if data_path is None:
        data_path = cfg.precip_gauge_csv

    # Load data with proper date parsing
    allprecip = pd.read_csv(data_path, parse_dates=['Date'], index_col='Date')
    
    # Ensure all columns are numeric (convert any string values to numeric, coercing errors to NaN)
    for col in allprecip.columns:
        allprecip[col] = pd.to_numeric(allprecip[col], errors='coerce')
    
    # Fill any NaN values with 0 (reasonable for precipitation data)
    allprecip = allprecip.fillna(0)
    
    return allprecip


def calculate_rolling_sums(allprecip, window=30):
    """
    Calculate rolling sum statistics across precipitation gauges.
    
    Parameters
    ----------
    allprecip : pd.DataFrame
        DataFrame containing precipitation data
    window : int, optional
        Rolling window size in days (default: 30)
        
    Returns
    -------
    pd.DataFrame
        DataFrame with rolling sum calculations
    """
    calcs = allprecip.copy()
    calcs['Average'] = allprecip.mean(axis=1)
    calcs['Med'] = allprecip.median(axis=1)
    calcs['Max'] = allprecip.max(axis=1)
    calcs['Min'] = allprecip.min(axis=1)
    
    df = calcs.copy()
    df['Daily_Med(mm)'] = df['Med'].rolling(window).sum()
    df['b_roll'] = df['Buckhorn'].rolling(window).sum()
    df['j_roll'] = df['Jackson'].rolling(window).sum()
    df['o_roll'] = df['Oneida'].rolling(window).sum()
    df['w_roll'] = df['Whitesburg'].rolling(window).sum()
    
    selected_columns = ['b_roll', 'j_roll', 'o_roll', 'w_roll']
    df['MaxRollingSum'] = df[selected_columns].max(axis=1)
    df['MinRollingSum'] = df[selected_columns].min(axis=1)
    
    return df


def plot_full_span_rolling_sum(df, output_path=None):
    """
    Generate full time span 30-day rolling sum precipitation plot.

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame with rolling sum calculations
    output_path : str or Path, optional
        Output path for the figure (defaults to
        ``config.figure_path / 'Month_AverageFullSpan.png'``)
    """
    if output_path is None:
        output_path = cfg.figure_path / 'Month_AverageFullSpan.png'

    plt.figure(figsize=(8.5, 4.8))
    
    plt.plot(df.index, df['MaxRollingSum'], label='Upper', 
             linestyle='--', linewidth=1, color='lightgray')
    plt.plot(df.index, df['MinRollingSum'], label='Lower', 
             linestyle='--', linewidth=1, color='lightgray')
    plt.fill_between(df.index, df['MaxRollingSum'], df['MinRollingSum'], 
                     color='lightsteelblue', alpha=0.5)
    plt.plot(df.index, df['Daily_Med(mm)'], linestyle='-', color='midnightblue')
    
    plt.ylabel('30-day Rolling Sum Precipitation (mm)', fontname='Arial', size=12)
    plt.xticks(fontname='Arial', size=12)
    plt.yticks(fontname='Arial', size=12)
    plt.gca().tick_params(axis='x', labeltop=False, top=True)
    plt.gca().tick_params(axis='y', which='both', labelright=False, right=True)
    plt.gca().tick_params(axis='both', direction='out')
    plt.gca().grid(True, alpha=0.3)
    plt.tight_layout()
    
    # Ensure output directory exists
    Path(output_path).parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(output_path, format='png', dpi=600, transparent=True)
    plt.show()
    print(f"Saved full span plot to {output_path}")


def plot_fourday_inset(allprecip, year=2022, months=[6, 7, 8],
                        output_path=None):
    """
    Generate 4-day rolling sum precipitation plot for specified period.

    Parameters
    ----------
    allprecip : pd.DataFrame
        DataFrame containing precipitation data
    year : int, optional
        Year to filter (default: 2022)
    months : list, optional
        List of months to include (default: [6, 7, 8] for June-August)
    output_path : str or Path, optional
        Output path for the figure (defaults to
        ``config.figure_path / 'FourDay_GaugeInset.png'``)
    """
    if output_path is None:
        output_path = cfg.figure_path / 'FourDay_GaugeInset.png'

    df_filtered = allprecip[allprecip.index.year == year]
    df_aroundJuly = df_filtered[df_filtered.index.month.isin(months)]
    df_July = df_aroundJuly.copy()
    
    df_July['RollFour_b'] = df_July['Buckhorn'].rolling(4).sum()
    df_July['RollFour_j'] = df_July['Jackson'].rolling(4).sum()
    df_July['RollFour_o'] = df_July['Oneida'].rolling(4).sum()
    df_July['RollFour_w'] = df_July['Whitesburg'].rolling(4).sum()
    
    plt.figure(figsize=(4.25, 2.25))
    
    plt.plot(df_July.index, df_July['RollFour_b'], linestyle='-', 
             color='red', alpha=0.5)
    plt.plot(df_July.index, df_July['RollFour_j'], linestyle='-', 
             color='yellow', alpha=0.5)
    plt.plot(df_July.index, df_July['RollFour_o'], linestyle='-', 
             color='purple', alpha=0.5)
    plt.plot(df_July.index, df_July['RollFour_w'], linestyle='-', 
             color='orange', alpha=0.5)
    
    plt.xticks(fontname='Arial', size=9)
    plt.yticks(fontname='Arial', size=10)
    plt.ylabel('Four-day Rolling\nSum Precipitation (mm)', fontname='Arial', size=10)
    plt.gca().tick_params(axis='x', labeltop=False, top=True)
    plt.gca().tick_params(axis='y', which='both', labelright=False, right=True)
    plt.gca().tick_params(axis='both', direction='out')
    plt.xticks(rotation=-35, ha='left')
    plt.gca().grid(True, alpha=0.3)
    plt.tight_layout()
    
    # Ensure output directory exists
    Path(output_path).parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(output_path, format='png', dpi=600, transparent=True)
    plt.show()
    print(f"Saved 4-day inset plot to {output_path}")


def calculate_daily_intensity_stats(allprecip, percentiles=(90, 95, 99),
                                    wet_threshold=1.0):
    """
    Summarize 1-day (daily) rainfall totals across the gauge network.

    Daily totals are the shortest duration resolvable from these gauge
    records, and are used to test whether any short-duration, high-intensity
    storm outside the July 2022 event approached its magnitude.

    Percentiles are reported two ways because the record is ~58% dry days:
    over all days (dry days included, which pulls low percentiles toward
    zero) and over wet days only (>= ``wet_threshold``), which is the
    conventional climatological definition (e.g. ETCCDI R95p/R99p).

    Parameters
    ----------
    allprecip : pd.DataFrame
        Daily precipitation totals, one column per gauge (mm)
    percentiles : tuple of int, optional
        Percentiles to report (default: 90th, 95th, 99th)
    wet_threshold : float, optional
        Minimum daily total (mm) for a day to count as wet (default: 1.0)

    Returns
    -------
    pd.DataFrame
        One row per gauge, plus a pooled row over all gauge-days and a
        "max across gauges" row representing the wettest reading anywhere
        in the network on each day
    """
    # Wettest reading anywhere in the network on each day
    max_anywhere = allprecip.max(axis=1)

    series = {gauge: allprecip[gauge] for gauge in allprecip.columns}
    # Pool every gauge-day observation into a single distribution
    series['Pooled (all gauge-days)'] = allprecip.stack()
    series['Max across gauges (anywhere)'] = max_anywhere

    rows = []
    for name, s in series.items():
        wet = s[s >= wet_threshold]
        row = {
            'Series': name,
            'N_obs': int(s.size),
            f'N_wet_days_ge_{wet_threshold:g}mm': int(wet.size),
            'Pct_wet_days': round(100 * wet.size / s.size, 1),
        }
        for p in percentiles:
            row[f'P{p}_all_days_mm'] = round(float(s.quantile(p / 100)), 2)
        for p in percentiles:
            row[f'P{p}_wet_days_mm'] = round(float(wet.quantile(p / 100)), 2)
        row['Max_mm'] = round(float(s.max()), 2)

        # Date (and gauge, where meaningful) of the record daily total
        if name == 'Pooled (all gauge-days)':
            date, gauge = s.idxmax()
        elif name == 'Max across gauges (anywhere)':
            date = s.idxmax()
            gauge = allprecip.loc[date].idxmax()
        else:
            date = s.idxmax()
            gauge = name
        row['Max_date'] = pd.Timestamp(date).date().isoformat()
        row['Max_gauge'] = gauge
        rows.append(row)

    return pd.DataFrame(rows)


def rank_daily_events(allprecip, top_n=20):
    """
    Rank individual days by the highest single-gauge total recorded that day.

    Answers the reviewer-facing question directly: if the July 2022 event
    dominates this ranking, no other short-duration storm in the record
    delivered comparable daily rainfall.

    Parameters
    ----------
    allprecip : pd.DataFrame
        Daily precipitation totals, one column per gauge (mm)
    top_n : int, optional
        Number of ranked days to return (default: 20)

    Returns
    -------
    pd.DataFrame
        Ranked days with each gauge's total, the network maximum, the gauge
        that recorded it, and the network mean for that day
    """
    ranked = allprecip.copy()
    ranked['Max_anywhere_mm'] = allprecip.max(axis=1)
    ranked['Max_gauge'] = allprecip.idxmax(axis=1)
    ranked['Network_mean_mm'] = allprecip.mean(axis=1)

    ranked = ranked.sort_values('Max_anywhere_mm', ascending=False).head(top_n)
    ranked = ranked.round(2)
    ranked.insert(0, 'Rank', range(1, len(ranked) + 1))
    ranked.index = ranked.index.date
    ranked.index.name = 'Date'

    return ranked


def count_threshold_exceedances(allprecip, thresholds=(75, 100, 150, 200),
                                wet_threshold=1.0):
    """
    Count days meeting or exceeding fixed daily-rainfall thresholds.

    All thresholds of interest sit far above ``wet_threshold``, so every
    qualifying day is by definition a wet day; the wet-day set matters only
    as the denominator for the reported percentages.

    Counts are given two ways: pooled gauge-days (each gauge counted
    separately, so one widespread storm can contribute up to four) and
    network-days (a calendar day counts once if any gauge qualifies).

    Parameters
    ----------
    allprecip : pd.DataFrame
        Daily precipitation totals, one column per gauge (mm)
    thresholds : tuple of float, optional
        Daily totals (mm) to test, using >= (default: 75, 100, 150, 200)
    wet_threshold : float, optional
        Minimum daily total (mm) for a day to count as wet (default: 1.0)

    Returns
    -------
    pd.DataFrame
        One row per threshold with per-gauge counts, pooled gauge-day and
        network-day counts, percentages of the respective wet-day totals,
        and network-day frequency per year
    """
    max_anywhere = allprecip.max(axis=1)
    n_years = len(allprecip) / 365.25

    n_wet_pooled = int((allprecip.stack() >= wet_threshold).sum())
    n_wet_network = int((max_anywhere >= wet_threshold).sum())

    rows = []
    for t in thresholds:
        per_gauge = (allprecip >= t).sum()
        n_pooled = int(per_gauge.sum())
        n_network = int((max_anywhere >= t).sum())

        row = {'Threshold_mm': t}
        row.update({f'{g}_days': int(per_gauge[g]) for g in allprecip.columns})
        row['Pooled_gauge_days'] = n_pooled
        row['Pct_of_wet_gauge_days'] = round(100 * n_pooled / n_wet_pooled, 3)
        row['Network_days_any_gauge'] = n_network
        row['Pct_of_wet_network_days'] = round(100 * n_network / n_wet_network, 3)
        row['Network_days_per_year'] = round(n_network / n_years, 2)
        rows.append(row)

    return pd.DataFrame(rows)


def report_daily_intensity(allprecip, output_dir=None, top_n=20):
    """
    Compute, print, and save the 1-day rainfall intensity tables.

    Parameters
    ----------
    allprecip : pd.DataFrame
        Daily precipitation totals, one column per gauge (mm)
    output_dir : str or Path, optional
        Directory for the CSV tables (defaults to ``config.figure_path``)
    top_n : int, optional
        Number of ranked days to save (default: 20)

    Returns
    -------
    tuple of pd.DataFrame
        (summary table, ranked daily events table)
    """
    if output_dir is None:
        output_dir = cfg.figure_path
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    summary = calculate_daily_intensity_stats(allprecip)
    ranked = rank_daily_events(allprecip, top_n=top_n)
    exceedance = count_threshold_exceedances(allprecip)

    with pd.option_context('display.width', 200,
                           'display.max_columns', None):
        print('\n1-day rainfall totals - distribution summary (mm)')
        print(summary.to_string(index=False))
        print(f'\nTop {top_n} days by highest single-gauge total (mm)')
        print(ranked.to_string())
        print('\nDaily threshold exceedances (>= threshold)')
        print(exceedance.to_string(index=False))

    summary_path = output_dir / 'precip_daily_percentiles.csv'
    ranked_path = output_dir / 'precip_daily_top_events.csv'
    exceedance_path = output_dir / 'precip_daily_exceedances.csv'
    summary.to_csv(summary_path, index=False)
    ranked.to_csv(ranked_path)
    exceedance.to_csv(exceedance_path, index=False)
    print(f'\nSaved daily percentile summary to {summary_path}')
    print(f'Saved ranked daily events to {ranked_path}')
    print(f'Saved threshold exceedances to {exceedance_path}')

    return summary, ranked, exceedance


def main():
    """
    Main function to run all precipitation analysis and plotting.
    """
    print("Loading precipitation data...")
    allprecip = load_precipitation_data()
    
    print("Calculating 30-day rolling sums...")
    df = calculate_rolling_sums(allprecip, window=30)
    
    print("Generating full span plot...")
    plot_full_span_rolling_sum(df)
    
    print("Generating 4-day rolling sum inset plot...")
    plot_fourday_inset(allprecip, year=2022, months=[6, 7, 8])

    print("Summarizing 1-day rainfall totals...")
    report_daily_intensity(allprecip)

    print("\nPrecipitation analysis complete!")


if __name__ == '__main__':
    main()

