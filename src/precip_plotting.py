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


def load_precipitation_data(data_path='data/eky-gauge-data.csv'):
    """
    Load precipitation data from CSV file.
    
    Parameters
    ----------
    data_path : str, optional
        Path to the precipitation data CSV file
        
    Returns
    -------
    pd.DataFrame
        DataFrame containing precipitation data from multiple gauges
    """
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


def plot_full_span_rolling_sum(df, output_path='outputs/figures/Month_AverageFullSpan.png'):
    """
    Generate full time span 30-day rolling sum precipitation plot.
    
    Parameters
    ----------
    df : pd.DataFrame
        DataFrame with rolling sum calculations
    output_path : str, optional
        Output path for the figure
    """
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
                        output_path='outputs/figures/FourDay_GaugeInset.png'):
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
    output_path : str, optional
        Output path for the figure
    """
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
    
    print("Precipitation analysis complete!")


if __name__ == '__main__':
    main()

