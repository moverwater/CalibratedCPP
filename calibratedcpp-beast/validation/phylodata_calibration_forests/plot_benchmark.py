#!/usr/bin/env python3
"""
Plot benchmark results from CalibratedCPP constraint tree benchmarks.

Usage:
    python plot_benchmark.py
"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path

# Use script directory for relative paths
SCRIPT_DIR = Path(__file__).parent
CSV_PATH = SCRIPT_DIR / "benchmark_results.csv"
PLOT_PATH = SCRIPT_DIR / "benchmark_analysis_plot.png"


def main():
    df = pd.read_csv(CSV_PATH)
    df = df[df['error'].isna() | (df['error'] == '')]

    print(f"Successful benchmarks: {len(df)}")

    # Fit multivariate: time ~ complexity + taxa
    X = np.column_stack([df['complexity'].values, df['taxa'].values, np.ones(len(df))])
    y = df['p50_us'].values
    coeffs, _, _, _ = np.linalg.lstsq(X, y, rcond=None)
    a, b, c = coeffs

    y_pred = a * df['complexity'].values + b * df['taxa'].values + c
    ss_res = np.sum((y - y_pred)**2)
    ss_tot = np.sum((y - np.mean(y))**2)
    r2 = 1 - ss_res / ss_tot

    print(f"Regression: time = {a:.5f}*complexity + {b:.3f}*taxa + {c:.1f}")
    print(f"R² = {r2:.4f}")

    # Create figure
    fig, axes = plt.subplots(1, 3, figsize=(16, 5))

    # Plot 1: Time vs Complexity (log-log)
    mask_nonzero = df['complexity'] > 0
    axes[0].scatter(df.loc[mask_nonzero, 'complexity'], df.loc[mask_nonzero, 'p50_us'],
                    alpha=0.6, s=40, label=f'complexity>0 (n={mask_nonzero.sum()})')
    axes[0].scatter(df.loc[~mask_nonzero, 'complexity'] + 0.5, df.loc[~mask_nonzero, 'p50_us'],
                    alpha=0.6, s=40, c='orange', label=f'complexity=0 (n={(~mask_nonzero).sum()})')
    axes[0].set_xscale('log')
    axes[0].set_yscale('log')
    axes[0].set_xlabel('Complexity Score')
    axes[0].set_ylabel('Median Time (μs)')
    axes[0].set_title('Time vs Complexity (Log-Log)')
    axes[0].legend()
    axes[0].grid(True, alpha=0.3)

    # Plot 2: Time vs Taxa
    axes[1].scatter(df['taxa'], df['p50_us'], alpha=0.6, s=40)
    axes[1].set_yscale('log')
    axes[1].set_xlabel('Number of Taxa')
    axes[1].set_ylabel('Median Time (μs)')
    axes[1].set_title('Time vs Taxa Count')
    axes[1].grid(True, alpha=0.3)

    # Plot 3: Actual vs Predicted
    axes[2].scatter(y_pred, y, alpha=0.6, s=40)
    axes[2].plot([0, max(y)], [0, max(y)], 'r--', label='Perfect fit', linewidth=2)
    axes[2].set_xscale('log')
    axes[2].set_yscale('log')
    axes[2].set_xlabel('Predicted Time (μs)')
    axes[2].set_ylabel('Actual Time (μs)')
    axes[2].set_title(f'Multivariate Fit (R²={r2:.4f})')
    axes[2].legend()
    axes[2].grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(PLOT_PATH, dpi=150)
    print(f"\nSaved: {PLOT_PATH}")


if __name__ == "__main__":
    main()
