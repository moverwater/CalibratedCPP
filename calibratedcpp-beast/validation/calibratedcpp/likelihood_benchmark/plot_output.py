#!/usr/bin/env python3
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path

# 1. Setup Relative Paths
SCRIPT_DIR = Path(__file__).parent
BENCH_CSV = SCRIPT_DIR / "benchmark_results.csv"
COMP_CSV = SCRIPT_DIR / "validation_results.csv"
README_PATH = SCRIPT_DIR / "README.md"
PLOT_PATH = SCRIPT_DIR / "combined_benchmark_validation.png"

# Style configuration
COLORS = {"Calibrated CPP": "black", "Heled & Drummond": "red"}

def update_readme(df):
    if not README_PATH.exists():
        return

    content = README_PATH.read_text(encoding='utf-8')

    # Prepare Clean Table
    display_df = df[['Benchmark', 'Param: nCalibrations', 'Score']].copy()
    display_df['Benchmark'] = display_df['Benchmark'].apply(
        lambda x: "Calibrated CPP" if "measureCPP" in x else "Heled & Drummond"
    )
    display_df.columns = ["Model", "Number of Calibrations", "Time (μs)"]
    md_table = display_df.to_markdown(index=False, tablefmt="github", floatfmt=".2f")

    # Slice and Stitch
    start_anchor = "<!-- start table -->"
    end_anchor = "<!-- end table -->"
    start_idx = content.find(start_anchor)
    end_idx = content.find(end_anchor)

    if start_idx != -1 and end_idx != -1:
        prefix = content[:start_idx + len(start_anchor)]
        suffix = content[end_idx:]
        new_content = f"{prefix}\n\n{md_table}\n\n{suffix}"
        README_PATH.write_text(new_content, encoding='utf-8')
        print("README.md updated.")

def main():
    if not BENCH_CSV.exists() or not COMP_CSV.exists():
        print("Error: One or both CSV files are missing.")
        return

    # Load Data
    df_bench = pd.read_csv(BENCH_CSV)
    df_comp = pd.read_csv(COMP_CSV)

    # Setup Plot (1 row, 2 columns)
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4.5))

    # --- Plot 1: Benchmark (Time) ---
    for label, group in df_bench.groupby('Benchmark'):
        name = "Calibrated CPP" if "measureCPP" in label else "Heled & Drummond"
        ax1.plot(group['Param: nCalibrations'], group['Score'], 
                 label=name, color=COLORS[name], marker='o', linewidth=1.5)
    
    ax1.set_yscale('log')
    ax1.set_xlabel('Number of calibrations (k)')
    ax1.set_ylabel('Time (μs)')
    ax1.set_title('Benchmark Performance')
    ax1.grid(True, alpha=0.3)

    # --- Plot 2: Likelihood (Validation) ---
    # Assuming CSV columns: BirthRate, HeledAndDrummond_LogLikelihood, CPP_LogLikelihood
    ax2.plot(df_comp['BirthRate'], df_comp['HeledAndDrummond_LogLikelihood'], 
             color=COLORS["Heled & Drummond"], label="Heled & Drummond", linewidth=1.5)
    ax2.scatter(df_comp['BirthRate'], df_comp['CPP_LogLikelihood'], 
                color=COLORS["Calibrated CPP"], label="Calibrated CPP", s=20, zorder=3)
    
    ax2.set_xlabel('Birth-rate')
    ax2.set_ylabel('Log-likelihood')
    ax2.set_title('Likelihood Validation')
    ax2.grid(True, alpha=0.3)

    # --- Global Formatting (Shared Legend) ---
    handles, labels = ax1.get_legend_handles_labels()
    fig.legend(handles, labels, loc='lower center', ncol=2, frameon=True)
    
    # Adjust layout to make room for the legend at the bottom
    plt.tight_layout(rect=[0, 0.08, 1, 1])
    
    plt.savefig(PLOT_PATH, dpi=150)
    print(f"Combined plot saved to {PLOT_PATH}")

    update_readme(df_bench)

if __name__ == "__main__":
    main()