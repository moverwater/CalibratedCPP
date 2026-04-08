import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os
from pathlib import Path

# --- 1. SETUP ---
# Set working directory to script location
SCRIPT_DIR = Path(__file__).parent
os.chdir(SCRIPT_DIR)

# --- 2. CONFIGURATION ---
target_columns = [
    'mrca.age(TaxonSet1)', 
    'mrca.age(TaxonSet2)', 
    'mrca.age(TaxonSet3)'
]

# Map specific files to their logical categories
# Format: (File Path, Analysis Type, Condition Type)
file_map = [
    # PRIOR FILES
    ('cats-0_sample_from_prior.log', 'Prior', 'Conditioned'),
    ('cats-0_not_conditioned_sample_from_prior.log', 'Prior', 'Not Conditioned'),
    
    # POSTERIOR FILES
    ('cats-0.log', 'Posterior', 'Conditioned'),
    ('cats-0_not_conditioned.log', 'Posterior', 'Not Conditioned')
]

burnin_fraction = 0.1

# --- 3. HELPER FUNCTION ---
def load_all_data(file_map):
    all_data = []
    
    for filename, analysis_type, condition_type in file_map:
        if os.path.exists(filename):
            try:
                # Read file
                df = pd.read_csv(filename, sep='\t', comment='#')
                
                # Check for columns
                cols = [c for c in target_columns if c in df.columns]
                if cols:
                    # Burn-in
                    drop_n = int(len(df) * burnin_fraction)
                    df = df.iloc[drop_n:].copy()
                    
                    # Melt to long format immediately for easier plotting
                    # Creates rows like: [Prior, Conditioned, TaxonSet1, 10.5]
                    melted = df[cols].melt(var_name='TaxonSet', value_name='Age')
                    melted['Analysis'] = analysis_type
                    melted['Condition'] = condition_type
                    
                    # Clean TaxonSet names
                    melted['TaxonSet'] = melted['TaxonSet'].str.replace('mrca.age(', '', regex=False).str.replace(')', '', regex=False)
                    
                    all_data.append(melted)
                    print(f"Loaded: {filename}")
                else:
                    print(f"Warning: Targets not found in {filename}")
            except Exception as e:
                print(f"Error reading {filename}: {e}")
        else:
            print(f"File not found: {filename}")

    if all_data:
        return pd.concat(all_data, ignore_index=True)
    return pd.DataFrame()

# --- 4. LOAD DATA ---
df = load_all_data(file_map)

if df.empty:
    print("No data loaded.")
    exit()

# --- 5. PLOTTING ---
# We create a figure with 2 Rows (Prior, Posterior) and 3 Columns (The Taxon Sets)
fig, axes = plt.subplots(2, 3, figsize=(18, 10), constrained_layout=True)

# Define unique Taxon Sets for iteration
taxon_sets = [c.replace('mrca.age(', '').replace(')', '') for c in target_columns]
analysis_types = ['Prior', 'Posterior']
colors = {'Conditioned': 'tab:blue', 'Not Conditioned': 'tab:orange'}

for row_idx, analysis in enumerate(analysis_types):
    # Filter for just Prior or just Posterior
    subset_analysis = df[df['Analysis'] == analysis]
    
    for col_idx, taxon in enumerate(taxon_sets):
        ax = axes[row_idx, col_idx]
        
        # Filter for specific TaxonSet
        data_to_plot = subset_analysis[subset_analysis['TaxonSet'] == taxon]
        
        if not data_to_plot.empty:
            # sns.kdeplot(
            #     data=data_to_plot,
            #     x='Age',
            #     hue='Condition',
            #     fill=True,
            #     palette=colors,
            #     alpha=0.4,
            #     linewidth=2,
            #     ax=ax,
            #     common_norm=False # Important: treats them as separate distributions
            # )
            sns.histplot(
                data=data_to_plot,
                x='Age',
                hue='Condition',
                fill=True,
                palette=colors,
                alpha=0.4,
                ax=ax,
                common_norm=False, # Normalize each group independently
                element="step",    # Use 'step' to draw outlines (cleaner for overlaps)
                stat="density",    # Normalize y-axis (crucial if file lengths differ)
                linewidth=0.5
            )
            
            # Labeling
            if row_idx == 0:
                ax.set_title(f"{taxon} (Prior)", fontsize=14, fontweight='bold')
            else:
                ax.set_title(f"{taxon} (Posterior)", fontsize=14, fontweight='bold')
                
            ax.set_xlabel("Age" if row_idx == 1 else "")
            
            # Handle Legend: Only show it on the very last plot to avoid clutter
            if (row_idx == 0 and col_idx == 2):
                sns.move_legend(ax, "upper right", title=None)
            else:
                if ax.get_legend(): ax.get_legend().remove()
        else:
            ax.set_visible(False)
plt.savefig("cats_comparison_plot.png", dpi=300, bbox_inches='tight')