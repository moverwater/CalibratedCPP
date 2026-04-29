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

# Map your 6 specific files to their logical categories
# Format: (File Path, Analysis Type, Condition Type)
file_map = [
    # PRIOR FILES
    ('cats-sample_from_prior.log', 'Prior', 'Conditioned'),
    ('cats-sample_from_prior_not_conditioned.log', 'Prior', 'Not Conditioned'),
    
    # HALF FILES
    ('cats-half.log', 'Half Alignment', 'Conditioned'),
    ('cats-half_not_conditioned.log', 'Half Alignment', 'Not Conditioned'),
    
    # FULL FILES 
    ('cats-full.log', 'Full Alignment', 'Conditioned'),
    ('cats-full_not_conditioned.log', 'Full Alignment', 'Not Conditioned')
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
# 1. Globally increase the font sizes for axes, ticks, and legend
plt.rcParams.update({
    'axes.labelsize': 18,   # Size of the "Age" x-axis label
    'xtick.labelsize': 14,  # Size of the numbers on the x-axis
    'ytick.labelsize': 14,  # Size of the numbers on the y-axis
    'legend.fontsize': 16   # Size of the text in the legend
})

fig, axes = plt.subplots(3, 3, figsize=(18, 14), constrained_layout=True)

taxon_sets = [c.replace('mrca.age(', '').replace(')', '') for c in target_columns]
analysis_types = ['Prior', 'Half Alignment', 'Full Alignment']
colors = {'Conditioned': 'tab:blue', 'Not Conditioned': 'tab:orange'}

title_map = {
    'TaxonSet1': 'Clade 1',
    'TaxonSet2': 'Clade 2',
    'TaxonSet3': 'Clade 3'
}

for row_idx, analysis in enumerate(analysis_types):
    subset_analysis = df[df['Analysis'] == analysis]
    
    for col_idx, taxon in enumerate(taxon_sets):
        ax = axes[row_idx, col_idx]
        data_to_plot = subset_analysis[subset_analysis['TaxonSet'] == taxon]
        
        if not data_to_plot.empty:
            sns.histplot(
                data=data_to_plot,
                x='Age',
                hue='Condition',
                fill=True,
                palette=colors,
                alpha=0.4,
                ax=ax,
                common_norm=False, 
                element="step",    
                stat="density",    
                linewidth=0.5
            )
            
            display_name = title_map.get(taxon, taxon)
            
            # 2. Increased title fontsize from 14 to 20
            ax.set_title(f"{display_name} ({analysis})", fontsize=20, fontweight='bold')
                
            ax.set_xlabel("Age" if row_idx == 2 else "")
            
            if (row_idx == 0 and col_idx == 2):
                sns.move_legend(ax, "upper right", title=None)
            else:
                if ax.get_legend(): ax.get_legend().remove()
        else:
            ax.set_visible(False)

# 3. Changed the output file extension to .pdf
plt.savefig("cats_comparison_plot.pdf", format='pdf', bbox_inches='tight')

import numpy as np
from scipy.stats import gaussian_kde
from scipy.spatial.distance import jensenshannon

# --- 6. JENSEN-SHANNON DIVERGENCE COMPUTATION ---
jsd_records = []

for analysis in analysis_types:
    subset_analysis = df[df['Analysis'] == analysis]
    
    for taxon in taxon_sets:
        display_name = title_map.get(taxon, taxon)
        data_taxon = subset_analysis[subset_analysis['TaxonSet'] == taxon]
        
        # Extract the continuous age traces for both conditions
        cond_data = data_taxon[data_taxon['Condition'] == 'Conditioned']['Age'].values
        not_cond_data = data_taxon[data_taxon['Condition'] == 'Not Conditioned']['Age'].values
        
        # Only compute if we actually have data for both
        if len(cond_data) > 1 and len(not_cond_data) > 1:
            
            # 1. Create a shared grid of values spanning the minimum and maximum ages
            min_val = min(cond_data.min(), not_cond_data.min())
            max_val = max(cond_data.max(), not_cond_data.max())
            grid = np.linspace(min_val, max_val, 1000)
            
            # 2. Fit KDEs to both datasets and evaluate them on the shared grid
            kde_cond = gaussian_kde(cond_data)(grid)
            kde_not_cond = gaussian_kde(not_cond_data)(grid)
            
            # 3. Normalize the KDE outputs so they sum to 1 (creating discrete probability distributions)
            p = kde_cond / np.sum(kde_cond)
            q = kde_not_cond / np.sum(kde_not_cond)
            
            # 4. Calculate Jensen-Shannon Distance and square it to get JS Divergence
            # scipy returns the distance (square root of divergence), so we square it
            js_dist = jensenshannon(p, q, base=2.0) 
            js_div = js_dist ** 2
            
            jsd_records.append({
                'Analysis': analysis,
                'Clade': display_name,
                'JSD': js_div
            })

# Convert our results into a new DataFrame
df_jsd = pd.DataFrame(jsd_records)

# --- 7. PLOT THE JSD BARPLOT ---
# Create a new figure specifically for the barplot
fig2, ax2 = plt.subplots(figsize=(12, 8))

sns.barplot(
    data=df_jsd,
    x='Analysis',
    y='JSD',
    hue='Clade',
    order=analysis_types,  # Keeps 'Prior', 'Half', 'Full' in chronological order
    palette='viridis',     # A clean, colorblind-friendly palette
    ax=ax2,
    edgecolor='black'      # Adds a crisp border to the bars
)

# Formatting the plot
ax2.set_title('Jensen-Shannon Divergence:\nConditioned vs. Not Conditioned', fontsize=20, fontweight='bold', pad=15)
ax2.set_ylabel('JS Divergence', fontsize=18)
ax2.set_xlabel('Analysis Type', fontsize=18)

# Customize the legend
sns.move_legend(ax2, "upper right", title="Taxon Set", title_fontsize=16, fontsize=14)

# Save the new plot as a PDF
plt.savefig("cats_jsd_barplot.pdf", format='pdf', bbox_inches='tight')

print("JSD Barplot saved as 'cats_jsd_barplot.pdf'")