"""
generate_k_comparison_chart.py

reads compiled connectivity summaries from the output directory and saves two
separate linear bar images reflecting independent spatial (0.01) and structural (0.10) limits.
"""

import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

# dynamic path routing: anchors output files inside the results folder
base_dir = os.path.dirname(os.path.abspath(__file__))
input_dir = os.path.join(base_dir, "results")
results_dir = os.path.join(base_dir, "results")
os.makedirs(results_dir, exist_ok=True)

output_img_a = os.path.join(results_dir, "connectivity_geographic_sharing.png")
output_img_b = os.path.join(results_dir, "connectivity_historical_sharing.png")

def load_and_clean_summary(file_name, domain_label):
    """loads a domain summary table and converts n/a string markers to numeric nans."""
    # look for summary files inside root output and nested subfolders
    paths_to_check = [
        os.path.join(input_dir, file_name),
        os.path.join(input_dir, "tables", file_name)
    ]

    file_path = None
    for p in paths_to_check:
        if os.path.exists(p):
            file_path = p
            break

    if file_path is None:
        print(f"Missing file: {file_name}")
        return pd.DataFrame()

    print(f"Found file: Loading {domain_label} dataset from: {os.path.relpath(file_path, base_dir)}")
    df = pd.read_csv(file_path)
    df.replace("N/A", np.nan, inplace=True)

    # standardize runtime case variations across columns
    for col in df.columns:
        if col.lower() == 'spatial_density':
            df['Spatial_Density'] = pd.to_numeric(df[col], errors='coerce')
        if col.lower() == 'structural_density':
            df['Structural_Density'] = pd.to_numeric(df[col], errors='coerce')

    df["Domain"] = domain_label
    return df

def generate_comparison_charts():
    """coordinates the multi-domain comparison graphics generation loops."""
    print("Gathering multi-domain datasets...")

    # load datasets using active summary names
    df_bio = load_and_clean_summary("biology_connectivity_summary.csv", "Biology")
    df_ling = load_and_clean_summary("linguistics_connectivity_summary.csv", "Linguistics")
    df_etc = load_and_clean_summary("culture_connectivity_summary.csv", "Culture")

    master_df = pd.concat([df_bio, df_ling, df_etc], ignore_index=True)

    if master_df.empty:
        print("Error: No data available to plot. Check that your connectivity summary csv files exist.")
        return

    # presentation style parameters
    plt.style.use('seaborn-v0_8-whitegrid' if 'seaborn-v0_8-whitegrid' in plt.style.available else 'default')

    # applied universal okabe-ito palette parameters to guarantee distinct cross-vision visibility
    palette_map = {"Biology": "#E69F00", "Linguistics": "#0072B2", "Culture": "#CC79A7"}

    # =========================================================================
    # IMAGE 1: GEOGRAPHIC TRAIT SHARING (0.01 SPATIAL GEODESIC BOUNDARY)
    # =========================================================================
    df_spatial = master_df.dropna(subset=["Spatial_Density"])
    if not df_spatial.empty:
        print("Generating Image 1: Geographic Sharing with 0.01 baseline floor...")
        fig, ax = plt.subplots(figsize=(7.5, 5.5), dpi=300)

        # explicit mapping to hue fixes deprecation warning gates
        sns.barplot(
            data=df_spatial, x="Domain", y="Spatial_Density", hue="Domain", ax=ax,
            palette=palette_map, legend=False, errorbar="sd", alpha=0.85,
            capsize=0.1, err_kws={'linewidth': 1.5}, edgecolor='black'
        )
        sns.stripplot(
            data=df_spatial, x="Domain", y="Spatial_Density", ax=ax,
            color="#222222", size=5, alpha=0.4, jitter=0.15, dodge=False
        )

        # map unique texture configurations to bars to guarantee greyscale scannability
        for idx, bar in enumerate(ax.patches):
            if idx == 0:
                bar.set_hatch('//')
            elif idx == 1:
                bar.set_hatch('')
            elif idx == 2:
                bar.set_hatch('..')

        # inject decoupled spatial boundary indicator line (0.01 floor) with dark contrast
        x_max_limit = 1.45
        ax.axhline(0.01, color='#000000', linestyle='--', linewidth=1.2, alpha=0.75)
        ax.text(x_max_limit, 0.01 + 0.005, 'Spatial Frontier Floor (0.01)',
                color='#000000', fontsize=8, fontweight='bold', ha='right', va='bottom')

        # capture scale window to show culture metrics and low spatial language bounds
        ax.set_ylim(-0.01, 0.25)

        ax.set_title("Geographic Trait Sharing\n(Do neighbors match each other?)", fontsize=12, fontweight="bold", pad=12)
        ax.set_xlabel("Evolutionary Domain", fontsize=10, fontweight="bold")
        ax.set_ylabel("Connection Strength (κ_spatial)", fontsize=10, fontweight="bold")

        plt.tight_layout()
        plt.savefig(output_img_a, bbox_inches="tight")
        plt.close()
        print(f"Saved geographic sharing chart to: '{os.path.relpath(output_img_a, base_dir)}'")

    # =========================================================================
    # IMAGE 2: HISTORICAL TREE SHARING (0.10 STRUCTURAL SAFETY BARRIER)
    # =========================================================================
    df_struct = master_df.dropna(subset=["Structural_Density"])
    if not df_struct.empty:
        print("Generating Image 2: Historical Tree Sharing with 0.10 safety floor...")
        fig, ax = plt.subplots(figsize=(7.5, 5.5), dpi=300)

        # explicit mapping to hue fixes deprecation warning gates
        sns.barplot(
            data=df_struct, x="Domain", y="Structural_Density", hue="Domain", ax=ax,
            palette=palette_map, legend=False, errorbar="sd", alpha=0.85,
            capsize=0.1, err_kws={'linewidth': 1.5}, edgecolor='black'
        )
        sns.stripplot(
            data=df_struct, x="Domain", y="Structural_Density", ax=ax,
            color="#222222", size=4, alpha=0.3, jitter=0.2, dodge=False
        )

        # map unique texture configurations to bars to guarantee greyscale scannability
        for idx, bar in enumerate(ax.patches):
            if idx == 0:
                bar.set_hatch('\\\\')
            elif idx == 1:
                bar.set_hatch('//')
            elif idx == 2:
                bar.set_hatch('')

        # inject decoupled structural safety threshold line (0.10 floor) with high-contrast color
        x_max_limit = 2.45
        ax.axhline(0.10, color='#000000', linestyle='--', linewidth=1.2, alpha=0.75)
        ax.text(x_max_limit, 0.10 + 0.015, 'Sparse Topology Ceiling (0.10)',
                color='#000000', fontsize=8, fontweight='bold', ha='right', va='bottom')

        # inject decoupled structural density buffer line (0.25 floor) with a distinct dot style
        x_max_limit = 2.45
        ax.axhline(0.25, color='#444444', linestyle=':', linewidth=1.5, alpha=0.85)
        ax.text(x_max_limit, 0.25 + 0.015, 'Dense Topology Floor (0.25)',
                color='#444444', fontsize=8, fontweight='bold', ha='right', va='bottom')

        # extend structural upper limit to clear high reptile benchmarks safely
        ax.set_ylim(-0.03, 1.05)

        ax.set_title("Historical Tree Sharing\n(Do family tree branches match each other?)", fontsize=12, fontweight="bold", pad=12)
        ax.set_xlabel("Evolutionary Domain", fontsize=10, fontweight="bold")
        ax.set_ylabel("Connection Strength (κ_structural)", fontsize=10, fontweight="bold")

        plt.tight_layout()
        plt.savefig(output_img_b, bbox_inches="tight")
        plt.close()
        print(f"Saved historical sharing chart to: '{os.path.relpath(output_img_b, base_dir)}'")

if __name__ == "__main__":
    generate_comparison_charts()
