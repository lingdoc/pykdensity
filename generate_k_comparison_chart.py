"""
generate_k_comparison_chart.py

Reads compiled connectivity summaries from the results directory and outputs two
separate linear bar images containing both the 0.005 and 0.15 pykdensity thresholds.
"""

import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

# 📍 DYNAMIC PATH ROUTING: Anchors directory locations relative to this script
base_dir = os.path.dirname(os.path.abspath(__file__))
results_dir = os.path.join(base_dir, "results")
output_img_a = os.path.join(results_dir, "connectivity_geographic_sharing.png")
output_img_b = os.path.join(results_dir, "connectivity_historical_sharing.png")

def load_and_clean_summary(file_name, domain_label):
    """Loads a domain summary table and converts N/A string markers to numeric NaNs."""
    file_path = os.path.join(results_dir, file_name)
    if not os.path.exists(file_path):
        return pd.DataFrame()

    df = pd.read_csv(file_path)
    df.replace("N/A", np.nan, inplace=True)
    df["Spatial_Density"] = pd.to_numeric(df["Spatial_Density"], errors='coerce')
    df["Structural_Density"] = pd.to_numeric(df["Structural_Density"], errors='coerce')
    df["Domain"] = domain_label
    return df

def add_unified_threshold_lines(ax, x_max_limit):
    """Draws both the 0.005 and 0.15 benchmark thresholds across the canvas panel."""
    # 📍 Threshold 1: Sparse Significance Boundary (0.005)
    ax.axhline(0.005, color='#e6550d', linestyle='--', linewidth=1.2, alpha=0.75)
    ax.text(x_max_limit, 0.005 + 0.003, 'Significance Boundary (0.005)',
            color='#a63603', fontsize=8, fontweight='bold', ha='right', va='bottom')

    # 📍 Threshold 2: Dense Saturation Floor (0.15)
    ax.axhline(0.15, color='#d62728', linestyle='--', linewidth=1.2, alpha=0.75)
    ax.text(x_max_limit, 0.15 - 0.003, 'Saturation Floor (0.15)',
            color='#b51a1a', fontsize=8, fontweight='bold', ha='right', va='top')

def generate_comparison_charts():
    print("Gathering multi-domain datasets...")

    df_bio = load_and_clean_summary("biology_connectivity_summary.csv", "Biology")
    df_ling = load_and_clean_summary("linguistics_connectivity_summary.csv", "Linguistics")
    df_etc = load_and_clean_summary("culture_connectivity_summary.csv", "Culture")

    master_df = pd.concat([df_bio, df_ling, df_etc], ignore_index=True)

    if master_df.empty:
        print("Error: No data available to plot. Run example_connectivity_pipeline.py first.")
        return

    # Presentation style parameters
    plt.style.use('seaborn-v0_8-whitegrid' if 'seaborn-v0_8-whitegrid' in plt.style.available else 'default')
    palette_map = {"Biology": "#2ca02c", "Linguistics": "#1f77b4", "Culture": "#e377c2"}

    # =========================================================================
    # IMAGE 1: GEOGRAPHIC TRAIT SHARING (LINEAR SCALE)
    # =========================================================================
    df_spatial = master_df.dropna(subset=["Spatial_Density"])
    if not df_spatial.empty:
        print("Generating Image 1: Geographic Sharing with dual thresholds...")
        fig, ax = plt.subplots(figsize=(7.5, 5.5), dpi=300)

        sns.barplot(
            data=df_spatial, x="Domain", y="Spatial_Density", ax=ax,
            palette=palette_map, errorbar="sd", alpha=0.85, capsize=0.1, err_kws={'linewidth': 1.5}
        )
        sns.stripplot(
            data=df_spatial, x="Domain", y="Spatial_Density", ax=ax,
            color="#333333", size=5, alpha=0.4, jitter=0.15, dodge=False
        )

        # Inject standard baseline limits
        add_unified_threshold_lines(ax, x_max_limit=1.45)

        # Pad top slightly past 0.15 so the text is fully legible
        ax.set_ylim(-0.005, 0.17)

        ax.set_title("Geographic Trait Sharing\n(Do neighbors match each other?)", fontsize=12, fontweight="bold", pad=12)
        ax.set_xlabel("Evolutionary Domain", fontsize=10, fontweight="bold")
        ax.set_ylabel("Connection Strength", fontsize=10, fontweight="bold")

        plt.tight_layout()
        plt.savefig(output_img_a, bbox_inches="tight")
        plt.close()
        print(f" Saved -> '{os.path.relpath(output_img_a, base_dir)}'")

    # =========================================================================
    # IMAGE 2: HISTORICAL TREE SHARING (LINEAR SCALE)
    # =========================================================================
    df_struct = master_df.dropna(subset=["Structural_Density"])
    if not df_struct.empty:
        print("Generating Image 2: Historical Tree Sharing with dual thresholds...")
        fig, ax = plt.subplots(figsize=(7.5, 5.5), dpi=300)

        sns.barplot(
            data=df_struct, x="Domain", y="Structural_Density", ax=ax,
            palette=palette_map, errorbar="sd", alpha=0.85, capsize=0.1, err_kws={'linewidth': 1.5}
        )
        sns.stripplot(
            data=df_struct, x="Domain", y="Structural_Density", ax=ax,
            color="#333333", size=4, alpha=0.3, jitter=0.2, dodge=False
        )

        # Inject standard baseline limits
        add_unified_threshold_lines(ax, x_max_limit=2.45)

        # Pad up past Biology's structural peaks (~0.65)
        ax.set_ylim(-0.02, 0.72)

        ax.set_title("Historical Tree Sharing\n(Do family tree branches match each other?)", fontsize=12, fontweight="bold", pad=12)
        ax.set_xlabel("Evolutionary Domain", fontsize=10, fontweight="bold")
        ax.set_ylabel("Connection Strength", fontsize=10, fontweight="bold")

        plt.tight_layout()
        plt.savefig(output_img_b, bbox_inches="tight")
        plt.close()
        print(f" Saved -> '{os.path.relpath(output_img_b, base_dir)}'")

if __name__ == "__main__":
    generate_comparison_charts()
