"""
generate_k_comparison_chart.py: Visualizes spatial and structural
connectivity metrics (\kappa) across biological, linguistic, and cultural domains.
"""
import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

def generate_publication_chart():
    # 1. Establish path registries
    RESULTS_DIR = "results"
    OUTPUT_IMAGE = os.path.join(RESULTS_DIR, "cross_domain_connectivity_comparison.png")

    bio_path = os.path.join(RESULTS_DIR, "biology_connectivity_summary.csv")
    ling_path = os.path.join(RESULTS_DIR, "linguistics_connectivity_summary.csv")
    cult_path = os.path.join(RESULTS_DIR, "culture_connectivity_summary.csv")

    # 2. Gather tracking logs
    plot_data = []

    # Ingest Biology
    if os.path.exists(bio_path):
        df = pd.read_csv(bio_path)
        for _, row in df.iterrows():
            plot_data.append({"Domain": "Biology", "Layer": "Spatial", "Density": 0.0}) # Spatial defaults to NaN/0
            plot_data.append({"Domain": "Biology", "Layer": "Structural", "Density": float(row["Kappa_Structural"])})

    # Ingest Linguistics (Aggregate means across features)
    if os.path.exists(ling_path):
        df = pd.read_csv(ling_path)
        plot_data.append({"Domain": "Linguistics", "Layer": "Spatial", "Density": df["Kappa_Spatial"].mean()})
        plot_data.append({"Domain": "Linguistics", "Layer": "Structural", "Density": df["Kappa_Structural"].mean()})

    # Ingest Culture
    if os.path.exists(cult_path):
        df = pd.read_csv(cult_path)
        for _, row in df.iterrows():
            plot_data.append({"Domain": "Culture", "Layer": "Spatial", "Density": float(row["Kappa_Spatial"])})
            plot_data.append({"Domain": "Culture", "Layer": "Structural", "Density": float(row["Kappa_Structural"])})

    if not plot_data:
        print("Error: No connectivity logs found inside the results directory.")
        return

    df_plot = pd.DataFrame(plot_data)

    # 3. Configure professional visual styling parameters
    sns.set_theme(style="whitegrid")
    plt.rcParams.update({
        "font.family": "serif",
        "font.size": 11,
        "axes.labelsize": 12,
        "axes.titlesize": 13,
        "xtick.labelsize": 11,
        "ytick.labelsize": 11
    })

    fig, ax = plt.subplots(figsize=(7, 5))

    # Render grouping bars cleanly side-by-side
    palette = {"Spatial": "#4A90E2", "Structural": "#E2844A"}
    sns.barplot(
        data=df_plot,
        x="Domain",
        y="Density",
        hue="Layer",
        palette=palette,
        edgecolor="0.2",
        linewidth=1.2,
        ax=ax
    )

    # 4. Refine grid axes boundaries
    ax.set_title("Empirical Adjacency Network Edge Densities ($\kappa$) Across Scientific Domains", pad=15)
    ax.set_xlabel("Research Domain Context", labelpad=10)
    ax.set_ylabel("Matrix Edge Density Scale ($\kappa$)", labelpad=10)
    ax.set_ylim(0.0, 0.7)

    # Add data value labels on top of bars
    for p in ax.patches:
        height = p.get_height()
        if height > 0.0001:  # Skip drawing label elements for zero values
            ax.annotate(f"{height:.4f}",
                        (p.get_x() + p.get_width() / 2., height),
                        ha='center', va='center',
                        xytext=(0, 8),
                        textcoords='offset points',
                        fontsize=10,
                        weight='bold')

    ax.legend(title="Network Layer Profile", loc="upper right", frameon=True)
    sns.despine(left=True, bottom=True)
    plt.tight_layout()

    # Export graphic asset
    plt.savefig(OUTPUT_IMAGE, dpi=300)
    plt.close()
    print(f"📊 Publication comparative visual saved successfully to: {OUTPUT_IMAGE}")

if __name__ == "__main__":
    generate_publication_chart()
