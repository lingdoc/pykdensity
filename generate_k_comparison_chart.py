"""
generate_k_comparison_chart.py

reads the summary files from the results folder and generates a bar chart
comparing spatial and structural data connectivity across domains.
"""
import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

def generate_comparison_chart():
    # setup paths for input files and the final image output
    results_dir = "results"
    output_image = os.path.join(results_dir, "cross_domain_connectivity_comparison.png")

    bio_path = os.path.join(results_dir, "biology_connectivity_summary.csv")
    ling_path = os.path.join(results_dir, "linguistics_connectivity_summary.csv")
    cult_path = os.path.join(results_dir, "culture_connectivity_summary.csv")

    # collect the rows from each summary table
    plot_data = []

    # load the biology summary rows
    if os.path.exists(bio_path):
        df = pd.read_csv(bio_path)
        for _, row in df.iterrows():
            plot_data.append({"Domain": "Biology", "Type": "Spatial", "Density": 0.0}) # spatial tracks as nan/0 here
            plot_data.append({"Domain": "Biology", "Type": "Structural", "Density": float(row["Structural_Density"])})

    # load the linguistics summary rows and calculate averages across features
    if os.path.exists(ling_path):
        df = pd.read_csv(ling_path)
        plot_data.append({"Domain": "Linguistics", "Type": "Spatial", "Density": df["Spatial_Density"].mean()})
        plot_data.append({"Domain": "Linguistics", "Type": "Structural", "Density": df["Structural_Density"].mean()})

    # load the culture summary rows
    if os.path.exists(cult_path):
        df = pd.read_csv(cult_path)
        for _, row in df.iterrows():
            plot_data.append({"Domain": "Culture", "Type": "Spatial", "Density": float(row["Spatial_Density"])})
            plot_data.append({"Domain": "Culture", "Type": "Structural", "Density": float(row["Structural_Density"])})

    if not plot_data:
        print("Error: no summary files found in the results folder.")
        return

    df_plot = pd.DataFrame(plot_data)

    # set standard chart text sizes and grid properties
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

    # draw the bars side by side using distinct colors
    palette = {"Spatial": "#4A90E2", "Structural": "#E2844A"}
    sns.barplot(
        data=df_plot,
        x="Domain",
        y="Density",
        hue="Type",
        palette=palette,
        edgecolor="0.2",
        linewidth=1.2,
        ax=ax
    )

    # configure title, axis labels, and chart limits
    ax.set_title("Data Connectivity Levels Across Fields of Study", pad=15)
    ax.set_xlabel("Field of Study", labelpad=10)
    ax.set_ylabel("Data Density Score", labelpad=10)
    ax.set_ylim(0.0, 0.7)

    # draw data value labels on top of each individual bar
    for p in ax.patches:
        height = p.get_height()
        if height > 0.0001:  # ignore labels for zero or empty values
            ax.annotate(f"{height:.4f}",
                        (p.get_x() + p.get_width() / 2., height),
                        ha='center', va='center',
                        xytext=(0, 8),
                        textcoords='offset points',
                        fontsize=10,
                        weight='bold')

    ax.legend(title="Connection Type", loc="upper right", frameon=True)
    sns.despine(left=True, bottom=True)
    plt.tight_layout()

    # save the image file to disk
    plt.savefig(output_image, dpi=300)
    plt.close()
    print(f"Comparison chart saved to: {output_image}")

if __name__ == "__main__":
    generate_comparison_chart()
