"""
example_connectivity_pipeline.py

This script shows how to run the connectivity package across three separate
fields of study. Each section represents a different data setup.
"""

import os
import glob
import pandas as pd
import numpy as np

# Import the core calculator function from our package
from pykdensity import calculate_densities

# 📍 DYNAMIC PATH ROUTING: anchors all folders relative to this script's location
base_dir = os.path.dirname(os.path.abspath(__file__))
data_dir = os.path.join(base_dir, "data")
results_dir = os.path.join(base_dir, "results")
os.makedirs(results_dir, exist_ok=True)

print("\nStarting cross-domain network data pipeline")

# 1. Biological data example
# This step handles family trees where you have matching species traits but no
# geographic location information. We tell the package to use the 'fixed' method,
# which counts connections as simple steps down from the root of the tree.
print("\n" + "-" * 80)
print("\n Processing domain: vertebrate ecology (biology)\n")

# Find all biological spreadsheet tables in the bio folder via dynamic routing
bio_csvs = sorted(glob.glob(os.path.join(data_dir, "bio", "*.csv")))
bio_results = []

for csv_path in bio_csvs:
    # Pull out the file name to use as a label in the final table
    label = str(os.path.basename(csv_path).split(".")[0])

    # Search for a matching tree file ending in .nex
    tree_path = csv_path.replace(".csv", "100trees.nex")
    if not os.path.exists(tree_path):
        tree_path = csv_path.replace(".csv", ".nex")

    local_csv_path = os.path.relpath(csv_path, base_dir)
    local_tree_path = os.path.relpath(tree_path, base_dir) if os.path.exists(tree_path) else "None"

    print(f"\n -> Processing file: {label}")
    print(f"    Spreadsheet: {local_csv_path}")
    print(f"    Tree file: {local_tree_path}")

    # load the table into python
    df_bio = pd.read_csv(csv_path)

    # run the analysis using tree steps while explicitly ignoring geographic settings
    k_spatial, k_structural = calculate_densities(
        data=df_bio,
        id_col='Species',
        tree=tree_path if os.path.exists(tree_path) else None,
        tree_type='fixed',  # measures relationships by counting shared family steps
        taxonomy_hierarchy_cols=['Order', 'Family'],
        structural_depth_threshold=5, # depth of phylogenetic/tree relations
        coord_cols=None,    # left blank because this dataset has no map coordinates
        theme_cols=None,
        verbose=True
    )

    bio_results.append({
        "Domain": "Biology",
        "Label": label,
        "Sample_Size_N": len(df_bio),
        "Kappa_Spatial": k_spatial if not pd.isna(k_spatial) else "N/A",
        "Kappa_Structural": k_structural if not pd.isna(k_structural) else "N/A"
    })

# save the biological summary table
if bio_results:
    bio_output_path = os.path.join(results_dir, "biology_connectivity_summary.csv")
    pd.DataFrame(bio_results).to_csv(bio_output_path, index=False)
    print(f"\nBiology results saved to: {os.path.relpath(bio_output_path, base_dir)}")

# 2. Linguistic data example
# This step handles individual language traits across folders. It combines real
# map positions (latitude and longitude numbers) with compressed family trees.
# We use the 'adaptive' setting so the package automatically handles uneven
# branch lengths across different language families.
print("\n" + "-" * 80)
print("\n Processing domain: Linguistic typology (linguistics)\n")

glottolog_path = os.path.join(data_dir, "tlu", "Glottolog_Languages.csv")
ling_results = []

# Ensure the main language map file is present before looping
if os.path.exists(glottolog_path):
    print(f" -> Loading language map files from: {os.path.relpath(glottolog_path, base_dir)}")
    gldf = pd.read_csv(glottolog_path)

    # Clean the data by keeping language codes and filtering out rows missing coordinates
    gldf = gldf.dropna(subset=['latitude', 'longitude'])
    gldf_shared = gldf[['glottocode', 'macroarea', 'Family_ID', 'longitude', 'latitude']].copy()
    gldf_shared['glottocode'] = gldf_shared['glottocode'].astype(str).str.strip().str.lower()

    # Locate text files inside separate nested subfolders
    feat_data_paths = sorted(glob.glob(os.path.join(data_dir, "tlu", "*", "*data.txt")))

    for feat_path in feat_data_paths:
        # Extract the folder name to identify the specific feature being tested
        feat_id = str(os.path.basename(os.path.dirname(feat_path)))

        # Look for a compressed tree file in the same folder
        feat_dir = os.path.dirname(feat_path)
        tree_matches = glob.glob(os.path.join(feat_dir, "*trees.gz"))
        tree_path = tree_matches[0] if tree_matches else None

        print(f"\n -> Assessing language feature: {feat_id}")

        # open the feature data file
        fdf = pd.read_csv(feat_path, delimiter="\t", header=None, names=["glottocode", "DV", "IV"])

        # skip folders that do not have enough data to be statistically meaningful
        if len(fdf) < 10:
            print(f"    Skipping: not enough data rows (found {len(fdf)}, needs at least 10)")
            continue

        # line up the feature text rows with the map coordinate master file
        fdf['glottocode'] = fdf['glottocode'].astype(str).str.strip().str.lower()
        ling_df = pd.merge(gldf_shared, fdf, on='glottocode', how='inner')

        local_tree_print = os.path.relpath(tree_path, base_dir) if tree_path else "None"
        print(f"    Matched {len(ling_df)} languages on the map.")
        print(f"    Tree file found: {local_tree_print}")

        # run the calculator using kilometer distances and dynamic tree branches
        k_spatial, k_structural = calculate_densities(
            data=ling_df,
            id_col='glottocode',
            tree=tree_path,
            tree_type='adaptive',  # uses percentages to handle varying tree lengths
            coord_cols=['latitude', 'longitude'], # 2D coordinates are converted to 3D on the backend
            spatial_threshold_km=500.0,    # flags items as connected if they are within 500km
            structural_depth_threshold=8, # sets the threshold level for historical branch similarity
            # A value of 15 pushes the adaptive ratio to its maximum 95% cap.
            # This ensures languages are only marked as connected if they share
            # almost their entire history, isolating tight local sub-families.
            # For a more realistic assessment, use a value like 7 or 8 (large families).
            verbose=True
        )

        ling_results.append({
            "Domain": "Linguistics",
            "Label": feat_id,
            "Sample_Size_N": len(ling_df),
            "Kappa_Spatial": k_spatial if not pd.isna(k_spatial) else "N/A",
            "Kappa_Structural": k_structural if not pd.isna(k_structural) else "N/A"
        })

    # save the linguistic findings
    if ling_results:
        ling_output_path = os.path.join(results_dir, "linguistics_connectivity_summary.csv")
        pd.DataFrame(ling_results).to_csv(ling_output_path, index=False)
        print(f"\nLinguistics results saved to: {os.path.relpath(ling_output_path, base_dir)}")
else:
    print(f"Cannot run linguistics loop: file missing at {os.path.relpath(glottolog_path, base_dir)}")

# 3. Cultural data example
# This step handles historical datasets that lack an explicit family tree file.
# Instead, we group data by matching named categories. Space is grouped by
# region names ('Zone'), and history is grouped by matching eras and formats.
print("\n" + "-" * 80)
print("\n Processing domain: Cultural evolution (culture)\n")

# Find excel spreadsheets in the etc folder
etc_excel_paths = sorted(glob.glob(os.path.join(data_dir, "etc", "*.xlsx")))
etc_results = []

for excel_path in etc_excel_paths:
    # use the spreadsheet filename as the tracker row title
    label = str(os.path.basename(excel_path).split(".")[0])
    local_excel_path = os.path.relpath(excel_path, base_dir)
    print(f" -> Opening spreadsheet: {label}")
    print(f"    Path: {local_excel_path}")

    # load the excel rows directly into python
    df_etc = pd.read_excel(excel_path)

    # run the analysis using text groups exclusively by passing tree=None
    k_spatial, k_structural = calculate_densities(
        data=df_etc,
        id_col='Name' if 'Name' in df_etc.columns else df_etc.columns[0],
        tree=None,              # skipped because no explicit tree file exists for this track
        tree_type='fixed',      # default setting ignored since no tree file is supplied
        coord_cols=['Zone'] if 'Zone' in df_etc.columns else None,  # groups items by shared region name
        theme_cols=['Century', 'Genre.1'] if 'Century' in df_etc.columns else None, # groups items by era and genre
        verbose=True
    )

    etc_results.append({
        "Domain": "Culture",
        "Label": label,
        "Sample_Size_N": len(df_etc),
        "Kappa_Spatial": k_spatial if not pd.isna(k_spatial) else "N/A",
        "Kappa_Structural": k_structural if not pd.isna(k_structural) else "N/A"
    })

# save the cultural metrics
if etc_results:
    etc_output_path = os.path.join(results_dir, "culture_connectivity_summary.csv")
    pd.DataFrame(etc_results).to_csv(etc_output_path, index=False)
    print(f"\nCulture results saved to: {os.path.relpath(etc_output_path, base_dir)}")

# 4. complete
print("\n" + "-" * 80)
print("Analysis complete. All files have been created in the results folder.\n")
