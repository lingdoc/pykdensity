"""
example_connectivity_pipeline.py

this script shows how to run the connectivity package across three separate
fields of study. each section represents a different data setup.
"""

import os
import glob
import pandas as pd
import numpy as np

# import the core calculator function from our package
from pykdensity import calculate_densities

# anchors all folders relative to this script's location
base_dir = os.path.dirname(os.path.abspath(__file__))
data_dir = os.path.join(base_dir, "data")
results_dir = os.path.join(base_dir, "results")
os.makedirs(results_dir, exist_ok=True)

print("\nStarting cross-domain network data pipeline")

# 1. biological data example
# this step handles family trees where you have matching species traits but no
# geographic location information. we tell the package to use the 'tree' mode,
# which looks only at the tree branches and ignores any text categories.
print("\n" + "-" * 80)
print("\n Processing domain: vertebrate ecology (biology)\n")

# find all biological spreadsheet tables in the bio folder via dynamic routing
bio_csvs = sorted(glob.glob(os.path.join(data_dir, "bio", "*.csv")))
bio_results = []

for csv_path in bio_csvs:
    # pull out the file name to use as a label in the final table
    label = str(os.path.basename(csv_path).split(".")[0])

    # search for a matching tree file ending in .nex
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

    # run the analysis using pure tree lines while explicitly ignoring geographic settings
    k_spatial, k_structural = calculate_densities(
        data=df_bio,
        id_col='Species',
        tree=tree_path if os.path.exists(tree_path) else None,
        tree_type='fixed',  # measures relationships by counting shared family steps
        struct_mode='tree', # tells the engine to look only at the tree file branches
        struct_depth=5, # depth of tree relations
        verbose=True
    )

    bio_results.append({
        "Domain": "Biology",
        "Label": label,
        "Sample_Size_N": len(df_bio),
        "Spatial_Density": k_spatial if not pd.isna(k_spatial) else "N/A",
        "Structural_Density": k_structural if not pd.isna(k_structural) else "N/A"
    })

# save the biological summary table
if bio_results:
    bio_output_path = os.path.join(results_dir, "biology_connectivity_summary.csv")
    pd.DataFrame(bio_results).to_csv(bio_output_path, index=False)
    print(f"\nBiology results saved to: {os.path.relpath(bio_output_path, base_dir)}")

# 2. linguistic data example
# this step handles individual language traits across folders. it combines real
# map positions (latitude and longitude numbers) with compressed family trees.
# we use the 'adaptive' setting so the package automatically handles uneven
# branch lengths across different language families.
print("\n" + "-" * 80)
print("\n Processing domain: linguistic typology (linguistics)\n")

glottolog_path = os.path.join(data_dir, "tlu", "Glottolog_Languages.csv")
ling_results = []

# ensure the main language map file is present before looping
if os.path.exists(glottolog_path):
    print(f" -> Loading language map files from: {os.path.relpath(glottolog_path, base_dir)}")
    gldf = pd.read_csv(glottolog_path)

    # clean the data by keeping language codes and filtering out rows missing coordinates
    gldf = gldf.dropna(subset=['latitude', 'longitude'])
    gldf_shared = gldf[['glottocode', 'macroarea', 'Family_ID', 'longitude', 'latitude']].copy()
    gldf_shared['glottocode'] = gldf_shared['glottocode'].astype(str).str.strip().str.lower()

    # locate text files inside separate nested subfolders
    feat_data_paths = sorted(glob.glob(os.path.join(data_dir, "tlu", "*", "*data.txt")))

    for feat_path in feat_data_paths:
        # extract the folder name to identify the specific feature being tested
        feat_id = str(os.path.basename(os.path.dirname(feat_path)))

        # look for a compressed tree file in the same folder
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

        # run the calculator using kilometer distances and dynamic tree branches, combine text categories
        k_spatial, k_structural = calculate_densities(
            data=ling_df,
            id_col='glottocode',
            tree=tree_path,
            tree_type='adaptive',  # uses percentages to handle varying tree lengths
            struct_depth=8, # sets the threshold level for historical branch similarity
            coord_cols=['latitude', 'longitude'], # coordinates are handled on the backend
            spatial_km=500.0,    # flags items as connected if they are within 500km
            spatial_cats='macroarea',      # filters space by matching global region text names
            struct_cats='Family_ID',       # links isolated items by matching language family name labels
            verbose=True
        )

        ling_results.append({
            "Domain": "Linguistics",
            "Label": feat_id,
            "Sample_Size_N": len(ling_df),
            "Spatial_Density": k_spatial if not pd.isna(k_spatial) else "N/A",
            "Structural_Density": k_structural if not pd.isna(k_structural) else "N/A"
        })

    # save the linguistic findings
    if ling_results:
        ling_output_path = os.path.join(results_dir, "linguistics_connectivity_summary.csv")
        pd.DataFrame(ling_results).to_csv(ling_output_path, index=False)
        print(f"\nLinguistics results saved to: {os.path.relpath(ling_output_path, base_dir)}")
else:
    print(f"Cannot run linguistics loop: file missing at {os.path.relpath(glottolog_path, base_dir)}")

# 3. cultural data example
# this step handles historical datasets that lack an explicit family tree file.
# instead, we group data by matching named categories. space is grouped by
# region names ('Zone'), and history is grouped by matching eras and formats.
print("\n" + "-" * 80)
print("\n Processing domain: cultural evolution (culture)\n")

# find excel spreadsheets in the etc folder
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
        spatial_cats='Zone' if 'Zone' in df_etc.columns else None,  # groups items by shared region name
        struct_cats=['Century', 'Genre.1'] if 'Century' in df_etc.columns else None, # groups items by era and genre lists
        verbose=True
    )

    etc_results.append({
        "Domain": "Culture",
        "Label": label,
        "Sample_Size_N": len(df_etc),
        "Spatial_Density": k_spatial if not pd.isna(k_spatial) else "N/A",
        "Structural_Density": k_structural if not pd.isna(k_structural) else "N/A"
    })

# save the cultural metrics
if etc_results:
    etc_output_path = os.path.join(results_dir, "culture_connectivity_summary.csv")
    pd.DataFrame(etc_results).to_csv(etc_output_path, index=False)
    print(f"\nCulture results saved to: {os.path.relpath(etc_output_path, base_dir)}")

# 4. complete
print("\n" + "-" * 80)
print("Analysis complete. all files have been created in the results folder.\n")
