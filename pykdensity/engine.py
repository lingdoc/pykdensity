import os
import io
import gzip
import re
import numpy as np
import pandas as pd
from scipy.spatial.distance import pdist, squareform
from Bio import Phylo

def _compute_edge_density(A):
    """Calculates the mathematical edge density ratio ignoring self-loops."""
    n = len(A)
    if n < 2:
        return 0.0
    np.fill_diagonal(A, 0.0)
    total_possible_edges = n * (n - 1)
    total_active_edges = np.sum(A)
    return float(total_active_edges / total_possible_edges)

def _parse_unified_tree(tree_path, df, id_col, tree_type='fixed', depth_threshold=3):
    """
    Unified tree ingestion engine.
    Extracts the Newick structural data and translation block from raw or gzipped files,
    maps tracking identifiers safely, and builds the network based on the tree_type mode.
    """
    n = len(df)
    df_clean = df.copy()

    # Standardize lookup keys to lowercase strings
    df_clean['match_key'] = df_clean[id_col].astype(str).str.strip().str.lower()
    species_list = df_clean['match_key'].tolist()

    # 1. Read the file stream, checking for compression format dynamically
    open_func = gzip.open if str(tree_path).endswith('.gz') else open

    translate_map = {}
    tree_string = ""
    in_translate = False

    with open_func(tree_path, 'rt', encoding='utf-8', errors='ignore') as f:
        for line in f:
            line_clean = line.strip()
            if not line_clean:
                continue

            # Handle embedded Nexus translation blocks if present
            if line_clean.lower().startswith("translate"):
                in_translate = True
                continue
            if in_translate and line_clean == ";":
                in_translate = False
                continue
            if in_translate:
                parts = re.split(r'\s+', line_clean.rstrip(',').rstrip(';'))
                if len(parts) >= 2:
                    idx_token = parts[0].strip()
                    raw_label = parts[1].strip().replace("'", "").replace('"', '')
                    # Isolate prefixes before trailing metadata underscores
                    clean_label = raw_label.split('_')[0] if '_' in raw_label else raw_label
                    translate_map[idx_token] = clean_label.lower()
                continue

            # Capture the Newick structural layout string
            if line_clean.lower().startswith("tree ") or (line_clean.startswith("(") and line_clean.endswith(";")):
                tree_string = line_clean.split("=", 1)[1].strip() if "=" in line_clean else line_clean
                break

    if not tree_string:
        raise ValueError(f"Could not extract a valid Newick topology from: {tree_path}")

    # Clean up standard metadata comment syntax blocks
    tree_string = re.sub(r'\[.*?\]', '', tree_string)

    # 2. Extract ancestral trace paths from the true tree layout
    lineage_paths_cache = {}

    try:
        # High-performance tree construction parsing via memory stream
        target_tree = Phylo.read(io.StringIO(tree_string), "newick")

        # Build node path lists down to every terminal leaf element
        for tip in target_tree.get_terminals():
            if not tip.name:
                continue
            tip_name = tip.name.strip().lower()

            # Resolve numeric token IDs if a translation block exists
            resolved_key = translate_map.get(tip_name, tip_name)
            # Remove spaces or underscores to maximize lookup matching rates
            resolved_key = resolved_key.replace(' ', '_')

            lineage_paths_cache[resolved_key] = [target_tree.root] + target_tree.get_path(tip)

    except Exception:
        # Fallback manual string parser if Bio.Phylo encounters complex layout variants
        node_counter = 10000
        stack = [0]

        for i, char in enumerate(tree_string):
            if char == '(':
                stack.append(node_counter)
                node_counter += 1
            elif char == ')':
                if len(stack) > 1:
                    stack.pop()
            elif char not in [',', ';', ':', ' ']:
                match = re.match(r'^([\d\w\.\-_]+)', tree_string[i:])
                if match:
                    token = match.group(1).strip().lower()
                    resolved_key = translate_map.get(token, token).replace(' ', '_')
                    lineage_paths_cache[resolved_key] = list(stack)

    # 3. Compute structural proximity adjacency connections
    A = np.zeros((n, n))
    mode = str(tree_type).strip().lower()

    # Pre-calculate paths for matching lookups
    cached_paths = []
    for sp in species_list:
        # Fallback match variants for text formatting quirks
        alt_sp = sp.replace(' ', '_')
        path = lineage_paths_cache.get(sp, lineage_paths_cache.get(alt_sp, [0]))
        cached_paths.append(path)

    for i in range(n):
        path_i = cached_paths[i]
        len_i = len(path_i)

        for j in range(i + 1, n):
            path_j = cached_paths[j]
            len_j = len(path_j)

            min_length = min(len_i, len_j)
            shared_depth = 0
            for k in range(min_length):
                if path_i[k] == path_j[k]:
                    shared_depth += 1
                else:
                    break

            # ROUTE A: Adaptive Proportional Ratio Calculation
            if mode == 'adaptive':
                ratio_threshold = min(max(depth_threshold * 0.10, 0.05), 0.95)
                rel_i = shared_depth / len_i if len_i > 0 else 0.0
                rel_j = shared_depth / len_j if len_j > 0 else 0.0

                if rel_i > ratio_threshold and rel_j > ratio_threshold:
                    A[i, j] = 1.0
                    A[j, i] = 1.0

            # ROUTE B: Fixed Direct Node Count Step Calculation
            else:
                if shared_depth > depth_threshold:
                    A[i, j] = 1.0
                    A[j, i] = 1.0

    return A

def calculate_densities(data, id_col, tree=None, tree_type='fixed', taxonomy_hierarchy_cols=None, coord_cols=None, theme_cols=None, spatial_threshold_km=500.0, structural_depth_threshold=3, verbose=True):
    """Unified engine function. Computes metrics cleanly without hiding structural errors."""
    df = data.copy()
    n = len(df)

    spatial_density = np.nan
    structural_density = np.nan

    if n < 2:
        return spatial_density, structural_density

    # 1. EVALUATE SPATIAL LAYER
    if coord_cols and not theme_cols:
        if len(coord_cols) == 3:
            coords = df[coord_cols].to_numpy().astype(float)
        else:
            R = 6371.0
            lat_rad = np.radians(df[coord_cols].to_numpy().astype(float))
            lon_rad = np.radians(df[coord_cols].to_numpy().astype(float))
            x = R * np.cos(lat_rad) * np.cos(lon_rad)
            y = R * np.cos(lat_rad) * np.sin(lon_rad)
            z = R * np.sin(lat_rad)
            coords = np.column_stack((x, y, z))

        pairwise_dist = squareform(pdist(coords, metric='euclidean'))
        A_space = (pairwise_dist <= spatial_threshold_km).astype(float)
        spatial_density = _compute_edge_density(A_space)

    elif coord_cols and len(coord_cols) == 1 and coord_cols[0] in df.columns and theme_cols:
        geo_vector = df[coord_cols[0]].to_numpy()
        A_space = (geo_vector[:, None] == geo_vector[None, :]).astype(float)
        spatial_density = _compute_edge_density(A_space)

    # 2. EVALUATE STRUCTURAL LAYER
    has_valid_tree = False

    if tree is not None and os.path.exists(tree):
        try:
            A_struct = _parse_unified_tree(
                tree_path=tree,
                df=df,
                id_col=id_col,
                tree_type=tree_type,
                depth_threshold=structural_depth_threshold
            )
            structural_density = _compute_edge_density(A_struct)
            has_valid_tree = True
        except Exception as e:
            print(f"Warning: Tree mapping parser failed: {e}")
            has_valid_tree = False

    # Taxonomy fallback routes are ONLY triggered if no physical tree asset was supplied
    if not has_valid_tree and taxonomy_hierarchy_cols is not None and tree is None:
        A_struct = np.zeros((n, n))
        for col in taxonomy_hierarchy_cols:
            col_vector = df[col].to_numpy()
            A_struct += (col_vector[:, None] == col_vector[None, :]).astype(float)
        A_struct = (A_struct > 0).astype(float)
        structural_density = _compute_edge_density(A_struct)

    elif not has_valid_tree and theme_cols is not None and tree is None:
        theme_matrix = df[theme_cols].to_numpy()
        A_struct = (theme_matrix[:, None, :] == theme_matrix[None, :, :]).all(axis=2).astype(float)
        structural_density = _compute_edge_density(A_struct)

    if verbose:
        print(f"📊 Derived Spatial Density   : {spatial_density if pd.isna(spatial_density) else f'{spatial_density:.4f}'}")
        print(f"🗂️  Derived Structural Density: {structural_density if pd.isna(structural_density) else f'{structural_density:.4f}'}")

    return spatial_density, structural_density
