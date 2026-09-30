import os
import io
import gzip
import re
import numpy as np
import pandas as pd
from scipy.spatial.distance import pdist, squareform
from Bio import Phylo

def _compute_edge_density(A):
    """calculates the ratio of active network connections, ignoring self-loops."""
    n = len(A)
    if n < 2:
        return 0.0
    np.fill_diagonal(A, 0.0)
    total_possible_edges = n * (n - 1)
    total_active_edges = np.sum(A)
    return float(total_active_edges / total_possible_edges)

def _parse_unified_tree(tree_path, df, id_col, tree_type='fixed', depth_threshold=3):
    """reads a tree file, extracts the newick structure, maps identifiers, and builds an adjacency matrix."""
    n = len(df)
    df_clean = df.copy()
    df_clean['match_key'] = df_clean[id_col].astype(str).str.strip().str.lower()
    species_list = df_clean['match_key'].tolist()

    open_func = gzip.open if str(tree_path).endswith('.gz') else open
    translate_map = {}
    tree_string = ""
    in_translate = False

    with open_func(tree_path, 'rt', encoding='utf-8', errors='ignore') as f:
        for line in f:
            line_clean = line.strip()
            if not line_clean:
                continue
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
                    clean_label = raw_label.split('_')[0] if '_' in raw_label else raw_label
                    translate_map[idx_token] = clean_label.lower()
                continue
            if line_clean.lower().startswith("tree ") or (line_clean.startswith("(") and line_clean.endswith(";")):
                tree_string = line_clean.split("=", 1)[1].strip() if "=" in line_clean else line_clean
                break

    if not tree_string:
        raise ValueError(f"Could not find a valid newick tree string in: {tree_path}")

    tree_string = re.sub(r'\[.*?\]', '', tree_string)
    lineage_paths_cache = {}

    try:
        target_tree = Phylo.read(io.StringIO(tree_string), "newick")
        for tip in target_tree.get_terminals():
            if not tip.name:
                continue
            tip_name = tip.name.strip().lower()
            resolved_key = translate_map.get(tip_name, tip_name).replace(' ', '_')
            lineage_paths_cache[resolved_key] = [target_tree.root] + target_tree.get_path(tip)
    except Exception:
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

    A = np.zeros((n, n))
    mode = str(tree_type).strip().lower()
    cached_paths = [lineage_paths_cache.get(sp, lineage_paths_cache.get(sp.replace(' ', '_'), [0])) for sp in species_list]

    for i in range(n):
        path_i = cached_paths[i]
        len_i = len(path_i)
        for j in range(i + 1, n):
            path_j = cached_paths[j]
            len_j = len(path_j)
            min_length = min(len_i, len_j)
            shared_depth = sum(1 for k in range(min_length) if path_i[k] == path_j[k])

            if mode == 'adaptive':
                ratio_threshold = min(max(depth_threshold * 0.10, 0.05), 0.95)
                if (shared_depth / len_i if len_i > 0 else 0) > ratio_threshold and (shared_depth / len_j if len_j > 0 else 0) > ratio_threshold:
                    A[i, j] = A[j, i] = 1.0
            else:
                if shared_depth > depth_threshold:
                    A[i, j] = A[j, i] = 1.0
    return A

def calculate_densities(data, id_col, tree=None, tree_type='fixed',
                        spatial_mode='auto', struct_mode='auto',
                        spatial_cats=None, struct_cats=None, coord_cols=None,
                        spatial_km=500.0, struct_depth=3, verbose=True):
    """
    Mode-driven density connectivity function. Accepts the following arguments:

    data                dataframe with observations [required]
    id_col              column for index identifiers (str) [required]
    tree                path to a Newick-style tree (str)
    tree_type           how to parse the tree, 'fixed' or 'adaptive', defaults to 'fixed' (str)
    spatial_mode        'coords', 'categorical', 'mixed', 'auto', or 'none' (str, will try to auto-detect)
    struct_mode         'tree', 'categorical', 'mixed', 'auto', or 'none' (str, will try to auto-detect)
    spatial_cats        categories for spatial density (str or list of str)
    struct_cats         categories for structural density (str or list of str)
    coord_cols          coordinate columns (2 or 3 columns containing Cartesian coordinates, float)
    spatial_km          distance (in km) to consider for proximity (float, defaults to 500.0)
    struct_depth        depth of tree connections to consider (int, defaults to 3)
    verbose             whether to print results to terminal (bool, defaults to True)

    """
    df = data.copy()
    n = len(df)
    spatial_density, structural_density = np.nan, np.nan
    if n < 2:
        return spatial_density, structural_density

    # 1. parameter type normalization
    s_cats = [spatial_cats] if isinstance(spatial_cats, str) else list(spatial_cats) if spatial_cats else []
    st_cats = [struct_cats] if isinstance(struct_cats, str) else list(struct_cats) if struct_cats else []

    # 2. automatic mode interpretation
    s_mode = str(spatial_mode).strip().lower()
    if s_mode == 'auto':
        if coord_cols and len(coord_cols) >= 2:
            s_mode = 'mixed' if s_cats else 'coords'
        else:
            s_mode = 'categorical' if s_cats else 'none'

    st_mode = str(struct_mode).strip().lower()
    if st_mode == 'auto':
        has_tree_file = tree is not None and os.path.exists(tree)
        if has_tree_file:
            st_mode = 'mixed' if st_cats else 'tree'
        else:
            st_mode = 'categorical' if st_cats else 'none'

    # 3. spatial density measure
    A_space = np.ones((n, n))
    calculated_space = False

    if s_mode in ['coords', 'mixed'] and coord_cols and len(coord_cols) >= 2:
        lat_rad = np.radians(df[coord_cols[0]].to_numpy().astype(float))
        lon_rad = np.radians(df[coord_cols[1]].to_numpy().astype(float))
        R = 6371.0
        coords = np.column_stack((R * np.cos(lat_rad) * np.cos(lon_rad), R * np.cos(lat_rad) * np.sin(lon_rad), R * np.sin(lat_rad)))
        A_space = (squareform(pdist(coords, metric='euclidean')) <= spatial_km).astype(float)
        calculated_space = True

    if s_mode in ['categorical', 'mixed'] and s_cats:
        A_space_cat = np.ones((n, n))
        valid_s_cols = [c for c in s_cats if c in df.columns]
        if valid_s_cols:
            for col in valid_s_cols:
                v = df[col].to_numpy()
                A_space_cat *= (v[:, None] == v[None, :]).astype(float)
            A_space = A_space * A_space_cat if calculated_space else A_space_cat
            calculated_space = True

    if calculated_space:
        spatial_density = _compute_edge_density(A_space)

    # 4. structural density measure
    A_tree = np.zeros((n, n))
    calculated_tree = False

    if st_mode in ['tree', 'mixed'] and tree and os.path.exists(tree):
        try:
            A_tree = _parse_unified_tree(tree, df, id_col, tree_type, struct_depth)
            calculated_tree = True
        except Exception as e:
            if verbose: print(f"Warning: Tree parsing failed: {e}")

    A_family = np.zeros((n, n))
    calculated_family = False
    if st_mode in ['categorical', 'mixed'] and st_cats:
        valid_st_cols = [c for c in st_cats if c in df.columns]
        if valid_st_cols:
            for col in valid_st_cols:
                v = df[col].to_numpy()
                A_family += (v[:, None] == v[None, :]).astype(float)
            A_family = (A_family > 0).astype(float)
            calculated_family = True

    if calculated_tree and calculated_family:
        structural_density = _compute_edge_density(((A_tree + A_family) > 0).astype(float))
    elif calculated_tree:
        structural_density = _compute_edge_density(A_tree)
    elif calculated_family:
        structural_density = _compute_edge_density(A_family)

    # data diagnostics printout
    if verbose:
        print("-" * 80)
        print(f"Spatial Mode       : {s_mode}")
        print(f"Structural Mode    : {st_mode}")
        print(f"Spatial density    : {spatial_density if pd.isna(spatial_density) else f'{spatial_density:.4f}'}")
        print(f"Structural density : {structural_density if pd.isna(structural_density) else f'{structural_density:.4f}'}")
        print("-" * 80)

    # evaluate spatial connectivity using continuous spatial boundary rules
    if not pd.isna(spatial_density):
        if spatial_density > 0.20:
            print("Spatial result     : High geographical clustering found.")
            print("Spatial rule       : High spatial autocorrelation risk. Spatial controls mandatory.")
        elif spatial_density <= 0.01:
            print("Spatial result     : Sparse geographical connections found.")
            print("Spatial rule       : High geographic dispersion. Standard models remain stable.")
        else:
            print("Spatial result     : Balanced, stable spatial baseline signal.")
        print("-" * 80)

    # evaluate structural connectivity using discrete tree boundary rules
    if not pd.isna(structural_density):
        if structural_density > 0.25:
            print("Structural result  : High structural data clustering found.")
            print("Structural rule    : Ideal for structure-aware models (PGLMM).")
        elif structural_density <= 0.10:
            print("Structural result  : Sparse Topology Danger Zone.")
            print("Structural rule    : High risk of parameter collapse. Continuous coordinate mapping (GPGLMM) required.")
        else:
            print("Structural result  : Balanced, stable structural baseline signal.")
            print("-" * 80)

    return spatial_density, structural_density
