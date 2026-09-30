# pykdensity

pykdensity is a Python package designed to measure connectivity (κ) within datasets used in biology, linguistics, and anthropology. Specifically, it calculates spatial density (how close your data points are to each other geographically) and structural density (how closely related your data points are evolutionarily or historically).

Measuring these attributes helps researchers identify when data points are too clustered, which can distort statistical results if left uncorrected.


## Installation

You can install this package directly from GitHub by adding it to your terminal or your project's requirements file.

```bash
pip install git+https://github.com/lingdoc/pykdensity
```


## Core Function Arguments

To use the package, import the `calculate_densities` function.

```python
from pykdensity import calculate_densities
```

It accepts the following arguments:

### Required Arguments
* **`data`** (pandas DataFrame): The table containing your research data.
* **`id_col`** (string): The name of the column that uniquely identifies each row (e.g., 'Species', 'glottocode', or 'Name').

### Optional Arguments
* **`tree`** (string): The file path to an evolutionary tree file (supports .nex, .csv, and compressed .trees.gz files).
* **`tree_type`** (string): Defines how structural node steps are calculated.
  * `"fixed"` (Default): Interprets the depth parameter as a flat integer count of branching events down from the master root node. Best for standard biological tree topologies.
  * `"adaptive"`: Dynamically converts the integer parameter into a percentage ratio (value * 10%) to normalize uneven branch lengths. Best for linguistic or high-variance tree lineages.
* **`spatial_mode`** (string): Controls spatial calculation. Accepts 'coords', 'categorical', 'mixed', 'auto', or 'none'.
* **`struct_mode`** (string): Controls structural calculation. Accepts 'tree', 'categorical', 'mixed', 'auto', or 'none'.
* **`spatial_cats`** (string or list of strings): Category columns used to group data by region names (e.g., 'macroarea' or 'Zone').
* **`struct_cats`** (string or list of strings): Category columns used to group data by flat lineages (e.g., 'Family_ID', 'Order', or 'Family').
* **`coord_cols`** (list of strings): Coordinates passed as a list of column headers containing points (e.g., `['latitude', 'longitude']` or `['X', 'Y', 'Z']`).
* **`spatial_km`** (float): The maximum distance in kilometers to consider two points geographically connected. Defaults to 500.0.
* **`struct_depth`** (integer): The minimum shared historical node depth or percentage ratio required to consider two points connected. Defaults to 3.
* **`verbose`** (boolean): Set to True to print progress logs and diagnostics to the terminal window. Defaults to True.


## Setting the Structural Depth Threshold

The `struct_depth` parameter controls how far back in historical or evolutionary time two observations must share a branch to be considered connected. Because biological and linguistic trees are formatted differently, the engine interprets this number in two distinct ways:

### 1. Biological Trees (Fixed Node Depth)
When parsing standard biological trees, this number represents a flat count of ancestral branching events starting from the root of the tree.
* **Low Values (e.g., 3 to 5):** Checks for broad, deep connections. Two species draw an edge if they simply belong to the same large taxonomic order or family.
* **High Values (e.g., 10+):** Restricts connections to highly specific, recent sub-clades. Only closely related sister species draw an edge.

### 2. Linguistic Trees (Adaptive Proportional Depth)
Linguistic trees often have wildly uneven branch lengths across language families. To prevent large, shallow language families from skewing your metrics, when `tree_type='adaptive'` the engine dynamically converts the integer value into a percentage threshold (multiplying the value by 10%).
* **Value of 3 (interprets as 30%):** A relaxed threshold. Two languages connect if they share even a minor portion of their historical path.
* **Value of 5 (interprets as 50%):** A balanced median benchmark. Requires languages to share at least half of their ancestral lineage history.
* **Value of 8 or higher (interprets as 80%+):** A restrictive threshold. Disconnects massive regional families from each other, only drawing an edge if the languages share a deep, specific local history.


## Examples

### 1. Using a Tree File (Biology Example)
If you have an explicit tree file showing how species are related:

```python
import pandas as pd
from pykdensity import calculate_densities

df = pd.read_csv("data/amphibians.csv")

spatial_k, structural_k = calculate_densities(
    data=df,
    id_col="Species",
    tree="data/amphibian100trees.nex",  
    tree_type="fixed",              
    struct_depth=5,
    verbose=True
)
```

### 2. Using Coordinates and Compressed Trees (Linguistics Example)
If you are tracking geographic points alongside a tree file:

```python
import pandas as pd
from pykdensity import calculate_densities

df = pd.read_csv("data/linguistics_dataset.csv")

spatial_k, structural_k = calculate_densities(
    data=df,
    id_col="glottocode",
    tree="tlu/0067_or_68KA/pruned_tree.trees.gz",
    tree_type="adaptive",
    coord_cols=["latitude", "longitude"],
    spatial_cats="macroarea",
    struct_cats="Family_ID",
    spatial_km=500.0,
    struct_depth=8,
    verbose=True
)
```

### 3. Using Categorical Groupings (Culture Example)
If your data does not use explicit tree paths but relies on regional and historical categories:

```python
import pandas as pd
from pykdensity import calculate_densities

df = pd.read_excel("data/Data_love_bywork_2019.xlsx")

spatial_k, structural_k = calculate_densities(
    data=df,
    id_col="Name",
    spatial_cats="Zone",
    struct_cats=["Century", "Genre.1"],
    verbose=True
)
```

## Advanced Multi-Layer Overlays (Mixed Mode)

A key feature of `pykdensity` is its ability to run **Mixed Mode** intersection metrics. This allows you to combine continuous topological layers (like geographic coordinates or branching tree structures) with flat categorical overlays (like macroareas or language family groupings) to get a more realistic picture of data clustering.

The engine handles these combined data spaces automatically using a logical intersection rule:

### 1. Spatial Mixed Overlays (`coord_cols` + `spatial_cats`)
When you provide both geographic coordinates and region category tags, the engine builds a dual-constraint network grid:
* **The Rule:** Two rows draw a connection **ONLY IF** they sit within the specified distance limit (e.g., `spatial_km=500.0`) **AND** belong to the exact same text category (e.g., `spatial_cats='macroarea'`).
* **Why use it:** This prevents distant geographic outliers that happen to share an artificial administrative boundary from bloating your connectivity scores.

### 2. Structural Mixed Overlays (`tree` + `struct_cats`)
When you pass a branching lineage file along with discrete flat family or order categories, the engine combines their signals:
* **The Rule:** Two rows draw a connection if they pass your ancestral branch path limit (`struct_depth`) **OR** if they match the categorical name labels exactly (e.g., `struct_cats='Family_ID'`).
* **Why use it:** This provides a crucial safety net for isolated linguistic or taxonomic data. If your tree file suffers from unanchored or missing deeper branches across family lines, the categorical overlay bridges the gap, allowing the engine to calculate a stable density score without dropping data.

### Example Configuration

To run a fully integrated, multi-layer mixed analysis, pass both continuous inputs and text category column strings simultaneously:

```python
import pandas as pd
from pykdensity import calculate_densities

df = pd.read_csv("data/advanced_dataset.csv")

# the engine auto-detects 'mixed' modes because both arrays are provided
spatial_k, structural_k = calculate_densities(
    data=df,
    id_col="glottocode",
    tree="data/pruned_tree.trees.gz",
    tree_type="adaptive",
    struct_depth=8,
    coord_cols=["latitude", "longitude"],
    spatial_km=500.0,
    spatial_cats="macroarea",     # intersects distance matrix with macro-regions
    struct_cats="Family_ID",      # unifies tree pathways with flat lineage IDs
    verbose=True
)
```

## Understanding the Outputs

The function returns two standard network adjacency metrics, both scaled strictly between `0.0` (no connectivity) and `1.0` (complete saturation):

* **Spatial Density:** The proportion of your data point pairs that sit close enough to each other geographically to cross your distance limit.
* **Structural Density:** The proportion of your data point pairs that share a deep historical or evolutionary branch segment.

### Matrix Interpretation Guidelines

Because these connectivity metrics reflect the strength of spatial and historical relationships in the data, they serve as a guide for selecting the right model architecture:

#### Spatial Density Tiers
* **Sparse Spatial Frontier (κ_spatial ≤ 0.01):**
  Indicates heavy geographic dispersion. While the continuous coordinate surface handles low density gracefully without collapsing your equations, it flags the need for specialized continuous models to fill extensive physical gaps. Standard models remain stable.
* **Dense Spatial Mesh (κ_spatial > 0.20):**
  Heavy geographic concentration. High risk of spatial autocorrelation (Galton's Problem), where proximity completely confounds historical independence. Spatial controls are mandatory.

#### Structural Density Tiers
* **Sparse Topology Danger Zone (κ_structural ≤ 0.10):**
  The structural matrix shatters into completely disconnected lineage islands. Traditional multi-level models (like `brms` unconstrained hierarchical setups) run out of shared branch traction, triggering severe parameter collapse or runaway errors. Continuous spatial cross-pooling architectures (like `GPGLMM` coordinate mapping) are required to anchor the parameters.
* **Dense Structural Mesh (κ_structural > 0.25):**
  The dataset contains a highly dense, robust historical signal. Ideal for structure-aware models like Phylogenetic GLMMs (PGLMMs). Simple regressions (like standard GLM) must be avoided, as high structural density violates row-independence assumptions, leading to artificially low p-values and high false-positive rates.

*Note: If your dataset lacks the columns needed for a specific calculation, the function returns `NaN` (Not a Number) for that metric instead of crashing, allowing your automation loops to finish running.*
