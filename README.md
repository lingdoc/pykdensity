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
* **`id_col`** (string): The name of the column that uniquely identifies each row (for example: 'Species', 'glottocode', or 'Name').

### Optional Arguments
* **`tree`** (string): The file path to an evolutionary tree file (supports .nex, .csv, and compressed .trees.gz files). If you provide a tree file, the package will use it to calculate structural relationships.
* **`tree_type`** (string): Defines how structural node steps are calculated.
  * `"fixed"` (Default): Interprets the threshold parameter as a flat integer count of branching events down from the master root node. Best for standard biological tree topologies.
  * `"adaptive"`: Dynamically converts the integer parameter into a percentage ratio (value * 10%) to normalize uneven branch lengths. Best for linguistic or high-variance tree lineages.
* **`taxonomy_hierarchy_cols`** (list of strings): If you do not have a tree file, you can pass a list of category columns (like `['Order', 'Family']`) to estimate structural relationships instead.
* **`coord_cols`** (list of strings): The geographic column names. If you pass two columns (like `['latitude', 'longitude']`), the package calculates real-world distances, converting Cartesian coordinates to 3D distance. If you pass a single column (like `['Zone']`), it groups items by matching region names.
* **`theme_cols`** (list of strings): Contextual categories (like `['Century', 'Genre']`) used to group cultural or historical data points.
* **`spatial_threshold_km`** (float): The maximum distance in kilometers to consider two points geographically connected. Defaults to 500.0.
* **`structural_depth_threshold`** (integer): The minimum number of shared historical steps required to consider two points structurally connected. Defaults to 3.
* **`verbose`** (boolean): Set to True to print progress updates to the terminal window. Defaults to True.


## Setting the Structural Depth Threshold

The `structural_depth_threshold` parameter controls how far back in historical or evolutionary time two observations must share a branch to be considered "connected." Because biological and linguistic trees are formatted differently, the engine interprets this number in two distinct ways:

### 1. Biological Trees (Fixed Node Depth)
When parsing standard biological trees (like `.nex` files), this number represents a flat count of ancestral branching events starting from the root of the tree.
* **Low Values (e.g., 3 to 5):** Checks for broad, deep connections. Two species will draw an edge if they simply belong to the same large taxonomic order or family.
* **High Values (e.g., 10+):** Restricts connections to highly specific, recent sub-clades. Only closely related sister species or animals in the exact same genus will draw an edge.

### 2. Linguistic Trees (Adaptive Proportional Depth)
Linguistic trees (like `.trees.gz` text strings) often have wildly uneven branch lengths across language families. To prevent large, shallow language families from skewing your metrics, when `tree_type='adaptive'` the engine dynamically converts the integer value into a percentage threshold (multiplying the value by 10%).
* **Value of 3 (interprets as 30%):** A relaxed threshold. Two languages draw a connection if they share even a minor portion of their historical path down from the root.
* **Value of 5 (interprets as 50%):** A balanced median benchmark. Requires languages to share at least half of their ancestral lineage history.
* **Value of 7 or higher (interprets as 70%+):** A more restrictive threshold. Disconnects massive regional families from each other, only drawing an edge if the languages share a deep, specific local history (like close dialects or sub-branches).


## Examples

### 1. Using a Tree File (Biology Example)
If you have an explicit tree file showing how species are related:

```python
import pandas as pd
from pykdensity import calculate_densities
# sample data from Munstermann et al (2022)
df = pd.read_csv("data/bio/mammals.csv")

# This will load/unzip the tree file and run the fixed-node biological calculation
spatial_k, structural_k = calculate_densities(
    data=df,
    id_col="Species",
    tree="data/bio/mammal100trees.nex",  
    tree_type="fixed",              
    structural_depth_threshold=5
)
```

### 2. Using Coordinates and Compressed Trees (Linguistics Example)
If you are tracking geographic points alongside a tree file:

```python
import pandas as pd
from pykdensity import calculate_densities
# sample data from Verkerk et al (2026)
df = pd.read_csv("data/tlu/Glottolog_Languages.csv")

# This reads a linguistic tracking tree formatted as a compressed .gz file
# using an adaptive tree method
spatial_k, structural_k = calculate_densities(
    data=df,
    id_col="glottocode",
    tree="data/tlu/0008KA/pruned_tree.trees.gz",
    tree_type="adaptive",
    structural_depth_threshold=8
)
```

### 3. Using Categorical Groupings (Culture Example)
If your data does not use explicit tree paths but relies on regional and historical categories:

```python
import pandas as pd
from pykdensity import calculate_densities
# sample data from Baumard et al (2022)
df = pd.read_excel("data/etc/Data_love_bywork_2019.xlsx")

# This loads the data and reads the categorical values
spatial_k, structural_k = calculate_densities(
    data=df,
    id_col="Name",
    coord_cols=["Zone"],
    theme_cols=["Century", "Genre"]
)
```

### Additional notes

More details can be found in the `calculate_k_densities.py` script in this repo, which runs the calculation over 3 example evolutionary datasets from biology, linguistics, and culture.


## Understanding the Outputs

The function returns two standard network adjacency metrics, both scaled strictly between `0.0` (no connectivity) and `1.0` (complete saturation):

* **Spatial Density:** The proportion of your data point pairs that sit close enough to each other geographically to cross your distance limit.
* **Structural Density:** The proportion of your data point pairs that share a deep historical or evolutionary branch segment.

### Matrix Interpretation Guidelines

Because these connectivity metrics reflect the strength of spatial and historical relationships in the data, they serve as a guide for selecting the right model architecture:

* **Sparse Matrix / Signal Deficit (κ ≤ 0.005):**
  The network matrix is extremely sparse (as seen in the linguistics track). When connectivity drops this low, standard variance-partitioning frameworks (like PGLMM or Gaussian Process tracks) often struggle or fail because there is almost no shared historical overlap between data points.
  * **Frequentist Models:** Typically flag this deficit by throwing optimization warnings, boundary constraints, or failing to converge entirely.
  * **Bayesian Models:** May successfully complete sampling chains and report technical convergence (stable trace plots and clean R-hat diagnostics) due to the smoothing influence of regularizing priors. However, the model may still suffer from hidden parameter explosion or structural variance collapse—where the posterior distribution simply mirrors the prior because the data signal is too weak. This calls the validity of the partitioned variance results into question.

  In this sparse tier, simpler flat regressions or strict categorical controls are often more stable, and in some cases the data may require special treatment.

* **Moderate Structural Signal (0.005 < κ ≤ 0.15):**
  The network possesses a mild, balanced signal (as seen in the culture track). There is enough shared historical overlap to separate background lineage history from your primary variables without overwhelming the model. Standard regressions with basic regional or family random effects usually perform well here.

* **Dense Matrix / Strong Covariance Signal (κ > 0.15):**
  The dataset contains a highly dense, robust historical or spatial signal (as seen in the biology track).
  * **Standard Architectures:** Simple regressions (like standard GLM) should be avoided here, as the high density violates basic row-independence assumptions, leading to artificially low p-values and high false-positive rates.
  * **Advanced Architectures:** This tier is **ideal for specialized, structure-aware models** like Phylogenetic GLMMs, continuous Gaussian Processes, or spatial autoregressive workflows. The high density provides a rich, strong signal that allows these architectures to map and control for background historical relationships.

*Note: If your dataset lacks the columns needed for a specific calculation, the function returns `NaN` (Not a Number) for that metric instead of crashing, allowing your automation loops to finish running.*
