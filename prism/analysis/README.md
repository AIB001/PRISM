# PRISM Analysis Module

A comprehensive toolkit for MD trajectory analysis with a focus on protein-ligand interactions.

---

## ⭐ Highlight: Interactive HTML Contact Visualization

**The most powerful feature of this module is the interactive HTML contact visualization generator** (`contact/htmlgen`), which creates beautiful, publication-ready visualizations of protein-ligand interactions with zero configuration.

### Quick Start - Generate Interactive HTML

```bash
# Command-line usage
python -m prism.analysis.contact.htmlgen trajectory.xtc topology.pdb ligand.mol2 -o contacts.html

# Python API
from prism.analysis.contact import generate_html
generate_html("trajectory.xtc", "topology.pdb", "ligand.mol2", "output.html")
```

**What you get:**
- 🎨 **2D/3D interactive visualization** - Switch between wheel layout and molecular structure
- 📊 **Color-coded contact frequencies** - Instantly identify key interactions (0-100% scale)
- 🖱️ **Draggable interface** - Rearrange residues and adjust layout in real-time
- 🔄 **Flexible display modes** - Choose between unique or duplicate residue modes
- 🎛️ **Customizable view** - Toggle grid lines, rotate 180°, hide/show elements
- 📸 **High-resolution export** - Export publication-quality PNG images (up to 8K)
- 📈 **Detailed statistics** - Contact frequencies, atom-level details, distance measurements
- 🔍 **Interactive tooltips** - Hover over atoms/residues for detailed information

**See the [Contact Analysis](#contact-analysis) section below for complete documentation.**

---

## Table of Contents

1. [Contact Analysis](#contact-analysis) ⭐ **Most Important**
   - [Interactive HTML Visualization](#interactive-html-visualization)
   - [Contact Detection Methods](#contact-detection-methods)
   - [Contact Frequency Analysis](#contact-frequency-analysis)
2. [Structural Analysis](#structural-analysis)
   - [RMSD (Root Mean Square Deviation)](#rmsd)
   - [RMSF (Root Mean Square Fluctuation)](#rmsf)
3. [Clustering Analysis](#clustering-analysis)
   - [K-means Clustering](#k-means-clustering)
   - [DBSCAN Clustering](#dbscan-clustering)
   - [Hierarchical Clustering](#hierarchical-clustering)
4. [Hydrogen Bond Analysis](#hydrogen-bond-analysis)
5. [Dihedral Angle Analysis](#dihedral-angle-analysis)
6. [Visualization Guidelines](#visualization-guidelines)
7. [Performance Optimization](#performance-optimization)

---

## Contact Analysis

### Interactive HTML Visualization

**The flagship feature of PRISM analysis module** - generates self-contained interactive HTML files for exploring protein-ligand contacts.

#### Features

##### Core Capabilities
- **Automatic ligand detection** - No manual specification needed
- **Fast vectorized analysis** - Handles trajectories with >1000 frames efficiently
- **Smart frame sampling** - Automatically samples large trajectories for performance
- **Two contact display modes** - Unique residue or allow duplicates for detailed analysis
- **Multiple visualization modes** - 2D wheel layout and 3D molecular structure
- **Real-time interactivity** - Drag, zoom, pan, rotate, and rearrange elements
- **Customizable display** - Toggle grid lines, hide/show connections and hydrogens

##### Visual Elements
- **Color-coded atoms** - Carbon (green), Oxygen (red), Nitrogen (blue), Sulfur (yellow)
- **Bond visualization** - Complete ligand chemical structure
- **Contact lines** - Thickness represents interaction strength
- **Frequency coloring** - Gradient from low (blue) to high (red) contact frequency

##### Analysis Output
- **Contact statistics** - Total contacts, high-frequency interactions (>50%), maximum frequency
- **Atom-level details** - Which ligand atoms interact with which residues
- **Distance measurements** - Average contact distances for each interaction
- **Automatic filtering** - Only shows significant contacts (>10% frequency)

#### Contact Display Modes

PRISM supports two contact visualization modes to handle different analysis scenarios:

##### Mode 1: Unique Residue Mode (Default)
**Best for:** Quick overview of main contacting residues

- Each residue appears **only once** in the visualization
- Shows the **best contact** (highest frequency) per residue
- Default: Top 20 contacts displayed
- Based on `residue_proportions` metric

**Example:** If ARG123 contacts both ligand atom C1 (80%) and C5 (75%), only the C1 contact is shown.

##### Mode 2: Allow Duplicate Residues Mode
**Best for:** Detailed analysis of multi-point interactions

- Same residue can appear **multiple times** with different ligand atoms
- Shows all **high-frequency atom-pair contacts**
- Default: Top 25 contacts displayed
- Based on `contact_frequencies` metric (atom-pair level)

**Example:** ARG123 appears twice - once for ARG123-C1 (80%) and once for ARG123-C5 (75%).

##### Comparison Table

| Feature | Unique Mode | Duplicate Mode |
|---------|-------------|----------------|
| **Parameter** | `allow_duplicate_residues=False` | `allow_duplicate_residues=True` |
| **Residue repetition** | One per residue | Multiple allowed |
| **Contact count** | Top 20 (default) | Top 25 (default) |
| **Data source** | Residue-level aggregation | Atom-pair level |
| **Use case** | Overview of key residues | Detailed multi-point contacts |

#### Command-Line Interface

```bash
# Mode 1: Unique residue mode (default)
python -m prism.analysis.contact.htmlgen \
    trajectory.xtc \
    topology.pdb \
    ligand.mol2 \
    -o contact_analysis.html

# Mode 2: Allow duplicate residues
python -m prism.analysis.contact.htmlgen \
    trajectory.xtc \
    topology.pdb \
    ligand.mol2 \
    -o contact_analysis.html \
    --allow-duplicates

# Custom maximum contacts
python -m prism.analysis.contact.htmlgen \
    trajectory.xtc \
    topology.pdb \
    ligand.mol2 \
    --allow-duplicates \
    --max-contacts 30
```

#### Python API

##### Quick Generation
```python
from prism.analysis.contact import generate_html

# Mode 1: Unique residue mode (default)
generate_html(
    trajectory="trajectory.xtc",
    topology="topology.pdb",
    ligand="ligand.mol2",
    output="contact_analysis.html"
)

# Mode 2: Allow duplicate residues
generate_html(
    trajectory="trajectory.xtc",
    topology="topology.pdb",
    ligand="ligand.mol2",
    output="contact_analysis.html",
    allow_duplicate_residues=True,
    max_contacts=25
)
```

##### Advanced Usage
```python
from prism.analysis.contact import HTMLGenerator

# Mode 1: Create generator with unique residue mode
generator = HTMLGenerator(
    trajectory_file="trajectory.xtc",
    topology_file="topology.pdb",
    ligand_file="ligand.mol2",
    allow_duplicate_residues=False,  # Default
    max_contacts=20                   # Default for unique mode
)

# Mode 2: Allow duplicate residues for detailed analysis
generator = HTMLGenerator(
    trajectory_file="trajectory.xtc",
    topology_file="topology.pdb",
    ligand_file="ligand.mol2",
    allow_duplicate_residues=True,
    max_contacts=25                   # Default for duplicate mode
)

# Run analysis
results = generator.analyze()
print(f"Total frames: {results['stats']['total_frames']}")
print(f"Total contacts: {results['stats']['total_contacts']}")
print(f"High-frequency contacts: {results['stats']['high_freq_contacts']}")

# Generate HTML
output_file = generator.generate(output_file="my_contacts.html")
print(f"Generated: {output_file}")
```

##### Integration with Analysis Workflow
```python
from prism.analysis.core import TrajectoryAnalyzer
from prism.analysis.contact import generate_html

# Run full PRISM analysis
analyzer = TrajectoryAnalyzer(...)
analyzer.run()

# Generate contact HTML
generate_html(
    trajectory=analyzer.trajectory_file,
    topology=analyzer.topology_file,
    ligand=analyzer.ligand_file
)
```

#### Configuration

##### Contact Detection Parameters
```python
class Config:
    contact_enter_threshold_nm = 0.35  # Distance to form contact (3.5 Å)
    contact_exit_threshold_nm = 0.4    # Distance to break contact (4.0 Å)
    distance_cutoff_nm = 0.5           # Maximum distance to consider (5.0 Å)
    min_frames_for_smoothing = 10      # Minimum frames for smoothing
    smooth_window = 11                 # Smoothing window size

# Use custom configuration
from prism.analysis.contact import FastContactAnalyzer
analyzer = FastContactAnalyzer(traj, config=custom_config)
```

##### Filtering Thresholds
- **Atom-residue contacts**: >10% frequency (line 148 in contact_analyzer.py)
- **Residue-level contacts**: >2% proportion (line 176 in contact_analyzer.py)

#### Module Components

##### `htmlgen.py` - HTML Generator
Main entry point for generating interactive visualizations.

**Key Classes:**
- `HTMLGenerator`: Main class for analysis and HTML generation
- `generate_html()`: Convenience function for one-line generation

##### `contact_analyzer.py` - FastContactAnalyzer
Fast vectorized contact analysis engine using MDTraj.

**Features:**
- Automatic ligand identification using PRISM utilities
- Efficient distance calculations with OpenMP parallelization
- Smart frame sampling for large trajectories (auto-sample to 1000 frames)
- Contact frequency and proportion calculations

**Usage:**
```python
from prism.analysis.contact import FastContactAnalyzer
import mdtraj as md

traj = md.load("trajectory.xtc", top="topology.pdb")
analyzer = FastContactAnalyzer(traj)
results = analyzer.calculate_contact_proportions()

print(f"Found {len(results['residue_proportions'])} contacting residues")
```

##### `data_processor.py` - DataProcessor
Process analysis results and ligand structures.

**Capabilities:**
- Load ligand files (.mol2, .sdf, .mol)
- Extract chemical structure and bonding information
- Calculate contact statistics
- Convert data for visualization

##### `visualization_generator.py` - VisualizationGenerator
Generate visualization data structures.

**Functions:**
- Create ligand atom data (positions, elements, bonds)
- Process contact data for rendering
- Calculate optimal layout coordinates

##### `html_builder.py` - HTMLBuilder
Build complete HTML with embedded JavaScript.

**Output:**
- Complete HTML5 document
- Embedded JSON data
- Included CSS styling and JavaScript code
- Support for high-DPI displays (Retina-ready)

##### `javascript_code.py` - JavaScript Engine
Complete visualization engine embedded in HTML.

**Features:**
- Canvas-based 2D/3D rendering
- Mouse interaction handlers (drag, zoom, pan)
- Layout algorithms (wheel, 3D coordinates)
- Animation and smooth transitions
- High-resolution image export functionality

#### Complete Example

```python
#!/usr/bin/env python3
from prism.analysis.contact import HTMLGenerator

# File paths
trajectory = "/path/to/md_trajectory.xtc"
topology = "/path/to/system.pdb"
ligand = "/path/to/ligand.mol2"

# Create generator
generator = HTMLGenerator(trajectory, topology, ligand)

# Run analysis (may take a few minutes for large trajectories)
print("Analyzing contacts...")
results = generator.analyze()

# Print summary
print(f"\nAnalysis Summary:")
print(f"  Total frames: {results['stats']['total_frames']}")
print(f"  Total contacts: {results['stats']['total_contacts']}")
print(f"  High-frequency contacts: {results['stats']['high_freq_contacts']}")
print(f"  Max frequency: {results['stats']['max_freq_percent']}%")

# Generate interactive HTML
print("\nGenerating HTML visualization...")
output_file = generator.generate("protein_ligand_contacts.html")
print(f"\n✓ Complete! Open {output_file} in a browser to explore.")
```

#### Use Cases

1. **MD Trajectory Analysis**
   - Visualize contact evolution over simulation time
   - Identify persistent vs transient interactions
   - Discover key binding site residues

2. **Binding Site Characterization**
   - Map residues involved in ligand recognition
   - Quantify interaction frequencies
   - Identify hotspot residues for mutagenesis studies

3. **Drug Design**
   - See which ligand atoms make important contacts
   - Guide chemical modifications
   - Compare different ligand analogs

4. **Publication Figures**
   - Export high-resolution images (up to 8K)
   - Interactive supplementary materials
   - Share with collaborators (single HTML file)

#### Performance Tips

**Large Trajectories (>1000 frames):**
- Automatic sampling to 1000 frames
- Maintains statistical accuracy
- Significantly faster analysis

**Memory Management:**
- Uses vectorized numpy operations
- Efficient MDTraj distance calculations
- Processes all atom pairs simultaneously

**Output File Size:**
- Typical HTML file: 500KB-2MB
- Self-contained (no external dependencies)
- Easy sharing (single file)
- Offline viewing (no internet required)

#### Troubleshooting

**Problem: Ligand not detected automatically**
```python
from prism.utils.ligand import identify_ligand_residue
import mdtraj as md

traj = md.load("trajectory.xtc", top="topology.pdb")
ligand_res = identify_ligand_residue(traj)
if not ligand_res:
    print("Ligand detection failed - check residue naming")
    print("Ensure ligand has a unique 3-letter code (not a standard amino acid)")
```

**Problem: HTML file too large**
- HTML embeds all data and JavaScript
- This is intentional for portability
- Typical size: 500KB-2MB (acceptable)
- For very large systems, consider frame sampling

**Problem: Slow analysis**
- Enable frame sampling (automatic for >1000 frames)
- Use `step` parameter to skip frames
- Check available CPU cores for parallelization

---

### Contact Detection Methods

#### Distance Calculation

PRISM contact analysis is based on **atom-pair distances**, not center-of-mass distances.

**Method:**
1. **Heavy atom selection:**
   - Ligand: All non-hydrogen atoms
   - Protein residues: Backbone and sidechain C, N, O, S atoms

2. **Distance calculation:**
   - For each (ligand atom, protein atom) pair, compute spatial distance
   - Uses MDTraj's `compute_distances()` with OpenMP optimization
   - **Closest distance**: For a given residue, take minimum distance across all atom pairs

3. **Contact criteria:**
   - Threshold: `contact_enter_threshold_nm` (default 0.35 nm = 3.5 Å)
   - Contact formed if distance < threshold

**Important Notes:**
- ✅ Uses actual atomic distances (not center-of-mass)
- ✅ Evaluates all atom pairs independently
- ✅ Takes minimum distance per residue

#### Contact Frequency Analysis

**1. Atom-Residue Contact Frequency (`contact_frequencies`)**

For each (ligand atom, residue) pair, calculate the proportion of frames where contact is formed:

$$
\text{Frequency} = \frac{\text{Number of frames with contact}}{\text{Total frames}}
$$

**2. Residue Contact Proportion (`residue_proportions`)**

For each residue, account for the number of ligand atoms it contacts:

$$
\text{Proportion} = \frac{\text{Total contacts}}{\text{Total frames} \times \text{Number of contacting ligand atoms}}
$$

This metric normalizes for residue size effects.

**3. Average Contact Distance (`residue_avg_distances`)**

For each residue, calculate mean distance across all contact frames:
- Only includes frames where distance < threshold
- Reports average distance in Ångströms

**Output Data Structure:**
```python
results = {
    'contact_frequencies': {
        (ligand_atom_idx, 'RES123'): 0.85,  # 85% contact frequency
        (ligand_atom_idx, 'RES456'): 0.60   # 60% contact frequency
    },
    'residue_proportions': {
        'RES123': 0.75,  # Residue contact proportion
        'RES456': 0.60
    },
    'residue_avg_distances': {
        'RES123': 3.2,  # Average distance in Å
        'RES456': 3.8
    },
    'residue_best_ligand_atoms': {
        'RES123': 5,  # Ligand atom with most contacts
        'RES456': 12
    },
    'ligand_residue': <MDTraj Residue>,  # Ligand residue object
    'total_frames': 1000                  # Total analyzed frames
}
```

**Visualization Options:**
- Bar charts: Key residue contact proportions
- Heatmaps: Ligand atom vs residue contact matrix
- Interactive HTML: 3D contact network visualization

#### Contact Number Timeseries

Track the number of contacts over simulation time.

**Method:**
1. Calculate distance matrix for all (ligand, protein) atom pairs
2. For each frame, count atom pairs with distance < threshold
3. Output timeseries data

**Output:**
```python
results = {
    'contact_numbers': np.array([45, 48, 52, ...]),  # Contacts per frame
    'times': np.array([0, 0.5, 1.0, ...]),           # Time points (ns)
    'total_pairs': 12500                              # Total atom pairs
}
```

**Visualization:**
- Timeseries plot: Contact number vs time
- Moving average overlay
- Statistical bands (mean ± std)

#### Residue-Ligand Distance Analysis

**Method 1: Contact Distance Distribution** (`analyze_contact_distances`)

- **Purpose**: Analyze distance distribution when contacts are formed
- **Method**: Calculate minimum distance between residue and ligand per frame
- **Filter**: No filtering (includes non-contact frames)
- **Use case**: Study dynamic contact/non-contact transitions

**Method 2: Distance Timeseries** (`analyze_residue_distance_timeseries`)

- **Purpose**: Track specific residue-ligand distance over time
- **Method**: Minimum distance per frame
- **Output**: Complete timeseries for each residue

**Output:**
```python
# Distance distribution
distances = {
    'ASP618': np.array([3.2, 3.5, 3.1, ...]),  # Distance array (Å)
    'ARG555': np.array([3.8, 4.2, 3.9, ...])
}

# Distance timeseries
timeseries = {
    'ASP618': np.array([3.2, 3.5, 3.1, ...]),  # Distance timeseries (Å)
    'times': np.array([0, 0.5, 1.0, ...])       # Time points (ns)
}
```

**Visualizations:**
- Histogram: Distance probability distribution
- Timeseries: Distance evolution over simulation
- Boxplot: Compare distance distributions across residues

---

## Structural Analysis

### RMSD

**Root Mean Square Deviation** - Measures structural deviation from reference.

#### Calculation Principle

$$
\text{RMSD} = \sqrt{\frac{1}{N} \sum_{i=1}^{N} (r_i - r_{0i})^2}
$$

Where:
- $N$ = number of atoms
- $r_i$ = coordinates of atom $i$ in current frame
- $r_{0i}$ = coordinates of atom $i$ in reference structure

#### Key Parameters

- **`align_selection`**: Atoms used for structural alignment (e.g., `"protein and name CA"`)
- **`calculate_selection`**: Atoms for RMSD calculation (can differ from alignment)
- **`ref_frame`**: Reference structure frame index (default: 0)
- **`step`**: Frame interval for analysis (speeds up calculation)

#### Calculation Workflow

1. Select reference structure (typically first frame or experimental structure)
2. Align trajectory to reference using `align_selection`
3. Calculate RMSD for `calculate_selection` atoms
4. Output RMSD value per frame (units: Å)

#### Typical Usage

```python
# Protein backbone RMSD (align and calculate on Cα)
rmsd = analyzer.calculate_rmsd(
    universe=topology_file,
    trajectory=trajectory_file,
    align_selection="protein and name CA",
    calculate_selection="protein and name CA"
)

# Ligand RMSD (align protein, calculate on ligand)
rmsd_ligand = analyzer.calculate_rmsd(
    universe=topology_file,
    trajectory=trajectory_file,
    align_selection="protein and name CA",  # Align protein first
    calculate_selection="resname LIG"       # Then calculate ligand RMSD
)
```

#### Visualization

**RMSD Timeseries Plot:**
- X-axis: Time (ns)
- Y-axis: RMSD (Å)
- Optional: Smoothed curve (Savitzky-Golay filter)
- Helps identify equilibration and conformational transitions

---

### RMSF

**Root Mean Square Fluctuation** - Measures atom/residue flexibility.

#### Calculation Principle

$$
\text{RMSF}_i = \sqrt{\frac{1}{T} \sum_{t=1}^{T} (r_i(t) - \bar{r}_i)^2}
$$

Where:
- $T$ = total frames
- $r_i(t)$ = coordinates of atom $i$ at time $t$
- $\bar{r}_i$ = average coordinates of atom $i$ across trajectory

#### Key Parameters

- **`align_selection`**: Atoms for trajectory alignment (ensures proper comparison)
- **`calculate_selection`**: Atoms for RMSF calculation
- **`auto_detect_chains`**: Automatically detect all protein/nucleic acid chains

#### Calculation Workflow

1. Align trajectory to average structure (computed automatically)
2. Calculate root mean square fluctuation per atom
3. For Cα/P atoms: One RMSF value per residue/nucleotide
4. Output RMSF values with corresponding residue numbers

#### Automatic Chain Detection

- **Protein chains**: Uses Cα atoms
- **Nucleic acid chains**: Uses P atoms (phosphate backbone)
- Generates independent RMSF data for each chain automatically

#### Visualization

**RMSF vs Residue Number Plot:**
- X-axis: Residue number
- Y-axis: RMSF (Å)
- Different chains shown in different colors/panels
- Highlights flexible regions (loops, termini) vs rigid regions (secondary structures)

---

## Clustering Analysis

Clustering identifies structurally similar conformations (states) in MD trajectories. PRISM supports three major clustering algorithms.

### K-means Clustering

#### Method Principle

K-means partitions conformations into k clusters by minimizing within-cluster variance.

#### Algorithm Steps

1. **Coordinate extraction**: Extract clustering atom coordinates (typically Cα or protein backbone)
2. **Standardization**: Normalize coordinates (zero mean, unit variance)
3. **Dimensionality reduction** (optional): Apply PCA (recommended)
4. **Clustering**: Assign frames to k clusters using K-means algorithm
5. **Centroid calculation**: Identify cluster centers (representative structures)

#### Key Parameters

- **`n_clusters`**: Number of clusters (must be specified)
- **`use_pca`**: Whether to use PCA dimensionality reduction (recommended: True)
- **`n_components`**: Number of PCA components (default: 10)
- **`align_selection`**: Atoms for alignment (typically Cα)
- **`cluster_selection`**: Atoms used for clustering (can be full protein or Cα)

#### Evaluation Metrics

- **Silhouette Score**: Range [-1, 1], closer to 1 indicates better clustering
- **Inertia**: Within-cluster sum of squares, lower is better

#### Output Data

```python
results = {
    'labels': np.array([0, 1, 1, 0, 2, ...]),        # Cluster label per frame
    'cluster_centers': np.array([[...], [...]]),     # Cluster center coordinates
    'silhouette_score': 0.65,                        # Silhouette coefficient
    'n_clusters': 5,                                 # Number of clusters
    'inertia': 1234.5                                # Inertia value
}
```

#### Visualizations

- **PCA projection scatter plot**: Clusters in first two principal component space
- **Cluster timeseries**: Which cluster each frame belongs to (state transitions)
- **Cluster size pie chart**: Proportion of frames in each cluster
- **Representative structures**: Cluster centroid structures saved as PDB files

---

### DBSCAN Clustering

#### Method Principle

**DBSCAN** (Density-Based Spatial Clustering of Applications with Noise) is a density-based method that:
- Automatically determines number of clusters
- Identifies noise points (outlier frames)
- Discovers arbitrary-shaped clusters

#### Algorithm Steps

1. Define neighborhood: Radius = `eps`
2. Core points: At least `min_samples` points within neighborhood
3. Border points: Within neighborhood of core point but not core itself
4. Noise points: Neither core nor border points

#### Key Parameters

- **`eps`**: Neighborhood radius
  - If `use_pca=True`: eps is distance in PCA space (dimensionless)
  - If `use_pca=False`: eps is RMSD distance (Å)
- **`min_samples`**: Minimum neighbors for core points (typically 5-10)
- **`use_pca`**: Whether to use PCA (recommended: True)

#### Parameter Selection Guidelines

- **Small eps** (0.3-0.5): Produces more small clusters
- **Large eps** (1.0-2.0): Produces fewer large clusters
- **`min_samples`**: Approximately $\ln(N)$, where $N$ is total frames

#### Output Data

```python
results = {
    'labels': np.array([0, 1, -1, 0, 2, ...]),   # -1 indicates noise
    'n_clusters': 4,                              # Auto-determined cluster count
    'n_noise': 15,                                # Number of noise points
    'silhouette_score': 0.58                      # Silhouette (excluding noise)
}
```

#### Visualizations

- **DBSCAN scatter plot**: Noise points marked with special color
- **Cluster statistics bar chart**: Frames per cluster and noise points

---

### Hierarchical Clustering

#### Method Principle

**Agglomerative Clustering** builds a hierarchical tree structure:
- **Agglomerative**: Start with each point as a cluster, progressively merge
- **Cuttable**: Cut tree at any level to obtain different numbers of clusters

#### Algorithm Steps

1. Each frame starts as individual cluster
2. Iteratively merge two closest clusters
3. Cut at specified level to obtain final clusters

#### Key Parameters

- **`n_clusters`**: Final number of clusters
- **`linkage`**: Inter-cluster distance calculation method
  - `ward`: Minimize within-cluster variance (recommended)
  - `complete`: Maximum distance
  - `average`: Average distance

#### Visualizations

- **Dendrogram**: Shows hierarchical structure
- **Heatmap**: Frame-to-frame RMSD matrix

---

### Finding Optimal Number of Clusters

PRISM provides `find_optimal_clusters()` to automatically determine optimal cluster count.

#### Evaluation Metrics

1. **Silhouette Score**: Measures cluster compactness and separation
2. **Elbow Method** (K-means only): Inertia vs cluster count curve

#### Usage

```python
results = analyzer.find_optimal_clusters(
    universe=topology_file,
    trajectory=trajectory_file,
    max_clusters=10,        # Test 2-10 clusters
    method="kmeans"
)
optimal_k = results['optimal_clusters']  # Recommended cluster count
```

#### Visualizations

- **Silhouette score curve**: Scores for different cluster counts
- **Elbow plot**: Inertia vs cluster count (for K-means)

---

## Hydrogen Bond Analysis

### Calculation Principle

Hydrogen bonds identified using geometric criteria:
- **Distance criterion**: Donor heavy atom - Acceptor heavy atom distance < 3.5 Å
- **Angle criterion**: D-H···A angle > 120°

### Analysis Content

1. Identify protein-ligand hydrogen bonds
2. Calculate hydrogen bond formation frequencies
3. Track hydrogen bonds over time

### Output Data

- List of hydrogen bond pairs with formation frequencies
- Hydrogen bond count timeseries
- Statistics for key residues

### Visualizations

- Hydrogen bond frequency bar chart
- Hydrogen bond network diagram
- Hydrogen bond count timeseries

---

## Dihedral Angle Analysis

### Analysis Types

#### 1. Backbone Dihedrals ($\phi$, $\psi$)

- $\phi$ (phi): C(-1) - N - Cα - C
- $\psi$ (psi): N - Cα - C - N(+1)
- Used for Ramachandran plot analysis

#### 2. Sidechain Dihedrals ($\chi_1$, $\chi_2$, ...)

- $\chi_1$: N - Cα - Cβ - Cγ
- $\chi_2$: Cα - Cβ - Cγ - Cδ
- Used for sidechain conformation analysis

### Visualizations

- **Ramachandran plot**: $\phi$-$\psi$ distribution
- **Dihedral timeseries**: Track conformational changes of specific residues
- **Dihedral distribution histograms**: Statistical conformational preferences

---

## Visualization Guidelines

### Universal Standards

**Font and Style:**
- **Font**: Times New Roman (global standard)
- **Font size**: Large enough for use as subplots
- **Resolution**: 300 DPI (publication quality)
- **File formats**: PNG (default), PDF (vector graphics)

**File Naming Convention:**
- Use descriptive names: `protein_rmsd_timeseries.png`
- Avoid generic names: `plot.png`

### Major Plot Types

#### Timeseries Plots
- **Usage**: RMSD, contact count, hydrogen bond count vs time
- **Features**: Optional smoothing curves, moving averages

#### Per-Residue Plots
- **Usage**: RMSF, B-factor, contact frequency vs residue number
- **Features**: Key residue annotations, region highlighting

#### Scatter Plots
- **Usage**: PCA projections for clustering
- **Features**: Different clusters in different colors

#### Heatmaps
- **Usage**: Contact matrices, RMSD matrices, distance matrices
- **Features**: Color mapping, colorbars

#### Distribution Plots
- **Usage**: Histograms, boxplots, violin plots
- **Features**: Statistical distribution characteristics

#### Bar Charts
- **Usage**: Contact frequencies, cluster sizes, hydrogen bond statistics
- **Features**: Sortable by value, highlighting

#### Multi-Panel Figures
- **Usage**: Comprehensive display of related analyses
- **Features**: Unified styling, shared axes

### Interactive Visualizations

**HTML Output:**
- Generated using Plotly for interactive plots
- Support zoom, hover information, export
- Suitable for exploratory analysis

**3D Visualizations:**
- Molecular structures overlaid with analysis data
- 3D contact network visualization

---

## Performance Optimization

### Parallel Computing

**Automatic Parallelization:**
- PRISM automatically detects available CPU cores
- Parallelizes distance calculations, RMSD operations via OpenMP
- No manual configuration required

**Manual Configuration:**
```python
import os
os.environ['OMP_NUM_THREADS'] = '8'  # Manually specify thread count
```

### Caching Mechanism

**Automatic Caching:**
- All calculation results automatically cached to `cache/` directory
- Uses pickle serialization for storage
- Cache keys automatically generated from parameters and filenames

**Cache Usage:**
```python
# First run: Execute calculation and cache
rmsd = analyzer.calculate_rmsd(...)

# Second run: Load from cache (nearly instant)
rmsd = analyzer.calculate_rmsd(...)  # Same parameters
```

**Clear Cache:**
```bash
rm -rf cache/  # Delete all cached data
```

### Memory Optimization

**Large Trajectory Processing:**
- Use `step` parameter for frame skipping
- Limit with `start_frame` and `end_frame`
- Automatic streaming processing

---

## Practical Guidelines

### Parameter Selection

**RMSD/RMSF Analysis:**
- Alignment selection: Use stable regions (e.g., Cα)
- Calculation selection: Choose based on research objective (ligand, binding site, etc.)

**Contact Analysis:**
- Threshold: 4.0-4.5 Å is commonly used
- Key residues: Use `analyze_key_residue_contacts()` to auto-filter top N

**Clustering Analysis:**
- K-means: First use `find_optimal_clusters()` to determine cluster count
- DBSCAN: Start with small `eps`, observe cluster count and noise ratio
- PCA dimensionality reduction significantly speeds up calculation

### Analysis Workflow

**Typical Workflow:**
1. **Preprocessing**: PBC correction, trajectory alignment
2. **Initial analysis**: RMSD, RMSF for overall stability assessment
3. **Detailed analysis**: Contact, hydrogen bond analysis for binding site
4. **Conformational analysis**: Clustering to identify major states
5. **Visualization**: Generate publication-quality figures

### Common Issues

**Q: RMSD continuously increasing?**
- Check if more stable alignment selection is needed
- May need to remove overall translation/rotation

**Q: Poor clustering results?**
- Try adjusting cluster count or algorithm parameters
- Check if PCA explained variance is sufficient (recommend >80%)
- May need longer trajectory or better sampling

**Q: How to choose contact analysis threshold?**
- 4.0-4.5 Å is common range (first coordination shell)
- Adjust based on specific system (larger ligands may use larger cutoff)

---

## API Reference

### Contact Analysis

#### `generate_html(trajectory, topology, ligand, output, allow_duplicate_residues, max_contacts)`
Quick HTML generation function.

**Parameters:**
- `trajectory` (str): Trajectory file path (.xtc, .dcd, etc.)
- `topology` (str): Topology file path (.pdb, .gro, etc.)
- `ligand` (str): Ligand structure file path (.mol2, .sdf, .mol)
- `output` (str): Output HTML file path (default: "contact_analysis.html")
- `allow_duplicate_residues` (bool): Allow same residue to appear multiple times (default: False)
- `max_contacts` (int): Maximum contacts to display (default: 20 for unique mode, 25 for duplicate mode)

**Returns:**
- `str`: Path to generated HTML file

#### `HTMLGenerator` Class

**Methods:**

##### `__init__(trajectory_file, topology_file, ligand_file, allow_duplicate_residues, max_contacts)`
Initialize generator with input files.

**Parameters:**
- `trajectory_file` (str): Trajectory file path
- `topology_file` (str): Topology file path
- `ligand_file` (str): Ligand structure file path
- `allow_duplicate_residues` (bool): Allow duplicate residues (default: False)
- `max_contacts` (int): Maximum contacts to display (default: auto-determined)

##### `analyze()`
Run complete contact analysis.

**Returns:**
- `dict`: Analysis results containing:
  - `contact_results`: Raw contact data
  - `ligand_data`: Ligand structure information
  - `contacts`: Processed contact data for visualization
  - `stats`: Summary statistics
  - `traj`: Loaded MDTraj trajectory

##### `generate(output_file="contact_analysis.html")`
Generate HTML visualization file.

**Parameters:**
- `output_file` (str): Output HTML file path

**Returns:**
- `str`: Path to generated HTML file

---

## Dependencies

### Required
- `mdtraj` - Trajectory analysis and manipulation
- `numpy` - Numerical operations and array processing
- `matplotlib` - Static plot generation
- `scikit-learn` - Clustering algorithms

### Optional
- `rdkit` - Enhanced ligand structure handling
- `plotly` - Interactive visualizations

### Included
- Complete JavaScript visualization engine (embedded in HTML output)
- No external JavaScript libraries required
- Works in any modern web browser (Chrome, Firefox, Safari, Edge)

---

## Citation

If you use PRISM analysis module in your research, please cite:

```
[PRISM citation information - to be added]
```

---

## Support

For issues, questions, or feature requests:
- GitHub Issues: [PRISM repository]
- Documentation: [PRISM docs]
- Examples: See `test/analysis/` directory

---

**Documentation Version**: v2.0
**Last Updated**: 2025-01-20
**Maintained by**: PRISM Development Team
