# Analysis Tools

PRISM provides comprehensive analysis and visualization tools for molecular dynamics trajectories, with special focus on protein-ligand interactions.

## Quick Start

### Basic Analysis

```python
from prism.analysis.visualization import generate_html

# Generate interactive visualization
generate_html(
    trajectory="output/GMX_PROLIG_MD/prod/md.xtc",
    topology="output/GMX_PROLIG_MD/prod/md.tpr", 
    ligand="ligand.sdf",
    output="contact_analysis.html"
)

# Open contact_analysis.html in browser for interactive visualization
```

### Using the Analysis Module

```python
from prism.analysis.visualization import HTMLGenerator

# Initialize analyzer
analyzer = HTMLGenerator(
    trajectory_file="trajectory.xtc",
    topology_file="topology.gro",
    ligand_file="ligand.sdf"
)

# Run analysis
results = analyzer.analyze()

# Generate visualization
analyzer.generate("my_analysis.html")
```

## Contact Analysis

### Protein-Ligand Contacts

PRISM automatically identifies and analyzes protein-ligand contacts:

```python
from prism.analysis.visualization import FastContactAnalyzer
import mdtraj as md

# Load trajectory
traj = md.load("trajectory.xtc", top="topology.gro")

# Analyze contacts
analyzer = FastContactAnalyzer(traj)
contact_results = analyzer.calculate_contact_proportions()

# Results include:
print(f"Contact frequencies: {contact_results['contact_frequencies']}")
print(f"Residue proportions: {contact_results['residue_proportions']}")
print(f"Average distances: {contact_results['residue_avg_distances']}")
```

### Understanding Contact Metrics

**Contact Frequency**: Fraction of frames where contact exists
- > 0.8: Very strong interaction
- 0.5-0.8: Strong interaction  
- 0.2-0.5: Moderate interaction
- < 0.2: Weak/transient interaction

**Distance Cutoffs**:
- Enter: 3.5 Å (0.35 nm)
- Exit: 4.0 Å (0.40 nm)
- Maximum: 5.0 Å (0.50 nm)

## Interactive HTML Visualization

### Features

The HTML visualization provides:
- **2D/3D molecular view** with real-time rotation
- **Contact frequency heatmap** with color coding
- **Interactive residue network** 
- **TOP 3 contacts** highlighted
- **Export capabilities** for publication

### Customization

```python
from prism.analysis.visualization import HTMLGenerator

# Custom configuration
generator = HTMLGenerator(
    trajectory_file="traj.xtc",
    topology_file="topol.gro",
    ligand_file="ligand.sdf"
)

# Modify configuration
generator.config.contact_enter_threshold_nm = 0.30  # Stricter cutoff
generator.config.distance_cutoff_nm = 0.45  # Shorter max distance

# Generate with custom settings
generator.generate("custom_analysis.html")
```

### Visualization Controls

In the generated HTML:
- **Drag residues**: Rearrange in 2D mode
- **Toggle 3D**: Switch between 2D/3D views
- **Show/Hide**: Connections, hydrogen atoms
- **Export**: High-resolution PNG (up to 8K)
- **Zoom**: Mouse wheel or buttons
- **Pan**: Middle-click or Shift+drag

## Trajectory Analysis with MDTraj

### Basic Trajectory Analysis

```python
import mdtraj as md
import numpy as np
import matplotlib.pyplot as plt

# Load trajectory
traj = md.load("output/GMX_PROLIG_MD/prod/md.xtc", 
               top="output/GMX_PROLIG_MD/solv_ions.gro")

print(f"Trajectory: {traj.n_frames} frames, {traj.n_atoms} atoms")
print(f"Time: {traj.time[0]} to {traj.time[-1]} ps")
```

### RMSD Analysis

```python
# Protein RMSD
protein_atoms = traj.topology.select("protein and name CA")
protein_traj = traj.atom_slice(protein_atoms)

# Align and calculate RMSD
protein_traj.superpose(protein_traj[0])
rmsd = md.rmsd(protein_traj, protein_traj[0])

# Plot
plt.figure(figsize=(10, 6))
plt.plot(traj.time, rmsd * 10)  # Convert to Angstroms
plt.xlabel("Time (ps)")
plt.ylabel("RMSD (Å)")
plt.title("Protein Cα RMSD")
plt.show()

# Ligand RMSD
ligand_atoms = traj.topology.select("resname LIG")
ligand_rmsd = md.rmsd(traj, traj[0], atom_indices=ligand_atoms)
```

### RMSF Analysis

```python
# Calculate RMSF per residue
rmsf = md.rmsf(protein_traj, protein_traj[0])

# Get residue numbers
residues = [protein_traj.topology.atom(i).residue.resSeq 
            for i in range(protein_traj.n_atoms)]

# Plot
plt.figure(figsize=(12, 6))
plt.plot(residues, rmsf * 10)
plt.xlabel("Residue Number")
plt.ylabel("RMSF (Å)")
plt.title("Residue Flexibility")
plt.show()
```

### Radius of Gyration

```python
# Calculate Rg over time
rg = md.compute_rg(traj.atom_slice(traj.topology.select("protein")))

plt.figure(figsize=(10, 6))
plt.plot(traj.time, rg)
plt.xlabel("Time (ps)")
plt.ylabel("Radius of Gyration (nm)")
plt.title("Protein Compactness")
plt.show()

# Check for unfolding
if rg[-1] > rg[0] * 1.2:
    print("Warning: Protein may be unfolding")
```

## Hydrogen Bond Analysis

```python
# Find hydrogen bonds
hbonds = md.baker_hubbard(traj, freq=0.5)  # Present in >50% of frames

print(f"Found {len(hbonds)} hydrogen bonds")

# Analyze protein-ligand H-bonds
ligand_atoms = set(traj.topology.select("resname LIG"))
protein_ligand_hbonds = []

for donor, hydrogen, acceptor in hbonds:
    if (donor in ligand_atoms) != (acceptor in ligand_atoms):
        # One is ligand, other is protein
        protein_ligand_hbonds.append((donor, hydrogen, acceptor))
        
print(f"Protein-ligand H-bonds: {len(protein_ligand_hbonds)}")

# Get details
for donor, hydrogen, acceptor in protein_ligand_hbonds[:5]:
    d_atom = traj.topology.atom(donor)
    a_atom = traj.topology.atom(acceptor)
    print(f"{d_atom.residue}-{d_atom.name} -> {a_atom.residue}-{a_atom.name}")
```

## Binding Site Analysis

### Contact Residues

```python
# Identify binding site residues
cutoff = 0.5  # nm
ligand_atoms = traj.topology.select("resname LIG")

# Find contacts
contacts = md.compute_neighbors(traj, cutoff, ligand_atoms)

# Get unique residues
binding_residues = set()
for frame_contacts in contacts:
    for atom_idx in frame_contacts:
        residue = traj.topology.atom(atom_idx).residue
        if residue.name not in ['LIG', 'HOH', 'WAT']:
            binding_residues.add(residue)

print(f"Binding site residues: {len(binding_residues)}")
for res in sorted(binding_residues, key=lambda x: x.resSeq):
    print(f"  {res.name}{res.resSeq}")
```

### Pocket Volume

```python
from scipy.spatial import ConvexHull

def calculate_pocket_volume(traj, pocket_residues):
    """Calculate binding pocket volume over time"""
    
    volumes = []
    for frame in traj:
        # Get pocket atom coordinates
        pocket_atoms = []
        for res in pocket_residues:
            for atom in res.atoms:
                pocket_atoms.append(atom.index)
        
        coords = frame.xyz[0, pocket_atoms]
        
        # Calculate convex hull volume
        try:
            hull = ConvexHull(coords)
            volumes.append(hull.volume)
        except:
            volumes.append(np.nan)
    
    return np.array(volumes)

# Calculate and plot
pocket_volume = calculate_pocket_volume(traj, binding_residues)

plt.figure(figsize=(10, 6))
plt.plot(traj.time, pocket_volume)
plt.xlabel("Time (ps)")
plt.ylabel("Pocket Volume (nm³)")
plt.title("Binding Pocket Volume")
plt.show()
```

## Ligand Dynamics

### Ligand Movement

```python
# Track ligand center of mass
ligand_atoms = traj.topology.select("resname LIG")
ligand_com = md.compute_center_of_mass(
    traj.atom_slice(ligand_atoms)
)

# Calculate displacement
displacement = np.linalg.norm(ligand_com - ligand_com[0], axis=1)

plt.figure(figsize=(10, 6))
plt.plot(traj.time, displacement * 10)
plt.xlabel("Time (ps)")
plt.ylabel("Displacement (Å)")
plt.title("Ligand Movement from Initial Position")
plt.show()
```

### Ligand Orientation

```python
# Define ligand vectors (customize for your ligand)
def get_ligand_orientation(traj):
    """Calculate ligand orientation angles"""
    
    # Example: vector between two specific atoms
    atom1 = traj.topology.select("resname LIG and name C1")[0]
    atom2 = traj.topology.select("resname LIG and name C10")[0]
    
    vectors = traj.xyz[:, atom2] - traj.xyz[:, atom1]
    
    # Calculate angles relative to initial
    initial = vectors[0]
    angles = []
    
    for vec in vectors:
        cos_angle = np.dot(vec, initial) / (
            np.linalg.norm(vec) * np.linalg.norm(initial)
        )
        angle = np.arccos(np.clip(cos_angle, -1, 1))
        angles.append(np.degrees(angle))
    
    return np.array(angles)

angles = get_ligand_orientation(traj)
plt.plot(traj.time, angles)
plt.xlabel("Time (ps)")
plt.ylabel("Rotation Angle (degrees)")
plt.show()
```

## Free Energy Calculations

### MM-PBSA (Simplified)

```python
# Basic MM-PBSA calculation (requires additional tools)
def calculate_binding_energy(traj, topology):
    """Simplified binding energy calculation"""
    
    # This is a simplified example
    # Real MM-PBSA requires specialized software
    
    # Select components
    protein = traj.topology.select("protein")
    ligand = traj.topology.select("resname LIG")
    complex_atoms = traj.topology.select("protein or resname LIG")
    
    energies = []
    
    for frame in traj[::10]:  # Every 10th frame
        # Would calculate:
        # E_binding = E_complex - E_protein - E_ligand
        # Including solvation effects
        pass
    
    return energies

# Note: Use gmx_MMPBSA or similar tools for accurate calculations
```

## Water Analysis

### Water Density

```python
# Analyze water distribution around ligand
def water_density_around_ligand(traj, cutoff=0.5):
    """Calculate water density around ligand"""
    
    ligand = traj.topology.select("resname LIG")
    water_oxygens = traj.topology.select("water and name O")
    
    water_counts = []
    
    for frame in traj:
        # Count waters within cutoff
        distances = md.compute_distances(
            frame,
            [[lig, wat] for lig in ligand for wat in water_oxygens]
        )
        
        within_cutoff = np.sum(distances < cutoff)
        water_counts.append(within_cutoff)
    
    return np.array(water_counts)

water_count = water_density_around_ligand(traj)

plt.figure(figsize=(10, 6))
plt.plot(traj.time, water_count)
plt.xlabel("Time (ps)")
plt.ylabel("Number of Water Molecules")
plt.title(f"Waters within 5 Å of Ligand")
plt.show()
```

## Clustering Analysis

### Conformational Clustering

```python
from sklearn.cluster import KMeans

# Cluster ligand conformations
ligand_atoms = traj.topology.select("resname LIG")
ligand_traj = traj.atom_slice(ligand_atoms)

# Align ligand trajectory
ligand_traj.superpose(ligand_traj[0])

# Reshape for clustering
X = ligand_traj.xyz.reshape((ligand_traj.n_frames, -1))

# Perform clustering
n_clusters = 5
kmeans = KMeans(n_clusters=n_clusters, random_state=42)
clusters = kmeans.fit_predict(X)

# Analyze clusters
for i in range(n_clusters):
    cluster_frames = np.where(clusters == i)[0]
    print(f"Cluster {i}: {len(cluster_frames)} frames "
          f"({100*len(cluster_frames)/len(clusters):.1f}%)")
    
    # Save representative structure
    representative = cluster_frames[len(cluster_frames)//2]
    ligand_traj[representative].save(f"cluster_{i}.pdb")
```

## Custom Analysis Functions

### Create Analysis Pipeline

```python
class TrajectoryAnalyzer:
    """Custom analysis pipeline"""
    
    def __init__(self, trajectory, topology):
        self.traj = md.load(trajectory, top=topology)
        self.results = {}
        
    def analyze_all(self):
        """Run all analyses"""
        self.calculate_rmsd()
        self.calculate_rmsf()
        self.analyze_contacts()
        self.analyze_hbonds()
        return self.results
    
    def calculate_rmsd(self):
        """Calculate RMSD"""
        protein = self.traj.topology.select("protein and name CA")
        rmsd = md.rmsd(self.traj, self.traj[0], atom_indices=protein)
        self.results['rmsd'] = {
            'values': rmsd,
            'mean': np.mean(rmsd),
            'std': np.std(rmsd)
        }
    
    def calculate_rmsf(self):
        """Calculate RMSF"""
        protein = self.traj.topology.select("protein and name CA")
        traj_ca = self.traj.atom_slice(protein)
        rmsf = md.rmsf(traj_ca, traj_ca[0])
        self.results['rmsf'] = rmsf
    
    def analyze_contacts(self):
        """Analyze contacts"""
        # Your contact analysis
        pass
    
    def analyze_hbonds(self):
        """Analyze hydrogen bonds"""
        hbonds = md.baker_hubbard(self.traj, freq=0.5)
        self.results['hbonds'] = hbonds
    
    def save_report(self, filename="analysis_report.txt"):
        """Save analysis report"""
        with open(filename, 'w') as f:
            f.write("=== Trajectory Analysis Report ===\n\n")
            
            if 'rmsd' in self.results:
                f.write(f"RMSD: {self.results['rmsd']['mean']:.3f} ± "
                       f"{self.results['rmsd']['std']:.3f} nm\n")
            
            if 'hbonds' in self.results:
                f.write(f"Hydrogen bonds: {len(self.results['hbonds'])}\n")

# Use the analyzer
analyzer = TrajectoryAnalyzer(
    "trajectory.xtc",
    "topology.gro"
)
results = analyzer.analyze_all()
analyzer.save_report()
```

## Visualization with NGLView

### Interactive 3D Visualization

```python
import nglview as nv
import mdtraj as md

# Load trajectory
traj = md.load("trajectory.xtc", top="topology.gro")

# Create viewer
view = nv.show_mdtraj(traj)

# Customize representation
view.clear_representations()
view.add_cartoon(selection="protein", color="secondary structure")
view.add_licorice(selection="LIG", color="element")
view.add_ball_and_stick(selection="protein and (sidechain and not hydrogen)")

# Add surfaces
view.add_surface(selection="protein", opacity=0.3)

# Display
view
```

## Performance Tips

### Efficient Loading

```python
# Load only specific atoms
protein_ligand = traj.topology.select("protein or resname LIG")
traj_subset = md.load("trajectory.xtc", 
                      top="topology.gro",
                      atom_indices=protein_ligand)

# Load specific time range
traj_partial = md.load("trajectory.xtc",
                      top="topology.gro", 
                      stride=10,  # Every 10th frame
                      frame=1000)  # Start from frame 1000
```

### Parallel Analysis

```python
from multiprocessing import Pool
import functools

def analyze_chunk(chunk_idx, traj_file, topology):
    """Analyze trajectory chunk"""
    # Load chunk
    traj = md.load(traj_file, top=topology,
                  frame=chunk_idx*1000, 
                  n_frames=1000)
    
    # Analyze
    rmsd = md.rmsd(traj, traj[0])
    
    return {'chunk': chunk_idx, 'rmsd': rmsd}

# Parallel analysis
with Pool(4) as pool:
    analyze_func = functools.partial(
        analyze_chunk,
        traj_file="trajectory.xtc",
        topology="topology.gro"
    )
    results = pool.map(analyze_func, range(4))
```

## Export Results

### Publication Figures

```python
import matplotlib.pyplot as plt
import seaborn as sns

# Set publication style
plt.style.use('seaborn-v0_8-paper')
sns.set_palette("husl")

fig, axes = plt.subplots(2, 2, figsize=(12, 10))

# RMSD
axes[0,0].plot(traj.time/1000, rmsd*10, linewidth=2)
axes[0,0].set_xlabel("Time (ns)")
axes[0,0].set_ylabel("RMSD (Å)")
axes[0,0].set_title("A. Protein Stability")

# RMSF
axes[0,1].plot(residues, rmsf*10, linewidth=2)
axes[0,1].set_xlabel("Residue")
axes[0,1].set_ylabel("RMSF (Å)")
axes[0,1].set_title("B. Residue Flexibility")

# Contacts
axes[1,0].hist(contact_frequencies, bins=20, edgecolor='black')
axes[1,0].set_xlabel("Contact Frequency")
axes[1,0].set_ylabel("Count")
axes[1,0].set_title("C. Contact Distribution")

# H-bonds
axes[1,1].plot(traj.time/1000, hbond_count, linewidth=2)
axes[1,1].set_xlabel("Time (ns)")
axes[1,1].set_ylabel("H-bonds")
axes[1,1].set_title("D. Hydrogen Bonds")

plt.tight_layout()
plt.savefig("analysis_figure.pdf", dpi=300, bbox_inches='tight')
plt.savefig("analysis_figure.png", dpi=300, bbox_inches='tight')
```

### Data Export

```python
import pandas as pd

# Create results dataframe
results_df = pd.DataFrame({
    'Time_ns': traj.time / 1000,
    'RMSD_A': rmsd * 10,
    'Rg_nm': rg,
    'Hbonds': hbond_count,
    'Ligand_displacement_A': displacement * 10
})

# Save to CSV
results_df.to_csv("analysis_results.csv", index=False)

# Save to Excel with multiple sheets
with pd.ExcelWriter("analysis_results.xlsx") as writer:
    results_df.to_excel(writer, sheet_name="Time Series", index=False)
    
    # Add contact data
    contacts_df = pd.DataFrame(contact_results['residue_proportions'].items(),
                              columns=['Residue', 'Contact_Frequency'])
    contacts_df.to_excel(writer, sheet_name="Contacts", index=False)
```

## Next Steps

- Review [Output Files](output-files.md) structure
- Learn [Advanced Usage](advanced-usage.md) techniques
- Troubleshoot issues with [Troubleshooting Guide](troubleshooting.md)