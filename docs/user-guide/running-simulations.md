# Running Simulations

After building your system with PRISM, you can run molecular dynamics simulations using either GROMACS or OpenMM. This guide covers both approaches and best practices.

## Quick Start

### Using GROMACS (Recommended)

```python
import prism

# Load the built system
sim = prism.model("output/GMX_PROLIG_MD")

# Run with GROMACS
results = sim.run(engine="gmx")
```

### Using OpenMM

```python
# Run with OpenMM (requires OpenMM installation)
results = sim.run(engine="openmm", platform="CUDA")
```

### Command-Line Execution

```bash
# Navigate to simulation directory
cd output/GMX_PROLIG_MD

# Run the provided script
bash localrun.sh
```

## Simulation Protocol

PRISM runs a standard 4-stage protocol:

1. **Energy Minimization (EM)** - Remove clashes
2. **NVT Equilibration** - Equilibrate temperature
3. **NPT Equilibration** - Equilibrate pressure
4. **Production MD** - Collect data

## Using the Simulation Module

### Basic Usage

```python
import prism

# Create simulation model from built system
sim = prism.model("output/GMX_PROLIG_MD")

# Check system information
sim.info()

# Run all stages
results = sim.run(engine="gmx")

# Access output files
print(f"Trajectory: {results['prod']['xtc']}")
print(f"Final structure: {results['prod']['gro']}")
```

### Stage Control

Run specific stages only:

```python
# Run only equilibration
results = sim.run(
    engine="gmx",
    stages=["em", "nvt", "npt"]  # Skip production
)

# Run only production (if equilibration done)
results = sim.run(
    engine="gmx",
    stages=["prod"]
)
```

### GPU Acceleration

```python
# GROMACS with GPU
results = sim.run(
    engine="gmx",
    gpu_id=0,  # GPU device ID
    ntomp=10,  # CPU threads
    ntmpi=1    # MPI ranks
)

# OpenMM with CUDA
results = sim.run(
    engine="openmm",
    platform="CUDA",
    device_index=0
)
```

## GROMACS Engine

### Configuration Options

```python
# Full control over GROMACS execution
results = sim.run(
    engine="gmx",
    gpu_id=0,        # GPU to use (-gpu_id)
    ntomp=15,        # OpenMP threads (-ntomp)
    ntmpi=1,         # Thread-MPI ranks (-ntmpi)
    continue_from=None  # Continue from checkpoint
)
```

### Generated Script

PRISM creates a `localrun.sh` script:

```bash
#!/bin/bash

# Energy Minimization
mkdir -p em
gmx grompp -f ../mdps/em.mdp -c solv_ions.gro -r solv_ions.gro \
           -p topol.top -o ./em/em.tpr -maxwarn 10
gmx mdrun -s ./em/em.tpr -deffnm ./em/em -ntmpi 1 -ntomp 10 -gpu_id 0 -v

# NVT Equilibration
mkdir -p nvt
gmx grompp -f ../mdps/nvt.mdp -c ./em/em.gro -r ./em/em.gro \
           -p topol.top -o ./nvt/nvt.tpr -maxwarn 10
gmx mdrun -ntmpi 1 -ntomp 15 -nb gpu -bonded gpu -pme gpu \
          -gpu_id 0 -s ./nvt/nvt.tpr -deffnm ./nvt/nvt -v

# NPT Equilibration
mkdir -p npt
gmx grompp -f ../mdps/npt.mdp -c ./nvt/nvt.gro -r ./nvt/nvt.gro \
           -t ./nvt/nvt.cpt -p topol.top -o ./npt/npt.tpr -maxwarn 10
gmx mdrun -ntmpi 1 -ntomp 15 -nb gpu -bonded gpu -pme gpu \
          -gpu_id 0 -s ./npt/npt.tpr -deffnm ./npt/npt -v

# Production
mkdir -p prod
gmx grompp -f ../mdps/md.mdp -c ./npt/npt.gro -r ./npt/npt.gro \
           -p topol.top -o ./prod/md.tpr -maxwarn 10
gmx mdrun -ntmpi 1 -ntomp 15 -nb gpu -bonded gpu -pme gpu \
          -gpu_id 0 -s ./prod/md.tpr -deffnm ./prod/md -v
```

### Manual GROMACS Commands

For custom control:

```bash
cd output/GMX_PROLIG_MD

# Energy minimization with custom settings
gmx grompp -f ../mdps/em.mdp -c solv_ions.gro -p topol.top -o em.tpr
gmx mdrun -v -deffnm em -nb gpu -gpu_id 0 -ntomp 8

# Production with specific options
gmx mdrun -deffnm prod/md \
  -nb gpu \           # Nonbonded on GPU
  -bonded gpu \       # Bonded on GPU  
  -pme gpu \          # PME on GPU
  -gpu_id 0 \         # GPU device
  -npme 0 \           # PME ranks (0=auto)
  -ntomp 10 \         # OpenMP threads
  -pin on \           # Pin threads
  -pinstride 1        # Pin stride
```

## OpenMM Engine

### Setup and Requirements

```python
# Check OpenMM availability
import prism

sim = prism.model("output/GMX_PROLIG_MD")

# Run with OpenMM
results = sim.run(
    engine="openmm",
    platform="CUDA",     # or "OpenCL", "CPU"
    device_index=0,      # GPU device
    num_threads=10       # For CPU platform
)
```

### Platform Selection

```python
# Auto-select best platform
results = sim.run(engine="openmm", platform="auto")

# Force CPU (for testing)
results = sim.run(engine="openmm", platform="CPU", num_threads=20)

# Multiple GPUs (OpenMM handles differently than GROMACS)
results = sim.run(
    engine="openmm",
    platform="CUDA",
    device_index="0,1"  # Use both GPU 0 and 1
)
```

### OpenMM Advantages

- Better Python integration
- Easier custom forces
- Direct trajectory analysis
- No file I/O overhead

### OpenMM Limitations

- May need manual topology fixes
- Less optimized than GROMACS for some systems
- Different performance characteristics

## Monitoring Simulations

### Check Progress

```python
import os
import time

def monitor_simulation(gmx_dir):
    """Monitor simulation progress"""
    
    log_files = {
        'em': 'em/em.log',
        'nvt': 'nvt/nvt.log',
        'npt': 'npt/npt.log',
        'prod': 'prod/md.log'
    }
    
    for stage, log_file in log_files.items():
        full_path = os.path.join(gmx_dir, log_file)
        if os.path.exists(full_path):
            # Get last line of log
            with open(full_path, 'r') as f:
                lines = f.readlines()
                if lines:
                    last_line = lines[-1]
                    if 'Step' in last_line:
                        print(f"{stage}: {last_line.strip()}")
        else:
            print(f"{stage}: Not started")

# Monitor every 60 seconds
while True:
    monitor_simulation("output/GMX_PROLIG_MD")
    time.sleep(60)
```

### Real-Time Analysis

```python
import mdtraj as md
import matplotlib.pyplot as plt

def plot_rmsd_realtime(traj_file, top_file):
    """Plot RMSD as simulation progresses"""
    
    # Load trajectory
    traj = md.load(traj_file, top=top_file)
    
    # Calculate RMSD
    rmsd = md.rmsd(traj, traj[0])
    
    # Plot
    plt.figure(figsize=(10, 6))
    plt.plot(traj.time, rmsd)
    plt.xlabel('Time (ps)')
    plt.ylabel('RMSD (nm)')
    plt.title('RMSD During Simulation')
    plt.show()

# Check production trajectory
plot_rmsd_realtime(
    "output/GMX_PROLIG_MD/prod/md.xtc",
    "output/GMX_PROLIG_MD/prod/md.tpr"
)
```

## Simulation Parameters

### Understanding MDP Files

PRISM generates optimized MDP files for each stage:

#### Energy Minimization (`em.mdp`)
```
integrator = steep        ; Steepest descent
emtol = 200.0            ; Convergence criterion
nsteps = 10000           ; Maximum steps
```

#### NVT Equilibration (`nvt.mdp`)
```
integrator = md          ; Molecular dynamics
dt = 0.002              ; 2 fs timestep
nsteps = 250000         ; 500 ps total
tcoupl = V-rescale      ; Temperature coupling
ref_t = 310             ; Target temperature (K)
```

#### NPT Equilibration (`npt.mdp`)
```
pcoupl = C-rescale      ; Pressure coupling
ref_p = 1.0            ; Target pressure (bar)
tau_p = 1.0            ; Pressure coupling time
```

#### Production (`md.mdp`)
```
nsteps = 250000000      ; 500 ns
nstxout-compressed = 250000  ; Save every 500 ps
```

### Modifying Parameters

```python
# Modify MDP files before running
def modify_mdp(mdp_file, parameters):
    """Modify MDP parameters"""
    
    with open(mdp_file, 'r') as f:
        lines = f.readlines()
    
    # Update parameters
    new_lines = []
    for line in lines:
        for param, value in parameters.items():
            if line.startswith(param):
                line = f"{param} = {value}\n"
        new_lines.append(line)
    
    with open(mdp_file, 'w') as f:
        f.writelines(new_lines)

# Example: Longer production run
modify_mdp(
    "output/mdps/md.mdp",
    {"nsteps": 500000000}  # 1 microsecond
)
```

## Restarting Simulations

### From Checkpoint

```bash
# GROMACS restart from checkpoint
cd output/GMX_PROLIG_MD/prod
gmx mdrun -s md.tpr -cpi md.cpt -deffnm md -append
```

```python
# Python restart
sim = prism.model("output/GMX_PROLIG_MD")
results = sim.run(
    engine="gmx",
    continue_from="npt"  # Continue from NPT checkpoint
)
```

### Extending Simulations

```bash
# Extend by 100 ns
gmx convert-tpr -s prod/md.tpr -extend 100000 -o prod/md_extended.tpr
gmx mdrun -s prod/md_extended.tpr -cpi prod/md.cpt -deffnm prod/md -append
```

## Multi-Simulation Workflows

### Batch Processing

```python
import prism
from pathlib import Path

# Run multiple systems
systems = Path("systems").glob("*/GMX_PROLIG_MD")

for system_dir in systems:
    try:
        sim = prism.model(str(system_dir))
        results = sim.run(engine="gmx")
        print(f"✓ Completed: {system_dir.parent.name}")
    except Exception as e:
        print(f"✗ Failed: {system_dir.parent.name}: {e}")
```

### Parallel Simulations

```python
from concurrent.futures import ProcessPoolExecutor
import prism

def run_simulation(gmx_dir):
    """Run single simulation"""
    sim = prism.model(gmx_dir)
    return sim.run(engine="gmx")

# Run 4 simulations in parallel
gmx_dirs = [
    "system1/GMX_PROLIG_MD",
    "system2/GMX_PROLIG_MD",
    "system3/GMX_PROLIG_MD",
    "system4/GMX_PROLIG_MD"
]

with ProcessPoolExecutor(max_workers=4) as executor:
    results = list(executor.map(run_simulation, gmx_dirs))
```

### Replica Exchange

```python
# Setup for REMD (manual process)
temperatures = [300, 310, 320, 330, 340, 350]

for i, temp in enumerate(temperatures):
    # Modify MDP for each temperature
    modify_mdp(
        f"replica_{i}/mdps/md.mdp",
        {"ref_t": temp}
    )
    
    # Build TPR files
    os.system(f"""
        cd replica_{i}/GMX_PROLIG_MD
        gmx grompp -f ../mdps/md.mdp -c npt.gro \
                   -p topol.top -o md_{i}.tpr
    """)

# Run REMD
os.system(f"""
    mpirun -np {len(temperatures)} gmx_mpi mdrun -multidir replica_* \
           -deffnm md -replex 1000
""")
```

## Performance Optimization

### Hardware Detection

```python
import subprocess
import torch  # If available

def detect_hardware():
    """Detect available hardware"""
    
    info = {}
    
    # CPU cores
    import multiprocessing
    info['cpu_cores'] = multiprocessing.cpu_count()
    
    # GPU detection
    try:
        result = subprocess.run(
            ['nvidia-smi', '--query-gpu=name,memory.total', 
             '--format=csv,noheader'],
            capture_output=True, text=True
        )
        if result.returncode == 0:
            gpus = result.stdout.strip().split('\n')
            info['gpus'] = gpus
    except:
        info['gpus'] = []
    
    # CUDA availability (if torch installed)
    try:
        info['cuda_available'] = torch.cuda.is_available()
        info['cuda_devices'] = torch.cuda.device_count()
    except:
        pass
    
    return info

hardware = detect_hardware()
print(f"Hardware: {hardware}")

# Optimize based on hardware
if hardware.get('gpus'):
    # Use GPU
    sim.run(engine="gmx", gpu_id=0, ntomp=hardware['cpu_cores']//2)
else:
    # CPU only
    sim.run(engine="gmx", ntomp=hardware['cpu_cores'])
```

### Performance Tuning

```python
# Benchmark different settings
import time

settings = [
    {"ntomp": 8, "gpu_id": 0},
    {"ntomp": 10, "gpu_id": 0},
    {"ntomp": 12, "gpu_id": 0},
    {"ntomp": 15, "gpu_id": 0},
]

for setting in settings:
    start = time.time()
    
    # Run short test
    sim.run(
        engine="gmx",
        stages=["em"],
        **setting
    )
    
    elapsed = time.time() - start
    print(f"Settings {setting}: {elapsed:.2f} seconds")
```

### Performance Metrics

```python
def parse_performance(log_file):
    """Extract performance metrics from GROMACS log"""
    
    with open(log_file) as f:
        lines = f.readlines()
    
    for line in lines:
        if 'Performance:' in line:
            # Extract ns/day
            parts = line.split()
            ns_per_day = float(parts[1])
            hours_per_ns = float(parts[3])
            return {
                'ns_per_day': ns_per_day,
                'hours_per_ns': hours_per_ns
            }
    
    return None

# Check performance
perf = parse_performance("output/GMX_PROLIG_MD/prod/md.log")
if perf:
    print(f"Performance: {perf['ns_per_day']:.2f} ns/day")
    
    # Estimate completion time
    total_ns = 500  # Target simulation length
    days_needed = total_ns / perf['ns_per_day']
    print(f"Estimated time for 500 ns: {days_needed:.1f} days")
```

## Troubleshooting

### Common Issues

#### "CUDA Error"

```python
# Fallback to CPU if GPU fails
try:
    results = sim.run(engine="openmm", platform="CUDA")
except:
    print("GPU failed, using CPU")
    results = sim.run(engine="openmm", platform="CPU")
```

#### "Atoms moving too fast"

```bash
# Reduce time step
modify_mdp("mdps/md.mdp", {"dt": 0.001})  # 1 fs instead of 2 fs

# Or increase equilibration
modify_mdp("mdps/nvt.mdp", {"nsteps": 500000})  # Longer NVT
```

#### "Checkpoint file corrupted"

```bash
# Start fresh from last good structure
gmx grompp -f ../mdps/md.mdp -c npt.gro -p topol.top -o md_new.tpr
gmx mdrun -s md_new.tpr -deffnm md_new
```

#### "Out of memory"

```python
# Reduce trajectory output frequency
modify_mdp(
    "mdps/md.mdp",
    {"nstxout-compressed": 500000}  # Save less frequently
)

# Or compress existing trajectory
os.system("gmx trjconv -f md.xtc -o md_compressed.xtc -dt 10")
```

## Analysis During Simulation

### On-the-fly Analysis

```python
import mdtraj as md
import numpy as np

def analyze_while_running(gmx_dir):
    """Analyze trajectory while simulation runs"""
    
    traj_file = f"{gmx_dir}/prod/md.xtc"
    tpr_file = f"{gmx_dir}/prod/md.tpr"
    
    if not os.path.exists(traj_file):
        print("No trajectory yet")
        return
    
    # Load current trajectory
    traj = md.load(traj_file, top=tpr_file)
    
    # Quick analyses
    print(f"Frames: {traj.n_frames}")
    print(f"Time: {traj.time[-1]} ps")
    
    # RMSD
    rmsd = md.rmsd(traj, traj[0])
    print(f"Current RMSD: {rmsd[-1]:.3f} nm")
    
    # Radius of gyration
    rg = md.compute_rg(traj)
    print(f"Current Rg: {rg[-1]:.3f} nm")
    
    return {
        'frames': traj.n_frames,
        'time': traj.time[-1],
        'rmsd': rmsd[-1],
        'rg': rg[-1]
    }
```

### Energy Monitoring

```bash
# Extract energies during run
gmx energy -f prod/md.edr -o energy.xvg << EOF
Potential
Kinetic-En.
Total-Energy
Temperature
Pressure
EOF

# Plot with xmgrace or Python
```

## Best Practices

### 1. Always Equilibrate
Never skip equilibration stages:
- EM removes clashes
- NVT stabilizes temperature
- NPT stabilizes density

### 2. Check Convergence
Monitor key properties:
```python
# Check if equilibrated
def check_equilibration(edr_file):
    """Check if system is equilibrated"""
    
    os.system(f"""
        echo "Temperature Pressure Density" | 
        gmx energy -f {edr_file} -o check.xvg
    """)
    
    # Load and analyze
    data = np.loadtxt("check.xvg", comments=['#', '@'])
    
    # Check last 20% is stable
    n = len(data)
    last_portion = data[int(0.8*n):]
    
    temp_std = np.std(last_portion[:, 1])
    press_std = np.std(last_portion[:, 2])
    
    print(f"Temperature StdDev: {temp_std:.2f} K")
    print(f"Pressure StdDev: {press_std:.2f} bar")
    
    if temp_std < 5 and press_std < 100:
        print("System appears equilibrated")
        return True
    else:
        print("System may need more equilibration")
        return False
```

### 3. Save Checkpoints
```python
# Ensure checkpoint saving
modify_mdp("mdps/md.mdp", {
    "nstcheckpoint": 1000000  # Every 2 ns
})
```

### 4. Use Appropriate Resources
- GPU for nonbonded calculations
- Multiple CPUs for PME
- Balance GPU and CPU load

### 5. Document Everything
```python
# Save simulation metadata
import json
from datetime import datetime

metadata = {
    "date": datetime.now().isoformat(),
    "system": "protein_ligand_complex",
    "engine": "gromacs",
    "gpu": "RTX_3080",
    "performance": "45 ns/day",
    "parameters": {
        "temperature": 310,
        "pressure": 1.0,
        "time": "500 ns"
    }
}

with open("simulation_metadata.json", "w") as f:
    json.dump(metadata, f, indent=2)
```

## Next Steps

- [Analyze trajectories](analysis-tools.md)
- [Understand outputs](output-files.md)
- [Advanced workflows](advanced-usage.md)