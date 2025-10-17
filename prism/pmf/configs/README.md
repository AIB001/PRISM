# PMF Configuration Files

This directory contains configuration files for PRISM PMF workflows.

## Configuration Files

### `equilibration_config.yaml`
Complete configuration template with all available options for PMF equilibration.

**Features:**
- Detailed Slurm job parameters
- Local execution settings
- Multiple predefined profiles (quick, standard, long, cpu_only)
- Complete simulation parameters
- Advanced options for optimization

### `equilibration_simple.yaml`
Simplified configuration for basic usage.

**Features:**
- Minimal required settings
- Easy-to-understand structure
- Basic profiles (quick, long)
- Suitable for most users

## Usage

### Using Default Configuration
```bash
# Uses equilibration_simple.yaml with default profile
python pmf_remodel_only.py --input ./gromacssim --output ./pmf_test
```

### Using Specific Profile
```bash
# Uses quick profile (6 hours, 8 CPUs)
python pmf_remodel_only.py --input ./gromacssim --output ./pmf_test \
    --equilibration-profile quick

# Uses long profile (48 hours, 16 CPUs)
python pmf_remodel_only.py --input ./gromacssim --output ./pmf_test \
    --equilibration-profile long
```

### Using Custom Configuration
```bash
# Uses your custom config file
python pmf_remodel_only.py --input ./gromacssim --output ./pmf_test \
    --equilibration-config ./my_custom_config.yaml
```

## Creating Custom Configuration

1. Copy `equilibration_simple.yaml` to your working directory
2. Modify the settings according to your cluster/system
3. Use the custom config with `--equilibration-config` option

### Important Settings to Customize

#### For Slurm Clusters
```yaml
slurm:
  partition: "your_gpu_partition"    # Your cluster's GPU partition name
  account: "your_account"            # If required by your cluster
  gres: "gpu:tesla:1"               # Adjust GPU resource format
  conda_env: "your_env_name"        # Your conda environment name
```

#### For Different Hardware
```yaml
local:
  cpus: 20           # Number of CPU cores
  gpu_id: 0          # GPU device ID (0, 1, 2, etc.)

slurm:
  cpus_per_task: 20  # CPUs per Slurm task
  memory: "64G"      # Memory allocation
```

## Configuration File Search Order

The system searches for configuration files in this order:

1. File specified with `--equilibration-config`
2. `./equilibration_config.yaml` (current directory)
3. `~/.prism/equilibration_config.yaml` (user home)
4. `prism/pmf/configs/equilibration_simple.yaml` (package default)

## Profiles

Profiles allow quick selection of predefined configurations:

- **default**: Uses the base configuration
- **quick**: Fast equilibration (6h, reduced resources)
- **standard**: Standard production equilibration (24h)
- **long**: Thorough equilibration (48h, more resources)
- **cpu_only**: For clusters without GPU access

## Example Workflow

1. Test with quick profile first:
   ```bash
   python pmf_remodel_only.py -i ./input -o ./test --equilibration-profile quick
   ```

2. If successful, run production with standard profile:
   ```bash
   python pmf_remodel_only.py -i ./input -o ./production --equilibration-profile standard
   ```

3. For critical calculations, use long profile:
   ```bash
   python pmf_remodel_only.py -i ./input -o ./final --equilibration-profile long
   ```