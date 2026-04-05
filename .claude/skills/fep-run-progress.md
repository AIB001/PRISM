# FEP Run Progress Check

Check FEP run progress and report complete status.

## Check Steps

### 1. Check Current mdrun Process Count
```bash
ps aux | grep "gmx mdrun" | grep -v grep | wc -l
```

### 2. Enter FEP Directory
```bash
cd /data2/gxf1212/work/PRISM/tests/gxf/FEP/unit_test/42-38/CHARMM-GUI/GMX_PROLIG_FEP
```

### 3. Report Bound Repeat Progress
Check repeat1, repeat2, repeat3 for:
- NPT completion count
- Production windows completion count (11 windows total)

### 4. Report Unbound Progress
If started, check each repeat's progress

### 5. Report GPU Usage
```bash
nvidia-smi --query-gpu=index,name,utilization.gpu,memory.used,memory.total --format=csv,noheader,nounits
```

### 6. Check GPU Process Allocation
For each running mdrun process, check its GPU allocation:
```bash
cat /proc/<PID>/environ | tr '\0' '\n' | grep CUDA
```

### 7. GPU Conflict Resolution
**Critical Rule**: If multiple GMX processes are running on the same GPU, you must terminate the extra processes!

**Check Method**:
- Check each process's CUDA_VISIBLE_DEVICES environment variable
- If multiple processes are on the same GPU, keep only the one with the longest runtime
- Use `kill -9 <PID>` to terminate conflicting processes

**Verify**: Ensure only 1 mdrun process per active GPU

### 8. Auto-start Condition
**If no mdrun processes are currently running**, automatically start the full FEP workflow:
```bash
cd /data2/gxf1212/work/PRISM/tests/gxf/FEP/unit_test/42-38/CHARMM-GUI/GMX_PROLIG_FEP && bash run_fep.sh
```

## Report Format

### Progress Report Includes:
- Current mdrun process count
- NPT and Production completion progress for each repeat
- GPU usage (utilization, memory usage)
- GPU process allocation details
- GPU conflict status and resolution results
- Overall progress percentage

### Sample Report:
```
## 📊 FEP Run Progress Report

### Current Run Status
- **mdrun processes**: X processes running
- **Run location**: Bound/Unbound repeatX
- **GPU conflicts**: ✅ Resolved / ⚠️ Conflicts exist

### GPU Usage
GPU 0: X% utilization, XXXMB/32768MB memory
GPU 1: X% utilization, XXXMB/32768MB memory
...

### Currently Running Windows
Window XX: 🔄 Running (X minutes)
...

### Bound Repeat Progress
**Repeat X**:
- ✅ NPT: Complete (1/1)
- ✅ Production: X/11 windows complete
...
```

## Important Reminders

1. **GPU conflicts must be resolved immediately**: Do not allow multiple processes to share the same GPU
2. **Check GPU allocation for all processes**: Use /proc/<PID>/environ to check CUDA_VISIBLE_DEVICES
3. **Prioritize long-running processes**: Terminate new processes, keep old processes
4. **Verify cleanup results**: After handling, verify remaining process count and GPU allocation
5. **Report in detail**: Include process ID, runtime, GPU allocation, and other detailed information

## FEP System Locations

Main test systems:
- 42-38 system: `/data2/gxf1212/work/PRISM/tests/gxf/FEP/unit_test/42-38/CHARMM-GUI/GMX_PROLIG_FEP`
- p38-19-24 system: `/data2/gxf1212/work/PRISM/tests/gxf/FEP/unit_test/p38-19-24/rtf_fep_output/GMX_PROLIG_FEP`
