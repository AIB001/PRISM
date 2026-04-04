# FEP运行进度检查

检查FEP运行进度并报告完整状态。

## 检查步骤

### 1. 查看当前mdrun进程数
```bash
ps aux | grep "gmx mdrun" | grep -v grep | wc -l
```

### 2. 进入FEP目录
```bash
cd /data2/gxf1212/work/PRISM/tests/gxf/FEP/unit_test/42-38/CHARMM-GUI/GMX_PROLIG_FEP
```

### 3. 报告Bound各repeat进度
检查repeat1、repeat2、repeat3的：
- NPT完成数
- Production windows完成数（共11个windows）

### 4. 报告Unbound进度
如果已开始，检查各repeat的进度

### 5. 报告GPU使用情况
```bash
nvidia-smi --query-gpu=index,name,utilization.gpu,memory.used,memory.total --format=csv,noheader,nounits
```

### 6. 检查GPU进程分配
对于每个运行的mdrun进程，检查其GPU分配：
```bash
cat /proc/<PID>/environ | tr '\0' '\n' | grep CUDA
```

### 7. GPU冲突处理
**关键规则**：如果有多个GMX进程运行在同一个GPU上，必须中断多余的进程！

**检查方法**：
- 查看每个进程的CUDA_VISIBLE_DEVICES环境变量
- 如果同一个GPU上有多个进程，只保留运行时间最长的
- 使用 `kill -9 <PID>` 终止冲突的进程

**验证**：确保每个活跃GPU上只有1个mdrun进程

### 8. 自动启动条件
**如果当前没有mdrun进程在运行**，自动启动完整FEP流程：
```bash
cd /data2/gxf1212/work/PRISM/tests/gxf/FEP/unit_test/42-38/CHARMM-GUI/GMX_PROLIG_FEP && bash run_fep.sh
```

## 报告格式

### 进度报告包含：
- 当前mdrun进程数
- 各repeat的NPT和Production完成进度
- GPU使用情况（利用率、显存占用）
- GPU进程分配详情
- GPU冲突状态和处理结果
- 总体进度百分比

### 示例报告：
```
## 📊 FEP运行进度报告

### 当前运行状态
- **mdrun进程数**: X个进程正在运行
- **运行位置**: Bound/Unbound repeatX
- **GPU冲突**: ✅ 已解决 / ⚠️ 存在冲突

### GPU使用情况
GPU 0: X% 利用率, XXXMB/32768MB 内存
GPU 1: X% 利用率, XXXMB/32768MB 内存
...

### 当前运行的Windows
Window XX: 🔄 运行中 (X分钟)
...

### Bound各Repeat进度
**Repeat X**:
- ✅ NPT: 完成 (1/1)
- ✅ Production: X/11 windows 完成
...
```

## 重要提醒

1. **GPU冲突必须立即处理**：不要让多个进程共享同一个GPU
2. **检查所有进程的GPU分配**：使用/proc/<PID>/environ查看CUDA_VISIBLE_DEVICES
3. **优先保留运行时间长的进程**：中断新进程，保留老进程
4. **验证清理结果**：处理后验证剩余进程数和GPU分配
5. **报告要详细**：包含进程ID、运行时间、GPU分配等详细信息

## FEP系统位置

主要测试系统：
- 42-38系统: `/data2/gxf1212/work/PRISM/tests/gxf/FEP/unit_test/42-38/CHARMM-GUI/GMX_PROLIG_FEP`
- p38-19-24系统: `/data2/gxf1212/work/PRISM/tests/gxf/FEP/unit_test/p38-19-24/rtf_fep_output/GMX_PROLIG_FEP`
