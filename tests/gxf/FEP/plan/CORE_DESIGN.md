# PRISM-FEP 核心算法设计

本文档详细描述原子映射和双拓扑构建的核心算法、数据结构和实现逻辑。

## 1. 原子映射模块 (`fep/core/mapping.py`)

### 1.1 功能需求
- **输入**: 两个配体的原子列表（包含坐标、元素、电荷、力场类型）
- **输出**: AtomMapping 对象，将原子分类为 common/transformed/surrounding
- **移植来源**: FEbuilder/setup/make_hybrid.py 的 `Ligand.get_correspondence()` 方法

### 1.2 核心数据结构

```python
from dataclasses import dataclass
from typing import List, Tuple, Optional
import numpy as np

@dataclass
class Atom:
    """
    原子对象
    
    设计说明:
    - name: 原子唯一标识符，如 'C1', 'H2', 'N3'
    - element: 元素符号，用于匹配时的类型检查 (C, H, N, O 等)
    - coord: numpy 数组 [x, y, z]，单位 Å
    - charge: 部分电荷，来自力场参数化
    - atom_type: 力场原子类型 (如 GAFF 的 'ca', 'ha', 'na')
    - index: 原子索引，用于后续引用
    """
    name: str
    element: str
    coord: np.ndarray
    charge: float
    atom_type: str
    index: int

@dataclass
class AtomMapping:
    """
    原子映射结果
    
    设计说明:
    这是整个 FEP 模块的核心数据结构，所有后续处理都基于这个分类
    
    分类逻辑:
    1. common: 两个配体都存在，且位置、元素、类型、电荷都接近
             在双拓扑中只保留一份，A/B 状态完全相同
             
    2. transformed_a/b: 只在一个配体中存在的原子
                        这些原子在状态 A 或 B 中会变为 dummy 类型（由力场提供，例如 `DUM_*`）
                        
    3. surrounding_a/b: 位置匹配但力场类型或电荷差异大的原子
                        这些是"过渡区域"，化学环境发生了变化
    """
    common: List[Tuple[Atom, Atom]]     # 共用原子对
    transformed_a: List[Atom]            # 配体A独有
    transformed_b: List[Atom]            # 配体B独有
    surrounding_a: List[Atom]            # 配体A（类型不同）
    surrounding_b: List[Atom]            # 配体B（类型不同）
```

### 1.3 距离匹配算法

**算法思路**: 贪心匹配 + 距离阈值

```python
class DistanceAtomMapper:
    """
    基于距离的原子映射器
    
    参数:
    - dist_cutoff: 距离阈值 (Å)，默认 0.6
      - 较小值 (0.3-0.5): 更严格，适合非常相似的分子
      - 较大值 (0.8-1.0): 更宽松，适合结构差异大的分子
    
    - charge_cutoff: 电荷差异阈值，默认 0.05
      - 超过此值的匹配原子归为 surrounding
    """
    
    def map(self, ligand_a: List[Atom], ligand_b: List[Atom]) -> AtomMapping:
        """
        主映射函数
        
        算法流程:
        
        Step 1: 距离匹配
        -----------
        遍历 ligand_a 的每个原子，在 ligand_b 中找到：
        - 距离 < dist_cutoff
        - 元素类型相同
        - 没有被其他原子匹配过
        的第一个原子，形成匹配对
        
        实现:
        - 使用 matched_b 集合记录已匹配的原子索引
        - 一旦找到匹配就 break，避免一对多匹配
        - 时间复杂度: O(N_a * N_b)，对小分子完全可接受
        
        Step 2: 识别 transformed atoms
        ----------------------------
        未出现在匹配对中的原子就是 transformed
        - 使用集合差集: 所有原子 - 匹配的原子
        
        Step 3: 识别 surrounding atoms
        -----------------------------
        遍历所有匹配对，检查：
        - atom_type 是否相同
        - charge 差异是否超过 charge_cutoff
        - 任一条件不满足 → surrounding
        
        Step 4: 手性检查（可选）
        -----------------------
        检查中心碳原子周围是否有 3 个共用原子 + 1 个变化原子
        可能需要标记手性翻转
        
        返回:
            AtomMapping: 映射结果
        """
        pass
```

### 1.4 边界情况处理

**问题1: 多个相同元素原子距离都 < cutoff**
- 解决: 选择第一个，避免重复匹配

**问题2: 手性中心的识别**
- 解决: 通过连接关系判断，需要原子的 bond 信息

**问题3: 浮点数比较的精度问题**
- 解决: 使用 charge_threshold 而不是直接 == 比较

---

### 1.5 pmx 蛋白突变输入（参考体系）

当前已跑通的参考体系使用 pmx 生成的蛋白突变拓扑（`tests/gxf/FEP/ref/system/newtop.top`）：
- 原子 A/B 状态已由 pmx 写入 `[ atoms ]`，无需再做距离映射
- `[ atoms ]` 包含 `typeB/chargeB/massB`，必须保留 `massB`
- dummy 类型由 `amber14sbmut.ff` 提供（`DUM_*`），不是固定字符串

因此在 PRISM-FEP 设计中应支持两条路径：
1. **配体映射路径**：使用本节的距离映射算法生成 A/B 拓扑  
2. **pmx 直读路径**：直接读取 A/B 拓扑（跳过映射），用于蛋白突变体系

---

## 2. 双拓扑构建模块 (`fep/core/dual_topology.py`)

### 2.1 功能需求
- **输入**: AtomMapping + 两个配体的力场参数（键、角、二面面角）
- **输出**: HybridAtom 列表（带 typeB/chargeB）
- **移植来源**: FEbuilder/setup/make_hybrid.py 的 `Ligand.combine()` 方法

### 2.2 核心数据结构

```python
from dataclasses import dataclass
from typing import Dict, List, Optional

@dataclass
class HybridAtom:
    """
    杂化原子 - 包含 A/B 两种状态
    
    设计说明:
    直接映射到 GROMACS ITP 文件的一行
    
    字段解释:
    - name: 原子名称
      共用原子: 原始名称 (如 C1, H2)
      变化原子: 加标签 (如 C1A, C1B) - 便于调试
      
    - index: 原子索引
      从 1 开始连续编号
      用于 bonds/angles/dihedrals 中的引用
      
    - state_a_type/state_b_type: A/B 状态的力场类型
      共用原子: 两者相同 (如都是 'ca')
      变化原子 A: state_a_type 正常，state_b_type = dummy 类型
      变化原子 B: state_a_type = dummy 类型，state_b_type 正常
      surrounding: 各自保留不同的类型
      
    - state_a_charge/state_b_charge: A/B 状态的电荷
      逻辑同 type
      
    - element: 元素符号
      用于质量查找
      
    - mass/mass_b: 原子质量
      pmx 蛋白突变拓扑中 `[ atoms ]` 含 `massB`，需要保留 A/B 两态质量
    """
    name: str
    index: int
    state_a_type: str
    state_a_charge: float
    state_b_type: Optional[str] = None      # None = 共用原子
    state_b_charge: Optional[float] = None
    element: str = ""
    mass: float = 0.0
    mass_b: Optional[float] = None  # pmx 拓扑需要 massB 时填充
```

### 2.3 构建流程

**原子顺序设计**（重要！）:
```
1. common atoms (共用原子)
   - 只保留一份
   - A/B 状态完全相同
   - 电荷根据 charge_strategy 处理
   
2. transformed atoms from ligand A
   - 加 'A' 后缀 (便于调试)
   - 状态 A: 正常参数
   - 状态 B: dummy 类型原子 (type=由力场提供, charge=0)
   
3. transformed atoms from ligand B
   - 加 'B' 后缀
   - 状态 A: 虚拟原子
   - 状态 B: 正常参数
   
4. surrounding atoms
   - 不加后缀
   - A/B 状态各自保留不同的类型和电荷
```

**伪代码**:
```python
class DualTopologyBuilder:
    def build(self, mapping: AtomMapping, 
                 params_a: Dict, params_b: Dict) -> List[HybridAtom]:
        """
        构建双拓扑的主流程
        
        流程:
        1. _build_atoms(): 构建杂化原子列表
        2. _map_parameters(): 映射参数到 hybrid atoms
        3. 返回 hybrid_atoms
        """
        self._build_atoms()
        self._map_parameters()
        return self.hybrid_atoms
    
    def _build_atoms(self):
        """
        构建杂化原子列表
        
        索引管理:
        - 从 1 开始 (GROMACS 约定)
        - 连续递增
        - 需要记录原始索引 → hybrid 索引的映射
        """
        index = 1
        dummy_type = "DUM_*"  # 实际应按原子类型映射（pmx/amber14sbmut 中为 DUM_*）
        
        # 1. 共用原子
        for a, b in mapping.common:
            charge = self._resolve_charge(a.charge, b.charge)
            hybrid_atoms.append(HybridAtom(
                name=a.name, index=index,
                state_a_type=a.atom_type, state_a_charge=charge,
                state_b_type=a.atom_type, state_b_charge=charge,  # 共用
                element=a.element, mass=params_a['masses'][a.atom_type]
            ))
            index += 1
        
        # 2. Transform atoms A
        for a in mapping.transformed_a:
            hybrid_atoms.append(HybridAtom(
                name=f"{a.name}A", index=index,
                state_a_type=a.atom_type, state_a_charge=a.charge,
                state_b_type=dummy_type, state_b_charge=0.0,  # B状态 dummy
                element=a.element, mass=params_a['masses'][a.atom_type]
            ))
            index += 1
        
        # 3. Transform atoms B
        for b in mapping.transformed_b:
            hybrid_atoms.append(HybridAtom(
                name=f"{b.name}B", index=index,
                state_a_type=dummy_type, state_a_charge=0.0,  # A状态 dummy
                state_b_type=b.atom_type, state_b_charge=b.charge,
                element=b.element, mass=params_b['masses'][b.atom_type]
            ))
            index += 1
        
        # 4. Surrounding atoms
        for a, b in zip(mapping.surrounding_a, mapping.surrounding_b):
            hybrid_atoms.append(HybridAtom(
                name=a.name, index=index,
                state_a_type=a.atom_type, state_a_charge=a.charge,
                state_b_type=b.atom_type, state_b_charge=b.charge,  # 各自保留
                element=a.element, mass=params_a['masses'][a.atom_type]
            ))
            index += 1
    
    def _resolve_charge(self, charge_a: float, charge_b: float) -> float:
        """
        处理共用原子的电荷
        
        三种策略:
        - mean: 取平均 (默认，平滑过渡)
        - ref: 使用参考配体的电荷 (保持参考静电环境)
        - mut: 使用突变配体的电荷 (保持突变静电环境)
        """
        if self.charge_strategy == 'ref':
            return charge_a
        elif self.charge_strategy == 'mut':
            return charge_b
        return (charge_a + charge_b) / 2.0
    
    def _map_parameters(self):
        """
        映射参数到 hybrid atoms
        
        关键点: GROMACS 不需要像 NAMD 那样合并参数值！
        
        只需将原始配体的 bonds/angles/dihedrals 原子索引
        映射到 hybrid_atoms 的索引即可
        
        GROMACS 会根据原子的 typeB/chargeB 自动选择 A 或 B 状态的参数
        
        实现:
        1. 建立原子名称 → hybrid_atom_index 的映射表
        2. 遍历配体A的所有参数
        3. 将原子名称映射到 hybrid_atoms 的索引
        4. 对配体B重复步骤
        5. GROMACS 会自动处理 A/B 状态的参数选择
        """
        # 建立映射表: atom_name -> hybrid_index
        name_to_index = {}
        for atom in self.hybrid_atoms:
            name_to_index[atom.name] = atom.index
        
        # TODO: 实现参数索引映射
        pass
```

### 2.4 关键设计决策

**为什么不需要参数合并？**
- NAMD: 两个独立的分子结构，需要合并参数
- GROMACS: 一个分子带 A/B 状态标记，自动选择状态

**为什么使用 dummy 原子？**
- 保持系统原子数一致（避免 GROMACS 报错）
- dummy 类型由力场定义（pmx/amber14sbmut 为 `DUM_*`），电荷通常为 0
- 质量不必为 0（pmx 拓扑含 `massB`），非键相互作用由 dummy 类型屏蔽

---

## 3. GROMACS vs NAMD 的关键差异

### 3.1 参数处理

**NAMD (FEbuilder)**:
- 需要参数合并 (merge_prm.py)
- 两个配体的键/角/二面角参数需要去重和合并
- 复杂的冲突处理逻辑

**GROMACS**:
- **不需要参数合并**
- 使用 typeB/chargeB 列，每个原子有两套状态
- bonds/angles/dihedrals 只需确保原子索引正确
- GROMACS 自动根据原子的 type/charge 选择状态

### 3.2 拓扑文件格式

**NAMD**: 
- .rtf (topology) + .prm (parameters)
- 需要合并到 hybrid.rtf/hybrid.prm
- 参数值需要合并处理

**GROMACS**:
- .itp (include topology)
- typeB/chargeB 格式，更简洁
- 不需要参数值合并，只需原子索引映射

### 3.3 实现影响

PRISM-FEP **不需要实现** FEbuilder 的 `merge_prm.py` 逻辑：
- ❌ 不需要复杂的参数值合并
- ✅ 只需要将原子索引映射到 hybrid_atoms
- ✅ GROMACS 会自动处理 A/B 状态的参数选择

---

## 4. 力场支持

### 4.1 支持的力场

| 力场 | CLI 参数 | 说明 | 示例文件 |
|------|----------|------|----------|
| GAFF | `--force-field gaff` | 通用小分子力场 | ligand.mol2 |
| GAFF2 | `--force-field gaff2` | GAFF 改进版 | ligand.mol2 |
| OpenFF | `--force-field openff` | 开放力场，基于 SMIRKS | ligand.sdf |
| **CGenFF** | `--force-field cgenff` | CHARMM 通用力场 | ligand.rtf + ligand.prm |
| OPLS-AA | `--force-field oplsaa` | OPLS 全原子力场 | ligand.mol2 |

### 4.2 CGenFF 支持

**特点**:
- PRISM 已有 `prism/forcefield/cgenff.py` 模块
- 支持 CHARMM-GUI 输出格式 (`.rtf` + `.prm`)
- 适合含药物分子的系统

**测试数据**:
- 位置: `tests/gxf/FEP/test/hif2a/42-38/`
- 格式: CHARMM-GUI 输出
- 包含: `toppar/lig.rtf`, `toppar/lig.prm`

---

## 5. 参考资料
- FEbuilder: `/home/gxf1212/data2/work/make_hybrid_top/FEbuilder/`
- GROMACS FEP: https://manual.gromacs.org/current/fep.html
- pmx: https://github.com/deGrootLab/pmx
- PRISM CGenFF: `/data2/gxf1212/work/PRISM/prism/forcefield/cgenff.py`
