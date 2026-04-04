# Skills 文档使用指南

## 概述

`.claude/skills/` 目录包含 PRISM 项目的技能文档，用于指导 Claude Code 完成特定任务。这些文件已包含在 git 仓库中，可以随时访问。

## Skills 文件列表

### FEP 相关
- **`fep-system-debug.md`** - FEP 体系搭建与运行调试
- **`fep-analysis.md`** - FEP 结果分析和可视化（多估算器）
- **`fep-forcefield-integration.md`** - FEP 力场集成测试
- **`fep-mapping-viz.md`** - FEP 原子映射可视化

### 通用调试
- **`gaff-debug.md`** - GAFF 力场生成调试

## 使用说明

### 如何使用 Skills 文档

1. **在 Claude Code 中**：
   - Skills 文档会自动加载
   - 开始对话时，Claude 会参考相关的 skills 文档
   - 可以直接询问："使用 fep-system-debug skill 来..."

2. **在命令行中**：
   ```bash
   # 查看可用的 skills
   ls .claude/skills/

   # 阅读特定 skill
   cat .claude/skills/fep-system-debug.md
   ```

### 测试数据引用

Skills 文档中的代码示例引用了 `examples/` 目录，但该目录不在 git 中。要使用这些示例：

#### 选项 1：使用你自己的数据

```bash
# 替换示例中的路径
mol2_a = "examples/ligands/oMeEtPh.mol2"  # 示例路径
mol2_a = "my_data/ligand.mol2"           # 你的数据
```

#### 选项 2：创建示例数据目录（可选）

如果你想保留完整的示例，可以在本地创建 `examples/` 目录：

```bash
# 从测试目录复制示例数据（不提交到 git）
mkdir -p examples/ligands
cp tests/gxf/FEP/unit_test/oMeEtPh-EtPh/*.mol2 examples/ligands/ 2>/dev/null || true

# 添加到 .gitignore（如果还没有）
echo "examples/" >> .gitignore
```

## LSP 工具说明

所有 skills 文档都包含 LSP 工具使用说明：

- **LSP 工具优先**：使用 `mcp__cclsp__` 工具进行代码导航
- **可用性**：仅在 Claude Code 工具上下文中可用
- **回退方案**：LSP 不可用时使用 `Grep` 工具

## 代码搜索工具对比

| 工具 | 优势 | 使用场景 |
|------|------|----------|
| `mcp__cclsp__find_definition` | 精确定位符号定义 | 查找类/函数定义 |
| `mcp__cclsp__find_references` | 查找所有引用 | 重构影响分析 |
| `mcp__cclsp__find_workspace_symbols` | 全局符号搜索 | 发现相关功能 |
| `Grep` | 文本搜索，广泛可用 | 快速关键词搜索 |

## 更新 Skills 文档

当更新 skills 文档时：

1. **避免硬编码测试路径**：使用 `examples/` 而不是 `tests/`
2. **保持代码示例简洁**：只展示核心逻辑
3. **添加必要说明**：解释如何替换为用户自己的数据
4. **更新 LSP 说明**：确保所有文档都有 LSP 工具章节

## Git 跟踪状态

Skills 文档现在已经被 git 跟踪：

```bash
# 查看跟踪状态
git status .claude/skills/

# 查看已跟踪的文件
git ls-files .claude/skills/
```

## 相关文档

- **项目指南**: `CLAUDE.md`
- **FEP 教程**: `/home/gxf1212/data/work/PRISM-Tutorial`
- **测试数据**: `tests/gxf/FEP/unit_test/`（本地，不在 git 中）
