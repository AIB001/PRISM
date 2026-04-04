# PR提交工作流程

记录PR提交和描述文件管理的完整流程和要求。

## PR描述文件管理

### 固定文件名规则
**唯一PR描述文件**：`/tmp/pr_body_expanded.md`

**严格规则**：
- ✅ 保留：`/tmp/pr_body_expanded.md` （固定使用，永不改名）
- ❌ 删除：所有其他PR相关的.md文件（pr_body_update.md, pr_description.md等）

**重要说明**：
- 这个文件包含PR #4的**完整描述**（整个FEP工作流程）
- 每次有新提交时，需要更新这个文件，包含**所有累积的更改**
- 用这个文件更新远程PR，让同事审阅时能看到完整的变更历史
- 写长一点没关系，重要的是**完整和详细**

### PR描述清理
```bash
# 删除所有其他PR文件
rm -f /tmp/pr_body_update.md /tmp/pr_description.md /tmp/prism_format_support.md

# 验证只剩下一个文件
ls -la /tmp/pr*.md
```

## PR提交流程

### 1. 检查未提交更改
```bash
cd /data2/gxf1212/work/PRISM
git status --short
```

### 2. 暂存所有更改
```bash
git add -A
git status --short
```

### 3. 提交更改
```bash
git commit -m "$(cat <<'EOF'
[module] Brief description of changes

- Main change 1
- Main change 2
- Main change 3

Co-Authored-By: Claude Sonnet 4.6 <noreply@anthropic.com>
EOF
)"
```

**Commit消息格式**：
- 使用 `[module]` 前缀（如 [fep], [builder], [analysis]）
- 简洁描述主要更改
- 列出关键变更点
- 包含Co-Authored-By标签

### 4. 推送到远程
```bash
git push origin gxf
```

### 5. 更新PR描述（关键步骤）
```bash
gh api repos/AIB001/PRISM/pulls/4 -X PATCH -F body=@/tmp/pr_body_expanded.md
```

**重要**：
- 每次推送新提交后，都要用 `/tmp/pr_body_expanded.md` 更新PR描述
- PR描述应该包含**完整的FEP工作流程描述**，不只是最新改动
- 体现从main分支到现在的**所有变更对比**
- 方便同事审阅时了解完整的工作内容

### 6. 验证PR更新
```bash
gh pr view 4 --json title,state,updatedAt --jq '{title, state, updated_at}'
```

## PR描述内容要求

### 必须包含的章节

```markdown
## Summary
简要概述（1-2句话）

## Main Features
### FEP workflow
- 完整的FEP工作流程特性
- 原子映射、混合拓扑等

### FEP analysis
- 分析方法、Bootstrap等

### Builder and force-field integration
- 力场集成、RTF支持等

### Analysis and architecture refactoring
- 架构重构、模块拆分等

## Command-line usage
实际的CLI使用示例（不是测试文件路径）

## Python API usage
Python接口使用示例

## Current validation status
### Validated end-to-end
- 完整验证的测试系统

### Validated for modeling / scaffold generation
- 建模验证的测试系统

### Additional validated paths
- 其他验证路径

## Review guidance
审阅顺序和建议

## Notes
注意事项和限制

## Files Changed
按类别列出的所有更改文件
```

### 内容规则
1. **完整详细**：包含所有主要功能、验证结果、使用方法
2. **方便审阅**：同事能看到完整的变更内容，不只是最新改动
3. **体现对比**：展示从main到现在的完整变更历史
4. **保持更新**：每次新提交都要更新PR描述，包含累积的所有更改
5. **实际代码**：根据实际代码描述，不虚构内容
6. **重点突出**：FEP模块是主题内容，要详细描述所有FEP特性
7. **多写features**：详细列出功能点，每个功能类别都要展开说明
8. **测试文件**：不要告诉同事测试文件在哪，只说有comprehensive test systems即可

### 必须保留的内容
- ✅ Main Features（所有主要功能类别）
- ✅ Command-line usage（实际CLI接口）
- ✅ Python API usage（Python接口）
- ✅ Current validation status（完整验证结果）
- ✅ Review guidance（审阅指南）
- ✅ Files Changed（按类别组织的文件列表）
- ✅ Notes（注意事项）

## 当前PR状态

### PR #4: "Integrate FEP workflow with analysis and builder refactors"
- **状态**: OPEN
- **分支**: gxf → main
- **URL**: https://github.com/AIB001/PRISM/pull/4
- **主题内容**: FEP工作流程（atom mapping, hybrid topology, scaffold generation, analysis, force-field integration）
- **最新添加**: RTF force field支持

### 常用命令
```bash
# 查看PR状态
gh pr view 4

# 查看PR提交数
gh pr view 4 --json commits --jq '.commits | length'

# 查看PR文件更改
gh pr view 4 --json changedFiles --jq '.changedFiles'

# 检查PR评论
gh pr view 4 --json comments --jq '.comments | length'

# 更新PR描述
gh api repos/AIB001/PRISM/pulls/4 -X PATCH -F body=@/tmp/pr_body_expanded.md
```

## 重要提醒

1. **文件名固定**：永远使用 `/tmp/pr_body_expanded.md`，不准改文件名
2. **清理旧文件**：每次更新PR前删除其他PR描述文件
3. **Commit格式**：使用 `[module] description` 格式，不要太长
4. **完整描述**：PR描述要完整详细，包含所有累积的更改
5. **每次更新**：每次推送新提交后，都要更新PR描述
6. **体现对比**：PR描述体现从main到现在的完整变更对比

## Git工作流程

### 分支管理
- **开发分支**: gxf
- **主分支**: main
- **PR目标**: main

### 提交流程
1. 在gxf分支上开发和测试
2. 提交更改到gxf分支
3. 推送到origin/gxf
4. **更新PR #4描述**（包含所有累积的更改）
5. 等待review和合并

### 不要做的事
- ❌ 不要直接推送到main分支
- ❌ 不要创建新的PR（使用PR #4）
- ❌ 不要随意更改PR描述文件名
- ❌ 不要在commit消息中写太长的描述
- ❌ 不要只描述最新改动，要包含所有累积更改
