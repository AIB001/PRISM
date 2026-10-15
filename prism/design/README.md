# PRISM-Design M0

将已审阅的 `prism-design-annotation/0.1` DesignSpec 转成真实局部生成任务，保留全部原始采样，再独立验证生成结果。第一轮支持单连接点片段替换和小片段生长。输入是 protein PDB、pocket PDB、带原子映射的参考配体 SDF、DesignSpec 和绑定的证据/审阅文件。

M0 的三维采样由已有 FLOWR.root 模型完成。新增的 [LLM 编辑规划器](#llm-编辑规划器) 可调用显式配置的在线/本地服务，输出并执行受限分子图编辑；该路线不需要微调。可选 M1 训练另见文末。通过结构与描述符检查不代表亲和力改善。

## 安装与模型

在 PRISM 的 Python 环境安装可选依赖：

```bash
pip install -e '.[design]'
python -m prism.design --help
# 安装后也可使用 prism design ... 或 prism-design ...
```

FLOWR.root 使用独立 Python 环境，避免改变 MD 环境中的 PyTorch/RDKit。当前适配的是 [FLOWR.root 官方源码](https://github.com/jule-c/flowr_root) 的精确版本：

- commit: `739a33c2aab74d074dd24605b9a2e34dc809ae9c`
- 本次测试权重：`flowr_root_v2.2.ckpt`
- 本地下载快照 SHA-256：`b818f41dc12ffb6bc558bb0ad997055581e07cd9e49dcac1b794ed9993c46e4c`。这是本次文件的校验值，并非发布者提供的签名。
- 模型环境按该版本的 `uv.lock` 安装：在 checkout 中执行 `uv sync --extra gpu --locked`。
- 模型权重的官方下载位置见该版本 README；权重不包含在 PRISM 安装包中。

复制 [配置模板](../configs/design.flowr_root.example.json)，填写源码、隔离环境 Python 和 checkpoint 的本地路径。路径可相对配置文件；Python 路径必须保留虚拟环境入口。适配器检查源码 commit、已跟踪文件是否改动及权重哈希，然后在每次运行中检查 checkpoint 的片段训练声明、原生固定掩码和实际预处理后的原子顺序/坐标。

```bash
python -m prism.design backend-check --backend-config /path/to/backend.json
```

`backend-check` 只校验绑定的文件；真正的运行还会进行原生能力探测。不支持的旧版生成模型会报错，不会自动降级成无约束生成。

## 一次运行

```bash
python -m prism.design validate \
  --spec /path/to/round01/complexes/184l/designspec.json --for-execution

python -m prism.design compile \
  --spec /path/to/round01/complexes/184l/designspec.json \
  --output /path/to/compiled_preview

python -m prism.design run \
  --spec /path/to/round01/complexes/184l/designspec.json \
  --backend-config /path/to/backend.json \
  --output /path/to/new_run --num-samples 16 --seed 20261013

python -m prism.design verify \
  --spec /path/to/round01/complexes/184l/designspec.json \
  --candidate /path/to/new_run/backend/raw_slots/slot_0001.sdf \
  --output /path/to/independent_check.json
```

`run` 自行编译，请使用新的输出目录；不会覆盖旧运行。`compile` 可单独用于检查编译结果。`verify` 只评估给定的原始候选，不调用生成模型。验证失败返回退出码 2；采样正常结束但没有合格候选，仍属于完成的运行，实际通过数见 `counts`。

DesignSpec 的相对引用以带哈希的 `selection_manifest.json` 所在标注根目录为基准。程序向上查找根目录，也支持 `--input-root`。源标注仍保持 draft 和原来的 `executable_in_current_prism: false`，编译器创建单独的 `prism-design-compiled/1` 产物，不篡改历史标注。

第一轮队列批处理：

```bash
python -m prism.design pilot \
  --queue /path/to/round01/pilot_queue.json \
  --backend-config /path/to/backend.json \
  --output /path/to/new_pilot --seed 20261013
```

只运行队列的 `first_smoke_cases`，预算来自每项的 `proposed_candidate_slots`。原第一轮队列是 6 × 16 = 96；不同案例使用基础 seed 加队列序号。输出包含逐案例报告及 `pilot_summary.json` / `pilot_summary.csv`。

## 执行和检查规则

1. 输入绑定：校验 DesignSpec schema、文件 SHA-256、化学状态、蛋白/口袋坐标系、原子映射、审阅决定与源质量标记。只有 `candidate_for_m0_pilot` 可编译执行。金属近接、明显原始碰撞、未解决源质量问题等需要先形成新的审阅版本。
2. 任务编译：将稳定 atom map 转成后端的重原子 0 起始索引。后端参考去掉 H 与 atom map 标签，并检查重原子化学结构、顺序和坐标不变。口袋只去掉 CONECT 行，原文件另存。原子映射标签在化学/手性识别前移除，避免标签制造假手性。
3. 条件采样：替换使用原生 `substructure_inpainting`，参数传入可编辑原子；生长固定全部参考重原子，默认新增 2 个原子。C–H 方向形成软生长先验，不能保证精确连接点。原生 `final_inpaint` 开启，默认价态修复、额外重试及最终有效性/去重过滤关闭。每个模型解码输出都进入原始记录。
4. 独立检查：在原受体坐标系寻找满足图和连接约束的核心映射；不对齐候选配体。检查元素、键级、固定手性、电荷、单一连接点、新片段大小、新增未指定手性及核心 RMSD；存在无法判断项时不会算通过。当前片段替换不支持在已指定手性中心切割，或跨越双键立体定义的边界。
5. 筛选：按 DesignSpec 记录硬约束和可计算软目标，并检查候选确实不同于参考配体、最近蛋白距离至少 2 Å、键长/共价半径和比例处于 0.60–1.45。再运行现有 PRISM 完整 QC；当前 M0 自动候选要求 QC 为 `PASS`。警告和失败结果仍保留，需另行复核。

原子子图匹配最多检查 2048 个映射。搜索被截断且没有找到合法核心时会记录不确定性。生长连接点和替换后电荷主要依赖后验筛选，不能将 native conditioning 解释为所有约束均已保证。

当前候选按每个案例内的 canonical isomeric SMILES 去重，保留首次通过的姿势；所有其他原始姿势仍在记录中。计数分开报告 `requested`、`raw_slots`、`decoded`、`not_attempted`、`hard_constraints_passed`、`design_selection_passed`、`design_and_qc_passed`、`unique_selected`。后端崩溃或超时后的未执行位置标记 `NOT_ATTEMPTED`，不会伪装成实际生成失败或成功。

## 产物与现有 PRISM 的衔接

每个案例输出：

| 文件/目录 | 内容 |
| --- | --- |
| `compiled_design.json`, `inputs/` | 编译任务、索引、输入副本、转换记录、哈希 |
| `backend/launch.json`, `process.json` | 实际命令、后端版本、进程状态 |
| `backend/capability_probe.json` | 实际原生掩码及预处理验证 |
| `backend/attempts.jsonl`, `raw_slots/` | 全部实际解码尝试和原始 SDF |
| `backend/runtime_summary.json`, `backend.log` | 后端统计、耗时、GPU 内存及错误 |
| `validation/slot_*.json` | 每个计划位置的独立检查、QC、拒绝原因 |
| `results.json`, `slots.csv`, `status.json` | 完整记录及统计 |
| `candidates/` | 通过筛选的原始分子及独立加氢副本 |
| `manifest.json` | 兼容现有 `prism prepare-md` 的已筛选候选清单 |

加氢使用现有 `hydrogenate_record`，保留重原子位置，不枚举 pH 相关质子化态/互变异构体。原始 SDF 始终单独保留。

```bash
prism prepare-md /path/to/new_run --output /path/to/md_inputs --only-pass
```

这一步导出 MD 输入。力场参数化、蛋白状态准备、最小化、MD 和 ABFE/FEP 应在后续阶段执行；本轮结构筛选没有运行这些计算。

## 测试与适用范围

`tests/test_design.py` 使用合成分子和合成坐标，不依赖 PDBbind 或 GPU。已有生成模块的单元测试可在隔离 GPU 的环境运行：

```bash
CUDA_VISIBLE_DEVICES='' python -m pytest \
  tests/test_design.py tests/test_generation.py tests/test_generation_quality.py \
  tests/test_generation_handoff.py tests/test_generation_weights.py \
  tests/test_generation_calibration.py tests/test_mcp_generation.py -q
```

真实 GPU 验证必须另行使用上述 `pilot` 命令。模型是已有 PDBbind 相关训练权重，本轮 PDBbind 案例属于工程可行性检查，未经训练集重叠排查，不能作为独立泛化基准。扩展数据前应先做靶点/序列/骨架拆分，再建立实测匹配优化对；当前标注不具备“优化后活性更好”的训练标签。

## M1 与数据整理

新增可移植的结构描述/片段重建数据格式、逐条审阅账本、肽断键证据和原生文本条件微调入口。详见 [M1 训练架构](M1_TRAINING_ARCHITECTURE.md) 及 [配置示例](../configs/design.m1.example.json)。

`prism design dataset-check` 检查数据并分组；`preflight` 在模型环境逐条验证；`train` 在同一环境训练/恢复。M1 增量权重可由现有 `run` 的 backend 配置加载，候选继续经过独立硬约束验证。

## LLM 编辑规划器

`plan-prepare → plan → plan-check` 将模型提出的原子/键编辑落实为可重放的分子图。它使用现有 DesignSpec 的源文件、审阅范围和硬约束，不需要 M1 checkpoint，也不改变 M1 的文本编码器。在线入口兼容 DeepSeek、Qwen 等 chat-completions 服务；离线入口可以导入本助手或其他模型根据同一证据包写出的 JSON。

### 准备结构证据

```bash
python -m prism.design plan-prepare \
  --spec /path/to/round01/complexes/184l/designspec.json \
  --goal '保留苯环核心，探索更紧凑的中性侧链；报告不确定性。' \
  --output /path/to/request
```

可以用 `--goal-file goal.txt` 代替命令行文字。省略目标时沿用 DesignSpec 的 prompt。输入必须通过现有 `for_execution=True` 检查，规划器只在该版本已审阅的固定核心、可编辑区域和连接点内工作。需要改变范围时应先建立相应的新 DesignSpec，规划器不会自行更改审阅级别。

输出 `context.json`、`messages.json`、`planning_request.json`。证据包括实际 atom maps、带编号 SMILES、元素/电荷/键/立体标签、计算描述符、每原子最近蛋白距离和 4 Å 内残基接触。每条事实有可引用的 evidence ID。距离由输入坐标重算；没有把距离邻近标成氢键、相互作用能或溶剂暴露。请求绑定源文件、证据和消息的哈希，执行前与重放时重新检查。

### 调用规划模型

设置 `DEEPSEEK_API_KEY` 环境变量，然后使用 [DeepSeek 配置](../configs/design.planner.deepseek.example.json)：

```bash
python -m prism.design plan \
  --request /path/to/request/planning_request.json \
  --llm-config prism/configs/design.planner.deepseek.example.json \
  --output /path/to/plan
```

也提供 [本地 Qwen 配置](../configs/design.planner.qwen_local.example.json)，指向用户已启动的兼容服务；模型名须与该服务的部署名一致。本模块不下载模型、不启动服务、不扫描或保存密钥。使用远程配置运行 `plan` 会向指定服务发送 `messages.json` 中的目标、分子和结构证据；可先检查该文件，或使用本地服务。

DeepSeek 示例采用官方当前的 `deepseek-flash` 名称。官方文档在 2026-10-13 查询时说明，旧 `deepseek-v4-flash` 别名已转向 V4.1 Flash；不能将旧别名当成冻结的模型版本。日志保存请求/返回的模型名、response ID 与 usage，模型别名本身不保证权重不变。参考 [模型说明](https://api-docs.deepseek.com/quick_start/pricing/) 和 [Chat Completions 协议](https://api-docs.deepseek.com/api/create-chat-completion/)。

默认最多 1 次调用，`max_attempts` 可设为 1–3。只有已完整返回但未通过格式/图检查的响应可在预算内修订；网络、凭证、截断和工具调用响应不会自动重试。每次 HTTP 请求超时上限为 60 秒。`request_options` 仅接受明确列出的思考模式参数，实际支持范围由服务决定；示例关闭 DeepSeek 思考模式以限制首轮预算，效果需要另行比较。任何响应都不能触发任意工具或 shell。

离线导入使用相同的检查路径：

```bash
python -m prism.design plan \
  --request /path/to/request/planning_request.json \
  --response /path/to/model_response.json \
  --output /path/to/plan
```

离线响应必须包含请求中的 `context_sha256`。其 schema 见 [edit-proposal-v1.json](schemas/edit-proposal-v1.json)，也包含在 `messages.json` 中。导入记录明确标记为 `imported_response`，不证明该文件来自哪个模型，不冒充一次在线调用。

### 编辑执行和验收

v1 支持 `mutate_atom`、`delete_atom`、`add_atom`、`add_bond`、`delete_bond`、`change_bond`。所有操作使用 atom maps；新增 ID 必须超过所有源重原子 ID，删除后也不能复用。普通隐式氢随新价态重算，不自动修复重原子图。暂不支持芳香图重写、显式氢操作或指定立体化学及其定义邻居的修改；新增未指定四面体/E-Z 立体化学会被拒绝。

模型还应提供简短设计解释、相关 evidence IDs、不确定性和未满足目标；可以明确 `abstain`。软件检查证据 ID 是否存在且涉及编辑原子，不能证明解释中每个语义判断正确。`goal_satisfaction_verified` 和 `semantic_claims_reviewed` 保持为 false，现有可计算约束逐项报告。

```bash
python -m prism.design plan-check \
  --plan /path/to/plan/edit_plan.json \
  --output /path/to/replay.json

python -m prism.design plan-check \
  --plan /path/to/plan/edit_plan.json \
  --candidate /path/to/posed_candidate.sdf \
  --output /path/to/pose_check.json
```

第一条重放实际编辑，检查引用文件、模型响应与保存结果一致。第二条先要求候选的化学/立体身份等于计划中的产物，再交给已有 `verify_candidate` 检查核心、连接、电荷、姿势和描述符条件。

`GRAPH_VALIDATED` 仅表示分子图和图级硬约束通过。软目标可失败，姿势约束标为 `not_evaluated_requires_3d`。`edited_graph.smi` 不包含三维坐标；源配体坐标不会被冒充为编辑产物的姿势。`ABSTAINED` 不产生分子，`REJECTED` 保留失败响应和原因。CLI 的拒绝结果退出码为 2。

| 输出 | 含义 |
| --- | --- |
| `attempt_*.txt`, `attempts.json` | 原始响应、逐次结果、服务元数据和失败原因 |
| `effective_llm_config.json` | 在线调用的实际参数快照，保存密钥环境变量名，不含密钥值 |
| `proposal.json` | 通过校验的可执行编辑操作或弃权 |
| `edited_graph.smi` | 实际执行得到的带编号分子图，无三维姿势 |
| `graph_validation.json` | 图级约束、计算描述符、尚未检查的姿势条件 |
| `edit_plan.json`, `status.json` | 绑定请求、产物和来源的状态记录 |

本次交付完成“证据 → LLM 方案 → 图编辑 → 独立重放检查”。现有 FLOWR 局部采样不保证精确生成这个图，本版没有新增 `run --plan` 来假装具备该能力。后续需用接受固定分子图的构象/对接流程获得姿势，再通过上述候选检查和完整 QC；当前计划不能直接进入 MD。逆合成评估、检索案例库和分轮设计反馈仍是后续模块。

测试：`python -m pytest tests/test_design_planner.py tests/test_design.py -q`。测试覆盖真实本地 HTTP 协议，但服务器响应是合成夹具，不能作为 DeepSeek/Qwen 的化学能力评测。
