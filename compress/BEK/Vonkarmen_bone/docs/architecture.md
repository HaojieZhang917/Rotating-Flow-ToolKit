# 代码结构与公共 API

## 1. 仓库边界

Git 根目录为 `Rotating-Flow-ToolKit/`。`compress/BEK/Vonkarmen_bone/` 使用自己的 `Project.toml` 和 `Manifest.toml`，因此可以独立复现实验环境，但不包含嵌套 `.git`。

## 2. 核心模块

`src/RotatingDiskFlow.jl` 是唯一包入口，公开以下模块和类型：

| API | 用途 |
|---|---|
| `CurveConfig` | 一条中性曲线的全部参数 |
| `compute_neutral_curve` | 单次延拓，不执行外层失败恢复 |
| `NeutralCurveRunner` | 稳健延拓、批处理、结果校验 |
| `NeutralContinuation` | 特征向量匹配和中性根校正 |
| `LopezBaseflow` | Lopez 广义 Boussinesq 基本流 |
| `LopezStability` | Lopez 空间/时间稳定性矩阵和特征值 |
| `CRD_BF` | 经典基本流、Chebyshev 网格等低层工具 |
| `SutherlandMarching` | Sutherland 模型径向 BDF2 推进 |

依赖关系为：

```text
RotatingDiskFlow
├── NeutralCurveRunner
│   ├── LopezBaseflow
│   ├── LopezStability
│   ├── NeutralContinuation
│   └── CRD_STA compressible operators
└── SutherlandMarching
```

## 3. 兼容层

根目录中的下列文件只负责载入 `src/` 实现：

```text
CRD_STA.jl
LopezBaseflow.jl
LopezStability.jl
NeutralContinuation.jl
NeutralCurveRunner.jl
SutherlandMarching.jl
```

因此现有 notebook 无需立即修改。不要同时在同一个 Julia session 中反复 `include` 兼容文件和 `using RotatingDiskFlow`；修改源码后应重启 kernel，避免 `Main` 中保留旧方法。

## 4. 数据流

一条中性曲线的计算顺序为：

1. 根据 `CurveConfig.model` 计算 Lopez 或可压缩基本流。
2. 将基本流插值到同一个 Chebyshev 相似坐标网格。
3. 在给定 `(R,beta)` 上组装空间二次特征值问题。
4. 用前一步特征值和特征向量作为 IAR 初值。
5. 通过特征向量重叠选择同一物理模态。
6. 调整 `R` 使 `imag(alpha)=0`。
7. 沿 `beta` 延拓；失败时缩小步长，并设置有限重试和终止条件。
8. 写入 Tecplot `.dat`，再执行单调性、残差和点数检查。

## 5. 生成文件

运行结果与日志写入项目根下的数据目录，不写入 `src/`。默认忽略：

```text
neutral_curve_batch/
neutral_curve_integrated/
grid_independence/
sutherland_marching_data/
```

论文所需的冻结数据应在确认后另行放入有版本号的归档目录，并记录求解器 commit、Julia 版本、网格和容差。
