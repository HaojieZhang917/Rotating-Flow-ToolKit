# NeutralCurveRunner 使用说明

`NeutralCurveRunner.jl` 用于计算旋转圆盘边界层的空间中性曲线，当前支持：

- Lopez 广义 Boussinesq 模型；
- 完全可压缩模型；
- 关闭扰动物性反馈的可压缩模型；
- 同时冻结基本流物性的可压缩模型；
- 基于特征向量重叠的模态连续追踪；
- 自适应 `beta` 步长、失败重试和 Type-I 回退搜索；
- 单工况、串行批处理和四进程并行批处理；
- 结果文件自动验证和断点续算。

## 1. 包依赖与载入

推荐通过 Julia 包入口载入：

```julia
using Pkg
Pkg.activate("/home/zhj/Rotating-Flow-ToolKit/compress/BEK/Vonkarmen_bone")
using RotatingDiskFlow
```

应当在项目目录下启动 Julia，以使用当前项目配置的依赖包：

```bash
cd /home/zhj/Rotating-Flow-ToolKit/compress/BEK/Vonkarmen_bone
julia --project=.
```

## 2. 在 Jupyter Notebook 中载入

新 notebook 的第一个单元格建议运行：

```julia
using Pkg
Pkg.activate(".")
using RotatingDiskFlow
```

已有 notebook 可以继续使用兼容入口：

```julia
include("NeutralCurveRunner.jl")
using .NeutralCurveRunner
```

同一个 Julia 会话中通常只需要执行一次 `include`。如果修改了源文件，建议重启
kernel 后重新载入，避免旧方法仍然留在 `Main` 中。

## 3. 计算单条 Lopez 中性曲线

```julia
config = CurveConfig(
    Tw = 1.10,
    omega = 0.0,
    R_initial = 500.0,
    beta_initial = 0.04,
    alpha_target = 0.1,
    num_modes = 2,
    model = :lopez,
    Mr = 0.3,
    Ro = -1.0,
    N_cheb = 69,
    beta_step = 8.0e-4,
    neutral_tol = 1.0e-7,
    output_dir = joinpath(pwd(), "neutral_curve_test"),
    keep_logs = false,
)

result = NeutralCurveRunner.run_case_with_retries(config)
```

`run_case_with_retries` 是当前批处理实际使用的稳健求解入口。该函数虽然没有被
`export`，仍可通过模块全名调用。它会依次执行：

1. 按给定初值搜索第一个中性点；
2. 用特征向量重叠追踪同一模态；
3. 在每个新 `beta` 上校正中性 Reynolds 数 `R`；
4. 求解困难时自动减小 `beta` 步长；
5. 普通设置失败后启用更稳健的双模态设置；
6. 必要时保存不连续的 Type-II 分支并从高 `beta` 搜索 Type-I 分支；
7. 成功后在 `keep_logs=false` 时删除单工况日志。

也可以直接调用导出的基础接口：

```julia
result = compute_neutral_curve(config)
```

该接口只执行一次延拓，不包含自动重试和 Type-I 回退，并且不会自动删除其日志。
研究新的分支结构时可以使用它，常规计算建议使用 `run_case_with_retries`。

## 4. 计算完全可压缩中性曲线

```julia
config = CurveConfig(
    Tw = 1.10,
    omega = 0.0,
    R_initial = 500.0,
    beta_initial = 0.04,
    alpha_target = 0.1,
    num_modes = 2,
    model = :compressible,
    Mr = 0.3,
    Ro = -1.0,
    N_cheb = 69,
    beta_step = 8.0e-4,
    neutral_tol = 1.0e-7,
    property_perturbations = true,
    base_property_variation = true,
    output_dir = joinpath(pwd(), "neutral_curve_test"),
    keep_logs = false,
)

result = NeutralCurveRunner.run_case_with_retries(config)
```

### 4.1 物性开关

| 工况 | `property_perturbations` | `base_property_variation` |
|---|---:|---:|
| 完整可压缩模型 | `true` | `true` |
| 关闭扰动物性反馈 | `false` | `true` |
| 关闭扰动物性反馈并冻结基本流物性 | `false` | `false` |

两个开关只影响 `model=:compressible`。Lopez 模型不读取这两个开关。

## 5. 主要参数

| 参数 | 默认值 | 说明 |
|---|---:|---|
| `Tw` | 必填 | 无量纲壁面温度 |
| `omega` | `0.0` | 给定无量纲频率 |
| `R_initial` | `500.0` | 初始中性点搜索所用 Reynolds 数 |
| `beta_initial` | `0.04` | 初始周向波数 |
| `alpha_target` | `0.1` | 初始特征值选择目标的实部 |
| `num_modes` | `1` | 请求追踪的模态数，只允许 `1` 或 `2` |
| `model` | `:lopez` | `:lopez` 或 `:compressible` |
| `Mr` | `0.3` | 当前可压缩基本流程序使用的马赫数参数 |
| `Ro` | `-1.0` | 旋转参数 |
| `N_cheb` | `69` | Chebyshev 配置点参数 |
| `beta_step` | `8e-4` | `beta` 延拓的名义步长，符号决定延拓方向 |
| `beta_bounds` | `(1e-3, 0.20)` | 允许的 `beta` 区间 |
| `min_beta_step` | `5e-5` | 自适应延拓允许的最小步长 |
| `neutral_tol` | `1e-7` | 中性条件残差容限 |
| `min_mode_overlap` | `0.60` | 连续模态的最小特征向量重叠 |
| `max_curve_points` | `500` | 单条曲线最大点数 |
| `output_dir` | `neutral_curve_batch` | 输出目录 |
| `keep_logs` | `false` | 稳健求解成功后是否保留日志 |

在 Type-I/Type-II 接近、耳状分支或模态相互作用区域，推荐设置：

```julia
num_modes = 2
min_mode_overlap = 0.5
```

正常精度计算建议先使用 `N_cheb=69`。网格收敛检查应分别使用多个 `N_cheb`
重复计算，而不是直接把生产计算的配置点数无限增大。

## 6. 返回结果

单工况函数返回一个命名元组：

```julia
result.data         # 中性曲线数值矩阵
result.path         # Tecplot DAT 文件路径
result.log_path     # 日志路径；日志可能已经被自动删除
result.stop_reason  # 延拓停止原因
result.validation   # 数据有效性检查结果
```

常见 `stop_reason` 包括：

| 状态 | 含义 |
|---|---|
| `:beta_limit` | 到达 `beta_bounds` |
| `:R_limit` | 曲线在足够多点后重新到达较大的 `R` |
| `:endpoint_no_root` | 减小步长后仍找不到下一中性点 |
| `:point_limit` | 到达 `max_curve_points` |

如果同时保存了不连续 Type-II 分支，稳健入口的返回结果还会包含：

```julia
result.typeII_path
```

## 7. 输出数据格式

每个 `.dat` 文件前两行是 Tecplot 变量和 Zone 信息，随后每行包含七列：

```text
omega  R  beta  alpha_r_1  alpha_i_1  alpha_r_2  alpha_i_2
```

中性条件为被追踪模态的空间增长率满足：

```math
\alpha_i = 0.
```

当只返回一个模态时，程序会将该模态同时写入第一组和最后一组特征值列。

典型文件名为：

```text
ome=0.0_Tw=1.1_model=lopez.dat
ome=0.0_Tw=1.1_model=compressible_Mr=0.3_propPert=on_baseProp=variable.dat
ome=0.0_Tw=1.1_model=compressible_Mr=0.3_propPert=off_baseProp=variable.dat
ome=0.0_Tw=1.1_model=compressible_Mr=0.3_propPert=off_baseProp=frozen.dat
```

## 8. 在 Notebook 中绘图

### 8.1 中性曲线

```julia
using Plots

data = result.data

plot(
    data[:, 2],
    data[:, 3],
    xlabel = "R",
    ylabel = "beta",
    label = "Lopez, Tw=1.10",
    linewidth = 2,
)
```

### 8.2 检查中性残差

```julia
neutral_residual = min.(abs.(data[:, 5]), abs.(data[:, 7]))

plot(
    data[:, 3],
    neutral_residual,
    yscale = :log10,
    xlabel = "beta",
    ylabel = "neutral residual",
    label = false,
)
```

### 8.3 读取已有输出文件

`load_curve` 是模块内部辅助函数，可通过模块全名调用：

```julia
path = joinpath(
    pwd(),
    "neutral_curve_test",
    "ome=0.0_Tw=1.1_model=lopez.dat",
)

data = NeutralCurveRunner.load_curve(path)
plot(data[:, 2], data[:, 3], xlabel="R", ylabel="beta", label="Lopez")
```

## 9. 标准批处理

当前 `standard_configs` 固定计算：

```math
T_w=1.06,1.08,\ldots,1.20.
```

每个壁温包含 Lopez 和三种可压缩物性组合，因此完整标准批次包含32条必需曲线。

### 9.1 串行批处理

```julia
results = run_standard_batch(
    output_dir = joinpath(pwd(), "neutral_curve_batch"),
    N_cheb = 69,
    keep_logs = false,
    resume = true,
    continue_on_error = true,
    worker = :all,
    prefer_typeI = false,
)
```

`worker` 可以选择：

```text
:all
:lopez
:compressible_full
:compressible_no_perturb
:compressible_frozen
```

`resume=true` 会读取状态文件、重新验证已有结果，并跳过已经完整完成的工况。
`prefer_typeI=true` 会直接采用高 `beta` 的 Type-I 初始化配置。

### 9.2 四进程并行批处理

```julia
result = run_parallel_standard_batch(
    output_dir = joinpath(pwd(), "neutral_curve_batch"),
    N_cheb = 69,
    keep_logs = false,
    resume = true,
)
```

该函数启动四个独立 Julia 进程，分别计算 Lopez 和三种可压缩物性工况。并行模式
会增加内存占用，使用前应确认机器能够同时容纳四套稳定性矩阵。

## 10. 直接从终端运行

默认命令会并行运行整个标准批次：

```bash
julia --project=. NeutralCurveRunner.jl
```

可以通过环境变量控制：

```bash
NEUTRAL_OUTPUT_DIR=neutral_curve_batch \
NEUTRAL_N_CHEB=69 \
NEUTRAL_KEEP_LOGS=false \
NEUTRAL_RESUME=true \
NEUTRAL_PARALLEL=true \
julia --project=. NeutralCurveRunner.jl
```

支持的环境变量如下：

| 环境变量 | 默认值 | 说明 |
|---|---|---|
| `NEUTRAL_OUTPUT_DIR` | `neutral_curve_batch` | 输出目录 |
| `NEUTRAL_N_CHEB` | `69` | Chebyshev 配置点参数 |
| `NEUTRAL_KEEP_LOGS` | `false` | 是否保留日志 |
| `NEUTRAL_RESUME` | `true` | 是否断点续算 |
| `NEUTRAL_PARALLEL` | `true` | 是否启动四进程并行批处理 |
| `NEUTRAL_WORKER` | `all` | 只运行指定模型组 |
| `NEUTRAL_PREFER_TYPEI` | `false` | 是否直接从 Type-I 配置开始 |
| `NEUTRAL_VALIDATE_ONLY` | `false` | 是否只验证已有批次 |

只运行 Lopez 模型的示例：

```bash
NEUTRAL_PARALLEL=false \
NEUTRAL_WORKER=lopez \
julia --project=. NeutralCurveRunner.jl
```

只验证已有文件：

```bash
NEUTRAL_VALIDATE_ONLY=true \
NEUTRAL_OUTPUT_DIR=neutral_curve_batch \
julia --project=. NeutralCurveRunner.jl
```

## 11. 验证输出文件

### 11.1 验证单条曲线

```julia
check = validate_curve_file(result.path, config)

check.ok
check.issues
size(check.data)
```

检查内容包括：

- 文件是否可以解析为七列数据；
- 点数是否达到最低要求；
- 是否包含 `NaN` 或 `Inf`；
- `R` 是否位于物理搜索区间；
- `beta` 是否严格单调；
- 中性残差是否满足容限。

### 11.2 验证标准批次

```julia
check = validate_standard_batch(
    output_dir = joinpath(pwd(), "neutral_curve_batch"),
    N_cheb = 69,
    cleanup_logs = true,
)

check.ok
check.invalid
check.summary_path
```

验证摘要会写入：

```text
neutral_curve_batch/batch_validation.tsv
```

只有全部必需曲线验证通过时，`cleanup_logs=true` 才会清理批处理日志。

## 12. 常见问题

### 12.1 初始中性点搜索失败

报错通常包含：

```text
No initial neutral crossing was found
```

首先检查 `R_initial`、`beta_initial` 和 `alpha_target` 是否靠近目标模态。Type-I
分支可以尝试：

```julia
config = CurveConfig(
    Tw = 1.10,
    model = :lopez,
    R_initial = 500.0,
    beta_initial = 0.12,
    alpha_target = 0.65,
    beta_step = -8e-4,
    num_modes = 2,
)
```

### 12.2 延拓在最后一个点无限寻找

当前稳健入口具有最小步长、最大扫描次数和 Type-I 回退控制。不要在外层自行无限
重复调用失败点。检查：

```julia
result.stop_reason
result.validation
```

并查看保留的日志：

```julia
config = CurveConfig(..., keep_logs=true)
```

### 12.3 曲线发生模态跳跃

建议使用：

```julia
num_modes = 2
min_mode_overlap = 0.5
beta_step = 4e-4
```

然后比较相邻点特征函数及日志中的 `overlap`。仅凭曲线形状不能判断是否发生模态
跳跃；耳状分支也可能是真实的物理解。

### 12.4 Notebook 中出现旧方法或关键词不匹配

如果修改过 `eigsol`、`CRD_STA.jl` 或 `NeutralCurveRunner.jl`，Notebook 可能仍保留
旧方法定义。最可靠的处理方式是：

1. 重启 Julia kernel；
2. 首先载入最新的 `NeutralCurveRunner.jl`；
3. 不再单独执行旧版 `eigsol` 或 `cur_legacy` 定义单元格。

## 13. 推荐工作流

正式计算建议按以下顺序进行：

1. 用 `N_cheb=69` 和单工况稳健入口试算；
2. 绘制中性曲线并检查残差；
3. 对关键点使用不同 `N_cheb` 进行收敛验证；
4. 确认初始模态和分支后运行批处理；
5. 使用 `validate_standard_batch` 生成最终验证摘要；
6. 再将通过验证的数据用于模型误差、临界点和物理机制分析。

## 14. Lopez Type-II 专用计算与分支合并

当通用入口在固定 `R=500` 下没有捕获低 `beta` 的 Type-II 模态时，可以使用
专用脚本。当前脚本为 `Tw=1.18` 和 `Tw=1.20` 提供了经过验证的低 `beta` 初值：

```bash
julia --project=. ComputeLopezTypeII.jl 1.18
julia --project=. ComputeLopezTypeII.jl 1.20
```

也可以在一次命令中依次计算：

```bash
julia --project=. ComputeLopezTypeII.jl 1.18 1.20
```

结果写入标准批处理目录，并使用独立的分支文件名：

```text
ome=0.0_Tw=1.18_model=lopez_branch=typeII.dat
ome=0.0_Tw=1.2_model=lopez_branch=typeII.dat
```

该脚本直接调用 `compute_neutral_curve`，不会在 Type-II 搜索失败后自动改追 Type I。
保存前还会检查点数、`beta` 范围、中性残差和活动特征值列是否发生切换。

将 Type-I 主曲线与 Type-II 分支合并为 Tecplot 多 Zone 文件：

```bash
julia --project=. MergeLopezNeutralBranches.jl
```

每个壁温会生成一个 `_allbranches.dat` 文件，其中 Type I 和 Type II 使用独立
Zone。这样可以在同一文件中绘图，同时不会在两个分支之间产生虚假的连接线。脚本
还会生成：

```text
neutral_curve_batch/lopez_allbranches_Tw1.08-1.2.dat
```

该汇总文件包含 `Tw=1.08:0.02:1.20` 的十四个 Zone，即每个温度各有一个
Type-I Zone 和一个 Type-II Zone。
