# RotatingDiskFlow

旋转圆盘边界层基本流、线性稳定性和中性曲线延拓工具。当前项目包含：

- Lopez 广义 Boussinesq 基本流与空间/时间稳定性算子；
- 完全可压缩旋转圆盘稳定性算子及物性反馈开关；
- 基于特征向量重叠的中性曲线延拓、失败恢复和批处理；
- Sutherland 物性律下的径向 BDF2 推进基本流求解器；
- Malik (1986) 中性点回归和 Lopez 算子有限差分线性化测试。

Git 仓库根目录是 `/home/zhj/Rotating-Flow-ToolKit`。本目录是其中可独立激活的 Julia 子项目。

## 快速开始

```bash
cd /home/zhj/Rotating-Flow-ToolKit/compress/BEK/Vonkarmen_bone
julia --project=.
```

在 Julia 或 Jupyter 中：

```julia
using RotatingDiskFlow

config = CurveConfig(
    Tw=1.10,
    model=:lopez,
    R_initial=500.0,
    beta_initial=0.04,
    alpha_target=0.10,
    num_modes=2,
    N_cheb=69,
)

result = RotatingDiskFlow.NeutralCurveRunner.run_case_with_retries(config)
```

旧 notebook 仍可使用：

```julia
include("NeutralCurveRunner.jl")
using .NeutralCurveRunner
```

根目录同名 `.jl` 文件只是兼容入口，新程序应优先 `using RotatingDiskFlow`。

## 基本流

Lopez 基本流：

```julia
profile = RotatingDiskFlow.LopezBaseflow.solve_lopez_baseflow(
    1.10;
    zmax=40.0,
    tol=1e-8,
    nout=2001,
)

profile.z
profile.F
profile.G
profile.H
profile.T
profile.metrics
```

Sutherland 径向推进模型中，目标位置 `R` 的局部马赫数满足 `Mr=Minf*R`：

```julia
Suth = RotatingDiskFlow.SutherlandMarching
R = 20.0
Mr = 0.3

params = Suth.Params(
    Tw=1.5,
    Minf=Mr/R,
    rmax=R,
    Nr=41,
    Nz=81,
    verbose=false,
)

r, z, stations = Suth.solve_marching(params)
profile_at_R = stations[end]
```

`solve_marching` 只返回内存数据；需要文件时再调用 `Suth.write_outputs(r,z,stations,params)`。

## 命令行入口

```bash
julia --project=. scripts/neutral_curves.jl
julia --project=. scripts/lopez_type_i.jl 1.08
julia --project=. scripts/lopez_type_ii.jl 1.18 1.20
julia --project=. scripts/merge_neutral_branches.jl
julia --project=. scripts/grid_independence.jl
julia --project=. scripts/sutherland_marching.jl --Tw=1.5 --Mr=0.3 --R=20
```

批处理参数也可通过 `NEUTRAL_OUTPUT_DIR`、`NEUTRAL_N_CHEB`、`NEUTRAL_PARALLEL`、`NEUTRAL_RESUME` 等环境变量设置，详见[中性曲线使用说明](docs/neutral-curves.md)。

## 测试

快速测试：

```bash
julia --project=. test/runtests.jl
```

包含 Malik 两个中性点的完整物理回归：

```bash
RUN_PHYSICS_REGRESSION=true julia --project=. test/runtests.jl
```

当前验证结果为快速测试 `48/48`、完整测试 `52/52`。测试内容和容差见[测试说明](docs/testing.md)。

## 项目结构

```text
src/                     核心 Julia 包
scripts/                 可执行入口
examples/                Notebook/绘图示例
test/                    单元、算子和物理回归测试
docs/                    使用、推导、矩阵审计和网格验证
neutral_curve_batch/     运行时生成，不纳入 Git
neutral_curve_integrated/运行时生成，不纳入 Git
```

详细模块关系见[代码结构](docs/architecture.md)。

## 文档

- [中性曲线使用说明](docs/neutral-curves.md)
- [代码结构与公共 API](docs/architecture.md)
- [测试与回归基准](docs/testing.md)
- [Lopez 稳定性方程推导](docs/lopez-stability-derivation.md)
- [可压缩稳定性矩阵审计](docs/crd-stability-audit.md)
- [双模态网格无关性验证](docs/grid-independence.md)

## Python 依赖

Lopez 基本流通过 `PyCall` 调用 SciPy。Python 依赖列在 `requirements.txt`：

```bash
python3 -m pip install -r requirements.txt
```

如 Julia 的 `PyCall` 未使用该 Python：

```julia
ENV["PYTHON"] = "/usr/bin/python3"
using Pkg
Pkg.build("PyCall")
```
