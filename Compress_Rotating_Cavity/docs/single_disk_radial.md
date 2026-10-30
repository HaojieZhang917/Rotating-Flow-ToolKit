# 旋转单盘：非相似稳态边界层的径向吸引性

这个算例检验一个具体问题：在稳态、轴对称、不可压、等温的单盘边界层模型中，给定偏离 von Kármán 主支的入口速度剖面，沿半径向外发展后，归一化剖面是否趋近该主支。外部流体静止、壁面无抽吸，外加径向压力梯度为零。

## 方程和边界条件

定义

\[
s=\log(r/r_0),\qquad \eta=z\sqrt{\Omega/\nu},
\]

\[
u_r=\Omega rF(s,\eta),\quad
u_\theta=\Omega rG(s,\eta),\quad
u_z=\sqrt{\nu\Omega}H(s,\eta).
\]

保留非相似径向导数，边界层方程为

\[
H_\eta+2F+F_s=0,
\]

\[
F_{\eta\eta}-HF_\eta-F^2+G^2-FF_s=0,
\]

\[
G_{\eta\eta}-HG_\eta-2FG-FG_s=0.
\]

壁面满足 \(H=F=0,G=1\)，远场满足 \(F=G=0\)。截断远场不额外指定 \(H\)，其值由卷吸和连续性确定。

在这套无量纲变量下，经典相似性主支是 \(F_s=G_s=0\) 的固定点。相似性解存在本身不保证固定点具有吸引性；径向谱和非线性入口实验用于检验后一个性质。

入口只指定 \(F(0,\eta),G(0,\eta)\)。非相似入口不能先用 \(H_\eta=-2F\) 生成轴向速度并作为额外入口条件，因为这样会漏掉 \(F_s\)。代码中该积分只用于第一站 Newton 初值；入口输出的 \(H\) 是 `NaN`，下游则联立求解完整连续性方程。

## 离散方法

- 法向使用 Chebyshev–Lobatto 配点；`degree=N` 表示 `N+1` 个节点。
- 径向使用隐式 BDF2，第一站用隐式 Euler 建立历史。步长略作调整，使末站恰好达到指定半径比。
- 每站通过解析 Jacobian 和回溯 Newton 法联立求解 `H,F,G`。
- 不把剖面投影回相似性解，不裁剪负速度，也不使用 FLARE 修改方程。

壁面和外缘的零径向速度不作为回流。若入口或已收敛站点的内部节点满足 `F < -reverse_tolerance`，返回 `:reverse_flow` 并停止推进；Newton 失败则返回 `:nonlinear_failure`。出现回流意味着单向向外推进失去适用性，不能把该结果解释成已求得回流区的正确解。

## 五组入口

令 \(F_*,G_*\) 为同一网格上的相似性参考解。入口设置为

\[
b(\eta)=\left(\frac{\eta}{c}\right)^p
\exp\left[p\left(1-\frac{\eta}{c}\right)\right],\quad
F_{\rm in}=F_*(1+a_Fb),\quad G_{\rm in}=G_*(1+a_Gb).
\]

| 算例 | \(a_F\) | \(a_G\) | \(c\) | \(p\) |
|---|---:|---:|---:|---:|
| `similarity` | 0 | 0 | 2.0 | 4 |
| `small_shape` | 0.02 | −0.02 | 1.5 | 4 |
| `small_opposite` | −0.02 | 0.02 | 3.0 | 6 |
| `large_radial` | 1.50 | −0.70 | 1.8 | 4 |
| `large_swirl` | −0.65 | 0.80 | 2.5 | 6 |

扰动在壁面高阶消失，保留低阶壁面相容条件；不同中心和幂次改变扰动形状。五个入口均保持径向外流。它们覆盖若干有限幅值变化，尚未穷尽所有可接受入口。

采用 Clenshaw–Curtis 积分权重计算相对误差

\[
E_{FG}(s)=\left[
\frac{\int[(F-F_*)^2+(G-G_*)^2]\,d\eta}
{\int(F_*^2+G_*^2)\,d\eta}
\right]^{1/2}.
\]

还记录各分量误差、壁面梯度、外缘卷吸、非线性残差和积分连续性误差。轴向速度不计入入口误差，因为它不是独立指定的入口数据。

## 运行与 API

在项目目录运行：

```bash
julia --project=. scripts/benchmark_single_disk_radial.jl
julia --project=. scripts/verify_single_disk_radial.jl
julia --project=. scripts/plot_single_disk_radial.jl
julia --project=. scripts/spectrum_single_disk_radial.jl
```

默认设置为 `N=80, eta_max=20, ds≈0.02, r/r0=100`。可通过命令行更改：

```bash
julia --project=. scripts/benchmark_single_disk_radial.jl --degree 80 --eta-max 24 --ds 0.01 --ratio 10 --output benchmark_results/single_disk_radial_eta24
julia --project=. scripts/plot_single_disk_radial.jl benchmark_results/single_disk_radial_eta24
```

`--tolerance 1e-11` 可调整 Newton 残差阈值，以检查末段平台的数值敏感性。阈值过低时，双精度配点残差可能无法达到要求，求解器会返回失败状态。

核心求解器只依赖 Julia 标准库 `LinearAlgebra`。绘图脚本用 Julia `PyCall` 调用已有的 matplotlib。可直接加载核心模块，也可通过项目包导入：

```julia
using CompressRotatingCavity.SingleDiskRadial
using LinearAlgebra
BLAS.set_num_threads(1)

cfg = RadialConfig(degree=80, eta_max=20.0, ds=0.02, s_end=log(10.0))
ref = reference_profile(cfg)
inlet = inlet_profile(ref; case=:large_radial)
run = march(ref, inlet)
@show run.status run.message run.metrics[end]

# 自定义入口：保持壁面/远场条件，并自行检查相容性。
custom = merge(inlet, (; F=copy(inlet.F), G=copy(inlet.G)))
# 修改 custom.F、custom.G 后调用 march(ref, custom)。
```

`run.F, run.G, run.H` 的行为法向节点、列为径向站点。`run.s` 是径向对数坐标；`exp.(run.s)` 才是半径比。`RadialConfig(method=:be)` 可改用一阶隐式 Euler。

## 输出与本次验证

输出目录 `benchmark_results/single_disk_radial/` 中：

- `history.csv`：全部径向站点的误差与诊断量。
- `profiles.csv`：若干半径处的完整法向剖面；`inlets.csv`：入口参数。
- `summary.csv, reference.txt, results.md`：参数和实际计算结果。
- `verification.csv/.txt`：独立的数值验证结果。
- `tolerance_check.txt`：小扰动末段对 Newton 停止阈值的敏感性检查。
- `error_decay.png/.pdf, profiles.png/.pdf`：误差衰减和入口/最终剖面对比。
- `spectrum.csv, spectrum.txt`：局部径向谱与分辨率检查。

验证脚本已通过 24 项检查，包括有限域参考壁面梯度与卷吸量、非零径向历史项下的解析 Jacobian、相似性入口保持、有限扰动推进、负径向入口停止，以及径向步长、法向节点数、外边界位置三项独立加密检查。

在 \(r/r_0=4\) 处，大径向扰动的 `ds≈0.04,0.02,0.01` 序列给出约 `2.06` 的观测收敛阶。`eta_max=20→24` 时，小扰动最终剖面在共同内区的相对差约 `2.1e-7`；参考解的外缘 \(H\) 差约 `6.7e-7`。因此相对同域离散参考解的极小误差，不代表无限域物理解具有同样的绝对精度。

本次五组入口全部达到 \(r/r_0=100\)。两组大扰动的 \(E_{FG}\) 分别从 `0.3903`、`0.2161` 降至约 `5.72e-9`、`9.95e-10`；相似性控制组的最大漂移约 `2.4e-16`。小扰动末段约 `1e-10` 的平台受 Newton 停止阈值影响，不能用于推断物理衰减饱和。

单独将 `small_shape` 的阈值从 `2e-10` 收紧到 `1e-11`，最终误差由 `6.55e-11` 变为 `5.21e-11`，末端恢复持续衰减。默认结果的最后两站均未进行 Newton 更新；严格阈值下最后一站进行了更新。该检查记录在 `tolerance_check.txt`。

`summary.csv` 中的斜率在 \(0.8\le s\le1.8\) 上拟合 \(\log E_{FG}\)，只是有限区间衰减率，不是特征值。局部谱脚本独立线性化上述完整非相似方程，采用 \(\delta q\propto e^{\kappa s}=(r/r_0)^\kappa\) 的约定，检验小扰动的径向增长或衰减。

在所检查的 \(|\kappa|\le50\) 窗口中，六个通过网格、计算域和模态形状匹配检查的候选模态均为实负值，最弱衰减模态为 \(\kappa\approx-3.8136\)。其 `N=48,64,80` 结果差约 `1e-11`，`eta_max=20→24` 的差约 `3.8e-6`。用该模态的 `1e-6` 幅值作为非线性推进入口，在短距离内拟合得到 `−3.813576`，与谱值及其符号一致。这是对所检查空间模态的数值证据；奇异质量矩阵和远场径向速度退化使高阶谱敏感，有限窗口计算不构成完整函数空间的谱证明。

这些结果支持所测入口在本模型内趋向 von Kármán 主支。它们不证明所有入口的全局吸引性，也不能直接推出真实有限圆盘流动的时间稳定性或三维稳定性。尚未指定 \(r_0\sqrt{\Omega/\nu}\)，因此默认半径比不是某一实际装置的雷诺数范围。
