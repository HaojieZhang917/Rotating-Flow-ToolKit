# 完整 fold 增广系统的一阶切向分析

## 本步骤范围

本步骤使用已保存的有限 `epsilon=-Hinf` fold 解，数值实现

```text
D_Z F_aug * Z_epsilon = -(0,0,0,1),
```

直接计算 fold locus 的

```text
dRo/d epsilon, dTw/d epsilon,
```

并检查固定-`Tw` 场 Jacobian 的左右零模和参数空间 solvability slope。这里不
新增 Ro 切片、不延拓到新 epsilon、不进入 outer 二阶或 adjoint Hessian 幅值
公式。读取的源数据是
`work/results/consistent_epsilon_fold_chain/convergence_N*/` 中已完成的完整自洽
fold 解。

## 完整增广系统

令 `U` 为离散的 `(H,F,G,T)`，`v` 为固定-`Tw` 场 Jacobian 的右零向量：

```text
F_aug(U,Tw,v,Ro;epsilon) = [
    R(U,Ro,Tw),
    J_U(U,Ro,Tw)*v,
    dot(v,v)-1,
    Hinf(U)+epsilon
] = 0.
```

在任一有限 `epsilon` fold 点，直接隐式微分：

```text
D_Z F_aug * Z_epsilon = [0,...,0,-1].
```

这个系统已经包含：基本流切向方程、fold 零向量方程的一阶导数、零向量规范，
以及 `(H_epsilon)_infinity=-1`。代码复用 `BEKConsistent.jl` 中解析 Jacobian、
方向 Hessian 和 `Ro` 导数；不对整个增广残差做低精度外层差分。

## 检验量

对每个既有小 `epsilon` 点记录：

- `dRo/d epsilon`、`dTw/d epsilon` 及其比值；
- 增广线性系统残差和 `dHinf/d epsilon+1`；
- 左零模投影
  `-<w,R_Ro>/<w,R_Tw>` 与 `(dTw/d epsilon)/(dRo/d epsilon)` 的一致性；
- 与既有四小点二次拟合导数
  `A+2B epsilon` 的差；
- 右零模各场分量的幅值以及在 `Y=epsilon*eta` 坐标下的保存剖面。

计算完成后将在本文追加结果，并在总体研究日志新增记录。

## 数值结果

### 完整增广切向系统可解且通过恒等检查

在每个已保存 fold 点，增广线性系统直接给出
`dHinf/d epsilon=-1`。同时，切向参数比

```text
(dTw/d epsilon)/(dRo/d epsilon)
```

与左零模投影

```text
-<w,R_Ro>/<w,R_Tw>
```

在所有九组 `N=240,280,320`、`epsilon=0.02,0.01,0.005` 计算中一致到
`O(10^-11)` 或更好。这验证了实现确实同时满足基本流切向、fold 条件导数、
零模规范和 `epsilon=-Hinf` 归一化，而不是对 branch table 做数值差分。

### 参数切向及分辨率

| `epsilon` | `N` | `dRo/d epsilon` | `dTw/d epsilon` | 比值 |
|---:|---:|---:|---:|---:|
| 0.02 | 240 | -0.43919661 | 0.53440714 | -1.21678338 |
| 0.02 | 280 | -0.43912223 | 0.53431672 | -1.21678357 |
| 0.02 | 320 | -0.43914114 | 0.53433971 | -1.21678354 |
| 0.01 | 240 | -0.42696196 | 0.51918248 | -1.21599234 |
| 0.01 | 280 | -0.42717724 | 0.51945113 | -1.21600843 |
| 0.01 | 320 | -0.42991433 | 0.52277931 | -1.21600810 |
| 0.005 | 240 | 0.02665453 | -0.03240236 | -1.21564145 |
| 0.005 | 280 | -0.54836632 | 0.66673254 | -1.21585246 |
| 0.005 | 320 | -0.46176837 | 0.56136318 | -1.21568133 |

`epsilon=0.02` 的绝对切向已随 `N` 收敛；`epsilon=0.01` 有约 `0.7%` 的
分辨率散布；`epsilon=0.005` 的绝对切向尚未收敛，不能用于极限幅值判断。
最后一点的 fold 二次非退化系数已经降到 `O(10^-11)`，其符号随 `N` 改变，
正是 L010 已识别的分辨率诱导 cusp-like 层。增广切向需要求逆，因此比 fold
坐标更早放大这类误差。尝试把既有 fold 从 `O(10^-9)` 残差继续校正到
`2e-11` 时停在离散/差分噪声底，不能据此制造虚假的高精度切向。

仅使用目前仍可解析的 `epsilon=0.02,0.01`，假设
`dRo/d epsilon=A_R+O(epsilon)`、`dTw/d epsilon=A_T+O(epsilon)`，`N=320`
两点线性截距为

```text
A_R = -0.4206875201,
A_T =  0.5112188991.
```

与既有 branch-coordinate 二次拟合
`(-0.4203424951,0.5106601141)` 分别相差 `0.082%` 和 `0.109%`。这是对目标幅值
很强的支持，但还不是空间分辨率完全收敛的第三个独立预测：相同截距在
`N=240,280` 为 `A_R=-0.41473,-0.41523`，说明导数极限仍比坐标极限敏感。

### Full-system fold 零模的 matched scaling

把右零模按 `max|v_G|=1` 归一化后，九组数据给出非常稳定的组合：

```text
max|v_H| / (epsilon max|v_G|)       = 3.063 -- 3.076,
max|v_F| / (epsilon^2 max|v_G|)     = 2.798 -- 2.823,
max|v_T| / max|v_G|                 = 0.567 -- 0.579,
|v_H(infinity)|/(epsilon max|v_G|)  = 2.756 -- 2.759.
```

因此原始 fixed-`Tw` fold eigenvector 的清晰极限缩放是

```text
v_G,v_T = O(1),    v_H = O(epsilon),    v_F = O(epsilon^2),
```

与 thermo-rotational outer 基本流的 distinguished scaling 完全一致。原始
欧氏归一化零模不会直接收敛为 leading outer shooting Jacobian 的零向量；经过
上述分量重标度后，它才具有稳定的 matched outer 结构。这正好解释了“有限
epsilon fold 奇异、归一化 leading outer BVP 正则”为什么不矛盾。

### 切线方向继续趋向解析相容曲线

尽管最小点的绝对速度受分辨率污染，参数比仍稳定，并从
`-1.2167835`、`-1.216008` 向解析极限

```text
Tcurve'(Ro_t) = -1.21543131
```

靠近。这是因为共同的近零非退化系数会放大 `dRo/d epsilon` 和
`dTw/d epsilon`，但在比值中大幅抵消。左零模投影与该比值的逐点一致又提供了
独立的 full-system solvability 检查。

## 当前判定

附件提出的两级机制得到部分完成且结构上被确认：

```text
leading normalized outer BVP  ->  selects (Ro_t,Tw,t),
full augmented fold tangent   ->  fixes the epsilon-parametrized departure speed.
```

完整一阶增广方程确实直接给出了有限-epsilon fold 的绝对切向，无需先推基本流
`O(epsilon^2)`。在可解析的 `epsilon=0.02,0.01` 上，它支持
`A_R≈-0.42034,A_T≈0.51066`；但由于 `epsilon=0.005` 已进入分辨率诱导退化层，
当前不能宣称幅值已经像极限位置和方向那样完全收敛。

最坚实的新结论是：full fold 零模的 matched scaling 已得到高度一致的三阶
分量结构，且增广 tangent/左零模 solvability/解析相容曲线三者相互闭合。

## 输出

- 共享完整切向实现：`work/src/BEKConsistent.jl` 中 `fold_hinf_tangent`；
- Julia 计算：`work/scripts/analyze_full_fold_tangent.jl`；
- Julia 汇总：`work/scripts/summarize_full_fold_tangent.jl`；
- 分辨率表：
  `work/results/full_fold_tangent_analysis/resolution_comparison.csv`；
- 零模 outer 缩放：
  `work/results/full_fold_tangent_analysis/right_mode_outer_scaling.csv`；
- 可解析两点极限估计：
  `work/results/full_fold_tangent_analysis/resolved_tangent_limit_estimates.csv`；
- 各 `N` 的 tangent 表与 mode 剖面：
  `work/results/full_fold_tangent_analysis/N240/`、`N280/`、`N320/`。
