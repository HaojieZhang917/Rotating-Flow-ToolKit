# 退化零阶 inner core 的极限一致性检验

## 本步骤范围

本步骤只检验一致模型的既有 `epsilon=-Hinf` fold 数据是否选择零阶退化核心

```text
F0 = G0 = H0 = 0,    T0 = Tw,t.
```

不新增基本流、fold continuation 或参数扫描，也不推导一阶 inner--outer
matching。这里的后处理不是冻结基本流或项开关模型，而是对完整自洽一致模型
已有解的极限兼容性诊断。

## 理论判据

取

```text
s(Ro)     = (2 - Ro - Ro^2)/2,
omega(Ro) = s(Ro) + Ro,
chi(T)    = 1 - gamma*(T - 1),       gamma = 1.
```

零方位速度的 leading radial balance 要求

```text
Delta_core = chi(Tw,t)*s(Ro_t)^2 - omega(Ro_t)^2 = 0.
```

等价的解析相容曲线为

```text
Tw,t = Tcurve(Ro_t)
     = 1 + [1 - (omega(Ro_t)/s(Ro_t))^2]/gamma.
```

诊断将同时报告 `Tw,t-Tcurve(Ro_t)`、`Delta_core`、
`Delta_core/omega^2`，以及由 outer slow manifold 反推的入口 `g_m`。

## 已有数据依赖与检验设计

- fold 点及其高阶离散半范围：
  `work/results/consistent_epsilon_fold_chain/recommended_epsilon_fold_table.csv`；
- 五点、去掉 `epsilon=0.1` 后四点的线性/二次外推：
  `work/results/consistent_fold_tail_scaling/polynomial_fits.csv`；
- 既有 `N=320` 剖面在固定 `eta=1,2,5` 的样本：
  `work/results/consistent_inner_leading_matching/fixed_eta_profiles.csv`。

将进行三项检查：

1. 比较两个数据窗口、线性/二次模型得到的极限点与解析相容曲线；
2. 在每个 fold 数据点的既有离散半范围内取全部符号角点，重新拟合并传播到
   `Delta_core`。这些半范围是数值敏感性包络，不是概率置信区间；
3. 对固定 `eta` 的 `H,F,G,Tw-T` 做五点及四点二次外推，仅检查零截距是否与
   数据趋势相容，不把它解释成一阶 matching 系数。

完成后的数值结果、证据边界和结论将追加在本文，并在
`docs/RESEARCH_PROGRESS_LOG.md` 新增记录；既有 L014 不改写。

## 结果

### 极限点与解析曲线

`gamma=1` 时，不同窗口和拟合阶数给出：

| 数据窗口 | 阶数 | `Ro_t` | 数值 `Tw_t` | 理论 `Tcurve(Ro_t)` | `Tw_t-Tcurve` | `Delta_core/omega_t^2` | `g_m` |
|---|---:|---:|---:|---:|---:|---:|---:|
| 五点 | 1 | -0.319987569 | 1.493851392 | 1.493895328 | -4.394e-5 | 8.681e-5 | 1.504e-4 |
| 五点 | 2 | -0.320760492 | 1.494843888 | 1.494835182 | 8.706e-6 | -1.723e-5 | -2.979e-5 |
| 去掉 `epsilon=0.1` | 1 | -0.320488961 | 1.494495573 | 1.494505108 | -9.535e-6 | 1.886e-5 | 3.263e-5 |
| 去掉 `epsilon=0.1` | 2 | -0.320719730 | 1.494786622 | 1.494785639 | **9.827e-7** | **-1.945e-6** | **-3.363e-6** |

最高信息量的四小点二次外推满足

```text
Tw,t - Tcurve(Ro_t) = 9.82746e-7,
Delta_core           = -1.20851e-6,
Delta_core/omega_t^2 = -1.94521e-6.
```

即只用外推 `Ro_t`，解析关系预测
`Tw,t=1.494785638928`，与独立数值外推
`Tw,t=1.494786621674` 相差不足 `10^-6`。二次模型还使曲线失配相对相应线性
外推减少约一个数量级；这表明小失配不是由指定解析条件强制得到的。

### 既有离散误差的传播

对每个 fold 点在保存的 `N=240--320` 半范围内取全部符号角点，再分别做二次
拟合。四小点窗口得到保守包络

```text
Tw,t - Tcurve in [-2.21594e-4,  2.23572e-4],
Delta_core       in [-2.74924e-4,  2.72508e-4],
g_m              in [-7.65343e-4,  7.57994e-4].
```

零值均位于包络内部。中心 `Delta_core` 的绝对值只约为该包络半宽的
`0.44%`。五点窗口也包含零。这里的包络是把已保存离散半范围作最保守独立组合
得到的确定性敏感性范围，不是统计置信区间；它可能显著高估真实的联合误差。

### 有限 `epsilon` 与固定 `eta` 的趋零证据

直接在五个有限 `epsilon` fold 点计算同一核心缺陷，其相对值从
`-4.254e-3` 依次降至
`-8.520e-4,-1.204e-4,-2.898e-5,-6.851e-6`。因此解析相容关系不是只在
外推截距偶然满足，而是沿数据链系统趋近；最后两点的缺陷有效衰减阶约为二阶。

既有 `N=320` 剖面在固定 `eta=1,2,5` 上，最后两点
`epsilon=0.01,0.005` 的有效衰减指数为：

| 场 | `eta=1` | `eta=2` | `eta=5` |
|---|---:|---:|---:|
| `H` | 1.974 | 2.012 | 2.000 |
| `F` | 2.016 | 2.007 | 1.993 |
| `G` | 0.997 | 0.996 | 0.990 |
| `Tw-T` | 1.005 | 1.005 | 1.005 |

所以在每个检查过的固定 `eta=O(1)` 位置，`F,H` 约按 `epsilon^2` 消失，
`G,Tw-T` 约按 `epsilon` 消失。它们共同外推至
`F0=G0=H0=0,T0=Tw,t`；这比单独观察一个很小的 `g_m` 更强。

## 判定与边界

在正则 inner/outer 结构和当前数据精度下，退化零阶核心检验**通过**：

```text
chi(Tw,t)*s_t^2 = omega_t^2
```

与独立外推误差相容，固定 `eta` 剖面也选择同一个零状态。因此当前数值证据支持

```text
F0 = G0 = H0 = 0,    T0 = Tw,t,    g_m = 0,
```

而 O(1) 的热/方位速度调整被排入 `Y=epsilon*eta` 慢外层。

这里的“通过”是数值分支选择在现有分辨率下的确认，不是严格数学证明：有限五点
不能证明极限唯一，也不能排除非正则 overlap；解析相容曲线本身只约束
`(Ro_t,Tw,t)` 的关系，不能单独选出 `Ro_t`、证明 fold，或给出一阶
solvability。上述问题均留待后续单独步骤，本次不进入一阶 matching。

## 输出

- Julia 后处理：`work/scripts/validate_degenerate_inner_core.jl`；
- 各拟合窗口相容性：
  `work/results/degenerate_inner_core_validation/fit_window_compatibility.csv`；
- 误差传播：
  `work/results/degenerate_inner_core_validation/uncertainty_envelopes.csv`；
- 有限 `epsilon` 趋势：
  `work/results/degenerate_inner_core_validation/finite_epsilon_compatibility.csv`；
- 固定 `eta` 趋零检验：
  `work/results/degenerate_inner_core_validation/fixed_eta_zero_intercepts.csv` 和
  `fixed_eta_effective_decay.csv`。
