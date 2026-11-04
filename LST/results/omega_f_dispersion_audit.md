# 自由涡波频率 `omega_f=-beta_f` 核查

## 代码事实

当前 `LSTParameters` 将 `omega` 作为输入参数。调用脚本中直接写入：

```julia
const BETA = 0.07759
const OMEGA = -BETA
p = LSTParameters(..., beta=BETA, omega=OMEGA, ...)
```

代表性入口包括 `scan_roughness_wavenumber.jl`、`scan_incident_roughness_wavenumbers.jl`、`reconstruct_gaussian_roughness.jl`、`scan_second_order_modal_peak.jl` 和 `run_spatial_lst_receptivity.jl`。

`farfield_matrix` 和 `farfield_state_pencil` 使用给定的 `p.omega`，求解的是法向波数 `kz`。`solve_farfield_state_modes` 返回 `kz` 和极化向量，没有求解 `omega` 的分支。

因此当前工作流是：

```text
(alpha, beta, omega) 给定  ->  求 far-field kz 和极化
```

而不是：

```text
(alpha, beta, kz) 给定  ->  求 omega
```

## 物理来源

当前基流在圆盘旋转坐标系中的远场为

\[
U_\infty=0,\qquad V_\infty=-1,\qquad W_\infty=0.
\]

在 Fourier 约定

\[
q'\sim e^{i(\alpha r+\beta\theta+k_z z-\omega t)}
\]

下，无黏均匀远场涡量波的对流关系为

\[
\omega=\alpha U_\infty+\beta V_\infty+k_zW_\infty.
\]

代入当前远场速度得到

\[
\omega=-\beta.
\]

所以 `omega=-beta` 有一个明确的物理解释：它对应于在惯性参考系中静止、在圆盘旋转坐标系中被远场旋转流以方位速度 `V_inf=-1` 对流的扰动。

但这个关系没有在当前程序中由色散方程自动推导或求解；它是根据该物理假设手动设定的输入。

## 有黏远场的严格关系

当前程序的严格远场条件是

\[
\det M_\infty(\alpha,\beta,\omega,k_z)=0,
\]

其中 `M_inf` 由 `farfield_matrix` 构造。当前实现固定 `alpha,beta,omega`，通过状态空间广义特征值问题求 `kz`。

在有黏性时，涡量分支满足对流项与黏性项共同平衡，`omega=-beta` 不再是所有有限黏性 `kz` 的完整色散关系；它只是无黏远场对流匹配的物理选取。当前算例在这个输入下得到的低衰减涡量根约为 `kz=-0.1491im`，该根是给定 `omega` 后求出的结果。

## 结论

1. `omega_f=-beta_f` 在代码中是人为指定的输入。
2. 它不是由当前远场求解器反解得到的频率。
3. 它具有无黏旋转坐标系远场对流色散关系的物理依据。
4. 当前代码真正求解的是给定 `(alpha,beta,omega)` 下的 `kz` 和极化。
5. 若要让频率由远场色散关系得到，需要固定 `(alpha,beta,kz)`，重新构造关于 `omega` 的远场特征方程并求解 `omega`。
6. 改变 `alpha_f`、`beta_f` 或参考系时，不能无条件继续使用 `omega=-beta_f`；必须重新说明所采用的远场物理模型。
