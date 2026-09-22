# 可压缩 BEK 空间稳定算子的 `Ro/Co` 推导与代码核对

日期：2026-09-12

## 1. 结论先行

保留的旧一般 `Ro` 可压缩空间稳定算子不能继续用于 Ekman/Bödewadt
中性根计算。写错的不是单独一个矩阵元，而是同一套错误同时存在于

- `RotatingDiskFlow/src/CRD_STA.jl::Spatial_mode_BEK1`（展开 QEP）；
- `RotatingDiskFlow/src/CRD_STA.jl::Spatial_mode_BEK`（张量 QEP）；
- `work/src/CompressibleBEKLinearStability.jl::assemble_qep` 对上述算子的
  `Ro=0,1` 仿射外推包装；
- `work/src/CompressibleBEKLinearStability.jl::source_profiles` 的一般 `Ro`
  热方程。

两套旧 QEP 在机器精度内一致，只说明它们是同一个公式的两种装配，不能作为
独立物理验证。严格的 `Ro=-1` von Kármán 切片可以在旧变量约定下复现
Turkyilmazoglu--Uygun (2006)，但把其中若干项直接替换为 `Ro` 和 `Co`
并不能得到一般 BEK 算子。

本审计之后已在 `work/src/CompressibleBEKLinearStability.jl` 实现统一的
Lingwood 直接变量候选算子；保留的 `RotatingDiskFlow/src/CRD_STA.jl` 未作
修改。实现和数值结果见第 8 节。该实现通过三个低马赫端点 Gate，但仍标记
为 `derived_bek_candidate`，不能把它说成两份 von Kármán 文献直接证明的
任意 `Ro` 可压缩模型。

## 2. 文献中已经确定的部分

### 2.1 可压缩 von Kármán 母方程

Turkyilmazoglu--Uygun (2006) 的式 (1)--(9) 在随盘旋转的坐标系中给出
可压缩基流；式 (10)--(15) 给出五个扰动方程。论文中的基流量是

\[
 (u_B,v_B,w_B)=(r\bar u_B,r\bar v_B,Re^{-1/2}\bar w_B),
\]

并以 Dorodnitsyn--Howarth 相似坐标离散。动量方程的判别性结构是

\[
 \mathcal U_p=\alpha\bar u_B+\beta\bar v_B-\omega,
 \qquad
 -2(\bar v_B+1)\hat v,
 \qquad
 +2(\bar v_B+1)\hat u.
\]

William 的学位论文第 2 章从旋转坐标系的原始连续性、动量、状态和能量
方程重新得到同一 disk limit：径向原方程中同时有曲率项
`-v^2/x`、Coriolis 项 `-2v` 和离心势项 `-x`；切向原方程有
`+2u`。第 2.30--2.32 式恢复相同的广义 von Kármán 基流方程。
该论文第 3.9--3.14 式又独立给出了线性化形式，可用于检查第 2 章推导的
扰动结果。

这两份文献只直接证明 `Ro=-1` 的可压缩旋转盘切片；它们没有证明当前代码
中 `Ro=0,1` 的插入方式。

### 2.2 BEK 的统一尺度

一般 BEK 扩展必须采用 Lingwood 约定

\[
 \Delta\Omega=\Omega_F-\Omega_D,
 \qquad Ro={\Delta\Omega\over\Omega},
 \qquad Co={2\Omega_D\over\Omega}=2-Ro-Ro^2.
\]

基流速度以**有符号的** `DeltaOmega` 缩放，而边界层长度以
`sqrt(nu/Omega)` 缩放。令 Lingwood 的直接基流变量为 `(U,V,W)`，则

\[
2U+W'=0,
\]

\[
Ro\{U^2+WU'-(V^2-1)\}-Co(V-1)-U''=0,
\]

\[
Ro\{2UV+WV'\}+CoU-V''=0.
\]

旧后端 `sol_baseflowODE` 内部求得 Lingwood 的 `(U,V,W)`，但返回
`(-U,-V,-W)`，随后 `baseflow_var` 又对 `Ro>=0` 做一次翻号。这就是旧代码
在 `Ro=-1` 使用论文变量、在 `Ro>=0` 使用直接变量的端点依赖约定。
新前端只在读取后端返回值时统一撤销第一次翻号，从而对所有 `Ro` 都使用
Lingwood 的直接 `(U,V,W)`；稳定算子内部不再根据端点翻剖面。

## 3. `Ro` 和 `Co` 应放在哪里

下面先用 Lingwood 的直接变量 `(U,V,W)` 写结构，以避免代码中历史翻号
造成歧义。令

\[
 \mathcal U=\bar\alpha U+\bar\beta V-\bar\omega,
 \qquad \Gamma=2RoV+Co.
\]

Newtonian 局部平行流极限的无粘骨架为

\[
(iR\bar\alpha+Ro)\hat u+iR\bar\beta\hat v+R\hat w'=0,
\]

\[
\{iR\mathcal U+RoU+RoW D\}\hat u
-\Gamma\hat v+RU'\hat w+iR\bar\alpha\hat p-\mathcal V_r=0,
\]

\[
\{iR\mathcal U+RoU+RoW D\}\hat v
+\Gamma\hat u+RV'\hat w+iR\bar\beta\hat p-\mathcal V_\theta=0,
\]

\[
\{iR\mathcal U+RoW'+RoW D\}\hat w
+R\hat p'-\mathcal V_z=0.
\]

由此得到不能混淆的放置规则：

| 项 | `Ro` | `Co` | 说明 |
|---|---:|---:|---|
| Doppler `alpha*U+beta*V-omega` | 否 | 否 | 局部波数、频率已随 BEK 局部尺度重标度 |
| 柱坐标连续性的 `u/r` | 是 | 否 | 归一化后给 `(Ro/R) u` |
| 动量中的 `U*u`、`W*D`、`W'` | 是 | 否 | 曲率/法向基流惯性 |
| `U' w`、`V' w` | 否 | 否 | 乘以 `R` 后不再有显式 `Ro` |
| 径向--切向速度耦合 | 通过 `2Ro*V` | 通过 `Co` | 合并为 `Gamma=2Ro*V+Co` |
| `rho_hat` 乘基流加速度 | 是 | 是 | `Co` 可在密度列出现，不能只检查速度交叉块 |
| 轴向动量、状态方程、热传导 | 按相应惯性/传导尺度 | 否 | `Co` 不直接进入 |
| 能量中的基流法向平流/压力功 | 是 | 否 | 与 `Ro*W*D` 同源 |
| 黏性耗散 | 取决于 Mach 定义 | 否 | 若 Mach 以 `Omega` 定义则出现 `Ro^2`; 若以 `DeltaOmega` 定义则吸收到 Mach 中 |

特别地，“Doppler 前整体乘 `Ro`”也是错误的。一般 BEK 的局部稳定尺度把
`Ro` 留在柱坐标曲率和法向相似平流中，而不留在
`alpha*U+beta*V-omega` 前。这一点与全局原始方程中看见的 `Ro` 不能直接
混用。

## 4. 对当前代码的逐块判定

### 4.1 明确错误：`Ro=-1 -> +1` 特殊分支

`Spatial_mode_BEK1` 第 396--398 行和 `Spatial_mode_BEK` 第 527--529 行
把 `Ro=-1` 改成 `+1`，但没有同时执行由尺度变化要求的完整
`(U,V,W,u,v,w,alpha,beta,omega)` 变换。这个分支是不连续的，也无法作为
一般 BEK 母算子的一部分保留。

### 4.2 明确错误：连续性曲率项没有 `Ro`

`D_11=rho`（第 626 行）把柱坐标 `u/r` 固定成了 von Kármán 切片的单位
系数。一般 BEK 局部方程要求该项随 `Ro`（符号还要与最终选择的扰动变量
变换一致）。同理，密度列中的基流法向输运 `C_14=rho*H`（第 587 行）
也不能在一般 `Ro` 下保持无系数；必须从压缩连续性整体重推。

这项在 Ekman 极限最具判别性：`Ro=0` 时局部柱坐标曲率贡献应按 BEK
尺度退化，而当前代码仍保留单位贡献。

### 4.3 部分结构正确，但输入剖面符号错误

当前动量块的以下组合与 BEK 骨架相符：

- `D_22=-rho*(2Ro*G+Co)`、`D_31=+rho*(2Ro*G+Co)`；
- `Ro*rho*F` 对角项；
- `Ro*rho^2*H*D` 法向平流项；
- `R*rho^2*DF`、`R*rho^2*DG` 基流剪切项不显式乘 `Ro`；
- `D_34=F*(2Ro*G+Co)+Ro*rho*H*DG` 具有“密度扰动乘基流切向加速度”的正确类型。

旧 `baseflow_var` 的端点依赖翻号使同一批矩阵元在不同 `Ro` 下代表不同
物理变量。新实现将剖面统一为 `(U,V,W)`；因此 `alpha*U`、`Ro*U`、
`Ro*W*D`、`U'`、`W'` 使用同一直接变量约定，`Gamma=2RoV+Co` 不再与
剖面符号冲突。

### 4.4 明确错误：能量行没有继承一般 BEK 的法向 `Ro`

第 617--621、643--646 行的 `H` 平流、压力功和热力耦合仍按
von Kármán 公式保留单位系数，而动量行第 590、596、602 行已经使用
`Ro*H`。同一个父方程中不可能出现这种尺度分裂。能量行至少需要对
`Ro*W*D`、相关压力功和柱坐标散度项整体重推；不能通过给现有几个
`H` 矩阵元机械乘 `Ro` 修补。

`Co` 不应直接进入能量方程，因为 Coriolis 力不做功；它只可通过基流
剖面间接影响能量方程。

### 4.5 明确错误：一般 `Ro` 热 BVP 的符号

当前 `source_profiles` 使用

```text
q'' =  Pr*Ro*H*q'
f'' =  Pr*Ro*H*f' + 2Pr*Ro*F*f + 2Pr*(F'^2+G'^2)
```

其中当前 `H=-W`。由 BEK 能量方程和 `Ro=-1` 的论文式 (9) 同时得到的
一致形式应为

\[
q''+Pr\,Ro\,Hq'=0,
\]

且 `f` 方程的两个线性输运项也应取相反号，才能在 `Ro=-1` 恢复
Turkyilmazoglu--Uygun 的

\[
f''+2Pr\,\psi f'-2Pr\,\psi'f
=2Pr\{(\psi'')^2+(\bar v_B')^2\}.
\]

黏性耗散项是否再显式带 `Ro^2`，取决于 `Mr` 是用 `Omega` 还是
`DeltaOmega` 定义。当前代码一方面令 `Ma=Mr/R`，另一方面又在温度重构
中使用 `(Ro*Mr)^2`，必须先写清 Mach 定义后再定系数。

### 4.6 仿射外推不是推导

`assemble_qep` 在 `Ro=0` 和 `Ro=1` 计算同一旧算子，再用

```text
C(Ro) = C(0) + Ro*(C(1)-C(0))
```

绕过 `Ro=-1 -> +1` 分支。这只是在代数上回收旧代码的显式线性项，既
没有修复连续性和能量行缺项，也没有修复 `(F,G,H)` 的物理含义。当前
`Ro=0,1` 根即使后向误差很小，也只是错误 QEP 的精确特征根。

## 5. 应采用的重构顺序

1. 固定直接 BEK 变量 `(U,V,W)`、扰动变量、正 `R` 约定，以及
   `(alpha,beta,omega)` 的尺度；禁止在底层函数内翻 `Ro`。
2. 从旋转系的可压缩守恒形式重新无量纲化，先写出含 `Ro/Co` 的连续性、
   三动量、总能量和状态方程。
3. 在线性化之前代入相似基流；保留 `rho_hat` 乘全部基流加速度，因此
   `Co` 同时进入交叉速度块和相应密度列。
4. 再作 Dorodnitsyn--Howarth 变换及压力消元；最后自动生成
   `(A0,A1,A2)`，不要手工维护两套公式。
5. 设置三个强制 Gate：`Ro=-1` 逐项恢复论文式 (10)--(15)；
   `Ro=0` 恢复 Ekman 的线性常系数远场；`Ro=1, Co=0` 不残留任何直接
   Coriolis 项。
6. 之后才做 `Mr -> 0`、网格收敛和基于特征向量重叠的连续模态跟踪。

## 6. 初始审计的证据范围

初始审计本身没有运行新的中性曲线计算，足以否定旧一般 `Ro` QEP，但不
等于验证替代模型。后续候选实现和计算作为新的证据层单独记录如下。

## 7. 来源

- M. Turkyilmazoglu and N. Uygun, *Direct Spatial Resonance in the
  Compressible Boundary Layer on a Rotating Disk* (2006)，PDF 第 4--7 页，
  尤其式 (1)--(15)。
- William 学位论文，PDF 第 32--50 页（第 2 章），尤其式 (2.4)--(2.32)；
  PDF 第 57--59 页的式 (3.9)--(3.14) 仅用于复核线性化。
- R. J. Lingwood, *Absolute instability of the Ekman layer and related
  rotating flows*, J. Fluid Mech. 331 (1997)，用于一般 BEK 的
  `Ro/Co` 与局部尺度约定。

## 8. 候选算子实现与数值 Gate

新前端实施了以下修正：

1. 所有端点统一使用直接 `(U,V,W)`，删除物理层面的 `Ro=-1 -> +1`
   特殊分支；调用旧张量函数时只用 `Ro=0,1` 两点回收其中显式仿射块。
2. 连续性中的 `u/r`、基流法向输运以及相应密度列使用 `+Ro`；能量行的
   `W` 平流与压力功伙伴也使用同一个 `+Ro`。早期试算误用了转换变量的
   `-Ro`，它在 `Ro=0` 不可见、却使 Bödewadt Gate 失败，现已纠正。
3. `Gamma=2RoV+Co` 保留在径向--切向速度耦合及相应密度--基流加速度
   项；`Ro=1` 时 `Co=0` 的直接 Gate 为零。
4. 热 BVP 改为满足两端 Dirichlet 条件的线性三对角离散；采用
   `Mr=r|DeltaOmega|/a_inf`，故 `Ma=Mr/R`，温度重构不再额外乘
   `(Ro*Mr)^2`。

在 `Mr=0.01`、`N=49/59` 的 Lingwood 舍入鼻点上，新候选分别给出

- von Kármán：`alpha=0.38061-0.00055i` (`N=49`)；
- Ekman：`alpha=0.529744-0.000010i`；
- Bödewadt：`alpha=0.482033-0.000354i`。

这三个结果均为低马赫、舍入参考点的近中性 Gate，不是严格
`Ma=0` 证明。尤其压力消元矩阵仍在低马赫下病态，不能把有限
`Mr=0.01` 称为数学极限。

从上述物理模态连续到 `Mr=0.3`，再连续到用户给定 `beta` 并对 `R` 求
`imag(alpha)=0`，得到 `N=69` 的定 `beta` 中性候选：

- Ekman (`Ro=0`, `beta=0.132`):
  `R=119.53588245`, `alpha=0.49836411827+1.59e-11i`；
- Bödewadt (`Ro=1`, `beta=0.127`):
  `R=32.94912489`, `alpha=0.61692247108+2.55e-11i`。

Ekman 的 `N=49/59/69` 中性 `R` 为
`119.62367869/119.48743987/119.53588245`，尚有约 `0.14` 的全跨度；
Bödewadt 为 `32.94912602/32.94912487/32.94912489`，已高度收敛。所有最终
QEP 后向误差均低于 `6e-18`，延拓末步重叠大于
`0.99999999999`。这些是固定 `beta` 的 similarity-subspace 空间中性点，
不是两参数中性曲线鼻点，也不是完整三维稳定性结论。
