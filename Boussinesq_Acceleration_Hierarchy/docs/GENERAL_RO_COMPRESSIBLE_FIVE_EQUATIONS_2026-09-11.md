# 一般 (Ro) 的可压缩 BEK 五扰动方程

> **2026-09-12 勘误：本文档已被
> `COMPRESSIBLE_BEK_RO_CO_OPERATOR_AUDIT_2026-09-12.md` 取代。**
> 本文正确识别了 `Gamma=2RoG+Co` 以及 Doppler 项不整体乘 `Ro`，但遗漏了
> 一般 BEK 连续性曲率项的 `Ro`、能量行法向输运的 `Ro`，并混用了
> `(F,G,H)=(-U,V,-W)` 与直接 Lingwood `(U,V,W)` 约定。以下内容只保留为
> 研究历史，不得再作为一般 `Ro` 数值计算的模型定义。

本文档固定当前代码所采用的变量、符号和验证范围。它不是把 (Ro=-1) 的 von Kármán 方程简单改名，而是说明哪些项来自旋转坐标中的基流惯性，哪些项来自曲率/Coriolis 耦合。

## 1. 扰动变量与状态方程

取局部波形

\[
 (\hat\rho,\hat u,\hat v,\hat w,\hat T)(y)\,
 e^{i(\alpha r+\beta\theta-\omega t)},
\]

并令

\[
 \mathcal U=\alpha F+\beta G-\omega,
 \qquad Co=2-Ro-Ro^2,
 \qquad \Gamma=2Ro\,G+Co .
\]

这里 (F,G,H) 是统一符号约定下的径向、切向和轴向基流系数，
\(F(\infty)=0, G(\infty)=1\)，而 (H) 是轴向质量流函数对应的剖面。Chapman 定律给出

\[
 \hat\mu=\hat T,\qquad \hat\lambda=-\frac23\hat T,
 \qquad \hat k=\frac{\hat T}{Pr},
 \qquad \gamma Ma^2\hat p=T\hat\rho+\rho\hat T .
\]

因此未知量只有五个：
\(\hat\rho,\hat u,\hat v,\hat w,\hat T\)。

## 2. 五个方程的结构

将旋转坐标中的连续性、径向动量、切向动量、轴向动量和能量方程线性化，并以 Dorodnitsyn--Howarth 坐标 (y) 表示，得到

\[
 \mathcal C(\hat\rho,\hat u,\hat v,\hat w)=0,
\]

\[
 \mathcal M_r(\hat\rho,\hat u,\hat v,\hat w,\hat T)=0,
 \qquad
 \mathcal M_\theta(\hat\rho,\hat u,\hat v,\hat w,\hat T)=0,
\]

\[
 \mathcal M_z(\hat\rho,\hat u,\hat v,\hat w,\hat T)=0,
 \qquad
 \mathcal E(\hat\rho,\hat u,\hat v,\hat w,\hat T)=0.
\]

在矩阵实现中，这五式分别对应 `Ta`, `A`, `B`, `C`, `D`, `Vxx`, `Vxy`, `Vxz`, `Vyz` 等块。最容易出错、也是一般 (Ro) 的判别性项如下：

* 连续性中的波数项为 (i\alpha R\rho\hat u+i\beta R\rho\hat v)，法向项为 (R\rho(D\rho+\rho D)\hat w)，基流密度扰动耦合为 (iR\mathcal U\hat\rho+2F\hat\rho+\rho DH\hat\rho+\rho H D\hat\rho)。
* 径向动量的基流法向惯性为 (Ro\,\rho F\hat u)，切向速度耦合为 (-\rho\Gamma\hat v)。
* 切向动量的径向速度耦合为 (+\rho\Gamma\hat u)，切向对角惯性为 (Ro\,\rho F\hat v)。
* 三个动量方程中的法向平流项均带 (Ro\,\rho^2H D)。
* 能量方程保留压力功、热传导、黏性耗散和基流温度梯度；不能把热方程替换成被动标量方程。

换言之，(Ro) 不应乘在 Doppler 项
\(\mathcal U\) 前面；它只出现在由旋转坐标基流惯性和法向相似缩放产生的项中。曲率/Coriolis 作用统一通过
\(\Gamma=2RoG+Co\) 出现。

## 3. 与 2006 年方程的对应

当 (Ro=-1)、(Co=2) 时，\(\Gamma=2(1-G)\)，在本文的 (G(\infty)=1) 约定下，经过基流符号转换后恢复 Turkyilmazoglu--Uygun 方程 (10)--(15) 的 von Kármán 形式。代码中的 `Spatial_mode_BEK` 仍承担逐项微分算子装配；新模块只负责统一 (Ro) 插值、基流坐标和边界消元。

## 4. 当前验证边界

已通过的只是：

1. (Ro=-1) 的可压缩 von Kármán 基准；
2. 张量装配与展开装配的数值一致性；
3. QEP 后向误差达到约 (10^{-12}\)--(10^{-14})；
4. (Ro\ne-1) 的低马赫诊断显示负 (Ro) 分支显著改善。

尚未宣称一般 (Ro) 理论验证完成。Lingwood 表中的五个鼻点是不可压缩、四舍五入的参考值，必须通过 (Ma\to0) 模态跟踪和网格收敛后才能作为正式中性曲线证据。
