# 间隙滤波：定转子腔中扰动穿过旋转核心的理论

本文推导"静盘侧扰动穿过核心到达转盘层"的传递规律，并说明如何用完整数值解检验它。数值验证结果见同目录下的 `report.md` 和本文第 8 节。

---

## 0. 记号与无量纲化

局部平行模型冻结在半径 $r$ 处，记 $R=r/\delta$，其中 $\delta=\sqrt{\nu/\Omega}$；间隙为 $h$，$Re_h=\Omega h^2/\nu=(h/\delta)^2$。算子 `cavity_pencil` 使用的局部单位是：

| 量 | 单位 |
|---|---|
| 长度 | $\delta$（轴向坐标 $\zeta=z^*/\delta=z\sqrt{Re_h}$，其中 $z=z^*/h\in[0,1]$） |
| 速度 | $\Omega r=R\,\Omega\delta$ |
| 时间 | $\delta/(\Omega r)$，因此频率单位为 $R\Omega$ |
| 波数 | $1/\delta$；周向波数 $\beta=n/R$ |

所以有量纲的角速度 $\Omega_*$ 换成局部单位是 $\Omega_*/(R\Omega)$。例如 $2K\Omega\to 2K/R$。

计算在**转盘坐标系**中进行，扰动写成 $\hat q(z)\,e^{i(\alpha r+\beta R\theta-\omega t)}$。外部扰动的频率约定为 $\omega_f=-c\,\beta_f$，$c$ 是图样在转盘系中向后转动的角速度，以 $\Omega$ 为单位。

---

## 1. 核心区的基流

在 Batchelor 流中，两层边界层之间的核心区满足以下条件（Re_h = 1000 的基流在 $0.15\lesssim z\lesssim0.75$ 范围内都近似成立）：

* 周向速度在实验室系中为 $K\Omega r$，即固体旋转。在转盘系中为 $(K-1)\Omega r$，在局部单位下就是 $G-1=K-1$。
* 径向速度 $F\approx0$，没有竖向剪切：$F'\approx G'\approx0$。
* 轴向速度 $W=\Omega\delta\,H$，在局部单位下为 $H/R$。$H\approx-0.77$，换算后只有约 0.002，先忽略，第 6 节再讨论它的影响。

所以，在随核心一起转动、角速度为 $K\Omega$ 的坐标系里，核心流体是**静止**的。扰动只受压力和 Coriolis 力作用，Coriolis 参数是

$$ f = 2K\Omega \quad(\text{局部单位：} f=2K/R). $$

---

## 2. 核心随动系中的线化方程

采用局部笛卡尔坐标：$x$ 沿径向，$y$ 沿周向，$z$ 沿轴向。忽略黏性和曲率，近似条件见第 6 节。离心力并入压力 $p$。

$$
\begin{aligned}
\partial_t u - f v &= -\partial_x p, \\
\partial_t v + f u &= -\partial_y p, \\
\partial_t w &= -\partial_z p, \\
\partial_x u+\partial_y v+\partial_z w &= 0 .
\end{aligned}
$$

### 2.1 正则模态

令 $(u,v,w,p)=(\hat u,\hat v,\hat w,P)(z)\,e^{i(kx+ly-\sigma t)}$，对应关系为 $k\leftrightarrow\alpha$、$l\leftrightarrow\beta$，$\sigma$ 是**核心系中**的频率：

$$
-i\sigma\hat u - f\hat v = -ikP,\qquad
f\hat u - i\sigma\hat v = -ilP,\qquad
-i\sigma\hat w = -P',\qquad
ik\hat u+il\hat v+\hat w'=0 .
$$

### 2.2 解出水平速度

前两式是关于 $\hat u,\hat v$ 的线性方程组，系数行列式为 $f^2-\sigma^2$。解得

$$
\hat u=\frac{\sigma k+ifl}{\sigma^2-f^2}\,P,\qquad
\hat v=\frac{\sigma l-ifk}{\sigma^2-f^2}\,P .
$$

水平散度中，两个 Coriolis 交叉项 $\pm ifkl$ 相互抵消：

$$
ik\hat u+il\hat v=\frac{i\sigma k_h^2}{\sigma^2-f^2}\,P,\qquad k_h^2=k^2+l^2 .
$$

### 2.3 竖向速度与连续性

由第三式得 $\hat w=-iP'/\sigma$，所以 $\hat w'=-iP''/\sigma$。代入连续性方程：

$$
\boxed{\;P''=m^2P,\qquad m^2=k_h^2\,\frac{\sigma^2}{\sigma^2-f^2}\;}
$$

由 2.2 的结果，$\hat u,\hat v$ 和 $P$ 只差一个常数因子，所以**水平速度满足同一个方程**。数值检验时直接用速度，就不必处理 $P_N$–$P_{N-2}$ 格式下的压力插值问题。

---

## 3. 转盘系频率到核心系频率：Doppler 频移

设坐标系 B 相对坐标系 A 以角速度 $\Omega_B$ 转动，则 $\theta_A=\theta_B+\Omega_B t$，于是

$$ e^{i(n\theta_A-\omega_A t)}=e^{i(n\theta_B-(\omega_A-n\Omega_B)t)}\;\Rightarrow\;\omega_B=\omega_A-n\Omega_B . $$

核心相对转盘的转速是 $\Omega_B=(K-1)\Omega$，换成局部单位（$n\Omega_B\to\beta(K-1)$），得到

$$
\boxed{\;\sigma=\omega+\beta(1-K)\;}\qquad\xrightarrow{\;\omega_f=-c\beta_f\;}\qquad \sigma_f=\beta_f\,(1-K-c).
$$

**直观理解**：σ 等于周向波数乘以"图样相对核心的转速"。图样在实验室系中以 $(1-c)\Omega$ 转动，核心以 $K\Omega$ 转动，两者之差为 $(1-c-K)\Omega$。

**和代码的对应**：`cavity_pencil` 中的 Doppler 项是 $-i[\omega-\beta(G-1)]$，在核心区 $G=K$；Coriolis 交叉项是 $\mp2(G+1)/R$，其中 $G$ 取转盘系的值，在核心区就是 $2K/R=f$。两者和上面的推导一致。

---

## 4. 三种区间

| 区间 | $m$ | 竖向结构 | 物理 |
|---|---|---|---|
| $\lvert\sigma\rvert\gg f$ | $m\approx k_h$ | $e^{\pm k_h\zeta}$ | 旋转可以忽略，扰动按势流（Laplace）衰减 |
| $f<\lvert\sigma\rvert\lesssim 2f$ | $m>k_h$ | 指数衰减 | 衰减比势流**更快**，$\lvert\sigma\rvert\to f^+$ 时 $m\to\infty$ |
| $\lvert\sigma\rvert<f$ | $m=i\mu$ 为纯虚数 | $\cos\mu\zeta$ | 惯性波，沿竖向传播，不衰减；群速度和转轴的夹角满足 $\cos\vartheta=\lvert\sigma\rvert/f$ |
| $\sigma\to0$ | $m\to0$ | 和 $\zeta$ 无关 | Taylor–Proudman 柱：扰动沿 z 不变，**完全透过** |

---

## 5. 传递函数

把核心看作 $\zeta_r<\zeta<\zeta_s$ 之间的一层流体，上下两端分别是转盘层外缘和静盘层外缘。到领头阶，转盘层对核心相当于一个**刚性盖子**：$\hat w(\zeta_r)=0$，也就是 $P'(\zeta_r)=0$。所以

$$ P(\zeta)=A\cosh\!\big(m(\zeta-\zeta_r)\big), $$

$$
\boxed{\;T\equiv\frac{|u_h(\zeta_r)|}{|u_h(\zeta_s)|}=\frac{1}{\lvert\cosh(m\,\Delta\zeta)\rvert},\qquad \Delta\zeta=(z_s-z_r)\sqrt{Re_h}\;}
$$

* 衰减区（$m$ 为实数）：$T\simeq 2e^{-m\Delta\zeta}$。**指数部分和"盖子"的具体模型无关**，前置系数 2 才依赖于盖子条件。
* 惯性波区（$m=i\mu$）：$T=1/\lvert\cos(\mu\Delta\zeta)\rvert$。当 $\mu\Delta\zeta=(j+\tfrac12)\pi$ 时出现间隙共振，无黏极限下发散，实际峰宽由黏性和 Ekman 阻尼决定。
* $\sigma=0$：$T=1$。

### 5.1 盖子条件的修正

实际上，转盘层会通过 Ekman 抽吸产生 $w_E\sim\tfrac12\delta_E\,\omega_z'$，其中 $\omega_z'\sim k_h u_h$，所以 $w_E/u_h\sim k_h\delta_E$（以 $\delta$ 为单位，$\delta_E=O(1)$）。衰减区中核心模态本身的 $w/u_h\sim m/k_h=O(1)$。因此，当 $k_h\delta_E$ 不小时，盖子条件会修正 $T$ 的 **O(1) 前置系数**，但不会改变指数 $e^{-m\Delta\zeta}$。

---

## 6. 被忽略的项及其量级（局部单位，Re_h = 1000，R = 330）

| 项 | 量级 | 和 $\sigma$ 相比 |
|---|---|---|
| 黏性 $(k_h^2-m^2)/R$ | $\lesssim k_h^2/330$ | 直接路径（$\sigma\approx0.057$）中只占 1e-3 到 1e-2；但 $\sigma\to0$ 时变得重要 |
| 轴向流输运 $(H/R)\,m$ | $0.0023\,m$ | $\sigma\to0$ 时重要，否则很小 |
| 曲率 $u/R$（连续性方程中） | $1/R\approx0.003$ | 和 $k_h$ 相比，在 $k_hR\gg1$ 时可以忽略 |
| 核心剪切 $F',G'$ | 在 $0.15\lesssim z\lesssim0.75$ 很小 | 靠近边界层外缘时增大 |

把这些项全部保留，核心区就是一个**常系数**问题。代入 $q\propto e^{\lambda\zeta}$，得到一个 4×4 二次特征值问题 $M(\lambda)=M_0+\lambda M_1+\lambda^2M_2$，可以算出它的**精确根**（`analyze_gap_filter.jl` 第 A 部分）。在无黏、$H=0$、不计曲率的极限下，模最小的两个根退化为 $\pm m$；其余根的模都很大，是只存在于薄层内的黏性根。

---

## 7. 对两条激发路径的含义

* **直接路径**：转盘模态的 $k_h=\sqrt{\alpha_n^2+\beta_t^2}\approx0.29\sim0.45$，频率比 $\sigma_t/f\gg1$，所以
  $$ \Lambda_{dir}\propto e^{-m\Delta\zeta}\approx e^{-k_h\sqrt{Re_h}\,(z_s-z_r)} .$$
  在 Re_h = 1000 下，Δζ ≈ 26，所以 $e^{-m\Delta\zeta}=10^{-3}\sim10^{-5}$。**$\ln\lvert\Lambda_{dir}\rvert$ 应该随 $\sqrt{Re_h}$ 线性下降，斜率为 $-k_h(z_s-z_r)$**，这是一个可以检验的定量预测。
* **散射路径**：入射扰动的 $k_h\Delta\zeta\lesssim1$（长波），或者 $\lvert\sigma_f\rvert<f$（与核心近似同步转动），此时 $T=O(1)$。粗糙元只需要提供波数差 $\kappa=\alpha_n-\alpha_f$。
* **结论**：两条路径之比的主导因子是 $e^{+m\Delta\zeta}$，主参数是 $\sqrt{Re_h}$，也就是间隙高度与边界层厚度之比。湍流模型、源项形状、基流修正这些因素，影响的都只是 O(1) 的前置系数。

---

## 8. 数值验证

见 `report.md`（由 `analyze_gap_filter.jl` 生成），对应关系如下：

* A：精确常系数根和无黏 $m$ 的比较；
* B：从完整数值解中拟合核心区指数；
* C：端到端传递比和 $1/\lvert\cosh m\Delta\zeta\rvert$ 的比较；
* D：Re_h 扫描；
* E：$\ln\lvert\Lambda_{dir}\rvert$ 和 $-m\Delta\zeta$ 的比较；
* F：识别共振峰。

验证结论汇总在本文末尾（第 9 节）。

---

## 9. 验证结论（2026-11-05）

脚本：`analyze_gap_filter.jl`（A–F 部分）、`analyze_gap_filter_core.jl`（核心区检验）、`export_gap_filter_profiles.jl`、`plot_gap_filter.py`。图：`gap_filter_verification.png`。

### 9.1 已证实的结论

1. **衰减区（$\lvert\sigma\rvert/f>2$）的竖向衰减律成立。** 在真正的核心区 $0.2\le z\le0.5$ 内，数值的 $\ln\lvert u(z_b)/u(z_a)\rvert$ 和带刚性盖的理论 $\ln\lvert\cosh m(\zeta_b-\zeta_r)/\cosh m(\zeta_a-\zeta_r)\rvert$ 相比，相对误差中位数为 8%，共 240 个工况。偏离较大的 16 个工况全部是 $n>0$、$\alpha\gtrsim0.3$，正好落在转盘 Type-I 模态的波数带内。我的推测是：在那里转盘层接近共振，不再是一个被动的"盖子"。这一点还没有单独验证。
2. **直接路径严格满足指数律。** 对 Type-I 目标点（R = 330，n = 16，ω̄ = 8），Re_h 从 1000 增大到 4000 时，$\ln\lvert\Lambda_{dir}\rvert$ 从 −10.8 变到 −25.0，而 $\ln\lvert\Lambda_{dir}\rvert+m\Delta\zeta=-0.48\pm0.01$ 保持不变：
   $$\Lambda_{dir}=0.62\,e^{-m\Delta\zeta},\qquad m=0.451\approx k_h .$$
   其中 Δζ 取两层外缘之间的距离（转盘侧取 G 的第一个极小，静盘侧取 G 的超调极大）。源项的位置和宽度按 δ 单位固定，以免源项本身随 Re_h 变化。
3. **常系数精确根和无黏 m 的关系**：在衰减区两者相差几个百分点，不对称的部分来自轴向流 H。在 $\sigma\to0$ 时，无黏的 $m=0$ 被黏性、H 和曲率共同改变为 $\pm0.03\text{–}0.04\,i$。

### 9.2 需要修正的结论

4. **"Taylor 柱完全透过"不成立。** 惯性区（$\lvert\sigma\rvert/f<0.8$）的 64 个工况里，数值解在核心区明显衰减，而且**每个工况在间隙内都存在临界层** $\sigma(z)=\beta(1-G(z)-c)=0$。原因有两个：
   * 在 Re_h = 1000 下，静盘 Bödewadt 层的超调振荡一直延伸到 z ≈ 0.55，"核心"只有 0.15 ≲ z ≲ 0.5，G 有 ±0.01 到 ±0.05 的变化；
   * σ → 0 时，σ 本身很小，G 的微小变化就会让 σ(z) 变号，产生临界层。
   
   临界层会吸收或散射波，类似于旋转流体中的 Booker–Bretherton 机制。所以惯性区的透过特性**不能用均匀核心理论描述**，需要考虑 G(z) 的变化（用 WKB 方法加临界层分析）。这是理论上的一个开放问题。
5. **对主线的影响**：散射路径的入射场之所以能穿过间隙，可靠的原因是**长波**（$k_h\Delta\zeta\lesssim1$，此时无论 σ 取何值 T = O(1)），而不是"和核心同步转动"。和核心接近同步转动的扰动反而会遇到临界层。

### 9.3 对研究主线的意义

* 直接路径的强度 $\propto e^{-m\Delta\zeta}$，其中 $m\simeq k_h(\alpha_n,\beta_t)$，主参数是 $\sqrt{Re_h}$。这是一条已经被验证的普适规律。
* 两条路径之比 $\sim\varepsilon_h\,e^{+m\Delta\zeta}\times O(1)$。湍流模型、源项形状等因素只影响 O(1) 的前置系数，这和前面的敏感性分析一致：改变它们时 Λ 最多变 2 倍。

### 9.4 补充与更正：σ/f 和临界层是两个独立的因素

* **衰减区的算例同样有临界层**，只是位于静盘层内（c = 0.75–1.0 时 z = 0.94–1.0，紧挨源项）。9.2 节中"0/240"统计的只是核心路径 0.12 < z < 0.9 上的临界层。以当地流速对流的扰动，源项必然位于自己的临界层附近。
* **对照算例**（ln|u(0.5)/u(0.2)|，Re_h = 1000，R = 330）：
  | 类别 | 参数 | σ/f | 数值 | 无黏 | 精确根 |
  |---|---|---|---|---|---|
  | 衰减区，核心有临界层（z = 0.80、0.88） | α = 0.2/0.3，n = −20，c = 0.62 | −2.1 | 2.47 / 2.51 | 1.99 / 3.14 | 2.19 / 3.23 |
  | 惯性区，核心无临界层 | α = 0.1/0.2/0.3，n = −5，c = 0.75 | 0.51 | 2.49 / 2.51 / 4.16 | −0.27 / −1.90 / −0.53 | 0.92 / 2.16 / 3.15 |
  | 惯性区，核心有临界层 | α = 0.2，n = −10，c = 0.70 | 0.23 | 2.29 | −0.17 | −0.38 |
* **解释**：当 $|\sigma|\gg f$ 时，$m^2/k_h^2=\sigma^2/(\sigma^2-f^2)\approx1$，对 σ 不敏感，扰动本质上是势流衰减，临界层只影响源项附近，进入 O(1) 前置系数。当 $|\sigma|\lesssim f$ 时，任何 O(f) 的频率偏移都会改变竖向结构。在局部单位下，$f=0.0019$，而核心 G 的 ±0.01 变化、核心轴向流输运 $W m$（$W=H/R\approx0.0023$）、黏性 $k_h^2/R$ 分别相当于 f 的约 16%、24%、6%。所以惯性区的失效**首先来自核心轴向流和黏性**（精确根理论能解释其中一部分），**路径上的临界层**（$\xi=w/(-i\sigma)\to\infty$，波与平均剪切强烈作用而被吸收）会使情况进一步恶化。9.2 节中"全部归因于临界层"的说法据此更正。

### 9.5 σ/f 刚刚大于 1 的情况：无黏奇异性被抹平

无黏理论给出 $m/k_h=g(\sigma/f)=x/\sqrt{x^2-1}$，在 $x\to1^+$ 时发散（σ/f = 1.02 时 g = 5.1）。数值检验（n = −10，c 取在核心速度之上，使核心路径上没有临界层，α = 0.03–0.12，数据见 `near_inertial.csv`）表明：

| σ/f | g（无黏） | α = 0.12 时：无黏 m / 精确根 Re λ |
|---|---|---|
| 1.02 | 5.07 | 0.63 / 0.16 |
| 1.10 | 2.40 | 0.30 / 0.16 |
| 1.50 | 1.34 | 0.17 / 0.14 |
| ≥ 2 | ≤ 1.15 | 两者一致 |

精确根始终只有 k_h 的约 1.3 倍，数值解也没有出现强衰减。原因是奇异性来自 $\hat u\propto P/(\sigma^2-f^2)$ 的分母为零，而核心轴向流的输运 $W m$（约为 f 的 12%）和黏性使有效频率成为复数，$\sigma_{eff}^2-f^2$ 不会等于零。

**结论**：在真实核心中，整个衰减区（σ/f ≳ 1）的透过率都可以写成 $T\approx1/\cosh(k_h\Delta\zeta)$，修正因子 $g_{eff}\lesssim1.5$。"长波透过"和"直接路径按指数衰减"是同一条规律的两个极限。σ/f 只在惯性区（σ/f < 1）中真正起决定作用，那里的行为不具普适性。
