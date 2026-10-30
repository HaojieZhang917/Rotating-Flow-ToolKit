# 自由来流涡扰动模型及其在动盘旋转系中的推导

日期：2026-10-26。研究背景：等温、常密度、低马赫数定转子系统；第一阶段没有实际静盘湍流数据，采用受控外部扰动。动盘与核心角速度分别记为 Ω_d、Ω_c，均为有符号常数，旋转轴为 z。以下先给出文献中的输入模型，再从 NS 推导旋转系模型。新的方程推导不应误认为下列论文已直接给出同一转静腔模型。

## 1. 文献中的几类模型

| 表述 | 如何定义输入 | 主要用途与限制 | 原始研究 |
|---|---|---|---|
| 单一 Fourier 涡模态 | 给波数、两个横向极化自由度、幅相；满足散度约束和外部流的色散关系 | 最适合受控感受性；单模态本身不是湍流 | Guo & Kloker 2026，§2.3.5 |
| 随机 Fourier 谱叠加 | 叠加许多散度为零的涡模态，指定能谱/谱张量和相位 | 表示弱外部湍流；各向同性不是腔体的默认事实 | Leib, Wundrow & Goldstein 1999 |
| OS/Squire 连续谱 | 用外部呈振荡、内部满足边界层动力学的连续谱特征函数构造输入 | 可同时描述外部扰动和剪切层响应；不同于离散失稳模态 | Jacobs & Durbin 2001；Zaki & Durbin 2005 |
| Kelvin/快速畸变形式 | 波矢与幅值随平均速度梯度演化 | 用于非均匀外流；下文直接推导其旋转系表达 | 本文由线性 NS 推导；不作为已完成的腔体模型验证 |

Leib 等研究表明，低频各向异性可以显著影响边界层响应，因此 Tu 和一条一维能谱不能唯一规定输入。Guo–Kloker 使用单 Fourier 分量和局部体力，研究外部涡扰动及其与粗糙的耦合；其几何是扫掠翼。连续谱模型则用于使输入与剪切层相容，不能把同名平板 OS/Squire 模型不经推导就替换为旋转盘算子。

## 2. 非旋转、均匀外流中的标准涡模态

令外流平均速度 U_e 为常向量，压力扰动用运动学压力 π=p'/ρ 表示。线性方程为

\[
(\partial_t+\boldsymbol U_e\cdot\nabla)\boldsymbol u'
=-\nabla\pi+\nu\nabla^2\boldsymbol u',\qquad
\nabla\cdot\boldsymbol u'=0.
\tag{1}
\]

取实波矢 k≠0，并写

\[
(\boldsymbol u',\pi)=\Re\{(\boldsymbol a,\hat\pi)
e^{i(\boldsymbol k\cdot\boldsymbol x-\omega t)}\}.
\tag{2}
\]

代入得

\[
\boldsymbol k\cdot\boldsymbol a=0,\qquad
[-i(\omega-\boldsymbol k\cdot\boldsymbol U_e)+\nu k^2]\boldsymbol a
=-i\boldsymbol k\hat\pi.
\tag{3}
\]

点乘 k 可得 π̂=0，非平凡横向幅值要求

\[
\boxed{\omega=\boldsymbol k\cdot\boldsymbol U_e-i\nu k^2.}
\tag{4}
\]

因此无黏极限为对流涡波 ω=k·U_e；保留黏性时幅值按 exp(−νk²t) 衰减。π̂=0 只对这一均匀、无外加强迫模型成立，不能推广到有剪切或旋转的全部涡扰动。

令 k_h=(k_x²+k_y²)^(1/2)>0，构造

\[
\boldsymbol e_1=\frac{(-k_y,k_x,0)}{k_h},\qquad
\boldsymbol e_2=\frac{\boldsymbol k}{k}\times\boldsymbol e_1,
\quad \boldsymbol a=A_1\boldsymbol e_1+A_2\boldsymbol e_2.
\tag{5}
\]

两个复数 A_1,A_2 给出独立极化；不能把三个速度幅值任意独立指定。当 k_h=0 时任选两个与 k 垂直的正交方向。

实 ω、实 k、ν>0 三者一般不能同时描述自由、恒幅时谐模态。若需要连续时谐输入，应明确它由源维持，或采用空间衰减的复波数问题。

## 3. 弱湍流的谱模型

线性范围可以叠加上述模态。采用 k_turb=(1/2)∫tr Φ(k)d³k 的约定时，各向同性谱张量可写为

\[
\Phi_{ij}(\boldsymbol k)=\frac{E(k)}{4\pi k^2}
\left(\delta_{ij}-\frac{k_i k_j}{k^2}\right),\qquad
\int_0^\infty E(k)dk=k_{\rm turb}.
\tag{6}
\]

这里的 k_turb 是湍动能，不是波数。对各向异性模型，更一般地可令 Φ=P_k C(k) P_k，其中 P_k=I−kkᵀ/k²，C 为半正定 Hermitian 矩阵。应施加 a(−k)=a(k)* 的实场约束。仅给初始随机相位并自由衰减不会得到统计定常湍流；定常输入需持续强迫或连续边界注入。

## 4. 坐标变换：动盘旋转不等于外部流体旋转

令 Q_d(t) 为绕 z 轴转过 Ω_dt 的正交矩阵，实验室系和动盘系分别以 I、d 标记：

\[
\boldsymbol x_I=Q_d\boldsymbol x_d,\qquad
\boldsymbol U_d=Q_d^T\boldsymbol U_I(Q_d\boldsymbol x_d,t)
-\boldsymbol\Omega_d\times\boldsymbol x_d.
\tag{7}
\]

定义约化运动学压力

\[
\Pi_d=\frac{p_I(Q_d\boldsymbol x_d,t)}{\rho}
-\frac12\Omega_d^2r^2.
\tag{8}
\]

恒速旋转系 NS 为

\[
\partial_t\boldsymbol U_d+(\boldsymbol U_d\cdot\nabla)\boldsymbol U_d
+2\boldsymbol\Omega_d\times\boldsymbol U_d
=-\nabla\Pi_d+\nu\nabla^2\boldsymbol U_d,
\quad \nabla\cdot\boldsymbol U_d=0.
\tag{9}
\]

离心项已包含于 Π_d，不能再次加入。微扰的离心势没有独立变化，π'=p'/ρ。

柱坐标 θ_d=θ_I−Ω_dt；径向和轴向总速度不变，周向总速度减去 Ω_dr。由于基流也减去同一项，对应点处的三个柱坐标扰动分量不需要再减 Ω_dr。纯角向模态的相位变换为

\[
e^{i(m\theta_I-\omega_I t)}
=e^{i[m\theta_d-(\omega_I-m\Omega_d)t]},
\quad
\boxed{\omega_d=\omega_I-m\Omega_d.}
\tag{10}
\]

Cartesian 单平面波在动盘系中的精确形式则是

\[
\boldsymbol u'_d=\Re\{Q_d^T\boldsymbol a_I
e^{i[(Q_d^T\boldsymbol k_I)\cdot\boldsymbol x_d-\omega_I t]}\}.
\tag{11}
\]

其波矢和矢量幅值都随坐标轴转动；通常不能将其视为一个固定 m、固定波矢的全盘模态。式(10)是纯柱坐标角向 Fourier 分量的频移，不是对任意 Cartesian 平面波直接加一个常数频移。

## 5. 任意轴对称基流上的旋转系线性方程

令动盘系基流为 U_0=(U(r,z),V(r,z),W(r,z))，扰动为 (u,v,w,π)。线性化式(9)：

\[
(\partial_t+\boldsymbol U_0\cdot\nabla)\boldsymbol u
+(\boldsymbol u\cdot\nabla)\boldsymbol U_0
+2\boldsymbol\Omega_d\times\boldsymbol u
=-\nabla\pi+\nu\nabla^2\boldsymbol u,
\quad\nabla\cdot\boldsymbol u=0.
\tag{12}
\]

以下给出精确柱坐标形式，避免将矢量 Laplacian 错写成逐分量标量 Laplacian。取 exp(imθ_d−iω_dt)，令

\[
\mathcal D_m=-i\omega_d+U\partial_r+W\partial_z+imV/r,
\quad
\Delta_m=\partial_{rr}+r^{-1}\partial_r-m^2/r^2+\partial_{zz}.
\]

则

\[
\mathcal D_m\hat u+U_r\hat u+U_z\hat w
-2(V/r+\Omega_d)\hat v
=-\partial_r\hat\pi
+\nu[(\Delta_m-r^{-2})\hat u-2im\hat v/r^2],
\tag{13a}
\]

\[
\mathcal D_m\hat v+(V_r+V/r+2\Omega_d)\hat u
+V_z\hat w+(U/r)\hat v
=-im\hat\pi/r
+\nu[(\Delta_m-r^{-2})\hat v+2im\hat u/r^2],
\tag{13b}
\]

\[
\mathcal D_m\hat w+W_r\hat u+W_z\hat w
=-\partial_z\hat\pi+\nu\Delta_m\hat w,
\tag{13c}
\]

\[
(\partial_r+r^{-1})\hat u+im\hat v/r+\partial_z\hat w=0.
\tag{13d}
\]

旋转核心、动盘边界层或非相似全局基态都可代入该式；它没有采用局部平行假设。若 U_0 是湍流平均场而非精确基本解，应另外说明冻结平均流和 Reynolds 应力闭合假设。

## 6. 近刚体旋转核心：为何恢复力由 Ω_c 决定

假定实验室系的理想核心流为 U_I,0=Ω_cr e_θ+W_c e_z，其中 W_c 为可选常轴向对流。动盘系中

\[
U=0,\qquad V=(\Omega_c-\Omega_d)r,\qquad W=W_c.
\tag{14}
\]

基压满足 ∂_rΠ_0=(Ω_c²−Ω_d²)r。定义

\[
\mathcal D_c=\partial_t+(\Omega_c-\Omega_d)\partial_\theta+W_c\partial_z,
\quad \Delta=\partial_{rr}+r^{-1}\partial_r+r^{-2}\partial_{\theta\theta}+\partial_{zz}.
\]

在式(13)中代入，得到精确核心方程

\[
\mathcal D_c u-2\Omega_c v=-\pi_r
+\nu(\Delta u-u/r^2-2v_\theta/r^2),
\tag{15a}
\]

\[
\mathcal D_c v+2\Omega_c u=-\pi_\theta/r
+\nu(\Delta v-v/r^2+2u_\theta/r^2),
\tag{15b}
\]

\[
\mathcal D_c w=-\pi_z+\nu\Delta w,
\qquad r^{-1}(ru)_r+r^{-1}v_\theta+w_z=0.
\tag{15c}
\]

例如径向的两项为 −2(V/r)v−2Ω_dv=−2Ω_cv；周向同理。这是平均旋转和科里奥利项共同作用的结果，不能只留下 2Ω_d。

两个必须满足的极限：

- Ω_c=0：惯性波恢复项消失；换坐标不会给物理静止外流制造惯性波。
- Ω_c=Ω_d：核心相对动盘静止，此时恢复项确实是 2Ω_d。

## 7. 刚体核心中的横向极化和色散关系

在随核心旋转的 Cartesian 系中，平均流仅为 W_ce_z。该系内可精确使用固定波矢平面波（没有曲率近似）。令频率 ω_c、内禀频率 σ=ω_c−k_zW_c，则

\[
[-i\sigma+\nu k^2]\boldsymbol a
+2\Omega_c\boldsymbol e_z\times\boldsymbol a
=-i\boldsymbol k\hat\pi,\qquad \boldsymbol k\cdot\boldsymbol a=0.
\tag{16}
\]

点乘 k 得

\[
\hat\pi=\frac{2i\Omega_c}{k^2}
\boldsymbol k\cdot(\boldsymbol e_z\times\boldsymbol a).
\tag{17}
\]

与普通均匀自由来流不同，旋转核心的涡性波一般带有压力扰动。

用式(5)的横向基写 a=A_1e_1+A_2e_2，记 λ=−iσ+νk² 和 f_k=2Ω_ck_z/k。投影后为

\[
\begin{pmatrix}\lambda&-f_k\\f_k&\lambda\end{pmatrix}
\begin{pmatrix}A_1\\A_2\end{pmatrix}=0.
\tag{18}
\]

因此 λ²+f_k²=0。可选择支号 s=±1 及单位复极化

\[
\boxed{\boldsymbol a_s=(\boldsymbol e_1-is\boldsymbol e_2)/\sqrt2,
\qquad \sigma_s=s\,2\Omega_c\frac{k_z}{k}-i\nu k^2.}
\tag{19}
\]

当 Ω_c=0 时两支退化为普通的两个横向涡扰动自由度；当 k_z=0 时出现零振荡频率的慢分量，不能将所有涡性扰动都解释成有非零频率的惯性波。

这与经典惯性波色散关系一致，旋转网格湍流实验也观察到该类波动，但本推导不意味着真实腔体输入全部由线性自由波组成。

## 8. 动盘系内的两种等价表示

### 8.1 纯柱坐标角向模态

对 exp(imθ_d+ik_zz−iω_dt)，其内禀频率为

\[
\boxed{\sigma=\omega_d-m(\Omega_c-\Omega_d)-k_zW_c
=\omega_I-m\Omega_c-k_zW_c.}
\tag{20}
\]

局部短波近似下，k=(k_r,m/r_0,k_z)，由式(19)得

\[
\boxed{\omega_{d,s}=m(\Omega_c-\Omega_d)+k_zW_c
+s\frac{2\Omega_c k_z}{\sqrt{k_r^2+(m/r_0)^2+k_z^2}}
-i\nu[k_r^2+(m/r_0)^2+k_z^2].}
\tag{21}
\]

式(21)是局部波包关系，要求波长相对曲率/基流变化尺度足够小。不能用常幅 exp(ik_rr+imθ+ik_zz) 在整个圆盘冒充精确 Cartesian 模态：精确连续性还包含 u/r。

对于无黏自由惯性波，必要频带条件是 |σ|≤2|Ω_c|。这里约束的是核心相对频率，不是直接约束实验室频率或动盘频率。真实强迫场、衰减场、剪切核心和有限腔模态应由相应方程决定，不能由这个必要条件宣称频带外完全没有扰动。

### 8.2 核心平面波的精确坐标变换

令 Q_Δ=Q[(Ω_c−Ω_d)t]，核心系波矢 k_0 和极化 a_s 固定，则动盘系 Cartesian 表示为

\[
\boldsymbol k_d(t)=Q_\Delta\boldsymbol k_0,\quad
\boldsymbol a_d(t)=Q_\Delta\boldsymbol a_s e^{-i\omega_ct},\quad
\pi_d(t)=\hat\pi_s e^{-i\omega_ct}.
\tag{22}
\]

把 a_d(t)exp(ik_d(t)·x_d) 的实部作为速度即可得到同一物理波，不需要局部平行假设。但该表示不等同于一个纯 m；需要纯 m 输入时，应使用式(13)/(15)的柱坐标径向解。

## 9. 精确柱坐标输入：Bessel 表示或矢量势

在无黏、W_c 常值、σ≠0 且 σ²≠4Ω_c² 时，以压力振幅 p_m(r) 表示式(15)的纯 m 模态：

\[
\hat u=\frac{i(\sigma p_m'-2\Omega_c m p_m/r)}{4\Omega_c^2-\sigma^2},\quad
\hat v=\frac{2\Omega_c p_m'-\sigma m p_m/r}{4\Omega_c^2-\sigma^2},\quad
\hat w=\frac{k_zp_m}{\sigma}.
\tag{23}
\]

连续性给出

\[
p_m''+\frac1r p_m'-\frac{m^2}{r^2}p_m+k_h^2p_m=0,
\quad
k_h^2=k_z^2\left(\frac{4\Omega_c^2}{\sigma^2}-1\right).
\tag{24}
\]

若 k_h 实数，轴上正则的解为 J_|m|(k_hr)，环形域还可使用 Y_|m|。有限容器边界会进一步选择允许的径向结构和频率；这里未施加容器条件。式(23)不用于 σ=0 或 |σ|=2|Ω_c| 等退化情形，应回到式(15)求解。

若只希望构造严格散度为零的受控输入，可先用矢量势

\[
\boldsymbol u=\nabla\times(\psi\boldsymbol e_z)
+\nabla\times\nabla\times(\phi\boldsymbol e_z).
\tag{25}
\]

令势函数含 exp(imθ+ik_zz)，则

\[
u=im\psi/r+ik_z\phi',\quad
v=-\psi'-mk_z\phi/r,\quad
w=-(\phi''+\phi'/r-m^2\phi/r^2).
\tag{26}
\]

该构造满足连续性，但还要满足动量方程或明确维持它的体力。局部包络宜放在势函数内；直接将一个散度为零的速度乘任意包络通常会破坏散度约束。

## 10. 非刚体核心的 Kelvin/快速畸变扩展

局部取平均流 U_0=A(t)x（跟随平均流轨迹，已消去均匀平移），令 u=a(t)exp[ik(t)·x]。在任意恒速 Ω_d 系中，从式(12)逐项比较 x 的系数得

\[
\dot{\boldsymbol k}=-A^T\boldsymbol k,
\quad
\dot{\boldsymbol a}=-A\boldsymbol a-2\boldsymbol\Omega_d\times\boldsymbol a
-i\boldsymbol k\hat\pi-\nu k^2\boldsymbol a.
\tag{27}
\]

由 d(k·a)/dt=0 得

\[
\hat\pi=\frac{2i}{k^2}\boldsymbol k\cdot
(A\boldsymbol a+\boldsymbol\Omega_d\times\boldsymbol a).
\tag{28}
\]

若保留均匀平移 b(t)，幅值方程还要加 −i(k·b)a。刚体核心对应 A=(Ω_c−Ω_d)J_z，式(27)给出 k 随相对角速度转动，式(28)恢复式(17)。对有剪切的核心，固定波矢和简单惯性波频率一般不再成立。

## 11. 无量纲式与建议的第一组算例

以任意指定长度 ℓ、时间 Ω_d^−1（Ω_d>0）无量纲化。令 χ=Ω_c/Ω_d、Re_ℓ=Ω_dℓ²/ν、κ=ℓk、W_c*=W_c/(Ω_dℓ)、ω*=ω_d/Ω_d，则式(21)为

\[
\omega_s^*=m(\chi-1)+\kappa_zW_c^*
+s\,2\chi\frac{\kappa_z}{|\boldsymbol\kappa|}
-i\frac{|\boldsymbol\kappa|^2}{Re_\ell}.
\tag{29}
\]

ℓ 若取盘间距，Re_ℓ 是间距 Reynolds 数；若取外半径则是整体旋转 Reynolds 数。应在与现有程序连接时明确长度和频率尺度，不能把局部稳定性 R 与上述 Re_ℓ 自动等同。

建议先取给定 χ 的近刚体核心、弱单模态入射。可用 m_in=0 与动盘非轴对称浅粗糙 h_n(r)e^{inθ_d} 的组合：混合项满足 m_out=m_in±n、ω_d,out=ω_d,in。这样可以有针对性地隔离粗糙对角向波数的转换；它是受控机制算例，不代表已经测得真实静盘输入为轴对称。

若输入从上方核心向动盘传播，选择入射支应看群速度而不是只看 k_z 符号。无黏理想核心中

\[
c_{g,z}=W_c+s\,2\Omega_c\frac{k_h^2}{k^3}.
\tag{30}
\]

当 z 由动盘指向核心时，入射能量应满足 c_g,z<0。有限壁面反射后会同时出现两支；初步无反射受控模型与真实有限腔必须区分。

## 12. 与边界层感受性模型的连接

输入应满足核心动力学或有明确强迫源。到达动盘边界层之后，用实际 U,V,W 代入式(13)求受迫响应，不预设它已经是某个离散横流模态。

可分解为 q=q_inc+q_sc。若 L_e q_inc=0，而实际边界层算子是 L，则散射场满足 Lq_sc=−(L−L_e)q_inc。光滑壁面给 q_sc,w=−q_inc,w；远场指定入射分量，对散射分量采用适当的出射/衰减条件。加入粗糙后还要构造混合边界源及体源。不能沿用离散 LST 的“所有扰动在外缘归零”，同时又声称已从外缘输入非零来流扰动。

若任意指定 m、ω、空间尺度和极化而不满足自由色散关系，可以作为有明确体力的受迫模型，但不能称为自由传播的核心波。若只追踪静止粗糙散射，旋转系 ω 保持不变；频率变换仍依式(10)。

## 13. Julia 数值核对

运行 `julia --project=. scripts/verify_rotating_vortical_models.jl`。脚本只使用 LinearAlgebra、Test、Printf，验证解析入射场，不是腔体/粗糙感受性求解器。本次 Julia 1.12.6 运行 497/497 检查通过：

- 平面波动量残差最大约 1.24e−16，散度约 5.55e−17。
- 不同 Ω_d 下同一核心波的幅值方程有限差分误差约 2.90e−10；波矢输运误差约 5.14e−11。
- 一般仿射剪切下 d(k·a)/dt 约 4.22e−17。
- 把精确旋转波转成柱坐标后，直接有限差分代入式(15)的四个方程，联合残差约 4.30e−8。
- 另对 m=−3,0,2、两支频率的 Bessel 入射场进行独立柱坐标有限差分残差检查，最大约 5.45e−8。
- 包含 Ω_c=0、非零 Ω_c=Ω_d、零竖直波数、轴向波矢和频移符号检查。以上不是对全部式(13)系数、有限容器边界或粗糙散射解的独立数值验证。

结果写入 `benchmark_results/rotating_vortical_models/verification.txt`。

## 14. 原始文献与阅读范围

1. [Leib, Wundrow & Goldstein 1999](https://doi.org/10.1017/S0022112098003504), *Effect of free-stream turbulence and other vortical disturbances on a laminar boundary layer*, JFM 380, 169–203。核对期刊摘要及 NASA 作者稿关于 Fourier 叠加的说明；作者稿年份为 1997，正式论文为 1999。
2. [Guo & Kloker 2026](https://doi.org/10.1017/jfm.2026.11937), *Swept-wing boundary-layer receptivity to free-stream vortical disturbances*, JFM 1041, A27。核对开放全文 §2.3.5 的连续性、受控单 Fourier 模态及 Taylor 对流模型；其三种速度组合并非三个独立极化。
3. [Jacobs & Durbin 2001](https://doi.org/10.1017/S0022112000002469), *Simulations of bypass transition*, JFM 428, 185–212。核对期刊摘要及作者原文入口，使用连续谱构造湍流入流。
4. [Zaki & Durbin 2005](https://doi.org/10.1017/S0022112005003800), *Mode interaction and the bypass route to transition*, JFM 531, 85–111。核对作者原文关于连续谱的段落及期刊摘要。
5. [Bewley et al. 2007](https://doi.org/10.1063/1.2747679), *Inertial waves in rotating grid turbulence*, PoF 19, 071701；[作者全文](https://users.ictp.it/~krs/pdf/2007_002.pdf)。核对色散关系与旋转湍流波动实验。
6. [Thomas 2024](https://doi.org/10.1063/5.0180608), *Receptivity of the rotating disk boundary layer to traveling disturbances*。此前已核对作者全文，用于后续旋转盘伴随感受性连接；不能以其孤立单盘基态替代实际核心流。

式(7)至式(30)是基于明确假设的本次推导，并通过上述有限范围的代数/坐标变换检查。真实静盘湍流的谱、非线性相互作用和有限腔反射尚未由本模型确定。
