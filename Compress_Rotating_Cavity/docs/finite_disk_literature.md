# 有限旋转圆盘的研究模型：文献与选型

调研日期：2026-10-26。目标是有限旋转盘嵌在静止平面内、外部没有强加主流的层流基本流。本文区分直接几何先例、自由圆盘边缘研究和旋转腔研究；不将检索结果视为穷尽性综述，也不据此声称研究空白。

## 可以确定的结论

有限盘研究没有统一的外边界模型。“有限半径”“有限计算域”“有限厚度实体圆盘”“容器内部分旋转底壁”是不同问题。求整体基本流的直接先例采用完整轴对称 Navier–Stokes 方程，保留径向、周向、轴向三个速度分量及全局耦合；转捩研究进一步采用三维非定常方程。相似性解可作为初值、内区对照或扰动分析基态，其使用方式需要逐篇核对。

| 模型 | 物理/数学闭合 | 代表文献 | 与目标的关系 |
|---|---|---|---|
| 部分旋转底壁、平坦无切应力顶面 | 有限圆柱容器；中央盘旋转，外环和侧壁静止 | Piva & Meiburg 2005 | 最直接的同类底壁几何；顶部仍有限 |
| 部分旋转端壁、封闭容器 | 中央盘旋转，其他固壁静止 | Yu et al. 2007 | 可作为全局基本流求解器的验证例 |
| 自由薄盘边缘与尾迹 | 盘内无滑移，盘外中面由对称性/无剪切条件闭合 | Leslie 1972 | 盘外无静止壁面，与嵌入式盘不同 |
| 有限域旋转边界层及显式盘缘 | 海绵区、零厚度盘缘、有限厚度盘缘作对比 | Appelquist et al. 2015, 2016 | 重点是全局扰动与非线性转捩 |
| 容器内独立实体圆盘 | 圆盘上下流体、盘缘、侧壁、顶底边界 | Liu & Ladd 2022 | 适合实际旋转盘反应装置；不是嵌入式盘 |
| 有限双盘/转静腔 | 完整稳态方程；外缘封闭或开放 | Dijkstra & van Heijst 1983; Brady & Durlofsky 1987 | 支持研究边界如何选择内区相似性分支 |
| 有限域振幅模型 | 空间变系数 Ginzburg–Landau 方程 | Healey 2010 | 用于研究扰动反馈，不能直接给出基本流 |

## 最接近嵌入式盘的先例

**Piva & Meiburg (2005)** 使用轴对称完整方程的流函数—周向涡量—周向速度形式，在有限圆柱容器内求稳态；顶部取平坦无切应力近似。参数为 Re=Ωa²/ν、H/a、容器半径/a。底壁交界以 tanh 过渡平滑，并检查过渡宽度影响。外环存在向内流动，不能用单向外推求整个流场。[作者全文，§II–III](https://sites.me.ucsb.edu/~meiburg/pubs/Piva_and_Meiburg_2005.pdf)

该文的壁面周向速度正则化，用有量纲变量表示为

\[
v_w(r)=\frac{\Omega r}{2}\left[1-\tanh\left(\frac{r-a}{\varepsilon}\right)\right].
\]

这是一种边缘模型；复现时应记录 ε，研究真实环隙时应显式建立几何。论文中的 flat stress-free 顶面不等于无限远场。

**Yu, Lee, Zeng & Low (2007)** 研究同类部分旋转端壁，但采用封闭圆柱容器。模型为稳态轴对称完整 NS，有限体积、二阶中心差分与 SIMPLEC。旋转壁给 u_r=u_z=0、u_θ=Ωr；静止固壁给三分量零速度。文章考察容器尺寸、回流与涡破裂，不能将其上盖条件直接移植到无界流体。[DOI](https://doi.org/10.1063/1.2731420)，[作者上传全文](https://www.researchgate.net/publication/252992425_Characterization_of_flow_behavior_in_an_enclosed_cylinder_with_a_partially_rotating_end_wall)。

## 显式盘缘、有限域与实际容器

**Leslie (1972)** 从自由有限薄盘的边缘排出流出发，用匹配渐近展开研究旋转径向射流中的尾迹。自由盘外没有静止固壁，不能用其盘外条件替代目标模型的无滑移外环。[原文入口](https://doi.org/10.1017/S1446788700013707)

**Appelquist et al. (2015)** 采用线性化不可压 NS，以径向阻尼区模拟下游湍流的影响，研究入射扰动和全局反馈。该处理定义的是扰动问题的下游闭合，不等于解析了实际盘缘。[期刊原文](https://www.cambridge.org/core/journals/journal-of-fluid-mechanics/article/abs/global-linear-instability-of-the-rotatingdisk-flow-investigated-through-simulations/DB6C29B986352A425A6B980C6C71F12F)

**Appelquist et al. (2016)** 在非线性 NS 模拟中比较海绵区、零厚度盘缘、有限厚度盘缘。结果区分线性全局稳定与非线性全局失稳，说明下游模型和外部湍流反馈必须随结论一起说明。它不是嵌入静止平面几何的直接验证。[期刊摘要](https://www.cambridge.org/core/journals/journal-of-fluid-mechanics/article/abs/on-the-global-nonlinear-instability-of-the-rotatingdisk-flow-over-a-finite-domain/1B578FF3C8880BC601EAB82D34BCC6EE)，[作者全文](https://bura.brunel.ac.uk/bitstream/2438/16980/1/FullText.pdf)

**Liu & Ladd (2022)** 建立有限实体盘及较大容器的三维非定常 NS 模型，比较自由表面/固壁顶部等配置。所研究装置可在 Re_a<1000 出现非定常状态，不能将其阈值套用于无限盘或嵌入式盘。其 OpenFOAM 算例有补充材料。[开放全文与补充材料](https://www.cambridge.org/core/journals/journal-of-fluid-mechanics/article/onset-of-instabilities-in-rotating-flows-by-direct-numerical-simulation/B0337EA81F64986E9A00FCC37AC9BD76)

**Healey (2010)** 的空间变系数复 Ginzburg–Landau 模型研究有限径向边缘如何影响全局模态。它适合机制分析，不能替代速度—压力基本流的边值问题。[期刊原文](https://www.cambridge.org/core/journals/journal-of-fluid-mechanics/article/model-for-unstable-global-modes-in-the-rotatingdisk-boundary-layer/9887B9BF7A344C293B8A23F17D78C8A9)

## 与相似性选择直接相关的双盘文献

**Dijkstra & van Heijst (1983)** 对圆柱包围的有限双盘流求完整稳态 NS，不施加相似性假设，并与实验比较。对其考察参数，盘内部分区域接近相似性形式；不应把该结果推广为所有几何均如此。[期刊原文](https://www.cambridge.org/core/journals/journal-of-fluid-mechanics/article/abs/flow-between-two-finite-rotating-disks-enclosed-by-a-cylinder/D37C9096B1B0D96DD89B5EE1DF0DD390)

**Brady & Durlofsky (1987)** 用渐近—数值方法比较有限双盘的封闭外缘和开放外缘。外缘条件可改变有限盘流接近哪个相似性分支；某些参数下近轴区也不呈相似性。该文是研究“相似性解是否代表有限装置基本流”的重要前例，但其对象是双盘。[期刊原文](https://www.cambridge.org/core/journals/journal-of-fluid-mechanics/article/abs/on-rotating-disk-flow/FA1EBEEA80AA8B483DBAD3152D1E7392)

还有一类 flush-mounted rotating-disc 文献研究通道/边界层减阻，例如 Ricco & Hahn (2013)。其中存在强加主流，常含圆盘阵列，不能直接用作静止外部流体中单盘的轴对称基本流依据。[原文](https://www.cambridge.org/core/journals/journal-of-fluid-mechanics/article/abs/turbulent-drag-reduction-through-rotating-discs/5D13FE3B1BE4E9E26BE9D92B34A3101F)

## 对本项目的建议（基于上述文献的判断）

先用 Piva–Meiburg 的有限容器几何验证 Julia 轴对称全局 NS 求解器，再研究目标的无限静止平面/无限上方流体。后者需更换为远场开放处理，独立扩大径向与轴向截断距离；不能只把容器半径加大而保留有限高度顶盖，就宣称达到无界流体极限。

层流基本流阶段仅限制 ∂t=∂θ=0，保留三个速度分量、两维空间变化及压力耦合。允许回流存在；相似性剖面只作初始猜测或对照。稳定性研究另对该全局基本流线性化，不能把径向吸引性等同于时间/非轴对称稳定性。

有价值的研究量是归一化剖面误差 E(r)、壁面剪切、卷吸和压力梯度随 Re_a、盘缘过渡/环隙尺度变化的规律，以及满足指定误差阈值的相似性有效半径区间。容器 H/a、R_c/a 是物理参数；无界流体的 L/a、H/a 是数值截断参数，应在报告中区分。

只实现有限盘基本流本身不构成新颖性证明；需要与已有的部分旋转端壁和有限盘相似性选择研究比较。本轮没有确认一篇与“无限静止平面内单个有限旋转盘、无限上方静止流体”完全同构且足以作为唯一标准基准的论文；这不等于该问题没有已有研究。

## 文献表与核验范围

1. Piva, M. & Meiburg, E. (2005). *Steady axisymmetric flow in an open cylindrical container with a partially rotating bottom wall*. Physics of Fluids **17**, 063603. DOI: [10.1063/1.1932664](https://doi.org/10.1063/1.1932664)。已核对作者全文的方法、边界条件与过渡层。
2. Yu, P., Lee, T. S., Zeng, Y. & Low, H. T. (2007). *Characterization of flow behavior in an enclosed cylinder with a partially rotating end wall*. Physics of Fluids **19**, 057104. DOI: [10.1063/1.2731420](https://doi.org/10.1063/1.2731420)。已核对作者上传全文的方法段；出版社链接未能读取。
3. Leslie, L. M. (1972). *The wake of a finite rotating disc*. DOI: [10.1017/S1446788700013707](https://doi.org/10.1017/S1446788700013707)。核对作者上传原文的引言、模型与方法，未复核渐近系数。
4. Appelquist, E., Schlatter, P., Alfredsson, P. H. & Lingwood, R. J. (2015). *Global linear instability of the rotating-disk flow investigated through simulations*. JFM **765**, 612–631. DOI: [10.1017/jfm.2015.2](https://doi.org/10.1017/jfm.2015.2)。核对期刊摘要。
5. 同作者 (2016). *On the global nonlinear instability of the rotating-disk flow over a finite domain*. JFM **803**, 332–355. DOI: [10.1017/jfm.2016.506](https://doi.org/10.1017/jfm.2016.506)。核对期刊摘要与作者全文的模型说明；未逐一复核外边界实现。
6. Liu, Z. & Ladd, A. J. C. (2022). *Onset of instabilities in rotating flows by direct numerical simulation*. JFM **945**, A31. DOI: [10.1017/jfm.2022.566](https://doi.org/10.1017/jfm.2022.566)。核对开放全文几何和数值方法。
7. Healey, J. J. (2010). *Model for unstable global modes in the rotating-disk boundary layer*. JFM **663**, 148–159. DOI: [10.1017/S0022112010003836](https://doi.org/10.1017/S0022112010003836)。核对期刊摘要。
8. Dijkstra, D. & van Heijst, G. J. F. (1983). *The flow between two finite rotating disks enclosed by a cylinder*. JFM **128**, 123–154. DOI: [10.1017/S0022112083000415](https://doi.org/10.1017/S0022112083000415)。核对期刊摘要和原文首页；年份按卷期，不按后来的电子上线日期。
9. Brady, J. F. & Durlofsky, L. (1987). *On rotating disk flow*. JFM **175**, 363–394. DOI: [10.1017/S0022112087000430](https://doi.org/10.1017/S0022112087000430)。核对期刊摘要。
10. Ricco, P. & Hahn, S. (2013). *Turbulent drag reduction through rotating discs*. JFM，原文链接见上。核对期刊摘要，仅用于区分有主流的嵌入式盘模型。

另检索了近年的有限半径、盘缘和部分旋转端壁工作。本笔记不把纳米流体/传热论文中单独使用的相似性 ODE、或标题中的 finite element 一词，视作已解析有限半径和真实盘缘的证据。核心模型判断依赖上列原始研究，不以年代新旧替代几何与方程核验。
