# 导师汇报前证据审计计划（2026-09-08）

## 加密 Rossby 扫描补充

- **目的**：检验原五个 Rossby 切片的拓扑结论是否因欠采样产生误导，并估计
  traditional 的 fold--tail 过渡和 consistent 的 smooth--fold、fold--tail
  过渡所在区间。
- **第一阶段统一扫描**：`Ro=-1.0,-0.9,-0.8,-0.7,-0.6,-0.5,-0.4,-0.3,
  -0.25,-0.2,-0.15,-0.1,-0.05`，两种热闭合均取 `N=120`，其余模型、映射、
  容差和分类判据与 L034 相同。
- **第二阶段边界复核**：根据第一阶段分类变化，只对相邻异类点的中间值补点，
  并至少用 `N=160` 复核最靠近每个过渡的代表点。
- **输出隔离**：所有新数据写入
  `work/results/presentation_evidence_audit_20260908/rossby_dense/`；不覆盖既有结果。
- **解释限制**：离散切片只能给出过渡区间，不能代替二维 fold continuation；
  `smooth` 仍表示在当前 `Tw<=1.99` 搜索窗内未见 fold/tail。

## 目标

在制作汇报幻灯片之前，独立复核学习手册中支撑主线的方程、数值数据和证据边界。
本次审计不使用历史未校正温度方程的结果，不把有限分辨率趋势升级为极限定理。

## 模型与参数

- acceleration-consistent BEK 相似基本流；恒物性、线性密度因子
  `chi=1-gamma*(T-1)`。
- 校正后的温度方程：`T'' + Pr*Ro*H*T' = 0`。
- `Pr=0.72`、`gamma=1`；有限问题采用 Pier 型真无穷映射
  `(a,b,c)=(2,0.6,0.5)`。
- fold 由完整增广系统 `R=0, J_U v=0, v^T v=1, Hinf+epsilon=0`
  定义；不得以 Newton 失败或曲线视觉转向替代。

## 新计算及输出目录

所有新增复核输出写入
`work/results/presentation_evidence_audit_20260908/`，不覆盖 `baselines/` 或既有
结果目录。

计划复核：

1. 用两个独立坐标算法重算 leading outer BVP 的
   `(Ro_t,Tw_t,p0)`，检查截断、步长和 shooting Jacobian。
2. 对 consistent residual 的解析 Jacobian、Ro 导数及方向 Hessian做独立中心差分
   检查。
3. 从保存的 corrected fold 状态重新校正代表性
   `epsilon=0.02,0.015,0.01`，检查增广残差、零模、切向和网格差异；重点区分
   `epsilon=0.01` 的趋势证据与收敛证据。
4. 复核 inner 相容曲线、解析斜率和 finite-fold 切向比值的算术一致性。
5. 审计手册中所有关键表格和结论是否指向 corrected 数据，必要时修改手册并重新
   编译、渲染检查。

## 预定结论等级

- **直接成立**：代数恒等式或在独立算法、网格/截断变化下稳定的数值量。
- **有限参数支持**：经完整增广 fold 求解确认，但尚无 `epsilon->0` 误差界。
- **候选机制**：与数据一致，仍缺 matched zero-mode/高阶可解性或可压缩交叉验证。
- **不得用于汇报结论**：依赖旧温度方程、未收敛最小 `epsilon`、或把相似子空间结果
  表述为完整三维稳定性的内容。

## 无 fold 约束的 generic-tail 极限图谱（2026-09-09）

- **问题**：检验普通 thermal-tail 解在 `epsilon=-Hinf->0+` 时是否保留非零
  inner core，并独立定位退化指标
  `D=chi(Tw)*s(Ro)^2-omega(Ro)^2=Ro*F''(0)` 的零点。
- **模型**：acceleration-consistent corrected-energy 基本流，`Pr=0.72`、`gamma=1`；
  只解 `R=0`，固定 `(Ro,Hinf)` 并释放 `Tw`，不调用 `J_U v=0` 或任何 fold
  求解器。
- **主采样**：`Ro=-.469,-.467,-.465,-.46,-.44,-.42,-.40,-.35,-.30`；
  `epsilon=.08,.05,.03,.02,.015,.01`。先用较低网格 pilot 验证分支路径，随后对
  小 epsilon 使用 `N=240,280,320,360` 及映射参数复核。
- **诊断**：保存 `Tw,D,F''(0)`，固定 `eta=1,2,5` 的场值，以及 `0<=eta<=8`
  的剖面。温度使用 `T-Tw`；core amplitude 先逐点外推
  `(F0,G0,H0)` 再积分，不能把有限-epsilon 振幅直接当作 `M0`。
- **证据限制**：多项式截距只有在窗口、阶数、网格和映射变化下稳定时才解释为
  `epsilon->0` 极限。固定-Ro generic tail 不满足 `Ro=Ro_t+O(epsilon)`，不得强迫
  使用 fold-limit 的参数展开。
- **输出隔离**：写入 `work/results/generic_tail_zero_core_20260909/`，不覆盖任何
  baseline 或现有审计结果。
