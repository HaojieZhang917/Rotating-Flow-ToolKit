# 总体研究进展日志

本文件是项目唯一的总体计算结论日志，用于在论文定稿时重建研究主线、
证据链和逻辑演化。详细模型定义、参数表、数据文件和专题报告仍保存在各自
文件中；这里记录每项计算回答了什么问题、得到什么结论、结论适用到哪里，
以及它如何改变下一步研究方向。

## 维护规则

1. 每次完成新的数值计算、收敛测试、分支延拓、模型交叉验证或机理诊断后，
   必须在本文件末尾追加一条记录。
2. 记录必须在结果被用于论文、汇报或下游计算之前完成。
3. 已有记录不删除、不静默改写。若新结果修正旧结论，新增一条记录并在两条
   记录中注明 `supersedes` / `superseded by`。
4. 明确区分：
   - 完全自洽的模型比较；
   - 冻结基本流或单项开关诊断；
   - 数值可解性范围；
   - Boussinesq 近似的物理精度范围。
5. “未发现”必须注明搜索范围，不能写成全局不存在；由等温解延拓得到的分支
   不能自动排除孤立解支。
6. 鞍结必须至少有分支回折和固定控制参数 Jacobian 奇异性证据；cusp 必须有
   两参数折叠延拓或等价的余维二条件。单个观测量的极值不称为鞍结。
7. 每条记录至少包含：日期、模型、参数和数值方法、问题、主要结果、验证、
   结论边界、主线意义、输出位置和下一步。

## 当前研究主线快照（2026-09-03）

1. 等温 BEK 方程和真无穷映射已经在 `-1 <= Ro <= 1` 上验证，可作为热基本流
   的统一数值底座。
2. 盘随动系传统固定离心 Boussinesq 闭合在 von Karman 端点产生两个非退化
   鞍结；有限折叠对只存在于非常窄的近 von Karman 区域。
3. 离开该区域后，负 `Ro` 分支的主要终止机制变成 `Hinf -> 0-` 导致的无限长
   热尾，而不是有限参数鞍结。
4. 传统模型中独立 `(T-1)` 项的正确盘系系数是
   `gamma*Co^2/(4Ro)`，其中 `gamma=beta_inf*T_inf`。当前空气/理想气体计算取
   `gamma=1`。该孤立项依赖参考系，因此全 BEK 传统结果必须标为“盘系固定
   离心闭合”，不能解释成参考系不变的物理预测。
5. 加速度一致的 von Karman 模型在等温连接加热分支上连续延拓到 `Tw=1.99`，
   没有固定-`Tw` 鞍结。`Hinf(Tw)` 在 `Tw≈1.647` 的极小值只是观测量极值。
6. `Ro=-1` 的传统鞍结具有明确的闭合选择性，但代表性一致模型扫描表明奇点
   并未在全家族中整体消失：一致模型在 `Ro=-0.5` 出现另一有限折叠，在
   `Ro=-0.25,-0.1` 保留热尾边界。当前更准确的主线是“闭合方式重新分配并
   转换 fold/tail 拓扑”，而不是“一致模型消灭全部鞍结”。
7. 这一拓扑重排仍需同一 `Ro` 网格上的两模型直接比较、两参数折叠延拓和匹配
   的低马赫数/可压缩计算确认其物理意义。

## 进展记录

### 2026-09-01 — L001：等温 BEK 真无穷映射基准

- **状态**：已验证的数值底座。
- **模型**：等温、不可压缩 Lingwood BEK 相似性方程。
- **参数与方法**：`Ro=-1,-0.75,...,1`；Pier 半无限映射；生产设置
  `N=120,(a,b,c)=(2,0.6,0.5)`。
- **问题**：现有 von Karman ODE 是否能一致扩展到完整 BEK 家族，并恢复三个
  经典端点。
- **主要结果**：所有采样 `Ro` 均收敛；`Ro=-1` 恢复 von Karman，`Ro=0`
  对应 Ekman 极限，`Ro=1` 对应 Bodewadt。`Ro=-1` 残差约 `4.1e-11`。
- **验证**：与保存的 von Karman 解及 Lingwood 壁面导数交叉比较；映射参数和
  阶数测试支持采用 `N=120,(2,0.6,0.5)`。
- **结论边界**：只验证等温基本流，不包含温度、密度闭合或分岔结论。
- **主线意义**：建立所有后续传统/一致模型必须共享的 ODE 与无穷远离散底座。
- **输出**：`work/results/gate1_isothermal/`；
  `docs/GATE1_ISOTHERMAL_BEK.md`。
- **下一步**：从同一有量纲母方程导出热闭合，不能经验添加通用温度源项。

### 2026-09-01 — L002：传统 von Karman 第一和第二鞍结复现

- **状态**：已验证的传统模型基线。
- **模型**：盘随动系传统固定离心 Boussinesq，`Ro=-1`，`Pr=0.72`，
  `gamma=1`。
- **参数与方法**：真无穷映射，`N=120,(a,b,c)=(2,0.6,0.5)`；以
  `Hinf` 和伪弧长延拓；固定-`Tw` Jacobian 检查。
- **问题**：旧有限域计算得到的折叠是否为远场截断或求解器伪影，是否存在
  第二个折叠。
- **主要结果**：第一折叠
  `(Hinf,Tw)≈(-0.53276,1.048021731)`；第二折叠
  `(Hinf,Tw)≈(-0.11332,1.030145)`。旧 `zmax=20` 的第二折叠位置具有明显
  有限域偏差。
- **验证**：真无穷映射、Jacobian 简单零方向、平方根律、非退化系数和
  相似性子空间时间谱均支持非退化 saddle-node 判定。
- **结论边界**：这是特定传统盘系闭合的数学分支拓扑，不代表所有 Boussinesq
  或可压缩模型都在该温度失去解。
- **主线意义**：提供需要由一致模型和可压缩模型解释的核心异常现象。
- **输出**：`work/results/gate1_ro_minus1_hinf_tw.png`；
  `baselines/saddle_node/`；`docs/CURRENT_STATUS.md`。
- **下一步**：追踪折叠随 `Ro` 的运动并区分有限折叠与远场热尾失效。

### 2026-09-01 至 2026-09-02 — L003：传统负 `Ro` 折叠对合并

- **状态**：折叠消失区间已解析；精确 cusp 坐标仍为暂定。
- **模型**：盘系传统固定离心 Boussinesq BEK。
- **参数与方法**：近 `Ro=-1` 加密扫描，`N=120`，真无穷映射。
- **问题**：为什么只有 von Karman 附近出现有限鞍结，折叠对在哪里消失。
- **主要结果**：两个折叠随 `Ro` 增大迅速靠近；`Ro=-0.9938` 仍检测到极近的
  折叠对，`Ro=-0.9936` 已未检测到两折叠。当前给出的 cusp 区间约为
  `-0.9938 < Ro_c < -0.9936`，但精确余维二坐标尚未由 bordered cusp 方程
  求出。
- **验证**：局部分支加密、固定控制参数条件数和折叠对连续靠近趋势一致。
- **结论边界**：可以说“折叠对在狭窄区间内消失”，暂不能把插值点称为精确
  cusp 坐标。
- **主线意义**：传统有限折叠不是全 BEK 普遍现象，而是近 von Karman 的
  狭窄拓扑区。
- **输出**：`work/results/traditional_bek_cusp/`；
  `docs/REFINED_FOLD_MERGER.md`。
- **下一步**：在统一闭合界面上做真正的两参数 fold/cusp continuation。

### 2026-09-02 — L004：传统全 BEK 加热分支适用性/可解性图

- **状态**：传统盘系闭合的生产基线。
- **模型**：盘系传统固定离心 Boussinesq，`-1 <= Ro <= 1`。
- **参数与方法**：正则热强度
  `B=Lambda_cf*(Tw-1)`；`N=120,(a,b,c)=(2,0.6,0.5)`；真无穷映射。
- **问题**：各 `Ro` 的等温连接加热分支在哪里终止，终止原因是什么。
- **主要结果**：
  - 只有近 von Karman 的窄区由有限 saddle-node 控制；
  - `-0.9936 <= Ro < 0` 的终止机制是 `Hinf -> 0-` 和无限长热尾；
  - `Ro->0-` 时传统模型允许的有限壁温区间收缩到等温点；
  - `Ro=0` 是差速尺度下的非一致系数极限，不分配非等温相似解；
  - `0<Ro<1` 在形式延拓至 `Tw=2` 的范围内没有检测到折叠；
  - `Ro=1` 在盘系传统闭合中温度对动量被动。
- **验证**：代表性点的 `N=80,120,160`、不同映射参数、固定-`B` Jacobian
  以及远场热衰减长度检查。
- **结论边界**：这是传统方程的数学可解性图，不是相对于可压缩模型的
  Boussinesq 精度范围；`Tw=2` 不是物理适用上限。
- **主线意义**：首次把“有限鞍结失效”和“无限热尾失效”分成两种机制，形成
  后续基本流数据层。
- **输出**：
  `work/results/full_bek_traditional_availability/production_v2_N120_a2.0_b0.6_c0.5/`；
  `docs/FULL_BEK_TRADITIONAL_RESULTS.md`。
- **下一步**：对同一 `Ro` 网格计算加速度一致端点，比较拓扑而不只比较剖面误差。

### 2026-09-02 — L005：`Ro=-0.5` 传统模型失效过程剖面

- **状态**：已完成的机理可视化。
- **模型**：盘系传统固定离心 Boussinesq，`Ro=-0.5`，`Pr=0.72`。
- **参数与方法**：`Tw=1` 至 `1.6`；生产真无穷映射。
- **问题**：传统模型接近适用性边界时，基本流和热尾如何演化。
- **主要结果**：加热削弱径向环流和轴向抽吸；温度尾随 `Hinf` 绝对值减小而
  显著伸长。该分支的主要问题是热尾逐步失去指数局域性，而不是在所考察温度
  内先遇到有限折叠。
- **验证**：剖面、远场衰减率和分支延拓结果相互一致。
- **结论边界**：描述的是传统盘系闭合的相似性解演化，不能直接代表真实空气
  在大温差下的定量行为。
- **主线意义**：给“无限热尾”提供了直接剖面证据，并解释其与 von Karman
  有限折叠的不同。
- **输出**：`docs/RO_M05_PROFILE_SWEEP.md` 及其记录的结果目录。
- **下一步**：在一致模型下检查抽吸是否仍趋零，以判断热尾是否同样存在。

### 2026-09-03 — L006：传统 `(T-1)` 项系数重新推导

- **状态**：已完成的模型定义修正；不是一次参数扫描。
- **母方程**：任意匀速旋转参考系中的有量纲 Navier--Stokes 方程，密度线性化
  `rho'/rho_inf=-gamma*(T-1)`。
- **问题**：BEK 径向方程中独立 `(T-1)` 项的系数和符号到底是什么。
- **主要结果**：在 Lingwood 盘随动系和差速尺度下
  `Lambda_D=gamma*Omega_D^2/(Omega*DeltaOmega)
  =gamma*Co^2/(4Ro)`。当前代码取 `gamma=1`，只对理想气体线性化成立；
  `T_inf=1` 本身不能推出 `gamma=1`。`Ro=-1` 时系数为 `-gamma`，准确恢复
  旧 von Karman 方程。
- **验证**：从有量纲旋转系方程、约化压力、Lingwood 相似变换逐项推导，并用
  Julia 验证 `Ro=-1,-0.5,0.5,1` 的端点恒等式。
- **结论边界**：该系数只属于“盘系、选择性保留参考离心密度扰动”的传统闭合。
  换到远场流体系或惯性系会得到不同的孤立系数；加速度一致模型中温度扰动
  乘完整惯性加速度，不存在唯一独立常系数源项。
- **主线意义**：把全 BEK 传统结果重新限定为参考系依赖的诊断基线，也解释了
  `Ro=0` 奇异和 `Ro=1` 被动温度为何不能直接当作普遍物理结论。
- **输出**：`work/src/BEKThermal.jl`；
  `docs/FULL_BEK_TRADITIONAL_AVAILABILITY.md`；本日志的模型限定。
- **下一步**：所有一致模型计算显式保留 `gamma`，并从完整惯性加速度构造残差。

### 2026-09-03 — L007：加速度一致 von Karman 加热分支

- **状态**：已验证的 `Ro=-1` 一致模型结果。
- **模型**：稳态加速度一致 Boussinesq；`Ro=-1`，`Pr=0.72`，`gamma=1`。
- **参数与方法**：真无穷映射，生产设置
  `N=120,(a,b,c)=(2,0.6,0.5)`；固定 `Tw` 从 `1` 延拓至 `1.99`；完整状态
  伪弧长和固定-`Tw` Jacobian 奇异值检查。
- **问题**：传统 von Karman 鞍结在密度与完整惯性加速度一致耦合后是否仍存在，
  分支上的 `Hinf` 回升是否为 saddle-node。
- **主要结果**：等温连接分支连续到 `Tw=1.99`，没有固定-`Tw` 鞍结。
  `Hinf(Tw)` 在
  `(Tw,Hinf)≈(1.64704657,-0.9120805930)` 有一个光滑极小值，但伪弧长
  `Tw` 始终增加，最小前向增量为 `6.84098e-3`，所以该点不是鞍结。
  `Hinf` 保持约 `-0.9`，没有无限热尾。
- **验证**：
  - 全部分支残差低于 `1e-10`；
  - 固定-`Tw` 行缩放奇异值比从 `7.70e-7` 平滑变化到 `5.77e-7`，未出现
    传统折叠的机器精度坍缩；
  - `N=80,120,160` 在 `Tw=1.2,1.8` 的 `Hinf` 和 `F'(0)` 收敛；
  - `Tw=1.2,1.5` 与保存的有限域 Lopez 结果在显示精度内一致。
- **结论边界**：只排除 `1<=Tw<=1.99` 内从等温解连接的加热分支上的折叠；
  未全局排除孤立解支。`Tw>1.2` 只作闭合模型形式诊断，不能解释成物理精度。
  固定 `Ro=-1,gamma=1` 是一参数切片，不能据此全局排除两参数 cusp。
- **主线意义**：最直接地支持“传统鞍结具有闭合选择性”。一致模型的光滑
  `Hinf` 极值也给出了一个重要警告：不能把任意分支投影极值误判为鞍结。
- **输出**：`work/results/vonkarman_consistent_baseflow/`；
  `work/src/BEKConsistent.jl`；
  `work/scripts/vonkarman_consistent_baseflow.jl`；
  `docs/VON_KARMAN_CONSISTENT_BASEFLOW.md`。
- **下一步**：在非零 `Ro` 上验证一致残差，然后计算完整一致 BEK 家族的
  `(Ro,Tw)` 解域和两参数 fold/cusp 图；同时在小温差区加入匹配可压缩基准。

### 2026-09-03 — L008：一致模型代表性 `Ro` 拓扑扫描

- **状态**：已验证的代表切片结果；不是完整二维图。
- **模型**：稳态加速度一致 Boussinesq，`Pr=0.72`，`gamma=1`。
- **参数与方法**：`Ro=-1,-0.75,-0.5,-0.25,-0.1,0.25,0.5,1`；
  `N=120,(a,b,c)=(2,0.6,0.5)` 真无穷映射；固定-`Tw` 自适应延拓到
  `Tw=1.99` 或分支边界；Jacobian、热尾长度和必要时的完整状态伪弧长检查。
- **问题**：传统模型的 fold/tail 拓扑在一致闭合下是整体消失、局部消失，
  还是转换为另一类边界。
- **主要结果**：
  - `Ro=-1,-0.75`：等温连接分支平滑到 `Tw=1.99`；
  - `Ro=-0.5`：存在有限折叠
    `(Tw_c,Hinf_c)≈(1.721704197,-0.29758392)`；
  - `Ro=-0.25,-0.1`：分别在 `Tw≈1.40512,1.18225` 进入
    `Hinf->0-` 无限热尾；
  - `Ro=0.25,0.5,1`：平滑到 `Tw=1.99`，且 `Hinf` 远离零；
  - 一致模型的 `Ro=1` 温度不再被动。
- **验证**：`Ro=-0.5` 伪弧长跨过 `dTw/ds=0` 并进入 `dTw<0` 返回支；
  `N=80,120,160` 给出
  `Tw_c=1.7217041979,1.7217041973,1.7217041945`。两个热尾切片在三种
  阶数下均保持 `Hinf->0-`，端点外推差异约 `2.5e-4` 和 `4.5e-5`。
- **结论边界**：只覆盖八个固定 `Ro` 的等温连接加热支。平滑到 `Tw=1.99`
  不代表大温差物理准确，也未排除孤立解支。尾端点尚未做专门渐近外推，不能
  按折叠坐标的精度引用。
- **主线意义**：否定了“一致模型在整个负 `Ro` 区域全部光滑”的简单假设。
  更强也更准确的结论是：密度-加速度闭合会把分支边界在 `Ro` 空间中迁移，
  并在有限折叠与无限热尾之间转换。
- **输出**：`work/results/consistent_representative_ro_sweep/`；
  `work/scripts/consistent_representative_ro_sweep.jl`；
  `docs/CONSISTENT_REPRESENTATIVE_RO_SWEEP.md`。
- **下一步**：在同一八个 `Ro` 上生成传统/一致统一拓扑表，量化“fold 消失、
  fold 产生、tail 保留、Bodewadt 激活”四类差异；随后对一致模型折叠区做
  两参数 continuation。

### 2026-09-03 — L009：传统/一致模型统一八点拓扑表

- **状态**：已完成并按共同判据复算。
- **模型**：盘系传统固定离心闭合与稳态加速度一致闭合；均为完全自洽基本流，
  不是冻结场或项开关诊断。
- **参数与方法**：八个共同 `Ro=-1,-0.75,-0.5,-0.25,-0.1,0.25,0.5,1`；
  `Pr=0.72,gamma=1`；`N=120,(a,b,c)=(2,0.6,0.5)` 真无穷映射；固定
  `Tw`、Jacobian 奇异值、必要时跨折叠延拓及 `Hinf->0-` 外推。
- **问题**：闭合改变究竟是统一移除奇异性，还是重新组织分支终止拓扑。
- **主要结果**：传统分类为
  `fold,tail,tail,tail,tail,smooth,smooth,passive`；一致分类为
  `smooth,smooth,fold,tail,tail,smooth,smooth,smooth`。传统 `Ro=-1` 折叠
  `Tw_c=1.0480217310`；一致 `Ro=-0.5` 折叠
  `Tw_c=1.7217041955`。闭合既能移除折叠，也能在其他 `Ro` 产生/迁移折叠，
  并能把传统被动端点变成主动耦合端点。
- **验证**：两个折叠均有固定-`Tw` Jacobian 奇异及跨折叠 `dTw/ds` 换号；
  传统温度与保存基线一致，一致温度与独立 `N=80,120,160` 结果一致。
- **结论边界**：tail 的 `Tw_c` 是 `Hinf=0` 外推；其 Jacobian 数值取最后
  有限热尾状态，不能冒充折叠奇异值。`Tw=1.99` 只表示形式延拓范围。
- **主线意义**：形成真正的 closure-dependent topology table，确立论文主张应为
  “全局终止拓扑重组”，而不是“一致闭合普遍消除折叠”。
- **输出**：`docs/TWO_MODEL_TOPOLOGY_COMPARISON.md`；
  `work/results/two_model_topology_comparison/`；
  `work/scripts/compare_traditional_consistent_topology.jl`。
- **下一步**：只延拓一致模型 `Ro=-0.5` fold locus，并与热尾边界对照。

### 2026-09-03 — L010：一致模型二维 fold--tail 拓扑

- **状态**：有限核心 fold 臂已验证；fold--tail 连接为强数值证据，精确奇异极限
  坐标仍暂定；有限 cusp 结论被收敛检查否定。
- **模型**：稳态加速度一致 BEK，完全自洽基本流。
- **参数与方法**：`Pr=0.72,gamma=1`；生产
  `N=120,(a,b,c)=(2,0.6,0.5)`；增广方程
  `R=0,Jv=0,v'v=1` 的完整状态伪弧长二维 fold continuation；固定
  `Hinf=-0.02,-0.01,-0.005,-0.0025` 热尾代理；`N=80,160` 拓扑收敛检查。
- **问题**：一致模型 fold 曲线由有限 cusp 结束，还是转入 `Hinf=0` 远场失局域
  边界。
- **主要结果**：种子修正为
  `(Ro,Tw,Hinf)=(-0.5,1.721704195473,-0.297584911495)`。fold 主臂在
  `Hinf=-0.1,-0.05,-0.02` 上的 `N=120/160` 坐标收敛。固定 `N=120`
  出现两个 cusp-like 非退化系数过零，但第一个表观转向的 `Hinf` 随
  `N=80,120,160` 从 `-1.73e-2` 移到 `-6.91e-3,-3.82e-3`，第二组漂移
  更大；没有收敛的有限 cusp。生产计算最靠近热尾处
  `Hinf=-1.475e-3`，fold 与外推 tail 的 `Tw` 差 `-7.93e-5`。
- **验证**：解析 `d(Jv)/du` 与 `R_Ro` 对有限差分的相对误差分别
  `7.5e-8,7.5e-10`；种子 fold 温度在 `N=80,120,160` 一致；有限负
  `Hinf` 层级的 fold 坐标收敛；表观 cusp 不随分辨率收敛而向热尾退缩。
- **结论边界**：当前支持“fold 臂进入远场失局域层”，不提供精确
  `Hinf=0` 正则解或严格连接证明。`N=160` 二次外推
  `(Ro,Tw)≈(-0.3072,1.4783)` 仅为定位提示，不能作为最终临界坐标。
- **主线意义**：最有价值的结果不是找到有限 cusp，而是显示有限核心折叠在
  接近热尾时转换为远场失局域控制；密度--加速度闭合确实重组全局终止机制。
- **输出**：`docs/CONSISTENT_FOLD_TAIL_CONTINUATION.md`；
  `work/results/consistent_fold_tail_continuation/`；
  `work/scripts/continue_consistent_fold_tail.jl`；
  `work/scripts/assess_fold_resolution_convergence.jl`。
- **下一步**：暂停新增物理，若需提高连接坐标精度，只做热尾适配坐标/专门奇异
  极限校正；在此之前不做稳定性、低 Mach 或新物性扫描。

### 2026-09-04 — L011：一致模型固定 `epsilon=-Hinf` 的 fold 数据链

- **状态**：五点数据链已完成；有限极限趋势已确认，未做渐近推导。
- **模型**：稳态加速度一致 BEK，完全自洽基本流。
- **参数与方法**：固定
  `epsilon=0.1,0.05,0.02,0.01,0.005`；增广方程
  `R=0,Jv=0,v'v=1,Hinf+epsilon=0`，释放 `Ro`；生产
  `N=120,(a,b,c)=(2,0.6,0.5),Pr=0.72,gamma=1`；长热尾点用
  `N=80,160,200,240,280,320` 做必要分辨率检查。
- **问题**：不用寻找 cusp，直接确认 fold 分支在 `epsilon->0+` 时如何接近热尾
  边界，以及 `Ro_f,Tw_f` 是否趋向有限值。
- **主要结果**：高阶 `N=240--320` 共识表为
  `(epsilon,Ro_f,Tw_f)`：
  `(0.1,-0.368218374,1.552811091)`、
  `(0.05,-0.342987270,1.521896650)`、
  `(0.02,-0.329321644,1.505245949)`、
  `(0.01,-0.324984263,1.499970045)`、
  `(0.005,-0.322827361,1.497347532)`。最后一点高阶半范围约
  `7.2e-5` (`Ro`) 和 `8.7e-5` (`Tw`)。序列随 `epsilon` 减小持续变平，支持
  有限极限。后三点线性截距诊断为
  `(Ro_t,Tw_t)≈(-0.32065867,1.49470958)`。
- **验证**：每点直接满足固定-`epsilon` fold 增广方程；行缩放固定-`Tw`
  Jacobian 最小奇异值为 `O(1e-15)`；生产残差低于 `4e-10`，高阶最大残差约
  `2e-9`；`N=320` 温度剖面单调且热尾随 `epsilon` 减小系统伸长。
- **结论边界**：`N=120` 只在 `epsilon>=0.02` 收敛；其 `0.01,0.005` 点不能
  用于极限判断。给出的截距只是有限性检查，不是渐近标度、严格连接坐标或
  matched-asymptotic 结果。
- **主线意义**：获得后续 fold-to-tail scaling 的第一条可靠控制参数数据链，
  并直接支持 fold 臂在接近热尾时趋向有限 `(Ro,Tw)`，而不是有限 cusp。
- **输出**：`docs/CONSISTENT_EPSILON_FOLD_CHAIN.md`；
  `work/results/consistent_epsilon_fold_chain/`；
  `work/scripts/consistent_epsilon_fold_chain.jl`；
  `work/scripts/summarize_epsilon_fold_convergence.jl`。
- **下一步**：按当前要求停止；不扩展 `Ro`、不做渐近推导。若后续启动
  matched asymptotics，应以本表及其分辨率误差为输入。

### 2026-09-04 — L012：一致模型 fold--tail 的 `epsilon` 标度拟合

- **状态**：已完成；只使用 L011 的五个既有数据点，没有新增基本流计算。
- **模型**：稳态加速度一致 BEK 的 fold 链；数值后处理而非新的闭合或冻结场诊断。
- **参数与方法**：输入
  `epsilon=0.1,0.05,0.02,0.01,0.005` 的 `N=240--320` 共识表；最后三个
  等比点的无极限差分指数；五点/去掉最大点的线性与二次最小二乘；输入半范围
  角点扰动；残差阶检查。全部用 Julia。
- **问题**：确定
  `Ro_f-Ro_t~A_R epsilon^p`、`Tw_f-Tw_t~A_T epsilon^q` 的指数和一阶系数，
  为后续 distinguished scaling 提供数值依据。
- **主要结果**：最后三个等比点给出
  `p_eff=1.007864,q_eff=1.008468`；输入角点范围分别为
  `[0.9551,1.0625]`、`[0.9557,1.0631]`，支持 `p=q=1`。去掉
  `epsilon=0.1` 的二次拟合为
  `Ro_t=-0.320719730,A_R=-0.420342495,B_R=-0.500001415` 和
  `Tw_t=1.494786622,A_T=0.510660114,B_T=0.630605042`。保留全部五点时
  `A_R=-0.415086,A_T=0.503275`，线性系数窗口变化约 `1.3--1.5%`。
- **验证**：四点输入角点扰动给出
  `A_R in [-0.42766,-0.41303]`、`A_T in [0.50178,0.51955]`；去掉最大
  `epsilon` 后截距变化仅 `O(5e-5)`。二次模型相对线性模型把 RMS 残差降低
  约 `13` 倍；全五点比较降低约 `28--30` 倍。归一化残差支持领先
  `O(epsilon)` 和次阶近似 `O(epsilon^2)`。
- **结论边界**：这是五点数值拟合，不是 matched asymptotics 或解析证明；
  最小两点误差除以 `epsilon^2` 后被放大，因此不宣称高精度二阶系数。
- **主线意义**：后续热尾匹配应采用
  `Ro-Ro_t=O(-Hinf), Tw-Tw_t=O(-Hinf)`，数值一阶系数约
  `A_R=-0.42,A_T=0.51`。
- **输出**：`docs/CONSISTENT_FOLD_TAIL_SCALING_FIT.md`；
  `work/results/consistent_fold_tail_scaling/`；
  `work/scripts/fit_consistent_fold_tail_scaling.jl`。
- **下一步**：按当前要求停止；不向更小 `epsilon` 计算，也不在本步骤开展
  matched asymptotics。

### 2026-09-04 — L013：一致模型 `epsilon->0+` 的 outer thermal-tail 方程

- **状态**：解析推导完成；没有新数值计算，也没有 inner--outer matching。
- **模型**：稳态加速度一致 BEK；从 `work/src/BEKConsistent.jl` 的完整自洽残差
  出发，不是冻结基本流诊断。
- **参数与方法**：`epsilon=-Hinf`，`Y=epsilon eta`；外层展开
  `H=epsilon h,F=epsilon^2 f,G=g,T=Theta`，`Ro=Ro_t+O(epsilon)`。
- **问题**：证明热尾长度、检查外层速度能否冻结为远场常数，并给出明确的 leading
  outer 温度问题。
- **主要结果**：若只取 far-field `H=-epsilon`，温度方程化为
  `Theta_YY+Pr Theta_Y=0`，故
  `Theta-1=C exp(-Pr epsilon eta)` 且
  `ell_T~1/(Pr epsilon)`。完整一致 outer 层满足
  `f=-h'/2`、`chi(Theta)(s+Ro_t g)^2=(s+Ro_t)^2`、
  `g''+chi[Ro_t h g'-(s+Ro_t g)h']=0`、
  `Theta''-Pr h Theta'=0`。
- **验证**：把 `F=0,G=1,H=-epsilon` 代回径向方程得到剩余
  `gamma(Theta-1)(s+Ro_t)^2/Ro_t`；因此在当前 `Ro_t≈-0.32` 且
  `Theta-1=O(1)` 时不能把常速度近似均匀用于整个 outer 层。随着
  `Y->infinity,Theta->1`，该剩余消失并恢复指数 far-far tail。
- **结论边界**：`Y->0` 的 `Theta_m,h,g` 数据必须由后续 matching 给出；本步骤
  没有推导 `A_R,A_T`，没有求解 coupled outer BVP，也没有新增 continuation。
- **主线意义**：确认 distinguished thermal coordinate 为
  `Pr epsilon eta`，同时揭示一致闭合会把 O(1) 温度尾耦合成慢 `G,H` 尾；简单
  标量指数解只在冻结轴向流或 far-far 子区间严格成立。
- **输出**：`docs/CONSISTENT_OUTER_THERMAL_TAIL.md`。
- **下一步**：按当前要求停止；inner--outer matching 留待单独任务。

### 2026-09-05 — L014：inner leading problem、正则匹配及退化核心检验

- **状态**：leading inner 方程及条件性匹配推导完成；仅对既有结果做代数与剖面
  后处理，未求解新的 BVP 或 continuation。
- **模型**：与 L013 相同的完整一致 BEK；不是冻结基本流模型。
- **参数与方法**：固定 eta 的 O(1) 展开，外层 Y=epsilon eta；假设正则外层入口
  与局域内层（integral H0 有限）。使用 L012 五点/四点二次拟合截距及 L011
  保存的五组 N=320 剖面，Julia 线性插值提取 eta=1,2,5 的场值。
- **问题**：内层无穷远应匹配什么状态；温度平台与方位速度平台是否独立；
  O(1) 内层核心是否可能退化。
- **主要结果**：inner leading ODE 保留完整动量项、壁面
  H0=F0=G0=0,T0=Tw_t；正则 overlap 为 H0,F0->0,G0->g_m,T0->Theta_m。
  积分温度方程并结合局域性给出 T0恒等于Tw_t，故 Theta_m=Tw_t，且
  g_m=[omega_t/sqrt(1-gamma(Tw_t-1))-s_t]/Ro_t。H0(infinity)=0只确定
  零阶匹配，不能确定外层系数h(0+)。
- **验证**：两组既有极限拟合分别给出 g_m=-2.9791e-5,-3.3629e-6；与
  零核心相容曲线 Tw=2-(omega/s)^2 的差为8.7059e-6,9.8275e-7。
  从epsilon=0.1到0.005，eta=1处G由0.074693降至0.0039318，
  Tw-T由0.024791降至0.0011558，H,F也趋近零；这些固定eta样本不是
  overlap平台的直接测量。
- **结论边界**：T0常数证明依赖正则匹配/内层积分有限；仅H0->0并不足够。
  H0~-a/eta的奇异重叠可能允许非恒温内层。本次代数与既有剖面支持
  H0=F0=G0=0,T0=Tw_t 的退化内层候选，但没有证明分支选择、唯一性，
  也没有由相容曲线确定Ro_t或fold条件。g_m微小负值不应解释为确定反转。
- **主线意义**：在正则结构下，内层温度平台就是极限壁温；慢外层承担到真实
  远场的O(1)温度与方位速度变化。O(1)核心并不必然保持非平凡流动。
- **输出**：docs/CONSISTENT_INNER_LEADING_MATCHING.md；
  work/scripts/check_consistent_inner_matching.jl；
  work/results/consistent_inner_leading_matching/。
- **下一步**：本步骤止于leading匹配条件。确定更高阶入口数据及分支选择需另做
  校正/可解性分析，尚未执行。

### 2026-09-05 — L015：退化零阶 inner core 的极限选择验证

- **状态**：在当前数据与正则匹配假设下通过；是数值极限选择验证，不是唯一性
  定理。
- **模型**：完整稳态加速度一致 BEK 既有 fold 解的代数/剖面后处理；不是冻结
  基本流、项开关或新 BVP。
- **参数与方法**：`Pr=0.72,gamma=1`；只用 L011 的五个
  `epsilon=0.1,...,0.005` 高阶共识点、L012 的两种拟合窗口和 L014 的
  `N=320` 固定-eta剖面。Julia 计算解析相容曲线、核心缺陷、全部数值半范围
  角点传播及固定-eta有效衰减指数。
- **问题**：fold--tail 极限是否真正选择
  `F0=G0=H0=0,T0=Tw,t`，亦即
  `chi(Tw,t)s_t^2=omega_t^2` 和 `g_m=0`。
- **主要结果**：四小点二次外推
  `(Ro_t,Tw,t)=(-0.320719730411,1.494786621674)`；只用 `Ro_t` 的解析曲线预测
  `Tw,t=1.494785638928`，差 `9.827e-7`。对应
  `Delta_core=-1.209e-6`、`Delta_core/omega_t^2=-1.945e-6`、
  `g_m=-3.363e-6`。五点二次外推温差为 `8.706e-6`；线性/二次及窗口变化均向
  零相容方向收敛。
- **验证**：四小点全部离散半范围角点传播给出
  `Delta_core in [-2.749e-4,2.725e-4]`，明确包含零，中心缺陷仅约包络半宽
  `0.44%`。有限-epsilon相对缺陷沿五点从 `-4.254e-3` 降至
  `-6.851e-6`。在固定 `eta=1,2,5`，最后两点的 `F,H` 有效衰减阶约 2，
  `G,Tw-T` 约 1，四场均趋零。
- **结论边界**：有限数据不能严格证明极限唯一或排除非正则 overlap；误差包络
  是确定性敏感性范围而非置信区间。相容曲线不独立选择 `Ro_t`，也不证明 fold
  条件或一阶 solvability。本步骤没有新基本流、continuation 或一阶 matching。
- **主线意义**：现有分支不只是给出很小的 `g_m`，而是同时满足解析径向平衡、
  误差包络和固定-eta剖面趋零；因此 O(1) inner BEK core 在 fold--tail 极限中
  塌缩，而 O(1) 热/方位调整被排入慢外层，是当前证据支持的核心理论结构。
- **输出**：`docs/DEGENERATE_INNER_CORE_VALIDATION.md`；
  `work/scripts/validate_degenerate_inner_core.jl`；
  `work/results/degenerate_inner_core_validation/`。
- **下一步**：按本任务要求停止；不在本记录中进入一阶展开或新增参数计算。

### 2026-09-05 — L016：leading outer 极限选择与增广 Jacobian 修正

- **状态**：leading outer BVP 已独立求解并收敛；附件提出的“基本流 augmented
  Jacobian 零模”假设被数值否定。
- **模型**：完整一致 BEK 的 leading thermo-rotational outer 极限；沿径向慢
  流形消去温度后的自洽 BVP，不是冻结场诊断。没有新增有限-epsilon基本流或
  continuation。
- **参数与方法**：`Pr=0.72,gamma=1`，`epsilon=-Hinf`，
  `Y=epsilon*eta`，`h(infinity)=-1`；Julia RK4 shooting，以
  `g_Y(0),Ro` 为两个未知量满足 `g(infinity)=1,h(infinity)=-1`；检查
  `Ymax=20,30,40`、步长 `0.01,0.005,0.0025` 和 Jacobian 差分步长。
- **问题**：leading outer BVP 能否独立选出 `(Ro_t,Tw,t)`；包含 `Ro` 的约化
  augmented Jacobian 是否奇异；`tau1/r1≈-1.21546` 是否来自其零模。
- **主要结果**：outer BVP 独立给出
  `Ro_t=-0.3207391663,Tw,t=1.4948092622,g_Y(0)=0.7881936107`，与四小点
  finite-epsilon外推只差 `1.94e-5,2.26e-5`，均在既有误差范围内。解析曲线
  斜率 `Tcurve'(Ro_t)=-1.21543131`，而拟合 `A_T/A_R=-1.21486673`，相差
  `5.65e-4`。
- **验证**：`Ymax=30->40` 的 `Ro_t` 变化约 `1.3e-10`；步长收敛至显示精度；
  末端残差 `<8e-15`。约化 shooting Jacobian
  `sigma_min=0.9550,sigma_max=3.4151,cond=3.576,det=3.2615`，差分步长
  `1e-4--1e-7` 下稳定，明显不奇异。
- **结论边界**：该 Jacobian 是消去 IVP 场自由度后的 leading outer基本流
  Jacobian，不是有限-epsilon固定-`Tw` Jacobian。结果不否定每个有限-epsilon
  点的 fold，而是否定把其零模直接等同于 leading outer基本流增广零模。
- **主线意义**：`epsilon` 已由 `-Hinf` 和 `h(infinity)=-1` 物理归一化；
  `tau1/r1` 是 `Tw=Tcurve(Ro)+O(epsilon^2)` 的链式法则切线斜率。若以后求
  `A_R`，必须展开包含固定-`Tw` 场 Jacobian 零向量条件的完整 fold 增广系统及 matching，不能直接
  套用当前基本流 BVP 的标准 fold Hessian 公式。
- **输出**：`docs/OUTER_LIMIT_SELECTION_AND_JACOBIAN.md`；
  `work/scripts/solve_outer_limit_and_jacobian.jl`；
  `work/results/outer_limit_selection_and_jacobian/`。
- **下一步**：按本步骤范围停止；没有执行 Hessian、adjoint solvability 或新的
  finite-epsilon计算。

### 2026-09-05 — L017：完整 fold 增广系统的一阶切向与零模缩放

- **状态**：一阶增广切向系统已实现并通过代数检查；切向方向和零模缩放已
  收敛，绝对极限幅值获强支持但最小 epsilon 尚未空间收敛。
- **模型**：完整稳态加速度一致 BEK 的有限-epsilon fold 系统
  `R=0,J_Uv=0,v'v=1,Hinf+epsilon=0`；完全自洽，不是 frozen/toggle 诊断。
- **参数与方法**：只读取既有 `epsilon=0.02,0.01,0.005` fold 解；
  `N=240,280,320`，`Pr=0.72,gamma=1`，真无穷映射。解析增广 Jacobian和方向
  Hessian，解 `D_Z F_aug Z_epsilon=(0,0,0,-1)`；左右零模投影和分辨率交叉检查。
- **问题**：`epsilon=-Hinf` 是否在完整 fold 系统一阶就固定
  `dRo/d epsilon,dTw/d epsilon`；finite-epsilon fold 零模如何进入 outer 极限。
- **主要结果**：`epsilon=0.02` 的三阶结果收敛至
  `dRo/d epsilon≈-0.43914,dTw/d epsilon≈0.53434`；`epsilon=0.01,N=320` 为
  `-0.429914,0.522779`。用这两个可解析点作导数线性截距给
  `(A_R,A_T)=(-0.420688,0.511219)`，与 branch-coordinate 拟合
  `(-0.420342,0.510660)` 相差 `0.082%,0.109%`。逐点切向比与左零模
  solvability slope 一致到 `O(10^-11)`。
- **验证**：所有点 `dHinf/d epsilon=-1`；增广线性残差约 `1e-10--4e-8`。
  按 `max|v_G|=1`，九组右零模统一满足
  `v_H=O(epsilon),v_F=O(epsilon^2),v_G,v_T=O(1)`；对应归一化常数跨
  epsilon和N只变化约百分之一级或更小。
- **结论边界**：`epsilon=0.005` 的绝对切向在 N 间从 `+0.0267` 到
  `-0.5484,-0.4618`，不可用；该点 fold 二次系数为 `O(1e-11)` 且符号随 N
  改变，属于既知分辨率诱导 cusp-like 层。故当前不宣称 `A_R,A_T` 已完全
  空间收敛，也未求解 epsilon=0 的解析一阶 matched BVP。
- **主线意义**：确认“outer BVP 选极限点、完整 fold 一阶系统选离开速度”的
  两级结构。有限-epsilon fold 零模经过分量重标度后具有稳定 slow-tail结构，
  解释了原系统奇异而 reduced outer BVP 正则。无需先做基本流二阶 Fredholm
  幅值选择。
- **输出**：`docs/FULL_FOLD_TANGENT_ANALYSIS.md`；
  `work/src/BEKConsistent.jl`；`work/scripts/analyze_full_fold_tangent.jl`；
  `work/scripts/summarize_full_fold_tangent.jl`；
  `work/results/full_fold_tangent_analysis/`。
- **下一步**：若要把 `A_R,A_T` 提升为独立收敛预测，应构造 epsilon=0 的解析
  matched tangent BVP 或采用 tail-adapted discretization；本步骤不进入二阶。

### 2026-09-05 — L018：resolved tangent window、零模 collapse 与数值 cusp 边界

- **状态**：三项数值证据已完成；给出保守收敛精度并修正 L017 对 `0.1%`
  匹配的引用方式。
- **模型**：完整自洽一致 BEK fold 系统及其 fixed-`Tw` 零模；不是冻结基本流
  或 term toggle。
- **参数与方法**：新增同一 fold locus 的 `epsilon=0.015`，
  `N=200,240,280,320`，并把已有 `N=320` 解升阶校正到 `N=360` 的
  `epsilon=0.02,0.015,0.01`；真无穷映射 `(2,0.6,0.5)`，`Pr=0.72,gamma=1`。
  完整增广隐式 tangent、`Y=epsilon*eta` mode collapse；另后处理既有
  `N=80,120,160` apparent-turn 数据。
- **问题**：是否存在可声明的 grid-resolved tangent window；fold eigenfunction
  是否具有真正 outer similarity profile；cusp-like 解析下限怎样随 N 后退。
- **主要结果**：`epsilon=0.02,0.015` 上 `N=240--360` 的 tangent 相对散布仅
  `1.69e-4,3.36e-4`；`epsilon=0.01` 增至 `6.89e-3`。两个 resolved 点的
  `N=360` 线性截距为 `(A_R,A_T)=(-0.4187145,0.5087672)`，在 `0.39%,0.37%`
  内确认 branch-fit `(-0.4203425,0.5106601)`。这是当前应引用的保守精度；
  L017 的 `0.08%,0.11%` 来自边缘 `N=320,epsilon=0.01`，不能标成完全收敛。
- **验证**：`N=360` 下三组 rescaled modes
  `v_G,v_T,v_H/epsilon,v_F/epsilon^2` 在 `Y` 上 collapse；0.015 对 0.01 的
  relative L2 差分别 `0.392%,0.686%,0.109%,4.90%`。左右零模 solvability与
  直接 tangent 比仍逐点一致。
- **数值边界**：apparent-turn 三点给有效幂 `2.267,2.064`，log fit
  `epsilon_min~N^-2.188`；`N^2 epsilon_min=110.95,99.55,97.74`。映射倒数
  第二节点解析满足 `eta_(N-1)~0.579N^2`，从而热尾解析要求直接解释
  `epsilon_min=O(N^-2)`。这是一条离散解析边界，不是物理 cusp law。
- **结论边界**：只有三个旧 apparent-turn 点，不能高精度拟合指数；
  `epsilon=0.01` 的绝对 tangent 尚未完全 N 收敛，`0.005` 仍排除。没有进行
  outer 二阶、稳定性或新物性计算。
- **主线意义**：full fold singularity 的 eigenvector 经 slow-tail 分量重标度后
  具有统一 outer similarity limit，而 normalized outer DAE 保持正则；这化解
  “原系统零特征值/极限 BVP 无零特征值”的表面悖论，并把 apparent cusp 定位为
  映射解析边界。
- **输出**：`docs/TANGENT_CONVERGENCE_AND_MODE_COLLAPSE.md`；
  `work/scripts/compute_intermediate_fold_tangent.jl`；
  `work/scripts/extend_fold_tangent_N360.jl`；
  `work/scripts/summarize_tangent_convergence_mode_collapse.jl`；
  `work/results/tangent_convergence_mode_collapse/`。
- **下一步**：当前 JFM 机制证据链已足够；除非需要更高精度独立 `A_R,A_T`，
  暂不进入 tail-adapted epsilon=0 matched tangent BVP 或 `O(epsilon^2)`。

### 2026-09-06 — L019：总体研究学习手册与代码--理论索引

- **状态**：文档综合完成；本条不包含新数值计算，也不改变 L001--L018 的证据等级。
- **范围**：依据现有 Julia 共享模块、执行脚本、专题报告、CSV 与图件，从等温
  BEK、两种热闭合和统一拓扑对照，整理到退化 inner core、outer DAE、完整 fold
  tangent、零模 collapse 及数值 cusp 分辨率边界。
- **主要结果**：生成 57 页中文循序学习手册，加入公式逐步推导、证据等级、常见
  混淆、四周学习路径、代码--方程映射、自测题与答案；明确采用 L018 的保守
  tangent 精度，并保留“outer BVP 正则、full fold 奇异”的修正解释。
- **结论边界**：手册只综合已有结果；没有新增基本流、延拓、稳定性、低马赫参考或
  matched-asymptotic 计算。未来结果若更新，仍须在本日志追加记录后再修订手册。
- **输出**：`docs/learning_manual/BEK_Boussinesq_Research_Learning_Manual.tex`；
  XeLaTeX 验证产物 `docs/learning_manual/BEK_Boussinesq_Research_Learning_Manual.pdf`。
- **下一步**：按手册四轮/四周路线重新手推和复现；科研计算本身仍按 L018 的范围
  与后续最小任务推进。

### 2026-09-06 — L020：热方程系数校正后的全链复算

- **状态**：已完成方程、导数和代表性拓扑复算；小 epsilon 极限仍有分辨率边界。
- **模型**：传统离心 Boussinesq 与一致加速度耦合模型；无穷映射；等温外部基流；不含温度物性项。
- **参数与方法**：Pr=0.72、gamma=1、N=80--280、a=2、b=0.6、c=0.5；Newton、固定 Hinf 增广 fold 和统一判据拓扑追踪。
- **问题**：检查从有量纲 BEK 缩放得到的温度对流系数是否改变既有 fold--tail 结论。
- **主要结果**：正确方程为 `T''+Pr*Ro*H*T'=0`。Ro=-1 传统 fold 保持 `Tw=1.048021731`；一致模型在 Ro=-0.5 和 -0.25 仍有 fold，在 Ro=-0.1 进入 thermal-tail。Ro>0 且 Hinf<0 时热模不局域。
- **验证**：一致模型方向 Hessian 相对误差 `1.68e-7`，Ro 导数误差 `4.43e-10`；N=280 epsilon 链给出 `(Ro_f,Tw_f)=(-0.5042793,1.6981602),(-0.4861967,1.6781872),(-0.4764253,1.6678403),(-0.4733199,1.6646192)`（epsilon=0.1,0.05,0.02,0.01）。修正 outer solve 在 Ymax=60 收敛到 `Ro_t=-0.4703674163,Tw_t=1.6615874195`。
- **结论边界**：旧的非 Ro=-1 热结果与旧一阶系数不再有效；epsilon=0.005 尚未空间收敛，不能据此拟合独立极限系数。定性 fold-to-tail 机制保留，但临界位置明显移动。
- **主线意义**：研究主线没有崩溃，反而暴露出必须修正的 BEK 缩放问题；后续论文只能使用 corrected_energy 输出。
- **输出**：`work/results/corrected_energy_*`、`docs/CORRECTED_ENERGY_RECOMPUTATION_2026-09-06.md`。
- **下一步**：完成 N=320/更高分辨率的小 epsilon 复核，并据此决定是否重做 fold tangent 与二维延拓；不修改学习手册。

### 2026-09-06 — L021：N=320 增广 tangent 与统一拓扑复核

- **状态**：已完成；学习手册未修改。
- **结果**：N=320 corrected fold tangent 在 epsilon=0.02,0.01 给出 `dRo/deps=-0.31217,-0.31033`、`dTw/deps=0.32544,0.32028`，比值分别 `-1.04252,-1.03205`；增广残差约 `4e-9`。
- **拓扑**：N=120 统一判据确认 traditional: `Ro=-1` fold、其余测试负 Ro 为 thermal-tail；consistent: `Ro=-1,-0.75` smooth，`Ro=-0.5,-0.25` fold，`Ro=-0.1` thermal-tail。
- **边界**：consistent `Ro=-0.25` 的 fold 位于 `Hinf≈-0.0381`，接近热尾，需更高分辨率后才可作精确临界值；正 Ro 由 corrected energy 的远场指数不局域，不纳入非平凡热解 continuation。
- **输出**：`work/results/corrected_energy_full_fold_tangent/N320`、`work/results/corrected_energy_two_model_topology/final_N120`。

### 2026-09-06 — L022：研究主线复核与学习手册同步

- **状态**：已完成逻辑复核并更新手册；未改变历史数据。
- **复核结论**：热方程系数修正改变了所有非 Ro=-1 的定量 fold 坐标、outer 极限坐标和 tangent 系数，但没有改变主线的四个结构性结论：闭合依赖的拓扑重组、一致模型 fold 向 thermal-tail 迁移、outer slow-manifold 与有限 epsilon fold 的两层结构，以及 fold-to-tail 的奇异重构。
- **手册更新**：同步修正温度方程、热尾长度 `1/[Pr*(-Ro)*epsilon]`（Ro<0）、outer reduced ODE、corrected epsilon 链、outer 极限候选和 tangent 表；明确旧系数已废止、epsilon=0.005 仍不用于独立极限系数。
- **验证边界**：XeLaTeX 已成功处理至图件加载阶段；当前工作目录缺少手册引用的若干历史 PNG，导致完整 PDF 覆盖编译无法完成，但本次修改未产生 LaTeX 语法错误。
- **输出**：`docs/learning_manual/BEK_Boussinesq_Research_Learning_Manual.tex`。

### 2026-09-07 — L023：学习手册 corrected 版与 WSL 构建修复

- **状态**：完成；本条取代 L022 中“PDF 因缺图未能完整编译”的构建状态。
- **主线复核**：解析相容曲线在 corrected outer 点满足 `Tcurve(Ro_t)-Tw_t=O(1e-13)`，且 `Tcurve'(Ro_t)=-1.021771561`；N=320 tangent 比值向该斜率逼近。退化核心、outer slow manifold 与两级选择机制因此保留。
- **手册修改**：删除所有未校正二维、scaling、tangent、mode-collapse 和 cusp 图表的定量引用；保留可由 corrected 计算支持的结论，并将尚待重算部分明确降级。
- **构建修复**：新增 `.vscode/settings.json` 和 `work/scripts/xelatex_windows_wsl.sh`，通过 Windows `pushd` 映射 WSL UNC 目录后调用现有 TeX Live，规避 LaTeX Workshop 找不到 Linux executable 及 `xdvipdfmx` 不接受 UNC 当前目录的问题。
- **验证**：recipe 实测退出码 0；两遍 XeLaTeX 生成 54 页 PDF 与 SyncTeX；抽查目录、方程、fold 链、inner 和 outer 页面，无裁切、重叠或乱码。
- **输出**：`docs/learning_manual/BEK_Boussinesq_Research_Learning_Manual.tex`、同名 PDF、`.vscode/settings.json`、`work/scripts/xelatex_windows_wsl.sh`。

### 2026-09-07 — L024：学习手册独立 outer 核验与证据审计

- **状态**：完成独立 leading outer 计算与手册算术核对；审计报告继续整理。此条不覆盖历史结果。
- **模型与范围**：corrected leading outer DAE，Pr=0.72、gamma=1，单调 g、q 和 chi 非零；以 g 为自变量的线性 p,h 系统在 [0,1] 积分，p(1)=0 选 r，h(1)=-1 归一化；不是完整有限 epsilon fold 重算或物理模型交叉验证。
- **参数与验证**：RK4 250、500、1000、2000 步；r∈[-0.48,-0.46] 二分。后两网格所报坐标在约 1e-13 内一致。得到 `Ro_t=-0.470367416464`、`Tw_t=1.661587419719`、`p0=0.654193538532`、`Tcurve'=-1.021771560631`。与原 corrected Ymax=60 shooting 坐标差约 2e-10，支持这一归一化约化问题的数值结果。
- **审计结论**：PDF 第 8 章仍引用旧 outer Jacobian 数字；第 9 章旧切向外推数值不能由当前表格推出。当前表格两点线性截距仅在算术上为 `-0.308496794732`、`0.315114319865`，不作为已收敛极限系数。有限 epsilon=.02、.01 的伴随斜率与解析极限斜率仍相差约 2.03%、1.01%；离散伴随恒等式的 1e-11 精度不等于渐近精度。
- **边界**：此独立计算不证明完整 fold 链收敛、统一匹配、全局唯一性、物理 cusp 的排除或三维稳定性。映射远端 dx 条件及新旧图表混用的审计结论见专题报告。
- **输出**：`work/scripts/audit_learning_manual_20260907.jl`、`work/results/learning_manual_audit_20260907/`、`docs/LEARNING_MANUAL_AUDIT_2026-09-07.md`。

### 2026-09-07 — L025：54 页学习手册方程与结论审计完成

- **状态**：完成全部编号方程的核验清单、主干推导重建、关键 PDF 页面渲染检查和证据追溯；原 PDF/TeX、共享求解模块和历史数据未修改。
- **明确更正**：L023 关于已清除旧图表的描述不完整：编译日志确认图 5.2 仍加载旧 consistent 曲线；式 (8.11) 沿用旧 outer Jacobian 数字；式 (9.10) 沿用旧 N360、epsilon=.02/.015 外推系数。L018 的 0.4% 精度不能迁移为 corrected tangent 精度。
- **模型/机制结论**：corrected 温度方程、加速度闭合、inner 相容、outer ODE 及伴随关系在所述假设下代数成立。L024 独立计算支持归一化 outer 候选点。程序实际采用 H_x(-1)=0 的映射端点正则性条件，手册不能将其视作与物理 H_eta(infinity)=0 无条件等价。
- **证据范围**：finite epsilon corrected 数据支持候选热尾机制和代表切片拓扑；不证明完整参数平面分类、统一匹配、绝对幅值极限、零模函数收敛、物理 cusp 排除或实际流动物理准确性。本次没有重跑全套高 N fold 和二维延拓。
- **输出**：`docs/LEARNING_MANUAL_AUDIT_2026-09-07.md`，其中逐式列明通过、条件成立、错误或证据不足；独立计算见 L024。

### 2026-09-07 — L026：corrected fold tangent 跨网格补算

- **状态**：完成；本条补充并收紧 L021 的单网格 tangent 证据。
- **模型与方法**：acceleration-consistent corrected-energy BEK；`Pr=0.72`、`gamma=1`、映射 `(a,b,c)=(2,0.6,0.5)`；`N=200,240,280,320`；`epsilon=0.02,0.015,0.01`；完整固定 `Hinf=-epsilon` 增广 fold 及其隐式切向，容差 `1e-9`。输入为 corrected fold profile；没有使用旧温度方程的 tangent 数据。
- **结果**：N280 与 N320 的 `(dRo/depsilon,dTw/depsilon)` 相对差在 epsilon=.02 为 `(0.1543%,0.1543%)`，epsilon=.015 为 `(0.1380%,0.1376%)`，epsilon=.01 增为 `(1.3669%,1.3704%)`。N320 的斜率依次为 `-1.0425152608,-1.0372555395,-1.0320453221`，相对 endpoint 解析斜率 `-1.0217715606` 的距离为 `2.030%,1.515%,1.005%`。直接切向与离散伴随斜率在所报位数内一致。
- **模式幅值**：N320 的 `max|vH|/max|vG|` 为 `.07183,.05381,.03584`，`max|vF|/max|vG|` 为 `.001748,.001001,.0004535`，与 `O(epsilon)`、`O(epsilon^2)` 分量幅值一致；这不是跨 epsilon 的函数空间 profile-collapse 证明。
- **结论边界**：epsilon=.02 和 .015 有两张高网格的有限 epsilon 一致性；epsilon=.01 的绝对切向未空间收敛。现有数据不足以给 `A_R,A_T` 的 epsilon-to-zero 极限或误差条。离散伴随恒等式的高精度不等于渐近极限精度。
- **输出**：`docs/CORRECTED_MANUAL_SUPPLEMENT_2026-09-07.md`、`work/scripts/corrected_manual_tangent_convergence_20260907.jl`、`work/scripts/summarize_corrected_manual_tangent_20260907.jl`、`work/results/corrected_manual_tangent_convergence_20260907/`。

### 2026-09-07 — L027：学习手册按审计结果修订并重新发布

- **状态**：完成；本条记录 L025 审计与 L026 补算之后的手册修订状态。
- **修改范围**：修正图 5.2 的证据版本、outer Jacobian 数字、有限 epsilon 方程写法和 tangent 表；撤下旧 `(-0.4187145,0.5087672)` 极限系数；澄清物理无穷远条件与程序映射端正则条件的区别；补入 simple-fold 非退化条件、Rossby/旋转量定义、传统闭合中 `gamma=1` 的实现限制和大温差物理适用性警告。
- **结论收紧**：不再把二维网格尺度解释为 corrected cusp 的证明或排除；不再声称 `A_R,A_T` 已收敛、零模剖面已经塌缩或模型具有全三维稳定性。代表性切片、归一化 outer 候选点和有限 epsilon tangent 分别按其实际证据范围陈述。
- **补算证据**：采用 L026 的 N=200--320、epsilon=.02/.015/.01 corrected full-fold tangent；手册仅报告跨网格可支持的有限 epsilon 数值和分量量级，不进行无误差条的 epsilon-to-zero 外推。
- **构建与版面验证**：两遍 XeLaTeX 成功生成 55 页 PDF；渲染检查全部页面，并对方程、图 5.2、新 tangent 表、结论和证据索引所在关键页再次放大核对，未发现裁切、重叠或乱码。文本核查确认过期 Jacobian 数字与畸形 `qquad` 已清除。
- **输出**：`docs/learning_manual/BEK_Boussinesq_Research_Learning_Manual.tex`、同名 PDF、`docs/LEARNING_MANUAL_AUDIT_2026-09-07.md`、`docs/CORRECTED_MANUAL_SUPPLEMENT_2026-09-07.md`。
- **仍开放**：corrected 二参数 cusp 搜索、有限 epsilon 到 outer 的统一匹配、`A_R,A_T` 极限及误差条、小温差/低马赫物理有效性检验和非相似三维稳定性；这些不能由本次修订推断。

### 2026-09-07 — L028：学习手册重构为单线递进式教程

- **状态**：完成教学结构重构；没有新增或改变数值模型及研究结论。
- **问题**：L027 版本虽然内容齐全，但有限 epsilon 数值问题、inner/outer 极限和计算方法交叉出现，初学者难以辨认当前坐标、量级、输入和输出。
- **结构修改**：新增“一条线读完全书”的六步研究路线、全书统一量级账本及阅读门槛；明确先完成物理定义、有限参数求解、fold 认证和已知数值事实，再由尾长推导 Y=epsilon*eta，最后连接 inner、outer DAE 与完整增广 fold 切向。
- **推导补充**：逐项从有量纲温度输运推出 T''+Pr*Ro*H*T'=0；列出 F,H,G,T 在 outer 坐标下的零阶、一阶和二阶导数量级；逐行解释径向黏性项为何为 epsilon^4、代数慢流形为何出现，以及怎样由温度方程和方位方程消元得到三个一阶 outer ODE。
- **计算方法补充**：用统一六步流程说明 Newton 基本解、伪弧长预测校正、候选 fold 检测、增广零模求解、非退化与跨网格验证、thermal-tail 分类之间的关系。
- **版面验证**：两遍 XeLaTeX 生成 58 页 PDF；按 10 页一组渲染检查全部页面，并放大核对路线图、量级表、fold 流程和 outer 推导页。检查点跨页造成的空白页已消除；未见裁切、重叠、乱码或 overfull 警告。
- **证据边界**：本次是教学表达和推导展开，不提升 L026--L027 的证据等级；corrected cusp、统一匹配、A_R,A_T 极限、物理有效性和非相似三维稳定性仍开放。
- **输出**：docs/learning_manual/BEK_Boussinesq_Research_Learning_Manual.tex 及同名 58 页 PDF。

### 2026-09-07 — L029：以“物理问题如何被解决”为中心重写学习手册

- **状态**：完成；本条取代 L028 的教学组织状态，但不改变 L026--L027 的数值证据与结论边界。
- **核心问题**：明确研究对象不是推导 ODE 本身，而是判定加热 BEK 旋转边界层的稳态支在改变 Ro 与 Tw 时以 smooth、fold 还是 thermal tail 终止，以及密度--加速度闭合如何重组该拓扑。
- **研究主线**：新增完整主线章，逐步说明相似性 BVP 的用途、固定 Ro 升温实验的分类决策树、每个拓扑表格条目的计算来源、为什么五条切片之后必须沿 fold 本身移动、固定 epsilon 增广系统的未知量和方程、长尾为何使有限网格外推失效，以及 inner、outer、完整 fold 切向分别填补的证据缺口。
- **5.4/5.5 补全**：以 consistent Ro=-0.5 为例列出从等温种子、固定 Tw 延拓、伪弧长穿越、固定 Hinf 局部参数、零模和非退化检查到 fold 坐标的全过程；明确 epsilon 链每一点都是对 (U,v,Ro,Tw) 的 2n+2 维增广 BVP 重新校正，不是插值或图上取极值。
- **尺度推导**：从远场温度指数严格得到 LT=O(epsilon^-1) 与 Y=epsilon*eta；由 G、T 的 O(1)边界跳跃和 Hinf=-epsilon 得到 outer 的 G,T,H 尺度；由连续性推出 F=epsilon^2 f；逐项代回温度、方位和径向方程验证 distinguished balance；再将光滑 outer 入口展开代回固定 eta，推出 inner 的 G、T空间变化、F、H 尺度。参数的 O(epsilon) 幂次被明确标为完整 fold 切向所验证的正则性，而不是由尾长单独推出。
- **详细推导保留**：保留并扩展温度无量纲化、inner 径向代数项 O(1)/O(epsilon) 展开、outer 导数账本、slow-manifold 退化和 DAE 到三个一阶 ODE 的消元。
- **构建与验证**：两遍 XeLaTeX 生成 66 页 PDF，无 overfull、未定义引用或编译错误；按 10 页一组渲染全部页面并放大检查主线章、尺度推导、inner/outer 和 fold 链页面。未发现裁切、重叠、乱码或检查点孤立占页。
- **证据边界**：本次仅重构研究逻辑并细化已有数学推导，不新增模型计算，不提升 cusp、AR/AT 极限、物理有效性或三维稳定性的证据等级。
- **输出**：docs/learning_manual/BEK_Boussinesq_Research_Learning_Manual.tex 及同名 66 页 PDF。

## 新记录模板

复制以下模板到本文件末尾，不要覆盖历史条目：

```markdown
### YYYY-MM-DD — LXXX：简短标题

- **状态**：已验证 / 暂定 / 探索性 / 被 LYYY 替代。
- **模型**：闭合、几何、是否自洽或冻结诊断。
- **参数与方法**：Ro、Tw、Pr、gamma、N、映射、容差、延拓方法。
- **问题**：本次计算要判定什么。
- **主要结果**：数值和机制。
- **验证**：收敛、残差、Jacobian、交叉模型或独立实现。
- **结论边界**：没有证明什么；参数和分支范围。
- **主线意义**：它如何支持、修正或否定当前论文叙事。
- **输出**：代码、数据、图和专题报告路径。
- **下一步**：由本结论直接导出的最小后续任务。
```
### 2026-09-07 — L030：补齐 inner--outer 匹配的逻辑前提

- **状态**：完成。针对学习手册在式 (7.13) 前直接使用 outer 入口值、却未先解释匹配原理的问题，重新组织了相关章节。
- **原有断层**：此前没有明确区分物理壁面 $\eta=0$ 与 outer 坐标的入口极限 $Y\to0$，也没有定义重叠区及“共同部分匹配”，因而 $f(0)$、$g(0)$、$h(0)$、$\theta_0(0)$ 的来源显得像未经说明的边界条件。
- **本次修改**：在 outer 入口展开和式 (7.13) 之前新增“先定义匹配：outer 的 $Y=0$ 不是物理壁面”一节；从同一个有限 $\epsilon$ 解的 inner/outer 两种表示出发，用 $\eta=\epsilon^{-\alpha}$（$0<\alpha<1$）显式得到 $1\ll\eta\ll\epsilon^{-1}$；逐项解释 $g(0)=0$、$\theta_0(0)=T_{w,t}$、$h(0)=0$ 及允许 $f(0)=f_0\ne0$ 的原因；并在式 (7.13) 后补出 $B_G=p_0$、$H_2\sim h_1\eta$、$B_T=\theta_1$ 与 $h_1+2f_0=0$ 等系数匹配关系。
- **研究逻辑边界**：该节把“候选 inner 主阶结构的提出”与“下一章代回方程验证”明确分开；$h(0)=0$ 使用了重叠区中的正则增长假设及有限 $\epsilon$ 解的支持，不把它表述成脱离假设的全局一致收敛定理。
- **构建与核验**：XeLaTeX 连续构建两次，得到 68 页 PDF；渲染并逐页检查全部 68 页，未见裁切、重叠或乱码；最终日志中无 overfull、未定义引用或 LaTeX 致命错误。
- **输出**：`docs/learning_manual/BEK_Boussinesq_Research_Learning_Manual.tex`；`docs/learning_manual/BEK_Boussinesq_Research_Learning_Manual.pdf`。

### 2026-09-07 — L031：渐近章节改为先 inner、后 outer、再匹配

- **状态**：完成；本条取代 L030 所述教学编排，未改变数值模型、参数或已有计算结论。
- **发现的问题**：L030 虽补入 matching 定义，仍在 inner 方程验证之前用 outer 入口展开推出 inner 幂次，随后又借用该 inner 幂次建立 outer 入口，存在前后依赖；旧式 (7.21d) 还把共同部分常数 $f_0$误写成固定 $\eta$下的整个 inner 函数。
- **本次重构**：第七章只保留远场热指数 $L_T=[\Pr(-\Ro)\epsilon]^{-1}$作为长尺度线索；第八章先由有限 $\epsilon$剖面提出候选 inner 幂次，再逐阶代回方程得到退化核心、解析相容曲线和线性增长，并由 $\epsilon\eta=O(1)$暴露失效尺度；第九章随后建立 outer scaling、DAE 与约化 ODE，最后才定义 overlap 并逐项匹配。
- **公式修正**：明确 $F^{\mathrm{in}}=\epsilon^2F_2(\eta)$、$F_2(0)=0$，匹配只要求 $F_2(\eta)\to f_0$（inner 远端的共同部分）；同时区分 leading outer 温度入口常数 $T_{w,t}$与下一阶参数修正 $\epsilon\tau_1$。
- **证据边界**：inner 幂次首先由有限 $\epsilon$数据提出，再由逐阶方程相容性检验；$h(0)=0$仍依赖 $H_2$在 overlap 至多线性增长的正则性假设。此次为推导逻辑与公式适用范围修订，没有建立统一误差界或新增数值证据。
- **构建与核验**：两遍 XeLaTeX 后生成 65 页 PDF；全部页面已渲染为 JPEG 并以联系表检查，另放大检查 inner、outer 和 matching 页面。最终交叉引用与版面检查见本次构建日志。
- **输出**：`docs/learning_manual/BEK_Boussinesq_Research_Learning_Manual.tex`；`docs/learning_manual/BEK_Boussinesq_Research_Learning_Manual.pdf`。

### 2026-09-08 — L032：导师汇报前的 outer 与解析导数独立复核

- **状态**：已验证；新结果写入独立审计目录，未覆盖 baseline 或既有结果。
- **模型与方法**：corrected-energy acceleration-consistent leading outer DAE，`Pr=0.72`、`gamma=1`。分别使用 `Y` 坐标 shooting 与以单调 `g` 为自变量的 RK4/bisection 独立实现；另对完整有限维残差的解析方向 Hessian 和 Rossby 导数作中心差分检查。
- **主要结果**：两套 outer 实现分别得到 `(Ro_t,Tw_t,p0)=(-0.470367416261,1.661587419512,0.654193538967)` 与 `(-0.470367416464,1.661587419719,0.654193538532)`，差值均小于 `4.4e-10`；相容曲线斜率为 `-1.021771560631`。`Ymax=60` 的 shooting 残差为 `2.14e-12`，入口 shooting Jacobian 最小奇异值约 `0.982`。解析方向 Hessian 的相对差分误差为 `1.5e-9`（步长 `1e-4`），Rossby 导数最低误差为 `4.1e-10`。
- **结论边界**：验证的是归一化 leading outer 问题及代码导数，不证明完整 finite-epsilon fold、全局唯一性、物理 Boussinesq 有效性或三维稳定性。
- **输出**：`work/scripts/presentation_evidence_audit_20260908.jl`；`work/results/presentation_evidence_audit_20260908/`。

### 2026-09-08 — L033：corrected fold、基本流与零模的分辨率复核

- **状态**：部分验证并收紧结论；`epsilon=0.005` 继续排除于定量结论之外。
- **模型与方法**：corrected-energy acceleration-consistent full fold 增广系统，`Pr=0.72`、`gamma=1`、映射 `(a,b,c)=(2,0.6,0.5)`；`N=240,280,320`，`epsilon=0.02,0.015,0.01`；复算 fold 切向，并将基本流与归一化 outer 解在 `0<=Y<=10` 比较。零模按 `vG,vT,vH/epsilon,vF/epsilon^2` 重标度，以 `epsilon=0.01` 为参考比较剖面。
- **主要结果**：`epsilon=0.02` 的 fold 坐标和绝对切向在三张高网格上稳定；`epsilon=0.015` 的 `N=280,320` 结果相近；`epsilon=0.01` 的 fold 位置仍接近，但绝对切向在 `N=240--320` 间明显非单调，不能称空间收敛。`N=320` 基本流从 `epsilon=.02` 到 `.01` 时，`G,T,H/epsilon` 对 leading outer 的 relative-L2 误差分别由 `0.407%,0.148%,0.815%` 降至 `0.175%,0.093%,0.406%`。零模 `G,T,H` 的跨 epsilon relative-L2 差约为 `0.99%,0.72%,0.31%`（`.02` 对 `.01`），但 `F/epsilon^2` 为 `12.0%`；`.015` 对 `.01` 的 `F` 差仍为 `6.75%`。
- **修正结论**：基本流 outer 极限已有强数值支持；零模只有 `G,T,H` 分量显示清楚的 emerging collapse，现阶段不能声称四分量 full eigenfunction 已统一收敛。`epsilon=.005` 的 inner `F/epsilon^2,H/epsilon^2` 已表现欠分辨，只能作为位置趋势提示。
- **输出**：`work/scripts/corrected_manual_tangent_convergence_20260907.jl`（新增可选 `--outdir`）；`work/scripts/summarize_asymptotic_profile_evidence_20260908.jl`；`work/results/presentation_evidence_audit_20260908/tangent/` 及同目录 CSV。

### 2026-09-08 — L034：代表性分支拓扑可重复性复算

- **状态**：已验证代表性切片；与 L021 的分类一致。
- **模型与方法**：traditional 与 acceleration-consistent corrected-energy 两模型，`Pr=0.72`、`gamma=1`、`N=120`，`Ro=-1,-0.75,-0.5,-0.25,-0.1`；固定 `Tw` continuation、必要时伪弧长穿越、统一 fold/tail/smooth 判据。
- **主要结果**：traditional 为 `fold, tail, tail, tail, tail`；consistent 为 `smooth, smooth, fold, fold, tail`。consistent 的两处 fold 复算为 `(Tw,Hinf)=(1.693338220638,-0.088697406043)` 与 `(1.586291890050,-0.038121687604)`；与旧表的温度坐标在约 `1e-11` 内复现，近奇异方向上的 `Hinf` 微小差异不改变分类。
- **结论边界**：这是五个 `Ro` 切片上的 N=120 拓扑复算，不是二维参数平面的完整分类或高精度临界曲线；`Ro=-0.25` 的 fold 靠近 tail，精确坐标仍需高网格确认。
- **输出**：`work/scripts/compare_traditional_consistent_topology.jl`（新增可选 `--outdir`）；`work/results/presentation_evidence_audit_20260908/topology/audit_N120/`。

### 2026-09-08 — L035：导师汇报前证据审计与学习手册结论收紧

- **状态**：完成；综合 L032--L034，不新增数值计算。
- **修订**：在第十章明确指出第八、九章基本流渐近方程不包含 `J_Uv=0`，因此只研究 tail 极限，fold 约束来自完整增广系统；fold--tail 联系目前是端点、切向和零模分量结构的联合数值--渐近证据，不是统一误差界下的存在定理。加入 finite-epsilon 基本流对 outer 的 relative-L2 误差，并把“完整零模已有 outer 相似极限”降为 `G,T,H` emerging collapse、`F` 分量尚未收敛。
- **证据分级**：outer 端点、代表切片拓扑、基本流剖面收敛和 resolved finite-epsilon folds 可用于导师汇报；`epsilon=.01` 绝对切向、`epsilon=.005`、四分量零模极限、cusp、三维稳定性和物理 Boussinesq 准确性不得提升为已验证结论。
- **构建与版面**：两遍 XeLaTeX 成功生成 67 页 PDF；渲染全部 67 页并以七张联系表逐页检查，另放大核对第九章新基本流误差表及第十章 fold 约束、切向和零模页，未见裁切、重叠或乱码。
- **输出**：`docs/PRESENTATION_EVIDENCE_AUDIT_2026-09-08.md`；`docs/learning_manual/BEK_Boussinesq_Research_Learning_Manual.tex` 及同名 PDF。

### 2026-09-08 — L036：13 点 Rossby 加密拓扑扫描（第一阶段）

- **状态**：探索性复算完成；出现非单调分类，必须经第二阶段升网格和边界加密后解释。
- **模型与方法**：traditional 与 acceleration-consistent corrected-energy 模型；`Pr=0.72`、`gamma=1`、`N=120`，`Ro=-1,-.9,-.8,-.7,-.6,-.5,-.4,-.3,-.25,-.2,-.15,-.1,-.05`；其余 continuation 和统一分类判据同 L034。
- **主要结果**：traditional 仅 `Ro=-1` 为 fold，其余 12 点均为 tail。consistent 在 `-1,-.9,-.8` 为 smooth，`-.7,-.6,-.5` 为 fold，`-.4,-.3` 为 tail，`-.25,-.2` 又为 fold，`-.15,-.1,-.05` 为 tail。
- **关键发现**：原五点表中每个采样点均复现，但 consistent 的分类不是由五点表可能暗示的单一 smooth--fold--tail 区间；`-.25,-.2` 的第二 fold 窗口可能是真实的额外拓扑结构，也可能是长尾欠分辨或 continuation 分支选择造成的表观窗口。
- **结论边界**：本条不得用于声称第二 fold 窗口真实存在；N=120 且靠近 `Hinf=0` 的结果只能用于定位复核区间。`smooth` 仅表示 `Tw<=1.99` 窗口内未见终止事件。
- **输出**：`work/results/presentation_evidence_audit_20260908/rossby_dense/stage1_N120/`。
- **下一步**：在四处分类转换附近加密 Ro，并用 `N=160` 复算第二窗口及相邻 tail 点；比较 fold 的网格漂移和 thermal-tail 解析长度。

### 2026-09-08 — L037：Rossby 分类转换区间加密（N=120）

- **状态**：探索性边界加密完成；确认需要升网格判断的窄转换区。
- **参数**：在第一阶段四处分类变化附近增加 `Ro=-.975,-.95,-.925,-.75,-.725,-.475,-.45,-.425,-.30,-.2875,-.275,-.2625,-.175,-.1625` 等点；统一 `N=120`。
- **主要结果**：consistent 的主 fold 链在 `Ro=-.475` 仍给 fold (`Hinf=-0.01779`) 而 `-.45` 为 tail，与独立 outer 候选 `Ro_t=-0.47037` 相符。第二窗口在 `-.30` 为 tail、`-.2875` 已为 fold，并持续到 `-.175`；`-.1625` 转为 tail。traditional 的 `-.975` 在 `Tw=1.99` 时 `Hinf=-.00173`，被当前规则记为 smooth，说明该处分类受 `Tw<=1.99` 截断影响，不能解释成新的稳态拓扑。
- **结论边界**：第二 fold 窗口在 N=120 上不是单个离群点，但其 `|Hinf|` 仅约 `0.03--0.05`，仍可能受热尾解析限制；主 fold--tail 边界和第二窗口均需 N=160 复算。此条不进入手册主结论。
- **输出**：`work/results/presentation_evidence_audit_20260908/rossby_dense/stage2_N120/`。

### 2026-09-08 — L038：分类边界的 N=160 复核

- **状态**：完成；第二 fold 窗口未通过网格位置一致性检验，暂判为长尾欠分辨产生的表观结构。
- **参数与方法**：两模型，`N=160`，`Ro=-.975,-.725,-.475,-.45,-.30,-.2875,-.20,-.175,-.1625`；方法同 L036--L037。
- **主要结果**：consistent 主边界保持为 `Ro=-.475` fold、`-.45` tail；前者的 `Hinf` 从 N120 的 `-0.01779` 漂移至 N160 的 `-0.01487`，符合接近 `Ro_t=-0.47037,Hinf=0` 时的高分辨率需求。第二窗口则不保持：`Ro=-.2875` 从 N120 fold 变为 N160 tail，`Ro=-.1625` 从 N120 tail 变为 N160 fold；`-.20,-.175` 虽仍为 fold，但 `Hinf` 分别漂移至 `-0.03371,-0.03819`。
- **解释**：分类窗口随 N 明显移动而不是坐标收敛，是 tail 长度超过映射有效解析范围的典型信号。它不能作为第二个物理 fold 岛的证据。
- **边界**：主边界稳定分类仍只由少数点夹逼；第二窗口需至少 N=200 再检查漂移方向，之后才能决定从 PPT 主结论中删除还是保留为数值警示。
- **输出**：`work/results/presentation_evidence_audit_20260908/rossby_dense/stage3_N160/`。

### 2026-09-08 — L039：表观第二 fold 窗口的 N=200--280 排除检验

- **状态**：完成；N=120 的第二 fold 窗口未通过网格收敛，判为长尾欠分辨造成的数值伪结构。
- **模型与方法**：acceleration-consistent corrected-energy，使用同一 continuation 与分类判据；N=200 检查 `Ro=-.30,-.2875,-.25,-.20,-.175,-.1625,-.15`，N=240 检查 `-.25,-.20,-.175,-.15,-.125`，N=280 检查 `-.25,-.20,-.15`。
- **主要结果**：`Ro=-.25` 在 N120 为 fold，到 N200、240、280 均为 tail；`Ro=-.20` 在 N120--200 为表观 fold，到 N240、280 为 tail；`Ro=-.15` 在 N200、240 为表观 fold，到 N280 为 tail。表观窗口随 N 增大持续向更小 `|Ro|` 移动并在固定 Ro 上消失。
- **结论**：没有证据支持第二个物理 fold 岛。该窗口是分辨率依赖的 apparent fold，不应进入物理拓扑图。
- **主线影响**：原五点中 `Ro=-.25` 的 consistent-fold 分类被否定；可靠主线应改为：存在一条从有限 `Hinf<0` fold 向 `Ro_t≈-0.47037,Hinf=0` 终止的主 fold 链，而较小 `|Ro|` 的粗网格 fold 判定属于数值伪结构。
- **输出**：`work/results/presentation_evidence_audit_20260908/rossby_dense/stage4c_N200/`、`stage5_N240/`、`stage6_N280/`。

### 2026-09-08 — L040：两处非长尾转换的局部加密

- **状态**：完成局部夹逼；转换机制仍需二维 fold continuation。
- **traditional 近 Ro=-1**：N=120/160 均确认 `Ro=-.995` 为有限 fold，坐标约 `(Tw,Hinf)=(1.05294,-0.49720)`；`Ro=-.99` 两张网格均在远离 tail 处停止但未得到可靠 fold crossing，标为 unresolved；`-.985,-.975` 在 `Tw=1.99` 上限前未见事件。因此“传统 fold 只存在于精确 Ro=-1”被否定，转换被夹在 `[-.995,-.985]`，但切片法尚不能判断是 cusp、分支连接还是需要更稳健的延拓。
- **consistent smooth--fold**：N=120 加密给 `Ro=-.710` smooth、`-.705` fold、`-.700` fold；N=160 复核同一分类。当前仅将左边界夹在 `[-.710,-.705]`；靠边界的弧长斜率很小，精确坐标不作定量声明。
- **主线影响**：加密采样支持 consistent 中等负 Ro 的主 fold 区和 `Ro_t≈-.47037` 的唯一网格收敛 fold--tail 终点；同时暴露 traditional 近 `-.99` 的独立开放问题。
- **输出**：`work/results/presentation_evidence_audit_20260908/rossby_dense/stage7_N120/`、`trad_m0985_N120/`、`stage8_N120/`、`stage9_N160/`；N160 `Ro=-.995,-.99` 的命令行复核结果见本条数值记录。

### 2026-09-08 — L041：加密拓扑结论写回手册与审计报告

- **状态**：完成；综合 L036--L040，不新增模型计算。
- **修正**：删除“consistent `Ro=-.25` 为物理 fold”的旧结论；将其与 `-.20,-.15` 的低网格 fold 统一标为随 N 移动并消失的 apparent fold。手册改为报告经过跨网格复核的主 consistent fold 区、`Ro_t≈-.470367` 的 fold--tail 端点，以及 traditional `Ro≈-.99` 的未解析窄转换。
- **图件**：新增 `second_window_grid_sensitivity.csv/png`，直接展示 N=120--280 分类不收敛；取代手册中可能误导的旧五点 N=120 拓扑图。
- **构建验证**：学习手册两遍 XeLaTeX 成功生成 68 页 PDF；最终日志无 overfull、未定义引用或致命错误，修改页已重新渲染检查，无裁切、重叠或乱码。
- **输出**：`docs/PRESENTATION_EVIDENCE_AUDIT_2026-09-08.md`；`docs/learning_manual/BEK_Boussinesq_Research_Learning_Manual.tex/pdf`；`work/scripts/summarize_dense_rossby_audit_20260908.jl`；`work/results/presentation_evidence_audit_20260908/rossby_dense/summary/`。

### 2026-09-08 — L042：生产五点表数据完整性恢复与最终范围确认

- **状态**：完成；数据恢复不改变科学结论。
- **数据核验**：重新按原 `N=120`、`Ro=(-1,-.75,-.5,-.25,-.1)`、两模型命令生成 `work/results/corrected_energy_two_model_topology/production_N120/topology_table.csv`；恢复为 10 行（含表头）并与 `final_N120` 的拓扑分类逐项一致。近奇异 fold 的坐标末位差异来自重复求解的收敛路径，不影响 `fold/smooth/tail` 判定。
- **最终范围**：加密扫描支持 consistent 主 fold 链及其 `Ro_t≈-0.470367` fold--tail 端点；N=120--280 不支持低 `|Ro|` 的第二 fold 岛。traditional 的 `Ro≈-.99` 转换仍为 unresolved，不能给出精确临界值。
- **证据边界**：以上是切片 continuation 与网格敏感性证据，不等同于二维参数平面完整拓扑定理；左侧 smooth--fold 边界和 traditional 窄转换仍需二维 fold continuation。

### 2026-09-08 — L043：代表性分支 Tecplot 数据导出

- **状态**：完成；仅转换 L036/L040 的既有 continuation 数据，没有重算或改变任何科学结论。
- **内容**：输出 `Ro=-1` 和 `Ro=-.5` 的 traditional/consistent 代表性分支，其中包括 fixed-$T_w$ continuation 段，以及可用的跨 fold pseudo-arclength 段。字段为 `point, Ro, Tw, Hinf, residual, thermal_length, segment_code`。
- **用途与边界**：用于在 Tecplot 中绘制 $H_\infty$--$T_w$ 分支图；`segment_code=1` 的 arc 文件原始数据未保存 thermal-length 诊断，导出为 `NaN`。所有有限-$T_w$分支及其分类仍须按照 L039--L042 的网格审计范围解读。
- **输出**：`work/scripts/export_representative_branch_tecplot_20260908.jl`；`work/results/presentation_evidence_audit_20260908/tecplot/representative_traditional_consistent_branches.dat`。

### 2026-09-09 — L044：无 fold generic-tail 极限的先导计算

- **状态**：pilot 完成；结果否定了“普通 tail 具有非零 $D_0$ 或非零 unscaled core”的首个工作假设，正式高网格扫描前必须更换辨识指标。
- **模型与路径**：acceleration-consistent corrected-energy 基本流，`Pr=.72,gamma=1,N=120`；固定 `Ro=-.469,-.40`，只解 `R=0`，以 `Hinf=-epsilon` 释放 `Tw`，`epsilon=.08,.03,.01`。全过程未调用 fold 或零模方程。
- **结果**：两条 generic-tail 切片的固定-$\eta$速度场均趋零；`Ro=-.40` 上 $G,F,H$主要呈 $O(\epsilon)$，而不是 fold-limit 的 $G=O(\epsilon),F,H=O(\epsilon^2)$。代数退化量 $D=\chi(T_w)s^2-\omega^2$ 在 `Ro=-.40` 从 `1.0855e-2` 降至 `1.1357e-3`，也趋向零；故 $D_0=0$ 和未缩放 $M_0=0$ 很可能是 regular tail 的共同性质，不能单独选出 $Ro_t$。
- **新机制假设**：特殊 junction 可能对应 generic-tail 中 $O(\epsilon)$ meridional inner core 系数消失，使 $F,H$的首个非零项从 $O(\epsilon)$提升为 $O(\epsilon^2)$。后续应外推 $D/\epsilon$ 及基于 $(F,H)/\epsilon$ 的 core amplitude，而不是只画 $D_0,M_0$。
- **数值警告**：collocation 方程的壁面行被边界条件替换，直接用谱矩阵计算 $F''(0)$在长尾欠分辨时不可靠；pilot 的 `D_from_wall` 在小 $\epsilon$未通过恒等式检查。主诊断必须使用参数代数式 $D=\chi(T_w)s^2-\omega^2$，壁面二阶导数只作升网格后的交叉检查。
- **输出**：`work/scripts/generic_tail_zero_core_scan_20260909.jl`；`work/results/generic_tail_zero_core_20260909/N120_a2p0000_b0p6000_c0p5000/`。

### 2026-09-09 — L045：无 fold generic-tail 一阶 meridional-core 退化扫描

- **状态**：主扫描、临界加密、网格与映射检查完成；支持一阶标度退化机制，`h0` plateau 仍未完全收敛。
- **模型与范围**：consistent corrected-energy basic-flow，仅解 `R=0`，固定 `Hinf=-epsilon`、释放 `Tw`，未调用 fold/零模方程。主扫描含 9 个 `Ro=-.469...-.30`、6 个 `epsilon=.08...01`、`N=240,320`；临界加密加入 `Ro=-.4700,-.4702`，用 `N=240,280,320,360`；`N=280` 比较三组映射。
- **主要结果**：普通 tail 的未缩放 core 与 `D` 均趋零，但 `D/epsilon`、`H/epsilon` 与 `Mmer/epsilon` 外推为非零；它们向 `Ro_t` 系统消失。近临界 `D1` 零点在三组 N280 映射上为 `[-.47038077,-.47036374]`，包含独立 outer 值 `-.4703674163`；综合报告为 `Ro_*=-.47038 +/- 3e-5`。N320 最接近的二次窗口结果为 `-.47036350`。
- **解释**：证据支持 generic tail 的 `F,H=O(epsilon)` 一阶 meridional core 在 `Ro_*` 消失，从而提升为临界类的 `F,H=O(epsilon^2)`；这比“nonzero core 到 zero core”更准确，并且完全独立于 fold 约束。
- **限制**：固定 `eta=1,2,5` 的 `H1` 零点接近同一区间，但 `eta=8` 仍受有限 epsilon 的 overlap 污染，故 `h0=H1(infinity)` 尚未独立高精度收敛；谱壁面 `F''(0)`未通过小 epsilon 交叉检查，主指标使用精确代数 `D`。
- **输出**：`docs/GENERIC_TAIL_FIRST_ORDER_DEGENERACY_2026-09-09.md`；`work/scripts/generic_tail_zero_core_scan_20260909.jl`；`work/scripts/analyze_generic_tail_zero_core_20260909.jl`；`work/results/generic_tail_zero_core_20260909/`。

### 2026-09-09 — L046：generic outer family 与 h0=0 成员

- **状态**：完成并通过步长/坐标检查；属于既有 leading outer 方程的结构推广，不作为独立方程级验证重复计数。
- **方法**：固定 `Ro`，从 `g=1,p=0,h=-1` 反向积分到 `g=0`，直接得到 `p0(Ro),h0(Ro)`；再以 `Y` 坐标正向积分检查远场。没有施加 `h0=0` 或 fold 条件。
- **结果**：generic family 连续覆盖 `Ro=-.50...-.30`，`h0` 从正值穿过零后变负；`h0=0` 在 `g_steps=250--16000` 上稳定为 `Ro=-.470367416464`、`Tw=1.661587419719`、`p0=.654193538532`。这精确重现 critical outer 点，并证明它是 generic family 的零入口输运成员。
- **解释边界**：该零点与旧 outer 值的一致是同一 ODE 的等价参数化，不是 independent basic-flow/fold evidence；新贡献是展示连续 family 及 matching freedom 的消失。
- **输出**：`work/scripts/generic_outer_family_20260909.jl`；`work/results/generic_outer_family_20260909/`。

### 2026-09-09 — L047：generic O(epsilon) inner BVP、剖面验证与局部退化律

- **状态**：完成；generic inner BVP 数值收敛，并由无-fold finite-epsilon profiles 支持。
- **方法**：在 outer `(p0,h0)` matching 数据下求解 generic first-order linear inner BVP；`L=12--24,deta=.01,.005`。将其与 `N=320` fixed-H basic-flow 的重标度剖面在 `eta=0:8` 比较，并对 `|Ro-Ro_t|<=.00164` 拟合局部系数。
- **结果**：`Ro=-.40` 时 inner `D1=.1110338`，与 finite-epsilon 外推 `.1110344` 一致；`epsilon=.03 -> .01` 时 `F1,H1,G1` relative-L2 error 分别约由 `21.3%,23.8%,5.87%` 降至 `7.4%,8.3%,1.97%`。局部律为 `h0~-7.78537 mu`、`D1~1.71942 mu`、`tau1~-1.35962 mu`、`Mmer,1~19.0614|mu|`，其中 `mu=Ro-Ro_t`。
- **结论**：generic tail 的 `O(epsilon)` meridional matching freedom 在 `Ro_t` 线性消失，critical `F,H=O(epsilon^2)`问题是其退化子类。靠近 `Ro_t` 时一阶目标趋零导致相对误差放大，尚不能用现有数据验证完整二阶 crossover profile。
- **输出**：`work/scripts/solve_generic_inner_first_order_20260909.jl`；`work/scripts/compare_generic_inner_finite_epsilon_20260909.jl`；`work/scripts/analyze_generic_outer_inner_degeneracy_20260909.jl`；`work/results/generic_inner_first_order_20260909/`；`work/results/generic_outer_inner_degeneracy_20260909/`。

### 2026-09-09 — L048：无 fold 的双参数 crossover 验证与 fold 后验比较

- **状态**：固定 $\Lambda=(Ro-Ro_t)/\epsilon$ 的基本流扫描完成；支持 $\mu=\epsilon\Lambda$ 为 generic 与 critical tail 的 distinguished crossover scaling。
- **模型与参数**：consistent corrected-energy basic flow，只解 `R=0`；规定 $H_\infty=-\epsilon$、释放 $T_w$，未使用 fold 或零模条件。取 $\Lambda=(-.4,-.3,0,.5,1,2)$、$\epsilon=(.03,.02,.015,.01)$、$N=240,320$，映射 $(a,b,c)=(2,.6,.5)$。
- **缩放结果**：在固定 $\eta=.5,1,2,5,8$ 上，$F/\epsilon^2,H/\epsilon^2,G/\epsilon,(T-T_t)/\epsilon$ 均随 $\epsilon$ 减小而收敛。N320 相对 $\epsilon=.01$ 的最坏 $F/\epsilon^2$ 差异由 `.03` 的 13.0% 降至 `.015` 的 3.55%，$H/\epsilon^2$ 由 19.2% 降至 5.34%；`.01` 上 N240--N320 的最坏差异分别仅 0.211% 与 0.0629%。最大非线性残差 $3.2\times10^{-9}$。
- **独立系数检查**：N320 的 $D/\epsilon^2=a+b\Lambda$ 全区间拟合斜率从 $b=1.6481$（$\epsilon=.03$）收敛至 $1.7054$（`.01`），逼近 generic inner 独立得到的 $c_D=1.719419$；拟合 RMS 从 0.0404 降至 0.00972。截距趋向约 0.255，对应 critical 二阶贡献。
- **fold 后验比较**：在完成无 fold 验证后，用修正 fold 链计算 $\Lambda_f$。N320 在 $\epsilon=.1,.05,.02,.01$ 分别为 `-.339119,-.316586,-.302897,-.298372`，证明已解析范围内 fold 进入 $\Lambda=O(1)$ 层；`.005` 的 N280/N320 值 `-.3654/-.3442` 未网格收敛，不用于极限外推。
- **结论边界**：已数值支持“first-order matching degeneracy 产生 distinguished crossover，finite-$\epsilon$ fold 进入该层”。尚未完成二阶 inner 与一阶 outer 修正的完整耦合 BVP，也未推导选择 fold 层内轨迹或极限 $\Lambda_f$ 的下一阶可解性条件。
- **输出**：`docs/GENERIC_TAIL_CROSSOVER_STUDY_2026-09-09.md`；`work/scripts/generic_tail_crossover_20260909.jl`；`work/scripts/analyze_generic_tail_crossover_20260909.jl`；`work/results/generic_tail_crossover_20260909/`。

### 2026-09-09 — L049：crossover fold 零模的分量尺度

- **状态**：完成已有修正 fold 零模的 crossover 重标度检查；为下一阶 Fredholm 问题确定了非一致分量尺度。
- **数据依赖**：读取 `presentation_evidence_audit_20260908/tangent/N200--320/tangent_table.csv`，不重解 fold；检查 $\epsilon=.02,.015,.01$。
- **结果**：跨 N240/280/320，$\max|v_H|/(\epsilon\max|v_G|)$ 与 $\max|v_F|/(\epsilon^2\max|v_G|)$ 收敛为有限非零量。N320 对 $\epsilon\to0$ 的线性外推分别为 `3.5766`、`4.6993`；$\max|v_T|/\max|v_G|\to0.7963$。因此 fold mode 满足 $v_G,v_T=O(1),v_H=O(\epsilon),v_F=O(\epsilon^2)$。
- **方向估计边界**：N320 三点切向外推给 $dRo/d\epsilon\to-0.3070$（RMS 0.0021），与有限 $\Lambda_f\approx-0.30$ 一致；但 N240 的 `.01` 切向为异常点，且坐标商与导数外推仍有有限-$\epsilon$差异，故不把该数值当成解析选择结果。
- **理论影响**：下一阶条件应写成 $\mathcal L_0v_1+(\mathcal L_1^{(0)}+\Lambda\mathcal L_1^{(\mu)})v_0=0$，由匹配伴随模给出 $\mathcal S_0+\Lambda\mathcal S_1=0$。在构造匹配权重前，不能用有限 collocation 向量的未加权点积冒充连续 Fredholm 内积。
- **输出**：`docs/CROSSOVER_FOLD_SOLVABILITY_PLAN_2026-09-09.md`；`work/scripts/analyze_crossover_fold_mode_20260909.jl`；`work/results/crossover_fold_solvability_20260909/`。

### 2026-09-09 — L050：连续 outer 零模算子、Green concomitant 与 graded-inner 警告

- **状态**：完成 matched-adjoint 推导的第一阶段；从连续方程得到 critical outer 的线性化 DAE、形式伴随与边界 concomitant，尚未计算 Fredholm 商。
- **outer 结果**：对 $v_H=\epsilon a,v_F=\epsilon^2b,v_G=c,v_T=d$，连续线性化包含连续性、线性化慢流形、方位和能量四式。逐项分部积分得到四个伴随方程及 $\mathscr B_o=\ell a+\sigma c'-\sigma'c+\chi Ro_th\sigma c+\tau d'-\tau'd+\Pran Ro_th\tau d$。
- **接口含义**：outer 远场可由 primal 自由度推出分离伴随条件；但 $Y\downarrow0$ 的 $c',d'$ 是 matching data，不能把 outer concomitant 单独置零，必须由 inner overlap 项抵消。
- **关键尺度修正**：全域 max-norm 尺度由 outer 主导。在固定 inner 坐标下实际为 $v_G,v_T=O(\epsilon)$、$v_H,v_F=O(\epsilon^2)$；因此仅缩放未知量会使 leading inner operator 丢失径向黏性行。正确对象必须是同时缩放变量列和残差行的 graded inner--outer operator。
- **结论边界**：本条建立了连续 outer Green identity，并识别出 naive leading-inner adjoint 的退化原因；尚未推导保留下一阶径向行的 inner concomitant，也未验证 interface cancellation、overlap invariance 或计算 $\mathcal S_0,\mathcal S_1$。
- **输出**：`docs/CROSSOVER_FOLD_SOLVABILITY_PLAN_2026-09-09.md`。

### 2026-09-09 — L051：保留径向黏性行的 graded inner 零模与伴随

- **状态**：完成 crossover inner 零模至 $O(\epsilon^2)$ 的连续分级展开，并得到第二级形式伴随和 concomitant；尚未完成 inner--outer 接口抵消。
- **primal hierarchy**：第一层为 $c''=d''=0$ 与径向代数关系；第二层同时出现 $a'+2b=0$、$b''+A_0c_2+C_0d_2=-A_1c-C_1d$、$c_2''+2\chi_ts_tb=0$、$d_2''=0$。其中 $A_1,C_1$ 通过 $q_1=s_t'\Lambda+Ro_tg_1$ 和 $\theta_1=\mathcal T'(Ro_t)\Lambda+B_T\eta$ 对 $\Lambda$ 呈仿射依赖。
- **adjoint 结果**：对第二层四式逐项分部积分，得到伴随系统 $-L'=0$、$2L+P''+2\chi_ts_tS=0$、$A_0P+S''=0$、$C_0P+Q''=0$，以及 $\mathscr B_i=La+Pb'-P'b+Sc_2'-S'c_2+Qd_2'-Q'd_2$。
- **选择项位置**：inner 候选贡献为 $-\int P(A_1c+C_1d)d\eta$，自然分成常数项和 $\Lambda$ 系数；但它不是完整 $\mathcal S_0,\mathcal S_1$，还必须加入 outer correction 和 overlap concomitant。
- **关键警告**：$c_2,d_2$ 在 overlap 含二次多项式共同部分。若不先做 common-part subtraction，所得积分会依赖任意接口位置，不能称为 matched Fredholm condition。
- **输出**：`docs/CROSSOVER_FOLD_SOLVABILITY_PLAN_2026-09-09.md`。

### 2026-09-09 — L052：同一 fold 零模的 inner/outer 双重尺度验证

- **状态**：完成现有修正 fold-mode profiles 的固定-$\eta$与固定-$Y$重采样；数值确认全域 max-norm scaling 不能直接作为 inner scaling。
- **inner 结果**：N320、固定 $\eta=1$，当 $\epsilon=.02\to.01$ 时，$v_G/(\epsilon\max|v_G|)$ 从 `2.0267` 到 `2.0403`，$v_T/(\epsilon\max|v_G|)$ 从 `-.5883` 到 `-.5848`，$v_H/(\epsilon^2\max|v_G|)$ 从 `-2.5556` 到 `-2.5829`，$v_F/(\epsilon^2\max|v_G|)$ 从 `2.4329` 到 `2.4561`。
- **outer 结果**：固定 $Y=.5$ 时，$v_H/(\epsilon\max|v_G|)$ 从 `-2.6869` 到 `-2.6891`，其余 outer-scaled 分量同样收敛；在较小 $Y=.1$，$v_F$ 收敛较慢，符合 overlap 污染预期。
- **结论**：同一全局零模在 outer 满足 $(v_H,v_F,v_G,v_T)=(O(\epsilon),O(\epsilon^2),O(1),O(1))$，在固定 inner 坐标则满足 $(O(\epsilon^2),O(\epsilon^2),O(\epsilon),O(\epsilon))$。这直接支持分区 graded operator 和接口 concomitant 的必要性。
- **输出**：`work/scripts/analyze_crossover_fold_mode_20260909.jl`；`work/results/crossover_fold_solvability_20260909/inner_fixed_eta_scaling.csv`；`outer_fixed_Y_scaling.csv`。

### 2026-09-09 — L053：first/second-grade 联合 inner adjoint block

- **状态**：将 inner 首层代数/线性增长方程与次层径向黏性方程合并为单一 graded block，修正只对次层子块取伴随的不完整性。
- **结果**：引入七个方程乘子 $(U,V,R,L,P,S,Q)$ 后，完整 inner concomitant 为 $Uc'-U'c+Vd'-V'd+La+Pb'-P'b+Sc_2'-S'c_2+Qd_2'-Q'd_2$；伴随方程显式包含 $A_1P,C_1P$，因此 $\Lambda$ 通过连续 inner operator 的系数进入。
- **意义**：fold 方向可等价地理解为使 graded inner--outer null BVP 存在非平凡解的参数，或理解为对不含 $\Lambda$ 的基准算子施加 Fredholm 条件；两种写法必须给同一 $\mathcal S_0+\Lambda\mathcal S_1$。
- **未完成项**：尚需推导 outer 下一阶零模方程、adjoint matching rules 与 common-part subtraction，才能验证接口 concomitant 抵消。
- **输出**：`docs/CROSSOVER_FOLD_SOLVABILITY_PLAN_2026-09-09.md`。

### 2026-09-09 — L054：outer 次阶零模方程与 primal common part

- **状态**：完成 outer $O(\epsilon)$ 零模修正的连续 Fréchet 表达，并列出进入 inner--outer concomitant 的 primal 共同部分；尚未反推出 adjoint matching rules。
- **outer forcing**：以规范化 DAE $\mathcal N(X;r)=0$ 表示，分解 $X_1=X_1^{(0)}+\Lambda X_1^{(\mu)}$，得到 $\mathcal L_o^{(0)}q_1=-\mathcal N_{XX}[X_1^{(0)},q_0]-\Lambda\{\mathcal N_{XX}[X_1^{(\mu)},q_0]+\mathcal N_{Xr}q_0\}$。
- **共同部分**：leading outer null mode 在入口满足 $a_o\sim A_1Y$、$b_o\to-A_1/2$、$c_o\sim C_1Y$、$d_o\sim D_1Y$，且 $2\chi_ts_tRo_tC_1-\gamma s_t^2D_1=0$。次阶 inner 的 $c_2,d_2$ 含与 outer Taylor 展开相同的二次和一次多项式。
- **匹配含义**：齐次 wall data 与 $c''=d''=0$ 强制 $c_{o1}(0)=d_{o1}(0)=0$。其余 adjoint 接口条件必须通过 inner/outer concomitant 对独立共同系数逐项抵消得到，不能在 $Y=0$ 人为补条件。
- **结论边界**：已把 outer forcing 分成常数与 $\Lambda$ 两部分，但尚未求 $X_1^{(0)},X_1^{(\mu)}$ 的完整剖面，也未计算 $\mathcal S_j^{(o)}$。
- **输出**：`docs/CROSSOVER_FOLD_SOLVABILITY_PLAN_2026-09-09.md`。

### 2026-09-09 — L055：四场连续伴随、区域权重与接口条件来源

- **canonical trace**：把 Green concomitant 改写为 $(v_H,v_F,v_F',v_G,v_G',v_T,v_T')$ 的共轭迹形式，明确匹配要求的是七个物理伴随迹的等式。它们必须随后在 overlap 中按 $\epsilon$ 展开；特别是不能把 inner 的 $P$ 直接等同于 outer 的 $\rho$。
- **状态**：完成 matched-adjoint 的连续理论锚定；未进行新的数值计算，未计算 Fredholm 系数或 fold 方向。
- **方法**：直接对 consistent corrected-energy 的四个连续残差行线性化，并在物理无权内积下逐项分部积分。该推导与 work/src/BEKConsistent.jl 的连续残差和 Jacobian 项逐项对应；不使用离散左零向量作为伴随输入。
- **结果**：得到四场形式伴随 (alpha, beta, zeta, tau) 与精确 Green concomitant。outer power counting 为 O(1, epsilon^2, 1, 1)：径向伴随 beta 的 epsilon^2 因子由径向代数项和 azimuthal-adjoint 二阶导数的平衡强制。fixed-inner 坐标下，四个物理伴随分量均为 O(1)，其首个非平凡系数给出既有 second-grade block 的 (L, P, S, Q)。seven-multiplier block 的 (U, V, R) 只能解释为分级行加权后的 correction/constraint bookkeeping，而不是额外物理伴随场。
- **接口结论**：在 1 << eta_m << epsilon^(-1)、Y_m=epsilon eta_m，正确条件为 B_phys^in(eta_m)-B_phys^out(Y_m) -> 0。它必须在 primal common part 的独立系数上逐项实现；不得在 outer DAE 的 Y=0 人为添加 sigma(0)=0、vartheta(0)=0 一类条件。outer concomitant 是物理边界型的 O(epsilon) 部分，径向项在该层为 O(epsilon^5)，说明其在 leading outer Green identity 中弱但未消失。
- **结论边界**：上述尺度来自连续方程的 power counting，尚未由独立离散左奇异向量测量。尚未解出 matched four-field adjoint BVP，尚未验证 interface cancellation、cutoff invariance、normalization invariance 或 grid/mapping convergence；因此不能把有限-epsilon 约 -0.30 的切向趋势提升为解析选择值。
- **下一步**：以 A1,C1,D1,A2,C2,D2 等共同部分逐项抵消反推出伴随入口匹配关系，再构造带连续 quadrature 权重的 inner--outer adjoint BVP。
- **输出**：docs/CROSSOVER_FOLD_SOLVABILITY_PLAN_2026-09-09.md。

### 2026-09-09 — L056：overlap 伴随迹的首阶关系与径向未闭合项

- **状态**：完成连续 physical-trace matching 的首阶推导；未解伴随 BVP，未计算 Fredholm 系数。
- **方法**：将精确 Green concomitant 规范写为 seven canonical primal traces $(v_H,v_F,v_F',v_G,v_G',v_T,v_T')$ 与其共轭伴随迹的双线性型；在 regular outer entrance 与 bounded leading inner $G,T$-adjoint 的条件下比较 overlap 首项。
- **已得接口关系**：$L_\infty=\ell(0)$、$S_\infty=\sigma(0)$、$Q_\infty=\vartheta(0)$、$S'_\infty=Q'_\infty=P_\infty=0$，并且 $L_\infty=-\chi_ts_tS_\infty$。这组关系由 physical trace cancellation 得到，不是把 inner/outer 方程变量按名字对应。
- **关键发现**：leading inner radial adjoint 满足 $P''''+4\chi_t^2s_t^2P=0$。本问题 $\chi_ts_t>0$，所以 bounded $P$ 在 large-inner limit 指数衰减；而 outer physical radial trace 为 $\beta_o=\epsilon^2\rho(0)+\cdots$。因此 $P$ 的首阶尾部不能在 algebraic overlap 直接匹配 $\rho(0)$。
- **结论边界**：这不证明 $\rho(0)=0$，也不单独证明存在第三层。它严格表明：要么更高阶 inner adjoint coefficient 提供该 $O(\epsilon^2)$ trace，要么需要 intermediate radial-adjoint transition layer；在分辨这两种可能性前，$\rho$ 的 entrance condition 和 $\mathcal S_0,\mathcal S_1$ 都未定义。
- **下一步**：将 exact physical $\beta$ 与 $\Pi_F$ traces 推到首个 algebraic overlap 阶，判定 higher-grade matching 与 intermediate-layer 机制。
- **输出**：docs/CROSSOVER_FOLD_SOLVABILITY_PLAN_2026-09-09.md。

### 2026-09-09 — L057：二阶径向伴随迹的局部匹配兼容性

- **状态**：完成 exact four-field adjoint 的内区 $O(\epsilon^2)$ 展开；这是理论推导，未进行新的数值计算，尚未解 matched-adjoint BVP 或计算 Fredholm 系数。
- **归一化修正**：outer DAE 的径向残差取为 $\mathscr A=\chi q^2-\omega^2$，而物理径向残差的 leading term 为 $R_F=-\mathscr A/Ro_t$。因此物理伴随迹应写成 $\beta=-\epsilon^2Ro_t\rho+O(\epsilon^3)$，不是 $\epsilon^2\rho$；相应地 $\Pi_{F,o}=\epsilon^3\{Ro_t\rho_Y-\chi Ro_t^2h\rho\}+\cdots$。
- **二阶方程**：定义 $q_2=s_2+Ro_tg_2+\Lambda g_1$、$A_2=-2(\chi_tq_2+\chi_1q_1+\chi_2s_t)$ 与 $C_2$ 为 $-\gamma A_r$ 的二阶系数，得到四场伴随的 $P_2,S_2,Q_2,L_2$ 方程。保留 $h_2,f_2$ 项后，$G$-伴随方程在大-$\eta$ common part 给出 $S_2''+A_0P_2\to\chi_tRo_t(h_1-2f_0)S_\infty$。
- **兼容性检查**：用 critical base matching $h_1+2f_0=0$ 及 $S_\infty=\sigma(0)$，该共同部分与 outer $G$-adjoint entrance equation 一致；物理径向迹要求 $P_2\to-Ro_t\rho(0)$。因此首个 algebraic overlap 检查不强制 $\rho(0)=0$，也不强制引入第三个 radial-adjoint layer。
- **边界**：这只是局部 common-part compatibility，不是全局存在性证明。仍需推导/求解 $P_3$ 及完整 inner--outer adjoint BVP，并检验 concomitant cancellation、interface-cutoff invariance、normalization/grid/mapping independence，之后才可定义 $\mathcal S_0,\mathcal S_1$ 和 $\Lambda_f$。
- **输出**：docs/CROSSOVER_FOLD_SOLVABILITY_PLAN_2026-09-09.md。

### 2026-09-09 — L058：七个物理共轭迹的 overlap 共同部分账本

- **状态**：完成连续 exact four-field Green concomitant 的七个物理共轭迹 \((\alpha,\beta,\Pi_F,\zeta,\Pi_G,\tau,\Pi_T)\) 在 critical crossover overlap 中的首个非零共同部分展开；这是解析一致性检查，未进行新的数值计算或 Fredholm 系数计算。
- **接口原则**：接口条件由 \(\mathscr B_{\rm phys}^{\rm in}-\mathscr B_{\rm phys}^{\rm out}\to0\) 对独立 primal canonical traces 的逐项消去产生，而不是按 inner/outer 乘子变量名称相等。导出的首阶关系为 \(L\to\ell_0\)、\(S\to\sigma_0\)、\(Q\to\vartheta_0\)，以及 \(P,P_1\to0\)、\(P_2\to-Ro_t\rho_0\)、\(P_2'\to0\)、\(P_3'\to-Ro_t\rho_{Y0}\)。
- **局部交叉验证**：\(S_1,Q_1\) 的线性共同部分和 \(S_2,Q_2\) 的二次共同部分逐项复现 outer \(G,T\) conjugate derivative traces。inner $G$、$T$ adjoint equations 分别复现 outer entrance 的 \(\sigma_{YY0}\)、\(\vartheta_{YY0}\) relations；L057 的 $P_2$ 径向迹兼容性因而保持成立。此前 Julia sanity-check 的失败来自 `2f0` 被 Julia 解析为 Float32 字面量而非 `2*f0`；改为显式乘法后差为 \(1.1\times10^{-16}\)，不构成推导修正。
- **结论边界**：这证明到所列 algebraic common-part 阶数的局部 Green-concomitant cancellation；尚未构造或求解全局 matched adjoint BVP，\(\rho_0,\rho_{Y0}\) 和 far-field data 尚未被选择，亦未证明 interface-cutoff independence。因此不得据此计算 \(\mathcal S_0,\mathcal S_1\) 或宣称得到解析 \(\Lambda_f\)。
- **下一步**：以该七迹 ledger 为 interface map，构造完整 inner--outer adjoint BVP；在多组 movable overlap cuts 上检验 \(\mathscr B_i-\mathscr B_o\to0\)、\(\mathcal S_j^{\rm in}+\mathcal S_j^{\rm out}\) 的 cutoff independence，以及 normalization/grid/mapping independence。
- **输出**：docs/CROSSOVER_FOLD_SOLVABILITY_PLAN_2026-09-09.md。

### 2026-09-09 — L059：径向导数迹的三阶进入位置

- **状态**：从 exact $G$ 与 $F$ physical-adjoint equations 写出 $O(\epsilon^3)$ 的 $S_3,P_3,L_3$ 项，并将 outer radial conjugate derivative trace 的首个非零入口阶与其配对；未求解该三阶系统。
- **结果**：$P_2$ 提供 outer radial value trace，而 $\Pi_F$ 的入口数据到 $O(\epsilon^3)$ 才出现，给出 $P_3'\to-Ro_t\rho_{Y0}$。因此 $\rho_{Y0}$ 不应被人为设为零；它对应 $P_3\sim-Ro_t\rho_{Y0}\eta+\mathrm{constant}$ 的线性共同部分。
- **结论边界**：三阶方程已定位 radial derivative matching 所在的首阶，但尚未检查 $S_3,P_3$ 是否确实复现同一多项式共同部分，更未证明 global BVP existence 或 cutoff independence。该项不能用于计算 $\mathcal S_0,\mathcal S_1$ 或 $\Lambda_f$。
- **下一步**：将三阶 $G$-adjoint common part 与 outer $\rho_Y$ 入口展开逐项比较；若局部兼容，再实施完整 matched-adjoint BVP 与 movable-cut tests。
- **输出**：docs/CROSSOVER_FOLD_SOLVABILITY_PLAN_2026-09-09.md。

### 2026-09-09 — L060：\(\rho_Y(0)\) 的显式三阶入口检验式

- **状态**：对 leading outer $G$-adjoint equation 在 \(Y=0\) 求一次导数，得到含 \(\rho_{Y0}\) 的显式入口关系；未进行数值 BVP 求解。
- **结果**：令 \(K=\chi q\)、\(M=\chi Ro_th\)、\(N=\chi Ro_tf\)，则 \(\sigma_{YYY0}+2Ro_t(K_{Y0}\rho_0+K_0\rho_{Y0})-(M_{YY0}\sigma_0+2M_{Y0}\sigma_{Y0})+2(N_{Y0}\sigma_0+N_0\sigma_{Y0})=0\)。它是 inner \(S_3''\) 与 \(P_3\sim-Ro_t\rho_{Y0}\eta+\mathrm{constant}\) 的线性 common part 必须满足的 outer target。
- **结论边界**：该式只定义了三阶局部匹配检验，尚未把所有 inner 三阶系数代入化简，未验证三阶 cancellation，亦不选择 \(\rho_{Y0}\)、\(\mathcal S_0,\mathcal S_1\) 或 \(\Lambda_f\)。
- **下一步**：完整展开 \(S_3\) 方程的线性 \(\eta\) 系数并对照该式；随后构造 global matched-adjoint BVP，执行 movable-cut、normalization、grid 和 mapping independence tests。
- **输出**：docs/CROSSOVER_FOLD_SOLVABILITY_PLAN_2026-09-09.md。

### 2026-09-09 — L061：非相似 heated-BEK 边界层 PDE 与 similarity Gate

- **状态**：完成轴对称 steady acceleration-consistent boundary-layer PDE 的有量纲和 BEK 无量纲推导；新增 Julia algebraic regression。没有求解二维 PDE，也没有产生 fold/tail 拓扑结论。
- **模型**：取 \(u_R=-\Omega RoRF(x,\eta)\)、\(u_\phi=\Omega R[s+RoG(x,\eta)]\)、\(u_z=-\ell\Omega RoH(x,\eta)\)，温度为标量 \(T(x,\eta)\)。第一版保留 \(xFF_x,xFG_x,xFT_x\)、壁法向对流和 \(F_{\eta\eta},G_{\eta\eta},T_{\eta\eta}\)，忽略径向扩散。consistent 密度因子从二维 material acceleration 重新施加，不是向 ODE 人工添加导数。
- **Gate 1**：令全部径向导数为零后严格恢复 work/src/BEKConsistent.jl 的四个 similarity residuals。随机局部状态回归覆盖 \(Ro=-0.50,-0.4703674,-0.40\)，每点 1000 组，最大差 \(4.441\times10^{-16}\)。
- **关键逻辑边界**：无限圆盘、均匀壁温和均匀刚体远场下，similarity 解是 PDE 的精确不变子空间；从精确 similarity profile 启动的径向 march 只能作为代码回归，不能检验 radial development 是否 regularize tail。真正 non-similar 计算还需有限圆盘/edge、\(T_w(x)\) 或非相似入口数据。
- **输出**：docs/NONSIMILAR_HEATED_BEK_BOUNDARY_LAYER_2026-09-09.md；work/scripts/verify_nonsimilar_similarity_gate_20260909.jl。

### 2026-09-09 — L062：三个关键 \(Ro\) 的径向推进方向审计

- **状态**：使用既有无-fold generic-tail 一阶 inner profiles 做 marching-admissibility diagnostic；未重新求解基本流。源数据为 work/results/generic_inner_first_order_20260909 下三个对应 profile。
- **结果**：在 \(0<\eta\le8\)，物理径向速度首阶系数 \(-RoF_1\) 在 \(Ro=-0.50\) 和 \(-0.40\) 都换号，零点约为 \(\eta=5.3315\) 和 \(4.6177\)。在 \(Ro=-0.470367\)，整体振幅仅 \(1.23\times10^{-10}\)，应解释为 first-order degeneracy，而不是可靠方向。
- **结论**：两个 generic benchmarks 在不同 wall-normal 区域具有相反径向传播方向；临界点的一阶 transport 又消失。单向抛物 marching 不适合作为三点物理验证的主模型。有限圆盘仍是首选 symmetry breaking，但必须用包含 radial diffusion/edge region 的 global axisymmetric model；抛物模型只适合作 regression 或受控 \(T_w(x)\) 响应实验。
- **证据边界**：换号数据属于 generic-tail first-order asymptotics，不等于 finite-\(\epsilon\) 全场特征分析，也尚未证明任何 finite-disk fold/tail 结果。
- **输出**：work/scripts/audit_similarity_radial_direction_20260909.jl；work/results/nonsimilar_heated_bek_20260909/radial_direction_audit.csv；docs/NONSIMILAR_HEATED_BEK_BOUNDARY_LAYER_2026-09-09.md。

### 2026-09-09 — L063：similarity 子空间的径向横向响应 pencil

- **状态**：完成 steady non-similar PDE 关于 similarity base 的连续 Mellin-mode 线性化，并实现共享 Julia pencil assembly；这是空间径向响应，不是 temporal 或 three-dimensional stability。
- **理论结果**：取扰动 \(e^{\lambda\log x}\widehat Q(\eta)\)，得到 \((J_{\rm sim}+\lambda M_r)\widehat Q=0\)。\(M_r\) 仅有 continuity-to-\(f\)、radial-to-\(f\)、azimuthal-to-\(g\)、energy-to-\(\theta\) 四个 interior blocks，系数分别为 \(1,\chi RoF_0,\chi RoF_0,PrRoF_0\)，boundary rows 为零。故 \(\lambda=0\) 严格恢复 fixed-\(Ro,T_w\) similarity Jacobian；simple similarity fold 必然对应该 reduced spatial model 的 zero-exponent singularity。
- **验证**：在 \(Ro=-0.50,-0.4703674,-0.40\) 的 degree-60 isothermal consistent bases 上，对 \(\lambda=-1.3,-0.2,0,0.4,2\) 比较 nonlinear non-similar residual 的 centered directional difference 与 analytic pencil。15 组 relative-infinity error 为 \(2.23\times10^{-12}\) 至 \(6.55\times10^{-11}\)。
- **证据边界**：isothermal bases 只用于 operator regression；尚未计算任何 heated-base generalized eigenvalue、resolvent norm 或 transverse conclusion。由于 \(M_r\) 奇异，必须区分 finite 与 infinite generalized eigenvalues，并检查 direct/adjoint residual 和 grid/map convergence。
- **下一步**：在每个 \(Ro\) 上明确选择基本流坐标 \(T_w\) 或 \(\epsilon\)，先计算远离 fold/tail 的 regular reference spectrum，再沿分支逼近 fold 或 tail，跟踪最接近零的 finite \(\lambda\) 及 forced-\(T_w(x)\) resolvent。
- **输出**：docs/SIMILARITY_TRANSVERSE_RADIAL_RESPONSE_2026-09-09.md；work/src/BEKRadialResponse.jl；work/scripts/verify_radial_response_pencil_20260909.jl；work/results/nonsimilar_heated_bek_20260909/radial_pencil_directional_check.csv。

### 2026-09-09 — L064：similarity fold 的径向调制正规形

- **状态**：完成 fold-centered radial amplitude equation 的形式推导；尚未计算投影系数。
- **结果**：对 \(U=U_f+A(\log x)v+\cdots\) 投影得到 \(c_rA_\xi+c_T\Delta T_w+c_R\Delta Ro+c_2A^2=\mathcal F_{\rm ns}+\cdots\)，其中 \(c_r=\langle w,M_rv\rangle\)、\(c_2=\langle w,R_{UU}[v,v]\rangle/2\)。固定点恢复 similarity saddle-node normal form，而附近分支的慢空间指数为 \(\lambda_{\rm slow}=-2c_2A_*/c_r+\cdots\)。
- **物理解释**：\(c_2\) 给 similarity branch 的局部弯曲，\(c_r\) 将其转换成 \(\log r\) 上的径向吸引/排斥；若 \(c_r=0\)，一阶 parabolized response 再次退化，必须保留 radial diffusion 或更高阶调制。若 \(c_r\to0\) 与 tail limit 同时发生，则存在新的 radial-nonsimilarity/tail distinguished limit。
- **证据边界**：这是形式局部理论，不证明 full finite-disk stability，也不排除 disconnected non-similar steady branches。连续 pairing、boundary terms、coefficient convergence 和 forced-response comparison 尚未完成。
- **下一步**：先在远离 tail 的收敛 simple fold 上计算 \(c_r,c_T,c_2\)，用 pencil 的近零 \(\lambda\) 和小径向壁温强迫验证正规形；再沿 finite-\(\epsilon\) fold chain 检查 \(c_r\) 是否趋零。
- **输出**：docs/SIMILARITY_TRANSVERSE_RADIAL_RESPONSE_2026-09-09.md。

### 2026-09-09 — L065：fold–tail 链的径向横截性 Gate

- **状态**：完成 corrected-energy finite-ε fold 链上的径向 Fredholm 投影、N=240/280/320 网格比较，以及六个 fold 邻近分支的直接 pencil 特征值验证；尚未完成 rational-map 参数变化。
- **规范修正**：直接强制 (w^Tv=1) 时，离散 (w^Tv) 在部分 ε 与网格上接近零并换号，使裸 (c_r,c_2) 产生假性极点。正式比较改用 ĉ_r=⟨w,M_rv⟩/⟨w,R_Tw⟩ 和 ĉ_2=⟨w,R_UU[v,v]⟩/(2⟨w,R_Tw⟩)，并固定 max|v_G|=1；该比值对伴随归一化不敏感。
- **主要结果**：在 ε=0.10,0.05,0.02,0.01,0.005 上，N=240/280/320 的 ĉ_r 分别约为 -0.02447,-0.01180,-0.004597,-0.002278,-0.001132。用 ε≥0.01 的四点拟合，N=280 与 N=320 均给出 |ĉ_r|≈0.260 ε^1.030。因此相对于普通壁温 fold unfolding，一阶径向横截性随 fold–tail 极限近似线性消失。
- **正规形交叉验证**：在 N=200、ε=0.10,0.05,0.02 的 fold 两侧各求一个 fixed-Tw 非线性分支，并用 bordered eigenpair corrector 直接求解 (J+λM_r)q=0。六个直接慢空间特征值与 λ_NF=-2c_2A*/c_r 的相对差为 1.2%–3.0%，且两支之间按预测换号。
- **次要结果边界**：N=320、ε≥0.01 的 provisional fits 为 |ĉ_2|≈3.96 ε^2.145、|c_2/c_r|≈15.2 ε^1.114；ε=0.005 的 c_2 未收敛，故暂不把二阶指数作为最终结论。当前最稳健结论仅是 ĉ_r=O(ε)。
- **物理范围**：这是稳态轴对称径向空间响应，不是 temporal stability 或完整三维稳定性。结果表明 fold–tail 极限同时呈现 similarity fold 与一阶 radial-response 的相对退化，但尚不能确定径向扩散提升后的 distinguished scale。
- **下一步**：先做 rational-map sensitivity；随后把完整二维方程中径向扩散项投影到 fold 模，比较其 ε 阶数与 ĉ_rA_ξ，推导需要保留径向扩散的 distinguished radial scale。
- **输出**：work/scripts/analyze_fold_radial_transversality_20260909.jl；work/scripts/validate_fold_radial_normal_form_20260909.jl；work/results/fold_radial_transversality_20260909/；docs/SIMILARITY_TRANSVERSE_RADIAL_RESPONSE_2026-09-09.md。

### 2026-09-10 — L066：Gate R4 径向扩散提升及修正尺度

- **状态**：从完整轴对称圆柱坐标 scalar/vector Laplacian 推导 reduced BEK 四场系统可自洽保留的径向扩散项，实现 Julia frozen-radius operators，并在 corrected-energy fold 链上完成 N=280/320 投影与幂律拟合。尚未做 rational-map variation 或 variable-radius amplitude BVP。
- **模型修正**：对 u_R=RA、u_φ=RB，vector Laplacian 分别给 R A_RR+3A_R+RA_zz 与 R B_RR+3B_R+RB_zz；温度给 T_RR+R^{-1}T_R+T_zz。在 ξ=log x 下得到 F_ξξ+2F_ξ、G_ξξ+2G_ξ、T_ξξ。H 行是连续性而非轴向动量，故没有向该行错误加入 H radial diffusion；若需 W 黏性必须恢复轴向动量及压力修正。
- **算子**：固定 R0 后 local pencil 为 J+λM_r+κ0(λ²D2+λD1)，κ0=(ℓ/R0)²；D2 的 FF、GG、TT blocks 为 I，D1 的 FF、GG blocks 为 2I。五个 λ 的独立代数回归误差为零。
- **主要结果**：以壁温 fold transversality 归一化，ε=0.10,0.05,0.02,0.01 时 N=320 的 ĉ_rr 为 -48.432,-97.545,-243.619,-479.476；N=280 为 -48.430,-97.541,-244.687,-505.996。拟合得到 |ĉ_rr|≈4.64ε^{-1.017} (N=280) 与 4.91ε^{-0.996} (N=320)。柱坐标一阶几何投影同样约为 O(ε^{-1})。
- **机制分解**：N=320 下 ĉ_rr 的发散主要来自温度块，G 块提供部分抵消，F 块很小；因此这不是单纯速度 Laplacian 的几何效应，而与 thermal-tail 的长程响应直接相关。
- **关键修正**：原假设 ĉ_rr=O(1) 被数据否定。由于 ĉ_r=O(ε) 而 ĉ_rr,ĉ_r,geo=O(ε^{-1})，对 O(1) 的 ξ 变化，平衡要求 κ0=O(ε²)，从而候选物理半径为 R0/ℓ=O(ε^{-1})，不是 ε^{-1/2}。
- **证据边界**：这是 steady axisymmetric frozen-radius projection，不是 finite-disk solution、temporal stability 或 3D stability。R0/ℓ=O(ε^{-1}) 还依赖 A 在 ξ 上 O(1) 变化；其他 radial modulation scales 必须在下一阶多尺度方程中重新计数。精确 prefactor 仍需 mapping sensitivity。
- **下一步**：完成 rational-map sensitivity 后，推导保留 e^{-2ξ}(A_ξξ+geometry A_ξ) 的 variable-radius fold amplitude equation，并据此决定 finite-disk 域尺度与边缘条件。
- **输出**：docs/RADIAL_DIFFUSION_PROMOTION_2026-09-10.md；work/src/BEKRadialResponse.jl；work/scripts/verify_radial_diffusion_operator_20260910.jl；work/scripts/analyze_fold_radial_transversality_20260909.jl；work/results/fold_radial_transversality_20260909/。

### 2026-09-10 — L067：Gate R4 的独立 rational-map 收敛

- **状态**：在 N=240 上以两套新 rational maps 重新求解 corrected-energy fixed-H∞ fold 链并重新计算径向横截/扩散投影；不是对默认映射旧剖面的重采样。
- **参数与源依赖**：保持 a=2、c=0.5、Pr=0.72、γ=1，取 b=0.55 与 0.65；每套均计算 ε=0.10,0.05,0.02,0.01。求解器为 work/scripts/consistent_epsilon_fold_chain.jl，新输出分别位于 work/results/corrected_energy_epsilon_fold_chain/map_b055_N240 与 map_b065_N240。
- **结果**：b=0.55 给 |ĉ_r| exponent 1.03044、|ĉ_rr| exponent -1.00086、|ĉ_r,geo| exponent -0.94336；b=0.65 分别给 1.03095、-1.01281、-0.94629。默认 b=0.60 的 N=320 结果为 1.03079、-0.99608、-0.94184。
- **结论**：ĉ_r=O(ε) 与 ĉ_rr=O(ε^{-1}) 对本次 mapping variation 稳健，Gate R4 通过网格与首轮独立映射检查。由此得到的 O(1)-ξ 候选平衡 R0/ℓ=O(ε^{-1}) 可以进入下一步 reduced theory；仍不能把它当作已经验证的 finite-disk 尺度。
- **下一步**：推导 variable-radius fold amplitude equation，保留 κ(R)=ℓ²/R²、二阶径向导数和 cylindrical first-derivative contribution，并重新检查允许的 ξ/R 慢尺度。
- **输出**：docs/RADIAL_DIFFUSION_PROMOTION_2026-09-10.md；work/results/fold_radial_transversality_20260909/fold_radial_transversality_map_b055_N240.csv；fold_radial_transversality_map_b065_N240.csv；scaling_fits_map_b055_N240.csv；scaling_fits_map_b065_N240.csv。

### 2026-09-11 — L068：\(X,Y\) 二维 thermo-radial outer 的首阶闭合

- **状态**：从含 radial diffusion 的 non-similar BEK residual 直接实施 \(X=\epsilon R/\ell,\ Y=\epsilon z/\ell\) 展开，得到首阶二维 outer PDE；尚未求解该 PDE 或 nonlinear amplitude BVP。
- **场尺度**：\(F=\epsilon^2f(X,Y)\)、\(H=\epsilon h(X,Y)\)、\(G=g(X,Y)\)、\(T=t(X,Y)\)。连续性为 \(h_Y+2f+Xf_X=0\)，径向首阶仍为 \(\chi(s_t+Ro_tg)^2=\omega_t^2\)。
- **二维闭合**：方位方程含 \(g_{YY}+g_{XX}+3g_X/X\)，能量方程含 \(t_{YY}+t_{XX}+t_X/X\)。令全部 \(X\) 导数为零严格恢复已有一维 critical outer DAE。
- **全方程边界**：径向 meridional inertia、\(f\) 的双向黏性、完整轴向动量、\(W\)-黏性和小压力修正在下一 meridional 阶共同出现。没有向连续性方程错误加入 \(H\)-diffusion。
- **回归**：对 \(\epsilon=0.2,0.05,0.01\) 各 1000 个随机局部状态，原 radially diffusive residual 与缩放后表达式的最大代数恒等误差为 \(7.105\times10^{-15}\)。
- **投影解释**：outer measure \(d\eta=dY/\epsilon\) 为 \(O(1)\) thermal direct/adjoint modes 提供 \(\epsilon^{-1}\) 积分长度，与 Gate R4 中温度块主导 \(\widehat c_{rr}\) 的数值结果一致。
- **振幅方程修正**：在 \(X=\epsilon e^\xi\) 下，线性 radial operator 为 \(\epsilon[C_{rr}A_{XX}+(C_rX+(C_{rr}+C_g)/X)A_X]\)。由于当前 \(\widehat c_2=O(\epsilon^2)\) 仅为 provisional result，非线性项除以共同 \(O(\epsilon)\) 后仍含 \(\epsilon A^2\)；未选择 amplitude/detuning scaling 前，不能把 \(A^2\) 与 radial operator 直接宣称同阶。
- **结论边界**：已建立“1D near-wall inner ↔ 2D thermo-radial outer”的 leading structure，但尚未得到完整二维速度/压力重构、边界条件、存在性或 finite-disk topology。
- **下一步**：由 inner/overlap matching 推导二维 outer 在 \(Y=0\) 的数据和轴/外缘条件，并选择 nonlinear amplitude、forcing 与 parameter detuning 的一致尺度。
- **输出**：docs/TWO_DIMENSIONAL_THERMO_RADIAL_OUTER_2026-09-11.md；work/scripts/verify_thermo_radial_outer_scaling_20260911.jl。

### 2026-09-11 — L069：full-axisymmetric hierarchy 与物理可实现性审计

- **状态**：从完整 steady axisymmetric continuity、三动量和能量方程重做 \(X=\epsilon R/\ell,\ Y=\epsilon z/\ell\) 量级审计；本项修正 L068 对“首阶闭合”的表述。尚未求解 full hierarchy 或调用 compressible solver。
- **保留结果**：\(G,T\) 的二维扩散方程、连续性和 leading centrifugal constraint \(\chi q^2=\omega_t^2\) 均通过 full-equation audit，故 thermo-radial promotion 本身保留。
- **关键修正**：径向/轴向 meridional momentum 虽比 leading centrifugal balance 低四个 \(\epsilon\) 幂，但它们是 \(f,h\) 的首个非零动量方程。以 \(P=P_0/\epsilon^2+\epsilon^2P_4+\cdots\) 展开后，\(P_{4X},P_{4Y}\) 的可积性增加一个 meridional-vorticity constraint。因此 L068 的四个 thermo-azimuthal 方程只是必要子系统，不是完整 full-NS closure。
- **重力 Gate**：若重力沿转轴，忽略自然对流相对于弱 meridional balance 要求 \(g|\beta\Delta T|/(\Omega^2\ell)\ll\epsilon^3\)，比简单比较 \(g\) 与 \(\Omega^2R\) 更严格；若重力不沿轴，axisymmetry 本身失效。
- **Boussinesq Gate**：当前 \(\gamma=1,\ T_{w,t}=1.6615874197\) 对应 \(|\beta\Delta T|=\gamma(T_w-1)=0.6615874197=O(1)\)。因此现有 critical fold--tail 不是定量 Boussinesq-valid prediction，只能先定位为 constant-property linear-density organising structure；必须由 matched low-Mach variable-density calculation 验证是否保留。
- **物理控制参数**：定义 \(\Pi_R=\epsilon\sqrt{Re_R}=\epsilon R/\ell\) 与 \(\Pi_H=\epsilon H_c/\ell\)。\(\Pi_R=O(1)\) 为 thermo-radial interaction，\(\Pi_H=O(1)\) 为 axial confinement。若 heated-BEK transition threshold 为 \(Re_{tr}\)，在转捩前到达 interaction 至少要求 \(\epsilon\sqrt{Re_{tr}}\gtrsim1\)；本项未假定任何普适临界数值。
- **可观测量**：后续二维验证必须输出 wall heat flux、radial/azimuthal wall shear、axial pumping 和 forcing/edge penetration length，而不能只验证 adjoint coefficients。
- **下一步**：先闭合 \(P_4\) meridional subsystem 和 inner/outer boundary data；并行构造 \((\Pi_R,|\beta\Delta T|,\Pi_H,\mathcal G)\) regime map。随后优先用现有 low-Mach solver 检查 \(O(1)\) density contrast 下 fold/tail 是否存在，再决定是否值得求 finite-disk field。
- **输出**：docs/FULL_AXISYMMETRIC_OUTER_HIERARCHY_AUDIT_2026-09-11.md；修订 docs/TWO_DIMENSIONAL_THERMO_RADIAL_OUTER_2026-09-11.md。

### 2026-09-11 — L070：可压缩 von Kármán 空间稳定性重构、算子审计与基准验证

- **状态**：完成一套分离“参数→基流→QEP 组装→目标模态→验证→输出”的 Julia 接口和命令行脚本；复核 Turkyilmazoglu & Uygun (2006) 的严格 von Kármán 切片，并完成有理映射/有限域五组网格交叉验证。尚未把任意 Ro 扩展认定为已验证模型，也未重算精确中性曲线。
- **关键修正 1**：原 notebook 在稳定性矩阵前调用 `interp(...,"phy")`，但论文式 (10)–(15) 使用 Dorodnitsyn–Howarth 相似坐标 y。新实现固定使用 `"sim"`；物理坐标仅允许作后处理。
- **关键修正 2**：旧 `f_q` 的边界回调错误地用数组首尾项代替端点求值，实测给出 f(40)=0.021989…、q(40)=0.085944…，未满足零远场条件。新实现使用 Julia `TwoPointBVProblem`，本次端点最大残差为 4.19e-17。
- **算子审计**：对五场未知量 (rho,u,v,w,T) 的 continuity、三动量和 energy 行核对了压力消元与 Chapman-law 物性扰动结构。紧凑张量组装与独立展开组装在五套网格上的相对差约为 A0=O(1e-16)、A1=O(1e-24)、A2=0；人工 continuity regularization 在新接口中默认为零。
- **基准参数**：Mr=0.1、Tw=1、Ro=-1、R=285.36、beta=0.07759、omega=0、gamma=1.4、Pr=0.72，内部 Ma=Mr/R=3.50434538828e-4。
- **数值结果**：有理映射 N=59/69/79 分别给 alpha=0.3851895715+0.0003923048i、0.3851895712+0.0003923072i、0.3851895692+0.0003923045i；有限域 ymax=40、N=79/99 给 0.3851895722+0.0003922959i、0.3851895702+0.0003923059i。后向误差为 7.54e-15 至 3.79e-14；两类网格约在 1e-8 内一致。
- **解释边界**：论文给出的 alpha_r=0.38482 是目标附近的舍入中性数据；固定其 R、beta 后，本计算得到小但非零 imag(alpha)=3.92e-4。因此本项验证的是收敛的 near-neutral spatial mode 与论文算子切片，不是新的精确中性点。精确中性点需以 imag(alpha)=0 做参数校正并用模态重叠连续跟踪。
- **Ro 证据边界**：继承代码的任意-Ro 版本含 Co=2-Ro-Ro^2、依 Ro 的基流翻号和 `Ro=-1 -> +1` 的内部特殊映射。论文不能独立验证这套约定；新 API 默认拒绝 Ro≠-1，只能通过显式开关产生带 `unverified_bek_extension` 标签的诊断结果。
- **输出**：work/src/CompressibleBEKLinearStability.jl；work/scripts/run_compressible_bek_lsa_20260911.jl；work/scripts/validate_refactored_compressible_bek_20260911.jl；work/results/compressible_bek_lsa_validation_20260911/；docs/COMPRESSIBLE_BEK_LSA_REFACTOR_AND_OPERATOR_AUDIT_2026-09-11.md。

### 2026-09-11 — L071：一般 (Ro) 五方程约定、低马赫诊断与独立 QEP 导数审计

- **状态**：将一般 (Ro) 的五扰动方程、状态方程消元和 (\Gamma=2RoG+Co) 的符号约定写入 `docs/GENERAL_RO_COMPRESSIBLE_FIVE_EQUATIONS_2026-09-11.md`；未把该扩展宣称为论文级验证模型。
- **低马赫诊断**：在 Lingwood (1997) 近似鼻点的 (Ro=1)、(R=27.4)、\(\beta=0.115\) 上，将 (Mr) 从 0.1 降至 0.01、\(N=49\)，目标模态得到 \(\alpha=0.4817132889-0.0017285274i\)，后向误差 \(7.26\times10^{-13}\)。它比 (Mr=0.1) 的错误分支更接近参考 \(\alpha_r\simeq0.487\)，但仍不是中性点，也不能替代连续模态跟踪。
- **算子检查**：新增 `directional_qep_test` 和 `work/scripts/audit_general_ro_operator_20260911.jl`，用 QEP 对 \(\alpha\) 的中心差分方向导数核对 \((A_1+2\alpha A_2)q\)。该审计运行受 Windows/WSL Julia 预编译缓存问题干扰，结果文件未作为证据使用，待修复运行环境后重跑。
- **证据边界**：当前一般 (Ro) 结果仍属于 `unverified_bek_extension`；Lingwood 表值是不可压缩且四舍五入的参考，正式中性曲线前必须做 (Mr\to0)、\(N\) 收敛和模态重叠跟踪。
- **输出**：docs/GENERAL_RO_COMPRESSIBLE_FIVE_EQUATIONS_2026-09-11.md；work/scripts/audit_general_ro_operator_20260911.jl；work/src/CompressibleBEKLinearStability.jl。

### 2026-09-11 — L072：低马赫压力重标度修正

- **发现**：当前压力消元形式含 `1/(gamma*Ma^2)`，直接令 (M_r\to0) 会造成矩阵表示病态，不能作为严格不可压缩极限验证。
- **修正**：有限压力要求 (T\hat\rho+\rho\hat T=O(Ma^2))。应定义
  (\hat\pi=\gamma Ma^2\hat p)，保留正则约束
  (\hat\pi-T\hat\rho-\rho\hat T=0)，再通过 Schur complement 或独立不可压缩 QEP 取极限。
- **结论边界**：(M_r=0.01) 仅是有限低马赫诊断，不能称为 (Ma\to0) 验证。后续 Lingwood 比较必须实现压力重标度/不可压缩算子。
- **输出**：docs/LOW_MACH_PRESSURE_SCALING_AUDIT_2026-09-11.md。

### 2026-09-11 — L073：(M_r=0.3) 参考表比较与模态选择警告

- **计算**：按用户给定参考表，取 von Kármán/Ekman/Bödewadt 对应 (Ro=-1,0,1)，固定 (T_w=1,\omega=0,\gamma=1.4,Pr=0.72)，(N=69)，比较 (M_r=0.3) 结果。
- **结果**：von Kármán 选中根为 (0.3900723-0.0030934i)，实部相对参考值误差 1.365%；Ekman 为 (0.2666017+0.3867959i)，实部误差 3.927%；Bödewadt 为 (0.5634678-0.1505188i)，实部误差 7.717%。各 QEP 后向误差均小于 (2\times10^{-13})。
- **关键判断**：Ekman/Bödewadt 候选根虚部很大，不是中性根；按实部最近选根会把模态选择误差误报为压缩性误差。当前一般 (Ro) 扩展仍未通过物理验证。
- **证据边界**：参考表是不可压缩且四舍五入的临界参数；本次固定 (R,\beta) 并未求解 (\operatorname{Im}\alpha=0) 的中性点，因此不能给出最终 (M_r=0.3) 压缩性误差。
- **输出**：docs/COMPRESSIBLE_MR03_REFERENCE_COMPARISON_2026-09-11.md；work/scripts/compare_compressible_bek_reference_Mr03_20260911.jl；work/scripts/list_compressible_bek_modes_Mr03_20260911.jl；work/results/compressible_bek_reference_Mr03_20260911/。

### 2026-09-11 — L074：旧 notebook 的 (Ro=0,1) 交叉比较

- **方法**：按旧 `CRD_STA_CACU.ipynb` 的实际函数链重跑 (M_r=0.3,T_w=1,N=69) 的 Ekman 与 Bödewadt 参考点。
- **结果**：旧程序的首个候选分别为 Ekman (0.318390-0.0959599i)、Bödewadt (0.550536-0.108147i)；新程序按目标实部选择的候选为 Ekman (0.266602+0.386796i)、Bödewadt (0.563468-0.150519i)。同一次 IAR 中还存在多个其他复根，返回顺序不稳定。
- **判断**：两套程序不能在 (Ro=0,1) 下逐根对应。差异来自旧 `phy` 坐标插值、旧 `f_q` 远场边界、基流翻号/`Ro=-1` 特殊分支及缺少模态重叠跟踪；旧结果不能作为一般 (Ro) 的独立真值。
- **输出**：docs/LEGACY_RO01_COMPARISON_2026-09-11.md；work/scripts/run_legacy_crd_sta_ro01_reference_20260911.jl。

### 2026-09-11 — L075：旧程序切换相似坐标后的 (Ro=0,1) 诊断

- **方法**：仅把旧 notebook 的 `interp(...,"phy")` 改为 `interp(...,"sim")`，其余旧算子和旧基流流程保持不变；(M_r=0.3,T_w=1,N=69)。
- **结果**：Ekman 首根 (0.3171723-0.0955872i)，Bödewadt 首根 (0.5619501-0.1140132i)。与旧物理坐标首根 (0.318390-0.0959599i)、(0.550536-0.1081474i) 变化很小。
- **判断**：坐标修正本身不能恢复给定临界参考值；两流动仍没有明确的近中性根。剩余差异来自旧热 BVP、(Ro) 符号/特殊分支和模态跟踪缺失。
- **输出**：docs/LEGACY_SIMCOORD_RO01_COMPARISON_2026-09-11.md；work/scripts/run_legacy_crd_sta_simcoord_ro01_20260911.jl；work/results/legacy_crd_sta_simcoord_ro01_20260911/。

### 2026-09-12 — L076：可压缩求解器的半无限网格修正与 Ro=0,1 连续模态复核

- **状态**：将新可压缩求解器的 rational Chebyshev 坐标改为真正的半无限表示；末节点为 y=Inf，全部有限节点严格递增。原 RotatingDiskFlow 后端和 baselines 未修改。
- **映射审计**：旧实现虽然按无穷远代数映射构造 D、D2，却把 y>40 的坐标全部截成 40；N=59,69,79,89,99 分别产生 13,15,17,20,22 个重复 y=40 节点。新实现保留相同 D、D2（相对变化均为 0），恢复真实有限节点并令唯一末节点为 Inf；y>40 的基流采用 F=0、G=1、T=rho=1、H=H∞，跨网格特征向量在 s=y/(1+y) 上插值。
- **热基流**：Mr=0.3 时 Ro=0 的 max|T-1|=3.33e-16；Ro=1 的 max(T)=1.1226852428、max|rho-1|=0.1092783962；热 BVP 端点误差小于 3.5e-17。
- **连续模态**：N=69 沿 Mr=0.01→0.30 的重叠跟踪得到 Ekman alpha=0.315366030307-0.088067397220i，Bödewadt alpha=0.563467793123-0.150518845867i；两者仍非中性。最终步重叠分别约为 1 与 0.998749。
- **网格 Gate**：Ekman 的 N=59/69/79 根收敛到约 0.3153660303-0.0880673972i。Bödewadt 在 N=59/69/79/89/99 的高重叠连续跟踪仍摆动，N=99 为 0.559223368023-0.160621751355i，相邻重叠均大于 0.993；因此半无限坐标修正没有解决该根簇的网格不收敛。
- **结论边界**：该修改消除了真实的坐标/采样缺陷，但 Ro=0,1 不一致仍不能归因于坐标截断。一般 Ro 算子和符号约定仍属 unverified_bek_extension；结果是 frozen-local-base-flow 的 similarity-subspace spatial diagnostic，不是自洽中性曲线或完整三维稳定性结论。
- **输出**：work/src/CompressibleBEKLinearStability.jl；work/scripts/audit_compressible_ro01_mode_tracking_20260912.jl；work/results/compressible_solver_ro01_audit_20260912/；docs/COMPRESSIBLE_RO01_SOLVER_AUDIT_2026-09-12.md。

### 2026-09-12 — L077：由两份可压缩 von Kármán 文献回溯的一般 BEK `Ro/Co` 算子审计

- **状态**：逐项核对 Turkyilmazoglu--Uygun (2006) 式 (1)--(15)、William 学位论文第 2 章旋转系原方程与 disk limit，并以论文第 3 章线性化及 Lingwood 一般 BEK 局部尺度检查当前代码；本项未运行新中性曲线计算。
- **母方程结论**：一般 BEK 使用 `Ro=DeltaOmega/Omega`、`Co=2-Ro-Ro^2` 和以有符号 `DeltaOmega` 缩放的速度。局部稳定方程的 Doppler `alpha*U+beta*V-omega` 不整体乘 `Ro`；`Ro` 出现在柱坐标连续性曲率、动量的 `U*u`/`W*D`/`W'`、能量的法向输运及对应密度--基流加速度项；`Co` 直接进入径向--切向 Coriolis 块，并因 `rho_hat` 乘基流加速度进入相应密度列，但不直接进入连续性、轴向动量、状态或能量方程。
- **代码判定**：`Spatial_mode_BEK1` 与 `Spatial_mode_BEK` 是同一错误一般化的两种装配。明确问题包括 `Ro=-1 -> +1` 的不完整特殊变换、连续性 `u/r` 固定为单位系数、能量行 `H` 平流未继承 `Ro`，以及将当前 `(F,G,H)=(-U,V,-W)` 当作直接 `(U,V,W)` 使用。两套算子的机器精度一致不能构成独立验证。
- **热方程勘误**：`source_profiles` 当前的 `q` 与 `f` 两个线性输运项符号与 `q''+Pr*Ro*H*q'=0` 及 `Ro=-1` 论文热方程相反；Mach 数究竟以 `Omega` 还是 `DeltaOmega` 定义尚未固定，因此耗散项显式 `Ro^2` 的位置也尚不能封板。
- **包装器判定**：`assemble_qep` 用 `C(0)+Ro[C(1)-C(0)]` 绕过特殊分支只是旧错误算子的仿射外推，不能修复缺项或变量约定；现有 Ekman/Bödewadt 根是低后向误差的错误-QEP 根，不可用于物理结论。
- **证据边界**：本项否定当前一般 `Ro` 可压缩 QEP，但没有实现替代算子。下一步必须从一个旋转系可压缩守恒母方程生成唯一 QEP，并通过 `Ro=-1` 论文逐项回归、`Ro=0` Ekman 远场 Gate 与 `Ro=1,Co=0` 无直接 Coriolis Gate 后再计算中性曲线。
- **输出**：docs/COMPRESSIBLE_BEK_RO_CO_OPERATOR_AUDIT_2026-09-12.md；勘误 docs/GENERAL_RO_COMPRESSIBLE_FIVE_EQUATIONS_2026-09-11.md。

### 2026-09-12 — L078：一般 BEK 候选算子修正、低马赫端点 Gate 与定 β 中性参数

- **实现**：在 `work/src/CompressibleBEKLinearStability.jl` 中统一采用 Lingwood 直接变量 `(U,V,W)`，只在读取 `sol_baseflowODE` 时撤销后端历史输出翻号；保留的 `RotatingDiskFlow/src/CRD_STA.jl` 与 `baselines/` 均未修改。连续性 `u/r`、基流法向输运和能量 `W` 平流/压力功伙伴统一使用 `+Ro`；`Gamma=2RoV+Co` 保留在径向--切向耦合及相应密度列。早期实现误用转换变量的 `-Ro`，其 Bödewadt 失败结果已由本条取代。
- **热模型**：以 `Mr=r|DeltaOmega|/a_inf`、`Ma=Mr/R` 定义 Mach 数，热 BVP 用满足双端 Dirichlet 条件的中心三对角离散；温度重构不再额外使用 `(Ro*Mr)^2`。半无限 rational grid 保持唯一 `Inf` 端点。
- **结构 Gate**：修正的连续性、`Gamma`、能量和 `Ro=1,Co=0` 块残差为零或舍入量；固定参考点的 QEP 后向误差均低于 `6e-11`。旧 compact/expanded 审计因两者使用不同变量约定而正式退役，不能再当物理验证。
- **低马赫端点 Gate**：在 `Mr=0.01`、`N=49/59` 的 Lingwood 舍入鼻点，von Kármán 得 `alpha=0.380609-0.000545i` (`N=49`)，Ekman 得 `0.529744-0.000010i`，Bödewadt 得 `0.482033-0.000354i`。三端点均恢复近中性物理分支；这是有限低马赫 Gate，不是严格不可压缩极限。
- **模态选择判定**：用户给定 `alpha_ref=0.2775` (Ekman) 与 `0.5231` (Bödewadt) 在原固定 `(R,beta)` 目标求解中选中了其他空间分支。以特征向量重叠从 Lingwood 模态沿 `Mr:0.01->0.3`、再沿 `beta` 连续后，两流动都能在扫描区间找到唯一固定-β 中性交点。
- **中性结果**：`Mr=0.3,Tw=1,omega=0,gamma=1.4,Pr=0.72,N=69` 下，Ekman (`Ro=0,beta=0.132`) 得 `R=119.5358824462`、`alpha=0.4983641183+1.59e-11i`；Bödewadt (`Ro=1,beta=0.127`) 得 `R=32.9491248879`、`alpha=0.6169224711+2.55e-11i`。对应后向误差分别 `2.08e-19` 与 `2.95e-18`。
- **网格范围**：Ekman 的 `N=49/59/69` 中性 `R=119.62367869/119.48743987/119.53588245`，全跨度约 `0.136`，故只报告当前 `N=69` 候选且仍需更高阶/映射检查；Bödewadt 为 `32.94912602/32.94912487/32.94912489`，已高度收敛。
- **证据边界**：这些是 derived-BEK candidate 的 frozen-local-base-flow、similarity-subspace、固定-β 空间中性点，不是两参数中性曲线鼻点、完整三维稳定性或两份 von Kármán 文献直接验证的一般-Ro 可压缩结论。
- **输出**：`docs/COMPRESSIBLE_BEK_RO_CO_OPERATOR_AUDIT_2026-09-12.md`；`docs/COMPRESSIBLE_BEK_CORRECTED_FIXED_BETA_NEUTRAL_STUDY_2026-09-12.md`；`work/scripts/check_corrected_lingwood_low_mach_20260912.jl`；`work/scripts/continue_corrected_compressible_neutral_ro01_20260912.jl`；`work/scripts/continue_lingwood_anchored_compressible_neutral_ro01_20260912.jl`；`work/scripts/refine_corrected_neutral_grid_ro01_20260912.jl`；`work/results/compressible_bek_corrected_operator_v8_20260912/`；`work/results/compressible_bek_corrected_lingwood_gate_v2_20260912/`；`work/results/compressible_bek_lingwood_anchored_neutral_20260912/`。

### 2026-09-15 — L080：一般 Ro 热 BVP 符号修正与 (M_r=0.3) 重测

- **修改**：将 `source_profiles` 的热输运和反应系数改为 `transport=Pr*Ro*W`、`reaction=2Pr*Ro*U`，对应统一原始能量方程的 `q''-Pr*Ro*W*q'=0` 和 `f''-Pr*Ro*W*f'-2Pr*Ro*U*f=2Pr*shear`。保留后端和 baselines 不变。
- **原始残差 Gate**：`N=120` 时 von Kármán/Ekman/Bödewadt 的文献兼容能量绝对残差分别为 `8.09e-8/3.68e-7/3.15e-7`，项平衡比为 `1.78e-6/2.59e-6/2.98e-6`；`y<=39` 扩展区结果相同。速度、状态残差继续在插值误差范围内。原错误符号的自残差变为 `1.87e-3` 和 `6.84e-2`，确认测试独立有效。
- **远场结果**：修正后 von Kármán 壁温 `q` 在 `y=10` 已降至 `2.82e-3`；Bödewadt 仍在 `y=39` 保持 `0.622`，过渡贴着人工端点；Ekman 仍是 `q=1-y/40` 的 `Ro=0` 退化解。因此微分方程通过，但 Bodewadt 壁温和严格 Ekman 半无限热极限尚未通过远场 Gate。
- **中性参数重测**：Ekman 物理锚定支在 `Mr=0.3,beta=0.132` 的 `N=49` 候选为 `R=119.6237`，`N=69` 同点虚部为 `2.27e-4`；Bödewadt `N=49` 中性候选为 `R=56.5816`，但 `N=59/69` 同点得到 `0.46241-0.09045i` 与 `0.52317-0.01498i`，网格不收敛。供应的不可压缩 Bodewadt 参考点在 `Mr=0.01` 给 `alpha=0.48170-0.00180i`，而 `Mr=0.3` 延拓后为 `0.54080-0.20592i`，不是中性。
- **结论边界**：热方程符号错误已修正；原始微分残差 Gate 通过，但完整一般-Ro 可压缩基本流尚未完全可信。Bodewadt 和严格 `Ro=0` 热远场问题解决、并完成独立高阶中性重算之前，所有 `Mr=0.3` 中性值仍为诊断候选。
- **输出**：`work/src/CompressibleBEKLinearStability.jl`；`work/scripts/validate_compressible_bek_baseflow_primitive_residuals_20260912.jl`；`work/scripts/continue_lingwood_anchored_compressible_neutral_ro01_20260912.jl`；`work/results/compressible_bek_baseflow_primitive_residual_gate_v3_20260915/`；`work/results/compressible_bek_lingwood_anchored_neutral_v3_20260915/`；`docs/COMPRESSIBLE_BEK_THERMAL_BVP_CORRECTION_2026-09-15.md`。

### 2026-09-15 — L081：热远场域长与 (Ro=0) 退化处理

- **实现**：`source_profiles` 新增 `ymax` 参数，热源网格与 `prepare_case` 的数值域同步；在精确 `Ro=0` 且 `Mr>0` 或 `Tw!=1` 时，将 `PreparedCase.evidence_scope` 标记为 `:ro0_thermal_degenerate`。
- **域长证据**：Bodewadt 壁温在 `L_y=40,60,80` 都有 `q(L_y-1)=0.621518`，而 `q(10)` 约为 1，说明过渡随人工端点移动。`Ro=0` 的 `q(1)` 从 `0.975`、`0.9833` 变为 `0.9875`，固定位置不收敛。
- **中性重定位**：Bodewadt `Mr=0.3,beta=0.127` 在各域长独立求得 `R=56.58165,56.61959,65.70576`；`L_y=80` 候选在 `N=59/69` 的虚部分别为 `-0.05744/-0.10136`，网格与域长均未通过。
- **结论**：完整可压缩 BEK 基本流仍不能封为全部 `Ro` 均可信。热微分方程已通过，但 Bodewadt 半无限热远场需要出流相容/非相似处理，`Ro=0` 需要奇异极限模型；当前 `Mr=0.3` Bodewadt 中性参数不报告为最终物理结果。
- **输出**：`work/src/CompressibleBEKLinearStability.jl`；`work/scripts/validate_bek_thermal_domain_and_ro0_20260915.jl`；`work/scripts/continue_bodewadt_neutral_domain_20260915.jl`；`work/results/compressible_bek_thermal_domain_ro0_20260915/`；`work/results/compressible_bek_bodewadt_neutral_domain_20260915/`；`docs/COMPRESSIBLE_BEK_THERMAL_BVP_CORRECTION_2026-09-15.md`。

### 2026-09-12 — L079：一般 Ro 可压缩基本流的稳态原始方程连续残差 Gate

- **模型与方法**：把当前 Lingwood 直接变量 `(U,V,W)`、温度 `T` 和 `rho=1/T` 代回 DH 变换前的统一轴对称可压缩 BEK 边界层方程；采用有符号差速速度尺度、`Co=2-Ro-Ro^2`、Chapman `mu=T`、`k=T/Pr`。分别计算连续、径向动量、切向动量、状态和能量逐点残差；`Ro=(-1,0,1)`、`N=(80,120,160,200)`，并同时检查 `y<=12` 与 `y<=39`。本项是局部相似边界层 Gate，不是完整有限半径 NS 残差。
- **通过部分**：`N=120` 时三端点连续/两动量方程的最大绝对残差分别约为 `9.65e-8`、`1.47e-7`、`5.38e-8`，状态残差为零；全部网格的速度残差不超过约 `3.1e-6`，高阶增长来自对插值源剖面的微分。因此当前 `(U,V,W)` 速度基本流在本 Gate 的精度范围内可信。
- **失败部分**：`(Mr,Tw)=(0.3,1)` 时，文献兼容原始能量残差在 `Ro=-1` 为 `4.29e-3`（项平衡比 `9.41e-2`），在 `Ro=1` 为 `1.50e-2`（项平衡比 `1.41e-1`），且不随 `N` 下降；与此同时当前生成方程的自残差仅为 `8.05e-8` 与 `4.16e-8`。这证明问题是热 BVP 的输运/反应符号，而非离散误差。
- **远场诊断**：`(Mr,Tw)=(0,0.8)` 时，von Karman 的 `q(30)=0.9983`、`q(35)=0.9586`、`q(39)=0.4710`，温度跃迁被推到人工端点，`y<=39` 能量残差为 `1.21e-1`；Bodewadt 近壁能量残差为 `1.53e-1`。`Ro=0` 退化为 `q''=0`，有限区间线性解不能在半无限域同时满足壁面与远场 Dirichlet 条件。
- **Gate 结论**：完整 `(U,V,W,T,rho)` 可压缩基本流 **FAIL**。速度与状态关系通过，但当前 `T` 及由它耦合得到的 `rho` 不能称为统一一般 `Ro` 的可信基本流；L078 的中性参数继续仅作 `derived_bek_candidate` 诊断，不能升级为自洽可压缩 BEK 结果。下一步应先修正/重构热 BVP（尤其处理 `Ro=0` 奇异极限），再重复本残差 Gate 和中性参数计算。
- **输出**：`work/src/CompressibleBEKBaseflowResiduals.jl`；`work/scripts/validate_compressible_bek_baseflow_primitive_residuals_20260912.jl`；`work/results/compressible_bek_baseflow_primitive_residual_gate_v2_20260912/`；`docs/COMPRESSIBLE_BEK_BASEFLOW_PRIMITIVE_RESIDUAL_GATE_2026-09-12.md`。
