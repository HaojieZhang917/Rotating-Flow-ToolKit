# Generic thermal-tail 一阶内区退化审计（2026-09-09）

## 研究问题

在完全不使用 fold 条件 `J_U v=0` 的情况下，检验 acceleration-consistent
basic-flow thermal-tail family 是否在 `Ro_t=-0.4703674163` 附近发生一阶
meridional-core 标度退化。

## 模型、路径与数据范围

- corrected-energy consistent model，`Pr=0.72, gamma=1`；
- 只求解 `R(U,Ro,Tw)=0`，固定 `Hinf=-epsilon` 并释放 `Tw`；
- 主扫描：`Ro=-.469,-.467,-.465,-.46,-.44,-.42,-.40,-.35,-.30`，
  `epsilon=.08,.05,.03,.02,.015,.01`，`N=240,320`；
- 临界邻域：增加 `Ro=-.4700,-.4702`，并以 `N=240,280,320,360` 复核；
- 映射检查（`N=280`）：`(a,b,c)=(2,.6,.5),(2,.55,.5),(2.5,.6,.5)`。

所有数据位于 `work/results/generic_tail_zero_core_20260909/`。没有调用 fold、零模或
增广系统。

## 指标

连续方程在壁面给出

```text
D = chi(Tw)*s(Ro)^2-omega(Ro)^2 = Ro*F''(0).
```

pilot 表明 `D->0` 和未缩放 core `->0` 是 generic tail 的共同零阶性质，故主指标改为

```text
D/epsilon = D1(Ro)+epsilon*D2(Ro)+...
```

并同步外推固定 eta 的 `H/epsilon`、`F/epsilon` 与

```text
Mmer/epsilon = [integral_0^8(F^2+H^2)deta]^(1/2)/epsilon.
```

## 结果

`N=320`、`epsilon<=.03` 二次截距给出的代表性一阶系数为：

| Ro | D1 | H1(eta=8) | Mmer,1 |
|---:|---:|---:|---:|
| -.469 | 2.3416e-3 | -1.1492e-2 | 2.5857e-2 |
| -.460 | 1.7668e-2 | -7.9305e-2 | 1.9137e-1 |
| -.440 | 5.0712e-2 | -2.1644e-1 | 5.3304e-1 |
| -.400 | 1.1103e-1 | -4.4139e-1 | 1.1026 |
| -.300 | 2.1371e-1 | -7.9171e-1 | 2.0170 |

这些量从普通 tail 向 `Ro_t` 系统趋零，支持 generic tail 具有非零 `O(epsilon)`
meridional core，而该一阶 core 在端点退化。

临界邻域对 `D/epsilon` 先作 epsilon 截距，再用
`Ro=-.469,-.4700,-.4702` 作局部线性零点，得到：

| N / mapping | epsilon window | D1 zero |
|---|---|---:|
| 240, default | <=.03 quadratic | -.47039488 |
| 320, default | <=.03 quadratic | -.47036350 |
| 280, default | <=.03 quadratic | -.47037685 |
| 280, b=.55 | <=.03 quadratic | -.47038077 |
| 280, a=2.5 | <=.03 quadratic | -.47036374 |

映射检查的范围为 `[-.47038077,-.47036374]`，包含独立 outer BVP 的
`Ro_t=-.4703674163`。综合网格和拟合窗口，当前安全表述为

```text
Ro_*^(basic-flow,D1) = -0.47038 +/- 0.00003,
```

并与 outer `Ro_t` 一致。特别地，`N=320` 近临界短距离外推与 outer 值只差
`3.9e-6`，但该单一最优差值不作为总体误差条。

固定 eta 的 `H1` 零点在 `eta=1,2,5` 已靠近同一区间；例如 `N=320`、
`epsilon<=.03` 时依次约为 `-.4703827,-.4703654,-.4703914`。`eta=8` 则漂移至
`-.4704763`，说明当前 epsilon 范围尚未形成可同时满足 `eta>>1` 与
`epsilon*eta<<1` 的干净 plateau。因此 `h0=H1(infinity)` 的零点仍属 emerging
evidence，不能宣称已独立高精度确定。

## 数值与逻辑边界

1. 壁面 collocation 行被边界条件替换，谱矩阵直接得到的 `F''(0)` 在小 epsilon
   下未收敛；主结果使用参数代数式 `D`。壁面导数只作未来的交叉检查。
2. `Mmer,1` 无符号，只支持幅值消失，不能单独证明穿零。
3. effective exponent 只显示 crossover；临界定位使用一阶系数截距，而不是
   `p_eff=2`。
4. 当前结果证明的是 basic-flow tail family 的强数值机制证据，不是统一
   inner--outer 存在/唯一性定理，也尚未推导完整 generic-tail leading BVP。
5. fold 未参与本次计算。它与 `Ro_*` 的重合仍是第二条、独立的已有 bifurcation
   证据，而不是本扫描的输入。

## 当前结论

现有数据不支持“generic tail 非零 core、临界点 zero core”的二分法，而支持更精细的
标度退化：整个 regular tail family 在零阶趋向 zero core；普通 tail 在一阶保留
`F1,H1` meridional core，而该一阶系数在 `Ro_*≈Ro_t` 消失，使临界类提升为
`F,H=O(epsilon^2)`。这为不借助 fold 条件解释 `Ro_t` 的特殊性提供了第一组独立证据。

## 下一计算：generic leading outer family

- 固定 `Ro`，将 outer 入口 `h(0)=h0` 与 `p0=g_Y(0)`同时作为未知量，满足
  `g(infinity)=1,h(infinity)=-1`；不施加 `h0=0`，也不调用 fold 条件。
- 首算法以单调 `g` 为坐标，从 `(g,p,h)=(1,0,-1)`反向积分至 `g=0`，直接得到
  `(p0,h0)`；用步数 `1000--16000` 检查离散收敛。
- 第二算法在 `Y` 坐标由所得 `(p0,h0)`正向积分，检查有限 `Ymax` 的远场残差，
  并与 fixed-epsilon basic-flow 外推的 `H1`比较。
- 扫描范围覆盖 `Ro=-.50...-.30`，并在 `-.4703674163`附近加密。只有当
  `h0(Ro)` 的零点对步长、截断和两种坐标实现稳定时，才将其解释为 generic
  outer family 的入口输运退化点。

## Generic outer family 结果

以 `g` 为自变量，从远场 `(g,p,h)=(1,0,-1)`反向积分至 `g=0`，直接得到每个
固定 `Ro` 的 `(p0,h0)`。再将该入口数据放回 `Y` 坐标正向积分；在 `Ymax=60`
上代表性远场误差约为 `1e-10--1e-7`，误差随 `Ro` 远离临界点而增大但不影响
`g` 坐标结果。

`h0=0` 的零点随 `g` 步数收敛如下：

| g steps | Ro at h0=0 | Tw | p0 |
|---:|---:|---:|---:|
| 250 | -.470367416465 | 1.661587419720 | .654193538533 |
| 1000 | -.470367416465 | 1.661587419720 | .654193538533 |
| 4000 | -.470367416464 | 1.661587419719 | .654193538532 |
| 16000 | -.470367416464 | 1.661587419719 | .654193538532 |

这证明 critical outer BVP 正是 generic outer family 的 `h0=0` 成员。不过它与旧
critical outer 计算使用同一 leading ODE，因此是结构重释和算法一致性，不是独立于
outer 方程的新物理验证。

## Generic first-order inner BVP 与 matching

对每个 fixed `Ro`，使用 outer 给出的 `(p0,h0)`求解线性一阶 inner 系统：

```text
H1'+2F1=0,
F1''-2*chi0*s*G1+(gamma*s^2/Ro)*Theta1=0,
G1''+chi0*Co*F1=0,
Theta1=tau1+Btheta*eta.
```

边界/匹配条件为 `F1(0)=G1(0)=H1(0)=0`、`F1(infinity)=0`、
`G1'(infinity)=p0`、`H1(infinity)=h0`。有限差分 BVP 在
`L=12,16,20,24` 与 `deta=.01,.005` 下收敛；代表性线性代数残差为
`1e-11--1e-10`。

例如：

| Ro | outer h0 | inner tau1 | inner D1 |
|---:|---:|---:|---:|
| -.400 | -.4407510534 | -.08851544 | .11103376 |
| -.469 | -.0106005914 | -.00185740 | .00234877 |
| Ro_t | approximately 0 | approximately 0 | approximately 0 |

`Ro=-.40` 的 inner `D1=.1110338` 与无-fold finite-epsilon 外推
`.1110344` 一致。将完整 `N=320` fixed-H 解的
`F/epsilon,H/epsilon,G/epsilon,(T-T0)/epsilon` 与 inner BVP 比较，普通 tail
上的 profile error 随 epsilon 系统下降；例如 `Ro=-.40`、`epsilon=.03 -> .01`
时 relative-L2 误差由 `21.3% -> 7.4%`（F1）、`23.8% -> 8.3%`（H1）、
`5.87% -> 1.97%`（G1）。靠近 `Ro_t` 时一阶目标本身趋零，相对误差被
`O(epsilon^2)` crossover 放大，只能使用绝对误差或更高阶展开。

## 局部退化律

在 `|Ro-Ro_t|<=1.64e-3` 内，对收敛的 generic outer/inner 解作过原点的局部
线性--二次拟合，得到

```text
h0     = -7.78537*(Ro-Ro_t) + O((Ro-Ro_t)^2),
D1     =  1.71942*(Ro-Ro_t) + O((Ro-Ro_t)^2),
tau1   = -1.35962*(Ro-Ro_t) + O((Ro-Ro_t)^2),
Mmer,1 = 19.0614*abs(Ro-Ro_t) + O((Ro-Ro_t)^2).
```

前三个拟合 RMS residual 分别约为 `1.8e-8,5.8e-9,4.1e-9`；无符号 norm 的
残差约为 `5.3e-5`。这支持 `h0,D1,tau1` 同时线性穿零，而整个 meridional
一阶场以线性幅值塌缩。

## 更新后的机制结论

目前得到的闭合图景是：generic tail 的零阶 core 均为零，但在一阶具有非平凡
`(F1,H1)`及入口输运 `h0`；`Ro->Ro_t` 时 `h0,D1,tau1` 与 meridional norm
同时消失，generic inner--outer matching 失去一个入口自由度，系统进入
`F,H=O(epsilon^2)`的 critical distinguished limit。critical outer BVP 因而是
generic outer family 的退化成员；finite-epsilon 无-fold basic-flow 外推和另一条
finite-epsilon fold chain 分别提供独立数值支持。
