# Leading outer 极限的参数选择与增广 Jacobian 检验

## 本步骤要回答的问题

附件提出：在 `(Ro_t,Tw,t)` 构造包含场变量和参数的 outer augmented Jacobian，
检查它是否具有 fold 零模，并用零模分量比解释 `tau1/r1`。本步骤直接检验这一
命题，同时明确当前渐近小参数的物理归一化。

本步骤不新增有限 `epsilon` 基本流或 continuation，不做 Hessian、adjoint
solvability 或一阶/二阶完整 matching。数值工作只求 leading outer 极限 BVP。

## 先验结构核对

当前小参数不是任意 bookkeeping parameter，而是原问题的物理输出/约束

```text
epsilon = -Hinf > 0.
```

在 `H=epsilon*h(Y), Y=epsilon*eta` 下，它给出精确归一化

```text
h(infinity) = -1.
```

因此重标度 `epsilon_tilde=c*epsilon` 会把远场条件改成
`h_tilde(infinity)=-1/c`；保持 `h(infinity)=-1` 后不能再任意重标度 epsilon。

此外，有限-epsilon fold 点满足的增广系统实际包含三部分：基本流残差、固定
`Tw` Jacobian 的零向量方程、以及 `Hinf+epsilon=0`。只把 `(Ro,Tw)` 加到
outer 基本流残差中，并不等于 fold 增广系统。

## Leading outer BVP

沿径向慢流形消去温度：

```text
q = s+r*g,
chi = omega^2/q^2,
theta(g;r) = 1 + [1-omega^2/q^2]/gamma.
```

令 `p=g'`。方位和温度方程可化成三个一阶方程

```text
g' = p,
p' = Pr*h*p + 3*r*p^2/q,
h' = [p' + chi*r*h*p]/(chi*q).
```

边界/归一化为

```text
g(0)=0, h(0)=0,
g(infinity)=1, h(infinity)=-1.
```

未知初始斜率 `p(0)` 和参数 `r` 由两个远场条件共同选择。壁温是输出
`Tw,t=theta(0;r)=Tcurve(r)`。

将用 Julia RK4 shooting 求解，并用中心差分构造约化增广 shooting Jacobian

```text
d(g(Ymax)-1, h(Ymax)+1) / d(p(0),r).
```

若它可逆，则归一化 leading outer BVP 是孤立解，没有附件所假定的基本流
augmented-Jacobian 零模；若它奇异，才支持该假设。将检查步长和 `Ymax` 收敛。

完成计算后在本文追加数值结果，并按项目规则新增研究日志记录。

## 数值结果

### Outer 极限独立选择临界坐标

shooting 未使用有限 `epsilon` 拟合坐标作为固定参数，只把它们用作 Newton 初值。
收敛结果为

```text
g_Y(0) =  0.7881936107,
Ro_t   = -0.3207391663,
Tw,t   =  1.4948092622.
```

其中 `Tw,t` 由径向慢流形在 `g(0)=0` 处输出。与四小点有限-epsilon二次外推

```text
Ro_t = -0.3207197304,
Tw,t =  1.4947866217
```

的差分别为 `-1.944e-5` 和 `2.264e-5`，均位于 L015 已传播的保守数值范围内。
这说明 leading outer BVP 本身能够选出同一个 fold--tail 极限，而不是只能在给定
`Ro_t` 后计算剖面。

`Ymax=20` 到 `30` 时 `Ro_t` 改变 `1.76e-7`；`Ymax=30` 到 `40` 只改变
`1.3e-10`。在 `Ymax=30`，RK4 步长从 `0.01` 减到 `0.0025` 后显示的
`Ro_t,Tw,t` 不变，末端残差低于 `8e-15`。

### 约化增广 Jacobian 不奇异

在最细结果上，约化 shooting Jacobian 为非奇异矩阵，其奇异值和行列式为

```text
sigma_max = 3.41513952,
sigma_min = 0.95500843,
condition = 3.57603076,
determinant = 3.26148703.
```

中心差分相对步长从 `1e-4` 改到 `1e-7` 时，`sigma_min` 始终为
`0.9550084`、条件数始终为 `3.57603`。因此不存在接近零的奇异值；这不是离散
步长掩盖的弱零模。

这里使用的是消去场初值自由度后的约化增广 Jacobian

```text
d(g(Ymax)-1,h(Ymax)+1) / d(g_Y(0),Ro).
```

它与归一化 leading outer BVP 的局部可解性等价，但不等同于有限-epsilon
固定-`Tw` 基本流 Jacobian。

### `tau1/r1` 的正确来源

outer 解给出的解析相容曲线斜率为

```text
dTcurve/dRo = -1.2154313135.
```

有限-epsilon四小点拟合则给出

```text
A_T/A_R = 0.5106601141/(-0.4203424951)
        = -1.2148667338.
```

二者只差 `5.646e-4`（相对约 `4.65e-4`）。这个比例不需要诉诸 outer
Jacobian 零向量。由于 L015 已显示

```text
Tw(epsilon)-Tcurve(Ro(epsilon)) = O(epsilon^2),
```

对它在 `epsilon=0` 求导便直接得到

```text
A_T = Tcurve'(Ro_t)*A_R,
```

所以 `tau1/r1` 是极限相容曲线的切线斜率，而不是当前 leading outer 基本流
BVP 的 null-vector 分量比。

## 对附件命题的判定

附件对 outer 标度、DAE 慢流形以及 `h_0(0)=0` 的识别是正确且有价值的；但把
leading outer 基本流增广 Jacobian 预期为 fold 奇异矩阵，需要修正：

1. `epsilon=-Hinf` 已通过 `h(infinity)=-1` 固定尺度，不需要另造
   `mu-mu_t=epsilon^2` 才能使 `A_R` 有绝对意义；
2. 归一化 leading outer BVP 是一个正则的特征参数问题，它独立选择 `Ro_t`，
   约化 Jacobian 明显可逆；
3. 每个有限-epsilon点的“fold”条件属于完整系统
   `R=0, J_U|Tw-fixed*v=0, v'*v=1, Hinf+epsilon=0`。只把 `Ro,Tw` 加入 outer
   基本流残差并没有包含 `J_Tw*v=0`，因此不能把两种 augmented Jacobian
   混同；
4. 固定-`Tw` fold 零模在 `epsilon->0` 时可能具有 inner/outer 不同缩放，仍需
   在完整 fold 增广系统中展开。leading outer 基本流 Jacobian 的可逆性并不
   否定有限-epsilon fold，而是否定“`A_R` 来自该基本流 Jacobian 零模幅值”的
   特定解释。

因此，现在不能使用附件中的

```text
A_R^2 = -2*mu2*<Psi,R_mu>/<Psi,R_XX[Phi,Phi]>
```

来解释 `A_R=-0.42034`：当前问题没有识别出这样的独立二次 detuning，且被测试
的 leading outer 基本流算子没有零模。若以后计算 `A_R`，正确对象应是经过
singular scaling 和 matching 后的**完整 fold 增广系统的一阶校正**，其中
`epsilon=-Hinf` 通过远场归一化直接进入，而不是当前基本流 BVP 的标准
saddle-node Hessian 公式。

## 输出

- Julia 程序：`work/scripts/solve_outer_limit_and_jacobian.jl`；
- 收敛与 Jacobian：
  `work/results/outer_limit_selection_and_jacobian/convergence_and_jacobian.csv`；
- Jacobian 差分步长检查：
  `work/results/outer_limit_selection_and_jacobian/jacobian_fd_sensitivity.csv`；
- leading outer 剖面：
  `work/results/outer_limit_selection_and_jacobian/outer_profile.csv`；
- 与有限-epsilon外推比较：
  `work/results/outer_limit_selection_and_jacobian/comparison_with_fold_extrapolation.txt`。
