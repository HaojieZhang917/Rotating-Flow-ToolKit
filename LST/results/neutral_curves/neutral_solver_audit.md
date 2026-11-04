# 空间 LST 中性曲线求解器核查

- 算子：`VonKarmanLST.jl` 的 `L0 + alpha*L1 + alpha^2*L2 = 0`，未重新构造 LST。
- 参考系：随旋转圆盘坐标系；基本流为 `U=F, V=G-1, W=H`。
- 全局到局部转换：`beta=m/R`、`omega=omega1/R`。算子中的 `beta*R` 与 `omega` 项均按现有工程原样使用；没有再次乘 R。
- 空间中性条件：`imag(alpha)=0`；增长率使用 `sigma=-imag(alpha)`，没有混用 temporal 特征值。
- 边界条件：有限 `zmax=30.0` 的齐次壁面/外缘截断，采用当前已验证的谱离散入口。
- 分支追踪：R 粗网格变号 + 一维二分；相邻点同时使用 alpha 连续性和离散速度特征函数重叠。
- Type-I 参考：R=285.36、beta=0.07759、alpha约为 `0.38505249+0.00034867im` 的现有回归点。
- Type-II 自动识别：当前没有将未有物理参考种子的分支强行命名为 Type-II；若没有独立 Type-II 参考模态，结果会标为未分类。

## 本次配置
`N=25`, `zmax=30.0`, `R=50.0:25.0:800.0`, `m=-80:10:150`, `omega1=[-5.0, 0.0, 4.0, 7.9, 10.0]`。

输出数据见 `neutral_type1_curves.csv`，中性点直接特征函数见 `neutral_direct_modes.jls`。
