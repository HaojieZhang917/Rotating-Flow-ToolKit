# Gaussian 局部粗糙元 Fourier 重构验证

参数：`R=285.36`, `alpha_f=0.38482`, `beta_f=0.07759`, `omega_f=-0.07759`, `A_inc=0.001 + 0.0im`, `epsilon=0.001`, `r_f=0.0`, `c=1.0`, `N=29`, `zmax=20.0`。每个单 Fourier 响应均使用六维 Modal 二阶求解器，目标波数为 `alpha_hf=alpha_f+kappa`。

## Fourier 变换验证
- 解析 Gaussian 峰值：`1`；数值逆变换峰值：`0.9999999980133817`，峰值位置：`0.0`。
- 最大绝对重构误差：`1.9866183142980276e-9`。逆变换使用 `1/(2*pi)`，正负 kappa 均参与。

## 单波数响应与完整局部响应
- [roughness_spectrum.csv](roughness_spectrum.csv) 中的 `G_single` 是单位粗糙度 Fourier 分量的速度能量响应，`abs_Hhat_sq_G_single` 是 Gaussian 谱权重后的传递谱，不是某个固定 r 处的局部能量。
- [reconstructed_response.csv](reconstructed_response.csv) 是先保留复数相位、再按 r 重构的完整响应；[local_energy.csv](local_energy.csv) 是重构后再计算的 `E_hf(r)`。
- 物理径向积分能量：`1.3980078155709509e-13`；离散 Parseval 能量：`1.3980078155709519e-13`；相对差：`7.222720412836265e-16`。
- 壁面强迫一致性最大相对误差：`1.6560117006689774`，比较的是重构 `q_hf(r,z=0)` 与 `-h(r) q_f,z(0)`。

## 重新解释单 kappa 峰值
此前 `|kappa|≈0.185` 的结果是单个单位 Fourier 粗糙度分量的传递响应峰值。完整 Gaussian 粗糙元的谱贡献由 `|Hhat(kappa)|^2 G_single(kappa)` 决定；因此峰值还要与 Gaussian 频谱宽度共同作用，不能直接把 `G_single` 的峰值当成局部粗糙元总能量峰值。
负、正 kappa 响应在固定复数入射涡波下不必共轭；实值 Gaussian 的几何谱具有共轭对称性，但耦合响应还包含固定 q_f、目标算子和相位因子。

## 重构收敛
见 [reconstruction_convergence.csv](reconstruction_convergence.csv)。计算采用均匀 kappa 网格，`L_r=2*pi/Delta-kappa` 的周期径向区间，因此物理径向积分与离散 Parseval 比较使用同一周期区间。
