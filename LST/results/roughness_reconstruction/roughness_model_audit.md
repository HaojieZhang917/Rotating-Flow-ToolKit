# 粗糙度模型核查

核查对象：`VonKarmanLST.jl`、六维二阶调用入口以及历史扫描脚本。

|检查对象|实际实现|源代码位置|
|---|---|---|
|粗糙度物理形式|单一复 Fourier 分量 `h= h_amplitude*exp(i*(kappa*r+n*theta))`；不是当前求解器内部直接求解的局部 Gaussian|`VonKarmanLST.jl:1037-1087`|
|kappa 定义|单一径向 Fourier 波数；定常畸变目标 `alpha_h=kappa`，二阶目标 `alpha_hf=alpha_f+kappa`|`VonKarmanLST.jl:1049-1067,1148-1174`|
|方位向波数|`n=0` 时 `beta_h=0`、`beta_hf=beta_f`；`n!=0` 必须显式提供 `target_beta` 或 `beta_shift`，代码不会自动把整数周向模态换算成局部 beta|`VonKarmanLST.jl:1059-1069,1003-1027`|
|qh 定义|单位粗糙度 Fourier 分量诱导的定常畸变；主扫描使用 `h_amplitude=1`|`VonKarmanLST.jl:1037-1087`；扫描调用见 `scan_roughness_wavenumber.jl:80-82`|
|体积源|单一 kappa 的交叉源 `-(B(qh,qf)+B(qf,qh))`；连续性和压力源为零|`VonKarmanLST.jl:1148-1174`|
|壁面强迫|`-h_amplitude*qf_z(0)`；二阶调用由 `wall_shear_response(qf)` 构造|`VonKarmanLST.jl:1203-1205`；六维调用见 `solve_six_state_lns.jl:340-342`|
|epsilon|单一 Fourier 分支求解中通过 `roughness_amplitude` 乘入；qh 本身不含主扫描的 `EPSROUGH`|`VonKarmanLST.jl:1188-1205`|
|A_inc|包含在 qf 的一阶响应中；二阶源因此线性含 `A_inc`，最终增益再除以 `epsilon^2*abs(A_inc)^2`|扫描脚本 `response_data` 和 `refine_negative_branch_peak.jl:33-37`|
|正负分支|调用端使用有符号 `kappa`；负分支为 `kappa_signed=-abs(kappa)`，qh、源项和目标算子共同使用该符号|`scan_incident_roughness_wavenumbers.jl:84-88`；`fixed_target_matching.jl:24-27`|
|已有 Gaussian 函数|已存在解析包络、解析 Fourier 谱、均匀 kappa 网格和复数场重构函数，但此前扫描脚本没有调用它们完成六维 Modal 全响应积分|`VonKarmanLST.jl:846-987`；仅有旧验证脚本 `validate_coupled_receptivity.jl:89-160`|
|当前 G_full|固定 `[0,10]` 上的速度平方能量除以 `epsilon^2*abs(A_inc)^2`；它是单 Fourier 分支传递增益，不是局部粗糙元总响应，也不是 Type I/II 伴随激发系数|`compare_second_order_outer_bc.jl:27-32`|
|当前粗糙度是否已完整局部化|没有。旧的 `reconstruct_roughness_field` 只对给定的 qhat 数组做谱叠加，历史主扫描逐个 kappa 求解并比较单分量增益|`VonKarmanLST.jl:913-987`；`scan_roughness_wavenumber.jl:78-108`|

## 结论

此前的 `|kappa_peak|≈0.185` 和 `G_full≈0.23051` 表示单位粗糙度 Fourier 分量、固定入射涡波和固定局部 Modal/Dirichlet 数值设置下的单分量速度能量传递响应。它可以作为传递函数的波数选择性结果，但不能直接解释为有限宽度局部粗糙元的总散射能量。
要得到局部 Gaussian 粗糙元，需要将单位 Fourier 响应乘以 `epsilon*Hhat(kappa)*Delta-kappa/(2*pi)`，保留每个分量的复相位和目标波数，然后在物理径向位置先叠加复数响应，再计算局部能量。

本次静态调用检查得到：`qh` 的单位幅值为 1.0 + 0.0im，二阶目标 `alpha_hf=0.19982`，源向量长度为 120。
