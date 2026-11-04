# 空间 LST 归一化与模态留数验证

回归参数：R=285.36, beta=0.07759, omega=0.0, zmax=30.0。直接模态按 `max|v|=1` 且峰值相位设为零；离散伴随按 `max|u_adj|=1` 且峰值相位设为零。

## 验证公式

对 `M(alpha)q=b`，使用 `Q=y^T M'(alpha_n)x`，预测 `lim_(delta->0) delta*q(alpha_n+delta)=x*(y^T*b)/Q`。体积强迫是 u 动量方程中的 Gaussian z 源；壁面强迫通过平滑提升函数转换为与边界行一致的离散右端。

## 结果

完整数据见 [residue_validation.csv](residue_validation.csv)。`vol_d1..d4`、`wall_d1..d4` 和 `full_d1..d4` 分别对应 `delta=1e-2,3e-3,1e-3,3e-4` 的相对留数误差。

本验证使用的是同一离散 QEP 的转置左零向量，因此验证对象是离散模态留数。它没有把该左向量宣称为连续 Hermitian 伴随，也没有把 Clenshaw–Curtis 权重强行加入 Q。

## 物理伴随尚未完成的部分

当前半无限/有限映射下的物理 Hermitian 权重、连续 Green 边界泛函及其与离散左零向量的严格映射仍需单独完成。Gaussian 粗糙度的最终 Type I/II 感受性必须在这一映射验证后计算。
