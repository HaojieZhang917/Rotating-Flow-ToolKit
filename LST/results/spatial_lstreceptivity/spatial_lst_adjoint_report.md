# von Kármán 空间 LST 与伴随模态验证
参数：R=285.36, beta=0.07759, omega_f=-0.07759。
目标 omega_f 下的候选记录在 `spatial_lstr_eigenvalues.csv`；当前候选没有被自动标记为 Type I/II，且其谱分支随 N 变化明显，不能作为最终目标模态。
omega=0 参考算例用于验证直接/离散伴随接口，N=89 的模态和伴随模态已经保存。参考 alpha 在 N=49,69,89 间稳定到约 1e-11，说明二次铅笔离散求解本身收敛。
Q、投影归一化和离散 Green 检查见 `spatial_mode_normalization.csv`。当前 Q 的 N=49,69,89 数值分别约为 `0.0644+0.0236i`、`0.0691+0.0672i`、`0.0673+0.1046i`，尚未收敛，原因是当前计算使用离散转置伴随和单位双线性度量，尚未提供与连续物理内积一致的权重归一化。
直接与伴随剖面见 `spatial_mode_profiles.csv`；二进制模态见 `spatial_direct_modes.jls` 和 `spatial_adjoint_modes.jls`。
直接/伴随离散残差约为 1e-17，投影算子对左右特征向量缩放不变，离散转置 Green 检查为零。有限外缘 `zmax=20,30,40,50` 的参考 alpha 已稳定到约 1e-13，数据见 `spatial_zmax_convergence.csv`。当前尚未完成 Type I/II 的物理分类、半无限物理内积的连续 Green 验证和 Gaussian 感受性系数；不能将候选模态直接用于最终感受性结论。
