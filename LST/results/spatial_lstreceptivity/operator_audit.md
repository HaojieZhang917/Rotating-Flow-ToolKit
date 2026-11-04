# 空间 LST 与离散伴随算子核查
当前变量为 `[u,v,w,p]`，空间算子为 `L0+alpha*L1+alpha^2*L2`，Fourier 约定为 `exp(i*(alpha*r+beta*s-omega*t))`。
当前离散伴随是施加相同边界约束后的二次铅笔逐项转置；`Q` 使用 `transpose(y)*(L1+2alpha*L2)*x`。
当前结果是离散双线性伴随框架，后续若采用物理 Hermitian 伴随，需要补充一致的权重矩阵和连续 Green 边界验证。
