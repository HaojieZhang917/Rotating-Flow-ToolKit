# 测试与回归基准

## 1. 快速测试

```bash
cd /home/zhj/Rotating-Flow-ToolKit/compress/BEK/Vonkarmen_bone
julia --project=. test/runtests.jl
```

快速测试当前包含 48 个断言：

- 包入口、配置检查和输出文件命名；
- Tecplot 曲线写入、读取和合法性检查；
- 自适应 `beta` 步长的缩小与恢复；
- Chebyshev 一阶、二阶导数的多项式精确性；
- Sutherland 黏性率、密度和 BDF2 公式；
- Lopez 算子对原非线性方程有限差分 Jacobian 的一致性；
- 静止/旋转坐标频率变换；
- 空间与时间算子的一致性；
- 壁面和远场边界自由度消元。
- 根目录旧 `include(...)` 入口的兼容性。

## 2. 物理回归

```bash
RUN_PHYSICS_REGRESSION=true julia --project=. test/runtests.jl
```

该测试额外调用 SciPy 基本流并验证 Malik (1986) 的两个等温中性点：

| 模态 | `R` | `beta` | 文献 `alpha_r` | 当前 `N=69` 结果 |
|---|---:|---:|---:|---:|
| Type I | 285.36 | 0.07759 | 0.38482 | `0.38505249 + 3.48666e-4im` |
| Type II | 440.88 | 0.04672 | 0.13228 | `0.13184955 + 9.07977e-4im` |

判据为：

```text
relative error(alpha_r) < 0.5%
abs(alpha_i) < 1.0e-3
```

这里的 `alpha_i` 是当前 rational 网格/远场截断在文献参数点上的离中性偏差，不是多项式特征值问题的代数残差。Type II 对远场和网格映射更敏感，因此使用项目原有的 `1e-3` 回归容差。

## 3. 网格无关性

```bash
julia --project=. scripts/grid_independence.jl
```

该脚本计算 `N=49,59,69,79,89` 下的 Type I 和 Type II 临界点，并比较插值到公共物理坐标后的特征函数。现有结果见[双模态网格无关性验证](grid-independence.md)。

## 4. 测试策略

普通开发先运行快速测试。修改以下内容时必须运行完整物理回归：

- 基本流方程或温度/密度定义；
- `CRD_STA.jl` 或 `LopezStability.jl` 任一矩阵块；
- 压力边界条件；
- Chebyshev 网格和相似坐标插值；
- 特征值筛选、特征向量追踪和中性根校正。
