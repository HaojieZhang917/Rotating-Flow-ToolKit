# Vonkarmen_bone Daily Workspace

该目录用于日常计算、Jupyter Notebook、结果检查和论文数据后处理。
可发行的 Julia 包位于：

```text
/home/zhj/Rotating-Flow-ToolKit/RotatingDiskFlow
```

## Notebook 推荐载入方式

```julia
using Pkg
Pkg.activate("/home/zhj/Rotating-Flow-ToolKit/RotatingDiskFlow")
using RotatingDiskFlow
```

已有 notebook 可以继续使用原来的兼容入口：

```julia
include("NeutralCurveRunner.jl")
using .NeutralCurveRunner
```

这些根目录兼容文件会自动加载 `../../../RotatingDiskFlow/src/` 中的发行版实现，因此旧 notebook 不需要批量改写。

## 日常输出

从本目录启动 Julia 时，中性曲线默认写入：

```text
Vonkarmen_bone/neutral_curve_batch/
```

发行包不会把运行数据写入自己的 `src/` 或仓库目录。该目录中的 notebook、历史分析脚本和计算结果不属于发行包 API。

完整使用文档：

```text
/home/zhj/Rotating-Flow-ToolKit/RotatingDiskFlow/README.md
/home/zhj/Rotating-Flow-ToolKit/RotatingDiskFlow/docs/neutral-curves.md
```
