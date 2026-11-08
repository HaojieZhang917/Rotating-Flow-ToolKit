# LST：定转子腔转盘边界层行进波的感受性

研究问题：定转子腔中，静盘侧已经湍流化而转盘层仍是层流。静盘侧的非定常扰动与转盘上的粗糙元共同作用，怎样激发转盘层中的**行进**失稳波？（粗糙元单独只能激发定常模态。）

完整的物理主线、推导和结果解读见学习讲义 [`notes/rotor_stator_lecture/lecture.pdf`](notes/rotor_stator_lecture/lecture.pdf)（源文件 `lecture.md`，HTML 版 `lecture.html`）。

## 主线与对应代码

| 环节 | 内容 | 代码 | 结果目录 `results/` |
|---|---|---|---|
| 工具 | 单盘 von Kármán LST（验证算例），已修正二阶双线性项的 1/R 错误 | `VonKarmanLST.jl`、`solve_six_state_lns.jl`、`lab_stationary_core.jl`、`traveling_map_core.jl`、`local_mode_tracking.jl` | — |
| 工具 | 腔体模块：Batchelor 基流（含湍流静盘层的混合长 / SA 闭合）、转盘系空间算子（P_N–P_{N−2} 压力格式）、留数、二阶散射 | `CavityLST.jl` | — |
| 单盘验证 | 复现 Thomas (2024) 图 2：行进波中性曲线与增长率图 | `compute_traveling_growth_maps.jl`、`plot_traveling_growth_maps.py` | `traveling_growth_maps/` |
| 单盘验证 | 涡波 × 粗糙元的行进波感受性（复数 α 留数） | `compute_traveling_receptivity.jl`、`plot_traveling_receptivity.py` | `traveling_receptivity/` |
| 腔体验证 | 基流与算子移植验证（对照论文表 1） | `validate_cavity_port.jl` | `cavity_port_validation/` |
| 腔体验证 | 二阶散射：双线性项恒等式、复数失谐、网格收敛 | `validate_cavity_scattering.jl` | `cavity_scattering_validation/` |
| 初步探索 | 静盘扰动模型与穿透 | `demo_stator_forcing.jl`、`plot_stator_forcing_demo.py` | `stator_forcing_demo/` |
| 初步探索 | 直接路径 vs 粗糙元散射路径（平价粗糙元高度） | `compute_cavity_scattering_scan.jl`、`plot_cavity_scattering_scan.py` | `cavity_scattering_scan/` |
| 敏感性 | 湍流静盘层基流（混合长）及其验证 | `validate_turbulent_baseflow.jl`、`plot_turbulent_baseflow.py` | `turbulent_baseflow/` |
| 敏感性 | 两种闭合模型（混合长、SA）的 K(R) 对比 | `compare_turbulence_closures.jl`、`plot_turbulence_closures.py` | `turbulent_baseflow/closures.*` |
| 敏感性 | 三种基流下目标模态与两条路径的变化 | `compute_baseflow_sensitivity.jl`、`plot_baseflow_sensitivity.py` | `baseflow_sensitivity/` |
| **间隙滤波** | 推导文档、核心区检验、Λ_dir = 0.62·exp(−mΔζ)、惯性区与临界层 | `analyze_gap_filter.jl`、`analyze_gap_filter_core.jl`、`export_gap_filter_profiles.jl`、`plot_gap_filter.py` | `gap_filter/`（含 `gap_filter_theory.md`） |
| 对接转盘失稳 | 转盘层增长率图与间隙透过窗口 | `compute_cavity_rotor_map.jl`、`plot_cavity_rotor_map.py` | `cavity_rotor_map/` |

下一步：长波转换系数 Λ_u、Λ_v（见讲义第 15 章）。

## 运行说明

* Julia 只用到标准库（`LinearAlgebra`、`Printf`、`Distributed`）；Python 需要 `numpy`、`matplotlib`。
* 计算脚本把结果写入 `results/<目录>/`，对应的 `plot_*.py` 从同一目录读数据画图。大多数脚本几分钟内跑完；`compute_cavity_rotor_map.jl` 和 `compute_traveling_growth_maps.jl` 使用 `Distributed`（多进程，避免 BLAS 线程争用），约 10 分钟。
* 注意：`.gitignore` 忽略 `*.csv`、`*.png`、`*.pdf`，所以结果数据和图片只保存在本地磁盘上，不进 git。
* 腔体中所有受迫问题都必须用 `cavity_pencil(...; pressure=:PN2)`。默认的 `:PN` 格式只保留来对照原代码，它带有伪压力模态。

## 整理记录（2026-11-08）

* 从工作区 `hello-ce41a1` 同步了上表中的全部代码和结果，从 `laughing-booth-264191` 同步了讲义。
* 删除了 10 月底单盘"定常涡波 × 粗糙元"阶段的约 30 个探索脚本和结果，以及定常方案 A、β 扫描、smoke 测试、日志、修正前备份等过程文件。
  * 原先已提交的文件可以从提交 `3f05324` 恢复，例如 `git checkout 3f05324 -- LST/run_spatial_lst_receptivity.jl`。
  * 未提交的过程文件已打包备份到 `~/Rotating-Flow-ToolKit_archive/LST_untracked_backup_2026-11-08.tar.gz`。
* `LST_BEK.jl`、`ALST_BEK.jl`、`Vonkarmen.txt`、`test.ipynb` 是原有文件，未改动。
