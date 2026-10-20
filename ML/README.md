# 转子边界层特征值代理模型

使用普通全连接神经网络学习指定目标模态的特征值：

`(Ts, R, alpha, beta) -> (omega_r, omega_i)`

`Ts` 表示壁面吹吸条件。时间模式采用 **`omega_i > 0` 为增长模态、`omega_i < 0` 为衰减模态、`omega_i = 0` 为中性** 的约定。直接使用文件中的实部、虚部作为标签，不交换两列、不改变符号，也没有加入物理残差约束。目标模态的预测不能直接解释为整个系统的稳定性或实际转捩。

此前 `Ts=0.0` 的三输入实验继续保留。旧报告曾按实部符号计算稳定性误差；该项不适用于已纠正的物理约定。此前仅排除成对 `(-1,-1)` 的全 Ts 数据备份在 `data/all_ts_paired_only`，对应 `runs/all_ts` 标记为 `SUPERSEDED`。这些旧实验不代表按最新规则清洗后的模型，不与新实验直接比较误差。

## 数据来源与读取规则

原始目录保持只读：

```text
/home/zhj/main/code/compress/compressible_bek/TwoDisk/data/stability/Zarf_neutral
```

- `--all-ts` 读取 `Ts=...` 子目录中的 `R=*.dat` 原始扫描文件；汇总和提取文件不重复计入。`Ts` 从目录读取，其余输入从数据行读取。
- 原始表头仍是 `R, alpha, beta, omega_r, omega_i` 五列，并验证 ZONE 点数。文件名中的 R 不一定等于行内 R，以行内值为准。
- **任何 `omega_i == -1` 的记录均为缺失，必须排除**，包括 `(-1,-1)` 和实部有效但虚部为 -1 的记录。仅 `omega_r == -1`、`omega_i != -1` 时仍保留。其他真实负值保留；非有限记录另行排除。
- 四参数实验以 `(Ts, R, alpha, beta)` 为去重键。相同输入、相同标签保留第一条；相同输入出现冲突标签则报错。不同 Ts 下相同三维参数是不同样本。
- 保存文件 SHA256、每条样本的源文件编号及原始行号，原始扫描文件不改写。

每个 Ts 内按确定的文件/行顺序组织有效样本，再独立使用种子 `42` 生成排列。每组训练数为 `floor(0.8N)`、验证数为 `floor(0.1N)`，其余为测试。各组索引合并后固定，训练时仅打乱训练样本的读取顺序。重复运行相同过滤和种子会得到相同划分；此次过滤条件改变后，旧实验的划分不能直接沿用。

`excluded_missing_sentinel` 是所有 `omega_i=-1` 缺失记录的总数，`excluded_paired_sentinel` 和 `excluded_omega_i_minus_one_only` 分别记录其中的成对缺失和仅虚部 -1 缺失；整体及每 Ts 均单独统计。

四个输入和两个输出均使用**合并训练集**的均值、总体标准差（`ddof=0`）标准化；验证和测试沿用这些统计量。数据保存为原始 `float64` 和标准化 `float32`，使用 `numpy.load(path, allow_pickle=False)` 读取。

这项划分评价已有 Ts 工况内的留出样本。未见 Ts 的泛化、参数空间外推和中性边界精度需要专门验证。

## 四输入主流程

数据准备和绘图可在 WSL 使用已有 `.venv-wsl`。实际 GPU 训练环境为 Windows `D:\python\python.exe`，PyTorch `1.8.1+cu111`、NumPy `1.20.1`，设备为两块 RTX 3080。配置和依赖版本记录在训练产物中。

```bash
cd /home/zhj/Rotating-Flow-ToolKit/ML
.venv-wsl/bin/python prepare_data.py --all-ts --seed 42
```

两个候选可分别在 Windows PowerShell 窗口运行：

```powershell
Set-Location '\\wsl.localhost\Ubuntu-24.04\home\zhj\Rotating-Flow-ToolKit\ML'
D:\python\python.exe train_mlp.py --data-dir data/all_ts --output-dir runs/all_ts_clean/baseline --seed 42 --hidden-sizes 128 128 128 --epochs 200 --device cuda:1 --resident-data --batch-size 8192 --eval-batch-size 32768 --threads 4 --lr 0.001 --lr-patience 8 --lr-factor 0.5 --min-lr 0.000001 --patience 40 --min-delta 0.0001 --target-rmse-r 0.001 --target-rmse-i 0.0001 --sign-tolerance 0.0001 --validation-only
D:\python\python.exe train_mlp.py --data-dir data/all_ts --output-dir runs/all_ts_clean/candidate_256 --seed 42 --hidden-sizes 256 256 256 --epochs 200 --device cuda:0 --resident-data --batch-size 8192 --eval-batch-size 32768 --threads 4 --lr 0.001 --lr-patience 8 --lr-factor 0.5 --min-lr 0.000001 --patience 40 --min-delta 0.0001 --target-rmse-r 0.001 --target-rmse-i 0.0001 --sign-tolerance 0.0001 --validation-only
```

两个 3080 各自训练一个候选模型，是独立实验，并非把同一个模型分配到两张卡。以上选择宽度 128 和 256 的三隐藏层 Tanh 网络，输出层线性。损失为两个标准化输出的平均 MSE，优化器为 Adam。`--resident-data` 将使用的集合放在对应 GPU 上。每个 epoch 遍历全部训练样本，最后不足一批的样本也使用；每轮参数更新次数为 `ceil(N_train / 8192)`。

每轮只使用验证集选择模型，候选以 `--validation-only` 运行。最终按验证表现选出模型，将其 checkpoint 用 `--evaluate-only --resume-model ...` 放到 `runs/all_ts_clean/final` 评估固定测试集；不要根据测试误差反过来选择候选。

```bash
.venv-wsl/bin/python plot_results.py --run-dir runs/all_ts_clean/final --data-dir data/all_ts
```

数据产物或运行目录已存在时，脚本拒绝覆盖。比较新模型时保留同一数据目录并指定新的输出目录；只有有意替换产物时才使用 `--overwrite`。四输入模型不能直接从旧三输入 checkpoint 恢复训练。

`train_mlp.py` 默认数据目录是 `data/all_ts`，默认输出目录是 `runs/all_ts_clean/baseline`。网络宽度、层数和训练轮数是实验设置，不代表已达到特定误差目标。原暂定的两个分量容差也不会因增长率约定修正而交换。

单点预测：

```powershell
D:\python\python.exe predict.py --model runs/all_ts_clean/final/best_model.pt --Ts 0.0 --R 300 --alpha 0.2 --beta -0.1
```

四输入 checkpoint 必须提供 `--Ts`。未指定模型时优先读取 `runs/all_ts_clean/final/best_model.pt`，尚无 final 时才回退到新 baseline。输出 JSON 包含四个输入、还原后的 `omega_r`、`omega_i` 及坐标范围提示。预测不依赖原始 `.dat` 文件。各列落在最小和最大值之间，也不能证明该参数组合附近已有样本。

## 原单 Ts 实验

`prepare_data.py --ts 0.0` 仍生成三输入 `(R, alpha, beta)` 数据，将 Ts 保存为固定工况。未指定 `--all-ts` 时保留这一行为。

```bash
.venv-wsl/bin/python prepare_data.py --ts 0.0 --seed 42 --output-dir data/ts_0.0
.venv-wsl/bin/python train_mlp.py --data-dir data/ts_0.0 --output-dir runs/ts_0.0/new_run --hidden-sizes 64 64
.venv-wsl/bin/python predict.py --model runs/ts_0.0/baseline/best_model.pt --R 300 --alpha 0.2 --beta -0.1
```

单 Ts 模型可省略 `--Ts`；如提供，必须与 checkpoint 的固定 Ts 一致。上述准备命令不会覆盖已存在的数据。

历史 `Ts=0.0` 数据采用旧缺失规则，旧样本数不再作为当前实验依据。新规则下的整体数量、缺失分项及每个 Ts 的划分见 `data/all_ts/metadata.json`。旧三输入 checkpoint 仍可读取，其结果只用于历史复现。

## 结果文件

| 文件 | 内容 |
|---|---|
| `data/all_ts/dataset.npz` | 原始/标准化样本、固定索引、样本来源 |
| `data/all_ts/scaler.npz` | 合并训练集的输入、输出均值和标准差 |
| `data/all_ts/metadata.json` | Ts 列表、整体和分组数量、划分规则、原始文件摘要 |
| `best_model.pt` | 网络结构和权重、输入顺序、适用 Ts、标准化参数、数据摘要 |
| `history.csv` | 每轮训练、验证损失和两分量验证 RMSE |
| `metrics.json` | 配置、验证和最终测试指标、参考基线 |
| `test_predictions.npz` | 原单位输入、真实值、预测值、误差、原样本索引 |
| `learning_curve.png`、`test_parity.png`、`report.md` | 由绘图脚本读取已保存结果生成 |

新候选结果位于 `runs/all_ts_clean/baseline` 和 `runs/all_ts_clean/candidate_256`，最终选定模型和测试结果位于 `runs/all_ts_clean/final`。旧 `runs/all_ts` 已作废但保留用于溯源。

本轮已完成：排除 1,061,227 条缺失记录，其中仅虚部为 -1 的额外记录为 123 条；有效样本 3,042,808 条。训练/验证/测试分别为 2,434,244 / 304,278 / 304,286，每轮训练 298 个 batch（最后一批 1,220 条）。两个候选各训练 200 轮，按验证指标选中 `4 -> 256 -> 256 -> 256 -> 2` 网络第 197 轮的 checkpoint。独立测试 RMSE：实部 **0.00304905**、虚部 **0.000363216**，尚未达到暂定的 0.001 / 0.0001 目标。

`runs/all_ts_clean/selection.json` 保存候选比较；`final/training_metrics.json` 和 `final/evaluation_metrics.json` 分别保留训练与最终评估的原始记录，`final/metrics.json` 汇总两者。最终 `best_model.pt` 是原选中 checkpoint 的逐字节副本，`history.csv` 保留该候选的完整 200 轮记录。五个已有 Ts 内的随机留出结果不能证明未见 Ts 的预测精度。

指标包括两分量的 MAE、RMSE、最大绝对误差、R²，以及按 `omega_i` 计算的增长率符号误差。均值预测和线性最小二乘基线只使用训练集拟合。整体误差需结合各 Ts 误差及局部参数覆盖判断。

`plot_results.py --data-dir` 可覆盖训练配置保存的数据路径，便于在 Windows 和 WSL 间读取同一运行结果。绘图脚本只读取已完成的测试预测，不会重新评估测试集。`diagnose_validation.py` 是原 `Ts=0.0` 三维数据的专用诊断，不适用于全 Ts 四维数据。

## 检查

```bash
.venv-wsl/bin/python -m unittest -v test_pipeline
```

测试使用临时合成数据，检查成对缺失和仅虚部 -1 的缺失过滤、仅实部 -1 和其他真实负值保留、重复与冲突标签、行内 R、固定分层划分、训练集标准化、三输入和四输入预测还原，以及虚部增长率判据；不改动原始扫描文件。

