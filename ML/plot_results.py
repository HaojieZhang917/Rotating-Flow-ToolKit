"""Plot the saved learning curve and held-out predictions of one MLP run."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def case_label(metrics: dict) -> str:
    if metrics.get("Ts") is not None:
        return f"Ts = {metrics['Ts']:g}"
    values = metrics.get("Ts_values", [])
    return "Multiple Ts" + (": " + ", ".join(f"{value:g}" for value in values) if values else "")


def write_report(run: Path, metrics: dict, data_dir: Path | None = None) -> None:
    """Summarize the saved run; no model selection or new test prediction occurs."""
    config = metrics["configuration"]
    continuation = metrics.get("training_continuation")
    previously_tested = bool(metrics.get("test_split_previously_evaluated", False))
    metadata = json.loads(((data_dir or Path(config["data_dir"])) / "metadata.json").read_text(encoding="utf-8"))
    features = metrics.get("feature_names", metadata["input_names"])
    counts, split = metadata["counts"], metrics["split_sizes"]
    if "excluded_paired_sentinel" in counts:
        missing_description = (
            f"排除 {counts['excluded_missing_sentinel']:,} 条 `omega_i=-1` 缺失记录："
            f"其中成对 `(-1,-1)` 为 {counts['excluded_paired_sentinel']:,} 条，"
            f"仅虚部为 -1 为 {counts['excluded_omega_i_minus_one_only']:,} 条。"
        )
        breakdown = ["", "| Ts | 原始记录 | 成对缺失 | 仅虚部 -1 缺失 | 有效样本 |",
                     "|---|---:|---:|---:|---:|"]
        for value, group in sorted(metadata.get("counts_by_ts", {}).items(), key=lambda item: float(item[0])):
            breakdown.append(f"| {value} | {group['raw_rows']:,} | {group['excluded_paired_sentinel']:,} | {group['excluded_omega_i_minus_one_only']:,} | {group['usable_rows']:,} |")
        breakdown.append("")
    else:
        missing_description = f"历史数据排除 {counts['excluded_missing_sentinel']:,} 条成对 `(-1,-1)` 缺失记录，未应用最新的全部 `omega_i=-1` 缺失过滤；本报告不将其作为清洗后的新实验。"
        breakdown = []
    batches = (split["train"] + config["batch_size"] - 1) // config["batch_size"]
    stop = {"early_stopping": "验证指标触发提前停止", "target_rmse_reached": "达到两分量验证误差目标", "epoch_limit": "达到设定的最大训练轮数", "evaluate_only": "仅评估已有模型"}.get(metrics.get("stop_reason"), metrics.get("stop_reason", "验证损失触发提前停止" if metrics["early_stopped"] else "达到设定的最大训练轮数"))
    lines = [
        f"# {case_label(metrics)}：特征值代理模型实验", "",
        f"输入是 `({', '.join(features)})`，输出是 `(omega_r, omega_i)`。采用时间模式约定：`omega_i > 0` 表示增长。", "",
        "## 数据与划分", "",
        f"读取 {counts['source_files']} 个原始扫描文件，共 {counts['raw_rows']:,} 条记录。{missing_description}",
        f"另排除非有限记录 {counts.get('excluded_nonfinite', 0):,} 条，去除重复记录 {counts.get('duplicate_rows_removed', 0):,} 条。有效样本 {counts['usable_rows']:,} 条。除已确认的缺失标记外，真实负值保留。", "",
        *breakdown,
        "| 集合 | 样本数 |", "|---|---:|",
        f"| 训练 | {split['train']:,} |", f"| 验证 | {split['validation']:,} |", f"| 测试 | {split['test']:,} |", "",
        f"划分种子为 {metadata['seed']}。固定集合索引保存在 `dataset.npz`，各原始文件 SHA256 与样本行号一并记录。原始 R 以数据行数值为准。", "",
        "标准化均值和标准差仅由训练集计算：", "",
        "| 分量 | 训练集均值 | 训练集标准差 |", "|---|---:|---:|",
    ]
    for prefix, names in (("x", metadata["input_names"]), ("y", metadata["target_names"])):
        for name, mean, std in zip(names, metadata["scaler"][prefix + "_mean"], metadata["scaler"][prefix + "_std"]):
            lines.append(f"| {name} | {mean:.9g} | {std:.9g} |")
    if metrics.get("selection_candidates"):
        lines.extend(["", "两块 RTX 3080 各训练一个候选，二者使用相同的数据划分和训练种子；各完成 200 轮。候选比较只使用验证集：", "",
                      "| 隐藏层宽度 | 最佳轮次 | 验证 omega_r RMSE | 验证 omega_i RMSE |",
                      "|---|---:|---:|---:|"])
        for candidate in metrics["selection_candidates"]:
            errors = candidate["validation_rmse"]
            lines.append(f"| {' × '.join(map(str, candidate['hidden_sizes']))} | {candidate['best_epoch']} | {errors['omega_r']:.7g} | {errors['omega_i']:.7g} |")
        lines.extend(["", "选择标准为 `max(RMSE_r / 0.001, RMSE_i / 0.0001)` 最小，随后仅对选定模型评估测试集。", ""])
    if metrics.get("selection_targets"):
        targets = metrics["selection_targets"]
        status = "已达到" if metrics.get("target_met_on_validation") else "尚未达到"
        lines.extend(["", f"此前暂定的验证目标为实部 RMSE ≤ {targets['omega_r']:g}、虚部 RMSE ≤ {targets['omega_i']:g}；当前模型**{status}这两个目标**。", ""])
    if continuation:
        resumed_epoch = continuation["resumed_from_epoch"]
        additional = continuation["additional_epochs_completed"]
        trajectory_total = continuation["trajectory_epochs_total"]
        discarded = ", ".join(map(str, continuation.get("discarded_source_epochs", [])))
        training_description = (
            f"本次从 `{continuation['source_run']}` 第 {resumed_epoch} 轮 checkpoint 恢复，"
            f"新增训练 **{additional}** 轮，训练轨迹到总编号 **{trajectory_total}**。"
            f"选取总编号 **{metrics['best_epoch']}** 的模型"
            f"（续训第 {metrics['best_epoch'] - resumed_epoch} 轮），"
            f"验证 MSE 为 **{metrics['best_validation_loss']:.7g}**。停止原因：{stop}。"
        )
        lines.extend([
            "", "## 续训记录", "",
            f"原候选实际执行 {continuation['source_epochs_executed']} 轮；本次曲线采用原前 {continuation['source_history_epochs_used']} 轮和本次新增 {additional} 轮。"
            f"原总编号 {discarded} 的试跑被回退，保留为历史记录，不混入本次曲线。", "",
            f"本次最多新增 {continuation['max_additional_epochs']} 轮；早停耐心为 {continuation['early_stopping_patience']} 轮。"
            f"选模评分为 `{config['selection_metric']}`；连续相应轮数没有超过 `{config['min_delta']:g}` 的评分改善时早停。"
            "任何更低评分仍可保存为最佳 checkpoint。两个分量的验证误差都达到目标时提前结束。", "",
        ])
        if continuation.get("initial_validation_rmse") and "validation" in metrics:
            initial = continuation["initial_validation_rmse"]
            current = metrics["validation"]["mlp"]["components"]
            lines.extend([
                "| 阶段 | 验证 omega_r RMSE | 验证 omega_i RMSE |",
                "|---|---:|---:|",
                f"| 续训前（第 {resumed_epoch} 轮） | {initial['omega_r']:.7g} | {initial['omega_i']:.7g} |",
                f"| 本次最佳（总编号 {metrics['best_epoch']}） | {current['omega_r']['rmse']:.7g} | {current['omega_i']['rmse']:.7g} |", "",
            ])
    else:
        training_description = f"完成 {metrics['epochs_completed']} 轮，停止原因：{stop}。选取第 **{metrics['best_epoch']}** 轮的模型，验证 MSE 为 **{metrics['best_validation_loss']:.7g}**。"
    test_description = "以下指标由最终选定的模型在固定测试集上计算；两个参考基线仅使用训练集拟合。MAE、RMSE 均已还原到原始数据尺度。"
    if previously_tested:
        test_description += " 本次使用上一阶段已经查看过的同一固定测试集再次评估。训练与选模仍只依据验证集；这些结果不作为首次完全未见数据的一次性终验。"
    lines.extend([
        "", "## 训练", "",
        f"网络结构：`{len(features)} -> {' -> '.join(map(str, config['hidden_sizes']))} -> 2`，共 {config['n_parameters']:,} 个可调参数。隐藏层 Tanh，输出层线性。",
        f"使用 Adam，学习率配置 `{config['lr']}`，损失为两个标准化输出的平均 MSE。训练种子 {config['seed']}，运行设备 `{config.get('resolved_device', config['device'])}`。", "",
        f"PyTorch `{config.get('torch_version', '未记录')}`，NumPy `{config.get('numpy_version', '未记录')}`。", "",
        f"每个 epoch 使用全部 {split['train']:,} 条训练样本；batch size={config['batch_size']:,}，每轮 {batches} 次参数更新，最后一批不丢弃。每轮另评估 {split['validation']:,} 条验证样本。",
        training_description,
        f"本报告记录的训练与评估耗时约 {metrics['runtime_seconds']/60:.2f} 分钟。", "",
        "![训练与验证损失](learning_curve.png)", "",
        "## 固定测试集再次评估" if previously_tested else "## 独立测试集结果", "",
        test_description, "",
        "| 模型 | 标准化 MSE | omega_r MAE | omega_r RMSE | omega_i MAE | omega_i RMSE |",
        "|---|---:|---:|---:|---:|---:|",
    ])
    for key, label in (("mlp", "MLP"), ("linear_baseline", "线性最小二乘"), ("mean_baseline", "训练集均值")):
        result = metrics["test"][key]
        r, i = result["components"]["omega_r"], result["components"]["omega_i"]
        lines.append(f"| {label} | {result['standardized_mse']:.7g} | {r['mae']:.7g} | {r['rmse']:.7g} | {i['mae']:.7g} | {i['rmse']:.7g} |")
    if metrics.get("by_ts"):
        lines.extend(["", "各 Ts 分组结果：", "",
                      "| Ts | 集合 | 样本数 | omega_r RMSE | omega_i RMSE | 虚部符号误差率 |",
                      "|---|---|---:|---:|---:|---:|"])
        for value, evaluations in sorted(metrics["by_ts"].items(), key=lambda item: float(item[0])):
            for split_name, label in (("validation", "验证"), ("test", "测试")):
                if split_name not in evaluations:
                    continue
                group = evaluations[split_name]
                real, imag = group["components"]["omega_r"], group["components"]["omega_i"]
                growth = group["growth_rate_sign"]
                sign_rate = f"{100 * growth['error_rate']:.3f}%" if growth.get("component") == "omega_i" or "omega_i" in growth.get("definition", "") else "旧判据无效"
                lines.append(f"| {value} | {label} | {group['n_samples']:,} | {real['rmse']:.7g} | {imag['rmse']:.7g} | {sign_rate} |")
    result = metrics["test"]["mlp"]
    lines.extend(["", "| MLP 分量 | R² | 最大绝对误差 |", "|---|---:|---:|"])
    for name in metrics["target_names"]:
        component = result["components"][name]
        r2 = f"{component['r2']:.7g}" if component["r2"] is not None else "未定义"
        lines.append(f"| {name} | {r2} | {component['max_abs_error']:.7g} |")
    sign = result["growth_rate_sign"]
    if "omega_i" in sign.get("definition", "") or sign.get("component") == "omega_i":
        sign_description = f"用虚部正负检查目标模态稳定侧/不稳定侧。测试集符号误差为 {sign['error_count']:,}/{result['n_samples']:,}，误差率 **{100*sign['error_rate']:.3f}%**；该数值包含接近零的标签。"
    else:
        sign_description = "历史指标按实部符号计算，不能用于当前 `omega_i` 增长率约定；此处不将旧百分比解释为增长率符号误差。特征值回归误差仍可使用。"
    command = f".venv-wsl/bin/python predict.py --model {run}/best_model.pt" + (" --Ts 0.0" if "Ts" in features else "") + " --R 300 --alpha 0.2 --beta -0.1"
    lines.extend([
        "", sign_description,
        "", "测试结果仍包含较大的局部误差，最大绝对误差见上表。需要结合参数位置检查这些样本，不能仅凭总体 RMSE 判断每个参数点的精度。",
        "", "![测试预测与参考值](test_parity.png)", "",
        "## 结果文件与复现", "",
        "- `best_model.pt`：最佳网络、输入列顺序、适用 Ts、标准化参数与数据版本摘要。",
        "- `history.csv`：每轮训练/验证损失。",
        "- `metrics.json`：模型配置、验证/测试完整指标及基线结果。",
        "- `test_predictions.npz`：原始参数、真实值、预测值、误差和数据集索引。", "",
        "在 WSL 项目目录中进行单点预测：", "", "```bash",
        command, "```", "",
        "准备阶段记录原始文件 SHA256 和样本来源；模型选择仅使用验证集。", "",
        "本次结果针对已有 Ts 工况内的随机留出样本。未参与训练的 Ts、完整参数区域以及中性边界附近的精度需要另行验证。", "",
    ])
    if not metrics["early_stopped"] and metrics["best_epoch"] == metrics.get("last_epoch", metrics["epochs_completed"]):
        lines.extend(["本轮最佳模型来自最后一个 epoch，验证曲线仍在下降，尚不能据此认为训练已经充分收敛。", ""])
    (run / "report.md").write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    run_root = Path(__file__).resolve().parent / "runs" / "all_ts_clean"
    continued = run_root / "final_continued"
    default_run = continued if (continued / "report.md").is_file() else run_root / "final"
    parser.add_argument("--run-dir", type=Path, default=default_run)
    parser.add_argument("--data-dir", type=Path, help="Override the prepared-data path saved in the run (e.g. when switching Windows/WSL)")
    args = parser.parse_args()
    run = args.run_dir.resolve()
    metrics = json.loads((run / "metrics.json").read_text())
    if "test" not in metrics or not (run / "test_predictions.npz").is_file():
        parser.error("This plot requires a finalized run with saved test predictions; it does not evaluate the test set itself")
    with (run / "history.csv").open() as handle:
        history = list(csv.DictReader(handle))
    epochs = [int(row["epoch"]) for row in history]
    train = [float(row["train_loss"]) for row in history]
    validation = [float(row["validation_loss"]) for row in history]
    fig, axis = plt.subplots(figsize=(7.0, 4.2), constrained_layout=True)
    axis.semilogy(epochs, train, label="Training", color="#2563eb")
    axis.semilogy(epochs, validation, label="Validation", color="#e87524")
    axis.axvline(metrics["best_epoch"], color="#555555", ls="--", lw=1, label=f"Selected epoch {metrics['best_epoch']}")
    if metrics.get("training_continuation"):
        resumed_epoch = metrics["training_continuation"]["resumed_from_epoch"]
        axis.axvline(resumed_epoch, color="#8b5cf6", ls=":", lw=1.2, label=f"Resumed from epoch {resumed_epoch}")
    axis.set(xlabel="Epoch", ylabel="MSE of standardized outputs", title=f"MLP: {case_label(metrics)}")
    axis.grid(alpha=0.2)
    axis.legend()
    fig.savefig(run / "learning_curve.png", dpi=170)
    plt.close(fig)

    with np.load(run / "test_predictions.npz", allow_pickle=False) as data:
        true, predicted = data["y_true"], data["y_pred"]
    # Keep plots readable. All held-out rows, not this subset, determine metrics.
    sample = np.random.default_rng(42).choice(len(true), min(10000, len(true)), replace=False)
    fig, axes = plt.subplots(1, 2, figsize=(9.5, 4.2), constrained_layout=True)
    labels = [r"$\omega_r$", r"$\omega_i$"]
    for column, axis in enumerate(axes):
        axis.scatter(true[sample, column], predicted[sample, column], s=5, alpha=0.2, edgecolors="none", rasterized=True)
        limits = [float(min(true[:, column].min(), predicted[:, column].min())),
                  float(max(true[:, column].max(), predicted[:, column].max()))]
        axis.plot(limits, limits, color="#d1495b", lw=1.2)
        component = metrics["test"]["mlp"]["components"][["omega_r", "omega_i"][column]]
        axis.set(xlabel=f"Reference {labels[column]}", ylabel=f"Predicted {labels[column]}",
                 title=f"Test RMSE = {component['rmse']:.3g}")
        axis.grid(alpha=0.2)
    fig.suptitle(f"{case_label(metrics)}; 10,000 displayed test points at most")
    fig.savefig(run / "test_parity.png", dpi=170)
    plt.close(fig)
    write_report(run, metrics, args.data_dir)
    print(f"Saved plots and report.md in {run}")


if __name__ == "__main__":
    main()
