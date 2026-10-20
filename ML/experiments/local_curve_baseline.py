"""Train-only local beta interpolation on each exact (R, alpha) grid pair.

This is a diagnostic baseline, not a neural model. It evaluates validation
only, never test labels, and does not modify prepared data or use the GPU.
"""

import argparse
import json
from pathlib import Path

import numpy as np


def metrics(prediction, target, mask):
    count = int(mask.sum())
    if count == 0:
        return {"count": 0}
    error = prediction[mask] - target[mask]
    return {
        "count": count,
        "rmse": np.sqrt(np.mean(error ** 2, axis=0)).tolist(),
        "mae": np.mean(np.abs(error), axis=0).tolist(),
        "max_absolute_error": np.max(np.abs(error), axis=0).tolist(),
    }


def run(path, degree):
    with np.load(path, allow_pickle=False) as dataset:
        x = dataset["x_raw"]
        y = dataset["y_raw"]
        train = dataset["train_idx"]
        validation = dataset["val_idx"]
    pairs, group = np.unique(x[:, :2], axis=0, return_inverse=True)
    group = group.reshape(-1)
    train = train[np.lexsort((x[train, 2], group[train]))]
    counts = np.bincount(group[train], minlength=len(pairs))
    offsets = np.concatenate(([0], np.cumsum(counts)))
    prediction = np.full((len(validation), 2), np.nan)
    interpolated = np.zeros(len(validation), dtype=bool)
    available = np.zeros(len(validation), dtype=bool)
    validation_order = np.argsort(group[validation], kind="stable")
    val_counts = np.bincount(group[validation], minlength=len(pairs))
    val_offsets = np.concatenate(([0], np.cumsum(val_counts)))
    nodes = degree + 1
    for gid in np.flatnonzero(val_counts):
        rows = validation_order[val_offsets[gid]:val_offsets[gid + 1]]
        source = train[offsets[gid]:offsets[gid + 1]]
        if len(source) < nodes:
            continue
        beta = x[source, 2]
        query = x[validation[rows], 2]
        insertion = np.searchsorted(beta, query)
        starts = np.clip(insertion - nodes // 2, 0, len(source) - nodes)
        selection = starts[:, None] + np.arange(nodes)[None, :]
        selected_beta = beta[selection]
        selected_y = y[source[selection]]
        weights = np.ones((len(rows), nodes))
        for i in range(nodes):
            for j in range(nodes):
                if i != j:
                    weights[:, i] *= (query - selected_beta[:, j]) / (
                        selected_beta[:, i] - selected_beta[:, j]
                    )
        prediction[rows] = np.sum(weights[:, :, None] * selected_y, axis=1)
        available[rows] = True
        interpolated[rows] = (query >= beta[0]) & (query <= beta[-1])
    extrapolated = available & ~interpolated
    return {
        "method": f"local degree-{degree} polynomial in beta at exact (R, alpha)",
        "train_count": len(train),
        "validation_count": len(validation),
        "unique_R_alpha_pairs": len(pairs),
        "available": metrics(prediction, y[validation], available),
        "interpolation_only": metrics(prediction, y[validation], interpolated),
        "endpoint_extrapolation": metrics(prediction, y[validation], extrapolated),
        "unavailable_count": int((~available).sum()),
        "notes": "Fixed train/validation indices; source targets are train-only; no test evaluation",
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dataset", type=Path, default=Path(__file__).resolve().parents[1] / "data/ts_0.0/dataset.npz")
    parser.add_argument("--degree", type=int, choices=(1, 2, 3, 4, 5), default=3)
    parser.add_argument("--report", type=Path)
    args = parser.parse_args()
    report = run(args.dataset, args.degree)
    result = json.dumps(report, indent=2)
    if args.report:
        args.report.write_text(result + "\n", encoding="utf-8")
    print(result)
