#!/usr/bin/env python3
"""Prepare one Ts slice or all Ts values for eigenvalue regression.

Run in WSL: python3 prepare_data.py --ts 0.0 --seed 42
Use --all-ts for four inputs (Ts, R, alpha, beta). Without it, --ts selects
three inputs (R, alpha, beta). Raw inputs are read-only. Outputs contain
fixed split indices and training-only scalers. Growth is determined by omega_i.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
from pathlib import Path
import re
import sys
import tempfile

import numpy as np


DEFAULT_SOURCE_ROOT = Path(
    "/home/zhj/main/code/compress/compressible_bek/TwoDisk/data/stability/Zarf_neutral"
)
INPUT_NAMES = ["R", "alpha", "beta"]
ALL_INPUT_NAMES = ["Ts", *INPUT_NAMES]
TARGET_NAMES = ["omega_r", "omega_i"]
FIELD_NAMES = INPUT_NAMES + TARGET_NAMES
OUTPUT_NAMES = ("dataset.npz", "scaler.npz", "metadata.json")
VARIABLES_RE = re.compile(
    r'VARIABLES\s*=\s*"R"\s*,\s*"alpha"\s*,\s*"beta"\s*,\s*'
    r'"omega_r"\s*,\s*"omega_i"\s*',
    re.IGNORECASE,
)
TITLE_RE = re.compile(r'TITLE\s*=\s*"[^\"]*"\s*', re.IGNORECASE)
ZONE_RE = re.compile(
    r'ZONE\s+T\s*=\s*"[^\"]*"\s*,\s*I\s*=\s*(\d+)\s*,\s*'
    r'J\s*=\s*(\d+)\s*,\s*F\s*=\s*POINT\s*',
    re.IGNORECASE,
)


def ts_label(value: float) -> str:
    """Use the directory naming convention, including Ts=0.0."""
    if not math.isfinite(value):
        raise ValueError("Ts must be finite")
    return str(0.0 if value == 0.0 else value)


def read_source(path: Path, file_id: int) -> tuple[np.ndarray, np.ndarray, dict]:
    """Read one Tecplot POINT file and retain physical source line numbers."""
    rows = []
    line_numbers = []
    digest = hashlib.sha256()
    header_stage = 0
    seen_title = False
    expected_rows = None
    with path.open("rb") as source:
        for line_number, raw_line in enumerate(source, 1):
            digest.update(raw_line)
            encoding = "utf-8-sig" if line_number == 1 else "utf-8"
            try:
                line = raw_line.decode(encoding).strip()
            except UnicodeDecodeError as exc:
                raise ValueError(f"{path}:{line_number}: invalid UTF-8") from exc
            if not line:
                continue
            if header_stage == 0:
                if TITLE_RE.fullmatch(line) and not seen_title:
                    seen_title = True
                    continue
                if not VARIABLES_RE.fullmatch(line):
                    raise ValueError(
                        f"{path}:{line_number}: expected VARIABLES header with "
                        f"columns {FIELD_NAMES} in this order"
                    )
                header_stage = 1
                continue
            if header_stage == 1:
                zone = ZONE_RE.fullmatch(line)
                if zone is None:
                    raise ValueError(
                        f"{path}:{line_number}: expected ZONE T=..., I=..., "
                        "J=..., F=POINT header"
                    )
                i_count, j_count = (int(value) for value in zone.groups())
                if i_count < 1 or j_count < 1:
                    raise ValueError(f"{path}:{line_number}: ZONE dimensions must be positive")
                expected_rows = i_count * j_count
                header_stage = 2
                continue
            fields = line.split()
            if len(fields) != 5:
                raise ValueError(f"{path}:{line_number}: expected 5 columns, got {len(fields)}")
            try:
                rows.append([float(value) for value in fields])
            except ValueError as exc:
                raise ValueError(f"{path}:{line_number}: nonnumeric data row") from exc
            line_numbers.append(line_number)
    if header_stage != 2:
        raise ValueError(f"{path}: missing required VARIABLES or ZONE header")
    if len(rows) != expected_rows:
        raise ValueError(
            f"{path}: ZONE declares {expected_rows} rows but file contains {len(rows)}"
        )

    raw = np.asarray(rows, dtype=np.float64).reshape(-1, 5)
    source_lines = np.asarray(line_numbers, dtype=np.int64)
    # Confirmed by the data owner: any omega_i=-1 marks a missing result.
    # This includes the original (-1,-1) pair. Other negative values remain.
    paired_missing = (raw[:, 3] == -1.0) & (raw[:, 4] == -1.0)
    missing = raw[:, 4] == -1.0
    # Report mutually exclusive exclusion counts, with the sentinel taking priority.
    nonfinite = (~np.isfinite(raw).all(axis=1)) & ~missing
    keep = ~(missing | nonfinite)
    summary = {
        "file_id": file_id,
        "name": path.name,
        "sha256": digest.hexdigest(),
        "zone_rows": expected_rows,
        "raw_rows": len(rows),
        "excluded_missing_sentinel": int(missing.sum()),
        "excluded_paired_sentinel": int(paired_missing.sum()),
        "excluded_omega_i_minus_one_only": int((missing & ~paired_missing).sum()),
        "excluded_nonfinite": int(nonfinite.sum()),
        "valid_rows_before_deduplication": int(keep.sum()),
    }
    return raw[keep], source_lines[keep], summary


def deduplicate(
    raw: np.ndarray, source_file_ids: np.ndarray, source_lines: np.ndarray,
    source_summaries: list[dict], input_dim: int = 3,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, int]:
    """Keep the first exact duplicate; reject contradictory target labels."""
    _, first, inverse = np.unique(raw[:, :input_dim], axis=0, return_index=True, return_inverse=True)
    inverse = inverse.reshape(-1)
    conflicts = np.any(raw[:, input_dim:] != raw[first[inverse], input_dim:], axis=1)
    if conflicts.any():
        row = int(np.flatnonzero(conflicts)[0])
        original = int(first[inverse[row]])
        locations = []
        for index in (original, row):
            name = source_summaries[int(source_file_ids[index])]["name"]
            locations.append(f"{name}:{int(source_lines[index])}")
        raise ValueError(
            "Conflicting eigenvalue labels for the same input "
            f"{raw[row, :input_dim].tolist()}: {locations[0]} has "
            f"{raw[original, input_dim:].tolist()}, {locations[1]} has "
            f"{raw[row, input_dim:].tolist()}. Resolve the conflict before splitting."
        )
    # Restore file/line order rather than np.unique's lexicographic input order.
    keep_indices = np.sort(first)
    kept_per_source = np.bincount(source_file_ids[keep_indices], minlength=len(source_summaries))
    for file_id, summary in enumerate(source_summaries):
        kept = int(kept_per_source[file_id])
        summary["usable_rows"] = kept
        summary["duplicate_rows_removed"] = summary["valid_rows_before_deduplication"] - kept
    return (
        raw[keep_indices], source_file_ids[keep_indices], source_lines[keep_indices],
        len(raw) - len(keep_indices),
    )


def statistics(values: np.ndarray, names: list[str]) -> dict:
    return {
        name: {"min": float(values[:, index].min()), "max": float(values[:, index].max())}
        for index, name in enumerate(names)
    }


def training_scaler(values: np.ndarray, names: list[str]) -> tuple[np.ndarray, np.ndarray]:
    mean = values.mean(axis=0, dtype=np.float64)
    std = values.std(axis=0, ddof=0, dtype=np.float64)
    if not np.isfinite(mean).all() or not np.isfinite(std).all():
        raise ValueError(f"Nonfinite training statistics for {names}")
    if (std == 0.0).any():
        constant = [name for name, value in zip(names, std) if value == 0.0]
        raise ValueError(f"Cannot standardize constant training components: {constant}")
    return mean, std


def save_outputs(output_dir: Path, arrays: dict, scalers: dict, metadata: dict) -> None:
    """Finish all temporary files before replacing the three final artifacts."""
    output_dir.mkdir(parents=True, exist_ok=True)
    temporary_paths = []
    try:
        for name, contents in (("dataset.npz", arrays), ("scaler.npz", scalers)):
            with tempfile.NamedTemporaryFile(
                mode="w+b", prefix=f".{name}.", suffix=".tmp", dir=output_dir, delete=False,
            ) as temporary:
                temporary_paths.append((Path(temporary.name), output_dir / name))
                # Pass an open stream so NumPy does not append '.npz' to the path.
                np.savez_compressed(temporary, **contents)
                temporary.flush()
                os.fsync(temporary.fileno())
        with tempfile.NamedTemporaryFile(
            mode="w", encoding="utf-8", newline="\n", prefix=".metadata.json.",
            suffix=".tmp", dir=output_dir, delete=False,
        ) as temporary:
            temporary_paths.append((Path(temporary.name), output_dir / "metadata.json"))
            json.dump(metadata, temporary, ensure_ascii=False, indent=2, allow_nan=False)
            temporary.write("\n")
            temporary.flush()
            os.fsync(temporary.fileno())
        # Metadata is replaced last. This is atomic per file, not across all files.
        for temporary, destination in temporary_paths:
            os.replace(temporary, destination)
    finally:
        for temporary, _ in temporary_paths:
            if temporary.exists():
                temporary.unlink()


def prepare(args: argparse.Namespace) -> dict:
    all_ts = bool(getattr(args, "all_ts", False))
    label = "all_ts" if all_ts else ts_label(args.ts)
    input_names = ALL_INPUT_NAMES if all_ts else INPUT_NAMES
    input_dim = len(input_names)
    if args.seed < 0:
        raise ValueError("The random seed must be nonnegative")
    source_root = args.source_root.expanduser().resolve()
    if all_ts:
        source_dirs = []
        for directory in source_root.glob("Ts=*"):
            if directory.is_dir():
                value = float(directory.name.removeprefix("Ts="))
                ts_label(value)
                source_dirs.append((value, directory))
        source_dirs.sort(key=lambda item: (item[0], item[1].name))
        if len({value for value, _ in source_dirs}) < 2:
            raise ValueError("--all-ts requires at least two distinct finite Ts values")
        source_dir = source_root
    else:
        source_dir = source_root / f"Ts={label}"
        if not source_dir.is_dir():
            raise ValueError(f"Ts directory does not exist: {source_dir}")
        source_dirs = [(float(args.ts), source_dir)]
    output_dir = (
        args.output_dir if args.output_dir is not None
        else Path(__file__).resolve().parent / "data" / (label if all_ts else f"ts_{label}")
    ).expanduser().resolve()
    try:
        output_dir.relative_to(source_root)
    except ValueError:
        pass
    else:
        raise ValueError("Output directory must be outside the raw source tree")
    if output_dir.exists() and not output_dir.is_dir():
        raise ValueError(f"Output path is not a directory: {output_dir}")
    existing = [name for name in OUTPUT_NAMES if (output_dir / name).exists()]
    if existing and not args.overwrite:
        raise ValueError(
            f"Output artifacts already exist in {output_dir}: {existing}. "
            "Use --overwrite explicitly to replace them."
        )
    # Deterministic ordering is part of the split definition. R comes from rows.
    paths = [(value, path) for value, directory in source_dirs
             for path in sorted(directory.glob("R=*.dat"), key=lambda path: path.name)]
    if not paths:
        raise ValueError(f"No R=*.dat files found in {source_dir}")
    all_rows, all_ids, all_lines, source_summaries = [], [], [], []
    for file_id, (value, path) in enumerate(paths):
        rows, lines, summary = read_source(path, file_id)
        if all_ts:
            rows = np.column_stack([np.full(len(rows), value), rows])
            summary["name"] = path.relative_to(source_root).as_posix()
        summary["Ts"] = value
        all_rows.append(rows)
        all_lines.append(lines)
        all_ids.append(np.full(len(rows), file_id, dtype=np.int32))
        source_summaries.append(summary)
    raw = np.concatenate(all_rows, axis=0)
    source_file_ids = np.concatenate(all_ids)
    source_lines = np.concatenate(all_lines)
    if len(raw) == 0:
        raise ValueError("No finite, nonsentinel samples found")
    raw, source_file_ids, source_lines, duplicate_count = deduplicate(
        raw, source_file_ids, source_lines, source_summaries, input_dim,
    )
    sample_count = len(raw)
    if sample_count < 10:
        raise ValueError("At least 10 unique valid samples are needed for nonempty 80/10/10 splits")
    x_raw, y_raw = raw[:, :input_dim].copy(), raw[:, input_dim:].copy()
    ts_values = sorted({value for value, _ in source_dirs})
    split_parts = {name: [] for name in ("train", "validation", "test")}
    split_by_ts, counts_by_ts = {}, {}
    for value in ts_values:
        members = np.flatnonzero(x_raw[:, 0] == value) if all_ts else np.arange(sample_count)
        if len(members) < 10:
            raise ValueError(f"Ts={value} needs at least 10 usable samples")
        # Resetting the same seed per Ts preserves the original single-Ts split.
        permutation = members[np.random.default_rng(args.seed).permutation(len(members))]
        train_size, val_size = math.floor(0.8 * len(members)), math.floor(0.1 * len(members))
        indices = {"train": permutation[:train_size],
                   "validation": permutation[train_size:train_size + val_size],
                   "test": permutation[train_size + val_size:]}
        key = ts_label(value)
        split_by_ts[key] = {name: len(part) for name, part in indices.items()}
        sources = [item for item in source_summaries if item["Ts"] == value]
        counts_by_ts[key] = {
            "source_files": len(sources), "usable_rows": len(members),
            **{name: sum(item[name] for item in sources) for name in
               ("raw_rows", "excluded_missing_sentinel", "excluded_paired_sentinel",
                "excluded_omega_i_minus_one_only", "excluded_nonfinite", "duplicate_rows_removed")},
            "omega_i_negative": int((y_raw[members, 1] < 0).sum()),
            "omega_i_zero": int((y_raw[members, 1] == 0).sum()),
            "omega_i_positive": int((y_raw[members, 1] > 0).sum()),
        }
        for name, part in indices.items():
            split_parts[name].append(part)
    train_idx, val_idx, test_idx = (np.concatenate(split_parts[name])
                                   for name in ("train", "validation", "test"))
    x_mean, x_std = training_scaler(x_raw[train_idx], input_names)
    y_mean, y_std = training_scaler(y_raw[train_idx], TARGET_NAMES)
    x = ((x_raw - x_mean) / x_std).astype(np.float32)
    y = ((y_raw - y_mean) / y_std).astype(np.float32)
    if not np.isfinite(x).all() or not np.isfinite(y).all():
        raise ValueError("Standardized samples cannot be represented as finite float32 values")
    counts = {
        "source_files": len(paths),
        "raw_rows": sum(item["raw_rows"] for item in source_summaries),
        "excluded_missing_sentinel": sum(item["excluded_missing_sentinel"] for item in source_summaries),
        "excluded_paired_sentinel": sum(item["excluded_paired_sentinel"] for item in source_summaries),
        "excluded_omega_i_minus_one_only": sum(item["excluded_omega_i_minus_one_only"] for item in source_summaries),
        "excluded_nonfinite": sum(item["excluded_nonfinite"] for item in source_summaries),
        "valid_rows_before_deduplication": sum(item["valid_rows_before_deduplication"] for item in source_summaries),
        "duplicate_rows_removed": duplicate_count,
        "usable_rows": sample_count,
        "omega_r_negative": int((y_raw[:, 0] < 0.0).sum()),
        "omega_r_zero": int((y_raw[:, 0] == 0.0).sum()),
        "omega_r_positive": int((y_raw[:, 0] > 0.0).sum()),
        "omega_i_negative": int((y_raw[:, 1] < 0.0).sum()),
        "omega_i_zero": int((y_raw[:, 1] == 0.0).sum()),
        "omega_i_positive": int((y_raw[:, 1] > 0.0).sum()),
    }
    scalers = {"x_mean": x_mean, "x_std": x_std, "y_mean": y_mean, "y_std": y_std}
    metadata = {
        "schema_version": 2,
        "Ts": None if all_ts else float(args.ts),
        "Ts_values": ts_values,
        "source_directory": str(source_dir),
        "input_names": input_names,
        "target_names": TARGET_NAMES,
        "growth_rate_component": "omega_i",
        "growth_convention": "omega_i > 0: growing; omega_i < 0: decaying; omega_i == 0: neutral",
        "counts": counts,
        "counts_by_ts": counts_by_ts,
        "split_by_ts": split_by_ts,
        "seed": args.seed,
        "split": {
            "train": len(train_idx), "validation": len(val_idx), "test": len(test_idx),
            "fractions": {"train": 0.8, "validation": 0.1, "test": 0.1},
            "method": "independent numpy.random.default_rng(seed).permutation(N_Ts) per Ts; floor(0.8N_Ts), floor(0.1N_Ts), remainder; concatenate splits in numeric Ts order",
        },
        "ranges": {"inputs": statistics(x_raw, input_names), "targets": statistics(y_raw, TARGET_NAMES)},
        "scaler": {
            "method": "standardization", "fit_split": "train", "std_ddof": 0,
            **{name: value.tolist() for name, value in scalers.items()},
            "transform": "normalized = (raw - mean) / std",
            "inverse_transform": "raw = normalized * std + mean",
        },
        "arrays": {
            "x_raw": f"float64, [N, {input_dim}]", "y_raw": "float64, [N, 2]",
            "x": f"float32, [N, {input_dim}]", "y": "float32, [N, 2]",
            "train_idx": "int64, indices into the N rows", "val_idx": "int64",
            "test_idx": "int64", "source_file_id": "int32, index into source_files",
            "source_line": "int64, original 1-based line number",
        },
        "filter": "exclude any omega_i=-1 as confirmed missing (including paired (-1,-1)), plus nonfinite rows; retain other negative values",
        "duplicates": "exactly equal inputs and targets keep first; conflicting targets abort",
        "file_order": "numeric Ts order, then lexicographic filename order; rows retain source line order",
        "source_files": source_summaries,
        "numpy_version": np.__version__,
        "load_npz": "numpy.load(path, allow_pickle=False)",
    }
    arrays = {
        "x_raw": x_raw, "y_raw": y_raw, "x": x, "y": y,
        "train_idx": train_idx, "val_idx": val_idx, "test_idx": test_idx,
        "source_file_id": source_file_ids, "source_line": source_lines,
    }
    save_outputs(output_dir, arrays, scalers, metadata)
    return {"output_directory": str(output_dir), "input_names": input_names,
            "Ts_values": ts_values, "counts": counts, "counts_by_ts": counts_by_ts,
            "split": metadata["split"], "split_by_ts": split_by_ts}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--data-root", "--source-root", dest="source_root", type=Path, default=DEFAULT_SOURCE_ROOT, help="Directory containing Ts=... subdirectories")
    selection = parser.add_mutually_exclusive_group()
    selection.add_argument("--ts", type=float, default=0.0, help="Fixed Ts slice (default: 0.0)")
    selection.add_argument("--all-ts", action="store_true", help="Combine all Ts directories using four input columns")
    parser.add_argument("--seed", type=int, default=42, help="Seed for the fixed 80/10/10 split (default: 42)")
    parser.add_argument("--output-dir", type=Path, help="Default: data/all_ts or data/ts_<Ts> alongside this script")
    parser.add_argument("--overwrite", action="store_true", help="Explicitly replace existing prepared artifacts; never delete the directory")
    args = parser.parse_args()
    try:
        summary = prepare(args)
    except (OSError, ValueError) as exc:
        print(f"Data preparation failed: {exc}", file=sys.stderr)
        return 1
    print(json.dumps(summary, ensure_ascii=False, indent=2, allow_nan=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
