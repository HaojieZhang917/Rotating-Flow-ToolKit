#!/usr/bin/env python3
"""Train a supervised MLP with fixed Ts or Ts as an input feature.

The data split and scalers are created by prepare_data.py. This script never
changes that split: only the order of training samples changes each epoch.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import inspect
import json
import math
import os
import random
import time
from pathlib import Path

import numpy as np

# Required by deterministic CUDA matrix multiplication; set before CUDA starts.
os.environ["CUBLAS_WORKSPACE_CONFIG"] = ":4096:8"
import torch
from torch import nn
from torch.utils.data import DataLoader, TensorDataset


ROOT = Path(__file__).resolve().parent
DEFAULT_DATA_DIR = ROOT / "data" / "all_ts"
DEFAULT_OUTPUT_DIR = ROOT / "runs" / "all_ts_clean" / "baseline"
FEATURE_NAMES = ["R", "alpha", "beta"]
ALL_FEATURE_NAMES = ["Ts", *FEATURE_NAMES]
TARGET_NAMES = ["omega_r", "omega_i"]


class MLP(nn.Module):
    """A small fully connected network with an unrestricted linear output."""

    def __init__(self, hidden_sizes: list[int], input_dim: int = 3):
        super().__init__()
        layers: list[nn.Module] = []
        if input_dim not in (3, 4):
            raise ValueError("Only three- or four-input eigenvalue models are supported")
        previous_size = input_dim
        for size in hidden_sizes:
            layers.extend([nn.Linear(previous_size, size), nn.Tanh()])
            previous_size = size
        layers.append(nn.Linear(previous_size, len(TARGET_NAMES)))
        self.network = nn.Sequential(*layers)

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        return self.network(x)


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def load_checkpoint(path: Path) -> dict:
    """Load this pipeline's local checkpoints on modern or legacy PyTorch.

    Old CUDA installations lack weights_only. Only locally generated, trusted
    checkpoints should be supplied when running with that legacy installation.
    """
    options = {"map_location": "cpu"}
    if "weights_only" in inspect.signature(torch.load).parameters:
        options["weights_only"] = True
    return torch.load(path, **options)


def load_prepared_data(data_dir: Path) -> tuple[dict, dict, dict]:
    """Load, validate and retain the original row indices of the fixed split."""
    with np.load(data_dir / "dataset.npz", allow_pickle=False) as archive:
        data = {name: archive[name].copy() for name in (
            "x_raw", "y_raw", "x", "y", "train_idx", "val_idx", "test_idx"
        )}
    with np.load(data_dir / "scaler.npz", allow_pickle=False) as archive:
        scalers = {name: np.asarray(archive[name], dtype=np.float64).copy()
                   for name in ("x_mean", "x_std", "y_mean", "y_std")}
    with (data_dir / "metadata.json").open(encoding="utf-8") as stream:
        metadata = json.load(stream)
    feature_names = metadata.get("input_names")
    if feature_names not in (FEATURE_NAMES, ALL_FEATURE_NAMES):
        raise ValueError("Metadata input_names must specify the three- or four-input column order")
    input_dim = len(feature_names)
    for key, expected_names in (("target_names", TARGET_NAMES),):
        if metadata.get(key) != expected_names:
            raise ValueError(f"Metadata {key} does not match the expected column order")

    n_samples = len(data["x_raw"])
    for name, width in (("x_raw", input_dim), ("x", input_dim), ("y_raw", 2), ("y", 2)):
        array = data[name]
        if array.shape != (n_samples, width):
            raise ValueError(f"{name} must have shape ({n_samples}, {width}); got {array.shape}")
        if not np.isfinite(array).all():
            raise ValueError(f"{name} contains non-finite values")
    if n_samples < 3:
        raise ValueError("At least three samples are required")
    for name, width in (("x_mean", input_dim), ("x_std", input_dim), ("y_mean", 2), ("y_std", 2)):
        if scalers[name].shape != (width,) or not np.isfinite(scalers[name]).all():
            raise ValueError(f"{name} must contain {width} finite values")
    if (scalers["x_std"] <= 0).any() or (scalers["y_std"] <= 0).any():
        raise ValueError("Standard deviations must be positive")

    splits = []
    for name in ("train_idx", "val_idx", "test_idx"):
        indices = data[name]
        if indices.ndim != 1 or indices.dtype.kind not in "iu" or len(indices) == 0:
            raise ValueError(f"{name} must be a nonempty one-dimensional integer array")
        if (indices < 0).any() or (indices >= n_samples).any():
            raise ValueError(f"{name} has out-of-bounds indices")
        data[name] = indices.astype(np.int64, copy=False)
        splits.append(data[name])
    combined = np.concatenate(splits)
    if len(combined) != n_samples or not np.array_equal(np.sort(combined), np.arange(n_samples)):
        raise ValueError("The splits must be disjoint and contain every sample exactly once")
    for raw_name, mean_name, std_name in (("x_raw", "x_mean", "x_std"), ("y_raw", "y_mean", "y_std")):
        train_values = data[raw_name][data["train_idx"]]
        if not np.allclose(scalers[mean_name], train_values.mean(axis=0, dtype=np.float64), rtol=1e-12, atol=1e-12):
            raise ValueError(f"{mean_name} was not fitted on this training split")
        if not np.allclose(scalers[std_name], train_values.std(axis=0, ddof=0, dtype=np.float64), rtol=1e-12, atol=1e-12):
            raise ValueError(f"{std_name} was not fitted on this training split with ddof=0")

    # Confirm that network inputs correspond to the stored physical values.
    for scaled_name, raw_name, mean_name, std_name in (
        ("x", "x_raw", "x_mean", "x_std"),
        ("y", "y_raw", "y_mean", "y_std"),
    ):
        expected = (data[raw_name] - scalers[mean_name]) / scalers[std_name]
        if not np.allclose(data[scaled_name], expected, rtol=1e-5, atol=5e-6):
            raise ValueError(f"{scaled_name} is inconsistent with the saved scalers")
        data[scaled_name] = np.asarray(data[scaled_name], dtype=np.float32)

    if input_dim == 3:
        ts = float(metadata.get("Ts", metadata.get("ts", 0.0)))
        if not math.isfinite(ts):
            raise ValueError("The fixed Ts must be finite")
        metadata["Ts"] = ts
        metadata["Ts_values"] = [ts]
    else:
        if metadata.get("Ts") is not None:
            raise ValueError("A four-input dataset must not declare a fixed Ts")
        actual_ts = np.unique(data["x_raw"][:, 0]).tolist()
        if len(actual_ts) < 2 or metadata.get("Ts_values") != actual_ts:
            raise ValueError("Ts_values does not match the variable Ts input column")
        metadata["Ts"] = None
    if metadata.get("growth_rate_component", "omega_i") != "omega_i":
        raise ValueError("The current physical convention requires omega_i as growth rate")
    return data, scalers, metadata


def make_loader(data: dict, indices: np.ndarray, batch_size: int, *,
                shuffle: bool = False, generator: torch.Generator | None = None) -> DataLoader:
    dataset = TensorDataset(torch.from_numpy(data["x"][indices]),
                            torch.from_numpy(data["y"][indices]))
    return DataLoader(dataset, batch_size=batch_size, shuffle=shuffle,
                      generator=generator, num_workers=0, drop_last=False)


class ResidentBatches:
    """Keep a split on the device and gather each batch without Python row reads."""

    def __init__(self, data: dict, indices: np.ndarray, batch_size: int,
                 device: torch.device, *, shuffle: bool = False,
                 generator: torch.Generator | None = None):
        self.x = torch.from_numpy(data["x"][indices]).to(device)
        self.y = torch.from_numpy(data["y"][indices]).to(device)
        self.batch_size = batch_size
        self.shuffle = shuffle
        self.generator = generator

    def __iter__(self):
        if self.shuffle:
            order = torch.randperm(len(self.x), device=self.x.device, generator=self.generator)
            for start in range(0, len(self.x), self.batch_size):
                indices = order[start:start + self.batch_size]
                yield self.x[indices], self.y[indices]
        else:
            for start in range(0, len(self.x), self.batch_size):
                yield self.x[start:start + self.batch_size], self.y[start:start + self.batch_size]

    def __len__(self) -> int:
        return math.ceil(len(self.x) / self.batch_size)


def batch_source(data: dict, indices: np.ndarray, batch_size: int, device: torch.device,
                 resident: bool, *, shuffle: bool = False,
                 generator: torch.Generator | None = None) -> DataLoader | ResidentBatches:
    if resident:
        return ResidentBatches(data, indices, batch_size, device, shuffle=shuffle, generator=generator)
    return make_loader(data, indices, batch_size, shuffle=shuffle, generator=generator)


def finite_batch(x: torch.Tensor, y: torch.Tensor) -> None:
    if not (torch.isfinite(x).all() & torch.isfinite(y).all()).item():
        raise FloatingPointError("A batch contains a non-finite input or label")


def train_epoch(model: MLP, loader: DataLoader | ResidentBatches, optimizer: torch.optim.Optimizer,
                device: torch.device) -> float:
    model.train()
    total_loss = 0.0
    sample_count = 0
    for x_batch, y_batch in loader:
        x_batch, y_batch = x_batch.to(device), y_batch.to(device)
        finite_batch(x_batch, y_batch)
        optimizer.zero_grad(set_to_none=True)
        prediction = model(x_batch)
        loss = nn.functional.mse_loss(prediction, y_batch, reduction="mean")
        if not torch.isfinite(loss).item():
            raise FloatingPointError("Training loss became non-finite")
        loss.backward()
        gradient_checks = [torch.isfinite(parameter.grad).all() for parameter in model.parameters()
                           if parameter.grad is not None]
        if gradient_checks and not torch.stack(gradient_checks).all().item():
            raise FloatingPointError("A network gradient became non-finite")
        optimizer.step()
        batch_size = len(x_batch)
        total_loss += loss.item() * batch_size
        sample_count += batch_size
    return total_loss / sample_count


@torch.no_grad()
def validation_statistics(model: MLP, loader: DataLoader | ResidentBatches,
                          device: torch.device, y_std: np.ndarray,
                          y_mean: np.ndarray | None = None,
                          y_raw: torch.Tensor | None = None) -> dict:
    model.eval()
    squared_error = torch.zeros(2, dtype=torch.float64, device=device)
    output_std = torch.as_tensor(y_std, dtype=torch.float64, device=device)
    output_mean = None if y_mean is None else torch.as_tensor(y_mean, dtype=torch.float64, device=device)
    sample_count = 0
    for x_batch, y_batch in loader:
        x_batch, y_batch = x_batch.to(device), y_batch.to(device)
        finite_batch(x_batch, y_batch)
        prediction = model(x_batch)
        if y_raw is None:
            error = (prediction - y_batch).to(torch.float64)
        else:
            # Evaluate against the original float64 labels rather than rounded
            # float32 standardized labels, matching final physical metrics.
            raw_batch = y_raw[sample_count:sample_count + len(x_batch)]
            error = (prediction.to(torch.float64) * output_std + output_mean - raw_batch) / output_std
        if not torch.isfinite(error).all().item():
            raise FloatingPointError("Validation loss became non-finite")
        squared_error += error.square().sum(dim=0)
        sample_count += len(x_batch)
    component_mse = squared_error.cpu().numpy() / sample_count
    return {"standardized_mse": float(component_mse.mean()),
            "physical_rmse": np.sqrt(component_mse) * y_std}


def validation_loss(model: MLP, loader: DataLoader | ResidentBatches, device: torch.device) -> float:
    return validation_statistics(model, loader, device, np.ones(2))["standardized_mse"]


@torch.no_grad()
def predict_loader(model: MLP, loader: DataLoader | ResidentBatches, device: torch.device,
                   scalers: dict) -> np.ndarray:
    model.eval()
    predictions = []
    for x_batch, y_batch in loader:
        finite_batch(x_batch, y_batch)
        prediction = model(x_batch.to(device)).cpu().numpy().astype(np.float64)
        if not np.isfinite(prediction).all():
            raise FloatingPointError("Predictions contain non-finite values")
        predictions.append(prediction)
    return np.concatenate(predictions) * scalers["y_std"] + scalers["y_mean"]


def cpu_snapshot(value):
    """Store optimizer and RNG state using tensors and simple containers only."""
    if isinstance(value, torch.Tensor):
        return value.detach().cpu().clone()
    if isinstance(value, dict):
        return {key: cpu_snapshot(item) for key, item in value.items()}
    if isinstance(value, list):
        return [cpu_snapshot(item) for item in value]
    if isinstance(value, tuple):
        return tuple(cpu_snapshot(item) for item in value)
    return value


def selection_score(statistics: dict, targets: np.ndarray | None) -> float:
    if targets is None:
        return statistics["standardized_mse"]
    return float(np.max(statistics["physical_rmse"] / targets))


def physical_metrics(y_true: np.ndarray, y_pred: np.ndarray,
                     y_std: np.ndarray, sign_tolerance: float) -> dict:
    """Report each component in its original units, including growth-rate sign."""
    error = y_pred - y_true
    if not np.isfinite(error).all():
        raise FloatingPointError("Physical prediction errors are non-finite")
    components = {}
    for column, name in enumerate(TARGET_NAMES):
        residual_sum = float(np.sum(error[:, column] ** 2))
        centered_sum = float(np.sum((y_true[:, column] - y_true[:, column].mean()) ** 2))
        components[name] = {
            "mae": float(np.mean(np.abs(error[:, column]))),
            "rmse": float(np.sqrt(np.mean(error[:, column] ** 2))),
            "max_abs_error": float(np.max(np.abs(error[:, column]))),
            "r2": 1.0 - residual_sum / centered_sum if centered_sum > 0 else None,
        }
    true_growth, pred_growth = y_true[:, 1], y_pred[:, 1]
    sign_errors = np.sign(true_growth) != np.sign(pred_growth)
    retained = np.abs(true_growth) > sign_tolerance
    return {
        "n_samples": len(y_true),
        "standardized_mse": float(np.mean((error / y_std) ** 2)),
        "components": components,
        "growth_rate_sign": {
            "component": "omega_i",
            "definition": "np.sign(predicted omega_i) != np.sign(true omega_i); zero is neutral",
            "error_count": int(sign_errors.sum()),
            "error_rate": float(sign_errors.mean()),
            "true_unstable_count": int((true_growth > 0).sum()),
            "true_stable_count": int((true_growth < 0).sum()),
            "true_neutral_count": int((true_growth == 0).sum()),
            "near_zero_exclusion_threshold": sign_tolerance,
            "n_samples_excluding_near_zero": int(retained.sum()),
            "error_rate_excluding_near_zero": float(sign_errors[retained].mean()) if retained.any() else None,
        },
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-dir", type=Path, default=DEFAULT_DATA_DIR)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--overwrite", action="store_true", help="Allow writing to an existing output directory; no directory is deleted")
    parser.add_argument("--hidden-sizes", type=int, nargs="+", help="Default: 64 64, or the resumed model architecture")
    parser.add_argument("--resume-model", type=Path, help="Continue from a compatible checkpoint; --epochs counts additional epochs")
    parser.add_argument("--reset-optimizer", action="store_true", help="Reinitialize Adam instead of restoring saved moments")
    parser.add_argument("--resident-data", action="store_true", help="Keep each used split on the device and batch with tensor indexing")
    parser.add_argument("--validation-only", action="store_true", help="Skip all test predictions and metrics during model selection")
    parser.add_argument("--evaluate-only", action="store_true", help="Evaluate --resume-model without updating any network parameter")
    parser.add_argument("--epochs", "--max-epochs", type=int, default=80)
    parser.add_argument("--patience", type=int, default=12, help="Early-stopping patience; 0 disables early stopping")
    parser.add_argument("--min-delta", type=float, default=1e-5,
                        help="Required selection-score improvement for resetting early-stopping patience")
    parser.add_argument("--batch-size", type=int, default=4096)
    parser.add_argument("--eval-batch-size", type=int, default=16384)
    parser.add_argument("--lr", type=float, help="Default: 1e-3, or the resumed optimizer learning rate")
    parser.add_argument("--lr-patience", type=int, default=0, help="ReduceLROnPlateau patience; 0 disables the scheduler")
    parser.add_argument("--lr-factor", type=float, default=0.5)
    parser.add_argument("--min-lr", type=float, default=1e-7)
    parser.add_argument("--target-rmse-r", type=float, help="Physical validation RMSE target for omega_r; specify both targets")
    parser.add_argument("--target-rmse-i", type=float, help="Physical validation RMSE target for omega_i; specify both targets")
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--threads", type=int, default=4)
    parser.add_argument("--device", default="cpu", help="cpu or an available CUDA device")
    parser.add_argument("--sign-tolerance", type=float, default=0.0,
                        help="Exclude |true omega_i| <= this value only from the additional sign-error metric")
    args = parser.parse_args()
    for name in ("epochs", "batch_size", "eval_batch_size", "threads"):
        if getattr(args, name) <= 0:
            parser.error(f"--{name.replace('_', '-')} must be positive")
    if args.patience < 0 or args.lr_patience < 0:
        parser.error("--patience and --lr-patience must be nonnegative")
    if args.hidden_sizes is not None and any(size <= 0 for size in args.hidden_sizes):
        parser.error("Hidden layer sizes must be positive")
    if not 0 <= args.seed < 2 ** 32:
        parser.error("--seed must be between 0 and 2**32 - 1")
    for name in ("lr", "min_lr", "min_delta", "sign_tolerance"):
        value = getattr(args, name)
        if value is not None and (not math.isfinite(value) or value < 0 or (name in ("lr", "min_lr") and value == 0)):
            parser.error(f"Invalid --{name.replace('_', '-')}")
    if not math.isfinite(args.lr_factor) or not 0 < args.lr_factor < 1:
        parser.error("--lr-factor must be between zero and one")
    if (args.target_rmse_r is None) != (args.target_rmse_i is None):
        parser.error("Specify both --target-rmse-r and --target-rmse-i")
    for name in ("target_rmse_r", "target_rmse_i"):
        value = getattr(args, name)
        if value is not None and (not math.isfinite(value) or value <= 0):
            parser.error(f"--{name.replace('_', '-')} must be finite and positive")
    if args.evaluate_only and args.resume_model is None:
        parser.error("--evaluate-only requires --resume-model")
    args.data_dir = args.data_dir.resolve()
    args.output_dir = args.output_dir.resolve()
    if args.resume_model is not None:
        args.resume_model = args.resume_model.resolve()
    if args.output_dir.exists() and not args.overwrite:
        parser.error(f"Output directory already exists: {args.output_dir}; choose another directory or use --overwrite")
    return args


def main() -> None:
    args = parse_args()
    start_time = time.perf_counter()
    random.seed(args.seed)
    np.random.seed(args.seed)
    torch.manual_seed(args.seed)
    torch.set_num_threads(args.threads)
    torch.use_deterministic_algorithms(True)
    torch.backends.cuda.matmul.allow_tf32 = False
    torch.backends.cudnn.allow_tf32 = False
    torch.backends.cudnn.benchmark = False
    torch.backends.cudnn.deterministic = True
    if hasattr(torch, "set_float32_matmul_precision"):
        torch.set_float32_matmul_precision("highest")
    device = torch.device(args.device)
    if device.type not in ("cpu", "cuda"):
        raise ValueError("Only cpu and CUDA devices are supported")
    if device.type == "cuda":
        if not torch.cuda.is_available():
            raise ValueError("CUDA was requested but is unavailable")
        if device.index is None:
            device = torch.device("cuda", torch.cuda.current_device())
        if device.index >= torch.cuda.device_count():
            raise ValueError(f"CUDA device index {device.index} is unavailable")
        torch.cuda.manual_seed_all(args.seed)

    data, scalers, metadata = load_prepared_data(args.data_dir)
    feature_names = metadata["input_names"]
    input_dim = len(feature_names)
    protected_directories = [args.data_dir]
    if metadata.get("source_directory"):
        source_directory = Path(metadata["source_directory"])
        # A Linux source path recorded in WSL metadata is not a native absolute
        # Windows path. Training only reads the prepared artifacts.
        if source_directory.is_absolute():
            protected_directories.append(source_directory.resolve())
    for protected_directory in protected_directories:
        if args.output_dir.is_relative_to(protected_directory):
            raise ValueError(f"Training output must not be inside a data directory: {protected_directory}")
    hashes = {name: sha256_file(args.data_dir / name)
              for name in ("dataset.npz", "scaler.npz", "metadata.json")}
    train_idx, val_idx, test_idx = (data[name] for name in ("train_idx", "val_idx", "test_idx"))
    targets = None if args.target_rmse_r is None else np.array([args.target_rmse_r, args.target_rmse_i])
    resume_checkpoint = None
    start_epoch = 0
    resume_details = {"checkpoint": None, "optimizer": "new", "scheduler": "disabled", "batch_rng": "new"}
    if args.resume_model is not None:
        resume_checkpoint = load_checkpoint(args.resume_model)
        if resume_checkpoint.get("data_hashes") != hashes:
            raise ValueError("Resume checkpoint data hashes differ from this fixed prepared dataset")
        if resume_checkpoint.get("Ts") != metadata["Ts"]:
            raise ValueError("Resume checkpoint has a different fixed Ts")
        if resume_checkpoint.get("feature_names") != feature_names or resume_checkpoint.get("target_names") != TARGET_NAMES:
            raise ValueError("Resume checkpoint has an incompatible column order")
        architecture = resume_checkpoint["architecture"]
        if architecture.get("input_dim") != input_dim or architecture.get("output_dim") != 2 or architecture.get("activation") != "tanh":
            raise ValueError(f"Resume requires a {input_dim}-input, 2-output Tanh MLP")
        resumed_sizes = architecture["hidden_sizes"]
        if args.hidden_sizes is not None and args.hidden_sizes != resumed_sizes:
            raise ValueError("--hidden-sizes cannot change the architecture of --resume-model")
        args.hidden_sizes = resumed_sizes
        for name, current_value in scalers.items():
            resumed_value = resume_checkpoint["scalers"][name].cpu().numpy()
            if not np.array_equal(resumed_value, current_value):
                raise ValueError(f"Resume checkpoint scaler {name} differs from the prepared data")
        start_epoch = int(resume_checkpoint.get("epoch", resume_checkpoint.get("best_epoch", 0)))
        resume_details["checkpoint"] = str(args.resume_model)
    elif args.hidden_sizes is None:
        args.hidden_sizes = [64, 64]

    generator_device = device if args.resident_data else torch.device("cpu")
    generator = torch.Generator(device=generator_device).manual_seed(args.seed)
    train_loader = None if args.evaluate_only else batch_source(
        data, train_idx, args.batch_size, device, args.resident_data, shuffle=True, generator=generator)
    val_loader = batch_source(data, val_idx, args.eval_batch_size, device, args.resident_data)
    val_raw = torch.from_numpy(data["y_raw"][val_idx]).to(device)
    model = MLP(args.hidden_sizes, input_dim=input_dim).to(device)
    if resume_checkpoint is not None:
        model.load_state_dict(resume_checkpoint["state_dict"])
    requested_lr = args.lr
    optimizer = torch.optim.Adam(model.parameters(), lr=requested_lr or 1e-3)
    if resume_checkpoint is not None:
        if not args.reset_optimizer and "optimizer_state_dict" in resume_checkpoint:
            optimizer.load_state_dict(resume_checkpoint["optimizer_state_dict"])
            resume_details["optimizer"] = "restored_saved_Adam_state"
            if requested_lr is not None:
                for group in optimizer.param_groups:
                    group["lr"] = requested_lr
        else:
            resume_details["optimizer"] = "reinitialized_by_request" if args.reset_optimizer else "reinitialized_legacy_checkpoint_has_no_Adam_state"
        if resume_checkpoint.get("seed") == args.seed and resume_checkpoint.get("batch_rng_device") == str(generator.device):
            if "batch_rng_state" in resume_checkpoint:
                generator.set_state(resume_checkpoint["batch_rng_state"])
                resume_details["batch_rng"] = "restored"
        else:
            resume_details["batch_rng"] = "reseeded_due_to_missing_state_seed_or_device_change"
        if resume_checkpoint.get("seed") == args.seed and "torch_rng_state" in resume_checkpoint:
            torch.set_rng_state(resume_checkpoint["torch_rng_state"])
            if device.type == "cuda" and len(resume_checkpoint.get("cuda_rng_states", [])) == torch.cuda.device_count():
                torch.cuda.set_rng_state_all(resume_checkpoint["cuda_rng_states"])
    args.lr = float(optimizer.param_groups[0]["lr"])
    scheduler = None
    if args.lr_patience:
        scheduler = torch.optim.lr_scheduler.ReduceLROnPlateau(
            optimizer, mode="min", factor=args.lr_factor, patience=args.lr_patience,
            threshold=args.min_delta, threshold_mode="abs", min_lr=args.min_lr)
        resume_details["scheduler"] = "new"
        if resume_checkpoint is not None and not args.reset_optimizer:
            old_config = resume_checkpoint.get("training_config", {})
            same_scheduler = all(old_config.get(name) == getattr(args, name)
                                 for name in ("lr_patience", "lr_factor", "min_lr", "min_delta"))
            same_targets = resume_checkpoint.get("selection_targets") == (None if targets is None else targets.tolist())
            if same_scheduler and same_targets and "scheduler_state_dict" in resume_checkpoint:
                scheduler.load_state_dict(resume_checkpoint["scheduler_state_dict"])
                resume_details["scheduler"] = "restored"

    args.output_dir.mkdir(parents=True, exist_ok=True)
    checkpoint_path = args.output_dir / "best_model.pt"
    configuration = {name: str(value) if isinstance(value, Path) else value
                     for name, value in vars(args).items()}
    configuration["loss"] = "MSE mean over samples and both standardized outputs"
    configuration["optimizer"] = "Adam"
    configuration["activation"] = "Tanh"
    configuration["n_parameters"] = sum(parameter.numel() for parameter in model.parameters())
    configuration["torch_version"] = str(torch.__version__)
    configuration["numpy_version"] = np.__version__
    configuration["resolved_device"] = str(device)
    configuration["precision"] = "float32 network and Adam; float64 validation accumulation; TF32 disabled"
    configuration["selection_metric"] = "standardized_MSE" if targets is None else "max(physical_RMSE_r/target_r, physical_RMSE_i/target_i)"
    history = []
    best_loss = float("inf")
    best_score = float("inf")
    best_epoch = start_epoch
    best_rmse = None
    patience_reference = float("inf")
    stale_epochs = 0
    stop_reason = "epoch_limit"
    training_min = torch.from_numpy(data["x_raw"][train_idx].min(axis=0))
    training_max = torch.from_numpy(data["x_raw"][train_idx].max(axis=0))

    def save_checkpoint(epoch: int, run_epoch: int) -> None:
        checkpoint = {
            "format_version": 2,
            "state_dict": cpu_snapshot(model.state_dict()),
            "optimizer_state_dict": cpu_snapshot(optimizer.state_dict()),
            "architecture": {"input_dim": input_dim, "hidden_sizes": args.hidden_sizes, "output_dim": 2, "activation": "tanh"},
            "Ts": metadata["Ts"], "Ts_values": metadata["Ts_values"],
            "feature_names": feature_names, "target_names": TARGET_NAMES,
            "growth_rate_component": "omega_i",
            "scalers": {name: torch.from_numpy(value.copy()) for name, value in scalers.items()},
            "training_input_min": training_min, "training_input_max": training_max,
            "seed": args.seed, "epoch": epoch, "run_epoch": run_epoch,
            "start_epoch": start_epoch, "best_epoch": best_epoch,
            "best_validation_loss": best_loss, "best_validation_score": best_score,
            "validation_rmse": best_rmse.tolist(),
            "selection_targets": None if targets is None else targets.tolist(),
            "data_dir": str(args.data_dir), "data_hashes": hashes,
            "training_config": configuration, "resume_details": resume_details,
            "batch_rng_state": generator.get_state(), "batch_rng_device": str(generator.device),
            "torch_rng_state": torch.get_rng_state(),
            "cuda_rng_states": torch.cuda.get_rng_state_all() if device.type == "cuda" else [],
            "early_stopping": {"reference": patience_reference, "stale_epochs": stale_epochs},
        }
        if scheduler is not None:
            checkpoint["scheduler_state_dict"] = scheduler.state_dict()
        torch.save(checkpoint, checkpoint_path)

    ts_description = str(metadata["Ts_values"]) if metadata["Ts"] is None else f"{metadata['Ts']:g} (fixed)"
    print(f"Ts={ts_description}; train={len(train_idx):,}, validation={len(val_idx):,}, test={len(test_idx):,}; device={device}", flush=True)
    print(f"MLP: {input_dim} -> {' -> '.join(map(str, args.hidden_sizes))} -> 2; parameters={configuration['n_parameters']:,}", flush=True)
    print(f"Batch source={'resident tensors' if args.resident_data else 'DataLoader'}; selection={configuration['selection_metric']}; test={'skipped' if args.validation_only else 'final evaluation only'}", flush=True)
    if resume_checkpoint is not None:
        statistics = validation_statistics(model, val_loader, device, scalers["y_std"], scalers["y_mean"], val_raw)
        best_loss = statistics["standardized_mse"]
        best_rmse = statistics["physical_rmse"]
        best_score = selection_score(statistics, targets)
        patience_reference = best_score
        save_checkpoint(start_epoch, 0)
        print(f"Resume epoch={start_epoch}; optimizer={resume_details['optimizer']}; validation RMSE_r={best_rmse[0]:.8g} RMSE_i={best_rmse[1]:.8g} score={best_score:.8g}", flush=True)
    if args.evaluate_only:
        stop_reason = "evaluate_only"
    elif targets is not None and best_score <= 1.0:
        stop_reason = "target_rmse_reached"
        print("The resumed model already meets both validation RMSE targets.", flush=True)
    with (args.output_dir / "history.csv").open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=["epoch", "run_epoch", "train_loss", "validation_loss", "validation_rmse_r", "validation_rmse_i", "validation_score", "best_validation_loss", "best_validation_score", "learning_rate", "elapsed_seconds"])
        writer.writeheader()
        epochs_to_run = args.epochs if stop_reason == "epoch_limit" else 0
        for run_epoch in range(1, epochs_to_run + 1):
            epoch = start_epoch + run_epoch
            train_loss = train_epoch(model, train_loader, optimizer, device)
            statistics = validation_statistics(model, val_loader, device, scalers["y_std"], scalers["y_mean"], val_raw)
            val_loss, rmse = statistics["standardized_mse"], statistics["physical_rmse"]
            score = selection_score(statistics, targets)
            new_best = score < best_score
            if score < patience_reference - args.min_delta:
                patience_reference, stale_epochs = score, 0
            else:
                stale_epochs += 1
            if scheduler is not None:
                scheduler.step(score)
            if new_best:
                best_loss, best_score, best_epoch, best_rmse = val_loss, score, epoch, rmse.copy()
                save_checkpoint(epoch, run_epoch)
            # Save every exact score improvement; min_delta controls patience only.
            row = {"epoch": epoch, "run_epoch": run_epoch, "train_loss": train_loss,
                   "validation_loss": val_loss, "validation_rmse_r": float(rmse[0]),
                   "validation_rmse_i": float(rmse[1]), "validation_score": score,
                   "best_validation_loss": best_loss, "best_validation_score": best_score,
                   "learning_rate": float(optimizer.param_groups[0]["lr"]),
                   "elapsed_seconds": time.perf_counter() - start_time}
            writer.writerow(row)
            stream.flush()
            history.append(row)
            print(f"epoch {epoch:4d} (+{run_epoch}/{args.epochs}): train={train_loss:.8g} val={val_loss:.8g} RMSE_r={rmse[0]:.8g} RMSE_i={rmse[1]:.8g} best-score={best_score:.8g} lr={row['learning_rate']:.3g} elapsed={row['elapsed_seconds']:.1f}s{' *' if new_best else ''}", flush=True)
            if targets is not None and score <= 1.0:
                stop_reason = "target_rmse_reached"
                print("Both physical validation RMSE targets have been reached.", flush=True)
                break
            if args.patience and stale_epochs >= args.patience:
                stop_reason = "early_stopping"
                print(f"Early stopping after {args.patience} epochs without improvement larger than {args.min_delta:g}.", flush=True)
                break

    # The test split is evaluated only after model selection is complete.
    checkpoint = load_checkpoint(checkpoint_path)
    model.load_state_dict(checkpoint["state_dict"])
    val_prediction = predict_loader(model, val_loader, device, scalers)
    evaluations_to_run = [("validation", val_idx, val_prediction)]
    test_prediction = None
    if not args.validation_only:
        test_loader = batch_source(data, test_idx, args.eval_batch_size, device, args.resident_data)
        test_prediction = predict_loader(model, test_loader, device, scalers)
        evaluations_to_run.append(("test", test_idx, test_prediction))

    # Fit two cheap reference models using training samples only.
    mean_prediction = data["y_raw"][train_idx].mean(axis=0)
    train_x = (data["x_raw"][train_idx] - scalers["x_mean"]) / scalers["x_std"]
    train_y = (data["y_raw"][train_idx] - scalers["y_mean"]) / scalers["y_std"]
    design = np.column_stack([np.ones(len(train_idx)), train_x])
    linear_coefficients, _, linear_rank, _ = np.linalg.lstsq(design, train_y, rcond=None)
    evaluation = {}
    by_ts = {str(value): {} for value in metadata["Ts_values"]}
    for split_name, indices, prediction in evaluations_to_run:
        true_y = data["y_raw"][indices]
        scaled_x = (data["x_raw"][indices] - scalers["x_mean"]) / scalers["x_std"]
        linear_prediction = np.column_stack([np.ones(len(indices)), scaled_x]) @ linear_coefficients
        linear_prediction = linear_prediction * scalers["y_std"] + scalers["y_mean"]
        evaluation[split_name] = {
            "mlp": physical_metrics(true_y, prediction, scalers["y_std"], args.sign_tolerance),
            "mean_baseline": physical_metrics(true_y, np.broadcast_to(mean_prediction, true_y.shape), scalers["y_std"], args.sign_tolerance),
            "linear_baseline": physical_metrics(true_y, linear_prediction, scalers["y_std"], args.sign_tolerance),
        }
        for value in metadata["Ts_values"]:
            mask = data["x_raw"][indices, 0] == value if input_dim == 4 else np.ones(len(indices), dtype=bool)
            if mask.any():
                by_ts[str(value)][split_name] = physical_metrics(
                    true_y[mask], prediction[mask], scalers["y_std"], args.sign_tolerance)
    metrics = {
        "Ts": metadata["Ts"], "Ts_values": metadata["Ts_values"],
        "feature_names": feature_names, "target_names": TARGET_NAMES,
        "growth_rate_component": "omega_i", "by_ts": by_ts,
        "configuration": configuration, "data_hashes": hashes,
        "split_sizes": {"train": len(train_idx), "validation": len(val_idx), "test": len(test_idx)},
        "best_epoch": best_epoch, "best_validation_loss": best_loss,
        "best_validation_score": best_score,
        "selection_targets": None if targets is None else {"omega_r": float(targets[0]), "omega_i": float(targets[1])},
        "target_met_on_validation": None if targets is None else all(
            evaluation["validation"]["mlp"]["components"][name]["rmse"] <= targets[column]
            for column, name in enumerate(TARGET_NAMES)),
        "epochs_completed": len(history), "start_epoch": start_epoch,
        "last_epoch": start_epoch + len(history), "early_stopped": stop_reason == "early_stopping",
        "stop_reason": stop_reason, "test_evaluated": not args.validation_only,
        "resume_details": resume_details,
        "runtime_seconds": time.perf_counter() - start_time,
        "baseline_details": {
            "mean": "Training-set mean of each physical output",
            "linear": "Ordinary least squares, fitted on standardized training inputs and outputs with an intercept",
            "linear_design_rank": int(linear_rank),
        },
        **evaluation,
    }
    with (args.output_dir / "metrics.json").open("w", encoding="utf-8") as stream:
        json.dump(metrics, stream, indent=2, ensure_ascii=False, allow_nan=False)
        stream.write("\n")
    if test_prediction is not None:
        np.savez_compressed(args.output_dir / "test_predictions.npz",
                            x_raw=data["x_raw"][test_idx], y_true=data["y_raw"][test_idx],
                            y_pred=test_prediction, error=test_prediction - data["y_raw"][test_idx],
                            dataset_indices=test_idx)
        print(f"Best epoch={best_epoch}; final test standardized MSE={evaluation['test']['mlp']['standardized_mse']:.8g}", flush=True)
    else:
        print(f"Best epoch={best_epoch}; validation-only run, test was not evaluated.", flush=True)
    print(f"Saved model and evaluation: {args.output_dir}", flush=True)


if __name__ == "__main__":
    main()
