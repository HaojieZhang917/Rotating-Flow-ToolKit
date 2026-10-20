"""Scientific data invariants and a small end-to-end pipeline check.

Run with: .venv-wsl/bin/python -m unittest -v test_pipeline
All fixtures and generated outputs live in temporary directories.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

import numpy as np


WORKSPACE = Path(__file__).resolve().parent


def run_script(name: str, *arguments: object) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        [sys.executable, str(WORKSPACE / name), *map(str, arguments)],
        text=True,
        capture_output=True,
        check=True,
    )


def write_fixture(source: Path, *, conflict: bool = False, ts: float = 0.0) -> Path:
    folder = source / f"Ts={ts}"
    folder.mkdir(parents=True)
    rows = [
        [100.0 + 5 * i, 0.05 + 0.005 * i, -0.25 + 0.002 * i,
         0.015 * i - 0.2 + ts, 0.002 * i - 0.04 + 0.001 * np.sin(i) - ts]
        for i in range(40)
    ]
    # Any omega_i=-1 is missing; omega_r=-1 with another omega_i remains valid.
    rows[3][3:] = [-1.0, 0.125]
    rows.append(rows[5].copy())
    rows.append([600.0, 0.2, -0.1, -1.0, -1.0])
    rows.append([605.0, 0.205, -0.1, 0.3, float("nan")])
    rows.append([610.0, 0.21, -0.1, 0.3, -1.0])
    if conflict:
        bad = rows[5].copy()
        bad[3] += 1.0
        rows.append(bad)
    # Deliberately different from the row values: R must come from the data.
    path = folder / "R=999.dat"
    header = (
        'TITLE = "Stability Analysis"\n'
        'VARIABLES = "R", "alpha", "beta", "omega_r", "omega_i"\n'
        f'ZONE T="R=999", I={len(rows)}, J=1, F=POINT\n'
    )
    body = "\n".join("\t".join(format(value, ".17g") for value in row) for row in rows)
    path.write_text(header + body + "\n", encoding="utf-8")
    return path


class DataPreparationTests(unittest.TestCase):
    def test_filtering_split_scaling_and_reproducibility(self) -> None:
        with tempfile.TemporaryDirectory(prefix="single-ts-ml-test-") as directory:
            root = Path(directory)
            source = root / "source"
            raw_file = write_fixture(source)
            before = hashlib.sha256(raw_file.read_bytes()).hexdigest()
            outputs = [root / "prepared_a", root / "prepared_b"]
            for output in outputs:
                run_script("prepare_data.py", "--data-root", source, "--ts", "0.0",
                           "--output-dir", output, "--seed", "42")
            self.assertEqual(before, hashlib.sha256(raw_file.read_bytes()).hexdigest())

            with np.load(outputs[0] / "dataset.npz", allow_pickle=False) as data, \
                 np.load(outputs[1] / "dataset.npz", allow_pickle=False) as repeat, \
                 np.load(outputs[0] / "scaler.npz", allow_pickle=False) as scaler:
                x_raw, y_raw = data["x_raw"], data["y_raw"]
                self.assertEqual(x_raw.shape, (40, 3))
                self.assertEqual(y_raw.shape, (40, 2))
                self.assertEqual(data["x"].dtype, np.float32)
                self.assertEqual(data["y"].dtype, np.float32)
                self.assertTrue(np.isfinite(data["x"]).all())
                self.assertTrue(np.isfinite(data["y"]).all())
                self.assertFalse(np.any(x_raw[:, 0] == 999))
                self.assertTrue(np.any(y_raw[:, 0] < 0))
                self.assertTrue(np.any((y_raw[:, 0] == -1) & (y_raw[:, 1] == 0.125)))
                self.assertFalse(np.any(np.all(y_raw == -1, axis=1)))
                self.assertFalse(np.any(y_raw[:, 1] == -1))

                indices = [data[key] for key in ("train_idx", "val_idx", "test_idx")]
                self.assertEqual([len(index) for index in indices], [32, 4, 4])
                np.testing.assert_array_equal(np.sort(np.concatenate(indices)), np.arange(40))
                for key in ("train_idx", "val_idx", "test_idx", "x_raw", "y_raw", "x", "y"):
                    np.testing.assert_array_equal(data[key], repeat[key])

                train = data["train_idx"]
                for raw, prefix, scaled in ((x_raw, "x", data["x"]), (y_raw, "y", data["y"])):
                    np.testing.assert_allclose(scaler[f"{prefix}_mean"], raw[train].mean(axis=0))
                    np.testing.assert_allclose(scaler[f"{prefix}_std"], raw[train].std(axis=0))
                    np.testing.assert_allclose(
                        scaled, (raw - scaler[f"{prefix}_mean"]) / scaler[f"{prefix}_std"],
                        rtol=2e-7, atol=2e-7,
                    )
                    np.testing.assert_allclose(
                        scaled.astype(np.float64) * scaler[f"{prefix}_std"] + scaler[f"{prefix}_mean"],
                        raw, rtol=2e-7, atol=2e-7,
                    )

            with self.assertRaises(subprocess.CalledProcessError) as existing:
                run_script("prepare_data.py", "--data-root", source, "--ts", "0.0",
                           "--output-dir", outputs[0], "--seed", "42")
            self.assertIn("already exist", existing.exception.stderr)

    def test_conflicting_duplicate_labels_are_rejected(self) -> None:
        with tempfile.TemporaryDirectory(prefix="single-ts-ml-conflict-") as directory:
            root = Path(directory)
            write_fixture(root / "source", conflict=True)
            with self.assertRaises(subprocess.CalledProcessError) as conflict:
                run_script("prepare_data.py", "--data-root", root / "source", "--ts", "0.0",
                           "--output-dir", root / "prepared")
            self.assertIn("Conflicting eigenvalue labels", conflict.exception.stderr)

    def test_training_checkpoint_predicts_in_original_units(self) -> None:
        with tempfile.TemporaryDirectory(prefix="single-ts-ml-training-") as directory:
            root = Path(directory)
            write_fixture(root / "source")
            prepared, run = root / "prepared", root / "run"
            run_script("prepare_data.py", "--data-root", root / "source", "--ts", "0.0",
                       "--output-dir", prepared)
            run_script("train_mlp.py", "--data-dir", prepared, "--output-dir", run,
                       "--epochs", "2", "--batch-size", "8", "--threads", "1",
                       "--patience", "2", "--hidden-sizes", "8", "8")
            metrics = json.loads((run / "metrics.json").read_text())
            self.assertEqual(metrics["Ts"], 0.0)
            self.assertEqual(metrics["split_sizes"], {"train": 32, "validation": 4, "test": 4})
            self.assertIn(metrics["best_epoch"], [1, 2])
            self.assertTrue(np.isfinite(metrics["test"]["mlp"]["standardized_mse"]))
            with np.load(run / "test_predictions.npz", allow_pickle=False) as predictions, \
                 np.load(prepared / "dataset.npz", allow_pickle=False) as data:
                indices = predictions["dataset_indices"]
                np.testing.assert_array_equal(indices, data["test_idx"])
                np.testing.assert_array_equal(predictions["y_true"], data["y_raw"][indices])
                np.testing.assert_array_equal(predictions["x_raw"], data["x_raw"][indices])
                np.testing.assert_allclose(predictions["error"], predictions["y_pred"] - predictions["y_true"])
                point = predictions["x_raw"][0]
                prediction = json.loads(run_script(
                    "predict.py", "--model", run / "best_model.pt", "--R", point[0],
                    "--alpha", point[1], "--beta", point[2],
                ).stdout)
                self.assertEqual(prediction["Ts"], 0.0)
                np.testing.assert_allclose(
                    [prediction["omega_r"], prediction["omega_i"]], predictions["y_pred"][0],
                    rtol=1e-5, atol=1e-7,
                )

    def test_scaler_fitted_on_all_samples_is_rejected(self) -> None:
        with tempfile.TemporaryDirectory(prefix="single-ts-ml-leakage-") as directory:
            root = Path(directory)
            write_fixture(root / "source")
            prepared = root / "prepared"
            run_script("prepare_data.py", "--data-root", root / "source", "--ts", "0.0",
                       "--output-dir", prepared)
            with np.load(prepared / "dataset.npz", allow_pickle=False) as archive:
                data = {key: archive[key].copy() for key in archive.files}
            # These transforms are internally consistent, but they used held-out rows.
            scalers = {}
            for prefix in ("x", "y"):
                raw = data[prefix + "_raw"]
                scalers[prefix + "_mean"] = raw.mean(axis=0)
                scalers[prefix + "_std"] = raw.std(axis=0)
                data[prefix] = ((raw - scalers[prefix + "_mean"]) / scalers[prefix + "_std"]).astype(np.float32)
            np.savez_compressed(prepared / "dataset.npz", **data)
            np.savez_compressed(prepared / "scaler.npz", **scalers)
            with self.assertRaises(subprocess.CalledProcessError) as leakage:
                run_script("train_mlp.py", "--data-dir", prepared, "--output-dir", root / "run",
                           "--epochs", "1", "--threads", "1")
            self.assertIn("training", leakage.exception.stderr.lower())

    def test_all_ts_filtering_split_scaling_and_prediction(self) -> None:
        """Identical triples at different Ts are distinct physical inputs."""
        with tempfile.TemporaryDirectory(prefix="all-ts-ml-test-") as directory:
            root = Path(directory)
            source = root / "source"
            files = [write_fixture(source, ts=value) for value in (0.0, 0.1)]
            before = [hashlib.sha256(path.read_bytes()).hexdigest() for path in files]
            prepared, repeated, single, run = (root / name for name in ("prepared", "repeat", "single", "run"))
            for output in (prepared, repeated):
                run_script("prepare_data.py", "--data-root", source, "--all-ts",
                           "--output-dir", output, "--seed", "42")
            run_script("prepare_data.py", "--data-root", source, "--ts", "0.0",
                       "--output-dir", single, "--seed", "42")
            self.assertEqual(before, [hashlib.sha256(path.read_bytes()).hexdigest() for path in files])
            metadata = json.loads((prepared / "metadata.json").read_text())
            self.assertEqual(metadata["input_names"], ["Ts", "R", "alpha", "beta"])
            self.assertIsNone(metadata["Ts"])
            self.assertEqual(metadata["Ts_values"], [0.0, 0.1])
            self.assertEqual(metadata["growth_rate_component"], "omega_i")
            self.assertEqual(metadata["counts"]["excluded_missing_sentinel"], 4)
            self.assertEqual(metadata["counts"]["excluded_paired_sentinel"], 2)
            self.assertEqual(metadata["counts"]["excluded_omega_i_minus_one_only"], 2)
            for counts in metadata["counts_by_ts"].values():
                self.assertEqual(counts["excluded_missing_sentinel"], 2)
                self.assertEqual(counts["excluded_paired_sentinel"], 1)
                self.assertEqual(counts["excluded_omega_i_minus_one_only"], 1)
            with np.load(prepared / "dataset.npz", allow_pickle=False) as data, \
                 np.load(repeated / "dataset.npz", allow_pickle=False) as repeat, \
                 np.load(single / "dataset.npz", allow_pickle=False) as old, \
                 np.load(prepared / "scaler.npz", allow_pickle=False) as scaler:
                x_raw, y_raw = data["x_raw"], data["y_raw"]
                self.assertEqual(x_raw.shape, (80, 4))
                self.assertEqual(y_raw.shape, (80, 2))
                self.assertEqual(len(np.unique(x_raw[:, 1:], axis=0)), 40)
                self.assertEqual(len(np.unique(x_raw, axis=0)), 80)
                self.assertFalse(np.any(x_raw[:, 1] == 999))
                self.assertFalse(np.any(np.all(y_raw == -1, axis=1)))
                self.assertFalse(np.any(y_raw[:, 1] == -1))
                self.assertTrue(np.any((y_raw[:, 0] == -1) & (y_raw[:, 1] == 0.125)))
                self.assertTrue(np.any(y_raw[:, 1] < 0))
                self.assertTrue(np.isfinite(data["x"]).all())
                self.assertTrue(np.isfinite(data["y"]).all())
                indices = [data[key] for key in ("train_idx", "val_idx", "test_idx")]
                self.assertEqual([len(index) for index in indices], [64, 8, 8])
                np.testing.assert_array_equal(np.sort(np.concatenate(indices)), np.arange(80))
                for key in ("train_idx", "val_idx", "test_idx", "x_raw", "y_raw", "x", "y"):
                    np.testing.assert_array_equal(data[key], repeat[key])
                for value in (0.0, 0.1):
                    self.assertEqual([int((x_raw[index, 0] == value).sum()) for index in indices], [32, 4, 4])
                # Per-Ts seeding preserves the original Ts=0 holdout membership.
                for key in ("train_idx", "val_idx", "test_idx"):
                    selection = data[key][x_raw[data[key], 0] == 0.0]
                    np.testing.assert_array_equal(x_raw[selection, 1:], old["x_raw"][old[key]])
                for raw, prefix in ((x_raw, "x"), (y_raw, "y")):
                    np.testing.assert_allclose(scaler[prefix + "_mean"], raw[data["train_idx"]].mean(axis=0))
                    np.testing.assert_allclose(scaler[prefix + "_std"], raw[data["train_idx"]].std(axis=0))
                    np.testing.assert_allclose(data[prefix], (raw - scaler[prefix + "_mean"]) / scaler[prefix + "_std"], rtol=2e-7, atol=2e-7)
            run_script("train_mlp.py", "--data-dir", prepared, "--output-dir", run,
                       "--epochs", "2", "--batch-size", "8", "--threads", "1",
                       "--patience", "2", "--hidden-sizes", "8", "8")
            metrics = json.loads((run / "metrics.json").read_text())
            self.assertEqual(metrics["feature_names"], ["Ts", "R", "alpha", "beta"])
            self.assertIsNone(metrics["Ts"])
            self.assertEqual(metrics["growth_rate_component"], "omega_i")
            with np.load(run / "test_predictions.npz", allow_pickle=False) as predictions:
                point = predictions["x_raw"][0]
                prediction = json.loads(run_script("predict.py", "--model", run / "best_model.pt",
                    "--Ts", point[0], "--R", point[1], "--alpha", point[2], "--beta", point[3]).stdout)
                self.assertEqual(prediction["input"], dict(zip(["Ts", "R", "alpha", "beta"], point)))
                np.testing.assert_allclose([prediction["omega_r"], prediction["omega_i"]],
                                           predictions["y_pred"][0], rtol=1e-5, atol=1e-7)
                with self.assertRaises(subprocess.CalledProcessError) as missing_ts:
                    run_script("predict.py", "--model", run / "best_model.pt",
                               "--R", point[1], "--alpha", point[2], "--beta", point[3])
                self.assertIn("Ts", missing_ts.exception.stderr)

    def test_growth_sign_uses_imaginary_component(self) -> None:
        from train_mlp import physical_metrics

        # Real-part signs are unchanged; exactly one imaginary-part sign is wrong.
        true = np.array([[-2.0, 0.2], [2.0, -0.2], [-3.0, 0.001], [3.0, 0.0]])
        pred = np.array([[-1.0, 0.3], [1.0, -0.3], [-1.0, -0.001], [1.0, 0.0]])
        result = physical_metrics(true, pred, np.ones(2), sign_tolerance=0.01)["growth_rate_sign"]
        self.assertIn("omega_i", result["definition"])
        self.assertEqual(result["error_count"], 1)
        self.assertEqual(result["error_rate"], 0.25)
        self.assertEqual(result["true_unstable_count"], 2)
        self.assertEqual(result["true_stable_count"], 1)
        self.assertEqual(result["true_neutral_count"], 1)
        self.assertEqual(result["error_rate_excluding_near_zero"], 0.0)


if __name__ == "__main__":
    unittest.main()
