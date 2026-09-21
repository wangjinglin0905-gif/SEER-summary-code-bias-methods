from __future__ import annotations

"""Strict QBA reanalysis after correction of the Cox sandwich variance.

This script retains every independent oracle-calibration run, regenerates the
global mechanism-prior mapping with a genuine truncated Normal draw for uhr,
and applies Normal bias distributions analytically. The latter removes
unnecessary inner Monte Carlo noise: if B~N(mu_B, sd_B^2) and the observed
log-HR is approximately Normal, the corrected distribution is available in
closed form.
"""

import argparse
import gc
import hashlib
import json
import math
import platform
import time
from multiprocessing import Pool
from pathlib import Path

import numpy as np
import pandas as pd

from strict_cox_rerun import generate_core, overlap_weight


MASTER_SEED = 20260817
CAL_SEED_BASE = 998877
N_CAL = 1_000_000
R_CAL = 100
N_EXTERNAL = 200

REPO_ROOT = Path(__file__).resolve().parents[2]
ROOT = REPO_ROOT
STRICT = ROOT / "analysis" / "strict_cox_rerun"
OUT = ROOT / "analysis" / "strict_qba"
CAL_DIR = OUT / "oracle_calibration_raw"

GRIDS = [
    {"scenario_id": 301, "se": 1.00, "sp": 1.00, "uhr": 1.00, "source": 3},
    {"scenario_id": 302, "se": 1.00, "sp": 1.00, "uhr": 0.50, "source": 7},
    {"scenario_id": 303, "se": 1.00, "sp": 1.00, "uhr": 0.33, "source": 11},
    {"scenario_id": 304, "se": 0.90, "sp": 0.90, "uhr": 1.00, "source": 15},
    {"scenario_id": 305, "se": 0.90, "sp": 0.90, "uhr": 0.50, "source": 19},
    {"scenario_id": 306, "se": 0.90, "sp": 0.90, "uhr": 0.33, "source": 23},
    {"scenario_id": 307, "se": 0.75, "sp": 0.85, "uhr": 1.00, "source": 27},
    {"scenario_id": 308, "se": 0.75, "sp": 0.85, "uhr": 0.50, "source": 31},
    {"scenario_id": 309, "se": 0.75, "sp": 0.85, "uhr": 0.33, "source": 35},
]


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest().upper()


def weighted_cox_beta(time_obs, event, x, weight) -> float:
    """Weighted one-covariate Cox coefficient; continuous times imply no ties."""
    order = np.argsort(-np.asarray(time_obs), kind="mergesort")
    d = np.asarray(event, float)[order]
    z = np.asarray(x, float)[order]
    w = np.asarray(weight, float)[order]
    beta = 0.0
    for _ in range(50):
        risk = w * np.exp(beta * z)
        s0 = np.maximum(np.cumsum(risk), 1e-300)
        s1 = np.cumsum(risk * z)
        s2 = np.cumsum(risk * z * z)
        zbar = s1 / s0
        m = d == 1
        score = np.sum(w[m] * (z[m] - zbar[m]))
        info = np.sum(w[m] * (s2[m] / s0[m] - zbar[m] * zbar[m]))
        if not np.isfinite(info) or info <= 1e-12:
            raise FloatingPointError("nonpositive Cox information")
        step = score / info
        beta += step
        if abs(step) < 1e-10:
            break
    if not np.isfinite(beta):
        raise FloatingPointError("nonfinite Cox coefficient")
    return float(beta)


def one_bias(seed_entropy, n, se, sp, uhr, eff=0.85) -> tuple[float, int, int]:
    rng = np.random.default_rng(np.random.SeedSequence(seed_entropy))
    _, x2, _, code, time_obs, event, cause = generate_core(
        rng, n, uhr, se, sp, eff, -0.8, -0.5
    )
    weights = overlap_weight(x2, code)
    event_c = event * cause
    beta = weighted_cox_beta(time_obs, event_c, code, weights)
    return beta - math.log(eff), int(event_c.sum()), int(np.sum(event * (1 - cause)))


def calibrate_grid(grid: dict) -> dict:
    CAL_DIR.mkdir(parents=True, exist_ok=True)
    path = CAL_DIR / f"qba_{grid['scenario_id']}_oracle_100_runs.csv"
    if path.exists():
        cached = pd.read_csv(path)
        if len(cached) == R_CAL and cached.status.eq("ok").all():
            b = cached.bias_loghr.to_numpy()
            print(f"[{grid['scenario_id']}] reused 100 oracle calibration runs", flush=True)
            return {**grid, "mu_B": b.mean(), "sd_B": b.std(ddof=1),
                    "mcse_mu_B": b.std(ddof=1) / math.sqrt(len(b)),
                    "successful": len(b), "failed": 0, "cache_status": "reused"}
    rows = []
    t0 = time.time()
    for rep in range(R_CAL):
        ent = [CAL_SEED_BASE, int(round(grid["se"] * 100)), int(round(grid["sp"] * 100)),
               int(round(grid["uhr"] * 100)), 85, rep]
        try:
            b, nc, no = one_bias(ent, N_CAL, grid["se"], grid["sp"], grid["uhr"])
            rows.append({"scenario_id": grid["scenario_id"], "rep": rep,
                         "n": N_CAL, "bias_loghr": b, "n_cancer": nc,
                         "n_other": no, "status": "ok", "failure_reason": ""})
        except Exception as exc:
            rows.append({"scenario_id": grid["scenario_id"], "rep": rep,
                         "n": N_CAL, "bias_loghr": np.nan, "n_cancer": np.nan,
                         "n_other": np.nan, "status": "failed", "failure_reason": repr(exc)[:300]})
        if rep % 10 == 9:
            gc.collect()
    raw = pd.DataFrame(rows)
    raw.to_csv(path, index=False)
    ok = raw[raw.status.eq("ok")].bias_loghr.to_numpy()
    print(f"[{grid['scenario_id']}] oracle mu={ok.mean():+.5f} sd={ok.std(ddof=1):.5f} "
          f"ok={len(ok)} t={time.time()-t0:.1f}s", flush=True)
    return {**grid, "mu_B": ok.mean(), "sd_B": ok.std(ddof=1),
            "mcse_mu_B": ok.std(ddof=1) / math.sqrt(len(ok)),
            "successful": len(ok), "failed": R_CAL - len(ok), "cache_status": "new"}


def beta_params(mean: float, sd: float) -> tuple[float, float]:
    variance = sd * sd
    total = mean * (1 - mean) / variance - 1
    return mean * total, (1 - mean) * total


def truncated_normal(rng, mean, sd, low, high, n):
    pieces = []
    count = 0
    while count < n:
        draws = rng.normal(mean, sd, max(64, 2 * (n - count)))
        keep = draws[(draws >= low) & (draws <= high)]
        if len(keep):
            pieces.append(keep)
            count += len(keep)
    return np.concatenate(pieces)[:n]


def external_mapping() -> dict:
    path = OUT / "global_mechanism_prior_200_raw_draws.csv"
    if path.exists():
        cached = pd.read_csv(path)
        if len(cached) == N_EXTERNAL and cached.status.eq("ok").all() and cached.uhr_draw_method.eq("accept-reject truncated Normal").all():
            b = cached.bias_loghr.to_numpy()
            print("reused 200 global mechanism-prior draws", flush=True)
            return {"mu_B": b.mean(), "sd_B": b.std(ddof=1),
                    "mcse_mu_B": b.std(ddof=1) / math.sqrt(len(b)),
                    "successful": len(b), "failed": 0, "cache_status": "reused"}
    rng = np.random.default_rng(12345)
    a_se, b_se = beta_params(0.85, 0.05)
    a_sp, b_sp = beta_params(0.90, 0.03)
    sensitivity = rng.beta(a_se, b_se, N_EXTERNAL)
    specificity = rng.beta(a_sp, b_sp, N_EXTERNAL)
    uhr = truncated_normal(rng, 0.40, 0.10, 0.15, 0.80, N_EXTERNAL)
    rows = []
    t0 = time.time()
    for rep in range(N_EXTERNAL):
        ent = [CAL_SEED_BASE + 1, rep, int(round(uhr[rep] * 100)),
               int(round(sensitivity[rep] * 100)), int(round(specificity[rep] * 100))]
        try:
            b, nc, no = one_bias(ent, 200_000, sensitivity[rep], specificity[rep], uhr[rep])
            rows.append({"rep": rep, "sensitivity": sensitivity[rep], "specificity": specificity[rep],
                         "uhr": uhr[rep], "bias_loghr": b, "n_cancer": nc, "n_other": no,
                         "uhr_draw_method": "accept-reject truncated Normal", "status": "ok",
                         "failure_reason": ""})
        except Exception as exc:
            rows.append({"rep": rep, "sensitivity": sensitivity[rep], "specificity": specificity[rep],
                         "uhr": uhr[rep], "bias_loghr": np.nan, "n_cancer": np.nan, "n_other": np.nan,
                         "uhr_draw_method": "accept-reject truncated Normal", "status": "failed",
                         "failure_reason": repr(exc)[:300]})
    raw = pd.DataFrame(rows)
    raw.to_csv(path, index=False)
    ok = raw[raw.status.eq("ok")].bias_loghr.to_numpy()
    print(f"global mechanism prior mu={ok.mean():+.5f} sd={ok.std(ddof=1):.5f} "
          f"ok={len(ok)} t={time.time()-t0:.1f}s", flush=True)
    return {"mu_B": ok.mean(), "sd_B": ok.std(ddof=1),
            "mcse_mu_B": ok.std(ddof=1) / math.sqrt(len(ok)),
            "successful": len(ok), "failed": N_EXTERNAL - len(ok), "cache_status": "new"}


def evaluate(calibration: pd.DataFrame, external: dict) -> tuple[pd.DataFrame, pd.DataFrame]:
    true_log = math.log(0.85)
    summary_rows = []
    per_rows = []
    for grid in GRIDS:
        cal = calibration.loc[calibration.scenario_id.eq(grid["scenario_id"])].iloc[0]
        rep_path = STRICT / "replicates" / f"scenario_{grid['source']:03d}.csv"
        reps = pd.read_csv(rep_path)
        reps = reps[reps.status.eq("ok") & np.isfinite(reps.logHR_c) & np.isfinite(reps.se_c)].copy()
        priors = [
            ("oracle_cell_specific", float(cal.mu_B), float(cal.sd_B),
             "DGM-informed benchmark"),
            ("global_mechanism_prior", float(external["mu_B"]), float(external["sd_B"]),
             "single prior induced by specified Se/Sp/uhr distributions"),
            ("oracle_mean_x1.25", 1.25 * float(cal.mu_B), float(cal.sd_B),
             "relative mean misspecification; weak when oracle mean is near zero"),
            ("oracle_mean_x1.50", 1.50 * float(cal.mu_B), float(cal.sd_B),
             "relative mean misspecification; weak when oracle mean is near zero"),
        ]
        unadj_cover = ((reps.logHR_c - 1.96 * reps.se_c <= true_log) &
                       (reps.logHR_c + 1.96 * reps.se_c >= true_log)).mean()
        unadj_bias = float((reps.logHR_c - true_log).mean())
        for prior_name, mu_b, sd_b, interpretation in priors:
            corrected = reps.logHR_c - mu_b
            total_se = np.sqrt(np.square(reps.se_c) + sd_b * sd_b)
            lo = corrected - 1.96 * total_se
            hi = corrected + 1.96 * total_se
            error = corrected - true_log
            cover = (lo <= true_log) & (hi >= true_log)
            summary_rows.append({
                "scenario_id": grid["scenario_id"], "replica_source": grid["source"],
                "se": grid["se"], "sp": grid["sp"], "uhr": grid["uhr"], "true_hr": 0.85,
                "prior": prior_name, "prior_interpretation": interpretation,
                "prior_mu_B": mu_b, "prior_sd_B": sd_b,
                "bias_loghr": float(error.mean()), "rmse_loghr": float(np.sqrt(np.mean(error * error))),
                "coverage": float(cover.mean()),
                "width_ratio": float(np.mean(total_se) / np.mean(reps.se_c)),
                "mcse_coverage": float(math.sqrt(cover.mean() * (1 - cover.mean()) / len(reps))),
                "n_eval": len(reps), "unadjusted_bias_loghr": unadj_bias,
                "unadjusted_coverage": float(unadj_cover),
            })
            per_rows.append(pd.DataFrame({
                "scenario_id": grid["scenario_id"], "replica_source": grid["source"],
                "iter": reps["iter"].to_numpy(), "prior": prior_name,
                "prior_mu_B": mu_b, "prior_sd_B": sd_b,
                "logHR_obs": reps.logHR_c.to_numpy(), "se_obs": reps.se_c.to_numpy(),
                "logHR_corrected": corrected.to_numpy(), "se_corrected": total_se.to_numpy(),
                "ci_low_log": lo.to_numpy(), "ci_high_log": hi.to_numpy(),
                "covered": cover.astype(int).to_numpy(),
            }))
    return pd.DataFrame(summary_rows), pd.concat(per_rows, ignore_index=True)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--workers", type=int, default=4)
    args = parser.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    CAL_DIR.mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    with Pool(min(args.workers, len(GRIDS))) as pool:
        calibration = pd.DataFrame(pool.map(calibrate_grid, GRIDS, chunksize=1)).sort_values("scenario_id")
    calibration.to_csv(OUT / "oracle_calibration_summary_9_rows.csv", index=False)
    external = external_mapping()
    (OUT / "global_mechanism_prior_summary.json").write_text(json.dumps(external, indent=2), encoding="utf-8")
    summary, per_rep = evaluate(calibration, external)
    summary.to_csv(OUT / "strict_qba_summary_36_rows.csv", index=False)
    per_rep.to_csv(OUT / "strict_qba_per_replicate_180000_rows.csv", index=False)

    old_path = ROOT / "data" / "legacy" / "qba_v2_summary.csv"
    old = pd.read_csv(old_path)
    old_prior_map = {"oracle": "oracle_cell_specific", "externally_informed": "global_mechanism_prior",
                     "moderate_mis": "oracle_mean_x1.25", "severe_mis": "oracle_mean_x1.50"}
    old["prior"] = old.prior.map(old_prior_map)
    comparison = summary.merge(old[["scenario_id", "prior", "coverage", "bias_logHR", "rmse_logHR", "width_ratio"]],
                               on=["scenario_id", "prior"], how="left", suffixes=("_strict", "_v2"))
    for metric in ("coverage", "bias_loghr", "rmse_loghr", "width_ratio"):
        old_col = {"bias_loghr": "bias_logHR_v2", "rmse_loghr": "rmse_logHR_v2"}.get(metric, f"{metric}_v2")
        new_col = f"{metric}_strict"
        if new_col in comparison and old_col in comparison:
            comparison[f"delta_{metric}"] = comparison[new_col] - comparison[old_col]
    comparison.to_csv(OUT / "strict_qba_vs_v2_comparison.csv", index=False)

    provenance = {
        "completed_at_local": time.strftime("%Y-%m-%d %H:%M:%S"),
        "elapsed_seconds": time.time() - t0, "oracle_n": N_CAL,
        "oracle_repetitions": R_CAL, "external_mapping_draws": N_EXTERNAL,
        "external_uhr_distribution": "Normal(mean=0.40, sd=0.10) truncated by accept-reject to [0.15,0.80]",
        "qba_application": "analytic Normal convolution; no inner Monte Carlo",
        "strict_cox_source": str(STRICT), "script_sha256": sha256(Path(__file__)),
        "python": platform.python_version(), "numpy": np.__version__, "pandas": pd.__version__,
    }
    (OUT / "strict_qba_provenance.json").write_text(json.dumps(provenance, indent=2), encoding="utf-8")
    print(summary.groupby("prior")[["bias_loghr", "rmse_loghr", "coverage", "width_ratio"]].mean().to_string())
    print(json.dumps(provenance, indent=2))


if __name__ == "__main__":
    main()
