from __future__ import annotations

"""Recompute the custom score-residual sandwich and compare with a fresh R run."""

import json
from pathlib import Path

import numpy as np
import pandas as pd


REPO_ROOT = Path(__file__).resolve().parents[2]
SOURCE = REPO_ROOT / "data" / "large" / "cox_validation_source"
OUT = REPO_ROOT / "verification" / "cox_reverification"


def cox_full_sandwich(time, event, x, weight):
    time = np.asarray(time, float)
    event = np.asarray(event, float)
    x = np.asarray(x, float)
    case_weight = np.asarray(weight, float)
    order = np.argsort(-time, kind="mergesort")
    d = event[order]
    z = x[order]
    w = case_weight[order]

    def components(beta):
        risk = w * np.exp(beta * z)
        s0 = np.maximum(np.cumsum(risk), 1e-300)
        s1 = np.cumsum(risk * z)
        s2 = np.cumsum(risk * z * z)
        zbar = s1 / s0
        variance = s2 / s0 - zbar * zbar
        events = d == 1
        score = np.sum(w[events] * (z[events] - zbar[events]))
        information = np.sum(w[events] * variance[events])
        return score, information, s0, zbar

    beta = 0.0
    for _ in range(50):
        score, information, _, _ = components(beta)
        step = score / information
        beta += step
        if abs(step) < 1e-10:
            break

    _, information, s0, zbar = components(beta)
    risk = w * np.exp(beta * z)
    event_weight = np.where(d == 1, w, 0.0)
    da0 = event_weight / s0
    da1 = event_weight * zbar / s0
    cumulative_a0 = np.cumsum(da0[::-1])[::-1]
    cumulative_a1 = np.cumsum(da1[::-1])[::-1]
    observed = np.where(d == 1, w * (z - zbar), 0.0)
    expected = risk * (z * cumulative_a0 - cumulative_a1)
    residual = observed - expected
    robust_se = np.sqrt(np.sum(np.square(residual))) / information
    return float(beta), float(robust_se), float(residual.sum())


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    r_estimates = pd.read_csv(OUT / "R_coxph_estimates_rerun.csv")
    rows = []
    for combo in range(24):
        data = pd.read_csv(SOURCE / f"data_combo_{combo:02d}.csv")
        for dataset_id, group in data.groupby("dataset_id", sort=True):
            beta, se, residual_sum = cox_full_sandwich(
                group.time, group.event_c, group.C, group.w
            )
            rows.append(
                {
                    "combo": combo,
                    "dataset_id": int(dataset_id),
                    "beta_python": beta,
                    "se_python": se,
                    "score_residual_sum": residual_sum,
                }
            )
    result = pd.DataFrame(rows).merge(
        r_estimates, on=["combo", "dataset_id"], validate="one_to_one"
    )
    result["abs_beta_difference"] = (result.beta_python - result.beta_r).abs()
    result["abs_se_difference"] = (result.se_python - result.se_r).abs()
    result.to_csv(OUT / "cox_python_vs_R_rerun_1200_fits.csv", index=False)

    summary = {
        "fits": len(result),
        "max_abs_beta_difference": float(result.abs_beta_difference.max()),
        "max_abs_robust_se_difference": float(result.abs_se_difference.max()),
        "median_abs_robust_se_difference": float(result.abs_se_difference.median()),
        "max_abs_score_residual_sum": float(result.score_residual_sum.abs().max()),
        "passed": bool(
            len(result) == 1200
            and result.abs_beta_difference.max() < 3e-6
            and result.abs_se_difference.max() < 3e-8
            and result.score_residual_sum.abs().max() < 1e-8
        ),
    }
    (OUT / "cox_reverification_summary.json").write_text(
        json.dumps(summary, indent=2), encoding="utf-8"
    )
    print(json.dumps(summary, indent=2))
    if not summary["passed"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
