from __future__ import annotations

"""Re-score the archived pRMST module against a coherent standard ATO truth.

The target treatment is the latent true strategy S*, not its misclassified code C.
The target population is defined by h(X2)=e_S(X2){1-e_S(X2)}, where
e_S(X2)=P(S*=1|X2). Because U is independent of X2 in the DGM and the
counterfactual outcome model depends on U but not X2, the target-standardized
truth reduces exactly to the population mixture P(U=1)=0.30. It therefore does
not vary with the observation-layer sensitivity/specificity.

The functional is the partial restricted mean time without cancer death,
integral_0^5 {1-CIF_c(t)} dt, in the presence of competing other-cause death.
It is not event-free RMST and does not treat other-cause death as the event.
"""

import hashlib
import json
import math
import platform
from pathlib import Path

import numpy as np
import pandas as pd


TAU = 5.0
SHAPE = 1.2
P_U1 = 0.30
EFF = 0.85

REPO_ROOT = Path(__file__).resolve().parents[2]
ROOT = REPO_ROOT
ARCHIVE = ROOT / "data" / "synthetic" / "rmst_archived_intervals"
OUT = ROOT / "analysis" / "rmst_standard_ato"


NONPH = {
    "early_benefit": {"cut": 1.0, "m1": 0.7, "m2": 1.0,
                      "final_label": "early benefit, then null"},
    # The archived label 'decay' was wrong: 0.6 -> 0.3 is stronger benefit.
    "decay": {"cut": 3.0, "m1": 0.6, "m2": 0.3,
              "final_label": "delayed intensification of benefit"},
    "crossing": {"cut": 2.0, "m1": 1.3, "m2": 0.8,
                 "final_label": "early harm, later benefit"},
}


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest().upper()


def cancer_multiplier(t: float | np.ndarray, treated: int, tag: str) -> np.ndarray:
    t = np.asarray(t)
    if treated == 0 or tag == "core":
        return np.ones_like(t, dtype=float) if treated == 0 else np.full_like(t, EFF, dtype=float)
    p = NONPH[tag]
    return np.where(t <= p["cut"], p["m1"], p["m2"]).astype(float)


def cum_hazard_cancer(t: float | np.ndarray, u: int, treated: int, uhr: float, tag: str) -> np.ndarray:
    t = np.asarray(t, dtype=float)
    base = math.exp(0.8) * (uhr ** u)
    if tag == "core":
        mult = EFF if treated else 1.0
        return base * mult * np.power(t, SHAPE)
    if not treated:
        return base * np.power(t, SHAPE)
    p = NONPH[tag]
    tk = np.power(t, SHAPE)
    ck = p["cut"] ** SHAPE
    return base * np.where(
        t <= p["cut"],
        p["m1"] * tk,
        p["m1"] * ck + p["m2"] * (tk - ck),
    )


def hazard_cancer(t: float | np.ndarray, u: int, treated: int, uhr: float, tag: str) -> np.ndarray:
    t = np.asarray(t, dtype=float)
    base = math.exp(0.8) * (uhr ** u)
    mult = cancer_multiplier(t, treated, tag)
    return SHAPE * base * mult * np.power(t, SHAPE - 1.0)


def cum_hazard_other(t: float | np.ndarray, u: int, uhr: float) -> np.ndarray:
    t = np.asarray(t, dtype=float)
    return math.exp(0.5) * (uhr ** u) * np.power(t, SHAPE)


def prmst_quad(u: int, treated: int, uhr: float, tag: str) -> tuple[float, float]:
    """Return pRMST and a nested Gauss-Legendre numerical-error estimate."""

    def loss_integrand(t: float) -> float:
        survival = math.exp(
            -float(cum_hazard_cancer(t, u, treated, uhr, tag))
            -float(cum_hazard_other(t, u, uhr))
        )
        return (TAU - t) * survival * float(hazard_cancer(t, u, treated, uhr, tag))

    def gauss_legendre(lo: float, hi: float, order: int) -> float:
        nodes, weights = np.polynomial.legendre.leggauss(order)
        x = 0.5 * (nodes + 1.0)
        # t=lo+(hi-lo)x^5 regularizes the t^(shape-1) behavior at zero.
        t = lo + (hi - lo) * np.power(x, 5)
        jacobian = 5.0 * (hi - lo) * np.power(x, 4)
        y = np.array([loss_integrand(float(v)) for v in t]) * jacobian
        return 0.5 * float(np.dot(weights, y))

    points = [] if tag == "core" or not treated else [NONPH[tag]["cut"]]
    intervals = [0.0] + points + [TAU]
    loss_fine = 0.0
    loss_coarse = 0.0
    for lo, hi in zip(intervals[:-1], intervals[1:]):
        loss_coarse += gauss_legendre(lo, hi, 256)
        loss_fine += gauss_legendre(lo, hi, 512)
    return TAU - loss_fine, abs(loss_fine - loss_coarse)


def prmst_dense(u: int, treated: int, uhr: float, tag: str, n_grid: int = 200_001) -> float:
    """Independent dense-grid check of the Fubini-integral implementation."""
    def segment(lo: float, hi: float) -> float:
        x = np.linspace(0.0, 1.0, n_grid)
        t = lo + (hi - lo) * np.power(x, 5)
        jacobian = 5.0 * (hi - lo) * np.power(x, 4)
        survival = np.exp(
            -cum_hazard_cancer(t, u, treated, uhr, tag)
            -cum_hazard_other(t, u, uhr)
        )
        y = (TAU - t) * survival * hazard_cancer(t, u, treated, uhr, tag) * jacobian
        h = 1.0 / (n_grid - 1)
        return float((h / 3.0) * (
            y[0] + y[-1] + 4.0 * np.sum(y[1:-1:2]) + 2.0 * np.sum(y[2:-1:2])
        ))

    points = [] if tag == "core" or not treated else [NONPH[tag]["cut"]]
    intervals = [0.0] + points + [TAU]
    return TAU - sum(segment(lo, hi) for lo, hi in zip(intervals[:-1], intervals[1:]))


def truth(uhr: float, tag: str) -> dict[str, float]:
    vals: dict[tuple[int, int], float] = {}
    errs: dict[tuple[int, int], float] = {}
    dense: dict[tuple[int, int], float] = {}
    for u in (0, 1):
        for treated in (0, 1):
            vals[(u, treated)], errs[(u, treated)] = prmst_quad(u, treated, uhr, tag)
            dense[(u, treated)] = prmst_dense(u, treated, uhr, tag)
    p0 = (1 - P_U1) * vals[(0, 0)] + P_U1 * vals[(1, 0)]
    p1 = (1 - P_U1) * vals[(0, 1)] + P_U1 * vals[(1, 1)]
    dense0 = (1 - P_U1) * dense[(0, 0)] + P_U1 * dense[(1, 0)]
    dense1 = (1 - P_U1) * dense[(0, 1)] + P_U1 * dense[(1, 1)]
    return {
        "truth_prmst_s0": p0,
        "truth_prmst_s1": p1,
        "truth_delta_prmst": p1 - p0,
        "quadrature_abs_error_bound": sum(errs.values()),
        "dense_grid_delta": dense1 - dense0,
        "quad_minus_dense": (p1 - p0) - (dense1 - dense0),
    }


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    old = pd.read_csv(ARCHIVE / "rmst_v2_summary.csv")

    truth_rows = []
    for uhr in sorted(old[old.tag.eq("core")].uhr.unique(), reverse=True):
        truth_rows.append({"uhr": uhr, "archived_tag": "core", "final_label": "proportional hazards", **truth(uhr, "core")})
    for tag in ("early_benefit", "decay", "crossing"):
        uhr = float(old.loc[old.tag.eq(tag), "uhr"].iloc[0])
        truth_rows.append({"uhr": uhr, "archived_tag": tag, "final_label": NONPH[tag]["final_label"], **truth(uhr, tag)})
    truth_df = pd.DataFrame(truth_rows)
    truth_df.to_csv(OUT / "rmst_standard_ato_truths.csv", index=False)

    rescored = []
    for row in old.itertuples(index=False):
        rep_path = ARCHIVE / f"rmst_scenario_{int(row.scenario_id)}.csv"
        rep = pd.read_csv(rep_path)
        tr = truth_df.loc[truth_df.archived_tag.eq(row.tag) & truth_df.uhr.eq(row.uhr)].iloc[0]
        target = float(tr.truth_delta_prmst)
        ok = rep[np.isfinite(rep.delta) & np.isfinite(rep.boot_lo) & np.isfinite(rep.boot_hi)].copy()
        errors = ok.delta - target
        coverage = float(((ok.boot_lo <= target) & (ok.boot_hi >= target)).mean())
        mean_delta = float(ok.delta.mean())
        bias = float(errors.mean())
        interval_width = ok.boot_hi - ok.boot_lo
        # Bias-shifted coverage is a descriptive interval-calibration diagnostic,
        # not an implementable correction and not a primary performance measure.
        centered_low = ok.boot_lo - bias
        centered_high = ok.boot_hi - bias
        centered_coverage = float(((centered_low <= target) & (centered_high >= target)).mean())
        rescored.append({
            "scenario_id": int(row.scenario_id), "se": row.se, "sp": row.sp,
            "uhr": row.uhr, "archived_tag": row.tag,
            "final_label": tr.final_label, "n": int(row.n), "R": int(row.R), "B": int(row.B),
            "truth_prmst_s0": tr.truth_prmst_s0, "truth_prmst_s1": tr.truth_prmst_s1,
            "truth_delta_prmst": target,
            "quadrature_abs_error_bound": tr.quadrature_abs_error_bound,
            "quad_minus_dense": tr.quad_minus_dense,
            "mean_delta": mean_delta, "bias": bias,
            "rmse": float(np.sqrt(np.mean(np.square(errors)))),
            "coverage": coverage,
            "coverage_mcse": math.sqrt(coverage * (1 - coverage) / len(ok)),
            "bias_shifted_coverage_descriptive": centered_coverage,
            "mean_interval_width": float(interval_width.mean()),
            "median_interval_width": float(interval_width.median()),
            "n_eval": len(ok), "min_boot_n": int(ok.boot_n.min()),
            "old_truth": row.truth_delta_prmst, "old_bias": row.bias,
            "old_coverage": row.coverage,
        })
    scored = pd.DataFrame(rescored).sort_values("scenario_id")
    scored["delta_truth_new_minus_old"] = scored.truth_delta_prmst - scored.old_truth
    scored["delta_bias_new_minus_old"] = scored.bias - scored.old_bias
    scored["delta_coverage_new_minus_old"] = scored.coverage - scored.old_coverage
    scored.to_csv(OUT / "rmst_standard_ato_rescored_12_scenarios.csv", index=False)

    provenance = {
        "estimand": "Difference in 5-y partial restricted mean time without cancer death under true strategies S*=1 vs S*=0 in the X2-defined true-strategy overlap population",
        "target_density": "h(X2)=e_S(X2){1-e_S(X2)}; U independent of X2 in DGM, so target truth averages U with P(U=1)=0.30",
        "competing_risk_functional": "integral_0^5 [1-CIF_c(t)] dt; other-cause death is competing, not part of the event definition",
        "truth_method": "512-node Gauss-Legendre quadrature of tau - integral_0^tau (tau-t) S_all(t) lambda_c(t) dt, with 256-node difference and dense-grid Simpson cross-check",
        "archived_interval_source": str(ARCHIVE),
        "archived_summary_sha256": sha256(ARCHIVE / "rmst_v2_summary.csv"),
        "script_sha256": sha256(Path(__file__)),
        "python": platform.python_version(), "numpy": np.__version__, "pandas": pd.__version__,
    }
    (OUT / "rmst_standard_ato_provenance.json").write_text(json.dumps(provenance, indent=2), encoding="utf-8")

    print(truth_df.to_string(index=False))
    print(scored[["scenario_id", "final_label", "se", "uhr", "truth_delta_prmst", "bias", "coverage", "coverage_mcse", "bias_shifted_coverage_descriptive"]].to_string(index=False))
    print(f"max |quad-dense| = {truth_df.quad_minus_dense.abs().max():.3e}")


if __name__ == "__main__":
    main()
