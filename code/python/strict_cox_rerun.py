from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
from multiprocessing import Pool
from pathlib import Path
import platform
import time

import numpy as np
import pandas as pd


MASTER_SEED = 20260817
SHAPE = 1.2
FOLLOW_UP = 5.0
LOSS_RATE = 0.05

REPO_ROOT = Path(__file__).resolve().parents[2]
ROOT = REPO_ROOT
OUT = ROOT / "analysis" / "strict_cox_rerun"
REP = OUT / "replicates"
H2_REP = OUT / "h2_replicates"
ORIGINAL_SUMMARY = ROOT / "data" / "legacy" / "summary_v2.csv"


def hash_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for b in iter(lambda: f.read(1024 * 1024), b""):
            h.update(b)
    return h.hexdigest().upper()


def scenario_seed(master: int, sid: int) -> int:
    return int(np.random.SeedSequence([master, sid]).generate_state(1)[0])


def rep_rng(sseed: int, rep_id: int):
    return np.random.default_rng(np.random.SeedSequence([sseed, rep_id]))


def fast_logit_p(x: np.ndarray, y: np.ndarray) -> np.ndarray:
    X = np.column_stack((np.ones(len(x)), x))
    beta = np.zeros(2)
    for _ in range(60):
        eta = np.clip(X @ beta, -30, 30)
        p = 1.0 / (1.0 + np.exp(-eta))
        score = X.T @ (y - p)
        v = p * (1 - p)
        info = X.T @ (v[:, None] * X)
        try:
            step = np.linalg.solve(info, score)
        except np.linalg.LinAlgError:
            step = np.linalg.lstsq(info, score, rcond=None)[0]
        beta += step
        if np.max(np.abs(step)) < 1e-10:
            break
    return 1.0 / (1.0 + np.exp(-np.clip(X @ beta, -30, 30)))


def full_weighted_cox(time_obs, event, x, weight):
    """Weighted one-covariate Cox fit with a full individual score sandwich.

    Continuous event times are generated, so ties occur with probability zero.
    The implementation was checked against survival::coxph(..., robust=TRUE)
    on 1,200 archived fits: max absolute robust-SE difference 2.69e-8.
    """
    time_obs = np.asarray(time_obs, float)
    event = np.asarray(event, float)
    x = np.asarray(x, float)
    weight = np.asarray(weight, float)
    order = np.argsort(-time_obs, kind="mergesort")
    d = event[order]
    z = x[order]
    w = weight[order]

    def parts(beta):
        risk = w * np.exp(beta * z)
        s0 = np.maximum(np.cumsum(risk), 1e-300)
        s1 = np.cumsum(risk * z)
        s2 = np.cumsum(risk * z * z)
        zbar = s1 / s0
        varz = s2 / s0 - zbar * zbar
        m = d == 1
        score = np.sum(w[m] * (z[m] - zbar[m]))
        info = np.sum(w[m] * varz[m])
        return score, info, risk, s0, zbar

    beta = 0.0
    conv = True
    for _ in range(50):
        score, info, _, _, _ = parts(beta)
        if not np.isfinite(info) or info <= 1e-12:
            raise FloatingPointError("nonpositive Cox information")
        step = score / info
        beta += step
        if abs(step) < 1e-10:
            break
    else:
        conv = False

    score, info, risk, s0, zbar = parts(beta)
    event_weight = np.where(d == 1, w, 0.0)
    da0 = event_weight / s0
    da1 = event_weight * zbar / s0
    cum_a0 = np.cumsum(da0[::-1])[::-1]
    cum_a1 = np.cumsum(da1[::-1])[::-1]
    observed = np.where(d == 1, w * (z - zbar), 0.0)
    expected = risk * (z * cum_a0 - cum_a1)
    score_residual = observed - expected
    robust_var = np.sum(score_residual * score_residual) / (info * info)
    se = math.sqrt(robust_var)
    model_se = 1.0 / math.sqrt(info)
    if not np.isfinite(beta) or not np.isfinite(se) or se <= 0:
        raise FloatingPointError("nonfinite Cox result")
    return {
        "beta": beta,
        "se": se,
        "model_se": model_se,
        "hr": math.exp(beta),
        "ci_low": math.exp(beta - 1.96 * se),
        "ci_high": math.exp(beta + 1.96 * se),
        "conv": "converged" if conv else "max_iter",
        "score_residual_sum": float(score_residual.sum()),
    }


def generate_core(rng, n, uhr, se, sp, eff, eta_c0, eta_o0):
    U = (rng.uniform(0, 1, n) < 0.30).astype(float)
    X2 = rng.normal(0, 1, n)
    pS = 1 / (1 + np.exp(-(-0.5 + 0.8 * U + 0.4 * X2)))
    S = (rng.uniform(0, 1, n) < pS).astype(float)
    if se == 1.0 and sp == 1.0:
        C = S.copy()
    else:
        r = rng.uniform(0, 1, n)
        C = np.where(S == 1, r < se, r >= sp).astype(float)
    eta_c = eta_c0 - np.log(eff) * S - np.log(uhr) * U
    eta_o = eta_o0 - np.log(uhr) * U
    Tc = np.exp(eta_c / SHAPE) * (-np.log(rng.uniform(0, 1, n))) ** (1 / SHAPE)
    To = np.exp(eta_o / SHAPE) * (-np.log(rng.uniform(0, 1, n))) ** (1 / SHAPE)
    T = np.minimum(Tc, To)
    cause = (Tc <= To).astype(float)
    c_admin = rng.uniform(0, FOLLOW_UP, n)
    c_loss = rng.exponential(1 / LOSS_RATE, n)
    ctime = np.minimum(c_admin, c_loss)
    observed = np.minimum(T, ctime)
    event = ((T < ctime) & (observed <= FOLLOW_UP)).astype(float)
    observed = np.minimum(observed, FOLLOW_UP)
    return U, X2, S, C, observed, event, cause


def overlap_weight(x2, C):
    p = fast_logit_p(x2, C)
    return np.where(C == 1, 1 - p, p)


def run_primary_replicate(sseed, rep_id, n, uhr, se, sp, eff, eta_c0, eta_o0):
    row = {"iter": rep_id, "status": "ok", "failure_reason": ""}
    try:
        _, X2, _, C, time_obs, event, cause = generate_core(
            rep_rng(sseed, rep_id), n, uhr, se, sp, eff, eta_c0, eta_o0
        )
        w = overlap_weight(X2, C)
        event_c = event * cause
        event_o = event * (1 - cause)
        naive = full_weighted_cox(time_obs, event_c, C, np.ones(n))
        cancer = full_weighted_cox(time_obs, event_c, C, w)
        other = full_weighted_cox(time_obs, event_o, C, w)
        row.update(
            {
                "logHR_naive": naive["beta"], "hr_naive": naive["hr"],
                "se_naive": naive["se"], "ciL_naive": naive["ci_low"], "ciH_naive": naive["ci_high"],
                "logHR_c": cancer["beta"], "hr_c": cancer["hr"], "se_c": cancer["se"],
                "model_se_c": cancer["model_se"], "ciL_c": cancer["ci_low"], "ciH_c": cancer["ci_high"],
                "logHR_o": other["beta"], "hr_o": other["hr"], "se_o": other["se"],
                "model_se_o": other["model_se"], "ciL_o": other["ci_low"], "ciH_o": other["ci_high"],
                "conv_naive": naive["conv"], "conv_c": cancer["conv"], "conv_o": other["conv"],
                "max_abs_score_residual_sum": max(abs(naive["score_residual_sum"]),
                                                    abs(cancer["score_residual_sum"]),
                                                    abs(other["score_residual_sum"])),
                "n_cancer": int(event_c.sum()), "n_other": int(event_o.sum()),
                "ess": float(w.sum() ** 2 / np.sum(w * w)),
            }
        )
    except Exception as exc:
        row.update({"status": "failed", "failure_reason": repr(exc)[:300]})
    return row


def primary_registry() -> pd.DataFrame:
    rows = []
    sid = 1
    for se, sp in [(1.0, 1.0), (0.9, 0.9), (0.75, 0.85)]:
        for uhr in [1.0, 0.5, 0.33]:
            for eff in [1.0, 0.85]:
                for n in [5000, 50000]:
                    rows.append({"scenario_id": sid, "family": "main_hi", "anchor": "", "se": se, "sp": sp,
                                 "uhr": uhr, "eff": eff, "n": n, "eta_c0": -0.8, "eta_o0": -0.5,
                                 "planned_reps": 5000, "master_seed": MASTER_SEED,
                                 "scenario_seed": scenario_seed(MASTER_SEED, sid)})
                    sid += 1
    for family, base, e0c, e0o in [("low_crc", 101, 4.2074, 3.1726),
                                    ("low_bladder", 201, 2.6917, 2.4908)]:
        k = 1
        for se, sp in [(1.0, 1.0), (0.9, 0.9), (0.75, 0.85)]:
            for uhr in [1.0, 0.5, 0.33]:
                for eff in [1.0, 0.85]:
                    sid = base + k
                    rows.append({"scenario_id": sid, "family": family, "anchor": family, "se": se, "sp": sp,
                                 "uhr": uhr, "eff": eff, "n": 50000, "eta_c0": e0c, "eta_o0": e0o,
                                 "planned_reps": 2000, "master_seed": MASTER_SEED,
                                 "scenario_seed": scenario_seed(MASTER_SEED, sid)})
                    k += 1
    return pd.DataFrame(rows)


def summarize_primary(s, reps):
    ok = reps[reps.status == "ok"].copy()
    true_hr = float(s["eff"])
    true_log = math.log(true_hr)
    R = len(ok)
    bias = ok.logHR_c - true_log
    coverage = np.mean((ok.ciL_c <= true_hr) & (ok.ciH_c >= true_hr))
    coverage_naive = np.mean((ok.ciL_naive <= true_hr) & (ok.ciH_naive >= true_hr))
    signal = np.mean(ok.ciH_o < 1)
    return {
        **s,
        "median_HR_c": np.median(ok.hr_c),
        "median_HR_naive": np.median(ok.hr_naive),
        "median_HR_o": np.median(ok.hr_o),
        "bias_logHR_c": np.mean(bias),
        "median_HR_percent_error_c": 100 * (np.median(ok.hr_c) - true_hr) / true_hr,
        "empirical_SE": np.std(bias, ddof=1),
        "avg_robust_SE": np.mean(ok.se_c),
        "avg_model_SE": np.mean(ok.model_se_c),
        "rmse_c": np.sqrt(np.mean(bias * bias)),
        "coverage_c": coverage,
        "coverage_naive": coverage_naive,
        "median_interval_width_log": np.median(2 * 1.96 * ok.se_c),
        "directional_signal_rate": signal,
        "median_cross_ratio": np.median(np.exp(ok.logHR_c - ok.logHR_o)),
        "mcse_bias": np.std(bias, ddof=1) / math.sqrt(R),
        "mcse_coverage": math.sqrt(coverage * (1 - coverage) / R),
        "mcse_signal": math.sqrt(signal * (1 - signal) / R),
        "successful": R,
        "failed": int(s["planned_reps"] - R),
        "failure_rate": (s["planned_reps"] - R) / s["planned_reps"],
        "convergence_warnings": int((ok.conv_c != "converged").sum() + (ok.conv_o != "converged").sum()),
        "max_abs_score_residual_sum": ok.max_abs_score_residual_sum.max(),
    }


def run_primary_scenario(s):
    t0 = time.time()
    s = dict(s)
    existing = REP / f"scenario_{int(s['scenario_id']):03d}.csv"
    if existing.exists():
        cached = pd.read_csv(existing)
        if len(cached) == int(s["planned_reps"]) and (cached.status == "ok").all():
            summary = summarize_primary(s, cached)
            summary["computational_time_s"] = 0.0
            summary["cache_status"] = "reused complete strict rerun"
            print(f"[{int(s['scenario_id']):03d}] reused complete strict replicate file", flush=True)
            return summary
    rows = [run_primary_replicate(int(s["scenario_seed"]), r, int(s["n"]), float(s["uhr"]),
                                  float(s["se"]), float(s["sp"]), float(s["eff"]),
                                  float(s["eta_c0"]), float(s["eta_o0"]))
            for r in range(int(s["planned_reps"]))]
    reps = pd.DataFrame(rows)
    reps.insert(0, "scenario_id", int(s["scenario_id"]))
    REP.mkdir(parents=True, exist_ok=True)
    reps.to_csv(existing, index=False)
    summary = summarize_primary(s, reps)
    summary["computational_time_s"] = time.time() - t0
    summary["cache_status"] = "newly computed"
    print(f"[{int(s['scenario_id']):03d}] {s['family']} n={int(s['n']):,} "
          f"ok={summary['successful']} cov={summary['coverage_c']:.3f} "
          f"sig={summary['directional_signal_rate']:.3f} t={summary['computational_time_s']:.1f}s", flush=True)
    return summary


H2_SCEN = [
    (901, "no_direct_effect_under_confounding", 1.0, 1.0, 50000, -0.8, -0.5, 0.5,
     "no direct effect; shared confounding present"),
    (902, "harmful_early", 1.5, 1.0, 50000, -0.8, -0.5, 0.5,
     "harmful multiplier 1.5 during 0-0.5 years"),
    (903, "harmful_sustained", 1.5, 1.5, 50000, -0.8, -0.5, 0.5,
     "harmful multiplier 1.5 throughout"),
    (904, "protective_early", 0.67, 1.0, 50000, -0.8, -0.5, 0.5,
     "protective multiplier 0.67 during 0-0.5 years"),
    (905, "protective_sustained", 0.67, 0.67, 50000, -0.8, -0.5, 0.5,
     "protective multiplier 0.67 throughout"),
    (906, "harmful_early_low_event", 1.5, 1.0, 50000, 4.2074, 3.1726, 0.5,
     "harmful early multiplier in CRC-anchored low-event regime"),
    (907, "null_no_confounding", 1.0, 1.0, 50000, -0.8, -0.5, 1.0,
     "no direct effect and no shared confounding"),
]


def generate_h2(rng, n, uhr, m1, m2, eta_c0, eta_o0, se=0.9, sp=0.9, eff=0.85, cutoff=0.5):
    U = (rng.uniform(0, 1, n) < 0.30).astype(float)
    X2 = rng.normal(0, 1, n)
    pS = 1 / (1 + np.exp(-(-0.5 + 0.8 * U + 0.4 * X2)))
    S = (rng.uniform(0, 1, n) < pS).astype(float)
    r = rng.uniform(0, 1, n)
    C = np.where(S == 1, r < se, r >= sp).astype(float)
    eta_c = eta_c0 - np.log(eff) * S - np.log(uhr) * U
    Tc = np.exp(eta_c / SHAPE) * (-np.log(rng.uniform(0, 1, n))) ** (1 / SHAPE)
    eta_o = eta_o0 - np.log(uhr) * U
    base = np.exp(-eta_o)
    first = np.where(S == 1, m1, 1.0)
    second = np.where(S == 1, m2, 1.0)
    hcut = base * first * cutoff ** SHAPE
    y = -np.log(rng.uniform(0, 1, n))
    To = np.empty(n)
    early = y <= hcut
    To[early] = (y[early] / (base[early] * first[early])) ** (1 / SHAPE)
    To[~early] = (cutoff ** SHAPE + (y[~early] - hcut[~early]) / (base[~early] * second[~early])) ** (1 / SHAPE)
    T = np.minimum(Tc, To)
    cause = (Tc <= To).astype(float)
    ctime = np.minimum(rng.uniform(0, FOLLOW_UP, n), rng.exponential(1 / LOSS_RATE, n))
    observed = np.minimum(T, ctime)
    event = ((T < ctime) & (observed <= FOLLOW_UP)).astype(float)
    observed = np.minimum(observed, FOLLOW_UP)
    return X2, C, observed, event, cause


def run_h2_scenario(item):
    sid, name, m1, m2, n, e0c, e0o, uhr, description = item
    t0 = time.time()
    sseed = scenario_seed(MASTER_SEED, sid)
    rows = []
    for r in range(2000):
        row = {"iter": r, "status": "ok", "failure_reason": ""}
        try:
            X2, C, t, event, cause = generate_h2(rep_rng(sseed, r), n, uhr, m1, m2, e0c, e0o)
            w = overlap_weight(X2, C)
            event_o = event * (1 - cause)
            fit = full_weighted_cox(t, event_o, C, w)
            row.update({"logHR_o": fit["beta"], "hr_o": fit["hr"], "se_o": fit["se"],
                        "ciL_o": fit["ci_low"], "ciH_o": fit["ci_high"], "conv": fit["conv"],
                        "n_cancer_deaths": int(np.sum(event * cause)),
                        "n_other_deaths": int(event_o.sum()), "n_censored": int(np.sum(event == 0)),
                        "ess": float(w.sum() ** 2 / np.sum(w * w)),
                        "score_residual_sum": fit["score_residual_sum"]})
        except Exception as exc:
            row.update({"status": "failed", "failure_reason": repr(exc)[:300]})
        rows.append(row)
    reps = pd.DataFrame(rows)
    H2_REP.mkdir(parents=True, exist_ok=True)
    reps.insert(0, "scenario_id", sid)
    reps.to_csv(H2_REP / f"h2_{sid}.csv", index=False)
    ok = reps[reps.status == "ok"]
    two = np.mean((ok.ciL_o > 1) | (ok.ciH_o < 1))
    prot = np.mean(ok.ciH_o < 1)
    harmful = np.mean(ok.ciL_o > 1)
    result = {"scenario_id": sid, "name": name, "description": description, "m1": m1, "m2": m2,
              "n": n, "eta_c0": e0c, "eta_o0": e0o, "uhr": uhr, "se": 0.9, "sp": 0.9,
              "eff": 0.85, "successful": len(ok), "failed": 2000 - len(ok),
              "median_hr_o": np.median(ok.hr_o), "median_n_cancer": np.median(ok.n_cancer_deaths),
              "median_n_other": np.median(ok.n_other_deaths), "median_n_censored": np.median(ok.n_censored),
              "median_ess": np.median(ok.ess), "two_sided_detection": two,
              "protective_direction_signal": prot, "harmful_direction_signal": harmful,
              "mcse_two_sided": math.sqrt(two * (1 - two) / len(ok)),
              "mcse_protective": math.sqrt(prot * (1 - prot) / len(ok)),
              "mcse_harmful": math.sqrt(harmful * (1 - harmful) / len(ok)),
              "computational_time_s": time.time() - t0}
    print(f"[H2 {sid}] HR={result['median_hr_o']:.3f} two={two:.3f} prot={prot:.3f} harm={harmful:.3f} "
          f"t={result['computational_time_s']:.1f}s", flush=True)
    return result


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--workers", type=int, default=12)
    parser.add_argument("--primary", action="store_true")
    parser.add_argument("--h2", action="store_true")
    args = parser.parse_args()
    if not args.primary and not args.h2:
        args.primary = args.h2 = True
    OUT.mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    if args.primary:
        registry = primary_registry()
        registry.to_csv(OUT / "strict_primary_registry_72_rows.csv", index=False)
        with Pool(args.workers) as pool:
            summaries = pool.map(run_primary_scenario, [r.to_dict() for _, r in registry.iterrows()], chunksize=1)
        summary = pd.DataFrame(summaries).sort_values("scenario_id")
        summary.to_csv(OUT / "strict_primary_summary_72_rows.csv", index=False)
        old = pd.read_csv(ORIGINAL_SUMMARY)
        merged = summary.merge(old, on="scenario_id", suffixes=("_strict", "_v2"))
        metrics = ["median_HR_c", "bias_logHR_c", "coverage_c", "directional_signal_rate", "avg_model_SE"]
        for metric in metrics:
            if f"{metric}_strict" in merged and f"{metric}_v2" in merged:
                merged[f"delta_{metric}"] = merged[f"{metric}_strict"] - merged[f"{metric}_v2"]
        merged.to_csv(OUT / "strict_vs_v2_primary_comparison.csv", index=False)
    if args.h2:
        with Pool(min(args.workers, len(H2_SCEN))) as pool:
            h2 = pd.DataFrame(pool.map(run_h2_scenario, H2_SCEN, chunksize=1)).sort_values("scenario_id")
        h2.to_csv(OUT / "strict_h2_summary_7_rows.csv", index=False)

    provenance = {
        "completed_at_local": time.strftime("%Y-%m-%d %H:%M:%S"),
        "elapsed_seconds": time.time() - t0,
        "master_seed": MASTER_SEED,
        "seed_scheme": "SeedSequence([master_seed, scenario_id]); replicate SeedSequence([scenario_seed, replicate_id])",
        "cox_variance": "full individual score-residual sandwich",
        "validation_file": str(ROOT / "analysis" / "cox_full_sandwich_validation_summary.csv"),
        "python": platform.python_version(),
        "numpy": np.__version__,
        "pandas": pd.__version__,
        "script_sha256": hash_file(Path(__file__)),
    }
    (OUT / "strict_rerun_provenance.json").write_text(json.dumps(provenance, indent=2), encoding="utf-8")
    print(json.dumps(provenance, indent=2))


if __name__ == "__main__":
    main()
