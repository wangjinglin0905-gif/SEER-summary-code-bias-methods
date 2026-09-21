from __future__ import annotations

"""Independent aggregation audit for the frozen SEER methods evidence.

This script does not regenerate the data-generating process. It rebuilds all
reported summaries from retained replicate-level outputs, verifies QBA
calibration distributions, checks the corrected pRMST truth definition, and
cross-checks the published SEER aggregate estimates. Patient-level SEER data
are neither read nor written.
"""

import hashlib
import json
import math
from pathlib import Path

import numpy as np
import pandas as pd


SOURCE = Path(__file__).resolve().parents[2]
QA = SOURCE / "verification" / "frozen_evidence"
QA.mkdir(parents=True, exist_ok=True)

STRICT = SOURCE / "analysis" / "strict_cox_rerun"
QBA = SOURCE / "analysis" / "strict_qba"
RMST = SOURCE / "analysis" / "rmst_standard_ato"
SEER = SOURCE / "analysis" / "seer_descriptive_sensitivity"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def scalar_close(a: float, b: float, atol: float = 2e-12) -> bool:
    return bool(np.isclose(float(a), float(b), rtol=0.0, atol=atol, equal_nan=True))


def add_check(checks: list[dict], name: str, passed: bool, observed, expected, note: str = "") -> None:
    checks.append(
        {
            "check": name,
            "passed": bool(passed),
            "observed": observed,
            "expected": expected,
            "note": note,
        }
    )


def audit_primary(checks: list[dict]) -> pd.DataFrame:
    registry = pd.read_csv(STRICT / "strict_primary_registry_72_rows.csv")
    published = pd.read_csv(STRICT / "strict_primary_summary_72_rows.csv").set_index("scenario_id")
    add_check(checks, "primary_registry_rows", len(registry) == 72, len(registry), 72)
    expected_eval = int((registry.planned_reps).sum())
    add_check(checks, "primary_planned_datasets", expected_eval == 252_000, expected_eval, 252_000)

    rows: list[dict] = []
    for reg in registry.itertuples(index=False):
        path = STRICT / "replicates" / f"scenario_{int(reg.scenario_id):03d}.csv"
        data = pd.read_csv(path)
        ok = data.loc[data.status.eq("ok")].copy()
        eff = float(reg.eff)
        true_log = math.log(eff)
        rebuilt = {
            "scenario_id": int(reg.scenario_id),
            "rows": len(data),
            "successful": len(ok),
            "failed": int((~data.status.eq("ok")).sum()),
            "median_HR_c": float(ok.hr_c.median()),
            "median_HR_naive": float(ok.hr_naive.median()),
            "median_HR_o": float(ok.hr_o.median()),
            "bias_logHR_c": float((ok.logHR_c - true_log).mean()),
            "median_HR_percent_error_c": float(100.0 * (ok.hr_c.median() - eff) / eff),
            "empirical_SE": float(ok.logHR_c.std(ddof=1)),
            "avg_robust_SE": float(ok.se_c.mean()),
            "avg_model_SE": float(ok.model_se_c.mean()),
            "rmse_c": float(np.sqrt(np.mean(np.square(ok.logHR_c - true_log)))),
            "coverage_c": float(((ok.ciL_c <= eff) & (ok.ciH_c >= eff)).mean()),
            "coverage_naive": float(((ok.ciL_naive <= eff) & (ok.ciH_naive >= eff)).mean()),
            "median_interval_width_log": float((2 * 1.96 * ok.se_c).median()),
            "directional_signal_rate": float((ok.ciH_o < 1).mean()),
            "median_cross_ratio": float(np.exp(ok.logHR_c - ok.logHR_o).median()),
            "max_abs_score_residual_sum": float(ok.max_abs_score_residual_sum.abs().max()),
        }
        rows.append(rebuilt)
        pub = published.loc[int(reg.scenario_id)]
        for field in [
            "successful", "failed", "median_HR_c", "median_HR_naive", "median_HR_o",
            "bias_logHR_c", "median_HR_percent_error_c", "empirical_SE", "avg_robust_SE",
            "avg_model_SE", "rmse_c", "coverage_c", "coverage_naive",
            "median_interval_width_log", "directional_signal_rate", "median_cross_ratio",
            "max_abs_score_residual_sum",
        ]:
            obs = rebuilt[field]
            exp = pub[field]
            exact = int(obs) == int(exp) if field in {"successful", "failed"} else scalar_close(obs, exp)
            add_check(checks, f"primary_{int(reg.scenario_id):03d}_{field}", exact, obs, exp)

    rebuilt_df = pd.DataFrame(rows)
    rebuilt_df.to_csv(QA / "primary_summary_independently_rebuilt.csv", index=False)
    add_check(
        checks,
        "primary_total_successful",
        int(rebuilt_df.successful.sum()) == 252_000,
        int(rebuilt_df.successful.sum()),
        252_000,
    )
    add_check(
        checks,
        "primary_total_failed",
        int(rebuilt_df.failed.sum()) == 0,
        int(rebuilt_df.failed.sum()),
        0,
    )
    max_score = float(rebuilt_df.max_abs_score_residual_sum.max())
    add_check(checks, "primary_score_residual_max_below_1e-8", max_score < 1e-8, max_score, "<1e-8")
    return rebuilt_df


def audit_h2(checks: list[dict]) -> pd.DataFrame:
    published = pd.read_csv(STRICT / "strict_h2_summary_7_rows.csv").set_index("scenario_id")
    add_check(checks, "h2_scenario_rows", len(published) == 7, len(published), 7)
    rows = []
    for sid, pub in published.iterrows():
        data = pd.read_csv(STRICT / "h2_replicates" / f"h2_{int(sid)}.csv")
        ok = data.loc[data.status.eq("ok")]
        rebuilt = {
            "scenario_id": int(sid),
            "successful": len(ok),
            "failed": int((~data.status.eq("ok")).sum()),
            "median_hr_o": float(ok.hr_o.median()),
            "two_sided_detection": float(((ok.ciH_o < 1) | (ok.ciL_o > 1)).mean()),
            "protective_direction_signal": float((ok.ciH_o < 1).mean()),
            "harmful_direction_signal": float((ok.ciL_o > 1).mean()),
        }
        rows.append(rebuilt)
        for field in ["successful", "failed", "median_hr_o", "two_sided_detection",
                      "protective_direction_signal", "harmful_direction_signal"]:
            obs, exp = rebuilt[field], pub[field]
            exact = int(obs) == int(exp) if field in {"successful", "failed"} else scalar_close(obs, exp)
            add_check(checks, f"h2_{int(sid)}_{field}", exact, obs, exp)
    rebuilt_df = pd.DataFrame(rows)
    rebuilt_df.to_csv(QA / "h2_summary_independently_rebuilt.csv", index=False)
    return rebuilt_df


def audit_qba(checks: list[dict]) -> pd.DataFrame:
    published = pd.read_csv(QBA / "strict_qba_summary_36_rows.csv")
    per = pd.read_csv(QBA / "strict_qba_per_replicate_180000_rows.csv")
    add_check(checks, "qba_summary_rows", len(published) == 36, len(published), 36)
    add_check(checks, "qba_per_replicate_rows", len(per) == 180_000, len(per), 180_000)

    true_log = math.log(0.85)
    rows = []
    for (sid, prior), group in per.groupby(["scenario_id", "prior"], sort=True):
        corrected_error = group.logHR_corrected - true_log
        rebuilt = {
            "scenario_id": int(sid),
            "prior": prior,
            "bias_loghr": float(corrected_error.mean()),
            "rmse_loghr": float(np.sqrt(np.mean(np.square(corrected_error)))),
            "coverage": float(group.covered.mean()),
            "width_ratio": float(group.se_corrected.mean() / group.se_obs.mean()),
            "n_eval": len(group),
        }
        rows.append(rebuilt)
    rebuilt_df = pd.DataFrame(rows)
    joined = rebuilt_df.merge(
        published[["scenario_id", "prior", "bias_loghr", "rmse_loghr", "coverage", "width_ratio", "n_eval"]],
        on=["scenario_id", "prior"], suffixes=("_rebuilt", "_published"), validate="one_to_one"
    )
    for row in joined.itertuples(index=False):
        key = f"qba_{int(row.scenario_id)}_{row.prior}"
        for field in ["bias_loghr", "rmse_loghr", "coverage", "width_ratio", "n_eval"]:
            obs = getattr(row, f"{field}_rebuilt")
            exp = getattr(row, f"{field}_published")
            exact = int(obs) == int(exp) if field == "n_eval" else scalar_close(obs, exp)
            add_check(checks, f"{key}_{field}", exact, obs, exp)
    rebuilt_df.to_csv(QA / "qba_summary_independently_rebuilt.csv", index=False)

    calibration = pd.read_csv(QBA / "oracle_calibration_summary_9_rows.csv").set_index("scenario_id")
    for sid, pub in calibration.iterrows():
        raw = pd.read_csv(QBA / "oracle_calibration_raw" / f"qba_{int(sid)}_oracle_100_runs.csv")
        bias = raw.loc[raw.status.eq("ok"), "bias_loghr"]
        for field, obs in {"mu_B": bias.mean(), "sd_B": bias.std(ddof=1), "successful": len(bias)}.items():
            exp = pub[field]
            exact = int(obs) == int(exp) if field == "successful" else scalar_close(obs, exp)
            add_check(checks, f"qba_calibration_{int(sid)}_{field}", exact, obs, exp)

    global_raw = pd.read_csv(QBA / "global_mechanism_prior_200_raw_draws.csv")
    global_summary = json.loads((QBA / "global_mechanism_prior_summary.json").read_text(encoding="utf-8"))
    bias = global_raw.loc[global_raw.status.eq("ok"), "bias_loghr"]
    for field, obs in {"mu_B": bias.mean(), "sd_B": bias.std(ddof=1), "successful": len(bias)}.items():
        exp = global_summary[field]
        exact = int(obs) == int(exp) if field == "successful" else scalar_close(obs, exp)
        add_check(checks, f"qba_global_{field}", exact, obs, exp)
    add_check(
        checks,
        "qba_global_uhr_bounds",
        bool(global_raw.uhr.between(0.15, 0.80, inclusive="both").all()),
        [float(global_raw.uhr.min()), float(global_raw.uhr.max())],
        "[0.15, 0.80]",
    )
    return rebuilt_df


def audit_rmst(checks: list[dict]) -> pd.DataFrame:
    data = pd.read_csv(RMST / "rmst_standard_ato_rescored_12_scenarios.csv")
    truths = pd.read_csv(RMST / "rmst_standard_ato_truths.csv")
    add_check(checks, "rmst_scenario_rows", len(data) == 12, len(data), 12)
    add_check(checks, "rmst_truth_rows", len(truths) == 6, len(truths), 6)
    core = data.loc[data.archived_tag.eq("core")]
    nunique = core.groupby("uhr").truth_delta_prmst.nunique(dropna=False)
    add_check(checks, "rmst_core_truth_invariant_to_se_sp", bool((nunique == 1).all()), nunique.to_dict(), "1 per uHR")
    add_check(checks, "rmst_dense_quadrature_agreement", float(truths.quad_minus_dense.abs().max()) <= 3e-14,
              float(truths.quad_minus_dense.abs().max()), "<=3e-14 years")
    add_check(checks, "rmst_bootstrap_count", int(data.min_boot_n.min()) == 1000,
              int(data.min_boot_n.min()), 1000)
    add_check(checks, "rmst_core_coverage_range", scalar_close(core.coverage.min(), 0.750, 1e-12)
              and scalar_close(core.coverage.max(), 0.956, 1e-12),
              [float(core.coverage.min()), float(core.coverage.max())], [0.750, 0.956])
    return data


def audit_cox_validation(checks: list[dict]) -> None:
    path = SOURCE / "analysis" / "cox_full_sandwich_validation_1200_fits.csv"
    data = pd.read_csv(path)
    add_check(checks, "cox_validation_fits", len(data) == 1200, len(data), 1200)
    max_beta = float((data.beta_full - data.beta_r).abs().max())
    max_se = float((data.se_full - data.se_r).abs().max())
    max_score = float(data.score_residual_sum.abs().max())
    add_check(checks, "cox_validation_max_beta_difference", max_beta < 3e-6, max_beta, "<3e-6")
    add_check(checks, "cox_validation_max_robust_se_difference", max_se < 3e-8, max_se, "<3e-8")
    add_check(checks, "cox_validation_score_residual_sum", max_score < 1e-8, max_score, "<1e-8")


def audit_seer_aggregates(checks: list[dict]) -> None:
    endpoint = pd.read_csv(SEER / "seer_endpoint_models.csv")
    diagnostics = pd.read_csv(SEER / "seer_diagnostics.csv")
    flow = pd.read_csv(SEER / "seer_cohort_flow.csv")
    targets = {
        ("CRC", "cancer death"): (54_035, 0.642, 0.578, 0.713),
        ("Bladder", "cancer death"): (36_500, 0.962, 0.845, 1.096),
    }
    endpoint_lower = endpoint.copy()
    endpoint_lower["cohort_key"] = endpoint_lower.cohort.astype(str).str.lower()
    endpoint_lower["endpoint_key"] = endpoint_lower.endpoint.astype(str).str.lower()
    for (cohort, endpoint_name), expected in targets.items():
        q = endpoint_lower.loc[
            endpoint_lower.cohort_key.str.contains(cohort.lower())
            & endpoint_lower.endpoint_key.str.contains(endpoint_name)
            & endpoint_lower.analysis.eq("primary unknown-category")
        ]
        if len(q) != 1:
            add_check(checks, f"seer_{cohort}_overall_unique", False, len(q), 1)
            continue
        row = q.iloc[0]
        observed = (int(row.n), round(float(row.hr), 3), round(float(row.lower), 3), round(float(row.upper), 3))
        add_check(checks, f"seer_{cohort}_cancer_hr", observed == expected, observed, expected)
    add_check(checks, "seer_flow_available", len(flow) > 0, len(flow), ">0")
    add_check(checks, "seer_diagnostics_available", len(diagnostics) > 0, len(diagnostics), ">0")


def main() -> None:
    checks: list[dict] = []
    primary = audit_primary(checks)
    h2 = audit_h2(checks)
    qba = audit_qba(checks)
    rmst = audit_rmst(checks)
    audit_cox_validation(checks)
    audit_seer_aggregates(checks)

    check_df = pd.DataFrame(checks)
    check_df.to_csv(QA / "verification_checks.csv", index=False)
    failed = check_df.loc[~check_df.passed]

    authoritative = [
        STRICT / "strict_primary_registry_72_rows.csv",
        STRICT / "strict_primary_summary_72_rows.csv",
        STRICT / "strict_h2_summary_7_rows.csv",
        QBA / "strict_qba_summary_36_rows.csv",
        QBA / "oracle_calibration_summary_9_rows.csv",
        QBA / "global_mechanism_prior_200_raw_draws.csv",
        RMST / "rmst_standard_ato_rescored_12_scenarios.csv",
        RMST / "rmst_standard_ato_truths.csv",
        SEER / "seer_endpoint_models.csv",
        SEER / "seer_joint_contrasts.csv",
        SOURCE / "analysis" / "cox_full_sandwich_validation_1200_fits.csv",
    ]
    manifest = [
        {"path": str(path), "bytes": path.stat().st_size, "sha256": sha256(path)}
        for path in authoritative
    ]
    (QA / "authoritative_file_manifest.json").write_text(
        json.dumps(manifest, indent=2, ensure_ascii=False), encoding="utf-8"
    )

    pooled = (
        pd.read_csv(QBA / "strict_qba_summary_36_rows.csv")
        .groupby("prior", sort=True)[["bias_loghr", "rmse_loghr", "coverage", "width_ratio"]]
        .mean()
        .reset_index()
    )
    pooled.to_csv(QA / "qba_pooled_means_rebuilt.csv", index=False)

    pooled_lines = [
        "| Prior | Bias (log HR) | RMSE (log HR) | Coverage | Width ratio |",
        "|---|---:|---:|---:|---:|",
    ]
    for row in pooled.itertuples(index=False):
        pooled_lines.append(
            f"| {row.prior} | {row.bias_loghr:.3f} | {row.rmse_loghr:.3f} | "
            f"{row.coverage:.3f} | {row.width_ratio:.3f} |"
        )

    report = [
        "# Frozen-evidence verification report",
        "",
        f"- Checks executed: {len(check_df):,}",
        f"- Checks passed: {int(check_df.passed.sum()):,}",
        f"- Checks failed: {len(failed):,}",
        f"- Primary retained datasets: {int(primary.successful.sum()):,}",
        f"- H2 retained datasets: {int(h2.successful.sum()):,}",
        f"- QBA evaluation rows: {int(qba.n_eval.sum()):,} across four priors",
        f"- pRMST scenario rows: {len(rmst):,}",
        "",
        "## Key rebuilt results",
        "",
        f"- Maximum primary score-residual sum: {primary.max_abs_score_residual_sum.max():.3e}",
        f"- Core pRMST coverage range: {rmst.loc[rmst.archived_tag.eq('core'), 'coverage'].min():.3f}–{rmst.loc[rmst.archived_tag.eq('core'), 'coverage'].max():.3f}",
        "- QBA pooled means:",
        "",
        *pooled_lines,
        "",
        "## Failed checks",
        "",
    ]
    if failed.empty:
        report.append("None.")
    else:
        report.append("See verification_checks.csv for the machine-readable failure ledger.")
    report.extend([
        "",
        "## Scope",
        "",
        "This audit independently re-aggregated frozen synthetic replicate-level outputs and checked frozen SEER aggregate tables. It did not re-read patient-level SEER files because the external G:/E: data volumes were not mounted during this audit.",
    ])
    (QA / "frozen_evidence_verification.md").write_text("\n".join(report) + "\n", encoding="utf-8")
    print(f"checks={len(check_df)} passed={int(check_df.passed.sum())} failed={len(failed)}")
    if not failed.empty:
        print(failed[["check", "observed", "expected"]].to_string(index=False))
        raise SystemExit(1)


if __name__ == "__main__":
    main()
