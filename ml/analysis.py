"""
Statistical Analysis for 5G-Radar Coexistence Experiments.

Functions for:
  - Descriptive statistics (mean, std, CI)
  - Significance tests (t-test, Wilcoxon, ANOVA)
  - Results table generation
  - Key findings extraction
"""

import numpy as np
import json
from pathlib import Path
from typing import Any

try:
    from scipy import stats as sp_stats
except ImportError:
    sp_stats = None


# ---------------------------------------------------------------------------
# 1. Descriptive statistics
# ---------------------------------------------------------------------------

def compute_statistics(values: list[float] | np.ndarray, confidence: float = 0.95) -> dict:
    """Compute mean, std, median, and confidence interval.

    Parameters
    ----------
    values : array-like
    confidence : float, default 0.95

    Returns
    -------
    dict with mean, std, median, ci_lower, ci_upper, n
    """
    arr = np.asarray(values, dtype=np.float64)
    n = len(arr)
    if n == 0:
        return {"mean": 0, "std": 0, "median": 0, "ci_lower": 0, "ci_upper": 0, "n": 0}

    mean = float(arr.mean())
    std = float(arr.std(ddof=1)) if n > 1 else 0.0
    median = float(np.median(arr))

    if n > 1 and sp_stats is not None:
        t_crit = sp_stats.t.ppf((1 + confidence) / 2, n - 1)
        margin = t_crit * std / np.sqrt(n)
    else:
        margin = 0.0

    return {
        "mean": mean,
        "std": std,
        "median": median,
        "ci_lower": mean - margin,
        "ci_upper": mean + margin,
        "n": n,
    }


# ---------------------------------------------------------------------------
# 2. Significance tests
# ---------------------------------------------------------------------------

def significance_test(
    baseline: list[float] | np.ndarray,
    proposed: list[float] | np.ndarray,
    test: str = "t-test",
    alpha: float = 0.05,
) -> dict:
    """Two-sample significance test.

    Parameters
    ----------
    baseline, proposed : array-like
    test : 't-test' | 'wilcoxon' | 'mannwhitney'
    alpha : significance level

    Returns
    -------
    dict with statistic, p_value, significant, effect_size (Cohen's d)
    """
    if sp_stats is None:
        return {"error": "scipy not installed"}

    b = np.asarray(baseline, dtype=np.float64)
    p = np.asarray(proposed, dtype=np.float64)

    if test == "t-test":
        stat, pval = sp_stats.ttest_ind(b, p, equal_var=False)
    elif test == "wilcoxon":
        if len(b) == len(p):
            stat, pval = sp_stats.wilcoxon(b, p)
        else:
            stat, pval = sp_stats.mannwhitneyu(b, p, alternative="two-sided")
    elif test == "mannwhitney":
        stat, pval = sp_stats.mannwhitneyu(b, p, alternative="two-sided")
    else:
        raise ValueError(f"Unknown test: {test}")

    # Cohen's d
    pooled_std = np.sqrt((b.var(ddof=1) + p.var(ddof=1)) / 2)
    cohens_d = (p.mean() - b.mean()) / pooled_std if pooled_std > 0 else 0.0

    return {
        "test": test,
        "statistic": float(stat),
        "p_value": float(pval),
        "significant": bool(pval < alpha),
        "effect_size_cohens_d": float(cohens_d),
        "alpha": alpha,
    }


def anova_test(groups: dict[str, list[float]], alpha: float = 0.05) -> dict:
    """One-way ANOVA across multiple groups.

    Parameters
    ----------
    groups : dict mapping group name to list of values
    alpha : significance level

    Returns
    -------
    dict with F-statistic, p_value, significant
    """
    if sp_stats is None:
        return {"error": "scipy not installed"}

    arrays = [np.asarray(v, dtype=np.float64) for v in groups.values()]
    f_stat, pval = sp_stats.f_oneway(*arrays)

    return {
        "test": "one-way ANOVA",
        "groups": list(groups.keys()),
        "F_statistic": float(f_stat),
        "p_value": float(pval),
        "significant": bool(pval < alpha),
        "alpha": alpha,
    }


# ---------------------------------------------------------------------------
# 3. Results table generation
# ---------------------------------------------------------------------------

def generate_results_table(
    results: dict[str, dict],
    metrics: list[str] | None = None,
) -> str:
    """Generate a markdown results table.

    Parameters
    ----------
    results : dict mapping method name to metric dict
        Each metric dict has keys like 'throughput', 'bler', etc.
        Values can be float or dict with 'mean'/'std'.
    metrics : list of metric names to include

    Returns
    -------
    str – markdown table
    """
    if metrics is None:
        # Gather all metric keys
        all_keys = set()
        for v in results.values():
            all_keys.update(v.keys())
        metrics = sorted(all_keys)

    # Header
    header = "| Method | " + " | ".join(metrics) + " |"
    sep = "|--------|" + "|".join(["--------"] * len(metrics)) + "|"
    rows = [header, sep]

    for method, vals in results.items():
        cells = [method]
        for m in metrics:
            v = vals.get(m, "-")
            if isinstance(v, dict) and "mean" in v:
                cells.append(f"{v['mean']:.4f} ± {v.get('std', 0):.4f}")
            elif isinstance(v, (int, float)):
                cells.append(f"{v:.4f}")
            else:
                cells.append(str(v))
        rows.append("| " + " | ".join(cells) + " |")

    return "\n".join(rows)


def generate_latex_table(
    results: dict[str, dict],
    metrics: list[str] | None = None,
    caption: str = "Performance comparison",
    label: str = "tab:results",
) -> str:
    """Generate a LaTeX table for paper inclusion."""
    if metrics is None:
        all_keys = set()
        for v in results.values():
            all_keys.update(v.keys())
        metrics = sorted(all_keys)

    ncols = len(metrics) + 1
    col_spec = "l" + "c" * len(metrics)

    lines = [
        f"\\begin{{table}}[htbp]",
        f"\\centering",
        f"\\caption{{{caption}}}",
        f"\\label{{{label}}}",
        f"\\begin{{tabular}}{{{col_spec}}}",
        "\\toprule",
        "Method & " + " & ".join(m.replace("_", " ").title() for m in metrics) + " \\\\",
        "\\midrule",
    ]

    for method, vals in results.items():
        cells = [method.replace("_", " ")]
        for m in metrics:
            v = vals.get(m, "-")
            if isinstance(v, dict) and "mean" in v:
                cells.append(f"${v['mean']:.3f} \\pm {v.get('std', 0):.3f}$")
            elif isinstance(v, (int, float)):
                cells.append(f"${v:.3f}$")
            else:
                cells.append(str(v))
        lines.append(" & ".join(cells) + " \\\\")

    lines += [
        "\\bottomrule",
        "\\end{tabular}",
        "\\end{table}",
    ]
    return "\n".join(lines)


# ---------------------------------------------------------------------------
# 4. Key findings extraction
# ---------------------------------------------------------------------------

def extract_key_findings(
    baseline_results: dict[str, dict],
    proposed_results: dict[str, dict],
    metric_key: str = "throughput",
) -> list[str]:
    """Extract key findings comparing baseline vs proposed methods.

    Returns list of finding strings.
    """
    findings = []

    # Best baseline
    best_baseline = None
    best_val = -np.inf
    for name, vals in baseline_results.items():
        v = vals.get(metric_key, {})
        mean = v.get("mean", v) if isinstance(v, dict) else v
        if isinstance(mean, (int, float)) and mean > best_val:
            best_val = mean
            best_baseline = name

    # Best proposed
    best_proposed = None
    best_prop_val = -np.inf
    for name, vals in proposed_results.items():
        v = vals.get(metric_key, {})
        mean = v.get("mean", v) if isinstance(v, dict) else v
        if isinstance(mean, (int, float)) and mean > best_prop_val:
            best_prop_val = mean
            best_proposed = name

    if best_baseline and best_proposed:
        improvement = ((best_prop_val - best_val) / max(abs(best_val), 1e-9)) * 100
        findings.append(
            f"Best proposed method ({best_proposed}) achieves "
            f"{improvement:+.1f}% {metric_key} vs best baseline ({best_baseline})"
        )

    return findings


# ---------------------------------------------------------------------------
# 5. Full analysis pipeline
# ---------------------------------------------------------------------------

def run_full_analysis(
    results_dir: str,
    output_dir: str | None = None,
) -> dict:
    """Run full statistical analysis on experiment results.

    Expects results_dir to contain JSON files with metrics per method.
    """
    results_path = Path(results_dir)
    if output_dir is None:
        output_dir = str(results_path / "analysis")
    Path(output_dir).mkdir(parents=True, exist_ok=True)

    # Load results
    all_results = {}
    for f in results_path.glob("*.json"):
        with open(f) as fh:
            all_results[f.stem] = json.load(fh)

    if not all_results:
        print(f"No JSON result files found in {results_dir}")
        return {}

    # Compute stats for each method
    summary = {}
    for method, data in all_results.items():
        method_stats = {}
        for metric, values in data.items():
            if isinstance(values, list):
                method_stats[metric] = compute_statistics(values)
            elif isinstance(values, (int, float)):
                method_stats[metric] = {"mean": values, "std": 0, "n": 1}
        summary[method] = method_stats

    # Generate tables
    md_table = generate_results_table(summary)
    with open(Path(output_dir) / "results_table.md", "w") as f:
        f.write(md_table)

    latex_table = generate_latex_table(summary)
    with open(Path(output_dir) / "results_table.tex", "w") as f:
        f.write(latex_table)

    # Save summary
    # Convert numpy types for JSON serialization
    def convert(obj):
        if isinstance(obj, (np.integer,)):
            return int(obj)
        if isinstance(obj, (np.floating,)):
            return float(obj)
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        return obj

    with open(Path(output_dir) / "summary.json", "w") as f:
        json.dump(summary, f, indent=2, default=convert)

    print(f"Analysis saved to {output_dir}")
    return summary


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------
if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Statistical analysis")
    parser.add_argument("--results-dir", type=str, default="data/results")
    parser.add_argument("--output-dir", type=str, default=None)
    args = parser.parse_args()

    run_full_analysis(args.results_dir, args.output_dir)
