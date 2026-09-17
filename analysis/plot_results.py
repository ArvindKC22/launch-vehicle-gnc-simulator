"""Post-process Monte Carlo outputs into dispersion and consistency figures."""
import argparse
from pathlib import Path
from statistics import NormalDist

import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
import numpy as np
import pandas as pd


def chi2_quantile(probability, dof):
    """Exact chi-square quantile when SciPy is available, otherwise Wilson-Hilferty."""
    try:
        from scipy.stats import chi2
        return float(chi2.ppf(probability, dof))
    except ImportError:
        z = NormalDist().inv_cdf(probability)
        return dof * max(0.0, 1.0 - 2.0 / (9.0 * dof) + z * np.sqrt(2.0 / (9.0 * dof))) ** 3


def aggregate_consistency(trace, metric):
    rows = trace[trace.metric == metric].copy()
    if rows.empty:
        return rows
    rows["epoch_s"] = rows.time_s.round(6)
    grouped = []
    for epoch, group in rows.groupby("epoch_s", sort=True):
        count = len(group)
        dof = int(group.dof.iloc[0])
        total_dof = count * dof
        is_nees = metric == "nees"
        grouped.append({
            "time_s": epoch,
            "mean": group.value.mean() if is_nees else group.value.sum() / total_dof,
            "lower": chi2_quantile(0.025, total_dof) / (count if is_nees else total_dof),
            "upper": chi2_quantile(0.975, total_dof) / (count if is_nees else total_dof),
            "count": count,
            "dof": dof,
        })
    return pd.DataFrame(grouped)


def plot_consistency(trace, output_path):
    panels = [("nees", "Ensemble NEES", 6.0),
              ("gps_nis", "GPS normalized NIS", 1.0),
              ("baro_nis", "Barometer normalized NIS", 1.0)]
    available = [(key, label, expected) for key, label, expected in panels
                 if (trace.metric == key).any()]
    if not available:
        return
    fig, axes = plt.subplots(len(available), 1, sharex=True,
                             figsize=(8, 2.8 * len(available)), squeeze=False)
    for ax, (key, label, expected) in zip(axes[:, 0], available):
        data = aggregate_consistency(trace, key)
        stride = max(1, int(np.ceil(len(data) / 2000)))
        shown = data.iloc[::stride]
        ax.plot(shown.time_s, shown["mean"], lw=0.9, label=label)
        ax.fill_between(shown.time_s, shown.lower, shown.upper,
                        color="tab:red", alpha=0.18, label="ensemble 95% bounds")
        ax.axhline(expected, color="black", ls="--", lw=1.2, label="expected")
        inside = ((data["mean"] >= data.lower) & (data["mean"] <= data.upper)).mean()
        ax.set_ylabel(label)
        ax.text(0.01, 0.96, f"epochs inside bounds: {inside:.1%}",
                transform=ax.transAxes, va="top")
        ax.legend(loc="best")
    axes[-1, 0].set_xlabel("simulation time (s)")
    fig.tight_layout()
    fig.savefig(output_path, dpi=180)
    plt.close(fig)


def standardized_sensitivity(data, xy):
    columns = {
        "wind_scale": "wind",
        "thrust_scale": "thrust",
        "dry_mass_scale": "mass",
        "aero_scale": "aero",
        "sensor_bias_scale": "sensor bias",
        "thrust_misalignment_deg": "thrust misalignment",
        "cg_offset_m": "CG offset",
    }
    columns = {column: label for column, label in columns.items() if column in data}
    x = data[list(columns)].to_numpy(dtype=float)
    x_std = x.std(axis=0, ddof=1)
    active = x_std > np.finfo(float).eps
    x = (x[:, active] - x[:, active].mean(axis=0)) / x_std[active]
    labels = np.asarray(list(columns.values()))[active]
    design = np.column_stack((np.ones(len(x)), x))
    scores = np.zeros(x.shape[1])
    r_squared = []
    for response in xy.T:
        response_std = response.std(ddof=1)
        if response_std <= np.finfo(float).eps:
            r_squared.append(0.0)
            continue
        y = (response - response.mean()) / response_std
        coefficients = np.linalg.lstsq(design, y, rcond=None)[0]
        with np.errstate(over="ignore", invalid="ignore", divide="ignore"):
            prediction = design @ coefficients
        r_squared.append(max(0.0, 1.0 - np.sum((y - prediction) ** 2) / np.sum(y ** 2)))
        scores += coefficients[1:] ** 2
    if scores.sum() > 0.0:
        scores /= scores.sum()
    return pd.DataFrame({"input": labels, "contribution": scores}), r_squared


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("csv")
    parser.add_argument("--out", default="figures")
    parser.add_argument("--consistency-csv")
    args = parser.parse_args()
    output = Path(args.out)
    output.mkdir(parents=True, exist_ok=True)
    all_runs = pd.read_csv(args.csv)
    data = all_runs[all_runs.impact_valid == 1].dropna(
        subset=["impact_x_m", "impact_y_m", "apogee_m"])
    numeric_columns = ["impact_x_m", "impact_y_m", "apogee_m"]
    data = data[np.isfinite(data[numeric_columns]).all(axis=1)]
    if len(data) < 3:
        raise SystemExit("Need at least 3 valid impact runs to compute an ellipse")

    xy = data[["impact_x_m", "impact_y_m"]].to_numpy()
    mean = xy.mean(axis=0)
    covariance = np.cov(xy, rowvar=False)
    values, vectors = np.linalg.eigh(covariance)
    order = values.argsort()[::-1]
    values, vectors = np.maximum(values[order], 0.0), vectors[:, order]
    angle = np.degrees(np.arctan2(vectors[1, 0], vectors[0, 0]))
    fig, ax = plt.subplots()
    ax.scatter(*xy.T, s=3, alpha=0.15)
    ax.scatter(*mean, c="red", label="Monte Carlo mean")
    ax.add_patch(Ellipse(mean, 6 * np.sqrt(values[0]), 6 * np.sqrt(values[1]),
                         angle=angle, fill=False, color="black", lw=2,
                         label="3-sigma ellipse"))
    ax.set(xlabel="Downrange from launch site (m)",
           ylabel="Crossrange from launch site (m)")
    ax.axis("equal")
    ax.legend()
    fig.tight_layout()
    fig.savefig(output / "landing_ellipse.png", dpi=180)
    plt.close(fig)

    fig, ax = plt.subplots()
    bins = min(40, max(8, int(np.sqrt(len(data)) * 2)))
    ax.hist(data.apogee_m, bins=bins)
    ax.set(xlabel="Apogee (m)", ylabel="runs")
    fig.tight_layout()
    fig.savefig(output / "apogee_distribution.png", dpi=180)
    plt.close(fig)

    if args.consistency_csv:
        plot_consistency(pd.read_csv(args.consistency_csv),
                         output / "nees_nis_consistency.png")

    sensitivity, r_squared = standardized_sensitivity(data, xy)
    sensitivity = sensitivity.sort_values("contribution")
    fig, ax = plt.subplots()
    ax.barh(sensitivity.input, sensitivity.contribution)
    ax.set(xlabel="first-order standardized contribution")
    ax.set_title(f"Linear surrogate: R² downrange={r_squared[0]:.2f}, "
                 f"crossrange={r_squared[1]:.2f}")
    if len(data) < 100:
        ax.text(0.99, 0.02, f"preliminary: only {len(data)} valid runs",
                transform=ax.transAxes, ha="right", color="darkred")
        print("Warning: fewer than 100 valid runs; sensitivity ranking is preliminary")
    fig.tight_layout()
    fig.savefig(output / "tornado.png", dpi=180)
    plt.close(fig)


if __name__ == "__main__":
    main()
