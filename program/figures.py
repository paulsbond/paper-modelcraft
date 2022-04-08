#!/usr/bin/python3
# -*- coding: utf-8 -*-

# Line weights should be between 0.35 and 1.5pt at final size.
# Monochrome and greyscale images should have a minimum resolution of 600 d.p.i.
# Single-column width (8.85 cm), part-page width (12 cm) or full-page width (18 cm).
# A page is approximately 18 x 24 cm
# At the final published size, the labelling on the figure should be approximately 8pt.
# Use Arial, Courier, Helvetica, Symbol, Times or Times New Roman

import datetime
import os
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.stats
from matplotlib.ticker import ScalarFormatter


_COLOURS = ["#377eb8", "#ff7f00", "#4daf4a"]
_DARKER = ["#2c6493", "#cc6500", "#3d8c3b"]
_MARKERS = ["o", "v", "s"]
_HATCHES = [None, "///", "..."]

plt.rc("axes", titlesize=8, labelsize=8, linewidth=0.6)
plt.rc("font", size=8, family="Arial")
plt.rc("legend", fontsize=8, numpoints=1)
plt.rc("xtick", labelsize=8)
plt.rc("ytick", labelsize=8)


def _make_figures():
    print("Making figures...")
    os.makedirs("figures", exist_ok=True)
    results_mr = pd.read_csv("results/results_mr.csv")
    results_ep = pd.read_csv("results/results_ep.csv")
    results_af = pd.read_csv("results/results_af.csv")
    for results in results_mr, results_ep, results_af:
        results.dropna(inplace=True)
        _add_extra_columns(results)
    _mrep(results_mr, results_ep)
    _time(results_mr)
    _ablation(results_mr)
    _af(results_af)
    _mr_stats(results_mr)
    _ep_stats(results_ep)
    _af_stats(results_af)


def _add_extra_columns(results):
    results["extra_completeness"] = (
        results["modelcraft_completeness"] - results["ccp4i_completeness"]
    )
    results["extra_seconds"] = results["modelcraft_seconds"] - results["ccp4i_seconds"]


def _mrep(results_mr, results_ep):
    _, ((ax1, ax2), (ax3, ax4), (ax5, ax6)) = plt.subplots(
        nrows=3,
        ncols=2,
        gridspec_kw={"height_ratios": [1.8, 1, 1]},
        figsize=(18 / 2.54, 20 / 2.54),
        dpi=600,
    )
    for ax, results, title in (
        (ax1, results_mr, "Molecular Replacement"),
        (ax2, results_ep, "Experimental Phasing"),
    ):
        _raw_completness(ax, results)
        ax.set_title(title)
    for ax, results, xkey, xlabel, xmin, xmax in (
        (ax3, results_mr, "resolution", "Resolution / Å", 3.5, 1.0),
        (ax4, results_ep, "resolution", "Resolution / Å", 3.5, 1.0),
        (ax5, results_mr, "f_map_correlation", "F-map Correlation", 0.2, 1.0),
        (ax6, results_ep, "f_map_correlation", "F-map Correlation", 0.2, 1.0),
    ):
        _binned_completeness(ax, results, xkey, xlabel, xmin, xmax)
    plt.tight_layout(pad=0.3)
    plt.savefig("figures/fig_mrep.png")
    plt.close()


def _time(results):
    fig = plt.figure(figsize=(8.85 / 2.54, 8.85 / 2.54), dpi=600)
    ax = fig.add_subplot(111)
    x = results["extra_seconds"] / (60 * 60)  # Convert seconds to hours
    y = results["extra_completeness"] * 100
    ax.plot(x, y, "kx", markersize=4, color=_COLOURS[0])
    ax.set_xscale("log")
    ax.xaxis.set_major_formatter(ScalarFormatter())
    ax.tick_params(direction="out", length=3, pad=3, top=False, right=False)
    ax.set_xlabel("Extra Time / h")
    ax.set_ylabel("Extra Completeness / p.p.")
    plt.tight_layout(pad=0.3)
    plt.savefig("figures/fig_time.png")
    plt.close()


def _ablation(results):
    labels = [
        "Sheetbend",
        "Pruning",
        "Parrot",
        "Dummy Atom",
        "Water",
        "Side Chain",
    ]
    keys = [
        "modelcraft_no_sheetbend",
        "modelcraft_no_pruning",
        "modelcraft_no_parrot",
        "modelcraft_no_dummy_atoms",
        "modelcraft_no_waters",
        "modelcraft_no_side_chain_fixing",
    ]
    metrics = ["completeness", "rwork", "rfree"]
    data = {label: {metric: [] for metric in metrics} for label in labels}
    for metric in metrics:
        for result in results.to_records():
            for label, key in zip(labels, keys):
                value = result[f"{key}_{metric}"] - result[f"modelcraft_{metric}"]
                data[label][metric].append(value * 100)
        for label in labels:
            data[label][metric + "_mean"] = np.mean(data[label][metric])
            data[label][metric + "_sem"] = scipy.stats.sem(data[label][metric])
    plot_data = {
        "Completeness": {
            "heights": [data[label]["completeness_mean"] for label in labels],
            "yerr": [data[label]["completeness_sem"] for label in labels],
        },
        "R-work": {
            "heights": [data[label]["rwork_mean"] for label in labels],
            "yerr": [data[label]["rwork_sem"] for label in labels],
        },
        "R-free": {
            "heights": [data[label]["rfree_mean"] for label in labels],
            "yerr": [data[label]["rfree_sem"] for label in labels],
        },
    }
    x = np.arange(len(labels))
    fig = plt.figure(figsize=(8.85 / 2.54, 8.85 / 2.54), dpi=600)
    ax = fig.add_subplot(111)
    width = 0.3
    for i, group in enumerate(plot_data):
        heights = plot_data[group]["heights"]
        yerr = plot_data[group]["yerr"]
        adjusted = x - width + i * width
        ax.bar(
            adjusted,
            heights,
            width,
            yerr=yerr,
            capsize=0,
            label=group,
            color=_COLOURS[i],
            hatch=_HATCHES[i],
            edgecolor=_DARKER[i],
        )
    plt.xticks(rotation=40)
    ax.legend(loc="lower right")
    ax.set_xticks(x)
    ax.set_xticklabels(labels)
    ax.margins(x=0.01)
    ax.tick_params(direction="out", length=3, pad=3, top=False, right=False)
    ax.set_xlabel("Step Removed")
    ax.set_ylabel("Mean Change / p.p.")
    plt.tight_layout(pad=0.3)
    plt.savefig("figures/fig_ablation.png")
    plt.close()


def _af(results):
    fig = plt.figure(figsize=(8.85 / 2.54, 8.85 / 2.54), dpi=600)
    ax = fig.add_subplot(111)
    _raw_completness(ax, results)
    plt.tight_layout(pad=0.3)
    plt.savefig("figures/fig_af.png")
    plt.close()


def _raw_completness(ax, results):
    x = results["ccp4i_completeness"] * 100
    y = results["modelcraft_completeness"] * 100
    min_, max_ = (0, 100)
    ax.plot([min_, max_], [min_, max_], "k--", alpha=0.5, linewidth=0.8)
    ax.plot(x, y, "kx", markersize=4, color=_COLOURS[0])
    ax.axis([min_, max_, min_, max_])
    ax.set_aspect("equal", "box")
    ax.tick_params(direction="out", length=3, pad=3, top=False, right=False)
    ax.set_xlabel("CCP4i Buccaneer Completeness / %")
    ax.set_ylabel("ModelCraft Completeness / %")


def _binned_completeness(ax, results, xkey, xlabel, xmin, xmax):
    data = {"ModelCraft": {"x": [], "y": []}, "CCP4i Buccaneer": {"x": [], "y": []}}
    for row in results.to_records():
        data["ModelCraft"]["x"].append(row[xkey])
        data["ModelCraft"]["y"].append(row["modelcraft_completeness"] * 100)
        data["CCP4i Buccaneer"]["x"].append(row[xkey])
        data["CCP4i Buccaneer"]["y"].append(row["ccp4i_completeness"] * 100)
    for i, key in enumerate(data):
        x = data[key]["x"]
        y = data[key]["y"]
        mean, _, se, bin_center = _bin_xy(x, y)
        ax.plot(bin_center, mean, label=key, color=_COLOURS[i], marker=_MARKERS[i])
        ax.fill_between(
            bin_center,
            mean - se,
            mean + se,
            alpha=0.5,
            color=_COLOURS[i],
            linewidth=0.0,
        )
    ax.axis([xmin, xmax, 0, 100])
    ax.legend(loc="lower right")
    ax.set_xlabel(xlabel)
    ax.set_ylabel("Completeness / %")


def _bin_xy(x, y, nbins=3):
    x = np.array(x)
    y = np.array(y)
    n, bin_edges = np.histogram(x, bins=nbins)
    sy, bin_edges = np.histogram(x, bins=nbins, weights=y)
    sy2, bin_edges = np.histogram(x, bins=nbins, weights=y * y)
    bin_center = (bin_edges[1:] + bin_edges[:-1]) / 2
    mean = sy / n
    std = np.sqrt(sy2 / (n - 1) - mean * mean)
    se = std / np.sqrt(n)
    return mean, std, se, bin_center


def _mr_stats(results):
    print("\n## MR stats")
    print(len(results), "samples")
    stat = sum(results["extra_completeness"] > 0) / len(results) * 100
    print(f"ModelCraft produced a more complete model in {stat}% of cases")
    stat = sum(results["modelcraft_completeness"] > 0.8) / len(results) * 100
    print(f"ModelCraft produced a >80% complete model in {stat}% of cases")
    stat = sum(results["ccp4i_completeness"] > 0.8) / len(results) * 100
    print(f"CCP4i produced a >80% complete model in {stat}% of cases")
    stat = (
        sum(
            (results["modelcraft_completeness"] < 0.2)
            & (results["ccp4i_completeness"] < 0.2)
        )
        / len(results)
        * 100
    )
    print(f"Both pipelines produced a <20% complete model in {stat}% of cases")
    print("Mean extra time:", _timestr(np.mean(results["extra_seconds"])))
    print("Mean CCP4i time:", _timestr(np.mean(results["ccp4i_seconds"])))
    print("Mean ModelCraft time:", _timestr(np.mean(results["modelcraft_seconds"])))
    print("Median extra time:", _timestr(np.median(results["extra_seconds"])))
    print("Mean extra completeness:", np.mean(results["extra_completeness"]))
    print("Mean CCP4i completeness:", np.mean(results["ccp4i_completeness"]))
    print("Mean ModelCraft completeness:", np.mean(results["modelcraft_completeness"]))
    print("Median extra completeness:", np.median(results["extra_completeness"]))


def _ep_stats(results):
    print("\n## EP stats")
    print(len(results), "samples")
    stat = (
        sum(
            (results["modelcraft_completeness"] > 0.8)
            & (results["ccp4i_completeness"] > 0.8)
        )
        / len(results)
        * 100
    )
    print(f"Both pipelines produced a >80% complete model in {stat}% of cases")
    stat = sum(results["extra_completeness"] > 0) / len(results) * 100
    print(f"ModelCraft produced a more complete model in {stat}% of cases")
    print("Mean CCP4i completeness:", np.mean(results["ccp4i_completeness"]))
    print("Mean ModelCraft completeness:", np.mean(results["modelcraft_completeness"]))


def _af_stats(results):
    print("\n## AF stats")
    print(len(results), "samples")
    stat = sum(results["extra_completeness"] > 0) / len(results) * 100
    print(f"ModelCraft produced a more complete model in {stat}% of cases")
    print("Mean CCP4i completeness:", np.mean(results["ccp4i_completeness"]))
    print("Mean ModelCraft completeness:", np.mean(results["modelcraft_completeness"]))
    rowidx = results["extra_completeness"].idxmax()
    print(results["id"][rowidx], "has the max extra completeness")
    rowidx = results["extra_completeness"].idxmin()
    print(results["id"][rowidx], "has the min extra completeness")


def _timestr(seconds):
    return str(datetime.timedelta(seconds=seconds))


if __name__ == "__main__":
    _make_figures()
