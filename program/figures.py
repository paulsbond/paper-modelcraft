#!/usr/bin/python3
# -*- coding: utf-8 -*-

import os
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# Line weights should be between 0.35 and 1.5pt at final size.
# Monochrome and greyscale images should have a minimum resolution of 600 d.p.i.
# Single-column width (8.85 cm), part-page width (12 cm) or full-page width (18 cm).
# A page is approximately 18 x 24 cm
# At the final published size, the labelling on the figure should be approximately 8pt.
# Use Arial, Courier, Helvetica, Symbol, Times or Times New Roman

_BLUE = "#377eb8"
_ORANGE = "#ff7f00"

plt.rc("axes", titlesize=8, labelsize=8, linewidth=0.6)
plt.rc("font", size=8, family="sans-serif")
plt.rc("legend", fontsize=8, numpoints=1)
plt.rc("xtick", labelsize=8)
plt.rc("ytick", labelsize=8)


def make_figures():
    print("Making figures...")
    results_mr = pd.read_csv("results/results_mr.csv")
    results_ep = pd.read_csv("results/results_ep.csv")
    results_af = pd.read_csv("results/results_af.csv")
    os.makedirs("figures", exist_ok=True)
    _completeness(results_mr, "mr")
    _completeness(results_ep, "ep")
    _completeness(results_af, "af")
    # mr_res
    # ep_res
    # mr_fmap
    # ep_fmap
    # mr_time
    # mr_ablation


def _completeness(results, type_):
    fig = plt.figure(figsize=(8.85 / 2.54, 8.85 / 2.54), dpi=600)
    ax = fig.add_subplot(111)
    sub = results[results["type"] == type_]
    sub = sub[["ccp4i_completeness", "modelcraft_completeness"]]
    sub = sub.dropna()
    x = sub["ccp4i_completeness"]
    y = sub["modelcraft_completeness"]
    min_, max_ = (0, 1)
    ax.plot([min_, max_], [min_, max_], "k-", alpha=0.5, linewidth=0.5)
    ax.plot(x, y, "kx", markersize=4, color=_BLUE)
    ax.axis([min_, max_, min_, max_])
    ax.set_aspect("equal", "box")
    ax.tick_params(direction="out", length=3, pad=3, top=False, right=False)
    ax.set_xlabel("CCP4i Completeness")
    ax.set_ylabel("ModelCraft Completeness")
    plt.tight_layout(pad=0.3)
    plt.savefig(f"figures/fig_{type_}.png")
    with open(f"figures/fig_{type_}.txt", "w") as stream:
        print("Points:", len(x), file=stream)
    plt.close()


def _binned(results, type_):
    fig = plt.figure(figsize=(8.85 / 2.54, 8.85 / 2.54), dpi=600)
    ax = fig.add_subplot(111)
    x, y = _get_xy(results, type_, "ccp4i_completeness", "modelcraft_completeness")
    af = results[results["type"] == "af"]
    ep = results[results["type"] == "ep"]
    ax1, ax2 = _axes()
    bin_centres = []
    for ax, subset in ((ax1, af), (ax2, ep)):
        data = {"ModelCraft": {"x": [], "y": []}, "CCP4i": {"x": [], "y": []}}
        for row in subset.to_records():
            data["ModelCraft"]["x"].append(row[xkey])
            data["ModelCraft"]["y"].append(row["modelcraft"])
            data["CCP4i"]["x"].append(row[xkey])
            data["CCP4i"]["y"].append(row["ccp4i"])
        for i, key in enumerate(data):
            x = data[key]["x"]
            y = data[key]["y"]
            mean, _, se, bin_center = _bin_xy(x, y)
            bin_centres.extend(bin_center)
            ax.plot(
                bin_center, mean, label=key, color=_COLOURS[i], marker=_MARKERS[i],
            )
            ax.fill_between(
                bin_center,
                mean - se,
                mean + se,
                alpha=0.5,
                color=_COLOURS[i],
                linewidth=0.0,
            )
    min_x, max_x = _min_max(bin_centres, pad=0.08)
    ax1.axis([min_x, max_x, 0, 100])
    ax2.axis([min_x, max_x, 0, 100])
    ax1.legend(loc=legend_loc)
    ax2.legend(loc=legend_loc)
    _save_fig(number, ax1, ax2, xlabel, "Completeness / %")


def _bin_xy(x, y, nbins=5):
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


def _min_max(x, pad=0.02):
    min_x = min(x)
    max_x = max(x)
    padding = (max_x - min_x) * pad
    min_x -= padding
    max_x += padding
    return min_x, max_x
