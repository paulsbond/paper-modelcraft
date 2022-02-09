#!/usr/bin/python3
# -*- coding: utf-8 -*-

# Line weights should be between 0.35 and 1.5pt at final size.
# Monochrome and greyscale images should have a minimum resolution of 600 d.p.i.
# Single-column width (8.85 cm), part-page width (12 cm) or full-page width (18 cm).
# A page is approximately 18 x 24 cm
# At the final published size, the labelling on the figure should be approximately 8pt.
# Use Arial, Courier, Helvetica, Symbol, Times or Times New Roman

import os
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


_COLOURS = ["#377eb8", "#ff7f00"]
_MARKERS = ["o", "v"]

plt.rc("axes", titlesize=8, labelsize=8, linewidth=0.6)
plt.rc("font", size=8, family="sans-serif")
plt.rc("legend", fontsize=8, numpoints=1)
plt.rc("xtick", labelsize=8)
plt.rc("ytick", labelsize=8)


def make_figures():
    print("Making figures...")
    os.makedirs("figures", exist_ok=True)
    results_mr = pd.read_csv("results/results_mr.csv")
    results_ep = pd.read_csv("results/results_ep.csv")
    results_af = pd.read_csv("results/results_af.csv")
    results_mr.dropna(inplace=True)
    results_ep.dropna(inplace=True)
    results_af.dropna(inplace=True)
    with open("figures/samples.txt", "w") as stream:
        print("MR:", len(results_mr), file=stream)
        print("EP:", len(results_ep), file=stream)
        print("AF:", len(results_af), file=stream)
    _completeness("mr", results_mr)
    _completeness("ep", results_ep)
    _completeness("af", results_af)
    _binned("mr_res", results_mr, 1, 3.5, "resolution", "Resolution / Å")
    _binned("ep_res", results_ep, 1, 3.5, "resolution", "Resolution / Å")
    _binned("mr_fmap", results_mr, 0.2, 1, "f_map_correlation", "F-map Correlation")
    _binned("ep_fmap", results_ep, 0.2, 1, "f_map_correlation", "F-map Correlation")
    # mr_time
    # mr_ablation


def _completeness(name, results):
    fig = plt.figure(figsize=(8.85 / 2.54, 8.85 / 2.54), dpi=600)
    ax = fig.add_subplot(111)
    x = results["ccp4i_completeness"]
    y = results["modelcraft_completeness"]
    min_, max_ = (0, 1)
    ax.plot([min_, max_], [min_, max_], "k-", alpha=0.5, linewidth=0.5)
    ax.plot(x, y, "kx", markersize=4, color=_COLOURS[0])
    ax.axis([min_, max_, min_, max_])
    ax.set_aspect("equal", "box")
    ax.tick_params(direction="out", length=3, pad=3, top=False, right=False)
    ax.set_xlabel("CCP4i Completeness")
    ax.set_ylabel("ModelCraft Completeness")
    plt.tight_layout(pad=0.3)
    plt.savefig(f"figures/fig_{name}.png")
    plt.close()


def _binned(name, results, xmin, xmax, xkey, xlabel):
    fig = plt.figure(figsize=(8.85 / 2.54, 8.85 / 2.54), dpi=600)
    ax = fig.add_subplot(111)
    data = {"ModelCraft": {"x": [], "y": []}, "CCP4i": {"x": [], "y": []}}
    for row in results.to_records():
        data["ModelCraft"]["x"].append(row[xkey])
        data["ModelCraft"]["y"].append(row["modelcraft_completeness"])
        data["CCP4i"]["x"].append(row[xkey])
        data["CCP4i"]["y"].append(row["ccp4i_completeness"])
    for i, key in enumerate(data):
        x = data[key]["x"]
        y = data[key]["y"]
        mean, _, se, bin_center = _bin_xy(x, y)
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
    ax.axis([xmin, xmax, 0, 1])
    ax.legend()
    ax.set_xlabel(xlabel)
    ax.set_ylabel("Completeness")
    plt.tight_layout(pad=0.3)
    plt.savefig(f"figures/fig_{name}.png")
    plt.close()


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

