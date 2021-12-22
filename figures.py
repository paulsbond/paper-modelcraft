#!/usr/bin/python3
# -*- coding: utf-8 -*-

import os
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


_COLOURS = ["#377eb8", "#ff7f00"]
_MARKERS = ["o", "v"]


plt.rc("axes", linewidth=0.6)
plt.rc("font", size=8, family="sans-serif")
plt.rc("legend", fontsize=8, numpoints=1)


def _main():
    results = pd.read_csv("results.csv")
    results.dropna(inplace=True)
    results["ccp4i"] *= 100
    results["modelcraft"] *= 100
    os.makedirs("figures", exist_ok=True)
    _comparison("1", results)
    _binned("2", results, "resolution", "Resolution / Å", "lower left")
    _binned("3", results, "f_map_correlation", "F-map Correlation", "lower right")


def _comparison(number, results):
    af = results[results["type"] == "af"]
    ep = results[results["type"] == "ep"]
    x1 = af["ccp4i"]
    x2 = ep["ccp4i"]
    y1 = af["modelcraft"]
    y2 = ep["modelcraft"]
    min_ = -3
    max_ = 103
    ax1, ax2 = _axes()
    ax1.plot([min_, max_], [min_, max_], "k-", alpha=0.5, linewidth=0.5)
    ax2.plot([min_, max_], [min_, max_], "k-", alpha=0.5, linewidth=0.5)
    ax1.plot(x1, y1, "kx", markersize=4, color=_COLOURS[0])
    ax2.plot(x2, y2, "kx", markersize=4, color=_COLOURS[0])
    ax1.axis([min_, max_, min_, max_])
    ax2.axis([min_, max_, min_, max_])
    ax1.set_aspect("equal", "box")
    ax2.set_aspect("equal", "box")
    _save_fig(number, ax1, ax2, "CCP4i Completeness / %", "ModelCraft Completeness / %")


def _binned(number, results, xkey, xlabel, legend_loc):
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


def _axes():
    fig = plt.figure(figsize=(12 / 2.54, 7 / 2.54), dpi=300)
    ax1 = fig.add_subplot(121)
    ax2 = fig.add_subplot(122, sharex=ax1, sharey=ax1)
    ax1.set_title("Molecular Replacement")
    ax2.set_title("Experimental Phasing")
    return ax1, ax2


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


def _save_fig(number, ax1, ax2, xlabel, ylabel):
    ax1.tick_params(direction="out", length=3, pad=3, top=False, right=False)
    ax2.tick_params(direction="out", length=3, pad=3, top=False, right=False)
    plt.setp(ax2.get_yticklabels(), visible=False)
    ax1.set_xlabel(xlabel)
    ax2.set_xlabel(xlabel)
    ax1.set_ylabel(ylabel)
    plt.tight_layout(pad=0.3)
    plt.savefig(f"figures/figure{number}.png")
    plt.close()


if __name__ == "__main__":
    _main()
