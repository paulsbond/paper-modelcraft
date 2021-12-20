#!/usr/bin/python3

import os
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


_COLOURS = [
    "#377eb8",
    "#ff7f00",
    "#4daf4a",
    "#f781bf",
    "#a65628",
    "#984ea3",
    "#999999",
    "#e41a1c",
    "#dede00",
]
_MARKERS = ["o", "v", "s", "p", "P", "*", "h", "X", "d"]


plt.rc("axes", linewidth=0.6)
plt.rc("font", size=8, family="sans-serif")
plt.rc("legend", fontsize=8, numpoints=1)


def _main():
    results = pd.read_csv("results.csv")
    results.dropna(inplace=True)
    os.makedirs("figures", exist_ok=True)
    af = results[results["type"] == "af"]
    ep = results[results["type"] == "ep"]
    _comparison("figures/01_af_raw.png", af, "ccp4i", "modelcraft")
    _comparison("figures/02_ep_raw.png", ep, "ccp4i", "modelcraft")
    _binned("figures/03_af_res.png", af, "resolution", (0, 4))
    _binned("figures/04_ep_res.png", ep, "resolution", (0, 4))
    _binned("figures/05_af_fmap.png", af, "f_map_correlation", (0, 1))
    _binned("figures/06_eo_fmap.png", ep, "f_map_correlation", (0, 1))
    _scatter("figures/07_af_ccp4i_res.png", af, "resolution", "ccp4i", (0, 4))
    _scatter("figures/08_ep_ccp4i_res.png", ep, "resolution", "ccp4i", (0, 4))
    _scatter("figures/09_af_ccp4i_fmap.png", af, "f_map_correlation", "ccp4i", (0, 1))
    _scatter("figures/10_ep_ccp4i_fmap.png", ep, "f_map_correlation", "ccp4i", (0, 1))
    _scatter("figures/11_af_modelcraft_res.png", af, "resolution", "modelcraft", (0, 4))
    _scatter("figures/12_ep_modelcraft_res.png", ep, "resolution", "modelcraft", (0, 4))
    _scatter("figures/13_af_modelcraft_fmap.png", af, "f_map_correlation", "modelcraft", (0, 1))
    _scatter("figures/14_ep_modelcraft_fmap.png", ep, "f_map_correlation", "modelcraft", (0, 1))


def _scatter(path, results, xkey, ykey, xlim):
    ax = _single_ax()
    x = results[xkey]
    y = results[ykey]
    ax.plot(x, y, "kx", markersize=3, color=_COLOURS[0], alpha=0.8)
    ax.set_xlim(*xlim)
    ax.set_ylim(0, 1)
    _save_fig(ax, path, xkey, ykey)


def _comparison(path, results, xkey, ykey):
    x = results[xkey]
    y = results[ykey]
    ax = _single_ax(8, 8)
    ax.plot([0, 1], [0, 1], "k--", alpha=0.5, linewidth=1)
    ax.plot(x, y, "kx", markersize=3, color=_COLOURS[0], alpha=0.8)
    ax.axis([0, 1, 0, 1])
    ax.set_aspect("equal", "box")
    _save_fig(ax, path, xkey, ykey)


def _binned(
    path,
    results,
    xkey,
    xlim,
    nbins=5,
    legend_loc="best",
):
    data = {"ModelCraft": {"x": [], "y": []}, "CCP4i": {"x": [], "y": []}}
    for row in results.to_records():
        data["ModelCraft"]["x"].append(row[xkey])
        data["CCP4i"]["x"].append(row[xkey])
        data["ModelCraft"]["y"].append(row["modelcraft"])
        data["CCP4i"]["y"].append(row["ccp4i"])
    ax = _single_ax()
    all_x = []
    for i, key in enumerate(data):
        x = data[key]["x"]
        y = data[key]["y"]
        all_x.extend(x)
        mean, _, se, bin_center = _bin_xy(x, y, nbins)
        ax.plot(
            bin_center,
            mean,
            label=key,
            color=_COLOURS[i],
            marker=_MARKERS[i],
        )
        ax.fill_between(
            bin_center,
            mean - se,
            mean + se,
            alpha=0.5,
            color=_COLOURS[i],
            linewidth=0.0,
        )
    min_x, max_x = _min_max(all_x)
    ax.set_xlim(*xlim)
    ax.set_ylim(0, 1)
    ax.legend(loc=legend_loc)
    _save_fig(ax, path, xkey, "Completeness")


def _single_ax(width_mm=12, height_mm=8):
    fig = plt.figure(figsize=(width_mm / 2.54, height_mm / 2.54), dpi=300)
    return fig.add_subplot(111)


def _save_fig(ax, path, xlabel, ylabel):
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.tick_params(direction="out", length=3, pad=3, top=False, right=False)
    plt.tight_layout(pad=0.3)
    plt.savefig(path)
    plt.close()


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
    min_x = min_x - (max_x - min_x) * pad
    max_x = max_x + (max_x - min_x) * pad
    return min_x, max_x


if __name__ == "__main__":
    _main()
