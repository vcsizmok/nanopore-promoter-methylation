#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Supplementary functional-characterization figure for promoter
methylation hotspots.

Panels
------
TAD organization
CTCF association
CpG / sequence features
Histone-mark enrichment
18-state ChromHMM enrichment/depletion
Cross-epigenome ChromHMM summary

All panels are also saved individually as PNG and PDF.
"""

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


# ============================================================
# BASE DIRECTORY
# ============================================================

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
BASE = os.path.join(SCRIPT_DIR, "results")


# ============================================================
# INPUT FILES
# ============================================================

TAD_SUMMARY_FILE = os.path.join(
    BASE,
    "TAD_overlap_analysis",
    "ALL_TAD_MAPS_permutation_summary.tsv"
)

CTCF_SUMMARY_FILE = os.path.join(
    BASE,
    "CTCF_analysis",
    "all_contexts",
    "ALL_CONTEXTS_CTCF_permutation_summary.tsv"
)

CPG_SUMMARY_FILE = os.path.join(
    BASE,
    "CpG_sequence_analysis",
    "CpG_sequence_permutation_summary.tsv"
)

HISTONE_SUMMARY_FILE = os.path.join(
    BASE,
    "functional_feature_analysis",
    "all_tissues_histone_comparison",
    "ALL_TISSUES_histone_permutation_summary.tsv"
)

CHROMHMM_SUMMARY_FILE = os.path.join(
    BASE,
    "ChromHMM_analysis",
    "ALL_EPIGENOMES_ChromHMM_permutation_summary.tsv"
)

CHROMHMM_CROSS_FILE = os.path.join(
    BASE,
    "ChromHMM_analysis",
    "cross_epigenome_summary",
    "ChromHMM_cross_epigenome_state_summary.tsv"
)


# ============================================================
# OUTPUT DIRECTORY
# ============================================================

OUTDIR = os.path.join(
    BASE,
    "reviewer_functional_feature_figures"
)

# Input summaries can be supplied without copying reference datasets into GitHub.
import argparse
_parser = argparse.ArgumentParser(description="Recreate final hotspot functional-characterization figure")
for _flag, _default in [("tad", TAD_SUMMARY_FILE), ("ctcf", CTCF_SUMMARY_FILE),
                        ("cpg", CPG_SUMMARY_FILE), ("histone", HISTONE_SUMMARY_FILE),
                        ("chromhmm", CHROMHMM_SUMMARY_FILE),
                        ("chromhmm-cross", CHROMHMM_CROSS_FILE)]:
    _parser.add_argument("--" + _flag, default=_default)
_parser.add_argument("--outdir", default=OUTDIR)
_args = _parser.parse_args()
TAD_SUMMARY_FILE = _args.tad
CTCF_SUMMARY_FILE = _args.ctcf
CPG_SUMMARY_FILE = _args.cpg
HISTONE_SUMMARY_FILE = _args.histone
CHROMHMM_SUMMARY_FILE = _args.chromhmm
CHROMHMM_CROSS_FILE = _args.chromhmm_cross
OUTDIR = _args.outdir

os.makedirs(
    OUTDIR,
    exist_ok=True
)


# ============================================================
# GLOBAL FONT SETTINGS
# ============================================================

TITLE_FONT = 15
AXIS_LABEL_FONT = 13
TICK_FONT = 11
SMALL_TICK_FONT = 10
LEGEND_FONT = 10
COLORBAR_FONT = 11
STAR_FONT = 13


plt.rcParams.update({
    "font.size": 11,
    "axes.titlesize": TITLE_FONT,
    "axes.labelsize": AXIS_LABEL_FONT,
    "xtick.labelsize": TICK_FONT,
    "ytick.labelsize": TICK_FONT,
    "legend.fontsize": LEGEND_FONT
})


# ============================================================
# DISPLAY ORDERS
# ============================================================

TAD_MAP_ORDER = [
    "GM12878",
    "mammary",
    "multi_cancer",
    "gonadal",
    "thyroid_prostate",
    "pancreas",
    "brain"
]

TAD_MAP_LABELS = {
    "GM12878": "GM12878",
    "mammary": "Mammary",
    "multi_cancer": "Multi-cancer",
    "gonadal": "Gonadal",
    "thyroid_prostate": "Thyroid/\nprostate",
    "pancreas": "Pancreas",
    "brain": "Brain"
}


CTCF_ORDER = [
    "MCF7",
    "Mammary_fibroblast",
    "Colon",
    "Liver",
    "Pancreas",
    "Lung",
    "Astrocytes"
]

CTCF_LABELS = {
    "MCF7": "MCF-7",
    "Mammary_fibroblast": "Mammary\nfibroblast",
    "Colon": "Colon",
    "Liver": "Liver",
    "Pancreas": "Pancreas",
    "Lung": "Lung",
    "Astrocytes": "Astrocytes"
}


HISTONE_ORDER = [
    "MCF7",
    "HMEC",
    "vHMEC",
    "Colon",
    "Liver",
    "Pancreas",
    "Lung",
    "Astrocytes"
]

HISTONE_LABELS = {
    "MCF7": "MCF-7",
    "HMEC": "HMEC",
    "vHMEC": "vHMEC",
    "Colon": "Colon",
    "Liver": "Liver",
    "Pancreas": "Pancreas",
    "Lung": "Lung",
    "Astrocytes": "Astrocytes"
}


EPIGENOME_ORDER = [
    "E119_HMEC",
    "E075_Colon",
    "E079_Esophagus",
    "E066_Liver",
    "E098_Pancreas",
    "E097_Ovary",
    "E096_Lung",
    "E125_Astrocytes",
    "E127_Keratinocytes",
    "E129_Osteoblasts"
]

EPIGENOME_LABELS = {
    "E119_HMEC": "HMEC",
    "E075_Colon": "Colon",
    "E079_Esophagus": "Esophagus",
    "E066_Liver": "Liver",
    "E098_Pancreas": "Pancreas",
    "E097_Ovary": "Ovary",
    "E096_Lung": "Lung",
    "E125_Astrocytes": "Astrocytes",
    "E127_Keratinocytes": "Keratinocytes",
    "E129_Osteoblasts": "Osteoblasts"
}


STATE_ORDER = [
    "1_TssA",
    "2_TssFlnk",
    "3_TssFlnkU",
    "4_TssFlnkD",
    "5_Tx",
    "6_TxWk",
    "7_EnhG1",
    "8_EnhG2",
    "9_EnhA1",
    "10_EnhA2",
    "11_EnhWk",
    "12_ZNF/Rpts",
    "13_Het",
    "14_TssBiv",
    "15_EnhBiv",
    "16_ReprPC",
    "17_ReprPCWk",
    "18_Quies"
]


# ============================================================
# HELPERS
# ============================================================

def check_file(filename):

    if not os.path.exists(filename):
        raise IOError(
            "Required file not found:\n{}".format(
                filename
            )
        )


def clean_axes(ax, remove_left=False):

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    if remove_left:
        ax.spines["left"].set_visible(False)

    ax.yaxis.set_ticks_position("left")
    ax.yaxis.set_label_position("left")

    ax.tick_params(
        axis="both",
        which="both",
        labelsize=TICK_FONT
    )

    ax.tick_params(
        axis="y",
        which="both",
        left=not remove_left,
        right=False,
        labelleft=True,
        labelright=False
    )


def symmetric_limit(values, minimum=2.0):

    values = np.asarray(
        values,
        dtype=float
    )

    values = values[
        np.isfinite(values)
    ]

    if len(values) == 0:
        return minimum

    limit = np.percentile(
        np.abs(values),
        97.5
    )

    return max(
        minimum,
        float(limit)
    )


def save_figure(fig, stem):

    fig.savefig(
        os.path.join(
            OUTDIR,
            stem + ".png"
        ),
        dpi=300,
        bbox_inches="tight"
    )

    fig.savefig(
        os.path.join(
            OUTDIR,
            stem + ".pdf"
        ),
        bbox_inches="tight"
    )

    plt.close(fig)


# ============================================================
# LOAD DATA
# ============================================================

def load_tad():

    check_file(TAD_SUMMARY_FILE)

    return pd.read_csv(
        TAD_SUMMARY_FILE,
        sep="\t"
    )


def load_ctcf():

    check_file(CTCF_SUMMARY_FILE)

    return pd.read_csv(
        CTCF_SUMMARY_FILE,
        sep="\t"
    )


def load_cpg():

    check_file(CPG_SUMMARY_FILE)

    return pd.read_csv(
        CPG_SUMMARY_FILE,
        sep="\t"
    )


def load_histone():

    check_file(HISTONE_SUMMARY_FILE)

    df = pd.read_csv(
        HISTONE_SUMMARY_FILE,
        sep="\t"
    )

    metric_labels = {
        "median_hotspot_H3K27me3":
            "H3K27me3",

        "median_hotspot_H3K4me3":
            "H3K4me3",

        "median_hotspot_H3K27me3_H3K4me3_log2ratio":
            "H3K27me3:H3K4me3"
    }

    df["metric_label"] = (
        df["metric"]
        .map(metric_labels)
    )

    return df


def load_chromhmm():

    check_file(CHROMHMM_SUMMARY_FILE)

    return pd.read_csv(
        CHROMHMM_SUMMARY_FILE,
        sep="\t"
    )


def load_cross_chromhmm():

    check_file(CHROMHMM_CROSS_FILE)

    return pd.read_csv(
        CHROMHMM_CROSS_FILE,
        sep="\t"
    )


# ============================================================
# TAD
# ============================================================

def plot_tad(ax, tad_df):

    wanted_tests = [
        "any_TAD_overlap",
        "complete_TAD_containment"
    ]

    sub = tad_df[
        tad_df["test"].isin(
            wanted_tests
        )
    ].copy()

    x = np.arange(
        len(TAD_MAP_ORDER)
    )

    width = 0.18

    overlap = (
        sub[
            sub["test"]
            ==
            "any_TAD_overlap"
        ]
        .set_index("TAD_map")
        .reindex(TAD_MAP_ORDER)
    )

    contained = (
        sub[
            sub["test"]
            ==
            "complete_TAD_containment"
        ]
        .set_index("TAD_map")
        .reindex(TAD_MAP_ORDER)
    )

    ax.bar(
        x - 0.27,
        overlap["observed"].values,
        width=width,
        label="Any overlap: observed"
    )

    ax.bar(
        x - 0.09,
        overlap["random_median"].values,
        width=width,
        alpha=0.45,
        label="Any overlap: null median"
    )

    ax.bar(
        x + 0.09,
        contained["observed"].values,
        width=width,
        label="Containment: observed"
    )

    ax.bar(
        x + 0.27,
        contained["random_median"].values,
        width=width,
        alpha=0.45,
        label="Containment: null median"
    )

    ax.set_xticks(x)

    ax.set_xticklabels(
        [
            TAD_MAP_LABELS[item]
            for item in TAD_MAP_ORDER
        ],
        rotation=45,
        ha="right",
        fontsize=TICK_FONT
    )

    ax.set_ylim(
        0,
        1.05
    )

    ax.set_ylabel(
        "Fraction of hotspots",
        fontsize=AXIS_LABEL_FONT
    )

    ax.set_title(
        "TAD organization",
        fontsize=TITLE_FONT,
        pad=6
    )

    ax.legend(
        frameon=False,
        fontsize=LEGEND_FONT,
        ncol=2,
        loc="upper right",
        bbox_to_anchor=(
            1.0,
            1.07
        ),
        borderaxespad=0.0
    )

    clean_axes(ax)


# ============================================================
# CTCF
# ============================================================

def plot_ctcf(ax, ctcf_df):

    metrics = [
        "fraction_CTCF_overlap",
        "median_CTCF_peaks_per_100kb",
        "median_CTCF_minimum_edge_distance_bp"
    ]

    metric_labels = [
        "Hotspot overlap",
        "Peak density",
        "Minimum-edge proximity"
    ]

    distance_metrics = [
        "median_CTCF_minimum_edge_distance_bp"
    ]

    x = np.arange(
        len(CTCF_ORDER)
    )

    width = 0.22

    offsets = (
        np.arange(len(metrics))
        -
        (len(metrics) - 1) / 2.0
    ) * width

    all_ratios = []

    for metric_index, metric in enumerate(
        metrics
    ):

        temp = (
            ctcf_df[
                ctcf_df["metric"]
                ==
                metric
            ]
            .set_index("dataset")
            .reindex(CTCF_ORDER)
        )

        observed = (
            temp["observed"]
            .values
            .astype(float)
        )

        null = (
            temp["random_median"]
            .values
            .astype(float)
        )

        if metric in distance_metrics:

            ratio = (
                null
                /
                observed
            )

        else:

            ratio = (
                observed
                /
                null
            )

        all_ratios.extend(
            ratio.tolist()
        )

        ax.bar(
            x + offsets[
                metric_index
            ],
            ratio,
            width=width,
            label=metric_labels[
                metric_index
            ],
            zorder=2
        )

    ax.axhline(
        1.0,
        color="black",
        linewidth=1,
        linestyle="--",
        zorder=6
    )

    ax.set_xticks(x)

    ax.set_xticklabels(
        [
            CTCF_LABELS[item]
            for item in CTCF_ORDER
        ],
        rotation=40,
        ha="right",
        fontsize=TICK_FONT
    )

    ax.set_ylim(
        0,
        max(
            1.15,
            np.nanmax(
                all_ratios
            )
            +
            0.08
        )
    )

    ax.set_ylabel(
        "CTCF association relative to matched null",
        fontsize=AXIS_LABEL_FONT
    )

    # closer to panel than before
    ax.set_title(
        "CTCF association",
        fontsize=TITLE_FONT,
        pad=3
    )

    # closer to x labels than before
    ax.legend(
        frameon=False,
        fontsize=LEGEND_FONT,
        ncol=3,
        loc="upper center",
        bbox_to_anchor=(
            0.5,
            -0.31
        ),
        borderaxespad=0.0,
        columnspacing=1.2,
        handlelength=1.8
    )

    clean_axes(ax)


# ============================================================
# CpG / SEQUENCE
# ============================================================

def plot_cpg(ax, cpg_df):

    metrics = [
        "median_CpG_per_kb",
        "median_CpG_observed_expected",
        "median_GC_fraction",
        "fraction_CpG_island_overlap",
        "median_CpG_island_fraction"
    ]

    labels = [
        "CpG density",
        "CpG O/E",
        "GC fraction",
        "CpG-island overlap",
        "CpG-island coverage"
    ]

    plot_df = (
        cpg_df[
            cpg_df["metric"].isin(
                metrics
            )
        ]
        .set_index("metric")
        .reindex(metrics)
    )

    ratios = (
        plot_df[
            "observed_to_random_median_ratio"
        ]
        .values
        .astype(float)
    )

    pvalues = (
        plot_df[
            "empirical_p_two_sided"
        ]
        .values
        .astype(float)
    )

    x = np.arange(
        len(metrics)
    )

    bars = ax.bar(
        x,
        ratios,
        width=0.65,
        zorder=2
    )

    ax.axhline(
        1.0,
        color="black",
        linewidth=1,
        linestyle="--",
        zorder=10
    )

    ax.set_ylim(
        0,
        max(
            1.15,
            np.nanmax(
                ratios
            )
            +
            0.08
        )
    )

    star_height = 1.035

    for bar, p in zip(
        bars,
        pvalues
    ):

        if (
            np.isfinite(p)
            and
            p < 0.05
        ):

            ax.text(
                bar.get_x()
                +
                bar.get_width()
                /
                2.0,
                star_height,
                "*",
                ha="center",
                va="bottom",
                fontsize=STAR_FONT,
                fontweight="bold",
                zorder=11
            )

    ax.set_xticks(x)

    ax.set_xticklabels(
        labels,
        rotation=35,
        ha="right",
        fontsize=TICK_FONT
    )

    ax.set_ylabel(
        "Observed / matched-null median",
        fontsize=AXIS_LABEL_FONT
    )

    # same title spacing as B
    ax.set_title(
        "CpG and sequence features",
        fontsize=TITLE_FONT,
        pad=3
    )

    clean_axes(ax)


# ============================================================
# HISTONE HEATMAP
# ============================================================

def plot_histone_heatmap(
    ax,
    histone_df
):

    metric_order = [
        "H3K27me3",
        "H3K4me3",
        "H3K27me3:H3K4me3"
    ]

    matrix = (
        histone_df.pivot(
            index="metric_label",
            columns="dataset",
            values="permutation_z"
        )
        .reindex(
            index=metric_order,
            columns=HISTONE_ORDER
        )
    )

    values = matrix.values

    limit = symmetric_limit(
        values
    )

    im = ax.imshow(
        values,
        aspect="auto",
        vmin=-limit,
        vmax=limit,
        cmap="RdBu_r"
    )

    ax.set_xticks(
        np.arange(
            len(HISTONE_ORDER)
        )
    )

    ax.set_xticklabels(
        [
            HISTONE_LABELS[item]
            for item in HISTONE_ORDER
        ],
        rotation=45,
        ha="right",
        fontsize=TICK_FONT
    )

    ax.set_yticks(
        np.arange(
            len(metric_order)
        )
    )

    ax.set_yticklabels(
        metric_order,
        fontsize=TICK_FONT
    )

    for i, metric_label in enumerate(
        metric_order
    ):

        for j, dataset in enumerate(
            HISTONE_ORDER
        ):

            row = histone_df[
                (
                    histone_df[
                        "metric_label"
                    ]
                    ==
                    metric_label
                )
                &
                (
                    histone_df[
                        "dataset"
                    ]
                    ==
                    dataset
                )
            ]

            if len(row) != 1:
                continue

            p = float(
                row[
                    "empirical_p"
                ].iloc[0]
            )

            if p < 0.05:

                ax.text(
                    j,
                    i,
                    "*",
                    ha="center",
                    va="center",
                    fontsize=STAR_FONT,
                    fontweight="bold"
                )

    ax.set_title(
        "Histone-mark enrichment",
        fontsize=TITLE_FONT,
        pad=5
    )

    ax.set_xlabel(
        "Reference context",
        fontsize=AXIS_LABEL_FONT
    )

    ax.set_ylabel("")

    clean_axes(
        ax,
        remove_left=True
    )

    return im


# ============================================================
# CHROMHMM HEATMAP
# ============================================================

def plot_chromhmm_heatmap(
    ax,
    chromhmm_df
):

    matrix = (
        chromhmm_df.pivot(
            index="state",
            columns="epigenome",
            values="permutation_z"
        )
        .reindex(
            index=STATE_ORDER,
            columns=EPIGENOME_ORDER
        )
    )

    values = matrix.values

    limit = symmetric_limit(
        values
    )

    im = ax.imshow(
        values,
        aspect="auto",
        vmin=-limit,
        vmax=limit,
        cmap="RdBu_r"
    )

    ax.set_xticks(
        np.arange(
            len(EPIGENOME_ORDER)
        )
    )

    ax.set_xticklabels(
        [
            EPIGENOME_LABELS[item]
            for item in EPIGENOME_ORDER
        ],
        rotation=45,
        ha="right",
        fontsize=SMALL_TICK_FONT
    )

    ax.set_yticks(
        np.arange(
            len(STATE_ORDER)
        )
    )

    ax.set_yticklabels(
        STATE_ORDER,
        fontsize=SMALL_TICK_FONT
    )

    if (
        "FDR_BH_all_tests"
        in chromhmm_df.columns
    ):

        for i, state in enumerate(
            STATE_ORDER
        ):

            for j, epigenome in enumerate(
                EPIGENOME_ORDER
            ):

                row = chromhmm_df[
                    (
                        chromhmm_df[
                            "state"
                        ]
                        ==
                        state
                    )
                    &
                    (
                        chromhmm_df[
                            "epigenome"
                        ]
                        ==
                        epigenome
                    )
                ]

                if len(row) != 1:
                    continue

                q = float(
                    row[
                        "FDR_BH_all_tests"
                    ].iloc[0]
                )

                if q < 0.05:

                    ax.text(
                        j,
                        i,
                        "*",
                        ha="center",
                        va="center",
                        fontsize=STAR_FONT,
                        fontweight="bold"
                    )

    ax.set_title(
        "ChromHMM state enrichment",
        fontsize=TITLE_FONT,
        pad=5
    )

    ax.set_xlabel(
        "Reference epigenome",
        fontsize=AXIS_LABEL_FONT
    )

    ax.set_ylabel(
        "ChromHMM state",
        fontsize=AXIS_LABEL_FONT
    )

    clean_axes(
        ax,
        remove_left=True
    )

    return im


# ============================================================
# CROSS-EPIGENOME CHROMHMM SUMMARY
# ============================================================

def plot_cross_chromhmm(
    ax,
    cross_df
):

    plot_df = (
        cross_df
        .set_index("state")
        .reindex(STATE_ORDER)
        .reset_index()
    )

    if (
        "median_permutation_z"
        in plot_df.columns
    ):

        median_z_col = (
            "median_permutation_z"
        )

    elif (
        "median_z"
        in plot_df.columns
    ):

        median_z_col = (
            "median_z"
        )

    else:

        raise ValueError(
            "No median Z-score column found "
            "in cross-epigenome summary."
        )

    if (
        "n_z_positive"
        in plot_df.columns
    ):

        n_pos_col = (
            "n_z_positive"
        )

    elif (
        "n_pos"
        in plot_df.columns
    ):

        n_pos_col = (
            "n_pos"
        )

    else:

        raise ValueError(
            "No positive-Z count column found."
        )

    if (
        "n_z_negative"
        in plot_df.columns
    ):

        n_neg_col = (
            "n_z_negative"
        )

    elif (
        "n_neg"
        in plot_df.columns
    ):

        n_neg_col = (
            "n_neg"
        )

    else:

        raise ValueError(
            "No negative-Z count column found."
        )

    y = np.arange(
        len(plot_df)
    )

    x = (
        plot_df[
            median_z_col
        ]
        .values
        .astype(float)
    )

    ax.axvline(
        0,
        linewidth=1,
        color="black"
    )

    for yi, value in zip(
        y,
        x
    ):

        ax.plot(
            [
                0,
                value
            ],
            [
                yi,
                yi
            ],
            linewidth=1,
            color="0.65"
        )

    ax.scatter(
        x,
        y,
        s=50,
        color="black",
        zorder=3
    )

    ax.set_yticks(y)

    ax.set_yticklabels(
        plot_df["state"],
        fontsize=SMALL_TICK_FONT
    )

    ax.invert_yaxis()

    ax.set_xlabel(
        "Median permutation Z-score",
        fontsize=AXIS_LABEL_FONT
    )

    ax.set_ylabel("")

    ax.set_title(
        "Cross-epigenome ChromHMM summary",
        fontsize=TITLE_FONT,
        pad=5
    )

    xmin = np.nanmin(x)
    xmax = np.nanmax(x)

    span = (
        xmax
        -
        xmin
    )

    if span <= 0:
        span = 1.0

    annotation_x = (
        xmax
        +
        0.18
        *
        span
    )

    for i, row in plot_df.iterrows():

        positive = int(
            row[
                n_pos_col
            ]
        )

        negative = int(
            row[
                n_neg_col
            ]
        )

        directional = (
            positive
            +
            negative
        )

        if positive > negative:

            text_value = (
                "{}/{} enriched"
            ).format(
                positive,
                directional
            )

        elif negative > positive:

            text_value = (
                "{}/{} depleted"
            ).format(
                negative,
                directional
            )

        else:

            text_value = "mixed"

        ax.text(
            annotation_x,
            i,
            text_value,
            va="center",
            fontsize=SMALL_TICK_FONT
        )

    ax.set_xlim(
        xmin
        -
        0.15
        *
        span,
        xmax
        +
        0.55
        *
        span
    )

    clean_axes(ax)


# ============================================================
# INDIVIDUAL FIGURES
# ============================================================

def make_tad_figure(
    tad_df
):

    fig, ax = plt.subplots(
        figsize=(
            10,
            4.5
        )
    )

    plot_tad(
        ax,
        tad_df
    )

    fig.tight_layout()

    save_figure(
        fig,
        "TAD_organization"
    )


def make_ctcf_figure(
    ctcf_df
):

    fig, ax = plt.subplots(
        figsize=(
            9,
            5.4
        )
    )

    plot_ctcf(
        ax,
        ctcf_df
    )

    # Less bottom space than before because legend moved closer
    fig.subplots_adjust(
        top=0.93,
        bottom=0.33,
        left=0.13,
        right=0.98
    )

    save_figure(
        fig,
        "CTCF_association"
    )


def make_cpg_figure(
    cpg_df
):

    fig, ax = plt.subplots(
        figsize=(
            7.5,
            4.8
        )
    )

    plot_cpg(
        ax,
        cpg_df
    )

    fig.tight_layout()

    save_figure(
        fig,
        "CpG_sequence_features"
    )


def make_histone_figure(
    histone_df
):

    fig, ax = plt.subplots(
        figsize=(
            10,
            3.6
        )
    )

    im = plot_histone_heatmap(
        ax,
        histone_df
    )

    cbar = fig.colorbar(
        im,
        ax=ax
    )

    cbar.set_label(
        "Permutation Z-score",
        fontsize=COLORBAR_FONT
    )

    cbar.ax.tick_params(
        labelsize=SMALL_TICK_FONT
    )

    fig.tight_layout()

    save_figure(
        fig,
        "histone_mark_enrichment"
    )


def make_chromhmm_heatmap_figure(
    chromhmm_df
):

    fig, ax = plt.subplots(
        figsize=(
            11,
            7.5
        )
    )

    im = plot_chromhmm_heatmap(
        ax,
        chromhmm_df
    )

    cbar = fig.colorbar(
        im,
        ax=ax
    )

    cbar.set_label(
        "Permutation Z-score",
        fontsize=COLORBAR_FONT
    )

    cbar.ax.tick_params(
        labelsize=SMALL_TICK_FONT
    )

    fig.tight_layout()

    save_figure(
        fig,
        "ChromHMM_state_enrichment"
    )


def make_cross_chromhmm_figure(
    cross_df
):

    fig, ax = plt.subplots(
        figsize=(
            7.5,
            7.5
        )
    )

    plot_cross_chromhmm(
        ax,
        cross_df
    )

    fig.tight_layout()

    save_figure(
        fig,
        "ChromHMM_cross_epigenome_summary"
    )


# ============================================================
# COMBINED FIGURE
# ============================================================

def make_combined_figure(
    tad_df,
    ctcf_df,
    cpg_df,
    histone_df,
    chromhmm_df,
    cross_df
):

    fig = plt.figure(
        figsize=(
            16,
            16
        ),
        constrained_layout=True
    )

    gs = fig.add_gridspec(
        4,
        2,
        height_ratios=[
            0.82,
            1.00,
            0.58,
            1.42
        ],
        hspace=0.08,
        wspace=0.16
    )

    # --------------------------------------------------------
    # TAD
    # --------------------------------------------------------

    ax_a = fig.add_subplot(
        gs[
            0,
            :
        ]
    )

    plot_tad(
        ax_a,
        tad_df
    )

    # --------------------------------------------------------
    # CTCF
    # --------------------------------------------------------

    ax_b = fig.add_subplot(
        gs[
            1,
            0
        ]
    )

    plot_ctcf(
        ax_b,
        ctcf_df
    )

    # --------------------------------------------------------
    # CpG
    # --------------------------------------------------------

    ax_c = fig.add_subplot(
        gs[
            1,
            1
        ]
    )

    plot_cpg(
        ax_c,
        cpg_df
    )

    # --------------------------------------------------------
    # Histone
    # --------------------------------------------------------

    ax_d = fig.add_subplot(
        gs[
            2,
            :
        ]
    )

    hist_im = plot_histone_heatmap(
        ax_d,
        histone_df
    )

    cbar_d = fig.colorbar(
        hist_im,
        ax=ax_d,
        fraction=0.025,
        pad=0.02
    )

    cbar_d.set_label(
        "Permutation Z-score",
        fontsize=COLORBAR_FONT
    )

    cbar_d.ax.tick_params(
        labelsize=SMALL_TICK_FONT
    )

    # --------------------------------------------------------
    # Bottom row
    # --------------------------------------------------------

    gs_bottom = gs[
        3,
        :
    ].subgridspec(
        1,
        2,
        width_ratios=[
            1.35,
            1.0
        ],
        wspace=0.16
    )

    # --------------------------------------------------------
    # ChromHMM heatmap
    # --------------------------------------------------------

    ax_e = fig.add_subplot(
        gs_bottom[
            0,
            0
        ]
    )

    chrom_im = plot_chromhmm_heatmap(
        ax_e,
        chromhmm_df
    )

    cbar_e = fig.colorbar(
        chrom_im,
        ax=ax_e,
        fraction=0.038,
        pad=0.03
    )

    cbar_e.set_label(
        "Permutation Z-score",
        fontsize=COLORBAR_FONT
    )

    cbar_e.ax.tick_params(
        labelsize=SMALL_TICK_FONT
    )

    # --------------------------------------------------------
    # Cross-epigenome summary
    # --------------------------------------------------------

    ax_f = fig.add_subplot(
        gs_bottom[
            0,
            1
        ]
    )

    plot_cross_chromhmm(
        ax_f,
        cross_df
    )

    save_figure(
        fig,
        "supplementary_functional_feature_analysis"
    )


# ============================================================
# MAIN
# ============================================================

def main():

    print(
        "Loading TAD results..."
    )

    tad_df = load_tad()

    print(
        "Loading CTCF results..."
    )

    ctcf_df = load_ctcf()

    print(
        "Loading CpG results..."
    )

    cpg_df = load_cpg()

    print(
        "Loading histone results..."
    )

    histone_df = load_histone()

    print(
        "Loading ChromHMM results..."
    )

    chromhmm_df = load_chromhmm()

    print(
        "Loading cross-epigenome "
        "ChromHMM results..."
    )

    cross_df = load_cross_chromhmm()

    print(
        "Saving TAD figure..."
    )

    make_tad_figure(
        tad_df
    )

    print(
        "Saving CTCF figure..."
    )

    make_ctcf_figure(
        ctcf_df
    )

    print(
        "Saving CpG figure..."
    )

    make_cpg_figure(
        cpg_df
    )

    print(
        "Saving histone figure..."
    )

    make_histone_figure(
        histone_df
    )

    print(
        "Saving ChromHMM heatmap..."
    )

    make_chromhmm_heatmap_figure(
        chromhmm_df
    )

    print(
        "Saving cross-epigenome "
        "ChromHMM summary..."
    )

    make_cross_chromhmm_figure(
        cross_df
    )

    print(
        "Saving combined figure..."
    )

    make_combined_figure(
        tad_df,
        ctcf_df,
        cpg_df,
        histone_df,
        chromhmm_df,
        cross_df
    )

    print(
        "\nFigures written to:"
    )

    print(
        OUTDIR
    )


if __name__ == "__main__":
    main()