#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
try:
    from scipy.stats import binomtest
except ImportError:
    from scipy.stats import binom_test as _binom_test
    class _BinomResult:
        def __init__(self, pvalue):
            self.pvalue = pvalue
    def binomtest(k, n, p=0.5, alternative="two-sided"):
        return _BinomResult(_binom_test(k, n, p=p, alternative=alternative))


# ============================================================
# PATHS
# ============================================================

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
BASE_DIR = os.path.join(SCRIPT_DIR, "results", "ChromHMM_analysis")

SUMMARY_FILE = os.path.join(
    BASE_DIR,
    "ALL_EPIGENOMES_ChromHMM_permutation_summary.tsv"
)

OUTDIR = os.path.join(
    BASE_DIR,
    "cross_epigenome_summary"
)

# Portable input/output overrides (applied before loading the summary).
import argparse
_parser = argparse.ArgumentParser(description="Summarize ChromHMM states across epigenomes")
_parser.add_argument("--input", default=SUMMARY_FILE, help="ALL_EPIGENOMES_ChromHMM_permutation_summary.tsv")
_parser.add_argument("--outdir", default=OUTDIR)
_args = _parser.parse_args()
SUMMARY_FILE = _args.input
OUTDIR = _args.outdir

os.makedirs(
    OUTDIR,
    exist_ok=True
)


# ============================================================
# STATE ORDER
# ============================================================

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
# FRIENDLY EPIGENOME NAMES
# ============================================================

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


# ============================================================
# LOAD
# ============================================================

df = pd.read_csv(
    SUMMARY_FILE,
    sep="\t"
)

print(
    "Rows loaded: {:,}".format(
        len(df)
    )
)

print(
    "Epigenomes: {}".format(
        df["epigenome"].nunique()
    )
)

print(
    "States: {}".format(
        df["state"].nunique()
    )
)


# ============================================================
# CROSS-EPIGENOME SUMMARY
# ============================================================

summary_rows = []

for state in STATE_ORDER:

    sub = df[
        df["state"]
        ==
        state
    ].copy()

    z = sub[
        "permutation_z"
    ].replace(
        [np.inf, -np.inf],
        np.nan
    ).dropna()

    log2_enrichment = sub[
        "log2_observed_to_random_mean"
    ].replace(
        [np.inf, -np.inf],
        np.nan
    ).dropna()

    n_total = len(
        sub
    )

    n_positive = int(
        np.sum(
            z > 0
        )
    )

    n_negative = int(
        np.sum(
            z < 0
        )
    )

    # --------------------------------------------------------
    # Sign consistency test
    #
    # Under null, positive/negative directions are equally
    # likely. Test whichever direction is more common.
    # --------------------------------------------------------

    n_directional = (
        n_positive
        +
        n_negative
    )

    if n_directional > 0:

        dominant_count = max(
            n_positive,
            n_negative
        )

        sign_consistency_p = (
            binomtest(
                dominant_count,
                n=n_directional,
                p=0.5,
                alternative="two-sided"
            ).pvalue
        )

    else:

        sign_consistency_p = np.nan

    if n_positive > n_negative:

        dominant_direction = "enriched"

    elif n_negative > n_positive:

        dominant_direction = "depleted"

    else:

        dominant_direction = "mixed"

    summary_rows.append({

        "state":
            state,

        "n_epigenomes":
            n_total,

        "median_permutation_z":
            np.nanmedian(
                z
            ),

        "mean_permutation_z":
            np.nanmean(
                z
            ),

        "min_permutation_z":
            np.nanmin(
                z
            ),

        "max_permutation_z":
            np.nanmax(
                z
            ),

        "n_z_positive":
            n_positive,

        "n_z_negative":
            n_negative,

        "dominant_direction":
            dominant_direction,

        "sign_consistency_p":
            sign_consistency_p,

        "median_log2_observed_expected":
            np.nanmedian(
                log2_enrichment
            ),

        "mean_log2_observed_expected":
            np.nanmean(
                log2_enrichment
            ),

        "n_nominal_p_lt_0.05":
            int(
                np.sum(
                    sub[
                        "empirical_p_two_sided"
                    ]
                    <
                    0.05
                )
            ),

        "n_within_epigenome_FDR_lt_0.05":
            int(
                np.sum(
                    sub[
                        "FDR_BH_within_epigenome"
                    ]
                    <
                    0.05
                )
            ),

        "n_global_FDR_lt_0.05":
            int(
                np.sum(
                    sub[
                        "FDR_BH_all_tests"
                    ]
                    <
                    0.05
                )
            )
    })


cross_summary = pd.DataFrame(
    summary_rows
)


# ============================================================
# BH FDR FOR SIGN-CONSISTENCY TESTS
# ============================================================

def benjamini_hochberg(pvalues):

    pvalues = np.asarray(
        pvalues,
        dtype=float
    )

    n = len(
        pvalues
    )

    order = np.argsort(
        pvalues
    )

    ranked = pvalues[
        order
    ]

    adjusted = np.empty(
        n,
        dtype=float
    )

    running_min = 1.0

    for i in range(
        n - 1,
        -1,
        -1
    ):

        rank = (
            i + 1
        )

        value = (
            ranked[i]
            *
            n
            /
            float(rank)
        )

        running_min = min(
            running_min,
            value
        )

        adjusted[
            order[i]
        ] = min(
            running_min,
            1.0
        )

    return adjusted


cross_summary[
    "sign_consistency_FDR"
] = benjamini_hochberg(
    cross_summary[
        "sign_consistency_p"
    ].values
)


# ============================================================
# SAVE SUMMARY
# ============================================================

cross_summary.to_csv(
    os.path.join(
        OUTDIR,
        "ChromHMM_cross_epigenome_state_summary.tsv"
    ),
    sep="\t",
    index=False
)


print(
    "\n=============================================="
)

print(
    "CROSS-EPIGENOME STATE SUMMARY"
)

print(
    "==============================================\n"
)

print(
    cross_summary[
        [
            "state",
            "median_permutation_z",
            "n_z_positive",
            "n_z_negative",
            "dominant_direction",
            "sign_consistency_p",
            "sign_consistency_FDR",
            "n_nominal_p_lt_0.05",
            "n_global_FDR_lt_0.05"
        ]
    ].to_string(
        index=False
    )
)


# ============================================================
# HEATMAP MATRIX
# ============================================================

z_matrix = df.pivot(
    index="state",
    columns="epigenome",
    values="permutation_z"
)

z_matrix = z_matrix.reindex(
    index=STATE_ORDER,
    columns=EPIGENOME_ORDER
)

z_matrix.columns = [
    EPIGENOME_LABELS[x]
    for x in z_matrix.columns
]

z_matrix.to_csv(
    os.path.join(
        OUTDIR,
        "ChromHMM_cross_epigenome_z_matrix.tsv"
    ),
    sep="\t"
)


# ============================================================
# HEATMAP
# ============================================================

fig, ax = plt.subplots(
    figsize=(11, 9)
)

matrix_values = (
    z_matrix.values
)

# Symmetric limits make enrichment/depletion visually comparable
finite_values = matrix_values[
    np.isfinite(
        matrix_values
    )
]

limit = np.percentile(
    np.abs(
        finite_values
    ),
    97.5
)

limit = max(
    limit,
    2.0
)

im = ax.imshow(
    matrix_values,
    aspect="auto",
    vmin=-limit,
    vmax=limit
)

ax.set_xticks(
    np.arange(
        len(
            z_matrix.columns
        )
    )
)

ax.set_xticklabels(
    z_matrix.columns,
    rotation=45,
    ha="right"
)

ax.set_yticks(
    np.arange(
        len(
            z_matrix.index
        )
    )
)

ax.set_yticklabels(
    z_matrix.index
)


# ------------------------------------------------------------
# Mark globally FDR-significant cells
# ------------------------------------------------------------

for i, state in enumerate(
    z_matrix.index
):

    for j, epigenome_label in enumerate(
        z_matrix.columns
    ):

        original_epigenome = (
            EPIGENOME_ORDER[j]
        )

        row = df[
            (
                df["state"]
                ==
                state
            )
            &
            (
                df["epigenome"]
                ==
                original_epigenome
            )
        ]

        if len(row) == 1:

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
                    fontsize=11
                )


cbar = fig.colorbar(
    im,
    ax=ax
)

cbar.set_label(
    "Permutation Z-score"
)

ax.set_title(
    "ChromHMM state enrichment across reference epigenomes"
)

ax.set_xlabel(
    "Reference epigenome"
)

ax.set_ylabel(
    "18-state ChromHMM state"
)

fig.tight_layout()

fig.savefig(
    os.path.join(
        OUTDIR,
        "ChromHMM_cross_epigenome_heatmap.png"
    ),
    dpi=300,
    bbox_inches="tight"
)

fig.savefig(
    os.path.join(
        OUTDIR,
        "ChromHMM_cross_epigenome_heatmap.pdf"
    ),
    bbox_inches="tight"
)

plt.close(
    fig
)


# ============================================================
# STATE-LEVEL SUMMARY PLOT
# ============================================================

plot_df = cross_summary.copy()

plot_df = plot_df.set_index(
    "state"
).loc[
    STATE_ORDER
].reset_index()


fig, ax = plt.subplots(
    figsize=(7, 8)
)

y = np.arange(
    len(
        plot_df
    )
)

ax.scatter(
    plot_df[
        "median_permutation_z"
    ],
    y,
    s=50
)

ax.axvline(
    0,
    linewidth=1
)

ax.set_yticks(
    y
)

ax.set_yticklabels(
    plot_df[
        "state"
    ]
)

ax.set_xlabel(
    "Median permutation Z-score across epigenomes"
)

ax.set_ylabel(
    ""
)

ax.set_title(
    "Cross-epigenome ChromHMM enrichment"
)

ax.invert_yaxis()

fig.tight_layout()

fig.savefig(
    os.path.join(
        OUTDIR,
        "ChromHMM_state_median_z.png"
    ),
    dpi=300,
    bbox_inches="tight"
)

fig.savefig(
    os.path.join(
        OUTDIR,
        "ChromHMM_state_median_z.pdf"
    ),
    bbox_inches="tight"
)

plt.close(
    fig
)


print(
    "\nResults written to:"
)

print(
    OUTDIR
)