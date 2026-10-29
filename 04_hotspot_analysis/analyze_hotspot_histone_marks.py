#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import numpy as np
import pandas as pd
import pyBigWig


# ============================================================
# USER SETTINGS
# ============================================================

# Repository inputs and shared GSC annotations (override via CLI if needed).
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
HOTSPOT_FILE = os.path.join(SCRIPT_DIR, "data", "methylated_promoter_hotspots.bed")
PROMOTER_FILE = os.path.join(SCRIPT_DIR, "data", "input_promoters_coordinates_methylation.tsv")
CHROMSIZES_FILE = "/gsc/www/bcgsc.ca/downloads/promoter_methylation_nanopore/annotations/genome/hg38.chrom.sizes"


H3K27_BASE = (
    "/projects/clininfo/vcsizmok/long_POG/aDMR/aDMR_analysis/"
    "promoter_clustering/H3K27me3"
)

H3K4_BASE = (
    "/projects/clininfo/vcsizmok/long_POG/aDMR/aDMR_analysis/"
    "promoter_clustering/H3K4me3"
)

OUTDIR = os.path.join(SCRIPT_DIR, "results", 'functional_feature_analysis', 'all_tissues_histone_comparison')

N_PERMUTATIONS = 10000
RANDOM_SEED = 12345

# Signal interval around promoter TSS
TSS_FLANK = 2000

# Pseudocount for H3K27me3/H3K4me3 ratio
RATIO_PSEUDOCOUNT = 1.0


# ============================================================
# HISTONE DATASETS
# ============================================================

DATASETS = {

    "MCF7": {
        "H3K27me3": os.path.join(
            H3K27_BASE,
            "MCF-7",
            "ENCFF210IBW.bigWig"
        ),
        "H3K4me3": os.path.join(
            H3K4_BASE,
            "MCF-7",
            "ENCFF267OQQ.bigWig"
        )
    },

    "HMEC": {
        "H3K27me3": os.path.join(
            H3K27_BASE,
            "HMEC",
            "ENCFF530KJI.bigWig"
        ),
        "H3K4me3": os.path.join(
            H3K4_BASE,
            "HMEC",
            "ENCFF977KZW.bigWig"
        )
    },

    "vHMEC": {
        "H3K27me3": os.path.join(
            H3K27_BASE,
            "vHMEC",
            "ENCFF121GTV.bigWig"
        ),
        "H3K4me3": os.path.join(
            H3K4_BASE,
            "vHMEC",
            "ENCFF526WUE.bigWig"
        )
    },

    "Colon": {
        "H3K27me3": os.path.join(
            H3K27_BASE,
            "colon",
            "ENCFF899SKW.bigWig"
        ),
        "H3K4me3": os.path.join(
            H3K4_BASE,
            "colon",
            "ENCFF487CTD.bigWig"
        )
    },

    "Liver": {
        "H3K27me3": os.path.join(
            H3K27_BASE,
            "liver",
            "ENCFF443GFE.bigWig"
        ),
        "H3K4me3": os.path.join(
            H3K4_BASE,
            "liver",
            "ENCFF280QYJ.bigWig"
        )
    },

    "Pancreas": {
        "H3K27me3": os.path.join(
            H3K27_BASE,
            "pancreas",
            "ENCFF209UPZ.bigWig"
        ),
        "H3K4me3": os.path.join(
            H3K4_BASE,
            "pancreas",
            "ENCFF285STS.bigWig"
        )
    },

    "Lung": {
        "H3K27me3": os.path.join(
            H3K27_BASE,
            "lung",
            "ENCFF893TFN.bigWig"
        ),
        "H3K4me3": os.path.join(
            H3K4_BASE,
            "lung",
            "ENCFF282VQS.bigWig"
        )
    },

    "Astrocytes": {
        "H3K27me3": os.path.join(
            H3K27_BASE,
            "astrocytes",
            "ENCFF665DLZ.bigWig"
        ),
        "H3K4me3": os.path.join(
            H3K4_BASE,
            "astrocytes",
            "ENCFF577BWJ.bigWig"
        )
    }
}


# ============================================================
# LOAD HOTSPOTS
# ============================================================

def load_hotspots(filename):

    df = pd.read_csv(
        filename,
        sep="\t"
    )

    df.columns = df.columns.str.strip()

    required = [
        "chr",
        "start",
        "end",
        "n_promoters"
    ]

    missing = [
        x for x in required
        if x not in df.columns
    ]

    if missing:
        raise ValueError(
            "Missing hotspot columns: "
            + ", ".join(missing)
        )

    df = df.rename(
        columns={
            "chr": "chrom"
        }
    )

    df["start"] = pd.to_numeric(
        df["start"],
        errors="raise"
    ).astype(int)

    df["end"] = pd.to_numeric(
        df["end"],
        errors="raise"
    ).astype(int)

    df["n_promoters"] = pd.to_numeric(
        df["n_promoters"],
        errors="raise"
    ).astype(int)

    df["length"] = (
        df["end"]
        -
        df["start"]
    )

    df["hotspot_id"] = [
        "hotspot_{:03d}".format(i + 1)
        for i in range(len(df))
    ]

    return df


# ============================================================
# LOAD PROMOTERS
# ============================================================

def load_promoters(filename):

    df = pd.read_csv(
        filename,
        sep="\t"
    )

    df.columns = df.columns.str.strip()

    required = [
        "symbol",
        "chrom",
        "tss_start",
        "promoter_start",
        "promoter_end"
    ]

    missing = [
        x for x in required
        if x not in df.columns
    ]

    if missing:
        raise ValueError(
            "Missing promoter columns: "
            + ", ".join(missing)
        )

    df["tss_start"] = pd.to_numeric(
        df["tss_start"],
        errors="raise"
    ).astype(int)

    df["promoter_start"] = pd.to_numeric(
        df["promoter_start"],
        errors="raise"
    )

    df["promoter_end"] = pd.to_numeric(
        df["promoter_end"],
        errors="raise"
    )

    # Important:
    # reproduce the promoter position used to define the hotspots
    df["midpoint"] = (
        (
            df["promoter_start"]
            +
            df["promoter_end"]
        ) // 2
    ).astype(int)

    df = df.reset_index(
        drop=True
    )

    df["promoter_index"] = df.index

    return df


# ============================================================
# LOAD CHROMOSOME SIZES
# ============================================================

def load_chromsizes(filename):

    df = pd.read_csv(
        filename,
        sep="\t",
        header=None,
        names=[
            "chrom",
            "size"
        ]
    )

    df["size"] = df["size"].astype(int)

    return dict(
        zip(
            df["chrom"],
            df["size"]
        )
    )


# ============================================================
# VALIDATE BIGWIG FILES
# ============================================================

def validate_bigwig_files():

    print(
        "\nChecking histone bigWig files..."
    )

    error = False

    for dataset, paths in DATASETS.items():

        for mark in [
            "H3K27me3",
            "H3K4me3"
        ]:

            filename = paths[
                mark
            ]

            exists = os.path.exists(
                filename
            )

            print(
                "{:<12s} {:<10s} {}".format(
                    dataset,
                    mark,
                    "OK" if exists else "MISSING"
                )
            )

            if not exists:
                print(
                    "    {}".format(
                        filename
                    )
                )

                error = True

    if error:

        raise IOError(
            "One or more bigWig files are missing."
        )


# ============================================================
# PROMOTER LOOKUP
# ============================================================

def build_promoter_lookup(promoters):

    lookup = {}

    for chrom, group in promoters.groupby(
        "chrom"
    ):

        group = group.sort_values(
            "midpoint"
        )

        lookup[
            chrom
        ] = {

            "positions":
                group[
                    "midpoint"
                ]
                .astype(int)
                .to_numpy(),

            "indices":
                group[
                    "promoter_index"
                ]
                .astype(int)
                .to_numpy()
        }

    return lookup


def get_promoter_indices_in_region(
    chrom,
    start,
    end,
    promoter_lookup
):

    if chrom not in promoter_lookup:

        return np.array(
            [],
            dtype=int
        )

    positions = promoter_lookup[
        chrom
    ]["positions"]

    indices = promoter_lookup[
        chrom
    ]["indices"]

    left = np.searchsorted(
        positions,
        start,
        side="left"
    )

    right = np.searchsorted(
        positions,
        end,
        side="right"
    )

    return indices[
        left:right
    ]


# ============================================================
# VALIDATE HOTSPOT PROMOTER COUNTS
# ============================================================

def validate_hotspot_promoter_counts(
    hotspots,
    promoter_lookup
):

    counts = []

    for _, row in hotspots.iterrows():

        idx = get_promoter_indices_in_region(
            row["chrom"],
            int(row["start"]),
            int(row["end"]),
            promoter_lookup
        )

        counts.append(
            len(idx)
        )

    hotspots = hotspots.copy()

    hotspots[
        "promoter_count_recalculated"
    ] = counts

    matches = (
        hotspots[
            "n_promoters"
        ]
        ==
        hotspots[
            "promoter_count_recalculated"
        ]
    )

    print(
        "\nPromoter-count validation:"
    )

    print(
        "{} / {} hotspots reproduce "
        "the n_promoters value".format(
            int(matches.sum()),
            len(hotspots)
        )
    )

    if not matches.all():

        print(
            hotspots.loc[
                ~matches,
                [
                    "chrom",
                    "start",
                    "end",
                    "n_promoters",
                    "promoter_count_recalculated"
                ]
            ].to_string(
                index=False
            )
        )

        raise ValueError(
            "Hotspot promoter-count validation failed."
        )

    return hotspots


# ============================================================
# VALID RANDOM INTERVALS
# ============================================================

def get_valid_start_ranges(
    chrom,
    length,
    target_n_promoters,
    chromsizes,
    promoter_lookup
):

    chrom_size = int(
        chromsizes[
            chrom
        ]
    )

    max_start = (
        chrom_size
        -
        length
    )

    positions = promoter_lookup[
        chrom
    ]["positions"]

    n = len(
        positions
    )

    k = int(
        target_n_promoters
    )

    ranges = []

    if k > n:
        return ranges

    for i in range(
        0,
        n - k + 1
    ):

        j = i + k - 1

        first = int(
            positions[i]
        )

        last = int(
            positions[j]
        )

        if (
            last - first
            >
            length
        ):
            continue

        lower = max(
            0,
            last - length
        )

        upper = min(
            max_start,
            first
        )

        if i > 0:

            prev_pos = int(
                positions[
                    i - 1
                ]
            )

            lower = max(
                lower,
                prev_pos + 1
            )

        if j + 1 < n:

            next_pos = int(
                positions[
                    j + 1
                ]
            )

            upper = min(
                upper,
                next_pos
                -
                length
                -
                1
            )

        if lower <= upper:

            ranges.append(
                (
                    int(lower),
                    int(upper)
                )
            )

    return ranges


def precompute_valid_start_ranges(
    hotspots,
    chromsizes,
    promoter_lookup
):

    cache = {}

    combinations = (
        hotspots[
            [
                "chrom",
                "length",
                "n_promoters"
            ]
        ]
        .drop_duplicates()
    )

    print(
        "\nPrecomputing promoter-matched "
        "random interval possibilities..."
    )

    for _, row in combinations.iterrows():

        key = (
            row["chrom"],
            int(row["length"]),
            int(row["n_promoters"])
        )

        ranges = get_valid_start_ranges(
            key[0],
            key[1],
            key[2],
            chromsizes,
            promoter_lookup
        )

        if len(ranges) == 0:

            raise RuntimeError(
                "No matched random intervals for "
                "{} length={} n_promoters={}".format(
                    key[0],
                    key[1],
                    key[2]
                )
            )

        sizes = np.array(
            [
                upper - lower + 1
                for lower, upper
                in ranges
            ],
            dtype=np.int64
        )

        cache[
            key
        ] = {
            "ranges":
                ranges,

            "cumulative":
                np.cumsum(
                    sizes
                ),

            "total":
                int(
                    sizes.sum()
                )
        }

    return cache


def sample_random_start(
    chrom,
    length,
    n_promoters,
    cache,
    rng
):

    key = (
        chrom,
        length,
        n_promoters
    )

    item = cache[
        key
    ]

    x = int(
        rng.randint(
            0,
            item[
                "total"
            ]
        )
    )

    range_idx = int(
        np.searchsorted(
            item[
                "cumulative"
            ],
            x,
            side="right"
        )
    )

    previous = (
        int(
            item[
                "cumulative"
            ][
                range_idx - 1
            ]
        )
        if range_idx > 0
        else 0
    )

    lower, upper = (
        item[
            "ranges"
        ][
            range_idx
        ]
    )

    return int(
        lower
        +
        (
            x - previous
        )
    )


# ============================================================
# GENERATE RANDOM HOTSPOT SET
# ============================================================

def generate_random_regions(
    hotspots,
    cache,
    rng
):

    rows = []

    for _, row in hotspots.iterrows():

        length = int(
            row[
                "length"
            ]
        )

        start = sample_random_start(
            row["chrom"],
            length,
            int(row["n_promoters"]),
            cache,
            rng
        )

        rows.append({
            "hotspot_id":
                row["hotspot_id"],

            "chrom":
                row["chrom"],

            "start":
                start,

            "end":
                start + length,

            "length":
                length,

            "n_promoters":
                int(
                    row[
                        "n_promoters"
                    ]
                )
        })

    return pd.DataFrame(
        rows
    )


# ============================================================
# PROMOTER INDICES PER REGION
# ============================================================

def get_region_promoter_indices(
    regions,
    promoter_lookup
):

    output = []

    for _, row in regions.iterrows():

        idx = get_promoter_indices_in_region(
            row["chrom"],
            int(row["start"]),
            int(row["end"]),
            promoter_lookup
        )

        output.append(
            idx
        )

    return output


# ============================================================
# BIGWIG SIGNAL
# ============================================================

def get_bigwig_mean_signal(
    bw,
    chrom,
    start,
    end
):

    chroms = bw.chroms()

    if chrom not in chroms:

        return np.nan

    chrom_length = int(
        chroms[
            chrom
        ]
    )

    start = max(
        0,
        int(start)
    )

    end = min(
        chrom_length,
        int(end)
    )

    if end <= start:
        return np.nan

    try:

        value = bw.stats(
            chrom,
            start,
            end,
            type="mean",
            exact=True
        )[0]

    except RuntimeError:

        return np.nan

    if value is None:
        return np.nan

    if not np.isfinite(
        value
    ):
        return np.nan

    return float(
        value
    )


# ============================================================
# PRECOMPUTE PROMOTER SIGNALS
# ============================================================

def calculate_promoter_signals(
    promoters,
    h3k27_file,
    h3k4_file
):

    bw27 = pyBigWig.open(
        h3k27_file
    )

    bw4 = pyBigWig.open(
        h3k4_file
    )

    h3k27 = np.full(
        len(promoters),
        np.nan
    )

    h3k4 = np.full(
        len(promoters),
        np.nan
    )

    ratio = np.full(
        len(promoters),
        np.nan
    )

    for i, row in promoters.iterrows():

        tss = int(
            row[
                "tss_start"
            ]
        )

        start = (
            tss
            -
            TSS_FLANK
        )

        end = (
            tss
            +
            TSS_FLANK
        )

        signal27 = get_bigwig_mean_signal(
            bw27,
            row["chrom"],
            start,
            end
        )

        signal4 = get_bigwig_mean_signal(
            bw4,
            row["chrom"],
            start,
            end
        )

        h3k27[i] = signal27
        h3k4[i] = signal4

        if (
            np.isfinite(signal27)
            and
            np.isfinite(signal4)
        ):

            ratio[i] = np.log2(
                (
                    signal27
                    +
                    RATIO_PSEUDOCOUNT
                )
                /
                (
                    signal4
                    +
                    RATIO_PSEUDOCOUNT
                )
            )

        if (
            (i + 1) % 5000
            ==
            0
        ):

            print(
                "  processed {:,} / {:,} promoters".format(
                    i + 1,
                    len(promoters)
                )
            )

    bw27.close()
    bw4.close()

    return {
        "H3K27me3":
            h3k27,

        "H3K4me3":
            h3k4,

        "log2ratio":
            ratio
    }


# ============================================================
# HOTSPOT-LEVEL STATISTICS
# ============================================================

def calculate_hotspot_metric(
    region_promoter_indices,
    promoter_values
):
    """
    For each hotspot:
        median signal across promoters.

    Final statistic:
        median across hotspots.

    This preserves equal weighting of hotspots rather than
    allowing promoter-dense hotspots to dominate.
    """

    hotspot_values = []

    for indices in region_promoter_indices:

        values = promoter_values[
            indices
        ]

        values = values[
            np.isfinite(
                values
            )
        ]

        if len(values) == 0:

            hotspot_values.append(
                np.nan
            )

        else:

            hotspot_values.append(
                np.median(
                    values
                )
            )

    hotspot_values = np.array(
        hotspot_values,
        dtype=float
    )

    if np.all(
        np.isnan(
            hotspot_values
        )
    ):

        overall = np.nan

    else:

        overall = np.nanmedian(
            hotspot_values
        )

    return (
        overall,
        hotspot_values
    )


# ============================================================
# EMPIRICAL P VALUES
# ============================================================

def empirical_p_greater(
    observed,
    random_values
):

    values = random_values[
        np.isfinite(
            random_values
        )
    ]

    return (
        np.sum(
            values >= observed
        )
        +
        1
    ) / float(
        len(values)
        +
        1
    )


def empirical_p_less(
    observed,
    random_values
):

    values = random_values[
        np.isfinite(
            random_values
        )
    ]

    return (
        np.sum(
            values <= observed
        )
        +
        1
    ) / float(
        len(values)
        +
        1
    )


def empirical_p_two_sided(
    observed,
    random_values
):

    p_low = empirical_p_less(
        observed,
        random_values
    )

    p_high = empirical_p_greater(
        observed,
        random_values
    )

    return min(
        1.0,
        2.0
        *
        min(
            p_low,
            p_high
        )
    )


# ============================================================
# ANALYZE ONE DATASET
# ============================================================

def analyze_dataset(
    dataset_name,
    dataset_paths,
    promoters,
    observed_region_indices,
    random_region_index_sets
):

    print(
        "\n"
        "============================================================"
    )

    print(
        "DATASET: {}".format(
            dataset_name
        )
    )

    print(
        "============================================================"
    )

    # --------------------------------------------------------
    # Promoter signals
    # --------------------------------------------------------

    signals = calculate_promoter_signals(
        promoters,
        dataset_paths[
            "H3K27me3"
        ],
        dataset_paths[
            "H3K4me3"
        ]
    )

    # --------------------------------------------------------
    # Coverage
    # --------------------------------------------------------

    coverage27 = np.mean(
        np.isfinite(
            signals[
                "H3K27me3"
            ]
        )
    )

    coverage4 = np.mean(
        np.isfinite(
            signals[
                "H3K4me3"
            ]
        )
    )

    coverage_ratio = np.mean(
        np.isfinite(
            signals[
                "log2ratio"
            ]
        )
    )

    print(
        "Promoter signal coverage:"
    )

    print(
        "  H3K27me3: {:.3f}".format(
            coverage27
        )
    )

    print(
        "  H3K4me3 : {:.3f}".format(
            coverage4
        )
    )

    print(
        "  Ratio   : {:.3f}".format(
            coverage_ratio
        )
    )

    # --------------------------------------------------------
    # Observed
    # --------------------------------------------------------

    observed27, observed27_hotspots = (
        calculate_hotspot_metric(
            observed_region_indices,
            signals[
                "H3K27me3"
            ]
        )
    )

    observed4, observed4_hotspots = (
        calculate_hotspot_metric(
            observed_region_indices,
            signals[
                "H3K4me3"
            ]
        )
    )

    observed_ratio, observed_ratio_hotspots = (
        calculate_hotspot_metric(
            observed_region_indices,
            signals[
                "log2ratio"
            ]
        )
    )

    # --------------------------------------------------------
    # Permutations
    # --------------------------------------------------------

    random27 = np.full(
        N_PERMUTATIONS,
        np.nan
    )

    random4 = np.full(
        N_PERMUTATIONS,
        np.nan
    )

    random_ratio = np.full(
        N_PERMUTATIONS,
        np.nan
    )

    print(
        "\nEvaluating {:,} permutations...".format(
            N_PERMUTATIONS
        )
    )

    for i in range(
        N_PERMUTATIONS
    ):

        region_indices = (
            random_region_index_sets[
                i
            ]
        )

        random27[i] = (
            calculate_hotspot_metric(
                region_indices,
                signals[
                    "H3K27me3"
                ]
            )[0]
        )

        random4[i] = (
            calculate_hotspot_metric(
                region_indices,
                signals[
                    "H3K4me3"
                ]
            )[0]
        )

        random_ratio[i] = (
            calculate_hotspot_metric(
                region_indices,
                signals[
                    "log2ratio"
                ]
            )[0]
        )

        if (
            (i + 1) % 1000
            ==
            0
        ):

            print(
                "  {:,} / {:,}".format(
                    i + 1,
                    N_PERMUTATIONS
                )
            )

    # --------------------------------------------------------
    # Dataset output directory
    # --------------------------------------------------------

    dataset_outdir = os.path.join(
        OUTDIR,
        dataset_name
    )

    os.makedirs(
        dataset_outdir,
        exist_ok=True
    )

    # --------------------------------------------------------
    # Save promoter signals
    # --------------------------------------------------------

    promoter_output = promoters[
        [
            "promoter_index",
            "symbol",
            "chrom",
            "tss_start"
        ]
    ].copy()

    promoter_output[
        "H3K27me3"
    ] = signals[
        "H3K27me3"
    ]

    promoter_output[
        "H3K4me3"
    ] = signals[
        "H3K4me3"
    ]

    promoter_output[
        "log2_H3K27me3_H3K4me3"
    ] = signals[
        "log2ratio"
    ]

    promoter_output.to_csv(
        os.path.join(
            dataset_outdir,
            "promoter_histone_signals.tsv"
        ),
        sep="\t",
        index=False
    )

    # --------------------------------------------------------
    # Save observed hotspot values
    # --------------------------------------------------------

    hotspot_output = pd.DataFrame({

        "hotspot_id":
            np.arange(
                1,
                len(
                    observed27_hotspots
                )
                +
                1
            ),

        "median_H3K27me3":
            observed27_hotspots,

        "median_H3K4me3":
            observed4_hotspots,

        "median_log2ratio":
            observed_ratio_hotspots
    })

    hotspot_output.to_csv(
        os.path.join(
            dataset_outdir,
            "observed_hotspot_histone_values.tsv"
        ),
        sep="\t",
        index=False
    )

    # --------------------------------------------------------
    # Save permutation distributions
    # --------------------------------------------------------

    permutation_df = pd.DataFrame({

        "permutation":
            np.arange(
                1,
                N_PERMUTATIONS
                +
                1
            ),

        "median_hotspot_H3K27me3":
            random27,

        "median_hotspot_H3K4me3":
            random4,

        "median_hotspot_H3K27me3_H3K4me3_log2ratio":
            random_ratio
    })

    permutation_df.to_csv(
        os.path.join(
            dataset_outdir,
            "histone_permutation_distributions.tsv"
        ),
        sep="\t",
        index=False
    )

    # --------------------------------------------------------
    # Summary
    # --------------------------------------------------------

    metrics = [

        (
            "median_hotspot_H3K27me3",
            observed27,
            random27,
            "greater"
        ),

        (
            "median_hotspot_H3K4me3",
            observed4,
            random4,
            "two-sided"
        ),

        (
            "median_hotspot_H3K27me3_H3K4me3_log2ratio",
            observed_ratio,
            random_ratio,
            "greater"
        )
    ]

    summary_rows = []

    for (
        metric,
        observed,
        random_values,
        alternative
    ) in metrics:

        random_values = random_values[
            np.isfinite(
                random_values
            )
        ]

        if alternative == "greater":

            p = empirical_p_greater(
                observed,
                random_values
            )

        else:

            p = empirical_p_two_sided(
                observed,
                random_values
            )

        random_mean = np.mean(
            random_values
        )

        random_sd = np.std(
            random_values,
            ddof=1
        )

        permutation_z = (
            (
                observed
                -
                random_mean
            )
            /
            random_sd
            if random_sd > 0
            else np.nan
        )

        row = {

            "dataset":
                dataset_name,

            "metric":
                metric,

            "alternative":
                alternative,

            "observed":
                observed,

            "random_mean":
                random_mean,

            "random_median":
                np.median(
                    random_values
                ),

            "random_2.5_percentile":
                np.percentile(
                    random_values,
                    2.5
                ),

            "random_97.5_percentile":
                np.percentile(
                    random_values,
                    97.5
                ),

            "permutation_z":
                permutation_z,

            "empirical_p":
                p
        }

        # Additional directional test for H3K4me3
        if (
            metric
            ==
            "median_hotspot_H3K4me3"
        ):

            row[
                "empirical_p_depletion"
            ] = empirical_p_less(
                observed,
                random_values
            )

        else:

            row[
                "empirical_p_depletion"
            ] = np.nan

        summary_rows.append(
            row
        )

    summary_df = pd.DataFrame(
        summary_rows
    )

    summary_df.to_csv(
        os.path.join(
            dataset_outdir,
            "histone_permutation_summary.tsv"
        ),
        sep="\t",
        index=False
    )

    return (
        summary_df,
        {
            "dataset":
                dataset_name,

            "H3K27me3_coverage":
                coverage27,

            "H3K4me3_coverage":
                coverage4,

            "ratio_coverage":
                coverage_ratio
        }
    )


# ============================================================
# MAIN
# ============================================================

def main():

    os.makedirs(
        OUTDIR,
        exist_ok=True
    )

    validate_bigwig_files()

    print(
        "\nLoading hotspots/promoters..."
    )

    hotspots = load_hotspots(
        HOTSPOT_FILE
    )

    promoters = load_promoters(
        PROMOTER_FILE
    )

    chromsizes = load_chromsizes(
        CHROMSIZES_FILE
    )

    print(
        "Hotspots: {}".format(
            len(hotspots)
        )
    )

    print(
        "Promoters: {:,}".format(
            len(promoters)
        )
    )

    promoter_lookup = (
        build_promoter_lookup(
            promoters
        )
    )

    hotspots = (
        validate_hotspot_promoter_counts(
            hotspots,
            promoter_lookup
        )
    )

    # ========================================================
    # Observed promoter indices per hotspot
    # ========================================================

    observed_region_indices = (
        get_region_promoter_indices(
            hotspots,
            promoter_lookup
        )
    )

    # ========================================================
    # Generate same permutations ONCE
    # ========================================================

    cache = (
        precompute_valid_start_ranges(
            hotspots,
            chromsizes,
            promoter_lookup
        )
    )

    rng = np.random.RandomState(
        RANDOM_SEED
    )

    random_region_index_sets = []

    print(
        "\nGenerating {:,} common matched "
        "random hotspot sets...".format(
            N_PERMUTATIONS
        )
    )

    for i in range(
        N_PERMUTATIONS
    ):

        random_regions = (
            generate_random_regions(
                hotspots,
                cache,
                rng
            )
        )

        region_indices = (
            get_region_promoter_indices(
                random_regions,
                promoter_lookup
            )
        )

        random_region_index_sets.append(
            region_indices
        )

        if (
            (i + 1) % 1000
            ==
            0
        ):

            print(
                "Generated {:,} / {:,}".format(
                    i + 1,
                    N_PERMUTATIONS
                )
            )

    # ========================================================
    # Analyze datasets
    # ========================================================

    combined_summaries = []
    coverage_rows = []

    for dataset_name in DATASETS:

        summary_df, coverage = (
            analyze_dataset(
                dataset_name,
                DATASETS[
                    dataset_name
                ],
                promoters,
                observed_region_indices,
                random_region_index_sets
            )
        )

        combined_summaries.append(
            summary_df
        )

        coverage_rows.append(
            coverage
        )

    # ========================================================
    # Combined summary
    # ========================================================

    combined = pd.concat(
        combined_summaries,
        ignore_index=True
    )

    combined.to_csv(
        os.path.join(
            OUTDIR,
            "ALL_TISSUES_histone_permutation_summary.tsv"
        ),
        sep="\t",
        index=False
    )

    coverage_df = pd.DataFrame(
        coverage_rows
    )

    coverage_df.to_csv(
        os.path.join(
            OUTDIR,
            "ALL_TISSUES_histone_signal_coverage.tsv"
        ),
        sep="\t",
        index=False
    )

    # ========================================================
    # Wide matrices for plotting
    # ========================================================

    z_matrix = combined.pivot(
        index="metric",
        columns="dataset",
        values="permutation_z"
    )

    z_matrix.to_csv(
        os.path.join(
            OUTDIR,
            "histone_permutation_z_matrix.tsv"
        ),
        sep="\t"
    )

    p_matrix = combined.pivot(
        index="metric",
        columns="dataset",
        values="empirical_p"
    )

    p_matrix.to_csv(
        os.path.join(
            OUTDIR,
            "histone_empirical_p_matrix.tsv"
        ),
        sep="\t"
    )

    # ========================================================
    # Print summary
    # ========================================================

    print(
        "\n"
        "============================================================"
    )

    print(
        "ALL-TISSUE HISTONE RESULTS"
    )

    print(
        "============================================================"
    )

    print(
        combined[
            [
                "dataset",
                "metric",
                "observed",
                "random_median",
                "permutation_z",
                "empirical_p",
                "empirical_p_depletion"
            ]
        ].to_string(
            index=False
        )
    )

    print(
        "\nResults written to:"
    )

    print(
        OUTDIR
    )



def _configure_cli():
    global HOTSPOT_FILE, PROMOTER_FILE, CHROMSIZES_FILE, OUTDIR, N_PERMUTATIONS, RANDOM_SEED
    import argparse
    parser = argparse.ArgumentParser(description=__doc__ or "Promoter methylation hotspot analysis")
    parser.add_argument("--hotspots", default=HOTSPOT_FILE)
    parser.add_argument("--promoters", default=PROMOTER_FILE)
    parser.add_argument("--chrom-sizes", default=CHROMSIZES_FILE)
    parser.add_argument("--outdir", default=OUTDIR)
    parser.add_argument("--permutations", type=int, default=N_PERMUTATIONS)
    parser.add_argument("--seed", type=int, default=RANDOM_SEED)
    parser.add_argument("--reference-root", default=None,
                        help="Directory corresponding to original promoter_clustering directory; no files copied")
    args = parser.parse_args()
    HOTSPOT_FILE = args.hotspots
    PROMOTER_FILE = args.promoters
    CHROMSIZES_FILE = args.chrom_sizes
    OUTDIR = args.outdir
    N_PERMUTATIONS = args.permutations
    RANDOM_SEED = args.seed
    if args.reference_root:
        _replace_reference_root(args.reference_root)

def _replace_reference_root(new_root):
    # Only rewrite original GSC reference paths, not outputs or user-selected input paths.
    import os
    original = "/projects/clininfo/vcsizmok/long_POG/aDMR/aDMR_analysis/promoter_clustering"
    for key, value in list(globals().items()):
        if key in ("OUTDIR", "HOTSPOT_FILE", "PROMOTER_FILE", "CHROMSIZES_FILE"):
            continue
        if key.endswith(("_DIR", "_BASE", "_FILE")) and isinstance(value, str) and value.startswith(original + os.sep):
            globals()[key] = os.path.join(new_root, os.path.relpath(value, original))
        if key in ("EPIGENOMES", "CTCF_DATASETS", "DATASETS", "TAD_MAPS") and isinstance(value, dict):
            def _remap(v):
                if isinstance(v, str) and v.startswith(original + os.sep):
                    return os.path.join(new_root, os.path.relpath(v, original))
                if isinstance(v, dict):
                    return {k: _remap(x) for k, x in v.items()}
                return v
            globals()[key] = _remap(value)

if __name__ == "__main__":
    _configure_cli()

    main()