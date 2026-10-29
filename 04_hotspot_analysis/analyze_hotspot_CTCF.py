#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Multi-context CTCF permutation analysis.

For each CTCF reference context, compare observed methylation hotspots
against the SAME 10,000 promoter-matched random hotspot sets.

Matching criteria:
    - chromosome
    - hotspot length
    - number of evaluable promoters

CTCF metrics:
    1. fraction of hotspots overlapping >=1 CTCF peak
    2. total overlapping CTCF peaks
    3. median CTCF peaks per 100 kb
    4. median distance from hotspot midpoint to nearest CTCF peak center
    5. median mean distance of hotspot edges to nearest CTCF peak center
    6. median minimum edge distance to nearest CTCF peak center

Compatible with Python 3.7.
"""

import os
import numpy as np
import pandas as pd


# ============================================================
# USER SETTINGS
# ============================================================

# Repository inputs and shared GSC annotations (override via CLI if needed).
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
HOTSPOT_FILE = os.path.join(SCRIPT_DIR, "data", "methylated_promoter_hotspots.bed")
PROMOTER_FILE = os.path.join(SCRIPT_DIR, "data", "input_promoters_coordinates_methylation.tsv")
CHROMSIZES_FILE = "/gsc/www/bcgsc.ca/downloads/promoter_methylation_nanopore/annotations/genome/hg38.chrom.sizes"

CTCF_BASE = "/projects/clininfo/vcsizmok/long_POG/aDMR/aDMR_analysis/promoter_clustering/CTCF_binding_sites"

OUTDIR = os.path.join(SCRIPT_DIR, "results", "CTCF_analysis", "all_contexts")

N_PERMUTATIONS = 10000
RANDOM_SEED = 12345


# ============================================================
# CTCF DATASETS
# ============================================================

CTCF_DATASETS = {

    "MCF7": os.path.join(
        CTCF_BASE,
        "MCF-7",
        "ENCFF278FNP.bed.gz"
    ),

    "Mammary_fibroblast": os.path.join(
        CTCF_BASE,
        "mammary_fibroblast",
        "ENCFF498KKL.bed.gz"
    ),

    "Colon": os.path.join(
        CTCF_BASE,
        "colon",
        "ENCFF115EJR.bed.gz"
    ),

    "Liver": os.path.join(
        CTCF_BASE,
        "liver",
        "ENCFF046XEZ.bed.gz"
    ),

    "Pancreas": os.path.join(
        CTCF_BASE,
        "pancreas",
        "ENCFF225IGV.bed.gz"
    ),

    "Lung": os.path.join(
        CTCF_BASE,
        "lung",
        "ENCFF629UDP.bed.gz"
    ),

    "Astrocytes": os.path.join(
        CTCF_BASE,
        "astrocytes",
        "ENCFF450EVR.bed.gz"
    )
}


DATASET_ORDER = [
    "MCF7",
    "Mammary_fibroblast",
    "Colon",
    "Liver",
    "Pancreas",
    "Lung",
    "Astrocytes"
]


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

    # Reproduce promoter coordinate used in hotspot discovery.
    df["midpoint"] = (
        (
            df["promoter_start"]
            +
            df["promoter_end"]
        )
        // 2
    ).astype(int)

    df = df.reset_index(
        drop=True
    )

    df["promoter_index"] = df.index

    return df


# ============================================================
# CHROMOSOME SIZES
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

    df["size"] = pd.to_numeric(
        df["size"],
        errors="raise"
    ).astype(int)

    return dict(
        zip(
            df["chrom"],
            df["size"]
        )
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

        lookup[chrom] = {

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


def count_promoters_in_region(
    chrom,
    start,
    end,
    promoter_lookup
):

    return len(
        get_promoter_indices_in_region(
            chrom,
            start,
            end,
            promoter_lookup
        )
    )


# ============================================================
# VALIDATE HOTSPOT PROMOTER COUNTS
# ============================================================

def validate_hotspot_promoter_counts(
    hotspots,
    promoter_lookup
):

    counts = []

    for _, row in hotspots.iterrows():

        counts.append(
            count_promoters_in_region(
                row["chrom"],
                int(row["start"]),
                int(row["end"]),
                promoter_lookup
            )
        )

    hotspots = hotspots.copy()

    hotspots[
        "promoter_count_recalculated"
    ] = counts

    matches = (
        hotspots["n_promoters"]
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
            "Promoter-count validation failed."
        )

    return hotspots


# ============================================================
# LOAD CTCF PEAKS
# ============================================================

def load_ctcf_peaks(filename):

    df = pd.read_csv(
        filename,
        sep="\t",
        compression="gzip",
        header=None,
        comment="#"
    )

    if df.shape[1] < 3:

        raise ValueError(
            "CTCF file contains fewer than 3 columns: "
            + filename
        )

    df = df.iloc[
        :,
        :3
    ].copy()

    df.columns = [
        "chrom",
        "start",
        "end"
    ]

    df["start"] = pd.to_numeric(
        df["start"],
        errors="raise"
    ).astype(int)

    df["end"] = pd.to_numeric(
        df["end"],
        errors="raise"
    ).astype(int)

    df["center"] = (
        (
            df["start"]
            +
            df["end"]
        )
        // 2
    ).astype(int)

    return df


# ============================================================
# CTCF LOOKUP
# ============================================================

def build_ctcf_lookup(ctcf):

    lookup = {}

    for chrom, group in ctcf.groupby(
        "chrom"
    ):

        group = group.sort_values(
            "start"
        )

        lookup[chrom] = {

            "starts":
                group[
                    "start"
                ]
                .astype(int)
                .to_numpy(),

            "ends":
                group[
                    "end"
                ]
                .astype(int)
                .to_numpy(),

            "centers":
                np.sort(
                    group[
                        "center"
                    ]
                    .astype(int)
                    .to_numpy()
                )
        }

    return lookup


# ============================================================
# CTCF REGION METRICS
# ============================================================

def count_ctcf_peaks_in_region(
    chrom,
    start,
    end,
    ctcf_lookup
):

    if chrom not in ctcf_lookup:

        return 0

    starts = ctcf_lookup[
        chrom
    ]["starts"]

    ends = ctcf_lookup[
        chrom
    ]["ends"]

    right = np.searchsorted(
        starts,
        end,
        side="left"
    )

    if right == 0:

        return 0

    count = np.sum(
        ends[:right]
        >
        start
    )

    return int(
        count
    )


def nearest_ctcf_center_distance(
    chrom,
    position,
    ctcf_lookup
):

    if chrom not in ctcf_lookup:

        return np.nan

    centers = ctcf_lookup[
        chrom
    ]["centers"]

    if len(centers) == 0:

        return np.nan

    idx = np.searchsorted(
        centers,
        position
    )

    distances = []

    if idx < len(centers):

        distances.append(
            abs(
                int(
                    centers[idx]
                )
                -
                position
            )
        )

    if idx > 0:

        distances.append(
            abs(
                int(
                    centers[
                        idx - 1
                    ]
                )
                -
                position
            )
        )

    return float(
        min(
            distances
        )
    )


# ============================================================
# ANALYZE ONE REGION
# ============================================================

def analyze_region(
    row,
    ctcf_lookup
):

    chrom = row[
        "chrom"
    ]

    start = int(
        row[
            "start"
        ]
    )

    end = int(
        row[
            "end"
        ]
    )

    length = (
        end
        -
        start
    )

    n_ctcf = count_ctcf_peaks_in_region(
        chrom,
        start,
        end,
        ctcf_lookup
    )

    peaks_per_100kb = (
        n_ctcf
        /
        float(length)
        *
        100000.0
    )

    midpoint = int(
        (
            start
            +
            end
        )
        // 2
    )

    midpoint_distance = nearest_ctcf_center_distance(
        chrom,
        midpoint,
        ctcf_lookup
    )

    start_distance = nearest_ctcf_center_distance(
        chrom,
        start,
        ctcf_lookup
    )

    end_distance = nearest_ctcf_center_distance(
        chrom,
        end,
        ctcf_lookup
    )

    edge_values = np.array(
        [
            start_distance,
            end_distance
        ],
        dtype=float
    )

    if np.all(
        np.isnan(
            edge_values
        )
    ):

        mean_edge_distance = np.nan
        min_edge_distance = np.nan

    else:

        mean_edge_distance = np.nanmean(
            edge_values
        )

        min_edge_distance = np.nanmin(
            edge_values
        )

    return {

        "CTCF_overlap":
            n_ctcf > 0,

        "CTCF_peak_count":
            n_ctcf,

        "CTCF_peaks_per_100kb":
            peaks_per_100kb,

        "CTCF_midpoint_distance_bp":
            midpoint_distance,

        "CTCF_mean_edge_distance_bp":
            mean_edge_distance,

        "CTCF_minimum_edge_distance_bp":
            min_edge_distance
    }


def analyze_regions(
    regions,
    ctcf_lookup
):

    rows = []

    for _, row in regions.iterrows():

        metrics = analyze_region(
            row,
            ctcf_lookup
        )

        output = {

            "chrom":
                row["chrom"],

            "start":
                int(row["start"]),

            "end":
                int(row["end"]),

            "length":
                int(
                    row["end"]
                    -
                    row["start"]
                )
        }

        for key, value in metrics.items():

            output[
                key
            ] = value

        if "hotspot_id" in row.index:

            output[
                "hotspot_id"
            ] = row[
                "hotspot_id"
            ]

        rows.append(
            output
        )

    return pd.DataFrame(
        rows
    )


# ============================================================
# COHORT STATISTICS
# ============================================================

def calculate_cohort_statistics(
    region_metrics
):

    return {

        "fraction_CTCF_overlap":
            region_metrics[
                "CTCF_overlap"
            ].mean(),

        "total_CTCF_peaks":
            region_metrics[
                "CTCF_peak_count"
            ].sum(),

        "median_CTCF_peaks_per_100kb":
            region_metrics[
                "CTCF_peaks_per_100kb"
            ].median(),

        "median_CTCF_midpoint_distance_bp":
            region_metrics[
                "CTCF_midpoint_distance_bp"
            ].median(),

        "median_CTCF_mean_edge_distance_bp":
            region_metrics[
                "CTCF_mean_edge_distance_bp"
            ].median(),

        "median_CTCF_minimum_edge_distance_bp":
            region_metrics[
                "CTCF_minimum_edge_distance_bp"
            ].median()
    }


# ============================================================
# PROMOTER-MATCHED RANDOMIZATION
# ============================================================

def get_valid_start_ranges(
    chrom,
    length,
    target_n_promoters,
    chromsizes,
    promoter_lookup
):

    if chrom not in chromsizes:

        raise ValueError(
            "Missing chromosome size: "
            + chrom
        )

    if chrom not in promoter_lookup:

        raise ValueError(
            "Missing promoter coordinates: "
            + chrom
        )

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

    valid_ranges = []

    if k > n:

        return valid_ranges

    # --------------------------------------------------------
    # ZERO PROMOTERS
    # --------------------------------------------------------

    if k == 0:

        previous_position = -1

        for position in positions:

            position = int(
                position
            )

            lower = (
                previous_position
                +
                1
            )

            upper = (
                position
                -
                length
                -
                1
            )

            lower = max(
                0,
                lower
            )

            upper = min(
                max_start,
                upper
            )

            if lower <= upper:

                valid_ranges.append(
                    (
                        int(lower),
                        int(upper)
                    )
                )

            previous_position = (
                position
            )

        lower = (
            previous_position
            +
            1
        )

        upper = (
            max_start
        )

        if lower <= upper:

            valid_ranges.append(
                (
                    int(lower),
                    int(upper)
                )
            )

        return valid_ranges

    # --------------------------------------------------------
    # ONE OR MORE PROMOTERS
    # --------------------------------------------------------

    for i in range(
        0,
        n - k + 1
    ):

        j = (
            i
            +
            k
            -
            1
        )

        first_promoter = int(
            positions[
                i
            ]
        )

        last_promoter = int(
            positions[
                j
            ]
        )

        if (
            last_promoter
            -
            first_promoter
            >
            length
        ):

            continue

        lower = max(
            0,
            last_promoter
            -
            length
        )

        upper = min(
            max_start,
            first_promoter
        )

        if i > 0:

            previous_promoter = int(
                positions[
                    i - 1
                ]
            )

            lower = max(
                lower,
                previous_promoter
                +
                1
            )

        if j + 1 < n:

            next_promoter = int(
                positions[
                    j + 1
                ]
            )

            upper = min(
                upper,
                next_promoter
                -
                length
                -
                1
            )

        if lower <= upper:

            valid_ranges.append(
                (
                    int(lower),
                    int(upper)
                )
            )

    return valid_ranges


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
        "\nPrecomputing valid promoter-matched intervals..."
    )

    for _, row in combinations.iterrows():

        chrom = row[
            "chrom"
        ]

        length = int(
            row[
                "length"
            ]
        )

        n_promoters = int(
            row[
                "n_promoters"
            ]
        )

        key = (
            chrom,
            length,
            n_promoters
        )

        ranges = get_valid_start_ranges(
            chrom,
            length,
            n_promoters,
            chromsizes,
            promoter_lookup
        )

        if len(ranges) == 0:

            raise RuntimeError(
                "No matched regions for "
                "{} length={} n_promoters={}".format(
                    chrom,
                    length,
                    n_promoters
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


def sample_one_random_region(
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

    random_index = int(
        rng.randint(
            0,
            item[
                "total"
            ]
        )
    )

    range_index = int(
        np.searchsorted(
            item[
                "cumulative"
            ],
            random_index,
            side="right"
        )
    )

    previous = 0

    if range_index > 0:

        previous = int(
            item[
                "cumulative"
            ][
                range_index - 1
            ]
        )

    lower, upper = (
        item[
            "ranges"
        ][
            range_index
        ]
    )

    start = int(
        lower
        +
        (
            random_index
            -
            previous
        )
    )

    return {

        "chrom":
            chrom,

        "start":
            start,

        "end":
            start
            +
            length,

        "length":
            length,

        "n_promoters":
            n_promoters
    }


def generate_random_regions(
    hotspots,
    cache,
    rng
):

    rows = []

    for _, row in hotspots.iterrows():

        rows.append(
            sample_one_random_region(
                row["chrom"],
                int(row["length"]),
                int(row["n_promoters"]),
                cache,
                rng
            )
        )

    return pd.DataFrame(
        rows
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
# BUILD SUMMARY
# ============================================================

def build_summary(
    dataset_name,
    observed_stats,
    distributions
):

    alternatives = {

        "fraction_CTCF_overlap":
            "greater",

        "total_CTCF_peaks":
            "greater",

        "median_CTCF_peaks_per_100kb":
            "greater",

        "median_CTCF_midpoint_distance_bp":
            "less",

        "median_CTCF_mean_edge_distance_bp":
            "less",

        "median_CTCF_minimum_edge_distance_bp":
            "less"
    }

    rows = []

    for metric in observed_stats:

        observed = float(
            observed_stats[
                metric
            ]
        )

        random_values = np.asarray(
            distributions[
                metric
            ],
            dtype=float
        )

        random_values = random_values[
            np.isfinite(
                random_values
            )
        ]

        random_mean = np.mean(
            random_values
        )

        random_median = np.median(
            random_values
        )

        random_sd = np.std(
            random_values,
            ddof=1
        )

        if random_sd > 0:

            permutation_z = (
                observed
                -
                random_mean
            ) / random_sd

        else:

            permutation_z = np.nan

        if random_median != 0:

            ratio = (
                observed
                /
                random_median
            )

        else:

            ratio = np.nan

        alternative = alternatives[
            metric
        ]

        if alternative == "greater":

            directional_p = empirical_p_greater(
                observed,
                random_values
            )

        else:

            directional_p = empirical_p_less(
                observed,
                random_values
            )

        two_sided_p = empirical_p_two_sided(
            observed,
            random_values
        )

        rows.append({

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
                random_median,

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

            "observed_to_random_median_ratio":
                ratio,

            "permutation_z":
                permutation_z,

            "empirical_p":
                directional_p,

            "empirical_two_sided_p":
                two_sided_p
        })

    return pd.DataFrame(
        rows
    )


# ============================================================
# MAIN
# ============================================================

def main():

    os.makedirs(
        OUTDIR,
        exist_ok=True
    )

    # --------------------------------------------------------
    # Check input files
    # --------------------------------------------------------

    for filename in [
        HOTSPOT_FILE,
        PROMOTER_FILE,
        CHROMSIZES_FILE
    ]:

        if not os.path.exists(
            filename
        ):

            raise IOError(
                "Missing input file:\n"
                + filename
            )

    print(
        "Checking CTCF files..."
    )

    for dataset in DATASET_ORDER:

        filename = CTCF_DATASETS[
            dataset
        ]

        if not os.path.exists(
            filename
        ):

            raise IOError(
                "Missing CTCF file for {}:\n{}".format(
                    dataset,
                    filename
                )
            )

        print(
            "{:<22s} {}".format(
                dataset,
                filename
            )
        )

    # --------------------------------------------------------
    # Load hotspot/promoter inputs
    # --------------------------------------------------------

    print(
        "\nLoading hotspot and promoter data..."
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

    promoter_lookup = build_promoter_lookup(
        promoters
    )

    print(
        "Hotspots: {}".format(
            len(hotspots)
        )
    )

    print(
        "Evaluable promoters: {:,}".format(
            len(promoters)
        )
    )

    hotspots = validate_hotspot_promoter_counts(
        hotspots,
        promoter_lookup
    )

    # --------------------------------------------------------
    # Precompute random interval possibilities
    # --------------------------------------------------------

    cache = precompute_valid_start_ranges(
        hotspots,
        chromsizes,
        promoter_lookup
    )

    # --------------------------------------------------------
    # Generate SAME random hotspot sets once
    # --------------------------------------------------------

    rng = np.random.RandomState(
        RANDOM_SEED
    )

    random_region_sets = []

    print(
        "\nGenerating {:,} common random hotspot sets...".format(
            N_PERMUTATIONS
        )
    )

    for i in range(
        N_PERMUTATIONS
    ):

        random_regions = generate_random_regions(
            hotspots,
            cache,
            rng
        )

        random_region_sets.append(
            random_regions
        )

        if (
            (i + 1)
            %
            1000
            ==
            0
        ):

            print(
                "Generated {:,} / {:,}".format(
                    i + 1,
                    N_PERMUTATIONS
                )
            )

    # --------------------------------------------------------
    # Analyze each CTCF context
    # --------------------------------------------------------

    all_summaries = []
    peak_count_rows = []

    for dataset in DATASET_ORDER:

        print(
            "\n"
            "============================================================"
        )

        print(
            "CTCF DATASET: {}".format(
                dataset
            )
        )

        print(
            "============================================================"
        )

        filename = CTCF_DATASETS[
            dataset
        ]

        ctcf = load_ctcf_peaks(
            filename
        )

        print(
            "CTCF peaks: {:,}".format(
                len(ctcf)
            )
        )

        peak_count_rows.append({

            "dataset":
                dataset,

            "file":
                filename,

            "n_CTCF_peaks":
                len(ctcf)
        })

        ctcf_lookup = build_ctcf_lookup(
            ctcf
        )

        # ----------------------------------------------------
        # Observed
        # ----------------------------------------------------

        observed_metrics = analyze_regions(
            hotspots,
            ctcf_lookup
        )

        observed_stats = calculate_cohort_statistics(
            observed_metrics
        )

        # ----------------------------------------------------
        # Random distributions
        # ----------------------------------------------------

        distributions = {

            metric:
                np.full(
                    N_PERMUTATIONS,
                    np.nan
                )

            for metric in observed_stats
        }

        print(
            "Evaluating common random hotspot sets..."
        )

        for i in range(
            N_PERMUTATIONS
        ):

            random_metrics = analyze_regions(
                random_region_sets[
                    i
                ],
                ctcf_lookup
            )

            random_stats = calculate_cohort_statistics(
                random_metrics
            )

            for metric in observed_stats:

                distributions[
                    metric
                ][
                    i
                ] = random_stats[
                    metric
                ]

            if (
                (i + 1)
                %
                1000
                ==
                0
            ):

                print(
                    "  {:,} / {:,}".format(
                        i + 1,
                        N_PERMUTATIONS
                    )
                )

        # ----------------------------------------------------
        # Output directory
        # ----------------------------------------------------

        dataset_outdir = os.path.join(
            OUTDIR,
            dataset
        )

        os.makedirs(
            dataset_outdir,
            exist_ok=True
        )

        # ----------------------------------------------------
        # Save observed hotspot metrics
        # ----------------------------------------------------

        observed_metrics.to_csv(
            os.path.join(
                dataset_outdir,
                "observed_hotspot_CTCF_metrics.tsv"
            ),
            sep="\t",
            index=False
        )

        # ----------------------------------------------------
        # Save permutation distributions
        # ----------------------------------------------------

        permutation_df = pd.DataFrame({

            "permutation":
                np.arange(
                    1,
                    N_PERMUTATIONS
                    +
                    1
                )
        })

        for metric in distributions:

            permutation_df[
                metric
            ] = distributions[
                metric
            ]

        permutation_df.to_csv(
            os.path.join(
                dataset_outdir,
                "CTCF_permutation_distributions.tsv"
            ),
            sep="\t",
            index=False
        )

        # ----------------------------------------------------
        # Summary
        # ----------------------------------------------------

        summary = build_summary(
            dataset,
            observed_stats,
            distributions
        )

        summary.to_csv(
            os.path.join(
                dataset_outdir,
                "CTCF_permutation_summary.tsv"
            ),
            sep="\t",
            index=False
        )

        all_summaries.append(
            summary
        )

        print(
            summary[
                [
                    "metric",
                    "observed",
                    "random_median",
                    "permutation_z",
                    "empirical_p",
                    "empirical_two_sided_p"
                ]
            ].to_string(
                index=False
            )
        )

    # ========================================================
    # COMBINED SUMMARY
    # ========================================================

    combined = pd.concat(
        all_summaries,
        ignore_index=True
    )

    combined.to_csv(
        os.path.join(
            OUTDIR,
            "ALL_CONTEXTS_CTCF_permutation_summary.tsv"
        ),
        sep="\t",
        index=False
    )

    # ========================================================
    # CTCF PEAK-COUNT METADATA
    # ========================================================

    peak_count_df = pd.DataFrame(
        peak_count_rows
    )

    peak_count_df.to_csv(
        os.path.join(
            OUTDIR,
            "CTCF_dataset_peak_counts.tsv"
        ),
        sep="\t",
        index=False
    )

    # ========================================================
    # Z-SCORE MATRIX
    # ========================================================

    z_matrix = combined.pivot(
        index="metric",
        columns="dataset",
        values="permutation_z"
    )

    z_matrix = z_matrix.reindex(
        columns=DATASET_ORDER
    )

    z_matrix.to_csv(
        os.path.join(
            OUTDIR,
            "CTCF_permutation_z_matrix.tsv"
        ),
        sep="\t"
    )

    # ========================================================
    # OBSERVED / NULL MATRIX
    # ========================================================

    ratio_matrix = combined.pivot(
        index="metric",
        columns="dataset",
        values="observed_to_random_median_ratio"
    )

    ratio_matrix = ratio_matrix.reindex(
        columns=DATASET_ORDER
    )

    ratio_matrix.to_csv(
        os.path.join(
            OUTDIR,
            "CTCF_observed_to_null_ratio_matrix.tsv"
        ),
        sep="\t"
    )

    # ========================================================
    # P-VALUE MATRIX
    # ========================================================

    p_matrix = combined.pivot(
        index="metric",
        columns="dataset",
        values="empirical_p"
    )

    p_matrix = p_matrix.reindex(
        columns=DATASET_ORDER
    )

    p_matrix.to_csv(
        os.path.join(
            OUTDIR,
            "CTCF_empirical_p_matrix.tsv"
        ),
        sep="\t"
    )

    # ========================================================
    # PRINT FINAL COMBINED SUMMARY
    # ========================================================

    print(
        "\n\n"
        "============================================================"
    )

    print(
        "ALL-CONTEXT CTCF RESULTS"
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
                "observed_to_random_median_ratio",
                "permutation_z",
                "empirical_p",
                "empirical_two_sided_p"
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