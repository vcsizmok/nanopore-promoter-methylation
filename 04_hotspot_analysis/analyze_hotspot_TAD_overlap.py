#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import pandas as pd
import numpy as np


# ============================================================
# USER SETTINGS
# ============================================================

# Repository inputs and shared GSC annotations (override via CLI if needed).
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
HOTSPOT_FILE = os.path.join(SCRIPT_DIR, "data", "methylated_promoter_hotspots.bed")
PROMOTER_FILE = os.path.join(SCRIPT_DIR, "data", "input_promoters_coordinates_methylation.tsv")
CHROMSIZES_FILE = "/gsc/www/bcgsc.ca/downloads/promoter_methylation_nanopore/annotations/genome/hg38.chrom.sizes"


TAD_DIR = (
    "/projects/clininfo/vcsizmok/long_POG/aDMR/aDMR_analysis/"
    "promoter_clustering/TAD_downloads"
)

OUTDIR = os.path.join(SCRIPT_DIR, "results", 'TAD_overlap_analysis')

TAD_MAPS = {

    "GM12878":
        os.path.join(
            TAD_DIR,
            "GM12878_TADs_GRCh38.bed"
        ),

    "mammary":
        os.path.join(
            TAD_DIR,
            "ENCFF840CVX_TADs_GRCh38.bed"
        ),

    "multi_cancer":
        os.path.join(
            TAD_DIR,
            "ENCFF549OBE_TADs_GRCh38.bed"
        ),

    "gonadal":
        os.path.join(
            TAD_DIR,
            "ENCFF514ILX_TADs_GRCh38.bed"
        ),

    "thyroid_prostate":
        os.path.join(
            TAD_DIR,
            "ENCFF572UTG_TADs_GRCh38.bed"
        ),

    "pancreas":
        os.path.join(
            TAD_DIR,
            "ENCFF499QXR_TADs_GRCh38.bed"
        ),

    "brain":
        os.path.join(
            TAD_DIR,
            "ENCFF204WPQ_TADs_GRCh38.bed"
        )
}


N_PERMUTATIONS = 10000
RANDOM_SEED = 12345


# ============================================================
# LOAD INPUT FILES
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
        col for col in required
        if col not in df.columns
    ]

    if missing:
        raise ValueError(
            "Missing hotspot columns: "
            + ", ".join(missing)
        )

    df = df.rename(
        columns={"chr": "chrom"}
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
        df["end"] - df["start"]
    )

    return df


def load_promoters(filename):

    df = pd.read_csv(
        filename,
        sep="\t"
    )

    df.columns = df.columns.str.strip()

    required = [
        "chrom",
        "promoter_start",
        "promoter_end"
    ]

    missing = [
        col for col in required
        if col not in df.columns
    ]

    if missing:
        raise ValueError(
            "Missing promoter columns: "
            + ", ".join(missing)
        )

    df["promoter_start"] = pd.to_numeric(
        df["promoter_start"],
        errors="raise"
    )

    df["promoter_end"] = pd.to_numeric(
        df["promoter_end"],
        errors="raise"
    )

    # Reproduce original hotspot analysis
    df["midpoint"] = (
        (
            df["promoter_start"]
            +
            df["promoter_end"]
        ) // 2
    ).astype(int)

    return df


def load_tads(filename):

    df = pd.read_csv(
        filename,
        sep="\t",
        header=None,
        comment="#"
    )

    if df.shape[1] < 3:

        raise ValueError(
            "TAD file must contain at least 3 columns: "
            + filename
        )

    df = df.iloc[:, :3].copy()

    df.columns = [
        "chrom",
        "start",
        "end"
    ]

    df["start"] = df["start"].astype(int)
    df["end"] = df["end"].astype(int)

    return df


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
# PROMOTER FUNCTIONS
# ============================================================

def build_promoter_lookup(promoters):

    promoter_dict = {}

    for chrom, group in promoters.groupby(
        "chrom"
    ):

        positions = np.sort(
            group[
                "midpoint"
            ]
            .astype(int)
            .to_numpy()
        )

        promoter_dict[
            chrom
        ] = positions

    return promoter_dict


def count_promoters_in_region(
    chrom,
    start,
    end,
    promoter_dict
):

    """
    Count promoter midpoints in inclusive interval [start, end].

    This matches the region-level definition used in the
    original hotspot analysis.
    """

    if chrom not in promoter_dict:
        return 0

    positions = promoter_dict[
        chrom
    ]

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

    return int(
        right - left
    )


def validate_hotspot_promoter_counts(
    hotspots,
    promoter_dict
):

    recalculated = []

    for _, row in hotspots.iterrows():

        count = count_promoters_in_region(
            row["chrom"],
            int(
                row["start"]
            ),
            int(
                row["end"]
            ),
            promoter_dict
        )

        recalculated.append(
            count
        )

    hotspots = hotspots.copy()

    hotspots[
        "promoter_count_recalculated"
    ] = recalculated

    hotspots[
        "promoter_count_matches"
    ] = (
        hotspots[
            "n_promoters"
        ]
        ==
        hotspots[
            "promoter_count_recalculated"
        ]
    )

    n_match = (
        hotspots[
            "promoter_count_matches"
        ].sum()
    )

    print(
        "\nPromoter-count validation:"
    )

    print(
        "{} / {} hotspots reproduce "
        "the n_promoters value".format(
            n_match,
            len(hotspots)
        )
    )

    if n_match != len(
        hotspots
    ):

        print(
            hotspots.loc[
                ~hotspots[
                    "promoter_count_matches"
                ],
                [
                    "chrom",
                    "start",
                    "end",
                    "n_promoters",
                    "promoter_count_recalculated"
                ]
            ]
        )

        raise ValueError(
            "Promoter-count validation failed."
        )

    return hotspots


# ============================================================
# PROMOTER-MATCHED RANDOMIZATION
# ============================================================

def get_valid_start_ranges(
    chrom,
    length,
    target_n_promoters,
    chromsizes,
    promoter_dict
):

    """
    Find all genomic start-coordinate ranges where an interval
    [start, start + length] contains exactly target_n_promoters
    promoter midpoints.
    """

    if chrom not in chromsizes:

        raise ValueError(
            "Missing chromosome size: "
            + chrom
        )

    if chrom not in promoter_dict:

        raise ValueError(
            "Missing promoters for: "
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

    positions = promoter_dict[
        chrom
    ]

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
    # Exactly zero promoters
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
    # >=1 promoters
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

        # Exclude previous promoter
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

        # Exclude next promoter
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
    promoter_dict
):

    cache = {}

    unique_combinations = (
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

    for _, row in unique_combinations.iterrows():

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

        valid_ranges = (
            get_valid_start_ranges(
                chrom,
                length,
                n_promoters,
                chromsizes,
                promoter_dict
            )
        )

        if len(
            valid_ranges
        ) == 0:

            raise RuntimeError(
                "No promoter-matched region exists: "
                "{} length={} promoters={}".format(
                    chrom,
                    length,
                    n_promoters
                )
            )

        range_sizes = np.array(
            [
                upper
                -
                lower
                +
                1
                for lower, upper
                in valid_ranges
            ],
            dtype=np.int64
        )

        total_valid_starts = int(
            range_sizes.sum()
        )

        cumulative_sizes = (
            np.cumsum(
                range_sizes
            )
        )

        cache[
            key
        ] = {
            "ranges":
                valid_ranges,

            "range_sizes":
                range_sizes,

            "cumulative_sizes":
                cumulative_sizes,

            "total_valid_starts":
                total_valid_starts
        }

        print(
            "{} length={} promoters={} : "
            "{} valid ranges, {:,} valid starts".format(
                chrom,
                length,
                n_promoters,
                len(
                    valid_ranges
                ),
                total_valid_starts
            )
        )

    return cache


def generate_one_promoter_matched_region(
    chrom,
    length,
    target_n_promoters,
    cache,
    rng
):

    key = (
        chrom,
        length,
        target_n_promoters
    )

    item = cache[
        key
    ]

    random_index = int(
        rng.randint(
            0,
            item[
                "total_valid_starts"
            ]
        )
    )

    range_index = int(
        np.searchsorted(
            item[
                "cumulative_sizes"
            ],
            random_index,
            side="right"
        )
    )

    previous_cumulative = 0

    if range_index > 0:

        previous_cumulative = int(
            item[
                "cumulative_sizes"
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

    offset = (
        random_index
        -
        previous_cumulative
    )

    selected_start = int(
        lower
        +
        offset
    )

    selected_end = (
        selected_start
        +
        length
    )

    return {
        "chrom":
            chrom,

        "start":
            selected_start,

        "end":
            selected_end,

        "length":
            length,

        "n_promoters":
            target_n_promoters
    }


def generate_promoter_matched_regions(
    hotspots,
    cache,
    rng
):

    regions = []

    for _, row in hotspots.iterrows():

        region = (
            generate_one_promoter_matched_region(
                chrom=
                    row[
                        "chrom"
                    ],

                length=
                    int(
                        row[
                            "length"
                        ]
                    ),

                target_n_promoters=
                    int(
                        row[
                            "n_promoters"
                        ]
                    ),

                cache=
                    cache,

                rng=
                    rng
            )
        )

        regions.append(
            region
        )

    return pd.DataFrame(
        regions
    )


# ============================================================
# PRE-GENERATE RANDOM REGION SETS
# ============================================================

def generate_all_random_region_sets(
    hotspots,
    cache,
    n_permutations,
    seed
):

    """
    Generate promoter-matched random sets ONCE.

    The exact same random sets are then evaluated against
    every TAD map, making comparison among TAD maps cleaner.
    """

    rng = np.random.RandomState(
        seed
    )

    random_sets = []

    print(
        "\nGenerating promoter-matched "
        "random region sets..."
    )

    for i in range(
        n_permutations
    ):

        random_regions = (
            generate_promoter_matched_regions(
                hotspots,
                cache,
                rng
            )
        )

        random_sets.append(
            random_regions
        )

        if (
            (i + 1)
            % 1000
            == 0
        ):

            print(
                "Generated {:,} / {:,}".format(
                    i + 1,
                    n_permutations
                )
            )

    return random_sets


# ============================================================
# TAD FUNCTIONS
# ============================================================

def build_tad_lookup(
    tads
):

    tad_dict = {}

    for chrom, group in tads.groupby(
        "chrom"
    ):

        group = group.sort_values(
            [
                "start",
                "end"
            ]
        )

        tad_dict[
            chrom
        ] = group[
            [
                "start",
                "end"
            ]
        ].to_numpy()

    return tad_dict


def build_boundary_lookup(
    tad_dict
):

    boundary_dict = {}

    for chrom, tads in tad_dict.items():

        boundaries = np.concatenate(
            [
                tads[:, 0],
                tads[:, 1]
            ]
        )

        boundary_dict[
            chrom
        ] = np.unique(
            np.sort(
                boundaries
            )
        )

    return boundary_dict


def classify_hotspot(
    chrom,
    start,
    end,
    tad_dict
):

    if chrom not in tad_dict:

        return (
            "outside",
            [],
            []
        )

    overlapping_tads = []
    containing_tads = []

    for tad_start, tad_end in tad_dict[
        chrom
    ]:

        if tad_end <= start:
            continue

        if tad_start >= end:
            break

        if (
            tad_start < end
            and
            tad_end > start
        ):

            overlapping_tads.append(
                (
                    int(
                        tad_start
                    ),
                    int(
                        tad_end
                    )
                )
            )

            if (
                tad_start
                <=
                start
                and
                tad_end
                >=
                end
            ):

                containing_tads.append(
                    (
                        int(
                            tad_start
                        ),
                        int(
                            tad_end
                        )
                    )
                )

    if containing_tads:

        return (
            "contained",
            containing_tads,
            overlapping_tads
        )

    if overlapping_tads:

        return (
            "partial",
            [],
            overlapping_tads
        )

    return (
        "outside",
        [],
        []
    )


def nearest_boundary_distance(
    chrom,
    position,
    boundary_dict
):

    if chrom not in boundary_dict:
        return np.nan

    boundaries = boundary_dict[
        chrom
    ]

    idx = np.searchsorted(
        boundaries,
        position
    )

    distances = []

    if idx < len(
        boundaries
    ):

        distances.append(
            abs(
                boundaries[
                    idx
                ]
                -
                position
            )
        )

    if idx > 0:

        distances.append(
            abs(
                boundaries[
                    idx - 1
                ]
                -
                position
            )
        )

    if not distances:
        return np.nan

    return min(
        distances
    )


# ============================================================
# ANALYZE REGIONS
# ============================================================

def analyze_regions(
    regions,
    tad_dict,
    boundary_dict
):

    results = []

    for _, row in regions.iterrows():

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

        (
            status,
            containing_tads,
            overlapping_tads
        ) = classify_hotspot(
            chrom,
            start,
            end,
            tad_dict
        )

        start_dist = (
            nearest_boundary_distance(
                chrom,
                start,
                boundary_dict
            )
        )

        end_dist = (
            nearest_boundary_distance(
                chrom,
                end,
                boundary_dict
            )
        )

        result = {
            "chrom":
                chrom,

            "start":
                start,

            "end":
                end,

            "length":
                end - start,

            "TAD_status":
                status,

            "overlaps_TAD":
                status in [
                    "contained",
                    "partial"
                ],

            "contained_in_TAD":
                status
                ==
                "contained",

            "n_overlapping_TADs":
                len(
                    overlapping_tads
                ),

            "n_containing_TADs":
                len(
                    containing_tads
                ),

            "start_boundary_distance":
                start_dist,

            "end_boundary_distance":
                end_dist,

            "mean_edge_boundary_distance":
                np.nanmean(
                    [
                        start_dist,
                        end_dist
                    ]
                ),

            "minimum_edge_boundary_distance":
                np.nanmin(
                    [
                        start_dist,
                        end_dist
                    ]
                )
        }

        for col in regions.columns:

            if col not in result:

                result[
                    col
                ] = row[
                    col
                ]

        results.append(
            result
        )

    return pd.DataFrame(
        results
    )


# ============================================================
# ANALYZE ONE TAD MAP
# ============================================================

def analyze_one_tad_map(
    map_name,
    tad_file,
    hotspots,
    random_sets
):

    print(
        "\n\n"
        "########################################"
    )

    print(
        "TAD MAP: {}".format(
            map_name
        )
    )

    print(
        "########################################"
    )

    print(
        tad_file
    )

    tads = load_tads(
        tad_file
    )

    print(
        "Number of TAD records: {}".format(
            len(
                tads
            )
        )
    )

    tad_dict = (
        build_tad_lookup(
            tads
        )
    )

    boundary_dict = (
        build_boundary_lookup(
            tad_dict
        )
    )

    # --------------------------------------------------------
    # Observed
    # --------------------------------------------------------

    observed = analyze_regions(
        hotspots,
        tad_dict,
        boundary_dict
    )

    observed_overlap = (
        observed[
            "overlaps_TAD"
        ].mean()
    )

    observed_containment = (
        observed[
            "contained_in_TAD"
        ].mean()
    )

    observed_mean_edge = (
        observed[
            "mean_edge_boundary_distance"
        ].median()
    )

    observed_min_edge = (
        observed[
            "minimum_edge_boundary_distance"
        ].median()
    )

    # --------------------------------------------------------
    # Random
    # --------------------------------------------------------

    random_overlap = []
    random_containment = []
    random_mean_edge = []
    random_min_edge = []

    for i, random_regions in enumerate(
        random_sets
    ):

        random_analysis = (
            analyze_regions(
                random_regions,
                tad_dict,
                boundary_dict
            )
        )

        random_overlap.append(
            random_analysis[
                "overlaps_TAD"
            ].mean()
        )

        random_containment.append(
            random_analysis[
                "contained_in_TAD"
            ].mean()
        )

        random_mean_edge.append(
            random_analysis[
                "mean_edge_boundary_distance"
            ].median()
        )

        random_min_edge.append(
            random_analysis[
                "minimum_edge_boundary_distance"
            ].median()
        )

        if (
            (i + 1)
            % 1000
            == 0
        ):

            print(
                "{}: analyzed {:,} / {:,} permutations".format(
                    map_name,
                    i + 1,
                    len(
                        random_sets
                    )
                )
            )

    random_overlap = np.array(
        random_overlap
    )

    random_containment = np.array(
        random_containment
    )

    random_mean_edge = np.array(
        random_mean_edge
    )

    random_min_edge = np.array(
        random_min_edge
    )

    # --------------------------------------------------------
    # Empirical one-sided P values
    # --------------------------------------------------------

    p_overlap = (
        np.sum(
            random_overlap
            >=
            observed_overlap
        )
        +
        1
    ) / float(
        len(
            random_overlap
        )
        +
        1
    )

    p_containment = (
        np.sum(
            random_containment
            >=
            observed_containment
        )
        +
        1
    ) / float(
        len(
            random_containment
        )
        +
        1
    )

    p_mean_edge = (
        np.sum(
            random_mean_edge
            <=
            observed_mean_edge
        )
        +
        1
    ) / float(
        len(
            random_mean_edge
        )
        +
        1
    )

    p_min_edge = (
        np.sum(
            random_min_edge
            <=
            observed_min_edge
        )
        +
        1
    ) / float(
        len(
            random_min_edge
        )
        +
        1
    )

    # --------------------------------------------------------
    # Per-map output directory
    # --------------------------------------------------------

    map_outdir = os.path.join(
        OUTDIR,
        map_name
    )

    os.makedirs(
        map_outdir,
        exist_ok=True
    )

    # Observed classifications
    observed.to_csv(
        os.path.join(
            map_outdir,
            "hotspot_TAD_classification.tsv"
        ),
        sep="\t",
        index=False
    )

    # Permutation distributions
    permutation_df = pd.DataFrame({

        "permutation":
            np.arange(
                1,
                len(
                    random_sets
                )
                +
                1
            ),

        "fraction_TAD_overlap":
            random_overlap,

        "fraction_TAD_contained":
            random_containment,

        "median_mean_edge_boundary_distance_bp":
            random_mean_edge,

        "median_minimum_edge_boundary_distance_bp":
            random_min_edge
    })

    permutation_df.to_csv(
        os.path.join(
            map_outdir,
            "TAD_permutation_distributions.tsv"
        ),
        sep="\t",
        index=False
    )

    # --------------------------------------------------------
    # Map-level summary
    # --------------------------------------------------------

    summary_rows = []

    tests = [

        (
            "any_TAD_overlap",
            observed_overlap,
            random_overlap,
            p_overlap
        ),

        (
            "complete_TAD_containment",
            observed_containment,
            random_containment,
            p_containment
        ),

        (
            "mean_edge_TAD_boundary_proximity_bp",
            observed_mean_edge,
            random_mean_edge,
            p_mean_edge
        ),

        (
            "minimum_edge_TAD_boundary_proximity_bp",
            observed_min_edge,
            random_min_edge,
            p_min_edge
        )
    ]

    for (
        test_name,
        observed_value,
        random_values,
        p_value
    ) in tests:

        summary_rows.append({

            "TAD_map":
                map_name,

            "TAD_file":
                tad_file,

            "test":
                test_name,

            "observed":
                observed_value,

            "random_mean":
                np.mean(
                    random_values
                ),

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

            "empirical_p":
                p_value
        })

    summary_df = pd.DataFrame(
        summary_rows
    )

    summary_df.to_csv(
        os.path.join(
            map_outdir,
            "TAD_permutation_summary.tsv"
        ),
        sep="\t",
        index=False
    )

    # --------------------------------------------------------
    # Console summary
    # --------------------------------------------------------

    print(
        "\nObserved classification:"
    )

    print(
        observed[
            "TAD_status"
        ].value_counts()
    )

    print(
        "\nANY TAD OVERLAP"
    )

    print(
        "Observed: {:.3f}".format(
            observed_overlap
        )
    )

    print(
        "Random mean: {:.3f}".format(
            np.mean(
                random_overlap
            )
        )
    )

    print(
        "P = {:.4g}".format(
            p_overlap
        )
    )

    print(
        "\nCOMPLETE CONTAINMENT"
    )

    print(
        "Observed: {:.3f}".format(
            observed_containment
        )
    )

    print(
        "Random mean: {:.3f}".format(
            np.mean(
                random_containment
            )
        )
    )

    print(
        "P = {:.4g}".format(
            p_containment
        )
    )

    print(
        "\nMEAN EDGE DISTANCE"
    )

    print(
        "Observed: {:,.0f} bp".format(
            observed_mean_edge
        )
    )

    print(
        "Random median: {:,.0f} bp".format(
            np.median(
                random_mean_edge
            )
        )
    )

    print(
        "P = {:.4g}".format(
            p_mean_edge
        )
    )

    print(
        "\nMINIMUM EDGE DISTANCE"
    )

    print(
        "Observed: {:,.0f} bp".format(
            observed_min_edge
        )
    )

    print(
        "Random median: {:,.0f} bp".format(
            np.median(
                random_min_edge
            )
        )
    )

    print(
        "P = {:.4g}".format(
            p_min_edge
        )
    )

    return summary_df


# ============================================================
# MAIN
# ============================================================

def main():

    os.makedirs(
        OUTDIR,
        exist_ok=True
    )

    print(
        "Loading hotspot and promoter data..."
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

    promoter_dict = (
        build_promoter_lookup(
            promoters
        )
    )

    print(
        "Number of hotspots: {}".format(
            len(
                hotspots
            )
        )
    )

    print(
        "Number of evaluable promoters: {}".format(
            len(
                promoters
            )
        )
    )

    # ========================================================
    # VALIDATE ORIGINAL PROMOTER COUNTS
    # ========================================================

    hotspots = (
        validate_hotspot_promoter_counts(
            hotspots,
            promoter_dict
        )
    )

    hotspots.to_csv(
        os.path.join(
            OUTDIR,
            "hotspot_promoter_count_validation.tsv"
        ),
        sep="\t",
        index=False
    )

    # ========================================================
    # PRECOMPUTE PROMOTER MATCHES
    # ========================================================

    promoter_match_cache = (
        precompute_valid_start_ranges(
            hotspots,
            chromsizes,
            promoter_dict
        )
    )

    # ========================================================
    # GENERATE RANDOM REGION SETS ONCE
    # ========================================================

    random_sets = (
        generate_all_random_region_sets(
            hotspots=
                hotspots,

            cache=
                promoter_match_cache,

            n_permutations=
                N_PERMUTATIONS,

            seed=
                RANDOM_SEED
        )
    )

    # ========================================================
    # RUN EACH TAD MAP
    # ========================================================

    all_summaries = []

    for (
        map_name,
        tad_file
    ) in TAD_MAPS.items():

        if not os.path.exists(
            tad_file
        ):

            print(
                "\nWARNING: TAD file not found, skipping:"
            )

            print(
                tad_file
            )

            continue

        map_summary = (
            analyze_one_tad_map(
                map_name=
                    map_name,

                tad_file=
                    tad_file,

                hotspots=
                    hotspots,

                random_sets=
                    random_sets
            )
        )

        all_summaries.append(
            map_summary
        )

    # ========================================================
    # COMBINED SUMMARY
    # ========================================================

    if len(
        all_summaries
    ) == 0:

        raise RuntimeError(
            "No TAD maps were analyzed."
        )

    combined_summary = pd.concat(
        all_summaries,
        ignore_index=True
    )

    combined_summary.to_csv(
        os.path.join(
            OUTDIR,
            "ALL_TAD_MAPS_permutation_summary.tsv"
        ),
        sep="\t",
        index=False
    )

    # ========================================================
    # WIDE SUMMARY FOR EASY COMPARISON
    # ========================================================

    wide_summary = (
        combined_summary[
            [
                "TAD_map",
                "test",
                "observed",
                "random_median",
                "random_2.5_percentile",
                "random_97.5_percentile",
                "empirical_p"
            ]
        ]
        .copy()
    )

    wide_summary.to_csv(
        os.path.join(
            OUTDIR,
            "ALL_TAD_MAPS_comparison.tsv"
        ),
        sep="\t",
        index=False
    )

    print(
        "\n\n======================================"
    )

    print(
        "ALL TAD MAPS COMPLETE"
    )

    print(
        "======================================"
    )

    print(
        "\nCombined summary:"
    )

    print(
        combined_summary[
            [
                "TAD_map",
                "test",
                "observed",
                "random_median",
                "empirical_p"
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