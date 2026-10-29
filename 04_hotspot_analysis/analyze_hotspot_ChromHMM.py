#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import gzip
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


CHROMHMM_DIR = (
    "/projects/clininfo/vcsizmok/long_POG/aDMR/aDMR_analysis/"
    "promoter_clustering/ChrmHMM_downloads"
)

OUTDIR = os.path.join(SCRIPT_DIR, "results", 'ChromHMM_analysis')

N_PERMUTATIONS = 10000
RANDOM_SEED = 12345


# ============================================================
# EPIGENOME PANEL
# ============================================================

EPIGENOMES = {
    "E119_HMEC": "E119_18_core_K27ac_hg38lift_mnemonics.bed.gz",
    "E075_Colon": "E075_18_core_K27ac_hg38lift_mnemonics.bed.gz",
    "E079_Esophagus": "E079_18_core_K27ac_hg38lift_mnemonics.bed.gz",
    "E066_Liver": "E066_18_core_K27ac_hg38lift_mnemonics.bed.gz",
    "E098_Pancreas": "E098_18_core_K27ac_hg38lift_mnemonics.bed.gz",
    "E097_Ovary": "E097_18_core_K27ac_hg38lift_mnemonics.bed.gz",
    "E096_Lung": "E096_18_core_K27ac_hg38lift_mnemonics.bed.gz",
    "E125_Astrocytes": "E125_18_core_K27ac_hg38lift_mnemonics.bed.gz",
    "E127_Keratinocytes": "E127_18_core_K27ac_hg38lift_mnemonics.bed.gz",
    "E129_Osteoblasts": "E129_18_core_K27ac_hg38lift_mnemonics.bed.gz"
}


# ============================================================
# 18-STATE ORDER
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

    # Reproduce original hotspot definition
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
    """
    Inclusive [start, end], matching original hotspot analysis.
    """

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
            promoter_lookup
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

    n_match = int(
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

    if n_match != len(hotspots):

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
            ].to_string(
                index=False
            )
        )

        raise ValueError(
            "Promoter-count validation failed."
        )

    return hotspots


# ============================================================
# LOAD CHROMHMM FILE
# ============================================================

def load_chromhmm(filename):

    df = pd.read_csv(
        filename,
        sep="\t",
        compression="gzip",
        header=None,
        names=[
            "chrom",
            "start",
            "end",
            "state"
        ]
    )

    df["start"] = df["start"].astype(int)
    df["end"] = df["end"].astype(int)

    return df


# ============================================================
# BUILD CHROMHMM LOOKUP
# ============================================================

def build_chromhmm_lookup(df):

    lookup = {}

    for chrom, group in df.groupby(
        "chrom"
    ):

        group = group.sort_values(
            "start"
        )

        lookup[
            chrom
        ] = {
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

            "states":
                group[
                    "state"
                ]
                .astype(str)
                .to_numpy()
        }

    return lookup


# ============================================================
# ASSIGN CHROMHMM STATE TO ONE TSS
# ============================================================

def get_state_at_tss(
    chrom,
    tss,
    chromhmm_lookup
):
    """
    Return one state at TSS.
    If TSS falls in a lift-over gap, return 'unmapped'.
    """

    if chrom not in chromhmm_lookup:

        return "unmapped"

    starts = chromhmm_lookup[
        chrom
    ]["starts"]

    ends = chromhmm_lookup[
        chrom
    ]["ends"]

    states = chromhmm_lookup[
        chrom
    ]["states"]

    idx = np.searchsorted(
        starts,
        tss,
        side="right"
    ) - 1

    if idx < 0:

        return "unmapped"

    if (
        tss >= starts[idx]
        and
        tss < ends[idx]
    ):

        return states[idx]

    return "unmapped"


# ============================================================
# ANNOTATE ALL PROMOTERS
# ============================================================

def annotate_promoters_with_chromhmm(
    promoters,
    chromhmm_lookup,
    epigenome_name
):

    states = []

    print(
        "\nAnnotating promoters for {}".format(
            epigenome_name
        )
    )

    for i, row in promoters.iterrows():

        state = get_state_at_tss(
            row["chrom"],
            int(
                row["tss_start"]
            ),
            chromhmm_lookup
        )

        states.append(
            state
        )

        if (
            (i + 1) % 5000
            == 0
        ):

            print(
                "  {:,} / {:,}".format(
                    i + 1,
                    len(promoters)
                )
            )

    return np.array(
        states,
        dtype=object
    )


# ============================================================
# GET UNIQUE PROMOTERS FROM REGION SET
# ============================================================

def get_region_promoter_indices(
    regions,
    promoter_lookup
):

    all_indices = []

    for _, region in regions.iterrows():

        indices = (
            get_promoter_indices_in_region(
                region[
                    "chrom"
                ],
                int(
                    region[
                        "start"
                    ]
                ),
                int(
                    region[
                        "end"
                    ]
                ),
                promoter_lookup
            )
        )

        all_indices.extend(
            list(
                indices
            )
        )

    return np.array(
        sorted(
            set(
                all_indices
            )
        ),
        dtype=int
    )


# ============================================================
# STATE FRACTIONS
# ============================================================

def calculate_state_fractions(
    promoter_indices,
    state_array,
    state_order
):
    """
    Exclude 'unmapped' from denominator.
    """

    states = state_array[
        promoter_indices
    ]

    mapped_mask = (
        states
        !=
        "unmapped"
    )

    mapped_states = states[
        mapped_mask
    ]

    n_mapped = len(
        mapped_states
    )

    n_total = len(
        states
    )

    result = {}

    if n_mapped == 0:

        for state in state_order:

            result[
                state
            ] = np.nan

    else:

        for state in state_order:

            result[
                state
            ] = (
                np.sum(
                    mapped_states
                    ==
                    state
                )
                /
                float(
                    n_mapped
                )
            )

    coverage = (
        n_mapped
        /
        float(
            n_total
        )
        if n_total > 0
        else np.nan
    )

    return (
        result,
        n_total,
        n_mapped,
        coverage
    )


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
    # >=1 PROMOTERS
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
        "\nPrecomputing promoter-matched "
        "random intervals..."
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

        ranges = (
            get_valid_start_ranges(
                chrom,
                length,
                n_promoters,
                chromsizes,
                promoter_lookup
            )
        )

        if len(
            ranges
        ) == 0:

            raise RuntimeError(
                "No matched regions for "
                "{} length={} promoters={}".format(
                    chrom,
                    length,
                    n_promoters
                )
            )

        sizes = np.array(
            [
                upper
                -
                lower
                +
                1
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

            "sizes":
                sizes,

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
                row[
                    "chrom"
                ],
                int(
                    row[
                        "length"
                    ]
                ),
                int(
                    row[
                        "n_promoters"
                    ]
                ),
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
        ~np.isnan(
            random_values
        )
    ]

    return (
        np.sum(
            values
            >=
            observed
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
        ~np.isnan(
            random_values
        )
    ]

    return (
        np.sum(
            values
            <=
            observed
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

    lower = empirical_p_less(
        observed,
        random_values
    )

    upper = empirical_p_greater(
        observed,
        random_values
    )

    return min(
        1.0,
        2.0
        *
        min(
            lower,
            upper
        )
    )


# ============================================================
# BENJAMINI-HOCHBERG FDR
# ============================================================

def benjamini_hochberg(
    pvalues
):

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

    promoter_lookup = (
        build_promoter_lookup(
            promoters
        )
    )

    print(
        "Hotspots: {}".format(
            len(
                hotspots
            )
        )
    )

    print(
        "Evaluable promoters: {:,}".format(
            len(
                promoters
            )
        )
    )

    hotspots = (
        validate_hotspot_promoter_counts(
            hotspots,
            promoter_lookup
        )
    )

    # ========================================================
    # OBSERVED PROMOTER SET
    # ========================================================

    observed_indices = (
        get_region_promoter_indices(
            hotspots,
            promoter_lookup
        )
    )

    print(
        "\nUnique promoters in observed hotspots: {}".format(
            len(
                observed_indices
            )
        )
    )

    # ========================================================
    # RANDOMIZATION CACHE
    # ========================================================

    cache = (
        precompute_valid_start_ranges(
            hotspots,
            chromsizes,
            promoter_lookup
        )
    )

    # ========================================================
    # GENERATE RANDOM PROMOTER INDEX SETS ONCE
    # ========================================================

    rng = np.random.RandomState(
        RANDOM_SEED
    )

    random_promoter_index_sets = []

    print(
        "\nGenerating {:,} matched random promoter sets...".format(
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

        random_indices = (
            get_region_promoter_indices(
                random_regions,
                promoter_lookup
            )
        )

        random_promoter_index_sets.append(
            random_indices
        )

        if (
            (i + 1) % 1000
            == 0
        ):

            print(
                "Generated {:,} / {:,}".format(
                    i + 1,
                    N_PERMUTATIONS
                )
            )

    # ========================================================
    # ANALYZE EACH EPIGENOME
    # ========================================================

    all_summary_rows = []
    coverage_rows = []

    for (
        epigenome_name,
        filename
    ) in EPIGENOMES.items():

        print(
            "\n\n=============================================="
        )

        print(
            "EPIGENOME: {}".format(
                epigenome_name
            )
        )

        print(
            "=============================================="
        )

        full_path = os.path.join(
            CHROMHMM_DIR,
            filename
        )

        if not os.path.exists(
            full_path
        ):

            raise IOError(
                "Missing ChromHMM file:\n"
                + full_path
            )

        chromhmm = load_chromhmm(
            full_path
        )

        chromhmm_lookup = (
            build_chromhmm_lookup(
                chromhmm
            )
        )

        state_array = (
            annotate_promoters_with_chromhmm(
                promoters,
                chromhmm_lookup,
                epigenome_name
            )
        )

        # ----------------------------------------------------
        # SAVE PROMOTER ANNOTATIONS
        # ----------------------------------------------------

        promoter_output = promoters[
            [
                "promoter_index",
                "symbol",
                "chrom",
                "tss_start",
                "midpoint"
            ]
        ].copy()

        promoter_output[
            "ChromHMM_state"
        ] = state_array

        epigenome_outdir = os.path.join(
            OUTDIR,
            epigenome_name
        )

        os.makedirs(
            epigenome_outdir,
            exist_ok=True
        )

        promoter_output.to_csv(
            os.path.join(
                epigenome_outdir,
                "all_promoters_ChromHMM_annotation.tsv"
            ),
            sep="\t",
            index=False
        )

        # ----------------------------------------------------
        # GLOBAL COVERAGE
        # ----------------------------------------------------

        global_mapped = np.sum(
            state_array
            !=
            "unmapped"
        )

        global_total = len(
            state_array
        )

        global_coverage = (
            global_mapped
            /
            float(
                global_total
            )
        )

        # ----------------------------------------------------
        # OBSERVED STATE FRACTIONS
        # ----------------------------------------------------

        (
            observed_fractions,
            observed_total,
            observed_mapped,
            observed_coverage
        ) = calculate_state_fractions(
            observed_indices,
            state_array,
            STATE_ORDER
        )

        coverage_rows.append({

            "epigenome":
                epigenome_name,

            "all_promoters_total":
                global_total,

            "all_promoters_mapped":
                global_mapped,

            "all_promoters_coverage":
                global_coverage,

            "hotspot_promoters_total":
                observed_total,

            "hotspot_promoters_mapped":
                observed_mapped,

            "hotspot_promoters_coverage":
                observed_coverage
        })

        print(
            "Global annotation coverage: {:.3f}".format(
                global_coverage
            )
        )

        print(
            "Hotspot promoter coverage: {:.3f}".format(
                observed_coverage
            )
        )

        # ----------------------------------------------------
        # PERMUTATION DISTRIBUTIONS
        # ----------------------------------------------------

        distributions = {
            state: []
            for state in STATE_ORDER
        }

        random_coverages = []

        for i, random_indices in enumerate(
            random_promoter_index_sets
        ):

            (
                random_fractions,
                random_total,
                random_mapped,
                random_coverage
            ) = calculate_state_fractions(
                random_indices,
                state_array,
                STATE_ORDER
            )

            for state in STATE_ORDER:

                distributions[
                    state
                ].append(
                    random_fractions[
                        state
                    ]
                )

            random_coverages.append(
                random_coverage
            )

        for state in STATE_ORDER:

            distributions[
                state
            ] = np.array(
                distributions[
                    state
                ],
                dtype=float
            )

        # ----------------------------------------------------
        # SAVE PERMUTATION DISTRIBUTIONS
        # ----------------------------------------------------

        permutation_df = pd.DataFrame({

            "permutation":
                np.arange(
                    1,
                    N_PERMUTATIONS
                    +
                    1
                ),

            "mapped_coverage":
                np.array(
                    random_coverages,
                    dtype=float
                )
        })

        for state in STATE_ORDER:

            permutation_df[
                state
            ] = distributions[
                state
            ]

        permutation_df.to_csv(
            os.path.join(
                epigenome_outdir,
                "ChromHMM_permutation_distributions.tsv"
            ),
            sep="\t",
            index=False
        )

        # ----------------------------------------------------
        # BUILD SUMMARY
        # ----------------------------------------------------

        epigenome_summary_rows = []

        for state in STATE_ORDER:

            observed = (
                observed_fractions[
                    state
                ]
            )

            random_values = (
                distributions[
                    state
                ]
            )

            random_mean = np.nanmean(
                random_values
            )

            random_median = np.nanmedian(
                random_values
            )

            random_sd = np.nanstd(
                random_values,
                ddof=1
            )

            if random_median > 0:

                ratio = (
                    observed
                    /
                    random_median
                )

            else:

                ratio = np.nan

            if random_mean > 0 and observed > 0:

                log2_ratio = np.log2(
                    observed
                    /
                    random_mean
                )

            elif observed == 0 and random_mean > 0:

                log2_ratio = -np.inf

            else:

                log2_ratio = np.nan

            if (
                random_sd > 0
                and
                not np.isnan(
                    random_sd
                )
            ):

                z_perm = (
                    observed
                    -
                    random_mean
                ) / random_sd

            else:

                z_perm = np.nan

            p_enrichment = (
                empirical_p_greater(
                    observed,
                    random_values
                )
            )

            p_depletion = (
                empirical_p_less(
                    observed,
                    random_values
                )
            )

            p_two_sided = (
                empirical_p_two_sided(
                    observed,
                    random_values
                )
            )

            row = {

                "epigenome":
                    epigenome_name,

                "state":
                    state,

                "observed_fraction":
                    observed,

                "random_mean":
                    random_mean,

                "random_median":
                    random_median,

                "random_2.5_percentile":
                    np.nanpercentile(
                        random_values,
                        2.5
                    ),

                "random_97.5_percentile":
                    np.nanpercentile(
                        random_values,
                        97.5
                    ),

                "observed_to_random_median_ratio":
                    ratio,

                "log2_observed_to_random_mean":
                    log2_ratio,

                "permutation_z":
                    z_perm,

                "empirical_p_enrichment":
                    p_enrichment,

                "empirical_p_depletion":
                    p_depletion,

                "empirical_p_two_sided":
                    p_two_sided
            }

            epigenome_summary_rows.append(
                row
            )

            all_summary_rows.append(
                row.copy()
            )

        epigenome_summary = pd.DataFrame(
            epigenome_summary_rows
        )

        epigenome_summary.to_csv(
            os.path.join(
                epigenome_outdir,
                "ChromHMM_permutation_summary.tsv"
            ),
            sep="\t",
            index=False
        )

        print(
            epigenome_summary[
                [
                    "state",
                    "observed_fraction",
                    "random_median",
                    "permutation_z",
                    "empirical_p_two_sided"
                ]
            ].to_string(
                index=False
            )
        )

    # ========================================================
    # COMBINED SUMMARY
    # ========================================================

    combined_summary = pd.DataFrame(
        all_summary_rows
    )

    # BH FDR across all state x epigenome tests
    combined_summary[
        "FDR_BH_all_tests"
    ] = benjamini_hochberg(
        combined_summary[
            "empirical_p_two_sided"
        ].values
    )

    # Also calculate FDR separately within each epigenome
    combined_summary[
        "FDR_BH_within_epigenome"
    ] = np.nan

    for epigenome_name in (
        combined_summary[
            "epigenome"
        ].unique()
    ):

        mask = (
            combined_summary[
                "epigenome"
            ]
            ==
            epigenome_name
        )

        adjusted = (
            benjamini_hochberg(
                combined_summary.loc[
                    mask,
                    "empirical_p_two_sided"
                ].values
            )
        )

        combined_summary.loc[
            mask,
            "FDR_BH_within_epigenome"
        ] = adjusted

    combined_summary.to_csv(
        os.path.join(
            OUTDIR,
            "ALL_EPIGENOMES_ChromHMM_permutation_summary.tsv"
        ),
        sep="\t",
        index=False
    )

    # ========================================================
    # COVERAGE SUMMARY
    # ========================================================

    coverage_df = pd.DataFrame(
        coverage_rows
    )

    coverage_df.to_csv(
        os.path.join(
            OUTDIR,
            "ALL_EPIGENOMES_ChromHMM_annotation_coverage.tsv"
        ),
        sep="\t",
        index=False
    )

    # ========================================================
    # WIDE MATRICES FOR HEATMAPS
    # ========================================================

    z_matrix = combined_summary.pivot(
        index="state",
        columns="epigenome",
        values="permutation_z"
    )

    z_matrix = z_matrix.reindex(
        STATE_ORDER
    )

    z_matrix.to_csv(
        os.path.join(
            OUTDIR,
            "ChromHMM_permutation_z_matrix.tsv"
        ),
        sep="\t"
    )

    log2_matrix = combined_summary.pivot(
        index="state",
        columns="epigenome",
        values="log2_observed_to_random_mean"
    )

    log2_matrix = log2_matrix.reindex(
        STATE_ORDER
    )

    log2_matrix.to_csv(
        os.path.join(
            OUTDIR,
            "ChromHMM_log2_enrichment_matrix.tsv"
        ),
        sep="\t"
    )

    fdr_matrix = combined_summary.pivot(
        index="state",
        columns="epigenome",
        values="FDR_BH_all_tests"
    )

    fdr_matrix = fdr_matrix.reindex(
        STATE_ORDER
    )

    fdr_matrix.to_csv(
        os.path.join(
            OUTDIR,
            "ChromHMM_FDR_matrix.tsv"
        ),
        sep="\t"
    )

    print(
        "\n\n=============================================="
    )

    print(
        "ALL CHROMHMM ANALYSES COMPLETE"
    )

    print(
        "=============================================="
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