"""Call transposon insertion sites from a SAM file and write BED / bedGraph tracks.

Ported verbatim from the CAPTn-seq Snakemake pipeline (pipeline/scripts/sam2sites.py);
only the entrypoint changed (argparse instead of the injected `snakemake` object).
"""

import argparse
import csv

import numpy as np
import pandas as pd


def parse_sam_file(samfile: str, aligned_length_threshold: int = 21, pos_offset: tuple = (3, -6)) -> pd.DataFrame:
    """
    Parses a SAM-like file and retains only first-in-pair reads
    with aligned read length above threshold.

    Determines strand from the SAM flag using bitwise logic.

    pos_offset: (forward_offset, reverse_offset) position-correction offset,
    empirically determined per transposon/library-prep system - passed in
    from config (sam2sites_params.pos_offset_fwd/pos_offset_rev) rather than
    hardcoded here, since a different experiment's transposon may need
    different values.
    """
    df = pd.read_csv(
        samfile,
        sep="\t",
        header=None,
        usecols=range(6),
        names=["name", "flag", "chr_name", "pos", "score", "cigar"],
        dtype={
            "name": str,
            "flag": int,
            "chr_name": str,
            "pos": int,
            "score": int,
            "cigar": str
        },
        engine="python",
        quoting=csv.QUOTE_NONE,
        on_bad_lines="warn"
    )

    # Determine strand and filter to first-in-pair reads only
    df["strand"] = df["flag"].apply(lambda x: "-" if x & 0x10 else "+") #	0x10 == True is mapped to the reverse strand
    df["first_in_pair"] = df["flag"].apply(lambda x: bool(x & 0x40)) #	0x40 == True is the first read in a pair
    df = df[df["first_in_pair"]].copy()

    # use regex to extract each int/marker and build a secondary df holding these info.
    #
    # typic input:
    # ------------
    #
    # 20      42M
    # 25    72M3S
    # Name: cigar
    #
    # typic output:
    # -------------
    #
    #           size marker
    #    match
    # 20 0        42      M
    # 25 0        72      M
    #    1         3      S

    # Extract cigar operations
    cigar_components = df["cigar"].str.extractall(r"(?P<size>\d+)(?P<op>[A-Z])")
    cigar_components["size"] = cigar_components["size"].astype(int)

    # Group cigar operations by index
    cigar_group = cigar_components.groupby(level=0)

    # First and last CIGAR op per read (used for validity check)
    cigar_first_op = cigar_group["op"].first()
    cigar_last_op = cigar_group["op"].last()

    # Aligned length = sum of sizes for relevant ops (M, D, N, =, X)
    aligned_ops = cigar_components[cigar_components["op"].isin(["M", "D", "N", "=", "X"])]
    aligned_length = aligned_ops.groupby(level=0)["size"].sum()

    # Inject back into main df using .map()
    df["cigar_first_op"] = df.index.map(cigar_first_op)
    df["cigar_last_op"] = df.index.map(cigar_last_op)
    df["aligned_read_length"] = df.index.map(aligned_length).fillna(0).astype(int)

    # Define valid reads based on strand and expected CIGAR structure
    valid_positive = (df["strand"] == "+") & (df["cigar_first_op"] == "M")
    valid_negative = (df["strand"] == "-") & (df["cigar_last_op"] == "M")
    df["valid"] = (valid_positive | valid_negative) & (df["aligned_read_length"] > aligned_length_threshold)

    # Position correction offset (empirically determined; not derived from a
    # documented formula), passed in from config. (+offset[0]) shifts a
    # forward-strand read's aligned start to the actual insertion site;
    # (offset[1]) does the same for a reverse-strand read's aligned end.

    # Compute adjusted start position
    valid_df = df[df["valid"]].copy()
    valid_df["start"] = np.where(
        valid_df["strand"] == "+",
        valid_df["pos"] + pos_offset[0],
        valid_df["pos"] + valid_df["aligned_read_length"] + pos_offset[1]
    )
    return valid_df


def filter_reads_by_site_support(sam_df: pd.DataFrame, score_threshold: int) -> pd.DataFrame:
    """
    Retain only reads whose (chr_name, start) occurs with at least `score_threshold` other reads.
    Keeps all columns including QNAMEs (name) and strand info.
    """
    group_keys = ["chr_name", "start"]
    # Count number of reads per position per strand
    site_counts = (
        sam_df.groupby(group_keys + ["strand"])["name"]
        .count()
        .unstack(fill_value=0)
        # CHANGE: unstack() only produces a column for strands actually
        #         present in sam_df. If every surviving read happens to be on
        #         one strand (plausible for a low-coverage/failing sample),
        #         "+" or "-" would be missing entirely and the score_pos/
        #         score_neg lookups below would raise a KeyError. Reindexing
        #         guarantees both columns always exist, filled with 0.
        .reindex(columns=["+", "-"], fill_value=0)
        .rename(columns={"+": "score_pos", "-": "score_neg"})
        .reset_index()
    )
    # Keep only positions with sufficient read support on either strand
    supported_sites = site_counts[
        (site_counts["score_pos"] > score_threshold) |
        (site_counts["score_neg"] > score_threshold)
    ][["chr_name", "start"]]
    # Inner join back with reads
    filtered_reads = sam_df.merge(supported_sites, on=["chr_name", "start"], how="inner")
    return filtered_reads


def sam_to_site(filtered_dedup_reads: pd.DataFrame) -> pd.DataFrame:
    """
    Convert a filtered and deduplicated read set into strand-specific insertion site counts.

    Parameters
    ----------
    filtered_dedup_reads : pd.DataFrame
        DataFrame with at least the following columns:
        - 'chr_name': reference sequence
        - 'start': insertion position (already corrected)
        - 'strand': '+' or '-'
        - 'name': QNAME (used for counting insertions)

    Returns
    -------
    pd.DataFrame
        A DataFrame with one row per (chr_name, start), containing:
        - 'chr_name'
        - 'start'
        - 'end': start + 1 (BED-style interval)
        - 'name': fixed to 'i' (placeholder for BED)
        - 'score_pos': insertion count on + strand
        - 'score_neg': insertion count on - strand
    """
    group_keys = ["chr_name", "start"]
    score_pos = (
        filtered_dedup_reads[filtered_dedup_reads["strand"] == "+"]
        .groupby(group_keys)["name"]
        .count()
        .rename("score_pos")
    )
    score_neg = (
        filtered_dedup_reads[filtered_dedup_reads["strand"] == "-"]
        .groupby(group_keys)["name"]
        .count()
        .rename("score_neg")
    )

    site_df = pd.concat([score_pos, score_neg], axis=1).reset_index()
    site_df["end"] = site_df["start"] + 1
    site_df["name"] = "i"
    site_df = site_df.fillna({"score_pos": 0, "score_neg": 0})
    site_df = site_df[["chr_name", "start", "end", "name", "score_pos", "score_neg"]]
    site_df = site_df.sort_values(by=["chr_name", "start"])
    return site_df


def site_to_interval(site_df: pd.DataFrame, *, normalization_value: int) -> pd.DataFrame:
    total_reads = site_df["score_pos"].sum() + site_df["score_neg"].sum()
    pos_interval_df = (
        site_df[["chr_name", "start", "end", "name", "score_pos"]]
        .rename(columns={"score_pos": "score"})
    )
    pos_interval_df["strand"] = "+"
    neg_interval_df = (
        site_df[["chr_name", "start", "end", "name", "score_neg"]]
        .rename(columns={"score_neg": "score"})
    )
    neg_interval_df["strand"] = "-"
    interval_df = pd.concat([pos_interval_df, neg_interval_df], ignore_index=True)
    # Normalize scores to target total, clip at 1 to ensure non-zero presence
    interval_df["normalized_score"] = (
        (interval_df["score"] / total_reads * normalization_value)
        .astype(int)
        .clip(lower=1)
    )
    interval_df["fake_score"] = 999 # Placeholder fixed score
    interval_df = interval_df.sort_values(by=["chr_name", "start"])
    interval_df = interval_df[interval_df["score"] != 0]
    return interval_df


def build_bed_file(interval_df: pd.DataFrame, filename):
    """
    bed_file
    --------

    "chr_name", "start", "end", "name", "fake_score", "strand"

    pVCR94deltaX3deltaacr2FRT	26	27	i	999	+
    pVCR94deltaX3deltaacr2FRT	47	48	i	999	+
    pVCR94deltaX3deltaacr2FRT	68	69	i	999	+
    pVCR94deltaX3deltaacr2FRT	97	98	i	999	-
    pVCR94deltaX3deltaacr2FRT	245	246	i	999	-
    """
    interval_df[["chr_name", "start", "end", "name", "fake_score", "strand"]].to_csv(
        filename, header=False, index=False, sep="\t"
    )


def build_stranded_unnormalized_bedgraph_file(
    interval_df: pd.DataFrame, filename):
    """
    stranded_unnormalized_bedgraph_file
    --------

    "chr_name", "start", "end", "score"

    pVCR94deltaX3deltaacr2FRT	26	27	2400
    pVCR94deltaX3deltaacr2FRT	47	48	60
    pVCR94deltaX3deltaacr2FRT	68	69	60
    pVCR94deltaX3deltaacr2FRT	97	98	1380
    pVCR94deltaX3deltaacr2FRT	245	246	3780
    """
    interval_df[["chr_name", "start", "end", "score"]].to_csv(
        filename, header=False, index=False, sep="\t"
    )


def build_stranded_bedgraph_file(interval_df: pd.DataFrame, filename):
    """
    stranded_bedgraph_file
    --------

    "chr_name", "start", "end", "normalized_score"

    pVCR94deltaX3deltaacr2FRT	26	27	2189
    pVCR94deltaX3deltaacr2FRT	47	48	54
    pVCR94deltaX3deltaacr2FRT	47	48	54
    pVCR94deltaX3deltaacr2FRT	68	69	54
    pVCR94deltaX3deltaacr2FRT	97	98	1258
    pVCR94deltaX3deltaacr2FRT	245	246	3448
    """
    interval_df[["chr_name", "start", "end", "normalized_score"]].to_csv(
        filename, header=False, index=False, sep="\t"
    )


def build_unstranded_bedgraph_file(interval_df: pd.DataFrame, filename):
    """
    unstranded_bedgraph_file
    --------

    "chr_name", "start", "end", "normalized_score"

    pVCR94deltaX3deltaacr2FRT	26	27	2189
    pVCR94deltaX3deltaacr2FRT	47	48	108<-(score for insertions on the same position, but different strands are merged)
    pVCR94deltaX3deltaacr2FRT	68	69	54
    pVCR94deltaX3deltaacr2FRT	97	98	1258
    pVCR94deltaX3deltaacr2FRT	245	246	3448
    """
    interval_df = interval_df.groupby(["chr_name", "start", "end"], sort=False).agg({'normalized_score': 'sum'}).reset_index()
    interval_df[["chr_name", "start", "end", "normalized_score"]].to_csv(
        filename, header=False, index=False, sep="\t"
    )


def build_unstranded_unnormalized_bedgraph_file(interval_df: pd.DataFrame, filename):
    """
    unstranded_bedgraph_file
    --------

    "chr_name", "start", "end", "score"

    pVCR94deltaX3deltaacr2FRT	26	27	2400
    pVCR94deltaX3deltaacr2FRT	47	48	60
    pVCR94deltaX3deltaacr2FRT	68	69	60<-(score for insertions on the same position, but different strands are merged)
    pVCR94deltaX3deltaacr2FRT	97	98	1380
    pVCR94deltaX3deltaacr2FRT	245	246	3780
    """

    interval_df = interval_df.groupby(["chr_name", "start", "end"], sort=False).agg({'score': 'sum'}).reset_index()
    interval_df[["chr_name", "start", "end", "score"]].to_csv(
        filename, header=False, index=False, sep="\t"
    )


def add_arguments(parser):
    """Register this subcommand's options on *parser*."""
    parser.add_argument("--sam", required=True)
    parser.add_argument("--norm-value", type=int, required=True)
    parser.add_argument("--read-len-threshold", type=int, required=True)
    parser.add_argument("--score-threshold", type=int, required=True)
    parser.add_argument("--pos-offset-fwd", type=int, required=True)
    parser.add_argument("--pos-offset-rev", type=int, required=True)
    parser.add_argument("--bed-output", required=True)
    parser.add_argument("--normalized-output", required=True)
    parser.add_argument("--unstranded-output", required=True)
    parser.add_argument("--unnormalized-output", required=True)
    parser.add_argument("--unstranded-unnormalized-output", required=True)


def run_args(args):
    """Run the subcommand from an ``argparse.Namespace``."""
    pos_offset = (args.pos_offset_fwd, args.pos_offset_rev)
    sam_df = parse_sam_file(args.sam, args.read_len_threshold, pos_offset)
    filtered_reads = filter_reads_by_site_support(sam_df, args.score_threshold)
    site_df = sam_to_site(filtered_reads)
    interval_df = site_to_interval(site_df, normalization_value=args.norm_value)
    build_bed_file(interval_df, args.bed_output)
    build_stranded_bedgraph_file(interval_df, args.normalized_output)
    build_unstranded_bedgraph_file(interval_df, args.unstranded_output)
    build_stranded_unnormalized_bedgraph_file(interval_df, args.unnormalized_output)
    build_unstranded_unnormalized_bedgraph_file(interval_df, args.unstranded_unnormalized_output)


def main(argv=None):
    parser = argparse.ArgumentParser(
        prog="captnseq-toolkit sam2sites",
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    add_arguments(parser)
    run_args(parser.parse_args(argv))


if __name__ == "__main__":
    main()
