#! /usr/bin/env python3

"""
Ribosome profiling data analysis functionality for Ribomala.
This module processes BAM files listed in a sample sheet.
Calculates E-, P- and A-site ribosome occupancy
"""

import logging
import os
from pathlib import Path
import polars as pl
import sys
from typing import Tuple, List, Dict

from ribomala import qc

# --------------------------------------------------


def parse_sample_sheet(
    sample_sheet_path: str,
) -> Tuple[List[str], Dict[str, pl.DataFrame]]:
    """
    Parse the sample sheet and extract file names and their corresponding read length, frame, and offset information.

    Parameters
    ----------
        sample_sheet_path: Path to the sample sheet CSV file

    Returns
    -------
        A tuple containing:
        - List of unique file names
        - Dictionary mapping file names to their read length, frame, and offset information as a polars DataFrame
    """

    logging.info(f"Parsing sample sheet: {sample_sheet_path}")

    # Convert to Path object
    sample_sheet = Path(sample_sheet_path)

    if not sample_sheet.exists():
        logging.error(f"Sample sheet not found: {sample_sheet}")
        raise FileNotFoundError(f"Sample sheet not found: {sample_sheet}")

    # Read the sample sheet
    try:
        df = pl.read_csv(sample_sheet, separator="\t")
        logging.debug(f"Successfully read sample sheet with {len(df)} entries")
    except Exception as e:
        logging.error(f"Failed to read sample sheet: {e}")
        raise

    # Validate required columns
    required_columns = ["file_name", "read_length", "frame", "offset"]
    missing_columns = [col for col in required_columns if col not in df.columns]
    if missing_columns:
        logging.error(f"Sample sheet missing required columns: {missing_columns}")
        raise ValueError(f"Sample sheet missing required columns: {missing_columns}")

    # Get unique file names
    unique_files = df["file_name"].unique().to_list()
    logging.info(f"Found {len(unique_files)} unique files in sample sheet")

    # Create a dictionary to store file-specific information
    file_info = {}

    # For each unique file, extract the corresponding rows and keep only read_length, frame, and offset columns
    for file_name in unique_files:
        file_df = (
            df.filter(pl.col("file_name") == file_name)
            .select(["read_length", "frame", "offset"])
            .unique()
        )
        file_info[file_name] = file_df
        logging.debug(
            f"File {file_name}: found {len(file_df)} unique read length/frame/offset combinations"
        )

    logging.info("Sample sheet parsing completed successfully")
    return unique_files, file_info


# --------------------------------------------------


def validate_input_files(
    input_dir: str, file_names: List[str]
) -> Tuple[List[str], List[str]]:
    """
    Validate that all files in the sample sheet exist in the input directory.

    Parameters
    ----------
        input_dir: Directory containing input BAM files
        file_names: List of file names from the sample sheet

    Returns
    -------
        A tuple containing:
        - List of valid file paths
        - List of missing files
    """

    logging.info(f"Validating {len(file_names)} input files in directory: {input_dir}")

    input_path = Path(input_dir)

    if not input_path.exists() or not input_path.is_dir():
        logging.error(f"Input directory not found or not a directory: {input_path}")
        raise NotADirectoryError(
            f"Input directory not found or not a directory: {input_path}"
        )

    valid_files = []
    missing_files = []

    for file_name in file_names:
        file_path = input_path / file_name
        if file_path.exists() and file_path.is_file():
            valid_files.append(str(file_path))
            logging.debug(f"File exists: {file_path}")
        else:
            missing_files.append(file_name)
            logging.error(f"File not found: {file_path}")
            sys.exit(1)

    logging.info(
        f"Validation complete: {len(valid_files)} valid files, {len(missing_files)} missing files"
    )

    return valid_files, missing_files


# --------------------------------------------------


def asite_pos(df: pl.DataFrame, offset_file: pl.DataFrame) -> pl.DataFrame:
    """
    Join the DataFrame with an offset file and compute the shifted position (pos + offset).

    Parameters
    ----------
    df : pl.DataFrame
        Input DataFrame that must include 'read_length', 'frame', 'pos', and 'transcript_id' columns.

    offset_file : str
        Path to a tab-delimited offset file with columns: read_length, frame, offset.

    Returns
    -------
    pl.DataFrame
        A DataFrame with "transcript_id" and 'a_site_pos' columns.
    """

    tmp_df = (
        df.join(offset_file, on=["read_length", "frame"], how="inner")
        .with_columns(a_site_pos=(pl.col("pos") + pl.col("offset")))
        .select(["transcript_id", "a_site_pos"])
    )

    return tmp_df


# --------------------------------------------------


def count_reads_on_pos(df: pl.DataFrame) -> pl.DataFrame:
    """
    Count reads on A-site positions by grouping by transcript_id and a_site_pos.

    Parameters
    ----------
    df : pl.DataFrame
        Input DataFrame that must include 'transcript_id' and 'a_site_pos' columns.

    Returns
    -------
    pl.DataFrame
        A DataFrame with "transcript_id", 'a_site_pos' and "reads" columns.
    """

    tmp_df = df.group_by("transcript_id", "a_site_pos").agg(pl.len().alias("reads"))

    return tmp_df


# --------------------------------------------------


def filter_and_comp_ep_pos(
    df: pl.DataFrame,
    fasta_index: pl.DataFrame,
    excl_start: int = 60,
    excl_end: int = 60,
) -> pl.DataFrame:
    """
    Filter transcripts based on length criteria and calculate ribosome site positions.

    This function expects CDS to be extended by 18 nt on both sides. It excludes the specified
    number of nucleotides from the start and end of transcripts, and keeps only positions that
    lie within the remaining CDS region. It also calculates E-site and P-site positions from the A-site.

    Parameters
    ----------
    df : pl.DataFrame
        Input DataFrame that must include 'transcript_id', 'a_site_pos', and 'reads' columns.
    fasta_index : pl.DataFrame or str
        Path to a CSV file or DataFrame containing transcript information with columns 'transcript_id',
        'length', and other metadata.
    excl_start : int, default 60
        Number of nucleotides to exclude from the start of the CDS (after the initial 18 nt extension).
    excl_end : int, default 60
        Number of nucleotides to exclude from the end of the CDS (before the final 18 nt extension).

    Returns
    -------
    pl.DataFrame
        DataFrame with 'transcript_id', 'length', 'e_site_pos', 'p_site_pos', 'a_site_pos',
        and 'reads' columns. Only includes transcripts with at least 3 nucleotides remaining
        after trimming and only positions within the valid range.
    """
    trim_nt = 18 + excl_start + excl_end + 18
    fasta_index = pl.read_csv(fasta_index)

    tmp_df = (
        fasta_index.select(["transcript_id", "length"])
        .with_columns(
            start=(18 + excl_start),
            end=(pl.col("length") - (18 + excl_end)),
            length=(pl.col("length") - trim_nt),
        )
        .filter(pl.col("length") >= 3)
        .join(df, on="transcript_id", how="inner")
        .filter(
            pl.col("a_site_pos") >= pl.col("start"),
            pl.col("a_site_pos") <= pl.col("end"),
        )
        .select(["transcript_id", "length", "a_site_pos", "reads"])
        .with_columns(
            p_site_pos=(pl.col("a_site_pos") - 3), e_site_pos=(pl.col("a_site_pos") - 6)
        )
    )

    return tmp_df


# --------------------------------------------------


def count_cds_reads(df: pl.DataFrame, fasta_index: str):
    """
    Aggregates reads on CDS of each transcripts.

    Parameters
    ----------
    df : pl.DataFrame
        Input DataFrame containing 'transcript_id', 'pos', and 'reads' columns.
    fasta_index : pl.DataFrame or str
        Path to a CSV file or DataFrame containing transcript sequences with columns
        'transcript_id', 'sequence' and 'length.

    Returns
    -------
    pl.DataFrame
        DataFrame with 'transcript_id', 'length', 'reads' and 'tpm' columns
    """

    fasta_index = pl.read_csv(fasta_index)

    filtered_df = df.group_by("transcript_id", "pos").agg(reads=pl.len())

    read_count = (
        fasta_index.select(["transcript_id", "length"])
        .join(filtered_df, on="transcript_id", how="inner")
        .filter(pl.col("pos") >= 18, pl.col("pos") <= (pl.col("length") - 18))
        .group_by("transcript_id", "length")
        .agg(total_reads=pl.sum("reads"))
        .with_columns(
            # Calculate Reads per Kilobase (RPK)
            rpk=pl.col("total_reads") / (pl.col("length") / 1_000)
        )
        .with_columns(
            # Calculate per million scaling factor (RPK / 10^6)
            scaling_factor=(pl.col("rpk").sum() / 10**6)
        )
        .with_columns(
            # Calculate TPM
            tpm=pl.col("rpk") / pl.col("scaling_factor")
        )
        .select(["transcript_id", "length", "total_reads", "tpm"])
    )

    return read_count


# --------------------------------------------------


def identify_codons(df: pl.DataFrame, fasta_index: str) -> pl.DataFrame:
    """
    Identify codons at E, P, and A-sites for each position in the transcripts.

    Parameters
    ----------
    df : pl.DataFrame
        Input DataFrame containing 'transcript_id', 'length', 'e_site_pos', 'p_site_pos',
        'a_site_pos', and 'reads' columns.
    fasta_index : pl.DataFrame or str
        Path to a CSV file or DataFrame containing transcript sequences with columns
        'transcript_id' and 'sequence'.

    Returns
    -------
    pl.DataFrame
        DataFrame with 'transcript_id', 'length', 'e_site_pos', 'e_codon', 'p_site_pos',
        'p_codon', 'a_site_pos', 'a_codon', and 'reads' columns, where the codon columns
        contain the 3-nucleotide sequences at each ribosome site.
    """

    fasta_index = pl.read_csv(fasta_index)

    tmp_df = (
        fasta_index.select(["transcript_id", "sequence"])
        .join(df, on="transcript_id", how="inner")
        .with_columns(
            e_codon=pl.col("sequence").str.slice(pl.col("e_site_pos"), 3),
            p_codon=pl.col("sequence").str.slice(pl.col("p_site_pos"), 3),
            a_codon=pl.col("sequence").str.slice(pl.col("a_site_pos"), 3),
        )
        .select(
            [
                "transcript_id",
                "length",
                "e_site_pos",
                "e_codon",
                "p_site_pos",
                "p_codon",
                "a_site_pos",
                "a_codon",
                "reads",
            ]
        )
    )

    return tmp_df


# --------------------------------------------------


def calc_enrichment_score(df: pl.DataFrame) -> pl.DataFrame:
    """
    Calculate enrichment score for each position in the transcripts.

    Implements the enrichment score calculation as described in Hussmann J et al,
    PLOS Genetics, 2015 (supplementary data equations 1 and 2). The enrichment
    score is calculated as the ratio of the observed read count to the mean read
    count per position for that transcript.

    Parameters
    ----------
    df : pl.DataFrame
        Input DataFrame containing codon positions and read count information.
        Must contain at minimum 'transcript_id', 'e_site_pos', 'e_codon', 'p_site_pos',
        'p_codon', 'a_site_pos', 'a_codon', 'reads', and 'length' columns.

    Returns
    -------
    pl.DataFrame
        DataFrame with "transcript_id", "e_site_pos", "e_codon", "p_site_pos", "p_codon",
        "a_site_pos", "a_codon", "mean_read_count", "reads", "length", and "enrichment_score"
        columns, where enrichment_score is calculated as reads/mean_read_count.

    Notes
    -----
    The function first calculates the mean read count per position for each transcript by
    dividing the total read count by the transcript length, then calculates the enrichment
    score for each position.
    """

    mean_read_counts = (
        df.group_by("transcript_id")
        .agg(
            total_reads=pl.sum("reads"),
            transcript_length=pl.first("length"),
        )
        .with_columns(
            mean_read_count=(pl.col("total_reads") / pl.col("transcript_length"))
        )
        .filter(pl.col("mean_read_count") > 0.1)
        .select(["transcript_id", "mean_read_count"])
    )

    enrichment_df = (
        df.join(mean_read_counts, on="transcript_id", how="inner")
        .with_columns(enrichment_score=(pl.col("reads") / pl.col("mean_read_count")))
        .select(
            [
                "transcript_id",
                "e_site_pos",
                "e_codon",
                "p_site_pos",
                "p_codon",
                "a_site_pos",
                "a_codon",
                "mean_read_count",
                "reads",
                "length",
                "enrichment_score",
            ]
        )
    )

    return enrichment_df


# --------------------------------------------------


def calc_enrichment_offset(
    df: pl.DataFrame,
    fasta_index: str,
    codon: str,
    excl_start: int,
    excl_end: int,
    offset_upstream: int,
    offset_downstream: int,
) -> pl.DataFrame:
    """
    Calculate median enrichment scores at codon-specific offsets for a given codon.

    This function takes a DataFrame of per-codon enrichment scores and a FASTA index file,
    and computes the median enrichment score for each codon offset (in codon units)
    relative to a reference codon (`codon`) at the A-site position. It filters out
    positions too close to transcript ends according to exclusion windows.

    Parameters
    ----------
    df : pl.DataFrame
        Input DataFrame containing, at minimum, the following columns:
          - "transcript_id" (str): Transcript identifier.
          - "a_site_pos" (int): Nucleotide position of the A-site codon.
          - "a_codon" (str): Codon at the A-site.
          - "enrichment_score" (float): Enrichment score for that codon position.
    fasta_index : str
        Path to a CSV file containing a FASTA index with at least:
          - "transcript_id" (str): Transcript identifier.
          - "sequence" (str): Full transcript nucleotide sequence.
          - "length" (int): Length of the transcript sequence.
    excl_start : int
        Number of nucleotides upstream of the A-site to exclude from analysis.
    excl_end : int
        Number of nucleotides downstream of the A-site to exclude from analysis.
    offset_upstream : int
        Maximum codon offset (in codons) upstream (negative direction) to include.
    offset_downstream : int
        Maximum codon offset (in codons) downstream (positive direction) to include.
    codon : str
        The reference codon at the A-site to center the offset analysis around.

    Returns
    -------
    pl.DataFrame
        A DataFrame with one row per codon offset, containing:
          - "offset" (int): Codon offset relative to the reference A-site codon.
          - "med_enrichment_score" (float): Median enrichment score across all transcripts
            at the given offset.
          - "pos0_codon" (str): The reference codon (same for all rows, equal to `codon`).

    Notes
    -----
    - Positions are measured in nucleotides; codon offsets are multiplied by 3 to
      translate to nucleotide offsets.
    - Transcripts for which the calculated offset position would fall within
      `excl_start` nt of the 5′ end or within `excl_end` nt of the 3′ end
      (after trimming) are automatically filtered out.
    - The FASTA index CSV is read into memory; ensure it contains valid sequences
      matching the transcript IDs in `df`.

    """

    trim_nt = 18 + excl_start + excl_end

    fasta_index = pl.read_csv(fasta_index)

    # Create a range of offsets from -offset_range to +offset_range
    offsets = pl.DataFrame(
        {"offset": list(range(-offset_upstream, offset_downstream + 1))}
    )

    subset_df = df.filter(pl.col("a_codon") == codon).select(
        ["transcript_id", "a_site_pos", "a_codon", "enrichment_score"]
    )

    offset_enrichment_df = (
        subset_df.join(offsets, how="cross")
        .with_columns(
            # Calculate offset position (in nucleotides: 3 * codon offset)
            offset_pos=(pl.col("a_site_pos") + (pl.col("offset") * 3))
        )
        .join(fasta_index, on="transcript_id", how="inner")
        # filter to ensure offset positions are within valid range
        .filter(
            (pl.col("offset_pos") >= (18 + excl_start)),
            (pl.col("offset_pos") <= (pl.col("length") + trim_nt - (18 + excl_end))),
        )
        .with_columns(
            offset_codon=pl.col("sequence").str.slice(pl.col("offset_pos"), 3)
        )
        .select(["transcript_id", "offset_pos", "offset", "offset_codon"])
        .rename({"offset_pos": "a_site_pos"})
        .join(df, how="inner", on=["transcript_id", "a_site_pos"])
        .select(
            [
                "transcript_id",
                "a_site_pos",
                "offset",
                "offset_codon",
                "enrichment_score",
            ]
        )
        .group_by("offset")
        .agg(med_enrichment_score=pl.median("enrichment_score"))
        .with_columns(pos0_codon=pl.lit(codon))
    )

    return offset_enrichment_df


# --------------------------------------------------


def run(args):
    """
    Execute the main analysis pipeline for Ribomala.

    Parameters
    ----------
    args : argparse.Namespace
        Parsed command-line arguments containing:
        - input (str): Path to input directory containing BAM files.
        - output (str): Path to output directory.
        - exclstart (int): Start position for exclusion in enrichment analysis.
        - exclend (int): End position for exclusion in enrichment analysis.
        - samples (str): Path to sample sheet CSV file.
        - txcsv (str): Path to transcriptome FASTA index CSV file.
        - upstream (int): Number of upstream positions for enrichment calculation.
        - downstream (int): Number of downstream positions for enrichment calculation.
        - codon (str): Comma-separated list of codons to analyze.

    Notes
    -----
    This function validates inputs, processes each BAM file for read alignment
    and enrichment, and writes the results to the output directory.
    Intermediate steps include read shifting, codon identification,
    and enrichment score calculations.
    """

    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s - %(levelname)s - %(message)s",
    )

    input_dir = Path(args.input)
    output_dir = Path(args.output)
    excl_start = args.exclstart
    excl_end = args.exclend
    sample_sheet = args.samples
    indexed_fasta = args.txcsv
    offset_upstream = args.upstream
    offset_downstream = args.downstream

    # Split codons by comma and strip whitespace
    codons = [cod.strip() for cod in args.codon.split(",") if cod.strip()]

    # Create output directory if it doesn't exist
    output_dir.mkdir(parents=True, exist_ok=True)

    try:
        unique_files, file_info = parse_sample_sheet(sample_sheet_path=sample_sheet)
        logging.info(f"Unique files to process: {unique_files}")

        for file_name, info_df in file_info.items():
            logging.info(f"Details for {file_name}:\n{info_df}")

    except Exception as e:
        logging.error(f"Error while parsing sample sheet: {e}")
        sys.exit(1)

    try:
        valid_files, missing_files = validate_input_files(input_dir, unique_files)
        if missing_files:
            logging.error(f"The following files are missing: {missing_files}")
            sys.exit(1)
        else:
            logging.info("All input files were found successfully.")

    except Exception as e:
        logging.error(f"Error while validating input files: {e}")
        sys.exit(1)

    # Analysis begins
    try:
        for bam_file in unique_files:
            sample_name = bam_file
            logging.info(f"Processing sample file: {bam_file}")
            bam_offsets = file_info[bam_file]

            output_name = Path(bam_file).stem

            bam_file = input_dir / bam_file
            logging.info("Processing BAM file")
            bam_file = qc.extract_bam_info(bam_file)

            logging.info("Calculating reads count on CDS")
            count_df = count_cds_reads(df=bam_file, fasta_index=indexed_fasta)

            logging.info("Calculating ribosome reading frame")
            bam_file = qc.calculate_frame(bam_file)
            logging.info("Shifting reads to the provided offset position")
            enrichment_df = asite_pos(df=bam_file, offset_file=bam_offsets)
            logging.info("Counting reads on offset codon")
            enrichment_df = count_reads_on_pos(enrichment_df)
            logging.info("Computing positions relative to provided offset")
            enrichment_df = filter_and_comp_ep_pos(
                df=enrichment_df,
                fasta_index=indexed_fasta,
                excl_start=excl_start,
                excl_end=excl_end,
            )
            logging.info("Identifying codons at each position from the offset")
            enrichment_df = identify_codons(df=enrichment_df, fasta_index=indexed_fasta)
            logging.info("Calculating enrichment scores")
            enrichment_df = calc_enrichment_score(enrichment_df)

            # Enrichment around the codons of choice
            for codon in codons:
                logging.info(f"Calculating offset enrichment for codon: {codon}")
                offset_enrichment_df = calc_enrichment_offset(
                    df=enrichment_df,
                    fasta_index=indexed_fasta,
                    excl_start=excl_start,
                    excl_end=excl_end,
                    offset_upstream=offset_upstream,
                    offset_downstream=offset_downstream,
                    codon=codon,
                )

                offset_out_path = Path(output_dir) / "offset_enrichment"
                offset_out_path.mkdir(parents=True, exist_ok=True)
                output_path = (
                    offset_out_path
                    / f"{output_name}_{codon}_offset_enrichment_scores.csv"
                )
                logging.info(f"Writing enrichment scores to {output_path}")
                offset_enrichment_df.write_csv(output_path)

            enrichment_out_path = Path(output_dir) / "enrichment_epa_site"
            enrichment_out_path.mkdir(parents=True, exist_ok=True)

            enrichment_output_path = (
                Path(enrichment_out_path) / f"{output_name}_enrichment_scores.csv"
            )

            logging.info(f"Writing enrichment scores to {enrichment_output_path}")
            enrichment_df.write_csv(enrichment_output_path)

            count_out_path = Path(output_dir) / "cds_counts"
            count_out_path.mkdir(parents=True, exist_ok=True)

            count_output_path = Path(count_out_path) / f"{output_name}_counts.csv"

            logging.info(f"Writing transcripts count to {count_output_path}")
            count_df.write_csv(count_output_path)

            logging.info(f"Finished processing sample: {sample_name}")

    except Exception as e:
        logging.error(f"Encountered error: {e}")
        sys.exit(1)
