#! /usr/bin/env python3

"""
Ribosome profiling data quality control functionality for Ribomala. The module performs the following tasks:

- processes BAM files listed in a sample list.
- extracts transcript ID, 5' position and read length
- calculates reading frames
- checks periodicity and ribosome reading frame
- Plots the data and saves them to the output directory.
"""

import logging
from pathlib import Path
import plotly.express as px
import polars as pl
import pysam
import sys

# --------------------------------------------------


def read_samples_file(sample_list_file: Path) -> list:
    """
    Read the sample list file or a single sample and return a list of sample filenames.

    Parameters
    ----------
    sample_list_file : Path
        Path to the sample list file or a single sample name.

    Returns
    -------
    list
        A list of sample file names.

    Notes
    -----
    If sample_list_file is an existing file, read sample names from it;
    otherwise, treat the provided value as a single sample name.
    """
    try:
        if sample_list_file.exists() and sample_list_file.is_file():
            with sample_list_file.open("r") as f:
                samples = [line.strip() for line in f if line.strip()]
            logging.info(f"Found {len(samples)} samples in {sample_list_file}")
        else:
            # Assume the argument itself is a single sample name.
            samples = [str(sample_list_file)]
            logging.info(f"Using single sample: {samples[0]}")
        return samples
    except Exception as e:
        logging.error(f"Error processing sample list input {sample_list_file}: {e}")
        raise


# --------------------------------------------------


def extract_bam_info(bam_file: str) -> pl.DataFrame:
    """
    Extract transcript ID, 5' position, and read length from a BAM file.

    Parameters
    ----------
    bam_file : str
        Path to the BAM file.

    Returns
    -------
    pl.DataFrame
        A Polars DataFrame with columns: transcript_id, pos, and read_length.

    Raises
    ------
    Exception
        If the BAM file cannot be processed.
    """
    bam_info = []
    try:
        with pysam.AlignmentFile(bam_file, "rb") as bam:
            logging.info(f"Reading: {bam_file}")
            for read in bam:
                transcript_id = read.reference_name
                pos = read.reference_start
                read_length = read.query_length
                bam_info.append((transcript_id, pos, read_length))
        logging.info(f"Extracted {len(bam_info)} records from {bam_file}")
    except Exception as e:
        logging.error(f"Failed to process {bam_file}: {e}")
        raise
    return pl.DataFrame(
        bam_info, schema=["transcript_id", "pos", "read_length"], orient="row"
    )


# --------------------------------------------------


def calculate_frame(df: pl.DataFrame) -> pl.DataFrame:
    """
    Calculate the reading frame as pos modulo 3.

    Parameters
    ----------
    df : pl.DataFrame
        Input DataFrame with a 'pos' column.

    Returns
    -------
    pl.DataFrame
        A DataFrame with an added 'frame' column (pos % 3).
    """
    return df.with_columns((pl.col("pos") % 3).alias("frame"))


# --------------------------------------------------


def check_periodicity(
    df: pl.DataFrame,
    output_dir: Path,
    sample_name: str,
    min_read_length: int = 28,
    max_read_length: int = 33,
) -> None:
    """
    Check periodicity in ribosome profiling data.

    Parameters
    ----------
    df : pl.DataFrame
        Input DataFrame with 'pos', 'read_length', and 'frame' columns.
    output_dir : Path
        Path to the output directory.
    sample_name : str
        Name of the sample (used in output filenames).
    min_read_length : int, optional
        Minimum read length to consider, by default 28.
    max_read_length : int, optional
        Maximum read length to consider, by default 33.

    Notes
    -----
    This function filters reads, groups by position, read_length, and frame,
    then generates and exports a Plotly bar plot and CSV data to the output directory.
    """
    periodicity_df = (
        df.filter(
            (pl.col("pos") <= 48)
            & (pl.col("read_length") >= min_read_length)
            & (pl.col("read_length") <= max_read_length)
        )
        .group_by(["pos", "read_length", "frame"])
        .agg(pl.count().alias("num_reads"))
        .with_columns((pl.col("pos") - 18).alias("pos"), pl.col("frame").cast(pl.Utf8))
    )

    read_length_order = list(range(min_read_length, max_read_length + 1))
    custom_colors = {"0": "#e7298a", "1": "#e6ab02", "2": "#1b9e77"}

    fig = px.bar(
        periodicity_df,
        x="pos",
        y="num_reads",
        color="frame",
        facet_col="read_length",
        facet_col_wrap=2,
        color_discrete_map=custom_colors,
        category_orders={"frame": ["0", "1", "2"], "read_length": read_length_order},
        labels={
            "pos": "Distance from start codon (nt)",
            "num_reads": "Number of reads",
            "frame": "Reading frame",
        },
        template="simple_white",
    )
    fig.update_xaxes(range=[-18, 30], tickvals=list(range(-18, 31, 3)))
    fig.update_yaxes(matches=None, showticklabels=True)

    html_file = output_dir / f"{sample_name}_periodicity.html"
    csv_file = output_dir / f"{sample_name}_periodicity_data.csv"
    fig.write_html(str(html_file))
    periodicity_df.write_csv(str(csv_file))
    logging.info(f"Periodicity plot saved to {html_file} and data to {csv_file}")


# --------------------------------------------------


def check_frame_dist(
    df: pl.DataFrame,
    output_dir: Path,
    sample_name: str,
    min_read_length: int = 28,
    max_read_length: int = 33,
) -> None:
    """
    Check ribosome reading frame distribution.

    Parameters
    ----------
    df : pl.DataFrame
        Input DataFrame with 'read_length' and 'frame' columns.
    output_dir : Path
        Path to the output directory.
    sample_name : str
        Name of the sample (used in output filenames).
    min_read_length : int, optional
        Minimum read length to consider, by default 28.
    max_read_length : int, optional
        Maximum read length to consider, by default 33.

    Notes
    -----
    This function groups data by read_length and frame,
    then generates and exports a Plotly bar plot and CSV data to the output directory.
    """
    frame_df = (
        df.filter(
            (pl.col("read_length") >= min_read_length)
            & (pl.col("read_length") <= max_read_length)
        )
        .group_by(["read_length", "frame"])
        .agg(pl.count().alias("num_reads"))
        .with_columns(pl.col("frame").cast(pl.Utf8))
    )

    custom_colors = {"0": "#e7298a", "1": "#e6ab02", "2": "#1b9e77"}
    fig = px.bar(
        frame_df,
        x="read_length",
        y="num_reads",
        color="frame",
        barmode="stack",
        color_discrete_map=custom_colors,
        category_orders={"frame": ["0", "1", "2"]},
        labels={
            "read_length": "Read length",
            "num_reads": "Number of reads",
            "frame": "Reading frame",
        },
        template="simple_white",
    )

    html_file = output_dir / f"{sample_name}_frame.html"
    csv_file = output_dir / f"{sample_name}_frame_data.csv"
    fig.write_html(str(html_file))
    frame_df.write_csv(str(csv_file))
    logging.info(f"Frame distribution plot saved to {html_file} and data to {csv_file}")


# --------------------------------------------------


def run(args):
    """
    Run the QC mode with the given arguments.

    Parameters
    ----------
    args : argparse.Namespace
        Parsed command-line arguments with the following attributes:
        - input: Path to the input directory containing BAM files
        - output: Path to the output directory for results
        - samples: Path to the samples list file or a single sample name
        - min: Minimum read length to consider
        - max: Maximum read length to consider

    Returns
    -------
    None

    Notes
    -----
    This function configures logging and executes the QC process.
    If QC fails, the program will exit with code 1.
    """

    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s - %(levelname)s - %(message)s",
    )

    input_dir = Path(args.input)
    output_dir = Path(args.output)
    min_rpf_len = args.min
    max_rpf_len = args.max

    # Create output directory if it doesn't exist
    output_dir.mkdir(parents=True, exist_ok=True)

    try:
        samples = read_samples_file(sample_list_file=Path(args.samples))

        for sample in samples:
            input_sample = input_dir / Path(sample)

            out_sample_name = output_dir / input_sample.stem

            logging.info(f"Processing sample: {sample}")
            bam_df = extract_bam_info(bam_file=input_sample)

            logging.info("Calculating reading frame.")
            bam_df = calculate_frame(df=bam_df)

            logging.info("Plotting periodicity.")
            check_periodicity(
                df=bam_df,
                output_dir=output_dir,
                sample_name=out_sample_name,
                min_read_length=min_rpf_len,
                max_read_length=max_rpf_len,
            )

            logging.info("Plotting reading frame distribution.")
            check_frame_dist(
                df=bam_df,
                output_dir=output_dir,
                sample_name=out_sample_name,
                min_read_length=min_rpf_len,
                max_read_length=max_rpf_len,
            )

            logging.info(f"Finished processing sample: {sample}")

    except Exception as e:
        logging.error(f"Processing BAM for QC failed: {e}")
        sys.exit(1)
