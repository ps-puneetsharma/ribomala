#! /usr/bin/env python3

"""
Indexing functionality for Ribomala.
This module handles indexing of the transcriptome fasta file.
"""

import logging
from pathlib import Path
import polars as pl
import pysam
import sys

# --------------------------------------------------


def index_fasta(fasta_file: str) -> pysam.FastaFile:
    """
    Indexes a transcriptome FASTA file if it is not already indexed and returns a `pysam.FastaFile` object.

    Parameters
    ----------
    fasta_file : str
        Path to the transcriptome FASTA file.

    Returns
    -------
    pysam.FastaFile
        An indexed `pysam.FastaFile` object for the given transcriptome.

    Notes
    -----
    This function checks if an index file (`.fai`) exists for the given FASTA file.
    If the index file is missing, it generates one using `pysam.faidx`.
    The function also creates a CSV file with transcript information.

    """

    fasta_path = Path(fasta_file)
    logging.info(f"Transcriptome to be indexed: {fasta_file}")
    index_path = fasta_path.with_name(fasta_path.name + ".fai")
    output_name = fasta_path.stem

    try:
        logging.info(f"Reading FASTA file: {fasta_file}")
        if not index_path.exists():
            pysam.faidx(str(fasta_path))

        indexed_fasta = pysam.FastaFile(str(fasta_path))
        transcripts = indexed_fasta.references
        sequences = [indexed_fasta.fetch(tid) for tid in transcripts]

        logging.info(f"Found {len(transcripts)} transcripts in FASTA file")

        # Create DataFrame
        fasta_df = pl.DataFrame({"transcript_id": transcripts, "sequence": sequences})

        # Add length column (+1 because we want 1-based counting)
        fasta_df = fasta_df.with_columns(
            length=(pl.col("sequence").str.len_chars() + 1),
        )

        # Write processed file as CSV
        output_path = fasta_path.parent / f"{output_name}.csv"
        logging.info(f"Writing CSV to {output_path}")
        fasta_df.write_csv(output_path)
        logging.info(f"{fasta_file} successfully indexed!")

    except Exception as e:
        logging.error(f"Error processing FASTA file: {e}")
        raise
    return indexed_fasta


# --------------------------------------------------


def run(args):
    """
    Run the index mode with the given arguments.

    Parameters
    ----------
    args : argparse.Namespace
        Command-line arguments for the index mode.
        Should contain a 'fasta' attribute with the path to the FASTA file.

    Returns
    -------
    None

    Notes
    -----
    This function configures logging and executes the indexing process.
    If indexing fails, the program will exit with code 1.
    """

    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s - %(levelname)s - %(message)s",
    )

    try:
        index_fasta(args.fasta)
    except Exception as e:
        logging.error(f"Indexing failed: {e}")
        sys.exit(1)
