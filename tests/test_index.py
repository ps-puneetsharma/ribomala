#! /usr/bin/env python3

import shutil
from pathlib import Path

import pysam
import polars as pl
import pytest

from ribomala.index import index_fasta


def test_index_fasta_creates_index_and_csv(tmp_path):
    # Prepare the FASTA file in the temporary directory
    src_file = Path(__file__).parent.parent / "data" / "toy_transcriptome.fa"
    fasta_file = tmp_path / "toy_transcriptome.fa"
    shutil.copy(src_file, fasta_file)

    # Run the indexing function on the temporary FASTA file.
    indexed_fasta = index_fasta(str(fasta_file))

    # Check that the index file (.fai) has been created.
    index_file = fasta_file.with_name(fasta_file.name + ".fai")
    assert index_file.exists(), "Index file was not created."

    # Check that the CSV file has been created.
    csv_file = fasta_file.parent / f"{fasta_file.stem}.csv"
    assert csv_file.exists(), "CSV file was not created."

    # Validate the CSV content using Polars.
    df = pl.read_csv(csv_file)
    expected_columns = {"transcript_id", "sequence", "length"}
    assert expected_columns.issubset(set(df.columns)), (
        "CSV file does not contain all expected columns."
    )
    assert len(df) >= 1, "CSV file is empty."

    # Confirm that the returned object is an instance of pysam.FastaFile.
    assert isinstance(indexed_fasta, pysam.FastaFile), (
        "Returned object is not a pysam.FastaFile instance."
    )
