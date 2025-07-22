# Ribosome profiling data analysis workflow

The workflow assumes that the reads look like:

`NNN <RPF> NNN NNN ADAPTOR`

Where `NNN` are random nucleotides around the ribosome protected fragment (RPF).

## Preparing genomes and annotations

### “Undesired RNA” annotation

Example:

Organism of interest: Homo sapiens (Human)

A FASTA file containing rRNA, snRNA and snoRNA sequences can be generated from [Ensembl biomart](https://www.ensembl.org/biomart/martview/) by clicking-on/selecting the following options in order:

- Database: Ensembl Genes (Version)
- Dataset: Human genes (Version)
- Click “Filters”
- Expand the “Gene” tab
- Select Transcript type: rRNA, snRNA, snoRNA
- Click “Attributes”
- Select “Sequences” (Radio button) and expand the “Sequences” tab
- Select “Unspliced (gene)”
- Expand the “Header information” tab
- Select “Gene stable ID”
- Click “Count”
- Click “Results”
- Download the FASTA file

tRNA FASTA sequences can be obtained from [Genomic tRNA database](http://gtrnadb.ucsc.edu/)

- Under “Links to Most Viewed Genomes”, select “Homo sapiens (GRCh38/hg38)”
- On the left hand side, select “FASTA Seqs”
- Download “High confidence tRNA sequences: hg38-tRNAs.fa”

Use [NCBI Nucleotide module](https://www.ncbi.nlm.nih.gov/nuccore/?term=) to procure 45S, 28S, 18S, 5.8S, 5S rRNA sequences.

- Search: biomol_rRNA[prop] AND “Homo sapiens”[Organism].
- On the top right, click “Send to”
- Select “File”
- Select Format: FASTA
- Click “Create File”

The above procedure will download all rRNA sequences (806 at the time of writing) that were shown in the result.

Once you have the above three FASTA files (for sno-, sn- and r-RNA; tRNA; and rRNA), concatenate them into one single FASTA file. This is your annotation for the “undesired” RNA.

```bash

# File from ENSEMBL: GRCh38_r_sno_sn_RNA_ENSEMBL.txt
# File from NCBI: GRCh38_rRNA_NCBI.txt
# File from tRNAdb: GRCh38_tRNA_tRNADB.txt

# Catenate the files
cat *.txt > non_target_rna.fa

```


### Transcriptome annotation

A FASTA file containing sequences of principle splice isoforms that are **extended by 18 nt in 5’- and 3’- ends** can be generated from [Ensembl biomart](https://www.ensembl.org/biomart/martview/) by clicking-on/selecting the following options in order:
 by clicking-on/selecting the following options in order:

- Database: Ensembl Genes (Version)
- Dataset: Human genes (Version)
- Click “Filters”
- Expand the “Gene” tab
- Select Transcript type: lncRNA (optional) , protein_coding
- Check Source (transcript): ensembl_havana
- Check Ensembl Canonical: Only
- Check APPRIS annotation
- Click “Attributes”
- Select “Sequences” (Radio button) and expand the “Sequences” tab
- Select “Coding sequence”
- Check Upstream flank and provide the value “18”
- Check the Downstream flank and provide the value “18”
- Expand the “Header information” tab
- Select “Transcript stable ID”
- Click “Count”
- Click “Results”
- Download the FASTA file (called GRCh38_ensembl_cds_plus18_annot.fa in this document)


## Indexing

```bash

dir_in="/dir/path"
gen_in="/dir/path"
THREADS=128

# Safety first
set -e

date +"%d-%m-%Y %T"
echo "-------------------------"
echo "Indexing non-target RNA fasta using Bowtie"
echo "-------------------------"

bowtie-build \
--threads $THREADS \
"$gen_in"/non_target_rna.fa \
"$gen_in"/non_target_rna


date +"%d-%m-%Y %T"
echo "-------------------------"
echo "Indexing CDS fasta using Bowtie"
echo "-------------------------"

bowtie-build \
--threads $THREADS \
"$gen_in"/GRCh38_ensembl_cds_plus18_annot.fa \
"$gen_in"/GRCh38_ensembl_cds_plus18_annot

date +"%d-%m-%Y %T"
echo "-------------------------"
echo "Indexing DONE!"
echo "-------------------------"




```

## Removing 5' random nt and 3' adaptors from fastq files

```bash

dir_in="/dir/path"
gen_in="/dir/path"
THREADS=128

for file in "$dir_in"/*.fastq.gz; do

    base=$(basename "$file" .fastq.gz)

    date +"%d-%m-%Y %T"
    echo "-------------------------"
    echo "Cutadapt working on "${base}""
    echo "-------------------------"

    # Remove 5' random nt and adaptor
    cutadapt \
    -j $THREADS \
    --quality-cutoff 30 \
    --minimum-length 25 \
    --cut 3 \
    --no-indels \
    --discard-untrimmed \
    --adapter CTGTAGGCACCATCAAT \
    --output "$dir_in"/${base}_no_adapter.fastq.gz \
    "$file" \
    1> "$dir_in"/${base}_adapter_removal_cutadapt_log.txt

    # Remove 3' random nt
    cutadapt \
    -j $THREADS \
    --quality-cutoff 30 \
    --minimum-length 25 \
    --cut -6 \
    --output "$dir_in"/${base}_clean.fastq.gz \
    "$dir_in"/${base}_no_adapter.fastq.gz \
    1> "$dir_in"/${base}_trim_cutadapt_log.txt

    # remove temporary file
    rm "$dir_in"/${base}_no_adapter.fastq.gz

    date +"%d-%m-%Y %T"
    echo "-------------------------"
    echo "Cutadapt done processing "${base}""
    echo "-------------------------"

done

```

## Mapping

```bash

dir_in="/dir/path"
gen_in="/dir/path"
THREADS=128

# Map to non-target
for file in "$dir_in"/*_clean.fastq.gz; do

    base=$(basename "$file" _clean.fastq.gz)

    echo "-------------------------"
    echo "Bowtie aligning "$base" to non-target RNA"
    echo "-------------------------"

    bowtie \
    --threads $THREADS \
    --chunkmbs 320 \
    --time \
    --best \
    --sam \
    --un "$dir_in"/${base}_no_ncrna.fastq \
    -x "$gen_in"/non_target_rna \
    -q "$file" \
    2> "$dir_in"/${base}_non_target_rna_map_log.txt \
    1> /dev/null

    echo "-------------------------"
    echo "Compressing fastq"
    echo "-------------------------"

    pigz -p64 "$dir_in"/${base}_no_ncrna.fastq

    echo "-------------------------"
    echo "Bowtie DONE aligning "$base" to non-target RNA"
    echo "-------------------------"


done

# Map to CDS
for file in "$dir_in"/*_no_ncrna.fastq.gz; do

    base=$(basename "$file" _no_ncrna.fastq.gz)

    echo "-------------------------"
    echo "Bowtie aligning "$base" to CDS"
    echo "-------------------------"

    bowtie \
    --threads "$THREADS" \
    --chunkmbs 320 \
    --time \
    --best \
    -v 1 \
    -m 1  \
    --norc \
    --strata \
    --sam \
    -x "$gen_in"/GRCh38_ensembl_cds_plus18_annot \
    -q "$file" \
    2> "$dir_in"/${base}_cds_map_log.txt | \
    samtools view \
    -@ 10 \
    -h -bS -F 4 > "$dir_in"/${base}_cds_unsorted.bam

    echo "-------------------------"
    echo "Samtool sorting: "${base}_cds_unsorted.bam""
    echo "-------------------------"

    # Sort the BAM file
    samtools sort \
    -@ $THREADS \
    -o "$dir_in"/${base}_cds_sorted.bam "$dir_in"/${base}_cds_unsorted.bam

    echo "-------------------------"
    echo "Samtool indexing: "${base}_cds_sorted.bam""
    echo "-------------------------"

    # Index the BAM file
    samtools index \
    -@ "$THREADS" \
    "$dir_in"/${base}_cds_sorted.bam

    echo "-------------------------"
    echo "Cleaning up: Removing unsorted BAM files"
    echo "-------------------------"

    rm "$dir_in"/${base}_cds_unsorted.bam

    
    echo "-------------------------"
    echo "Bowtie DONE aligning "$base""
    echo "-------------------------"

done

```


## Ribomala usage

### Indexing

```py linenums="1"

ribomala index \
--fasta /path/GRCh38_ensembl_cds_plus18_annot.fa

```

The above code will generate `GRCh38_ensembl_cds_plus18_annot.csv` file. This file will be required by the `analysis` module.

### QC

```python linenums="1"

ribomala qc \
--samples /path/samples.txt \
--min 27 \
--max 34 \
--input /path/ \
--output /path/output/

```

where `samples.txt` contains:

```
> head sample.txt
sample_1.bam
sample_2.bam
sample_3.bam

```

### Analysis

```py linenums="1"

ribomala analysis \
--txcsv /path/GRCh38_ensembl_cds_plus18_annot.csv \
--samples /path/sample_sheet.txt \
--input /path/bams/ \
--output /path/output/ \
--codon CAT,CAC

```

where `sample_sheet.csv` contains:

```

file_name	condition	read_length	frame	offset	rep
sample_1.bam	X	30	2	16	A
sample_1.bam	X	31	2	16	A
sample_2.bam	X	29	0	15	B
sample_2.bam	X	30	0	15	B
sample_3.bam	Y	29	0	15	A
sample_3.bam	Y	30	0	15	A
sample_3.bam	Y	30	2	16	A
sample_3.bam	Y	31	2	16	A
sample_4.bam	Y	29	0	15	B
sample_4.bam	Y	30	0	15	B


```