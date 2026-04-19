# Bulk RNA-seq Pipeline

A Snakemake pipeline for bulk RNA-seq analysis: FastQC → fastp → FastQC → MultiQC → Salmon → MultiQC → DESeq2 (PCA + DEG).

## Features

- Automatic paired-end FASTQ detection (`_R1/_R2`, `_1/_2`, `_f/_r`, etc.)
- Decoy-aware Salmon index construction
- fastp adapter trimming and quality filtering
- MultiQC reports at each stage
- DESeq2 differential expression with PCA, Volcano plots, and Heatmaps
- Conda environment management per rule

---

## Directory Structure

```
bulk-rnaseq-pipeline/
├── Snakefile               # Main pipeline
├── config.yaml             # All parameters
├── samples.tsv             # Sample metadata
├── envs/
│   ├── qc.yaml             # FastQC, fastp, MultiQC
│   ├── salmon.yaml         # Salmon
│   └── deseq2.yaml         # R, DESeq2, tximport
├── scripts/
│   └── deseq2_analysis.R   # DESeq2 script
├── refs/                   # Reference files (user-provided)
│   ├── genome.fa
│   ├── transcriptome.fa
│   └── annotation.gtf
└── data/
    └── raw/                # Input FASTQ files
```

---

## Quick Start

### 1. Install Snakemake

```bash
conda create -n snakemake -c conda-forge -c bioconda snakemake=8.5.3 python=3.11
conda activate snakemake
```

### 2. Clone the repository

```bash
git clone https://github.com/your-username/bulk-rnaseq-pipeline.git
cd bulk-rnaseq-pipeline
```

### 3. Prepare reference files

Download reference files for your species. Example for human (GRCh38):

```bash
mkdir -p refs

# Genome FASTA
wget -O refs/genome.fa.gz \
  https://ftp.ensembl.org/pub/release-111/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
gunzip refs/genome.fa.gz

# Transcriptome FASTA
wget -O refs/transcriptome.fa.gz \
  https://ftp.ensembl.org/pub/release-111/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz
gunzip refs/transcriptome.fa.gz

# GTF annotation
wget -O refs/annotation.gtf.gz \
  https://ftp.ensembl.org/pub/release-111/gtf/homo_sapiens/Homo_sapiens.GRCh38.111.gtf.gz
gunzip refs/annotation.gtf.gz
```

For mouse (GRCm39), replace `homo_sapiens` with `mus_musculus` and update `org.Mm.eg.db` in `envs/deseq2.yaml`.

### 4. Place FASTQ files

```bash
mkdir -p data/raw
# Copy or symlink your FASTQ files
cp /path/to/your/data/*.fastq.gz data/raw/
```

Supported naming formats:

| R1 | R2 |
|---|---|
| `sample_R1.fastq.gz` | `sample_R2.fastq.gz` |
| `sample_1.fastq.gz` | `sample_2.fastq.gz` |
| `sample_f.fastq.gz` | `sample_r.fastq.gz` |
| `sample_F.fastq.gz` | `sample_R.fastq.gz` |

### 5. Edit `samples.tsv`

```tsv
sample    condition    replicate
ctrl_1    control      1
ctrl_2    control      2
ctrl_3    control      3
treat_1   treatment    1
treat_2   treatment    2
treat_3   treatment    3
```

- `sample`: must match the prefix of your FASTQ filenames
- `condition`: group label used in DESeq2 (e.g., `control`, `treatment`)
- `replicate`: biological replicate number

### 6. Edit `config.yaml`

Key settings to update:

```yaml
fastq_dir: "data/raw"

genome:
  genome_fasta: "refs/genome.fa"
  transcriptome_fasta: "refs/transcriptome.fa"
  gtf: "refs/annotation.gtf"
  salmon_index: "refs/salmon_index"

deseq2:
  contrasts:
    - ["treatment", "control"]   # [treatment_group, reference_group]
  padj_cutoff: 0.05
  log2fc_cutoff: 1.0
```

### 7. Run the pipeline

**Dry run (check workflow without executing):**

```bash
snakemake --use-conda -n
```

**Local execution:**

```bash
snakemake --use-conda --cores 16
```

**HPC with SLURM:**

```bash
snakemake --use-conda \
  --executor slurm \
  --default-resources slurm_account=your_account mem_mb=32000 \
  --jobs 50
```

**HPC with PBS/Torque:**

```bash
snakemake --use-conda \
  --cluster "qsub -l nodes=1:ppn={threads} -l mem={resources.mem_mb}mb" \
  --jobs 50
```

---

## Output Structure

```
results/
├── fastqc/
│   ├── raw/                    # Pre-trim FastQC reports
│   └── trimmed/                # Post-trim FastQC reports
├── fastp/                      # fastp HTML + JSON reports
├── trimmed/                    # Trimmed FASTQ files
├── multiqc/
│   ├── pre_trim/               # MultiQC: raw FastQC
│   ├── post_trim/              # MultiQC: fastp + trimmed FastQC
│   └── salmon/                 # MultiQC: Salmon mapping stats
├── salmon/
│   └── {sample}/
│       ├── quant.sf            # Transcript-level quantification
│       └── lib_format_counts.json
└── deseq2/
    ├── pca_plot.pdf            # PCA of all samples
    ├── {contrast}_DEG_results.csv    # Full DEG table
    ├── {contrast}_volcano.pdf        # Volcano plot
    └── {contrast}_heatmap.pdf        # Top 50 DEGs heatmap
```

### DEG results table columns

| Column | Description |
|---|---|
| `gene_id` | Ensembl gene ID |
| `baseMean` | Mean normalized count |
| `log2FoldChange` | log2 fold change (treatment / control) |
| `lfcSE` | Standard error of LFC |
| `stat` | Wald statistic |
| `pvalue` | Raw p-value |
| `padj` | BH-adjusted p-value |
| `diffexpressed` | `UP`, `DOWN`, or `NO` |

---

## Configuration Reference

| Parameter | Default | Description |
|---|---|---|
| `fastp.min_length` | 36 | Minimum read length after trimming |
| `fastp.qualified_quality_phred` | 20 | Minimum base quality |
| `salmon.lib_type` | `A` | Library type (`A` = auto-detect) |
| `salmon.extra` | `--gcBias --seqBias` | Extra Salmon flags |
| `deseq2.padj_cutoff` | 0.05 | Adjusted p-value threshold |
| `deseq2.log2fc_cutoff` | 1.0 | |log2FC| threshold |
| `deseq2.count_type` | `lengthScaledTPM` | tximport count scaling method |

---

## Multiple Comparisons

To run multiple pairwise comparisons, add entries to `contrasts` in `config.yaml`:

```yaml
deseq2:
  contrasts:
    - ["treatment_A", "control"]
    - ["treatment_B", "control"]
    - ["treatment_A", "treatment_B"]
```

Each contrast produces its own CSV, volcano plot, and heatmap.

---

## Requirements

- Snakemake >= 8.0
- conda / mamba
- ~32 GB RAM recommended for human genome indexing
- ~100 GB disk space for intermediate files (per 6-sample experiment)

---

## Citation

If you use this pipeline, please cite the underlying tools:

- **fastp**: Chen et al., *Bioinformatics* 2018
- **Salmon**: Patro et al., *Nature Methods* 2017
- **DESeq2**: Love et al., *Genome Biology* 2014
- **tximport**: Soneson et al., *F1000Research* 2015
- **MultiQC**: Ewels et al., *Bioinformatics* 2016
