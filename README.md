# STAMMP

STAMMP (Shotgun Taxonomic and Metabolic Metagenomic Pipeline) is the Urtecho Lab's pipeline for processing paired-end metagenomic reads through human depletion, taxonomic classification, AGORA genome abundance estimation, and MICOM community metabolic modeling.

---

## How the pipeline works

The pipeline runs in two sequential stages. The first stage uses Snakemake to process raw reads through a series of filtering and alignment steps that ultimately produce a per-taxon read count for each sample. The second stage uses those counts to build and simulate community metabolic models with MICOM.

Because large sample sets can be memory-intensive to run all at once, the pipeline is designed around batching: raw FASTQ files are split into small groups using `batch_samples.sh`, and the Snakemake workflow is run on one batch at a time by the orchestration script `run_batched_snakemakes.sh`. After all batches complete, their per-sample count files are merged into a single matrix before MICOM runs.

### Stage 1: AGORA read processing (Snakemake)

Each sample passes through five Snakemake rules in order:

**`qc_reads`** — fastp trims adapters and filters low-quality reads. The cleaned paired-end reads are written to `outputs/agora_outputs/cleaned/`.

**`bowtie_filter`** — Bowtie2 aligns reads to the CHM13+Y human reference genome. Only read pairs that do not map (i.e., non-human reads) are kept. Output goes to `outputs/agora_outputs/bowtie_filtered/`.

**`kraken_classify`** — Kraken2 assigns taxonomic labels to the human-depleted reads using the 8 GB standard database. This produces a `.kraken` classification file and a `.report` summary per sample.

**`kraken_extract`** — A custom script (`agora_filter_scripts/extract_kraken_reads.py`) reads the Kraken2 output and discards all reads classified as eukaryotic (taxid 2759, including children). The surviving reads are compressed with pigz and written to `outputs/agora_outputs/kraken_filtered/`. Intermediate uncompressed files are cleaned up automatically.

**`minimap_agora_count`** — minimap2 aligns the eukaryote-filtered reads to a pre-built AGORA genome index using short-read mode (`-x sr`). SAMtools filters for properly paired primary alignments, and a Perl one-liner counts hits per AGORA taxon ID. The result is a tab-separated `<sample>.counts` file in `outputs/agora_outputs/AGORA_mapped/`. Cleaned and Bowtie-filtered intermediate files are deleted at this point to save disk space.

After all batches finish, `agora_filter_scripts/count_files_to_table.py` merges all `.counts` files into a single CSV matrix (`agora_mapped_counts.csv`) where rows are taxon IDs and columns are sample IDs.

### Stage 2: MICOM metabolic modeling

`pipeline_scripts/MICOM_pipeline.py` takes the merged count matrix and runs three steps:

**Normalization** — read counts are divided by ungapped genome size (from `micom_references/genome_size_ungapped.csv`) to produce relative abundances. The normalized table is saved to `micom_outputs/processed/`.

**Model building** — MICOM's `build` function constructs a community metabolic model for each sample using the AGORA species database (`micom_references/agora_species_103.qza`). Taxa below a relative abundance cutoff of 0.0001 are excluded. Models are saved as pickle files under `micom_outputs/models/`.

**Growth simulation** — MICOM's `grow` function simulates community growth on a Western diet medium (`micom_references/western_diet_gut_agora.qza`) at a cooperative tradeoff value of 0.7. Results are saved as `micom_outputs/models/res_western.zip`.

---

## Dependencies

The pipeline uses two separate conda environments that cannot be merged due to conflicting dependencies.

**`ut-agora_pipeline_clone`** is required for the Snakemake stage and must contain: Snakemake, fastp, Bowtie2, Kraken2, minimap2, samtools, pigz, and the Python dependencies for `extract_kraken_reads.py`.

**`micomPy310`** is required for the MICOM stage and must contain: Python 3.10, micom, cobra, and pandas.

The pipeline also depends on three external databases that are stored on the Zarrinpar cluster. These paths are hard-coded in `pipeline_scripts/capscan_batch_snakemake_v2` and can be overridden via `--config` arguments if running on a different system:

- `kraken_db` — Kraken2 standard 8 GB database (`k2_standard_08_GB_20250714`)
- `bowtie_idx` — Bowtie2 CHM13+Y human genome index
- `agora_db` — minimap2 AGORA genome index (`agora_genomes.idx`)

The MICOM reference files (`agora_species_103.qza`, `western_diet_gut_agora.qza`, `genome_size_ungapped.csv`, `reference_taxonomy.csv`) ship inside the `micom_references/` directory of this repository.

---

## Repository structure

```
STAMMP/
├── pipeline_scripts/
│   ├── run_batched_snakemakes.sh       Orchestration script — runs all batches, merges, launches MICOM
│   ├── capscan_batch_snakemake_v2      Snakemake workflow (qc → bowtie → kraken → minimap → counts)
│   ├── batch_samples.sh                Splits a flat FASTQ directory into batch_01/, batch_02/, ...
│   └── MICOM_pipeline.py               Normalizes counts, builds models, runs growth simulations
├── agora_filter_scripts/
│   ├── extract_kraken_reads.py         Filters reads by Kraken2 taxon, excluding eukaryotes
│   └── count_files_to_table.py         Merges per-sample .counts files into a single CSV matrix
├── micom_references/
│   ├── agora_species_103.qza           AGORA species model database (MICOM input)
│   ├── western_diet_gut_agora.qza      Western diet growth medium (MICOM input)
│   ├── genome_size_ungapped.csv        Genome sizes for read count normalization
│   └── reference_taxonomy.csv          Taxonomy lookup for AGORA taxon IDs
├── utils/
│   └── analysis.py                     Shared analysis utilities (imported by MICOM_pipeline.py)
└── batched_pipeline_guide.ipynb        Full walkthrough with commands for each pipeline stage
```

---

## How to run

> **Note on `--skip-micom`:** `run_batched_snakemakes.sh` calls `MICOM_pipeline.py` with whatever `python` is on `PATH` at the time — it does not switch conda environments partway through. Since the AGORA stage and the MICOM stage require two different, non-mergeable environments (see [Dependencies](#dependencies)), running the script *without* `--skip-micom` will fail at the MICOM step unless your active environment happens to have both dependency sets installed. The worked example below runs the two stages separately, each in its correct environment, which is the combination that actually works.

### Worked example

This example assumes you have already cloned STAMMP and have paired-end FASTQ files sitting in a directory — substitute your own paths for `~/data/example_project`.

**Step 0 — Enter the repository**

```bash
cd ~/STAMMP
```

**Step 1 — Batch the raw FASTQ files**

Your raw reads live in a flat directory, e.g. `~/data/example_project/raw_fastqs/`, containing pairs like `sample01_1.fastq.gz` / `sample01_2.fastq.gz`. `batch_samples.sh` moves (not copies) them into numbered `batch_XX/` subdirectories in place:

```bash
bash pipeline_scripts/batch_samples.sh ~/data/example_project/raw_fastqs '*_1.fastq.gz' 5
```

This creates `~/data/example_project/raw_fastqs/batch_01/`, `batch_02/`, ... each with up to 5 sample pairs. Omit the trailing `5` to use the default batch size of 6.

**Step 2 — Run the AGORA (Snakemake) stage**

```bash
conda activate /home/maholton/.conda/envs/maholton/envs/ut-agora_pipeline_clone

bash pipeline_scripts/run_batched_snakemakes.sh \
    --input  ~/data/example_project/raw_fastqs \
    --output ~/data/example_project \
    --cores 32 \
    --skip-micom
```

This creates the output directory structure, auto-generates a sample list, runs each batch through Snakemake sequentially, and merges the per-sample `.counts` files into `~/data/example_project/outputs/agora_outputs/AGORA_mapped/agora_mapped_counts.csv`.

**Step 3 — Run the MICOM stage**

Switch environments, then call `MICOM_pipeline.py` directly:

```bash
conda activate /home/maholton/.conda/envs/maholton/envs/micomPy310

/home/maholton/.conda/envs/maholton/envs/micomPy310/bin/python \
    pipeline_scripts/MICOM_pipeline.py \
    --agora     ~/data/example_project/outputs/agora_outputs/AGORA_mapped/agora_mapped_counts.csv \
    --output    ~/data/example_project/outputs \
    --reference micom_references
```

Growth simulation results land at `~/data/example_project/outputs/micom_outputs/models/res_western.zip`.

### Arguments (`run_batched_snakemakes.sh`)

| Argument | Required | Default | Description |
|---|---|---|---|
| `--input <dir>` | Yes | | Directory containing `batch_01/`, `batch_02/`, ... |
| `--output <dir>` | Yes | | Project root; all outputs written under `<dir>/outputs/` |
| `--cores <N>` | No | 32 | Total thread budget passed to Snakemake |
| `--skip-micom` | No | off | Stop after the AGORA merge step; do not run MICOM |

### Arguments (`MICOM_pipeline.py`)

| Argument | Required | Description |
|---|---|---|
| `-a, --agora <file>` | Yes | Merged AGORA count matrix CSV (output of Step 2) |
| `-o, --output <dir>` | Yes | Project `outputs/` directory; `micom_outputs/` is created underneath |
| `-r, --reference <dir>` | Yes | Path to `micom_references/` (ships with this repo) |

### Running a single batch manually

To debug or rerun one batch directly with `snakemake`, without going through the orchestration script:

```bash
conda activate /home/maholton/.conda/envs/maholton/envs/ut-agora_pipeline_clone

snakemake -s pipeline_scripts/capscan_batch_snakemake_v2 \
    --cores 32 \
    --rerun-incomplete \
    --config \
        raw_dir=~/data/example_project/raw_fastqs/batch_01 \
        output_dir=~/data/example_project \
        sample_list=~/data/example_project/outputs/agora_outputs/sample_list.txt \
        threads=8
```

`sample_list.txt` must already exist — it's normally auto-generated by `run_batched_snakemakes.sh` (one sample ID per line, no extension). `threads=8` with `--cores 32` lets Snakemake schedule up to 4 samples in parallel; see [Parallelism](#parallelism) below.

---

## Output structure

```
<project_dir>/outputs/
├── agora_outputs/
│   ├── sample_list.txt               Master list of all sample IDs (one per line)
│   ├── cleaned/                      fastp-trimmed reads (*_1/2.fastq.gz)
│   ├── bowtie_filtered/              Human-depleted reads (*_1/2.fastq)
│   ├── kraken_filtered/              Eukaryote-filtered reads (*_1/2.fastq.gz) + Kraken2 output
│   └── AGORA_mapped/
│       ├── <sample>.counts           Per-sample tab-separated taxon hit counts
│       └── agora_mapped_counts.csv   Merged count matrix (taxon x sample)
├── micom_outputs/
│   ├── processed/
│   │   ├── abundance_normalized_table.csv
│   │   └── abundance_normalized.csv
│   └── models/
│       ├── <sample>.pickle           Per-sample MICOM community models
│       └── res_western.zip           GrowthResults object (Western diet simulation)
└── humann_outputs/                   Reserved for HUMAnN (run separately if needed)
```

Intermediate files (cleaned reads, Bowtie-filtered reads) are removed by Snakemake after each sample's `.counts` file is written to conserve disk space.

---

## Parallelism

Within each Snakemake batch, jobs run in parallel subject to the `--cores` budget. The `threads` config parameter (default: 8) sets how many threads each rule is allocated. Snakemake schedules as many concurrent jobs as fit within the budget: `concurrent_samples = floor(--cores / threads)`. With `--cores 32` and `threads=8`, four samples process in parallel. The `kraken_extract` rule is capped at 4 threads regardless of the `threads` setting, since the Python extraction step is single-threaded and pigz does not benefit meaningfully from more than 4 threads at typical metagenomic file sizes.

Batches themselves run sequentially — one batch completes before the next starts. This is intentional: intermediate files from completed batches are cleaned up before the next batch begins, keeping total disk usage bounded.
