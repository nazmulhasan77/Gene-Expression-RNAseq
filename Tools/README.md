# Installation Guide for RNA-Seq Analysis Tools

## 1. Prerequisites

Make sure you have the following installed:

* **Conda** (recommended for bioinformatics tools management)
  Install Miniconda (lightweight Conda):

  ```bash
  wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
  bash Miniconda3-latest-Linux-x86_64.sh
  # Follow prompts, then restart terminal or run:
  source ~/.bashrc
  ```

* Basic tools: `wget`, `curl`, `git`

  ```bash
  sudo apt update && sudo apt install wget curl git -y
  ```

---

## 2. FastQC — Quality Check

**Installation:**

Via Conda (recommended):

```bash
conda install -c bioconda fastqc
```

Or download standalone from:
[https://www.bioinformatics.babraham.ac.uk/projects/fastqc/](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)

---

## 3. MultiQC — Aggregates QC reports

**Installation:**

Via Conda:

```bash
conda install -c bioconda multiqc
```

Or via pip:

```bash
pip install multiqc
```

---

## 4. Trimmomatic — Adapter trimming and quality filtering

**Installation:**

Via Conda:

```bash
conda install -c bioconda trimmomatic
```

Or download jar from:
[http://www.usadellab.org/cms/?page=trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic)

---

## 5. Cutadapt — Adapter trimming (alternative)

**Installation:**

Via Conda:

```bash
conda install -c bioconda cutadapt
```

Or via pip:

```bash
pip install --user cutadapt
```

---

## 6. STAR — RNA-Seq read aligner (splice-aware)

**Installation:**

Via Conda:

```bash
conda install -c bioconda star
```

Or download from:
[https://github.com/alexdobin/STAR/releases](https://github.com/alexdobin/STAR/releases)

---

## 7. HISAT2 — Fast splice-aware aligner

**Installation:**

Via Conda:

```bash
conda install -c bioconda hisat2
```

Or download from:
[https://daehwankimlab.github.io/hisat2/](https://daehwankimlab.github.io/hisat2/)

---

## 8. FeatureCounts — Read counting from alignments

**Installation:**

Via Conda (part of Subread package):

```bash
conda install -c bioconda subread
```

Then use `featureCounts` command.

---

## 9. HTSeq-count — Counting reads per gene

**Installation:**

Via Conda:

```bash
conda install -c bioconda htseq
```

Or via pip:

```bash
pip install htseq
```

---

## 10. Kallisto — Pseudo-alignment and quantification

**Installation:**

Via Conda:

```bash
conda install -c bioconda kallisto
```

Or download precompiled binaries:
[https://pachterlab.github.io/kallisto/download](https://pachterlab.github.io/kallisto/download)

---

## 11. Salmon — Fast transcript quantification

**Installation:**

Via Conda:

```bash
conda install -c bioconda salmon
```

Or download binaries from:
[https://combine-lab.github.io/salmon/](https://combine-lab.github.io/salmon/)

---

## 12. R and Bioconductor — Statistical analysis and visualization

**Installation:**

### Install R:

Ubuntu:

```bash
sudo apt install r-base
```

Windows/Mac: Download from [https://cran.r-project.org/](https://cran.r-project.org/)

### Install R packages:

Launch R (type `R` in terminal) and run:

```r
install.packages("BiocManager")
BiocManager::install(c("DESeq2", "edgeR", "tximport", "apeglm", "EnhancedVolcano"))
```

---

## 13. Optional Utilities

* **Samtools** (for BAM file processing)

  ```bash
  conda install -c bioconda samtools
  ```
* **Picard tools** (BAM file QC)

  ```bash
  conda install -c bioconda picard
  ```

---

## 14. Example: Creating a Conda Environment with All Tools

You can create a dedicated environment for RNA-Seq analysis:

```bash
conda create -n rnaseq-env -c bioconda -c conda-forge fastqc multiqc trimmomatic cutadapt star hisat2 subread htseq kallisto salmon samtools picard r-base
conda activate rnaseq-env
```

Then install R packages inside R as shown above.

---

# Notes

* Always check the versions with commands like `fastqc --version`, `star --version`, `R --version`.
* Use `conda update --all` inside your environment to keep tools updated.
* For Mac, most tools can also be installed via **Homebrew**.
* For Windows, using Conda (Anaconda/Miniconda) is usually easiest. Some tools may require WSL (Windows Subsystem for Linux).

