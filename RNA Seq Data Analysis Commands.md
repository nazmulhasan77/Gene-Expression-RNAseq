# **RNA-Seq Data Analysis Pipeline Documentation**  
*From Raw Reads to Gene Counts*  

This document outlines the step-by-step process performed in the terminal for RNA-Seq data analysis, including quality control, read trimming, alignment, and quantification.  

---

## **1. Environment Setup**  
### **1.1 Create a Python Virtual Environment**  
```bash
python3 -m venv myenv           # Create a virtual environment
source myenv/bin/activate       # Activate it
```

### **1.2 Install Required Packages**  
```bash
pip install biopython pandas matplotlib seaborn cutadapt
```

---

## **2. Data Acquisition**  
### **2.1 Download Example Dataset**  
```bash
mkdir -p ~/rnaseq_practice && cd ~/rnaseq_practice
curl -O http://data.biostarhandbook.com/rnaseq/projects/griffith/griffith-data.tar.gz
tar -xvf griffith-data.tar.gz
```yes
**Output Structure**:  
```
rnaseq_practice/
├── reads/       # FASTQ files (e.g., HBR_1_R1.fq, UHR_1_R1.fq)
└── refs/        # Reference files (22.fa, 22.gtf)
```

---

## **3. Quality Control (FastQC)**  
### **3.1 Run FastQC on Raw Reads**  
```bash
cd ~/rnaseq_practice/reads
fastqc *.fq                      # Generates HTML/zip reports for each file
xdg-open HBR_1_R1_fastqc.html    # View report in browser
```
**Purpose**: Check per-base quality, adapter contamination, and sequence duplication.  

---

## **4. Read Trimming (Cutadapt)**  
### **4.1 Trim Low-Quality Bases and Short Reads**  
```bash
cutadapt -q 20 -m 30 -o UHR_1_R1_trimmed.fq UHR_1_R1.fq
```
**Parameters**:  
- `-q 20`: Trim bases with quality < 20.  
- `-m 30`: Discard reads shorter than 30 bp.  

**Output**:  
- `UHR_1_R1_trimmed.fq` (trimmed reads).  

### **4.2 Re-run FastQC on Trimmed Reads**  
```bash
fastqc UHR_1_R1_trimmed.fq
```

---

## **5. Alignment (HISAT2)**  
### **5.1 Build HISAT2 Index**  
```bash
cd ~/rnaseq_practice/refs
hisat2-build 22.fa 22_index      # Creates index files (e.g., 22_index.1.ht2)
```

### **5.2 Align Trimmed Reads to Reference**  
```bash
hisat2 -x 22_index -U ../reads/UHR_1_R1_trimmed.fq -S ../reads/UHR_1_aligned.sam
```
**Alignment Stats**:  
```
47.85% overall alignment rate
```

---

## **6. Convert SAM to BAM (Samtools)**  
```bash
samtools view -bS ../reads/UHR_1_aligned.sam | samtools sort -o ../reads/UHR_1_aligned.bam
samtools index ../reads/UHR_1_aligned.bam
```
**Output**:  
- `UHR_1_aligned.bam` (sorted BAM file).  
- `UHR_1_aligned.bam.bai` (BAM index).  

---

## **7. Quantification (featureCounts)**  
### **7.1 Install Subread**  
```bash
sudo apt install subread
```

### **7.2 Count Reads per Gene**  
```bash
featureCounts -a 22.gtf -o ../reads/UHR_1_counts.txt ../reads/UHR_1_aligned.bam
```
**Output**:  
- `UHR_1_counts.txt`: Gene-level counts.  
- `UHR_1_counts.txt.summary`: Summary statistics.  

### **7.3 View Counts**  
```bash
head ../reads/UHR_1_counts.txt
```
**Example Output**:  
```
Geneid              Chr     Start   End     Strand  Length  Counts
ENSG00000277248.1   chr22   10736171 10736283 -      113     0
ENSG00000274237.1   chr22   10936023 10936161 -      139     0
...
```

---

## **8. Next Steps**  
1. **Differential Expression**: Use `DESeq2` or `edgeR` in R.  
2. **Visualization**: Plot counts with `matplotlib/seaborn`.  
3. **Multi-Sample Analysis**: Repeat for all samples (e.g., `HBR_2_R1.fq`, `UHR_2_R1.fq`).  

---

## **Key Commands Summary**  
| Step               | Command                                                                 |
|--------------------|-------------------------------------------------------------------------|
| Quality Control    | `fastqc *.fq`                                                          |
| Trimming           | `cutadapt -q 20 -m 30 -o trimmed.fq input.fq`                          |
| Alignment          | `hisat2 -x index -U trimmed.fq -S aligned.sam`                         |
| SAM → BAM          | `samtools view -bS aligned.sam | samtools sort -o aligned.bam`       |
| Quantification     | `featureCounts -a annotation.gtf -o counts.txt aligned.bam`            |

---

## **Troubleshooting**  
- **Error**: `Command 'python' not found` → Use `python3`.  
- **Missing Files**: Verify paths (e.g., `UHR_1_R2.fq` was missing in paired-end trimming).  
- **Low Alignment Rate**: Check quality reports or adjust trimming parameters.  

This pipeline processes **single-end RNA-Seq data** from raw reads to gene counts. Modify for paired-end data by adding `-1`/`-2` flags in `hisat2` and `cutadapt`.  
