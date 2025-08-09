## **Differential Expression Analysis and Visualization** 
for RNA-Seq data, using **DESeq2** in R and visualizing with **ggplot2** (R) or **matplotlib/seaborn** (Python).

## Step 1: Load Count Data in R

Assuming you have a gene counts file like `HBR_1_counts.txt` with raw counts per gene.

```r
# Load necessary libraries
library(DESeq2)
library(ggplot2)
library(pheatmap)

# Read count data
counts <- read.table("HBR_1_counts.txt", header=TRUE, row.names=1)

# Check the data
head(counts)
```

---

## Step 2: Create Sample Metadata Table

Create a metadata table describing samples and conditions (e.g. treated vs control).

```r
sampleTable <- data.frame(
  sampleName = colnames(counts),
  condition = c("control", "treatment")  # change accordingly
)
```

---

## Step 3: Prepare DESeq2 Dataset

```r
dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = sampleTable,
                              design = ~ condition)
```

---

## Step 4: Pre-filter Low Counts (Optional but recommended)

```r
dds <- dds[rowSums(counts(dds)) > 10, ]
```

---

## Step 5: Run DESeq2 Analysis

```r
dds <- DESeq(dds)
res <- results(dds)

# View summary
summary(res)
```

---

## Step 6: Extract Significant Differentially Expressed Genes

```r
resSig <- subset(res, padj < 0.05 & abs(log2FoldChange) > 1)
head(resSig)
```

---

## Step 7: Plot MA Plot

```r
plotMA(res, main="DESeq2 MA-plot", ylim=c(-5,5))
```

---

## Step 8: Visualize Counts of a Specific Gene

```r
gene <- rownames(resSig)[1]
plotCounts(dds, gene=gene, intgroup="condition")
```

---

## Step 9: Heatmap of Top Variable Genes

```r
library(pheatmap)
vsd <- vst(dds, blind=FALSE)
select <- order(rowVars(assay(vsd)), decreasing=TRUE)[1:20]
pheatmap(assay(vsd)[select, ], cluster_rows=TRUE, show_rownames=TRUE,
         cluster_cols=TRUE, annotation_col=sampleTable)
```

---

# Optional: Visualization in Python with matplotlib/seaborn

If you want to export your DE results and plot in Python:

```python
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# Load DE results CSV
df = pd.read_csv("DESeq2_results.csv")

# Volcano plot
df['significant'] = (df['padj'] < 0.05) & (abs(df['log2FoldChange']) > 1)
plt.figure(figsize=(10,6))
sns.scatterplot(data=df, x='log2FoldChange', y=-np.log10(df['padj']),
                hue='significant', palette={True:'red', False:'grey'})
plt.xlabel('Log2 Fold Change')
plt.ylabel('-Log10 Adjusted P-value')
plt.title('Volcano Plot')
plt.show()
```

