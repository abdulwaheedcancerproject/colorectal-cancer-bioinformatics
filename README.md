Colorectal Cancer RNA-seq Differential Expression Analysis

Project Overview

This project investigates gene-expression differences between colorectal cancer tissue and adjacent normal tissue using RNA-sequencing data.

The aim is to develop a reproducible bioinformatics workflow for differential gene-expression analysis, quality control, and biological interpretation. A second dataset is being considered for cross-dataset comparison, subject to verification of study independence and compatibility.

This is an ongoing independent bioinformatics portfolio project.

Research Questions

- Which genes show evidence of differential expression between colorectal cancer and normal tissue?
- Are the observed expression patterns consistent across independent datasets?
- Which candidate genes and biological pathways warrant further investigation?

Datasets

The project currently includes two candidate datasets from the NCBI Gene Expression Omnibus (GEO):

- GSE251845 — primary dataset currently undergoing analysis.
- GSE196006 — candidate dataset for subsequent comparison.

The GSE251845 dataset contains 22 colorectal tumour samples and 22 matched adjacent normal samples, representing 22 patient pairs, according to the current sample metadata.

The relationship between the two studies must be verified before describing GSE196006 as an independent validation cohort.

Raw count matrices are stored locally and excluded from this repository.

Analysis Workflow

1. Data preparation

- Inspected the raw count matrix.
- Removed five HTSeq summary rows.
- Checked sample identifiers and corrected a case-sensitivity mismatch.
- Matched count-matrix samples with sample metadata.
- Checked patient pairing and count-data validity.

2. Differential expression

- Used DESeq2 for normalization and differential-expression analysis.
- Applied a paired experimental design: `~ patient_id + condition`.
- Compared tumour samples against matched normal samples.
- Applied multiple-testing correction using adjusted p-values.

3. Quality control

The current quality-control workflow includes:

- Library-size assessment.
- DESeq2 size-factor inspection.
- Normalized expression distribution plots.
- Principal component analysis (PCA).
- Paired PCA visualization.
- Sample-to-sample Pearson correlation.
- Investigation of potentially unusual samples.

All 44 samples have been retained so far. The samples investigated during QC did not provide sufficient evidence to justify exclusion.

Preliminary Results

For GSE251845:

- 44 samples: 22 tumour and 22 matched normal.
- 46,569 genes retained after removing genes with zero total counts.
- 15,158 genes with adjusted p-value below 0.05.
- 16,690 genes with adjusted p-value below 0.10.

These are preliminary computational results. Their biological interpretation and robustness require further investigation.

Tools and Technologies

- Python
- pandas
- R
- DESeq2
- ggplot2
- pheatmap
- Git and GitHub
- NCBI GEO

Repository Structure

text
colorectal-cancer-bioinformatics/
├── .gitignore
├── README.md
├── deseq2_analysis.R
├── results/
│   ├── GSE251845_deseq2_results.csv
│   ├── GSE251845_PCA_condition.png
│   ├── GSE251845_expression_boxplot.png
│   ├── GSE251845_library_sizes.csv
│   ├── GSE251845_paired_PCA.png
│   ├── GSE251845_paired_PCA_distances.csv
│   └── GSE251845_sample_correlation.png
└── Python analysis scripts
```

Planned Next Steps

1. Investigate the top differentially expressed genes.
2. Generate and interpret a volcano plot.
3. Annotate Ensembl gene identifiers.
4. Perform functional and pathway enrichment analysis.
5. Analyse GSE196006 separately.
6. Verify whether cross-dataset validation is scientifically appropriate.
7. Prepare a final report documenting the findings and limitations.

Reproducibility

Analysis scripts and selected results are tracked with Git. Raw datasets and the local Python virtual environment are excluded from version control.

The dataset accessions, processing steps, statistical methods, and analysis decisions will be documented throughout the project.

Limitations

The results are preliminary. Differential-expression findings require biological interpretation, and independent validation has not yet been established.

Author

Independent bioinformatics portfolio project.


