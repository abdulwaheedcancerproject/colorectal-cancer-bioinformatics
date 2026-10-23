Reproducible Cross-Dataset Analysis of Colorectal Cancer Gene Expression



Project Overview



This project aims to investigate gene-expression differences between colorectal cancer tissue and adjacent normal tissue using RNA-sequencing count data. The goal is to build a reproducible bioinformatics workflow in Python, with appropriate statistical methods and independent dataset validation where supported by the study designs.



This is an ongoing independent portfolio project.



Research Questions



\- Which genes show evidence of differential expression between colorectal cancer and normal tissue?

\- Are the observed expression patterns consistent across independent datasets?

\- Which candidate genes and biological pathways warrant further investigation?



Datasets



The project currently includes two candidate datasets from the NCBI Gene Expression Omnibus (GEO):



\- GSE251845

\- GSE196006



Their sample metadata, study designs, and relationship to one another must be verified before making claims about patient pairing or independent validation.



Raw count matrices are kept locally and are not committed to this repository.



Current Progress



\- Downloaded the candidate RNA-seq count matrices.

\- Inspected the GSE251845 count matrix and sample names.

\- Identified and corrected a sample-name case-sensitivity issue during sample matching.

\- Removed five HTSeq summary rows from the working count matrix during quality-control checks.

\- Calculated per-sample library sizes and compared cancer and normal groups.

\- Initialized Git version control and published the project repository.



Planned Workflow



1\. Verify dataset metadata, sample labels, and patient pairing.

2\. Assess data quality and sample-level variation.

3\. Apply an appropriate normalization and differential-expression method.

4\. Correct for multiple hypothesis testing.

5\. Compare findings across datasets if their independence and compatibility are established.

6\. Investigate candidate genes and biological pathways.

7\. Document methods, limitations, and reproducible results.



Repository Structure



text

colorectal-cancer-bioinformatics/

├── .gitignore

├── check\_files.py

├── check\_library\_sizes.py

├── check\_pairs.py

├── check\_pairs\_v2.py

├── check\_samples.py

├── clean\_counts.py

├── compare\_library\_sizes.py

├── inspect\_counts.py

└── inspect\_matrix.py





Tools



\- Python

\- pandas

\- Git and GitHub

\- NCBI GEO



Reproducibility



The analysis scripts and project documentation are tracked with Git. Raw datasets and the local Python virtual environment are excluded from version control. Dataset accessions and processing steps will be documented as the workflow develops.



Limitations



This project is in progress. Differential-expression findings, biological interpretations, and independent validation conclusions have not yet been established.



Author



Independent bioinformatics portfolio project.



