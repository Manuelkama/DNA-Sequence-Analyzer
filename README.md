
# DNA Sequence Analyzer

## 1. Project Objective

Analyze a human mitochondrial DNA sequence using Python and
Biopython. Investigate nucleotide composition, GC content,
and the COX1 gene. Translate COX1 into its protein sequence
and validate the result against the NCBI reference annotation.

## 2. Dataset

- Database: NCBI Nucleotide
- Accession: OR327029.1
- Organism: Homo sapiens
- Sample: Voucher BuKa3C
- Sequence type: Complete mitochondrial genome
- Genome length: 16,567 bp
- Gene analyzed: COX1

The dataset represents an ancient human specimen.

## 3. Tools and Libraries

- Python
- Biopython
- Pandas
- Matplotlib
- JupyterLab

## 4. Methodology

1. Imported the mitochondrial genome in FASTA format.
2. Counted A, T, C and G nucleotides.
3. Calculated genome GC content.
4. Visualized nucleotide composition.
5. Extracted COX1 using GenBank coordinates.
6. Calculated COX1 nucleotide composition and GC content.
7. Translated COX1 using vertebrate mitochondrial genetic code (table 2).
8. Analyzed codon frequencies.
9. Compared the translated protein with the NCBI reference.
10. Exported the extracted gene and analysis results.

## 5. Results

| Measurement | Result |
|---|---:|
| Genome length | 16,567 bp |
| Genome GC content | 44.40% |
| COX1 length | 1,542 bp |
| COX1 GC content | 46.17% |
| COX1 protein length | 513 amino acids |
| Total codons | 514 |
| Unique codons | 55 |
| NCBI protein match | True |

## 6. Validation

The translated COX1 protein was compared with the
NCBI annotated protein WND63996.1.

Both sequences contained 513 amino acids and matched
exactly.

## 7. Scientific Reference

Cock PJA et al. (2009).
Biopython: freely available Python tools for computational
molecular biology and bioinformatics.
Bioinformatics.

https://doi.org/10.1093/bioinformatics/btp163

## 8. Data Source

NCBI Nucleotide:
https://www.ncbi.nlm.nih.gov/nuccore/OR327029.1
