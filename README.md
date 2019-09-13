### 

This repository contains the code used in the analyses for the paper "Expert Curation of the Human and Mouse Olfactory Receptor Gene Repertoires Identifies Conserved Coding Regions Split Across Two Exons", by If Barnes, Ximena Ibarra-Soria and collaborators.

#### Transcriptome of human olfactory mucosa

The markdown document `human_ORexpression` contains the code used to analyse the bulk RNA-seq data from human biopsies of olfactory mucosa. 

- It produces normalised counts for all OR genes.
- Compares transcript length as a function of amount of sequencing data available.
- Analyses the split OR genes presented in the paper (Figure 5).

Data is available at the EGA under study [EGAS00001001486](https://www.ebi.ac.uk/ega/studies/EGAS00001001486).

#### Single-cell RNA-seq of mouse OSNs

The markdown document `singleCell_splitORexpression` contains the code used to analyse the single-cell RNA-seq data from 34 manually picked OSNs, obtained from heterozygous OMP-GFP mice, at postnatal day 3.

- It performs QC and normalisation.
- Analyses the expression of OR genes in each cell, and investigates mismapping events.
- Compares the two cells expressing split ORs to the rest that express intronless ORs.

Data is available at ArrayExpress under accession number [E-MTAB-8285](https://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-8285).