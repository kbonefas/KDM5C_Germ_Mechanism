# The X-linked intellectual disability gene KDM5C is a sex-biased brake against germline programs in somatic lineages

* PubMed Link: N/A

## Abstract

Mutations in numerous chromatin-modifying enzymes cause neurodevelopmental disorders (NDDs) with unknown mechanisms. Loss of repressive chromatin regulators can lead to the aberrant transcription of tissue-specific genes outside of their intended context, however the mechanisms and consequences of their dysregulation are largely unknown. Here, we examine how the X-linked intellectual disability gene lysine demethylase 5c (KDM5C), an eraser of histone 3 lysine 4 di and tri-methylation (H3K4me2/3), contributes to tissue identity. We found male _Kdm5c_ knockout (-KO) mice, which recapitulate key human neurological phenotypes, aberrantly express many liver, muscle, ovary, and testis genes within the amygdala and hippocampus. Gonad-enriched genes misexpressed in the _Kdm5c_-KO brain are unique to germ cells, indicating an erosion of the soma-germline boundary. Germline genes are typically decommissioned in somatic lineages in the post-implantation epiblast, yet _Kdm5c_-KO epiblast-like cells (EpiLCs) aberrantly expressed key regulators of germline identity and meiosis, including _Dazl_ and _Stra8_. Characterizing germline gene misexpression in males and female mutants revealed germline gene repression is sexually dimorphic, with female EpiLCs requiring a higher dose of KDM5C to maintain germline gene suppression. Using a comprehensive list of mouse germline-enriched genes, we found KDM5C is selectively recruited to a subset of germline gene promoters that contain CpG islands (CGIs) to facilitate DNA CpG methylation (CpGme) during ESC to EpiLC differentiation. However, late-stage spermatogenesis genes devoid of promoter CGIs can become expressed in _Kdm5c_-KO cells via ectopic activation by RFX transcription factors. Together, these data demonstrate KDM5C's fundamental role in tissue identity and indicate that KDM5C acts as a brake against runaway activation of germline developmental programs in somatic lineages.

## Repository Structure

This repository contains code run on a linux server for large sequencing data analysis (code/serverscripts), as well as scripts used to generate figures locally. All locally run code was performed through a snakemake pipeline, and information on how files were generated can be found in the Snakefile. After individual figures were generated (results/figure_pieces), figures were combined together and adjusted for publication, such as altering colors or sizing (submission/compiled_figs). Supplementary tables included in this study can be found in results/ or data/processed, which were then combined together into multi-sheet excel files for publication (submission/tables).

## Sequencing Data Availability

* Previously published RNA-seq datasets
  * Male wild-type and _Kdm5c_-KO adult amygdala and hippocampus - GEO: GSE127722.
  * Male and female wild-type, _Kdm5c_-KO, and _Kdm5c_-HET EpiLCs - GEO: GSE96797.
* Previously published ChIP-seq datasets
  * KDM5C ChIP-seq in wild-type and _Kdm5c_-KO EpiLCs - GEO: GSE96797
  * KDM5C ChIP-seq in wild-type and _Kdm5c_-KO mouse primary neuron cultures (PNCs) from the cortex and hippocampus - GEO: GSE61036.
  * H3K4me2 ChIP-seq in male wild-type and _Kdm5c_-KO EpiLCs - GEO: GSE96797
  * H3K4me3 ChIP-seq ChIP-seq in male wild-type and _Kdm5c_-KO amygdala - GEO: GSE127817
* Whole genome bisulfite sequencing (WGBS)
  * Wild-type and _Kdm5c_-KO naive emrbyonic stem cells (nESCs) and 96 hour extended EpiLCs (exEpiLCs) - SRA bioProject PRJNA1165148.
