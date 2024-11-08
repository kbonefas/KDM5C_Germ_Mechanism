# KDM5C is a sex-biased brake against germline gene expression programs in somatic lineages

* PubMed Link: N/A
* BioRxiv Link:

## Abstract

The division of labor among cellular lineages is a pivotal step in the evolution of multicellularity. In mammals, the soma-germline boundary is formed during early embryogenesis, when genes that drive germline identity are repressed in somatic lineages through DNA and histone modifications at promoter CpG islands (CGIs). Somatic misexpression of germline genes is a signature of cancer and observed in select neurodevelopmental disorders. However, it is currently unclear if all germline genes use the same repressive mechanisms and if factors like development and sex influence their dysregulation. Here, we examine how cellular context influences the formation of somatic tissue identity in mice lacking lysine demethylase 5c (KDM5C), an X chromosome eraser of histone 3 lysine 4 di and tri-methylation (H3K4me2/3). We found male _Kdm5c_ knockout (-KO) mice aberrantly express many tissue-specific genes within the brain, the majority of which are unique to the germline. By developing a comprehensive list of mouse germline-enriched genes, we observed _Kdm5c_-KO cells aberrantly express key drivers of germline fate during early embryogenesis but late-stage spermatogenesis genes within the mature brain. KDM5C binds CGIs within germline gene promoters to facilitate DNA CpG methylation as embryonic stem cells differentiate into epiblast-like cells (EpiLCs). However, the majority of late-stage spermatogenesis genes expressed within the _Kdm5c_-KO brain did not harbor promoter CGIs. These CGI-free germline genes were not bound by KDM5C and instead expressed through ectopic activation by RFX transcription factors. Furthermore, germline gene repression is sexually dimorphic, as female EpiLCs require a higher dose of KDM5C to maintain germline silencing. Altogether, these data revealed distinct regulatory classes of germline genes and sex-biased silencing mechanisms in somatic cells.

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
