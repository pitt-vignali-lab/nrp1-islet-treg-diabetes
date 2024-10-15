# Nrp1 loss on intra-islet Tregs in autoimmune diabetes

Analysis code for:

> **Treg insufficiency in autoimmune diabetes is driven by selective loss of
> Neuropilin-1 on intra-islet regulatory T cells**
> *The Journal of Immunology* (2024).
> PMID [39109924](https://pubmed.ncbi.nlm.nih.gov/39109924/)
> · PMC[11371503](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC11371503/)

This repository contains the single-cell analysis code behind **Figure 4**,
**Figure 5** and **Supplemental Table 1**. The remaining figures are flow
cytometry, diabetes incidence, histology and adoptive transfer experiments,
which are described in the paper and have no code here.

---

## What the analysis does

Tregs and Tconv cells were sorted from the non-draining lymph node and the
islets of six `RosaNrp1` and six Cre-control NOD mice at 6 weeks of age, and
profiled by 5′ single-cell RNA-seq with paired single-cell TCR-seq and Nrp1
CITE-seq.

- Tregs are separated from the Tconv mixture by **NMF**, then clustered and
  annotated.
- Pseudobulk profiles per genotype and tissue drive **SPIA** pathway analysis
  (Fig 4a) and **GSEA** cross-verification, showing increased Treg function and
  TCR stimulation when Nrp1 is restored.
- The Nrp1 CITE-seq channel splits Tregs into Nrp1<sup>+</sup> and
  Nrp1<sup>−</sup> populations, which are then compared within and across
  genotypes and tissues (Fig 5).
- Paired TCR libraries are analysed for repertoire diversity and clonality.

---

## Repository layout

```
1.Data/
  revision/   main dataset and derived objects
  initial/    earlier processing of the same experiment, used by Fig 5d-e
2.Script/     analysis scripts, numbered in run order
3.Output/     figures and tables produced by the scripts
```

Scripts expect the **repository root** as the R working directory:

```r
setwd("/path/to/nrp1-islet-treg-diabetes")
source("2.Script/01_sample_info_combine.R")
```

Each script carries a header giving its figure panel, inputs, outputs, the
upstream script it depends on, and its original filename in the lab archive.

> **Why two data namespaces.** The dataset was processed twice over the life of
> the project. Figures 4 and 5a–c and the TCR analyses come from the later
> processing (`revision`); Figures 5d and 5e were rendered from the earlier one
> (`initial`). Both rounds produce objects with the same filenames but different
> contents, so they are kept apart rather than merged.

---

## Figure → script map

| Figure | Analysis | Script |
| --- | --- | --- |
| — | Combine per-sample metadata with CellRanger matrices | [`01_sample_info_combine.R`](2.Script/01_sample_info_combine.R) |
| — | QC, filtering, normalisation, dimensionality reduction | [`02_qc_and_cleaning.R`](2.Script/02_qc_and_cleaning.R) |
| — | Identify Tregs by NMF | [`03_nmf_treg_identification.R`](2.Script/03_nmf_treg_identification.R) |
| — | NMF clustering including Foxp3<sup>+</sup> cells | [`04_nmf_cluster_with_foxp3_cells.R`](2.Script/04_nmf_cluster_with_foxp3_cells.R) |
| — | Build integrated Treg and Tconv objects | [`05_scvi_integration.R`](2.Script/05_scvi_integration.R) |
| — | Transfer cell-type annotations | [`06_annotation_prediction.R`](2.Script/06_annotation_prediction.R) |
| — | Assign final cluster identities | [`07_cluster_labeling.R`](2.Script/07_cluster_labeling.R) |
| — | Remove APC contamination from the Tconv object | [`08_tconv_remove_apcs.R`](2.Script/08_tconv_remove_apcs.R) |
| **4** | UMAPs, cluster statistics, pseudobulk aggregation | [`09_umap_statistics_pseudobulk.R`](2.Script/09_umap_statistics_pseudobulk.R) |
| **4a** | Unbiased SPIA on the KEGG database | [`10_fig4a_spia_kegg.R`](2.Script/10_fig4a_spia_kegg.R) |
| **4** | GSEA cross-verification of the SPIA pathways | [`11_fig4_treg_gsea.R`](2.Script/11_fig4_treg_gsea.R) |
| **5a** | Nrp1<sup>+</sup> / Nrp1<sup>−</sup> split from CITE-seq | [`12_fig5a_citeseq_nrp1_split.R`](2.Script/12_fig5a_citeseq_nrp1_split.R) |
| **5b–c** | DE between Nrp1<sup>+</sup> and Nrp1<sup>−</sup> Tregs; Venn and volcano | [`13_fig5bc_de_nrp1_pos_neg.R`](2.Script/13_fig5bc_de_nrp1_pos_neg.R) |
| **5d** | GSEA against TCR signalling gene sets | [`14_fig5d_gsea_tcr_signaling.R`](2.Script/14_fig5d_gsea_tcr_signaling.R) |
| **5e** | Heat map of Treg stability and function markers | [`15_fig5e_treg_stability_heatmap.R`](2.Script/15_fig5e_treg_stability_heatmap.R) |
| — | Treg TCR repertoire diversity | [`16_tcr_treg_diversity.R`](2.Script/16_tcr_treg_diversity.R) |
| — | Treg clonotype quantification | [`17_tcr_treg_clonotypes.R`](2.Script/17_tcr_treg_clonotypes.R) |
| — | Tconv repertoire diversity, clonotypes, clonal homeostasis | [`18_tcr_tconv_repertoire.R`](2.Script/18_tcr_tconv_repertoire.R) |
| **Supp Table 1** | Collect significant pseudobulk DEGs | [`19_supp_table1_pseudobulk_degs.R`](2.Script/19_supp_table1_pseudobulk_degs.R) |

Run `01` → `09` in order; `09` writes the pseudobulk tables that `10`, `11` and
`19` read.

---

## Data availability

**No sequencing data or intermediate objects are included in this repository.**

| Dataset | Source |
| --- | --- |
| 5′ scRNA-seq, scTCR-seq and Nrp1 CITE-seq generated for this study | GEO [**GSE262874**](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE262874) |

Download the CellRanger output into `1.Data/revision/` (gene expression under
`filtered_feature_bc_matrix/`, hashtag and Nrp1 counts under `umi_count/` and
`Nrp1umi_count/`, and TCR contigs under `TCR/`).

### Reference and curated inputs

| File | Status |
| --- | --- |
| `1.Data/initial/Nrp1 altered genes.csv` | **Included.** Curated list of Nrp1-altered genes with functional annotation, used for Fig 5d–e. |
| `1.Data/revision/RefData/` | Not included. Local KEGG data directory for `SPIA`, generated with `SPIA::makeSPIAdata()`, plus a mouse Entrez-to-symbol table. |
| `1.Data/initial/Nrp1 Fl vs WT Til Tregs abbys paper.xlsx`, `1.Data/initial/ITRhi_vs_ITRlo.xlsx` | Not included. Published supplementary tables from earlier studies, read by `14_fig5d_gsea_tcr_signaling.R` as comparison gene sets. |

---

## Software

| Tool | Version |
| --- | --- |
| bcl2fastq | 2.20.2 |
| CellRanger (`count` and `vdj`) | 6.0.1, mm10 reference |
| CITE-seq-Count | 1.4.3 |
| R | 4.1.3 |
| Seurat | 4.1.0 |
| NMF | 0.23.0 |
| presto (`wilcoxauc`) | 1.0.0 |
| SPIA | 2.44.0 |
| fgsea | 1.18.0 |
| scRepertoire | 1.3.5 |
| ComplexHeatmap | — |
| scvi-tools | via `reticulate` / `sceasy` |

SPIA pathways were called significant at Bonferroni-adjusted *pGFWER* < 0.05;
GSEA pathways at adjusted *P* < 0.05.

---

## Contact

**Jian Cui** — Senior Bioinformatics Scientist  
Vignali Lab, Department of Immunology  
School of Medicine, University of Pittsburgh  
<cuij@pitt.edu>

Questions about the analysis or this code are welcome.

## License

[MIT](LICENSE)
