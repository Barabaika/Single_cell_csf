<h1 align="center">
  scRNA-seq for drug repurposing in multiple sclerosis
  <br>
</h1>


[![](https://img.shields.io/github/languages/code-size/Barabaika/Single_cell_csf)](https://img.shields.io/github/languages/code-size/Barabaika/Single_cell_csf)
[![](https://img.shields.io/github/languages/top/Barabaika/Single_cell_csf)](https://img.shields.io/github/languages/top/Barabaika/Single_cell_csf)
[![](https://img.shields.io/github/issues/Barabaika/Single_cell_csf)](https://img.shields.io/github/issues/Barabaika/Single_cell_csf)

This repository contains primary source code for *"scRNA-seq for drug repurposing in multiple sclerosis"* manuscript. 

## Introduction

Multiple sclerosis (MS) is an autoimmune disease of a central nervous system still lacking a cure. Treatment typically focuses on slowing the progression and managing MS symptoms. Single cell transcriptomics study the immune system - the key player in the MS onset and development - in great detail increasing our understanding of MS mechanisms and stimulating the discovery of the targets for potential therapies. Still, de novo drug development takes decades, which can be reduced by drug repositioning. A promising approach is to select potential drugs based on activated or inhibited genes and pathways. In this study, we explored the public data of  peripheral blood mononuclear cells  and cerebrospinal fluid Cerebrospinal fluid (CSF) cells of patients with MS and idiopathic intracranial hypertension (IIH). We demonstrate that AIM2 inflammasome, complement activation pathways are activated in MS in different CSF and PBMC immune cells. Using genes from top activated pathways, we detected several promising small molecules, including *AG14361*, *FGIN-1-27*, *CA-074*, *ARP 101*, *Flunisolide* and *JAK3 Inhibitor VI for granulocytes*, to reverse MS immune cells transcriptomic signatures. Among these molecules we also detected an FDA approved MS drug Mitoxantrone supporting the reliability of our approach. 

## Requirements

Main dependencies are:
* *Seurat* 
* *dplyr*
* *harmony*
* *preprocessCore* >= 1.50.0
* *ggplot2*


## Data Availability

Data was taken from the study of Schafflick *et al.* (2020) ["Integrated single cell analysis of blood and cerebrospinal fluid leukocytes in multiple sclerosis"](https://doi.org/10.1038/s41467-019-14118-w)

## Scripts

Core analysis scripts can be found at `./core`

Overall pipeline of the analysis:

1. **Preprocessing, clustering and markers**

+ `CSF_analysis_umap_annot_de.R`
+ `PBMC_analysis_umap_annot_de.R`


2. **Cluster annotation with marker genes**

+ `python search_marker_genes.py`
+ `Search_for_cell_types.ipynb`

3. **Estimation of pathway activation by Oncobox**

```bash
source ~/venv/bin/activate

./oncobox_calc.sh \
	/data/ashevtsov/MS_data/CSF/geom_mean \
	/data/ashevtsov/MS_data/CSF/oncobox_res

```

Making a table with genes from top 10 biopathways (min gene number = 1)

```bash
python3 genes_from_top_pathways.py > genes_from_top10_pathways.log
```

4. **Connectivity Mapping**

+ `cmap_run.py`
+ `python cmap_results_merge.py`
+ `MS_state_plots.ipynb`

## 🆘 Help

Please feel free to contact **Andrey Shevtsov** (Andredeeandredee@gmail.com) and **Mikhail Raevskiy** (raevskii.mm@phystech.edu) if you have any questions about the software.
