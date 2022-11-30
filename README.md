# Single_cell_csf

Data was taken from the state Integrated single cell analysis of blood and cerebrospinal fluid leukocytes in multiple sclerosis (https://doi.org/10.1038/s41467-019-14118-w)

Core scripts for analysing data are:
  CSF_analysis_umap_annot_de.R
  PBMC_analysis_umap_annot_de.R
  
Overall pipeline of the research:

1) ***Preprocessing, clustering, markers, umap plots***

CSF_analysis_umap_annot_de.R
PBMC_analysis_umap_annot_de.R


1.5)***(Annotating clusters ) Make a table with marker genes***

python search_marker_genes.py
Search_for_cell_types.ipynb


2)***Making geom mean of every cell type***

/data/ashevtsov/Single_cell_csf/for_inregnet/geom_mean_by_cluster_name.R

3)***Finding biopathways by oncobox***

source ~/venv/bin/activate

./oncobox_calc.sh \
	/data/ashevtsov/MS_data/CSF/geom_mean \
	/data/ashevtsov/MS_data/CSF/oncobox_res

7)*** Making a table with genes from top 10 biopathways (min gene number = 1)***

python genes_from_top_pathways.py > /data/ashevtsov/MS_data/CSF/genes_from_top10_pathways.log

8)***running connectivity mapping***

cmap_run.py

9)***merging tables to final molecules table***

python cmap_results_merge.py

10)***making cmap plots***

MS_state_plots.ipunb