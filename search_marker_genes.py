#!/usr/bin/python3
# -*- coding: utf-8 -*-

import os
import pandas as pd
import numpy as np

markers_list = ["CD3E", "TRAC", "LCK", "IL7R", "CD4", "CD8B", "CD8A", "CCL5", "CD8na", "CCR7" , "FOXP3", "CTLA4" , "TRDC", "GNLY", "NKG7", "FCGR3A", "CD16", "PRF1", "SELL", "CD62L", "XCL1", "CD79A", "IGHD", "CD37", "CD27", "IGHM", "IGHG", "CD38", "TNFRSF17", "CD269", "LYZ", "WDFY4", "XCR1", "BATF3", "FCER1A", "CD1C", "CLEC10A", "S100A8", "S100A9", "FCGR3A", "CD16", "CD14", "TCF4", "E2-2","TNFRSF21", "DR6", "GNG11", "GATA3", "CCR3", "CCR4", "CXCR3", "CD68", "Cd7", "GZMB"]

SAMPLE = "CSF"
PATH_TO_CLUSTER_MARKERS = f"/data/ashevtsov/MS_data/{SAMPLE}/{SAMPLE}_integrated_markers.csv"
OUTPATH = f"/data/ashevtsov/MS_data/{SAMPLE}/{SAMPLE}_needed_markers_in_clusters.csv"


clust_markers_df = pd.read_csv(PATH_TO_CLUSTER_MARKERS)
# filtr by P value
clust_markers_df = clust_markers_df[clust_markers_df.p_val < 0.001]

res_dict = {}
clusters_n_list = list(range(clust_markers_df.cluster.max()))
for marker in markers_list:
    res_dict[marker] = dict(zip(clusters_n_list, [np.nan for _ in range(clust_markers_df.cluster.max() + 1)]))
    clusters = clust_markers_df[clust_markers_df.gene == marker].cluster.tolist()
    if len(clusters) > 0:
        for cluster in clusters:
            res_dict[marker][int(cluster)] = clust_markers_df[
              (clust_markers_df.gene == marker) & (clust_markers_df.cluster == cluster)
              ].avg_log2FC.values[0]

res_df = pd.DataFrame(res_dict)
res_df = res_df.T
        
res_df.to_csv(OUTPATH)  
