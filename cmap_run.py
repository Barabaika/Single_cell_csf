#!/usr/bin/python3
# -*- coding: utf-8 -*-

import requests
import json
import pandas as pd
import numpy as np
import os
from time import sleep

REVERSED_MOD = True
sample = "CSF"

PATH = f"/data/ashevtsov/MS_data/{sample}/oncobox_res/merged/{sample}_genes_from_top10_all_cells_bipathways_min1genes.csv"
OUTPATH = f"/data/ashevtsov/MS_data/{sample}/cmap_res/"

# CSF
if sample == "CSF":
  SAMPLES = ["T cell", "ab T cell", "a CD8+ T cell",
                     "na CD8+ T cell", "monocyte",
                     "mDC2", "mDC", "NK cell", "gd T cell", "naive B cell",
                     "granulocyte", "Plasmablast", "GrB+ pDC", "mDC1"]

elif sample == "PBMC":           
  SAMPLES = ["na CD8+ T cell", "a CD8+ T cell",
                     "naive B cell", "monocyte",
                     "gd T cell", "ab T cell", "mDC",
                     "granulocyte", "NK cell", "megakaryocyte", "GrB+ pDC", "mDC2"]


SAMPLES = [ "_" + "".join(i.split(" ")) + f"_{sample}" for i in SAMPLES ]


url = "https://maayanlab.cloud/L1000CDS2/query"
def upperGenes(genes):
    # The app uses uppercase gene symbols. So it is crucial to perform upperGenes() step.
    return [gene.upper() for gene in genes]

if REVERSED_MOD:
  agregative = False
else:
  agregative = True
  
path_genes_df = pd.read_csv(PATH)

for sample in SAMPLES:
  print(sample)
  
  pos_genes = path_genes_df[sample + "_pos"].unique().tolist()
  neg_genes = path_genes_df[sample + "_neg"].unique().tolist()
  print(f"initial unqiue genes number: pos- {len(pos_genes)}, neg- {len(neg_genes)}")
  
  intersected = list(set(pos_genes) & set(neg_genes))
  intersected.append(np.nan)
  pos_genes = [g for g in pos_genes if g not in intersected]
  neg_genes = [g for g in neg_genes if g not in intersected]
  print(f"Nonintersected unqiue genes number: pos- {len(pos_genes)}, neg- {len(neg_genes)}")
  
  # gene-set search example
  data = {"upGenes": pos_genes,
          "dnGenes": neg_genes}
  data["upGenes"] = upperGenes(data["upGenes"])
  data["dnGenes"] = upperGenes(data["dnGenes"])
  config = {"aggravate":agregative,"searchMethod":"geneSet","share":True,"combination":False,"db-version":"latest"}
  # metadata = [{"key":"Tag","value":"gene-set python example"}]
  payload = {"data":data,"config":config} # ,"meta":metadata
  headers = {"content-type":"application/json"}
  r = requests.post(url,data=json.dumps(payload),headers=headers)
  resGeneSet = r.json()
  res = pd.DataFrame(resGeneSet["topMeta"])
  filename = f"{sample}_cmap_res.csv"
  res.to_csv(os.path.join(OUTPATH, filename))
  
  sleep(10)
  
