import pandas as pd
import os
import numpy as np

PBMC = False
if PBMC:
  print('merging results for PBMC')
  os.chdir('/data/ashevtsov/MS_data/PBMC/cmap_res')
  merged_file_name = '/data/ashevtsov/MS_data/PBMC/PBMC_cmap_res_all.csv'
  name_stop = -18
else:
  print('merging results for CSF')
  os.chdir('/data/ashevtsov/MS_data/CSF/cmap_res')
  merged_file_name = '/data/ashevtsov/MS_data/CSF/CSF_cmap_res_all.csv'
  name_stop = -17
  
  
# merge molecules per all cell types
merged_df = pd.DataFrame(index = range(50))

for file in os.listdir():
  df = pd.read_csv(file)
  colname = file[1:name_stop]
  merged_df[colname] = df['pert_desc']

merged_df.to_csv(merged_file_name)
