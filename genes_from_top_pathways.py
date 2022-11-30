import pandas as pd
import os

MIN_GENES_IN_PATH = 1
sample = 'CSF'

folder_paths = [f"/data/ashevtsov/MS_data/{sample}/oncobox_res"]
out_top10_path_each_celltype_path = f'/data/ashevtsov/MS_data/{sample}/oncobox_res/merged/{sample}_top10_all_cells_bipathways_min{MIN_GENES_IN_PATH}genes.csv'
out_genes_top10_path_each_celltype_path = f'/data/ashevtsov/MS_data/{sample}/oncobox_res/merged/{sample}_genes_from_top10_all_cells_bipathways_min{MIN_GENES_IN_PATH}genes.csv'

def find_genes_by_path(path):
    
        genes_path_i = []
        # need to rename reactome and biocarta databases foldes to lowercase and KEGG adjasted to KEGG!
        if path.split('_')[0] in ['NCI', 'biocarta', 'reactome', 'KEGG']:
            df = pd.read_csv(os.path.join("/data/ashevtsov/oncoboxlib/databases/" + path.split('_')[0] + ' 1.123/arr.csv'))
            a = df[df[path] > 0]['gene'].to_numpy()
            for i in a:
                genes_path_i.append(i)

        else:
            try:
                df = pd.read_csv(os.path.join("/data/ashevtsov/oncoboxlib/databases/" + 'Metabolism' + ' 1.123/arr.csv'))
                a = df[df[path] > 0]['gene'].to_numpy()
                for i in a:
                    genes_path_i.append(i)
            except:
                df = pd.read_csv(os.path.join("/data/ashevtsov/oncoboxlib/databases/" + 'Qiagen 1.123' + '/arr.csv'))
                a = df[df[path] > 0]['gene'].to_numpy()
                for i in a:
                    genes_path_i.append(i)
        return genes_path_i
      
def make_list_genes_of_minNgene_pathways(path, N, ascending):
#     print(path)
    check = 0
    genes_for_each_path = {}
    
    df  = pd.read_csv(path)
    pathways = []
    df = df.sort_values(by='Tumour_geomean', ascending= ascending)
    i = 0
    while check < 10:
        biopath = df.iloc[i, 0]
        # print(biopath)
        i += 1
        genes_for_path = find_genes_by_path(biopath)
        if len(genes_for_path) >= N:
            pathways.append(biopath)
            print(f'Added {biopath}, number genes - {len(genes_for_path)}')
            genes_for_each_path[biopath] = genes_for_path
            check += 1
        else:
            print(biopath, 'NOT enough genes', 'got-', len(genes_for_path), 'need-', N)
    genes = []
    for i in genes_for_each_path.values():
        genes +=  i
    return genes_for_each_path, genes, pathways
# d = {}
# 
# 
# for folder_path in folder_paths:
#     for file in os.listdir(folder_path):
#         if 'merged' not in file:
#             file_path = os.path.join(folder_path, file)
#             results = pd.read_csv( file_path, sep = ',')
#             # print(results)
#             top_pos_pathways = list((results.sort_values('Tumour_geomean',ascending = False).head(10))['pathway'])
#             top_neg_pathways = list((results.sort_values('Tumour_geomean',ascending = True).head(10))['pathway'])
#             
#             d[file[:-4] + '_pos'] = top_pos_pathways
#             d[file[:-4] + '_neg'] = top_neg_pathways
#         
# df = pd.DataFrame(d)
# df.to_csv(out_top10_path_each_celltype_path)

genes_csf_all = {}
pathways_all = {}

for folder_path in folder_paths:
    for file in os.listdir(folder_path):
        if 'merged' not in file:
              print('-----------',file, '--------------')
              path = os.path.join(folder_path, file)
              
              genes_for_each_path, genes_pos, pathways_pos = make_list_genes_of_minNgene_pathways(path, MIN_GENES_IN_PATH, False)
              genes_csf_all[file[:-4] + '_pos'] = genes_pos
              pathways_all[file[:-4] + '_pos'] = pathways_pos
              
              genes_for_each_path, genes_neg, pathways_neg = make_list_genes_of_minNgene_pathways(path, MIN_GENES_IN_PATH, True)
              genes_csf_all[file[:-4] + '_neg'] = genes_neg
              pathways_all[file[:-4] + '_neg'] = pathways_neg

df_genes = pd.DataFrame(dict([ (k,pd.Series(v)) for k,v in genes_csf_all.items() ]))
df_genes.to_csv(out_genes_top10_path_each_celltype_path)

df_paths = pd.DataFrame(pathways_all)
df_paths.to_csv(out_top10_path_each_celltype_path)
