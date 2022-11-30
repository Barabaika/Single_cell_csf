import os
import subprocess


'''
1- untar folder with raw data in PATH (tar -xvf GSE138266_RAW.tar)
2- run this code, it with split .gz files with correct names to folders with names from file names 0:3  

'''
PATH = '/data/ashevtsov/MS_data/raw_data'
os.chdir(PATH)

batch_names = []

for file in os.listdir('./'):
  if '.gz' in file:
    print(file)
    batch_name = file.split('_')[:3]
    batch_name = '_'.join(batch_name)
    
    if batch_name not in batch_names:
      batch_names.append(batch_name)
      
      subprocess.run(["mkdir", batch_name])
      # bashCommand = "mkdir $batch_name"
      # process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE)
      # output, error = process.communicate()
      # !mkdir $batch_name
    
    if 'matrix' in file:
      mv_and_rename = batch_name + '/matrix.mtx.gz'
    elif 'barcodes' in file:
      mv_and_rename = batch_name + '/barcodes.tsv.gz'
    elif 'gene' in file:
      mv_and_rename = batch_name + '/features.tsv.gz'
      
    subprocess.run(["mv", file, mv_and_rename])
    # bashCommand = "mv $file $mv_and_rename"
    # process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE)
    # output, error = process.communicate()
    # !mv $file $mv_and_rename
    
