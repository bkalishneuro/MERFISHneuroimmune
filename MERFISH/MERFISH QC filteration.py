##This code was used to perform initial QC filtration, normalisation and doublet removal
 
import numpy as np
import pandas as pd
import scanpy as sc
import sys
import scrublet as scr
from matplotlib import pyplot as plt
from observable_jupyter import embed
from clustergrammer2 import net, Network, CGM2
from scipy.stats.stats import pearsonr   
from copy import deepcopy
from glob import glob
import anndata


def find_hdf5_files(main_directory):
    return [os.path.abspath(os.path.join(root, file)) for root, _, files in os.walk(main_directory) for file in fnmatch.filter(files, '*.hdf5')]

def filter_date_time(file_list):
    pattern = r'\d{2}-\d{2}-\d{2}_\d{2}-\d{2}'
    return [file for file in file_list if re.search(pattern, file)]

hdf5_files = filter_date_time(find_hdf5_files(MAIN_DIRECTORY))
hdf5_files.sort()
datasets = [hdf5_file for hdf5_file in hdf5_files]
print('Total Number of Files You will be Working on Today:', len(datasets))
print('File are sorted alphabetically, They are:')
for i in range(len(datasets)):
    line = datasets[i] + '                   file index: ' + str(i)
    print(line)
	
	
# load the file
anndata_filename = [x for x in glob('*.hdf5') if '_clustered' not in x][0]
print(anndata_filename)

# filter out blanks before clusering

keep_genes = [x for x in ad_viz.var.index.tolist() if 'Blank' not in x]

ad_viz = ad_viz[:, keep_genes]

#remove unclustered cells
ROI = ad_viz.obs['Custom cells groups']

keep_cells = ROI[(ROI != 'Unclustered')].index.tolist()

ad_viz = ad_viz[keep_cells]

ad_viz

#plot histogram
ser_volume = ad_viz.obs['volume']
ser_volume.hist(range = (0, 3000), bins=100)

# filter cells based on volume less than 50 and greater than 3X median
min_volume = 50
median = ad_viz.obs['volume'].median()
median
max_volume = median * 3
keep_cells = ser_volume[(ser_volume > min_volume) & (ser_volume < max_volume)].index.tolist()
ad_viz = ad_viz[keep_cells]
ad_viz

#Normalize by cell volume
ad_viz.to_df()
cell_by_gene = deepcopy(ad_viz.to_df())
cols_to_divide = [col for col in cell_by_gene.columns]

#add total raw transcript count for each cell in obs
ad_viz.obs['BarcodeCountOrig'] = cell_by_gene.sum(axis=1)
ad_viz.obs
cell_by_gene[cols_to_divide]= cell_by_gene[cols_to_divide].div(ad_viz.obs['volume'], axis=0)
cell_by_gene

#multiply Normalised count by 1000
cell_by_gene = cell_by_gene *1000
cell_by_gene

#add total normalised transcript count for each cell in obs
ad_viz.obs['BarcodeCountNor'] = cell_by_gene.sum(axis=1)
ad_viz.obs

ad_viz2 = anndata.AnnData(X=cell_by_gene, obs=ad_viz.obs , var = ad_viz.var) 
ad_viz2
#ad_viz.to_df = pd.DataFrame(cell_by_gene)
ad_viz2.to_df()

#filter by expression: filter out cells with gene expression zero
ser_exp = ad_viz2.obs['BarcodeCountOrig']
min_exp = 0
keep_cells = ser_exp[ser_exp > min_exp].index.tolist()
ad_viz2 = ad_viz2[keep_cells]
ad_viz2


#filter by expression values more than 98% quantile
 
quantilexpo98 = ad_viz2.obs['BarcodeCountNor'].quantile(0.98)

selectcells = ad_viz2.obs[(ad_viz2.obs.BarcodeCountNor <= quantilexpo98) & (ad_viz2.obs.BarcodeCountNor >= quantilexpo2)].index
print(quantilexpo2, quantilexpo98)
len(selectcells)

ad_viz2 = ad_viz2[selectcells].copy()
ad_viz2

#Remove doublets using Scrublet
counts_matrix = deepcopy(ad_viz2.to_df())
counts_matrix
scrub = scr.Scrublet(counts_matrix, expected_doublet_rate=0.18)
doublet_scores, predicted_doublets = scrub.scrub_doublets(min_counts=2, 
                                                          min_cells=3, 
                                                          min_gene_variability_pctl=85, 
                                                          n_prin_comps=30)
														  
scrub.plot_histogram();

#select threshold to remove 15-20% doublets and use below, here I have set it to 0.18
predicted_doublets18 = scrub.call_doublets(threshold=0.20)
scrub.plot_histogram();

#select singlets
singlets18 = np.where(predicted_doublets18 == False)[0]
len(singlets18)
ad_viz18_singlets = ad_viz2[singlets18, :]

#save count matrix as csv
df = ad_viz18_singlets.to_df()
#df.index = df.index.str[5:]
df.index = df.index.str

filtered_filecsv = anndata_filename.split('.hdf5')[0] + '_singlet_cellbygene.csv'
df.to_csv(f"{filtered_filecsv}", index=True)

#save spatial file as csv
dfobs = ad_viz18_singlets.obs
filtered_filecsv = anndata_filename.split('.hdf5')[0] + '_singlet_spatial.csv'
dfobs.to_csv(f"{filtered_filecsv}", index=True)


#save filtered anndata
filtered_filename = anndata_filename.split('.hdf5')[0] + '_filtered.hdf5'
ad_viz.write_h5ad(filtered_filename)