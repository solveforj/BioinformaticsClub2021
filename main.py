import scanpy as sc
import pandas as pd

# Set up data loading

#Data files
sample_strings = ['Duo_M1', 'Duo_M2', 'Jej_M1', 'Jej_M2', 'Il_M1', 'Il_M2']
sample_id_strings = ['3', '4', '5', '6', '7', '8']
file_base = 'GSM283657'
exp_string = '_Regional_'
data_file_end = '_matrix.mtx'
barcode_file_end = '_barcodes.tsv'
gene_file_end = '_genes.tsv'
cc_genes_file = '../Macosko_cell_cycle_genes.txt'

sample = sample_strings.pop(0)
sample_id = sample_id_strings.pop(0)
data_file = file_base+sample_id+exp_string+sample+data_file_end
barcode_file = file_base+sample_id+exp_string+sample+barcode_file_end
gene_file = file_base+sample_id+exp_string+sample+gene_file_end

######################################################################
adata = sc.read("GSM2836573_Regional_Duo_M1_matrix.mtx.gz", cache=True)
adata = adata.transpose()
adata.X = adata.X.toarray()

barcodes = pd.read_csv("GSM2836573_Regional_Duo_M1_barcodes.tsv.gz", header=None, sep='\t')
genes = pd.read_csv("GSM2836573_Regional_Duo_M1_genes.tsv.gz", header=None, sep='\t')
barcodes.rename(columns={0:'barcode'}, inplace=True)
barcodes.set_index('barcode', inplace=True)
adata.obs = barcodes
adata.obs['sample'] = [sample]*adata.n_obs
adata.obs['region'] = [sample.split("_")[0]]*adata.n_obs
adata.obs['donor'] = [sample.split("_")[1]]*adata.n_obs
genes.rename(columns={0:'gene_id', 1:'gene_symbol'}, inplace=True)
genes.set_index('gene_symbol', inplace=True)
adata.var = genes

for i in range(len(sample_strings)):
    # Parse Filenames
    sample = sample_strings[i]
    sample_id = sample_id_strings[i]
    data_file = file_base + sample_id + exp_string + sample + data_file_end
    barcode_file = file_base + sample_id + exp_string + sample + barcode_file_end
    gene_file = file_base + sample_id + exp_string + sample + gene_file_end

    print(sample, sample_id, data_file, barcode_file, gene_file)


