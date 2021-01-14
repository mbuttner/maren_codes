import scanpy as sc
import numpy as np
import matplotlib.pyplot as pl

def check_barcodes(adata, batch_key='sample', plot_show=True):
    #pairwise comparison of all barcodes per sample and visualise the fraction of overlapping barcodes as heatmap
    barcodes_list = {}
    for sample in adata.obs[batch_key].cat.categories:
        tmp = adata.obs_names[np.in1d(adata.obs[batch_key],sample)].values
        tmp_2 = [idx.split('-')[0] for idx in tmp]
        barcodes_list[sample] = tmp_2
    
    keys = barcodes_list.keys()
    corr_mat = np.zeros((len(keys), len(keys)))
    for key in enumerate(keys):
        for key2 in enumerate(keys):
            corr_mat[key[0], key2[0]] = len(np.intersect1d(barcodes_list[key[1]], 
                                                       barcodes_list[key2[1]]))/len(barcodes_list[key[1]])
    
    if plot_show:
        fig, ax = pl.subplots()
        im = pl.imshow(corr_mat)
        ax.set_xticks(np.arange(len(keys)))
        ax.set_yticks(np.arange(len(keys)))
        ax.set_xticklabels(keys)
        ax.set_yticklabels(keys)
        pl.ylabel('Query')
        pl.xlabel('Reference')
        pl.setp(ax.get_xticklabels(), rotation=45, ha="right",
                 rotation_mode="anchor")
        pl.show()
    
    return corr_mat
