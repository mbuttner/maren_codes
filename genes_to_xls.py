

def genes_to_xls(adata=None, filename='genes_hvg_regressed.xlsx'):
    import xlsxwriter
    import pandas as pd
    writer = pd.ExcelWriter(filename, engine='xlsxwriter')
    for i in adata.uns[ 'rank_genes_groups_gene_names'].dtype.names:
        d = {'gene':adata.uns['rank_genes_groups_gene_names'][i],
             'score':adata.uns['rank_genes_groups_gene_scores'][i], 
             'mean':adata.var['mean'][adata.uns['rank_genes_groups_gene_names'][i]]} #note the mean comes with an index.
        df = pd.DataFrame(data=d)
        df.to_excel(writer,sheet_name=i)
    writer.save()
    return
