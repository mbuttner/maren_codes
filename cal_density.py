def cal_density(adata, mapping='X_tsne'):

    from scipy.stats import gaussian_kde

    x = adata.obsm[mapping][:,0]
    y = adata.obsm[mapping][:,1]


    # Calculate the point density
    xy = np.vstack([x,y])
    z = gaussian_kde(xy)(xy)
    dens_name = mapping[2:] + '_density'
    adata.obs[dens_name]=z

    return(adata)


def plot_density(adata,adata_full, mapping='X_tsne'):
    import matplotlib.pyplot as pl

    from scipy.stats import gaussian_kde
    dens_name = mapping[2:] + '_density'
    z = adata.obs[dens_name]
    x = adata.obsm[mapping][:,0]
    y = adata.obsm[mapping][:,1]
    idx = z.argsort()
    x, y, z = x[idx], y[idx], z[idx]

    fig, ax = pl.subplots()
    ax.scatter(adata_full.obsm[mapping][:,0],
               adata_full.obsm[mapping][:,1], c='grey', edgecolor='')
    ax.scatter(x, y, c=z, s=10, edgecolor='',cmap='Spectral_r')
    pl.grid(b=False)
    pl.show()


    return()
