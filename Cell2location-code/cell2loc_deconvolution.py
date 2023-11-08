# -*- coding: utf-8 -*-
"""
Created on Wed Oct 12 23:09:59 2022

@author: alban
"""

#%%
import scanpy as sc

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt 
import matplotlib as mpl
import os

os.environ["THEANO_FLAGS"] = 'device=cuda,floatX=float32,force_device=True'
import cell2location


from matplotlib import rcParams
rcParams['pdf.fonttype'] = 42 # enables correct plotting of text for PDFs


#import sys


#%%


wd = r"C:\Users\alban\OneDrive\Documents\GitHub\Deconvolution-for-ST-in-Cardiorenal-Disease\Cell2location-code"


#results_folder = 'results/2022_10_12_atherosub4000/'
results_folder = os.path.join(wd, 'C2L_results')
try: os.mkdir(results_folder)
except: print("Directory already exists")

os.chdir(results_folder)
# create paths and names to results folders for reference regression and cell2location models
ref_run_name = f'{results_folder}\\reference_signatures'

ref_path = r"C:\Users\alban\OneDrive\Documents\GitHub\Deconvolution-for-ST-in-Cardiorenal-Disease\R-code\athero_subset3000.h5ad"
annotation_key = 'Final_annotation'

#%%

# Load in single cell subsampled data
adata_ref = sc.read(ref_path)


adata_ref.var['SYMBOL'] = adata_ref.var.index 
# rename 'GeneID-2' as necessary for your data
#adata_ref.var.set_index('features', drop=True, inplace=True)

# delete unnecessary raw slot (to be removed in a future version of the tutorial)
#del adata_ref.raw


from cell2location.utils.filtering import filter_genes
selected = filter_genes(adata_ref, cell_count_cutoff=5, cell_percentage_cutoff2=0.03, nonz_mean_cutoff=1.12)


#%%
# filter the object
adata_ref = adata_ref[:, selected].copy()

# prepare anndata for the regression model
cell2location.models.RegressionModel.setup_anndata(adata=adata_ref, 
                        # 10X reaction / sample / batch
                        batch_key='orig.ident', 
                        # cell type, covariate used for constructing signatures
                        labels_key=annotation_key, 
                        # multiplicative technical effects (platform, 3' vs 5', donor effect)
                        categorical_covariate_keys=None
                       )


# create the regression model
from cell2location.models import RegressionModel
mod = RegressionModel(adata_ref) 

# view anndata_setup as a sanity check
mod.view_anndata_setup()

mod.train(max_epochs=500, batch_size = 2500, use_gpu=True)

mod.plot_history(20)

adata_ref = mod.export_posterior(
    adata_ref, sample_kwargs={'num_samples': 1000, 'batch_size': 2500, 'use_gpu': True}
)


# Save model
mod.save(f"{ref_run_name}", overwrite=True)

# Save anndata object with results
adata_file = f"{ref_run_name}/sc.h5ad"
adata_ref.write(adata_file)
adata_file


#%%

# export estimated expression in each cluster
if 'means_per_cluster_mu_fg' in adata_ref.varm.keys():
    inf_aver = adata_ref.varm['means_per_cluster_mu_fg'][[f'means_per_cluster_mu_fg_{i}' 
                                    for i in adata_ref.uns['mod']['factor_names']]].copy()
else:
    inf_aver = adata_ref.var[[f'means_per_cluster_mu_fg_{i}' 
                                    for i in adata_ref.uns['mod']['factor_names']]].copy()
inf_aver.columns = adata_ref.uns['mod']['factor_names']
inf_aver.iloc[0:5, 0:5]

#%%

### Update with appropriate sample names
# CVD
samples_list = ['CVD1', 'CVD2', 'CVD3', 'CVD4', 'CVD5',
                'CVD6', 'CVD7', 'CVD8', 'CVD9', 'CVD10']

# Data directory, one folder per sample
data_dir = r"C:\Users\alban\Data\CVD"


    
for sample in samples_list:
    
    print(sample)
    folder = results_folder
    print(folder)
    run_name = os.path.join(folder, sample)
    print(run_name)
    

    adata_vis = sc.read_visium(path = os.path.join(wd, data_dir, sample),
                           count_file = os.path.join(sample + ' filtered_feature_bc_matrix.h5'),
                           load_images = True,
                           source_image_path =  os.path.join(wd, data_dir,
                                                             sample + '/spatial/tissue_hires_image.png/'))



    adata_vis.obs['n_counts'] = adata_vis.X.sum(axis=1).A1
    # but plots well
    plt.rcParams["figure.figsize"] = (12, 12)
    sc.pl.spatial(adata_vis, img_key="hires", color=["n_counts"],
                  save=os.path.join(sample + '_vis_1_nCounts.png'))

    #plt.savefig(os.path.join(results_folder, sample_name + '_vis_1_nCounts.png'))
    
    
    # Add the sample name to every cell observation metadata
    adata_vis.obs['sample'] = sample
    # Duplicate the genes names (rownames) from the metadata to a column
    adata_vis.var['SYMBOL'] = adata_vis.var_names


    # Rename gene_ids column to ENSEMBL
    adata_vis.var.rename(columns={'gene_ids': 'ENSEMBL'}, inplace=True)
    # Replace the row names with the ENSEMBL column
    adata_vis.var_names = adata_vis.var['ENSEMBL']
    # Drop the ENSEMBL column
    adata_vis.var.drop(columns='ENSEMBL', inplace=True)
    
    
    
    # find mitochondria-encoded (MT) genes
    adata_vis.var['MT_gene'] = [gene.startswith('MT') for gene in adata_vis.var['SYMBOL']]
        
    # remove MT genes for spatial mapping (keeping their counts in the object)
    adata_vis.obsm['MT'] = adata_vis[:, adata_vis.var['MT_gene'].values].X.toarray()
    adata_vis = adata_vis[:, ~adata_vis.var['MT_gene'].values]


    #mod.plot_QC()

    
    # find shared genes and subset both anndata and reference signatures
    # exchange ensembl ids adata_vis.var_names for adata_vis.var.SYMBOL as ensembl ids are missing in sc data (alsaigh)
    adata_vis.var_names = adata_vis.var.SYMBOL


    # Remove duplicate genes to allow intersection, 2 genes ['TBCE', 'HSPA14']    
    
    dup = adata_vis.var.index.duplicated(keep='first')
    
    adata_vis.var['dup'] = dup
    
    adata_vis = adata_vis[:, ~adata_vis.var['dup'].values]


    intersect = np.intersect1d(adata_vis.var_names, inf_aver.index)
    adata_vis = adata_vis[:, intersect].copy()
    inf_aver = inf_aver.loc[intersect, :].copy()


    # prepare anndata for cell2location model
    cell2location.models.Cell2location.setup_anndata(adata=adata_vis, batch_key="sample")


    # create and train the model
    mod = cell2location.models.Cell2location(
        adata_vis, cell_state_df=inf_aver, 
        # the expected average cell abundance: tissue-dependent 
        # hyper-prior which can be estimated from paired histology:
        N_cells_per_location=8,
        # hyperparameter controlling normalisation of
        # within-experiment variation in RNA detection:
        detection_alpha=20
    ) 
    mod.view_anndata_setup()

    #max_epochs=16000
    mod.train(max_epochs=20000, 
              # train using full data (batch_size=None)
              batch_size=None, 
              # use all data points in training because 
              # we need to estimate cell abundance at all locations
              train_size=1,
              use_gpu=True)
    
    # plot ELBO loss history during training, removing first 100 epochs from the plot
    mod.plot_history(100)
    plt.legend(labels=['full data training']);

    
    # In this section, we export the estimated cell abundance (summary of the posterior distribution).
    adata_vis = mod.export_posterior(
        adata_vis, sample_kwargs={'num_samples': 1000, 'batch_size': mod.adata.n_obs, 'use_gpu': True}
    )
    
    # Save model
    mod.save(f"{run_name}", overwrite=True)
    
    # mod = cell2location.models.Cell2location.load(f"{run_name}", adata_vis)
    
    # Save anndata object with results
    adata_file = f"{run_name}/sp.h5ad"
    adata_vis.write(adata_file)
    adata_file


    mod.plot_QC()
    

    # add 5% quantile, representing confident cell abundance, 'at least this amount is present',
    # to adata.obs with nice names for plotting
    adata_vis.obs[adata_vis.uns['mod']['factor_names']] = adata_vis.obsm['q05_cell_abundance_w_sf']
    
    factor_names = adata_vis.obs[adata_vis.uns['mod']['factor_names']]
    adata_vis.obs[adata_vis.uns['mod']['factor_names']] = factor_names
    
    obs = adata_vis.obs
    adata_vis2 = adata_vis    
    adata_vis2.obs = obs

    # select one slide
   # from cell2location.utils import select_slide
   # slide = select_slide(adata_vis, 'FW106008')

    names  = [s.strip('q05_cell_abundance_w_sf') for s in adata_vis.obsm['q05_cell_abundance_w_sf']]

    # plot in spatial coordinates
    with mpl.rc_context({'axes.facecolor':  'black',
                         'figure.figsize': [4.5, 5]}):
    
        sc.pl.spatial(adata_vis, cmap='magma',
                      # show 9 cell types
                      color=names,
                      #color=["B_cells", "Endothelial_cells", "Fibroblasts", "Macrophages_and_MDSCs", "NK_and_T_cells", "SMCs_Fibroblasts_Hybrids", "Mast_cells", "SMCs_maintype", "SMCs_subtype"], 
                      ncols=9, size=1.3, 
                      img_key='hires',
                      # limit color scale at 99.2% quantile of cell abundance
                      vmin=0, vmax='p99.2',
                      save=os.path.join(sample + '_vis_all_celltypes.png')
                      )


    # Now we use cell2location plotter that allows showing multiple cell types in one panel
    from cell2location.plt import plot_spatial
    
    # select up to 6 clusters 
    clust_labels = ["Endothelial_cells", "Fibroblasts", "Macrophages_and_MDSCs", "SMCs_Fibroblasts_Hybrids", "SMCs_maintype", "SMCs_subtype"]
    clust_col = ['' + str(i) for i in clust_labels] # in case column names differ from labels
    
    #slide = select_slide(adata_vis, 'V1_Human_Lymph_Node')
    
    with mpl.rc_context({'figure.figsize': (16, 16)}):
        fig = plot_spatial(
            adata=adata_vis, 
            # labels to show on a plot
            color=clust_col, labels=clust_labels, 
            show_img=True,
            # 'fast' (white background) or 'dark_background'
            style='fast', 
            # limit color scale at 99.2% quantile of cell abundance
            max_color_quantile=0.992,
            # size of locations (adjust depending on figure size)
            circle_diameter=6, 
            colorbar_position='right'
        )
        
        plt.savefig(os.path.join(results_folder, sample + '_multiplot' + '.png'), dpi=150)
        
    with mpl.rc_context({'figure.figsize': (16, 16)}):
        fig = plot_spatial(
            adata=adata_vis, 
            # labels to show on a plot
            color=clust_col, labels=clust_labels, 
            show_img=True,
            # 'fast' (white background) or 'dark_background'
            style='fast', 
            # limit color scale at 99.2% quantile of cell abundance
            max_color_quantile=0.9,
            # size of locations (adjust depending on figure size)
            circle_diameter=6, 
            colorbar_position='right'
        )

    adata_vis = mod.export_posterior(
        adata_vis, sample_kwargs={'num_samples': 1000, 'batch_size': mod.adata.n_obs, 'use_gpu': True}
    )
    print(adata_vis)
    adata_vis.obsm['q05_cell_abundance_w_sf'].to_csv(os.path.join(results_folder, sample + '_Cell2location_result.txt'))

#%%