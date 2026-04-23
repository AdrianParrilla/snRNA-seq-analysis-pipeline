#!/opt/conda/bin/python

import os
import numpy as np
import pandas as pd
import anndata as ad
import harmonypy as hm
import scanpy as sc
import scanorama
import scanpy.external as sce
from scvi.model import SCVI, SCANVI
from scib_metrics.benchmark import Benchmarker
import scib
import argparse
import torch  
import anndata2ri
import rpy2.rinterface_lib.callbacks as rcb
from rpy2.robjects.conversion import localconverter
import rpy2.robjects as ro


# Load R libraries
ro.r('''
suppressMessages({
    library(Seurat)
    library(SingleCellExperiment)
    library(STACAS)
})

# Increase max size for parallelization to ~100 GB to handle large datasets
options(future.globals.maxSize = 100 * 1024^3)
''')


# Code adapted from https://github.com/openproblems-bio/task_batch_integration/tree/main/src/methods

def get_hvg_batch(adata, layer, batch_key, n_genes):

    if n_genes > adata.n_vars or n_genes <= 0:
        print(f"\n>>> Using all {adata.n_vars} features as batch-aware HVGs...", flush=True)
        adata.var["highly_variable"] = True 
        adata.var["highly_variable_nbatches"] = adata.obs[batch_key].nunique()

    else:
        print(f"\n>>> Computing {n_genes} batch-aware HVGs...", flush=True)

        # seurat_v3 requires raw counts!!
        sc.pp.highly_variable_genes(
            adata, n_top_genes=n_genes, flavor='seurat_v3', layer=layer, batch_key=batch_key, subset=False
            )

        var_genes = adata.var.index[adata.var["highly_variable"]].tolist()

        print(f"Number of highly variable genes = {len(var_genes)}")

    return adata
 

def preprocess_adata(adata, layer, batch_key, n_genes):

    adata.layers[layer] = adata.X.copy() #save raw counts

    # Exclude outliers
    adata = adata[~adata.obs['outlier']].copy()

    adata = get_hvg_batch(adata, layer, batch_key, n_genes)

    sc.pp.normalize_total(adata, target_sum=None, inplace=True) # normalize from counts in case adata.X is scaled
    sc.pp.log1p(adata)

    sc.pp.pca(adata, use_highly_variable=True, key_added= 'Unintegrated')
    
    return adata


def merge_adata(*adata_list, **kwargs):
    """
    Merge adatas from list while remove duplicated ``obs`` and ``var`` columns.
    Based on scib -> https://github.com/theislab/scib/blob/59ae6eee5e611d9d3db067685ec96c28804e9127/scib/utils.py#L51C1-L72C62

    :param adata_list: ``anndata`` objects to be concatenated
    :param kwargs: arguments to be passed to ``anndata.AnnData.concatenate``
    """

    if len(adata_list) == 1:
        return adata_list[0]

    # Make sure that adatas do not contain duplicate columns
    for _adata in adata_list:
        for attr in ("obs", "var"):
            df = getattr(_adata, attr)
            dup_mask = df.columns.duplicated()
            if dup_mask.any():
                print(
                    f"Deleting duplicated keys `{list(df.columns[dup_mask].unique())}` from `adata.{attr}`."
                )
                setattr(_adata, attr, df.loc[:, ~dup_mask])

    return ad.AnnData.concatenate(*adata_list, **kwargs) 


def dim_reduction_plotting(adata, filename, batch_key, reduction:str):

    sc.pp.neighbors(adata, use_rep=reduction)
    sc.tl.umap(adata, min_dist=0.2, spread=5, random_state=42)

    axes = sc.pl.umap(adata, color='scDblFinder_class', show=False)
    fig = axes.get_figure()
    fig.savefig(f'{filename}_UMAP_{reduction}.png')

    print(f'\nUMAP plot for {reduction} succesfully generated!', flush=True)


def plot_doublets(adata, filename):

    sc.pp.neighbors(adata, use_rep='Unintegrated')
    sc.tl.umap(adata, min_dist=0.2, spread=5, random_state=42)

    axes = sc.pl.umap(adata, color=batch_key, show=False)
    fig = axes.get_figure()
    fig.savefig(f'{filename}_UMAP_{reduction}.png')

    # exclude doublets
    adata = adata[adata.obs['scDblFinder_class'] == 'singlet']

    print(f'\nUMAP plot for {reduction} succesfully generated!', flush=True)



def harmony_integration(adata, adata_integ, batch_key):

    print('\n>>> Running Harmony', flush=True)

    harmony_out = hm.run_harmony(
        adata_integ.obsm['Unintegrated'], 
        adata_integ.obs, 
        batch_key
    )

    adata.obsm['Harmony'] = harmony_out.Z_corr 

    print('\n>>> Harmony integration done!', flush=True)

    return adata


def scanorama_integration(adata, adata_integ, batch_key):

    print('\n>>> Running Scanorama', flush=True)

    if not isinstance(adata_integ.obs[batch_key].dtype, pd.CategoricalDtype):
        adata_integ.obs[batch_key] = adata_integ.obs[batch_key].astype('category')

    split = []
    for i in batch_categories:
        split.append(adata_integ[adata_integ.obs[batch_key] == i].copy())
    
    scanorama.integrate_scanpy(split, dimred = 50)

    scanorama_int = [ad.obsm['X_scanorama'] for ad in split]

    all_s = np.concatenate(scanorama_int)

    adata.obsm["Scanorama"] = all_s

    print('\n>>> Scanorama integration done!', flush=True)
    
    return adata


def seurat_RPCA(adata, adata_integ, layer, batch_key):
    """
    Seurat RPCA returns the integrated matrix instead of the embbeding
    """
    print('\n>>> Running Seurat RPCA integration', flush=True)

    adata_integ.obs[batch_key] = adata_integ.obs[batch_key].astype(str)

    adata_integ.X = adata_integ.layers[layer] # setting back raw counts for seurat pipeline

    if adata_integ.uns:
        del adata_integ.uns

    with localconverter(anndata2ri.converter):
        ro.globalenv["seurat"] = adata_integ
        ro.globalenv["batch_key"] = batch_key

    ro.r('''
        assayNames(seurat)[1] <- "counts"
        
        seurat_obj <- as.Seurat(seurat, counts = "counts", data = NULL)
         
        # Get the active assay name dynamically to be safe.
        assay_name <- DefaultAssay(seurat_obj)


        # Seurat v5 requires layers to be split by batch BEFORE integration
        seurat_obj[[assay_name]] <- split(seurat_obj[[assay_name]], f = seurat_obj@meta.data[[batch_key]])
        
        # RPCA requires PCA to be run first)
        seurat_obj <- NormalizeData(seurat_obj, verbose = FALSE)
        seurat_obj <- FindVariableFeatures(seurat_obj, verbose = FALSE)
        seurat_obj <- ScaleData(seurat_obj, verbose = FALSE)
        seurat_obj <- RunPCA(seurat_obj, verbose = FALSE)
        
        # Run RPCA Integration
        seurat_obj <- IntegrateLayers(
            object = seurat_obj,
            method = RPCAIntegration,
            orig.reduction = "pca",
            new.reduction = "integrated.rpca",
            verbose = FALSE
        )
        
        # Extract the integrated cell embeddings
        rpca_embeddings <- Embeddings(seurat_obj, "integrated.rpca")
        ''')

    # Fetch the embeddings back to Python
    embeddings = ro.r('rpca_embeddings')
    
    adata.obsm["Seurat_RPCA"] = np.array(embeddings)

    # Clean up R environment to free memory
    ro.r('rm(seurat, seurat_obj, rpca_embeddings); gc()')

    print('\n>>> Seurat RPCA integration done!', flush=True)

    return adata


def STACAS(adata, adata_integ, layer, batch_key):

    print('\n>>> Running STACAS integration', flush=True)

    adata_integ.obs[batch_key] = adata_integ.obs[batch_key].astype(str)

    adata_integ.X = adata_integ.layers[layer] # setting back raw counts for seurat pipeline

    if adata_integ.uns:
        del adata_integ.uns

    with localconverter(anndata2ri.converter):
        ro.globalenv["seurat"] = adata_integ
        ro.globalenv["batch_key"] = batch_key

        ro.r('''

        assayNames(seurat)[1] <- "counts"

        seurat_obj <- as.Seurat(seurat, counts = "counts", data = NULL)
                 
        seurat_obj <- NormalizeData(seurat_obj, verbose = FALSE)
             
        seurat_list <- SplitObject(seurat_obj, split.by = batch_key)
             
        stacas_obj <- Run.STACAS(seurat_list, verbose = FALSE)
  
        stacas_embeddings <- Embeddings(stacas_obj, "pca")
             

        ''')

    embeddings = ro.r('stacas_embeddings')
    adata.obsm["STACAS"] = np.array(embeddings)

    ro.r('rm(seurat, seurat_obj, seurat_list, stacas_obj, stacas_embeddings); gc()')

    print('\n>>> STACAS integration done!', flush=True)

    return adata


def scVI_integration(adata, adata_integ, layer, batch_key):
    '''
    For this method, raw counts must be used
    '''

    print('\n>>> Running scVI', flush=True)

    # detect device
    device = "gpu" if torch.cuda.is_available() else "cpu"
    print(f"Using device {device}", flush=True)

    SCVI.setup_anndata(adata_integ, layer= layer, batch_key= batch_key) # set up model
    model_scvi = SCVI(adata_integ, n_layers=2, n_latent=30)

    model_scvi.train(accelerator= device, early_stopping=True)

    adata.obsm["scVI"] = model_scvi.get_latent_representation() 

    print('\n>>> scVI integration done!', flush=True)

    return adata, model_scvi


def scANVI_integration(adata, adata_integ, model_scvi, label_key, batch_key):

    print('\n>>> Running scANVI', flush=True)

    # detect device
    device = "gpu" if torch.cuda.is_available() else "cpu"
    print(f"Using device {device}", flush=True)

    model_scanvi = SCANVI.from_scvi_model(model_scvi,
    adata=adata_integ,
    labels_key= label_key,
    unlabeled_category="Unknown"
    )

    model_scanvi.train(max_epochs=20, n_samples_per_label=100)

    adata.obsm["scANVI"] = model_scanvi.get_latent_representation()

    print('\n>>> scANVI integration done!', flush=True)

    return adata


def integration_metrics(adata, batch_key, label_key, filename):

    print('\n>>> Calculating integration metrics', flush=True)

    valid_embeds = [emb for emb in adata.obsm.keys() if 'umap' not in emb.lower()]

    print(f"Detected embeddings: {valid_embeds}")

    bm = Benchmarker(
        adata,
        batch_key=batch_key,
        label_key=label_key,
        embedding_obsm_keys=valid_embeds,
        n_jobs=1
        )

    bm.benchmark()

    metrics = bm.get_results(min_max_scale=False)

    metrics.to_csv(f"{filename}_benchmark_results.csv")

    bm.plot_results_table(
        min_max_scale=False, 
        show=False,
        save_dir=".")

    print('\n>>> Metrics calculation done!', flush=True)



def main(adata_dir, n_genes= 2000, layer='counts', label_key='cell_type', batch_key='batch', filename = 'adata'):

    print('\n>>> Loading adata...', flush=True)
    adata = sc.read_h5ad(adata_dir)
    print('\n>>> adata succesfully loaded', flush=True)


    if batch_key not in adata.obs.columns:
        raise ValueError(f"Error: batch_key '{batch_key}' not present in adata.obs")
    else:
        adata.obs[batch_key] = adata.obs[batch_key].astype('category')
    
    if label_key in adata.obs.columns:
        adata.obs[label_key] = adata.obs[label_key].astype('category')


    # ---------- Preprocess adata --------------
    adata = preprocess_adata(adata, layer, batch_key, n_genes)

    adata_integ = adata[:, adata.var["highly_variable"]].copy() # subset adata to HVG

    # ---------- Harmony integration ------------
    adata = harmony_integration(adata, adata_integ, batch_key)

    # ---------- Seurat RPCA integration ---------
    adata = seurat_RPCA(adata, adata_integ, layer, batch_key)

    # ----------- STACAS integration ---------
    adata = STACAS(adata, adata_integ, layer, batch_key)
    
    # ---------- Scanorama integrarion ----------
    adata = scanorama_integration(adata, adata_integ, batch_key)


    # ---------- scVI integrarion ---------------
    adata, model_scvi = scVI_integration(adata, adata_integ, layer, batch_key)


    lk_present = label_key in adata.obs.columns

    if not lk_present:
        print(f"Label key: {lk_present} not present in adata.obs, skipping scANVI integration and metrics calculation.", flush=True)
    else:
        adata = scANVI_integration(adata, adata_integ, model_scvi, label_key, batch_key)
        integration_metrics(adata, batch_key, label_key, filename)

    
    adata.write(f"{filename}_integration_benchmark.h5ad")
    print(f'Adata saved!', flush=True)

    plot_methods = list(adata.obsm.keys())

    for reduction in plot_methods:
        dim_reduction_plotting(adata, filename, batch_key, reduction)



    print("\n>>> Integration benchmark completed!", flush=True)



if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='Integrate scRNA seq data from different batches using scVI')
    parser.add_argument('-i', '--adata_dir', type=str, required=True, help='Input directory of the h5ad file containing the anndata object.')
    parser.add_argument('-lk', '--label_key', type=str, required=False, default='cell_type', help='Column from adata.obs where the cell types are specified. Default: cell_type')
    parser.add_argument('-bk', '--batch_key', type=str, required=False, default='batch', help='Column from adata.obs where the sample batch is specified. Default: batch')
    parser.add_argument('-n', '--n_genes', type=int, required=False, default=2000, help='Number of high variable genes to select. Default: 2000')
    parser.add_argument('-f', '--filename', type=str, required=False, default='adata', help='Filename for the output files')

    args = parser.parse_args()

    print(
        f"adata_dir: {args.adata_dir}\n"
        f"n_genes: {args.n_genes}\n"
        f"label_key: {args.label_key}\n"
        f"batch_key: {args.batch_key}\n"
        f"filename: {args.filename}"
    )

    main(
     adata_dir= args.adata_dir,
     n_genes= args.n_genes,
     label_key= args.label_key,
     batch_key= args.batch_key, 
     filename=args.filename
     )

