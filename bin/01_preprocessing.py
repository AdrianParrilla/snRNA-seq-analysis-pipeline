#!/opt/env/bin/python

import os
import argparse
from glob import glob
import warnings
import numpy as np
import pandas as pd
import anndata as ad
import scanpy as sc
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.stats import binom
from itertools import product
import pickle
import logging 
import anndata2ri
from scipy.stats import median_abs_deviation as mad
import rpy2.rinterface_lib.callbacks as rcb
from rpy2.robjects.conversion import localconverter
import rpy2.robjects as ro

warnings.simplefilter("ignore", FutureWarning)

rcb.logger.setLevel(logging.ERROR)


# Load R libraries
ro.r('''
suppressMessages({
    library(Seurat)
    library(scater)
    library(scDblFinder)
    library(SingleCellExperiment)
    library(BiocParallel)
})
''')



def load_adatas(input_pkl):

    with open(input_pkl, 'rb') as f:
        adatas = pickle.load(f)

    return adatas


def mad_outlier(adata, metric: str, nmads: int, upper_only: bool = False):
    M = adata.obs[metric]

    if not upper_only:
        return (M < np.median(M) - nmads * mad(M)) | (M > np.median(M) + nmads * mad(M))

    return (M > np.median(M) + nmads * mad(M))


def flag_outliers(adata):

    bool_vector = (
        mad_outlier(adata, "log1p_total_counts", 5)
        | mad_outlier(adata, "log1p_n_genes_by_counts", 5)
        | mad_outlier(adata, "pct_counts_in_top_20_genes", 5)
        | mad_outlier(adata, "pct_counts_mt", 3, upper_only=True)
        | (adata.obs["pct_counts_mt"] > 5) # setting maximum tolerated mitochondrial percentage
        | (adata.obs["n_genes_by_counts"] < 100) # flagging cells with less than 100 genes
    )

    adata.obs['outlier'] = bool_vector.astype(bool)
    

    return adata


def detect_doublets(adata):
    '''
    Perform doublet detection with scDblFinder

    returns: The input AnnData object with two new .obs columns:
        - 'scDblFinder_score'
        - 'scDblFinder_class'
    '''

    data_mat = adata.X.T

    with localconverter(anndata2ri.converter):
        ro.globalenv["data_mat"] = data_mat
        data_mat = ro.r('as(data_mat, "dgCMatrix")')

    

    ro.r(f'''
    set.seed({123})
    sce <- scDblFinder(
        SingleCellExperiment(list(counts=data_mat))
    )
    ''')

    doublet_score = np.array(ro.r('sce$scDblFinder.score'))
    doublet_class = np.array(ro.r('sce$scDblFinder.class'))

    # Add to AnnData.obs
    adata.obs['scDblFinder_score'] = doublet_score
    adata.obs['scDblFinder_class'] = doublet_class

    adata.obs['scDblFinder_class'] = pd.Categorical(adata.obs['scDblFinder_class']).rename_categories({1: 'singlet', 2: 'doublet'})

    return adata


def qc(adata):
    """
    Perform QC and detect doublets on a given Anndata
    """

    adata.var_names_make_unique()

    # mitochondrial genes
    adata.var["mt"] = adata.var_names.str.startswith("MT-")
    # ribosomal genes
    adata.var["ribo"] = adata.var_names.str.startswith(("RPS", "RPL"))

    adata.var["hb"] = adata.var_names.str.contains(r"^HB[ABDEGMQZ]\d*(?!\w)")

    sc.pp.calculate_qc_metrics(
        adata, qc_vars=["mt", "ribo", "hb"], inplace=True, percent_top=[20], log1p=True
        )

    remove = ['total_counts_mt', 'log1p_total_counts_mt', 'total_counts_ribo', 'log1p_total_counts_ribo', 'total_counts_hb', 'log1p_total_counts_hb']   

    adata.obs.drop(columns=remove, inplace=True)

    adata = flag_outliers(adata)

    adata = detect_doublets(adata)

    return adata


def parse_metadata(adata, metadata_path, sample_key="sample"):
    """
    Merges metadata from a CSV into an AnnData object.
    """
    print(f"Parsing metadata from {metadata_path}...")
    
    try:
        # Load metadata
        samples_metadata = pd.read_csv(metadata_path)
        
        # Ensure the join key exists in both objects
        if sample_key not in adata.obs.columns:
            raise KeyError(f"'{sample_key}' not found in adata.obs")
        if sample_key not in samples_metadata.columns:
            raise KeyError(f"'{sample_key}' not found in metadata file")

        # Setting index on metadata allows for a cleaner join
        samples_metadata = samples_metadata.set_index(sample_key)
        
        # Efficiently map the new columns onto the observations
        adata.obs = adata.obs.join(samples_metadata, on=sample_key, how="left")

        for col in adata.obs.columns:
            if adata.obs[col].dtype == 'object':
                adata.obs[col] = adata.obs[col].astype('category')
        
        print("Metadata successfully integrated.")
        return adata

    except FileNotFoundError:
        print(f"Error: The file at {metadata_path} was not found.")
    except Exception as e:
        print(f"An error occurred: {e}")
        raise


def main(input_pkl,metadata_path, sample_key):

    adatas = load_adatas(input_pkl)

    adatas_dict = {}

    for sample,adt in adatas.items(): 
        print(f"Processing {sample}", flush=True)
        adatas_dict[sample] = qc(adt)

    # concatenate adatas
    adatas_concat = ad.concat(adatas_dict, label="sample") 
    adatas_concat.obs_names_make_unique()

    # remove genes express in less than 3 cells
    sc.pp.filter_genes(adatas_concat, min_cells=3)

    # remove deprecated genes
    adatas_concat = adatas_concat[:, ~adatas_concat.var.index.str.contains('DEPRECATED_')] 
    
    adatas_concat.obs["sample"] = adatas_concat.obs["sample"].str.replace(r"_P\d+$", "", regex=True) #remove _P0* tail from sample names 
    
    # setting cell names as column
    adatas_concat.obs['cell_name'] = adatas_concat.obs.index.values.copy()
    
    # set cell id as row names
    adatas_concat.obs = adatas_concat.obs.set_index("cell_name")

    #Add metadata
    adata = parse_metadata(adatas_concat, metadata_path, sample_key=sample_key)


    print("Saving adata after QC...")
    adatas_concat.write('01_adata_concat_QC.h5ad')
    print("QC completed!")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Perform QC of scRNA seq data')
    parser.add_argument('--data_pkl', required=True, help='Path to the pkl file containing a dictionary sample:Anndata')
    parser.add_argument('--metadata_path', required=True, help='Path to the sample metadata')
    parser.add_argument('--sample_key', required=True, help='Metadata field containing sample names')

    
    args = parser.parse_args()
    
    main(args.data_pkl,
        args.metadata_path,
        args.sample_key)