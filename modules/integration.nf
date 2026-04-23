process integration { 
    label "integration"
    stageInMode 'copy'
    
    input:
        path adata
        val filename
        val n_genes
        val label_key
        val batch_key

    output:
        path "${filename}_integration_benchmark.h5ad", emit: adata_integrated
        path "*_UMAP_*.png"                        , emit: umap_plots
        path "*_benchmark_results.csv"             , emit: metrics_csv, optional: true
        path "*.svg"                               , emit: scib_plot, optional: true

    script:
        """
        03_integration.py \\
            --adata_dir ${adata} \\
            --filename ${filename} \\
            --n_genes ${n_genes} \\
            --label_key ${label_key} \\
            --batch_key ${batch_key}
        """
        }