process plot_QC { 
    label "plot_QC" 

    input:
        path adata
        val sample_key
        val color_by

    output:
        path "QC_violin_plots.png",      emit: qc_violins
        path "QC_doublet_count.png",      emit: qc_doublets
        path 'QC_summary.csv',           emit: qc_summary
        path 'QC_thresholds.csv'       , emit: qc_thresholds
    

    script:
        """
        02_plot_QC.py --adata_dir ${adata} --sample_key ${sample_key} --color_by ${color_by}
        """
        }