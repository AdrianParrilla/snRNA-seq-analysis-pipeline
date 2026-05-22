#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// Include modules
include { check_cb_samplesheet } from './modules/check_cb_samplesheet'
include { cellbender } from './modules/cellbender'
include { demultiplex } from './modules/demultiplex'
include { run_QC } from './modules/run_QC'
include { merge_adatas } from './modules/merge_adatas'
include { plot_QC } from './modules/plot_QC'
include { integration } from './modules/integration'
include { integration_metrics } from './modules/integration_metrics'
include { integration_plots } from './modules/integration_plots'



workflow CELLBENDER {
    main:

        def out_dir   = "${params.outdir}/sc_processing"
        
        
        check_cb_samplesheet(params.samplesheet)
        
        ch_data_dirs = check_cb_samplesheet.out.validated_csv
            .splitCsv(header: true)
            .map { row -> 
                def alignment_dir = file(row.alignment_dir, type: 'dir', checkIfExists: true)
                
                def lr    = row.learning_rate ?: null
                def exp   = row.expected_cells ?: null
                def total = row.total_droplets_included ?: null
                
                return tuple(row.dataset, alignment_dir, out_dir, lr, exp, total) 
            }
      
        cellbender(ch_data_dirs)
        
       
        ch_data_dirs = channel
            .fromPath(params.samplesheet)
            .splitCsv(header: true)
            .map { row ->
                def alignment_dir = file(row.alignment_dir, type: 'dir')
                
                return tuple(row.dataset, alignment_dir, out_dir)
            }
        
        ch_demux_input = ch_data_dirs.join(cellbender.out.cb_matrix_raw)

        

        demultiplex(ch_demux_input)
        
}



workflow QC_WORKFLOW {
    main:

        def out_dir   = "${params.outdir}/sc_processing"

        ch_h5ad_files = channel
                    .fromPath("${params.outdir}/sc_processing/checkpoints/samples_demultiplexed/*.h5ad", checkIfExists: true)
                    .map { it -> 
                        return tuple(it.baseName, it, out_dir) 
                    }

        run_QC(ch_h5ad_files)
        merge_adatas(run_QC.out.qc_h5ad.collect(), params.metadata, out_dir, params.merge_adatas.sample_key_merge, params.filename)
        plot_QC(merge_adatas.out.merged_adata, out_dir, params.plot_QC.sample_key_plot, params.plot_QC.color_by, params.filename)
        
}


workflow INTEGRATION_WORKFLOW {
    main:
        def out_dir = "${params.outdir}/sc_processing"

        adata_annotated = channel.fromPath(params.integration.adata_annotated)

        integration(adata_annotated, params.filename, params.integration.n_genes, params.integration.label_key, params.integration.batch_key, out_dir)
        integration_metrics(integration.out.adata_integrated, params.filename, params.integration.label_key, params.integration.batch_key, out_dir)
        integration_plots(integration_metrics.out.metrics_csv, out_dir)
}


workflow {
    if (params.run_workflow == 'QC') {
        log.info "Executing QC Workflow..."
        QC_WORKFLOW()
        
    } else if (params.run_workflow == 'integration') {
        log.info "Executing Integration Workflow..."
        INTEGRATION_WORKFLOW()

    } else if (params.run_workflow == 'cellbender') {
        log.info "Executing Cellbender Workflow..."
        CELLBENDER()
        
    } else {
        exit 1, "ERROR: Invalid workflow specified. Set params.run_workflow to 'QC', 'integration' or 'cellbender' in your config."
    }
}
