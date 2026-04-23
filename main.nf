#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// Include modules
include { parse_cb_samplesheet } from './modules/parse_cb_samplesheet'
include { cellbender } from './modules/cellbender'
include { QC } from './modules/QC'
include { plot_QC } from './modules/plot_QC'
include { integration } from './modules/integration'

workflow {
    main:
    input_mat = channel.fromPath(params.cellbender.input_mat)
    data_pkl = channel.fromPath(params.QC.data_pkl)
    metadata = channel.fromPath(params.metadata)
    adata = channel.fromPath(params.adata_QC)

    //cellbender_results = cellbender(input_mat)

    //QC(data_pkl, metadata, params.sample_key)
    //plot_QC(QC.out.adata_QC, params.sample_key, params.plot_QC.color_by)
    integration(adata, params.filename, params.integration.n_genes, params.integration.label_key, params.integration.batch_key)

    }


// Make 2 different workflows

