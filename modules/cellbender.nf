process cellbender { 
    label "cellbender"

    debug true

    input:
        path input_mat

    output:
        path "*_cellbender_filtered.h5", emit: cb_matrix
        path "*_cellbender_metrics.csv", emit: metrics
        path "*_cellbender.pdf"        , emit: report
        path "*_cellbender_report.html", emit: html_report, optional: true
        path "*_cellbender.log"        , emit: log

    script:
        def out_prefix = "${input_mat.baseName}_cellbender"
        def expected_cells = params.cellbender.expected_cells ? "--expected-cells ${params.cellbender.expected_cells}" : ""
        def total_droplets = params.cellbender.total_droplets_included ? "--total-droplets-included ${params.cellbender.total_droplets_included}" : ""
        """
        cellbender remove-background \\
            --cuda \\
            --input ${input_mat} \\
            --output ${out_prefix}.h5 \\
            --learning-rate ${params.cellbender.learning_rate} \\
            ${expected_cells} \\
            ${total_droplets}
        """
        } 