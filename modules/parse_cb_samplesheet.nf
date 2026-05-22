process parse_cb_samplesheet { 
    label "samplesheet"

    debug true

    input:
        path samplesheet

    output:
        path "01_adata_concat_QC.h5ad", emit: samplesheet

    script:
        """
        parse_cb_samplesheet.py --samplesheet ${samplesheet}
        """
        }