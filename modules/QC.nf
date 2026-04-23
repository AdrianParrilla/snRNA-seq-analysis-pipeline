process QC { 
    label "QC"

    input:
        path data_pkl
        path metadata 
        val sample_key

    output:
        path "01_adata_concat_QC.h5ad", emit: adata_QC

    script:
        """
        01_preprocessing.py --data_pkl ${data_pkl} --metadata_path ${metadata} --sample_key ${sample_key}
        """
        }