rule species_id:
    input:
        nb = f"{workflow.basedir}/notebooks/ag-vampir/species-id.ipynb",
        kernel = "results/.kernel.set",
        vcf = expand("results/vcfs/targets/{dataset}.annot.vcf", dataset=dataset),
        bed = config['targets'],
        metadata = config["metadata"],
    output:
        nb = "results/notebooks/ag-vampir/species-id.ipynb",
        docs_nb = "docs/ampseeker-results/notebooks/ag-vampir/species-id.ipynb"
    conda:
        "../envs/AmpSeeker-python.yaml"
    log:
        "logs/notebooks/ag-vampir/species-id.log"
    params:
        dataset = dataset,
    shell:
        """
        papermill {input.nb} {output.nb} -k AmpSeq_python -p metadata_path {input.metadata} -p vcf_path {input.vcf} -p bed_targets_path {input.bed} 2> {log}
        cp {output.nb} {output.docs_nb} 2>> {log}
        """