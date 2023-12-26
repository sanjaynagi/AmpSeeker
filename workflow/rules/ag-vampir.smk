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


rule process_toc:
    input:
        input_nb = f"{workflow.basedir}/notebooks/ag-vampir/process-toc.ipynb",
        toc = "docs/ampseeker-results/_toc.yml",
    output:
        out_nb = "results/notebooks/ag-vampir/process-toc.ipynb",
        toc_complete = touch("results/.toc.ag-vampir.complete"),
    conda:
        f'{workflow.basedir}/envs/AmpSeeker-python.yaml'
    log:
        "logs/notebooks/ag-vampir/process-toc.log"
    params:
        wkdir = wkdir
    shell:
        """
        papermill -k AmpSeq_python {input.input_nb} {output.out_nb} -p wkdir {params.wkdir} 2> {log}
        """
