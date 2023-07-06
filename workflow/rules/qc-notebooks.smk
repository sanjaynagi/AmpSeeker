
rule reads_per_well:
    input:
        nb = f"{workflow.basedir}/notebooks/reads-per-well.ipynb",
        kernel = "results/.kernel.set",
        bam = expand("results/alignments/{sample}.bam", sample=samples),
        bai = expand("results/alignments/{sample}.bam.bai", sample=samples),
        stats = expand("results/alignments/bamStats/{sample}.flagstat", sample=samples),
        metadata = config["metadata"],
    output:
        nb = "results/notebooks/reads-per-well.ipynb",
        docs_nb = "docs/ampseeker-results/notebooks/reads-per-well.ipynb"
    conda:
        "../envs/AmpSeeker-python.yaml"
    log:
        "logs/notebooks/reads-per-well.log"
    shell:
        """
        papermill {input.nb} {output.nb} -k AmpSeq_python -p metadata_path {input.metadata} 2> {log}
        cp {output.nb} {output.docs_nb} 2>> {log}
        """

rule read_qc:
    input:
        nb = f"{workflow.basedir}/notebooks/read-quality.ipynb",
        kernel = "results/.kernel.set",
        fastp = expand("results/fastp_reports/{sample}.json", sample=samples),
        index_qc = rules.index_read_fastqc.output,
        metadata = config["metadata"],
    output:
        nb = "results/notebooks/read-quality.ipynb",
        docs_nb = "docs/ampseeker-results/notebooks/read-quality.ipynb"
    conda:
        "../envs/AmpSeeker-python.yaml"
    log:
        "logs/notebooks/reads-quality.log"
    params:
        index_qc = config['bcl-convert'],
        wd = wkdir
    shell:
        """
        papermill {input.nb} {output.nb} -k AmpSeq_python -p metadata_path {input.metadata} -p index_read_qc {params.index_qc} -p wkdir {params.wd} 2> {log}
        cp {output.nb} {output.docs_nb} 2>> {log}
        """
