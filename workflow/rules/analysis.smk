
rule igv_notebook:
    input:
        nb=f"{workflow.basedir}/notebooks/IGV-explore.ipynb",
        kernel="results/.kernel.set",
        alignments=expand("results/alignments/{sample}.bam", sample=samples),
        genome=config["reference-fasta"],
        index=config["reference-fasta"] + ".fai",
        gff3=config["reference-gff3"],
        metadata=config["metadata"],
    output:
        nb="results/notebooks/IGV-explore.ipynb",
        docs_nb="docs/ampseeker-results/notebooks/IGV-explore.ipynb",
    conda:
        "../envs/AmpSeeker-python.yaml"
    log:
        "logs/notebooks/IGV-explore.log",
    params:
        reference_name=config["reference-name"],
    shell:
        """
        papermill {input.nb} {output.nb} -k AmpSeq_python -p metadata_path {input.metadata} -p genome_name {params.reference_name} -p reference_fasta {input.genome} -p reference_gff3 {input.gff3} 2> {log}
        cp {output.nb} {output.docs_nb} 2>> {log}
        """


rule pca:
    input:
        nb=f"{workflow.basedir}/notebooks/principal-component-analysis.ipynb",
        kernel="results/.kernel.set",
        vcf=expand("results/vcfs/targets/{dataset}.annot.vcf", dataset=dataset),
        metadata="results/config/metadata.qcpass.tsv",
    output:
        nb="results/notebooks/principal-component-analysis.ipynb",
        docs_nb="docs/ampseeker-results/notebooks/principal-component-analysis.ipynb",
    conda:
        "../envs/AmpSeeker-python.yaml"
    log:
        "logs/notebooks/principal-component-analysis.log",
    params:
        dataset=dataset,
        column=config["analysis"]["pca"]["colour-column"],
    shell:
        """
        papermill {input.nb} {output.nb} -k AmpSeq_python -p metadata_path {input.metadata} -p dataset {params.dataset} -p vcf_path {input.vcf} -p colour_column {params.column} 2> {log}
        cp {output.nb} {output.docs_nb} 2>> {log}
        """


rule sample_map:
    input:
        nb=f"{workflow.basedir}/notebooks/sample-map.ipynb",
        kernel="results/.kernel.set",
        metadata="results/config/metadata.qcpass.tsv",
    output:
        nb="results/notebooks/sample-map.ipynb",
        docs_nb="docs/ampseeker-results/notebooks/sample-map.ipynb",
    conda:
        "../envs/AmpSeeker-python.yaml"
    log:
        "logs/notebooks/sample-map.log",
    params:
        dataset=dataset,
    shell:
        """
        papermill {input.nb} {output.nb} -k AmpSeq_python -p metadata_path {input.metadata} -p dataset {params.dataset} 2> {log}
        cp {output.nb} {output.docs_nb} 2>> {log}
        """


rule allele_frequencies:
    input:
        nb=f"{workflow.basedir}/notebooks/allele-frequencies.ipynb",
        kernel="results/.kernel.set",
        metadata="results/config/metadata.qcpass.tsv",
        bed=config["targets"],
        vcf=expand("results/vcfs/targets/{dataset}.annot.vcf", dataset=dataset),
    output:
        nb="results/notebooks/allele-frequencies.ipynb",
        docs_nb="docs/ampseeker-results/notebooks/allele-frequencies.ipynb",
    conda:
        "../envs/AmpSeeker-python.yaml"
    log:
        "logs/notebooks/allele-frequencies.log",
    params:
        dataset=dataset,
        cohort_column=config["analysis"]["allele-frequencies"]["cohort-column"],
    shell:
        """
        papermill {input.nb} {output.nb} -k AmpSeq_python -p metadata_path {input.metadata} -p dataset {params.dataset} -p bed_path {input.bed} -p cohort_column {params.cohort_column} -p vcf_path {input.vcf} 2> {log}
        cp {output.nb} {output.docs_nb} 2>> {log}
        """


rule snp_dataframe:
    input:
        nb=f"{workflow.basedir}/notebooks/snp-dataframe.ipynb",
        kernel="results/.kernel.set",
        vcf=expand(
            "results/vcfs/{call_type}/{dataset}.annot.vcf",
            dataset=dataset,
            call_type=call_types,
        ),
    output:
        nb="results/notebooks/snp-dataframe.ipynb",
        docs_nb="docs/ampseeker-results/notebooks/snp-dataframe.ipynb",
    conda:
        "../envs/AmpSeeker-python.yaml"
    log:
        "logs/notebooks/snp-dataframe.log",
    params:
        dataset=dataset,
        wkdir=wkdir,
    shell:
        """
        papermill {input.nb} {output.nb} -k AmpSeq_python -p dataset {params.dataset} -p wkdir {params.wkdir} 2> {log}
        cp {output.nb} {output.docs_nb} 2>> {log}
        """

rule genetic_diversity:
    input:
        nb=f"{workflow.basedir}/notebooks/genetic-diversity.ipynb",
        kernel="results/.kernel.set",
        vcf=expand("results/vcfs/amplicons/{dataset}.annot.vcf", dataset=dataset),
        metadata="results/config/metadata.qcpass.tsv",
    output:
        nb="results/notebooks/genetic-diversity.ipynb",
        docs_nb="docs/ampseeker-results/notebooks/genetic-diversity.ipynb",
    conda:
        "../envs/AmpSeeker-python.yaml"
    log:
        "logs/notebooks/genetic-diversity.log",
    params:
        dataset=dataset,
        wkdir=wkdir,
        cohort_column=config["analysis"]["genetic-diversity"]["cohort-column"],
    shell:
        """
        papermill {input.nb} {output.nb} -k AmpSeq_python -p dataset {params.dataset} -p vcf_path {input.vcf} -p wkdir {params.wkdir} -p cohort_column {params.cohort_column} -p metadata_path {input.metadata} 2> {log}
        cp {output.nb} {output.docs_nb} 2>> {log}
        """
