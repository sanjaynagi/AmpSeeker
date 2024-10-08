rule species_id:
    input:
        nb=f"{workflow.basedir}/notebooks/ag-vampir/species-id.ipynb",
        kernel="results/.kernel.set",
        vcf=expand("results/vcfs/targets/{dataset}.annot.vcf", dataset=dataset),
        bed=config["targets"],
        metadata="results/config/metadata.qcpass.tsv",
    output:
        nb="results/notebooks/ag-vampir/species-id.ipynb",
        docs_nb="docs/ampseeker-results/notebooks/ag-vampir/species-id.ipynb",
        aims = "results/ag-vampir/aims/taxon_aims.tsv",
        taxon_complete=touch("results/.taxon.complete"),
    conda:
        "../envs/AmpSeeker-python.yaml"
    log:
        "logs/notebooks/ag-vampir/species-id.log",
    params:
        dataset=dataset,
        cohort_cols=cohort_cols,
        wkdir=wkdir,
    shell:
        """
        papermill {input.nb} {output.nb} -k AmpSeq_python -p metadata_path {input.metadata} -p wkdir {params.wkdir} -p vcf_path {input.vcf} -p bed_targets_path {input.bed} -p cohort_cols {params.cohort_cols} 2> {log}
        cp {output.nb} {output.docs_nb} 2>> {log}
        """


rule kdr_origin:
    input:
        nb=f"{workflow.basedir}/notebooks/ag-vampir/kdr-origins.ipynb",
        kernel="results/.kernel.set",
        vcf=expand("results/vcfs/targets/{dataset}.annot.vcf", dataset=dataset),
        metadata="results/config/metadata.qcpass.tsv",
        kdr_origin_SNPs="resources/ag-vampir/Kdr_marker_SNPs.csv",
        taxon_complete="results/.taxon.complete",
    output:
        nb="results/notebooks/ag-vampir/kdr-origins.ipynb",
        docs_nb="docs/ampseeker-results/notebooks/ag-vampir/kdr-origins.ipynb",
        kdr_origins="results/kdr-origins/kdr_origins.csv",
        kdr_genhap_origins="results/kdr-origins/kdr_genhap_origins.csv",
    conda:
        "../envs/AmpSeeker-python.yaml"
    log:
        "logs/notebooks/ag-vampir/kdr-origins.log",
    params:
        dataset=dataset,
        cohort_cols=cohort_cols,
        wkdir=wkdir,
    shell:
        """
        papermill {input.nb} {output.nb} -k AmpSeq_python -p metadata_path {input.metadata} -p vcf_path {input.vcf} -p cohort_cols {params.cohort_cols} -p wkdir {params.wkdir} -p kdr_marker_snps_path {input.kdr_origin_SNPs} 2> {log}
        cp {output.nb} {output.docs_nb} 2>> {log}
        """


