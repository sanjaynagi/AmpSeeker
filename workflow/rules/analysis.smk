rule population_structure:
    input:
        nb=f"{workflow.basedir}/notebooks/population-structure.ipynb",
        kernel="results/.kernel.set",
        vcf=expand("results/vcfs/amplicons/{dataset}.annot.vcf", dataset=dataset),
        metadata="results/config/metadata.qcpass.tsv",
        taxon_complete = "results/.taxon.complete" if panel == "ag-vampir" else [],
        taxon_colours = "results/.taxon.colour.complete" if panel == "ag-vampir" else [],
    output:
        nb="results/notebooks/population-structure.ipynb",
        docs_nb="docs/ampseeker-results/notebooks/population-structure.ipynb",
    conda:
        "../envs/AmpSeeker-python.yaml"
    log:
        "logs/notebooks/population-structure.log",
    params:
        dataset=dataset,
        cohort_cols=cohort_cols,
        wkdir=wkdir,
        platform = config['platform']
    shell:
        """
        papermill {input.nb} {output.nb} \
            -k AmpSeq_python \
            -p wkdir {params.wkdir} \
            -p platform {params.platform} \
            -p metadata_path {input.metadata} \
            -p dataset {params.dataset} \
            -p vcf_path {input.vcf} \
            -p cohort_cols {params.cohort_cols} 2> {log}
        
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
        papermill {input.nb} {output.nb} \
            -k AmpSeq_python \
            -p metadata_path {input.metadata} \
            -p dataset {params.dataset} 2> {log}
        
        cp {output.nb} {output.docs_nb} 2>> {log}
        """


rule allele_frequencies:
    input:
        nb=f"{workflow.basedir}/notebooks/allele-frequencies.ipynb",
        kernel="results/.kernel.set",
        metadata="results/config/metadata.qcpass.tsv",
        bed=config["targets"],
        vcf=expand("results/vcfs/targets/{dataset}.annot.vcf", dataset=dataset),
        whole_amp_vcf=expand("results/vcfs/amplicons/{dataset}.annot.vcf", dataset=dataset),
        taxon_complete = "results/.taxon.complete" if panel == "ag-vampir" else [],
        taxon_colours = "results/.taxon.colour.complete" if panel == "ag-vampir" else [],
    output:
        nb="results/notebooks/allele-frequencies.ipynb",
        docs_nb="docs/ampseeker-results/notebooks/allele-frequencies.ipynb",
    conda:
        "../envs/AmpSeeker-python.yaml"
    log:
        "logs/notebooks/allele-frequencies.log",
    params:
        dataset=dataset,
        cohort_cols=cohort_cols,
        wkdir = wkdir,
        platform = config['platform']
    shell:
        """
        papermill {input.nb} {output.nb} \
            -k AmpSeq_python \
            -p metadata_path {input.metadata} \
            -p platform {params.platform} \
            -p dataset {params.dataset} \
            -p wkdir {params.wkdir} \
            -p bed_path {input.bed} \
            -p cohort_cols {params.cohort_cols} \
            -p vcf_path {input.vcf} 2> {log}
        
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
        platform = config['platform']
    shell:
        """
        papermill {input.nb} {output.nb} \
            -k AmpSeq_python \
            -p dataset {params.dataset} \
            -p platform {params.platform} \
            -p wkdir {params.wkdir} 2> {log}
        
        cp {output.nb} {output.docs_nb} 2>> {log}
        """

rule genetic_diversity:
    input:
        nb=f"{workflow.basedir}/notebooks/genetic-diversity.ipynb",
        kernel="results/.kernel.set",
        vcf=expand("results/vcfs/amplicons/{dataset}.annot.vcf", dataset=dataset),
        metadata="results/config/metadata.qcpass.tsv",
        taxon_complete = "results/.taxon.complete" if panel == "ag-vampir" else [],
        taxon_colours = "results/.taxon.colour.complete" if panel == "ag-vampir" else [],
    output:
        nb="results/notebooks/genetic-diversity.ipynb",
        docs_nb="docs/ampseeker-results/notebooks/genetic-diversity.ipynb",
        table = expand("results/genetic-diversity/{coh}.pi.tsv", coh=cohort_cols.split(",")),
        samples_table = "results/genetic-diversity/samples.pi.tsv"
    conda:
        "../envs/AmpSeeker-python.yaml"
    log:
        "logs/notebooks/genetic-diversity.log",
    params:
        dataset=dataset,
        wkdir=wkdir,
        cohort_cols=cohort_cols,
        platform = config['platform']
    shell:
        """
        papermill {input.nb} {output.nb} \
            -k AmpSeq_python \
            -p dataset {params.dataset} \
            -p platform {params.platform} \
            -p vcf_path {input.vcf} \
            -p wkdir {params.wkdir} \
            -p cohort_cols {params.cohort_cols} \
            -p metadata_path {input.metadata} 2> {log}
        
        cp {output.nb} {output.docs_nb} 2>> {log}
        """











# rule igv_notebook:
#     input:
#         nb=f"{workflow.basedir}/notebooks/IGV-explore.ipynb",
#         kernel="results/.kernel.set",
#         alignments=expand("results/alignments/{sample}.bam", sample=samples),
#         genome=config["reference-fasta"],
#         index=config["reference-fasta"] + ".fai",
#         gff3=config["reference-gff3"],
#         metadata="results/config/metadata.tsv",
#     output:
#         nb="results/notebooks/IGV-explore.ipynb",
#         docs_nb="docs/ampseeker-results/notebooks/IGV-explore.ipynb",
#     conda:
#         "../envs/AmpSeeker-python.yaml"
#     log:
#         "logs/notebooks/IGV-explore.log",
#     params:
#         reference_name=config["reference-name"],
#         wkdir = wkdir
#     shell:
#         """
#         papermill {input.nb} {output.nb} -k AmpSeq_python -p metadata_path {input.metadata} -p wkdir {params.wkdir} -p genome_name {params.reference_name} -p reference_fasta {input.genome} -p reference_gff3 {input.gff3} 2> {log}
#         cp {output.nb} {output.docs_nb} 2>> {log}
#         """
