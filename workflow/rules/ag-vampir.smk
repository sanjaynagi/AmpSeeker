rule species_id:
    input:
        nb=f"{workflow.basedir}/notebooks/ag-vampir/species-id.ipynb",
        kernel="results/.kernel.set",
        vcf=expand("results/vcfs/targets/{dataset}.annot.vcf", dataset=dataset),
        bed=config["targets"],
        qc_nb="results/notebooks/sample-quality-control.ipynb",
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
        platform='illumina'
    shell:
        """
        papermill {input.nb} {output.nb} \
            -k AmpSeq_python \
            -p platform {params.platform} \
            -p dataset {params.dataset} \
            -p metadata_path {input.metadata} \
            -p wkdir {params.wkdir} \
            -p vcf_path {input.vcf} \
            -p bed_targets_path {input.bed} \
            -p cohort_cols {params.cohort_cols} 2> {log}

        cp {output.nb} {output.docs_nb} 2>> {log}
        """


rule kdr_analysis:
    input:
        nb=f"{workflow.basedir}/notebooks/ag-vampir/kdr-analysis.ipynb",
        kernel="results/.kernel.set",
        vcf=expand("results/vcfs/targets/{dataset}.annot.vcf", dataset=dataset),
        vcf2=expand("results/vcfs/amplicons/{dataset}.annot.vcf", dataset=dataset),
        metadata="results/config/metadata.qcpass.tsv",
        kdr_origin_SNPs="resources/ag-vampir/Kdr_marker_SNPs.csv",
        taxon_complete="results/.taxon.complete",
        taxon_colours="results/.taxon.colour.complete",
    output:
        nb="results/notebooks/ag-vampir/kdr-analysis.ipynb",
        docs_nb="docs/ampseeker-results/notebooks/ag-vampir/kdr-analysis.ipynb",
        kdr_origins="results/ag-vampir/kdr-origins/kdr_origins.tsv",
        kdr_genhap_origins="results/ag-vampir/kdr-origins/kdr_genhap_origins.tsv",
        kdr_origins_counts="results/ag-vampir/kdr-origins/kdr_origin_counts.xlsx",
        kdr_origin_freqs="results/ag-vampir/kdr-origins/kdr_origin_freqs.xlsx"
    conda:
        "../envs/AmpSeeker-python.yaml"
    log:
        "logs/notebooks/ag-vampir/kdr-analysis.log",
    params:
        dataset=dataset,
        cohort_cols=cohort_cols,
        wkdir=wkdir,
        platform='illumina'
    shell:
        """
        papermill {input.nb} {output.nb} \
            -k AmpSeq_python \
            -p platform {params.platform} \
            -p dataset {params.dataset} \
            -p metadata_path {input.metadata} \
            -p vcf_path {input.vcf} \
            -p cohort_cols {params.cohort_cols} \
            -p wkdir {params.wkdir} \
            -p kdr_marker_snps_path {input.kdr_origin_SNPs} 2> {log}
        
        cp {output.nb} {output.docs_nb} 2>> {log}
        """


rule overwrite_metadata_colors:
    input:
        metadata = "results/config/metadata.qcpass.tsv",
        taxon="results/.taxon.complete"
    output:
        colors = "results/config/metadata_colours.json",
        taxon_colours=touch("results/.taxon.colour.complete")
    params:
        cohort_columns = config["cohort-columns"]
    run:
        import json
        import os
        import pandas as pd
        
        color_mappings = create_color_mapping(
            pd.read_csv(input.metadata, sep="\t"), 
            columns=params.cohort_columns
        )
        
        with open(output.colors, 'w') as f:
            json.dump(color_mappings, f, indent=2)