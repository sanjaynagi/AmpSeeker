from snakemake.utils import validate
import pandas as pd

# Split into two sample sets as bcftools merge cant take over 1000 files
# So we must do two rounds of merging
if len(metadata) > 1000:
    large_sample_size = True
    n_samples = len(metadata)
    half = int(n_samples/2)
    samples1 = metadata['sampleID'][:half]
    samples2 = metadata['sampleID'][half:]
else:
    large_sample_size = False
    samples1 = []
    samples2 = []

rule set_kernel:
    input:
        f'{workflow.basedir}/envs/AmpSeeker-python.yaml'
    output:
        touch("results/.kernel.set")
    conda: f'{workflow.basedir}/envs/AmpSeeker-python.yaml'
    log:
        "logs/set_kernel.log"
    shell: 
        """
        python -m ipykernel install --user --name=AmpSeq_python 2> {log}
        """

rule process_notebooks:
    input:
        input_nb = f"{workflow.basedir}/notebooks/process-notebooks.ipynb",
        snp_df = "docs/ampseeker-results/notebooks/snp-dataframe.ipynb",
        igv = "docs/ampseeker-results/notebooks/IGV-explore.ipynb" if config["analysis"]["igv"] else [],
        coverage = "docs/ampseeker-results/notebooks/coverage.ipynb" if config["quality-control"]["coverage"] else [],
        pca = "docs/ampseeker-results/notebooks/principal-component-analysis.ipynb" if config["analysis"]["pca"] else [],
        af = "docs/ampseeker-results/notebooks/allele-frequencies.ipynb" if config["analysis"]["allele-frequencies"] else [],
        sample_map = "docs/ampseeker-results/notebooks/sample-map.ipynb" if config["analysis"]["sample-map"] else [],
        read_quality = "docs/ampseeker-results/notebooks/read-quality.ipynb" if config["quality-control"]["fastp"] else [],
        reads_per_well = "docs/ampseeker-results/notebooks/reads-per-well.ipynb" if plate_info else [],
    output:
        out_nb = "results/notebooks/process-notebooks.ipynb",
    conda:
        f'{workflow.basedir}/envs/pythonGenomics.yaml'
    log:
        "logs/notebooks/process_notebooks.log"
    params:
        wkdir = wkdir
    shell:
        """
        papermill -k AmpSeq_python {input.input_nb} {output.out_nb} -p wkdir {params.wkdir} 2> {log}
        """


def AmpSeekerOutputs(wildcards):
    inputs = []
    
    inputs.extend(
        expand(
            [
                "results/alignments/{sample}.bam",
                "results/alignments/{sample}.bam.bai",
            ],                sample=samples,
                ),
    )
    
            
    if config['bcl-convert']:
        inputs.extend(["results/qc/index-read-qc/I1.html", "results/qc/index-read-qc/I2.html"])

    if plate_info:
        inputs.extend(["results/notebooks/reads-per-well.ipynb", 
                        "docs/ampseeker-results/notebooks/reads-per-well.ipynb"])
 
    if large_sample_size:
        inputs.extend(
            expand(
                [
                    "results/vcfs/{call_type}/{dataset}.complete.merge_vcfs",
                    "results/vcfs/{call_type}/{dataset}.annot.vcf"
                ],
                dataset=config['dataset'], 
                call_type=call_types
                )
            )
    else:
        inputs.extend(expand("results/vcfs/{call_type}/{dataset}.annot.vcf", dataset=config['dataset'], call_type=call_types))

    if config['quality-control']['coverage']:
            inputs.extend(
                expand(
                    [
                        "results/coverage/{sample}.per-base.bed.gz",
                        "results/notebooks/coverage.ipynb",
                        "docs/ampseeker-results/notebooks/coverage.ipynb"
                    ],
                sample=samples)
                )

    if config['quality-control']['fastp']:
        inputs.extend(
            expand(
                [
                    "results/qc/fastp_reports/{sample}.html",
                    "results/notebooks/read-quality.ipynb",
                    "docs/ampseeker-results/notebooks/read-quality.ipynb"
                ],
                sample=samples,
            )
        )


    if config['quality-control']['qualimap']:
        inputs.extend(               
            expand(
                [
                    "results/qc/qualimap/{sample}",
                ],
                sample=samples,
            )
        )
    
    if config['quality-control']['multiqc']:
        inputs.extend(
            expand(
                [
                    "results/qc/multiqc/multiqc_report.html",
                ],
            )
        )
    
    if config['quality-control']['stats']:
        inputs.extend(
            expand(
                [
                    "results/alignments/bamStats/{sample}.flagstat",
                    "results/qc/{dataset}.merged.vcf.txt",
                ],
                sample=samples, 
                dataset=config['dataset'],
            )
        )

    if config['analysis']['pca']['activate']:
        inputs.extend(
            expand(
                [
                    "results/notebooks/principal-component-analysis.ipynb",
                    "docs/ampseeker-results/notebooks/principal-component-analysis.ipynb"
                ],
            )
        )

    if config['analysis']['allele-frequencies']['activate']:
        inputs.extend(
            expand(
                [
                    "results/notebooks/allele-frequencies.ipynb",
                    "docs/ampseeker-results/notebooks/allele-frequencies.ipynb"
                ],
            )
        )

    if config['analysis']['sample-map']:
        inputs.extend(
            expand(
                [
                    "results/notebooks/sample-map.ipynb",
                    "docs/ampseeker-results/notebooks/sample-map.ipynb"
                ],
            )
        )

    if config['analysis']['igv']:
        inputs.extend(
            [
                "results/notebooks/IGV-explore.ipynb",
                "docs/ampseeker-results/notebooks/IGV-explore.ipynb"
            ]
            )

    if config['build-jupyter-book']:
        inputs.extend(["results/ampseeker-results/_build/html/index.html"])
    
    return inputs








def welcome(version):
    import datetime

    print("---------------------------- AmpSeeker ----------------------------")
    print(f"Running AmpSeeker snakemake workflow in {workflow.basedir}\n")
    print(f"Authors:   Sanjay Curtis Nagi, Trevor Mugoya, Edward Lukyamezi")
    print(f"Workflow Version: {version}")
    print("Execution time: ", datetime.datetime.now())
    print(f"Dataset: {config['dataset']}", "\n")
