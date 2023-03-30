from snakemake.utils import validate
import pandas as pd

# Split into two sample sets as bcftools merge cant take over 1000 files
# So we must do two rounds of merging
large_sample_size = False
if len(metadata) > 1000:
    large_sample_size = True
    n_samples = len(metadata)
    half = int(n_samples/2)
    samples1 = metadata['sampleID'][:half]
    samples2 = metadata['sampleID'][half:]
else:
    samples1 = []
    samples2 = []


rule set_kernel:
    input:
        f'{workflow.basedir}/envs/AmpSeq_python.yaml'
    output:
        touch("results/.kernel.set")
    conda: f'{workflow.basedir}/envs/AmpSeq_python.yaml'
    log:
        "logs/set_kernel.log"
    shell: 
        """
        python -m ipykernel install --user --name=AmpSeq_python 2> {log}
        """

def AmpSeekerOutputs(wildcards):
    inputs = []
    if config['BCL_conversion']['activate']:
        inputs.extend(expand("results/reads/{sample}_{n}.fastq.gz", sample=samples, n=[1,2]))
        
    if config['Mapping']['activate']:
        if large_sample_size:
            inputs.extend(expand("results/vcfs/{dataset}.complete.merge_vcfs", dataset=config['dataset']) if large_sample_size else [])
        else:
            inputs.extend(expand("results/vcfs/{dataset}.merged.vcf", dataset=config['dataset']))

    if config['Stats']['activate']:
        inputs.extend(expand("results/alignments/bamStats/{sample}.flagstat", sample=samples))

    if config['Coverage']['activate']:
        if reference_type == "amplicon":
            inputs.extend(expand("results/coverage/{sample}.per-base.bed.gz", sample=samples))
        if reference_type == "wholegenome":
            inputs.extend(expand("results/wholegenome/coverage/windowed/{sample}.regions.bed.gz", sample=samples))

    inputs.extend(["results/notebooks/IGV-explore.ipynb"])
    
    return inputs

