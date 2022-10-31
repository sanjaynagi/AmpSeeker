from snakemake.utils import validate
import pandas as pd

# this container defines the underlying OS for each job when using the workflow
# with --use-conda --use-singularity
singularity: "docker://continuumio/miniconda3"

##### load config and conditional logic for splitting sample files if greater than 1000 #####

configfile: "config/config.yaml"
metadata = pd.read_csv(config['metadata'], sep="\t")
samples = metadata['sampleID']
sequence_data = config['sequence_data']
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

# Broken Validation Code
# validate(config, schema="../schemas/config.schema.yaml")
# samples = pd.read_csv(config["samples"], sep="\t").set_index("sample", drop=False)
# samples.index.names = ["sample_id"]
# validate(samples, schema="../schemas/samples.schema.yaml")

def RefBed(sequence_data):
    if sequence_data == "amplicon":
        ref=config['ref']['amplicon']
        bed=config['bed']['amplicon']
        return ref, bed
    if sequence_data == "wholegenome":
        ref=config['ref']['wholegenome']
        bed=config['bed']['wholegenome']
        return ref, bed

REF, BED = RefBed(sequence_data)

def coverage(wildcards):
    if sequence_data == "amplicon":
        return expand("results/coverage/{sample}.per-base.bed.gz", sample=samples)
    if sequence_data == "wholegenome":
        return expand("results/wholegenome/coverage/windowed/{sample}.regions.bed.gz", sample=samples)

def stats(wildcards):
    return expand("results/alignments/bamStats/{sample}.flagstat", sample=samples,ref=REF)

def vcf(wildcards):
    return expand("results/vcfs/{dataset}_LSTM_merged.vcf", ref=REF, dataset=config['dataset'])

def merge_vcfs(metadata):
    return "results/vcfs/.complete.merge_vcfs" if large_sample_size else []
