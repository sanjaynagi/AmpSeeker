from snakemake.utils import validate

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

def AmpSeekerOutputs(wildcards):
    inputs = []
    if config['Mapping']['activate']:
        if large_sample_size:
            inputs.extend("results/vcfs/{dataset}.complete.merge_vcfs" if large_sample_size else [])
        else:
            inputs.extend(expand("results/vcfs/{dataset}.merged.vcf", ref=REF, regions=BED, dataset=config['dataset']))

    if config['Stats']['activate']:
        inputs.extend(expand("results/alignments/bamStats/{sample}.flagstat", sample=samples,ref=REF))

    if config['Coverage']['activate']:
        if sequence_data == "amplicon":
            inputs.extend(expand("results/coverage/{sample}.per-base.bed.gz", sample=samples))
        if sequence_data == "wholegenome":
            inputs.extend(expand("results/wholegenome/coverage/windowed/{sample}.regions.bed.gz", sample=samples))
    
    return inputs
