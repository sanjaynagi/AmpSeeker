import pandas as pd
import numpy as np
import allel


def load_vcf(vcf_path, metadata):
    """
    Load VCF and filter poor-quality samples
    """
        
    # load vcf and get genotypes and positions
    vcf = allel.read_vcf(vcf_path, fields="*")
    samples = vcf['samples']    # keep only samples in qcpass metadata 
    sample_mask = np.isin(vcf['samples'], metadata.sample_id)
    
    # remove low quality samples 
    geno = allel.GenotypeArray(vcf['calldata/GT'])
    geno = geno.compress(sample_mask, axis=1)
    pos = vcf['variants/POS']
    contig = vcf['variants/CHROM']
    indel = vcf['variants/INDEL']
    
    # remove any indels 
    geno = geno.compress(~indel, axis=0)
    pos = pos[~indel]
    contig = contig[~indel]
    
    metadata = metadata.set_index('sample_id')
    samples = samples[sample_mask]
    
    return geno, pos, contig, metadata.loc[samples, :].reset_index() ,vcf['variants/REF'][~indel], vcf['variants/ALT'][~indel]