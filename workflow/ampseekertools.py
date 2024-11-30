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

    ref = vcf['variants/REF'][~indel]
    alt = vcf['variants/ALT'][~indel]
    ann = read_ANN_field(vcf_path)[~indel]

    metadata = metadata.set_index('sample_id')
    samples = samples[sample_mask]

    return geno, pos, contig, metadata.loc[samples, :].reset_index() , ref, alt, ann


def read_ANN_field(vcf_file):
    anns = []
    with open(vcf_file, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue  # Skip header lines
            fields = line.strip().split('\t')
            info_field = fields[7]
            info_pairs = info_field.split(';')
            ann_value = None
            for pair in info_pairs:
                if pair.startswith('ANN='):
                    ann_value = pair.split('=')[1]
                    break
            anns.append(ann_value)

    return np.array(anns)


def pca(geno, metadata, n_components = 3, query=None, missing_threshold=0.05):
    """
    Load genotype data and run PCA 
    """

    if query:
        mask = metadata.eval(query)
        metadata = metadata[mask]
        geno = geno.compress(mask, axis=1)

    # find segregating sites
    print("removing any invariant and highly missing sites")
    ac = geno.count_alleles()
    seg = ac.is_segregating()
    gn_seg = geno.compress(seg, axis=0)

    # remove highly missing sites
    missing_mask = gn_seg.is_missing().sum(axis=1) > gn_seg.shape[1] * missing_threshold
    gn_seg = gn_seg.compress(~missing_mask, axis=0)
    gn_alt = gn_seg.to_n_alt()
    
    coords, model = allel.pca(gn_alt, n_components=n_components)
    # flip axes back so PC1 is same orientation in each window 
    for i in range(n_components):
        c = coords[:, i]
    if np.abs(c.min()) > np.abs(c.max()):
        coords[:, i] = c * -1
    
    pca_df = pd.DataFrame(coords)
    pca_df.columns = [f"PC{pc+1}" for pc in range(n_components)]
    pca_df = pd.concat([metadata.reset_index(drop=True), pca_df], axis=1)
    
    return pca_df, model
