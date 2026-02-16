import pandas as pd
import numpy as np
import allel
import re

def load_bed(bed_path: str, expected_cols = ['contig', 'start', 'end', 'amplicon_id', 'mutation', 'ref', 'alt']) -> pd.DataFrame:
    """
    Reads a BED-like file with optional REF and ALT columns.
    Returns a dataframe with consistent column naming.
    """
    # Read the file without predefined column names
    df = pd.read_csv(bed_path, sep="\t", header=None)
    
    # Assign only the first N column names that exist
    df.columns = expected_cols[: df.shape[1]]
    
    # Add any missing columns as NaN
    for col in expected_cols[df.shape[1] :]:
        df[col] = pd.NA
    
    return df[expected_cols]

def load_vcf(vcf_path, metadata, platform, filter_missing=None):
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

    ref = vcf['variants/REF']
    alt = vcf['variants/ALT']
    ann = read_ANN_field(vcf_path)

    if platform == "illumina": # remove any indels 
        indel = vcf['variants/INDEL']
    elif platform == "nanopore": # remove any indels > 1bp
        indel = ~vcf['variants/is_snp']

    geno = geno.compress(~indel, axis=0)
    pos = pos[~indel]
    contig = contig[~indel]
    ref = ref[~indel]
    alt = alt[~indel]
    ann = ann[~indel]

    if filter_missing:
        missing_mask = geno.is_missing().sum(axis=1) > geno.shape[1] * filter_missing
        geno = geno.compress(~missing_mask, axis=0)
        pos = pos[~missing_mask]
        contig = contig[~missing_mask]
        ref = ref[~missing_mask]
        alt = alt[~missing_mask]
        ann = ann[~missing_mask]    

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

    # remove invariant sites
    loc_var = np.any(gn_alt != gn_alt[:, 0, np.newaxis], axis=1)
    gn_var = np.compress(loc_var, gn_alt, axis=0)
    
    coords, model = allel.pca(gn_var, n_components=n_components)

    # flip axes back so PC1 is same orientation in each window 
    for i in range(n_components):
        c = coords[:, i]
    if np.abs(c.min()) > np.abs(c.max()):
        coords[:, i] = c * -1
    
    pca_df = pd.DataFrame(coords)
    pca_df.columns = [f"PC{pc+1}" for pc in range(n_components)]
    pca_df = pd.concat([metadata.reset_index(drop=True), pca_df], axis=1)
    
    return pca_df, model




import numba
@numba.njit(parallel=True)
def multiallelic_diplotype_pdist(X, metric):
    """Optimised implementation of pairwise distance between diplotypes.

    N.B., here we assume the array X provides diplotypes as genotype allele
    counts, with axes in the order (n_samples, n_sites, n_alleles).

    Computation will be faster if X is a contiguous (C order) array.

    The metric argument is the function to compute distance for a pair of
    diplotypes. This can be a numba jitted function.

    """
    n_samples = X.shape[0]
    n_pairs = (n_samples * (n_samples - 1)) // 2
    out = np.zeros(n_pairs, dtype=np.float32)

    # Loop over samples, first in pair.
    for i in range(n_samples):
        x = X[i, :, :]

        # Loop over observations again, second in pair.
        for j in numba.prange(i + 1, n_samples):
            y = X[j, :, :]

            # Compute distance for the current pair.
            d = metric(x, y)

            # Store result for the current pair.
            k = square_to_condensed(i, j, n_samples)
            out[k] = d

    return out


@numba.njit
def square_to_condensed(i, j, n):
    """Convert distance matrix coordinates from square form (i, j) to condensed form."""

    assert i != j, "no diagonal elements in condensed matrix"
    if i < j:
        i, j = j, i
    return n * j - j * (j + 1) // 2 + i - 1 - j


@numba.njit
def multiallelic_diplotype_mean_cityblock(x, y):
    """Compute the mean cityblock distance between two diplotypes x and y. The
    diplotype vectors are expected as genotype allele counts, i.e., x and y
    should have the same shape (n_sites, n_alleles).

    N.B., here we compute the mean value of the distance over sites where
    both individuals have a called genotype. This avoids computing distance
    at missing sites.

    """
    n_sites = x.shape[0]
    n_alleles = x.shape[1]
    distance = np.float32(0)
    n_sites_called = np.float32(0)

    # Loop over sites.
    for i in range(n_sites):
        x_is_called = False
        y_is_called = False
        d = np.float32(0)

        # Loop over alleles.
        for j in range(n_alleles):
            # Access allele counts.
            xc = np.float32(x[i, j])
            yc = np.float32(y[i, j])

            # Check if any alleles observed.
            x_is_called = x_is_called or (xc > 0)
            y_is_called = y_is_called or (yc > 0)

            # Compute cityblock distance (absolute difference).
            d += np.fabs(xc - yc)

        # Accumulate distance for the current pair, but only if both samples
        # have a called genotype.
        if x_is_called and y_is_called:
            distance += d
            n_sites_called += np.float32(1)

    # Compute the mean distance over sites with called genotypes.
    if n_sites_called > 0:
        mean_distance = distance / n_sites_called
    else:
        mean_distance = np.nan

    return mean_distance
