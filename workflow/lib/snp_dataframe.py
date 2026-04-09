import numpy as np
import pandas as pd
import shared

import warnings
warnings.filterwarnings('ignore')


def _is_missing_genotype(genotype):
    """Return True if a VCF genotype string represents a missing call.

    Works for any ploidy: haploid ``.``, diploid ``./.``, tetraploid
    ``./././.`` etc. A genotype is treated as missing if every allele
    position is ``.``.
    """
    if genotype is None:
        return True
    alleles = str(genotype).replace('|', '/').split('/')
    return all(a == '.' for a in alleles)

def vcf_to_excel(vcf_path, excel_path, convert_genotypes=False, split_multiallelic=False):
    """Load a VCF and optionally export it to Excel.

    Parameters
    ----------
    vcf_path : str or path-like
        Path to the input VCF file.
    excel_path : str or path-like or None
        Output path for the Excel workbook. If `None`, no file is written.
    convert_genotypes : bool, default=False
        Whether to convert genotype strings to alternate-allele counts.
    split_multiallelic : bool, default=False
        Whether to expand multiallelic records into one row per alternate
        allele.

    Returns
    -------
    pandas.DataFrame
        DataFrame representation of the VCF, optionally transformed.
    """
    # Read VCF and create a dictionary
    vcf_df = vcf_to_df(vcf_path)
    samples = vcf_df.columns[6:]

    # Create a DataFrame for variants
    if split_multiallelic:
        vcf_df = split_rows_with_multiple_alleles(vcf_df, samples)

    # Convert genotypes to 0/1/2
    if convert_genotypes:
        vcf_df = convert_genotype_to_alt_allele_count(vcf_df, samples)
    
    if excel_path:
        # Write to Excel
        vcf_df.to_excel(excel_path, index=False)

    return vcf_df

def vcf_to_df(vcf_path, seg=True):
    """Convert a VCF file into a tabular dataframe.

    Parameters
    ----------
    vcf_path : str or path-like
        Path to the input VCF file.
    seg : bool, default=True
        Whether to retain only segregating sites.

    Returns
    -------
    pandas.DataFrame
        Variant-level metadata joined with per-sample genotype strings.
    """
    core = shared.load_vcf(vcf_path)
    contig = core["contig"]
    pos = core["pos"]
    filter_pass = core["filter_pass"]
    ref = core["ref"]
    alt = np.array([",".join([a for a in row if a != ""]) for row in core["alt"]], dtype=object)
    ann = core["ann"]
    geno = core["geno"]

    if seg:
        ac = geno.count_alleles(max_allele=3)
        mask_seg = ac.is_segregating()
        geno = geno.compress(mask_seg, axis=0)
        contig = contig[mask_seg]
        pos = pos[mask_seg]
        filter_pass = filter_pass[mask_seg]
        ref = ref[mask_seg]
        alt = alt[mask_seg]
        ann = ann[mask_seg]

    ### make pd dataframe version of vcf
    vcf_df = pd.DataFrame({'CHROM': contig, 'POS': pos, 'FILTER_PASS': filter_pass, 'REF': ref, 'ALT': alt, 'ANN': ann})
    # make pd dataframe version of genotypes
    geno_df = pd.DataFrame(geno.to_gt().astype(str), columns=core["samples"])
    vcf = pd.concat([vcf_df, geno_df], axis=1)
    return vcf

def split_rows_with_multiple_alleles(df, samples):
    """Expand multiallelic rows into one row per alternate allele.

    Parameters
    ----------
    df : pandas.DataFrame
        VCF-derived dataframe containing an `ALT` column and genotype columns.
    samples : sequence of str
        Names of the genotype columns to rewrite for each split allele.

    Returns
    -------
    pandas.DataFrame
        Dataframe with multiallelic variants expanded to separate rows.
    """
    # Create an empty list to store the new rows
    new_rows = []
    # Iterate through each row
    for index, row in df.iterrows():
        alt_alleles = row['ALT'].split(',')
        # Check if there are multiple alleles in the ALT field
        if len(alt_alleles) > 1:
            for allele_num, allele in enumerate(alt_alleles):
                # Create a new row for each allele
                new_row = row.copy()
                new_row['ALT'] = allele
                # Update genotype fields
                for col in samples:
                    genotype = row[col]
                    # Split the genotype and process it
                    if not _is_missing_genotype(genotype):
                        gt_alleles = str(genotype).replace('|', '/').split('/')
                        new_gt = ['0' if (gt != '.' and int(gt) != allele_num + 1 and gt != '0') else gt for gt in gt_alleles]
                        new_row[col] = '/'.join(new_gt)
                new_rows.append(new_row)
        else:
            new_rows.append(row)
    
    # Create a new DataFrame from the new rows
    new_df = pd.DataFrame(new_rows).reset_index(drop=True)
    return new_df

def convert_genotype_to_alt_allele_count(df, samples):
    """Replace genotype strings with alternate-allele counts.

    Parameters
    ----------
    df : pandas.DataFrame
        VCF-derived dataframe containing genotype strings.
    samples : sequence of str
        Names of genotype columns to convert.

    Returns
    -------
    pandas.DataFrame
        Copy of `df` with genotype columns converted to counts or `NaN` for
        missing calls.
    """
    nalt_df = df.copy()
    # Keep genotype columns as object so integer/NaN assignments are portable
    # across pandas versions (string extension dtype is stricter in 3.11 CI).
    nalt_df = nalt_df.astype({col: "object" for col in samples})
    # Iterate through each row
    for index, row in df.iterrows():
        # Update genotype fields
        for col in samples:
                genotype = row[col]
                if not _is_missing_genotype(genotype):
                    # Split the genotype and count alt alleles. This works
                    # for any ploidy: haploid "0"/"1", diploid "0/1",
                    # tetraploid "0/0/1/1" etc.
                    alleles = str(genotype).replace('|', '/').split('/')
                    alt_allele_count = sum(1 for allele in alleles if allele != '0' and allele != '.')
                    nalt_df.at[index, col] = alt_allele_count
                else:
                    nalt_df.at[index, col] = np.nan

    return nalt_df
