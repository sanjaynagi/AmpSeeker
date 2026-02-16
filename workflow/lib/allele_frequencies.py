import numpy as np
import pandas as pd
import plotly.express as px

import warnings
warnings.filterwarnings('ignore')

import shared as amp

def vcf_to_snp_dataframe(vcf_path, metadata, platform, filter_missing=None):
    geno, pos, contig, metadata, ref, alt, ann = amp.load_variants(
        vcf_path=vcf_path,
        metadata=metadata,
        platform=platform,
        filter_missing=filter_missing,
        filter_indel=True,
    )
    
    # make dataframe of variant positions and merge with bed
    snp_df = pd.DataFrame({'contig':contig, 'pos':pos, 'ref':ref, 'alt':[list(a[a != ""]) for a in alt], 'ann':ann})

    snp_df = snp_df.explode('alt').reset_index().rename(columns={'index':'variant_index'})
    snp_df = snp_df.assign(alt_index=snp_df.groupby(['contig','pos']).cumcount() + 1) 
    snp_df = snp_df.assign(label=lambda x: x.pos.astype(str) + " | " +  x.alt.fillna('NA'))
    snp_df.head(2)

    # split and find correct annotation 
    df = snp_df.assign(ann=lambda x: x.ann.str.split(","))
    anns = []
    for i, row in df.iterrows():
        alt = row['alt']
        if row['ann'] == None:
            ann = ""
        else:
            # keep only RD Vgsc annotations
            if 'AGAP004707' in ','.join(row['ann']):
                row['ann'] = [a for a in row['ann'] if "AGAP004707-RD" in a]

            ann = ','.join([a for a in row['ann'] if a.startswith(alt)])
        anns.append(ann)

    snp_df = snp_df.assign(ann=anns)
    
    return snp_df, geno

def calculate_frequencies_cohort(snp_df, metadata, geno, cohort_col, af_filter, missense_filter):
    np.seterr(all="ignore")
    
    df = snp_df.copy()
    
    # get indices of each cohort
    coh_dict = {}
    cohs = metadata[cohort_col].unique()
    cohs = cohs[~pd.isnull(cohs)]
    for coh in cohs:
        coh_dict[coh] = np.where(metadata[cohort_col] == coh)[0]
    
    # get allele counts for each population
    ac = geno.count_alleles_subpops(coh_dict, max_allele=3)
    
    for coh in cohs:
        total_counts = []
        alt_counts = []
        for i, row in df.iterrows():
            var_idx = row['variant_index']
            alt_idx = row['alt_index']
            total_counts.append(ac[coh][var_idx,:].sum())
            alt_counts.append(ac[coh][var_idx, alt_idx])

        df.loc[:, f'count_{coh}'] = np.array(alt_counts)
        df.loc[:, f'frq_{coh}'] = np.round(np.array(alt_counts)/np.array(total_counts), 3)
    
    freq_df = df.set_index('label').filter(like='frq')
    
    ann_df = snp_df.ann.str.split("|", expand=True).iloc[:, :11].drop(columns=[0,7,8])
    ann_df.columns = ['type', 'effect', 'gene', 'geneID', 'modifier', 'transcript', 'base_change', 'aa_change']
    snp_df = pd.concat([snp_df[['contig', 'pos', 'ref', 'alt']], ann_df], axis=1)
    snp_freq_df = pd.concat([snp_df, freq_df.reset_index()], axis=1)

    snp_freq_df = snp_freq_df.assign(label=
                  lambda x: x.contig + " | " + x.gene + " | " + x.pos.astype(str) + " | " + x.aa_change.str.replace("p.", "") + " | " + x.alt.fillna(" ")
                 )
    
    if af_filter:
        af_pass = (snp_freq_df.filter(like='frq') > 0.05).any(axis=1)
        snp_freq_df = snp_freq_df[af_pass]
    
    if missense_filter:
        snp_freq_df = snp_freq_df.query("type == 'missense_variant'")
    
    return snp_freq_df.set_index('label')

def plot_allele_frequencies(df, cohort_col, colscale="Reds"):
        
    fig = px.imshow(
            img=df,
            zmin=0,
            zmax=1,
            width=np.max([800, df.shape[1] * 100]),
            height=200 + (df.shape[0] * 18),
            text_auto=True,
            aspect=1,
            color_continuous_scale=colscale,
            title=f"Allele frequencies | by {cohort_col}",
        template='simple_white'
        )
    fig.update(layout_coloraxis_showscale=False)

    return fig
