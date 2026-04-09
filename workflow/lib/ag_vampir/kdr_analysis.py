import numpy as np
import pandas as pd
import plotly.express as px
import re

def _f_kdr_origin_gen(genotypes, clean = True):
    """Infer 995F kdr haplotype origin calls for one sample.

    Parameters
    ----------
    genotypes : pandas.Series
        Per-sample genotype calls containing the required F-origin markers.
    clean : bool, default=True
        Whether to post-process the raw origin string into a simplified call.

    Returns
    -------
    str
        Raw or cleaned 995F origin assignment for the sample.
    """
    if 'sample_name' in genotypes.index:
        sample_name = genotypes['sample_name']
    else:
        sample_name = genotypes.name
    # Check for the 995F mutations
    if pd.isnull(genotypes['kdr-995F']):
        kdr_F_origins = 'F:unknown'
    elif genotypes['kdr-995F'] == 'AA':
        kdr_F_origins = 'F:wt_hom'
    elif genotypes['kdr-995F'] == 'AT':
        kdr_F_origins = 'F:het'
    elif genotypes['kdr-995F'] == 'TT':
        kdr_F_origins = 'F:hom'
    else:
        print(f'Unexpected kdr F genotype. {sample_name} {genotypes["kdr-995F"]}')
        kdr_F_origins = 'Fail. Unexpected kdr F genotype'
    # If the individual has Fkdr, find out its origins
    # For F homozygotes
    if kdr_F_origins == 'F:hom':
        if pd.isnull(genotypes['Def-F1']):
            kdr_F_origins = f'{kdr_F_origins},F1?'
        elif genotypes['Def-F1'] == 'AA':
            kdr_F_origins = f'{kdr_F_origins},F1_hom'
        elif genotypes['Def-F1'] == 'AG':
            kdr_F_origins = f'{kdr_F_origins},F1_het'
        #
        if pd.isnull(genotypes['Def-F2']):
            kdr_F_origins = f'{kdr_F_origins},F2?'
        elif genotypes['Def-F2'] == 'AA':
            kdr_F_origins = f'{kdr_F_origins},F2_hom'
        elif genotypes['Def-F2'] == 'AG':
            kdr_F_origins = f'{kdr_F_origins},F2_het'
        #
        if pd.isnull(genotypes['Def-F3F4-2']):
            kdr_F_origins = f'{kdr_F_origins},F3F4?'
        elif genotypes['Def-F3F4-2'] == 'TT':
            if pd.isnull(genotypes['Def-F3']):
                kdr_F_origins = f'{kdr_F_origins},(F3F4)_hom'
            elif genotypes['Def-F3'] == 'CC':
                kdr_F_origins = f'{kdr_F_origins},F3_hom'
            elif genotypes['Def-F3'] == 'CG':
                kdr_F_origins = f'{kdr_F_origins},F3_het,F4_het'
            elif genotypes['Def-F3'] == 'GG':
                kdr_F_origins = f'{kdr_F_origins},F4_hom'
        elif genotypes['Def-F3F4-2'] == 'AT':
            if pd.isnull(genotypes['Def-F3']):
                kdr_F_origins = f'{kdr_F_origins},(F3F4)_het'
            elif genotypes['Def-F3'] == 'CC':
                kdr_F_origins = f'Fail. Genotypes suggest that sample {sample_name} is heterozygote for F3F4, but homozygote for F3.'
            elif genotypes['Def-F3'] == 'CG':
                kdr_F_origins = f'{kdr_F_origins},F3_het'
            elif genotypes['Def-F3'] == 'GG':
                kdr_F_origins = f'{kdr_F_origins},F4_het'
        #
        if pd.isnull(genotypes['Def-F5-2']):
            kdr_F_origins = f'{kdr_F_origins},F5?'
        elif genotypes['Def-F5-2'] == 'GG':
            kdr_F_origins = f'{kdr_F_origins},F5_hom'
        elif genotypes['Def-F5-2'] == 'AG':
            kdr_F_origins = f'{kdr_F_origins},F5_het'
    # for F heterozygotes
    elif kdr_F_origins == 'F:het':
        if pd.isnull(genotypes['Def-F1']):
            kdr_F_origins = f'{kdr_F_origins},F1?'
        elif genotypes['Def-F1'] == 'AA':
            kdr_F_origins = f'Fail. Genotypes suggest that sample {sample_name} is heterozygote for F kdr, but homozygote for F1.'
        elif genotypes['Def-F1'] == 'AG':
            kdr_F_origins = f'{kdr_F_origins},F1_het'
        #
        if pd.isnull(genotypes['Def-F2']):
            kdr_F_origins = f'{kdr_F_origins},F2?'
        elif genotypes['Def-F2'] == 'AA':
            kdr_F_origins = f'Fail. Genotypes suggest that sample {sample_name} is heterozygote for F kdr, but homozygote for F2.'
        elif genotypes['Def-F2'] == 'AG':
            kdr_F_origins = f'{kdr_F_origins},F2_het'
        #
        if pd.isnull(genotypes['Def-F3F4-2']):
            kdr_F_origins = f'{kdr_F_origins},F3F4?'
        elif genotypes['Def-F3F4-2'] == 'TT':
            kdr_F_origins = f'Fail. Genotypes suggest that sample {sample_name} is heterozygote for F kdr, but homozygote for F3F4.'
        elif genotypes['Def-F3F4-2'] == 'AT':
            if pd.isnull(genotypes['Def-F3']):
                kdr_F_origins = f'{kdr_F_origins},(F3F4)_het'
            elif genotypes['Def-F3'] == 'CC':
                kdr_F_origins = f'Fail. Genotypes suggest that sample {sample_name} is heterozygote for F kdr and F3F4, but homozygote for F3.'
            elif genotypes['Def-F3'] == 'CG':
                kdr_F_origins = f'{kdr_F_origins},F3_het'
            elif genotypes['Def-F3'] == 'GG':
                kdr_F_origins = f'{kdr_F_origins},F4_het'
        #
        if pd.isnull(genotypes['Def-F5-2']):
            kdr_F_origins = f'{kdr_F_origins},F5?'
        elif genotypes['Def-F5-2'] == 'GG':
            kdr_F_origins = f'Fail. Genotypes suggest that sample {sample_name} is heterozygote for F kdr, but homozygote for F5.'
        elif genotypes['Def-F5-2'] == 'AG':
            kdr_F_origins = f'{kdr_F_origins},F5_het'
    if clean:
        return(_kdr_gen_cleanup(kdr_F_origins))
    else:
        return(kdr_F_origins)

def _s_kdr_origin_gen(genotypes, clean = True, alternate_S4S5 = False):
    """Infer 995S kdr haplotype origin calls for one sample.

    Parameters
    ----------
    genotypes : pandas.Series
        Per-sample genotype calls containing the required S-origin markers.
    clean : bool, default=True
        Whether to post-process the raw origin string into a simplified call.
    alternate_S4S5 : bool, default=False
        Whether to use the alternate S4/S5 marker naming scheme.

    Returns
    -------
    str
        Raw or cleaned 995S origin assignment for the sample.
    """
    if 'sample_name' in genotypes.index:
        sample_name = genotypes['sample_name']
    else:
        sample_name = genotypes.name
    # Check for the 995S mutations
    if pd.isnull(genotypes['kdr-995S']):
        kdr_S_origins = 'S:unknown'
    elif genotypes['kdr-995S'] == 'TT':
        kdr_S_origins = 'S:wt_hom'
    elif genotypes['kdr-995S'] == 'CT':
        kdr_S_origins = 'S:het'
    elif genotypes['kdr-995S'] == 'CC':
        kdr_S_origins = 'S:hom'
    else:
        print(f'Unexpected kdr S genotype. {sample_name} {genotypes["kdr-995S"]}')
        kdr_S_origins = 'Fail. Unexpected kdr S genotype'
    # If the individual has Skdr, find out its origins
    # For S homozygotes
    if kdr_S_origins == 'S:hom':
        if pd.isnull(genotypes['Def-S1-3']):
            kdr_S_origins = f'{kdr_S_origins},S1?'
        elif genotypes['Def-S1-3'] == 'CC':
            kdr_S_origins = f'{kdr_S_origins},S1_hom'
        elif genotypes['Def-S1-3'] == 'CT':
            kdr_S_origins = f'{kdr_S_origins},S1_het'
        #
        if pd.isnull(genotypes['Def-S2S4']):
            kdr_S_origins = f'{kdr_S_origins},S2S4?'
        elif genotypes['Def-S2S4'] == 'TT':
            if pd.isnull(genotypes['Def-S2-4']):
                kdr_S_origins = f'{kdr_S_origins},(S2S4)_hom'
            elif genotypes['Def-S2-4'] == 'AA':
                kdr_S_origins = f'{kdr_S_origins},S2_hom'
            elif genotypes['Def-S2-4'] == 'AT':
                kdr_S_origins = f'{kdr_S_origins},S2_het,S4_het'
            elif genotypes['Def-S2-4'] == 'TT':
                kdr_S_origins = f'{kdr_S_origins},S4_hom'
        elif genotypes['Def-S2S4'] == 'CT':
            if pd.isnull(genotypes['Def-S2-4']):
                kdr_S_origins = f'{kdr_S_origins},(S2S4)_het'
            elif genotypes['Def-S2-4'] == 'AA':
                kdr_S_origins = f'Fail. Genotypes suggest that sample {sample_name} is heterozygote for S2S4, but homozygote for S2.'
            elif genotypes['Def-S2-4'] == 'AT':
                kdr_S_origins = f'{kdr_S_origins},S2_het'
            elif genotypes['Def-S2-4'] == 'TT':
                kdr_S_origins = f'{kdr_S_origins},S4_het'
        #
        if pd.isnull(genotypes['Def-S3']):
            kdr_S_origins = f'{kdr_S_origins},S3?'
        elif genotypes['Def-S3'] == 'GG':
            kdr_S_origins = f'{kdr_S_origins},S3_hom'
        elif genotypes['Def-S3'] == 'GT':
            kdr_S_origins = f'{kdr_S_origins},S3_het'
        # 
        if alternate_S4S5:
            if pd.isnull(genotypes['Def-S4S5-2']):
                kdr_S_origins = f'{kdr_S_origins},S4S5?'
            elif genotypes['Def-S4S5-2'] == 'TT':
                if pd.isnull(genotypes['Def-S5']):
                    kdr_S_origins = f'{kdr_S_origins},(S4S5)_hom'
                elif genotypes['Def-S5'] == 'CC':
                    kdr_S_origins = f'{kdr_S_origins},S5_hom'
                elif genotypes['Def-S5'] == 'AC':
                    kdr_S_origins = f'{kdr_S_origins},S5_het,S4_het'
                elif genotypes['Def-S5'] == 'AA':
                    kdr_S_origins = f'{kdr_S_origins},S4_hom'
            elif genotypes['Def-S4S5-2'] == 'GT':
                if pd.isnull(genotypes['Def-S5']):
                    kdr_S_origins = f'{kdr_S_origins},(S4S5)_het'
                elif genotypes['Def-S5'] == 'CC':
                    kdr_S_origins = f'Fail. Genotypes suggest that sample {sample_name} is heterozygote for S4S5, but homozygote for S5.'
                elif genotypes['Def-S5'] == 'AC':
                    kdr_S_origins = f'{kdr_S_origins},S5_het'
                elif genotypes['Def-S5'] == 'AA':
                    kdr_S_origins = f'{kdr_S_origins},S4_het'
        else :
            if pd.isnull(genotypes['Def-S4S5']):
                kdr_S_origins = f'{kdr_S_origins},S4S5?'
            elif genotypes['Def-S4S5'] == 'CC':
                if pd.isnull(genotypes['Def-S5']):
                    kdr_S_origins = f'{kdr_S_origins},(S4S5)_hom'
                elif genotypes['Def-S5'] == 'CC':
                    kdr_S_origins = f'{kdr_S_origins},S5_hom'
                elif genotypes['Def-S5'] == 'AC':
                    kdr_S_origins = f'{kdr_S_origins},S4_het,S5_het'
                elif genotypes['Def-S5'] == 'AA':
                    kdr_S_origins = f'{kdr_S_origins},S4_hom'
            elif genotypes['Def-S4S5'] == 'CT':
                if pd.isnull(genotypes['Def-S5']):
                    kdr_S_origins = f'{kdr_S_origins},(S4S5)_het'
                elif genotypes['Def-S5'] == 'CC':
                    kdr_S_origins = f'Fail. Genotypes suggest that sample {sample_name} is heterozygote for S4S5, but homozygote for S5.'
                elif genotypes['Def-S5'] == 'AC':
                    kdr_S_origins = f'{kdr_S_origins},S5_het'
                elif genotypes['Def-S5'] == 'AA':
                    kdr_S_origins = f'{kdr_S_origins},S4_het'
    # for S heterozygotes
    elif kdr_S_origins == 'S:het':
        if pd.isnull(genotypes['Def-S1-3']):
            kdr_S_origins = f'{kdr_S_origins},S1?'
        elif genotypes['Def-S1-3'] == 'CC':
            kdr_S_origins = f'Fail. Genotypes suggest that sample {sample_name} is heterozygote for S kdr, but homozygote for S1.'
        elif genotypes['Def-S1-3'] == 'CT':
            kdr_S_origins = f'{kdr_S_origins},S1_het'
        #
        if pd.isnull(genotypes['Def-S2S4']):
            kdr_S_origins = f'{kdr_S_origins},S2S4?'
        elif genotypes['Def-S2S4'] == 'TT':
            kdr_S_origins = f'Fail. Genotypes suggest that sample {sample_name} is heterozygote for S kdr, but homozygote for S2S4.'
        elif genotypes['Def-S2S4'] == 'CT':
            if pd.isnull(genotypes['Def-S2-4']):
                kdr_S_origins = f'{kdr_S_origins},(S2S4)_het'
            elif genotypes['Def-S2-4'] == 'AA':
                kdr_S_origins = f'Fail. Genotypes suggest that sample {sample_name} is heterozygote for S kdr and S2S4, but homozygote for S2.'
            elif genotypes['Def-S2-4'] == 'AT':
                kdr_S_origins = f'{kdr_S_origins},S2_het'
            elif genotypes['Def-S2-4'] == 'TT':
                kdr_S_origins = f'{kdr_S_origins},S4_het'
        #
        if pd.isnull(genotypes['Def-S3']):
            kdr_S_origins = f'{kdr_S_origins},S3?'
        elif genotypes['Def-S3'] == 'GG':
            kdr_S_origins = f'Fail. Genotypes suggest that sample {sample_name} is heterozygote for S kdr, but homozygote for S3.'
        elif genotypes['Def-S3'] == 'GT':
            kdr_S_origins = f'{kdr_S_origins},S3_het'
        # 
        if alternate_S4S5:
            if pd.isnull(genotypes['Def-S4S5_2']):
                kdr_S_origins = f'{kdr_S_origins},S4S5?'
            elif genotypes['Def-S4S5-2'] == 'TT':
                kdr_S_origins = f'Fail. Genotypes suggest that sample {sample_name} is heterozygote for S kdr, but homozygote for S4S5.'
            elif genotypes['Def-S4S5-2'] == 'GT':
                if pd.isnull(genotypes['Def-S5']):
                    kdr_S_origins = f'{kdr_S_origins},(S4S5)_het'
                elif genotypes['Def-S5'] == 'CC':
                    kdr_S_origins = f'Fail. Genotypes suggest that sample {sample_name} is heterozygote for S kdr and S4S5, but homozygote for S5.'
                elif genotypes['Def-S5'] == 'AC':
                    kdr_S_origins = f'{kdr_S_origins},S5_het'
                elif genotypes['Def-S5'] == 'AA':
                    kdr_S_origins = f'{kdr_S_origins},S4_het'
        else :
            if pd.isnull(genotypes['Def-S4S5']):
                kdr_S_origins = f'{kdr_S_origins},S4S5?'
            elif genotypes['Def-S4S5'] == 'CC':
                kdr_S_origins = f'Fail. Genotypes suggest that sample {sample_name} is heterozygote for S kdr, but homozygote for S4S5.'
            elif genotypes['Def-S4S5'] == 'CT':
                if pd.isnull(genotypes['Def-S5']):
                    kdr_S_origins = f'{kdr_S_origins},(S4S5)_het'
                elif genotypes['Def-S5'] == 'CC':
                    return(f'Fail. Genotypes suggest that sample {sample_name} is heterozygote for S kdr and S4S5, but homozygote for S5.')
                elif genotypes['Def-S5'] == 'AC':
                    kdr_S_origins = f'{kdr_S_origins},S5_het'
                elif genotypes['Def-S5'] == 'AA':
                    kdr_S_origins = f'{kdr_S_origins},S4_het'
    if clean:
        return(_kdr_gen_cleanup(kdr_S_origins))
    else:
        return(kdr_S_origins)

def _402_kdr_origin_gen(genotypes, clean = True):
    """Infer 402L kdr haplotype origin calls for one sample.

    Parameters
    ----------
    genotypes : pandas.Series
        Per-sample genotype calls containing the `kdr-402L` marker.
    clean : bool, default=True
        Whether to post-process the raw origin string into a simplified call.

    Returns
    -------
    str
        Raw or cleaned 402L origin assignment for the sample.
    """
    if 'sample_name' in genotypes.index:
        sample_name = genotypes['sample_name']
    else:
        sample_name = genotypes.name
    # Check for the 995F mutations
    if pd.isnull(genotypes['kdr-402L']):
        kdr_402_origins = '402:unknown'
    elif genotypes['kdr-402L'] == 'GG':
        kdr_402_origins = '402:wt_hom'
    elif genotypes['kdr-402L'] == 'CG':
        kdr_402_origins = '402:het,402LC_het'
    elif genotypes['kdr-402L'] == 'GT':
        kdr_402_origins = '402:het,402LT_het'
    elif genotypes['kdr-402L'] == 'CT':
        kdr_402_origins = '402:hom,402LC_het,402LT_het'
    elif genotypes['kdr-402L'] == 'CC':
        kdr_402_origins = '402:hom,402LC_hom'
    elif genotypes['kdr-402L'] == 'TT':
        kdr_402_origins = '402:hom,402LT_hom'
    else:
        print(f'Unexpected kdr 402 genotype. {sample_name} {genotypes["kdr-402L"]}')
        kdr_402_origins = 'Fail. Unexpected kdr 402 genotype'
        
    if clean:
        return(_kdr_gen_cleanup(kdr_402_origins))
    else:
        return(kdr_402_origins)

def _kdr_gen_cleanup(kdr_origin_str):
    """Simplify raw kdr origin strings into diploid haplotype calls.

    Parameters
    ----------
    kdr_origin_str : str
        Raw origin string produced by one of the kdr origin helper functions.

    Returns
    -------
    str
        Comma-separated pair of simplified haplotype calls.
    """
    if re.search('Fail', kdr_origin_str):
        return('?,?')
    if re.search('wt', kdr_origin_str):
        return('wt,wt')
    kdr_type = re.findall('.*(?=:)', kdr_origin_str)[0]
    outcomes = kdr_origin_str.split(',')
    origins = np.unique(outcomes[1:])
    established_origins = [o for o in origins if not re.search(r'\?', o)]
    # Remove "_het" and "_hom" text
    established_origins = [re.sub('_het', '', o) for o in established_origins]
    # Now double up each _hom entry. this is cludgy, but coudn't find a more elegant way
    for i in range(len(established_origins)):
        o = established_origins[i]
        if re.search('_hom', o):
            o = re.sub('_hom', '', o)
            established_origins[i] = o
            established_origins.append(o)
    if re.search('hom', outcomes[0]):
        if len(established_origins) == 1:
            return(f'{kdr_type},{established_origins[0]}')
        elif len(established_origins) == 2:
            return(','.join(established_origins))
        else:
            return(f'{kdr_type},{kdr_type}')
    if re.search('het', outcomes[0]):
        if len(established_origins) == 1:
            return(f'wt,{established_origins[0]}')
        else:
            return(f'wt,{kdr_type}')
    else:
        return('?,?')

def kdr_origin(genotypes, alternate_S4S5 = False, clean = True, include_402 = None):
    """Assemble kdr origin calls for one sample into a dataframe.

    Parameters
    ----------
    genotypes : pandas.Series
        Per-sample genotype calls used to infer kdr origin backgrounds.
    alternate_S4S5 : bool, default=False
        Whether to use the alternate S4/S5 marker naming scheme.
    clean : bool, default=True
        Whether to simplify raw origin strings before returning them.
    include_402 : bool, optional
        Whether to include 402L origin calls. If `None`, this is inferred from
        the available genotype columns.

    Returns
    -------
    pandas.DataFrame
        Single-row dataframe containing the inferred origin calls.
    """
    if 'sample_name' in genotypes.index:
        sample_name = genotypes['sample_name']
    else:
        sample_name = genotypes.name
    if include_402 == None:
        if 'kdr-402L' in genotypes.index:
            include_402 = True
        else:
            include_402 = False
    if include_402 == False:
        kdr_origins = pd.DataFrame({'kdr_F_origin': [_f_kdr_origin_gen(genotypes, clean)], 
                                    'kdr_S_origin': [_s_kdr_origin_gen(genotypes, clean, alternate_S4S5)]
                                    }, index = [sample_name]
        
        )
    else:
        kdr_origins = pd.DataFrame({'kdr_F_origin': [_f_kdr_origin_gen(genotypes, clean)], 
                                    'kdr_S_origin': [_s_kdr_origin_gen(genotypes, clean, alternate_S4S5)],
                                    'kdr_402_origin': [_402_kdr_origin_gen(genotypes, clean)]
                                    }, index = [sample_name]
        )
    return(kdr_origins)

def get_single_gen_call(x):  
    """Collapse separate kdr origin calls into a single genotype summary.

    Parameters
    ----------
    x : pandas.Series
        Row containing cleaned kdr origin calls.

    Returns
    -------
    str
        Comma-separated summary call for the sample.
    """
    if 'kdr_402_origin' in x.index:
        return(_get_single_gen_call_with_402(x))
    else:
        return(_get_single_gen_call_no_402(x))

def _get_single_gen_call_no_402(x): 
    """Summarise F and S origin calls when 402L is absent.

    Parameters
    ----------
    x : pandas.Series
        Row containing `kdr_F_origin` and `kdr_S_origin`.

    Returns
    -------
    str
        Two-part haplotype summary or `"?,?"` when the calls are inconsistent.
    """
    if 'sample_name' in x.index:
        sample_name = x['sample_name']
    else:
        sample_name = x.name
    joined_calls = np.array(x['kdr_F_origin'].split(',') + x['kdr_S_origin'].split(','))
    # There should be at least two 'wt' calls
    if np.sum(joined_calls == 'wt') < 2:
        print(f'Too many different mutant haplotype backgrounds in sample {sample_name}')
        return('?,?')
    # Otherwise, drop two wildtype calls
    else:
        which_drop = np.where(joined_calls == 'wt')[0][:2]
        return(','.join(np.delete(joined_calls, which_drop)))

def _get_single_gen_call_with_402(x): 
    """Summarise F, S, and 402L origin calls into one genotype label.

    Parameters
    ----------
    x : pandas.Series
        Row containing `kdr_F_origin`, `kdr_S_origin`, and `kdr_402_origin`.

    Returns
    -------
    str
        Two-part haplotype summary or `"?,?"` when the calls are inconsistent.
    """
    if 'sample_name' in x.index:
        sample_name = x['sample_name']
    else:
        sample_name = x.name
    joined_calls = np.array(x['kdr_F_origin'].split(',') + 
                            x['kdr_S_origin'].split(',') +
                            x['kdr_402_origin'].split(',')
    )
    # There should be at least four 'wt' calls
    if np.sum(np.isin(joined_calls,  ['wt', '?'])) < 4:
        print(f'Too many different mutant haplotype backgrounds in sample {sample_name}')
        return('?,?')
    # Otherwise, drop four wildtype calls
    else:
        which_drop = np.concatenate([
            np.where(joined_calls == 'wt')[0],
            np.where(joined_calls == '?')[0]
        ])[:4]
        return(','.join(np.delete(joined_calls, which_drop)))

def signif(x, n_figs):
    """Round values to a fixed number of significant figures.

    Parameters
    ----------
    x : array-like
        Numeric input values.
    n_figs : int
        Number of significant figures to retain.

    Returns
    -------
    numpy.ndarray or scalar
        Rounded values with the same broadcasted shape as `x`.
    """
    power = 10 ** np.floor(np.log10(np.abs(x).clip(1e-200)))
    rounded = np.round(x / power, n_figs - 1) * power
    return rounded

def _dipclust_concat_subplots(
    figures,
    width,
    height,
    row_heights,
    title,
    xaxis_range,
):
    """Stack dendrogram-related traces into a shared subplot figure.

    Parameters
    ----------
    figures : sequence
        Plotly figures or traces to append vertically.
    width : int
        Output figure width in pixels.
    height : int
        Output figure height in pixels.
    row_heights : sequence of float
        Relative heights for each subplot row.
    title : str
        Figure title.
    xaxis_range : tuple or list
        Shared x-axis range applied to the composed figure.

    Returns
    -------
    plotly.graph_objects.Figure
        Combined subplot figure.
    """
    from plotly.subplots import make_subplots  # type: ignore
    import plotly.graph_objects as go  # type: ignore

    # make subplots
    fig = make_subplots(
        rows=len(figures),
        cols=1,
        shared_xaxes=True,
        vertical_spacing=0.02,
        row_heights=row_heights,
    )

    for i, figure in enumerate(figures):
        if isinstance(figure, go.Figure):
            # This is a figure, access the traces within it.
            for trace in range(len(figure["data"])):
                fig.append_trace(figure["data"][trace], row=i + 1, col=1)
        else:
            # Assume this is a trace, add directly.
            fig.append_trace(figure, row=i + 1, col=1)

    fig.update_xaxes(visible=False)
    fig.update_layout(
        title=title,
        width=width,
        height=height,
        hovermode="closest",
        plot_bgcolor="white",
        xaxis_range=xaxis_range,
    )

    return fig

def plot_dendrogram(
    dist,
    linkage_method,
    count_sort,
    distance_sort,
    render_mode,
    width,
    height,
    title,
    line_width,
    line_color,
    marker_size,
    leaf_data,
    leaf_hover_name,
    leaf_hover_data,
    leaf_color,
    leaf_symbol,
    leaf_y,
    leaf_color_discrete_map,
    leaf_category_orders,
    template,
    y_axis_title,
    y_axis_buffer,
):
    """Plot a hierarchical clustering dendrogram with annotated leaves.

    Parameters
    ----------
    dist : array-like
        Condensed distance vector passed to SciPy linkage.
    linkage_method : str
        Linkage method used by `scipy.cluster.hierarchy.linkage`.
    count_sort : str or bool
        Leaf ordering option passed to `scipy.cluster.hierarchy.dendrogram`.
    distance_sort : str or bool
        Distance ordering option passed to `scipy.cluster.hierarchy.dendrogram`.
    render_mode : str
        Plotly render mode for lines and points.
    width : int
        Figure width in pixels.
    height : int
        Figure height in pixels.
    title : str
        Figure title.
    line_width : float
        Width applied to dendrogram line traces.
    line_color : str
        Colour applied to dendrogram line traces.
    marker_size : float
        Marker size for leaf points.
    leaf_data : pandas.DataFrame
        Metadata for leaves in the same order as the samples in `dist`.
    leaf_hover_name : str
        Column shown as the primary hover label for leaves.
    leaf_hover_data : list or dict
        Additional columns shown in leaf hover labels.
    leaf_color : str
        Column used to colour leaf markers.
    leaf_symbol : str
        Column used to set marker symbols.
    leaf_y : float
        Y position used to place leaf markers.
    leaf_color_discrete_map : dict
        Plotly colour map for leaf categories.
    leaf_category_orders : dict
        Plotly category ordering for leaf annotations.
    template : str
        Plotly template name.
    y_axis_title : str
        Y-axis title.
    y_axis_buffer : float
        Padding added to the y-axis range.

    Returns
    -------
    tuple
        Tuple of `(fig, leaf_data)` where `leaf_data` has been reordered to
        match the plotted dendrogram leaves.
    """
    import scipy.cluster.hierarchy as sch
    # Hierarchical clustering.
    Z = sch.linkage(dist, method=linkage_method)

    # Compute the dendrogram but don't plot it.
    dend = sch.dendrogram(
        Z,
        count_sort=count_sort,
        distance_sort=distance_sort,
        no_plot=True,
    )

    # Compile the line coordinates into a single dataframe.
    icoord = dend["icoord"]
    dcoord = dend["dcoord"]
    line_segments_x = []
    line_segments_y = []
    for ik, dk in zip(icoord, dcoord):
        # Adding None here breaks up the lines.
        line_segments_x += ik + [None]
        line_segments_y += dk + [None]
    df_line_segments = pd.DataFrame({"x": line_segments_x, "y": line_segments_y})

    # Convert X coordinates to haplotype indices (scipy multiplies coordinates by 10).
    df_line_segments["x"] = (df_line_segments["x"] - 5) / 10

    # Plot the lines.
    fig = px.line(
        df_line_segments,
        x="x",
        y="y",
        render_mode=render_mode,
        template=template,
    )

    # Reorder leaf data to align with dendrogram.
    leaves = dend["leaves"]
    n_leaves = len(leaves)
    leaf_data = leaf_data.iloc[leaves]

    # Add scatter plot to draw the leaves.
    fig.add_traces(
        list(
            px.scatter(
                data_frame=leaf_data,
                x=np.arange(n_leaves),
                y=np.repeat(leaf_y, n_leaves),
                color=leaf_color,
                symbol=leaf_symbol,
                render_mode=render_mode,
                hover_name=leaf_hover_name,
                hover_data=leaf_hover_data,
                template=template,
                color_discrete_map=leaf_color_discrete_map,
                category_orders=leaf_category_orders,
            ).select_traces()
        )
    )

    # Style the lines and markers.
    line_props = dict(
        width=line_width,
        color=line_color,
    )
    marker_props = dict(
        size=marker_size,
    )
    fig.update_traces(line=line_props, marker=marker_props)

    # Style the figure.
    fig.update_layout(
        width=width,
        height=height,
        title=title,
        autosize=True,
        hovermode="closest",
        # I cannot get the xaxis title to appear below the plot, and when
        # it's above the plot it often overlaps the title, so hiding it
        # for now.
        xaxis_title=None,
        yaxis_title=y_axis_title,
        showlegend=True,
    )

    # Style axes.
    fig.update_xaxes(
        mirror=False,
        showgrid=False,
        showline=False,
        showticklabels=False,
        ticks="",
        range=(-2, n_leaves + 2),
    )
    fig.update_yaxes(
        mirror=False,
        showgrid=False,
        showline=False,
        showticklabels=True,
        ticks="outside",
        range=(leaf_y - y_axis_buffer, np.max(dcoord) + y_axis_buffer),
    )

    return fig, leaf_data
