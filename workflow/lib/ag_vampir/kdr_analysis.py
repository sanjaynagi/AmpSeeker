import numpy as np
import pandas as pd
import plotly.express as px
import re

def _f_kdr_origin_gen(genotypes, clean = True):
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
    if 'kdr_402_origin' in x.index:
        return(_get_single_gen_call_with_402(x))
    else:
        return(_get_single_gen_call_no_402(x))

def _get_single_gen_call_no_402(x): 
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
