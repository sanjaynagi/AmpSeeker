import numpy as np
import pandas as pd
import plotly.express as px
import warnings

warnings.filterwarnings('ignore')

import allel
from scipy.spatial.distance import squareform

import shared

def plot_pca(pca_df, colour_column, cohort_columns, dataset,  x='PC1',y='PC2',z='PC3', color_mapping=None, height=500, width=750):
    fig= px.scatter_3d(
        pca_df, 
        x=x, 
        y=y,
        z=z, 
        title=f"PCA {dataset} | PC1 vs PC2 vs PC3 coloured by {colour_column}",
        color=colour_column, 
        hover_data=cohort_columns + ['sample_id'],
        color_discrete_map=color_mapping[colour_column], 
        template='simple_white',
        height=height,
        width=width
    )

    return fig


def compute_njt_inputs(geno, metadata, cohort_col):
    # Find segregating sites and remove highly missing sites.
    ac = geno.count_alleles()
    seg = ac.is_segregating()
    gn_seg = geno.compress(seg, axis=0)
    missing_mask = gn_seg.is_missing().sum(axis=1) > gn_seg.shape[1] * 0.4
    gn_seg = gn_seg.compress(~missing_mask, axis=0)

    ac = allel.GenotypeArray(gn_seg).to_allele_counts(max_allele=3)
    x = np.ascontiguousarray(np.swapaxes(ac.values, 0, 1))

    dists = shared.multiallelic_diplotype_pdist(x, metric=shared.multiallelic_diplotype_mean_cityblock)
    dists = squareform(dists)
    df_samples = metadata.set_index('sample_id')
    df = df_samples[[cohort_col]]

    df_dist_matrix = pd.DataFrame(dists, index=df_samples.index.to_list(), columns=df_samples.index.to_list())
    df_dists = df_dist_matrix.stack().reset_index().set_axis('sample_id_x sample_id_y distance'.split(), axis=1)
    df_dists = df_dists.merge(df, left_on='sample_id_x', right_index=True).merge(
        df, left_on='sample_id_y', right_index=True, suffixes=('_x', '_y')
    )
    df_dists = df_dists[df_dists['sample_id_x'] != df_dists['sample_id_y']]
    df_dists = df_dists.assign(
        dedup=np.array([''.join(sorted([a, b])) for a, b in zip(df_dists.sample_id_x, df_dists.sample_id_y)]).astype(str)
    )
    df_dists = df_dists.sort_values('sample_id_x').drop_duplicates('dedup').drop('dedup', axis=1)
    df_dists[cohort_col] = df_dists[f'{cohort_col}_x'] + ' | ' + df_dists[f'{cohort_col}_y']
    df_dists = df_dists.drop([f'{cohort_col}_x', f'{cohort_col}_y'], axis=1)
    df_grp_dists = df_dists.groupby(cohort_col).agg({'distance': 'mean'}).sort_values('distance').rename(
        columns={'distance': 'mean_distance'}
    ).reset_index()
    df_dists = df_dists.merge(df_grp_dists, on=cohort_col).assign(
        normalised_dist=lambda x: x.distance - x.mean_distance
    ).sort_values('normalised_dist')

    far_samples = df_dists.sort_values('normalised_dist', ascending=False)[: int(df_dists.shape[0] * 0.005)][
        ['sample_id_x', 'sample_id_y']
    ].values.flatten()
    far_samples, far_counts = np.unique(far_samples, return_counts=True)
    exclude_outliers = far_samples[far_counts > int(df_samples.shape[0] * 0.1)]

    dists = df_dist_matrix.drop(exclude_outliers, axis=0).drop(exclude_outliers, axis=1).values
    leaf_data = df_samples.query('sample_id not in @exclude_outliers').reset_index()

    nan_mask = np.isnan(dists).sum(axis=0) > 0
    dists = dists[~nan_mask][:, ~nan_mask]

    return dists, leaf_data, exclude_outliers


def run_njt_analysis(geno, metadata, cohort_cols, cohort_col, color_mapping, wkdir):
    import anjl

    dists, leaf_data, exclude_outliers = compute_njt_inputs(geno=geno, metadata=metadata, cohort_col=cohort_col)
    print(f'excluding extreme outliers from NJT {exclude_outliers}')

    z = anjl.dynamic_nj(dists)
    figures = []
    for col in cohort_cols:
        fig = anjl.plot(
            z,
            leaf_data=leaf_data,
            color=col,
            hover_name='sample_id',
            hover_data=cohort_cols,
            color_discrete_map=color_mapping[col],
            marker_size=8,
        )
        fig.write_image(f"{wkdir}/results/njt_{col}.png", scale=2)
        figures.append(fig)

    return figures
