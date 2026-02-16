import allel
import numba
import numpy as np
import pandas as pd
import plotly.express as px
from sklearn.neighbors import KNeighborsClassifier

def _aims_n_alt(gt, aim_alts, data_alts):
    n_sites = gt.shape[0]
    n_samples = gt.shape[1]
    # create empty array
    aim_n_alt = np.empty((n_sites, n_samples), dtype=np.int8)

    # for every site
    for i in range(n_sites):
        # find the index of the correct tag snp allele
        tagsnp_index = np.where(aim_alts[i] == data_alts[i])[0]
        for j in range(n_samples):
            if tagsnp_index.shape[0] == 0:
                aim_n_alt[i, j] = 0
                continue

            n_tag_alleles = np.sum(gt[i, j] == tagsnp_index[0])
            n_missing = np.sum(gt[i, j] == -1)
            if n_missing != 0:
                aim_n_alt[i,j] = -1
            else:
                aim_n_alt[i, j] = n_tag_alleles

    return aim_n_alt

def pca_all_samples(gn_wgs, gn_amp, df_wgs_samples, df_amp_samples, n_components=6, cohort_cols='sample_id'):
    """
    Perform PCA on both WGS and amplicon samples combined
    
    Parameters:
    -----------
    gn_wgs : GenotypeArray
        Genotypes of samples with known taxon
    gn_amp : GenotypeArray
        Genotypes of samples to be assigned
    df_wgs_samples : DataFrame
        Metadata for samples with known taxon, including 'taxon' column
    n_components : int
        Number of principal components to compute
        
    Returns:
    --------
    pca_df : DataFrame
        DataFrame containing PCA coordinates and metadata for all samples
    """
    # Create metadata for amplicon samples
    df_amp_samples = df_amp_samples.assign(sample_type='amplicon')
    
    # Add sample type to WGS metadata
    df_wgs_samples = df_wgs_samples.copy()
    df_wgs_samples = df_wgs_samples.assign(sample_type='WGS')
    # Combine metadata
    df_all_samples = pd.concat([df_wgs_samples.reset_index(), df_amp_samples[['sample_id', 'sample_type'] + cohort_cols.split(",")]], ignore_index=True)
    
    # Combine genotype data
    print(f"Combining {gn_wgs.shape[1]} WGS samples and {gn_amp.shape[1]} amplicon samples")
    combined_geno = allel.GenotypeArray(np.concatenate([gn_wgs[:], gn_amp[:]], axis=1))
    
    # Perform PCA
    print("Performing PCA on combined data")
    gn_alt = combined_geno.to_n_alt()
    print("Removing any invariant sites")
    loc_var = np.any(gn_alt != gn_alt[:, 0, np.newaxis], axis=1)
    gn_var = np.compress(loc_var, gn_alt, axis=0)
    
    print(f"Running PCA with {gn_var.shape[0]} variable sites and {gn_var.shape[1]} samples")
    coords, model = allel.pca(gn_var, n_components=n_components)
    
    # Flip axes back so PC1 is same orientation in each window
    for i in range(n_components):
        c = coords[:, i]
        if np.abs(c.min()) > np.abs(c.max()):
            coords[:, i] = c * -1
    
    # Create PCA DataFrame
    pca_df = pd.DataFrame(coords)
    pca_df.columns = [f"PC{pc+1}" for pc in range(n_components)]
    pca_df = pd.concat([df_all_samples.reset_index(drop=True), pca_df], axis=1)
    
    print("PCA completed successfully")
    return pca_df

def assign_taxa(pca_df, method='knn', n_neighbors=5, probability_threshold=0.8, **kwargs):
    """
    Assign taxa to amplicon samples based on PC1-4 coordinates
    
    Parameters:
    -----------
    pca_df : DataFrame
        DataFrame with PCA coordinates and metadata for all samples
    method : str
        Classification method to use ('knn' or 'svm')
    n_neighbors : int
        Number of neighbors to use for KNN classification (ignored if method='svm')
    probability_threshold : float
        Minimum probability required for taxon assignment
    **kwargs : dict
        Additional parameters for the classifier
    
    Returns:
    --------
    assignment_df : DataFrame
        DataFrame with taxon assignments for amplicon samples
    """
    from sklearn.svm import SVC
    
    # Separate training (WGS) and test (amplicon) data
    wgs_samples = pca_df[pca_df['sample_type'] == 'WGS']
    amp_samples = pca_df[pca_df['sample_type'] == 'amplicon']
    
    print(f"Training data: {len(wgs_samples)} WGS samples")
    print(f"Test data: {len(amp_samples)} amplicon samples")
    
    # Extract features for training (PC1-4)
    pc_features = ['PC1', 'PC2', 'PC3', 'PC4']
    X_train = wgs_samples[pc_features].values
    y_train = wgs_samples['taxon'].values
    
    # Extract features for testing
    X_test = amp_samples[pc_features].values
    
    # Initialize classifier based on method
    if method.lower() == 'knn':
        print(f"Training KNN classifier with {n_neighbors} neighbors")
        classifier = KNeighborsClassifier(n_neighbors=n_neighbors, **kwargs)
    elif method.lower() == 'svm':
        print("Training SVM classifier")
        # Ensure probability=True for SVM to get prediction probabilities
        svm_kwargs = {'probability': True}
        svm_kwargs.update(kwargs)
        classifier = SVC(**svm_kwargs)
    else:
        raise ValueError(f"Unsupported method: {method}. Choose 'knn' or 'svm'")
    
    # Train classifier
    classifier.fit(X_train, y_train)
    
    # Predict taxa
    predictions = classifier.predict(X_test)
    probabilities = classifier.predict_proba(X_test)
    max_probs = np.max(probabilities, axis=1)
    
    # Create assignment DataFrame
    assignment_df = amp_samples[['sample_id']].copy()
    
    # Apply probability threshold for assignment
    assignment_df['predicted_taxon'] = ['unassigned'] * len(amp_samples)
    assignment_df['probability'] = 0.0
    assignment_df['classifier'] = method.lower()
    
    for i, (pred, prob) in enumerate(zip(predictions, max_probs)):
        if prob >= probability_threshold:
            assignment_df.loc[assignment_df.index[i], 'predicted_taxon'] = pred
            assignment_df.loc[assignment_df.index[i], 'probability'] = prob
        else:
            # For low confidence assignments, still store the prediction but mark as unassigned
            assignment_df.loc[assignment_df.index[i], 'low_confidence_prediction'] = pred
            assignment_df.loc[assignment_df.index[i], 'low_confidence_probability'] = prob
    
    # Count assignments
    assigned_count = sum(assignment_df['predicted_taxon'] != 'unassigned')
    print(f"Assigned {assigned_count} out of {len(amp_samples)} samples with confidence ≥ {probability_threshold}")
    
    if assigned_count < len(amp_samples):
        print(f"{len(amp_samples) - assigned_count} samples were below the confidence threshold")
    
    # Print summary of assignments
    print("\nTaxon assignment summary:")
    print(assignment_df['predicted_taxon'].value_counts())
    
    return assignment_df

def plot_pca_3d_with_assignments(pca_df, assignment_df, title="PCA with Taxon Assignments", 
                              height=800, width=800):
    """
    Create a 3D plot of PCA results with taxon assignments
    
    Parameters:
    -----------
    pca_df : DataFrame
        DataFrame with PCA coordinates and metadata
    assignment_df : DataFrame
        DataFrame with taxon assignments for amplicon samples
    title : str
        Plot title
    height, width : int
        Plot dimensions
    
    Returns:
    --------
    plotly.graph_objects.Figure
    """
    import plotly.graph_objects as go
    
    # Create a copy of the PCA DataFrame for visualization
    plot_df = pca_df.copy()
    
    # Update amplicon samples with assigned taxa
    for i, row in assignment_df.iterrows():
        sample_id = row['sample_id']
        if row['predicted_taxon'] != 'unassigned':
            # Update the taxon for this sample
            plot_df.loc[plot_df['sample_id'] == sample_id, 'display_taxon'] = row['predicted_taxon'] + " (predicted)"
        else:
            # Mark as unassigned
            plot_df.loc[plot_df['sample_id'] == sample_id, 'display_taxon'] = 'unassigned'
    
    # For WGS samples, use the known taxon but add a "reference" label
    plot_df.loc[plot_df['sample_type'] == 'WGS', 'display_taxon'] = plot_df.loc[plot_df['sample_type'] == 'WGS', 'taxon'] + " (VObs reference)"
    
    # Split the data by sample type for different marker styles
    wgs_samples = plot_df[plot_df['sample_type'] == 'WGS']
    amp_samples = plot_df[plot_df['sample_type'] == 'amplicon']
    
    # Create an empty figure
    fig = go.Figure()
    
    # Add reference (WGS) samples - use circles
    for taxon in wgs_samples['display_taxon'].unique():
        subset = wgs_samples[wgs_samples['display_taxon'] == taxon]
        fig.add_trace(go.Scatter3d(
            x=subset['PC1'],
            y=subset['PC2'],
            z=subset['PC3'],
            mode='markers',
            marker=dict(
                size=5,
                symbol='circle',
                opacity=0.7
            ),
            name=taxon,
            hovertemplate="<b>%{text}</b><br>PC1: %{x:.2f}<br>PC2: %{y:.2f}<br>PC3: %{z:.2f}<br>Sample Type: Reference<br>Taxon: " + taxon.replace(" (reference)", "") + "<extra></extra>",
            text=subset['sample_id']
        ))
    
    # Add amplicon samples - use diamonds/cross for greater visibility
    for taxon in amp_samples['display_taxon'].unique():
        subset = amp_samples[amp_samples['display_taxon'] == taxon]
        
        # Skip unassigned for a moment
        if taxon == 'unassigned':
            continue
        
        # For assigned samples, display probability
        probs = []
        for sid in subset['sample_id']:
            prob = assignment_df.loc[assignment_df['sample_id'] == sid, 'probability'].values[0]
            probs.append(f"{prob:.2f}")
        
        fig.add_trace(go.Scatter3d(
            x=subset['PC1'],
            y=subset['PC2'],
            z=subset['PC3'],
            mode='markers',
            marker=dict(
                size=7,
                symbol='diamond',
                opacity=0.9
            ),
            name=taxon,
            hovertemplate="<b>%{text}</b><br>PC1: %{x:.2f}<br>PC2: %{y:.2f}<br>PC3: %{z:.2f}<br>Sample Type: Amplicon<br>Assigned Taxon: " + taxon.replace(" (predicted)", "") + "<br>Probability: %{customdata}<extra></extra>",
            text=subset['sample_id'],
            customdata=probs
        ))
    
    # Add unassigned amplicon samples
    unassigned = amp_samples[amp_samples['display_taxon'] == 'unassigned']
    if len(unassigned) > 0:
        # For unassigned, show low confidence prediction if available
        hover_data = []
        for sid in unassigned['sample_id']:
            row = assignment_df.loc[assignment_df['sample_id'] == sid].iloc[0]
            if 'low_confidence_prediction' in row and not pd.isna(row['low_confidence_prediction']):
                hover_data.append(f"{row['low_confidence_prediction']} ({row['low_confidence_probability']:.2f})")
            else:
                hover_data.append("No prediction")
        
        fig.add_trace(go.Scatter3d(
            x=unassigned['PC1'],
            y=unassigned['PC2'],
            z=unassigned['PC3'],
            mode='markers',
            marker=dict(
                size=7,
                symbol='x',
                color='gray',
                opacity=0.7
            ),
            name='unassigned',
            hovertemplate="<b>%{text}</b><br>PC1: %{x:.2f}<br>PC2: %{y:.2f}<br>PC3: %{z:.2f}<br>Sample Type: Amplicon<br>Status: Unassigned<br>Low confidence: %{customdata}<extra></extra>",
            text=unassigned['sample_id'],
            customdata=hover_data
        ))
    
    # Improve layout
    fig.update_layout(
        title=title,
        scene=dict(
            xaxis_title='PC1',
            yaxis_title='PC2',
            zaxis_title='PC3'
        ),
        height=height,
        width=width,
        margin=dict(l=0, r=0, b=0, t=40),
        legend_title_text='Taxa'
    )
    
    return fig

def _melt_gt_counts(gt_counts):
    n_snps, n_samples, n_alleles = gt_counts.shape
    # Use a float array to allow NaN values
    melted_counts = np.full((n_snps * (n_alleles - 1), n_samples), np.nan, dtype=np.float64)

    for i in range(n_snps):
        for j in range(n_samples):
            for k in range(n_alleles - 1):
                # Check if the genotype count is valid (== 2)
                if gt_counts[i][j].sum() == 2:
                    melted_counts[(i * (n_alleles - 1)) + k][j] = gt_counts[i][j][k + 1]
                else:
                    # Assign NaN for missing or invalid data
                    melted_counts[(i * (n_alleles - 1)) + k][j] = np.nan

    return melted_counts

def get_consensus_taxon(row):
    # Extract the taxon values
    aim_taxon = row['aim_taxon']
    pca_taxon = row['pca_taxon']
    tree_taxon = row['tree_taxon']
    
    # Count occurrences of each taxon
    taxon_counts = {}
    for taxon in [aim_taxon, pca_taxon, tree_taxon]:
        if pd.notna(taxon):  # Skip NaN values
            taxon_counts[taxon] = taxon_counts.get(taxon, 0) + 1
    
    # Find the most common taxon
    max_count = 0
    max_taxon = None
    for taxon, count in taxon_counts.items():
        if count > max_count:
            max_count = count
            max_taxon = taxon
    
    # If at least two columns agree, return that taxon, otherwise 'unassigned'
    if max_count >= 2:
        return max_taxon
    else:
        return 'unassigned'
