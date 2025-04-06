def load_metadata(metadata_path, from_sample_sheet=False, write=False):
    import os
    import pandas as pd
    import re
    import warnings

    if metadata_path.endswith("SampleSheet.csv"): 
        from_sample_sheet = True

    if not from_sample_sheet:
        metadata = pd.read_csv(metadata_path, sep="\t")
        
        # Check that metadata contains the column 'sample_id'
        assert 'sample_id' in metadata.columns, "Metadata must contain 'sample_id' column"
        
        # Check for duplicates in sample_id
        duplicates = metadata['sample_id'].duplicated()
        assert not duplicates.any(), f"Duplicate sample_ids found: {metadata['sample_id'][duplicates].tolist()}"
        
        # Ensure sample_id is a string
        metadata['sample_id'] = metadata['sample_id'].astype(str)
        
        # Check for NaN or empty values in sample_id
        assert not metadata['sample_id'].isna().any(), "NaN values found in 'sample_id' column"
        assert not (metadata['sample_id'] == '').any(), "Empty values found in 'sample_id' column"
        
    else:
        metadata = metadata_from_sample_sheet(metadata_path)
        
        # Check required columns for sample sheet
        required_columns = ['sample_id', 'sample_name', 'index', 'index2']
        for col in required_columns:
            assert col in metadata.columns, f"Sample sheet must contain '{col}' column"
        
        # Check for duplicates in sample_id
        duplicates = metadata['sample_id'].duplicated()
        assert not duplicates.any(), f"Duplicate sample_ids found: {metadata['sample_id'][duplicates].tolist()}"
        
        # Ensure sample_id == sample_name
        assert (metadata['sample_id'] == metadata['sample_name']).all(), "sample_id must match sample_name in sample sheet"
        
        # Ensure sample_id and sample_name are strings
        metadata['sample_id'] = metadata['sample_id'].astype(str)
        metadata['sample_name'] = metadata['sample_name'].astype(str)
        
        # Check for NaN or empty values in critical columns
        for col in ['sample_id', 'sample_name']:
            assert not metadata[col].isna().any(), f"NaN values found in '{col}' column"
            assert not (metadata[col] == '').any(), f"Empty values found in '{col}' column"
        
        # Validate that index sequences only contain valid nucleotides (A, C, G, T)
        for idx_col in ['index', 'index2']:
            if idx_col in metadata.columns:
                # Skip validation for empty index2 (some protocols only use one index)
                if idx_col == 'index2' and metadata[idx_col].isna().all():
                    continue
                    
                # Fill NaN values with empty string for validation
                indices = metadata[idx_col].fillna('')
                invalid_indices = [i for i in indices if i and not re.match(r'^[ACGT]+$', i)]
                assert len(invalid_indices) == 0, f"Invalid nucleotides found in {idx_col}: {invalid_indices}"
    
    if write:
        os.makedirs("results/config", exist_ok=True)
        if 'sample_ID' in metadata.columns:
            metadata = metadata.rename(columns={'sample_ID':'sample_id'})
        metadata.to_csv("results/config/metadata.tsv", sep="\t", index=False)

    return metadata


def metadata_from_sample_sheet(sample_sheet_path):
    with open(sample_sheet_path, 'r') as f:
        content = f.read()
    
    data_section = re.search(r'\[(BCLConvert_)?Data\],+\n(.*?)($|\[)', content, re.DOTALL)
    if not data_section:
        raise ValueError("Could not find [Data] or [BCLConvert_Data] section in the sample sheet")
    
    data_content = data_section.group(2)
    lines = data_content.strip().split('\n')

    df = pd.DataFrame([line.split(',') for line in lines[1:]], columns=lines[0].split(','))
    df = df.rename(columns={'sample_ID': 'sample_id'})
    
    # Split well column into well_letter and well_number
    df['well_letter'] = df['well'].str[0]
    df['well_number'] = df['well'].str[1:].astype(str)

    return df.drop(columns=['well']).rename(columns={'plate_name':'plate'})