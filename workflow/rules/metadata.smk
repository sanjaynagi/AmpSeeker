def load_metadata(metadata_path, from_sample_sheet=False, write=False):
    import os

    if isinstance(metadata_path, str) and metadata_path.endswith("SampleSheet.csv"):
        from_sample_sheet = True

    if not from_sample_sheet:
        if isinstance(metadata_path, list):
            raise ValueError("metadata must be a single file path when from-bcl is False.")
        if os.path.exists(metadata_path):
            metadata = pd.read_csv(metadata_path, sep="\t")
        else:
            raise ValueError(f"Metadata file '{metadata_path}' does not exist at the path provided in the config. Is the path correct, or do you mean to run from_bcl=True instead?")
    else:
        if isinstance(metadata_path, list):
            if len(metadata_path) == 0:
                raise ValueError("No SampleSheet.csv paths were provided for from-bcl mode.")
            metadata = metadata_from_sample_sheets(metadata_path)
        else:
            metadata = metadata_from_sample_sheet(metadata_path)

    assert all(col in metadata.columns for col in config['cohort-columns']), "Not all provided cohort columns are present in the metadata file" 
    
    if write:
        os.makedirs("results/config", exist_ok=True)
        if 'sample_ID' in metadata.columns:
            metadata = metadata.rename(columns={'sample_ID':'sample_id'})
        metadata.to_csv("results/config/metadata.tsv", sep="\t", index=False)

    return metadata


def metadata_from_sample_sheets(sample_sheet_paths):
    import os

    metadata = []
    for sample_sheet_path in sample_sheet_paths:
        if not os.path.exists(sample_sheet_path):
            raise ValueError(f"SampleSheet file '{sample_sheet_path}' does not exist.")
        df = metadata_from_sample_sheet(sample_sheet_path).copy()
        df["illumina_run_dir"] = os.path.dirname(sample_sheet_path)
        metadata.append(df)

    metadata = pd.concat(metadata, ignore_index=True)

    duplicate_ids = metadata["sample_id"].duplicated(keep=False)
    if duplicate_ids.any():
        duplicated = metadata.loc[duplicate_ids].copy()
        metadata_cols = [c for c in duplicated.columns if c not in ["illumina_run_dir"]]
        conflict_samples = []

        for sample_id, group in duplicated.groupby("sample_id"):
            group = group[metadata_cols]
            # Ignore run-level duplication if metadata values are identical.
            if any(group[col].dropna().astype(str).nunique() > 1 for col in group.columns if col != "sample_id"):
                conflict_samples.append(sample_id)

        if len(conflict_samples) > 0:
            raise ValueError(
                "Duplicate sample_id values with conflicting metadata were found across SampleSheets: "
                + ", ".join(map(str, conflict_samples))
            )

        metadata = metadata.drop_duplicates(subset=["sample_id"], keep="first")

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
    
    if 'well' in df.columns:
        # Split well column into well_letter and well_number
        df['well_letter'] = df['well'].str[0]
        df['well_number'] = df['well'].str[1:].astype(str)
        df = df.drop(columns=['well'])
    
    if 'plate_name' in df.columns:
        df = df.rename(columns={'plate_name':'plate'})

    return df
