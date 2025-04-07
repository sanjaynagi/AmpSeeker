def load_metadata(metadata_path, from_sample_sheet=False, write=False):
    import os

    if metadata_path.endswith("SampleSheet.csv"): 
        from_sample_sheet = True

    if not from_sample_sheet:
        metadata = pd.read_csv(metadata_path, sep="\t")
    else:
        metadata = metadata_from_sample_sheet(metadata_path)
    
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
    
    if 'well' in df.columns:
        # Split well column into well_letter and well_number
        df['well_letter'] = df['well'].str[0]
        df['well_number'] = df['well'].str[1:].astype(str)
        df = df.drop(columns=['well'])
    
    if 'plate_name' in df.columns:
        df = df.rename(columns={'plate_name':'plate'})

    return df