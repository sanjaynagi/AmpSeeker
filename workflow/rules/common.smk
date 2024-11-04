import pandas as pd


# Split into two sample sets as bcftools merge cant take over 1000 files
# So we must do two rounds of merging
if len(metadata) > 1000:
    large_sample_size = True
    n_samples = len(metadata)
    half = int(n_samples / 2)
    samples1 = metadata["sample_id"][:half]
    samples2 = metadata["sample_id"][half:]
else:
    large_sample_size = False
    samples1 = []
    samples2 = []


def load_metadata(metadata_path):
    # load panel metadata
    if metadata_path.endswith(".xlsx"):
        metadata = pd.read_excel(metadata_path, engine="openpyxl")
    elif metadata_path.endswith(".tsv"):
        metadata = pd.read_csv(metadata_path, sep="\t")
    elif metadata_path.endswith(".csv"):
        metadata = pd.read_csv(metadata_path, sep=",")
    else:
        raise ValueError("Metadata file must be .xlsx or .csv")
    return metadata


rule set_kernel:
    input:
        f"{workflow.basedir}/envs/AmpSeeker-python.yaml",
    output:
        touch("results/.kernel.set"),
    conda:
        f"{workflow.basedir}/envs/AmpSeeker-python.yaml"
    log:
        "logs/set_kernel.log",
    shell:
        """
        python -m ipykernel install --user --name=AmpSeq_python 2> {log}
        """

def generate_distinct_colors(n):
    """
    Generate n visually distinct colors using HSV color space.
    """
    import colorsys
    colors = []
    for i in range(n):
        # Use golden ratio to space hues evenly
        hue = i * 0.618033988749895 % 1
        # Use fixed saturation and value for consistency
        saturation = 0.7
        value = 0.95
        # Convert HSV to RGB
        rgb = colorsys.hsv_to_rgb(hue, saturation, value)
        # Convert RGB to hex
        hex_color = "#{:02x}{:02x}{:02x}".format(
            int(rgb[0] * 255),
            int(rgb[1] * 255),
            int(rgb[2] * 255)
        )
        colors.append(hex_color)
    return colors

def create_plotly_color_mapping(
    df_samples,
    columns,
    color_sequence='Plotly',
):
    """
    Create color mappings for categorical columns using Plotly color sequences.
    """
    import plotly.express as px
    # Get the color sequence
    if isinstance(color_sequence, str):
        if hasattr(px.colors.qualitative, color_sequence):
            colors = getattr(px.colors.qualitative, color_sequence)
        else:
            raise ValueError(f"Unknown color sequence: {color_sequence}. "
                           f"Available sequences: {dir(px.colors.qualitative)}")
    else:
        colors = color_sequence
    
    color_mappings = {}
    for col in columns:
        unique_values = df_samples[col].unique()
        n_categories = len(unique_values)
        n_colors = len(colors)
        
        if n_categories > n_colors:
            raise ValueError(
                f"Column '{col}' has {n_categories} categories but color sequence "
                f"only has {n_colors} colors. Either use a different color sequence "
                "or set repeated_colors=True"
            )
        
        # Create color mapping for this column
        column_colors = {}
        for i, value in enumerate(unique_values):
            color_idx = i % n_colors
            column_colors[value] = colors[color_idx]
            
        color_mappings[col] = column_colors
    
    return color_mappings

def get_fastqs(wildcards):
    """
    Get FASTQ files from metadata sheet.
    """
    metadata = load_metadata(config["metadata"])
    fastq_cols = ["fq1", "fq2"]

    if config["fastq"]["auto"]:
        for i, col in enumerate(fastq_cols):
            metadata = metadata.assign(
                **{col: f"resources/reads/" + metadata["sample_id"] + f"_{i+1}.fastq.gz"}
            )
        metadata = metadata.set_index("sample_id")
    else:
        assert (
            "fq1" in metadata.columns
        ), f"The fq1 column in the metadata does not seem to exist. Please create one, or use the 'auto' option and name the fastq files as specified in the config/README.md"
        assert (
            "fq2" in metadata.columns
        ), f"The fq2 column in the metadata does not seem to exist. Please create one, or use the 'auto' option and name the fastq files as specified in the config/README.md"

        metadata = metadata.set_index("sample_id")

    u = metadata.loc[wildcards.sample, fastq_cols].dropna()
    return [u.fq1, u.fq2]


def AmpSeekerOutputs(wildcards):
    inputs = []

    inputs.extend(
        expand(
            [
                "results/alignments/{sample}.bam",
                "results/alignments/{sample}.bam.bai",
            ],
            sample=samples,
        ),
    )

    if config["bcl-convert"]:
        inputs.extend(
            [
            "results/qc/index-read-qc/I1.html", 
            "results/qc/index-read-qc/I2.html",
            "results/notebooks/run-statistics.ipynb",
            "docs/ampseeker-results/notebooks/run-statistics.ipynb",            ]
        )

    if plate_info:
        inputs.extend(
            [
                "results/notebooks/reads-per-well.ipynb",
                "docs/ampseeker-results/notebooks/reads-per-well.ipynb",
            ]
        )

    if large_sample_size:
        inputs.extend(
            expand(
                [
                    "results/vcfs/{call_type}/{dataset}.complete.merge_vcfs",
                    "results/vcfs/{call_type}/{dataset}.annot.vcf",
                ],
                dataset=config["dataset"],
                call_type=call_types,
            )
        )
    else:
        inputs.extend(
            expand(
                "results/vcfs/{call_type}/{dataset}.annot.vcf",
                dataset=config["dataset"],
                call_type=call_types,
            )
        )

    if config["quality-control"]["coverage"]:
        inputs.extend(
            expand(
                [
                    "results/coverage/{sample}.per-base.bed.gz",
                    "results/notebooks/coverage.ipynb",
                    "docs/ampseeker-results/notebooks/coverage.ipynb",
                ],
                sample=samples,
            )
        )

    if config["quality-control"]["fastp"]:
        inputs.extend(
            expand(
                [
                    "results/qc/fastp_reports/{sample}.html",
                    "results/notebooks/read-quality.ipynb",
                    "docs/ampseeker-results/notebooks/read-quality.ipynb",
                ],
                sample=samples,
            )
        )

    if config["quality-control"]["qualimap"]:
        inputs.extend(
            expand(
                [
                    "results/qc/qualimap/{sample}",
                ],
                sample=samples,
            )
        )

    if config["quality-control"]["multiqc"]:
        inputs.extend(
            expand(
                [
                    "results/qc/multiqc/multiqc_report.html",
                ],
            )
        )

    if config["quality-control"]["stats"]:
        inputs.extend(
            expand(
                [
                    "results/alignments/bamStats/{sample}.flagstat",
                    "results/qc/{dataset}.merged.vcf.txt",
                ],
                sample=samples,
                dataset=config["dataset"],
            )
        )

    if config["analysis"]["population-structure"]:
        inputs.extend(
            expand(
                [
                    "results/notebooks/population-structure.ipynb",
                    "docs/ampseeker-results/notebooks/population-structure.ipynb",
                ],
            )
        )

    if config["analysis"]["genetic-diversity"]:
        inputs.extend(
            expand(
                [
                    "results/notebooks/genetic-diversity.ipynb",
                    "docs/ampseeker-results/notebooks/genetic-diversity.ipynb",
                ],
            )
        )

    if config["analysis"]["allele-frequencies"]:
        inputs.extend(
            expand(
                [
                    "results/notebooks/allele-frequencies.ipynb",
                    "docs/ampseeker-results/notebooks/allele-frequencies.ipynb",
                ],
            )
        )

    if config["analysis"]["sample-map"]:
        inputs.extend(
            expand(
                [
                    "results/notebooks/sample-map.ipynb",
                    "docs/ampseeker-results/notebooks/sample-map.ipynb",
                ],
            )
        )

    if config["analysis"]["igv"]:
        inputs.extend(
            [
                "results/notebooks/IGV-explore.ipynb",
                "docs/ampseeker-results/notebooks/IGV-explore.ipynb",
            ]
        )

    if config["build-jupyter-book"]:
        inputs.extend(["results/ampseeker-results/_build/html/index.html"])

    # ag-vampir/species-id notebook
    if config["panel"] == "ag-vampir":
        inputs.extend(
            [
                "results/notebooks/ag-vampir/species-id.ipynb",
                "docs/ampseeker-results/notebooks/ag-vampir/species-id.ipynb",
                "results/notebooks/ag-vampir/kdr-origins.ipynb",
                "docs/ampseeker-results/notebooks/ag-vampir/kdr-origins.ipynb",
            ],
        )

    return inputs


def welcome(version):
    import datetime

    print("---------------------------- AmpSeeker ----------------------------")
    print(f"Running AmpSeeker snakemake workflow in {workflow.basedir}\n")
    print(f"Workflow Version: {version}")
    print("Execution time: ", datetime.datetime.now())
    print(f"Dataset: {config['dataset']}", "\n")
