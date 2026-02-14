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

import pandas as pd

# Predefined color palettes
PALETTES = {
    'Plotly': ['#636EFA', '#EF553B', '#00CC96', '#AB63FA', '#FFA15A', '#19D3F3', '#FF6692', '#B6E880', '#FF97FF', '#FECB52'],
    'Safe': ['#88CCEE', '#CC6677', '#DDCC77', '#117733', '#332288', '#AA4499', '#44AA99', '#999933', '#882255', '#661100', '#6699CC', '#888888'],
    'Vivid': ['#E58606', '#5D69B1', '#52BCA3', '#99C945', '#CC61B0', '#24796C', '#DAA51B', '#2F8AC4', '#764E9F', '#ED645A', '#CC3A8E', '#A5AA99'],
    'Pastel': ['#B6B6F2', '#FFB6B6', '#BFFFBF', '#FFE2BF', '#BFE2FF', '#FFB6FF', '#E2FFB6', '#FFD9B6', '#B6FFE2', '#D9B6FF'],
    'Dark2': ['#1B9E77', '#D95F02', '#7570B3', '#E7298A', '#66A61E', '#E6AB02', '#A6761D', '#666666']
}

def create_color_mapping(df, columns, repeated_colors=True, unassigned_value='unassigned'):
    """
    Create color mappings using different palettes for each column.
    
    Args:
        df: Input dataframe
        columns: List of column names to create color mappings for
        repeated_colors: If True, colors will repeat if there are more categories than colors
        unassigned_value: Value that should be colored grey (default: 'unassigned')
            
    Returns:
        Dictionary of color mappings for each column
    """
    color_mappings = {}
    grey_color = "#808080"  # Standard grey color
    
    for i, col in enumerate(columns):
        # Cycle through palettes
        palette_name = list(PALETTES.keys())[i % len(PALETTES)]
        colors = PALETTES[palette_name]
        
        unique_values = df[col].unique()
        n_categories = len(unique_values)
        n_colors = len(colors)
        
        # Create color mapping for this column
        column_colors = {}
        for j, value in enumerate(unique_values):
            if value == unassigned_value:
                column_colors[value] = grey_color
            else:
                color_idx = j % n_colors if repeated_colors else j
                column_colors[value] = colors[color_idx]
            
        color_mappings[col] = column_colors
    
    return color_mappings


def get_fastqs(wildcards, platform):
    """
    Get FASTQ files from metadata sheet.
    """
    metadata = load_metadata("results/config/metadata.tsv")

    if platform == "illumina":
        fastq_cols = ["fq1", "fq2"]
    elif platform in ["nanopore", "ont"]:
        fastq_cols = ["fq1"]

    if fastq_auto:
        for i, col in enumerate(fastq_cols):
            metadata = metadata.assign(
                **{col: f"resources/reads/" + metadata["sample_id"] + f"_{i+1}.fastq.gz"}
            )
        metadata = metadata.set_index("sample_id")
    else:
        for col in fastq_cols:
            assert (
                col in metadata.columns
            ), f"The {col} column in the metadata does not seem to exist. Please create one, or use the 'auto' option."

        metadata = metadata.set_index("sample_id")

    u = metadata.loc[wildcards.sample, fastq_cols].dropna()
    return u.values


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

    if config["from-bcl"]:
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
                    "results/vcfs/{call_type}/.complete.{dataset}.merge_vcfs",
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
                    "results/qc/{dataset}-{call_type}.merged.vcf.txt",
                ],
                sample=samples,
                dataset=config["dataset"],
                call_type=call_types,
            )
        )

    if config["analysis"]["population-structure"]:
        inputs.extend(
                [
                    "results/notebooks/population-structure.ipynb",
                    "docs/ampseeker-results/notebooks/population-structure.ipynb",
                ],
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

    if config["build-jupyter-book"]:
        inputs.extend(["results/ampseeker-results/_build/html/index.html"])

    # ag-vampir/species-id notebook
    if config["panel"] == "ag-vampir":
        inputs.extend(
            [
                "results/notebooks/ag-vampir/species-id.ipynb",
                "docs/ampseeker-results/notebooks/ag-vampir/species-id.ipynb",
                "results/notebooks/ag-vampir/kdr-analysis.ipynb",
                "docs/ampseeker-results/notebooks/ag-vampir/kdr-analysis.ipynb",
            ],
        )

    return inputs


def welcome(version, using_user_metadata_colours):
    import datetime

    # Detect rich availability once
    try:
        from rich.console import Console
        from rich.panel import Panel
        from rich.table import Table
        use_rich = True
    except ImportError:
        use_rich = False

    now = datetime.datetime.now().replace(microsecond=0)

    # Determine input description
    if config['platform'] == 'illumina' and config["from-bcl"]:
        input_str = f"Illumina BCL folders: {', '.join(illumina_dirs)}"
    elif config['platform'] == 'illumina' and not config["from-bcl"] and fastq_auto:
        input_str = "FASTQs auto-detected in resources/reads/"
    elif config['platform'] == 'illumina':
        input_str = "FASTQ paths provided in metadata (fq1/fq2)"
    else:
        input_str = "Nanopore FASTQ paths provided in metadata (fq1)"

    if use_rich:
        console = Console()

        table = Table.grid(padding=(0, 2))
        table.add_row("Working dir:", workflow.basedir)
        table.add_row("Workflow version:", version)
        table.add_row("Execution time:", str(now))
        table.add_row("Dataset:", config["dataset"])
        table.add_row("Platform:", config["platform"])
        table.add_row("Input:", input_str)
        if using_user_metadata_colours:
            table.add_row("Metadata colours: user-provided schema in config/metadata_colours.json")

        console.print(
            Panel(
                table,
                title="[bold cyan3]AmpSeeker[/bold cyan3]",
                border_style="cyan3",
            )
        )

    else:
        # Plain fallback (no dependencies)
        sep = "-" * 64
        print(f"\n{sep}")
        print(" AmpSeeker ")
        print(f"{sep}")
        print(f"{'Working dir:':16} {workflow.basedir}")
        print(f"{'Workflow version:':16} {version}")
        print(f"{'Execution time:':16} {now}")
        print(f"{'Dataset:':16} {config['dataset']}")
        print(f"{'Platform:':16} {config['platform']}")
        print(f"{'Input:':16} {input_str}")
        if using_user_metadata_colours:
            print(f"{'Metadata colours:':16} user-provided schema in config/metadata_colours.json")
        print(f"{sep}\n")
