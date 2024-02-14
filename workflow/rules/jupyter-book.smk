rule jupyterbook:
    input:
        toc="docs/ampseeker-results/_toc.yml",
        toc_update="results/.toc.update.complete",
        pages="docs/ampseeker-results",
        run_info="docs/ampseeker-results/notebooks/run-information.ipynb",
        run_statistics="docs/ampseeker-results/notebooks/run-statistics.ipynb" if config['bcl-convert'] else [],
        process_notebooks="results/notebooks/process-notebooks.ipynb",
        snp_df="docs/ampseeker-results/notebooks/snp-dataframe.ipynb",
        sample_quality_control="docs/ampseeker-results/notebooks/sample-quality-control.ipynb",
        igv=(
            "docs/ampseeker-results/notebooks/IGV-explore.ipynb"
            if config["analysis"]["igv"]
            else []
        ),
        coverage=(
            "docs/ampseeker-results/notebooks/coverage.ipynb"
            if config["quality-control"]["coverage"]
            else []
        ),
        pca=(
            "docs/ampseeker-results/notebooks/principal-component-analysis.ipynb"
            if config["analysis"]["pca"]
            else []
        ),
        af=(
            "docs/ampseeker-results/notebooks/allele-frequencies.ipynb"
            if config["analysis"]["allele-frequencies"]
            else []
        ),
        sample_map=(
            "docs/ampseeker-results/notebooks/sample-map.ipynb"
            if config["analysis"]["sample-map"]
            else []
        ),
        read_quality=(
            "docs/ampseeker-results/notebooks/read-quality.ipynb"
            if config["quality-control"]["fastp"]
            else []
        ),
        reads_per_well=(
            "docs/ampseeker-results/notebooks/reads-per-well.ipynb"
            if plate_info
            else []
        ),
    output:
        directory("results/ampseeker-results/_build/html/"),
        home_page="results/ampseeker-results/_build/html/index.html",
    log:
        "logs/jupyterbook/jupyterbook.log",
    conda:
        "../envs/AmpSeeker-jupyterbook.yaml"
    shell:
        """
        jupyter-book build --all {input.pages} --path-output results/ampseeker-results &&
        ln -sf results/ampseeker-results/_build/html/index.html AmpSeeker-results.html
        """

rule process_toc:
    input:
        input_nb=f"{workflow.basedir}/notebooks/process-toc.ipynb",
        toc="docs/ampseeker-results/_toc.yml",
        kernel="results/.kernel.set",
    output:
        out_nb="results/notebooks/process-toc.ipynb",
        toc_complete=touch("results/.toc.update.complete"),
    conda:
        f"{workflow.basedir}/envs/AmpSeeker-python.yaml"
    log:
        "logs/notebooks/process-toc.log",
    params:
        wkdir=wkdir,
        bcl_convert = config['bcl-convert'],
        panel=config["panel"],
    shell:
        """
        papermill -k AmpSeq_python {input.input_nb} {output.out_nb} -p wkdir {params.wkdir} -p bcl_convert {params.bcl_convert} -p panel {params.panel} 2> {log}
        """


rule process_notebooks:
    input:
        input_nb=f"{workflow.basedir}/notebooks/process-notebooks.ipynb",
        snp_df="docs/ampseeker-results/notebooks/snp-dataframe.ipynb",
        igv=(
            "docs/ampseeker-results/notebooks/IGV-explore.ipynb"
            if config["analysis"]["igv"]
            else []
        ),
        coverage=(
            "docs/ampseeker-results/notebooks/coverage.ipynb"
            if config["quality-control"]["coverage"]
            else []
        ),
        pca=(
            "docs/ampseeker-results/notebooks/principal-component-analysis.ipynb"
            if config["analysis"]["pca"]
            else []
        ),
        af=(
            "docs/ampseeker-results/notebooks/allele-frequencies.ipynb"
            if config["analysis"]["allele-frequencies"]
            else []
        ),
        sample_map=(
            "docs/ampseeker-results/notebooks/sample-map.ipynb"
            if config["analysis"]["sample-map"]
            else []
        ),
        genetic_diversity=(
            "docs/ampseeker-results/notebooks/genetic-diversity.ipynb"
            if config["analysis"]["genetic-diversity"]
            else []
        ),
        read_quality=(
            "docs/ampseeker-results/notebooks/read-quality.ipynb"
            if config["quality-control"]["fastp"]
            else []
        ),
        reads_per_well=(
            "docs/ampseeker-results/notebooks/reads-per-well.ipynb"
            if plate_info
            else []
        ),
    output:
        out_nb="results/notebooks/process-notebooks.ipynb",
    conda:
        f"{workflow.basedir}/envs/AmpSeeker-python.yaml"
    log:
        "logs/notebooks/process_notebooks.log",
    params:
        wkdir=wkdir,
        panel=config["panel"],
    shell:
        """
        papermill -k AmpSeq_python {input.input_nb} {output.out_nb} -p wkdir {params.wkdir} -p panel {params.panel} 2> {log}
        """
