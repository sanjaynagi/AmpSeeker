rule jupyterbook:
    input:
        pages = "docs/ampseeker-results",
        igv = "docs/ampseeker-results/notebooks/IGV-explore.ipynb" if config["analysis"]["igv"] else [],
        coverage = "docs/ampseeker-results/notebooks/coverage.ipynb" if config["quality-control"]["coverage"] else [],
        pca = "docs/ampseeker-results/notebooks/principal-component-analysis.ipynb" if config["analysis"]["pca"] else [],
        af = "docs/ampseeker-results/notebooks/allele-frequencies.ipynb" if config["analysis"]["allele-frequencies"] else [],
        sample_map = "docs/ampseeker-results/notebooks/sample-map.ipynb" if config["analysis"]["sample-map"] else [],
    output:
        directory("results/ampseeker-results/_build/html/"),
        home_page = "results/ampseeker-results/_build/html/index.html"
    log:
        "logs/jupyterbook/jupyterbook.log"
    conda:
        "../envs/AmpSeeker-jupyterbook.yaml"
    shell:
        """
        jupyter-book build --all {input.pages} --path-output results/ampseeker-results &&
        ln -sf results/ampseeker-results/_build/html/index.html AmpSeeker-results.html
        """