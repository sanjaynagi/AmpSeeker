rule jupyterbook:
    input:
        pages = "docs/ampseeker-results",
        igv = "docs/ampseeker-results/notebooks/IGV-explore.ipynb",
        coverage = "docs/ampseeker-results/notebooks/coverage.ipynb",
    output:
        directory("docs/ampseeker-results/_build/html/"),
        home_page = "docs/ampseeker-results/_build/html/index.html"
    log:
        "logs/jupyterbook/jupyterbook.log"
    conda:
        "../envs/AmpSeeker-jupyterbook.yaml"
    params:
        dataset = dataset
    shell:
        """
        jupyter-book build --all {input.pages} &&
        ln -sf docs/ampseeker-results/_build/html/index.html AmpSeeker-{params.dataset}-results.html
        """