rule jupyterbook:
    input:
        pages = "docs/ampseeker-book",
        igv = "docs/ampseeker-book/notebooks/IGV-explore.ipynb",
        coverage = "docs/ampseeker-book/notebooks/coverage.ipynb",
    output:
        directory("docs/ampseeker-book/_build/html/"),
        home_page = "docs/ampseeker-book/_build/html/index.html"
    log:
        "logs/jupyterbook/jupyterbook.log"
    conda:
        "../envs/AmpSeeker-jupyterbook.yaml"
    shell:
        """
        jupyter-book build --all {input.pages} &&
        ln -sf docs/ampseeker-book/_build/html/index.html AmpSeeker-results.html
        """