rule jupyterbook:
    input:
        page = f"{workflow.basedir}/../docs/ampseeker-book",
        alignment = f"{workflow.basedir}/../docs/ampseeker-book/notebooks/IGV-explore.ipynb",
    output:
        directory("docs/ampseeker-book/_build/html/"),
        home_page = "docs/ampseeker-book/_build/html/index.html"
    log:
        "logs/jupyterbook/jupyterbook.log"
    conda:
        "../envs/AmpSeeker-jupyterbook.yaml"
    shell:
        """
        jupyter-book build --all {input.page}
        """