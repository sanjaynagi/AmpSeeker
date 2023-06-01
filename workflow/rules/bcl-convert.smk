rule bcl_convert:
    input:
        sample_csv = config["illumina-dir"] + "SampleSheet.csv",
        illumina_in_dir = config["illumina-dir"]
    output: 
        temp(directory("resources/reads/"))
    log:
        "logs/bcl_convert.log"
    shell:
        "bcl-convert --bcl-input-directory {input.illumina_in_dir} --force --output-directory {output} --sample-sheet {input.sample_csv} 2> {log}"

rule rename_fastq:
    input:
        reads = "resources/reads/"
    output:
        output_reads = expand("results/reads/{sample}_{n}.fastq.gz", n=[1,2], sample=samples)
    log:
        "logs/rename_fastq.log"
    conda:
        "../envs/AmpSeeker-cli.yaml"
    shell:
        """
        rename s/S[[:digit:]]\+_L001_R// {input.reads}/*.gz 
        rename s/_001// {input.reads}/*.gz 
        cp -r resources/reads/* results/reads/ 
        """
