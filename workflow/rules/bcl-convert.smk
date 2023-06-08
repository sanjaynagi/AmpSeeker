rule bcl_convert:
    input:
        sample_csv = config["illumina-dir"] + "SampleSheet.csv",
        illumina_in_dir = config["illumina-dir"]
    output: 
        reads_dir = temp(directory("resources/reads")),
        fastq_list = "resources/reads/Reports/fastq_list.csv"
    log:
        "logs/bcl_convert.log"
    shell:
        "bcl-convert --bcl-input-directory {input.illumina_in_dir} --force --output-directory {output.reads_dir} --sample-sheet {input.sample_csv} 2> {log}"


rule rename_fastq:
    """
    If users demultiplex from BCL
    """
    input:
        reads = "resources/reads/",
        fastq_list = "resources/reads/Reports/fastq_list.csv"
    output:
        output_reads = expand("results/reads/{sample}_{n}.fastq.gz", n=[1,2], sample=samples)
    log:
        "logs/rename_fastq.log"
    conda:
        "../envs/AmpSeeker-cli.yaml"
    shell:
        """
        while IFS="," read -r _ sample_id _ _ read1 read2; do
            echo renaming $read1 and $read2 to ${{sample_id}}_1.fastq.gz and ${{sample_id}}_2.fastq.gz 2>> {log}
            mv $read1 results/reads/${{sample_id}}_1.fastq.gz &&
            mv $read2 results/reads/${{sample_id}}_2.fastq.gz
        done < {input.fastq_list}
        """

rule symlink_fastq:
    """
    If reads are provided by user and bcl-convert not used
    """
    input:
        input_reads = "resources/reads/{sample}_{n}.fastq.gz"
    output:
        output_reads = "results/reads/{sample}_{n}.fastq.gz"
    log:
        "logs/symlink_fastq/{sample}_{n}.log"
    conda:
        "../envs/AmpSeeker-cli.yaml"
    shell:
        """
        ln -sf {input.input_reads} {output.output_reads}
        """