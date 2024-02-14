rule bcl_convert:
    input:
        sample_csv=config["illumina-dir"] + "SampleSheet.csv",
        illumina_in_dir=config["illumina-dir"],
    output:
        reads_dir=directory("resources/bcl_output"),
        fastq_list="resources/bcl_output/Reports/fastq_list.csv",
        demultiplex_stats = "resources/bcl_output/Reports/Demultiplex_Stats.csv",
    singularity:
        "docker://nfcore/bclconvert"
    log:
        "logs/bcl_convert.log",
    shell:
        "bcl-convert --bcl-input-directory {input.illumina_in_dir} --output-directory {output.reads_dir} --sample-sheet {input.sample_csv} --force 2> {log}"


rule rename_fastq:
    """
    If users demultiplex from BCL than rename, otherwise symlink to results
    """
    input:
        reads="resources/bcl_output/",
        read_dir=rules.bcl_convert.output,
        fastq_list="resources/bcl_output/Reports/fastq_list.csv",
    output:
        output_reads=expand(
            "resources/reads/{sample}_{n}.fastq.gz", n=[1, 2], sample=samples
        ),
    log:
        "logs/rename_fastq.log",
    conda:
        "../envs/AmpSeeker-cli.yaml"
    params:
        wd=wkdir,
    shell:
        """
        while IFS="," read -r _ sample_id _ _ read1 read2; do
            
            if [[ $sample_id == "RGSM" ]]; then
               continue
            fi
            
            echo renaming $read1 and $read2 to ${{sample_id}}_1.fastq.gz and ${{sample_id}}_2.fastq.gz 
            mv $read1 resources/reads/${{sample_id}}_1.fastq.gz
            mv $read2 resources/reads/${{sample_id}}_2.fastq.gz
        done < {input.fastq_list}
        """
