import os 

rule bcl_convert:
    input:
        sample_csv=os.path.join(config["illumina-dir"], "SampleSheet.csv"),
        illumina_in_dir=config["illumina-dir"],
    output:
        output_reads=expand(
            "resources/reads/{sample}_{n}.fastq.gz", n=[1, 2], sample=samples
        ),
        demultiplex_stats = "resources/reads/Stats/DemultiplexingStats.xml",
    conda:
        "../envs/AmpSeeker-bcl2fastq.yaml"
    log:
        "logs/bcl2fastq.log",
    shell:
        """
        bcl2fastq --runfolder-dir {input.illumina_in_dir} --output-dir resources/reads/ --sample-sheet {input.sample_csv} 2> {log}
        rename 's/L001_//; s/_001//; s/S\d+_R//; s/S\d+_//' resources/reads/*fastq.gz 2>> {log}
        """

