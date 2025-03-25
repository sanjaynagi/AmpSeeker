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
       r"""
        bcl2fastq --runfolder-dir {input.illumina_in_dir} --output-dir resources/reads/ --sample-sheet {input.sample_csv} 2> {log}
        rename 's/_S\d+_L001_//; s/R1_001/_1/' resources/reads/*R1_*fastq.gz 2>> {log}
        rename 's/_S\d+_L001_//; s/R2_001/_2/' resources/reads/*R2_*fastq.gz 2>> {log}
        rename 's/_S\d+_L001_//; s/I1_001/_I1/' resources/reads/*I1_*fastq.gz 2>> {log}
        rename 's/_S\d+_L001_//; s/I2_001/_I2/' resources/reads/*I2_*fastq.gz 2>> {log}
        """

