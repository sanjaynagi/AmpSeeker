rule bcl_convert:
    input:
        sample_csv = config["illumina_dir"] + "SampleSheet.csv",
        illumina_in_dir = config["illumina_dir"]
    output: 
        directory ("resources/reads/")
        # prim_reads = expand("resources/reads/{sample}_S{sample_number}_L001_R{n}_001.fastq.gz", n=[1,2], sample=samples, sample_number=samples.str.extract('(\d+)'))
    log:
        "logs/bcl_convert.log"
    params:
        # output_reads = config["read_dir"],
    shell:
        "bcl-convert --bcl-input-directory {input.illumina_in_dir} --force --output-directory {output.output_reads} --sample-sheet {input.sample_csv} 2> {log}"

rule rename_fastq:
    input:
        # reads=read_dir + "*.fastq.gz"
        # reads=read_dir
        reads = "resources/reads/"
    output:
        # store = directory (expand("{reads}", reads=read_dir)),
        output_reads = expand("results/reads/{sample}_{n}.fastq.gz", n=[1,2], sample=samples)
       
    log:
        "logs/rename_fastq.log"
    params:
        # workflow_reads = "results/reads/",
        # metadata = pd.read_csv(config['metadata'], sep="\t"),
        # present_samples = metadata['sampleID']
    conda:
        "../envs/AmpSeq_cli.yaml"
    shell:
        """
        rename s/S[[:digit:]]\+_L001_R// {input.reads}/*.gz 
        rename s/_001// {input.reads}/*.gz 
        cp -r resources/reads/* results/reads/ 
        """
    #    cp {input.reads}/sample* {params.workflow_reads} 2>>{log}