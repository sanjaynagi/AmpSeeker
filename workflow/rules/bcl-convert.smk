import os


def illumina_dir_from_run_id(run_id):
    return illumina_dirs[int(run_id)]


rule bcl_convert_run:
    input:
        sample_csv=lambda w: os.path.join(illumina_dir_from_run_id(w.run_id), "SampleSheet.csv"),
        illumina_in_dir=lambda w: illumina_dir_from_run_id(w.run_id),
    output:
        done=touch("resources/reads/bcl_runs/{run_id}/.complete"),
        stats_json="resources/reads/bcl_runs/{run_id}/Stats/Stats.json",
        demultiplex_stats="resources/reads/bcl_runs/{run_id}/Stats/DemultiplexingStats.xml",
    conda:
        "../envs/AmpSeeker-bcl2fastq.yaml"
    log:
        "logs/bcl2fastq/{run_id}.log",
    params:
        outdir="resources/reads/bcl_runs/{run_id}",
    shell:
        """
        rm -rf {params.outdir}
        mkdir -p {params.outdir}

        bcl2fastq \
          --runfolder-dir {input.illumina_in_dir} \
          --output-dir {params.outdir} \
          --sample-sheet {input.sample_csv} 2> {log}
        """


rule bcl_convert:
    input:
        per_run_done=expand("resources/reads/bcl_runs/{run_id}/.complete", run_id=illumina_run_ids),
        per_run_stats=expand("resources/reads/bcl_runs/{run_id}/Stats/Stats.json", run_id=illumina_run_ids),
        per_run_demultiplex=expand(
            "resources/reads/bcl_runs/{run_id}/Stats/DemultiplexingStats.xml",
            run_id=illumina_run_ids,
        ),
    output:
        output_reads=expand("resources/reads/{sample}_{n}.fastq.gz", n=[1, 2], sample=samples),
        index1_reads="resources/reads/I1_combined.fastq.gz",
        index2_reads="resources/reads/I2_combined.fastq.gz",
        stats_json="resources/reads/Stats/Stats.json",
        demultiplex_stats="resources/reads/Stats/DemultiplexingStats.xml",
    run:
        import glob
        import json
        import shutil

        os.makedirs("resources/reads/Stats", exist_ok=True)

        for sample in samples:
            sample = str(sample)
            r1_files = sorted(glob.glob(f"resources/reads/bcl_runs/*/{sample}_*_R1_001.fastq.gz"))
            r2_files = sorted(glob.glob(f"resources/reads/bcl_runs/*/{sample}_*_R2_001.fastq.gz"))

            if len(r1_files) == 0 or len(r2_files) == 0:
                raise ValueError(f"Missing one or more run FASTQs for sample: {sample}")

            with open(f"resources/reads/{sample}_1.fastq.gz", "wb") as out_fh:
                for f in r1_files:
                    with open(f, "rb") as in_fh:
                        shutil.copyfileobj(in_fh, out_fh)

            with open(f"resources/reads/{sample}_2.fastq.gz", "wb") as out_fh:
                for f in r2_files:
                    with open(f, "rb") as in_fh:
                        shutil.copyfileobj(in_fh, out_fh)

        i1_files = sorted(glob.glob("resources/reads/bcl_runs/*/*_I1_001.fastq.gz"))
        i2_files = sorted(glob.glob("resources/reads/bcl_runs/*/*_I2_001.fastq.gz"))
        if len(i1_files) > 0:
            with open(output.index1_reads, "wb") as out_fh:
                for f in i1_files:
                    with open(f, "rb") as in_fh:
                        shutil.copyfileobj(in_fh, out_fh)
        if len(i2_files) > 0:
            with open(output.index2_reads, "wb") as out_fh:
                for f in i2_files:
                    with open(f, "rb") as in_fh:
                        shutil.copyfileobj(in_fh, out_fh)

        shutil.copyfile(input.per_run_demultiplex[0], output.demultiplex_stats)

        merged_stats = {"ConversionResults": []}
        for stats_path in input.per_run_stats:
            with open(stats_path) as fh:
                data = json.load(fh)
            merged_stats["ConversionResults"].extend(data.get("ConversionResults", []))

        with open(output.stats_json, "w") as fh:
            json.dump(merged_stats, fh)
