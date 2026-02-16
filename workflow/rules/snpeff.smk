
rule snp_eff_db_download:
    """
    Download the snpEff database for your species
    """
    output:
        touch("results/vcfs/annotations/mysnpeffdb/.db.dl"),
    log:
        "logs/snpEff/snpEffDbDownload.log",
    conda:
        "../envs/AmpSeeker-snpeff.yaml"
    retries: 3
    params:
        ref=config["reference-snpeffdb"],
        dataDir="results/vcfs/annotations/mysnpeffdb",
    shell:
        "snpEff download {params.ref} -dataDir {params.dataDir} 2> {log}"


rule create_custom_snp_eff_db:
    """
    Create a custom SnpEff database from a reference genome and GFF file
    """
    input:
        fa=config['reference-fasta'],
        gff=config['reference-gff3'],
    output:
        "results/vcfs/annotations/mysnpeffdb/sequences.fa",
        "results/vcfs/annotations/mysnpeffdb/genes.gff",
        "results/vcfs/annotations/mysnpeffdb/snpEffectPredictor.bin"
    log:
        "logs/snpEff/createCustomSnpEffDb.log",
    conda:
        "../envs/AmpSeeker-snpeff.yaml"
    params:
        dataDir=lambda x: wkdir + "/results/vcfs/annotations/",
        wkdir=wkdir
    shell:
        """
        ln -s {params.wkdir}/{input.fa} {params.dataDir}/mysnpeffdb/sequences.fa 2> {log}
        ln -s {params.wkdir}/{input.gff} {params.dataDir}/mysnpeffdb/genes.gff 2> {log}
        snpEff build -gff3 -v -dataDir {params.dataDir} -configOption mysnpeffdb.genome=mysnpeffdb mysnpeffdb -noCheckCds -noCheckProtein 2>> {log}
        """


rule snp_eff:
    """
    Run snpEff on the VCFs 
    """
    input:
        calls="results/vcfs/{call_type}/{dataset}.merged.vcf",
        db="results/vcfs/annotations/mysnpeffdb/snpEffectPredictor.bin" if config['custom-snpeffdb'] else "results/vcfs/annotations/mysnpeffdb/.db.dl",
    output:
        calls="results/vcfs/{call_type}/{dataset}.annot.vcf",
        stats="results/vcfs/{call_type}/{dataset}.summary.html",
        csvStats="results/vcfs/{call_type}/{dataset}.summary.csv",
    log:
        "logs/snpEff/{dataset}_{call_type}.log",
    conda:
        "../envs/AmpSeeker-snpeff.yaml"
    retries:
        3
    params:
        db=config["reference-snpeffdb"] if not config['custom-snpeffdb'] else "mysnpeffdb",
        prefix=lambda w, output: os.path.splitext(output[0])[0],
        dataDir=lambda x: wkdir + "/results/vcfs/annotations/",
    shell:
        """
        snpEff eff {params.db} -dataDir {params.dataDir} -configOption mysnpeffdb.genome=mysnpeffdb -stats {output.stats} -csvStats {output.csvStats} -ud 0 {input.calls} > {output.calls} 2> {log}
        """
