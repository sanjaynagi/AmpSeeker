"""Check that masked clair3 gVCFs merge so no-data samples stay missing, not wild type.

Mirrors the bcftools commands in mask_clair3_gvcf and bcftools_merge (workflow/rules).
"""
import shutil
import subprocess

import pytest

pytestmark = pytest.mark.skipif(
    shutil.which("bcftools") is None or shutil.which("samtools") is None,
    reason="bcftools and samtools are required",
)

REF = "ACGT" * 100  # 400 bp; position 100 is T
HEADER = """##fileformat=VCFv4.2
##FILTER=<ID=PASS,Description="p">
##ALT=<ID=NON_REF,Description="x">
##INFO=<ID=END,Number=1,Type=Integer,Description="end">
##FORMAT=<ID=GT,Number=1,Type=String,Description="g">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="q">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="d">
##FORMAT=<ID=AD,Number=R,Type=Integer,Description="a">
##FORMAT=<ID=MIN_DP,Number=1,Type=Integer,Description="m">
##FORMAT=<ID=PL,Number=G,Type=Integer,Description="p">
##contig=<ID=c1,length=400>
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{sample}
"""


def block(start, end, gt, gq, min_dp):
    ref = REF[start - 1]
    return f"c1\t{start}\t.\t{ref}\t<NON_REF>\t0\t.\tEND={end}\tGT:GQ:MIN_DP:PL\t{gt}:{gq}:{min_dp}:0,60,600\n"


def variant(pos, gt, dp):
    return (
        f"c1\t{pos}\t.\t{REF[pos - 1]}\tA,<NON_REF>\t60\tPASS\t.\tGT:GQ:DP:AD:PL\t"
        f"{gt}:60:{dp}:0,{dp},0:81,80,0,990,990,990\n"
    )


def run(cmd, **kw):
    return subprocess.run(cmd, shell=True, check=True, capture_output=True, text=True, **kw).stdout


def test_masked_gvcf_merge_keeps_no_data_samples_missing(tmp_path):
    ref = tmp_path / "ref.fa"
    ref.write_text(">c1\n" + REF + "\n")
    run(f"samtools faidx {ref}")

    gvcfs = {
        # alt homozygous at 100
        "s1": block(1, 99, "0/0", 50, 30) + variant(100, "1/1", 20) + block(101, 400, "0/0", 50, 30),
        # covered, wild type
        "s2": block(1, 400, "0/0", 50, 30),
        # clair3 labels a zero-depth block 0/0: must be masked
        "s3": block(1, 90, "0/0", 50, 30) + block(91, 120, "0/0", 1, 0) + block(121, 400, "0/0", 50, 30),
        # no record at all across 91-120
        "s4": block(1, 90, "0/0", 50, 30) + block(121, 400, "0/0", 50, 30),
    }
    masked = []
    for name, body in gvcfs.items():
        raw = tmp_path / f"{name}.g.vcf"
        raw.write_text(HEADER.format(sample=name) + body)
        out = tmp_path / f"{name}.calls.vcf.gz"
        run(
            f"bcftools +setGT {raw} -Oz -o {out} -- -t q -n . "
            "-i '(FMT/MIN_DP<10) | (FMT/DP<10)'"
        )
        run(f"bcftools index -t {out}")
        masked.append(str(out))

    merged = run(
        f"bcftools merge -g {ref} -Ou {' '.join(masked)} | "
        "bcftools view -a -i 'ALT!=\"<NON_REF>\"' -H"
    ).splitlines()

    assert len(merged) == 1
    fields = merged[0].split("\t")
    assert fields[1] == "100" and fields[4] == "A"  # <NON_REF> removed
    gts = [f.split(":")[0] for f in fields[9:]]
    assert gts == ["1/1", "0/0", "./.", "./."]
