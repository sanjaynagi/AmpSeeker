import numpy as np
import pandas as pd

import ampseeker as amp


DATA_DIR = "tests/data"
VCF_PATH = f"{DATA_DIR}/test.annot.vcf"
METADATA_PATH = f"{DATA_DIR}/test.metadata.tsv"
PANEL_PATH = f"{DATA_DIR}/test.panel.bed"

EXPECTED = {
    "n_samples": 30,
    "core_n_variants": 95,
    "core_geno_shape": (95, 30, 2),
    "core_pos_min": 34264,
    "core_pos_max": 49438586,
    "core_filter_pass_true": 0,
    "core_indel_true": 2,
    "core_is_snp_true": 83,
    "core_ann_missing": 11,
    "bed_shape": (93, 7),
    "bed_contigs": ["2L", "2R", "3L", "3R", "X"],
    "illumina_shape": (93, 30, 2),
    "nanopore_shape": (83, 30, 2),
    "illumina_missing_shape": (89, 30, 2),
    "illumina_segregating_shape": (41, 30, 2),
    "illumina_maf005_shape": (41, 30, 2),
    "vcf_df_all_shape": (95, 36),
    "vcf_df_seg_shape": (41, 36),
    "vcf_df_multiallelic_rows": 16,
    "split_df_shape": (112, 36),
    "nalt_non_null": 3140,
    "nalt_counts": {0: 2534, 1: 104, 2: 502},
    "snp_df_shape": (110, 8),
    "snp_geno_shape": (93, 30, 2),
    "snp_ann_empty": 12,
    "snp_alt_nan": 10,
    "freq_shape": (110, 15),
    "freq_columns": ["frq_Gaoua", "frq_Kisumu", "frq_Tiefora"],
    "freq_min": 0.0,
    "freq_max": 1.0,
    "freq_first_label": "2L | AGAP004679 | 209536 | Val305Ile | A",
    "njt_shape": (30, 30),
    "njt_outliers": [],
}


def _variant_count(path):
    with open(path, "r") as f:
        return sum(1 for line in f if not line.startswith("#"))


def test_load_bed_has_expected_columns():
    df = amp.load_bed(PANEL_PATH)
    assert list(df.columns) == ["contig", "start", "end", "amplicon_id", "mutation", "ref", "alt"]
    assert df.shape == EXPECTED["bed_shape"]
    assert sorted(df.contig.unique().tolist()) == EXPECTED["bed_contigs"]
    assert df["mutation"].nunique() == EXPECTED["bed_shape"][0]


def test_read_ann_field_matches_variant_rows():
    anns = amp.read_ANN_field(VCF_PATH)
    assert len(anns) == _variant_count(VCF_PATH)
    assert isinstance(anns, np.ndarray)
    assert int(np.sum(pd.isna(anns))) == EXPECTED["core_ann_missing"]


def test_load_vcf_returns_aligned_outputs():
    metadata = pd.read_csv(METADATA_PATH, sep="\t")
    geno, pos, contig, md_out, ref, alt, ann = amp.load_variants(
        vcf_path=VCF_PATH, metadata=metadata, platform="illumina", filter_indel=True
    )

    n_variants = geno.shape[0]
    assert n_variants == len(pos) == len(contig) == len(ref) == len(alt) == len(ann)
    assert geno.shape == EXPECTED["illumina_shape"]
    assert md_out.shape == metadata.shape
    assert list(md_out["sample_id"]) == list(metadata["sample_id"])


def test_load_vcf_core_mode_returns_raw_arrays():
    core = amp.load_vcf(VCF_PATH)
    expected_keys = {
        "samples",
        "geno",
        "pos",
        "contig",
        "filter_pass",
        "ref",
        "alt",
        "ann",
        "indel",
        "is_snp",
    }
    assert expected_keys.issubset(set(core.keys()))
    assert core["geno"].shape[0] == len(core["pos"]) == len(core["contig"]) == len(core["ref"]) == len(core["alt"]) == len(core["ann"])
    assert core["geno"].shape == EXPECTED["core_geno_shape"]
    assert len(core["samples"]) == EXPECTED["n_samples"]
    assert core["samples"][0] == "G1"
    assert core["samples"][-1] == "KIS_3"
    assert int(core["pos"].min()) == EXPECTED["core_pos_min"]
    assert int(core["pos"].max()) == EXPECTED["core_pos_max"]
    assert int(np.sum(core["filter_pass"])) == EXPECTED["core_filter_pass_true"]
    assert int(np.sum(core["indel"])) == EXPECTED["core_indel_true"]
    assert int(np.sum(core["is_snp"])) == EXPECTED["core_is_snp_true"]


def test_vcf_to_df_contains_core_columns():
    df = amp.vcf_to_df(VCF_PATH, seg=False)
    df_seg = amp.vcf_to_df(VCF_PATH, seg=True)
    expected = {"CHROM", "POS", "FILTER_PASS", "REF", "ALT", "ANN"}
    assert expected.issubset(set(df.columns))
    assert df.shape == EXPECTED["vcf_df_all_shape"]
    assert df_seg.shape == EXPECTED["vcf_df_seg_shape"]


def test_split_rows_with_multiple_alleles_expands_rows():
    df = amp.vcf_to_df(VCF_PATH, seg=False)
    samples = df.columns[6:]
    assert int(df["ALT"].astype(str).str.contains(",").sum()) == EXPECTED["vcf_df_multiallelic_rows"]

    split_df = amp.split_rows_with_multiple_alleles(df, samples)
    assert split_df.shape == EXPECTED["split_df_shape"]
    assert int(split_df["ALT"].astype(str).str.contains(",").sum()) == 0


def test_convert_genotype_to_alt_allele_count_values():
    df = amp.vcf_to_df(VCF_PATH, seg=False)
    samples = df.columns[6:]
    split_df = amp.split_rows_with_multiple_alleles(df, samples)
    nalt_df = amp.convert_genotype_to_alt_allele_count(split_df, samples)

    vals = pd.unique(nalt_df.loc[:, samples].values.ravel())
    vals = [v for v in vals if not pd.isna(v)]
    assert all(v in (0, 1, 2) for v in vals)
    series_vals = pd.Series(nalt_df.loc[:, samples].values.ravel())
    assert int(series_vals.notna().sum()) == EXPECTED["nalt_non_null"]
    counts = series_vals.value_counts(dropna=True).sort_index().to_dict()
    assert counts == EXPECTED["nalt_counts"]


def test_vcf_to_snp_dataframe_basic_contract():
    metadata = pd.read_csv(METADATA_PATH, sep="\t")
    snp_df, geno = amp.vcf_to_snp_dataframe(VCF_PATH, metadata=metadata, platform="illumina")

    assert {"contig", "pos", "ref", "alt", "ann", "variant_index", "alt_index", "label"}.issubset(snp_df.columns)
    assert geno.shape == EXPECTED["snp_geno_shape"]
    assert snp_df.shape == EXPECTED["snp_df_shape"]
    assert snp_df["label"].nunique() == EXPECTED["snp_df_shape"][0]
    assert int((snp_df["ann"] == "").sum()) == EXPECTED["snp_ann_empty"]
    assert int(snp_df["alt"].isna().sum()) == EXPECTED["snp_alt_nan"]


def test_calculate_frequencies_cohort_has_frequency_columns():
    metadata = pd.read_csv(METADATA_PATH, sep="\t")
    snp_df, geno = amp.vcf_to_snp_dataframe(VCF_PATH, metadata=metadata, platform="illumina")
    out = amp.calculate_frequencies_cohort(
        snp_df=snp_df,
        metadata=metadata,
        geno=geno,
        cohort_col="location",
        af_filter=False,
        missense_filter=False,
    )

    freq_cols = [c for c in out.columns if c.startswith("frq_")]
    assert out.shape == EXPECTED["freq_shape"]
    assert sorted(freq_cols) == EXPECTED["freq_columns"]
    assert float(out[freq_cols].min().min()) == EXPECTED["freq_min"]
    assert float(out[freq_cols].max().max()) == EXPECTED["freq_max"]
    assert out.index[0] == EXPECTED["freq_first_label"]


def test_plot_allele_frequencies_returns_figure():
    df = pd.DataFrame({"frq_A": [0.1, 0.5], "frq_B": [0.2, 0.9]}, index=["v1", "v2"])
    fig = amp.plot_allele_frequencies(df=df, cohort_col="location")
    assert hasattr(fig, "to_plotly_json")


def test_plot_pca_returns_figure():
    pca_df = pd.DataFrame(
        {
            "sample_id": ["s1", "s2"],
            "location": ["A", "B"],
            "PC1": [0.1, -0.1],
            "PC2": [0.2, -0.2],
            "PC3": [0.3, -0.3],
        }
    )
    color_mapping = {"location": {"A": "#000000", "B": "#ffffff"}}
    fig = amp.plot_pca(
        pca_df=pca_df,
        colour_column="location",
        cohort_columns=["location"],
        dataset="unit",
        color_mapping=color_mapping,
    )
    assert hasattr(fig, "to_plotly_json")


def test_compute_njt_inputs_returns_square_distance_matrix():
    metadata = pd.read_csv(METADATA_PATH, sep="\t")
    geno, *_ = amp.load_variants(vcf_path=VCF_PATH, metadata=metadata, platform="illumina", filter_indel=True)
    dists, leaf_data, exclude_outliers = amp.compute_njt_inputs(geno, metadata, cohort_col="location")

    assert dists.ndim == 2
    assert dists.shape == EXPECTED["njt_shape"]
    assert leaf_data.shape[0] == dists.shape[0] == EXPECTED["n_samples"]
    assert isinstance(exclude_outliers, np.ndarray)
    assert exclude_outliers.tolist() == EXPECTED["njt_outliers"]
    assert int(np.isnan(dists).sum()) == 0


def test_load_variants_optional_filters_match_expected_shapes():
    metadata = pd.read_csv(METADATA_PATH, sep="\t")
    g_i, *_ = amp.load_variants(VCF_PATH, metadata=metadata, platform="illumina", filter_indel=True)
    g_n, *_ = amp.load_variants(VCF_PATH, metadata=metadata, platform="nanopore", filter_indel=True)
    g_miss, *_ = amp.load_variants(
        VCF_PATH, metadata=metadata, platform="illumina", filter_indel=True, filter_missing=0.4
    )
    g_seg, *_ = amp.load_variants(
        VCF_PATH, metadata=metadata, platform="illumina", filter_indel=True, segregating_only=True
    )
    g_maf, *_ = amp.load_variants(
        VCF_PATH, metadata=metadata, platform="illumina", filter_indel=True, filter_maf=0.05
    )

    assert g_i.shape == EXPECTED["illumina_shape"]
    assert g_n.shape == EXPECTED["nanopore_shape"]
    assert g_miss.shape == EXPECTED["illumina_missing_shape"]
    assert g_seg.shape == EXPECTED["illumina_segregating_shape"]
    assert g_maf.shape == EXPECTED["illumina_maf005_shape"]


def test_get_consensus_taxon_majority_rule():
    row = pd.Series({"aim_taxon": "gambiae", "pca_taxon": "gambiae", "tree_taxon": "coluzzii"})
    assert amp.get_consensus_taxon(row) == "gambiae"

    row2 = pd.Series({"aim_taxon": "gambiae", "pca_taxon": "coluzzii", "tree_taxon": "arabiensis"})
    assert amp.get_consensus_taxon(row2) == "unassigned"


def test_signif_rounding_behaviour():
    arr = np.array([1234.0, 0.01234])
    out = amp.signif(arr, 2)
    assert np.allclose(out, np.array([1200.0, 0.012]))


def test_get_single_gen_call_no_402():
    row = pd.Series({"kdr_F_origin": "wt,F1", "kdr_S_origin": "wt,S2"})
    call = amp.get_single_gen_call(row)
    assert call == "F1,S2"


def test_melt_gt_counts_shape_and_missing():
    gt_counts = np.array(
        [
            [[2, 0, 0], [1, 1, 0]],
            [[0, 0, 0], [0, 2, 0]],
        ],
        dtype=np.int64,
    )
    out = amp._melt_gt_counts(gt_counts)
    assert out.shape == (4, 2)
    assert np.isnan(out[2, 0])
