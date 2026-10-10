import pandas as pd

from captnseq_toolkit.sites2genes import compute_gene_insertions


def test_counts_insertions_inside_gene(tmp_path):
    features = tmp_path / "cds.bed"
    features.write_text("chr\t0\t1000\tgeneA\t0\t1\nchr\t2000\t3000\tgeneB\t0\t1\n")
    # 10 columns: gene BED6 (chr,start,end,name,score,strand) + insertion (chr,start,end,score)
    sample = tmp_path / "intersect.bed"
    sample.write_text(
        "chr\t0\t1000\tgeneA\t0\t1\tchr\t500\t501\t10\n"
        "chr\t0\t1000\tgeneA\t0\t1\tchr\t600\t601\t5\n"
    )
    out = tmp_path / "out.tsv"
    compute_gene_insertions(features, sample, 0, 100, out)
    df = pd.read_csv(out, sep="\t").set_index("gene_name")
    assert df.loc["geneA", "ins_count"] == 2
    assert df.loc["geneA", "read_count"] == 15
    assert df.loc["geneB", "ins_count"] == 0
