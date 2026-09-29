#!/usr/bin/env python3
"""Sum Salmon NumReads per gene over all libraries (complements the TPM matrix).

TPM >= 0.5 is a threshold; NumReads says whether *any* read was assigned to the gene at all.
Output columns: salmon_reads_total, salmon_n_samples_ge10reads.

Usage: salmon_gene_numreads.py --quants-dir <...>/expression_validation/quants --out salmon_numreads.tsv
"""
import argparse
import glob
import os

import pandas as pd


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--quants-dir", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    tot = nsamp = None
    files = sorted(glob.glob(os.path.join(a.quants_dir, "*_quant", "quant.sf")))
    for f in files:
        q = pd.read_csv(f, sep="\t", index_col=0, usecols=["Name", "NumReads"])["NumReads"]
        q.index = q.index.str.replace(r"_t\d+$", "", regex=True)
        q = q.groupby(level=0).sum()
        tot = q if tot is None else tot.add(q, fill_value=0)
        s = (q >= 10).astype(int)
        nsamp = s if nsamp is None else nsamp.add(s, fill_value=0)
    pd.DataFrame({"salmon_reads_total": tot, "salmon_n_samples_ge10reads": nsamp}).to_csv(a.out, sep="\t")
    print("%d libraries -> %s" % (len(files), a.out))


if __name__ == "__main__":
    main()
