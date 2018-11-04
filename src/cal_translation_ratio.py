import pandas as pd
import numpy as np
import os, sys, glob, os.path
pd.set_option("display.max_column", 100)




def main():
    ribo_file = "/archive2/tmhyxb9/FBL/cuffdiff_FBL/genes.fpkm_tracking"
    ribo_seq = pd.read_csv(ribo_file, sep="\t", index_col="gene_id")
    ribo_seq = ribo_seq.loc[:,["control_FPKM", "FBL_FPKM"]]
    RNAC_file = "/archive2/tmhyxb9/FBL/legancy/FBL_kd/genes.fpkm_tracking"
    RNA_seq = pd.read_csv(RNAC_file, sep="\t", index_col="gene_id")
    RNA_seq = RNA_seq.loc[:, ["C_FPKM", "T_FPKM"]]
    ribo_RNA_seq = ribo_seq.join(RNA_seq)

    ribo_RNA_seq["control_ratio"] = ribo_RNA_seq.control_FPKM/ribo_RNA_seq.C_FPKM
    ribo_RNA_seq["FBL_ratio"] = ribo_RNA_seq.FBL_FPKM / ribo_RNA_seq.T_FPKM

    ribo_RNA_seq["FBL_ratio_diff"] = ribo_RNA_seq.FBL_ratio / ribo_RNA_seq.control_ratio

    ribo_RNA_seq.to_csv("FBL_translation_ratio.xls", sep="\t")

    print ribo_seq.head()
    print RNA_seq.head()
    print ribo_RNA_seq.head()
    pass



if __name__ == "__main__":
    main()


## translation ratio = ribo-seq fpkm / RNA-seq fpkm