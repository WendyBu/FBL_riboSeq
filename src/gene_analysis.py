"""
up and down regulated gene pathway analysis.
generate venn diagram, heatmap, gsea enrichment, kegg pathway, NCBI-query important genes
generate snoRNA?
"""


import pandas as pd
import numpy as np
import os, sys, glob, os.path
import Venn_diagram as vd
pd.set_option("display.max_column", 100)



def split_up_down(filename):
    df = pd.read_csv(filename, sep="\t", index_col=0)
    # remove all the 0 in ribo-seq and RNA-seq
    filter0 = (df.iloc[:, 2] != 0) & (df.iloc[:, 3] != 0) & (df.iloc[:, 4] != 0) & (df.iloc[:, 5] != 0)
    filter_df = df[filter0]
    up_df = filter_df[filter_df.iloc[:, -1] >= 1.5]
    down_df = filter_df[filter_df.iloc[:, -1] < 0.67]
    return up_df.index.tolist(), down_df.index.tolist()


def draw_venn_diagram(lst1, lst2,lst3, name1, name2, name3, title):
    vd.draw_venn3(lst1, lst2, lst3, name1, name2, name3, title=title)
    pass


def main():
    FBL_filename = sys.argv[1]
    FBL_up, FBL_down = split_up_down(FBL_filename)
    print "FBL:", len(FBL_up), len(FBL_down)
    EZH2sh1_filename = sys.argv[2]
    EZH2sh1_up, EZH2sh1_down = split_up_down(EZH2sh1_filename)
    print "EZH2 sh1:", len(EZH2sh1_up), len(EZH2sh1_down)
    EZH2sh2_filename = sys.argv[3]
    EZH2sh2_up, EZH2sh2_down = split_up_down(EZH2sh2_filename)
    print "EZH2 sh2:", len(EZH2sh2_up), len(EZH2sh2_down)
    draw_venn_diagram(FBL_up, EZH2sh1_up, EZH2sh2_up, 'FBL_up', 'EZH2sh1_up', 'EZH2sh2_up', 'Translation Efficiency Increased')
    draw_venn_diagram(FBL_down, EZH2sh1_down, EZH2sh2_down, 'FBL_down', 'EZH2sh1_down', 'EZH2sh2_down', 'Translation Efficiency Decreased')
    with open("../results/geneset.txt", 'w+') as f:
        f.writelines("FBL_UP:\n")
        f.writelines(" ".join(FBL_up))
        f.writelines("\n")
        f.writelines("FBL_DOWN:\n")
        f.writelines(" ".join(FBL_down))
        f.writelines("\n")
        f.writelines("EZH2sh1_UP:\n")
        f.writelines(" ".join(EZH2sh1_up))
        f.writelines("\n")
        f.writelines("EZH2sh1_DOWN:\n")
        f.writelines(" ".join(EZH2sh1_down[0:2999]))
        f.writelines("\n")
        f.writelines("EZH2sh2_UP:\n")
        f.writelines(" ".join(EZH2sh2_up[0:2999]))
        f.writelines("\n")
        f.writelines("EZH2sh2_DOWN:\n")
        f.writelines(" ".join(EZH2sh2_down[0:2999]))

    pass



if __name__ == "__main__":
    main()


# python gene_analysis.py ../data/FBL_translation_ratio.xls ../data/EZH2sh1_translation_ratio.xls ../data/EZH2sh2_translation_ratio.xls
#