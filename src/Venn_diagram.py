"""input:  2 or more gene_exp.diff from Cuffdiff
return:  VennDiagram
two different types: python builtin; one from github(import pyvenn)
"""

import pandas as pd
import numpy as np
from matplotlib_venn import venn2, venn3
from matplotlib import pyplot as plt
import sys
import pyvenn
pd.set_option('display.max_columns', 50)



def draw_venn2(list1, list2, name1="", name2="", title=""):
    v = venn2([set(list1), set(list2)], set_labels=(name1, name2))
    v.get_patch_by_id('100').set_color('blue')
    v.get_patch_by_id('110').set_color('grey')
    v.get_patch_by_id('010').set_color('lightgreen')
    plt.title(title+'\n')
    pass


def draw_venn3(lst1, lst2, lst3, name1="", name2="", name3="", title=""):
    labels = pyvenn.get_labels([lst1, lst2, lst3], fill=['number', 'logic'])
    fig, ax = pyvenn.venn3(labels, names=[name1, name2, name3])
    fig.savefig("../results/" + title, bbox_inches='tight')  # can save files
    plt.close()
    pass


def draw_venn4(lst1, lst2, lst3, lst4, name1="", name2="", name3="", name4="", title="venn4.png"):
    labels = pyvenn.get_labels([lst1,lst2,lst3,lst4], fill=['number', 'logic'])
    fig, ax = pyvenn.venn4(labels, names=[name1, name2, name3, name4])
    fig.savefig("../results/"+title, bbox_inches='tight')  # can save files
    plt.close()
    pass



def getList_InRange(df, cutoff1=1, cutoff2=-1):
    """
    :param df: gene.diff file convert df from cuffdiff
    :param cutoff1: upregulated 2 fold
    :param cutoff2: downregulated 2 fold
    :return: genelist
    """
    df = df[~df["log2(fold_change)"].isin(["#NAME?"])]
    df["log2(fold_change)"] = df["log2(fold_change)"].apply(lambda x: np.float(x))
    df_c = df [(df["log2(fold_change)"]>cutoff1) | (df["log2(fold_change)"]<cutoff2)]
    df_list = df_c["gene_id"].tolist()
    return df_list


def geneList_up(df, cutoff1):
    df = df[~df["log2(fold_change)"].isin(["#NAME?"])]
    df["log2(fold_change)"] = df["log2(fold_change)"].apply(lambda x: np.float(x))
    df_c = df [(df["log2(fold_change)"]>cutoff1)]
    df_list = df_c["gene_id"].tolist()
    return df_list


def geneList_down(df, cutoff2):
    df = df[~df["log2(fold_change)"].isin(["#NAME?"])]
    df["log2(fold_change)"] = df["log2(fold_change)"].apply(lambda x: np.float(x))
    df_c = df [(df["log2(fold_change)"]<cutoff2)]
    df_list = df_c["gene_id"].tolist()
    return df_list


def genelist_Pvalue(df, p_value=0.05):
    df_c = df [(df["p_value"]<=p_value)]
    df_list = df_c["gene_id"].tolist()
    return df_list


def main(file1, file2):
    df1 = pd.read_csv(file1, sep="\t")
    df2 = pd.read_csv(file2, sep="\t")
    lst1 = getList_InRange(df1)
    lst2 = getList_InRange(df2)
    draw_venn2(lst1, lst2)
    pass


if __name__ == '__main__':
    main()









