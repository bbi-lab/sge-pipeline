import pandas as pd
import numpy as np
import pysam
import glob

import sge_counts
import sge_target



def calcPairwisePearsonR(target, snvfile1, snvfile2):
    '''calculates the pairwise Pearson r between two named sample ids
    
    '''
    s1counts = sge_counts.getSNVCounts(snvfile1, augment=True, pseudocount=0)
    s2counts = sge_counts.getSNVCounts(snvfile2, augment=True, pseudocount=0)
    s1df = s1counts.melt(id_vars=["chrom", "pos", "target", "repl", "day"],
                     value_vars=["A", "C", "G", "T"],
                     var_name="mutant_allele", value_name="count_s1")
    s2df = s2counts.melt(id_vars=["chrom", "pos", "target", "repl", "day"],
                     value_vars=["A", "C", "G", "T"],
                     var_name="mutant_allele", value_name="count_s2")
    df = s1df.merge(s2df, on=["chrom", "pos", "mutant_allele"])
    df = df[(df["pos"] >= target.editstartpos) & (df["pos"] <= target.editendpos)]
    df = df.merge(target.refdf, on="pos")
    df = df[df["ref"] != df["mutant_allele"]]   

    # compute SNV correlation
    mycorr = df[["count_s1", "count_s2"]].corr().iloc[0,1]
    return mycorr


# we need the reference sequence here so that we can eliminate the ref positions
# and the skip positions from the dataframe so that they don't impact the correlation
def calcMeanPearsonR(targetfile, targetname, countsdir):
    '''calculates the mean Pearson r between all samples for a specific
    target, stratified by timepoint (day)

    if only a single observation/replicate is availble for a given timepoint, 
    no value is returned for that timepoint

    returns: dictionary of (string) day-->(float) correlation
    '''
    mean_corrs = {}
    target = sge_target.Target(targetname, targetfile)
    samples = target.getSamplelist(countsdir, include_neg=False)
    for day, repllist in samples.items():
        if len(repllist) > 1:
            corrs = []
            for index, repl1 in enumerate(repllist):
                for repl2 in repllist[index+1:]:
                    mycorr = calcPairwisePearsonR(target, repl1, repl2)
                    corrs.append(mycorr)
            mean_corrs[day] = np.mean(corrs)
    return mean_corrs



def getVEPdf(vepfile):
    '''reads the output of Variant Effect Predictor files, converts to 
    pandas df, and returns it

    '''
    vepdf = pd.read_csv(vepfile, sep="\t", header=10)
    vepdf[["chrom", "pos"]] = vepdf["Location"].str.split(":", expand=True)
    vepdf["pos"] = vepdf["pos"].astype(int)
    return vepdf
