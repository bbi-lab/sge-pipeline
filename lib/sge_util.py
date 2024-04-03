import pandas as pd
import numpy as np
import pysam
import glob

import sge_counts

def getTargetEditRegion(targetfile, targetname):
    '''returns a 3-tuple (chrom, startpos, endpos) representing the edited region of 
    a specific target name
    
    coordinates are 1-based and inclusive
    '''
    targetdf = pd.read_csv(targetfile,
                           header=0, sep="\t",
                           dtype={'skip_pos': 'string'})
    chrom = targetdf.loc[targetdf["exonname"] == targetname, "chrom"].values[0]
    startpos = targetdf.loc[targetdf["exonname"] == targetname, "editstart"].values[0]
    endpos = targetdf.loc[targetdf["exonname"] == targetname, "editstop"].values[0]
    return (chrom, startpos, endpos)


def getTargetAmpliconRegion(targetfile, targetname):
    '''returns a 3-tuple (chrom, startpos, endpos) representing the coordinates of the 
    amplicon sequenced for a specific target name

    coordinates are 1-based and inclusive
    '''
    targetdf = pd.read_csv(targetfile,
                           header=0, sep="\t",
                           dtype={'skip_pos': 'string'})
    chrom = targetdf.loc[targetdf["exonname"] == targetname, "chrom"].values[0]
    startpos = targetdf.loc[targetdf["exonname"] == targetname, "ampstart"].values[0]
    endpos = targetdf.loc[targetdf["exonname"] == targetname, "ampstop"].values[0]

    return (chrom, startpos, endpos)


def getTargetRequiredEdits(targetfile, targetname):
    '''returns a list of any required SNV edits (integer pos) for a given target name

    '''
    targetdf = pd.read_csv(targetfile,
                           header=0, sep="\t",
                           dtype={'skip_pos': 'string'})
    always_edited = []
    editstring = targetdf.loc[targetdf["exonname"] == targetname, "required_edits"].values[0]
    tedits = editstring.split(",")
    for edit in tedits:
        pos = int(edit[:-1])
        always_edited.append(pos)
    return always_edited


def getTargetSkipPositions(targetfile, targetname):
    '''returns a list of any positions (integer values) to be skipped in the
    generation of counts/frequencies for a given target name
    
    '''
    targetdf = pd.read_csv(targetfile,
                           header=0, sep="\t",
                           dtype={'skip_pos': 'string'})
    skip_pos = []
    skipstring = targetdf.loc[targetdf["exonname"] == targetname, "skip_pos"].values[0]
    if not pd.isna(skipstring):
        skip_parts = skipstring.split(",")
        for p in skip_parts:
            skip_pos.append(int(p))
    return skip_pos


def getTargetSampleList(targetname, countsdir, include_neg=False):
    '''performs a directory lookup of all files matching a specific target name in a 
    given counts directory, and returns a dictionary of filenames:
    (string) day --> [ (string) full_filename1, (string) full_filename2, ... ]

    if include_neg is True, the negative control sample is included; otherwise it is
    excluded from the list of samples

    '''
    samples = {}
    if countsdir.endswith("/"):
        countsdir = countsdir[:-1]
    for fullfn in glob.glob("%s/%s_*.snvs.tsv" % (countsdir, targetname)):
        fn = fullfn.split("/")[-1]
        parts = fn.split(".")[0].split("_")
        repl = parts[2]
        if not include_neg:
            if repl == 'NC': # throw out the negative control
                continue
        day = parts[3]
        if day not in samples:
            samples[day] = []
        samples[day].append(fullfn)
    return samples


# read reference sequence, return dataframe
def getReferenceSequence(chrom, startpos, endpos, 
                         refpath="/net/shendure/vol10/nobackup/genome/fasta/Homo_sapiens.GRCh38.dna.toplevel.fa"):
    '''extracts the reference sequence within the given 1-based start and 1-based 
    end coordinates, inclusive, and returns the sequence as a pandas dataframe with
    columns "pos" (integer position) and "ref" (reference base)
    
    '''
    regionstring = chrom + ":" + str(startpos) + "-" + str(endpos)
    with pysam.FastaFile(refpath) as fafh:
        if regionstring.startswith("chr"):
            regionstring = regionstring[3:]
        refseq = fafh.fetch(region=regionstring)
    refdh = pd.DataFrame([a for a in zip(range(startpos, endpos+1, 1), refseq)], 
            columns=["pos", "ref"])
    return refdh



def calcPairwisePearsonR(snvfile1, snvfile2, targetfile, targetname):
    '''calculates the pairwise Pearson r between two named sample ids
    
    '''
    chrom, startpos, endpos = getTargetEditRegion(targetfile, targetname)
    refdh = getReferenceSequence(chrom, startpos, endpos)
    s1counts = sge_counts.getSNVCounts(snvfile1, augment=True, pseudocount=0)
    s2counts = sge_counts.getSNVCounts(snvfile2, augment=True, pseudocount=0)
    s1df = s1counts.melt(id_vars=["chrom", "pos", "target", "repl", "day"],
                     value_vars=["A", "C", "G", "T"],
                     var_name="mutant_allele", value_name="count_s1")
    s2df = s2counts.melt(id_vars=["chrom", "pos", "target", "repl", "day"],
                     value_vars=["A", "C", "G", "T"],
                     var_name="mutant_allele", value_name="count_s2")
    df = s1df.merge(s2df, on=["chrom", "pos", "mutant_allele"])
    df = df[(df["pos"] >= startpos) & (df["pos"] <= endpos)]
    df = df.merge(refdh, on="pos")
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
    samples = getTargetSampleList(targetname, countsdir)
    for day, repllist in samples.items():
        if len(repllist) > 1:
            corrs = []
            for index, repl1 in enumerate(repllist):
                for repl2 in repllist[index+1:]:
                    mycorr = calcPairwisePearsonR(repl1, repl2, targetfile,
                                                  targetname)
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
