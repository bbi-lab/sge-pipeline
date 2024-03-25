import os
import sys
import glob
import argparse
from collections import defaultdict

import pysam
import altair as alt
import pandas as pd
import numpy as np


def getSNVCounts(filename, augment=True):
    '''read in a counts file for SNVs
    
    args: filename: full path to the file with counts
          augment: if True, adds metadata to the dataframe as extra columns
    '''
    df = pd.read_csv(filename, header=0, sep="\t")
    if augment:
        sampleid = df["sampleid"][0]
        parts = sampleid.split("_")
        day = parts[3]
        repl = parts[2]
        gene = parts[0]
        df["gene"] = gene
        df["repl"] = repl
        df["day"] = day
    df = df.rename(columns={'n_A': 'A', 
                            'n_C': 'C',
                            'n_G': 'G', 
                            'n_T': 'T'})
    df = df.reset_index(drop=True)
    return df


def getReadStats(filename, augment=True):
    '''read in the stats file for SNVs and dels for a given sample

    args: filename: full path to the file with stats
          augment: if True, adds metadata to the dataframe as extra columns
    
    '''    
    df = pd.read_csv(filename, header=0, sep="\t")
    if augment:
        sampleid = df["sampleid"][0]
        parts = sampleid.split("_")
        day = parts[3]
        repl = parts[2]
        gene = parts[0]
        df["gene"] = gene
        df["repl"] = repl
        df["day"] = day

    # compute percentages from counts
    df["pct_bad_cigar"] = df["bad_cigar"] / df["total_reads"] * 100.0
    df["pct_wild_type"] = df["wild_type"] / df["total_reads"] * 100.0
    df["pct_missing_req_edit"] = df["missing_req_edit"] / df["total_reads"] * 100.0
    df["pct_too_many_snvs"] = df["too_many_snvs"] / df["total_reads"] * 100.0
    df["pct_del_plus_errors"] = df["del_plus_errors"] / df["total_reads"] * 100.0
    df["pct_no_snv_edit"] = df["no_snv_edit"] / df["total_reads"] * 100.0
    df["pct_snv_reads"] = df["snv_reads"] / df["total_reads"] * 100.0
    df["pct_deletion_reads"] = df["deletion_reads"] / df["total_reads"] * 100.0
    df = df.reset_index(drop=True)
    
    return df

