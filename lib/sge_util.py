import pandas as pd
import pysam


def getTargetEditRegion(targetfile, targetname):
    targetdf = pd.read_csv(targetfile,
                           header=0, sep="\t",
                           dtype={'skip_pos': 'string'})
    chrom = targetdf.loc[targetdf["exonname"] == targetname, "chrom"].values[0]
    startpos = targetdf.loc[targetdf["exonname"] == targetname, "editstart"].values[0]
    endpos = targetdf.loc[targetdf["exonname"] == targetname, "editstop"].values[0]
    return (chrom, startpos, endpos)


def getTargetAmpliconRegion(targetfile, targetname):
    targetdf = pd.read_csv(targetfile,
                           header=0, sep="\t",
                           dtype={'skip_pos': 'string'})
    chrom = targetdf.loc[targetdf["exonname"] == targetname, "chrom"].values[0]
    startpos = targetdf.loc[targetdf["exonname"] == targetname, "ampstart"].values[0]
    endpos = targetdf.loc[targetdf["exonname"] == targetname, "ampstop"].values[0]

    return (chrom, startpos, endpos)


def getTargetRequiredEdits(targetfile, targetname):
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


# read reference sequence, return dataframe
def getReferenceSequence(chrom, startpos, endpos, 
                         refpath="/net/shendure/vol10/nobackup/genome/fasta/Homo_sapiens.GRCh38.dna.toplevel.fa"):
    regionstring = chrom + ":" + str(startpos) + "-" + str(endpos)
    with pysam.FastaFile(refpath) as fafh:
        if regionstring.startswith("chr"):
            regionstring = regionstring[3:]
        refseq = fafh.fetch(region=regionstring)
    refdh = pd.DataFrame([a for a in zip(range(startpos, endpos+1, 1), refseq)], 
            columns=["pos", "ref"])
    return refdh



