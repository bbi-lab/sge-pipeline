import string
import pandas as pd
import pysam
import glob
import re

class Target:

    def __init__(self, targetname, targetfile):
        self.targetname = targetname
        self.targetfile = targetfile

        parts = targetname.split("_")
        self.gene = parts[0]
        if parts[1][-1] in (string.ascii_letters):
            self.exon = parts[1][:-1]
        else:
            self.exon = parts[1]

        targetdf = pd.read_csv(targetfile, header=0, 
                               sep="\t", dtype={'skip_pos': 'string'})

        # extract the features from the target dataframe
        self.chrom = targetdf.loc[targetdf["target"] == targetname, "chrom"].values[0]
        self.editstartpos = targetdf.loc[targetdf["target"] == targetname, "editstart"].values[0]
        self.editendpos = targetdf.loc[targetdf["target"] == targetname, "editstop"].values[0]
        self.editregion = self.chrom + "+"
        self.ampstartpos = targetdf.loc[targetdf["target"] == targetname, "ampstart"].values[0]
        self.ampendpos = targetdf.loc[targetdf["target"] == targetname, "ampstop"].values[0]

        # region strings
        self.editregionstring = self.chrom + ":" + str(self.editstartpos) + "-" + str(self.editendpos)
        self.ampregionstring = self.chrom + ":" + str(self.ampstartpos) + "-" + str(self.ampendpos)

        self.required_edits = self.getRequiredEdits(targetdf)
        self.skip_pos = self.getSkipPos(targetdf)
        self.cigar = targetdf.loc[targetdf["target"] == targetname, "cigar"].values[0]
        
        self.refdf = self.getReferenceSequence()
        self.homopolymer_pos = self.findHomopolymers()

    def getRequiredEdits(self, targetdf):
        required_edits = {}
        editstring = targetdf.loc[targetdf["target"] == self.targetname, "required_edits"].values[0]
        tedits = editstring.split(",")
        for edit in tedits:
            pos = int(edit[:-1])
            base = edit[-1]
            required_edits[pos] = base
        return required_edits
    
    def getSkipPos(self, targetdf):
        skip_pos = {}
        skipstring = targetdf.loc[targetdf["target"] == self.targetname, "skip_pos"].values[0]
        if pd.isna(skipstring) or not skipstring:
            return skip_pos
        skip_parts = skipstring.split(",")
        for p in skip_parts:
            skip_pos[int(p)] = True
        return skip_pos
    

    def getReferenceSequence(self, refpath="/net/shendure/vol10/nobackup/genome/fasta/Homo_sapiens.GRCh38.dna.toplevel.fa"):
        '''extracts the reference sequence within the given 1-based start and 1-based 
        end coordinates, inclusive, and returns the sequence as a pandas dataframe with
        columns "pos" (integer position) and "ref" (reference base)
        
        '''
        with pysam.FastaFile(refpath) as fafh:
            if self.ampregionstring.startswith("chr"):
                rs = self.ampregionstring[3:]
            else:
                rs = self.ampregionstring
            refseq = fafh.fetch(region=rs)
        refdf = pd.DataFrame([a for a in zip(range(self.ampstartpos, self.ampendpos+1, 1), refseq)], 
                columns=["pos", "ref"])
        return refdf


    def findHomopolymers(self, min_length=4):
        '''find all positions within the reference sequence that comprise homopolymer runs with 
        length at least <min_length> bases, and return a list of offsets of such positions
        '''
        if self.refdf.empty:
            return []
        poslist = []
        refseq = ''.join(self.refdf["ref"].to_list())
        for base in ("A", "C", "G", "T"):
            groups = re.finditer(r'%s{%d,}' % (base, min_length), refseq)
            for g in groups:
                for pos in range(g.start(), g.end(), 1):
                    poslist.append(pos)
        return poslist

        
 
    def getSNVSampleList(self, countsdir, include_neg=False):
        '''performs a directory lookup of all SNV counts files matching a 
        specific target name in a 
        given counts directory, and returns a dictionary of filenames:
        (string) day --> [ (string) full_filename1, (string) full_filename2, ... ]

        if include_neg is True, the negative control sample is included; otherwise it is
        excluded from the list of samples

        '''
        samples = {}
        if countsdir.endswith("/"):
            countsdir = countsdir[:-1]
        for fullfn in glob.glob("%s/%s_*.snvs.tsv" % (countsdir, self.targetname)):
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
        self.snvsamples = samples
        days = sorted(samples.keys())

        return samples
    
    
    def getDelSampleList(self, countsdir, include_neg=False):
        '''performs a directory lookup of all del counts files matching a 
        specific target name in a 
        given counts directory, and returns a dictionary of filenames:
        (string) day --> [ (string) full_filename1, (string) full_filename2, ... ]

        if include_neg is True, the negative control sample is included; otherwise it is
        excluded from the list of samples

        '''
        samples = {}
        if countsdir.endswith("/"):
            countsdir = countsdir[:-1]
        for fullfn in glob.glob("%s/%s_*.dels.tsv" % (countsdir, self.targetname)):
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
        self.delsamples = samples
        return samples
