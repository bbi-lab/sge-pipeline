# sge-pipeline

This repository holds the analysis code and supporting files necessary for computing variant-level scores for Saturation Genome Editing (SGE) datasets.

* `bin` contains executable scripts 
* `lib` contains python library code required for the analysis
* `etc` contains per-gene configuration files, defining the assay targets and amplicons
* `meta` contains additional supporting files, such as a conda environment definition and a `seqspec` specification file. 

## Installation
1. create a sge conda environment:
```
$ mamba env create -f /net/bbi/vol1/data/sge-analysis/etc/sge.environment.yaml
```

2. stick the python library code from `lib` into your site-python somewhere
3. create a top-level directory, `$SGEDIR`, to hold the code
4. create `$SGEDIR/bin` and put the scripts from `bin` inside
5. create `$SGEDIR/etc` to hold per-gene configuration files


## Generate per-variant counts files
To run the pipeline, you need to be able to submit cluster jobs – i.e., with qsub  This document does not cover logging into the cluster or basic cluster job control. 

* Step 1: `mamba activate sge`

* Step 2: Make a TSV file with one line for each dataset you want to process (i.e., one line per pair of FASTQ files).  The file must have five columns.  There should not be a header row.  The columns, in order, are:
  * Gene (e.g., SFPQ)
  * Exon name  (e.g., X2A)
  * Library name, usually following a specific format: Gene_target_replicate_day (e.g., SFPQ_X2A_R1R4_D05)
  * Full path to the read1 FASTQ file (e.g., /net/bbi/vol1/data/SGE/demux/241026_xxx/example.r1.fq.gz)
  * Full path to the read2 FASTQ file

  Each pair of FASTQ files typically follows a specific naming convention: `${GENE}_${EXON}_${REPL}_${DAY}.r1.fastq.gz` for read 1, and `${GENE}_${EXON}_${REPL}_${DAY}.r2.fastq.gz` for read 2.  
  * `${REPL}` is typically  one of: `R1R4`, `R2R5`, or `R3R6` for three experimental replicates; `lib` for the SNV library at day 0, or `NC` for the negative (no-editing) control.  
  * `${DAY}` is generally one of `D00`, `D05`, `D11`, `D13`, or `D17`.  `D00` is used only for the SNV library. 
  * `${GENE}_${EXON}` should correspond to an entry in the target file (see the `etc/` folder for more information).

  * Example: `BARD1_X4A_R3R6_D05.r1.fastq.gz` 

    The library name can be any arbitrary string, but there are several downstream steps (scoring, normalizing) that expect names in a certain format, in order to be able to automatically determine how many replicates are available and what timepoints were used.  So, for initial QC, you could use a name like SFPQ_X2A_testPCRconditionA_R1R4_D13 but eventually you will want to run the pipeline with the naming scheme shown above for production data.

* Step 3: Confirm that each target in the file above has a definition in the gene-specific target files
  A target consists of the gene name and the exon name, separated by an underscore (e.g., `SFPQ_X2A`)
  Gene-specific targets are in `$SGEDIR/<gene>/targets.tsv`
  If there is not a target entry in this file (for example, if there is a new target), you will need to add one before you can run the pipeline.
	
* Step 4: Build the pipeline script
  Note that this step has not been made to be generic -- it still relies on some local cluster configuration, including queue names.  You may need to edit the `makePipeline` script to modify the skeleton of the pipeline files to get this to run with in a different environment.

```
$ $SGEDIR/bin/makePipeline \
  -s <TSV file from step 2> \
  -o <output file for pipeline script> \
  -f <directory to put the merged fastq files in> \
  -b <directory to put the bam files in> \
  -c <directory to put the variant counts files in>
```

As an Example:
```
$ bin/makePipeline \
  -s EXAMPLE_FILELIST.tsv \
  -o $HOME/SGEpipeline.Nov1.sh \
  -f /path/to/merged_fastq_dir/ \
  -b /path/to/bam_dir/ \
  -c /path/to/counts_dir/ \  
  -n 20  # number of tasks to run at once 
```

* Step 5: Run the pipeline
  First, login to a cluster submission node, then:
	```
	$ qsub <pipeline_script.sh>
	```

Typically, the pipeline takes anywhere from five to 60 minutes per sample, depending mostly on the size of the input FASTQ files and the other load on the cluster.  Keep in mind that you will typically be analyzing many samples at the same time in parallel (remember the -n option when building the pipeline to set the concurrency), so the total runtime for a flowcell of data can be less than an hour.

Warning: the pipeline will overwrite any existing files that have the same names; please make sure you don’t clobber your old analysis when re-running the pipeline with the same output directories unless that’s what you intend!  


## Adding a new target
Each SGE target – an exon, portion of an exon, or other region in which SNVs and deletions are programmed – has to have a target definition in order to be analyzed with the SGE pipeline.

* Step 1: Make a backup copy of the targets file before you modify it
```
$ cp $SGEDIR/etc//${GENE}/targets.tsv /some/where/safe/
```

* Step 2: Gather what you’ll need
For each target, you need to know:
  * Amplicon coordinates: to get these, find the primer sequences used for the sequencing amplicon.  Pop them into the UCSC Genome Browser’s “In-Silico PCR” software to get the coordinates.  Make sure you have selected the right genome assembly: “Dec 2013 (GRCh38/hg38)”
  * Edited region coordinates: to get these, find the SGE oligo used for this target, and pop it into the UCSC Genome Browser’s BLAT search – again, making sure the right genome assembly is selected.  
  * The locations of the “fixed edits” – the edits that should be found on every valid read.
  * Any positions that you need to skip in the analysis.  Typically these will not be known in advance, but rather will emerge after running the pipeline and identifying issues.

* Step 3: Edit the targets.tsv file
Please ensure that the fields are separated by tabs and not spaces, and add a new row to the file:
  * In the first field, add the name of the target
  * In the second field, put the chromosome of the target, starting with ‘chr’
  * In the third and fourth fields, add the start and end coordinates of the edited region, as taken directly from the UCSC Genome Browser blat search.  Do not include commas.
  * In the fifth and sixth fields, put the start and end coordinates of the sequencing amplicon, as taken from the UCSC Genome Browser “in silico PCR” output. Do not include commas.
  * In the seventh column, add the fixed edits.  The format for this field is a comma-separated list with no spaces between items. Each item in the list is the coordinate of one of the fixed edits (no commas!) followed by the mutant allele.  
  * In the eighth column, you need to add the “valid CIGAR string” for the amplicon.  In general, you should just use the length of the amplicon – ampstop - ampstart + 1 – followed by the letter M.  Example: for BARD1 X1A, the amplicon start coordinate is 214809432 and the end coordinate is 214809629, making the length 214809629 - 214809432 + 1 = 198. Adding an M produces the valid cigar string “198M”.
  * In the final column, you can optionally add any positions, as a comma-separated list, that need to be skipped in the analysis. 

* Step 4: Extract the reference sequence
Once you have modified the targets file, you will run a program that generates a mini reference genome in fasta format for all targets for a specific gene.
```
$ $SGEDIR/bin/extractReferenceFastas \
  -t $SGEDIR/etc/${GENE}/targets.tsv \
  -o $SGEDIR/etc/${GENE}/
```

Step 5: Get variant annotations using vep
  The procedure (behind the scenes) is to create a VCF file containing all possible variants in the edited region, and (separately) a VCF file for the 3bp deletions.  These VCF files are used as input into the Variant Effect Predictor (vep) which marks up each variant with its predicted effect (nonsense, missense, etc) and outputs a TSV file of the results.  That TSV file then gets used by downstream programs.
```
$ $SGEDIR/bin/getVariantAnnotations \
  -t $SGEDIR/etc/${GENE}/targets.tsv \
  -n <new_target_name> \
  -o $SGEDIR/etc/${GENE}/
```


## Compute scores per variant
Once you have generated all of the per-variant counts files, you will want to compute scores

```
$SGEDIR/bin/scoreSNVs \
  -v
  -o ${GENE}.snvscores.tsv \
  -z ${GENE}.thresholds.tsv \
  -c /path/to/counts/${GENE} \
  -g ${GENE} \
```

For a complete list of options, see `scoreSNVs --help`




