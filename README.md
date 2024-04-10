# sge-pipeline

## Organization 
SGE files are generally stored in two locations.
* `/net/bbi/vol1/data/sge-seq/` contains only the raw sequencing data in FASTQ format corresponding to the various genes, exons, timepoints, and replicates.
  * `fastq/` contains one directory per gene.  Within the gene-specific folder, each pair of FASTQ files follows a specific naming convention:
  `${GENE}_${EXON}_${REPL}_${DAY}.r1.fastq.gz` for read 1, and `${GENE}_${EXON}_${REPL}_${DAY}.r2.fastq.gz` for read 2.  
    * `${REPL}` is currently one of: `R1R4`, `R2R5`, or `R3R6` for three experimental replicates; `lib` for the SNV library at day 0, or `NC` for the negative (no-editing) control.  
    * `${DAY}` is generally one of `D00`, `D05`, `D11`, `D13`, or `D17`.  `D00` is used only for the SNV library. 
    * `${GENE}_${EXON}` should correspond to an entry in the target file (see the `etc/` folder for more information).

    Example: `BARD1_X4A_R3R6_D05.r1.fastq.gz` 

    Note that these files _are_ backed up.


  * `nobackup/merged/` contains the same set of gene-specific folders, each of which holds the results of merging the overlapping paired-end read files from the `fastq/` directory.  This step is performed automatically when running the pipeline; these files do not need to be manually generated. 
  The files are named according the format `${GENE}_${EXON}_${REPL}_${DAY}.merged.fq.gz`.  
  These files are then processed automatically to remove any reads containing `N` bases, yielding filtered files named according to the format `${GENE}_${EXON}_${REPL}_${DAY}.merged.noNs.fq.gz`

    Note that these files _are not_ backed up as they can be trivially regenerated from the underlying raw data.

* `/net/bbi/vol1/data/sge-analysis/` holds the pipeline output files and its plumbing.  
  * `bin/` contains executable pipeline scripts.  See the `bin/` directory in this repository for more details.
  * `etc/` contains supporting files for each gene and target, including mini reference genomes for each target, annotations of each possible variant in VEP output format, and target definitions.   
  * `lib/` contains python library code for running the pipeline and should not need to be edited.
  * `nobackup/bam/` contains BAM files resulting from mapping the merged, filtered, target-specific FASTQ files against the target-specific reference sequence.  Each BAM file is in a gene-specific directory.  
  * `nobackup/counts/` contains up to three TSV files per target and sample: 
    * `${GENE}_${EXON}_${REPL}_${DAY}.dels.tsv` contains counts of the identified 3bp deletions
    * `${GENE}_${EXON}_${REPL}_${DAY}.snvs.tsv` contains counts of the identified SNVs
    * `${GENE}_${EXON}_${REPL}_${DAY}.readstats.tsv` contains summary statistics about the number of reads going into the analysis, and the numbers of those reads falling into different categories

    Example: `SFPQ_X3B_R1R4_D05.dels.tsv `


## common workflows
This section provides an overview of the steps required to accomplish various tasks in the SGE pipeline.  


### Adding a new gene
A typical workflow is adding a new gene to the pipeline -- specifically, a gene for which no exons have been previously processed.  This step requires a bit of plumbing -- creating directories and some speciality files -- that's not necessary if you just want to add a new exon to an _existing_ gene.  

1. ssh to the appropriate server and establish a `qlogin` session with minimal resource requirements (1 slot, 1 GB RAM).  

1. Activate the `sge` conda enviromment if this is not done automatically for you: 

    `$ conda activate sge`

1. Create the necessary directories:

```
$ cd /net/bbi/vol1/data/
$ GENE=YFG && mkdir -m 0775 sge-seq/fastq/${GENE} sge-seq/nobackup/merged/${GENE}
$ GENE=YFG && mkdir -m 0775 sge-analysis/etc/${GENE} \
  sge-analysis/nobackup/bam/${GENE}/ \
  sge-analysis/nobackup/counts/${GENE}/
```

4. You will need a new targets file for _YFG_.  You may find it easier to copy an existing one and edit it, rather than starting from scratch.
```
$ cd /net/bbi/vol1/data/sge-analysis
$ GENE=YFG && cp etc/BARD1/targets.tsv etc/${GENE}/targets.tsv
```

5. Create a Manifest for the FASTQ files for this gene.  See "Creating a Manifest" below.  

### Adding a new target to an existing gene
In this example, we will use _YFG_ as the name of the existing gene.  It is important to substitute the name of the actual gene when running the commands below.  

1. ssh to the appropriate server and establish a `qlogin` session with minimal resource requirements (1 slot, 1 GB RAM).  

1. Activate the `sge` conda enviromment if this is not done automatically for you: 

    `$ conda activate sge`

1. Update the targets file for _YFG_:
   * Make a backup of the existing targets file:
   
     ```
     $ cp /net/bbi/vol1/data/sge-analysis/etc/YFG/targets.tsv /net/bbi/vol1/data/sge-analysis/etc/YFG/targets.tsv.20240408
     ```
   * Obtain the coordinates and metadata of the new target.  You will need a number of pieces of information:
     * amplicon start and stop coordinates
     * edited region start and stop coordinates
     * positions of and expected nucleotide changes to required edits
     * any germline SNPs within the amplicon coordinates in the HAP1 genome, and other positions that need to be skipped in the analysis
   * Edit the targets file (`emacs` or `vi` or your favorite editor) and add a new row with the metadata and coordinates you obtained
   * Save the targets file

1. Create a reference genome for the new target by running `extractReferenceFastas`.  This command will automatically read the targets file you just updated and will extract the reference sequence corresponding to your new target, and write it to the `--outdir` (output directory) below.  The reference genome will also be indexed for use with `bwa`.  There should be no reason to modify the `outdir` from what is shown, other than to change the gene name.  
    ```
    $ /net/bbi/vol1/data/sge-analysis/bin/extractReferenceFastas \
    --targets /net/bbi/vol1/data/sge-analyis/etc/YFG/targets.tsv \
    --outdir /net/bbi/vol1/data/sge-analysis/etc/YFG/
    ```

1. You will likely want variant annotations.  The default is to use `vep`, but in theory, other annotations can easily be supported (not covered here).  The process works by creating a skeleton VCF file containing all possible SNV variants, then annotating that VCF with `vep` and outputting the results in standard `vep` format.  The process is then repeated for 3-bp deletions.  
    ```
    $ /net/bbi/vol1/data/sge-analysis/bin/getVariantAnnotations \
    --verbose \
    --targetname name_of_new_target \
    --targets /net/bbi/vol1/data/sge-analyis/etc/YFG/targets.tsv \
    --outdir /net/bbi/vol1/data/sge-analyis/etc/YFG/ \
    ```
    The unannotated VCF files are not strictly needed after this command runs and can be deleted, but they take up virtually zero space, so there isn't any reason to worry about them.

1. Creating a Manifest: If you do not have FASTQ files to process, you're done for now.  If you do, you'll want to run them through the pipeline to extract the variant counts.  You'll create a Manifest that shows the source of the files.  Documentation on how to produce this Manifest is forthcoming.  For now, please see one of the existing Manifest files for the current SGE genes, e.g. `/net/bbi/vol1/data/sge-seq/fastq/RAD51D/MANIFEST.txt`

```
$ /net/bbi/vol1/data/sge-analysis/bin/manifestToPipelineInput -h
usage: Convert MANIFEST to pipeline input files [-h] -o OUTPUT [-f FILTER] [-v] N [N ...]

positional arguments:
  N                     path to one or more MANIFEST files

options:
  -h, --help            show this help message and exit
  -o OUTPUT, --output OUTPUT
                        Path to output file for use with pipeline
  -f FILTER, --filter FILTER
                        Only process sampels with MANIFEST date matching <filter>
  -v, --verbose         Turn on verbose output
  
```

### Adding new fastq files for an existing gene
You may want to add a new replicate or timepoint to an existing target.  The process here is very similar to the process described above in "Adding a new target to an existing gene," except that certain steps can be skipped.  Start from step 6 above.

### Score a target

Once the various timepoints, replicates, and libraries have been processed by the pipeline to produce counts files, the scoring process brings those files together to compute and output scores for each variant.  This scoring process is done separately for the SNVs and the deletions, but a single script computes both sets of scores.  Scores are output in TSV format for ease of downstream visualization (etc), but in future work, scores may be output in VCF format.  The script can also produce figures showing the distributions of scores broken down by variant annotation (stop gain, missense, etc).  

The script to perform the analysis is called `scoreTarget`.  It takes a lot of arguments.  To see a list, run `/net/bbi/vol1/data/sge-analysis/bin/scoreTarget --help` on the command line.  Here's an example command for scoring the X4D target of BARD1, producing TSVs and figures for both deletions and SNVs:

```
/net/bbi/vol1/data/sge-analysis/bin/scoreTarget \
  -v \
  -n BARD1_X4D \
  -t /net/bbi/vol1/data/sge-analysis/etc/BARD1/targets.tsv \
  -s /net/bbi/vol1/data/sge-analysis/work/BARD1_X4D.snvscores.tsv \
  -c /net/bbi/vol1/data/sge-analysis/nobackup/counts/BARD1 \
  -V /net/bbi/vol1/data/sge-analysis/tmp/myvcf/BARD1_X4D.snvs.vep.tsv \
  -S /net/bbi/vol1/data/sge-analysis/work/BARD1_X4D.snvfig.html \
  -U /net/bbi/vol1/data/sge-analysis/tmp/myvcf/BARD1_X4D.dels.vep.tsv \
  -d /net/bbi/vol1/data/sge-analysis/work/BARD1_X4D.delscores.tsv \
  -D /net/bbi/vol1/data/sge-analysis/work/BARD1_X4D.delfig.html
```

Output figures can be in either html or png format, auto-detected based on the file extension you specify.  

Scores are currently based on the median log2 ratios between the per-variant frequency at time _t_ and in the library at day 0.  In most cases, _t_ is Day 13, although in some cases, Day 11 may be used.  The timepoint should be auto-detected.  In the future, other scoring strategies may be supported.

The output TSV files include the per-replicate, per-timepoint scores as well as the overall score described above.
