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
to be written

### Adding new fastq files for an existing gene
You may want to add a new exon region, or to add a new, later timepoint or additional replicate to existing exon target.  

### Adding a new gene
