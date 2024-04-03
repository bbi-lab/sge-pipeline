## Supporting files for the pipeline targets

This `etc` folder contains files that provide the target-specific configuration for the pipeline.  Each gene should have a file called `${GENE}.targets.tsv` that defines the list of targets for that gene.  This file should be created before the pipeline is run for the first time, or no output will be generated and the pipeline will not run properly.  

The target file is a tab-delimited file with a fixed list of required fields, and must be manually created. To create a target file for a new gene, it may be easiest to copy an existing file to use as a template, and to update the values in the file.  Target files have the following columns:

* `exonname`: the name of the target region to be edited, typically `${GENE}_${EXON}` (e.g., `BARD1_X1A`).  Not all target exons comprise multiple sub-regions (e.g., `BARD1_X2`).
* `chrom`: the chromosome on which the target resides.  Chromosome names should begin with `chr`
(e.g., `chr20` instead of `20`).
* `editstart`: the leftmost, 1-based, inclusive position of the edited target region
* `editstop`: the rightmost, 1-based, inclusive position of the edited target region.
* `ampstart`: the leftmost, 1-based, inclusive position of the sequenced amplicon
* `ampstop`: the rightmost, 1-based, inclusive position of the sequenced amplicon
* `required_edits`: a comma-separated list of all required edits that must be present within an amplicon. Each element of the list has the format `<Pos><Allele>` (e.g. `1234321G` representing a required G at position 1234321). No spaces between the elements of the list!
* `cigar`: a CIGAR string representing the valid alignment between the reference sequence and the product of the merged paired-end reads of the amplicon.  Typically this is an integer followed by `M` (e.g., `231M`).  Only used for SNVs; ignored for programmed 3-bp deletions.  
* `skip_pos`: a comma-sepated list of all positions to be skipped in the analysis -- for example, because of a germline SNP, or some kind of recurrent sequencing issue.  May be blank. Example: `1234321,1234567`.  No spaces!  Only used for SNVs; ignored for deletions.  

