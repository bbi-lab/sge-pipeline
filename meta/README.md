## meta files for SGE pipeline

* `conda-sge.yaml` contains the definition for a conda environment that will allow you to run the various scripts and programs in the pipeline.  Assuming you have a working `conda` installation already (not covered here), run `conda env create -f conda-sge.yaml` to create the environment named `sge`.  Once the environment has been created, run `conda activate sge` to activate and use the environment.
* `seqspec-sge.yaml` contains a seqspec definition for the SGE assay's DNA modality.  The two "onlist" files are supporting files for this spec.  
