## Python SGE library code

This directory contains python library code to support the SGE pipeline and scripts.  Must be added to your `PYTHONPATH` before running any custom scripts.  For supplied scripts and code, the installed location (currently `/net/bbi/vol1/data/sge-analysis/lib/`) should be automatically recognized, so no specific action is needed.  

Figures in the pipeline are made with [Altair](https://altair-viz.github.io/index.html), and the theme used to (loosely) customize the figures is defined in the `sge_altair.py` module.  Importing this module will automatically set `altair` defaults and enable the theme/customizations.  

