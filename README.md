# MAWQ (Microbiome Analysis using Workflow of QIIME) v0.1.0

Background
------

[QIIME 1.8.0 (stable public release)](http://qiime.org/) is a bioinformatics pipeline developed to facilitate the analysis of 16S data generated from NGS platforms.
The process of moving from actual sequencing data to meaningful PCoA plots involves a lot of manual steps and invocation of various QIIME scripts. For this pupose, we have developed a standard analysis procedure and wrapped all of the essential steps in one script that can take information directly from the sequencer and perform all of the downstream analyses. This will greatly expedite the process of standard analyses and saves a lot and time energy for the biologists.
Since the analysis is based on QIIME, we call our wrapper script MAWQ (canonically pronounced mock).


Required Packages
------

**Python:**

- [FLASh](http://ccb.jhu.edu/software/FLASH/)
- [BIOM 1.3.1](http://biom-format.org/)
- [QIIME 1.8.0 (stable public release)](https://github.com/qiime/qiime-deploy)

The script has been tested on Ubuntu 12.04.3 LTS.

Input files
------

The input file for the wrapper is a single tab-delimited file, with the first column containing 

1) **GenomeA_kaas.txt**

2) **GenomeA_kaas.txt**


How to use
------

There are two scripts in the src folder. They are called:

- ```mawq_miseq_localhost.py```
- ```mawq_miseq_amazon.py```
