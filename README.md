# MAWQ
Microbiome Analysis using Workflow of QIIME

Background
------

[QIIME 1.8.0 (stable public release)](https://github.com/qiime/qiime-deploy) is a bioinformatics pipeline developed to facilitate the analysis of 16S data generated from NGS platforms.
The process of moving from actual sequencing data to meaningful PCoA plots involves a lot of manual steps and invocation of various QIIME scripts. For this pupose, we have decided to wrap all the essential steps in one platform 


Required Packages
------

**Python:**

- [FLASh](http://ccb.jhu.edu/software/FLASH/)
- [BIOM 1.3.1](http://biom-format.org/)
- [QIIME 1.8.0 (stable public release)](https://github.com/qiime/qiime-deploy)

The script has been tested on Ubuntu 12.04.3 LTS.

Input files
------


1) **GenomeA_kaas.txt**

2) **GenomeA_kaas.txt**


How to use
------

There are two scripts in the src folder. They are called:

- ```mawq_miseq_localhost.py```
- ```mawq_miseq_amazon.py```
