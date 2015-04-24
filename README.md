# MAWQ (Microbiome Analysis using Workflow of QIIME) v0.1.0

Background
------

[QIIME 1.8.0 (stable public release)](http://qiime.org/) is a bioinformatics pipeline developed to facilitate the analysis of 16S data generated from NGS platforms.
The process of moving from actual sequencing data to meaningful PCoA plots involves a lot of manual steps and invocation of various QIIME scripts. For this pupose, we have developed a standard analysis procedure and wrapped all of the essential steps in one script that can take information directly from the sequencer and perform all of the downstream analyses. This will greatly expedite the process of standard analyses and save a lot and time energy for biologists.
Since the analysis is based on QIIME, we call our wrapper script MAWQ (canonically pronounced mock) as an abbreviation of Microbiome Analysis using Workflow of QIIME.


Required Packages
------

**Python:**

- [FLASh](http://ccb.jhu.edu/software/FLASH/)
- [BIOM 1.3.1](http://biom-format.org/)
- [QIIME 1.8.0 (stable public release)](https://github.com/qiime/qiime-deploy)

The script has been tested on Ubuntu 12.04.3 LTS. Currently, the pipeline is only exclusive to handling [Illumina MiSeq](http://www.illumina.com/systems/miseq.html) data.

Input files
------

The input file for the wrapper is a single tab-delimited file, with the first column containing the name of the output folder generated from Illumina's MiSeq platform. The second column contains the name of the QIIME-based mapping file of that corresponding output folder.

How to use
------

There are two scripts in the src folder. They are called:

- ```mawq_miseq_localhost.py```
- ```mawq_miseq_amazon.py```

The scripts can be run in the terminal simply by typing ```mawq_miseq_localhost.py``` or ```mawq_miseq_amazon.py```. The difference between localhost and amazon script is that the former can be run on an in-house server having QIIME installation whereas the latter can be invoked from an inhouse machine but the analyses will run on the Amazon EC2 instance of QIIME. Additionally, the Amazon EC2 instance key will be needed if using ```mawq_miseq_amazon.py```. 

The scripts run interactively, by asking a series of questions about the parameters needed to be specified for analysis. Once all the questions are answered, the script will run and give meaningful results from raw sequencing data.

Screenshot
------

![mawq_screenshot](https://github.com/alifar76/MAWQ/blob/master/img/mawq_screenshot.png "")
