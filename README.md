# TEspeX

TEspeX (Transposable Elements SPEific eXpression) is a tool for the TE expression quantification from RNA-seq data. The rationale of this pipeline is to map reads against TE consensus sequences, coding transcripts and non-coding transcripts and to select and count reads mapping with best alignment score only against TE consensus sequences. This should avoid the quantification of reads that may be generated from TE-fragments embedded in coding and non-coding annotated transcripts (39% of human protein coding genes contain TE-fragments). 

If you are curious abou the use of this pipeline, you can take a look here where it has been tested on a C. elegans embyonic dataset (Ansaloni F., et al - BMC Bioinformatics, 2019).

The pipeline has been written in python3 so YOU MUST use python3 and it has been tested on Ubuntu, CentOS and Mac OS X systems.

# How to install

open a Terminal and type 'git clone https://github.com/fansalon/TEspeX

Within the downoloaded folder you should have a folder called 'bin/' that contains all the executables file of the required programs (STAR, samtools, Picard and bedtools). The pipeline is written in order to refer to these executables and it is not advisible to change the location of these folders and files. The use of different version of these softwares may generate different and unpredictable results.

TEspeX takes also advantage of the python3 libraries: sys, time, os, argparse, gzip, subprocess, math, pysam and pandas.
All these libraries except for pysam and pandas are python standard libraries and should not require installation while pysam and pandas do require an installation.
To install pysam and pandas please open a terminal and type:

pip3 install pandas==0.23.0
pip3 install pysam==0.14.1

Please install pandas and pysam even if they are already installed in your machine. TEspeX has been tested on these specific versions.






