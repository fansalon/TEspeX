# TEspeX

TEspeX (Transposable Elements SPEific eXpression) is a tool for the TE expression quantification from RNA-seq data. The rationale of this pipeline is to map reads against TE consensus sequences, coding transcripts and non-coding transcripts and to select and count reads mapping with best alignment score only against TE consensus sequences. This should avoid the quantification of reads that may be generated from TE-fragments embedded in coding and non-coding annotated transcripts (39% of human protein coding genes contain TE-fragments). 

If you are curious abou the use of this pipeline, you can take a look here where it has been tested on a *C. elegans* embyonic dataset (Ansaloni F., et al - BMC Bioinformatics, 2019).

The pipeline has been written in python3 so **YOU MUST** use python3 and it has been tested on Ubuntu, CentOS and Mac OS X systems.

# How to install TEspeX

**Unix**

open a Terminal and type:
```
git clone https://github.com/fansalon/TEspeX
```

This should install TEspeX.

Copy the downoloaded folder in the directory in which you wish to install TEspeX, move to that directory and type:
```
cd TEspex/
tespex=$PWD
cd  bin
ls
ls picard
```
If the installation was succesfully a file called 'picard.jar' should be contained in the 'picard' directory.\
To check java is properly installed on your machine and picard is working, type:
```
cd picard
java -jar picard.jar
```
If the picard help is printed everything is fine, if an error rises it may be that java is not installed on your machine. Install Java and retry.

Now all the dependencies (STAR2.6.0c, samtools-1.3.1, pandas 0.23.0 and pysam 0.14.1) should be installed in the bin/ directory within the TEspeX/ directory.\Please install STAR, samtools, pandas and pysam even if they are  already installed on your machina. TEspeX has been tested on these specific versions  and the use of different version of these softwares may generate different and unpredictable results.

install STAR2.6.0c: \
```
cd $tespex/bin
wget -O STAR-2.6.0c.tar.gz "https://github.com/alexdobin/STAR/archive/2.6.0c.tar.gz"
tar -zxvf STAR-2.6.0c.tar.gz
cd STAR-2.6.0c/bin/Linux_x86_64_static/
./STAR --version
```
'STAR_2.6.0c' should be printed to screen.

install samtools-1.3.1\
cd $tespex/bin






Within the downoloaded folder you should have a folder called 'bin/' that contains all the executables file of the required programs (STAR, samtools, Picard and bedtools). The pipeline is written in order to refer to these executables and it is not advisible to change the location of these folders and files. The use of different version of these softwares may generate different and unpredictable results.

TEspeX takes also advantage of the python3 libraries: sys, time, os, argparse, gzip, subprocess, math, pysam and pandas.
All these libraries except for pysam and pandas are python standard libraries and should not require installation while pysam and pandas do require an installation.

To install pysam and pandas please open a terminal and type:

pip3 install pandas==0.23.0\
pip3 install pysam==0.14.1



# How to run TEspeX
TEspeX can be run calling directly the script or through a wrapper (wrapper.py file, contained in the 'master' folder). You should use the wrapper ONLY IF you dispose of a queue managment system.
  
# TEspeX in standard mode
If the installation was successful you should have locally a copy of the folder containing TEspeX. Typing 'ls' you should see a folder called 'TEspeX' o something similar. If not, your installation failed for some reasons. (Are you sure be in the directory in which you call the command 'git clone https://github.com/fansalon/TEspeX' ?)

If you can see the TEspeX/ directory type:

cd TEspeX/\
python3 ./TEspeX.py --help\

This command shows the help that should be something very simalr to:

usage: TExspec_v0.1.py [-h] --TE TE --cdna CDNA --ncrna NCRNA --sample SAMPLE
                      --paired PAIRED --length LENGTH --out OUT
                       [--num_threads NUM_THREADS] [--remove REMOVE]

arguments:\
  **-h**, --help            show this help message and exit\
  **--TE TE**               fa/fa.gz file containing TE consensus sequences
                        [required]\
  **--cdna CDNA**           fa/fa.gz file containing cdna Ensembl sequences
                        [required]\
  **--ncrna NCRNA**         fa/fa.gz file containing ncrna Ensembl sequences
                        [required]\
  **--sample SAMPLE**       txt file containing fq/fq.gz FULL PATHS. If reads are
                        single end, one path should be written in each line.
                        If reads are paired end the two mates should be
                        written in the same line separated by \t [required]\
  **--paired PAIRED**       T (true) or F (false) [required]\
  **--length LENGTH**       length of the read given as input. This is used to
                        calculate STAR index parameters [required]\
  **--out OUT**             directory where the output files will be written\
  **--num_threads NUM_THREADS**
                        number of threads used by STAR and samtools [2]\
  **--remove REMOVE**       if this parameter is set to T all the bam files are
                        removed. If it is F they are not removed [T]
                        
All the arguments, except fot --num_threads and --remove, are required. We suggest to use as argument of --TE file a fasta file containing TE consensus sequences (from RepBase?) and as arguments of the --cdna and --ncrna the transcriptome files containing cdna and ncrna from ensembl (or genecode if working with human or mouse).\
--num_threads if not specified is set to 2 while --remove is set to T by  default (meaning all the bam files are removed).

In the folder 'example' you can find a copy of the files used to perform the TE expression analysis in a sample of *C. elegans*

To test the pipeline launch the following command:



# TEspeX in wrapper mode






