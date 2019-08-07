# TEspeX

TEspeX (Transposable Elements SPEific eXpression) is a tool for the TE expression quantification from RNA-seq data. The rationale of this pipeline is to map reads against TE consensus sequences, coding transcripts and non-coding transcripts and to select and count reads mapping with best alignment score only against TE consensus sequences. This should avoid the quantification of reads that may be generated from TE-fragments embedded in coding and non-coding annotated transcripts (39% of human protein coding genes contain TE-fragments). 


The pipeline has been written in python3 so **YOU MUST** use python3 and it has been tested on Ubuntu, CentOS and Mac OS X systems.

# How to install TEspeX

## **Prerequisites**
* Python3 (>=3.4.6)
* pip3
* java


## **Unix**

open a Terminal and type:
```
git clone https://github.com/fansalon/TEspeX
```

This should download locally TEspeX.

Copy the downoloaded folder in the directory in which you wish to install TEspeX, move to that directory and type:
```
cd TEspex/
tespex=$PWD
```
A file called 'picard.jar' should be contained in the 'bin/picard' directory.\
To check whether java is properly installed on your machine and picard is working, type:
```
java -jar $tespex/bin/picard/picard.jar
```
If the picard help is printed everything is fine, if an error rises it may be that java is not installed on your machine. Install Java and retry.

Now all the dependencies (STAR2.6.0c, samtools-1.3.1, pandas 0.23.0 and pysam 0.14.1) should be installed in the bin/ directory within the TEspeX/ directory.\
Please install STAR, samtools, pandas and pysam even if they are  already installed on your machine. TEspeX has been tested on these specific versions  and the use of different versions of these softwares may generate different and unpredictable results.

install STAR2.6.0c:
```
cd $tespex/bin
wget -O STAR-2.6.0c.tar.gz https://github.com/alexdobin/STAR/archive/2.6.0c.tar.gz
tar -zxvf STAR-2.6.0c.tar.gz
cd STAR-2.6.0c/bin/
mkdir tespex/
cp Linux_x86_64_static/STAR tespex/
cd tespex/
./STAR --version
```
this should return:
```STAR_2.6.0c``` 

install samtools-1.3.1
```
cd $tespex/bin
wget -O samtools-1.3.1.tar.bz2 https://github.com/samtools/samtools/releases/download/1.3.1/samtools-1.3.1.tar.bz2
tar xjf samtools-1.3.1.tar.bz2
cd samtools-1.3.1
./configure --prefix=$PWD/
(if during the --configure step you encounter an error like: "configure: error: curses development files not found" please relaunch the command adding the "--without-curses" flag)
make
make install
bin/samtools --version
```
this should return something like:
```
samtools 1.3.1
Using htslib 1.3.1
Copyright (C) 2016 Genome Research Ltd.
```


TEspeX takes also advantage of the python3 libraries: sys, time, os, argparse, gzip, subprocess, math, pysam and pandas.
All these libraries except for pysam and pandas are python standard libraries and should not require installation while pysam and pandas do require an installation.

To install pysam and pandas please open a terminal and type:
```
pip3 install --user pandas==0.23.0
pip3 install --user pysam==0.14.1
```

## **Mac OS**

open a Terminal and type:
```
git clone https://github.com/fansalon/TEspeX
```

This should download locally TEspeX.

Copy the downoloaded folder in the directory in which you wish to install TEspeX, move to that directory and type:
```
cd TEspex/
tespex=$PWD
```
A file called 'picard.jar' should be contained in the 'bin/picard' directory.\
To check whether java is properly installed on your machine and picard is working, type:
```
java -jar $tespex/bin/picard/picard.jar
```
If the picard help is printed everything is fine, if an error rises it may be that java is not installed on your machine. Install Java and retry.

Now all the dependencies (STAR2.6.0c, samtools-1.3.1, pandas 0.23.0 and pysam 0.15.1) should be installed in the bin/ directory within the TEspeX/ directory.\
Please install STAR, samtools, pandas and pysam even if they are  already installed on your machine. TEspeX has been tested on these specific versions  and the use of different versions of these softwares may generate different and unpredictable results.

In Mac OS X it seems that pysam 0.14.1 cannot be installed through pip3 (Python3.7), to overcome this issue we  reccomend to download pysam 0.15.1. The pipeline has been tested on Mac OS X with pysam 0.15.1 and generates the same results obtained with pysam 0.14.1 on Unix systems.

install STAR2.6.0c:
```
cd $tespex/bin
curl -L -o STAR-2.6.0c.tar.gz https://github.com/alexdobin/STAR/archive/2.6.0c.tar.gz
tar -zxvf STAR-2.6.0c.tar.gz
cd STAR-2.6.0c/bin/
mkdir tespex/
cp MacOSX_x86_64/STAR tespex/
cd tespex/
./STAR --version
```
this should return:
```STAR_2.6.0c``` 

install samtools-1.3.1
```
cd $tespex/bin
curl -L -o samtools-1.3.1.tar.bz2 https://github.com/samtools/samtools/releases/download/1.3.1/samtools-1.3.1.tar.bz2
tar xjf samtools-1.3.1.tar.bz2
cd samtools-1.3.1
./configure --prefix=$PWD/
(if during the --configure step you encounter an error like: "configure: error: curses development files not found" please relaunch the command adding the "--without-curses" flag)
make
make install
bin/samtools --version
```
this should return something like:
```
samtools 1.3.1
Using htslib 1.3.1
Copyright (C) 2016 Genome Research Ltd.
```

TEspeX takes also advantage of the python3 libraries: sys, time, os, argparse, gzip, subprocess, math, pysam and pandas.
All these libraries except for pysam and pandas are python standard libraries and should not require installation while pysam and pandas do require an installation.

To install pysam and pandas please open a terminal and type:
```
pip3 install --user pandas==0.23.0
pip3 install --user pysam==0.15.1
```



# Input files
* TE consensus sequences in fasta format
* coding transcripts in fasta format (we suggest the 'cdna.fa.gz' fasta file downloaded from ensembl)
* non coding trasncripts in fasta format (we suggest the 'ncrna.fa.gz' fasta file downloaded from ensembl)
* RN-seq data in fastq/fastq.gz format

# Output files
* outfile.txt: txt file containing the raw counts of reads specifically mapping on TEs. The first column contains the TE names as they are in the TE consensus fasta file, the other columns contain the read counts for each fq input file
* mapping_stats.txt: txt file containing mapping statistics 

# How to run TEspeX
TEspeX can be run calling directly the script from the command line. However, for those users disposing a queue managment system we suggest to use a wrapper script that can launch several jobs in parallel saving consistent amount of time. We provide  the wrapper file we use on our PBS queue managment system (see the section TEspeX in  wrapper mode for more info).

  
# TEspeX in standard mode
All the dependencies should have been installed properly. To test this, type:
```
cd $tespex
python3 TEspeX_v0.1.py --help
```

This command shows the help that should be something very simalr to:
```
usage: TExspec_v0.1.py [-h] --TE TE --cdna CDNA --ncrna NCRNA --sample SAMPLE 
                        --paired PAIRED --length LENGTH --out OUT
                       [--num_threads NUM_THREADS] [--remove REMOVE]

arguments:
  -h, --help                  show this help message and exit
  --TE TE                     fa/fa.gz file containing TE consensus sequences in fasta format [required]
  --cdna CDNA                 fa/fa.gz file containing cdna Ensembl sequences in fasta format [required]
  --ncrna NCRNA               fa/fa.gz file containing ncrna Ensembl sequences in fasta format [required]
  --sample SAMPLE             txt file containing fq/fq.gz FULL PATHS. If reads are single end, one path should be written in each line. If reads are paired end the two mates should be written in the same line separated by \t [required]
  --paired PAIRED             T (true) or F (false). T means the reads are paired and consequently the sample file is expected to contain 2 columns. F means the reads are not paired, sample file is expected to contain  1 single column [required]
  --length LENGTH             length of the read given as input. This is used to calculate STAR index parameters. If your fq/fq.gz file contains reads with different length specify the shorter length [required]
  --out OUT                   directory where the output files will be written. This directory is created by the pipeline, specificy a non-yet-existing directory [required]
  --num_threads NUM_THREADS   number of threads used by STAR and samtools [2]
  --remove REMOVE             T (true) or F (false). If this parameter is set to T all the bam files are removed. If it is F they are not removed [T]
```

All the arguments, except fot ```--num_threads``` and ```--remove```, are required. We suggest to use as argument of ```--TE``` argument a fasta file containing TE consensus sequences (from RepBase?) and as arguments of the ```--cdna``` and ```--ncrna``` arguments the transcriptome files containing cdna and ncrna from ensembl (or genecode if working with human or mouse data).\
```--num_threads``` if not specified is set to 2 while ```--remove``` is set to T by default (meaning all the bam files are removed).

In the folder 'example' you can find a copy of the files used to perform the TE expression analysis on 2 *C. elegans* embryonic fastq files. To test whether the pipeline is working properly, please launch it using the input files in the 'example' folder as explained below.

1. first create the sample file typing:
```
ls $tespex/example/*.fastq.gz > $tespex/example/reads.txt
```
  Please notice that the sample file can contain as many file as you want (one per raw - 1 column if SE, 2 column if PE). They will be analized one-by-one by the pipeline.
  
2. launch the pipeline typing the following command:
```
python3 TEspeX_v0.1.py --TE example/RepBase_single_line.fa.gz \
--cdna example/Caenorhabditis_elegans.WBcel235.cdna.all.fa.gz \
--ncrna example/Caenorhabditis_elegans.WBcel235.ncrna.fa.gz \
--sample example/reads.txt --paired F --length 50 --out test
```
Launching this command the pipeline will first merge together the three fasta file creating a reference transcriptome (TE_transc_reference.fa) and then it will create a STAR index of this file using the ```--length``` parameter for the calculation of genomeSAindexNbase and genomeChrBinNbits. The TE_transc_reference.fa is written in the directory indicated with ```--out```  while the index files are contained in the  ```index``` folder within the  ```--out``` folder.\
Then reads of the first fastq sample (SRR3170296_partial.fastq.gz, in this case) are mapped, filtered and counted. In this example, the ```--paired``` parameter is set to F, and so the pipeline is expecting one fastq/fastq.gz file per row in the  ```reads.txt``` file. If you have paired-end data please write the fastq_1 and fastq_2 on the same raw separating them with \t and set ```--paired``` to T.\
All the output files generated during this step are written in the test/SRR3170296_partial folder.\
Then the second sample (SRR3170297_partial.fastq.gz, in this case) will be analyzed and the output files are written in the test/SRR3170297_partial folder.\
When all the samples contained in the ```--sample``` file are analyzed, the  output files are merged together. The raw read counts for each TE, for each sample, are written in a file called outfile.txt in the ```--out``` directory. The file contains in the first column the names of the TEs (as they are in the fasta file) and a column with the raw read counts for each fastq analyzed.\
Moreover a file called mapping_stats.txt containing i) total number of reads, ii) number of mapped reads, iii) number of reads mapping with best alignment score against TEs contained in the ```--TE``` file (please beaware: for each read there could be more than 1 best alignment), iv) number of TE specific reads (reads mapping with best alignment score only on TEs) and v) number of TE aspecific reads (reads mapping with best alignment score on both TEs and coding/noncoding transcripts) is provided.\
The pipeline prints in the Log.final.out all the commands that are launched in real-time, the user can read it to follow all the opearation the pipeline is doing.

The pipeline launched on the example files should take less than 5 minutes and it should create 5 files (TE_transc_reference.fa, TE_transc_reference.fai, Log.file.out, outfile.txt and mapping_stats.txt) and 3 directories (index/, SRR3170296_partial and SRR3170297_partial) within the ```--out``` directory.\
To check the pipeline run correctly, please test there are no differences between the 2 .txt files contained in your ```--out``` folder and the ones conteined in the example one typing:
```
cd $tespex
diff test/outfile.txt example/outfile.txt
diff test/mapping_stats.txt example/mapping_stats.txt
```
If nothing is printed it means all went fine.




# TEspeX in wrapper mode

The wrapper.py is the wrapper file **we use on our PBS queue managment system** thus it is written to fit the requirements of our system. If you dispose a PBS queue managment machine you should modify according to your machine settings:
* queue name (#PBS -q)
* walltime (#PBS -l)








The TEspeX pipeline has been developed by Federico Ansaloni, PhD student in the Computational Genomics lab (SISSA/ISAS - Trieste - Italy) of pof. Remo Sanges.\
To report bugs or suggestions please feel free to write an email to federico.ansaloni@gmail.com



