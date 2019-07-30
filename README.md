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

This should download locally TEspeX.

Copy the downoloaded folder in the directory in which you wish to install TEspeX, move to that directory and type:
```
cd TEspex/
tespex=$PWD
cd bin
ls
ls picard
```
A file called 'picard.jar' should be contained in the 'picard' directory.\
To check whether java is properly installed on your machine and picard is working, type:
```
cd picard
java -jar picard.jar
```
If the picard help is printed everything is fine, if an error rises it may be that java is not installed on your machine. Install Java and retry.

Now all the dependencies (STAR2.6.0c, samtools-1.3.1, pandas 0.23.0 and pysam 0.14.1) should be installed in the bin/ directory within the TEspeX/ directory.\
Please install STAR, samtools, pandas and pysam even if they are  already installed on your machine. TEspeX has been tested on these specific versions  and the use of different versions of these softwares may generate different and unpredictable results.

install STAR2.6.0c:
```
cd $tespex/bin
wget -O STAR-2.6.0c.tar.gz https://github.com/alexdobin/STAR/archive/2.6.0c.tar.gz
tar -zxvf STAR-2.6.0c.tar.gz
cd STAR-2.6.0c/bin/Linux_x86_64_static/
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
\
**Mac OS**

open a Terminal and type:
```
git clone https://github.com/fansalon/TEspeX
```

This should download locally TEspeX.

Copy the downoloaded folder in the directory in which you wish to install TEspeX, move to that directory and type:
```
cd TEspex/
tespex=$PWD
cd bin
ls
ls picard
```
A file called 'picard.jar' should be contained in the 'picard' directory.\
To check whether java is properly installed on your machine and picard is working, type:
```
cd picard
java -jar picard.jar
```
If the picard help is printed everything is fine, if an error rises it may be that java is not installed on your machine. Install Java and retry.

Now all the dependencies (STAR2.6.0c, samtools-1.3.1, pandas 0.23.0 and pysam 0.14.1) should be installed in the bin/ directory within the TEspeX/ directory.\
Please install STAR, samtools, pandas and pysam even if they are  already installed on your machine. TEspeX has been tested on these specific versions  and the use of different versions of these softwares may generate different and unpredictable results.

install STAR2.6.0c:
```
cd $tespex/bin
curl -L -o STAR-2.6.0c.tar.gz https://github.com/alexdobin/STAR/archive/2.6.0c.tar.gz
tar -zxvf STAR-2.6.0c.tar.gz
cd STAR-2.6.0c/bin/MacOSX_x86_64/
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
\
**Windows**

Use an Unix system.
\

# How to run TEspeX
TEspeX can be run calling directly the script or through a wrapper (wrapper.py file, contained in the 'master' folder). You should use the wrapper ONLY IF you dispose of a PBS queue managment system.
  
# TEspeX in standard mode
All the dependencies should have been insalled properly. To test this, type:
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

All the arguments, except fot ```--num_threads``` and ```--remove```, are required. We suggest to use as argument of ```--TE``` argument a fasta file containing TE consensus sequences (from RepBase?) and as arguments of the ```--cdna``` and ```--ncrna``` arguments the transcriptome files containing cdna and ncrna from ensembl (or genecode if working with human or mouse).\
```--num_threads``` if not specified is set to 2 while ```--remove``` is set to T by default (meaning all the bam files are removed).

In the folder 'example' you can find a copy of the files used to perform the TE expression analysis in a sample of *C. elegans*. To test whether the pipeline is working properly, please launch it using the input files in the 'example' folder as explained below.

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
Then reads of the first fastq sample (SRR3170296_partial.fastq.gz) are mapped, filtered and counted. In this example, the ```--paired``` parameter is set to F, and so the pipeline is expecting one fastq/fastq.gz file per row in the  ```reads.txt``` file. If you have paired-end data please write the fastq_1 and fastq_2 on the same raw separating them with \t and set ```--paired``` to T.\
All the output files generated during this step are written in the test/SRR3170296_partial folder.\
Then the second sample (SRR3170297_partial.fastq.gz) will be analyzed and the output files are written in the test/SRR3170297_partial folder.\
When all the samples contained in the ```--sample``` file are analyzed, the raw read counts for each TE, for each sample, are written in a file called outfile.txt in the ```--out``` directory. The file contains in the first column the names of the TEs (as they are in the fasta file) and a column with the raw read counts for each fastq analyzed.\
Moreover a file called mapping_stats.txt containing i) total number of reads, ii) number of mapped reads, iii) number of reads mapping with best alignment score against TEs contained in the ```--TE``` file (please beaware: for each read there could be more than 1 best alignment), iv) number of TE specific reads (reads mapping with best alignment score only on TEs) and v) number of TE aspecific reads (reads mapping with best alignment score on both TEs and coding/noncoding transcripts) is provided.\
The pipeline prints in the Log.final.out all the commands that are launched in real-time, the user can read it to follow all the opearation the pipeline is doing.

The pipeline launched with 20 threads (```--num_threads 20```) should take 4 minutes and it creates 5 files (TE_transc_reference.fa, TE_transc_reference.fai, Log.file.out, outfile.txt and mapping_stats.txt) and 3 directories (index/, SRR3170296_partial and SRR3170297_partial) within the ```--out``` directory.\
To check the pipeline run correctly, please test there are no differences between the 2 .txt files contained in your ```--out``` folder and the ones conteined in the example one typing:
```
cd $tespex
diff test/outfile.txt example/outfile.txt
diff test/mapping_stats.txt example/mapping_stats.txt
```
If nothing is printed it means all went fine.\




# TEspeX in wrapper mode






