# TEspeX

TEspeX (Transposable Elements SPEcific eXpression) is a tool for the TE expression quantification from Illummina short RNA-seq reads. The rationale of the pipeline is to map the reads against a reference transcriptome composed by i) TE consensus sequences, ii) coding transcripts and iii) non-coding transcripts and to quantify the expression of TEs avoiding counting reads deriving from TE fragments embedded in coding/non-coding non-TE annotated transcripts.

Techincally, TEspeX:
* merges the 3 fasta files generating a reference transcriptome
* builds the transcriptome index (STAR)
* maps the reads to the transcriptome assigning the best alignment score (AS) to all the alignments showing the best AS
* selects only the best scoring alignments (samtools)
* discards all the reads mapping with best scoring alignments on coding transcripts/non-coding transcripts
* counts only the TE-specific reads


Possible scenarios:
* A read maps with best alignment score to TE consensus sequences but not to coding/non-coding transcripts --> counted as TE-specific
* A read maps with best alignment score to TE consensus sequences but not to coding/non-coding transcripts. However, the reads multi-map to >10 loci and it is aligned to all of them with an alignment score (AS) ranging between the maxAS and the maxAS-1 --> discarded as not assignable to a specific TE subfamily
* A read maps to *BOTH* TE consensus sequences and coding/non-coding transcripts. However, the alignment(s) on TEs is flagged as best whereas the one(s) on non-TE transcripts is not -->  counted as TE-specific
* A read maps with best alignment score to *BOTH* TE consensus sequences and coding/non-coding transcripts --> discarded as it may be generated from TE fragments embedded in coding/non-coding transcripts
* A read maps with best alignment score to coding/non-coding transcripts --> discarded as it is (probably) transcribed from coding/non-coding transcripts.



The pipeline has been written in python3 and it has been tested on Ubuntu, CentOS and Mac OS X systems.

# How to install TEspeX

## **Prerequisites (installed and in $PATH)**
* conda
* java
* zlib (properly installed and configured as in: http://jianghao.wang/post/2017-11-07-install-packages-on-hpc/#install-zlib)


## **Unix**

open a Terminal and type:
```
git clone https://github.com/fansalon/TEspeX
```

This should download locally TEspeX.

Copy the downloaded folder in the directory you wish TEspeX to be installed, move to that directory and type:
```
cd TEspeX/
tespex=$PWD
```
A file called 'picard.jar' is contained in the 'bin/picard' directory.\
To check whether java is properly installed on your machine and picard properly works, type:
```
java -jar $tespex/bin/picard/picard.jar
```
If the picard help is printed everything is fine. If an error rises java is not (properly) installed on your machine. Install Java and retry.

All the dependencies (STAR2.6.0c, samtools-1.3.1, pandas 0.23.0 and pysam 0.15.1) have now to be installed in the bin/ directory within the TEspeX/ directory.\
Please install STAR, samtools, pandas and pysam even if they are already installed on your machine. TEspeX has been tested on these specific versions and the use of different versions of these softwares may generate different and unpredictable results.

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


In order to ensure that TEspeX is used with the python3 version/python3 library versions it has been tested with, a conda environment is created and within the environment all the required libraries are installed.\
**The conda environment needs to be activated every time TEspeX is used.**

To create the conda environment and install the required libraries type:
```
# create the environment using python 3.6
conda create -n TEspeX_deps python=3.6
## --> you will be asked to let conda download and install new packages: type Y

# activate the environment - to be done every time TEspeX is used
source activate TEspeX_deps

# install the required version of pandas and pysam
pip3 install --user pandas==0.23.0
pip3 install --user pysam==0.15.1

# to check everything properly worked
which python3
## --> /path/to/envs/TEspeX_deps/bin/python3
which pip3
## --> /path/to/envs/TEspeX_deps/bin/pip3
python3 --version
## --> Python 3.6.13 :: Anaconda, Inc.
pip3 --version
## --> pip 21.1.3 from /path/to/envs/TEspeX_deps/lib/python3.6/site-packages/pip (python 3.6)

# deactivate the environment
conda deactivate
```



## **Mac OS**

open a Terminal and type:
```
git clone https://github.com/fansalon/TEspeX
```

This should download locally TEspeX.

Copy the downoloaded folder in the directory you wish TEspeX to be installed, move to that directory and type:
```
cd TEspeX/
tespex=$PWD
```
A file called 'picard.jar' should be contained in the 'bin/picard' directory.\
To check whether java is properly installed on your machine and picard properly works, type:
```
java -jar $tespex/bin/picard/picard.jar
```
If the picard help is printed everything is fine. If an error rises java is not (properly) installed on your machine. Install Java and retry.

All the dependencies (STAR2.6.0c, samtools-1.3.1, pandas 0.23.0 and pysam 0.15.1) have now to be installed in the bin/ directory within the TEspeX/ directory.\
Please install STAR, samtools, pandas and pysam even if they are already installed on your machine. TEspeX has been tested on these specific versions and the use of different versions of these softwares may generate different and unpredictable results.


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

In order to ensure that TEspeX is used with the python3 version/python3 library versions it has been tested with, a conda environment is created and within the environment all the required libraries are installed.\
**The conda environment needs to be activated every time TEspeX is used.**

To create the conda environment and install the required libraries type:
```
# create the environment using python 3.6
conda create -n TEspeX_deps python=3.6
## --> you will be asked to let conda download and install new packages: type Y

# activate the environment - to be done every time TEspeX is used
source activate TEspeX_deps

# install the required version of pandas and pysam
pip3 install --user pandas==0.23.0
pip3 install --user pysam==0.15.1

# to check everything properly worked
which python3
## --> /path/to/envs/TEspeX_deps/bin/python3
which pip3
## --> /path/to/envs/TEspeX_deps/bin/pip3
python3 --version
## --> Python 3.6.13 :: Anaconda, Inc.
pip3 --version
## --> pip 21.1.3 from /path/to/envs/TEspeX_deps/lib/python3.6/site-packages/pip (python 3.6)

# deactivate the environment
conda deactivate
```



# Input files
* TE consensus sequences in fasta (fa/fa.gz) format (--TE argument). Follow these instructions to generate a proper TE consensus sequences input file (https://github.com/fansalon/TEconsensus)
* coding transcripts in fasta format (we suggest the 'cdna.fa.gz' fasta file downloaded from ensembl) (--cdna argument)
* non coding trasncripts in fasta format (we suggest the 'ncrna.fa.gz' fasta file downloaded from ensembl) (--ncrna argument)
* RNA-seq data in fastq (fq/fq.gz) format. TEspeX expects the full path of the fq/fq.gz files to be written in a plain txt file (1 file per row. If paired-end files, TEspeX expects the fq/fq.gz to be listed in two tab-separated columns - PE1 in column 1, PE2 of the same fq/fq.gz in column 2)


# Output files
* outfile.txt: txt file containing the raw counts of reads mapping specifically on TEs. The first column contains the TE names as they are in the TE consensus fasta file, the other columns contain the read counts for each fq input file
* mapping_stats.txt: txt file containing mapping statistics

# How to run TEspeX
TEspeX can be run calling directly the script from the command line. However, for those users disposing a queue managment system we suggest to use a wrapper script that can launch several jobs in parallel saving consistent amount of time. Wrapper file we use on our slurm queue managment system is provided (see the section TEspeX in  wrapper mode for more info).

  
# TEspeX in standard mode
All the dependencies should have been installed properly. To test this, type:
```
cd $tespex
source activate TEspeX_deps
python3 TEspeX.py --help
```

This command shows the help that should be something very similar to:
```
usage: TEspeX.py [-h] --TE TE --cdna CDNA --ncrna NCRNA --sample SAMPLE
                 --paired PAIRED --length LENGTH --out OUT --strand STRAND
                 [--num_threads NUM_THREADS] [--remove REMOVE] [--index INDEX]
                 [--version]

optional arguments:
  -h, --help            show this help message and exit
  --TE TE               fa/fa.gz file containing TE consensus sequences in
                        fasta format [required]
  --cdna CDNA           fa/fa.gz file containing cdna Ensembl sequences in
                        fasta format [required]
  --ncrna NCRNA         fa/fa.gz file containing ncrna Ensembl sequences in
                        fasta format [required]
  --sample SAMPLE       txt file containing fq/fq.gz FULL PATHS. If reads are
                        single end, one path should be written in each line.
                        If reads are paired end the two mates should be
                        written in the same line separated by \t [required]
  --paired PAIRED       T (true) or F (false). T means the reads are paired
                        and consequently the sample file is expected to
                        contain 2 columns. F means the reads are not paired,
                        sample file is expected to contain 1 single column
                        [required]
  --length LENGTH       length of the read given as input. This is used to
                        calculate STAR index parameters. If your fq/fq.gz file
                        contains reads with different length specify the most
                        frequent read length [required]
  --out OUT             directory where the output files will be written. This
                        directory is created by the pipeline, specificy a non-
                        yet-existing directory
  --strand STRAND       strandeness of the RNAseq library. no =
                        unstranded/htseqcount 'no', yes = htseqcount 'yes',
                        reverse = htseqcount 'reverse'
  --num_threads NUM_THREADS
                        number of threads used by STAR and samtools [2]
  --remove REMOVE       T (true) or F (false). If this parameter is set to T
                        all the bam files are removed. If it is F they are not
                        removed [T]
  --index INDEX         If you want TEspeX to build the index for you, leave
                        the default value [recommended]. Otherwise provide
                        FULL path to a directoray containing STAR indexes
                        [not_recommended] [F]
  --version             show the version number and exit
```

All the arguments, except fot ```--num_threads```, ```--remove``` and ```--index``` are required.\
We suggest to use as argument of ```--TE``` argument a fasta file containing TE consensus sequences and as arguments of the ```--cdna``` and ```--ncrna``` arguments the transcriptome files containing cdna and ncrna from ensembl (or genecode if working with human or mouse data).\
The ```--length``` of the read is only needed to build the index of the reference transcriptome - in case of trimmed reads just provide the most frequent read lenght.\
As ```--strand``` TEspeX expects the same nomenclature as the one used by htseq-count (if you are unsure about the strandedness of your data please take a look at: https://chipster.csc.fi/manual/library-type-summary.html) .\
By default ```--num_threads``` is set to 2,  ```--remove``` is set to T by default (meaning all the bam files are removed)  and ```--index``` is set to F (meaning TEspeX will take care about index building).

# Test TEspeX

In the folder 'example' a copy of the files used to perform the TE expression analysis on 2 *C. elegans* embryonic fastq files (Tintori SC, et al. - Dev. Cell - 2016 - https://www.ncbi.nlm.nih.gov/pubmed/27554860) is deposited. Please, **please**, **DO** test whether the pipeline is working properly following the steps below:

1. create the list of fq/fq.gz to be analysed:
```
ls $tespex/example/*.fastq.gz > $tespex/example/reads.txt
```
  Please note that the sample file can contain as many file as you want (one per raw - 1 column if SE, 2 columns if PE). They will be analized one-by-one by the pipeline.
  
2. launch the pipeline typing the following command:
```
python3 TEspeX.py --TE example/ce.Dfam.fa.gz \
--cdna example/Caenorhabditis_elegans.WBcel235.cdna.all.fa.gz \
--ncrna example/Caenorhabditis_elegans.WBcel235.ncrna.fa.gz \
--sample example/reads.txt --paired F --length 50 --out test --strand no
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
If nothing is printed it means all went fine. Otherwise, error has rose.




# TEspeX in wrapper mode

The file while wrapper_slurm.py is the **SLURM** wrapper file. 

To see the help, type:
```
cd $tespex
source activate TEspeX_deps
python3 wrapper_slurm.py -h
```

This should print to screen:
```
usage: wrapper_slurm.py [-h] --script SCRIPT --TE TE --cdna CDNA --ncrna NCRNA
                        --sample SAMPLE --paired PAIRED --length LENGTH --out
                        OUT --strand STRAND --job JOB --q Q --walltime
                        WALLTIME [--num_threads NUM_THREADS] [--remove REMOVE]
                        [--version]

optional arguments:
  -h, --help            show this help message and exit
  --script SCRIPT       path to the TEspeX*.py pipeline for calculation of TE
                        expression [required]
  --TE TE               fa/fa.gz file containing TE consensus sequences in
                        fasta format [required]
  --cdna CDNA           fa/fa.gz file containing cdna Ensembl sequences in
                        fasta format [required]
  --ncrna NCRNA         fa/fa.gz file containing ncrna Ensembl sequences in
                        fasta format [required]
  --sample SAMPLE       txt file containing fq/fq.gz FULL PATHS. If reads are
                        single end, one path should be written in each line.
                        If reads are paired end the two mates should be
                        written in the same line separated by \t [required]
  --paired PAIRED       T (true) or F (false). T means the reads are paired
                        and consequently the sample file is expected to
                        contain 2 columns. F means the reads are not paired,
                        sample file is expected to contain 1 single column
                        [required]
  --length LENGTH       length of the read given as input. This is used to
                        calculate STAR index parameters. If your fq/fq.gz file
                        contains reads with different length specify the most
                        frequent read length [required]
  --out OUT             directory where the output files will be written. This
                        directory is created by the pipeline, specificy a non-
                        yet-existing directory
  --strand STRAND       strandeness of the RNAseq library. no =
                        unstranded/htseqcount 'no', yes = htseqcount 'yes',
                        reverse = htseqcount 'reverse'
  --job JOB             number of jobs that can be run at the same time
  --q Q                 name of the the queue of your SLURM system you want
                        TEspeX to be run (SBATCH -p parameter) [required]
  --walltime WALLTIME   walltime of each job TEspeX will launch on your SLURM
                        system reported in hh:mm:ss format (SBATCH -t
                        parameter) [required]
  --num_threads NUM_THREADS
                        number of threads used by STAR and samtools. Minimum
                        number of threads: 4 [4]
  --remove REMOVE       T (true) or F (false). If this parameter is set to T
                        all the bam files are removed. If it is F they are not
                        removed [T]
  --version             show the version number and exit
```

The parameters are exactly the same of TEspeX.py script except for 3 parameters:
* --script: it requires the path to TEspeX.py script
* --q: the name of the SLURM queue to be used for the TEspeX job
* --walltime: walltime of each job TEspeX will launch on your SLURM system reported in hh:mm:ss format (SBATCH -t parameter)
* --job: it requires the number of jobs you want to run at the same time. This depends on the settings of your system. If you can run 40 jobs at the same time and you set ```--job 40``` and you have 80 fq/fq.gz written in the txt file given as input to ```--sample``` the wrapper_slurm.py script will: 
    * subset the ```--sample``` file in 40 sub-files containing 2 (80/40) fq/fq.gz each (named: sample0, sample1, .., sample40)
    * create 40 folders (named: 0, 1, .., 40)
    * launch 40 different jobs (named: job_0, job_1, .., job_40)
    * when all the jobs have finished the cleanup.py job is automatically launched and all the output files are merged together

When all is done you should have in your ```--out``` folder: slurm output files and 4 output files (TE_transc_reference.fa, TE_transc_reference.fai and the 2 output files mapping_stats_total.txt and outfile_total.txt) and 3 directories (index/, mappings/ and tmp/). In the mappings/ directory there is one directory for each fq/fq.gz analyzed containing all the temporary output files.


# How to test for differentially expressed TEs
TEspeX provides a raw-count output file that can be potentially used as input for any of the several tools developed to detect differential expression of genes among two biological conditions (e.g. DESeq2 and edgeR). However, the large majority (all?) of such tools estimate the library size (i.e., sequencing depth) of each analysed sample by summing together the reads mapped on all the genes. While this assumption is perfectly working when analysing gene expression data, this might be unprecise when analysing TE expression data (i.e. the sum of the reads mapping on TE consensus is not a good aproximation of the total library size).

Thus, we **strongly** recommend not to allow the tool to automatically calculate the library size, providing instead as library size the total number of reads TEspeX has succesfully mapped on the reference transcriptome (i.e., TE consensus+cdna+ncrna). This information is contained in the 3rd column of the mapping_stats.txt/mapping_stats_total.txt TEspeX outputs in the working directory.

In edgeR this can be easily done with the following commands:
```
norm <- read.table("mapping_stats_total.txt",header=T,sep='\t')
y <- DGEList(counts=counts,group=group,lib.size = norm$mapped)  # assuming conuts is the matrix containing the raw counts and group the factor containing metadata information
```

It is likely that **the perfect tool** for testing for DE TEs does not exist yet as most of the assumptions made on gene expression data are not necessarily true in the TE scenario. We have internally tested the DE analysis downstream to TEspeX with i) edgeR, ii) DESeq2 and iii) applying a Welch t.test to normalised counts (RPM: raw counts / tot mapped reads \*1M).

Our tests suggest that edgeR works slightly better than both DESeq2 (as it better handles the TEs showing no expression in multiple samples) and t.test (as it is less sensible to the sample size) and we thus suggest to use edgeR. DE testing downstream to TEspeX has been, nevertheless, tested with all the 3 methods providing consistent results in all the tested scenarios.

If you are interested in the development of **the perfect tool** for testing for DE TEs - or at least to try - feel free to contact Federico at (federico.ansaloni@gmail.com).



# Development and help
The TEspeX pipeline has been developed by Federico Ansaloni, former PhD student in the Computational Genomics lab (SISSA/ISAS - Trieste - Italy) of prof. Remo Sanges (https://www.sangeslab.eu). \
Nicolo' Gualandi and Mauro Esposito contributed to the pipeline testing.\
To report bugs or suggestions please feel free to write to the TEspeX supporting group https://groups.google.com/g/tespex-help


