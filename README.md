# TEspeX

Read about TEspeX in:

Ansaloni *et al.*, Bioinformatics, 2022: https://doi.org/10.1093/bioinformatics/btac526<br />
Ansaloni *et al.*, BMC Bioinformatics, 2019: https://doi.org/10.1186/s12859-019-3088-7

# Overview

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

# How to install TEspeX on Unix systems

## **Install TEspeX and TEspeX prerequisites**

**0. Prerequisites**

Please note how this manual assumes you have already installed and in $PATH:
  * gcc
  * make
  * bash

You can easily check this by typing ```which <name>```. If they are not installed, please install them before installing TEspeX.

**1. conda**
  * first check if conda, activate and deactivate are already installed and in $PATH typing:
    * ```which conda```
    * ```which activate```
    * ```which deactivate```
   
If paths are printed, you can directly go to 2., otherwise follow the steps below:

  * download and install miniconda
    * download the installer ```wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh```
    * install it ```bash Miniconda3-latest-Linux-x86_64.sh```
    * follow the prompts on the installer screens
    * there is no need to run conda init. What is, however, required is the add /path/to/miniconda3/bin to the $PATH
    * this depends on where you have installed miniconda3 so first declare where the miniconda3 is  ```conda_dir=</path/to/miniconda3/>```
    * then add it to the .bashrc   ```echo "export PATH=\${PATH}:${conda_dir}/bin" >> ~/.bashrc```
    * refresh it ```source ~/.bashrc```
    * check conda is in $PATH ```which conda```
    


**2. java (JDK)**
  * test if java JDK is already installed and in path:
    * ```java -version```
  
If this prints something like:
```
openjdk version "11.0.11" 2021-04-20
OpenJDK Runtime Environment (build 11.0.11+9-Ubuntu-0ubuntu2.18.04)
OpenJDK 64-Bit Server VM (build 11.0.11+9-Ubuntu-0ubuntu2.18.04, mixed mode, sharing)
```

Skip to 3.; otherwise:
  * ```conda install -c anaconda openjdk```
  
**3. git**
  * test if git is already installed:
    * ```which git```

If the git version is printed you can directly go to 4., otherwise type:

  * ```conda install git```

**4. zlib**

It is likely that you already have zlib installed and properly configured somewhere. To check this type:
  * ```echo $LD_LIBRARY_PATH```
  * ```echo $CFLAGS```
  * ```echo $LDFLAGS```

If paths to zlib lib and include folders are returned, everything should be OK - skip to 5.

Otherwise:
  * ```cd </path/to/miniconda3/>```
  * ```wget https://www.zlib.net/fossils/zlib-1.2.11.tar.gz```
  * ```tar -zxvf zlib-1.2.11.tar.gz```
  * ```cd zlib-1.2.11/```
  * ```./configure --prefix=$PWD/packages```
  * ```make```
  * ```make install```
  * ```export LD_LIBRARY_PATH=$PWD/packages/lib/:$LD_LIBRARY_PATH ```
  * ```export CFLAGS="-I$PWD/packages/include"```
  * ```export LDFLAGS="-L$PWD/packages/lib"```

**5. clone TEspeX and create working directory variable**

  * ```cd </path/where/install/TEspeX>```
  * ```git clone https://github.com/fansalon/TEspeX```
  * ```cd TEspeX/```
  * ```tespex=$PWD```

**6. Picard**

A file called 'picard.jar' is contained in the 'bin/picard' directory.\
To check whether java is properly installed on your machine and picard properly works, type:

  * ```java -jar $tespex/bin/picard/picard.jar```
  
If the picard help is printed everything is fine. If an error rises java is not (properly) installed on your machine. Possible solutions: i) go back to 2. and check conda has successfully installed java, ii) check that the java installed by conda is in $PATH and iii) check that the invoked java is really the one installed by conda.


**7. STAR**

STAR executables should be in the 'bin/STAR-2.6.0c' directory.\
To check whether STAR is properly working, type:

 * ```$tespex/bin/STAR-2.6.0c/STAR --version```

this should return:
```STAR_2.6.0c``` 

**8. samtools**

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

**9. create a conda environment where the required version of python, pysam and pandas are installed**

To create the conda environment and install the required libraries type:
```
conda create -n TEspeX_deps --override-channels -c bioconda -c defaults python=3.6 pandas=0.23.0 pysam'>=0.15.0,<=0.15.1'
## --> you will be asked to let conda download and install new packages: type Y

# activate the environment - to be done every time TEspeX is used
# **PLEASE DO NOT CHANGE THE ENVIRONMENT NAME AS IT WILL COMPROMISE THE FUNCTIONING OF THE PIPELINE**
source activate TEspeX_deps

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

**10. Test TEspeX has been successfully installed**

In the folder 'example' a copy of the files used to perform the TE expression analysis on 2 *C. elegans* embryonic fastq files (Tintori SC, et al. - Dev. Cell - 2016 - https://www.ncbi.nlm.nih.gov/pubmed/27554860) is deposited. Please, **please**, **DO** test whether the pipeline is working properly following the steps below:

  * change to TEspeX  dir and create the list of fq/fq.gz to be analysed:
    * ```cd $tespex```
    * ```ls $tespex/example/*.fastq.gz > $tespex/example/reads.txt```
  * activate TEspeX, **this has to be done every time TEspeX is launched**
    * ```source activate TEspeX_deps```
  * launch the pipeline typing the following command:
    ```
    python3 TEspeX.py --TE example/ce.Dfam.fa.gz \
    --cdna example/Caenorhabditis_elegans.WBcel235.cdna.all.fa.gz \
    --ncrna example/Caenorhabditis_elegans.WBcel235.ncrna.fa.gz \
    --sample example/reads.txt --paired F --length 50 --out test --strand no
    ```
  * check the output files are as expected:
    * ```cd $tespex```
    * ```diff test/outfile.txt example/outfile.txt```
    * ```diff test/mapping_stats.txt example/mapping_stats.txt```
    
    If nothing is printed it means all went fine. Otherwise, error has rose.
    
    
    
# Input files
* TE consensus sequences in fasta (fa/fa.gz) format (--TE argument). Follow these instructions to generate a proper TE consensus sequences input file (https://github.com/fansalon/TEconsensus)
* coding transcripts in fasta format (we suggest the '*pc_transcripts.fa' fasta file downloaded from gencode) (--cdna argument)
* non coding trasncripts in fasta format (we suggest the '*.lncRNA_transcripts.fa' fasta file downloaded from gencode) (--ncrna argument)
* please, pay specific attention to the selection of the input files. In gencode, the *transcripts.fa and *.lncRNA_transcripts.fa contain duplicated sequences. TEspeX requires no duplicated sequence headers are found across the 3 input fasta files. To this end, prior to start the analysys, TEspeX will check for duplicated sequence headers across the 3 input fasta files and stop the run if found
* RNA-seq data in fastq (fq/fq.gz) format. TEspeX expects the full path of the fq/fq.gz files to be written in a plain txt file (1 file per row. If paired-end files, TEspeX expects the fq/fq.gz to be listed in two tab-separated columns - PE1 in column 1, PE2 of the same fq/fq.gz in column 2)
```
# If single-ended
/path/to/fq1
/path/to/fq2
/path/to/fq3
/path/to/fq4

# If paired-ended
/path/to/fq1_1  /path/to/fq1_2
/path/to/fq2_1  /path/to/fq2_2
/path/to/fq3_1  /path/to/fq3_2
/path/to/fq4_1  /path/to/fq4_2
```


# Output files
* outfile.txt: txt file containing the raw counts of reads mapping specifically on TEs. The first column contains the TE names as they are in the TE consensus fasta file, the other columns contain the read counts for each fq input file
* mapping_stats.txt: txt file containing mapping statistics

# How to run TEspeX
TEspeX can be run calling directly the script from the command line. However, for those users disposing a queue managment system we suggest to use a wrapper script that can launch several jobs in parallel saving consistent amount of time. Wrapper file we use on our slurm queue managment system is provided (see the section TEspeX in  wrapper mode for more info).

  
# TEspeX in standard mode

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
                 [--mask MASK] [--version]

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
  --mask MASK           fasta file containing sequences to be masked. If this
                        file is provided, the sequences contained within it
                        are considered as coding/non-coding transcripts and
                        are added to the --cdna and --ncrna fasta files. This
                        might be of help if the users wish to consider some
                        specific regions as belonging to coding/non-coding
                        transcripts even though they are not reported in
                        --cdna and --ncrna fasta files. (e.g., N kb downstream
                        to the transcript TTS for a better handling of
                        readthrough process or non-genic TE-derived sequences
                        known to be passively transcribed from criptic
                        promoters). [F]
  --multimap MULTIMAP   maximum number of loci a read/read pair is allowed to
                        map to. This value will be provided to STAR
                        --outFilterMultimapNmax and --winAnchorMultimapNmax
                        parameters. It is warmly suggested not to change the
                        default value [10]
  --version             show the version number and exit
```

All the arguments, except fot ```--num_threads```, ```--remove```, ```--index```, ```--mask``` and ```--multimap``` are required.\
We suggest to use as argument of ```--TE``` argument a fasta file containing TE consensus sequences and as arguments of the ```--cdna``` and ```--ncrna``` arguments the transcriptome files containing cdna and ncrna from ensembl (or genecode if working with human or mouse data).\
The ```--length``` of the read is only needed to build the index of the reference transcriptome - in case of trimmed reads just provide the most frequent read lenght.\
As ```--strand``` TEspeX expects the same nomenclature as the one used by htseq-count (if you are unsure about the strandedness of your data please take a look at: https://chipster.csc.fi/manual/library-type-summary.html) .\
By default ```--num_threads``` is set to 2,  ```--remove``` is set to T by default (meaning all the bam files are removed)  and ```--index``` is set to F (meaning TEspeX will take care about index building).



# TEspeX in wrapper mode

The file wrapper_slurm.py is the **SLURM** wrapper file. 

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
                        [--mask MASK] [--module MODULE] [--version]

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
                        calculateSTAR index parameters. If your fq/fq.gz file
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
  --mask MASK           fasta file containing sequences to be masked. If this
                        file is provided, the sequences contained within it
                        are considered as coding/non-coding transcripts and
                        are added to the --cdna and --ncrna fasta files. This
                        might be of help if the users wish to consider some
                        specific regions as belonging to coding/non-coding
                        transcripts even though they are not reported in
                        --cdna and --ncrna fasta files. (e.g., N kb downstream
                        to the transcript TTS for a better handling of
                        readthrough process or non-genic TE-derived sequences
                        known to be passively transcribed from criptic
                        promoters). [F]
  --module MODULE       list of modules to be uploaded within each job by
                        'module load'. If more than one modules need to be
                        loaded provide a comma separated list (e.g.,
                        java,samtools) Reported modules should be available
                        through the 'module av' command. [F]
  --multimap MULTIMAP   maximum number of loci a read/read pair is allowed to
                        map to. This value will be provided to STAR
                        --outFilterMultimapNmax and --winAnchorMultimapNmax
                        parameters. It is warmly suggested not to change the
                        default value [10]
  --version             show the version number and exit

```

The parameters are exactly the same of TEspeX.py script except for 3 parameters:
* --module: it might be that while running job some installed modules need to be upoloaded. If this is the case, provide the list of modules to this parameter. TEspeX will then take care to upload the listed modules in each job by ```module load <module name>```
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


# Overview
If you found TEspeX usefull for your research, please cite:

Ansaloni *et al.*, Bioinformatics, 2022: https://doi.org/10.1093/bioinformatics/btac526<br />
Ansaloni *et al.*, BMC Bioinformatics, 2019: https://doi.org/10.1186/s12859-019-3088-7


