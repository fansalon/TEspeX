# Copyrigth (C) 2019 Federico Ansaloni

# This file is part of TEspeX.

# TEspeX is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# TEspeX is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with TEspeX.  If not, see <http://www.gnu.org/licenses/>.


import sys
import time
import os
from os import listdir
import argparse
import gzip
import subprocess
import math
import pysam
import pandas
from functools import reduce
import csv

# 1.
# define the help function
def help():
  # define 2 global variables because they will be used by more than 2 functions
  global dir
  global num_threads
  global bin_path

  # dir from which the script has been launched. This will be usefull to call all the other softwares that should be in the bin/ folder
  bin_path = os.path.dirname(os.path.realpath(__file__)) + "/bin/"

  parser = argparse.ArgumentParser()

  # create argument list
  parser.add_argument('--TE', type=str, help='fa/fa.gz file containing TE consensus sequences in fasta format [required]', required=True)
  parser.add_argument('--cdna', type=str, help='fa/fa.gz file containing cdna Ensembl sequences in fasta format [required]', required=True)
  parser.add_argument('--ncrna', type=str, help='fa/fa.gz file containing ncrna Ensembl sequences in fasta format [required]', required=True)
  parser.add_argument('--sample', type=str, help='txt file containing fq/fq.gz FULL PATHS. If reads are single end, one path should be written in each line. If reads are paired end the two mates should be written in the same line separated by \\t [required]', required=True)
  parser.add_argument('--paired', type=str, help='T (true) or F (false). T means the reads are paired and consequently the sample file is expected to contain 2 columns. F means the reads are not paired, sample file is expected to contain  1 single column [required]', required=True)
  parser.add_argument('--length', type=int, help='length of the read given as input. This is used to calculate STAR index parameters. If your fq/fq.gz file contains reads with different length specify the shorter length [required]', required=True)
  parser.add_argument('--out', type=str, help='directory where the output files will be written. This directory is created by the pipeline, specificy a non-yet-existing directory', required=True)
  parser.add_argument('--strand', type=str, help='strandeness of the RNAseq library. no = unstranded/htseqcount \'no\', yes = htseqcount \'yes\', reverse = htseqcount \'reverse\'', required=True)
  parser.add_argument('--num_threads', type=int, default=2, help='number of threads used by STAR and samtools [2]', required=False)
  parser.add_argument('--remove', type=str, default='T', help='T (true) or F (false). If this parameter is set to T all the bam files are removed. If it is F they are not removed [T]', required=False)
  parser.add_argument('--index', type=str, default='F', help='If you want TEspeX to build the index for you, leave the default value [reccomended]. Otherwise provide FULL path to a directoray containing STAR indexes [not_reccomended] [F]', required=False)

  # create arguments
  arg = parser.parse_args()
  te = os.path.abspath(arg.TE)
  cDNA = os.path.abspath(arg.cdna)
  ncRNA = os.path.abspath(arg.ncrna)
  sample_file = os.path.abspath(arg.sample)
  prd = arg.paired
  rl = arg.length
  dir = os.path.abspath(arg.out)
  strandeness = arg.strand
  num_threads = arg.num_threads
  rm = arg.remove
  index = arg.index
  if index != "F":
    index = os.path.abspath(arg.index)

#  global dir

  # create the outDir
  while True:
    try:
      os.mkdir(dir)
      break
    except FileExistsError:
      print("ERROR: "+dir+" directory already exists")
      sys.exit(1)


  # create a list with the arguments that are files
  argList = []
  argList.append(te)
  argList.append(cDNA)
  argList.append(ncRNA)
  argList.append(sample_file)
  # check that the input files exist
  for i in range(0, len(argList)):
    if os.path.isfile(argList[i]):
      continue
    else:
      print("ERROR!\n%s: no such file or directory" % (argList[i]))
      sys.exit(1)

  # check strand argument is 0, 1 or 2
  strand_list = [ "no", "yes", "reverse"]
  if strandeness not in strand_list:
    print("ERROR!\nunrecognized --strand parameter. Please specify no, yes or reverse")
    sys.exit(1)
  # check T F args
  truefalse = [ 'T', 'F' ]
  if prd not in truefalse or rm not in truefalse:
    print("ERROR!\nunrecognized --paired or --remove parameters. Please specify T or F")
    sys.exit(1)
  # check index arg
  if index == 'F':
    tmp = True
  else:
    if os.path.isdir(index):
      tmp = True
      prev_dir = '/'.join(index.split("/")[:-1])
      if os.path.isfile(prev_dir+"/TE_transc_reference.fa"):
        tmp2 = True
      else:
        print("ERROR: %s no such file or directory" % (prev_dir+"/TE_transc_reference.fa"))
        print("It seems you are using --index parameter. This is not reccomended. However, if you really want to use it, provide %s/TE_transc_reference.fa file" % (prev_dir))
        sys.exit(1)
    else:
      print("ERROR: %s no such file or directory" % (index))
      print("Please specify --index F [default] of an existing directory containing STAR indexes")
      sys.exit(1)

  return te, cDNA, ncRNA, sample_file, prd, rl, dir, strandeness, num_threads, rm, bin_path, index


# 2.
# this function writes the message to the log file in the output directory
def writeLog(message):
  print(message)
  with open(dir+"/Log.file.out", 'a') as logfile:
    logfile.write("[%s] " % (time.asctime()))
    logfile.write("%s\n" % (message))

# 3.
# this function takes as input a string containing a shell command and executes it
def bash(*command):
  def riseError(popen_var):
    popen_var.wait()
    if int(popen_var.returncode) != 0:
      out, err = popen_var.communicate()
      err_msg = "ERROR!\nexit code: " + str(popen_var.returncode) + "\n" + err.decode("UTF-8") + "\n" + out.decode("UTF-8")
      print(err_msg)
      writeLog(err_msg)
      sys.exit(1)
  # iterate though the arguments
  count = 0
  for arg in command:
     count += 1
     writeLog("executing "+arg)
     # if it is the 1st command, launch it normally
     if int(count) == 1:
       cmd = subprocess.Popen(arg, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
       riseError(cmd)
     # if not launch it using as stdin the stdout of the previous command
     else:
       cmd = subprocess.Popen(arg, shell=True, stdin=cmd.stdout, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
       riseError(cmd)

# 4.
# this function takes as input 3 fasta files: TE, ensembl-cdna, ensembl-ncrna, adds _transc and _transp to fasta names
# and merge the 3 files together creating the reference
def createReference(fasta, tag):
  # this function takes a line of a fasta file and return the name+tag or the nt sequence
  def lineSplitting(riga):
    if '>' in riga:
      riga_new = riga.split()[0] + tag
    else:
      riga_new = riga.split()[0]
    return riga_new

  with open(dir+"/TE_transc_reference.fa",'a') as output:
    # for each arguments of the function, unzip the file if it is zipped
    path_filename, file_extension = os.path.splitext(fasta)
    if file_extension == ".gz":
      writeLog("detected zipped file: %s" % (fasta))
      with gzip.open(fasta, 'rb') as f:
        for lineZIP in f:
          line = lineZIP.decode('utf8').strip()
          line_new = lineSplitting(line)
          output.write("%s\n" % (line_new))
    else:
      with open(fasta) as f:
        for line in f:
          line_new = lineSplitting(line)
          output.write("%s\n" % (line_new))

  fastaRef = (os.path.abspath(dir+"/TE_transc_reference.fa"))

  return fastaRef

# 7.
# this function creates the index of the reference file
def star_ind(genome, r_length):
  # to avoid STAR segmentation fault calculate genomeSAindexNbases and genomeChrBinNbits parameters
  readL = int(r_length)

  faidx_com = bin_path + "samtools-1.3.1/bin/samtools faidx " + genome
  bash(faidx_com)
  bed = pandas.read_table(genome+".fai", sep='\t', header=None)
  genome_length = sum(bed.iloc[:,1])
  chrom = len(bed.iloc[:,0])

  # now calculate the parameters for STAR
  genomeSAindexNbase = int(min(14, (math.log2(genome_length)) / 2 - 1))
  genomeChrBinNbits = int(min(18, math.log2(max(genome_length/chrom,readL))))

  # create dir where write indexes
  os.mkdir("index")
  os.chdir("index")
  # then we can call the STAR index function using the number of threads that is passed from command line
  starCmd = bin_path + "STAR-2.6.0c/bin/tespex/STAR --runThreadN " +str(num_threads)+ " --runMode genomeGenerate --genomeDir " +os.path.abspath(".")+ " --genomeFastaFiles " +genome+ " --genomeSAindexNbases " +str(genomeSAindexNbase)+ " --genomeChrBinNbits " +str(genomeChrBinNbits)
  bash(starCmd)

  os.chdir(dir)

# 8.
# map the reads to the reference. The argument of this function is a file with the full path to the reads
# if the reads are paired they are written on the same line separated by \t
#def star_aln(fq_list, bedReference, pair, rm):
def star_aln(fq_list, strandn, fastaReference, pair, rm, *index_dir):
  output_names = []				# this is the list that will contain the names of the bedtools coverage output files
  statOut = []					# this is the list that will contain mapping statistics
  statOut.append("SRR\ttot\tmapped\tTE-best\tspecificTE\tnot_specificTE")
  
  # define the general command (no reads and no zcat)
  # if  optional arg index_dir is passed it means that indexes have aady been generated and are in index_dir directory
  # otherwise each directory has its own index
  if index_dir:
    index = index_dir[0]
    command = bin_path + "STAR-2.6.0c/bin/tespex/STAR --outSAMunmapped None --outSAMprimaryFlag AllBestScore --outFilterMismatchNoverLmax 0.04 --outMultimapperOrder Random --outSAMtype BAM Unsorted --outStd BAM_Unsorted --runThreadN " +str(num_threads)+ " --genomeDir " +index
  else:
    command = bin_path + "STAR-2.6.0c/bin/tespex/STAR --outSAMunmapped None --outSAMprimaryFlag AllBestScore --outFilterMismatchNoverLmax 0.04 --outMultimapperOrder Random --outSAMtype BAM Unsorted --outStd BAM_Unsorted --runThreadN " +str(num_threads)+ " --genomeDir " +os.path.abspath("index")
  # for every line of the file launch the analysis
  with open(fq_list) as reads:
    for line in reads:
      writeLog("\n\nI am working with %s" % (line[:-1]))
      os.chdir(dir)
      lin = line.split()

      path_filename, file_extension = os.path.splitext(lin[0])	# separate full_path+file and extension
      filenam = os.path.basename(path_filename)			# extrapolate filename without full_path and extension
      filename = os.path.splitext(filenam)[0]			# (if the file is fq.gz .fq will remain in filenam)
      os.mkdir(filename)					# (if the file is fq.gz .fq will remain in filenam)
      os.chdir(filename)

      # single end
      if len(lin) == 1:
        writeLog("single end reads detected")
        if len(lin) == 1 and pair == "T":
          writeLog("ERROR: %s file contains SE reads but you specify PE reads from command line. Exiting.." % (fq_list))
          print("ERROR: %s file contains SE reads but you specify PE reads from command line. Exiting.." % (fq_list))
          sys.exit(1)
        if file_extension == ".gz":
          gzipped = True
          command_final = command + " --readFilesIn " +lin[0]+ " --readFilesCommand gunzip -c > " +filename+ ".bam"
        else:
          gzipped = False
          command_final = command + " --readFilesIn " +lin[0]+ " > " +filename+ ".bam"
      # paired end
      elif len(lin) == 2:
        writeLog("paired end reads detected")
        if len(lin) == 2 and pair == "F":
          writeLog("ERROR: %s file contains PE reads but you specify SE reads from command line. Exiting.." % (fq_list))
          print("ERROR: %s file contains PE reads but you specify SE reads from command line. Exiting.." % (fq_list))
          sys.exit(1)
        if file_extension == ".gz":
          gzipped = True
          command_final = command + " --readFilesIn " +lin[0]+ " " +lin[1]+ " --readFilesCommand gunzip -c > "+filename+ ".bam"
        else:
          gzipped = False
          command_final = command + " --readFilesIn " +lin[0]+ " " +lin[1]+ " > " +filename+ ".bam"
    # 8.1
    # map reads to reference
      bash(command_final)

    # 8.2
    # extract primary alignments (best score alignments)
      if pair == "F":
        if strandn == "no":
          prim_cmd = bin_path + "samtools-1.3.1/bin/samtools view -@ " +str(num_threads)+ " -b -F 0x100 -o " +filename+ "_mappedPrim.bam " +filename+ ".bam"
        elif strandn == "yes":
          prim_cmd = bin_path + "samtools-1.3.1/bin/samtools view -@ " +str(num_threads)+ " -b -F 0x10 -F 0x100 -o " +filename+ "_mappedPrim.bam " +filename+ ".bam"
        elif strandn == "reverse":
          prim_cmd = bin_path + "samtools-1.3.1/bin/samtools view -@ " +str(num_threads)+ " -b -f 0x10 -F 0x100 -o " +filename+ "_mappedPrim.bam " +filename+ ".bam"
      elif pair == "T":
        if strandn == "no":
          prim_cmd = bin_path + "samtools-1.3.1/bin/samtools view -@ " +str(num_threads)+ " -b -f 0x40 -F 0x100 -o " +filename+ "_mappedPrim.bam " +filename+ ".bam"
        elif strandn == "yes":
          prim_cmd = bin_path + "samtools-1.3.1/bin/samtools view -@ " +str(num_threads)+ " -b -f 0x40 -f 0x20 -F 0x100 -o " +filename+ "_mappedPrim.bam " +filename+ ".bam"
        elif strandn == "reverse":
          prim_cmd = bin_path + "samtools-1.3.1/bin/samtools view -@ " +str(num_threads)+ " -b -f 0x40 -f 0x10 -F 0x100 -o " +filename+ "_mappedPrim.bam " +filename+ ".bam"

      bash(prim_cmd)

    # 8.3
    # create list containing name of the reads mapping with best score alignmets only on TEs. These reads are mapping specifically on TEs
      writeLog("selecting reads mapping specifically on TEs")
      TE = []
      mrna = []
      bamfile = pysam.AlignmentFile(filename+"_mappedPrim.bam", "rb")
      for aln in bamfile.fetch(until_eof=True):
        if "_transc" in aln.reference_name:
          mrna.append(aln.query_name)
        elif "_transp" in aln.reference_name:
          TE.append(aln.query_name)
      bamfile.close()
      final = list( set(TE) - set(mrna) ) 			# these reads map with best score only on TEs and not on transcripts
      not_specific = list( set(TE) - set(final) )		# these reads map with best score on both TEs and transcripts

      # write the 2 lists in 2 output files
      with open("specificTE.txt", 'w') as out1:
        for i in range(0, len(final)):
          out1.write("%s\n" % (final[i]))
      with open("not-specificTE.txt", 'w') as out2:
        for j in range(0, len(not_specific)):
          out2.write("%s\n" % (not_specific[j]))

    # 8.4
    # usem picard to extract alignmets corresponing to reads mapping specifically on TEs
      if os.stat("specificTE.txt").st_size != 0:		# if specific read file not empty
        picard = "java -jar " + bin_path + "picard/picard.jar FilterSamReads I="+filename+"_mappedPrim.bam O="+filename+"_specificTE.bam FILTER=includeReadList RLF=specificTE.txt"
        bash(picard)
      else:							# if file empty, create a bam containig only the header
        header_bam = bin_path + "samtools-1.3.1/bin/samtools view -@ " +str(num_threads)+ " -H " +filename+ "_mappedPrim.bam -b -o " +filename+"_specificTE.bam"
        bash(header_bam)

    # 8.5
    # count the reads mapping specifically on TEs using htseqcount
      writeLog("counting TE expression levels considering TE-specific reads containded in " + filename + "_specificTE.bam")
      def counts(bam,fa):
        name = filename
        bam_chr = []
        bamfile = pysam.AlignmentFile(bam, "rb")
        for aln in bamfile.fetch(until_eof=True):
          bam_chr.append(aln.reference_name)
        bamfile.close()

        fa_chr = []
        with open(fa) as fa_f:
          for line in fa_f:
            if line.startswith(">"):
              if "_transp" in line:
                fa_chr.append((line.split()[0]).split(">")[1])

        with open(name+"_counts",'w') as output:
          output.write("TE\t%s\n" % (name))
          for chr in fa_chr:
            output.write("%s\t%s\n" % (chr, bam_chr.count(chr)))
      counts(filename+"_specificTE.bam", fastaReference)
#      # append the name of the bedtools coverage output in the list
      output_names.append(os.path.abspath(".")+"/"+filename+ "_counts")

    # 8.6
    # create a file with statistics
      writeLog("calculating the mapping statistics...")
      # total reads
      def starFinalHandling(star_final):
        tot_map = 0
        mapped = 0
        with open(star_final) as s_stat:
          for line in s_stat:
            if "Number of input reads" in line:
              numb = line.split("\t")[1]
              tot_map = str(numb[:-1])
            elif "Uniquely mapped reads number" in line:
              numb = line.split("\t")[1]
              mapped = str(int(mapped) + int(numb))
            elif "Number of reads mapped to multiple loci" in line:
              numb = line.split("\t")[1]
              mapped = str(int(mapped) + int(numb))
        return tot_map, mapped
      tot, map = starFinalHandling("Log.final.out")

      # reads mapped on TEs
      mapTE = len(final) + len(not_specific)
      # reads specifically mapped against TE
      specific = len(final)
      # reads not specifically mapped against TE
      not_spec = len(not_specific)
      # write into output list
      riga_stat = str(filename)+"\t"+str(tot)+"\t"+str(map)+"\t"+str(mapTE)+"\t"+str(specific)+"\t"+str(not_spec)
      statOut.append(riga_stat)

      # remove the bam files
      if rm == 'T':
        os.remove(filename+".bam")
        os.remove(filename+ "_mappedPrim.bam")
        os.remove(filename+"_specificTE.bam")

  return output_names, statOut

# 9.
# this function takes as input the list 'out' containing the full path to bedtools coverage output files
# and the list 'stat' containing the mapping statitistics for each fq analyzed and write the 2 lists
# in 2 output files
def createOut(out, stat):
  pd = []
  for count in out:
    pdFile = pandas.read_table(count,sep='\t',header=0)
    pd.append(pdFile)
  count_final = reduce(lambda left,right: pandas.merge(left,right,on='TE'), pd)
  count_final.to_csv(dir+"/outfile.txt",sep='\t',index=False,float_format='%.2f')

  # create the output file with mapping statistics
  with open(dir+"/mapping_stats.txt", 'w') as mapS:
    for i in range(0, len(stat)):
      mapS.write("%s\n" % (stat[i]))

  writeLog("DONE")
  writeLog("output files "+dir+"/outfile.txt and "+dir+"/mapping_stats.txt have been correctly created")

# 10.
# main
def main():
  TE, cdna, ncrna, sample, paired, read_length, dir, strand, num_threads, remove, bin_path, indici = help()
  os.chdir(dir)
  writeLog("\nuser command line arguments:\nTE file = %s\ncdna file = %s\nncrna file = %s\nsampleFile file = %s\npaired = %s\nreadLength = %s\noutDir = %s\nstrand = %s\nnum_threads = %s \nremove = %s\nindex = %s\n" % (TE, cdna, ncrna, sample, paired, read_length, dir, strand, num_threads, remove, indici))

  # create reference transcriptome and STAR index
  if indici == 'F':
    writeLog("creating reference file %s/TE_transc_reference.fa" % (dir))
    createReference(TE, "_transp")
    createReference(cdna, "_transc")
    reference = createReference(ncrna, "_transc")
    star_ind(reference, read_length)
    outfile, statfile = star_aln(sample, strand, reference, paired, remove)
  else:
    writeLog("reading index in %s" % (indici))
    prev_dir = '/'.join(indici.split("/")[:-1])
    reference = prev_dir+"/TE_transc_reference.fa"
    writeLog("reading reference in %s" % (reference))
    outfile, statfile = star_aln(sample, strand, reference, paired, remove, indici)
  createOut(outfile, statfile)


if __name__ == "__main__":
  main()
