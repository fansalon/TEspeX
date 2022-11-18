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

try:
  import argparse
  import sys
  # ensure that only the modules installed within the TEspeX_deps env are loaded - this basically deletes from sys.path all the paths not containing TEspeX_deps
  new_path = []
  for path in sys.path:
    l = path.find("TEspeX_deps")
    if l != -1:
      new_path.append(path)
  # sys.path is now equal to new_path --> if TEspeX_deps env has not been activated sys.path will be an empty list
  sys.path = new_path
  # now import other paths
  import time
  import os
  from os import listdir
  import gzip
  import subprocess
  import math
  import pysam
  import pandas
  from functools import reduce
  import csv
except ModuleNotFoundError:
  print("ERROR: it seems like none of your sys.path paths contains the TEspeX_deps one...")
  print("Did you forget to activate TEspeX_deps environment through source activate TEspeX_deps?")
  sys.exit(1)

__version__ = 'part of TEspeX v2.0.0'

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
  parser.add_argument('--length', type=int, help='length of the read given as input. This is used to calculate STAR index parameters. If your fq/fq.gz file contains reads with different length specify the shorter length [required]', required=True)
  parser.add_argument('--out', type=str, help='directory where the output files will be written. This directory is created by the pipeline, specificy a non-yet-existing directory', required=True)
  parser.add_argument('--num_threads', type=int, default=2, help='number of threads used by STAR and samtools [2]', required=False)
  parser.add_argument('--mask', type=str, default='F', help='fasta file containing sequences to be masked. If this file is provided, the sequences contained within it are considered as coding/non-coding transcripts and are added to the --cdna and --ncrna fasta files. This might be of help if the users wish to consider some specific regions as belonging to coding/non-coding transcripts even though they are not reported in --cdna and --ncrna fasta files. (e.g., N kb downstream to the transcript TTS for a better handling of readthrough process or non-genic TE-derived sequences known to be passively transcribed from criptic promoters). [F]', required=False)
  parser.add_argument('--version', action='version', version='%(prog)s ' + __version__, help='show the version number and exit')

  # create arguments
  arg = parser.parse_args()
  te = os.path.abspath(arg.TE)
  cDNA = os.path.abspath(arg.cdna)
  ncRNA = os.path.abspath(arg.ncrna)
  rl = arg.length
  dir = os.path.abspath(arg.out)
  num_threads = arg.num_threads
  mask = arg.mask
  if mask != "F":
    mask = os.path.abspath(arg.mask)

  # create a list with the arguments that are files
  argList = []
  argList.append(te)
  argList.append(cDNA)
  argList.append(ncRNA)
  if mask != "F":
    argList.append(mask)
  # check that the input files exist
  for i in range(0, len(argList)):
    if os.path.isfile(argList[i]):
      continue
    else:
      print("ERROR!\n%s: no such file or directory" % (argList[i]))
      sys.exit(1)

  return te, cDNA, ncRNA, rl, dir, num_threads, bin_path, mask


# 2.
# this function writes the message to the log file in the output directory
def writeLog(message):
  print(message)
  with open(dir+"/index.log.out", 'a') as logfile:
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

# 5.
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
  starCmd = bin_path + "STAR-2.6.0c/STAR --runThreadN " +str(num_threads)+ " --runMode genomeGenerate --genomeDir " +os.path.abspath(".")+ " --genomeFastaFiles " +genome+ " --genomeSAindexNbases " +str(genomeSAindexNbase)+ " --genomeChrBinNbits " +str(genomeChrBinNbits)
#  starCmd = starCmd + " --limitGenomeGenerateRAM 500000000000" #delete
  bash(starCmd)

  os.chdir(dir)


# check there are no duplicated seq in the reference file
def checkReference(file_name):
  names = []
  with open(file_name) as inpf:
    for line in inpf:
      if line.startswith(">"):
        names.append(line[:-1])
  if len(names) != len(list(set(names))):
    writeLog("ERROR: duplicated sequences found in fasta reference file %s" % (file_name))
    writeLog("This is likely due to the selection of the wrong files from Ensembl/Gencode")
    writeLog("Be sure to have downloaded the *pc_transcripts.fa and the *.lncRNA_transcripts.fa, and NOT the *transcripts.fa, as *transcripts.fa and *.lncRNA_transcripts.fa contain duplicated sequences\n\n")
    sys.exit(1)
  else:
    writeLog("index file is OK")


# main
def main():
  TE, cdna, ncrna, read_length, dir, num_threads, bin_path, maskfile = help()
  os.chdir(dir)
  writeLog("\nuser command line arguments:\nTE file = %s\ncdna file = %s\nncrna file = %s\nreadLength = %s\noutDir = %s\nnum_threads = %s\nmask file = %s " % (TE, cdna, ncrna, read_length, dir, num_threads, maskfile))
  if maskfile == "F":
    writeLog("creating reference file %s/TE_transc_reference.fa concatenating %s, %s and %s" % (dir,TE,cdna,ncrna))
  else:
    writeLog("creating reference file %s/TE_transc_reference.fa concatenating %s, %s, %s and %s" % (dir,TE,cdna,ncrna,maskfile))
  reference = createReference(TE, "_transp")
  reference = createReference(cdna, "_transc")
  reference = createReference(ncrna, "_transc")
  if maskfile == "F":
    writeLog("\nNo masked fasta sequence provided.\n")
  else:
    reference = createReference(maskfile, "_transc")
  # check the reference file does not contain any duplicataed sequences
  checkReference(reference)
  star_ind(reference, read_length)
  writeLog("index job has done!")


if __name__ == "__main__":
  main()
