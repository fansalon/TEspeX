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
import math

__version__ = 'part of TEspeX v1.0.2'

######################################################## Functions used to check parsed arguments fullfil the TEspeX expectations
# Strand
def checkStrand(strand_param):
  strand_list = [ "no", "yes", "reverse"]
  if strand_param not in strand_list:
    print("ERROR!\nunrecognized --strand parameter. Please specify no, yes or reverse")
    sys.exit(1)
# Paired, remove
def checkPrdRm(paired,remove):
  truefalse = [ 'T', 'F' ]
  if paired not in truefalse or remove not in truefalse:
    print("ERROR!\nunrecognized --paired or --remove parameters. Please specify T or F")
    sys.exit(1)
# Sample file
def checkSampleFile(sfile,paired):
  ct = 0
  with open(sfile) as inpt:
    for line in inpt:
      ct += 1
      line = line.strip("\n")
      line = line.split("\t") # line is a list, if SE only 1 element is in the list, if PE 2
      # check the expected number of fq is reported (2 if PE,1 if SE)
      if paired == "T":
        if len(line) != 2:
          print("ERROR: --paired T is specified therefore 2 columns are expected in the --sample file. However, %s are detected at line %s" % (str(len(line)),str(ct)))
          sys.exit()
      else:
        if len(line) != 1:
          print("ERROR: --paired F is specified therefore 1 column is expected in the --sample file. However, %s are detected at line %s" % (str(len(line)),str(ct)))
          sys.exit()
      for ln in line:
        # check is full path
        if not os.path.isabs(ln):
          print("ERROR: the file provided in --sample does not contain absolute paths to fastq/gz files. Please provide absolute paths to the fastq/gz and re-run TEspeX")
          print("Exiting....")
          sys.exit(1)
        # check file exist
        else:
          if not os.path.isfile(ln):
            print("ERROR: file %s does not exist" % (ln))
            print("Exiting....")
            sys.exit(1)
# Check pysam and pandas version
def checkPy():
  pandas_ver = pandas.__version__
  pysam_ver = pysam.__version__
  if pandas_ver != "0.23.0":
    print("ERROR: 0.23.0 pandas version is required, %s detected" % str(pandas_ver))
    print("Please, install the correct version of pandas (pip3 install --user pandas==0.23.0) and re-run TEspeX")
    print("Exiting....")
    sys.exit(1)
  if pysam_ver != "0.15.0" and pysam_ver != "0.15.1":
    print("ERROR: 0.15.0 or 0.15.1 pysam version is required, %s detected" % str(pysam_ver))
    print("Please, install the correct version of pysam and re-run TEspeX")
    print("Exiting....")
    sys.exit(1)
# Check walltime format
def checkWallTime(wwtime):
  if len(wwtime.split(":"))!= 3:
    print("ERROR: value provided as --walltime is not formatted as expected. TEspeX expects the walltime to be reported in hh:mm:ss format (e.g. 12:00:00) whereas you provided: %s" % str(wwtime))
    print("Please, provide --walltime as required and re-run TEspeX")
    print("Exiting....")
    sys.exit(1)
# check input files/dirs exist
def checkInputs(dir,ppipeline,tte,ccDNA,nncRNA,ssample_file):
  # create a list with the arguments that are files
  argList = []
  argList.append(ppipeline)
  argList.append(tte)
  argList.append(ccDNA)
  argList.append(nncRNA)
  argList.append(ssample_file)
  # check that the input files exist
  for i in range(0, len(argList)):
    if os.path.isfile(argList[i]):
      continue
    else:
      print("ERROR!\n%s: no such file or directory" % (argList[i]))
      sys.exit()
  # check the basename of the outdir exists - if yes, create outdir
  dname = os.path.dirname(dir)
  if os.path.isdir(dname):
    try:
      os.mkdir(dir)
    except FileExistsError:
      print("ERROR: "+dir+" directory already exists")
      sys.exit()
  else:
    print("ERROR: dirname of the --out parameter does not exist: %s"  % (dname))
    print("Please, specify the name of a non-existing directory, in an already existing path and retry")
    sys.exit()
######################################################## End checks

# 1.
# define the help function - parse input args
def help():
  # define 2 global variables because they will be used by more than 2 functions
  global dir
  global num_threads

  # this script requires the launcher.py script
  clean = os.path.dirname(os.path.realpath(__file__)) + "/cleanup.py"
  job_index = os.path.dirname(os.path.realpath(__file__)) + "/index.py"
  def fileTrue(file):
    if os.path.isfile(file):
      True
    else:
      print("\nERROR:")
      print("wrapper.py requires index.py and cleanup.py scripts.\nThe 2 scripts should be in the same directory of wrapper.py that seems to be %s" % (os.path.dirname(os.path.realpath(__file__))))
      print("\n")
      sys.exit()
  fileTrue(clean)
  fileTrue(job_index)

  parser = argparse.ArgumentParser()

  # create argument list
  parser.add_argument('--script', type=str, help='path to the TEspeX*.py pipeline for calculation of TE expression [required]', required=True)
  parser.add_argument('--TE', type=str, help='fa/fa.gz file containing TE consensus sequences in fasta format [required]', required=True)
  parser.add_argument('--cdna', type=str, help='fa/fa.gz file containing cdna Ensembl sequences in fasta format [required]', required=True)
  parser.add_argument('--ncrna', type=str, help='fa/fa.gz file containing ncrna Ensembl sequences in fasta format [required]', required=True)
  parser.add_argument('--sample', type=str, help='txt file containing fq/fq.gz FULL PATHS. If reads are single end, one path should be written in each line. If reads are paired end the two mates should be written in the same line separated by \\t [required]', required=True)
  parser.add_argument('--paired', type=str, help='T (true) or F (false). T means the reads are paired and consequently the sample file is expected to contain 2 columns. F means the reads are not paired, sample file is expected to contain  1 single column [required]', required=True)
  parser.add_argument('--length', type=int, help='length of the read given as input. This is used to calculate STAR index parameters. If your fq/fq.gz file contains reads with different length specify the most frequent read length [required]', required=True)
  parser.add_argument('--out', type=str, help='directory where the output files will be written. This directory is created by the pipeline, specificy a non-yet-existing directory', required=True)
  parser.add_argument('--strand', type=str, help='strandeness of the RNAseq library. no = unstranded/htseqcount \'no\', yes = htseqcount \'yes\', reverse = htseqcount \'reverse\'', required=True)
  parser.add_argument('--job', type=str, help='number of jobs that can be run at the same time', required=True)
  parser.add_argument('--q', type=str,help='name of the the queue of your SLURM system you want TEspeX to be run (SBATCH -p parameter) [required]',required=True)
  parser.add_argument('--walltime', type=str,help='walltime of each job TEspeX will launch on your SLURM system reported in hh:mm:ss format (SBATCH -t parameter) [required]',required=True)
  parser.add_argument('--num_threads', type=int, default=4, help='number of threads used by STAR and samtools. Minimum number of threads: 4 [4]', required=False)
  parser.add_argument('--remove', type=str, default='T', help='T (true) or F (false). If this parameter is set to T all the bam files are removed. If it is F they are not removed [T]', required=False)
  parser.add_argument('--version', action='version', version='%(prog)s ' + __version__, help='show the version number and exit')

  # create arguments
  arg = parser.parse_args()
  pipeline = os.path.abspath(arg.script)
  te = os.path.abspath(arg.TE)
  cDNA = os.path.abspath(arg.cdna)
  ncRNA = os.path.abspath(arg.ncrna)
  sample_file = os.path.abspath(arg.sample)
  prd = arg.paired
  rl = arg.length
  dir = os.path.abspath(arg.out)
  strandeness = arg.strand
  njob = arg.job
  num_threads = arg.num_threads
  rm = arg.remove
  queue_slurm = arg.q
  wtime_slurm = arg.walltime

  #### Check the parsed arguments fullfill the TEspeX requirments
  # check number of threads
  if num_threads < 4:
    print("ERROR: when run in wrapper mode at lest 4 threads should be selected. %s selected." % (str(num_threads)))
    sys.exit(1)
  # check strand argument is no, yes or reverse
  checkStrand(strandeness)
  # check paired and rm  are T/F 
  checkPrdRm(prd,rm)
  # check that the file containing the fq path i) really contains fullpath and ii) the fq exist
  checkSampleFile(sample_file,prd)
  # check Pandas and pysam versions
  checkPy()
  # check wall time formatting is the correct one hh:mm:ss
  checkWallTime(wtime_slurm)
  # check input args and create the wd
  checkInputs(dir,pipeline,te,cDNA,ncRNA,sample_file)

  return job_index, clean, pipeline, te, cDNA, ncRNA, sample_file, prd, rl, dir, strandeness, njob, num_threads, rm, queue_slurm, wtime_slurm

# 2.
# this function writes the message to the log file in the output directory
#def writeLog(message):
#  with open(dir+"/Log.launcher.out", 'a') as logfile:
#    logfile.write("[%s] " % (time.asctime()))
#    logfile.write("%s\n" % (message))

# 3.
# this function takes as input a string containing a shell command and executes it
def bash(command):
  cmd = subprocess.Popen(command, shell=True, stderr=subprocess.PIPE, stdout=subprocess.PIPE)
  IDtmp = (cmd.stdout.read()).decode('utf8')		# return the name of the job ID
  print("executing: %s" % (command))
  err, out = cmd.communicate()
  if int(cmd.returncode) != 0:
    print("Error in", command)
    print(err.decode("UTF-8"))
    print(out.decode("UTF-8"))
    sys.exit()
  return IDtmp

# create a file with fq samples to be analyzed by each job
def createSample(fq_file, jobs):
  fq_list = []
  count = 0
  with open(fq_file) as fq:
    for line in fq:
      fq_list.append(line[:-1])
      count += 1

  # check the number of jobs is not greater than the number of samples
  if int(jobs) > int(count):
    print("ERROR: number of jobs (--job parameter) is greater than the number of fq in %s file" % (fq_file))
    sys.exit()

  # split the list 'l' in 'n' lists where 'n' is the number of jobs
  def splitter(l, n):
    dictio = {}
    dictio_val = []
    missing = []
    ratio = math.trunc(len(l) / float(n))
    for i in range(0, int(n)):
      a = (i * ratio)
      b = (a + ratio)
      dictio[i] = l[a:b]
    # create list with dictio values
    for k in dictio.keys():
      for val in dictio[k]:
        dictio_val.append(val)
    # search every element of the input list in the dictio values in order to understand who is missing
    for j in range(0, len(l)):
      if l[j] not in dictio_val:
        missing.append(l[j])
    # add missing fq to the dictio (the 1st missing to the first key, the 2nd to the 2nd key..)
    for z in range(0, len(missing)):
      dictio[z].append(missing[z])
    # create n output file where n = number of keys
    for key in dictio.keys():
      with open("sample"+str(key)+".txt",'w') as outf:
        for v in dictio[key]:
          outf.write("%s\n" % (v))
  # launch splitter function dividing the fastq list for the number of jobs
  splitter(fq_list, jobs)


# create index
def createIndex(indexpy,te_fa,cdna_fa,ncrna_fa,read_lg,cpu,queue,wallt):
  with open("index.job",'w') as ind:
    ind.write("#!/bin/bash\n#\n#SBATCH -p "+ queue +"\n#SBATCH -N 1\n#SBATCH --sockets-per-node=2\n#SBATCH --threads-per-core=2\n#SBATCH --cores-per-socket="+str(int(cpu/4))+ "\n#SBATCH -t "+ wallt +"\ncd $SLURM_SUBMIT_DIR\n\n")
    ind.write("set -e\n\n")
    ind.write("source activate TEspeX_deps\n\n")
    ind.write("time python3 "+indexpy+" --TE "+te_fa+" --cdna "+cdna_fa+" --ncrna "+ncrna_fa+" --length "+str(read_lg)+" --out "+dir+" --num_threads "+str(cpu)+"\n\n")
  cmd_index = "sbatch index.job"
  index_job_id = bash(cmd_index)
  index_job_id_list = [ index_job_id.split(" ")[3].strip("\n") ]

  return index_job_id_list


# create the jobs
def createJob(jobs, py, te, cDna, ncRna, prd, read_leng, strandn, threads, remove, indicix,queue,wallt):
  job_list = []
  for i in range(0, int(jobs)):
    with open("job"+str(i), 'w') as outj:
      job_list.append("job"+str(i))
      outj.write("#!/bin/bash\n#\n#SBATCH -p "+ queue +"\n#SBATCH -N 1\n#SBATCH --sockets-per-node=2\n#SBATCH --threads-per-core=2\n#SBATCH --cores-per-socket="+str(int(threads/4))+ "\n#SBATCH -t "+ wallt +"\ncd $SLURM_SUBMIT_DIR\n\n")
      outj.write("set -e\n\n")
      outj.write("source activate TEspeX_deps\n\n")
      outj.write("time python3 "+py+" --TE "+te+" --cdna "+cDna+" --ncrna "+ncRna+" --sample sample"+str(i)+".txt --paired "+prd+" --length "+str(read_leng)+" --out "+str(i)+ " --strand "+strandn+ " --num_threads "+str(threads)+" --remove "+remove+ " --index "+indicix+"\n\n")

  return job_list

# launch the job
def launchJob(jl,depend_index):
  jobs = []
  for j in jl:
    command = "sbatch --dependency=afterok:" + depend_index[0] +" " +j
    jobID = bash(command)
    jobs.append(jobID.split(" ")[3].strip("\n"))
    job_number = j.split("job")[1]
    full_path = os.getcwd()+"/"+str(job_number)+"/"
    print("please refer to %sLog.file.out" % (full_path))

  return jobs

# clean up. This function must wait all the job to finish, then it should clean and merge everything
def cleanUP(job_id, cleanupy,queue,wallt):
  depend = ''
  for i in range(0, len(job_id)):
    if i != (len(job_id) - 1):
      depend = depend + job_id[i] + ":"
    else:
      depend = depend + job_id[i]
  with open("cleanup.job", 'w') as out:
      out.write("#!/bin/bash\n#\n#SBATCH -p "+ queue +"\n#SBATCH -N 1\n#SBATCH --sockets-per-node=2\n#SBATCH --threads-per-core=2\n#SBATCH --cores-per-socket=1"+ "\n#SBATCH -t "+ wallt +"\ncd $SLURM_SUBMIT_DIR\n\n")
      out.write("set -e\n\n")
      out.write("source activate TEspeX_deps\n\n")
      out.write("time python3 " + cleanupy + " --wd " + dir + " --job " + str(len(job_id)) + "\n" )
  cmd1 = "sbatch --dependency=afterok:" + depend +" cleanup.job"
  bash(cmd1)
  print("launched the cleanup.job job that will wait all the jobs to finish and then make the clean up of everything. Cross the fingers!")

# main
def main():
  index_script, clean_script, pyscript, TE, cdna, ncrna, sample, paired, read_length, dir, strand, num_job, num_threads, remove, slurm_q, slurm_wtime = help()
  os.chdir(dir)
  print("\nuser command line arguments:\nTE file = %s\ncdna file = %s\nncrna file = %s\nsampleFile file = %s\npaired = %s\nreadLength = %s\noutDir = %s\nstrand = %s\nnum_job = %s\nnum_threads = %s \nremove = %s\nq = %s\nwalltime = %s\n" % (TE, cdna, ncrna, sample, paired, read_length, dir, strand, num_job, num_threads, remove, slurm_q, slurm_wtime))
  createSample(sample, num_job)
  print("You have asked to split the analysis in %s different jobs:" % (num_job))
  # create index
  index_job_id = createIndex(index_script,TE,cdna,ncrna,read_length,num_threads,slurm_q,slurm_wtime)
  #  create mapping jobs
  jlist = createJob(num_job, pyscript, TE, cdna, ncrna, paired, read_length, strand, num_threads, remove, dir+"/index",slurm_q,slurm_wtime)
  jid = launchJob(jlist,index_job_id)
  cleanUP(jid, clean_script,slurm_q,slurm_wtime)


if __name__ == "__main__":
  main()
