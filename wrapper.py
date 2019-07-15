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

# 1.
# define the help function
def help():
  # define 2 global variables because they will be used by more than 2 functions
  global dir
  global num_threads

  # this script requires the launcher.py script
  clean = os.path.dirname(os.path.realpath(__file__)) + "/cleanup.py"
  if os.path.isfile(clean): 
    True
  else:
    print("\nERROR:")
    print("wrapper.py requires cleanup.py script.\nclean.py should be in the same directory of wrapper.py that seems to be %s" % (os.path.dirname(os.path.realpath(__file__))))
    print("\n")
    sys.exit()

  parser = argparse.ArgumentParser()
  
  # create argument list
  parser.add_argument('--script', type=str, help='path to the python pipeline for calculation of TE expression [required]', required=True)
  parser.add_argument('--TE', type=str, help='fa/fa.gz file containing TE consensus sequences [required]', required=True)
  parser.add_argument('--cdna', type=str, help='fa/fa.gz file containing cdna Ensembl sequences [required]', required=True)
  parser.add_argument('--ncrna', type=str, help='fa/fa.gz file containing ncrna Ensembl sequences [required]', required=True)
  parser.add_argument('--sample', type=str, help='txt file containing fq/fq.gz FULL PATHS. If reads are single end, one path should be written in each line. If reads are paired end the two mates should be written in the same line separated by \\t [required]', required=True)
  parser.add_argument('--paired', type=str, help='T (true) or F (false) [required]', required=True)
  parser.add_argument('--length', type=int, help='length of the read given as input. This is used to calculate STAR index parameters [required]', required=True)
  parser.add_argument('--out', type=str, help='directory where the output files will be written', required=True)
  parser.add_argument('--job', type=str, default=1, help='number of jobs that can be run at the same time', required=True)
  parser.add_argument('--num_threads', type=int, default=2, help='number of threads used by STAR and samtools [2]', required=False)
  parser.add_argument('--remove', type=str, default='T', help='if this parameter is set to T all the bam files are removed. If it is F they are not removed [T]', required=False)

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
  njob = arg.job
  num_threads = arg.num_threads
  rm = arg.remove

#  global dir

  # create the outDir
  while True:
    try:
      os.mkdir(dir)
      break
    except FileExistsError:
      print("ERROR: "+dir+" directory already exists")
      sys.exit()
  

  # create a list with the arguments that are files
  argList = []
  argList.append(pipeline)
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
      sys.exit()

  return clean, pipeline, te, cDNA, ncRNA, sample_file, prd, rl, dir, njob, num_threads, rm

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

# create the jobs
def createJob(jobs, py, te, cDna, ncRna, prd, read_leng, threads, remove):
  job_list = []
  for i in range(0, int(jobs)):
    with open("job"+str(i), 'w') as outj:
      job_list.append("job"+str(i))
      outj.write("#!/bin/bash\n#\n#PBS -q regular\n#PBS -j oe\n#PBS -l nodes=1:ppn="+str(threads)+"\n#PBS -l walltime=12:00:00\ncd $PBS_O_WORKDIR\n\n")
      outj.write("time python3 "+py+" --TE "+te+" --cdna "+cDna+" --ncrna "+ncRna+" --sample sample"+str(i)+".txt --paired "+prd+" --length "+str(read_leng)+" --out "+str(i)+ " --num_threads "+str(threads)+" --remove "+remove+"\n\n")

  return job_list 

# launch the job
def launchJob(jl):
  jobs = []
  for j in jl:
    command = "qsub "+j
    jobID = bash(command)
    jobs.append(jobID.split(".master")[0])
    job_number = j.split("job")[1]
    full_path = os.getcwd()+"/"+str(job_number)+"/"
    print("please refer to %sLog.file.out" % (full_path))

  return jobs

# clean up. This function must wait all the job to finish, then it should clean and merge everything
def cleanUP(job_id, cleanupy):
  depend = ''
  for i in range(0, len(job_id)):
    if i != (len(job_id) - 1):
      depend = depend + job_id[i] + ":"
    else:
      depend = depend + job_id[i]
  with open("cleanup.job", 'w') as out:
    out.write("#!/bin/bash\n#\n#PBS -q regular\n#PBS -j oe\n#PBS -l nodes=1:ppn=1\n#PBS -l walltime=01:00:00\n#PBS -W depend=afterok:"+depend+"\ncd $PBS_O_WORKDIR\n\n")
    out.write("time python3 " + cleanupy + " --wd " + dir + " --job " + str(len(job_id)) + "\n" )
  cmd1 = "qsub cleanup.job"
  bash(cmd1)

  print("launched the cleanup.job job that will wait all the jobs to finish and then make the clean up of everything. Cross the fingers!")

# main
def main():
  clean_script, pyscript, TE, cdna, ncrna, sample, paired, read_length, dir, num_job, num_threads, remove = help()
  os.chdir(dir)
  print("\nuser command line arguments:\nTE file = %s\ncdna file = %s\nncrna file = %s\nsampleFile file = %s\npaired = %s\nreadLength = %s\noutDir = %s\nnum_job = %s\nnum_threads = %s \nremove = %s\n" % (TE, cdna, ncrna, sample, paired, read_length, dir, num_job, num_threads, remove))
  createSample(sample, num_job)
  print("You have asked to split the analysis in %s different jobs:" % (num_job))
  jlist = createJob(num_job, pyscript, TE, cdna, ncrna, paired, read_length, num_threads, remove)
  jid = launchJob(jlist)
  cleanUP(jid, clean_script)


if __name__ == "__main__":
  main()
