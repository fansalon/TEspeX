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
  import subprocess
  import pandas as pd
except ModuleNotFoundError:
  print("ERROR: it seems like none of your sys.path paths contains the TEspeX_deps one...")
  print("Did you forget to activate TEspeX_deps environment through source activate TEspeX_deps?")
  sys.exit(1)

__version__ = 'part of TEspeX v1.0.3'

# 1.
# define the help function
def help():
  # define 2 global variables because they will be used by more than 2 functions
  global dir

  parser = argparse.ArgumentParser()
  
  # create argument list
  parser.add_argument('--wd', type=str, help='wrapper.py working directory (--out parameter of wrapper.py)', required=True)
  parser.add_argument('--job', type=int, help='number of jobs (--job parameter of wrapper.py', required=True)
  parser.add_argument('--version', action='version', version='%(prog)s ' + __version__, help='show the version number and exit')

  # create arguments
  arg = parser.parse_args()
  dir = os.path.abspath(arg.wd)
  njob = arg.job

  # check that the input files exist
  if os.path.exists(dir):
    True
  else:
    print("ERROR!\n%s: no such file or directory" % (dir))
    sys.exit()

  return dir, njob

# 2.
# this function takes as input a string containing a shell command and executes it
def bash(command):
  cmd = subprocess.Popen(command, shell=True, stderr=subprocess.PIPE, stdout=subprocess.PIPE)
  print("executing: %s" % (command))
  err, out = cmd.communicate()
  if int(cmd.returncode) != 0:
    print("Error in", command)
    print(err.decode("UTF-8"))
    print(out.decode("UTF-8"))
    sys.exit()

# 3. start to clean!
def clean(dir, jobs):
  os.chdir(dir)
  bash("mkdir tmp/")
  bash("mkdir mappings/")
  # create 2 empty pandas dataframe. These 2 dfs will be used to merge together the results
  out = pd.DataFrame()
  tot = pd.DataFrame()
  for i in range(0, int(jobs)):
    os.chdir(str(i))
    # remove all the directories with index except one that is moved to wd. Do the same for annotation files
  #  if i == 0:
  #    c1 = "mv index/ " + dir
  #    bash(c1)
  #    c2 = "mv TE_transc_reference* " + dir
  #    bash(c2)
  #  else:
  #    bash("rm -r index/")
  #    bash("rm TE_transc_reference*")
    # mv all the other files to tmp/ adding the number of the job to the  end of the file
    c3 = "mv outfile.txt " + dir + "/tmp/outfile_" + str(i) + ".txt"
    bash(c3)
    c4 = "mv mapping_stats.txt " + dir + "/tmp/mapping_stats_" + str(i) + ".txt"
    bash(c4)
    c5 = "mv Log.file.out " + dir + "/tmp/Log.file_" + str(i) + ".out"
    bash(c5)
    # now in each directory there should be only the dir with mappings
    c6 = "mv * " + dir + "/mappings"
    bash(c6)
    # go back to wd and delete the empty dir
    os.chdir(dir)
    c7 = "rm -r " + str(i)
    bash(c7)
    
    # now merge together all the outfile and mapping_stat files
    # mapping_stat
    os.chdir("tmp/")
    file = "mapping_stats_" + str(i) + ".txt"
    df = pd.read_csv(file,header=0,sep='\t')	# read the mapping file
    frames = [ out, df ]			
    out = pd.concat(frames, sort = False)	# cat with the previous one. If it is the 1st cycle the prev. dataframe will be empty
    # outfile
    file2 = "outfile_" + str(i) + ".txt"	# read the output file
    df2 = pd.read_csv(file2,header=0,sep='\t')
    if i == 0:					# if it is the 1st cycle add the column TE to the empty df 'tot'
      tot["TE"] = df2["TE"]
    tot = pd.merge(tot, df2, on = "TE")		# merge together tot and df2 according to TE column

    # return to wd
    os.chdir(dir)
  
  # move the job job.o and sample file to tmp
  c8 = "mv *job* sample*txt " + dir + "/tmp"
  bash(c8)
  # sort mapping file and write to output
  out_sorted = out.sort_values(by=["SRR"])
  out_sorted.to_csv(dir+"/mapping_stats_total.txt", sep = '\t', header = True, index = False)
  # sort outfile
  tot = tot.reindex(sorted(tot.columns), axis = 1)
  te = tot["TE"]					# TE column should always be the first so
  tot.drop(["TE"], axis=1,inplace=True)			# I delete it and insert in 1st position
  tot.insert(0,"TE",te)
  #tot["TE"] = tot["TE"].str.replace("_transp","")
  tot.to_csv(dir+"/outfile_total.txt", sep = '\t', header = True, index = False)

# main
def main():
  dir, num_job = help()
  os.chdir(dir)
  clean(dir, num_job)

###################################### MAIN
if __name__ == "__main__":
  main()
