#!/usr/bin/python3
# Python Script to run binomial.R 
#
# **Column variables will need to be set below to tell program which ones to analyze

import sys, re, glob, subprocess

def main(arg1, arg2):
   colint_xobs = 5  #R is 1 indexed!! (acc = 6, core = 5)
   colint_xtotal = 4  #Linux is 0 indexed !! for the portals.tsv file (acc = 5, core = 4)

   colint_yobs = 4  #R is 1 indexed!! (mixed = 4, total = 3)
   colint_ytotal = 3 #Linux is 0 indexed !! for the portals.tsv file (mixed = 3, total = 2)

   portal_dict = file2dict(arg1)
   list_files = listDirfiles(arg2)
   for portal, values in portal_dict.items():
      if len(portal) == 0:
         continue
      matching_file = getmatchingfile(list_files, portal)
      if not matching_file:
         print("error: no matching file for ", portal, file=sys.stderr)
         continue
      x_total = values.split("\t")[colint_xtotal]
      y_total = values.split("\t")[colint_ytotal]
      cmd = " ".join(("Rscript", "binomial.R", str(matching_file), str(portal), str(colint_xobs), str(colint_yobs), str(x_total), str(y_total)))
      for line in run_cmd(cmd).split("\n"):
         if not len(line) == 0:
            print("\t".join(line.split()))
   return

def run_cmd(cmd):
   p = subprocess.Popen(cmd, shell=True, executable='/bin/bash', stdout=subprocess.PIPE)
   (output, err) = p.communicate()
   p_status = p.wait()
   return output.decode('utf-8').strip()

def getmatchingfile(file_list, id):
   r = re.compile(".*" + id + ".*")
   matches = list(filter(r.match, file_list))
   if matches:
      return matches[0]
   else:
      return

def listDirfiles(dir):
   return glob.glob(str(dir + "/*"))

def file2dict(filename):
   new_dict = {}
   with open(filename) as f:
      for line in f:
         elements = line.strip().split('\t')
         portal = elements[0]
         new_dict[portal] = "\t".join(elements)
   return new_dict

if len(sys.argv) > 2:
   main(sys.argv[1], sys.argv[2])
else:
   print("improper usage: python3 tool.py ./portals_orthogroup_counts.tsv ./dir_indv_FREQ_files/")
