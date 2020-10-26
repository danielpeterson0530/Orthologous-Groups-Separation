#!/usr/bin/python3
# Python Script to concatenate sequences from each single-copy-orthogroup for each portal
#

import sys, glob

def main(arg1):
   fasta_dict = {}
   for file in listDirfiles(arg1):
      fastafile2dict(file, fasta_dict)
   for portal, seq in fasta_dict.items():
      print(">"+str(portal))
      print(seq)
   return

def fastafile2dict(filename, fasta_dict):
   strings = open(filename).read().strip().split('>')
   for line in strings:
      if len(line) == 0:
         continue
      parts = line.split()
      label = parts[0]
      portal = label.split("|")[1]
      fasta_dict[portal] = fasta_dict.get(portal, "") + ''.join(parts[1:])
   return fasta_dict

def listDirfiles(dir):
   return glob.glob(str(dir + "/*"))

if len(sys.argv) > 1:
   main(sys.argv[1])
else:
   print("improper usage: python3 tool.py ./single_copy_orthogroups_dir")
