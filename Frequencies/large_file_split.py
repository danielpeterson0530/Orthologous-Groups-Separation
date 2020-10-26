#!/usr/bin/python3
# Python script to split large .tsv file into individual [id].tsv files based on column
# 
#  ** Used to split OUTPUT.tsv into individual files for enrichment analysis
import sys,re

def main(arg):
   portals = get_portals(arg, 0)
   for portal in portals:
      returnmatchinglines(arg, portal)
   return

def returnmatchinglines(filename, id):
   l = open(str(id +".kegg3.tsv"), "a")
   r = re.compile(".*" + id + ".*")
   with open(filename) as f:
      for line in f:
         if re.search(r, line):
            match = re.findall(r, line)[0]
            l.write(str(match + "\n"))
   l.close()
   return

def get_portals(filename, col):
   portals = []
   with open(filename) as f:
      for line in f:
         elements = line.split()
         portalid = elements[col]
#         qid = elements[col]
#         portalid = qid.split("|")[1]
         if not portalid in portals:
            portals.append(portalid)
   return portals

if len(sys.argv) > 1:
   main(sys.argv[1])
else:
   print("improper usage: python3 tool.py ./filename")
