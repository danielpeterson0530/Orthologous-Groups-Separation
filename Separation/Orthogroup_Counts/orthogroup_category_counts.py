#!/usr/bin/python3
# Python script that returns numbers of orthogroups per catefory for each species
# 
# ** Must be placed in the OrthoFinder/results/Orthogroups/ directory to access files: Orthogroups_UnassignedGenes.tsv and Orthogroups.tsv
# *** OrthoFinder must have been run with file names appended with [filename].P or [filename].NP to recognize pathogenic/non-pathogenic species

import sys

def main():
   core_ortho_dict = sort_core_P_NP()
   acc_ortho_dict = sort_acc_P_NP()

   core_portal_dict = build_core_portal_dict()
   acc_portal_dict = build_acc_portal_dict()

   ortho_dict = merge(core_ortho_dict, acc_ortho_dict)
#   f = open(acc_outfile, "a")
   print("portal", "total_ct", "mixed_ct", "pcore_ct", "npcore_ct", "acc_ct", sep="\t")
   for portal, acc_ogs in acc_portal_dict.items():
      acc_ct = 0
      for acc_og in acc_ogs.split("\t"):
         group = str(ortho_dict.get(acc_og, ""))
         if group == "Acc":
            acc_ct += 1

      core_ogs = core_portal_dict.get(portal,"")
      pcore_ct = 0
      npcore_ct = 0
      mixed_ct = 0
      for core_og in core_ogs.split("\t"):
         group = str(ortho_dict.get(core_og, ""))
         if group == "Pcore":
            pcore_ct += 1
         elif group == "NPcore":
            npcore_ct += 1
         elif group == "Mixed":
            mixed_ct += 1
      total_ct = int(mixed_ct) + int(pcore_ct) + int(npcore_ct) + int(acc_ct)
      print(portal, total_ct, mixed_ct, pcore_ct, npcore_ct, acc_ct, sep="\t")


#      f.write(str(id + "\n"))
#   f.close()
   return


def build_acc_portal_dict():
   portal_dict = {}
   arg2 = "./Orthogroups_UnassignedGenes.tsv"
   for line in open(arg2).read().strip().split("\n")[1:]:
      og = line.split("\t")[0]
      portals = []
      for item in line.split("\t")[1:]:
         if not len(item) == 0:
            for id in item.strip().split(", "):
               portal = id.split("|")[1]
               if not portal in portals:
                  portals.append(portal)
      for portal_id in portals:
         portal_dict[portal_id] = portal_dict.get(portal_id, "") + str(og) + "\t"
   return portal_dict

def build_core_portal_dict():
   portal_dict = {}
   arg1 = "./Orthogroups.tsv"
   for line in open(arg1).read().strip().split("\n")[1:]:
      og = line.split("\t")[0]
      portals = []
      for item in line.split("\t")[1:]:
         if not len(item) == 0:
            for id in item.strip().split(", "):
               portal = id.split("|")[1]
               if not portal in portals:
                  portals.append(portal)
      for portal_id in portals:
         portal_dict[portal_id] = portal_dict.get(portal_id, "") + str(og) + "\t"
   return portal_dict

def sort_acc_P_NP():
   ortho_dict = {}
   arg2 = "./Orthogroups_UnassignedGenes.tsv"
   ortho_list = file2list(arg2)
   names_list = ortho_list[0]
   ortho_list = ortho_list[1:]
   ortho_dict = {}
##build index lists for which columns are UN,P,NP organisms
   p_spp_index = []
   np_spp_index = []
   for name in names_list:
      if name.split(".")[-1] == "P":
         p_spp_index.append(str(names_list.index(name)))
      elif name.split(".")[-1] == "NP":
         np_spp_index.append(str(names_list.index(name)))

##for each orthogroup row
   for line in ortho_list:
      i = 1
      p_ct = 0
      np_ct = 0
##for each organism in orthogroup
      while i < len(names_list):
         if len(line[i]) > 0:
##if 1+ P organism matches, add one to P
            if str(i) in p_spp_index:
               p_ct = p_ct + 1
##if 1+ NP organism matches, add one to NP
            elif str(i) in np_spp_index:
               np_ct = np_ct + 1
         i += 1

##if orthogroup is P organism it is Pacc
      if p_ct > 0 and not np_ct > 0:
         ortho_dict[str(line[0])] = "Acc"

##if orthogroup is NP organism it is NPacc
      elif np_ct > 0 and not p_ct > 0:
         ortho_dict[str(line[0])] = "Acc"

   return ortho_dict


def sort_core_P_NP():
   arg1 = "./Orthogroups.GeneCount.tsv"
   ortho_list = file2list(arg1)
   names_list = ortho_list[0]
   ortho_list = ortho_list[1:]
   ortho_dict = {}

##build index lists for which columns are P,NP organisms
   p_spp_index = []
   np_spp_index = []
   for name in names_list:
      if name.split(".")[-1] == "P":
         p_spp_index.append(str(names_list.index(name)))
      elif name.split(".")[-1] == "NP":
         np_spp_index.append(str(names_list.index(name)))

##for each orthogroup row
   for line in ortho_list:
      i = 1
      p_ct = 0
      np_ct = 0
##for each organism in orthogroup
      while i < len(names_list):
         if int(line[i]) > 0:
##if 1+ P organism matches, add one to P
            if str(i) in p_spp_index:
               p_ct = p_ct + 1
##if 1+ NP organism matches, add one to NP
            elif str(i) in np_spp_index:
               np_ct = np_ct + 1
         i += 1

##if orthogroup has 2+ P organisms and no NP organisms, it is Pcore
      if p_ct > 1 and not np_ct > 0:
         ortho_dict[str(line[0])] = "Pcore"

##if orthogroup has 2+ organisms and no P organisms, it is NPcore
      elif np_ct > 1 and not p_ct > 0:
         ortho_dict[str(line[0])] = "NPcore"

      else:
         ortho_dict[str(line[0])] = "Mixed"

   return ortho_dict

def merge(x, y):
   z = x.copy()
   z.update(y)
   return z

def file2list(filename):
   new_list = []
   strings = open(filename).read().strip().split("\n")
   for line in strings:
      if len(line) == 0:
         continue
      parts = line.split('\t')
      new_list.append(list(parts))
   return new_list

main()
