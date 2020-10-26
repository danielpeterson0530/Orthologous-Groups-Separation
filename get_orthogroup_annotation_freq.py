#!/usr/bin/python3
# Python Script to determine frequencies for concatentated GO annotation file
# Requires:    2_portal_orthogroup_group_id.tsv
#              functional_annotation_[go/ipr/kegg/etc].tsv
# **Change go_file_col to change counted column and out_extention as needed. 

import sys

def main(arg1, arg2):
   go_file_col = 6
   out_extention = ".GO.FREQ"

   portal_og_dict, og_id_dict, og_group_dict = get_portal_og_dict(arg1)
   go_dict = gofile2dict(arg2, go_file_col)

   for portal, ogs in portal_og_dict.items():
      num_ogs = len(unique(ogs.split("\t")[1:]))
      all_goterms = {}
      accog_goterm = {}
      mixedog_goterm = {}
      coreog_goterm = {}
      og_acc_ct = 0
      og_mixed_ct = 0
      og_core_ct = 0
      for og in unique(ogs.split("\t")):
         if len(og) == 0:
            continue
         group = og_group_dict.get(og,"")
         if not group:
            print("group error", file= sys.stderr)
         if group == "acc":
            og_acc_ct += 1
         elif group == "mixed":
            og_mixed_ct += 1
         elif group == "pcore":
            og_core_ct += 1
         elif group == "npcore":
            og_core_ct += 1

         ids = og_id_dict.get(og,"")
         if not ids:
            print("id error", file= sys.stderr)
         id_term_store = []
         for id in ids.split("\t"):
            if len(id) == 0:
               continue
            if not portal == str(id.split("|")[1]):
               continue
            shortid = "|".join(id.split("|")[0:3])
            go_terms = go_dict.get(shortid, "")
            if go_terms:
               for term in go_terms.split("\t"):
                  if len(term) == 0:
                     continue
                  if term in id_term_store:
                     continue
                  else:
                     id_term_store.append(term)
                  if not term in all_goterms.keys():
                     all_goterms[term] = ""
                  if group == "acc":
                     accog_goterm[term] = str(int(accog_goterm.get(term, 0)) + 1)
                  elif group == "mixed":
                     mixedog_goterm[term] = str(int(mixedog_goterm.get(term, 0)) + 1)
                  elif group == "pcore":
                     coreog_goterm[term] = str(int(coreog_goterm.get(term, 0)) + 1)
                  elif group == "npcore":
                     coreog_goterm[term] = str(int(coreog_goterm.get(term, 0)) + 1)
                  else:
                     print("group error 2", file= sys.stderr)

      for term in all_goterms.keys():
         af = accog_goterm.get(term, 0)
         mf = mixedog_goterm.get(term, 0)
         cf = coreog_goterm.get(term, 0)
         tf = int(af) + int(mf) + int(cf)
#         print(portal, term, num_ogs, og_mixed_ct, og_core_ct, og_acc_ct, sep="\t")
         print(portal, term, tf, mf, cf, af, sep="\t")
   return 

def unique(list1):
   list_set = set(list1)
   return list(list_set)

def gofile2dict(filename, col_num):
   go_dict = {}
   with open(filename) as f:
      for line in f:
         elements = line.strip().split('\t')
         protid = "|".join(line.split("|")[0:3])
         goid = elements[col_num]
         go_dict[protid] = str(go_dict.get(protid, "")) + "\t" + str(goid)
   return go_dict

def get_portal_og_dict(filename):
   portal_og_dict = {}
   og_id_dict = {}
   og_group_dict = {}
   with open(filename) as f:
      for line in f:
         if len(line) == 0:
            continue
         parts = line.strip().split("\t")
         portal = parts[0]
         group = parts[1]
         og = parts[2]
         shortprotid = "|".join(parts[3].split("|")[0:3])
         if portal == "portal":
            continue
         portal_og_dict[portal] =str(portal_og_dict.get(portal, "")) + "\t" +  str(og)
         og_id_dict[og] =str(og_id_dict.get(og, "")) + "\t" +  str(shortprotid)
         og_group_dict[og] = str(group)
   return (portal_og_dict, og_id_dict, og_group_dict)

## requires sys package
if len(sys.argv) > 2:
   main(sys.argv[1], sys.argv[2])
else:
   print("improper usage: python3 tool.py ./results/2_portal_orthogroup_group_id.tsv ./all_functional_annotations.GO.tsv")
