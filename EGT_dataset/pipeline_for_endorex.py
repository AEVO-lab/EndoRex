#!/usr/bin/env python3
# -*- coding: utf-8 -*-
###
###   Goal:
###      From MitoCOGs data (Kannan et al., 2014) to produce a dataset
###      for reconstructing Emdosymbiotic Gene Transfer evolutionary history
###      with EndoRex (Manuel Lafond -> https://github.com/AEVO-lab/EndoRex)
###
###   INPUT:
###      1- Gene/Orthologous families file
###         (dataset/11_Plantae_dataset/11_Plantae_MitoCOGs_nuc+mito_sorted_by_MitoCOGs.txt)
###      2- Species file
###         (dataset/11_Plantae_dataset/11_Plantae.txt)
###      3- Species tree file
###         (dataset/11_Plantae_dataset/11_Plantae_species_tree.nwk)
###      4- Sequences file for mitochondrial-encoded proteins of mitochondrial origins
###         (dataset/11_Plantae_dataset/MitoCOGs_sequences.fa)
###      5- Sequences file for nuclear-encoded proteins of mitochondrial origins
###         (dataset/11_Plantae_dataset/nuclear_encoded_orthologs_sequences_corrected.fa)
###      6- OUTPUT directory
###         (dataset/11_Plantae_dataset/RESULTS)
###
###   OUTPUT:
###      - INPUT gene trees for EndoRex
###        -> format example: ((((A__0__a1, A__1__a2), A__0__a3), A__0__a4), B__1__b1);
###           A,B             -> species name
###           0,1             -> gene encoding location: 0->mitochondria / 1->nucleus
###           a1,a2,a3,a4,b1  -> gene id
###
###   Name: pipeline_for_endorex.py           Author: Yoann Anselmetti
###   Creation date: 2021/01/05               Last modification: 2021/02/17
###


from sys import argv, path, exit, stdout
from os import close, path, makedirs
from shutil import rmtree
from datetime import datetime
import errno
import subprocess

from ete3 import Tree
### BIOPYTHON
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def mkdir_p(dir_path):
   try:
      makedirs(dir_path)
   except OSError as exc: # Python >2.5
      if exc.errno == errno.EEXIST and path.isdir(dir_path):
         pass
      else:
         raise


def set_record(ID,sep,location,dict_seq_record):
   tmpID=ID
   if not ID in dict_seq_record:
      tmpID="NA"+sep+ID.split(sep)[1]
   record=dict_seq_record[tmpID]

   SEQ=record.seq
   ID=ID.split(sep)[0]+sep+location+sep+ID.split(sep)[1]
   record=SeqRecord(SEQ,id=ID,description="")

   return record



if __name__ == '__main__':

   start_time = datetime.now()

   families_file=argv[1]
   species_file=argv[2]
   species_tree_file=argv[3]  
   mitocogs_seq=argv[4]
   nuclear_seq=argv[5]
   OUTPUT=argv[6]

   sep="__"
   verbose=1

   NOTUNG="../bin/Notung-2.9.1.5/Notung-2.9.1.5.jar"
   RAXML="../bin/raxmlHPC"
   MUSCLE="../bin/muscle"
   EndoRex="../../endorex.py"

   mkdir_p(OUTPUT)
   OUTPUT=path.abspath(OUTPUT)
   print("OUTPUT: "+OUTPUT)
   stdout.flush()

   output_fasta=OUTPUT+"/FASTA"
   output_msa=OUTPUT+"/MSA"
   raxml_matrix_name="PROTGAMMAAUTO"
   output_unrooted_tree=OUTPUT+"/TREE/RAxML/"+raxml_matrix_name
   output_rooted_tree=OUTPUT+"/TREE/NOTUNG"
   output_EGTrec_tree=OUTPUT+"/TREE/EndoRex"
   mkdir_p(output_fasta)
   mkdir_p(output_msa)
   if path.isdir(output_unrooted_tree):
      rmtree(output_unrooted_tree)
   mkdir_p(output_unrooted_tree)
   mkdir_p(output_rooted_tree)
   mkdir_p(output_EGTrec_tree)


   ### Get association between species ID and name
   in_species=open(species_file,"r")
   dict_species_id_name,dict_species_name_id=dict(),dict()
   for line in in_species:
      name=line.split(",")[0].replace(" ","_")
      ID=line.split(",")[1]
      dict_species_id_name[ID]=name
      dict_species_name_id[name]=ID


   ### Store MitoCOGs info in dict() "dict_ID"
   input_families=open(families_file,"r")
   dict_ID=dict()
   for line in input_families:
      seqID,ncbi_taxID,protein_length,MitoCOG_ID,loc=line.rstrip().split(",")
      ### Replace ncbi_taxID by species name
      species=dict_species_id_name[ncbi_taxID]
      if not MitoCOG_ID in dict_ID:
         dict_ID[MitoCOG_ID]=dict()
      dict_ID[MitoCOG_ID][species+sep+seqID]=loc





   ### GET SEQUENCES OF MITOCOGS FROM MITOCHONDRION
   mitocogs_dict=dict()
   for record in SeqIO.parse(mitocogs_seq, "fasta"):
      protein_id=record.id.split("|")[1].split()[0]
      ### Discard gene sequences that have been added to NCBI genes
      ### -> (lcl i.e. unamed/hypothetical protein)
      if "[" in record.description:
         species=record.description.split("[")[1].split("]")[0].replace(" ","_")
         if species in dict_species_name_id:
            mitocogs_dict[species+sep+protein_id]=record
      else:
         mitocogs_dict["NA"+sep+protein_id]=record


   ### GET SEQUENCES OF MITOCOGS FROM MITOCHONDRION
   nuclear_encoded_dict=dict()
   for record in SeqIO.parse(nuclear_seq, "fasta"):
      protein_id=record.id.split("|")[1].split()[0]
      ### Discard gene sequences that have been added to NCBI genes
      ### -> (lcl i.e. unamed/hypothetical protein)
      species=""
      if "[" in record.description:
         species=record.description.split("[")[1].split("]")[0].replace(" ","_")
         if species in dict_species_name_id:
            nuclear_encoded_dict[species+sep+protein_id]=record
      else:
         nuclear_encoded_dict["NA"+sep+protein_id]=record





   for MitoCOG_ID in dict_ID:
      if verbose:
         print("\nSTATS on mitochondrial protein-coding genes loation (mt or nuc):")
         print("- "+MitoCOG_ID+":")
      
      fasta_file=output_fasta+"/"+MitoCOG_ID+".fasta"
      out_fasta=open(fasta_file,"w")
      mt,nuc=0,0
      species_list,record_list=list(),list()
      for ID in dict_ID[MitoCOG_ID]:
         loc=dict_ID[MitoCOG_ID][ID]
         if loc=="mitochondrion":
            record_list.append(set_record(ID,sep,"0",mitocogs_dict))
            mt+=1
         elif loc=="nucleus":
            record_list.append(set_record(ID,sep,"1",nuclear_encoded_dict))
            nuc+=1
         else:
            exit("ERROR, location should be \"mitochondrion\" or \"nucleus\" and not: \""+loc+"\"")


         if verbose:
            species=ID.split(sep)[0]
            if not species in species_list:
               species_list.append(species)
      if verbose:
         print("\t+ #species: "+str(len(species_list)))
         print("\t+ #sequences: "+str(len(dict_ID[MitoCOG_ID])))
         print("\t\t* Mitochondria: "+str(mt))
         print("\t\t* Nucleus: "+str(nuc))

      ### Write FASTA file of the current MitoCOGs
      SeqIO.write(record_list,out_fasta,'fasta')
      out_fasta.close()

      if verbose:
         print("\n### START - Multiple Sequence Alignment (MSA) with MUSCLE")
      ### Run MUSCLE
      msa_file=output_msa+"/"+MitoCOG_ID+".afa"
      arg_list=[MUSCLE,"-in",fasta_file,"-out",msa_file]
      subprocess.check_call(arg_list)
      if verbose:
         print("### END - Multiple Sequence Alignment (MSA) with MUSCLE")


      ### Transform MSA file from FASTA to PHYLIP format
      if verbose:
         print("\n### START - Transform MSA FASTA format file in PHYLIP format file")
      records=SeqIO.parse(msa_file,"fasta")
      output_phylip=output_msa+"/"+MitoCOG_ID+".phylip"
      count=SeqIO.write(records,output_phylip,"phylip-relaxed") 
      if verbose:
         print("### END - Transform MSA FASTA format file in PHYLIP format file")

      ### Gene tree inference with raxmlHPC
      if verbose:
         print("\n### START - Protein-coding gene tree inference with RAxML")
      arg_list=[RAXML,"-m",raxml_matrix_name,"-p","123456","-s",output_phylip,"-w",output_unrooted_tree,"-n",MitoCOG_ID]
      subprocess.check_call(arg_list)
      if verbose:
         print("### END - Protein-coding gene tree inference with RAxML")

      ### Gene tree rooting with NOTUNG
      if verbose:
         print("\n### START - Protein-coding gene tree rooting with NOTUNG")
      gene_tree=output_unrooted_tree+"/RAxML_bestTree."+MitoCOG_ID
      arg_list=["java","-jar",NOTUNG,"-g",gene_tree,"-s",species_tree_file,"--root","--nolosses","--treeoutput","newick","--outputdir",output_rooted_tree]
      subprocess.check_call(arg_list)
      if verbose:
         print("### END - Protein-coding gene tree rooting with NOTUNG")


      ### EGT-reconciliation with EndoRex
      if verbose:
         print("\n### START - EGT reconciliation with EndoRex")
      gene_tree_file=output_rooted_tree+"/RAxML_bestTree."+MitoCOG_ID+".rooting.0"
      gene_tree_str=""
      in_gene_tree=open(gene_tree_file,'r')
      gene_tree_str=str(in_gene_tree.read().strip())
      in_gene_tree.close()
      # print("GENE_TREE: "+gene_tree_str)
      # print("TYPE GENE_TREE: ",type(gene_tree_str))
      in_species_tree=open(species_tree_file,'r')
      species_tree_str=str(in_species_tree.read().strip())
      in_species_tree.close()
      # print("SPECIES_TREE: "+species_tree_str)
      # print("TYPE SPECIES_TREE: ",type(species_tree_str))

      endorex_settings_list=list()
      # S1: Default setting
      endorex_settings_list.append(["S1","1.0","1.0","1.0","1.0","1.0","1.0"])
      # S2: Duplication & Loss costs used in NOTUNG rooting 
      endorex_settings_list.append(["S2","1.5","1.0","1.0","1.0","1.0","1.0"])
      # S3: S2 + higher costs for EGT (+1.0) 
      endorex_settings_list.append(["S3","1.5","1.0","2.0","2.0","2.0","2.0"])
      # S4: S3 + higher costs for EGTr (+1.0)
      endorex_settings_list.append(["S4","1.5","1.0","3.0","3.0","2.0","2.0"])
      # S5: S4 + higher costs for EGT orientation 1->0 (+1.0)
      endorex_settings_list.append(["S5","1.5","1.0","3.0","4.0","2.0","3.0"])
      # S6: S2 + higher costs for EGT orientation 1->0 (+1.0)
      endorex_settings_list.append(["S6","1.5","1.0","2.0","3.0","2.0","3.0"])


      for endorex_setting in endorex_settings_list:
         ID,dup,loss,fzero,fone,pzero,pone=endorex_setting         
         setting_dir="EndoRex_"+ID
         arg_list=["python3",EndoRex,"-s",species_tree_str,"-g",gene_tree_str,"-lm","none", \
                   "-d",dup,"-l",loss, \
                   "-f0",fzero,"-f1",fone, \
                   "-p0",pzero,"-p1",pone,]
         result=subprocess.run(arg_list,capture_output=True)

         endorex_dir=output_EGTrec_tree+"/no_loss_mode/"+setting_dir
         mkdir_p(endorex_dir)
         output_egt_tree_file=endorex_dir+"/"+MitoCOG_ID+"_DELrec.tree"
         if verbose:
            print("ENDOREX command line for setting "+ID+":")
            print("\t",end="")
            print(*arg_list)
            print("ENDOREX tree:")
            print("\t"+result.stdout.decode("utf-8"))
         with open(output_egt_tree_file,'w') as f:
            f.write(result.stdout.decode("utf-8"))
      if verbose:
         print("### END - EGT reconciliation with EndoRex")


   end_time = datetime.now()
   print('\nDuration: {}'.format(end_time - start_time))
