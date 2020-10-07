#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct  1 17:32:28 2020

@author: caitlincherryh
"""

# Version 3 from Suha

# Import libraries
import os
import numpy as np
from Bio.Nexus import Nexus
from Bio import AlignIO
from Bio.Alphabet import generic_protein
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment

## For AA alignments
if __name__ == '__main__': 
    # Assign rootDir <- the directory containing all directories to climb
    rootDir = '/Users/caitlincherryh/Documents/C1_EmpiricalTreelikeness/01_Data_Misof2014/'
    # Begin iterating through directories
    for DirName, subdirList, fileList in os.walk(rootDir):
        # If you can see a file called alignment.nex - i.e. if you can see an alignment
        if 'alignment.nex' in fileList:
            head_DirName, datas = os.path.split(DirName) # Break the last part of the path off (part after last slash)
            dset = os.path.basename(head_DirName) # Take the tail of the rootDir (file at end of directory name)
            print(DirName) # print the directory you're in
            aln_path = os.path.join(DirName,'alignment.nex') # create the path to a specific aligment
            aln = AlignIO.read(open(aln_path), "nexus") # open the nexus file as a multiple sequence alignment
            dat = Nexus.Nexus() # create a nexus object
            dat.read(aln_path) # read in the data as a nexus object
            aln_matrix = np.array([list(rec) for rec in aln]) # Iterate through each sequence in the aln and add them into a new matrix, turn it into an array
            taxa_list = dat.get_original_taxon_order() # extract the list of taxa
            unique_taxa_list = list(set(taxa_list)) # this is an unordered list of the taxa
            # Iterate through the loci/the charsets to separate into new files 
            for n in dat.charsets.keys():
                print(n)
                new_matrix = aln_matrix[:,dat.charsets[n]] # make a new matrix with only the characters for that charset
                new_aln = [] # initialise empty list
                for i in range(0, dat.ntax):
                    # For AA sequence:
                    new_aln.append(SeqRecord(Seq("".join(new_matrix[i,:]), generic_protein), id=taxa_list[i])) # append each row in the matrix as a new sequence
                new_msa = MultipleSeqAlignment(new_aln, generic_protein) # turn the list of SeqRecords into a MultipleSequenceAlignment object
                # AlignIO.write(MultipleSeqAlignment( new_aln[i] for i in range(0, len(taxa_list))), os.path.join(DirName,n+'.nex'), "nexus") #output nexus
                output_name = DirName+"/loci/"+n+".nex" # create a new output name based on the charset key
                AlignIO.write(new_msa,output_name,"nexus") # output this subset of alignment: write the msa as a nexus file
                
                
