#!/bin/python

# This script will read through output files under a specified directory,
# parse the output files to read in P values and output a results report with 
# raw results

import os
import re

# Folder with gsea output folders to parse
indir = '/home/vanessa/Documents/Work/GENE_EXPRESSION/nextbio/ASD/gsea_output'
indir = '/home/vanessa/Documents/Work/GENE_EXPRESSION/nextbio/SCHIZO/gsea_output'
indir = '/home/vanessa/Documents/Work/GENE_EXPRESSION/nextbio/ALZ/gsea_output'
indir = '/home/vanessa/Documents/Work/GENE_EXPRESSION/nextbio/PTSD/gsea_output'
indir = '/home/vanessa/Documents/Work/GENE_EXPRESSION/nextbio/BPD/gsea_output'

# Significance threshold
thresh = 0.05

# Output folder for results report
outdir = '/home/vanessa/Documents/Work/GENE_EXPRESSION/nextbio/BPD/gsea_report/'

# Make the output directory if it doesn't exist
if not os.path.isdir(outdir):
  os.makedirs(outdir)

# Get list of folders in indir
folders = os.walk(indir).next()[1]

# Regular expression to find RESULTS report files
expression = re.compile("RESULTS")

# We will keep a huge compilation of all results
matrix = []
matrix.append("FOLDER\tGS\tSIZE\tSOURCE\tES\tNES\tNOM p-val\tFDR q-val\tFWER p-val\tTag /%\tGene /%\tSignal\tFDR (median)\tglob.p.val\n")

# For each folder, read in output SUMMARY files
for f in folders:
  # Get all files in folder
  files = os.listdir(indir + "/" + f)
  results = []
  for i in files:
    match = expression.search(i)
    if match:
      results.append(i)
  # Read in results files
  for r in results:
    filey = open(indir + '/' + f + '/' + r,'r')
    for m in filey.readlines():
      matrix.append(f+ "\t" + m) 


# Write to output file
outfile = open(outdir + "unthresh_matrix_all" + str(len(matrix)) + ".txt",'w')
outfile.writelines(matrix)
outfile.close()


# Now use findSigResults.R to threshold!
















  

