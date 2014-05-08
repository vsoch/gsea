# Gene Set Enrichment Analysis - R
These scripts will be use to run gsea ,jar on Sherlock using R.  (currently in development!)

## Input Preparation

### parseCEL.R
NOTE: This script is intended to fix normalization for OLD analysis that have already been prepared (eg, the .gct files and .cls files and .chip files created), but not properly normalized! Raw data will be read from CEL files, background correction (PM and MM correction), and then quantile normalization.  Original probe names are in rows, and samples in columns.  The sample names should be the GSE IDs (extracted from the names of the cel files), and probe names will be the original probe names, and we run gsea with a "chip" file that has the gene--> probe lookup.  Again - this is ONLY for runs for which a .gct, .cls, and .chip already exist - a new script will be made to accomplish this task for completely new raw data. 
Is run by:

### run_parseCEL.R
in which you should specify the following input variables:
 - jobname: will be the name of the job submission
 - celdir: is the directory with input .CEL files
 - outdir: is the "norm" directory that normalized data, formatted for gsea, will go
 - olddir: directory with OLD formatted data (.gct) files are.  
 - rundir: the directory to output the inputParams.txt file, should also contain the "chip" and "cls" and "norm" folders.

The script will create a job file in the .jobs folder in the PWD, and both error and terminal output will go to .out.  The script will also submit the job to run on the sherlock cluster, with a max time of 2 days and 8GB memory.

### createTermGeneDatabase.R
Create a gene database from lists saved in Rda files.  This is customized to work with Vanessa's brain gene terms and would need to be modified for someone else's use.

### combineGeneSetsforPCA.R
Creates a compiled data frame for a specific subset of terms for performing PCA / visualization.  Again is customized for Vanessa's brain terms and would need customization for other use.

## Run GSEA

### run_gsea.R
This script will create .job files in the /.jobs folder and submit them to run gsea with the gsea java .jar file on sherlock.

This script takes as input:
- gseadir: the full path to the java .jar gsea file
- inputdata: read in from NORM column in inputParams.txt
- inputchip: read in from CHIP column in inputParams.txt
- inputcls: read in from CLASS column in inputParams.txt
- inputdb: a list specified at the top of the file, full paths to 1+ gene set databases
- outdir: the full path of top level folder (where gsea will create subdirectory)

You do not need to specify
- inputprefix: Is extracted from the input data path - this is the name of the subdirectory that will be created.

### TO DO

- parse other input types

## Output Parsing

### findSigOutput.py
NOTE: This is for GSEA R package only, not to be used with .java result
Should be run first after all gsea runs are complete.  A Python script that reads through files under a specific directory, parses the gsea R package output files, and reads in P values.  Does not correct for multiple comparisons, but just outputs all results into a table.

### findSigResults.R
NOTE: This is for GSEA R package only, not to be used with .java result
Should be run second, after findSigOutput.py, to threshold results and report significant under some p value.

### GOenrichment.R
Looks at network of gene functions for a gene term subset

### amiGOImages.R
Produces RamiGO images for a list of GO terms.  Currently does not function because the RamiGO package needs to be updated for changes to Gene Ontology.

### tanimotoAssess.R
Calculates tanimoto scores for gene sets in a database to assess overlap.  Tanimoto is set intersection / union

## Gene Conversion

### geneToKegg.R
Has code to convert from KO identifiers to swiss, or swiss to KO.
 
