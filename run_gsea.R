# run_gsea.R will read in parameters from an input file "inputParams.txt," create
# submission scripts to run gsea.R on sherlock, and submit them.
# USAGE RSCRIPT gsea.R gseadir inputprefix inputdata inputchip inputcls inputdb outdir

gseadir = "/share/PI/dpwall/SOFTWARE/GSEA-P-R/GSEA.1.0.R"                         # Path to main GSEA program
inputfile = "/scratch/PI/dpwall/DATA/GENE_EXPRESSION/gsea/ASD/inputParams.txt"    # Path to input parameter file
inputdb = c("/share/PI/dpwall/SOFTWARE/GSEA-P-R/GeneSetDatabases/brainTerms.gmt","/share/PI/dpwall/SOFTWARE/GSEA-P-R/GeneSetDatabases/ASD.gmt")
outdirtop = "/scratch/PI/dpwall/DATA/GENE_EXPRESSION/gsea/ASD/output"             # Top level output directory - subdirectories will be made inside

# Read in input parameter file - create job script and submit for each entry
inputfile = read.csv(inputfile,sep="\t",head=TRUE)
for (i in 1:dim(inputfile)[1]){
  normdata = inputfile$NORMDATA[i]
  folder = strsplit(normdata,"/")[[1]]
  folder = folder[length(folder)]
  folder = gsub(".gct","",folder)
  inputprefix = gsub(".gct","",strsplit(i,"_")[[1]][2])
  inputchip = inputfile$CHIP[i]
  inputcls = inputfile$CLASS[i]
  outdir = paste(outdirtop,"/",folder)
  dir.create(outdir, showWarnings = FALSE)
  
  for (db in inputdb){
    
    jobby = paste(i,".job",sep="")
    sink(paste(".jobs/",jobby,sep=""))
    cat("#!/bin/bash\n")
    cat("#SBATCH --job-name=",jobby,"\n",sep="")  
    cat("#SBATCH --output=.out/",jobby,".out\n",sep="")  
    cat("#SBATCH --error=.out/",jobby,".err\n",sep="")  
    cat("#SBATCH --time=2-00:00\n",sep="")
    cat("#SBATCH --mem=8000\n",sep="")
    cat("Rscript /scratch/PI/dpwall/SCRIPT/R/gsea/gsea.R",gseadir,normdata,inputchip,inputcls,db,outdir,"\n")
    sink()
    
    # SUBMIT R SCRIPT TO RUN ON CLUSTER  
    system(paste("sbatch",paste(".jobs/",jobby,sep="")))
  }  
}
