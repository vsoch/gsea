# This is the cel file directory
# run_parseCEL.R
# This script will, for a folder of .CEL files, read in the files, sort 
# for interactive node: srun -n 13 -N 1 --mem=64000 --pty bash --time=24:00:00

jobname = "CEL.PARSE"
celdir = '/scratch/PI/dpwall/DATA/GENE_EXPRESSION/gsea/BPD/cel'
outdir = '/scratch/PI/dpwall/DATA/GENE_EXPRESSION/gsea/BPD/norm'
olddir = '/scratch/PI/dpwall/DATA/GENE_EXPRESSION/gsea/BPD/OLD/input'
topdir = '/scratch/PI/dpwall/DATA/GENE_EXPRESSION/gsea/BPD'

jobby = paste(jobname,".job",sep="")
sink(paste(".jobs/",jobby,sep=""))
cat("#!/bin/bash\n")
cat("#SBATCH --job-name=",jobby,"\n",sep="")  
cat("#SBATCH --output=.out/",jobby,".out\n",sep="")  
cat("#SBATCH --error=.out/",jobby,".err\n",sep="")  
cat("#SBATCH --time=2-00:00\n",sep="")
cat("#SBATCH --mem=8000\n",sep="")
cat("Rscript /scratch/PI/dpwall/SCRIPT/R/gsea/parseCEL.R",celdir,outdir,olddir,topdir,"\n")
sink()

# SUBMIT R SCRIPT TO RUN ON CLUSTER  
system(paste("sbatch",paste(".jobs/",jobby,sep="")))
