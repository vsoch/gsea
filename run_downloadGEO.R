# This is the cel file directory
# run_parseCEL.R
# This script will, for a folder of .CEL files, read in the files, sort 
# for interactive node: srun -n 13 -N 1 --mem=64000 --pty bash --time=24:00:00

jobname = "GEO.PARSE"
# These are the GEO IDs
GEOSASD = c("GSE6575","GSE7329","GSE18123","GSE25507","GSE28475","GSE28521","GSE28475","GSE37772","GSE38322","GSE39447")
PTSDGEOS = c("GSE860","GSE42002")
ALZGEOS = c("GSE29378","GSE28146","GSE18309","GSE16759","GSE15222","GSE12685","GSE1297","GSE5281","GSE4757")
SZOGEOS = c("GSE53987","GSE12649","GSE12654","GSE12679","GSE17612","GSE21138","GSE21935","GSE25673","GSE46509","GSE38485","GSE27383","GSE37981","GSE35978","GSE26927")
BPDGEOS = c("GSE5388","GSE5329","GSE7036","GSE11767","GSE23848","GSE46449")

outdir = '/scratch/PI/dpwall/DATA/GENE_EXPRESSION/gsea/BPD/quant'
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
cat("Rscript /scratch/PI/dpwall/SCRIPT/R/gsea/downloadGEO.R",GEOS,outdir,olddir,topdir,"\n")
sink()

# SUBMIT R SCRIPT TO RUN ON CLUSTER  
system(paste("sbatch",paste(".jobs/",jobby,sep="")))
