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
BRDGEOS = c("GSE52222")
DOWNGEOS = c("GSE5390","GSE55426","GSE50586","GSE42114","GSE47014","GSE48611","GSE42956","GSE42772","GSE35665","GSE16176","GSE9762","GSE5390","GSE6408","GSE6283","GSE1789")
TURNERGEOS = c("GSE46687","GSE22551")
SLEEPGEOS = c("GSE48113","GSE49800","GSE38792","GSE37667","GSE21592")
ALCGEOS = c("GSE53808","GSE49393","GSE49376","GSE44456","GSE25999","GSE20568","GSE10356","GSE3632")
MDDGEOS = c("GSE12654","GSE24095","GSE58105","GSE54572","GSE54571","GSE54570","GSE54575","GSE54568")
FXGEOS = c("GSE48903","GSE48902","GSE48873","GSE41273")
INTDIS = c("GSE23358")
ADRENOGEOS = c("GSE34309")
ALSGEOS = c("GSE3307","GSE33855","GSE56808","GSE51684","GSE51741","GSE39543","GSE39642","GSE28253","GSE26927","GSE21450","GSE18920","GSE7950","GSE4595","GSE3307")
disorder = "LD"

outdir = paste('/scratch/PI/dpwall/DATA/GENE_EXPRESSION/gsea/',disorder,'/quant',sep="")
#olddir = '/scratch/PI/dpwall/DATA/GENE_EXPRESSION/gsea/BPD/OLD/input'
topdir = paste("/scratch/PI/dpwall/DATA/GENE_EXPRESSION/gsea/",disorder,sep="")

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
