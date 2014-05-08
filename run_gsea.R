# run_gsea.R will read in parameters from an input file "inputParams.txt," create
# submission scripts to run the java version of GSEA on sherlock, and submit them.
# USAGE RSCRIPT gsea.R gseadir inputprefix inputdata inputchip inputcls inputdb outdir

gseadir = "/share/PI/dpwall/SOFTWARE/GSEA-P-R/gsea2-2.0.14.jar"                  # Path to main GSEA program
inputfile = "/scratch/PI/dpwall/DATA/GENE_EXPRESSION/gsea/ASD/inputParam.txt"    # Path to input parameter file
inputdb = c("/scratch/PI/dpwall/DATA/GENE_EXPRESSION/gsea/GENE_DATABASE/brainTerms.gmt","/scratch/PI/dpwall/DATA/GENE_EXPRESSION/gsea/GENE_DATABASE/ASD.gmt")
outdirtop = "/scratch/PI/dpwall/DATA/GENE_EXPRESSION/gsea/ASD/output"             # Top level output directory - subdirectories will be made inside
setwd('/scratch/PI/dpwall/SCRIPT/R/gsea')

# Read in input parameter file - create job script and submit for each entry
inputfile = read.csv(inputfile,sep="\t",head=TRUE)
for (i in 1:dim(inputfile)[1]){
  normdata = as.character(inputfile$NORMDATA[i])
  inputchip = inputfile$CHIP[i]
  inputcls = inputfile$CLASS[i]
  
  for (db in inputdb){
    
    dbname = strsplit(db,'/')[[1]]
    dbname = gsub('.gmt','',dbname[length(dbname)])
  
    folder = strsplit(as.character(normdata),"/")[[1]]
    folder = folder[length(folder)]
    folder = paste(gsub(".gct","",folder),"_",dbname,sep="")
    inputprefix = folder
    #outdir = paste(outdirtop,"/",folder,"/",sep="")
    #dir.create(outdir, showWarnings = FALSE)
    
    jobby = paste(folder,".job",sep="")
    sink(file=paste(".jobs/",jobby,sep=""))
    cat("#!/bin/bash\n")
    cat("#SBATCH --job-name=",jobby,"\n",sep="")  
    cat("#SBATCH --output=.out/",jobby,".out\n",sep="")  
    cat("#SBATCH --error=.out/",jobby,".err\n",sep="")  
    cat("#SBATCH --time=2-00:00\n",sep="")
    cat("#SBATCH --mem=8000\n",sep="")
    cat("java -cp",gseadir,"xtools.gsea.Gsea -res",normdata,"-cls",as.character(inputcls),"-gmx",db,"-chip",as.character(inputchip),"-collapse true -mode Max_probe -norm meandiv -nperm 1000 -permute phenotype -rnd_type no_balance -scoring_scheme weighted -rpt_label",inputprefix,"-metric Signal2Noise -sort real -order descending -include_only_symbols true -make_sets true -median false -num 100 -plot_top_x 20 -rnd_seed timestamp -save_rnd_lists false -set_max 500 -set_min 15 -zip_report false -out",outdirtop,"-gui false\n")
    #cat("Rscript /scratch/PI/dpwall/SCRIPT/R/gsea/gsea.R",gseadir,inputprefix,normdata,as.character(inputchip),as.character(inputcls),db,outdir,"\n")
    sink()
    
    # SUBMIT R SCRIPT TO RUN ON CLUSTER  
    #system(paste("sbatch",paste(".jobs/",jobby,sep="")))
  }  
}
