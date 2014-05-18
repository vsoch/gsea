# cmapGSEA will run gsea with the brainterms set for all sets of medications in
# the connectivity map.  The top section runs "cmapGSEA.R" for each unique drug in the
# connectivity Map database, which prepares the input data for GSEA.  The bottom half of
# the script runs the GSEA using the gsea .jar file=, also by way of submitting jobs
# to the Sherlock cluster.

# VSochat May 2014
require(gdata)
library('stringr')

# Read in file with instance information
cmap = '/scratch/PI/dpwall/DATA/DRUG/CONNECTIVITY_MAP'
instances = paste(cmap,'/instances/',sep="")
meta = read.xls(paste(cmap,"/cmap_instances_02.xls",sep=""),head=TRUE,sep=",")
chipdir = paste(cmap,"/chip",sep="")
setwd('/scratch/PI/dpwall/SCRIPT/R/gsea')

# GSEA Paths
gseadir = "/share/PI/dpwall/SOFTWARE/GSEA-P-R/gsea2-2.0.14.jar"                  # Path to main GSEA program
inputdb = "/scratch/PI/dpwall/DATA/GENE_EXPRESSION/gsea/GENE_DATABASE/brainTerms.gmt"

# Output directory will be here
outdir = paste(cmap,"/gsea",sep="")
datadir = paste(cmap,"/input",sep="")
clsdir = paste(cmap,"/cls",sep="")
dir.create(outdir, showWarnings = FALSE)
dir.create(datadir, showWarnings = FALSE)
dir.create(clsdir, showWarnings = FALSE)

# Get rid of empty ones (silly excel sheets)
meta = meta[-which(meta$cmap_name==""),]

# For each unique drug, create and run a script to perform gsea
drugs = as.character(unique(meta$cmap_name))

# We will first parse all data
for (d in drugs){
  drug = meta[which(meta$cmap_name==d),]
  
  # Find the unique chips
  counts = table(drug$array3)
  chips = names(counts[which(counts != 0)])
  chips = gsub("-","",chips)

  for (c in chips){
    tmp = drug[drug$array3==c,]
    # These are the drug exposed cell lines
    pscan = as.character(tmp$perturbation_scan_id)
    # These are unique verhicles (the controls)
    vscan = as.character(tmp$vehicle_scan_id4)
    # Loop through vehicle scans, find the multiple files
    multipleIDX = c()
    for (v in 1:length(vscan)){
      vehicle = vscan[v]
      if (substr(vehicle,1,1) == "."){
        multipleIDX = c(multipleIDX,v)
      }      
    }
    if (!is.null(multipleIDX)){
      vscan2 = c()
      for (m in multipleIDX){
        tosplit = vscan[m]
        extension = strsplit(tosplit,"[.]")[[1]]
        extension = extension[extension!=""]
        probe = strsplit(pscan[m],"[.]")[[1]][1]
        for (e in extension){
          vscan2 = c(vscan2,paste(probe,".",e,sep=""))        
        }
      }
      # Get rid of the multiple IDS
      vscan = vscan[-multipleIDX]
      # Add the controls
      vscan = c(vscan,vscan2)
      
      # Form the full paths
      pscan = sapply(instances,paste,sapply(pscan,paste,".CEL.bz2",sep=""),sep="")
      vscan = sapply(instances,paste,sapply(vscan,paste,".CEL.bz2",sep=""),sep="")
      
      # Make sure that every file exists
      ps = c()
      vs = c()
      for (p in pscan){
        if (file.exists(p)){
          ps = c(ps,p)
        }
      }
      for (v in vscan){
        if (file.exists(v)){
          vs = c(vs,v)
        }
      }
    } else {
      ps = sapply(instances,paste,sapply(pscan,paste,".CEL.bz2",sep=""),sep="")
      vs = sapply(instances,paste,sapply(vscan,paste,".CEL.bz2",sep=""),sep="")
    }    
    nump = length(ps) 
    numc = length(vs)
    pscan = paste(ps,collapse=",")
    vscan = paste(vs,collapse=",")
    # Combine into one list
    is = paste(pscan,vscan,sep=",")
    
    # Format drug name for file and script submission
    name = as.character(d)
    name = str_replace_all(name,"[[:punct:]]","")
    name = gsub(" ","",name)
    
    # Create class file
    classfile = paste(clsdir,"/",name,"_",c,".cls",sep="")
    cat(nump + numc,"2 1\n",file=classfile)
    cat("# DRUG CON\n",file=classfile,append=TRUE)
    cat(rep("1",nump),rep("2",numc),"\n",file=classfile,append=TRUE)
    
    # Create job file and submit
    if (!file.exists(paste(datadir,"/",name,"_",c,".gct",sep=""))) {
      jobby = paste(substr(name,1,3),c,"pre.job",sep="")
      sink(file=paste(".jobs/",jobby,sep=""))
      cat("#!/bin/bash\n")
      cat("#SBATCH --job-name=",jobby,"\n",sep="")  
      cat("#SBATCH --output=.out/",jobby,".out\n",sep="")  
      cat("#SBATCH --error=.out/",jobby,".err\n",sep="")  
      cat("#SBATCH --time=2-00:00\n",sep="")
      cat("#SBATCH --mem=8000\n",sep="")
      cat("Rscript /scratch/PI/dpwall/SCRIPT/R/gsea/cmapPrep.R",name,is,datadir,c)  
      sink()
      system(paste("sbatch",paste(".jobs/",jobby,sep="")))
    }
  }
}

drugs = list.files("/scratch/PI/dpwall/DATA/DRUG/CONNECTIVITY_MAP/input", pattern=".gct")

# Now (when the above data prep is done) we will run GSEA
topdir = "/scratch/PI/dpwall/DATA/DRUG/CONNECTIVITY_MAP/"
for (d in drugs){
  name = as.character(d)
  chip = strsplit(d,"_")[[1]]
  chip = gsub(".gct","",paste(chip[2:length(chip)],collapse=""))
  chip = paste(topdir,"chip/",chip,".chip",sep="")
  input = paste(topdir,"input/",name,sep="")
  cls = paste(topdir,"cls/",gsub(".gct",".cls",name),sep="")
  
  jobby = paste(d,"_gsea.job",sep="")
  sink(file=paste(".jobs/",jobby,sep=""))
  cat("#!/bin/bash\n")
  cat("#SBATCH --job-name=",jobby,"\n",sep="")  
  cat("#SBATCH --output=.out/",jobby,".out\n",sep="")  
  cat("#SBATCH --error=.out/",jobby,".err\n",sep="")  
  cat("#SBATCH --time=2-00:00\n",sep="")
  cat("#SBATCH --mem=8000\n",sep="")
  cat("java -cp",gseadir,"xtools.gsea.Gsea -res",input,"-cls",as.character(cls),"-gmx",inputdb,"-chip",as.character(chip),"-collapse true -mode Max_probe -norm meandiv -nperm 1000 -permute genes -rnd_type no_balance -scoring_scheme weighted -rpt_label",gsub(".gct","",d),"-metric Signal2Noise -sort real -order descending -include_only_symbols true -make_sets true -median false -num 100 -plot_top_x 20 -rnd_seed timestamp -save_rnd_lists false -set_max 500 -set_min 15 -zip_report false -out",outdir,"-gui false\n")
  sink()
  
  system(paste("sbatch",paste(".jobs/",jobby,sep="")))
}
