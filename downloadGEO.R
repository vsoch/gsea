# GEOParse: Script 1
# This script will, for a folder of .CEL files, read in the files, sort 
# for interactive node: srun -n 12 -N 1 --mem=64000 --pty bash --time=24:00:00
# USAGE: RSCRIPT parseCEL.R celdir outdir rundir olddir
# RSCRIPT parseCEL.R /scratch/PI/dpwall/DATA/GENE_EXPRESSION/gsea/ASD/cel /scratch/PI/dpwall/DATA/GENE_EXPRESSION/gsea/ASD/norm /scratch/PI/dpwall/DATA/GENE_EXPRESSION/gsea/ASD/OLD/input /scratch/PI/dpwall/DATA/GENE_EXPRESSION/gsea/ASD'

library('oligo')
library('affy')
library('GEOquery')
library(Biobase)
library('preprocessCore')  # for quantile normalization
library("hgu95av2.db")     # for annotation


# VSochat April 2014

args <- commandArgs(TRUE)
GEOS = args[1]
outdir = args[2]
topdir = args[3]
olddir = args[4]

# Create output directory, if it doesn't exist
dir.create(outdir, showWarnings = FALSE)

cat(celdir,":cels\n")
cat(outdir,":output\n")
cat(topdir,":topdir\n")
cat(olddir,":oldgctdir\n\n")

# This is a log file, used to keep track of parameters to run GSEA
# ONLY RUN BELOW IF LOGFILE IS NOT YET CREATED
fileConn<-file(paste(topdir,"/inputParam.txt",sep=""))
writeLines("NORMDATA\tCHIP\tCLASS", fileConn)
close(fileConn)

# For each GEO, download data
for (g in GEOS){
  data = getGEO(g)
  e = data[[1]]
  dat = exprs(e) # Here is the expression matrix
  pheno = pData(e) # Here is the phenotypic data
  genes = fData(e) # This would be data on probes

  # 1. NORMALIZATION --------------------------------------------
  # Look at MA plot
  ma.plot( rowMeans(log2(dat)), log2(dat[, 1])-log2(dat[, 2]), cex=1 )
  
  # IF WE NEED TO QUANTILE NORMALIZE
  norm = normalize.quantiles(dat)    
  ma.plot( rowMeans(log2(norm)), log2(norm[, 1])-log2(norm[, 2]), cex=1 )
  
  # IF DATA ALREADY NORMALIZED (look at pheno$data_processing)
  norm = dat
  # -------------------------------------------------------------
  
  rownames(norm) = rownames(e)
  colnames(norm) = colnames(e)
  
  # 2. MAKE CHIP FILE --------------------------------------------
  # A) If chip file already exists, just specify path
  chipfile = paste(topdir,"/chip/",g,"_filt.chip",sep="")
  file.exists(chipfile)
  
  # B) If chip file already exists and we need to "filter"
  chipfile = paste(topdir,"/chip/",g,"_series_matrix_filt.chip",sep="")
  chipdata = read.csv(chipfile,head=TRUE,sep="\t")
  # Get rid of empty genesdi
  chipdata = chipdata[!is.na(chipdata[,2]),]
  chipdata = chipdata[-which(chipdata[,2]==""),]
  # Write to file
  chipfile = gsub("[.]chip","_filt.chip",chipfile)
  colnames(chipdata) = c("Probe Set ID","Gene Symbol","Gene Title")
  write.table(chipdata,file=chipfile,sep="\t",col.names=TRUE,row.names=FALSE,quote=FALSE)  
  # C) If not, then specify what new path will be
  # If files do not contain gene probe names, need to look up
  chipfile = paste(topdir,"/chip/",g,"_filt.chip",sep="")
  # NO LOOKUP
  tmp=cbind(as.character(genes$ID),as.character(genes$Symbol),as.character(genes$Definition))
  # IF NEED TO LOOKUP Map from accession numbers to gene symbols
  sym = select(hgu95av2.db, as.character(genes$GB_ACC), "SYMBOL", "ACCNUM")
  # Also get description of gene
  name = select(hgu95av2.db, as.character(genes$GB_ACC), "GENENAME", "ACCNUM")
  # Match probe IDs to ACCNUM
  idx = match(sym$ACCNUM,genes$GB_ACC)
  probes = genes$ID[idx]
  # Merge to make chip file fields
  tmp = cbind(as.character(probes),sym$SYMBOL,name$GENENAME)
  colnames(tmp) = c("Probe Set ID","Gene Symbol","Gene Title")
  write.table(tmp,file=chipfile,sep="\t",col.names=TRUE,row.names=FALSE,quote=FALSE)
  # -------------------------------------------------------------
  
  # 3. MAKE CLASS FILE --------------------------------------------
  classfile = paste(topdir,"/cls/",g,".cls",sep="")
  labels = as.character(pheno$characteristics_ch1.2)
  unilabels = unique(labels)
  ASD = c(1)  
  CON = c(2)
  labelidx = which(labels %in% unilabels[c(ASD,CON)])
  subset = labels[labelidx]
  binlabels = array(dim=length(subset))
  binlabels[subset %in% unilabels[ASD]] = 1
  binlabels[subset %in% unilabels[CON]] = 2
  # Print to class file
  cat(c(length(binlabels),2,1,collapse=" "),"\n",file=classfile)
  cat("# BPD HC\n",append=TRUE,file=classfile)
  cat(binlabels,file=classfile,append=TRUE)
  # Filter data to include these
  norm = norm[,labelidx]
  # --------------------------------------------------------------
  
  # 4. MAKE DATA FILE --------------------------------------------
  outfile = paste(outdir,"/",g,"_BPDvsHC.gct",sep="")
  DESCRIPTION = rep('na',dim(norm)[1])
  NAME = rownames(norm)
  filtered = cbind(NAME,DESCRIPTION,norm)
  # Now print new file
  dimprint = paste(dim(filtered)[1],dim(filtered)[2]-2,sep="\t")
  fileConn = file(outfile)
  writeLines(c("#1.2",dimprint),sep="\n",fileConn)
  close(fileConn)  
  write.table(filtered,file=outfile,append=TRUE,col.names=TRUE,row.names=FALSE,quote=FALSE,sep="\t")   
  # -------------------------------------------------------------
  
  # Lastly, we want to document the input files for each, so we can run programatically
  cat(paste(outfile,chipfile,classfile,sep="\t"),file=paste(topdir,"/inputParam.txt",sep=""),append=TRUE,sep="\t")
  cat("\n",file=paste(topdir,"/inputParam.txt",sep=""),append=TRUE)
  
  # FOR OLD ANALYSES ONLY (to get groups)  
  # For each of the old analyses, get GSMIDs from the files
  oldies = list.files(olddir,pattern=paste(g,"*",sep=""))
  for (o in oldies) {
    if (substring(o,nchar(o)-3,nchar(o)) == ".gct"){
      # Read in the data file
      # Get first two rows
      firsttwo = readLines(paste(olddir,"/",o,sep=""),n=2)
      raw = read.csv(paste(olddir,"/",o,sep=""),skip=2,sep="\t")
      subid = colnames(raw)[3:ncol(raw)]
      filtered = norm[,subid]
      o = gsub("-","",o)
      classfile = paste(topdir,"/cls/",gsub(".gct",".cls",o),sep="")
      
      # Add probe names, description(with na)
      DESCRIPTION = rep('na',dim(filtered)[1])
      NAME = rownames(filtered)
      filtered = cbind(NAME,DESCRIPTION,filtered)
      
      # Now print new file
      outfile = paste(outdir,"/",gsub("series_matrix_","",o),sep="")
      dimprint = paste(dim(filtered)[1],dim(filtered)[2]-2,sep="\t")
      fileConn = file(outfile)
      writeLines(c("#1.2",dimprint),sep="\n",fileConn)
      close(fileConn)  
      write.table(filtered,file=outfile,append=TRUE,col.names=TRUE,row.names=FALSE,quote=FALSE,sep="\t")   
    
      # Lastly, we want to document the input files for each, so we can run programatically
      cat(paste(outfile,chipfile,classfile,sep="\t"),file=paste(topdir,"/inputParam.txt",sep=""),append=TRUE,sep="\t")
      cat("\n",file=paste(topdir,"/inputParam.txt",sep=""),append=TRUE)
    } 
  }
  
  # We already have old chip files that we can use - moved into .chip folder
}