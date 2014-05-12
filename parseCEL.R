# GEOParse: Script 1
# This script will, for a folder of .CEL files, read in the files, sort 
# for interactive node: srun -n 12 -N 1 --mem=64000 --pty bash --time=24:00:00
# USAGE: RSCRIPT parseCEL.R celdir outdir rundir olddir
# RSCRIPT parseCEL.R /scratch/PI/dpwall/DATA/GENE_EXPRESSION/gsea/ASD/cel /scratch/PI/dpwall/DATA/GENE_EXPRESSION/gsea/ASD/norm /scratch/PI/dpwall/DATA/GENE_EXPRESSION/gsea/ASD/OLD/input /scratch/PI/dpwall/DATA/GENE_EXPRESSION/gsea/ASD'

library('oligo')
library('affy')

# VSochat April 2014

# SZO GSE25673 is NA, GSE38485 is bxg
# BPD is ready
# SZO, ALZ are done for cel
args <- commandArgs(TRUE)
celdir = args[1]
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
fileConn<-file(paste(topdir,"/inputParam.txt",sep=""))
writeLines("NORMDATA\tCHIP\tCLASS", fileConn)
close(fileConn)

# Get list of folders in ASD/cel directory
folders = list.files(celdir,pattern="GSE*")

# For each folder, format files
for (f in folders){
  folder = paste(celdir,'/',f,sep='')
  setwd(folder)
  files = list.files(paste(celdir,'/',f,sep=''),pattern="*.CEL")
  # Read in cel files
  affy.data = ReadAffy(filenames = files)
  # Summarize and normalize with MAS5
  eset.mas5 = mas5(affy.data)
  exprSet.nologs = exprs(eset.mas5)
  # Another option
  #affy.data <- read.celfiles(files)
  #exprSet.nologs = rma(affy.data)
  #expr = exprs(exprSet.nologs)
  
  # Get the sample names
  samp <- sampleNames(affy.data)
  samp = gsub(".cel.gz","",samp)
  samp = gsub(".CEL.gz","",samp)
  colnames(exprSet.nologs) = samp
  
  # Read in the chip file
  chipfile = paste(topdir,"/chip/",f,".chip",sep="")
  chipdata = read.csv(chipfile,head=TRUE,sep="\t")
  # Get rid of empty genes
  chipdata = chipdata[!is.na(chipdata[,2]),]
  chipdata = chipdata[-which(chipdata[,2]==""),]
  # Write to file
  chipfile = gsub("[.]chip","_filt.chip",chipfile)
  colnames(chipdata) = c("Probe Set ID","Gene Symbol","Gene Title")
  write.table(chipdata,file=chipfile,sep="\t",col.names=TRUE,row.names=FALSE,quote=FALSE)
  
  # For each of the old analyses, get GSMIDs from the files
  oldies = list.files(olddir,pattern=paste(f,"*",sep=""))
  for (o in oldies) {
    if (substring(o,nchar(o)-3,nchar(o)) == ".gct"){
      # Read in the data file
      # Get first two rows
      #firsttwo = readLines(paste(olddir,"/",o,sep=""),n=2)
      raw = read.csv(paste(olddir,"/",o,sep=""),skip=2,sep="\t")
      subid = colnames(raw)[3:ncol(raw)]
      filtered = exprSet.nologs[,subid]
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