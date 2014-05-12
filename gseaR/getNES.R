# This script will read in a report file (generated with findSigResults.R
# and a specified threshold, and output another file in the report
# directory with the NES score and genes for each significant gene set

# This is the disorder folder
disorder = "SZO"
# This is the output result file
result = paste("/scratch/PI/dpwall/DATA/GENE_EXPRESSION/gsea/",disorder,"/report/",disorder,".txt",sep="")
# This is the file we will write to parse for gephi input
gephi = paste("/scratch/PI/dpwall/DATA/GENE_EXPRESSION/gsea/",disorder,"/report/",disorder,"_gephiLookup.txt",sep="")
# This is the output directory for gsea
outdir = paste("/scratch/PI/dpwall/DATA/GENE_EXPRESSION/gsea/",disorder,"/output/",sep="")

# FDR corrected threshold
threshold = .05

# Read in the file
result = read.csv(result,head=TRUE,sep="\t")
result = result[result$FDR.q.val <= threshold,]

# These are the unique reports that we need to go through
reports = unique(result$NAME)

output = c()
# For each result above the threshold, save the FDRq value, gene set name, and NES
for (r in 1:length(reports)[1]) {
  folder = as.character(reports[r])
  tmp = result[result$NAME == folder,]
  for (g in 1:dim(tmp)[1]){
    geneset = tmp[g,]    
    # Read in each gene set file
    listfile = paste(outdir,folder,"/",as.character(geneset$NAME.1),".xls",sep="")
    listfile = read.table(listfile,head=TRUE,sep="\t")
    if (dim(listfile)[2] == 9){
      listfile = cbind(rep(folder,dim(listfile)[1]),rep(geneset$NAME.1,dim(listfile)[1]),rep(geneset$SIZE,dim(listfile)[1]),rep(geneset$NES,dim(listfile)[1]),rep(geneset$FDR.q.val,dim(listfile)[1]),listfile[,c(2,5,6,7,8)])
    }else {
      listfile = cbind(rep(folder,dim(listfile)[1]),rep(geneset$NAME.1,dim(listfile)[1]),rep(geneset$SIZE,dim(listfile)[1]),rep(geneset$NES,dim(listfile)[1]),rep(geneset$FDR.q.val,dim(listfile)[1]),listfile[,c(2,6,7,8,9)])
    }
    output = rbind(output,listfile)
  }  
}
output = as.data.frame(output)
colnames(output) = c("FOLDER","GENE_SET","SIZE","NES","FDR_Q_VALUE","GENE","RANK_IN_GENE_LIST","RANK_METRIC_SCORE","RUNNING_ES","CORE_ENRICHMENT")

# Print to file
write.table(output,file=gephi,row.names=FALSE,sep="\t")