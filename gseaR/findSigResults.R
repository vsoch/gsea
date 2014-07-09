# Get list of files in the output directory
disorder = "LYME"
outdir = paste("/scratch/PI/dpwall/DATA/GENE_EXPRESSION/gsea/",disorder,"/gsea",sep="")
folders = list.files(outdir)
# Threshold results at FDR q value:
threshold = 1
# This is the report directory to create and write results to
writedir = paste("/scratch/PI/dpwall/DATA/GENE_EXPRESSION/gsea/",disorder,"/report",sep="") 

reportFinal = c()
for (f in folders){
  reportdir = paste(outdir,"/",f,"/",sep="")
  visual = cat(list.files(reportdir,pattern="gsea_report_for"),sep="\n")
  report = list.files(reportdir,pattern=paste("gsea_report_for_MENIN",sep=""))
  if (length(report)!=0){
    idx = grep("*.xls",report)
    report = report[idx]
    report = read.csv(paste(reportdir,"/",report,sep=""),sep="\t",head=TRUE)
    report = report[report$FDR.q.val < threshold,]
    tmp = cbind(rep(f,dim(report)[1]),report)
    reportFinal = rbind(reportFinal,tmp)
   }
}
colnames(reportFinal)[1] = "NAME"
sort(reportFinal$FDR.q.val,decreasing=TRUE)

# Save to report folder
dir.create(writedir, showWarnings = FALSE)
write.table(reportFinal,file=paste(writedir,"/",disorder,"nothresh.txt",sep=""),sep="\t",row.names=FALSE)
