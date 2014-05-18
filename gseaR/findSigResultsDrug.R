# Get list of files in the output directory
outdir = "/scratch/PI/dpwall/DATA/DRUG/gsea"
folders = list.files(outdir)
# Threshold results at FDR q value:
threshold = .25
# This is the report directory to create and write results to
writedir = "/scratch/PI/dpwall/DATA/DRUG/report/" 

reportFinal = c()
for (f in folders){
  reportdir = paste(outdir,"/",f,"/",sep="")
  report = list.files(reportdir,pattern="gsea_report_for_DRUG")
  idx = grep("*.xls",report)
  report = report[idx]
  report = read.csv(paste(reportdir,"/",report,sep=""),sep="\t",head=TRUE)
  report = report[report$FDR.q.val < threshold,]
  tmp = cbind(rep(f,dim(report)[1]),report)
  reportFinal = rbind(reportFinal,tmp)
}
colnames(reportFinal)[1] = "NAME"

# Save to report folder
dir.create(writedir, showWarnings = FALSE)
write.table(reportFinal,file=paste(writedir,"/DRUGS_GSEA.txt",sep=""),sep="\t",row.names=FALSE)