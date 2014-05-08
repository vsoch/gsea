# Should be run after findSigOutput.py to threshold results

# Read in data file (output from findSigOutput.py)
data = read.csv('/home/vanessa/Documents/Work/GENE_EXPRESSION/nextbio/BPD/gsea_report/unthresh_matrix_all4559.txt',sep="\t",head=TRUE)
data = data[2:dim(data)[1],]

data$FDR.q.val = as.numeric(as.character(data$FDR.q.val))

# Which are significant?
sig = data[which(data$FDR.q.val <= .05),]

# Write results to file
write.table(sig, file = '/home/vanessa/Documents/Work/GENE_EXPRESSION/nextbio/ALZ/gsea_report/sig_matrix_all.txt')

reports = c()
for (s in 1:dim(sig)[1]){
  result = sig[s,]
  folder =   paste("/home/vanessa/Documents/Work/GENE_EXPRESSION/nextbio/ASD/gsea_output/",as.character(result$FOLDER),"/",sep="")
  # Get files that match the result pattern
  g = strsplit(as.character(result$GS),"_"); g=paste(g[[1]][1:2],sep="_")
  
  # Get all the files
  files =  list.files(folder)
  files.1 = list.files(folder)
  files = strsplit(files,"\\.")
  
  # Get rid of starter string
  #files = strsplit(sub("\\.", "*", files), "\\*") 
  
  for (f in 1:length(files)){
    if (files[[f]][2] == result$GS) {
      report = cbind(result,paste(folder,files.1[f],sep=""))
      reports = rbind(reports,report)
    }    
  }
    
}

colnames(reports)[15] = "filepath"
save(reports,file="/home/vanessa/Documents/Work/GENE_EXPRESSION/nextbio/gsea_report/sig_data_all.Rda")

###################
load("/home/vanessa/Documents/Work/GENE_EXPRESSION/nextbio/gsea_report/sig_data_all.Rda")
idx = seq(from=1,to=dim(reports)[1],by=2)
library(gridExtra)

# Now let's write a nice PDF report!
# This is the list of PDF files that we will merge, in this order
order = c()
count=1
for (r in idx) {
  pdfy=paste("/home/vanessa/Documents/Work/GENE_EXPRESSION/nextbio/gsea_report/table",count,".pdf",sep="")
  pdf(file=pdfy, height=2, width=20)  
  grid.table(reports[r,1:14])
  dev.off()
  count = count+1
  order = c(order,pdfy,as.character(reports[r,15]))
}

# Here is the commpand to run in the terminal to combine files into one PDF
cat("pdftk",order,"cat output /home/vanessa/Documents/Work/GENE_EXPRESSION/nextbio/gsea_report/gsea_report_diabolical.pdf")
