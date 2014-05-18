# Format GSE (array download from GEO) for gsea chip

# Read in file (look at visually to determine how many to skip)
filey = 'GPL3921-25447.txt'
data = read.csv(filey,head=TRUE,sep="\t",skip=16)
tmp =data.frame(as.character(data$ID),as.character(data$Gene.Symbol),as.character(data$Gene.Title))
tmp = tmp[-which(tmp[,2]==""),]
colnames(tmp) = c("Probe Set ID","Gene Symbol","Gene Title")
write.table(tmp,file="/home/vanessa/Documents/Work/GENE_EXPRESSION/cmap/probes/HTHGU133A.chip",col.names=TRUE,row.names=FALSE,quote=FALSE,sep="\t")