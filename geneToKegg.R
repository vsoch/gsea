# This script will convert a list of gene symbols to kegg (KO) identifiers (TOP) 
# and Kegg to Swiss prot (BOTTOM)

# GENE SYMBOL TO KEGG --------------------------------------------------------------------
# We need libraries to read and parse URL
library(RCurl)
options(RCurlOptions = list(cainfo = system.file("CurlSSL", "cacert.pem", package = "RCurl")))

# Read in file with Gene Symbols
# This is tab separated file with Id, Label, color from Geshi
infile = '/home/vanessa/Desktop/group5[Nodes].csv'
data = read.table(infile,sep="\t",head=TRUE)


# For each gene, look up the hsa identifiers
hsa = c()
for (g in 1:length(data$Id)) {
  query = paste("http://rest.kegg.jp/find/hsa/",data$Id[g],sep="")
  result = getURL(query)
  hsa = c(hsa,strsplit(result,"\t")[[1]][1])  
}

# Squash them together
data = cbind(data,hsa)

# Filter out genes missing hsa ID
missing = data[-grep("hsa:",hsa),]
missing$hsa = rep("<NA>",dim(missing)[2])
missing$ko = rep("<NA>",dim(missing)[2])
data = data[grep("hsa:",hsa),]

# Now look up KO ID
ko = c()
for (g in 1:length(data$Id)) {
  query = paste("http://rest.kegg.jp/link/ko/",as.character(data$hsa[g]),sep="")
  result = getURL(query)
  ko = c(ko,as.character(gsub("\n","",strsplit(result,"\t")[[1]][2])))  
}

# More squashing
ko = as.character(ko)
data = cbind(data,ko)

# Add missing back
data = rbind(data,missing)

# Write to file
outfile = gsub(".csv","_ko.csv",infile)
write.table(data,file=outfile,row.names=FALSE,col.names=TRUE,sep="\t")


# KEGG TO SWISS PROT -(one to many)--------------------------------------------------------
infile = '/home/vanessa/Downloads/FOAM-onto_rel1(1).tsv'
data = read.table(infile,sep="\t",head=TRUE)

# For each gene, look up swiss prot
swiss = list()
for (g in 5679:length(data$KO)) {
  cat("processing",g,"of",length(data$KO),"\n")
  tmp = c()
  query = paste('http://rest.genome.jp/link/uniprot/ko:',data$KO[g],sep="")
  result = getURL(query)
  entries = strsplit(result,"\n")
  for (e in entries[[1]]){
    tmp = c(tmp,strsplit(e,"\t")[[1]][2])   
  }
  swiss = c(swiss,list(tmp))  
}

# Save into Rda object - a list with 1) the KO id, and 2) a list of identifiers
# EG, res$ko[1] corresponds to IDs as res$uniprot[[1]]
res = list(ko = data$KO, uniprot = swiss)


# SWISS PROT TO KEGG -(one to one)--------------------------------------------------------
infile = 'uniprot-seq-ids.dat'
data = readLines(infile)
  

# For each gene, look up swiss prot
ko = list()
alls = c()
gene = c()
for (g in 702:length(data)) {
  cat("processing",g,"of",length(data),"\n")
  tmp = strsplit(data[g],'"'); tmp = tmp[[1]]
  tmp = gsub(" ","",tmp,fixed=TRUE)
  gene = c(gene,substring(tmp[1],3,nchar(tmp[1])-1))
  tmp = tmp[2:length(tmp)]
  idx = which(tmp == "")
  tmp = tmp[-idx]
  kos = c()
  for (t in tmp){
    query = paste('http://rest.genome.jp/link/ko/uniprot:',t,sep="")
    result = getURL(query)
    entries = strsplit(result,"\n")
    for (e in entries[[1]]){
      kos = c(kos,strsplit(e,"\t")[[1]][2])   
    }
  }
  ko = c(ko,unique(kos))
  alls = c(alls,unique(kos))
}
  
result = list(id=gene,singleList=alls,ko=ko)
save(result,file="/scratch/users/vsochat/DATA/GUTMAP/koFromSwiss.Rda")