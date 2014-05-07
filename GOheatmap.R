# This script will read in all GO annotations, and create
# a heatmap to show distribution of each

library(org.Hs.eg.db)
library(GO.db)
library("GOstats")

# Read in list of terms
uids = read.csv('/home/vanessa/Documents/Dropbox/Code/Python/neurosynth/data/featureUID.txt',sep=',',head=FALSE)
terms = uids[,2]

input = '/home/vanessa/Documents/Work/GENE_EXPRESSION/GO/list'

file=paste(outdir,"/",term,"_list.trm",sep=""))

# Create dict with counts
dict = list()

for (f in 1:length(terms)) {
  term = terms[f]
  file = read.csv(file=paste(input,"/",term,"_list.trm",sep=""),head=FALSE)
  file = file$V1
  dict = c(dict,paste(file))
}

# Get unique IDs
unis = unique(dict)
dict = c()
idx = seq(1,length(unis))

for (u in 1:length(unis)){
  dict[unis[[u]]] = 0
}

HM = c()

for (f in 1:length(terms)) {
  term = terms[f]
  file = read.csv(file=paste(input,"/",term,"_list.trm",sep=""),head=FALSE)
  file = file$V1
  tmp = dict
  for (l in 1:length(file)){
    tmp[as.character(file[l])] = tmp[as.character(file[l])] + 1
  }
  HM = rbind(HM,tmp)
}

HM[is.na(HM)] = 0
dim(HM)
rownames(HM) = terms
heatmap(HM)

# Distance matrix
disty = dist(as.matrix(HM))
hc = hclust(disty)
plot(hc,main="Behavioral Terms Grouped by GO Annotation Similarity")
groups = cutree(hc,)

for (g in 1:25){
  cat(as.character(terms[which(groups==g)]))
  cat("\n")
}
