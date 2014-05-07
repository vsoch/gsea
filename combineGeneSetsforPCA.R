# Create compiled data frame for specific subset of terms
# For performing PCA / visualization

# Set working directory
setwd("/home/vanessa/Documents/Work/GENE_EXPRESSION/neurosynth/results/sigresults")

# Define term sets
termsets = c("TRM397_angry_probeSet_MEdata_up.Rda","TRM261_noun_probeSet_MEdata_up.Rda","TRM312_skills_probeSet_MEdata_up.Rda","TRM515_covert_probeSet_MEdata_up.Rda")
datapath = "/home/vanessa/Documents/Work/GENE_EXPRESSION/neurosynth/probeSets/9mmsq/data/"

# Read in the term sets
terms = list()
suids = list()
for (t in 1:length(termsets)){
  load(paste(datapath,termsets[t],sep=""))
  terms = c(terms,list(mat))
  suids = c(suids,list(rownames(mat)))
}

# Look at overlaps:
overlaps = array(dim=c(4,4))
for (i in 1:4) {
  for (j in 1:4) {
    overlaps[i,j] = length(which(suids[[i]] %in% suids[[j]]))
  }
}

rownames(overlaps) = c("angry","noun","skills","covert")
colnames(overlaps) = c("angry","noun","skills","covert")

# Number of samples for each one?
samps = c()
for (t in 1:length(termsets)){
   samps = c(samps,length(suids[[t]]))
}

# Combine sample data - make unique labels for each
for (i in 1:4) {
  tmp = suids[[i]]
  term = colnames(overlaps)[i]
  for (t in 1:length(tmp)){
    tmp[t] = paste(tmp[t],"_",term,sep="")
  }
  rownames(terms[[i]]) = tmp
}

# Get all unique gene names
allgenes = c()
for (i in 1:4) {
  allgenes = c(allgenes,as.character(colnames(terms[[i]])))
}
allgenes = unique(allgenes)

# Get the total number of samples
nsamp = 0
for (i in c(1:4)) {
  nsamp = nsamp + nrow(terms[[i]])
}

# Create data matrix for result
matrix = array(data=0,dim=c(1,length(allgenes)))
colnames(matrix) = allgenes

# Now fill in values
for (c in 1:length(suids)){
  idx = which(colnames(matrix) %in% colnames(terms[[c]]))
  tmp = array(data=0,dim=c(nrow(terms[[c]]),length(allgenes)))
  tmp[,idx] = terms[[c]]
  rownames(tmp) = rownames(terms[[c]])  
  matrix = rbind(matrix,tmp)
}

# Get rid of first row
matrix = matrix[2:dim(matrix)[1],]

# Save to file
save(matrix,file="sigsamples_n4.Rda")
write.table(matrix,file="sigSamplesTable.csv",sep=",")

# Now try PCA! (BrainTermPCA.R)