# This will calculate tanimoto scores for the C1 gene set in GSEA
infile = "/home/vanessa/Documents/Dropbox/Code/R/GSEA-P-R/GeneSetDatabases/C1.gmt"
raw <- scan(infile, what="", sep="\n")

geneID = c()
genes = list()
desc = c()
for (r in 1:length(raw)){
  dat = strsplit(raw[r],"\t") 
  geneID = c(geneID,dat[[1]][1])
  genes = c(genes,list(dat[[1]][-c(1,2)]))
  desc = c(desc,dat[[1]][2])
}

TS = array(dim=c(length(genes),length(genes)))

# Now calculate tanimoto scores
for (u in 1:length(genes)) {
  cat ("Calculating Tanimotos for Term",u,"of",length(genes),"\n")
  # Now we assess overlap of sets - calculate tanimoto score, intersection / union
  g1 = as.character(genes[[u]])
  tanimoto = c()  # This will be vector of scores for one value
  for (o in 1:length(genes)) {
    g2 = as.character(genes[[o]])
    inter = length(which(g1 %in% g2))
    uni = length(unique(c(g1,g2)))
    score = inter/uni
    tanimoto = c(tanimoto,score)
  }
  TS[u,] = tanimoto
}

rownames(TS) = geneID
colnames(TS) = geneID
write.table(TS,file="/scratch/users/vsochat/DATA/GENE_EXPRESSION/neurosynth/tanimoto/C1tani.csv")
save(TS,file="/scratch/users/vsochat/DATA/GENE_EXPRESSION/neurosynth/tanimoto/C1tani.Rda")
