# This script will read in gene sets, and create a "database" of gene sets to be used with GSEA.
# Read in file with gene_id lookup - we need to overlap genes in allen brain atlas with this data

# Read in list of terms
uids = read.csv('/home/vanessa/Documents/Dropbox/Code/Python/neurosynth/data/featureUID.txt',sep=',',head=FALSE)
terms = uids[,2]

# Write to this file
outfile = paste("/home/vanessa/Documents/Dropbox/Code/R/GSEA-P-R/GeneSetDatabases/brainTerms.gmt",sep="")

# Read in file with probes
probes = read.csv("/home/vanessa/Documents/Work/ALLEN/Probes.csv",sep=",",header=FALSE)

# UP PROBES
sink(outfile)
for (t in 1:length(terms)) {
  term = terms[t]
  tid = uids[t,1]
    
  if (file.exists(paste("/home/vanessa/Documents/Work/GENE_EXPRESSION/neurosynth/probeSets/9mmsq/list/",tid,"_",term,"_probeSet_up.Rda",sep=""))) {
    
    ptmp = probes
    
    # Load the up set
    load(file=paste("/home/vanessa/Documents/Work/GENE_EXPRESSION/neurosynth/probeSets/9mmsq/list/",tid,"_",term,"_probeSet_up.Rda",sep=""))
    # Sort by pvalue significance
    resultup$pval = as.numeric(as.character(resultup$pval))
    resultup = resultup[with(resultup, order(pval)),]
    pid = as.numeric(gsub("pid_","",resultup$probes))
    # Find probes in probes table
    ptmp = ptmp[which(ptmp$V1 %in% pid),]
    # Get the gene list
    genes = unique(ptmp$V4)
    label = paste(as.character(term),"_up",sep="")
    cat(label,"Allen Brain Atlas expression up regulated subset for neurosynth fdr .05 corrected brain map sorted by pval more sig first",as.character(genes),"\n",sep="\t") 
    rm(resultup)
  }
  if (file.exists(paste("/home/vanessa/Documents/Work/GENE_EXPRESSION/neurosynth/probeSets/9mmsq/list/",tid,"_",term,"_probeSet_down.Rda",sep=""))) {

    # Load the down set
    load(file=paste("/home/vanessa/Documents/Work/GENE_EXPRESSION/neurosynth/probeSets/9mmsq/list/",tid,"_",term,"_probeSet_down.Rda",sep=""))
    ptmp = probes
    # Sort by pvalue significance
    resultdown$pval = as.numeric(as.character(resultdown$pval))
    resultdown = resultdown[with(resultdown, order(pval)),]
    pid = as.numeric(gsub("pid_","",resultdown$probes))
    # Find probes in probes table
    ptmp = ptmp[which(ptmp$V1 %in% pid),]
    # Get the gene list
    genes = unique(ptmp$V4)    
    label = paste(as.character(term),"_down",sep="")
    cat(label,"Allen Brain Atlas expression down regulated subset for neurosynth fdr .05 corrected brain map sorted by pval more sig first",as.character(genes),"\n",sep="\t") 
    rm(resultdown)
  }
}
sink()