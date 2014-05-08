# This script will look at network of gene functions for each of term subsets

source("http://bioconductor.org/biocLite.R")
setwd('/home/vanessa/Documents/Work/GENE_EXPRESSION/neurosynth/probeSets/9mmsq/list')
#A set of annotation maps describing the entire Gene Ontology
biocLite("GO.db")
n
biocLite("topGO")
n
biocLite("GOstats")
n
biocLite("org.Hs.eg.db")
n
library(org.Hs.eg.db)
library(GO.db)
library("GOstats")

# Read in list of terms
uids = read.csv('/home/vanessa/Documents/Dropbox/Code/Python/neurosynth/data/featureUID.txt',sep=',',head=FALSE)
terms = uids[,2]

x <- org.Hs.egGO
entrez_object <- org.Hs.egGO
mapped_genes <- mappedkeys(entrez_object)
entrez_to_go <- as.list(x[mapped_genes])
go_object <- as.list(org.Hs.egGO2EG)

# Read in file with entrez-affy lookup
lookup = read.csv(file='/home/vanessa/Documents/Work/GENE_EXPRESSION/probes_entrez_affy_lookup.csv',head=TRUE,sep=",")

# UPDATED FOR SHAPLEY UP AND DOWN SETS 4/1/2014
# Here is if we want to save annotatoin objects
for (t in 293:length(terms)) {
  cat("Processing",t,"of",length(terms),"\n")
  term = terms[t]
  tid = uids[t,1]
  cat("Term",t,"of",length(terms),"\n")

  # UP REGULATED LIST
  if (file.exists(paste("/home/vanessa/Documents/Work/GENE_EXPRESSION/neurosynth/probeSets/9mmsq/list/",tid,"_",term,"_probeSet_up.Rda",sep=""))) {
    load(file=paste("/home/vanessa/Documents/Work/GENE_EXPRESSION/neurosynth/probeSets/9mmsq/list/",tid,"_",term,"_probeSet_up.Rda",sep=""))
  
    # Get entrez ids
    idx = which(lookup$probe_id %in% gsub("pid_","",resultup$probes))
    probes = lookup$entrez_id[idx]
  
    # Find the indices for entrez gene ids that we have in our subset
    idx = which(mapped_genes %in% probes)
    go_subset = entrez_to_go[idx]
    idx = which(probes %in% mapped_genes)
    genes = as.character(probes[idx])
    genes = unique(genes)
    universe = mapped_genes
  
    # Now perform enrichment analysis
    params <- new('GOHyperGParams',
                geneIds=genes,
                universeGeneIds=universe,
                ontology='BP',
                pvalueCutoff=0.001,
                conditional=F,
                testDirection='over',
                annotation="org.Hs.eg.db"
    )
    hgOver <- hyperGTest(params)
    # Print list to file
    tabley = summary(hgOver)
    if (dim(tabley)[1] != 0) {
      cat(tabley$GOBPID,sep="\n",file=paste("/home/vanessa/Documents/Work/GENE_EXPRESSION/GO/list/shapleyFilter/",term,"_up_list.trm",sep=""))  
    }
    rm(resultup)
  }
  
  # DOWN REGULATED LIST
  if (file.exists(paste("/home/vanessa/Documents/Work/GENE_EXPRESSION/neurosynth/probeSets/9mmsq/list/",tid,"_",term,"_probeSet_down.Rda",sep=""))) {
    load(file=paste("/home/vanessa/Documents/Work/GENE_EXPRESSION/neurosynth/probeSets/9mmsq/list/",tid,"_",term,"_probeSet_down.Rda",sep=""))
    
    # Get entrez ids
    idx = which(lookup$probe_id %in% gsub("pid_","",resultdown$probes))
    probes = lookup$entrez_id[idx]
    
    # Find the indices for entrez gene ids that we have in our subset
    idx = which(mapped_genes %in% probes)
    go_subset = entrez_to_go[idx]
    idx = which(probes %in% mapped_genes)
    genes = as.character(probes[idx])
    genes = unique(genes)
    universe = mapped_genes
    
    # Now perform enrichment analysis
    params <- new('GOHyperGParams',
                  geneIds=genes,
                  universeGeneIds=universe,
                  ontology='BP',
                  pvalueCutoff=0.001,
                  conditional=F,
                  testDirection='over',
                  annotation="org.Hs.eg.db"
    )
    hgOver <- hyperGTest(params)
    # Print list to file
    tabley = summary(hgOver)
    if (dim(tabley)[1] != 0) {
      save(tabley,file=paste("/home/vanessa/Documents/Work/GENE_EXPRESSION/GO/list/shapleyFilter/",term,"_down_list.Rda",sep=""))  
    }
    rm(resultdown)
  }
}