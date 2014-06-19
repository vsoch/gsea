library(RCurl)

# import the gene symbol list - should be a csv file with the following columns
# DISORDER   TERM         RESULT    GENE      METRIC    SCORE       RUNNING.ES CORE
# AD         AUTOMATIC_UP GSE12685  UBOX5     95        0.16748269  0.1515688  Yes
listID_genes_info = read.csv("/home/vanessa/Documents/Work/GENE_EXPRESSION/neurosynth/results/sigresults/AllGenesNoFilter.csv", header=T) 

# For each term
terms = as.character(unique(listID_genes_info$TERM))
for (t in terms){

  # Get the unique genes
  genes = listID_genes_info[which(listID_genes_info$TERM %in% t),]
  genes = as.character(sort(unique(genes$GENE)))

  # For each gene, look up the hsa identifiers
  hsa = c()
  for (g in 1:length(genes)) {
    query = paste("http://rest.kegg.jp/find/hsa/",genes[g],sep="")
    result = getURL(query)
    id = strsplit(result,"\t")[[1]][1]
    hsa = c(hsa,id)  
    cat(g,"-",id,"\n")
  }
  
  # These are the genes we are missing
  missing = genes[which(hsa=="\n")]
  cat("Missing:",missing,"\n")
  
  # OR just get rid of
  if (length(missing) > 0){
    genes = genes[-which(hsa=="\n")]
    hsa = hsa[-which(hsa=="\n")]
  }
  # Now we need to look up brite ids
  
  bridsok = list() # A list with br: bride IDs for each of the hsa
  emptybid = c() # List of empty indices, to remove from hsa and genes after
  for (i in 1:length(hsa)){
    brids = getURL(paste("http://rest.kegg.jp/link/br/",hsa[i], sep=""))
    # Get rid of newlines
    brids = strsplit(brids,"\t|\n")[[1]]
    # Get the ones that begin with br:
    brids = brids[grep("br:hsa",brids)]
    if (length(brids)>0){
      bridsok = c(bridsok,list(brids))
      cat(hsa[i],"-",brids,"\n")
    } else {
      cat("Missing,",hsa[i],i,"\n")
      emptybid = c(emptybid,i)
    }
  }
  
  # Get rid of ones with missing
  if (length(emptybid)>0){
    genes = genes[-emptybid]
    emptyhsa = hsa[emptybid]
    hsa = hsa[-emptybid]
  }
  
  # Now obtain entire pathway
  Pathway = list()
  genelist = c()
  hsalist = c()
  britelist = c()
  
  for (i in 1:length(bridsok)) { 
    cat("Processing",i,"of",length(bridsok),"\n")
    # This will hold our result, with A,B,C,D, etc.
    entry = c()
    
    # We need to do query for each bridsok
    for (j in 1:length(bridsok[[i]])) {
      # Download the whole BRITE associated with the KEGG ID
      tmpl = getURL(paste("http://rest.kegg.jp/get/",bridsok[[i]][j], sep=""))
      tmpl = strsplit(tmpl,"\n")[[1]]
      # Find line with our id
      id = gsub("hsa:","",hsa[[i]])  
    
      # Find the gene in the hierarchy, could be level B,C,D,E...
      match = grep(paste("[[:blank:]]",id,"[[:blank:]]",sep="+"),tmpl,val=TRUE)
      matchindex = grep(paste("[[:blank:]]",id,"[[:blank:]]",sep="+"), tmpl)
      
      # For each match:
      for (m in 1:length(match)) {
        tmp = tmpl
        mat = match[m]
        midx = matchindex[m]
        
        letter = substr(mat,start=1,stop=1)
        entry[letter] = mat
    
        # Get rid of all the hierarchy including and after our match
        tmp = tmp[-seq(from=midx,to=length(tmp))]
    
        # Now reverse the list
        tmp = rev(tmp)
    
        # While we still have letters left
        while(letter != "A") {
          # Get rid of the letter we already found
          nix = grep(paste(letter,"[[:blank:]]",sep="+"), tmp)
          tmp = tmp[-nix]
    
          # The top is our next letter
          midx = 1
          mat = tmp[midx]  
          letter = substr(mat,start=1,stop=1)
          entry[letter] = mat
        }
    
        # Fill in Zero for the letters we don't have
        letters = c("A","B","C","D","E")
        entry[letters[-which(letters %in% names(entry))]] = 0
    
        # Now add object to pathway list, geneid, and pathway id
        Pathway = c(Pathway,list(entry))
        genelist = c(genelist,genes[i])
        hsalist = c(hsalist,hsa[i])
        britelist = c(britelist,bridsok[[i]][j])
      } # End loop through match list
      
    } # End loop through brite IDS
    
  } # End loop through all hsa

  # Save to term output file 
  PATHWAY = list(ABCDE = Pathway,genes=genelist,hsaid=hsalist,brite=britelist,term=t,missinghsa=missing,missingbid=emptyhsa)
  save(PATHWAY,file=paste(t,"_pathway.Rda",sep=""))
} # End loop through term