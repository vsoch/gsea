# This script will produce RamiGO images for lists of GO terms
run <- function(term) {
  
  
  # CAN"T RUN THIS - RAMIGO IS BROKEN
  
library('RamiGO')

# Read in terms and uids
uids = read.csv('C:/Users/Vanessa/Documents/Dropbox/Code/Python/neurosynth/data/featureUID.txt',sep=',',head=FALSE)
terms = uids[,2]
inputdir = "C:/Users/Vanessa/Documents/Work/DATA/GENE_EXPRESSION/GO/shapleyFilter"

for (t in 1:length(terms)) {
  term = terms[t]
  cat("Processing",t,"of",length(terms),"\n")

  # DOWN
  file = paste(inputdir,"/",term,"_down_list.trm",sep="")
  if file.exists(file){
    GO = read.csv(file,sep="\n",head=FALSE); GO = as.character(GO$V1)
    color = rainbow(length(GO))[rank(GO)]
    pngRes <- getAmigoTree(goIDs=GO, color=color, filename=paste("C:/Users/Vanessa/Documents/Work/DATA/GENE_EXPRESSION/GO/img/updown",term,"_down_amiGO"), picType="png", saveResult=TRUE)
  }  
  # UP
  file = paste(inputdir,"/",term,"_up_list.trm",sep="")
  if file.exists(file){
    GO = read.csv(file,sep="\n",head=FALSE); GO = as.character(GO$V1)
    color = rainbow(length(GO))[rank(GO)]
    pngRes <- getAmigoTree(goIDs=GO, color=color, filename=paste("C:/Users/Vanessa/Documents/Work/DATA/GENE_EXPRESSION/GO/img/updown/",term,"_down_amiGO"), picType="png", saveResult=TRUE)
    }
}

# This will print HTML for images
for (i in 1:length(terms)){
  term = terms[i]
  cat('<option value="',as.character(term),'" data-picture="img/',as.character(term),' _amiGO.png">',as.character(term),'</option>\n',sep="")    
}
  

