library(topGO)
library(RColorBrewer)
library(tidyverse)

ab <- "sa"
setwd(paste0("/home/ajperez/Nextcloud/ncbidatasets/", ab))

# Files
FOLDER <- "enrichment"
file_genes <- "freq"
#file_genes <- "percent_nocrispr_crispr_02_07.id" #spacers|repeats|flanks
#file_genes <- "spacers_sinthie.id" #spacers|repeats|flanks
### grep -E "ab00001_04005|ab00004_00550" pangenome_annot.tsv OXA-23
file_background <- paste0("roary/", ab, "2/pangenome_references_", ab , "_uniprot_bacteria_go.tsv")
Nodes <- 10 # number of processes to show

if (!dir.exists(FOLDER)){
  dir.create(FOLDER)
}

# Iterate through the two ontologies
#for (ONT in c("BP", "CC", "MF")) {
for (ONT in c("BP")) {
  letter <- substr(ONT, 2, 2)
  Ontology <- paste0("GO.", letter, ".ID") #PFC (BP MF CC)
    
  #Create temp file
  data <- read.csv(file_background, sep = "\t", header = TRUE, row.names = NULL)[,(c('X.ID', Ontology))]
  data[[Ontology]] <- as.character(gsub(';', ', ', data[[Ontology]]))
  file_temp <- paste0(file_background,"2")
  write.table(data, file = file_temp, sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
  genes <- read.csv(file_genes, header=F)$V1

  # Get background annotation
  GOesByID <- readMappings(file = file_temp)
  bg_genes <- names(GOesByID)

  compared_genes <- factor(as.integer(bg_genes %in% genes))
  names(compared_genes) <- bg_genes
  
  # Create topGO object
  GOdata <- new("topGOdata", ontology = ONT, allGenes = compared_genes,
                annot = annFUN.gene2GO, gene2GO = GOesByID)
  asd <- unlist(Term(GOTERM))
  
  # Run Fisher test
  resultFisher <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
  
  # Create and print table with enrichment result
  allRes <- GenTable(GOdata, classicFisher = resultFisher, topNodes = Nodes)
  
  # Different palettes
  palette <- c("#F52A2A", "#D561EA", "#61B0EA", "green", "#E89B57", "#E4EA61", "white") # alternative palette
  myPalette <- colorRampPalette(rev(brewer.pal(9, "YlOrRd")))
  
  # Figure
  ########
  layout(t(1:2), widths=c(8,3))
  par(mar=c(4, .5, .7, .7), oma=c(3, 15, 3, 4), las=1)
  
  pvalue <- as.numeric(gsub("<", "", allRes$classicFisher)) # remove '<' symbols
  allRes$classicFisher <- pvalue
  max_value <- as.integer(max(-log(pvalue)))+1
  pv_range <- exp(-seq(max_value, 0, -1))
  allRes <- mutate(allRes, plot_id = paste(GO.ID, Term, sep = " - "))
  
  mylabels <- paste (allRes$GO.ID, "-",  asd[allRes$GO.ID])
  mybreaks <- 10^-(0:30)
  
  p <- ggplot(data=allRes, aes(x=reorder(plot_id, Significant), y = Significant)) +
    geom_bar(stat="identity", color="black", aes(fill=as.numeric(log(classicFisher))), size = 0.3)+
    geom_text(aes(label=mylabels), position=position_fill(vjust=0), hjust=0, fontface="bold", size = 5) +
    coord_flip() + 
    theme(panel.background = element_blank(), panel.grid.major.x = element_line(colour = "darkgrey", size=0.75),
          panel.grid.minor.x = element_line(colour = "grey", size=0.75), axis.title.y=element_blank(), 
          axis.text.y=element_blank(), axis.ticks.y=element_blank(), 
          axis.ticks.x =element_blank(), axis.line.y=element_blank(), axis.text=element_text(size=12)) +
    ylab("Number of genes") +
    guides(fill = guide_colourbar(barheight = 25, reverse=T)) +
    scale_fill_gradientn(name = "p-value", colours = myPalette(4), breaks = log(mybreaks), 
                         guide = guide_colourbar(reverse = TRUE), labels = mybreaks) +
    scale_y_continuous(breaks = seq(0, max(allRes$Significant), by = 10))
  print(p)
  
  print(allRes)
  next
  
  # Save results
  if (exists("allRes") & nrow(allRes) > 4) {
    # Save table and figure 
    write.table(allRes, file = paste0(FOLDER, "/enrichment_", ONT, ".tsv"), quote = F, row.names = F, sep = "\t")
    pdf(paste0(FOLDER, "/enrichment_", ONT, ".pdf"), width=16, height=8, paper='special',)
    print(p)
    dev.off()
  }
  
  #List of genes by enriched GO
  GOnames <- as.vector(allRes$GO.ID)
  allGenes <- genesInTerm(GOdata, GOnames)
  significantGenes <- list() 
  for(x in 1:Nodes){
    significantGenes[[x]] <- allGenes[[x]][allGenes[[x]] %in% as.vector(genes)]
  }
  names(significantGenes) <- allRes$Term
  significantGenes
  new_folder <- (paste0(FOLDER, "/", ONT, "/"))
  if (!file.exists(new_folder)) { dir.create(new_folder) }
  for(x in 1:Nodes){
    write.table(significantGenes[x], quote = FALSE, sep = "\t", 
                file = paste(new_folder, x, "-", str_replace_all(names(significantGenes)[x], "/", "_"), ".tsv", sep=""),
                col.names = F, row.names = F)
  }
  print(allRes)
}

