# Functional enrichments from: id_file and Sma3s annotation
# change lines in block INPUT (###HERE)
# UPOBioinfo, 2019
library(topGO)
library(ggplot2)
library(RColorBrewer)

# INPUT (###HERE)
setwd("/home/alumno/Descargas")                    # working directory
file_genes <- "core_proteins.id"                   # file with identifiers
file_background <- "pangenome_uniprot_BsEc_go.tsv" # Annotation file from Sma3s
Nodes <- 15                                        # number of annotations to show

# Result folder
FOLDER <- "enrichment"
if (!dir.exists(FOLDER)){
  dir.create(FOLDER)
}

# Iterate through the two ontologies
for (ONT in c("BP", "CC", "MF")) {
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
  #allRes <- mutate(allRes, plot_id = paste(GO.ID, Term, sep = " - "))
  
  mylabels <- paste (allRes$GO.ID, "-",  asd[allRes$GO.ID])
  mybreaks <- 10^-(0:30)
  
  p <- ggplot(data=allRes, aes(x=reorder(Term, Significant), y = Significant)) +
    geom_bar(stat="identity", color="black", aes(fill=as.numeric(log(classicFisher))), size = 0.3)+
    geom_text(aes(label=mylabels), position=position_fill(vjust=0), hjust=0, fontface="bold", size = 5) +
    coord_flip() + 
    theme(panel.background = element_blank(), panel.grid.major.x = element_line(colour = "darkgrey", size=0.75),
          panel.grid.minor.x = element_line(colour = "grey", size=0.75), axis.title.y=element_blank(), 
          axis.text.y=element_blank(), axis.ticks.y=element_blank(), axis.text.x = element_text(angle = 45, hjust = 0),
          axis.ticks.x =element_blank(), axis.line.y=element_blank(), axis.text=element_text(size=12)) +
    ylab("Number of genes") +
    guides(fill = guide_colourbar(barheight = 25, reverse=T)) +
    scale_fill_gradientn(name = "p-value", colours = myPalette(4), breaks = log(mybreaks), 
                         guide = guide_colourbar(reverse = TRUE), labels = mybreaks) +
    scale_y_continuous(breaks = seq(0, max(allRes$Significant), by = 20))
  print(p)
  
  print(allRes)
  
  # Save results
  if (exists("allRes") & nrow(allRes) > 4) {
    # Save table and figure 
    write.table(allRes, file = paste0(FOLDER, "/enrichment_", ONT, ".tsv"), quote = F, row.names = F, sep = "\t")
    pdf(paste0(FOLDER, "/enrichment_", ONT, ".pdf"), width=16, height=8, paper='special')
    png(paste0(FOLDER, "/enrichment_", ONT, ".png"), width = 900, height = 600)
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
}

