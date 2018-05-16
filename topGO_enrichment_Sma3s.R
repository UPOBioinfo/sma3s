# User folder
setwd("/home/ajperez/R/topGO")

# Mandatory library
library(topGO)

# Input files and (gene IDs and background annotation with Sma3s format)
file_genes <- "pratense.id"
file_background <- "transcriptoma_phleum_pratense_UniSprPla_a123_r20_t0.1.tsv"
Nodes <- 40 # number of processes to show

file_temp <- paste0(file_background,"2")
linux <- paste("grep -v '#'", file_background, "| cut -f1,5 | sed 's/;/, /g' >", file_temp)
linux_rm <- paste("rm",file_temp)
system(linux) # create temp file

# Get gene IDs for the enrichment
genes <- read.csv(file_genes, header=F)$V1

# Get background annotation
GOesByID <- readMappings(file = file_temp)
bg_genes <- names(GOesByID)

# Compare genes vs bg_genes
compared_genes <- factor(as.integer(bg_genes %in% genes))
names(compared_genes) <- bg_genes

# Create topGO object
GOdata <- new("topGOdata", ontology = "BP", allGenes = compared_genes,
              annot = annFUN.gene2GO, gene2GO = GOesByID)

# Run Fisher test
resultFisher <- runTest(GOdata, algorithm = "classic", statistic = "fisher")

# Create and print table with enrichment result
allRes <- GenTable(GOdata, classicFisher = resultFisher, topNodes = Nodes)
allRes

# Graphic
#########
layout(t(1:2), widths=c(8,1))
par(mar=c(4, .5, .7, .7), oma=c(3, 15, 3, 4), las=1)

#Generate a continuous color palette
rbPal <- colorRampPalette(c('red','white','blue'))
pvalue <- as.numeric(gsub("<", "", allRes$classicFisher)) # remove '<' symbols
max_value <- as.integer(max(-log(pvalue)))+1
pv_range <- exp(-seq(max_value, 0, -1))
allRes$Color <- rbPal(max_value)[cut(pvalue, pv_range)]

# Print figure
o <- order(allRes$Significant, decreasing = TRUE)
barplot(allRes$Significant[o],  names.arg=allRes$Term[o], las=2, horiz= TRUE, col=allRes$Color[o],
        xlab="Number of sequences", main=file_genes, sub=file_background, cex.names=0.85)

# Legend
image(0, seq(1, max_value), t(seq_along(seq(1, max_value))), col=rev(rbPal(max_value)), axes=FALSE, ann=FALSE)
pv_label <- exp(-seq(log(1), -log(min(pvalue)), l=6))
pv_label <- formatC(pv_label, format = "e", digits = 2)
axis(4, at=seq(1, max_value, length=6), labels=c(1, pv_label[2:6]), cex.axis=0.85)
title("p-value", cex.main = 0.8)

# Remove temp file
system(linux_rm)

