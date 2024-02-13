#Load the libraries
library(phyloseq)
library(ggplot2)
library(readxl)
library(DESeq2)
library(patchwork)
library(apeglm)
library(ggfortify)
library(cluster)
library(tidyverse)
library(rstatix)
library(ggpubr)
library(devtools)
library(ggbiplot)
library(RColorBrewer)
library(EnhancedVolcano)
library(pals)
library(pheatmap)
library(grid)
library(gridExtra)
library(PCAtools)
library(vsn)
library(hexbin)
library(ggplot2)
library(reshape2)
library(Hmisc)
library(stats)
library(corrplot)
library(ggplot2)
library(dplyr)
library(pheatmap)
library(RColorBrewer)
library(janitor)


#Read the counts table for the grouped dataset
dataset <- read.delim("../Input_data/new_singlem_table_v1.txt", header = T,check.names = F)
wd <- dataset[,-c(1:10, 47:52)]
rownames(wd) <- dataset$ConsensusLineage
dim(wd)

data <- dplyr::filter(wd, !grepl('d__Eukaryota', rownames(wd)))
dim(data)
#No need to remove Euks, such low numbers but tidy up the data

###Load the metadata file

metadata <- read.delim("../Input_data/metadata.txt", header = T, row.names = 1, check.names = F)
metadata <- metadata[1:54,]
dim(metadata)

###Prepare a DESeq object by creating a count matrix input with the count data
#For this, metadata inputs need to be factorised before running the matrix
colData <- metadata

# It is absolutely critical that the columns of the count matrix and the rows 
# of the column data (information about samples) are in the same order. 
# DESeq2 will not make guesses as to which column of the count matrix belongs to 
# which row of the column data, these must be provided to DESeq2 already in consistent order.
names(data) = gsub(pattern = "_S.*", replacement = "", x = names(data))
names(data) = gsub(pattern = "_COM.*", replacement = "", x = names(data))

all(rownames(colData) == colnames(data))

###Create the count matrix with DESeqDataSetFromMatrix
count.data.set <- DESeqDataSetFromMatrix(countData= data, 
                                         colData=colData, design= ~ Site+Season)
count.data.set

###Run the differential expression analyses with DESeq2
#The results give you the base mean, the log2fold changes and pvalue/adjusted pvalues
dds <- DESeq(count.data.set)
dim(dds)
dds

####Keep those OTUs with hits >=1 since we have low number of OTUs.
keep <- rowSums(counts(dds)) >=1
dds <- dds[keep,]
dim(dds)

#Normalize with variance stabilizing transformation with blind dispersion
vsd <- varianceStabilizingTransformation(dds)
vsd
head(assay(vsd))
dim(vsd)
vsd_tax <- rownames(vsd)
vsd_norm <- cbind(vsd_tax, assay(vsd))
#write.table(vsd_norm, file = "../Results/normalized_deseq2_wd2.tsv", sep = ",", quote = F, row.names = F)
library(vsn)
library(matrixStats)
rv <- rowVars(assay(vsd))
o <- order(rv,decreasing=TRUE)
dists <- dist(t(assay(vsd)))
hc <- hclust(dists, method = "ward.D2")
plot(hc, labels=vsd$Sample_Season)

#How to get PCA scree plot?
pca_res <- prcomp(t(assay(vsd)), scale. = TRUE)
## calculate the variance for each gene
rv <- rowVars(assay(vsd))
## select the ntop genes by variance
select <- order(rv, decreasing=TRUE)[seq_len(min(500, length(rv)))]
## perform a PCA on the data in assay(x) for the selected genes
pca <- prcomp(t(assay(vsd)[select,]))
## the contribution to the total variance for each component
percentVar <- pca$sdev^2 / sum( pca$sdev^2 )
##plot the "percentVar"
scree_plot=data.frame(percentVar)
scree_plot[,2]<- c(1:54)
colnames(scree_plot)<-c("variance","component_number")
ggplot(scree_plot, mapping=aes(x=component_number, y=variance))+geom_bar(stat="identity")


###########Plot the PCA plot#############
data.PCA <- plotPCA(vsd, intgroup=c("Site", "Season"), returnData=TRUE)
percentVar <- round(100 * attr(data.PCA, "percentVar"))
#data.PCA$Month <- factor(data.PCA$Month, levels = c("February", "May", "July", "August", "October", "December", "Feb-21"))

all_PCA <- ggplot(data.PCA, aes(PC1, PC2, color=Site, shape=Season)) +
  geom_point(size=2)+ 
  ggtitle("PCA_Community") +
  theme(panel.background = element_blank(),
        axis.line = element_line()) +
  scale_color_manual(values = c("#0072B2", "#009E73", "#D55E00")) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  theme(panel.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent",
                                       colour = NA_character_),
        legend.background = element_rect(fill = "transparent"))
#ggforce::geom_mark_ellipse(aes(color = data.PCA$Season), show.legend = NA)

all_PCA


ggsave("../Figures/PCA-genus-counts.png", bg = "transparent", plot = all_PCA, dpi = 600, units = c("in"), width = 5, height = 5)

###########################Biplot#####################
p <- pca(assay(vsd), metadata = colData, removeVar = 0.1)
screeplot(p, axisLabSize = 18, titleLabSize = 22)

#--------------------------------------------------------------------
####Only necessary if having other details####

# filter the expression data to match the samples in our pdata
mat <- assay(vsd)
mat <- mat[,which(colnames(mat) %in% rownames(metadata))]

dim(mat)
# check that sample names match exactly between pdata and expression data 
all(colnames(mat) == rownames(metadata))

###Run biplot again

p <- pca(mat, metadata = colData, removeVar = 0.1)
screeplot(p, axisLabSize = 18, titleLabSize = 22)

#-----------------------------------------------------------------------------------

Site <- biplot(p, colby = 'Site',
               showLoadings = F,
               shape = 'Season',
               colLegendTitle = 'Site', 
               shapekey = c('Spring'=15, 'Summer'=17, 'Winter'=16, 'Autumn' = 18),
               #colkey = c('BR1' = "#0072B2", 'BR1' =  "#009E73",  'BR1'  = "#D55E00"),
               #colkey = c('BR1' = "#bccce3", 'BR1' =  "#b9dcca",  'BR1'  = "#e7c0c1"),
               colkey = c('BR1' = "#BFD3E6FF", 'BR1' =  "#0570B0FF",  'BR1'  = "#08519CFF"),
               # encircle config
               encircle = TRUE,
               #ellipseFill = TRUE,
               #hline = 0, 
               #vline = c(-40, 0, 40),
               legendPosition = 'right', 
               legendLabSize = 15, 
               legendIconSize = 7.0,
               lab = p$metadata$Site,
               labSize = 7.0,
               gridlines.major = F,
               gridlines.minor = F,
               drawConnectors = FALSE,
               title = "PCA variation between sites") + 
  theme(panel.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent",
                                       colour = NA_character_),
        legend.background = element_rect(fill = "transparent"))
Site

Season <- biplot(p, colby = 'Season',
                 shape = 'Season',
                 #showLoadings = T,
                 colLegendTitle = 'Season', 
                 shapekey = c('Spring'=15, 'Summer'=17, 'Winter'=16, 'Autumn' = 18),
                 colkey = c( 'Spring'="#41AE76FF", 'Summer'="#006D2CFF", 'Winter'="#99D8C9FF", 'Autumn' = "#CCECE6FF"),
                 # encircle config
                 encircle = TRUE,
                 encircleFill = TRUE,
                 #hline = 0, vline = c(-30, 0, 30),
                 legendPosition = 'right', 
                 legendLabSize = 15, 
                 legendIconSize = 7.0,
                 lab = p$metadata$Site,
                 labSize = 7.0,
                 gridlines.major = F,
                 gridlines.minor = F,
                 drawConnectors = FALSE,
                 title = "PCA variation between seasons") + 
  theme(panel.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent",
                                       colour = NA_character_),
        legend.background = element_rect(fill = "transparent"))
Season

gg <- grid.arrange(Site, Season, nrow = 1)

ggsave("../Figures/Biplot-genus-counts_v2.png", plot = gg , dpi = 600, units = c("in"), width = 16, height = 8)
