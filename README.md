# Metagenomics and machine learning approach 

## Data interpretation


```
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

```
#Read the counts table for the grouped dataset

```
dataset <- read.delim("../Input_data/new_singlem_table_v1.txt", header = T,check.names = F)
wd <- dataset[,-c(1:10, 47:52)]
rownames(wd) <- dataset$ConsensusLineage
dim(wd)

data <- dplyr::filter(wd, !grepl('d__Eukaryota', rownames(wd)))
dim(data)

```
###Load the metadata file

```

metadata <- read.delim("../Input_data/metadata.txt", header = T, row.names = 1, check.names = F)
metadata <- metadata[1:54,]
dim(metadata)

```
###Prepare a DESeq object by creating a count matrix input with the count data
#For this, metadata inputs need to be factorised before running the matrix
```
colData <- metadata

```
# It is absolutely critical that the columns of the count matrix and the rows 
# of the column data (information about samples) are in the same order. 
# DESeq2 will not make guesses as to which column of the count matrix belongs to 
# which row of the column data, these must be provided to DESeq2 already in consistent order.
names(data) = gsub(pattern = "_S.*", replacement = "", x = names(data))
names(data) = gsub(pattern = "_COM.*", replacement = "", x = names(data))
```

all(rownames(colData) == colnames(data))

```
###Create the count matrix with DESeqDataSetFromMatrix

```
count.data.set <- DESeqDataSetFromMatrix(countData= data, 
                                         colData=colData, design= ~ Site+Season)
count.data.set

```
###Run the differential expression analyses with DESeq2
#The results give you the base mean, the log2fold changes and pvalue/adjusted pvalues

```
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

################################NMDS###########
sampleDists <- dist(t(assay(vsd)))
library(vegan)
library(ggplot2)
df = metadata

library(lubridate)
res <- hms(df$`Daylight hours`)        # format to 'hours:minutes:seconds'
df$Daylight <- hour(res)*60 + minute(res) + seconds(res)/60  
df

env = df %>% select(6:10,12:15,19,20,23,27)
env_scaled = scale(env, center = TRUE, scale = TRUE)

br.mds <- metaMDS(sampleDists, distance = "bray", trace = FALSE, autotransform = FALSE)
plot(br.mds)
en = envfit(br.mds$points, env, permutations = 999, na.rm = TRUE)
plot(en)

data.scores = as.data.frame(br.mds$points)
#data.scores <- as.data.frame(scores(br.mds, "species"))
data.scores$Site = metadata$Site
data.scores$Season = factor(metadata$Season, levels = c("Autumn", "Winter", "Spring", "Summer"))  

en_coord_cont = as.data.frame(scores(en, "vectors")) * ordiArrowMul(en)
en_coord_cat = as.data.frame(scores(en, "factors")) * ordiArrowMul(en)


#########Another NMDS plot + variation

xx = ggplot(data.scores, aes(x = MDS1, y = MDS2)) + 
  geom_point(size = 4, aes(shape = Season, colour = Site))+ 
  geom_segment(aes(x = 0, y = 0, xend = MDS1, yend = MDS2), 
               data = en_coord_cont, size =0.5, alpha = 0.5, colour = "grey30") +
  geom_text(data = en_coord_cont, aes(x = MDS1+1, y = MDS2), colour = "grey30", 
            fontface = "bold", label = row.names(en_coord_cont)) + 
  theme(axis.text.y = element_text(colour = "black", size = 12, face = "bold"), 
        axis.text.x = element_text(colour = "black", face = "bold", size = 12), 
        legend.text = element_text(size = 12, face ="bold", colour ="black"), 
        legend.position = "right", axis.title.y = element_text(face = "bold", size = 14), 
        axis.title.x = element_text(face = "bold", size = 14, colour = "black"), 
        legend.title = element_text(size = 14, colour = "black", face = "bold"), 
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1.2),
        legend.key=element_blank()) + 
  scale_color_manual(values = c("#0072B2", "#009E73", "#D55E00")) +
  labs(x = "NMDS1", colour = "Site", y = "NMDS2", shape = "Season") 

xx

#########Make it pretty

xy = ggplot(data.scores, aes(x = MDS1, y = MDS2)) + 
  geom_point(size = 4, aes(shape = Season, colour = Site))+ 
  coord_fixed(ratio=2) +
  geom_segment(aes(x = 0, y = 0, xend = MDS1, yend = MDS2), 
               arrow = arrow(length = unit(0.5, "cm")), colour = "grey60",
               data = en_coord_cont, size =0.75) +
  #geom_text(data = en_coord_cont, aes(x = MDS1, y = MDS2), hjust=0.5,vjust = 0.9, colour = "black", size = 5, position=position_dodge(width=3),label = row.names(en_coord_cont))+ 
  theme(axis.text.y = element_text(colour = "black", size = 12), 
        axis.text.x = element_text(colour = "black", size = 12), 
        #axis.line = element_line(size = 1, colour = "black"),
        legend.text = element_text(size = 12, colour ="black"), 
        legend.position = "right", axis.title.y = element_text(size = 14), 
        axis.title.x = element_text(size = 14, colour = "black"), 
        legend.title = element_text(size = 14, colour = "black"), 
        panel.border = element_rect(colour = "black", fill = NA, size = 1),
        panel.background = element_rect(fill = "transparent",
                                        colour = NA_character_), # necessary to avoid drawing panel outline
        panel.grid.major = element_blank(), # get rid of major grid
        panel.grid.minor = element_blank(), # get rid of minor grid
        #plot.background = element_rect(fill = "transparent",colour = NA_character_), # necessary to avoid drawing plot outline
        #legend.background = element_rect(fill = "transparent"),
        #legend.box.background = element_rect(fill = "transparent"),
        #legend.key = element_rect(fill = "transparent")
  ) + 
  xlim(-75,100)+
  ylim(-50,50) +
  scale_color_manual(values = c("#BFD3E6FF", "#0570B0FF", "#08519CFF")) +
  scale_shape_manual(values = c('Spring'=15, 'Summer'=17, 'Winter'=16, 'Autumn' = 18)) +
  labs(x = "NMDS1", colour = "Site", y = "NMDS2", shape = "Season") 

xy


ggsave("../Figures/NMDS_v2.png", plot = xy, bg = "transparent", width = 8, height = 8, units = "in", dpi = 600)
###############Corrplot##############
#define the columns that contain your abundance data. Change the number after the ":" to subset your data
df = read.delim("../Input_data/metadata.txt", header = T)
df = df[-c(55:69),]
dim(df)

library(lubridate)
res <- hms(df$Daylight.hours)        # format to 'hours:minutes:seconds'
df$Daylight <- hour(res)*60 + minute(res) + seconds(res)/60
df$Ammonium[is.na(df$Ammonium)] = 0
df$Ammonium[is.na(df$Ammonia)] = 0
df$DOC[is.na(df$DOC)] = 0
df
env = df[,c(6,7,9,10,11,13,14,15,16,20,21,24,28)]
env2 <- env %>% group_by(Sample_Season) %>% 
  dplyr::summarise(across(everything(), mean)) %>% t %>% data.frame() %>% row_to_names(1)
dim(env2)
dim(data.phm5)
cor_topOTUs <- data.frame(cbind(t(data.phm5), t(env2)), check.names = F)
str(cor_topOTUs)
cor_topOTUs$rowname <- rownames(cor_topOTUs)
cor_topOTUs[,-63] <- sapply(cor_topOTUs[,-63], as.numeric)
colnames(cor_topOTUs)
cc = cor(cor_topOTUs[,c(51:62)], method = "spearman", use = "complete.obs")
dim(cc)
###############################
png(filename = "../Figures/corrplot-v1.png", width = 6, height = 6, units = "in", res = 400)

corrplot(cc, method = 'square', diag = FALSE, order = 'hclust', tl.col = "black",
         rect.col = 'black', rect.lwd = 2, tl.offset=0.5,tl.cex = 1.25, addrect = 3,
         col=colorRampPalette(c("#CCEBC5FF", "#A8DDB5FF",  "#4EB3D3FF", "#0868ACFF"))(200))
dev.off()

####################################################
library("pheatmap")
select <- order(rowMeans(counts(dds,normalized=TRUE)),
                decreasing=TRUE)
df <- as.data.frame(colData(dds)[,c("Site", "Month")]) 
df <- as.data.frame(rownames(colData))
DF <- assay(vsd)[select,]
colnames(DF) <- paste0(colData$Site)
new_df_order <- DF[,c(1,2,3,10,11,12,19,20,21,28,29,30,37,38,39,46,47,48,
                          4,5,6,13,14,15,22,23,24,31,32,33,40,41,42,49,50,51,
                          7,8,9,16,17,18,25,26,27,34,35,36,43,44,45,52,53,54)]
colData_ordered <- colData[c(1,2,3,10,11,12,19,20,21,28,29,30,37,38,39,46,47,48,
                             4,5,6,13,14,15,22,23,24,31,32,33,40,41,42,49,50,51,
                             7,8,9,16,17,18,25,26,27,34,35,36,43,44,45,52,53,54),]
data.hm <- new_df_order
dim(data.hm)
data.phm2 <- data.hm
colnames(data.phm2) <- colData_ordered$Sample_Season
data.phm3 <- data.phm2 %>% t() %>% as_tibble() %>%  rownames_to_column()
data.phm3$type <- colData_ordered$Sample_Season
data.phm4 <- data.phm3 %>% group_by(type) %>% 
  dplyr::summarise(across(everything(), mean)) %>% 
  select(-rowname) %>% t %>% as.data.frame() %>% row_to_names(1) 
str(data.phm4)

data.phm4[,-13] <- sapply(data.phm4[,-13], as.numeric)
data.phm5 <- (data.phm4[1:51,-13])

data.phm6 <- data.phm5[order(rownames(data.phm5)), ]

df_hm_s2 <- data.frame(data.phm6)
df_hm_s2$taxonomy <- rownames(data.phm6)
df_hm_s2[,-13] <- sapply(df_hm_s2[,-13], as.numeric)
df_hm_s3 <- (df_hm_s2)
df_hm_s3$taxonomy <- gsub("Unassigned;", "", df_hm_s3$taxonomy, fixed=TRUE)
#df_hm_s3$taxonomy <- gsub("d__Archaea;", "", df_hm_s3$taxonomy, fixed=TRUE)
df_hm_s4 <- separate(df_hm_s3,taxonomy,into = c("Root", "Domain", "Phylum","Class", "Order", "Family", "Genus"),sep = ";",remove = FALSE,extra = "merge")
df_hm_s4$Phylum <- gsub("unclassified", "Root", df_hm_s4$Phylum , fixed=TRUE)
df_hm_s4$Genus[is.na(df_hm_s4$Genus)] <- "g__"
df_hm_s4$Genus <- gsub("Unassigned", "g__", df_hm_s4$Genus , fixed=TRUE)
df_hm_s4$Genus <- gsub("g__", "", df_hm_s4$Genus , fixed=TRUE)
df_hm_s4$newtax <- factor(paste0("(",df_hm_s4$Order,")", ";", df_hm_s4$Genus))
df_hm_s4$newtax<- gsub(" ", "", df_hm_s4$newtax , fixed=TRUE)
df_hm_s5 <- df_hm_s4[,c(1:12)]
rownames(df_hm_s5) <- df_hm_s4$newtax


p <- pheatmap(df_hm_s5[-51,], cluster_rows=FALSE, show_rownames=TRUE,
              fontsize_col = 12, show_colnames = TRUE,
              cluster_cols=FALSE, 
              annotation_legend = TRUE,
              legend = FALSE,
              annotation_names_col = TRUE,
              gaps_col = c(4,8), angle_col = 90,
              gaps_row = c(2,3,7,13,14,15,16,17,18,23,26,30,34,42,48,49),
              cellwidth = 15, cellheight = 13, 
              #filename = "../Figures/Top50OTUs_Mar23_v1.png",
              border_color = "black",
              color = hcl.colors(6, "GnBu"))

p

dev.off()

################Correlations
df = read.delim("../Input_data/metadata.txt", header = T)
df = df[-c(55:69),]
dim(df)

library(lubridate)
res <- hms(df$Daylight.hours)        # format to 'hours:minutes:seconds'
df$Daylight <- hour(res)*60 + minute(res) + seconds(res)/60
df$Ammonium[is.na(df$Ammonium)] = 0
df$Ammonium[is.na(df$Ammonia)] = 0
df$DOC[is.na(df$DOC)] = 0
df
env = df[,c(6,7,9,10,11,13,14,15,16,20,21,24,28)]
env2 <- env %>% group_by(Sample_Season) %>% 
  dplyr::summarise(across(everything(), mean)) %>% t %>% data.frame() %>% row_to_names(1)
dim(env2)

cor_topOTUs <- data.frame(cbind(t(data.phm5), t(env2)), check.names = F)
str(cor_topOTUs)
cor_topOTUs$rowname <- rownames(cor_topOTUs)
cor_topOTUs[,-64] <- sapply(cor_topOTUs[,-64], as.numeric)
colnames(cor_topOTUs)
cc = cor(cor_topOTUs[,-64], method = "spearman", use = "complete.obs")
dim(cc)

rr <- rcorr(cc, type=c("spearman"))
rr_p <- data.frame(rr$P, check.names = F)
rrpc <- round(rr_p[-c(52:63), 52:63], digit = 3)
rrpc[rrpc < 0.001] <- "*"
str(rrpc)
rrpc[rrpc >= 0.001] <- ""
rrpc<- rrpc[order(rownames(rrpc)), ]


cc_01 <- cc[-c(52:63),c(52:63)]
dim(cc_01)
dim(rrpc)
cc_01 <- cc_01[order(rownames(cc_01)), ]

df_cc_s2 <- data.frame(cc_01)
df_cc_s2$taxonomy <- rownames(cc_01)
df_cc_s2[,-13] <- sapply(df_cc_s2[,-13], as.numeric)
df_cc_s3 <- (df_cc_s2)
df_cc_s3$taxonomy <- gsub("Unassigned;", "", df_cc_s3$taxonomy, fixed=TRUE)
df_cc_s3$taxonomy <- gsub("d__Archaea;", "", df_cc_s3$taxonomy, fixed=TRUE)
df_cc_s4 <- separate(df_cc_s3,taxonomy,into = c("Root", "Domain", "Phylum","Class", "Order", "Family", "Genus"),sep = ";",remove = FALSE,extra = "merge")
df_cc_s4$Phylum <- gsub("unclassified", "Root", df_cc_s4$Phylum , fixed=TRUE)
df_cc_s4$Genus[is.na(df_cc_s4$Genus)] <- "g__"
df_cc_s4$Genus <- gsub("Unassigned", "g__", df_cc_s4$Genus , fixed=TRUE)
df_cc_s4$Genus <- gsub("g__", "", df_cc_s4$Genus , fixed=TRUE)
df_cc_s4$newtax <- factor(paste0("(",df_cc_s4$Order,")", ";", df_cc_s4$Genus))
df_cc_s4$newtax<- gsub(" ", "", df_cc_s4$newtax , fixed=TRUE)
df_cc_s5 <- df_cc_s4[,c(1:12)]
rownames(df_cc_s5) <- df_cc_s4$newtax


p_corrplot_all <- pheatmap(df_cc_s5[-51,], cluster_rows=FALSE, show_rownames= TRUE, 
                           show_colnames = TRUE, cluster_cols=TRUE, fontsize_col = 12,
                           annotation_legend = TRUE, fontsize_number = 10,
                           legend = FALSE, 
                           display_numbers = rrpc[-51,],  number_color = "white",
                           gaps_col = c(4,8), angle_col = 90,
                           gaps_row = c(2,3,7,13,14,15,16,17,18,23,26,30,34,42,48,49),
                           cellwidth = 15, cellheight = 13, 
                           border_color = "black", treeheight_col = 0,
                           legend_breaks = c(1, 0.4, 0.7, 0, -0.7, -0.4, -1),
                           #filename = "../Figures/Correlations_nutrients_v2.png",
                           color = c( "#d04d34", "#e59e8c" ,"white", "white", "#a3dc9d", "#448631"))

p_corrplot_all

dev.off()

plot_list=list()
plot_list[['p']]=p[[4]]
plot_list[['p_corrplot_all']]= p_corrplot_all[[4]]

p2 <- grid.arrange(grobs=plot_list, nrow = 1, widths=c(4,4))

ggsave("../Figures/FigureS2_v3.png",plot = p2, dpi = 1200, height = 12, width = 11, units = c("in"))

###Calculate relative abundance
wd.rl <- data.frame(wd/colSums(wd)*100)
wd.rl2 <- wd.rl[rowSums(wd.rl[])>0,]

####Filter out data

###RT = 0.1%

RT <- filter_all(wd.rl2, all_vars(. < 0.01))
AT <- filter_all(wd.rl2, all_vars(. >= 0.1))
AT_any <- AT <- filter_all(wd.rl2, any_vars(. >= 0.1))

Pelagis <- dplyr::filter(AT_any, grepl('f__Pelagibacteraceae', rownames(AT_any)))
Flavos <- dplyr::filter(AT_any, grepl('f__Flavobacteriaceae', rownames(AT_any)))
Rhodos <- dplyr::filter(AT_any, grepl('f__Rhodobacteraceae', rownames(AT_any)))
Pseudos<- dplyr::filter(AT_any, grepl('f__HTCC2089', rownames(AT_any)))
SAR86 <- dplyr::filter(AT_any, grepl('f__D2472', rownames(AT_any)))

Selected_taxa <- rbind(Pelagis, Flavos, Rhodos, Pseudos,SAR86)

df_selected <- Selected_taxa[,c(1,2,3,10,11,12,19,20,21,28,29,30,37,38,39,46,47,48,
                      4,5,6,13,14,15,22,23,24,31,32,33,40,41,42,49,50,51,
                      7,8,9,16,17,18,25,26,27,34,35,36,43,44,45,52,53,54)]
colData_ordered <- colData[c(1,2,3,10,11,12,19,20,21,28,29,30,37,38,39,46,47,48,
                             4,5,6,13,14,15,22,23,24,31,32,33,40,41,42,49,50,51,
                             7,8,9,16,17,18,25,26,27,34,35,36,43,44,45,52,53,54),]
colnames(df_selected) <- colData_ordered$Sample_Season
rownames(df_selected) <- rownames(df_selected)
df_s1 <- df_selected %>% t() %>% as_tibble() %>%  rownames_to_column()
df_s1$type <- colData_ordered$Sample_Season
df_s2<- df_s1 %>% group_by(type) %>% 
  dplyr::summarise(across(everything(), mean)) %>% 
  select(-rowname) %>% t %>% as.data.frame() %>% row_to_names(1) 
str(df_s2)
df_s2$taxonomy <- rownames(df_s2)
df_s2[,-13] <- sapply(df_s2[,-13], as.numeric)
df_s3 <- (df_s2)
df_s3$taxonomy <- gsub("Unassigned;", "", df_s3$taxonomy, fixed=TRUE)
df_s3$taxonomy <- gsub("d__Archaea;", "", df_s3$taxonomy, fixed=TRUE)
df_s4 <- separate(df_s3,taxonomy,into = c("Root", "Domain", "Phylum","Class", "Order", "Family", "Genus"),sep = ";",remove = FALSE,extra = "merge")
df_s4$Phylum <- gsub("unclassified", "Root", df_s4$Phylum , fixed=TRUE)
df_s4$Genus[is.na(df_s4$Genus)] <- "g__"
df_s4$Genus <- gsub("Unassigned", "g__", df_s4$Genus , fixed=TRUE)
df_s4$Genus <- gsub("g__", "", df_s4$Genus , fixed=TRUE)
df_s4$newtax <- factor(paste0("(",df_s4$Order,")", ";", df_s4$Genus))
df_s4$newtax<- gsub(" ", "", df_s4$newtax , fixed=TRUE)
df_s5 <- df_s4[,c(1:12)]
rownames(df_s5) <- df_s4$newtax

df_s6<- filter_all(df_s2[,-13], any_vars(. >= 0.1))


data.reordered <- data[order(rowSums(data),decreasing=T),]
data.reordered2 <- data.reordered[rownames(df_s6),]
colnames(data.reordered2)
new_df_o <- data.reordered2[,c(1,2,3,10,11,12,19,20,21,28,29,30,37,38,39,46,47,48,
                      4,5,6,13,14,15,22,23,24,31,32,33,40,41,42,49,50,51,
                      7,8,9,16,17,18,25,26,27,34,35,36,43,44,45,52,53,54)]
colnames(new_df_o) <- colData_ordered$Site
data.phm3 <- new_df_o %>% t() %>% as_tibble() 
data.phm3$type <- colData_ordered$Site
data.phm4 <- data.phm3 %>% group_by(type) %>% 
  dplyr::summarise(across(everything(), mean)) %>% 
  t %>% as.data.frame() %>% row_to_names(1) 
str(data.phm4)
data.phm4$taxonomy <- rownames(data.phm4)
data.phm4[,-4] <- sapply(data.phm4[,-4], as.numeric)


df_s3 <- (data.phm4)
df_s3$taxonomy <- gsub("Unassigned;", "", df_s3$taxonomy, fixed=TRUE)
df_s3$taxonomy <- gsub("d__Archaea;", "", df_s3$taxonomy, fixed=TRUE)
df_s4 <- separate(df_s3,taxonomy,into = c("Root", "Domain", "Phylum","Class", "Order", "Family", "Genus"),sep = ";",remove = FALSE,extra = "merge")
df_s4$Phylum <- gsub("unclassified", "Root", df_s4$Phylum , fixed=TRUE)
df_s4$Genus[is.na(df_s4$Genus)] <- "g__"
df_s4$Genus <- gsub("Unassigned", "g__", df_s4$Genus , fixed=TRUE)
df_s4$Genus <- gsub("g__", "", df_s4$Genus , fixed=TRUE)
df_s4$newtax <- factor(paste0("(",df_s4$Order,")", ";", df_s4$Genus))
df_s4$newtax<- gsub(" ", "", df_s4$newtax , fixed=TRUE)


df_plot_s1 <- df_s4[,c(1:3)]
rownames(df_plot_s1) <- df_s4$newtax
df_plot_s1 <-  df_plot_s1[order(row.names(df_plot_s1)), ]
############Sites##################
data.m <- as.matrix(df_plot_s1[nrow(df_plot_s1):1,])
data_percentage <- apply(data.m, 1, function(x){x*100/sum(x,na.rm=T)})
col <- brewer.pal(3, "Spectral")
mdat = melt(data_percentage, id.vars=c("BR1", "BR1", "BR1"))

mdat.S <- mdat[1:24,]
mdat.S <- arrange(mdat.S, Var1, desc(value))
mdat.R <- mdat[25:51,]
mdat.R <- arrange(mdat.R, Var1, desc(value))
mdat.H <- mdat[52:69,]
mdat.H <- arrange(mdat.H, Var1, desc(value))
mdat.P <- mdat[70:93,]
mdat.P <- arrange(mdat.P, Var1, desc(value))
mdat.F <- mdat[94:150,]
mdat.F <- arrange(mdat.F, Var1, desc(value))



mdat.o <- rbind(mdat.F, mdat.P, mdat.R,mdat.H,mdat.S)
mdat.o$Var2 <- factor(mdat.o$Var2)
mdat.o$Var1 <- factor(mdat.o$Var1, levels = c("BR1", "BR1", "BR1"))

plot_1 = ggplot(mdat.o, aes(x=Var2, y=value, fill=Var1)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), legend.position = "bottom", legend.text = element_text(size = 12),
        #axis.text.y = element_text(hjust=0, size = 12), axis.ticks.y = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_blank()) +
  geom_bar(position = "stack", stat="identity")+
  geom_text(aes(label = ifelse(value > 2, paste0(format(round(value, digits = 1)),"%" ), " "),
                hjust = "centre"), 
            position = position_stack(vjust = 0.5), size = 3)+
  #geom_col(position = position_stack(reverse = TRUE)) +
  #scale_fill_brewer(palette ="Blues") +
  scale_fill_manual(values = c('BAY' = "#BFD3E6FF", 'BR1' =  "#0570B0FF",  'BR2'  = "#08519CFF"),
                    breaks = c("BAY", "BR1", "BR2"))+
  ylab("Mean relative proportion of OTUs") + 
  #xlab("Top50 OTUs") +
  guides(fill=guide_legend(title="Sites", reverse = TRUE)) + 
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank())+
  #theme(panel.background = element_rect(fill = "transparent"),
  #plot.background = element_rect(fill = "transparent",
  #colour = NA_character_),
  #legend.background = element_rect(fill = "transparent")) +
  coord_flip()
#theme(plot.margin = unit(c(1, 5, 1, 1), "lines"))
#scale_y_discrete(position = "top")
plot_1


###############Seasons##################
data.reordered <- data[order(rowSums(data),decreasing=T),]
#data.reordered2 <- data.reordered[rownames(df_s6),]
colnames(data.reordered2)
new_df_o <- data.reordered2[,c(1,2,3,10,11,12,19,20,21,28,29,30,37,38,39,46,47,48,
                               4,5,6,13,14,15,22,23,24,31,32,33,40,41,42,49,50,51,
                               7,8,9,16,17,18,25,26,27,34,35,36,43,44,45,52,53,54)]
colnames(new_df_o) <- colData_ordered$Season
data.phm3 <- new_df_o %>% t() %>% as_tibble() 
data.phm3$type <- colData_ordered$Season
data.phm4 <- data.phm3 %>% group_by(type) %>% 
  dplyr::summarise(across(everything(), mean)) %>% 
  t %>% as.data.frame() %>% row_to_names(1) 
str(data.phm4)
data.phm4$taxonomy <- rownames(data.phm4)
data.phm4[,-5] <- sapply(data.phm4[,-5], as.numeric)


df_s3 <- (data.phm4)
df_s3$taxonomy <- gsub("Unassigned;", "", df_s3$taxonomy, fixed=TRUE)
df_s3$taxonomy <- gsub("d__Archaea;", "", df_s3$taxonomy, fixed=TRUE)
df_s4 <- separate(df_s3,taxonomy,into = c("Root", "Domain", "Phylum","Class", "Order", "Family", "Genus"),sep = ";",remove = FALSE,extra = "merge")
df_s4$Phylum <- gsub("unclassified", "Root", df_s4$Phylum , fixed=TRUE)
df_s4$Genus[is.na(df_s4$Genus)] <- "g__"
df_s4$Genus <- gsub("Unassigned", "g__", df_s4$Genus , fixed=TRUE)
df_s4$Genus <- gsub("g__", "", df_s4$Genus , fixed=TRUE)
df_s4$newtax <- factor(paste0("(",df_s4$Order,")", ";", df_s4$Genus))
df_s4$newtax<- gsub(" ", "", df_s4$newtax , fixed=TRUE)


df_plot_s2 <- df_s4[,c(1:4)]
rownames(df_plot_s2) <- df_s4$newtax
df_plot_s2 <-  df_plot_s2[order(row.names(df_plot_s2)), ]
data.m <- as.matrix(df_plot_s2[nrow(df_plot_s2):1,])
data_percentage <- apply(data.m, 1, function(x){x*100/sum(x,na.rm=T)})
col <- brewer.pal(3, "Spectral")

mdat = melt(data_percentage, id.vars=c("BR1", "BR1", "BR1"))
#measure.vars = rownames(data_percentage))
mdat.o <- mdat[order(mdat$Var2, decreasing = T),]
mdat.o$Var2 <- factor(mdat.o$Var2)
mdat.o$Var1 <- factor(mdat.o$Var1, levels = c("Autumn", "Summer", "Spring", "Winter"))

plot_2 = ggplot(mdat.o, aes(x=Var2, y=value, fill=Var1)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), legend.position = "bottom", legend.text = element_text(size = 12),
        axis.text.y = element_text(hjust=0, size = 12), axis.ticks.y = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_blank()) +
  geom_bar(position = "stack", stat="identity")+
  geom_text(aes(label = ifelse(value > 2, paste0(format(round(value, digits = 1)),"%" ), " "),
                hjust = "centre"), 
            position = position_stack(vjust = 0.5), size = 3)+
  #geom_col(position = position_stack(reverse = TRUE)) +
  #scale_fill_brewer(palette ="Greens") +
  scale_fill_manual(values= c('Spring'="#41AE76FF", 'Summer'="#006D2CFF", 'Winter'="#99D8C9FF", 'Autumn' = "#CCECE6FF"),
                    breaks = c("Autumn", "Summer", "Spring", "Winter"))+
  ylab("Mean relative proportion of OTUs") + 
  #xlab("Top50 OTUs") +
  guides(fill=guide_legend(title="Seasons", reverse = TRUE)) +
  #theme(panel.background = element_rect(fill = "transparent"),
  #plot.background = element_rect(fill = "transparent",
  #colour = NA_character_),
  #legend.background = element_rect(fill = "transparent")) +
  coord_flip()+
  #theme(plot.margin = unit(c(1, 5, 1, 1), "lines"))
  scale_x_discrete(position = "top")
plot_2


#################Sampleseasons

#data.reordered <- data[order(rowSums(data),decreasing=T),]
select <- order(rowMeans(counts(dds,normalized=TRUE)),
                decreasing=TRUE)
data.reordered <- assay(vsd)[select,][1:51,]
data.reordered2 <- data[order(rowSums(data),decreasing=T),]
data.reordered3 <- data.reordered2[rownames(data.reordered),]
colnames(data.reordered3)
new_df_o <- data.reordered3[,c(1,2,3,10,11,12,19,20,21,28,29,30,37,38,39,46,47,48,
                               4,5,6,13,14,15,22,23,24,31,32,33,40,41,42,49,50,51,
                               7,8,9,16,17,18,25,26,27,34,35,36,43,44,45,52,53,54)]
colnames(new_df_o) <- colData_ordered$Sample_Season
data.phm3 <- new_df_o %>% t() %>% as_tibble() 
data.phm3$type <- colData_ordered$Sample_Season
data.phm4 <- data.phm3 %>% group_by(type) %>% 
  dplyr::summarise(across(everything(), mean)) %>% 
  t %>% as.data.frame() %>% row_to_names(1) 
str(data.phm4)
data.phm4$taxonomy <- rownames(data.phm4)
data.phm4[,-13] <- sapply(data.phm4[,-13], as.numeric)


df_s3 <- data.phm4[order(row.names(data.phm4)), ]
df_s3$taxonomy <- gsub("Unassigned;", "", df_s3$taxonomy, fixed=TRUE)
#df_s3$taxonomy <- gsub("d__Archaea;", "", df_s3$taxonomy, fixed=TRUE)
df_s4 <- separate(df_s3,taxonomy,into = c("Root", "Domain", "Phylum","Class", "Order", "Family", "Genus"),sep = ";",remove = FALSE,extra = "merge")
df_s4$Phylum <- gsub("unclassified", "Root", df_s4$Phylum , fixed=TRUE)
df_s4$Genus[is.na(df_s4$Genus)] <- "g__"
df_s4$Genus <- gsub("Unassigned", "g__", df_s4$Genus , fixed=TRUE)
df_s4$Genus <- gsub("g__", "", df_s4$Genus , fixed=TRUE)
df_s4$newtax <- factor(paste0("(",df_s4$Order,")", ";", df_s4$Genus))
df_s4$newtax<- gsub(" ", "", df_s4$newtax , fixed=TRUE)


df_plot_s2 <- df_s4[,c(5:8)]
rownames(df_plot_s2) <- df_s4$newtax
data.m <- as.matrix(df_plot_s2[nrow(df_plot_s2):1,])
data_percentage <- apply(data.m[-1,], 1, function(x){x*100/sum(x,na.rm=T)})
col <- brewer.pal(3, "Spectral")

mdat = melt(data_percentage)
mdat.o <- mdat[order(mdat$Var2, decreasing = T),]
mdat.o$Var2 <- factor(mdat.o$Var2)
mdat.o$Var1 <- factor(mdat.o$Var1)

plot_3_BR1 = ggplot(mdat.o, aes(x=Var2, y=value, fill=Var1)) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        legend.position = "bottom", 
        legend.text = element_text(size = 12),
        #axis.text.y = element_text(hjust=0, size = 12), 
        axis.ticks.y = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        ) +
  geom_bar(position = "stack", stat="identity")+
  geom_text(aes(label = ifelse(value > 2, paste0(format(round(value, digits = 1)),"%" ), " "),
                hjust = "centre"), 
            position = position_stack(vjust = 0.5), size = 3)+
  scale_fill_manual(values= c('BR1Spring'="#41AE76FF", 'BR1Summer'="#006D2CFF", 'BR1Winter'="#99D8C9FF", 'BR1Autumn' = "#CCECE6FF"),
                    breaks = c("BR1Autumn", "BR1Summer", "BR1Spring", "BR1Winter")
                    )+
  ylab("Mean relative proportion of OTUs") + 
  guides(fill=guide_legend(title="Seasons", reverse = TRUE)) +
  coord_flip()+
  scale_x_discrete(position = "top")
plot_3_BR1

library(cowplot)
plot3 <- cowplot::plot_grid(plot_3_BAY,NULL,plot_3_BR1,NULL,plot_3_BR2 +
                              theme(#axis.text.y = element_blank(),
                                #axis.ticks.y = element_blank(),
                                #axis.title.y = element_blank()), 
                              ),
                            rel_widths = c(3, -1, 3, -1, 3),
                            align = "v", nrow = 1)
plot3


ggsave("../Figures/Top50OTUs_seasons_v2.png", plot = plot3, dpi = 600, units = c("in"), width = 20, height = 10)

```

![pca-ALL](https://github.com/aprabhu90/Brisbane-river-microbiome/assets/80237948/f76cf2ac-b10f-4d29-8646-8b84026595eb)
