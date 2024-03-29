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


################Plot the abundance of top 50 OTUs##########
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

#########################Plot genus level relative abundance for top 50 OTUs ###########################
data <- dataset
data.phm <- data_phm
data.reordered <- data[order(rowSums(data),decreasing=T),]
data.ro.50 <- data.reordered[rownames(data_phm),]
new_df <- data.ro.50[order(row.names(data.ro.50)),]
#new_df <- new_df[-51,]
colnames(new_df)
new_df_o <- new_df[,c(1,2,3,10,11,12,19,20,21,28,29,30,37,38,39,46,47,48,
                      4,5,6,13,14,15,22,23,24,31,32,33,40,41,42,49,50,51,
                      7,8,9,16,17,18,25,26,27,34,35,36,43,44,45,52,53,54)]
colnames(new_df_o) <- colData_ordered$Site
rownames(new_df_o) <- rownames(data.phm2)
data.phm3 <- new_df_o %>% t() %>% as_tibble() #%>%  rownames_to_column()
data.phm3$type <- colData_ordered$Site
data.phm4 <- data.phm3 %>% group_by(type) %>% 
  dplyr::summarise(across(everything(), mean)) %>% 
  #select(-rowname) %>% 
  t %>% as.data.frame() %>% row_to_names(1) 
str(data.phm4)
data.phm4$taxonomy <- rownames(data.phm4)
data.phm4[,-4] <- sapply(data.phm4[,-4], as.numeric)
data.phm5 <- (data.phm4[,-4])
dim(data.phm5)


data.phm4$taxonomy <-rownames(data.phm4)
df_s3 <- (data.phm4)
df_s4 <- separate(df_s3,taxonomy,into = c("Domain", "Phylum","Class", "Order", "Family", "Genus"),sep = ";",remove = FALSE,extra = "merge")
df_s4$Genus[is.na(df_s4$Genus)] <- "g__"
df_s4$newtax <- factor(paste0("(",df_s4$Order,")", " ", df_s4$Genus))
rownames(df_s4) <- df_s4$newtax
df_s5 <- df_s4[,1:3]

############Sites##################
data.m <- as.matrix(df_s5[nrow(df_s5):1,])
data_percentage <- apply(data.m, 1, function(x){x*100/sum(x,na.rm=T)})
col <- brewer.pal(3, "Spectral")
mdat = melt(data_percentage, id.vars=c("BAY", "BR1", "BR2"))
#measure.vars = rownames(data_percentage))
mdat.o <- mdat[order(mdat$Var2, decreasing = T),]
mdat.o$Var2 <- factor(mdat.o$Var2)
mdat.o$Var1 <- factor(mdat.o$Var1, levels = c("BR2", "BR1", "BAY"))


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
                    breaks = c("BR2", "BR1", "BAY"))+
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
colnames(new_df_o) <- colData_ordered$Season
rownames(new_df_o) <- rownames(data.phm2)
data.phm3 <- new_df_o %>% t() %>% as_tibble() #%>%  rownames_to_column()
data.phm3$type <- colData_ordered$Season
data.phm4 <- data.phm3 %>% group_by(type) %>% 
  dplyr::summarise(across(everything(), mean)) %>% 
  #select(-rowname) %>% 
  t %>% as.data.frame() %>% row_to_names(1) 
str(data.phm4)
data.phm4$taxonomy <- rownames(data.phm4)
data.phm4[,-5] <- sapply(data.phm4[,-5], as.numeric)
data.phm5 <- (data.phm4[,-5])
dim(data.phm5)

data.phm5$taxonomy <-rownames(data.phm5)
df_s3 <- (data.phm5)
df_s4 <- separate(df_s3,taxonomy,into = c("Domain", "Phylum","Class", "Order", "Family", "Genus"),sep = ";",remove = FALSE,extra = "merge")
df_s4$Genus[is.na(df_s4$Genus)] <- "g__"
df_s4$newtax <- factor(paste0("(",df_s4$Order,")", " ", df_s4$Genus))
rownames(df_s4) <- df_s4$newtax
df_s5 <- df_s4[,1:4]

data.m <- as.matrix(df_s5[nrow(df_s5):1,])
data_percentage <- apply(data.m, 1, function(x){x*100/sum(x,na.rm=T)})
col <- brewer.pal(3, "Spectral")

mdat = melt(data_percentage, id.vars=c("BAY", "BR1", "BR2"))
mdat.o <- mdat[order(mdat$Var2, decreasing = T),]
mdat.o$Var2 <- factor(mdat.o$Var2)
mdat.o$Var1 <- factor(mdat.o$Var1, levels = c("Autumn", "Summer", "Spring", "Winter"))


plot_2 = ggplot(mdat.o, aes(x=Var2, y=value, fill=Var1)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), legend.position = "bottom", legend.text = element_text(size = 12),
        axis.text.y = element_text(hjust=0, size = 14), axis.ticks.y = element_blank(),
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


library(cowplot)
plot3 <- cowplot::plot_grid(plot_1, NULL, plot_2 +
                              theme(#axis.text.y = element_blank(),
                                #axis.ticks.y = element_blank(),
                                #axis.title.y = element_blank()), 
                              ),
                            rel_widths = c(1, -0.45, 1),
                            align = "v", nrow = 1)
plot3




ggsave("../Figures/Sites_and_Seasons_v1.png", 
       #bg = "transparent", 
       plot = plot3, dpi = 1200, height = 13, width = 13, units = c("in"))

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

##########################Plot Archaea and Bacteria relative abundance 
####################Archaea
data <-wd
data.reordered <- data[order(rowSums(data),decreasing=T),]
A.data <- dplyr::filter(data.reordered, grepl('d__Archaea', rownames(data.reordered)))
dim(A.data)
A.data$Taxonomy <- rownames(A.data)
A.o <- data.frame(t(colSums(A.data[10:nrow(A.data), -55])))#/colSums(archaea.data.reordered[,-55]))) 
A.o$Taxonomy <- "Others"
A.10 <- A.data[1:10,]
colnames(A.o) <- colnames(A.10)
A.m.o <- rbind(A.10, A.o)
new_df <- A.m.o[,-55]
colnames(new_df)
rownames(new_df) <- A.m.o$Taxonomy
new_df_o <- new_df[,c(1,2,3,10,11,12,19,20,21,28,29,30,37,38,39,46,47,48,
                      4,5,6,13,14,15,22,23,24,31,32,33,40,41,42,49,50,51,
                      7,8,9,16,17,18,25,26,27,34,35,36,43,44,45,52,53,54)]
colnames(new_df_o) <- colData_ordered$Sample_Season
rownames(new_df_o) <- rownames(new_df)
data.phm3 <- new_df_o %>% t() %>% as_tibble() %>%  rownames_to_column()
data.phm3$type <- colData_ordered$Sample_Season
data.phm4 <- data.phm3 %>% group_by(type) %>% 
  dplyr::summarise(across(everything(), mean)) %>% 
  select(-rowname) %>% t %>% as.data.frame() %>% row_to_names(1) 
str(data.phm4)
data.phm4$taxonomy <- rownames(data.phm4)
data.phm4[,-13] <- sapply(data.phm4[,-13], as.numeric)
data.phm5 <- (data.phm4)
data.phm5$taxonomy <- gsub("Unassigned;", "", data.phm5$taxonomy, fixed=TRUE)
data.phm5$taxonomy <- gsub("d__Archaea;", "", data.phm5$taxonomy, fixed=TRUE)
data.phm6 <- separate(data.phm5,taxonomy,into = c("Domain", "Phylum","Class", "Order", "Family", "Genus"),sep = ";",remove = FALSE,extra = "merge")
data.phm6$Phylum <- gsub("unclassified", "Root", data.phm6$Phylum , fixed=TRUE)
data.phm6$Genus[is.na(data.phm6$Genus)] <- "g__"
data.phm6$Genus <- gsub("Unassigned", "g__", data.phm6$Genus , fixed=TRUE)
data.phm6$newtax <- factor(paste0("(",data.phm6$Phylum,")", ";", data.phm6$Genus))
data.phm6$newtax<- gsub(" ", "", data.phm6$newtax , fixed=TRUE)
data.phm7 <- data.phm6[,c(1:12,20)]
rownames(data.phm7) <- data.phm7$newtax

data.phm8 <- read.delim("../Results/for_plot_RL.txt", header = T)
data.phm9 <- data.phm8[,-c(1,14:19)] %>% relocate(Taxa)

A <- data.phm9[c(21:40),]
B<-data.phm9[c(1:20),]
A.m <- melt(A)
colourCount = length(unique(A.m$Taxa))
getPalette = colorRampPalette(rev(brewer.pal(6, "Spectral")))


plot_A1 = ggplot(A.m, aes(x=variable, y=value, fill=Taxa)) +
  theme_bw()+
  theme(panel.border = element_blank(),
        legend.position = "right",
        axis.text.x = element_text(angle = 45, size = 14, hjust = 1, vjust = 1.5),
        axis.ticks.x = element_blank(),
        legend.text = element_text(size=15),
        panel.background = element_rect(fill = "white",
                                        colour = NA_character_), # necessary to avoid drawing panel outline
        panel.grid.major = element_blank(), # get rid of major grid
        panel.grid.minor = element_blank(), # get rid of minor grid
        plot.background = element_rect(fill = "white",
                                       colour = NA_character_))+ # necessary to avoid drawing plot outline
  #legend.background = element_rect(fill = "transparent"),
  #legend.box.background = element_rect(fill = "transparent"),
  #legend.key = element_rect(fill = "transparent"))+
  #axis.text.x = element_text(size = 10, angle = 90)) +
  geom_bar(position = "fill", stat="identity", color = "black")+
  #scale_fill_viridis(option = "turbo", discrete = T) + 
  paletteer::scale_fill_paletteer_d("pals::tol")+ 
  ylab("Mean relative proportion of OTUs") + 
  xlab("BAY") +
  guides(fill=guide_legend(title="OTUs")) 
#scale_y_discrete(position = "top")
plot_A1

ggsave("../Figures/Archaea_v1.png", plot = plot_A1, dpi = 600, units = c("in"), width = 12, height = 10)


#################Bacteria############
data <-wd
data.reordered <- data[order(rowSums(data),decreasing=T),]
B.data <- dplyr::filter(data.reordered, grepl('d__Bacteria', rownames(data.reordered)))
dim(B.data)
B.data$Taxonomy <- rownames(B.data)
B.o <- data.frame(t(colSums(B.data[21:nrow(B.data),-55])))#/colSums(archaea.data.reordered[,-55]))) 
B.o$Taxonomy <- "Others"
B.20 <- B.data[1:20,]
colnames(B.o) <- colnames(B.20)
B.m.o <- rbind(B.20,B.o)
new_df <- B.m.o[,-55]
colnames(new_df)
rownames(new_df) <- B.m.o$Taxonomy
new_df_o <- new_df[,c(1,2,3,10,11,12,19,20,21,28,29,30,37,38,39,46,47,48,
                      4,5,6,13,14,15,22,23,24,31,32,33,40,41,42,49,50,51,
                      7,8,9,16,17,18,25,26,27,34,35,36,43,44,45,52,53,54)]
colnames(new_df_o) <- colData_ordered$Sample_Season
rownames(new_df_o) <- rownames(new_df)
data.phm3 <- new_df_o %>% t() %>% as_tibble() %>%  rownames_to_column()
data.phm3$type <- colData_ordered$Sample_Season
data.phm4 <- data.phm3 %>% group_by(type) %>% 
  dplyr::summarise(across(everything(), mean)) %>% 
  select(-rowname) %>% t %>% as.data.frame() %>% row_to_names(1) 
str(data.phm4)
data.phm4$taxonomy <- rownames(data.phm4)
data.phm4[,-13] <- sapply(data.phm4[,-13], as.numeric)
data.phm5 <- (data.phm4)
data.phm5$taxonomy <- gsub("Unassigned;", "", data.phm5$taxonomy, fixed=TRUE)
data.phm5$taxonomy <- gsub("d__Bacteria;", "", data.phm5$taxonomy, fixed=TRUE)
data.phm6 <- separate(data.phm5,taxonomy,into = c("Root", "Phylum","Class", "Order", "Family", "Genus"),sep = ";",remove = FALSE,extra = "merge")
data.phm6$Phylum <- gsub("Unassigned", "Root", data.phm6$Phylum , fixed=TRUE)
data.phm6$Genus[is.na(data.phm6$Genus)] <- "g__"
data.phm6$Genus <- gsub("Unassigned", "g__", data.phm6$Genus , fixed=TRUE)
data.phm6$Family[is.na(data.phm6$Family)] <- "f__"
data.phm6$newtax <- factor(paste0("(", data.phm6$Phylum,")", ";", data.phm6$Family, ";", data.phm6$Genus))
data.phm6$newtax<- gsub(" ", "", data.phm6$newtax , fixed=TRUE)
data.phm7 <- data.phm6[,c(1:12,20)]
rownames(data.phm7) <- data.phm7$newtax

B.m <- melt(data.phm7)
colourCount = length(unique(B.m$newtax))
getPalette = colorRampPalette(rev(brewer.pal(6, "Spectral")))

plot_B = ggplot(B.m, aes(x=variable, y=value, fill=newtax)) +
  theme_bw()+
  theme(panel.border = element_blank(),
        legend.position = "right",
        axis.text.x = element_text(angle = 45, size = 14, hjust = 1, vjust = 1.5),
        axis.ticks.x = element_blank(),
        legend.text = element_text(size=15),
        panel.background = element_rect(fill = "white",
                                        colour = NA_character_), # necessary to avoid drawing panel outline
        panel.grid.major = element_blank(), # get rid of major grid
        panel.grid.minor = element_blank(), # get rid of minor grid
        plot.background = element_rect(fill = "white",
                                       colour = NA_character_))+ # necessary to avoid drawing plot outline
  #legend.background = element_rect(fill = "transparent"),
  #legend.box.background = element_rect(fill = "transparent"),
  #legend.key = element_rect(fill = "transparent"))+
  #axis.text.x = element_text(size = 10, angle = 90)) +
  geom_bar(position = "fill", stat="identity", color = "black")+
  #scale_fill_viridis(option = "turbo", discrete = T) + 
  paletteer::scale_fill_paletteer_d("pals::stepped")+ 
  ylab("Mean relative proportion of OTUs") + 
  xlab("BAY") +
  guides(fill=guide_legend(title="OTUs", ncol = 1)) 
#scale_y_discrete(position = "top")
plot_B

ggsave("../Figures/Bacteria_v1.png", plot = plot_B, dpi = 600, units = c("in"), width = 12, height = 10)


##############################################################
plotAB <- ggpubr::ggarrange(plot_A1, plot_B, heights = c(3, 4), nrow = 1, align = "v")
plotAB

ggsave("../Figures/AB_v1.png", plot = plotAB, dpi = 600, units = c("in"), width = 25, height = 10)
