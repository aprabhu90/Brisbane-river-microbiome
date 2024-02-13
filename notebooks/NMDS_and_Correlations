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

################Correlations with physico-chemical parameters ########################

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











