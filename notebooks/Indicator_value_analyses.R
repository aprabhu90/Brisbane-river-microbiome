library(indicspecies)
set.seed(1234)
#Read the counts table for the grouped dataset
dataset <- read.csv("../Input_data/singlem_addnames.csv", header = T,check.names = F)
dataset$OTU_ID <- format(paste0( "OTU", 1:nrow(dataset)))
wd <- dataset[,-c(1,56)]
rownames(wd) <- dataset$OTU_ID
dim(wd)
data <- wd
str(data)
data2 <- data[rowSums(data) >=10,]
###Load the metadata file

metadata <- read.delim("../Input_data/metadata.txt", header = T, row.names = 1, check.names = F)
metadata <- metadata[1:54,]
dim(metadata)


colData <- metadata

# It is absolutely critical that the columns of the count matrix and the rows 
# of the column data (information about samples) are in the same order. 
# DESeq2 will not make guesses as to which column of the count matrix belongs to 
# which row of the column data, these must be provided to DESeq2 already in consistent order.


all(rownames(colData) == colnames(data2))

###Create the count matrix with DESeqDataSetFromMatrix
count.data.set <- DESeqDataSetFromMatrix(countData= data2, 
                                         colData=colData, design= ~ Site+Season)
count.data.set

###Run the differential expression analyses with DESeq2
#The results give you the base mean, the log2fold changes and pvalue/adjusted pvalues
dds <- DESeq(count.data.set)
dim(dds)
dds

####Keep those OTUs with hits >=1 since we have low number of OTUs.
#keep <- rowSums(counts(dds)) >=1
#dds <- dds[keep,]
#dim(dds)

#Normalize with variance stabilizing transformation with blind dispersion
vsd <- varianceStabilizingTransformation(dds)
vsd
head(assay(vsd))
dim(vsd)
vsd_norm <- data.frame(assay(vsd))

set.seed(1)
df <- cbind(t(vsd_norm),metadata)
df[,1:1107] <- sapply(df[,1:1107], as.numeric)
######Nitrates
df$N <- ifelse(df$Nitrate <= 0.02, 'Normal', ifelse(df$Nitrate <= 0.2, 'High(10X)', 'Severe(50X)'))


indval_analysis_N1 = multipatt(df[,1:1107], df$N, 
                            func = 'IndVal.g', duleg = TRUE, control=how(nperm=1000))

summary(indval_analysis_N1)

#Extract data in a table like format
indval_table_N1 = indval_analysis_N1$sign
#Significant p value trimming
indval_table_N2 = indval_table_N1[which(indval_table_N1$p.value < 0.001),]
indval_table_N3 <- left_join(indval_table_N2 %>% mutate(OTU_ID = rownames(indval_table_N2)), dataset[,c(1,56)], by=c("OTU_ID"))
write.csv(indval_table_N3, file="../Output/Indval_analyses_N2_withtax.txt", quote = F)
#Indval analysis summary
summary(indval_analysis_N1, indvalcomp = TRUE)
#RENAME FILE TO SPECIFIC TREATMENT ANALYZED
capture.output(summary(indval_analysis_N1, indvalcomp = TRUE, alpha = 0.001, At = 0.7, Bt = 0.7), file="../Output/Indval_analyses_N4.csv")
write.csv(indval_table_N2, file="../Output/Indval_analyses_N2_v4.csv")

data_N <- data2 %>% filter(rownames(data2) %in% rownames(indval_table_N2))
#write.csv(data_N, file="../Output/Indval_analyses_N_subset1.csv")

#########Phosphorus
df$P <- ifelse(df$Phosphorus <= 0.1, 'Moderately High(3x)', ifelse(df$Phosphorus <= 0.15, 'High(5X)', 'Severe(15X)'))


indval_analysis_P1 = multipatt(df[,1:1107], df$P, 
                               func = 'IndVal.g', duleg = TRUE, control=how(nperm=1000))

summary(indval_analysis_P1)

#Extract data in a table like format
indval_table_P1 = indval_analysis_P1$sign
#Significant p value trimming
indval_table_P2 = indval_table_P1[which(indval_table_P1$p.value < 0.05),]
indval_table_P3 <- left_join(indval_table_P2 %>% mutate(OTU_ID = rownames(indval_table_P2)), dataset[,c(1,56)], by=c("OTU_ID"))
write.csv(indval_table_P3, file="../Output/Indval_analyses_P2_withtax.txt", quote = F)
#Indval analysis summary
summary(indval_analysis_P1, indvalcomp = TRUE)
#RENAME FILE TO SPECIFIC TREATMENT ANALYZED
capture.output(summary(indval_analysis_P1, indvalcomp = TRUE, alpha = 0.05, At = 0.7, Bt = 0.7), file="../Output/Indval_analyses_P4.csv")
write.csv(indval_table_P2, file="../Output/Indval_analyses_P2_v4.csv")


data_P <- data2 %>% filter(rownames(data2) %in% rownames(indval_table_P2))
#write.csv(data_P, file="../Output/Indval_analyses_P_subset1.csv")




#write.csv(wd, file="../Output/Data_for_ML.csv")
#write.csv(vsd_norm, file="../Output/Data_for_ML_normalized.csv")


################Correlations

#define the columns that contain your abundance data. Change the number after the ":" to subset your data
select <- order(rowMeans(counts(dds,normalized=TRUE)),
                decreasing=TRUE)
DF <- assay(vsd)[select,]
com = DF
dim(com)

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

cor_topOTUs <- cbind(t(com), env)

str(cor_topOTUs)
cor_topOTUs$rowname <- rownames(cor_topOTUs)
cor_topOTUs[,-c(1108,1121)] <- sapply(cor_topOTUs[,-c(1108,1121)], as.numeric)
colnames(cor_topOTUs)
cc = cor(cor_topOTUs[,-c(1108,1121)], method = "spearman", use = "complete.obs")
dim(cc)


rr <- rcorr(cc, type=c("spearman"))
rr_p <- data.frame(rr$P, check.names = F)
rrpc <- round(rr_p[-c(1108:1119), 1108:1119], digit = 3)
rrpc[rrpc < 0.001] <- "*"
str(rrpc)
rrpc[rrpc >= 0.001] <- ""

rrpc$OTU_ID <- rownames(rrpc)
rrpc_02 <- left_join(rrpc, dataset[,c(1,56)],
                   by = "OTU_ID")
rrpc_03 <- rrpc %>% filter(rownames(rrpc) %in% rownames(indval_table_P2))
rownames(rrpc_03) <- rrpc_03$Taxonomy
rrpc_04 <- rrpc_03[ order(row.names(rrpc_03)), ]

#write.table(rrpc, "../Output/correlations_rrpc_top50_v2.csv", sep = ",", quote = FALSE)

cc_01 <- data.frame(cc[-c(1108:1119),c(1108:1119)])
dim(cc_01)
dim(rrpc)

cc_01$OTU_ID <- rownames(cc_01)
cc_02 <- left_join(cc_01, dataset[,c(1,56)],
                   by = "OTU_ID")

cc_03 <- cc_02 %>% filter(cc_02$OTU_ID %in% rownames(indval_table_P2))
rownames(cc_03) <- cc_03$Taxonomy
cc_04 <- cc_03[ order(row.names(cc_03)), ]

###Phosphorus 

p_corrplot_all <- pheatmap(cc_04[,-c(13:14)], cluster_rows=FALSE, show_rownames= TRUE, 
                           show_colnames = TRUE, cluster_cols=TRUE, fontsize_col = 12,
                           annotation_legend = TRUE, fontsize_number = 5,
                           legend = FALSE,  #annotation_col = annotation_col,
                           display_numbers = rrpc_04[,-13],  number_color = "white",
                           gaps_col = c(4,8), angle_col = 90,
                           #gaps_row = c(2,3,7,13,14,15,16,17,18,23,26,30,34,42,48,49),
                           cellwidth = 15, cellheight = 15, 
                           border_color = "black", treeheight_col = 0,
                           legend_breaks = c(1, 0.4, 0.7, 0, -0.7, -0.4, -1),
                           #filename = "../Figures/Correlations_nutrients_v1.png",
                           #color = colorRampPalette((brewer.pal(n =3, name ="RdBu")))(5))
                           color = c( "#d04d34", "#e59e8c" ,"white", "white", "#a3dc9d", "#448631"))

p_corrplot_all

dev.off()

plot_list=list()
plot_list[['p']]=p[[4]]
plot_list[['p_corrplot_all']]= p_corrplot_all[[4]]

p2 <- grid.arrange(grobs=plot_list, nrow = 1, widths=c(4,4))

ggsave("../Figures/FigureS2_v2.png",plot = p2, dpi = 1200, height = 12, width = 11, units = c("in"))












###############XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX##############################
#########Temp
df$Temp <- ifelse(df$Temperature <= 20, 'Low(<20)', ifelse(df$Temperature <= 25, 'Ambient(20-25C)', 'High(25-30C))'))


indval_analysis_T1 = multipatt(df[,1:3098], df$Temp, 
                               func = 'IndVal.g', duleg = TRUE, control=how(nperm=1000))

summary(indval_analysis_T1)

#Extract data in a table like format
indval_table_T1 = indval_analysis_T1$sign
#Significant p value trimming
indval_table_T2 = indval_table_T1[which(indval_table_T1$p.value <= 0.001),]
#Indval analysis summary
summary(indval_analysis_T1, indvalcomp = TRUE)
#RENAME FILE TO SPECIFIC TREATMENT ANALYZED
capture.output(summary(indval_analysis_T1, indvalcomp = TRUE, alpha = 0.05, At = 0.7, Bt = 0.7), file="../Results/Indval_analyses_T2.csv")
write.csv(indval_table_T2, file="../Results/Indval_analyses_T2_v2.csv")

##############pH
df$pHInd <- ifelse(df$pH <= 8, 'Non-saline(pH<=8)', 'Saline(pH>8)')


indval_analysis_pH1 = multipatt(df[,1:3098], df$pHInd, 
                               func = 'IndVal.g', duleg = TRUE, control=how(nperm=1000))

summary(indval_analysis_pH1)

#Extract data in a table like format
indval_table_pH1 = indval_analysis_pH1$sign
#Significant p value trimming
indval_table_pH2 = indval_table_pH1[which(indval_table_pH1$p.value <= 0.001),]
#Indval analysis summary
summary(indval_analysis_pH1, indvalcomp = TRUE)
#RENAME FILE TO SPECIFIC TREATMENT ANALYZED
capture.output(summary(indval_analysis_pH1, indvalcomp = TRUE, alpha = 0.05, At = 0.7, Bt = 0.7), file="../Results/Indval_analyses_T2.csv")
write.csv(indval_table_pH2, file="../Results/Indval_analyses_pH2_v2.csv")

#######DO
df$DOInd <- ifelse(df$DO <= 7, 'Low DO (DO<7)', 'Normal DO (DO>=7)')


indval_analysis_DO1 = multipatt(df[,1:3098], df$DOInd, 
                                func = 'IndVal.g', duleg = TRUE, control=how(nperm=1000))

summary(indval_analysis_DO1)

#Extract data in a table like format
indval_table_DO1 = indval_analysis_DO1$sign
#Significant p value trimming
indval_table_DO2 = indval_table_DO1[which(indval_table_DO1$p.value <= 0.001),]
#Indval analysis summary
summary(indval_analysis_DO1, indvalcomp = TRUE)
#RENAME FILE TO SPECIFIC TREATMENT ANALYZED
capture.output(summary(indval_analysis_DO1, indvalcomp = TRUE, alDOa = 0.05, At = 0.7, Bt = 0.7), file="../Results/Indval_analyses_T2.csv")
write.csv(indval_table_DO2, file="../Results/Indval_analyses_DO2_v2.csv")


#########Salinity
df$SalInd <- ifelse(df$Salinity <= 30, 'Brackish( <= 30 PPT)', 'Marine(<= 30PPT')


indval_analysis_Sal1 = multipatt(df[,1:3098], df$SalInd, 
                                func = 'IndVal.g', duleg = TRUE, control=how(nperm=1000))

summary(indval_analysis_Sal1)

#Extract data in a table like format
indval_table_Sal1 = indval_analysis_Sal1$sign
#Significant p value trimming
indval_table_Sal2 = indval_table_Sal1[which(indval_table_Sal1$p.value <= 0.001),]
#Indval analysis summary
summary(indval_analysis_Sal1, indvalcomp = TRUE)
#RENAME FILE TO SPECIFIC TREATMENT ANALYZED
capture.output(summary(indval_analysis_Sal1, indvalcomp = TRUE, alSala = 0.05, At = 0.7, Bt = 0.7), file="../Results/Indval_analyses_T2.csv")
write.csv(indval_table_Sal2, file="../Results/Indval_analyses_Sal2_v2.csv")


#########Turbidity
df$TurbInd <- ifelse(df$Turbidity > 10, 'Brackish( > 10 NTU)', 'Coastal(<= 10 NTU)')


indval_analysis_Turb1 = multipatt(df[,1:3098], df$TurbInd, 
                                 func = 'IndVal.g', duleg = TRUE, control=how(nperm=1000))

summary(indval_analysis_Turb1)

#Extract data in a table like format
indval_table_Turb1 = indval_analysis_Turb1$sign
#Significant p value trimming
indval_table_Turb2 = indval_table_Turb1[which(indval_table_Turb1$p.value <= 0.001),]
#Indval analysis summary
summary(indval_analysis_Turb1, indvalcomp = TRUE)
#RENAME FILE TO SPECIFIC TREATMENT ANALYZED
capture.output(summary(indval_analysis_Turb1, indvalcomp = TRUE, alTurba = 0.05, At = 0.7, Bt = 0.7), file="../Results/Indval_analyses_T2.csv")
write.csv(indval_table_Turb2, file="../Results/Indval_analyses_Turb2_v2.csv")

#########TDS
df$TDSInd <- ifelse(df$TDS < 20, 'Brackish( < 20 mg/L)', 'Coastal( > 20 mg/L)')


indval_analysis_TDS1 = multipatt(df[,1:3098], df$TDSInd, 
                                  func = 'IndVal.g', duleg = TRUE, control=how(nperm=1000))

summary(indval_analysis_TDS1)

#Extract data in a table like format
indval_table_TDS1 = indval_analysis_TDS1$sign
#Significant p value trimming
indval_table_TDS2 = indval_table_TDS1[which(indval_table_TDS1$p.value <= 0.001),]
#Indval analysis summary
summary(indval_analysis_TDS1, indvalcomp = TRUE)
#RENAME FILE TO SPECIFIC TREATMENT ANALYZED
capture.output(summary(indval_analysis_TDS1, indvalcomp = TRUE, alTDSa = 0.05, At = 0.7, Bt = 0.7), file="../Results/Indval_analyses_T2.csv")
write.csv(indval_table_TDS2, file="../Results/Indval_analyses_TDS2_v2.csv")


#########Silicon
df$SilInd <- ifelse(df$Silicon > 1, 'Brackish( >1 mg/L)', 'Coastal( <=1 mg/L)')


indval_analysis_Sil1 = multipatt(df[,1:3098], df$SilInd, 
                                 func = 'IndVal.g', duleg = TRUE, control=how(nperm=1000))

summary(indval_analysis_Sil1)

#Extract data in a table like format
indval_table_Sil1 = indval_analysis_Sil1$sign
#Significant p value trimming
indval_table_Sil2 = indval_table_Sil1[which(indval_table_Sil1$p.value <= 0.001),]
#Indval analysis summary
summary(indval_analysis_Sil1, indvalcomp = TRUE)
#RENAME FILE TO SPECIFIC TREATMENT ANALYZED
capture.output(summary(indval_analysis_Sil1, indvalcomp = TRUE, alSila = 0.05, At = 0.7, Bt = 0.7), file="../Results/Indval_analyses_T2.csv")
write.csv(indval_table_Sil2, file="../Results/Indval_analyses_Sil2_v2.csv")

#########################################################################

metadata <- read.delim("../Input_data/metadata_with microgAMBI.txt", header = T)

microgAMBI <- metadata[,c(11,27,13,28,29,30,31)]

microgAMBI$P <- ifelse(microgAMBI$Phosphorus <= 0.1, 'Elevated', ifelse(microgAMBI$Phosphorus <= 0.15, 'High', 'Severe'))
microgAMBI$N <- ifelse(microgAMBI$Nitrate <= 0.02, 'Normal', ifelse(microgAMBI$Nitrate <= 0.2, 'Elevated', 'High'))

coeff <- 10

microgAMBI.allColor <- "lightblue"
PhosphorusColor <- "black"

pcm <- melt(microgAMBI)
m <- # Value used to transform the data
ggplot(microgAMBI, aes(x=seasongroup, group = 1)) +
  
  #geom_line( aes(y=microgAMBI$microgAMBI.all, color = microgAMBI.allColor)) + 
  geom_col(aes(y=microgAMBI$microgAMBI.all, fill = microgAMBI.allColor)) +
  geom_line( aes(y=microgAMBI$Phosphorus*coeff, fill = PhosphorusColor)) + # Divide by 10 to get the same range than the temperature
  geom_hline(yintercept = 0.2, color = "red") +
  scale_y_continuous(
    
    # Features of the first axis
    name = "First Axis",
    
    # Add a second axis and specify its features
    sec.axis = sec_axis(~.*coeff, name="Second Axis")
  ) +
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90))
m

