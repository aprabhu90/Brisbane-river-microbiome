library(indicspecies)

#Read the counts table for the grouped dataset
dataset <- read.delim("../Input_data/new_singlem_table_v1.txt", header = T,check.names = F)
dataset$OTU_ID <- format(paste0( "OTU", 1:nrow(dataset)))
wd <- dataset[,-c(1:10, 47:52, 71)]
rownames(wd) <- dataset$OTU_ID
dim(wd)
data <- wd

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
names(data) = gsub(pattern = "_S.*", replacement = "", x = names(data))
names(data) = gsub(pattern = "_COM.*", replacement = "", x = names(data))

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

df <- cbind(t(vsd_norm),metadata)
df[,1:2020] <- sapply(df[,1:2020], as.numeric)
######Nitrates
df$N <- ifelse(df$Nitrate <= 0.02, 'Normal', ifelse(df$Nitrate <= 0.2, 'High(10X)', 'Severe(50X)'))


indval_analysis_N1 = multipatt(df[,1:2020], df$N, 
                            func = 'IndVal.g', duleg = TRUE, control=how(nperm=1000))

summary(indval_analysis_N1)

#Extract data in a table like format
indval_table_N1 = indval_analysis_N1$sign
#Significant p value trimming
indval_table_N2 = indval_table_N1[which(indval_table_N1$p.value <= 0.001),]
#Indval analysis summary
summary(indval_analysis_N1, indvalcomp = TRUE)
#RENAME FILE TO SPECIFIC TREATMENT ANALYZED
capture.output(summary(indval_analysis_N1, indvalcomp = TRUE, alpha = 0.05, At = 0.7, Bt = 0.7), file="../Results/Indval_analyses_N3.csv")
write.csv(indval_table_N2, file="../Results/Indval_analyses_N2_v3.csv")

data_N <- data2 %>% filter(rownames(data2) %in% rownames(indval_table_N2))
write.csv(data_N, file="../Results/Indval_analyses_N_subset1.csv")

#########Phosphorus
df$P <- ifelse(df$Phosphorus <= 0.1, 'Moderately High(3x)', ifelse(df$Phosphorus <= 0.15, 'High(5X)', 'Severe(15X)'))


indval_analysis_P1 = multipatt(df[,1:2020], df$P, 
                               func = 'IndVal.g', duleg = TRUE, control=how(nperm=1000))

summary(indval_analysis_P1)

#Extract data in a table like format
indval_table_P1 = indval_analysis_P1$sign
#Significant p value trimming
indval_table_P2 = indval_table_P1[which(indval_table_P1$p.value <= 0.01),]
#Indval analysis summary
summary(indval_analysis_P1, indvalcomp = TRUE)
#RENAME FILE TO SPECIFIC TREATMENT ANALYZED
capture.output(summary(indval_analysis_P1, indvalcomp = TRUE, alpha = 0.05, At = 0.7, Bt = 0.7), file="../Results/Indval_analyses_P3.csv")
write.csv(indval_table_P2, file="../Results/Indval_analyses_P2_v3.csv")


data_P <- data2 %>% filter(rownames(data2) %in% rownames(indval_table_P2))
write.csv(data_P, file="../Results/Indval_analyses_P_subset1.csv")




write.csv(wd, file="../Results/Data_for_ML.csv")
write.csv(vsd_norm, file="../Results/Data_for_ML_normalized.csv")

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


