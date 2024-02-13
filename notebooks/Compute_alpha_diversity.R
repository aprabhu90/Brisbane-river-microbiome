###########Load the libraries
library(phyloseq)
library(ape)
library(DESeq2)
library(dplyr)
library(ggplot2)
library(gplots)
library(lme4)
library(miLineage)
library(phangorn)
library(plotly)
library(tidyr)
library(vegan)
library(VennDiagram)
library(venneuler)
library(viridis)

data <- read.csv("../Input data/All_species.csv", header = T)
dim(data)
data <- data[-nrow(data),]
tail(data)
data <- data[,-c(2:10,47:52)]
metadata <- read.delim("metadata.txt", header = T)
metadata <- metadata[-c(1:9, 46:51),]
dim(metadata)

#########Euks, Proks seperately
data_Euks <- dplyr::filter(data, grepl('d__Eukaryota', data$Taxonomy))
dim(data_Euks)

archaea.data <- dplyr::filter(data, grepl("d__Archaea", data$Taxonomy))
dim(archaea.data)

bacteria.data <- dplyr::filter(data, grepl("d__Bacteria", data$Taxonomy))
dim(bacteria.data)

#keep <- rowSums(data_Euks[,-1]) >= 10
data <- as_tibble(data_Euks) 
dim(data)
#keep <- rowSums(archaea.data[,-1]) >= 10
data <- as_tibble(archaea.data)
dim(data)
#keep <- rowSums(bacteria.data[,-1]) >= 10
data <- as_tibble(bacteria.data)
dim(data)

data_prok <- dplyr::filter(data, !grepl('d__Eukaryota', data$Taxonomy))
data <- as_tibble(data_prok)

#Prepare a phyloseq object
otu_mat<- data.frame(data)
colnames(otu_mat)[1] <- "Taxonomy"
tax_mat<- data.frame(otu_mat$Taxonomy)
colnames(tax_mat)[1] <- "Taxonomy"
tax <- separate(
  tax_mat,
  col = 1,
  into = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"),
  sep = ";",
  remove = TRUE,
  convert = FALSE,
  extra = "warn",
  fill = "warn")
tax_mat <- cbind(tax_mat, tax)
samples_df <- as_tibble(metadata)
otu_mat <- otu_mat %>%
  tibble::column_to_rownames("Taxonomy") 
tax_mat <- tax_mat %>% 
  tibble::column_to_rownames("Taxonomy")
samples_df <- samples_df %>% 
  tibble::column_to_rownames("SampleID") 

####Convert the dataframes to matrices#####
otu_mat <- as.matrix(otu_mat)
tax_mat <- as.matrix(tax_mat)
dim(otu_mat)
dim(tax_mat)

#####Perform the rarefaction#####
rarecurve(t(otu_mat[,colSums(otu_mat) > 1]),step = 1000, label = F,xlim = c(0,40000),sample = 1000)
rarefy_threshold <- min(colSums((otu_mat))) 
nrow(otu_mat)
ncol(otu_mat)
species_rare <- t(rrarefy(t(otu_mat[,colSums(otu_mat) >= rarefy_threshold]), rarefy_threshold))
nrow(species_rare)
ncol(species_rare)

####Designate names for rarefied objects in phyloseq#####
OTU = otu_table(species_rare, taxa_are_rows = TRUE)
TAX = tax_table(tax_mat)
samples = sample_data(samples_df)

####Run phyloseq#####
phy_seq <- phyloseq(OTU, TAX, samples)
phy_seq
sample_names(phy_seq)
rank_names(phy_seq)
sample_variables(phy_seq)
SampleID <- row.names(samples_df)
metadata <- cbind(SampleID, samples_df)
####Load information into breakaway

#devtools::install_github("bryandmartin/corncob")
library(corncob)
library(breakaway)
#data("soil_phylo")
#remotes::install_github("adw96/DivNet")
library(DivNet)
#remotes::install_github("mikemc/speedyseq")
library(speedyseq)

phy_seq %>% sample_data %>% head

richness_water <- phy_seq %>% breakaway

plot(richness_water, physeq=phy_seq, color="Site", shape = "Season")

summary(richness_water) %>% as_tibble

meta <- phy_seq %>%
  sample_data %>%
  as_tibble %>%
  mutate("sample_names" = phy_seq %>% sample_names )

combined_richness_water <- meta %>%
  left_join(summary(richness_water),
            by = "sample_names")

water_phylum <- phy_seq %>%
  tax_glom(taxrank="Phylum")

water_family <- phy_seq %>%
  speedyseq::tax_glom(taxrank="Family")

water_species <- phy_seq %>%
  speedyseq::tax_glom(taxrank="Species")

###Run breakaway with Divnet

dv_phylum <- DivNet::divnet(water_phylum, X = NULL)
dv_family <- DivNet::divnet(water_family, X = NULL)
dv_species <- DivNet::divnet(water_species, X = NULL)

dv_species

combined_shannon_water_species_shannon <- meta %>%
  left_join(dv_species$shannon %>% summary,
            by = "sample_names")
combined_shannon_water_species_shannon

combined_shannon_water_species_simpson <- meta %>%
  left_join(dv_species$simpson %>% summary,
            by = "sample_names")


write.csv(combined_shannon_water_species_shannon, file="../Results/prokaryotic-Richness-breakaway-species-shannon.csv")
write.csv(combined_shannon_water_species_simpson, file="../Results/prokaryotic-Richness-breakaway-species-simpson.csv")
###Plot richness

combined_shannon_water_species <- read.csv("../Results/prokaryotic-Richness-breakaway-species-shannon.csv", header = T, row.names = 1)

combined_shannon_water_species2 <- combined_shannon_water_species %>%
                                   group_by(Site, Season) %>%
                                   dplyr::summarise(mean(combined_shannon_water_species$estimate),sd(combined_shannon_water_species$estimate))
combined_shannon_water_species3<- left_join(combined_shannon_water_species, combined_shannon_water_species2, by = c("Site" = "Site"))

my_comparisons_sites <- list( c("BAY", "BR1"), c("BR1", "BR2"), c("BAY", "BR2") )


breakaway_site<-ggplot(combined_shannon_water_species3, aes(x=Site,y=estimate, fill = Site))+
  facet_wrap(~Season.x, ncol = 4) + 
  geom_boxplot(coef = 10) +
  #geom_errorbar(aes(ymin=`mean(combined_shannon_water_species$estimate)`-`sd(combined_shannon_water_species$estimate)`, 
                    #ymax=`mean(combined_shannon_water_species$estimate)` + `sd(combined_shannon_water_species$estimate)`),
                    #width=.6, position = position_dodge(width=2))+
  #geom_line(position = position_dodge(width=2))+
  #geom_point(position = position_dodge(width=2), size=3)+
  theme_classic()+
  theme(legend.position = "right",
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        #axis.title.y = element_text(colour = "black", size = 12)
        axis.title = element_blank()
  ) +
  scale_y_continuous()+
  scale_fill_manual(values = c('BAY' = "#BFD3E6FF", 'BR1' =  "#0570B0FF",  'BR2'  = "#08519CFF"))+
  stat_compare_means(comparisons = my_comparisons_sites, label = "p.signif", 
                     hide.ns = T, symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1, 2), 
                                                     symbols = c("*", "*", "*", "*", ".", " "),
                                                     vjust = 0.1))

breakaway_site

ggsave("../Figures/breakaway_shannon_site_species3.png", plot = breakaway_site, dpi = 600, units = c("mm"), width = 150, height = 60)


breakaway_season<-ggplot(combined_shannon_water_species, aes(x=Season, y=estimate, fill = Site))+
  facet_wrap(~Site, ncol = 4) + 
  geom_boxplot(coef = 20) +
  #geom_errorbar(aes(ymin=estimate-lower, ymax=estimate+upper), width=.6, position = position_dodge(width=2))+
  #geom_line(position = position_dodge(width=2))+
  #geom_point(position = position_dodge(width=2), size=3)+
  theme_bw()+
  theme(legend.position = "right",
        axis.ticks.x = element_blank(),
        axis.text.x = element_text(),
        axis.title.y = element_text(colour = "black", size = 12),
  ) +
  scale_y_continuous()+
  scale_fill_manual(values = c("#0072B2", "#009E73", "#D55E00")) +
  ggtitle("Breakaway alpha estimate")

breakaway_season
ggsave("../Figures/breakaway_shannon_season2_species.png", plot = breakaway_season, dpi = 600, units = c("in"), width = 12, height = 6)
write.csv(combined_shannon_water_species, file="../Results/prokaryotic-Richness-breakaway-species.csv")

###Do statistical testing 

### Need a non parametric test for a non-normalized data
hist(combined_shannon_water_species$estimate)
anova <- aov(estimate ~ Site*Season , data = combined_shannon_water_species)
summary(anova)
TukeyHSD(anova)
kruskal.test(estimate ~ Site , data = combined_shannon_water_species)
capture.output(TukeyHSD(anova), file = "../Results/Anova/breakaway_site+season_anova.txt")
#capture.output(TukeyHSD(anova), file = "../Results/Anova/breakaway_site_TukeyHSD.txt")

anova <- aov(estimate ~ Season , data = combined_shannon_water_species)
summary(anova)
TukeyHSD(anova)
kruskal.test(estimate ~ Site , data = combined_shannon_water_species)
capture.output(TukeyHSD(anova), file = "../Results/Anova/BA_season_anova.txt")
#capture.output(TukeyHSD(anova), file = "../Results/Anova/BD_season_TukeyHSD.txt")

anova <- aov(estimate ~ Site+Season , data = combined_shannon_water_species)
summary(anova)
TukeyHSD(anova)
kruskal.test(estimate ~ Site , data = combined_shannon_water_species)
capture.output(TukeyHSD(anova), file = "../Results/Anova/BD_siteplusseason_anova.txt")
#capture.output(TukeyHSD(anova), file = "../Results/Anova/BD_site+season_TukeyHSD.txt")
