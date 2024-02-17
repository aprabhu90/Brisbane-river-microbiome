library(reshape2)
library(dplyr)
library(tidyverse)
library(ggridges)
library(ggpattern)

###Prepare a dotplot of completeness contamination

mags <- read.delim("../Input_data/MAGs_orginal.txt", header = T, check.names = F)


###Add a barplot contam to left and completeness to bottom
library(ggplot2)
cc <- ggplot(mags, aes(x = Completeness, y = Contamination)) +
      geom_point(aes(color = factor(type_checkm1)), position = "jitter",size=2, alpha = 0.7) +  
      scale_colour_manual(values=c("Finished" = "#6c93ae","Medium"="#8476c3","High"="#a353d7")) +
      theme_bw() +
      theme(legend.position = "bottom",
            panel.grid = element_blank()) +
      xlab("Completeness(%)") +
      ylab("Contamination(%)") +
      labs(color= "Genome quality") + 
      guides(colour = guide_legend(override.aes = list(size=5)))
cc


library(ggplot2)
cg <- ggplot(mags, aes(x = mags$Est_genome_size/1000000, y = mags$GC)) +
  facet_wrap(~Domain)+
  geom_point(aes(color = factor(Phylum)),position = "jitter",size=2, alpha = 0.7) +  
  geom_smooth(aes(group = Phylum))+
  scale_colour_manual(values=c("#a64557",
                                        "#de475f",
                                        "#dc8576",
                                        "#b64327",
                                        "#e37639",
                                        "#976c30",
                                        "#d0a134",
                                        "#c7ad68",
                                        "#6f782d",
                                        "#93b937",
                                        "#468730",
                                        "#82b773",
                                        "#50c264",
                                        "#307a52",
                                        "#3da08b",
                                        "#59cbb4",
                                        "#4dbcdf",
                                        "#6091d0",
                                        "#5c66b6",
                                        "#775dcf",
                                        "#b18bdc",
                                        "#c65ac9",
                                        "#a2438c",
                                        "#94537b",
                                        "#dd85b5",
                                        "#db4690")) +
  theme_bw() +
  theme(legend.position = "bottom",
        panel.grid = element_blank()) +
  xlab("Estimated genome size") +
  ylab("GC content(%)") +
  labs(color= "Phylum") + 
  xlim(0, 7.5)+
  guides(colour = guide_legend(override.aes = list(size=5)))
cg


library(VennDiagram)

# Generate 3 sets of 200 words
set1 <- paste(rep("word_" , 200) , sample(c(1:1000) , 200 , replace=F) , sep="")
set2 <- paste(rep("word_" , 200) , sample(c(1:1000) , 200 , replace=F) , sep="")
set3 <- paste(rep("word_" , 200) , sample(c(1:1000) , 200 , replace=F) , sep="")

# Chart
venn.diagram(
  x = set_list,
  category.names = c("Archaea ", "Bacteria", 
                     "Unknown species" , "Known species",
                     "Unknown genus" , "Known genus"
                     ),
  filename = '#14_venn_diagramm.png',
  output=TRUE
)

set2 <- mags %>% dplyr::filter(mags$Domain == "d__Archaea") %>% select(`Bin Id`)
set3 <- mags %>% dplyr::filter(mags$Domain == "d__Bacteria") %>% select(`Bin Id`)
set4 <- mags %>% dplyr::filter(mags$Species == "s__") %>% select(`Bin Id`)
set5 <- mags %>% dplyr::filter(mags$Species != "s__") %>% select(`Bin Id`)
set6 <- mags %>% dplyr::filter(mags$Genus == "g__") %>% select(`Bin Id`)
set7 <- mags %>% dplyr::filter(mags$Genus != "g__") %>% select(`Bin Id`)


set_list <- list(set2$`Bin Id`, set3$`Bin Id`, 
                 set4$`Bin Id`, set5$`Bin Id`,
                 set6$`Bin Id`, set7$`Bin Id`
                 )

lab = mags %>% group_by(Sequencing) %>% tally()

compl_mags <- mags %>% dplyr::filter(mags$Completeness >= 70) %>% tally()

cg <- ggplot(mags, aes(x = round(Completeness), fill = Sequencing))+
            geom_bar() +
            theme_bw() +
            theme(plot.background = element_blank(),
                  panel.grid = element_blank(),
                  legend.position = "bottom")+
            xlab("Completeness (%)") +
            ylab("Genomes")+
            guides(colour = guide_legend(override.aes = list(size=5))) +
            scale_fill_manual(values = c("Illumina" = "#869339", "Nanopore" = "#d8434d"))
            
cg
ggsave("../Figures/plotcomp.png", plot = cg, dpi = 600, units = c("in"), width = 5, height = 5)


###Plot mags
mags_4 <- mags
mags_4$Est_genome_size <- mags_4$Est_genome_size/1000000
mags_4$Phylum <- gsub("p__*", "", mags_4$Phylum)

est <- ggplot(mags_4, aes(x = GC,  y= reorder(Phylum, desc(Phylum)), fill = Phylum)) +
  geom_boxplot(coef = 5, alpha = 0.5) + 
  guides(fill=guide_legend(ncol =1)) +
  theme(legend.position = "none", 
        panel.grid = element_blank(),
        axis.text.x = element_text(size = 10,vjust = 0.7),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks = element_blank(),
        panel.background = element_blank()) +
  xlab("GC (%)")  +
  scale_fill_manual(values = c("#e386a3",
                                        "#e14180",
                                        "#d8434d",
                                        "#b05252",
                                        "#d1552b",
                                        "#df976d",
                                        "#9c622b",
                                        "#db9332",
                                        "#c6ad45",
                                        "#7a6e2b",
                                        "#869339",
                                        "#9cb934",
                                        "#90b772",
                                        "#4b7a2e",
                                        "#56b94f",
                                        "#37835c",
                                        "#51c49e",
                                        "#49b9d2",
                                        "#6194d3",
                                        "#626edc",
                                        "#6160a5",
                                        "#be8fd8",
                                        "#9d53c4",
                                        "#d65ebf",
                                        "#985180",
                                        "#b53d7d"))
                                        
est
#ggsave("../Figures/GC2.png", plot = est, dpi = 600, units ="in", width = 5, height = 5)




est_2 <- ggplot(mags_4, aes(x = Est_genome_size,  y= reorder(Phylum, desc(Phylum)), fill = Phylum)) +
  xlim(0.5,8)+
  geom_boxplot(coef = 5, alpha = 0.5) + 
  guides(fill=guide_legend(ncol =1)) +
  theme(legend.position = "none", 
        panel.grid = element_blank(),
        axis.text.x = element_text(size = 10,vjust = 0.7),
        axis.text.y = element_text(size = 10),
        axis.title.y = element_blank(),
        axis.ticks = element_blank(), panel.background = element_blank()) +
  xlab("Estimated genome size (Mb)")  +
  scale_fill_manual(values = c("#e386a3",
                                        "#e14180",
                                        "#d8434d",
                                        "#b05252",
                                        "#d1552b",
                                        "#df976d",
                                        "#9c622b",
                                        "#db9332",
                                        "#c6ad45",
                                        "#7a6e2b",
                                        "#869339",
                                        "#9cb934",
                                        "#90b772",
                                        "#4b7a2e",
                                        "#56b94f",
                                        "#37835c",
                                        "#51c49e",
                                        "#49b9d2",
                                        "#6194d3",
                                        "#626edc",
                                        "#6160a5",
                                        "#be8fd8",
                                        "#9d53c4",
                                        "#d65ebf",
                                        "#985180",
                                        "#b53d7d"))

est_2
#ggsave("../Figures/est_size1.png", plot = est_2, dpi = 600, units ="in", width = 6, height = 10)

library(cowplot)
est3 <- cowplot::plot_grid(est_2, NULL, est,
                  rel_widths = c(1, -0.1, 1),
                  align = "v", nrow = 1)

est3

#ggsave("../Figures/est1.png", plot = est3, dpi = 600, units ="in", width = 8, height = 8)


########Plot sankey###########

library(ggsankey)
library(dplyr)
library(ggplot2)
library(ggforce)

mags <- read.delim("../Input_data/MAGs_orginal.txt", header = T, check.names = F)
P1 <- mags 
P1$Domain <- gsub("d__", "", P1$Domain)
P1$Phylum <- gsub("p__", "", P1$Phylum)
P1$Class <- gsub("c__", "", P1$Class)
P1$Order <- gsub("o__", "", P1$Order)

P2 <- 
  P1%>% 
  dplyr::mutate(Domain_colour = case_when(Domain == "Bacteria"~"#Af4b4e",
                                          T~"#7d579a")) %>%
  dplyr::mutate(Phylum_colour = case_when(Domain == "Bacteria"~"#Ca4fa5",
                                          T~"#E8eb32"),
                Class_colour = case_when(Domain == "Bacteria"~"#4baf75",
                                         T~"#Afad4b"),
                Order_colour = case_when(Domain == "Bacteria"~"#4ba6af",
                                         T~"#B98a2b")
  )

df.AB <- P2[,37:40] # P1
df.AB1 <- df.AB %>% select(Domain, Phylum, Class, Order) %>% make_long(Domain, Phylum, Class, Order)
df.AB1$node <- factor(df.AB1$node)
reorder1 <- data.frame(levels(df.AB1$node))
reorder2 <- data.frame(reorder1[c(11,112,74,73,111,92,91,70,69,80,115,13,50,124,
                                  17,137,136,85,79,105,
                                  93,44,125,123,118,106,95,94,63,53,52,43,45,42,38,31,26,12,1,
                                  10,132,131,128,126,110,114,113,109,108,107,104,101,100,99,98,96,86,82,65,66,51,39,28,25, 
                                  90,89,88,87, #Plancto
                                  83,64,119,81,134, #Pastesci
                                  77,76,75, #Nitrospinota
                                  68,127,67,133,140, #Myxo
                                  58,59,46,129, #Marghulis
                                  62,61,60, #Marini
                                  54,122,130, #3Latesci
                                  49,48,47,57, #3Gemma
                                  40,16,2,97, ##Firmicutes
                                  37,116,
                                  36,14,15,
                                  32,33,84,135,27, ##Cyanos
                                  30,35,117,121, 56,55, #Chloroflex
                                  24,120, ##Bdello_C
                                  23,22,21, ##Bdello
                                  19,103,20,102,18,29,34,41,78, #Bacteroidota
                                  6,9,8,71,72,4,3,7, ##Actino
                                  5,139,138),]) ##Acido
df.AB1$node <- factor(df.AB1$node, levels = rev(reorder2[,1]))

dagg <- df.AB1%>%
  dplyr::group_by(node)%>%
  tally()

df2 <- merge(df.AB1, dagg, by.x = 'node', by.y = 'node', all.x = TRUE)

tax_palette.l <- list()
for (tax_level in c("Domain", "Phylum", "Class", "Order")){
  tax_palette.l <- c(tax_palette.l, colour_list_from_dataframe(P2, tax_level))
}


if(!require("pacman")){
  install.packages("pacman")
}

pacman::p_load("tidyverse", "ggnewscale", "colorspace", "ggplot2", 
               "reshape2", "openxlsx", "viridis", "ggtext",
               "ggsankey", "colorspace")



# Assign colours to a dataframe for specified discrete columns.
assign_colours_to_df <- function(dataframe.df, 
                                 columns, 
                                 my_palette = NULL, 
                                 auto_assign =T, 
                                 my_default_palette=NULL){
  # TODO - make option for default continuous palette
  dataframe.df <- as.data.frame(dataframe.df)
  # Palettes
  colour_palette_soft_8 <- c("#8b90bc","#76cc5d","#9558b7","#d2c351","#cd5f88","#89cab7","#d06842","#858658")
  colour_palette_10 <- c("#8eec45","#0265e8","#f6a800","#bf6549","#486900","#c655a0","#00d1b6","#ff4431","#aeb85c","#7e7fc8")
  colour_palette_15 <- c("#77b642","#7166d9","#cfa240","#b351bb","#4fac7f","#d44891","#79843a","#c68ad4","#d15a2c","#5ba7d9","#ce4355","#6570ba",
                         "#b67249","#9b4a6f","#df8398")
  colour_palette_20 <- c("#0057b4","#7fff56","#d600bc","#d8d500","#e76eff","#019932","#9f8fff","#ffc730","#007fac","#a20019","#06fefd","#ff6782",
                         "#00774c","#e0c8ff","#717a00","#4b2952","#e2ed7d","#46321e","#ffbd76","#ffb4c6")
  colour_palette_30 <- c("#009348","#f579fe","#4fe16e","#b40085","#4d7e00","#4742b4","#f0c031","#016dd9","#d45200","#7499ff","#ef4d2d","#01c9c8",
                         "#f8394b","#88d7a6","#d20063","#c8cc5d","#882986","#fdb95d","#404f8f","#917300","#f3aefc","#5c5800","#ff75c3","#00674a",
                         "#ba001c","#979760","#8b354c","#ff875f","#943105","#cf9478")
  colour_palette_206 <- c("#cfefb4","#7d8b00","#a70079","#552155","#632900","#ffb173","#fbdcf2","#015a6a","#43fdf7","#ff443a","#008186","#3b8aff",
                          "#8b5fff","#ff9777","#4200a9","#85f6fd","#c96000","#36218a","#d28900","#0137d7","#30325b","#ff836b","#008b4f","#21ff9d",
                          "#00794d","#870052","#e9ec4b","#ce006b","#6e0044","#8a6500","#006971","#432e4b","#ca8dff","#f20059","#44ffe2","#00be5c",
                          "#a0d2ff","#1914ab","#4d284e","#59d7ff","#ab9aff","#0151d9","#1de740","#e24500","#9fc400","#610769","#0a4600","#1e365b",
                          "#018f3f","#b15fff","#009c5e","#005290","#506100","#f49aff","#0187c1","#ffb5f4","#daf100","#70081d","#ff9890","#c1baff",
                          "#ffbe5a","#1b3466","#ff2a7f","#ff5d3c","#e47800","#ac6bff","#1f6000","#006627","#4f4000","#dcd6ff","#ffd7c1","#ed2de4",
                          "#a50038","#a5a8ff","#0f2f7f","#b11700","#00e06b","#ffabb8","#015780","#82eaff","#1b2a88","#6f1600","#d3ef9c","#746e00",
                          "#01d851","#625300","#01d799","#96fd6c","#ff5ca1","#7b0017","#004c2b","#baf678","#f8aaff","#007c1b","#01a88a","#a71ed8",
                          "#fb8cff","#840079","#276d00","#556655","#02b0de","#c0efd7","#63193e","#8e9984","#017ac9","#ff925f","#ff63d7","#294100",
                          "#28baff","#5b2523","#35ab00","#69132e","#8a3b00","#a67700","#7fff6a","#002f96","#681a0b","#4d3003","#ff7de6","#0190d8",
                          "#a69700","#ff6282","#d3f266","#ffc4cf","#ffac3c","#d064ff","#d07aff","#c3005d","#9d0067","#0167c1","#8cfe82","#ffd68f",
                          "#8cfcaf","#f50096","#00c2a2","#aa5e00","#02c16d","#4e4bf6","#ffd962","#004793","#93d800","#462a58","#323a03","#4f9eff",
                          "#2b3a25","#2defff","#02edd6","#864e00","#ffc59f","#e7e9ab","#014cc4","#437bff","#00afba","#ff7d82","#8a1ed4","#ff48b3",
                          "#acf7ab","#005550","#7600a6","#bc0028","#00adab","#02dfbf","#ba004c","#004760","#ebc5ff","#0162d7","#9b3900","#5869ff",
                          "#ff6160","#87b6ff","#ff6796","#ff8422","#ff8440","#b500a8","#937fff","#0132bd","#f48e00","#1e8800","#462370","#3e3614",
                          "#9ca800","#efe5bf","#aeb6a0","#d9aaff","#d8ef89","#cec800","#ffb8b3","#4a2c42","#01715b","#b8ebff","#ff9ec0","#ff93ec",
                          "#ffe0aa","#65b300","#6a8b00","#f6e77c","#ff85c0","#5de522","#a5f6ca","#c70077","#5a4149","#a3b700","#ff63c4","#63fecd",
                          "#93f6e7","#01b4a4")
  # Create palette list for each column
  palette_list <- list()
  for (col in columns){
    column_values <- as.character(unique(dataframe.df[,col]))
    column_values <- column_values[!is.na(column_values)]
    n_values <- length(column_values)
    if (is.null(my_default_palette)){
      if (n_values <= 10){
        # default_colour_palette <- setNames(colorRampPalette(colour_palette_10)(n_values), column_values)
        default_colour_palette <- colour_palette_10
      } else if (n_values > 10 & n_values <= 15){
        # default_colour_palette <- setNames(colorRampPalette(colour_palette_15)(n_values), column_values)
        default_colour_palette <- colour_palette_15
      } else if (n_values > 15 & n_values <= 20){
        # default_colour_palette <- setNames(colorRampPalette(colour_palette_20)(n_values), column_values)
        default_colour_palette <- colour_palette_20
      } else if (n_values > 20 & n_values <= 30){
        # default_colour_palette <- setNames(colorRampPalette(colour_palette_30)(n_values), column_values)
        default_colour_palette <- colour_palette_30
      } else{
        # default_colour_palette <- setNames(colorRampPalette(colour_palette_206)(n_values), column_values)
        default_colour_palette <- colour_palette_206
      }
    } else{
      default_colour_palette <- my_default_palette
    }
    
    if (typeof(my_palette) == "character"){
      if (length(my_palette) >= n_values){
        palette_list[[col]] <- my_palette  
      } else{
        stop(paste0("Provided palette has ", length(my_palette), " colours, but column ", col, " has ", n_values, ". Provide a longer palette."))
      }
    } 
    else if (is.null(my_palette) & auto_assign == T){
      palette_list[[col]] <- default_colour_palette
    }
    else if (typeof(my_palette) == "list"){
      if (! col %in% names(my_palette)){
        if (auto_assign == T){
          warning(paste0("Column ", col, " does not have a specified palette. Assigning a default palette"))
          palette_list[[col]] <- default_colour_palette
        } else{
          warning(paste0("Column ", col, " does not have a specified palette. Skipping."))          
        }
      } else{
        palette_list[[col]] <- my_palette[[col]]
      }
    }
    else {
      stop("Provided palette should be a character vector of colours or a list of separate character vectors, or set auto_assign=T")
    }
    # if (typeof(palette_list[[col]]) == "list"){
    #   col_colours <- palette_list[[col]]
    # } else{
    #   col_colours <- setNames(palette_list[[col]][1:n_values], column_values)
    # }
    if (!is.null(names(palette_list[[col]]))){
      col_colours <- palette_list[[col]]
    } else{
      col_colours <- setNames(palette_list[[col]][1:n_values], column_values)
    }
    dataframe.df[,paste0(col, "_colour")] <- as.character(lapply(as.character(dataframe.df[,col]), function(x) as.character(col_colours[x])))      
  }
  dataframe.df
}

#' Return list of colours for variable values from metadata dataframe
#' assumes _colour column present
#' @param variable column (variable) to get colours for
colour_list_from_dataframe <- function(dataframe.df, variable){
  temp <- unique(dataframe.df[,c(variable, paste0(variable,"_colour"))])
  colours.l <- setNames(temp[,2],temp[,1])
  if (is.null(levels(dataframe.df[,variable])) == F){
    colours.l <- colours.l[levels(dataframe.df[,variable])]
  }
  colours.l <- colours.l[!is.na(names(colours.l))]
  colours.l
}

# Load bin data
bin_data.df <-  read.delim("../Input_data/MAGs_orginal.txt", sep = "\t", header = T)

# Clean up taxonomy level names
bin_data.df <-
  bin_data.df %>%
  # dplyr::mutate(across(c(Domain, Phylum, Class, Order, Genus, Species), ~gsub("^[a-z]__", "", .x))) %>%
  dplyr::select(Bin.Id, Domain, Phylum, Class, Order, Genus, Species)

# Some 'Orders' are undefined, make them "Unassigned"
bin_data.df <- 
  bin_data.df %>%
  dplyr::mutate(Order = ifelse(Order == "o__", "Unassigned", Order))

# --------------------------------------------------------------------------------
# Assign colours
# bin_data.df <- assign_colours_to_df(bin_data.df, 
#                      columns = c("Domain", "Phylum", "Class", "Order")
#                      )

# Colours based on domain colour. Gradually get lighter
# bin_data.df <-
#   bin_data.df %>%
#   dplyr::mutate(Domain_colour = case_when(Domain == "d__Bacteria"~"#Af4b4e",
#                                           T~"#7d579a")) %>%
#   dplyr::mutate(Phylum_colour = colorspace::lighten(Domain_colour, .5),
#                 Class_colour = colorspace::lighten(Phylum_colour, .5),
#                 Order_colour = colorspace::lighten(Class_colour, .5))

# Unique colour per level. Different colours for domains
bin_data.df <- 
  bin_data.df %>% 
  dplyr::mutate(Domain_colour = case_when(Domain == "d__Bacteria"~"#Af4b4e",
                                          T~"#7d579a")) %>%
  dplyr::mutate(Phylum_colour = case_when(Domain == "d__Bacteria"~"#Ca4fa5",
                                          T~"#E8eb32"),
                Class_colour = case_when(Domain == "d__Bacteria"~"#4baf75",
                                         T~"#Afad4b"),
                Order_colour = case_when(Domain == "d__Bacteria"~"#4ba6af",
                                         T~"#B98a2b")
  )
# --------------------------------------------------------------------------------

# Make Unassigned orders "grey"
# bin_data.df <-
#   bin_data.df %>% dplyr::mutate(Order_colour = ifelse(Order == "Unassigned", "grey", Order_colour))

# Separate the archaeal and bacterial bin data
# archaea_bin_data.df <- bin_data.df %>% dplyr::filter(Domain == "Archaea")
# bacteria_bin_data.df <- bin_data.df %>% dplyr::filter(Domain == "Bacteria")
archaea_bin_data.df <- bin_data.df %>% dplyr::filter(Domain == "d__Archaea")
bacteria_bin_data.df <- bin_data.df %>% dplyr::filter(Domain == "d__Bacteria")

# bacteria_bin_data.df <-
  # bacteria_bin_data.df %>%
  # dplyr::arrange(Phylum) %>%
#   dplyr::group_by(Phylum) %>%
#   dplyr::arrange(Class, Order,.by_group = T) %>%
#   dplyr::group_by(Phylum, Class) %>%
#   dplyr::arrange(Order, .by_group = T)
#   

# To use geom_sankey the aestethics x, next_x, node and next_node are required. 
# The last stage should point to NA. The aestethics fill and color will affect both nodes and flows.

# Process the bin data to work with ggsankey. Use the provided 'make_long' function
archaea_bin_data.df <-
  archaea_bin_data.df %>% 
  make_long(Domain, Phylum, Class, Order) %>%
  dplyr::group_by(node) %>%
  # dplyr::mutate(Count = n(), Label = paste0(node, " (",Count,")")) 
  dplyr::mutate(Count = n(), Label = gsub("[a-z]__", "", paste0(node, " (",Count,")")))

bacteria_bin_data.df <-
  bacteria_bin_data.df %>% 
  make_long(Domain, Phylum, Class, Order) %>%
  dplyr::group_by(node) %>%
  # dplyr::mutate(Count = n(), Label = paste0(node, " (",Count,")"))
  dplyr::mutate(Count = n(), Label = gsub("[a-z]__", "", paste0(node, " (",Count,")")))

# Combine the separate bin dataframes back together
sankey_data.df <-
  rbind(archaea_bin_data.df,
        bacteria_bin_data.df)

# Set the order of the nodes based on the order in the separate
# bin dataframes
sankey_data.df$node <- 
  factor(sankey_data.df$node, levels = c(archaea_bin_data.df$node %>% unique(),
                                       bacteria_bin_data.df$node %>% unique())
       )

# Create a combined palette for all lineages
tax_palette.l <- list()
for (tax_level in c("Domain", "Phylum", "Class", "Order")){
  tax_palette.l <- c(tax_palette.l, colour_list_from_dataframe(bin_data.df, tax_level))
}

# Create a darker pallete for flow outlines
tax_palette_darker.l <- sapply(tax_palette.l, colorspace::darken, .1)

# Create plot
myplot <- 
  ggplot(sankey_data.df, aes(x = x, 
                             next_x = next_x,
                             node = node, 
                             next_node = next_node, 
                             fill = node, 
                             label = Label,
                             colour = factor(node)
                             )
         ) +
  geom_sankey(flow.alpha = .8, width = 0.2,
              node.color = "black")+
  geom_sankey_label(size = 8, color = "black", label.size = .1) +
  scale_fill_manual(values = tax_palette.l) +
  scale_color_manual(values = tax_palette.l) +
  theme_sankey(base_size = 0) +
  labs(x = NULL) +
  theme(legend.position = "none",
        plot.title = element_text(hjust = .5))
ggsave(plot = myplot, filename = "../Figures/sankey_plot-1.pdf", height = 40, width = 30, device = "pdf")

#
# df2
myplot <- 
  ggplot(df2, aes(x = x, 
                             next_x = next_x,
                             node = node, 
                             next_node = next_node, 
                             fill = node, 
                             label = paste0(node," (", n, ")"),
                             colour = factor(node)
  )
  ) +
  geom_sankey(flow.alpha = .8, width = 0.2,
              node.color = "black")+
  geom_sankey_label(size = 8, color = "black", label.size = .1) +
  #scale_fill_manual(values = tax_palette.l) +
  scale_color_manual(values = tax_palette.l) +
  theme_sankey(base_size = 0) +
  labs(x = NULL) +
  theme(legend.position = "none",
        plot.title = element_text(hjust = .5))
ggsave(plot = myplot, filename = "../Figures/sankey_plot-2.pdf", height = 40, width = 30, device = "pdf")
