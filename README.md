# Metagenomics and machine learning approach 

## Data interpretation
```
vibrant <- read.table("../Data/VIBRANT_all.tsv", quote="\"", comment.char="")
vibrant$V2 <- str_remove(vibrant$V2, "_VIBRANT.tsv")
vibrant$seq_name <- paste0(vibrant$V2, "_", vibrant$V1)
sum(duplicated(vibrant$seq_name))
```

![pca-ALL](https://github.com/aprabhu90/Brisbane-river-microbiome/assets/80237948/f76cf2ac-b10f-4d29-8646-8b84026595eb)
