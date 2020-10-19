#' Downstream analysis
#' 
#' This a script used for the analysis of the tool's
#' output that has been followed.
#' 
#' @param k value indicating the maximum limit of the range investigated
#' 
#' @param count.table path to the output k-mer count table
#' 
#' @param meta path to meta-data file to be used
#' 
#' @param Ncounts path to file containing the N counts of the sequences
#' 
#' @param Acounts path to file containing the A counts of the sequnces
#' 
#' @seq.clustering path to sequence clustering file
#' 
#' 


################### Loading libraries #############

library(data.table)
library(tidyverse)
library(caret)
library(reshape2)
library(ggbiplot)
library(stringdist)
library(ggdendro)
library(dendextend)
library(ggpubr)
library(matrixStats)
library(ComplexHeatmap)
library(gplots)

################### Reading Data #############

k = 80

# Read count table
count.table = "results_k_4-10/ClusteringData/clustData.csv"

kfeatures = fread(count.table)

kmers.map = data.table(k_ID = paste("k", 1:(ncol(kfeatures) - 1), sep = ""),
                       kmer = colnames(kfeatures)[2:ncol(kfeatures)])

colnames(kfeatures) = c("seqIndex", kmers.map$k_ID)
kfeatures = as.data.frame(kfeatures)
row.names(kfeatures) = kfeatures$seqIndex


# Read meta data file
meta = "SARS-meta-data-file_modified.txt"

meta = fread(meta, fill = TRUE, na.strings = c("","NA"))

meta = meta %>% drop_na(ID) %>% drop_na(Region) %>% drop_na(collectionDate)
list = meta$ID

meta$month = str_split(meta$collectionDate, "-", simplify = TRUE)[,2]
meta$month = as.numeric(meta$month)
meta[which(meta$month == 12), ]$month = 0

kfeatures = kfeatures[which(kfeatures$seqIndex %in% list),]


# clean N
Ncounts = "SARS-N_counts.csv"
Ncounts = fread(Ncounts)
Ncounts = Ncounts %>% filter(N_counts == 0)
list = Ncounts$ID

kfeatures = kfeatures[which(kfeatures$seqIndex %in% list),]

meta = meta[which(meta$ID %in% list),]


# clean A
Acounts = "SARS-A_counts.csv"
Acounts = fread(Acounts)
Acounts = Acounts %>%  filter(A_counts == "[6, 5]" | A_counts == "[5]" | A_counts == "[7, 6, 5]" )
list = Acounts$ID

kfeatures = kfeatures[which(kfeatures$seqIndex %in% list),]

meta = meta[which(meta$ID %in% list),]

# keep only counts
kfeatures = kfeatures[,2:ncol(kfeatures)]


rm(list, Acounts, Ncounts)



################### Zero- and Near Zero-Variance Predictors ############

nzv = nearZeroVar(kfeatures)

kfeatures = kfeatures[, -nzv]

rm(nzv)

################### PCA analysis #################################


# Perform PCA
kfeatures.pca = prcomp(kfeatures, center = TRUE, scale. = TRUE)

# get proportion
stats = summary(kfeatures.pca)
importance = stats$importance
importance = as.data.frame(t(importance))

ggplot(data = importance, 
       mapping = aes(x = 1:length(`Proportion of Variance`), 
                     y = `Proportion of Variance`)) + 
  # geom_line() +
  geom_point(color = "red", size = 1.3) + 
  labs(x = "PC")


# PC1 ~ PC2 
x = as.data.table(kfeatures.pca$x[,1:4])
x$ID = row.names(kfeatures.pca$x)

who = match(x$ID, meta$ID)

x$Region = meta[who, ]$Region
x$Day = meta[who, ]$`Day Number`

ggplot(data = x, aes(x = PC1, y = PC2)) + 
  geom_point(aes(color = Day), alpha = 0.4) 

# time correlation
pca_data = x[,1:4]
pca_data$Region = meta$Region
pca_data$Day = meta$`Day Number`
pca.data$source = meta$source



cor.test(pca_data[,2], pca_data$Day, method ="pearson")

ggscatter(pca_data, x = "PC1", y = "Day",fill = "Region", shape = 21,
          add = "reg.line",  
          add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
          conf.int = TRUE) + 
  stat_cor(method = "pearson", label.x = 13, label.y = 1) + geom_smooth()

pca.full = reshape2::melt(pca.data, id.vars = c("ID", "Region", "Day", "source"))

ggscatter(pca.full, x = "PC", y = "Day", fill = "Region", shape = 21,size = 5, alpha = 0.5,
          palette = c("yellow", "red", "blue", "purple", "green", "mediumvioletred"),
          add = "reg.line", 
          add.params = list(color = "black", fill = "gray"), # Customize reg. line
          conf.int = TRUE) + 
  stat_cor(method = "pearson", label.x = 7, label.y = 1) + geom_smooth()+ facet_wrap(~category)


# per country
x = x[,c("PC1", "PC2", "ID", "Region", "Day")]

ggplot(data = x, aes(x = PC1, y = PC2)) +
  geom_point(aes(color = Day), alpha = 0.4) +
  facet_wrap(vars(Region), ncol = 3)

# PCA with directions

gr = ggbiplot(kfeatures.pca, 
              choices = c(1, 2),
              ellipse = FALSE,
              groups = as.factor(meta[who, ]$Region),
              obs.scale = 1,
              var.axes=TRUE, 
              var.scale = 1,
              alpha = 0.4) +
  theme_minimal()+
  theme(legend.position = "bottom")


# directions in every principal component
x = as.data.table(kfeatures.pca$x[,1:4])
rotation = as.data.table(kfeatures.pca$rotation[, 1:4])

colnames(x) = paste(colnames(x), " ", importance[colnames(x), ]$`Proportion of Variance` * 100, "%", sep = "")
colnames(rotation) = colnames(x)

x$ID = row.names(kfeatures.pca$x)
rotation$kmer_ID = row.names(kfeatures.pca$rotation)

who = match(x$ID, meta$ID)
x$Region = meta[who, ]$Region
x$Day = meta[who, ]$`Day Number`

x.melt = reshape2::melt(x, id.vars = c("ID", "Region", "Day"))
rotation.melt = reshape2::melt(rotation, id.vars = c("kmer_ID"))

# gr = ggplot(data = x.melt, mapping = aes(x = Region, y = value)) + 
#   # geom_point(aes(color = Day)) +
#   geom_jitter(aes(color = Day), width = 0.2, alpha = 0.4) +
#   facet_wrap(vars(variable), ncol = 2) +
#   # stat_summary(fun = mean, fun.min = mean, fun.max = mean,
#   #              geom = "crossbar", width = 0.75) +
#   labs(x = "", y = "")


rotation.melt$value = 100 * rotation.melt$value
rotation.melt$JitterVariable = ave(as.numeric(rotation.melt$variable), rotation.melt$variable, 
                                   FUN = function(x) x + rnorm(length(x), sd = .15))


gr = ggplot(data = x.melt, mapping = aes(x = variable, y = value)) + 
  geom_jitter(aes(fill = Day), shape = 21, colour = "NA", size = 2, stroke = 0.1, alpha = 0.4) +
  geom_segment(data = rotation.melt, 
               mapping = aes(x = JitterVariable, y = 0, 
                             xend = JitterVariable, yend = value),
               arrow = arrow(length = unit(0.2, "cm")),
               size = 1) +
  scale_color_manual(values = c("darkolivegreen", "red")) +
  labs(x = "", y = "")

################### Heatmap - Hierarchical clustering ###################

df = t(kfeatures)

# split rows into clusters
row_dend = hclust(dist(df, method = "euclidean"), method = "ward.D") # row clustering

# split columns
col_dend = hclust(dist(t(df), method="euclidean"), method = "ward.D") # column clustering

clusters_rows = dendextend::cutree(row_dend, k = 4)
clusters_col = dendextend::cutree(col_dend, k = 4)

# write.table(as.data.frame(clusters_rows), file = "clusters_rows.csv", col.names = T, row.names = T, quote = F, sep = ";")
# write.table(as.data.frame(clusters_col), file = "clusters_columns.csv", col.names = T, row.names = T, quote = F, sep = ";")

column_ha3 = HeatmapAnnotation(Region = meta$Region, 
                               Day = anno_lines(meta$`Day Number`),
                               col = list(Region = c("Asia" = "orange", "North America" = "gray", "Oceania" = "steelblue2", 
                                                     "Europe" = "yellow", "South America" = "red", "Africa" = "purple")),
                               na_col = "white")

Heatmap(df, top_annotation = column_ha3, 
        show_row_names = FALSE, show_column_names = FALSE, 
        split = clusters_rows, column_split = clusters_col)

################### Statistical tests ###################

cut.info = as.data.frame(clusters_col)
cut.info$ID = row.names(cut.info)
meta.all = merge(meta,cut.info, by.x = "ID", by.y = "ID")

# analysis for each country

names.meta = colnames(meta.all)
meta.countries = meta.all[!is.na(meta.all$Country),]
names = unique(meta.countries$Country)

analysis_country = matrix(0, nrow = length(names), ncol = 2)


for (i in 1:length(names)){
  meta.test = meta.countries %>% select(Country, clusters_col)
  inter = names[i]
  meta.test = as.data.frame(meta.test)
  meta.test = meta.test[!is.na(meta.test[,1]),]
  colnames(meta.test) = c("factor","clusters_col")
  
  
  test_data = meta.test %>% 
    group_by(clusters_col) %>% 
    summarise(count = n(),
              with = length(which(factor == inter))) %>% 
    mutate(non = count - with)
  
  test_data = as.matrix(test_data[,c(3,4)])
  
  x = fisher.test(test_data, alternative = "two.sided")
  
  
  analysis_country[i,1] = paste(names[i])
  analysis_country[i,2] = paste(x$p.value)
}

# analysis for each Region
names.meta = colnames(meta.all)
meta.countries = meta.all[!is.na(meta.all$Region),]
names = unique(meta.countries$Region)

analysis_table_Region = matrix(0, nrow = (length(names)), ncol = 2)

for (i in 1:length(names)){
  meta.test = meta.countries %>% select(Region, clusters_col)
  inter = names[i]
  meta.test = as.data.frame(meta.test)
  meta.test = meta.test[!is.na(meta.test[,1]),]
  colnames(meta.test) = c("factor","clusters_col")
  
  
  test_data = meta.test %>% 
    group_by(clusters_col) %>% 
    summarise(count = n(),
              with = length(which(factor == inter))) %>% 
    mutate(non = count - with)
  
  test_data = as.matrix(test_data[,c(3,4)])
  
  x = fisher.test(test_data, alternative = "two.sided", workspace = 2e9)
  
  
  analysis_table_Region[i,1] = paste(names[i])
  analysis_table_Region[i,2] = paste(x$p.value)
}


xx = meta.all %>% group_by(Region, clusters_col) %>% summarise(count = n())

ggplot(xx, aes(x = Region, y = count, fill = count)) +
  geom_bar(stat = "identity", width = 0.8) + 
  facet_grid(~clusters_col) +
  coord_flip()

ggplot(meta.all, aes(x = as.character(clusters_col), y = `Day Number`)) + 
  geom_jitter(position = position_jitter(0.1), shape = 21, color = "gray61", fill = "blue", size = 5)+
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median,
               geom = "crossbar", width = 0.5, color="grey31") + 
  geom_text(x = 3, y = 3, label = "p-value < 2.2e-16")


t.test(meta.all %>% filter(clusters_col == 1) %>% select(`Day Number`), meta.all %>% filter(clusters_col == 2) %>% select(`Day Number`))
t.test(meta.all %>% filter(clusters_col == 1) %>% select(`Day Number`), meta.all %>% filter(clusters_col == 3) %>% select(`Day Number`))
t.test(meta.all %>% filter(clusters_col == 3) %>% select(`Day Number`), meta.all %>% filter(clusters_col == 2) %>% select(`Day Number`))


ggplot(meta.all, aes(x = Region, y = `Day Number`, fill = Region)) +
  geom_jitter(position = position_jitter(0.1), shape = 21, color = "gray61", size = 3)+ 
  facet_grid(~clusters_col) 

ggplot(meta.all, aes(x = as.character(clusters_col), y = `Day Number`, fill = Region)) +
  geom_jitter(position = position_jitter(0.1), shape = 21, color = "gray61", size = 3)+ 
  facet_grid(~Region) 



################### Features comparison ###################

seq.clustering = "clusters_rows.csv"

seq.clustering = fread(seq.clustering)
colnames(seq.clustering) = c("kID", "cut")

seq.clustering = seq.clustering[which(seq.clustering$kID %in% colnames(kfeatures)), ]
who = match(seq.clustering$kID, kmers.map$k_ID)

seq.clustering$kmer = kmers.map[who, ]$kmer
seq.clustering$length = str_length(seq.clustering$kmer)

kmer.distances = stringdistmatrix(a = seq.clustering$kmer, method = "lv")
kmer.distances = as.matrix(kmer.distances)

colnames(kmer.distances) = seq.clustering$kID
row.names(kmer.distances) = seq.clustering$kID

kmer.distances = as.dist(kmer.distances)

hc = hclust(kmer.distances, method = "ward.D2")

dhc = as.dendrogram(hc)

# Careful: the function below works for four clusters only
# you need to update colors in order to use it

colLab <- function(n){
  if(is.leaf(n)){

    #I take the current attributes
    a = attributes(n)

    #I deduce the line in the original data, and so the treatment and the specie.
    ligne = match(attributes(n)$label, seq.clustering$kID)

    treatment = seq.clustering[ligne,]$cut;

    if(treatment == 1){ col_treatment = "blue"}
    if(treatment == 2){ col_treatment = "red"}
    if(treatment == 3){ col_treatment = "darkolivegreen"}
    if(treatment == 4){ col_treatment = "black"}

    #Modification of leaf attribute
    attr(n,"nodePar") = c(a$nodePar, list(cex = 1.5,
                                          lab.cex = 1,
                                          pch = 20,
                                          col = col_treatment,
                                          lab.col = col_treatment,
                                          lab.font = 1,
                                          lab.cex = 1))
  }

  return(n)
}


dL = dendrapply(dhc, colLab)
dL = as.dendro(dL)

gr = ggdendrogram(dL, rotate = TRUE, theme_dendro = FALSE)





