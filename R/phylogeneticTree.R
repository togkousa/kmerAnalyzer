# clear
cat("\014")
rm(list = ls())

# libraries
library(ape)
library(ade4)
library(sets)
library(tseries)
library(seqinr)
library(rlist)
#library(stats)
#library(phyloseq)
#library(phyclust)

# filepath
path <- "phylogeneticTree/clustalo/SARSCoV2.aln"
sequences <- read.alignment(file = path, format="clustal")

# loading metadata
metadata_path <-  "clustalo/final_metafile.txt"
metadata <- read.table(metadata_path, header = TRUE, sep = "\t", skip = 0)

# Reading clusters file
filepath = 'clustalo/clusters_columns.csv'
clustersfile = read.csv(file = filepath)
original_clusters = clustersfile$cluster

# Keeping only the sequences from metadata
to_be_deleted = c()
for (i in 1:length(sequences$nam)){
  id <- sequences$nam[i]
  if (!id %in% metadata$accession){
    to_be_deleted <- c(to_be_deleted, i)
  }
  else{
    sequences$nam[i] <- as.character(metadata$ID[which(id == metadata$accession)])
  }
}
sequences$nb <- length(metadata$accession)
sequences$nam <- sequences$nam[-c(to_be_deleted)]
sequences$seq <- sequences$seq[-c(to_be_deleted)]

####################### SAMPLING ###############################################
# IF you want to take the total data set, set the sample_size variable
# equal to sequences$nb
# sample_size <- sequences$nb
sample_size <- length(sequences$nam)
names <- sequences$nam[c(1:sample_size)]
seqs <- sequences$seq[c(1:sample_size)]
original_clusters <- original_clusters[c(1:sample_size)]
rm(sequences)

for (i in 1:sample_size){
  
  seqs[i] <- strsplit(as.character(seqs[i]), "")
  
}

########### The followning commands are merged in line 73 ######################
#x_dna <- as.DNAbin(seqs, labels = names)

# Calculating Genetic Distances
#D <- dist.dna(x_dna, model = "K80", pairwise.deletion = FALSE)
#temp <- t(as.matrix(D))
#temp <- temp[, ncol(temp):1]
#par(mar = c(1, 5, 5, 1))
#image(x = 1:sample_size, y = 1:sample_size, temp, col = rev(heat.colors(100)),xaxt = "n", yaxt = "n", xlab = "", ylab = "")
#axis(side = 2, at = 1:sample_size, lab = names, las = 2, cex.axis = 0.5)
#axis(side = 3, at = 1:sample_size, lab = names, las = 3, cex.axis = 0.5)
################################################################################

################### TREE CONSTRUCTION ##########################################
# Argument model = "logdet" determines the method that is used to calcilate distance
# between aligned dna sequences. You can check different mehtods in the following link:: 
# https://www.rdocumentation.org/packages/ape/versions/5.4-1/topics/dist.dna

tre <- nj(dist.dna(as.DNAbin(seqs, labels = names), model = "logdet", pairwise.deletion = FALSE))
tre$tip.label = names
write.tree(tre, file = "clustalo/phylotree.nex", append = FALSE,
           digits = 10, tree.names = FALSE)

# load the tree as follows
mytree <- read.tree("phylogeneticTree/clustalo/SARSCoV2_tree.nex")

######### PLOTTING
# Classic plot
plot(mytree, cex = 0.4, tip.color = original_clusters+1)
title("SARS-CoV-2 Genome Sequences Phylogenetic Tree")

# Second plot 
# By changing the type parameter. you can choose the type of the plot you want
# to produce. More info in the following link:
# https://www.rdocumentation.org/packages/ape/versions/5.4-1/topics/plot.phylo

original_clusters = as.character(original_clusters)
original_clusters[which(original_clusters == "1")] = "slateblue1"
original_clusters[which(original_clusters == "2")] = "red"
original_clusters[which(original_clusters == "3")] = "darkgreen"

plot(mytree, type = "phylogram", use.edge.length = FALSE,
     node.pos = 2, show.tip.label = TRUE, show.node.label = FALSE,
     edge.color = "black", edge.width = 1, edge.lty = 1, font = 1,
     cex = 0.4, adj = NULL, srt = 0, no.margin = FALSE,
     root.edge = FALSE, label.offset = 0, underscore = FALSE,
     x.lim = NULL, y.lim = NULL, direction = "rightwards",
     lab4ut = NULL, tip.color = original_clusters, plot = TRUE,
     rotate.tree = 0, open.angle = 0, node.depth = NULL,
     align.tip.label = FALSE)
title("SARS-CoV-2 Genome Sequences Phylogenetic Tree ")