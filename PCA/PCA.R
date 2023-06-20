#!/usr/bin/env Rscript

library(mixOmics)
library(tidyverse)
library(optparse)


list_options = list(
  make_option(c("-f", "--file"), type = "character", help = "file_name with the median for all species"),
  make_option(c("-o", "--out"), type = "character", default = getwd(), help = "Directory to export results"),
  make_option(c("-t", "--taxo"), type = "character", help = "Path to taxonomy info"),
  make_option(c("-s", "--subset_level"), type = "character",default = NULL, help = "Level of taxonomy to extract species from (Default = None)"),
  make_option(c("-n", "--subset_name"),type="character", default = NULL, help = "Name of taxonomic group of interest"),
  make_option(c("-g", "--groups"), type = "character", default = "phylum", help ="Taxonomic level to colorize the species")
)

opt_parser = OptionParser(option_list=list_options)
opt = parse_args(opt_parser)


data <- read.csv(opt$file, sep = ",", dec = ".")
taxo <- read.csv(opt$taxo, sep ="\t")
names = taxo[taxo[[opt$subset_level]] == opt$subset_name,]
species = unlist(names$Species)
subset_data = data[data$Species %in% species,]

labels <- unlist(names[[opt$groups]])
metrics <- subset_data[-c(1)]
# metrics <- data[-c(1)]
col.pal <- c("#3e505b", "#000000","#0bb4ff", "#50e991", "#9b19f5", "#00bfa0" ,
             "#b30000","#dc0ab4", "#132c6f","#e6d800" , "#006400", "#D36060",
             "#EA7317", "#2364AA", "#4287f5")

col.pal.sacc <- c("#3e505b", "#000000","#0bb4ff", "#50e991", "#9b19f5", "#00bfa0" ,
                      "#b30000","#dc0ab4", "#132c6f")   
col.pal.3 <- c("#b30000","#e6d800","#006400")

result <- pca(metrics, scale = T)
pdf(file = paste(opt$out,"/ACP_",opt$subset_level,".pdf",sep=""))
result_plot <- pca(metrics, scale = T, ncomp = 5)
plot(result_plot)
plotIndiv(result, style = "ggplot2", group = labels,ellipse.level = 0.8,
          title = paste("PCA ",opt$subset_name,"\nGroups = ",opt$groups,sep=""), legend = TRUE, size.legend = rel(0.8),
          ellipse = T)
plotVar(result, var.names = T, title = paste("PCA_",opt$subset_name,sep=""))
dev.off()

#sPCA
test.keepX = c(seq(1,26,1))
tuning_spca <- tune.spca(metrics, ncomp = 2, nrepeat = 5, folds = floor(nrow(metrics)/3), test.keepX = test.keepX)
res_spca <- spca(metrics,ncomp = 2, center = T, scale = T, keepX = tuning_spca$choice.keepX, max.iter = 500)
pdf(file = paste(opt$out,"/sACP_",opt$subset_level,".pdf",sep=""))
plotIndiv(res_spca, style = "ggplot2", group = labels,
          title = paste("sPCA_",opt$subset_name,"\nGroups = ",opt$groups,sep=""), legend = T, size.legend = rel(0.8), ellipse = T)
plotLoadings(res_spca, method = 'mean', contrib = 'max')
plotVar(res_spca, var.names = T, title = paste("sPCA_",opt$subset_name,sep=""))
dev.off()

