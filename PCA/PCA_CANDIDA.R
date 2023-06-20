
library(mixOmics)
library(tidyverse)
library(ggplot2)
library(optparse)


list_options <- list(
  make_option(c("-f", "--file"), type = "character", help = "file_name with the median for all species"),
  make_option(c("-o", "--out"), type = "character", default = getwd(), help = "Directory to export results"),
  make_option(c("-c", "--candida"), type = "character", help = "Path to candida info")
)

opt_parser <- OptionParser(option_list=list_options)
opt <- parse_args(opt_parser)


data <- read.csv(opt$file, sep = ",", dec = ".")
candida <- read.csv(opt$candida, sep =",")

rownames(data) <- data$Species
candida_metrics <- data[data$Species %in% unlist(candida$species),]

species_names <- candida_metrics$Species
candida_clean <- candida[candida$species %in% unlist(species_names),]
labels = candida_clean$Clade
candida_metrics <- candida_metrics[-c(1)]
result <- pca(candida_metrics, scale = T)
#PCA
pdf(file = paste(opt$out,"/ACP_Candida.pdf",sep=""))
result_plot <- pca(candida_metrics, scale = T, ncomp = 5)
plot(result_plot)
plotIndiv(result, style = "ggplot2", group = labels,ellipse.level = 0.8,
          title = "PCA Candida\nGroups = clades", legend = TRUE, size.legend = rel(0.8), ellipse = T)
plotVar(result, var.names = T, title = "PCA_Candida")
dev.off()

#sPCA
test.keepX = c(seq(1,26,1))
tuning_spca <- tune.spca(candida_metrics, ncomp = 2, nrepeat = 5, folds = floor(nrow(candida_metrics)/3), test.keepX = test.keepX)
res_spca <- spca(candida_metrics,ncomp = 2, center = T, scale = T, keepX = tuning_spca$choice.keepX, max.iter = 500)
pdf(file = paste(opt$out,"/sACP_Candida.pdf",sep=""))
plotIndiv(res_spca, style = "ggplot2", group = labels,
          title = "sPCA Candida\nGroups = clades", legend = T, size.legend = rel(0.8), ellipse = T)
plotLoadings(res_spca, method = 'mean', contrib = 'max')
plotVar(res_spca, var.names = T, title = "PCA_Candida")
dev.off()

