# loading packages
library(ape)
library(dartR)
library(stringr)
library(dplyr)
library(tictoc)
library(data.table)

# args <- commandArgs(trailingOnly = TRUE)
# 
# .libPaths("/home/549/jm9807/.R")
# 
# setwd("/g/data/te48")

# installing the developing version of dartR
# gl.install.vanilla.dartR(flavour = "dev")
########################################################################
###########################  loading datasets ##########################
########################################################################
# Objects loaded with function load_platy are: 
# genes 
# gff 
# gl_26 
# gl_ox  
# gl_dart  
# gl_dart_26  
source_load <- paste0(getwd(),"/scripts/load_platy_martin.R")
source(source_load)

chrom_platypus <-
  read.csv(paste0(getwd(), "/info/chrom_platypus.csv"))

for(i in chrom_platypus$Name_b){
# the only parameter of the function is the chromosome name
  gene_ox <- load_platy_3(i)
########################################################################
###########################  analyses ##########################
########################################################################
# gff <- load_list$gff
# genes <- load_list$genes
# gene_ox <- load_list$gl_ox
# gene_dart <- load_list$gl_dart
# gene_26 <- load_list$gl_26
# gene_dart_26 <- load_list$gl_dart_26

# gene_26 <- gl.drop.pop(gene_26,pop.list = "TENTERFIELD" )
  # res <- gl.outflank(gene_ox,RightTrimFraction=0.25)
  # res <- res$outflank
  # saveRDS(res,paste0("/Users/s441489/final_platy/outflank/ox_outflank_",i,".rda"))
  #
  # pcoa <- gl.pcoa(gene_ox)
  # gl.pcoa.plot(pcoa,gene_ox)
  #  ggsave(paste0("ox_plot_",i,".pdf"),  width = 6, height = 6, units = "in", dpi="retina",bg = "transparent" )
  #
  #
  # saveRDS(pcoa,paste0("/Users/s441489/final_platy/outflank/ox_pca_",i,".rda"))

# 
# MHC_name <- "MHC"
# mhc <- gff[gff$type=="gene",]
# mhc <- mhc[grep(MHC_name,mhc$attributes),]
# mhc <- cbind(mhc, as.data.frame(str_split(mhc$attributes, pattern = ";", simplify = TRUE)))
# 
# loci <- gene_26$position
# mhc_genes <- NULL
# for(i in 1:nrow(mhc)){
#   genes <- mhc[i,]
#   y  <- list(genes$start,genes$end)
#   mhc_genes_tmp <- which(loci %inrange% y)
#   mhc_genes <- c(mhc_genes,mhc_genes_tmp)
# }
# 
# gene_26_mhc <- gl.keep.loc(gene_26,loc.list = locNames(gene_26)[mhc_genes])
# 
# dataset_list <- list(dart=gene_dart_26,jenna=gene_26_mhc)
# names_datset <- names(dataset_list) 
# 
# for(i in 1:length(dataset_list)){
#   db <- dataset_list[[i]]
#   gl.report.heterozygosity(db,verbose = 0)
#   ggsave(paste0("het_mhc_",names_datset[i],".pdf"),  width = 6, height = 6, units = "in", dpi="retina",bg = "transparent" )
# }
# 
# 
# lapply(dataset_list,function(x){
#   pcoa <- gl.pcoa(x,verbose = 0,plot.out = FALSE)
#   gl.pcoa.plot(pcoa,x,verbose = 0)
# })
# 
# 
# res <- gl.outflank(gene_26)





}















