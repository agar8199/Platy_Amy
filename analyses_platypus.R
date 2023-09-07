#installing the developing version of dartR
 library(devtools)
 install_github("green-striped-gecko/dartR@dev_luis")
# loading packages
library(ape)
library(dartR)
library(stringr)
library(dplyr)
# library(tictoc)
library(data.table)

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
source_load <- paste0(getwd(),"/scripts/load_platy.R")
source(source_load)
# the only parameter of the function is the chromosome name
load_list <- load_platy(chrom="X2")
########################################################################
###########################  analyses ##########################
########################################################################
gff <- load_list$gff
genes <- load_list$genes
gene_ox <- load_list$gl_ox
gene_ox_newcc <- load_list$gl_ox_new
gene_dart <- load_list$gl_dart
gene_26 <- load_list$gl_26
gene_dart_26 <- load_list$gl_dart_26

ifn_name <- "IFN"
ifn <- gff[gff$type=="gene",]
ifn <- ifn[grep(ifn_name,ifn$attributes),]
ifn <- cbind(ifn, as.data.frame(str_split(ifn$attributes, pattern = ";", 
                                          simplify = TRUE)))

loci <- gene_26$position
ifn_genes <- NULL
for(i in 1:nrow(ifn)){
  genes <- ifn[i,]
  y  <- list(genes$start,genes$end)
  ifn_genes_tmp <- which(loci %inrange% y)
  ifn_genes <- c(ifn_genes,ifn_genes_tmp)
}

gene_26_ifn <- gl.keep.loc(gene_26,loc.list = locNames(gene_26)[ifn_genes])

dataset_list <- list(dart = gene_dart_26, gl_26 = gene_26_ifn)
names_datset <- names(dataset_list)

for(i in 1:length(dataset_list)){
  db <- dataset_list[[i]]
  if(is.null(db)){
    next()
  }else{
    gl.report.heterozygosity(db,verbose = 0)
    ggsave(paste0("het_ifn_",names_datset[i],".pdf"),  width = 6, height = 6,
           units = "in", dpi="retina",bg = "transparent" )
  }
 

}

