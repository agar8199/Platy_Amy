
load_platy_3 <- function(chrom){
  
  tic("Loading process took")
  
  cat("LOADING CHROMOSOME", chrom,"...\n\n")
  # working directory
  path_dir <- getwd()
  # reading chromosome information
  chrom_platypus <-
    read.csv(paste0(path_dir, "/info/chrom_platypus.csv"))
  chrom_platypus <- chrom_platypus[1:38, ]
  # Chromosome to analyse
  chrom_row <- which(chrom_platypus$Name_b == chrom)
  chrom_ox <- paste0("Martin_chr_", chrom)
  chrom_jenna <- paste0("chr_",chrom)
  chrom <- paste0("chr_", chrom_platypus[chrom_row,"Name"])
  chrom_dir <- paste0(path_dir, "/", chrom, "/")
  
  ########################################################################
  ################### loading GFF file ####################################
  ########################################################################
  # cat("LOADING GFF FILE...\n")
  # 
  # # decompressing GFF
  # system(paste0("gzip -dk ", chrom_dir, chrom, "_mOrnAna1.pri.v4.gff.gz"))
  # gff <- read.gff(paste0(chrom_dir, chrom, "_mOrnAna1.pri.v4.gff"))
  # genes <- gff[gff$type == "gene", ]
  # genes_meta_tmp <- genes$attributes
  # genes_meta <-
  #   as.data.frame(str_split(genes_meta_tmp, pattern = "=|;|-", simplify = TRUE))
  # genes_meta <- genes_meta[, c(3, 5)]
  # colnames(genes_meta) <- c("name", "ID")
  # genes$attributes <- genes_meta$name
  # genes$length <- genes$end - genes$start
  # 
  # cat("  Objects 'genes'and 'gff' created\n\n")
  # print(table(gff$type))
  
  ########################################################################
  ################### loading Jenna samples vcf ##########################
  ########################################################################
  cat("\nLOADING Oxford GENOMES DATASET...\n")
  
  # samples information
  samples_info_ox <-
    read.csv(paste0(path_dir, "/info/", "locations_oxford.csv"))
  # decompressing vcf
  system(paste0("gzip -dk ", chrom_dir, chrom_ox, ".vcf.gz"))
  gl_vcf <- gl.read.vcf(paste0(chrom_dir, chrom_ox, ".vcf"),verbose = 0)
  # filtering
  # INFO codes in vcf file
  # RFGQ_ALL "Empirical quality score (phred scaled) for the call - the geometric
  # mean of all sample RFGQ probabilities"
  # AC "Allele count in genotypes for each ALT allele, in the same order as listed"
  # AN "Total number of alleles in called genotypes"
  # DP "Combined depth across samples"
  # MQ "RMS mapping quality"
  # MQ0 "Number of MAPQ == 0 reads covering this record"
  # NS "Number of samples with data"
  # END "End position on CHROM"
  # MP "Model posterior"
  # names(gl_vcf$other$loc.metrics)
  # gl.report.locmetric(gl_vcf,metric ="RFGQ_ALL" )
  gl_vcf <-
    gl.filter.locmetric(
      gl_vcf,
      metric = "RFGQ_ALL",
      lower = 20,
      upper = max(gl_vcf$other$loc.metrics$RFGQ_ALL),
      verbose = 0
    )
  # assigning information to individuals and populations
  gl_vcf$other$ind.metrics <- samples_info_ox
  indNames(gl_vcf) <- gl_vcf$other$ind.metrics$Sample
  pop(gl_vcf) <- gl_vcf$other$ind.metrics$Region
  
  # ind_keep <- which(gl_vcf$other$ind.metrics$analysis=="YES")
  # gl_vcf <- gl.keep.ind(gl_vcf,ind.list = indNames(gl_vcf)[ind_keep])
  # gl_vcf <- gl.filter.monomorphs(gl_vcf,verbose = 0)
  
  # summary(gl_vcf$other$loc.metrics$RFGQ_ALL)
  
  # gl_vcf <- gl.filter.locmetric(gl_vcf,metric = "RFGQ_ALL",upper = max(gl_vcf$other$loc.metrics$RFGQ_ALL),lower = 20)
  gl_ox <- gl_vcf
  
  # popNames(gl_26)
  # 
  # gl_26 <- gl.keep.loc(gl_26,loc.list = locNames(gl_26)[1:10000])
  # gl_26_2 <- gl.subsample.loci(gl_26,n=10000)
  # 
  # pcoa <- gl.pcoa(gl_26_2)
  # gl.pcoa.plot(pcoa,gl_26_2)
  # 
  # 
  # res <- gl.report.ld.map(gl_26_2,ld_max_pairwise = 2000000,ind.limit = 4)
  # 
  # # saveRDS(res, file="ld_2.rds")
  # 
  # res_2 <- gl.ld.distance(ld_report=res,ld_resolution = 50000)
  
  cat("  Object 'gl_ox' created with", nLoc(gl_ox),"loci\n\n")
  
  #######################################################################
  ################## loading Oxford samples #############################
  #######################################################################
  # cat("LOADING OXFORD DATASET...\n")
  # 
  # # samples information
  # samples_info_ox <-
  #   read.csv(paste0(path_dir, "/info/", "locations_oxford.csv"))
  # # decompressing vcf
  # system(paste0("gzip -dk ", chrom_dir, chrom, "_oxford.vcf.gz"))
  # file_info <- file.info(paste0(chrom_dir, chrom, "_oxford.vcf"))
  # if(file_info$size > 14000){
  # 
  #   gl_vcf_ox <- gl.read.vcf(paste0(chrom_dir, chrom, "_oxford.vcf"),verbose = 0)
  #   # filtering
  #   # names(gl_vcf_ox$other$loc.metrics)
  #   gl_vcf_ox$other$loc.metrics$QUAL <-
  #     as.numeric(gl_vcf_ox$other$loc.metrics$QUAL)
  #   #gl.report.locmetric(gl_vcf_ox,metric ="QUAL")
  # 
  #   gl_vcf_ox <- gl_vcf_ox[order(gl_vcf_ox$ind.names), ]
  # 
  #   # assigning information to individuals and populations
  #   gl_vcf_ox$other$ind.metrics <- samples_info_ox
  #   pop(gl_vcf_ox) <- gl_vcf_ox$other$ind.metrics$Region
  # 
  #   #changing name of chromosome names to the used in genbank
  #   contigs_ox <-
  #     read.csv(paste0(path_dir, "/info/contigs_Oxford.csv"))
  #   t1 <- as.data.frame(as.character(gl_vcf_ox$chromosome))
  #   colnames(t1) <- "Sequence_Name"
  #   t2 <- left_join(t1, contigs_ox, by = 'Sequence_Name')
  #   gl_vcf_ox$chromosome <- as.factor(t2$GenBank_Accn)
  # 
  #   cat("  Remapping Oxford dataset...\n")
  # 
  #   # the system used by the BED format is zero-based for the coordinate start and
  #   # one-based for the coordinate end. Thus, the nucleotide with the coordinate 1 in
  #   # a genome will have a value of 0 in column 2 and a value of 1 in column 3.
  #   # A thousand-base BED interval with the following start and end:
  #   # chr7    0    1000
  #   # writing bed file to convert coordinates between genomes using liftOver
  #   bed_ox <-
  #     cbind(as.character(gl_vcf_ox$chromosome),
  #           (gl_vcf_ox$position - 1),
  #           gl_vcf_ox$position)
  #   bed_file <- paste0(chrom_dir, chrom , "_preLift_ox.bed")
  #   write.table(
  #     bed_ox,
  #     file = bed_file,
  #     col.names = FALSE,
  #     row.names = FALSE,
  #     quote = FALSE,
  #     sep = "\t"
  #   )
  #   liftOver_exe <- paste0(path_dir, "/programs/liftOver")
  #   chain_file_ox <-
  #     paste0(path_dir,
  #            "/info/GCA_002966995.1ToGCF_004115215.2.over.chain.gz")
  #   conversion_file <- paste0(chrom_dir, chrom, "_mOrnAna1_v4.bed")
  #   unMapped_file <- paste0(chrom_dir, chrom, "_unMapped")
  # 
  #   # Granting access to the executable
  #   # and making it executable
  #   # This should be done just once
  #   # system(paste0("chmod +x ",liftOver_exe))
  #   # system(paste0("chmod 755 ",liftOver_exe))
  # 
  #   # running liftOver
  #   system(paste(
  #     liftOver_exe,
  #     bed_file,
  #     chain_file_ox,
  #     conversion_file,
  #     unMapped_file
  #   ))
  # 
  #   #loading mapped SNPs
  #   map <- read.table(conversion_file)
  # 
  #   # loading unmapped SNPs
  #   if(length(readLines(unMapped_file))>0){
  #     unmap <- read.table(unMapped_file)
  #   }else{
  #     unmap <- NULL
  #   }
  # 
  #   loc_unmap <- unlist(lapply(unmap$V2, function(x) {
  #     which(gl_vcf_ox$position == x)
  #   }))
  # 
  #   gl_vcf_ox <-
  #     gl.drop.loc(gl_vcf_ox, loc.list = locNames(gl_vcf_ox)[loc_unmap],verbose = 0)
  # 
  #   gl_vcf_ox$chromosome <- as.factor(map$V1)
  #   gl_vcf_ox$position <- map$V2
  # 
  #   loc_chr <-
  #     which(as.character(gl_vcf_ox$chromosome) != chrom_platypus[chrom_row, "RefSeq_mOrnAna_v4"])
  # 
  #   gl_vcf_ox <-
  #     gl.drop.loc(gl_vcf_ox, loc.list = locNames(gl_vcf_ox)[loc_chr],verbose = 0)
  # 
  #   gl_ox <- gl_vcf_ox
  #   
  #   ind_keep_ox <- which(gl_ox$other$ind.metrics$analysis=="YES")
  #   gl_ox <- gl.keep.ind(gl_ox,ind.list = indNames(gl_ox)[ind_keep_ox])
  #   gl_ox <- gl.filter.monomorphs(gl_ox,verbose = 0)
  # 
  #   cat("  Object 'gl_ox' created with", nLoc(gl_ox),"loci\n")
  # 
  # }else{
  #   cat(" No loci in the Oxford dataset for this chromosome\n\n")
  #   gl_ox <- NULL
  # }
  # 
  
  ########################################################################
  ################### loading Jenna samples vcf ##########################
  ########################################################################
  cat("\nLOADING 26 GENOMES DATASET...\n")
  
  # samples information
  samples_info <-
    read.csv(paste0(path_dir, "/info/", "samples_jenna.csv"))
  # decompressing vcf
  system(paste0("gzip -dk ", chrom_dir, chrom_jenna, ".vcf.gz"))
  gl_jenna <- gl.read.vcf(paste0(chrom_dir, chrom_jenna, ".vcf"),verbose = 0)
  # filtering
  # INFO codes in vcf file
  # RFGQ_ALL "Empirical quality score (phred scaled) for the call - the geometric
  # mean of all sample RFGQ probabilities"
  # AC "Allele count in genotypes for each ALT allele, in the same order as listed"
  # AN "Total number of alleles in called genotypes"
  # DP "Combined depth across samples"
  # MQ "RMS mapping quality"
  # MQ0 "Number of MAPQ == 0 reads covering this record"
  # NS "Number of samples with data"
  # END "End position on CHROM"
  # MP "Model posterior"
  # names(gl_vcf$other$loc.metrics)
  # gl_vcf <- gl.filter.monomorphs(gl_vcf,verbose = 0)
  #gl.report.locmetric(gl_vcf,metric ="RFGQ_ALL" )
  gl_jenna <-
    gl.filter.locmetric(
      gl_jenna,
      metric = "RFGQ_ALL",
      lower = 20,
      upper = max(gl_jenna$other$loc.metrics$RFGQ_ALL),
      verbose = 0
    )
  # assigning information to individuals and populations
  gl_jenna$other$ind.metrics <- samples_info
  indNames(gl_jenna) <-
    c(
      "S11",
      "S13",
      "S14",
      "S17",
      "S20",
      "S22",
      "S23",
      "S27",
      "S3",
      "S31",
      "S35",
      "S36",
      "S5",
      "S8",
      "S9",
      "T1",
      "T16",
      "T18",
      "T21",
      "T23",
      "T27",
      "T29",
      "T36",
      "T38",
      "T5",
      "T8"
    )
  pop(gl_jenna) <- gl_jenna$other$ind.metrics$Group
  
  gl_26 <- gl_jenna
  
  cat("  Object 'gl_26' created with", nLoc(gl_26),"loci\n\n")
  
  ########################################################################
  ################### loading DArT samples ###############################
  ########################################################################
  # cat("LOADING DArT DATASET...\n")
  # 
  # gl_vcf_dart <- gl.load(paste0(chrom_dir, chrom, "_dart.rds"),verbose = 0)
  # gl_vcf_dart$chromosome <-
  #   as.factor(gl_vcf_dart$other$loc.metrics$Chrom_Platypus_Chrom_NCBIv1)
  # gl_vcf_dart$position <-
  #   gl_vcf_dart$other$loc.metrics$ChromPos_Platypus_Chrom_NCBIv1
  # 
  # # decompressing fasta
  # ref_gen_file_zip <- paste0(chrom_dir, chrom, "_mOrnAna1.pri.v4.fa.gz")
  # system(paste0("gzip -dk ", ref_gen_file_zip))
  # ref_gen_file_unzip <- paste0(chrom_dir, chrom, "_mOrnAna1.pri.v4.fa")
  # cat("  Remapping DArT dataset...\n")
  
  # remapping using BLAST
  # Installing BLAST
  #'  You can download the BLAST installs from:
  #'  https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/
  #'  It is important to install BLAST in a path that does not contain spaces for
  #'  this function to work.
  
  # gl_vcf_dart <- gl.blast(gl_vcf_dart, ref_genome = ref_gen_file_unzip,verbose = 0)
  # pos <-  gl_vcf_dart$other$loc.metrics$SnpPosition + 
  #   gl_vcf_dart$other$loc.metrics$sstart
  # 
  # gl_vcf_dart$position <- as.integer(pos)
  # 
  # gl_dart <- gl_vcf_dart
  # 
  # cat("  Object 'gl_dart' created with", nLoc(gl_dart),"loci\n\n")
  
  ########################################################################
  ################### separating samples from Jenna's rivers #############
  ########################################################################
  # cat("SUBSETTING gl_dart TO COINCIDE WITH gl_26 ...\n")
  # 
  # gl_vcf_dart_jenna <-
  #   gl.keep.pop(gl_vcf_dart,
  #               pop.list = c("TENTERFIELD",
  #                            "SEVERN_BELOW",
  #                            "SEVERN_ABOVE"),verbose = 0)
  # gl_vcf_dart_jenna <-
  #   gl.filter.monomorphs(gl_vcf_dart_jenna, verbose = 0)
  # 
  # #changing the names so they coincide with Jenna's names
  # indNames(gl_vcf_dart_jenna) <-
  #   c(
  #     "T27",
  #     "T35",
  #     "S4",
  #     "S12",
  #     "S20",
  #     "S28",
  #     "S36",
  #     "T11",
  #     "T19",
  #     "T28",
  #     "T36",
  #     "S5",
  #     "S13",
  #     "S21",
  #     "S29",
  #     "S37",
  #     "T4",
  #     "T12",
  #     "T20",
  #     "T29",
  #     "S6",
  #     "S14",
  #     "S22",
  #     "S30",
  #     "S38",
  #     "T5",
  #     "T13",
  #     "T21",
  #     "T30",
  #     "T38",
  #     "S7",
  #     "S15",
  #     "S23",
  #     "S31",
  #     "S39",
  #     "T6",
  #     "T14",
  #     "T22",
  #     "T31",
  #     "T39",
  #     "S8",
  #     "S16",
  #     "S24",
  #     "S32",
  #     "S40",
  #     "T7",
  #     "T15",
  #     "T23",
  #     "T32",
  #     "T40",
  #     "S9",
  #     "S17",
  #     "S25",
  #     "S33",
  #     "S41",
  #     "T8",
  #     "T16",
  #     "T24",
  #     "S2",
  #     "T33",
  #     "T41",
  #     "S10",
  #     "S18",
  #     "S26",
  #     "S34",
  #     "T1",
  #     "T9",
  #     "T17",
  #     "T25",
  #     "S3",
  #     "T34",
  #     "S11",
  #     "S19",
  #     "S27",
  #     "S35",
  #     "T2",
  #     "T10",
  #     "T18",
  #     "T26"
  #   )
  # 
  # gl_dart_26 <- gl_vcf_dart_jenna
  # 
  # cat("  Object 'gl_dart_26' created with", nLoc(gl_dart_26),"loci\n\n")
  # 
  ########################################################################
  # testing percentage of correct genotypes between DArT and new dataset##
  ########################################################################
  # cat("Testing percentage of correct genotypes between DArT and 26 genomes\n")
  gl_26$loc.names <- as.character(gl_26$position)
  gl_ox$loc.names <- as.character(gl_ox$position)
  res_1 <- gl.report.pa(gl_26,loc_names = TRUE)
  res_2 <- gl.report.pa(gl_26,loc_names = TRUE,  method ="one2rest")  
  res_3 <- gl.report.pa(gl_ox,loc_names = TRUE)
  res_4 <- gl.report.pa(gl_ox,loc_names = TRUE,  method ="one2rest")
  gl.save(gl_26,file = paste0("gl_26_",chrom,".rds"))
  gl.save(gl_ox,file = paste0("gl_ox_",chrom,".rds"))
  saveRDS(res_1,file = paste0("pa_pairwise_26_",chrom,".rds"))
  saveRDS(res_2,file = paste0("pa_one2rest_26_",chrom,".rds"))
  saveRDS(res_3,file = paste0("pa_pairwise_ox_",chrom,".rds"))
  saveRDS(res_4,file = paste0("pa_one2rest_ox_",chrom,".rds"))
  
  toc()
  #   
}

# genotypes_2_datasets <- function(x,y){
#   
#   cat("Testing percentage of correct genotypes between two datasets\n")
#   
#   x <- gl.keep.ind(x, ind.list = indNames(y),verbose = 0)
#   x <- gl.filter.callrate(x, threshold = 1,verbose=0)
#   x <- gl.filter.monomorphs(x,verbose = 0)
#   
#   y <- gl.filter.callrate(y, threshold = 1,verbose = 0)
#   
#   t2 <- which(x$loc.names %in% y$loc.names)
#   t3 <- which(y$loc.names %in% x$loc.names)
#   
#   x <- gl.keep.loc(x, loc.list = locNames(x)[t2],verbose = 0)
#   x <- x[, order(x$loc.names)]
#   x <- x[order(indNames(x)), ]
#   
#   y <- gl.keep.loc(y, loc.list = locNames(y)[t3],verbose = 0)
#   y <- y[, order(y$loc.names)]
#   y <- y[order(indNames(y)), ]
#   
#   s1 <- as.matrix(x)
#   s2 <- as.matrix(y)
#   
#   s1[s1 == 2] <- 0
#   s2[s2 == 2] <- 0
#   
#   s2 <- s2[,(1:ncol(s2))-1]
#   
#   # percentage of incorrect called genotypes between DArT and new dataset
#   genotypes <- nInd(y) * nLoc(y)
#   correct_genotypes <- sum(colSums(s1 == s2))
#   
#   incorrect_percentage <- (correct_genotypes / genotypes) * 100
#   
#   incorrect_percentage <-  round(incorrect_percentage,2)
#   cat(incorrect_percentage,"percentage of", genotypes,"genotypes were the same in dataset 1 and dataset 2\n\n")
# }
# 
# gl_26_2 <- gl.subsample.loci(x_martin,n=10000)
# pcoa <- gl.pcoa(gl_26_2)
# gl.pcoa.plot(pcoa,gl_26_2)
# 
# 
# 
# 
# 