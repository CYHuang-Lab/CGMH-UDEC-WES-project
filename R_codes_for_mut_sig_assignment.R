### R codes for mutational signature assignment (SBS96 and ID)
rm(list = ls())

# Install dependency packages and load in data
library(remotes)
library(dplyr)
library(data.table)

# to install mSigAct R package, please use: 
# remotes::install_github(repo = "steverozen/mSigAct", ref = "v2.3.2-branch")
library(mSigAct) 

# to install ICAMS R package, please use: 
# remotes::install_github(repo = "steverozen/ICAMS", ref = "v3.0.6-branch")
library(ICAMS) 

# to install cosmicsig R package, please use: 
# remotes::install_github(repo = "Rozen-Lab/cosmicsig", ref = "v1.0.7-branch")
library(cosmicsig) 

download.file("https://github.com/CYHuang-Lab/CGMH-UDEC-WES-project/blob/main/CGMH-UDEC-catSBS96.csv")
download.file("https://github.com/CYHuang-Lab/CGMH-UDEC-WES-project/blob/main/CGMH-UDEC-catSBS96.csv")

catSBS96 <- ReadCatalog("CGMH-UDEC-catSBS96.csv") %>% 
  as.catalog(ref.genome = "GRCh38", region = "genome", catalog.type = "counts")

catID <- ReadCatalog("CGMH-UDEC-catID.csv") %>% 
  as.catalog(ref.genome = "GRCh38", region = "genome", catalog.type = "counts")

sigs_sbs96_genome <- cosmicsig::COSMIC_v3.2$signature$GRCh38$SBS96

sigs_sbs96_exome <- ICAMS::TransformCatalog(
  catalog = sigs_sbs96_genome, target.region = "exome", target.ref.genome = "GRCh37")

sigs_prop_sbs96 <- mSigAct::ExposureProportions( 
  mutation.type = "SBS96", cancer.type = "Uterus-AdenoCA") #PCAWG7::CancerTypes()

sbs96.base <- paste0("SBS",c(1,2,5,13,28,40))
sbs96.pole <- paste0("SBS",c("10a","10b"))
sbs96.dmmr <- paste0("SBS",c(6,14,15,26,44))
sbs96.pold1 <- "SBS20"
sbs96.full <- c(sbs96.base, sbs96.pole, sbs96.dmmr, sbs96.pold1)

sigs_sbs96_exome <- sigs_sbs96_exome[, sbs96.full, drop = FALSE]
my_opts <- mSigAct::DefaultManyOpts(likelihood.dist = "multinom")

## set ID/SBS threshold to decide SBS universe for the reconstruction
sig.table <- data.frame(
  ID = colnames(catSBS96), 
  SBS = colSums(catSBS96), 
  CCTA = catSBS96["CCTA",],
  TCTA = catSBS96["TCTA",],
  Homopolymer_T = colSums(catID[grep("DEL:T",rownames(catID)),])
)

sigs_to_use <- list()
for(i in 1:nrow(sig.table)){
  sigs_to_use[[i]] <- sbs96.base
  if(sig.table$CCTA[i] > 3000){
    sigs_to_use[[i]] <- c(sigs_to_use[[i]], sbs96.pold1)
  } 
  
  if(sig.table$TCTA[i] > 1000){
    sigs_to_use[[i]] <- c(sigs_to_use[[i]], sbs96.pole)
  } 
  
  if(sig.table$Homopolymer_T[i] > 100){
    sigs_to_use[[i]] <- c(sigs_to_use[[i]], sbs96.dmmr)
  } 
}

## mutational signature assignment for SBS catalogs
retval <- lapply(1:nrow(sig.table), function(x){
  mSigAct::SparseAssignActivity(
    spectra = catSBS96[,x,drop=FALSE],
    sigs = sigs_sbs96_exome[,sigs_to_use[[x]]],
    output.dir = "./mutsig/SBS_sparse/",
    max.level = length(sigs_to_use[[x]]) - 1,
    p.thresh = 0.05 / length(sigs_to_use[[x]]),
    m.opts = my_opts,
    num.parallel.samples = 1,
    mc.cores.per.sample = 1,
    seed = 3838,
    max.subsets = 1e15,
    drop.low.mut.samples = FALSE
  )
})

saveRDS(retval, "./mutsig/SBS_sparse/sparse.out.SBS96.rds")
