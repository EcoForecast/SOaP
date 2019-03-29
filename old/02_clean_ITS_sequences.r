# clean up ITS sequences, match them to covariate sample IDs.
# adapted from Colin Averill's script at: NEFI_microbe/ITS/data_construction/NEON_ITS/4._sequence_processing_raw_fastq/5._map_tax_fun_cleanup.r

# using sequences from 2013-2014 that Colin Averill received from Lee and had already processed. 
# this does not include the 2016-2017 validation data, which came directly from Lee; 
# we are waiting for a response from her concerning the state of processing
rm(list=ls())
library(data.table)
library(RCurl)
source("functions/common_group_quantification.r")

# these files are on the SCC, with universal read/write/execute permissions
# read in sequence variant (SV) table - contains the number of SVs per sample
sv <- readRDS("/projectnb/talbot-lab-data/NEFI_data/ITS/scc_gen/NEON_fastq_SV.table.rds")
# read in functional group table - contains the taxonomic/functional group assignment for each SV
fun <- readRDS("/projectnb/talbot-lab-data/NEFI_data/ITS/scc_gen/NEON_fastq_fun.rds")

#remove samples with less than 1000 reads from sv.----
sv <- as.data.frame(sv)
sv$seq.depth <- rowSums(sv)
sv <- sv[sv$seq.depth > 1000,]
sv$seq.depth <- NULL

#kill SVs that no longer have any sequences or are not fungi.----
to_remove <- colnames(sv[,colSums(sv) == 0])
fungi.check <- rownames(fun[fun$kingdom != 'Fungi' | is.na(fun$kingdom),])
to_remove <- c(to_remove, fungi.check)
sv <-  sv[,!(colnames(sv) %in% to_remove) ]
fun <- fun[!(rownames(fun) %in% to_remove),]

#rarefy SV table to 1k reads/sample.----
set.seed(5)
sv.1k <- vegan::rrarefy(sv, 1000)
sv.1k <- sv.1k[,colSums(sv.1k) > 0]

#Build phylo-functional group taxonomy tables.----
fg <- data.table(fun)
fg[grep('Arbuscular'     , guild), fg := 'Arbuscular'     ]
#fg[grep('Pathogen'       , guild), fg := 'Pathogen'       ]
fg[grep('Animal Pathogen', guild), fg := 'Animal_Pathogen']
fg[grep('Plant Pathogen' , guild), fg := 'Plant_Pathogen' ]
fg[grep('Saprotroph'     , guild), fg := 'Saprotroph'     ]
fg[grep('Wood Saprotroph', guild), fg := 'Wood_Saprotroph']
fg[grep('Ectomycorrhizal', guild), fg := 'Ectomycorrhizal']
fg <- fg[,.(kingdom, phylum, class, order, family, genus, species, fg)]
fg <- as.data.frame(fg)
rownames(fg) <- rownames(fun)
fg.1k <- fg[rownames(fg) %in% colnames(sv.1k),]

# calculate absolute and relative abundances
groups <- unique(fg.1k$fg)
groups <- groups[!is.na(groups)]
quant <- common_group_quantification(sv = sv.1k, tax = fg.1k, groups = groups, tax_level = "fg")

# save output
saveRDS(sv.1k, "data/ITS_sv.rds")
saveRDS(fg.1k, "data/ITS_fg_assignments.rds")
saveRDS(quant, "data/ITS_fg_abundances.rds")
