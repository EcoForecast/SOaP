#load data.
dat <- readRDS("data/NEON_soil_tree_merge.rds")
unique(dat$siteID)
unique(dat[,c("plotID", "dateID")])
p.codes <- read.csv("data/NEON_tree/NEON_DP1.10098.plantStatus_decode.csv")
myc.spp <- readRDS("data/NEON_tree/myc_species_assignments.rds")
myc.gen <- read.csv("data/NEON_tree/tedersoo_2017_em_genera.csv")
myc.spp$genus_spp <- myc.spp$Species
poa.gen <- readRDS("data/NEON_tree/poaceae_genera_wikipedia.rds")

am_genera <- c('Thuja','Fraxinus','Nyssa','Celtis','Cornus','Diospyros','Ilex','Lonicera','Magnolia','Viburnum',poa.gen)
erm_genera <- c('Rhododendron','Vaccinium')

#Everything is ID'd at least to genus. Species has lots of qualifiers. Lets clear these up.
dat$species <- sub("^(\\S*\\s+\\S+).*", "\\1", dat$scientificName)
dat$genus   <- sub(" .*", ""                 , dat$scientificName)

#assign mycorrhizal status
dat <- merge(dat,myc.spp[,c('Species','MYCO_ASSO')], by.x = 'species', by.y = 'Species', all.x=T)
dat[dat$genus %in% myc.gen$genus,]$MYCO_ASSO <- 'ECM'
dat[dat$genus %in% am_genera,    ]$MYCO_ASSO <- 'AM'
dat[dat$genus %in% erm_genera,   ]$MYCO_ASSO <- 'ERM'

#subset to trees that have a stem diameter measurement.
dat <- dat[grep('tree',dat$growthForm),]
dat <- dat[!(is.na(dat$stemDiameter)),]

#how much of the basal area is assigned a mycorrhizal association? 96%. We good.
dat$basal <- pi * (dat$stemDiameter/2)^2
metric <- round((1 - (sum(dat[is.na(dat$MYCO_ASSO),]$basal) / sum(dat$basal)))*100, 2)
cat(paste0(metric,'% of trees assigned a mycorrhizal association.'))

#We now need to account for dead trees, insect damaged trees.
#deal with legacy codes in data using key.
for(i in 1:nrow(dat)){
  if(dat$plantStatus[i] %in% p.codes$lovElementName){
    dat$plantStatus[i] <- as.character(p.codes[p.codes$lovElementName == dat$plantStatus[i],]$lovElementCode)
  }
}
#2. aggregate to site scale - this is where uncertainty comes in.----
output.list <- list()
n.reps <- 1000
for(i in 1:n.reps){
  sim.dat <- dat
  sim.dat$stemDiameter <- rnorm(nrow(sim.dat),log(sim.dat$stemDiameter), 0.0316)
  sim.dat$stemDiameter <- exp(sim.dat$stemDiameter)
  sim.dat$basal <- pi * (sim.dat$stemDiameter/2)^2
  
  #Assign live, live_ECM and dead basal area.
  sim.dat$basal_live <- ifelse(grepl('Live', sim.dat$plantStatus) == T, sim.dat$basal, 0)
  sim.dat$basal_dead <- ifelse(grepl('Dead', sim.dat$plantStatus) == T | grepl('dead', sim.dat$plantStatus) == T, sim.dat$basal, 0)
  sim.dat$basal_ECM  <- ifelse(sim.dat$MYCO_ASSO == 'ECM', sim.dat$basal_live, 0)
  
  #aggregate.
  site.level            <- aggregate(basal_live ~ siteID, data = sim.dat, FUN = sum, na.rm=T, na.action = na.pass)
  site.level$basal_ECM  <- aggregate(basal_ECM  ~ siteID, data = sim.dat, FUN = sum, na.rm=T, na.action = na.pass)[,2]
  site.level$basal_dead <- aggregate(basal_dead ~ siteID, data = sim.dat, FUN = sum, na.rm=T, na.action = na.pass)[,2]
  site.level <- site.level[!(site.level$basal_live == 0),]
  site.level$relEM <- site.level$basal_ECM / site.level$basal_live
  site.level$live_fraction <- site.level$basal_live / (site.level$basal_live + site.level$basal_dead)
  
  #save to list
  output.list[[i]] <- site.level
}
for(i in 1:length(output.list)){
  if(i == 1){
    basal_live <- output.list[[i]]$basal_live
    basal_ECM  <- output.list[[i]]$basal_ECM
    basal_dead <- output.list[[i]]$basal_dead
    relEM      <- output.list[[i]]$relEM
  }
  if(i > 1){
    basal_live <- cbind(basal_live, output.list[[i]]$basal_live)
    basal_ECM  <- cbind(basal_ECM , output.list[[i]]$basal_ECM )
    basal_dead <- cbind(basal_dead, output.list[[i]]$basal_dead)
    relEM      <- cbind(relEM     , output.list[[i]]$relEM     )
  }
}
desc.list <- list(basal_live,basal_ECM,basal_dead,relEM)
site.mean <- list()
site.sd   <- list()
for(i in 1:length(desc.list)){
  site.mean[[i]] <- rowMeans(desc.list[[i]])
  site.sd  [[i]] <- matrixStats::rowSds(desc.list[[i]])
}
site.mean <- data.frame(site.mean)
site.sd   <- data.frame(site.sd  )
colnames(site.mean) <- c('basal_live','basal_ECM','basal_dead','relEM')
colnames(site.sd  ) <- c('basal_live_sd','basal_ECM_sd','basal_dead_sd','relEM_sd')
plotID <- output.list[[1]]$siteID
siteID <- unique(substring(siteID,1,4))
site.level <- data.frame(siteID, site.mean,site.sd)
site.level


#3. save plot-level EM relative abundance.----
saveRDS(site.level, "data/site_level_trees.rds")