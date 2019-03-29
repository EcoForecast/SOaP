# visualize plot and site level means of NEON functional groups.

rm(list=ls())
library(ggplot2)
library(reshape2)
library(gridExtra)
library(stringr)

#load data and format.----
d <- readRDS("data/ITS_fg_abundances.rds")

# our 5 sites
sites <- c("DSNY", "HARV", "OSBS", "CPER", "STER")

# subset to abundances
abundances <- d$abundances

# convert to relative abundances
# have to add 1 so that there are no zero values (relevant for later analyses)
y <- as.matrix((abundances + 1) / rowSums(abundances + 1)) #use this for stanadard dirichlet.
y <- as.data.frame(y)

# fix rownames and add categorization
rownames(y) <- substr(rownames(y), 1, nchar(rownames(y))-4)
y$site <- substr(rownames(y), 1, 4)
y$plot <- substr(rownames(y), 1, 8)
y$dateID <- substr(rownames(y), nchar(rownames(y))-7, nchar(rownames(y)))
y$dateID <- paste(substr(y$dateID, 1, 4), substr(y$dateID, 5, 6), sep="-")

# remove any rows with weird dateIDs
y <- y[str_detect(y$dateID, "\\d\\d\\d\\d\\-\\d\\d"),]
y <- y[which(y$dateID != "2104-05"), ] # regex didn't pick up this one...

#subset to our 5 sites
y <- y[which(y$site %in% sites),]

# aggregate by site/date, melt into long format
all_data <- list()
for (s in 1:5){
df <- y[which(y$site %in% sites[s]),]
df.no.plot <- df[,-c(9)]
df.by.date <- aggregate(. ~site+dateID, data=df.no.plot, mean, na.rm=TRUE)

df_long <- melt(df.by.date, 
                id.vars = c("site", "dateID"), 
                variable.name = "Functional_Group")
all_data[[s]] <- df_long
}
all_data_to_plot <- do.call(rbind, all_data) # combine data from 5 sites

# make pretty facet plots
p <- ggplot(all_data_to_plot, aes(x = dateID, y = value, fill = Functional_Group)) + 
  geom_bar(stat = "identity") + theme(axis.text.x = element_text(angle = 320, hjust=-.1)) + 
  ggtitle("Site-level averages of functional groups over time") + facet_grid(site ~ .)
p

# save plot
ggsave("Functional_groups_by_site_and_month.png")
