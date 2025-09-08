## Conduct sensitivity analysis to affect if results are robust to salmon
# size

# when making maps, sub >= for copied symbol to get nice >= sign

## 0. Script description ----

# Script uses data (summer-salish-Chinook-before-clustering.RDS and winter-
# salish-Chinook-before-clustering.RDS) output from summer-grid-analysis.R and
# winter-grid-analysis.R
# These are the data frames that were used in cluster analyses (e.g., stomachs
# that were empty or contain only unidentified prey were removed, minor prey
# taxa were removed).

# 1. Visualize size of Chinook salmon by PFMA and season.

# 2. Visualize diet of Chinook salmon by 5 cm length bin, separately for summer
# and winter.

# 3. Conducts cluster analysis with length subsets (>= 55 cm and >= 62 cm),
# separately for summer and winter, to assess if regionalization results are
# robust to regionally-varying size regulations.

# load data and scripts
library(tidyverse)
library(sf)
library(reshape2)
library(vegan)
library(labdsv)
library(ggpubr)

summer.salish.ch <- readRDS("data/out/summer-diet-for-size-sensitivity.RDS")
winter.salish.ch <- readRDS("data/out/winter-diet-for-size-sensitivity.RDS")

nrow(summer.salish.ch) # 968
nrow(winter.salish.ch) # 723

bc.coast <- read_sf("data/BC-coast-shapefile/bc-coast.shp")
grid <- readRDS("data/out/grid_17_5.RDS")
attr(grid, "grid") # 17_5

## 1. Size by PFMA ----

both.seasons <- rbind(select(summer.salish.ch, Season, PFMA, `Final Length (cm)`),
                      select(winter.salish.ch, Season, PFMA, `Final Length (cm)`))
both.seasons.length <- subset(both.seasons, !is.na(`Final Length (cm)`))
nrow(both.seasons.length)

# by season and PFMA
length.summary <- both.seasons.length %>% group_by(Season, PFMA) %>% 
  summarize(median_length = median(`Final Length (cm)`))

# by season
both.seasons.length %>% group_by(Season) %>% 
  summarize(median_length = median(`Final Length (cm)`),
            min_length = min(`Final Length (cm)`),
            max_length = max(`Final Length (cm)`))

both.seasons.length.N <- both.seasons.length %>% group_by(Season, PFMA) %>% tally()

# size by PFMA
ggplot(data = both.seasons.length,
       aes(x = factor(PFMA), y = `Final Length (cm)`, colour = Season)) +
  geom_jitter(width = 0.2, height = 0, alpha = 0.3, size = 1.5) +
  geom_boxplot(outliers = FALSE, fill = "transparent", size = 0.8) +
  geom_text(data = both.seasons.length.N,
            aes(x = factor(PFMA), y = 100, label = n), colour = "black") +
  scale_colour_manual(values = c("firebrick3", "dodgerblue3")) +
  xlab("Pacific Fishery Management Area") + ylab("Fork length (cm)") + 
  facet_wrap(~Season, nrow = 2) + 
  theme_bw() +
  theme(axis.title = element_text(size = 12),
        axis.text = element_text(size = 11),
        strip.background = element_blank(),
        strip.text.x = element_text(size = 12),
        panel.border = element_rect(colour = "black", linewidth = 0.9),
        legend.position = "none")
ggsave("figures/supplement/PFMA-size.PNG", 
       width = 20, height = 17, units = "cm", dpi = 600)

## 2. Calculate diet by size bins in summer and winter ----

## 2A. Summer ----
range(summer.salish.ch$`Final Length (cm)`, na.rm = TRUE)
summer.salish.ch$length_bin <- cut(summer.salish.ch$`Final Length (cm)`,
  breaks = seq(45, 100, 5), right = FALSE)
summer.salish.ch$length_bin <- as.character(summer.salish.ch$length_bin)
summer.salish.ch <- relocate(summer.salish.ch, length_bin, .after = "Final Length (cm)")

summer.with.length <- subset(summer.salish.ch, !is.na(`Final Length (cm)`))
nrow(summer.with.length) # 863

summer.length.bins <- sort(unique(summer.with.length$length_bin))
summer.length.bins
names(summer.with.length)
ncol(summer.salish.ch[20:35]) # 16
summer.length.diet <- data.frame(matrix(NA, nrow = 16, ncol = 0))

for (i in 1:length(summer.length.bins)) {
  bin <- subset(summer.with.length, length_bin == summer.length.bins[i])
  mean.diet <- colMeans(bin[20:35])
  summer.length.diet <- cbind(summer.length.diet, mean.diet)
  colnames(summer.length.diet)[i] <- summer.length.bins[i]
}

colSums(summer.length.diet) # expected since minor prey removed already

# adjust for minor prey removed
summer.length.rare.prey.adjust <- data.frame(t(data.frame(1 - colSums(summer.length.diet))))
rownames(summer.length.rare.prey.adjust) <- "Rare prey"
names(summer.length.rare.prey.adjust) <- names(summer.length.diet)
summer.length.rare.prey.adjust

summer.length.diet <- rbind(summer.length.diet, summer.length.rare.prey.adjust)
summer.length.diet <- summer.length.diet * 100
colSums(summer.length.diet)
summer.length.diet$Group <- rownames(summer.length.diet)

summer.length.diet.long <- melt(summer.length.diet)
names(summer.length.diet.long) <- c("Group", "Length bin (cm)", "Percentage (%)")
summer.length.diet.long$Season <- "Summer"

# 2B. Winter ----
range(winter.salish.ch$`Final Length (cm)`, na.rm = TRUE)
winter.salish.ch$length_bin <- cut(winter.salish.ch$`Final Length (cm)`,
                                   breaks = seq(45, 100, 5), right = FALSE)
winter.salish.ch$length_bin <- as.character(winter.salish.ch$length_bin)
winter.salish.ch <- relocate(winter.salish.ch, length_bin, .after = "Final Length (cm)")

winter.with.length <- subset(winter.salish.ch, !is.na(`Final Length (cm)`))
nrow(winter.with.length) # 695

winter.length.bins <- sort(unique(winter.with.length$length_bin))
winter.length.bins
names(winter.with.length)
winter.length.diet <- data.frame(matrix(NA, nrow = 18, ncol = 0))

for (i in 1:length(winter.length.bins)) {
  bin <- subset(winter.with.length, length_bin == winter.length.bins[i])
  mean.diet <- colMeans(bin[20:37])
  winter.length.diet <- cbind(winter.length.diet, mean.diet)
  colnames(winter.length.diet)[i] <- winter.length.bins[i]
}

colSums(winter.length.diet) # expected because minor prey removed

# account for rare prey
winter.length.rare.prey.adjust <- data.frame(t(data.frame(1 - colSums(winter.length.diet))))
rownames(winter.length.rare.prey.adjust) <- "Rare prey"
names(winter.length.rare.prey.adjust) <- names(winter.length.diet)
winter.length.rare.prey.adjust

winter.length.diet <- rbind(winter.length.diet, winter.length.rare.prey.adjust)
winter.length.diet <- winter.length.diet * 100
colSums(winter.length.diet)
winter.length.diet$Group <- rownames(winter.length.diet)

winter.length.diet.long <- melt(winter.length.diet)
names(winter.length.diet.long) <- c("Group", "Length bin (cm)", "Percentage (%)")
winter.length.diet.long$Season <- "Winter"

## plot both seasons together
both.seasons.length.diet <- rbind(summer.length.diet.long, winter.length.diet.long)

sort(unique(both.seasons.length.diet$Group))
both.seasons.length.diet$PlotGroup <- both.seasons.length.diet$Group
both.seasons.length.diet$PlotGroup[both.seasons.length.diet$Group == "Pacific herring"] <- "Herring"
both.seasons.length.diet$PlotGroup[both.seasons.length.diet$Group == "Pacific sand lance"] <- "Sand lance"
both.seasons.length.diet$PlotGroup[both.seasons.length.diet$Group == "Northern anchovy"] <- "Anchovy"

both.seasons.length.diet$PlotGroup[both.seasons.length.diet$Group %in% c("Barracudina",
  "Deep-sea smelt", "Eulachon", "Rockfish", "Salmon", "Sculpin", 
  "Threespine stickleback", "Tubesnout")] <- "Other fish"
both.seasons.length.diet$PlotGroup[both.seasons.length.diet$Group %in% c("Amphipod",
  "Crustacean larvae", "Planktonic shrimp")] <- "Other crustaceans"
both.seasons.length.diet$PlotGroup[both.seasons.length.diet$Group == "Euphausiid"] <- "Krill"
both.seasons.length.diet$PlotGroup[both.seasons.length.diet$Group == "Gadiformes"] <- "Cod and hake"

greater.5 <- subset(both.seasons.length.diet, `Percentage (%)` > 5)
sort(unique(greater.5$PlotGroup))

sort(unique(both.seasons.length.diet$PlotGroup))
both.seasons.length.diet$PlotGroup <- factor(both.seasons.length.diet$PlotGroup,
  levels = c("Herring", "Sand lance", "Anchovy", "Cod and hake", "Myctophid", "Surfperch",
             "Other fish", "Krill", "Mysid", "Pandalid", "Other crustaceans",
             "Polychaete", "Squid", "Rare prey"))

summer.length.bin.n <- data.frame(table(summer.with.length$length_bin))
summer.length.bin.n$season <- "Summer"
winter.length.bin.n <- data.frame(table(winter.with.length$length_bin))
winter.length.bin.n$season <- "Winter"

both.seasons.length.bin.n <- rbind(summer.length.bin.n, winter.length.bin.n)
names(both.seasons.length.bin.n) <- c("Length bin (cm)", "N", "Season")

ggplot() +
  geom_bar(data = both.seasons.length.diet,
           aes(x = `Length bin (cm)`, weight = `Percentage (%)`, fill = PlotGroup)) +
  scale_fill_manual(values = c("#08306b", "#08519c", "#2171b5", "#4292c6",
    "#6baed6", "#9ecae1", "#c6dbef", "#67000d", "#a50f15", "#cb181d", "#ef3b2c",
    "#fb6a4a", "#fc9272", "grey60")) +
  geom_text(data = both.seasons.length.bin.n, 
            aes(x = `Length bin (cm)`, y = 95, label = N), 
            colour = "white", size = 3.5) +
  facet_wrap(~ Season, ncol = 1) +
  ylab("Percentage (%)") + labs(fill = "Group") +
  theme_bw() +
  theme(legend.title = element_text(hjust = 0.5, size = 11),
        legend.text = element_text(size = 10),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 11),
        strip.background = element_blank(),
        strip.text.x = element_text(size = 12),
        panel.border = element_rect(colour = "black", linewidth = 0.9),
        legend.position = "right")
ggsave("figures/supplement/size-diet.PNG", 
       width = 21, height = 14, units = "cm", dpi = 600)

## 3. Cluster analysis ----
summer55 <- subset(summer.salish.ch, `Final Length (cm)` >= 55) # 847
summer62 <- subset(summer.salish.ch, `Final Length (cm)` >= 62) # 822
winter55 <- subset(winter.salish.ch, `Final Length (cm)` >= 55) # 678
winter62 <- subset(winter.salish.ch, `Final Length (cm)` >= 62) # 621

## summer >= 55 cm
#summer55 %>% group_by(GridID) %>% tally() %>% view() # all have at least 10 stomachs

summer55.grid.levels <- sort(unique(summer55$GridID))
names(summer55)
summer55.prop <- data.frame(matrix(NA, nrow = 16, ncol = 0)) # need to change

for (i in 1:length(summer55.grid.levels)) {
  grid.subset <- subset(summer55, GridID == summer55.grid.levels[i])
  mean.diet <- colMeans(grid.subset[20:35])
  summer55.prop <- cbind(summer55.prop, mean.diet)
  colnames(summer55.prop)[i] <- summer55.grid.levels[i]
}

colSums(summer55.prop) # expected

# calculate Bray-Curtis dissimilarity matrix
summer55.prop.t <- as.data.frame(t(summer55.prop))
summer55.bray <- vegdist(summer55.prop.t, method = "bray")

# hierarchical cluster analysis
summer55.cluster <- hclust(summer55.bray, "average")
plot(summer55.cluster, main = "", sub = "", xlab = "",
     ylab = "Bray-Curtis dissimilarity", hang = -1)

# determine optimal number of clusters using indicator species analysis
summer55.indval.fill <- data.frame(k = 2:10, sum.sig.indicators = NA, 
                                prop.sig.indicators = NA)

for (z in 2:10) {
  z.clusters <- data.frame(Cluster = cutree(summer55.cluster, k = z))
  z.clusters$GridID <- rownames(z.clusters)
  
  indiv.diet <- merge(x = summer55, y = z.clusters, by = "GridID",
                      all.x = FALSE, all.y = FALSE)

  set.seed(530)
  z.indval <- indval(indiv.diet[20:35], indiv.diet$Cluster, numitr = 1000)
  
  # calculate sum of significant indicators
  significant.indicators <- z.indval$indcls[z.indval$pval < 0.05]
  sum.sig.indicators <- sum(significant.indicators)
  
  # calculate prop with sig indicators
  clusters.with.sig.indval <- factor(z.indval$maxcls[z.indval$pval < 0.05])
  prop.sig.indicators <- length(levels(clusters.with.sig.indval)) / z
  
  summer55.indval.fill$sum.sig.indicators[summer55.indval.fill$k == z] <- sum.sig.indicators
  summer55.indval.fill$prop.sig.indicators[summer55.indval.fill$k == z] <- prop.sig.indicators
}

summer55.indval.fill[which.max(summer55.indval.fill$sum.sig.indicators), ]
# k = 3 is supported

summer55.n10.grid <- subset(grid, GridID %in% summer55.grid.levels)
summer55.cluster.cut <- data.frame(Cluster = cutree(summer55.cluster, k = 3))
summer55.cluster.cut$GridID <- rownames(summer55.cluster.cut)
summer55.n10.grid <- merge(x = summer55.n10.grid, y = summer55.cluster.cut,
                           by = "GridID", all.x = FALSE, all.y = FALSE)
summer55.n10.grid$Cluster <- as.character(summer55.n10.grid$Cluster)

summer55.n10.grid$Cluster[summer55.n10.grid$GridID %in% c(15, 16, 23)] <- "4"
summer55.n10.grid$GridID[summer55.n10.grid$GridID == 103] <- 99

summer55.map <- ggplot() +
  geom_sf(data = bc.coast) +
  geom_sf(data = summer55.n10.grid, aes(fill = Cluster), alpha = 0.8,
          show.legend = FALSE) +
  geom_sf_text(data = summer55.n10.grid, aes(label = GridID)) +
  scale_fill_manual(values = c("#0072B2", "#D55E00", "#009E73", "#56b4e9")) +
  coord_sf(datum = 4326,
           crs = st_crs(32610),
           xlim = c(300*10^3, 540*10^3), ylim = c(5325*10^3, 5600*10^3)) +
  theme_bw() + 
  theme(axis.title = element_blank(),
        panel.grid = element_blank(),
        axis.text = element_text(size = 10),
        legend.position = "inside",
        legend.position.inside = c(0.88, 0.82),
        legend.background = element_blank()) +
  ggtitle(paste0("Summer >= 55 cm"))
summer55.map

## 3B summer >= 62 cm ----
#summer62 %>% group_by(GridID) %>% tally() %>% view() # all have at least 10

summer62.grid.levels <- sort(unique(summer62$GridID))
names(summer62)
summer62.prop <- data.frame(matrix(NA, nrow = 16, ncol = 0)) # need to change

for (i in 1:length(summer62.grid.levels)) {
  grid.subset <- subset(summer62, GridID == summer62.grid.levels[i])
  mean.diet <- colMeans(grid.subset[20:35])
  summer62.prop <- cbind(summer62.prop, mean.diet)
  colnames(summer62.prop)[i] <- summer62.grid.levels[i]
}

colSums(summer62.prop) # expected

# calculate Bray-Curtis dissimilarity matrix
summer62.prop.t <- as.data.frame(t(summer62.prop))
summer62.bray <- vegdist(summer62.prop.t, method = "bray")

# hierarchical cluster analysis
summer62.cluster <- hclust(summer62.bray, "average")
plot(summer62.cluster, main = "", sub = "", xlab = "",
     ylab = "Bray-Curtis dissimilarity", hang = -1)

summer62.indval.fill <- data.frame(k = 2:10, sum.sig.indicators = NA, 
                                   prop.sig.indicators = NA)

for (z in 2:10) {
  z.clusters <- data.frame(Cluster = cutree(summer62.cluster, k = z))
  z.clusters$GridID <- rownames(z.clusters)
  
  indiv.diet <- merge(x = summer62, y = z.clusters, by = "GridID",
                      all.x = FALSE, all.y = FALSE)
  nrow(indiv.diet) # should always be 968
  
  set.seed(530)
  z.indval <- indval(indiv.diet[20:35], indiv.diet$Cluster, numitr = 1000)
  
  # calculate sum of significant indicators
  significant.indicators <- z.indval$indcls[z.indval$pval < 0.05]
  sum.sig.indicators <- sum(significant.indicators)
  
  # calculate prop with sig indicators
  clusters.with.sig.indval <- factor(z.indval$maxcls[z.indval$pval < 0.05])
  prop.sig.indicators <- length(levels(clusters.with.sig.indval)) / z
  
  summer62.indval.fill$sum.sig.indicators[summer62.indval.fill$k == z] <- sum.sig.indicators
  summer62.indval.fill$prop.sig.indicators[summer62.indval.fill$k == z] <- prop.sig.indicators
}

summer62.indval.fill[which.max(summer62.indval.fill$sum.sig.indicators), ]
# k = 3 is supported

summer62.n10.grid <- subset(grid, GridID %in% summer62.grid.levels)
summer62.cluster.cut <- data.frame(Cluster = cutree(summer62.cluster, k = 3))
summer62.cluster.cut$GridID <- rownames(summer62.cluster.cut)
summer62.n10.grid <- merge(x = summer62.n10.grid, y = summer62.cluster.cut,
                           by = "GridID", all.x = FALSE, all.y = FALSE)
summer62.n10.grid$Cluster <- as.character(summer62.n10.grid$Cluster)

summer62.n10.grid$Cluster[summer62.n10.grid$GridID %in% c(15, 16, 23, 17)] <- "4"

summer62.n10.grid$GridID[summer62.n10.grid$GridID == 103] <- 99

summer62.map <- ggplot() +
  geom_sf(data = bc.coast) +
  geom_sf(data = summer62.n10.grid, aes(fill = Cluster), alpha = 0.8,
          show.legend = FALSE) +
  geom_sf_text(data = summer62.n10.grid, aes(label = GridID)) +
  scale_fill_manual(values = c("#0072B2", "#D55E00", "#009E73", "#56b4e9")) +
  coord_sf(datum = 4326,
           crs = st_crs(32610),
           xlim = c(300*10^3, 540*10^3), ylim = c(5325*10^3, 5600*10^3)) +
  theme_bw() + 
  theme(axis.title = element_blank(),
        panel.grid = element_blank(),
        axis.text = element_text(size = 10),
        legend.position = "inside",
        legend.position.inside = c(0.88, 0.82),
        legend.background = element_blank()) +
  ggtitle(paste0("Summer >= 62 cm"))
summer62.map

ggarrange(summer55.map, summer62.map, labels = c("A", "B"))
ggsave("figures/summer/summer-size-maps.PNG", width = 25, height = 13, units = "cm",
       bg = "white", dpi = 600)

## 3C winter >= 55 cm ----
winter55 %>% group_by(GridID) %>% tally() %>% view() # grid 54 only has 8 stomachs now

winter55.grid.levels <- sort(unique(winter55$GridID))
names(winter55)
winter55.prop <- data.frame(matrix(NA, nrow = 18, ncol = 0)) # need to change

for (i in 1:length(winter55.grid.levels)) {
  grid.subset <- subset(winter55, GridID == winter55.grid.levels[i])
  mean.diet <- colMeans(grid.subset[20:37])
  winter55.prop <- cbind(winter55.prop, mean.diet)
  colnames(winter55.prop)[i] <- winter55.grid.levels[i]
}

colSums(winter55.prop) # expected

# calculate Bray-Curtis dissimilarity matrix
winter55.prop.t <- as.data.frame(t(winter55.prop))
winter55.bray <- vegdist(winter55.prop.t, method = "bray")

# hierarchical cluster analysis
winter55.cluster <- hclust(winter55.bray, "average")
plot(winter55.cluster, main = "", sub = "", xlab = "",
     ylab = "Bray-Curtis dissimilarity", hang = -1)

winter55.indval.fill <- data.frame(k = 2:10, sum.sig.indicators = NA, 
                                   prop.sig.indicators = NA)

for (z in 2:10) {
  z.clusters <- data.frame(Cluster = cutree(winter55.cluster, k = z))
  z.clusters$GridID <- rownames(z.clusters)
  
  indiv.diet <- merge(x = winter55, y = z.clusters, by = "GridID",
                      all.x = FALSE, all.y = FALSE)
  nrow(indiv.diet) # should always be 968
  
  set.seed(530)
  z.indval <- indval(indiv.diet[20:37], indiv.diet$Cluster, numitr = 1000)
  
  # calculate sum of significant indicators
  significant.indicators <- z.indval$indcls[z.indval$pval < 0.05]
  sum.sig.indicators <- sum(significant.indicators)
  
  # calculate prop with sig indicators
  clusters.with.sig.indval <- factor(z.indval$maxcls[z.indval$pval < 0.05])
  prop.sig.indicators <- length(levels(clusters.with.sig.indval)) / z
  
  winter55.indval.fill$sum.sig.indicators[winter55.indval.fill$k == z] <- sum.sig.indicators
  winter55.indval.fill$prop.sig.indicators[winter55.indval.fill$k == z] <- prop.sig.indicators
}

winter55.indval.fill[which.max(winter55.indval.fill$sum.sig.indicators), ]

# make map
cutree(winter55.cluster, k = 5)
winter55.n10.grid <- subset(grid, GridID %in% winter55.grid.levels)
winter55.cluster.cut <- data.frame(Cluster = cutree(winter55.cluster, k = 5))
winter55.cluster.cut$GridID <- rownames(winter55.cluster.cut)
winter55.n10.grid <- merge(x = winter55.n10.grid, y = winter55.cluster.cut,
                           by = "GridID", all.x = FALSE, all.y = FALSE)
winter55.n10.grid$Cluster <- as.character(winter55.n10.grid$Cluster)

winter55.n10.grid$Cluster[winter55.n10.grid$GridID == 16] <- "2b"
winter55.n10.grid$Cluster[winter55.n10.grid$GridID == 17] <- "3b"

ggplot() +
  geom_sf(data = bc.coast) +
  geom_sf(data = winter55.n10.grid, aes(fill = Cluster), alpha = 0.8,
          show.legend = FALSE) +
  geom_sf_text(data = winter55.n10.grid, aes(label = GridID)) +
  scale_fill_manual(values = c("#deb887", "#0072B2", "#56b4e9", "#5823AD", "#976EDB", 
                               "#D55E00", "#009E73")) +
  coord_sf(datum = 4326,
           crs = st_crs(32610),
           xlim = c(300*10^3, 540*10^3), ylim = c(5325*10^3, 5600*10^3)) +
  theme_bw() + 
  theme(axis.title = element_blank(),
        panel.grid = element_blank(),
        axis.text = element_text(size = 10)) +
  ggtitle(paste0("Winter >= 55 cm"))
ggsave("figures/winter/winter55-map.PNG", width = 17, height = 17, 
       units = "cm", bg = "white", dpi = 600)

## 3D winter >= 62 cm ----
winter62 %>% group_by(GridID) %>% tally() %>% view() 
# grid 54 only has 8 stomachs now, grid 17 only has 6, grid 15 only has 4

# first run cluster analysis including all 18 grid cells
winter62.grid.levels <- sort(unique(winter62$GridID))

colSums(winter62[20:37])
winter62 <- select(winter62, -Tubesnout) # no longer any tubesnout at min 62 cm

names(winter62)
winter62.prop <- data.frame(matrix(NA, nrow = 17, ncol = 0)) # need to change

for (i in 1:length(winter62.grid.levels)) {
  grid.subset <- subset(winter62, GridID == winter62.grid.levels[i])
  mean.diet <- colMeans(grid.subset[20:36])
  winter62.prop <- cbind(winter62.prop, mean.diet)
  colnames(winter62.prop)[i] <- winter62.grid.levels[i]
}

colSums(winter62.prop) # expected

# calculate Bray-Curtis dissimilarity matrix
winter62.prop.t <- as.data.frame(t(winter62.prop))
winter62.bray <- vegdist(winter62.prop.t, method = "bray")

# hierarchical cluster analysis
winter62.cluster <- hclust(winter62.bray, "average")
plot(winter62.cluster, main = "", sub = "", xlab = "",
     ylab = "Bray-Curtis dissimilarity", hang = -1)

winter62.indval.fill <- data.frame(k = 2:10, sum.sig.indicators = NA, 
                                   prop.sig.indicators = NA)

for (z in 2:10) {
  z.clusters <- data.frame(Cluster = cutree(winter62.cluster, k = z))
  z.clusters$GridID <- rownames(z.clusters)
  
  indiv.diet <- merge(x = winter62, y = z.clusters, by = "GridID",
                      all.x = FALSE, all.y = FALSE)
  nrow(indiv.diet) # should always be 968
  
  set.seed(530)
  z.indval <- indval(indiv.diet[20:36], indiv.diet$Cluster, numitr = 1000)
  # 20:36 since no tubesnout
  
  # calculate sum of significant indicators
  significant.indicators <- z.indval$indcls[z.indval$pval < 0.05]
  sum.sig.indicators <- sum(significant.indicators)
  
  # calculate prop with sig indicators
  clusters.with.sig.indval <- factor(z.indval$maxcls[z.indval$pval < 0.05])
  prop.sig.indicators <- length(levels(clusters.with.sig.indval)) / z
  
  winter62.indval.fill$sum.sig.indicators[winter62.indval.fill$k == z] <- sum.sig.indicators
  winter62.indval.fill$prop.sig.indicators[winter62.indval.fill$k == z] <- prop.sig.indicators
}

winter62.indval.fill[which.max(winter62.indval.fill$sum.sig.indicators), ]
# k = 4 is supported

# make min 62 cm, all grid cell map
cutree(winter62.cluster, k = 4)
winter62.n10.grid <- subset(grid, GridID %in% winter62.grid.levels)
winter62.cluster.cut <- data.frame(Cluster = cutree(winter62.cluster, k = 4))
winter62.cluster.cut$GridID <- rownames(winter62.cluster.cut)
winter62.n10.grid <- merge(x = winter62.n10.grid, y = winter62.cluster.cut,
                           by = "GridID", all.x = FALSE, all.y = FALSE)
winter62.n10.grid$Cluster <- as.character(winter62.n10.grid$Cluster)

# # divide into regions
# winter62.n10.grid$Cluster[winter62.n10.grid$GridID %in% c(16, 17)] <- "2B"
# winter62.n10.grid$Cluster[winter62.n10.grid$GridID == 15] <- "1A"

winter62.map.all <- ggplot() +
  geom_sf(data = bc.coast) +
  geom_sf(data = winter62.n10.grid, aes(fill = Cluster), alpha = 0.8) +
  geom_sf_text(data = winter62.n10.grid, aes(label = GridID)) +
  scale_fill_manual(values = c("#009E73", "#0072B2", "#D55E00", "#5823AD")) +
  
  coord_sf(datum = 4326,
           crs = st_crs(32610),
           xlim = c(300*10^3, 540*10^3), ylim = c(5325*10^3, 5600*10^3)) +
  theme_bw() + 
  theme(axis.title = element_blank(),
        panel.grid = element_blank(),
        axis.text = element_text(size = 10),
        legend.position = "inside",
        legend.position.inside = c(0.88, 0.86),
        legend.background = element_blank()) +
  ggtitle(paste0("Winter >= 62 cm, all grid cells"))
winter62.map.all

## repeat cluster analysis without grid cells with <10 stomachs
winter62.grid.tally <- data.frame(table(winter62$GridID))
winter62.grid.tally10 <- subset(winter62.grid.tally, Freq >= 10)

winter62.prop.t$GridID <- rownames(winter62.prop.t)
winter62.prop.t10 <- subset(winter62.prop.t, GridID %in% winter62.grid.tally10$Var1) # 15
winter62.prop.t10 <- select(winter62.prop.t10, -GridID)

winter62.bray10 <- vegdist(winter62.prop.t10, method = "bray")

# hierarchical cluster analysis
winter62.cluster10 <- hclust(winter62.bray10, "average")
plot(winter62.cluster10, main = "", sub = "", xlab = "",
     ylab = "Bray-Curtis dissimilarity", hang = -1)

winter62.indval.fill10 <- data.frame(k = 2:10, sum.sig.indicators = NA, 
                                   prop.sig.indicators = NA)
names(winter62)

for (z in 2:10) {
  z.clusters <- data.frame(Cluster = cutree(winter62.cluster10, k = z))
  z.clusters$GridID <- rownames(z.clusters)
  
  indiv.diet <- merge(x = winter62, y = z.clusters, by = "GridID",
                      all.x = FALSE, all.y = FALSE)
  nrow(indiv.diet) # should always be 968
  
  set.seed(530)
  z.indval <- indval(indiv.diet[20:36], indiv.diet$Cluster, numitr = 1000)
  
  # calculate sum of significant indicators
  significant.indicators <- z.indval$indcls[z.indval$pval < 0.05]
  sum.sig.indicators <- sum(significant.indicators)
  
  # calculate prop with sig indicators
  clusters.with.sig.indval <- factor(z.indval$maxcls[z.indval$pval < 0.05])
  prop.sig.indicators <- length(levels(clusters.with.sig.indval)) / z
  
  winter62.indval.fill10$sum.sig.indicators[winter62.indval.fill10$k == z] <- sum.sig.indicators
  winter62.indval.fill10$prop.sig.indicators[winter62.indval.fill10$k == z] <- prop.sig.indicators
}


winter62.indval.fill10[which.max(winter62.indval.fill10$sum.sig.indicators), ]

# make min 62 cm, grid cells with >= 10 stomachs map
winter62.indval.fill10[which.max(winter62.indval.fill10$sum.sig.indicators), ]
cutree(winter62.cluster10, k = 4)
winter62.n10.grid10 <- subset(grid, GridID %in% winter62.grid.tally10$Var1)
winter62.cluster.cut10 <- data.frame(Cluster = cutree(winter62.cluster10, k = 4))
winter62.cluster.cut10$GridID <- rownames(winter62.cluster.cut10)
winter62.n10.grid10 <- merge(x = winter62.n10.grid10, y = winter62.cluster.cut10,
                           by = "GridID", all.x = FALSE, all.y = FALSE)
winter62.n10.grid10$Cluster <- as.character(winter62.n10.grid10$Cluster)

winter62.map10 <- ggplot() +
  geom_sf(data = bc.coast) +
  geom_sf(data = winter62.n10.grid10, aes(fill = Cluster), alpha = 0.8) +
  geom_sf_text(data = winter62.n10.grid10, aes(label = GridID)) +
  scale_fill_manual(values = c("#0072B2", "#D55E00", "#009E73", "#5823AD")) +
  coord_sf(datum = 4326,
           crs = st_crs(32610),
           xlim = c(300*10^3, 540*10^3), ylim = c(5325*10^3, 5600*10^3)) +
  theme_bw() + 
  theme(axis.title = element_blank(),
        panel.grid = element_blank(),
        axis.text = element_text(size = 10),
        legend.position = "inside",
        legend.position.inside = c(0.88, 0.86),
        legend.background = element_blank()) +
  ggtitle(paste0("Winter >= 62 cm, only grid cells with >= 10 stomachs"))
winter62.map10

ggarrange(winter62.map.all, winter62.map10, labels = c("A", "B"))
ggsave("figures/winter/winter-62-maps.PNG", 
       width = 26, height = 13, units = "cm", bg = "white", dpi = 600)

## 5. Map of PFMAs given size regulations ----
# PFMA shapefile available at Pacific Salmon Foundation Marine Data Centre
# https://soggy2.zoology.ubc.ca/geonetwork/srv/eng/catalog.search#/metadata/25ebea3f-f5df-4304-baae-d39065aacee3

pfma <- read_sf("data/DFO-PFMAs/DFO_BC_PFMA_SUBAREAS_CHS_V3_G.shp")
head(pfma)

salish.sea.pfma <- subset(pfma, MGNT_AREA %in% c(13, 14, 15, 16, 17, 18, 19, 20,
                                                 28, 29))
salish.sea.pfma$size_regulation <- "Minimum 62 cm"
salish.sea.pfma$size_regulation[salish.sea.pfma$MGNT_AREA == 20] <- "Minimum 45 cm"
salish.sea.pfma$size_regulation[salish.sea.pfma$MGNT_AREA == 19 &
  salish.sea.pfma$SUBAREA %in% c(1, 2, 3, 4)] <- "Minimum 45 cm"

ggplot() +
  geom_sf(data = bc.coast) +
  geom_sf(data = salish.sea.pfma, aes(fill = size_regulation)) +
  scale_fill_manual(values = c("#50C96A", "#6eb0ff")) +
  coord_sf(datum = 4326,
           crs = st_crs(32610),
           xlim = c(300*10^3, 540*10^3), ylim = c(5325*10^3, 5600*10^3)) +
  labs(fill = "Minimum size regulation") + 
  theme_bw() + 
  theme(axis.title = element_blank(),
        panel.grid = element_blank(),
        axis.text = element_text(size = 10),
        legend.title = element_text(size = 11),
        legend.text = element_text(size = 10),
        legend.position = "inside",
        legend.position.inside = c(0.78, 0.92),
        legend.background = element_blank())
ggsave("figures/maps/size-regulations.PNG", 
       width = 14, height = 14, units = "cm", dpi = 600)
