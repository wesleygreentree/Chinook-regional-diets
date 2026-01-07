# Winter analysis 

## 0. Script description ----

# 1. Spatial processing: transforms longitude-latitude coordinates into UTM
# zone 10, so that they can be assigned to 25 x 25 km grid cells. The grid
# used is the one determined to be most appropriate for the sampling coverage
# in scripts/select-grid.R

# 2. Process diet composition: stomachs that were empty or only contained
# unidentified material are removed, diet proportions for analysis are 
# calculated as % identified weight. 

# 3. Remove rare prey: rare prey taxa are removed prior to multivariate analyses, using
# a vector of summer rare prey taxa identified in scripts/data-processing.R
# Rare prey were defined as taxa contributing <0.5% mean identified weight
# in all Pacific Fishery Management Areas. 

# 4. Cluster analysis: Mean diet composition in each grid cell is calculated,
# a Bray-Curtis dissimilarity matrix is calculated, and hierarchical cluster
# analysis is conducted with the Bray-Curtis dissimilarity matrix

# 5. Optimal number of clusters: To interpret results, we need to define the
# appropriate number of clusters. Two methods are used: 
# (A) silhouette width identifies the number of clusters that optimizes 
# within-cluster dissimilarity to be low and among-cluster dissimilarity to be
# high (i.e. clusters separate well). The mean silhouette width calculation is
# modified from the original (Rousseeuw 1987) to avoid penalizing clusters with
# only one grid cell, as those are still ecologically meaningful.
# (B) Indicator species analysis, which identifies the optimal number of
# clusters as the typology where each cluster has a significant indicator 
# species (p < 0.05) and where the sum of significant indicator species among
# clusters is greater than at any other number of clusters.

# 6-8: visualization of cluster analysis results. Map, dendrogram, bar charts
# of mean diet composition in each region identified.

# 9. Principal coordinate analysis (PCoA): PCoA is a multivariate ordination
# method that is used here to provide insight into the strength of
# differentiation among regions, as the variation explained by each axis (R^2)
# can be calculated. A Lingoes correction was used to correct for negative
# eigenvalues resulting from the Bray-Curtis dissimilarity matrix (Legendre
# and Legendre 2012).

# 10. Indicator species are determined for each region, following Dufrene and 
# Legendre (1997).

# 11. Produce distribution maps for select prey.

# 12. Seasonal sensitivity analyses: determine if results robust to the 
# designation of winter as October to March.

# load packages
library(tidyverse)
library(sf)
library(reshape2)
library(vegan)
library(cluster)
library(ape)
library(labdsv)
library(ggdendro)
library(ggpubr)

# read data
bc.coast <- read_sf("data/BC-coast-shapefile/bc-coast.shp")
grid <- readRDS("data/out/grid_17_5.RDS")
attr(grid, "grid") # 17_5
salish.ch <- readRDS("data/out/Salish-Sea-Chinook-salmon-diet.RDS")
nrow(salish.ch) # 2655
table(salish.ch$Season) # 1589 summer, 1066 winter

## 1. Spatial processing ----

# remove stomachs without a creel site
table(salish.ch$`Site Resolution`) # 2495 creel site, 160 PFMA or PFMA subarea
salish.ch <- subset(salish.ch, `Site Resolution` == "creel site")
nrow(salish.ch) # 2495

# convert to UTM
bc.coastUTM <- st_transform(bc.coast, crs = st_crs(32610))

salish.ch.sf <- st_as_sf(salish.ch, coords = c("Longitude","Latitude"),
                         crs = st_crs(4326), remove = FALSE)
salish.chUTM <- st_transform(salish.ch.sf, crs = st_crs(32610))

# join to grid
joined <- st_join(x = salish.chUTM, y = grid)

# pull out UTM coordinates and remove geometry column
joined.coords <- do.call(rbind, st_geometry(joined)) %>% 
  as_tibble() %>% setNames(c("utm_x", "utm_y"))
joined.without.geometry <- st_drop_geometry(joined)
joined <- cbind(joined.without.geometry, joined.coords)

summer.salish.ch <- subset(joined, Season == "Summer") # 1479
winter.salish.ch <- subset(joined, Season == "Winter") # 1016

winter.months <- winter.salish.ch %>% group_by(Year, Month) %>% tally()
winter.months$EcologicalYear <- winter.months$Year
winter.months$EcologicalYear[winter.months$Month %in% c(1, 2, 3)] <- winter.months$Year[winter.months$Month %in% c(1, 2, 3)] - 1

winter.months$Month <- factor(winter.months$Month, levels = c(10, 11, 12, 1, 2, 3))
winter.months$Winter <- NA
winter.months$Winter[winter.months$EcologicalYear == 2017] <- "Winter 2017-18"
winter.months$Winter[winter.months$EcologicalYear == 2018] <- "Winter 2018-19"
winter.months$Winter[winter.months$EcologicalYear == 2019] <- "Winter 2019-20"
winter.months$Winter[winter.months$EcologicalYear == 2020] <- "Winter 2020-21"
winter.months$Winter[winter.months$EcologicalYear == 2021] <- "Winter 2021-22"

ggplot() +
  geom_col(data = winter.months, aes(x = Month, y = n, fill = Month),
           show.legend = FALSE) +
  scale_fill_viridis_d(option = "plasma", direction = -1) +
  ylab("Number of stomachs") +
  theme_bw() +
  theme(axis.title = element_text(size = 12),
        axis.text = element_text(size = 11),
        strip.text = element_text(size = 12),
        strip.background = element_blank()) +
  facet_wrap(~ Winter, ncol = 1, scales = "free_y")
ggsave("figures/winter/winter-month-n.PNG", width = 17, height = 14, units = "cm",
       dpi = 1600)

# 2. Process diet composition ----

# save coordinates of all samples for maps
all.winter <- winter.salish.ch

# remove empty stomachs
nrow(winter.salish.ch) # 1016
nrow(subset(winter.salish.ch, TotalProp == 0)) # 174
winter.salish.ch <- subset(winter.salish.ch, TotalProp != 0)
winter.salish.ch <- winter.salish.ch[, !names(winter.salish.ch) == "Empty"]
nrow(winter.salish.ch) # 842
174/1016 # 17.1%

names(winter.salish.ch) # 16:54
sum(colMeans(winter.salish.ch[16:54])) # 1 # good

# remove stomachs with only unidentified material
nrow(subset(winter.salish.ch, `Unidentified material` == 1))
winter.salish.ch <- subset(winter.salish.ch, `Unidentified material` != 1)
nrow(winter.salish.ch) # 762

# calculate % identified weight
names(winter.salish.ch)
winter.salish.ch$IdentifiedTotal <- rowSums(winter.salish.ch[16:53])

winter.salish.ch[16:53] <- winter.salish.ch[16:53] / winter.salish.ch$IdentifiedTotal
unique(rowSums(winter.salish.ch[16:53])) # 1 # good

# Select grids with n >= 10
winter.tally <- winter.salish.ch %>% group_by(GridID) %>% tally()
winter.tally <- subset(winter.tally, n >= 10) # 18 grid cells
winter.n10.grids <- winter.tally$GridID
winter.n10.grids

# save data for species accumulation curves before rare prey are removed
saveRDS(winter.salish.ch, "data/out/winter-Salish-Chinook-for-specaccum.RDS")

winter.salish.ch <- subset(winter.salish.ch, GridID %in% winter.n10.grids)
nrow(winter.salish.ch) # 723

## 3. Remove rare prey ----
winter.major.prey <- readRDS("data/out/winter-major-prey-vector.RDS") # from scripts/data-processing.R
names(winter.salish.ch)

# exclude rare prey for analysis
winter.salish.ch <- select(winter.salish.ch, `Fish Code`, Species, Season, 
  `Date Resolution`, Year, Month, Day, Region, PFMA, `Site Resolution`,
  Longitude, Latitude, utm_x, utm_y, `Final Length (cm)`, Project, `Data flag`,
  GridID, all_of(winter.major.prey))

names(winter.salish.ch)
major.prey.row.sum <- rowSums(winter.salish.ch[19:36])
hist(major.prey.row.sum)
view(major.prey.row.sum) # only 1 < 99%, so that's good
# I don't normalize % identified to sum to 1 after removing major prey

# save for size sensitivity analysis
saveRDS(winter.salish.ch, "data/out/winter-diet-for-size-sensitivity.RDS")

# 4. Cluster analysis ----

ncol(winter.salish.ch[19:36]) # 16
winter.grid.prop <- data.frame(matrix(NA, nrow = 18, ncol = 0)) 

for (i in 1:length(winter.n10.grids)) {
  grid.subset <- subset(winter.salish.ch, GridID == winter.n10.grids[i])
  mean.diet <- colMeans(grid.subset[19:36])
  winter.grid.prop <- cbind(winter.grid.prop, mean.diet)
  colnames(winter.grid.prop)[i] <- winter.n10.grids[i]
}

colSums(winter.grid.prop) # only cell 16 < 1, expected

# calculate Bray-Curtis dissimilarity matrix
winter.grid.prop.t <- as.data.frame(t(winter.grid.prop))
grid.bray <- vegdist(winter.grid.prop.t, method = "bray")

# hierarchical cluster analysis
grid.cluster <- hclust(grid.bray, "average")
plot(grid.cluster, main = "", sub = "", xlab = "",
     ylab = "Bray-Curtis dissimilarity", hang = -1)

# cophenetic correlation - does the dendrogram reflect the dissimilarity matrix
coph <- cophenetic(grid.cluster)
cor(grid.bray, coph) # 0.92, good fit

## 5. Optimal number of clusters ----

## 5A. Silhouette width ----
# Calculate silhouette width to determine optimal number of clusters
# A corrected version is used that removes the penalty against single clusters.
# In the original method (Rousseeuw 1987), single-member clusters receive
# silhouette = 0. Here, they are excluded from the calculation of mean
# silhouette width. 
# The uncorrected silhouette width is also calculated

sil.fill <- data.frame(k = 2:10, uncorrected_silhouette = NA,
                       corrected_silhouette = NA)
for (b in 2:10) {
  sil <- silhouette(cutree(grid.cluster, k = b), grid.bray)
  sil.fill$uncorrected_silhouette[sil.fill$k == b] <- summary(sil)$avg.width
  
  sil <- data.frame(sil[, 1:3])
  
  freq <- data.frame(table(sil$cluster))
  
  n1 <- subset(freq, Freq == 1)
  n1.clusters <- n1$Var1
  
  corrected.sil <- subset(sil, !cluster %in% n1.clusters)
  sil.fill$corrected_silhouette[sil.fill$k == b] <- mean(corrected.sil$sil_width)
  
}

corrected.sil.plot <- ggplot() +
  geom_col(data = sil.fill, aes(x = factor(k), y = corrected_silhouette),
           fill = "grey50") +
  xlab("Number of clusters") + ylab("Corrected mean silhouette width") +
  scale_y_continuous(expand = c(0, 0)) +
  theme_minimal() +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12))
corrected.sil.plot

## 5B. Indicator species ----
# While the cluster analysis is conducted at the grid-level, I conducted
# indicator species analysis at the individual stomach level, which mirrors
# the calculation of indicator species analysis for each region identified by
# the cluster analysis.

# Potentially, one might argue that the indicator species should be calculated
# using average grid diet composition. I disagree, as this method is likely to
# over-inflate indicator values. Still, I wanted to check that the number of 
# clusters suggested is the same (line 314).

# I follow the approach of Borcard et al. 2018. Numerical Ecology with R.

# calculate IndVal at stomach level
indiv.indval.fill <- data.frame(k = 2:10, sum.sig.indicators = NA, 
                                prop.sig.indicators = NA)

for (z in 2:10) {
  z.clusters <- data.frame(Cluster = cutree(grid.cluster, k = z))
  z.clusters$GridID <- rownames(z.clusters)
  
  indiv.diet <- merge(x = winter.salish.ch, y = z.clusters, by = "GridID",
                      all.x = FALSE, all.y = FALSE)
  nrow(indiv.diet) # should always be 723
  
  set.seed(530)
  z.indval <- indval(indiv.diet[19:36], indiv.diet$Cluster, numitr = 1000)
  
  # calculate sum of significant indicators
  significant.indicators <- z.indval$indcls[z.indval$pval < 0.05]
  sum.sig.indicators <- sum(significant.indicators)
  
  # calculate prop with sig indicators
  clusters.with.sig.indval <- factor(z.indval$maxcls[z.indval$pval < 0.05])
  prop.sig.indicators <- length(levels(clusters.with.sig.indval)) / z
  
  indiv.indval.fill$sum.sig.indicators[indiv.indval.fill$k == z] <- sum.sig.indicators
  indiv.indval.fill$prop.sig.indicators[indiv.indval.fill$k == z] <- prop.sig.indicators
}

indiv.indval.fill$prop_binary <- "<1"
indiv.indval.fill$prop_binary[indiv.indval.fill$prop.sig.indicators == 1] <- "1"
indiv.indval.fill$prop_binary <- factor(indiv.indval.fill$prop_binary, 
                                        levels = c("1", "<1"))

indiv.indval.plot <- ggplot() +
  geom_col(data = indiv.indval.fill, 
           aes(x = factor(k), y = sum.sig.indicators, 
               fill = prop_binary)) +
  scale_fill_manual(values = c("#467f7b", "#cae6e4")) +
  xlab("Number of clusters") + ylab("Sum of significant indicator values") +
  labs(fill = "Proportion of clusters with\nsignificant indicator species") +
  scale_y_continuous(expand = c(0, 0)) +
  theme_minimal() +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        legend.position = "inside",
        legend.position.inside = c(0.84, 0.86),
        legend.title = element_text(size = 13, hjust = 1),
        legend.text = element_text(size = 12))
indiv.indval.plot

ggarrange(corrected.sil.plot, indiv.indval.plot, ncol = 1, labels = c("A", "B"))
ggsave("figures/winter/cluster-validation-winter.PNG", 
       width = 18.5, height = 20, units = "cm", bg = "white", dpi = 600)

## calculate IndVal at grid cell level, expect inflated values
grid.indval.fill <- data.frame(k = 2:10, sum.sig.indicators = NA, 
                               prop.sig.indicators = NA)

for (p in 2:10) {
  set.seed(256)
  p.indval <- indval(winter.grid.prop.t, cutree(grid.cluster, k = p))
  
  significant.indicators <- p.indval$indcls[p.indval$pval < 0.05]
  sum.sig.indicators <- sum(significant.indicators)
  
  # calculate prop with significant indicators
  clusters.with.sig.indval <- factor(p.indval$maxcls[p.indval$pval < 0.05])
  prop.sig.indicators <- length(levels(clusters.with.sig.indval)) / p
  
  grid.indval.fill$sum.sig.indicators[grid.indval.fill$k == p] <- sum.sig.indicators
  grid.indval.fill$prop.sig.indicators[grid.indval.fill$k == p] <- prop.sig.indicators
}

grid.indval.fill # no level have a significant indicator for each cluster, 5 has highest sum

# 6. Map clusters ----
grid.cluster5 <- data.frame(Cluster = cutree(grid.cluster, k = 5))
grid.cluster5$GridID <- rownames(grid.cluster5)

n10.grid <- subset(grid, GridID %in% winter.n10.grids)
n10.grid <- merge(x = n10.grid, y = grid.cluster5, by = "GridID",
                  all.x = FALSE, all.y = FALSE)
n10.grid$Cluster <- factor(n10.grid$Cluster)

n10.grid$WinterRegion[n10.grid$Cluster == 1] <- "Western JDF"
n10.grid$WinterRegion[n10.grid$Grid == 16] <- "Central JDF"
n10.grid$WinterRegion[n10.grid$GridID == 17] <- "Eastern JDF"
n10.grid$WinterRegion[n10.grid$GridID == 27] <- "Haro Strait"
n10.grid$WinterRegion[n10.grid$Cluster == 2 & !n10.grid$GridID == 16] <- "Strait of Georgia"
n10.grid$WinterRegion[n10.grid$GridID == 81] <- "Northern SOG"
n10.grid$WinterRegion[n10.grid$GridID == 56] <- "Howe Sound"

region.palette <- c(`Haro Strait` = "#D55E00", `Howe Sound` = "#009E73",
  `Western JDF` = "#deb887", `Eastern JDF` = "#976EDB", `Northern SOG` = "#5823AD",
  `Strait of Georgia` = "#0072B2", `Central JDF` = "#56b4e9")

n10.grid$WinterRegion <- factor(n10.grid$WinterRegion, 
  levels = c("Northern SOG", "Strait of Georgia", "Howe Sound",
             "Haro Strait", "Eastern JDF", "Central JDF", "Western JDF"))
winter.map <- ggplot() +
  geom_sf(data = bc.coastUTM) +
  geom_point(data = all.winter, aes(x = utm_x, y = utm_y), alpha = 0.3) +
  geom_sf(data = n10.grid, aes(fill = WinterRegion), alpha = 0.8) +

  geom_sf_text(data = n10.grid, aes(label = GridID), colour = "black", 
               size = 4.5) +
  scale_fill_manual(values = region.palette) +
  
  coord_sf(datum = 4326,
           crs = st_crs(32610),
           xlim = c(300*10^3, 540*10^3), ylim = c(5325*10^3, 5600*10^3)) +
  theme_bw() + 
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank(),
        panel.border = element_rect(colour = "black", linewidth = 0.8),
        legend.position = "inside",
        legend.position.inside = c(0.88, 0.88),
        legend.background = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(size = 11.5))
winter.map

## 7. Make dendrogram ----
df <- data.frame(x = 1:18, y = rep(-0.0085, 18),
  region = c("Haro Strait", "Howe Sound", "Western JDF", "Eastern JDF", 
  "Northern SOG", rep("Strait of Georgia", 2), "Central JDF", 
  rep("Strait of Georgia", 10)))

ggdendrogram(grid.cluster)
grid.cell.order <- grid.cluster$labels[grid.cluster$order]
grid.cell.order

winter.dendrogram <- ggdendrogram(grid.cluster) +
  geom_point(data = df, aes(x = x, y = y, colour = region),
             size = 7, shape = 15, show.legend = FALSE) +
  scale_colour_manual(values = region.palette) +
  
  annotate(geom = "rect", xmin = 0.55, xmax = 1.45, ymin = -0.045, ymax = -0.039,
           fill = "#D55E00") +
  annotate(geom = "rect", xmin = 1.55, xmax = 2.45, ymin = -0.045, ymax = -0.039,
           fill = "#009E73") +
  annotate(geom = "rect", xmin = 2.55, xmax = 3.45, ymin = -0.045, ymax = -0.039,
           fill = "#deb887") +
  annotate(geom = "rect", xmin = 3.55, xmax = 5.45, ymin = -0.045, ymax = -0.039,
           fill = "#5823AD") +
  annotate(geom = "rect", xmin = 5.55, xmax = 18.45, ymin = -0.045, ymax = -0.039,
           fill = "#0072B2") +
  annotate(geom = "text", x = 1, y = -0.06, label = "1", size = 4) +
  annotate(geom = "text", x = 2, y = -0.06, label = "2", size = 4) +
  annotate(geom = "text", x = 3, y = -0.06, label = "3", size = 4) +
  annotate(geom = "text", x = 4.5, y = -0.06, label = "4", size = 4) +
  annotate(geom = "text", x = 12, y = -0.06, label = "5", size = 4) +
  
  annotate(geom = "segment", x = 0.001, y = -0.02, yend = 0.52,
           colour = "black", linewidth = 1.2) +
  scale_x_continuous(limits = c(0, 19), expand = c(0, 0), breaks = 1:18,
                     labels = grid.cell.order) +
  
  theme(axis.title.y = element_text(size = 14),
        axis.title.x = element_blank(),
        axis.ticks.y.left = element_line(),
        axis.text.y = element_text(size = 11.5, angle = 0, colour = "black"),
        axis.text.x = element_text(size = 11.5, angle = 0, hjust = 0.5,
          vjust = 15.3, colour = "black", margin = margin(b = -20)),
        legend.title = element_blank(),
        legend.text = element_text(size = 11.5),
        legend.position = "inside", 
        legend.position.inside = c(0.83, 0.8),
        plot.background = element_blank(),
        plot.margin = margin(t = -5, b = 2)) +
  ylab("Bray-Curtis dissimilarity") 
winter.dendrogram

## 8. Stacked bar plot ----
winter.salish.ch <- merge(x = winter.salish.ch, y = grid.cluster5, by = "GridID",
                          all.x = FALSE, all.y = FALSE)
winter.salish.ch <- relocate(winter.salish.ch, GridID, .after = `Data flag`)
winter.salish.ch <- relocate(winter.salish.ch, Cluster, .after = GridID)

winter.salish.ch$WinterRegion[winter.salish.ch$Cluster == 1] <- "Western JDF"
winter.salish.ch$WinterRegion[winter.salish.ch$Grid == 16] <- "Central JDF"
winter.salish.ch$WinterRegion[winter.salish.ch$GridID == 17] <- "Eastern JDF"
winter.salish.ch$WinterRegion[winter.salish.ch$GridID == 27] <- "Haro Strait"
winter.salish.ch$WinterRegion[winter.salish.ch$Cluster == 2 &
                              !winter.salish.ch$GridID == 16] <- "Strait of Georgia"
winter.salish.ch$WinterRegion[winter.salish.ch$GridID == 81] <- "Northern SOG"
winter.salish.ch$WinterRegion[winter.salish.ch$GridID == 56] <- "Howe Sound"

winter.region.levels <- sort(unique(winter.salish.ch$WinterRegion))
names(winter.salish.ch)
ncol(winter.salish.ch[20:37])
region.prop <- data.frame(matrix(NA, nrow = 18, ncol = 0))

for (j in 1:length(winter.region.levels)) {
  region.subset <- subset(winter.salish.ch, WinterRegion == winter.region.levels[j])
  mean.diet <- colMeans(region.subset[20:37])
  region.prop <- cbind(region.prop, mean.diet)
  colnames(region.prop)[j] <- winter.region.levels[j]
}

# account for the rare prey removed:
colSums(region.prop)

rare.prey.adjust <- data.frame(t(data.frame(1 - colSums(region.prop))))
rownames(rare.prey.adjust) <- "Rare prey"
names(rare.prey.adjust) <- names(region.prop)
rare.prey.adjust

region.prop <- rbind(region.prop, rare.prey.adjust)
region.prop <- region.prop * 100 # 
colSums(region.prop) # good
region.prop$Group <- rownames(region.prop)

region.prop.long <- melt(region.prop)
names(region.prop.long) <- c("Group", "WinterRegion", "Percentage")

winter5 <- subset(region.prop.long, Percentage >= 5)
sort(unique(winter5$Group))

region.prop.long$PlotGroup <- region.prop.long$Group
region.prop.long$PlotGroup[region.prop.long$PlotGroup == "Pacific herring"] <- "Herring"
region.prop.long$PlotGroup[region.prop.long$PlotGroup == "Pacific sand lance"] <- "Sand lance"
region.prop.long$PlotGroup[region.prop.long$PlotGroup == "Northern anchovy"] <- "Anchovy"
region.prop.long$PlotGroup[region.prop.long$PlotGroup == "Gadiformes"] <- "Cod and hake"
region.prop.long$PlotGroup[region.prop.long$PlotGroup == "Planktonic shrimp"] <- "Plank. shrimp"
region.prop.long$PlotGroup[region.prop.long$PlotGroup == "Euphausiid"] <- "Krill"
region.prop.long$PlotGroup[region.prop.long$PlotGroup == "Threespine stickleback"] <- "Stickleback"

region.prop.long$PlotGroup[!region.prop.long$PlotGroup %in% c("Herring",
  "Sand lance", "Anchovy", "Surfperch", "Myctophid", "Mysid",
  "Krill", "Cod and hake", "Stickleback",  "Pandalid", "Plank. shrimp")] <- "Other"

region.prop.long$PlotGroup <- factor(region.prop.long$PlotGroup,
  levels = c("Herring", "Sand lance", "Anchovy", "Surfperch", "Myctophid", "Cod and hake",
             "Stickleback", "Krill", "Mysid", "Pandalid", "Plank. shrimp", "Other"))

region.prop.long$region_label[region.prop.long$WinterRegion == "Central JDF"] <- "Central\nJDF"
region.prop.long$region_label[region.prop.long$WinterRegion == "Eastern JDF"] <- "Eastern\nJDF"
region.prop.long$region_label[region.prop.long$WinterRegion == "Haro Strait"] <- "Haro\nStrait"
region.prop.long$region_label[region.prop.long$WinterRegion == "Howe Sound"] <- "Howe\nSound"
region.prop.long$region_label[region.prop.long$WinterRegion == "Northern SOG"] <- "Northern\nSOG"
region.prop.long$region_label[region.prop.long$WinterRegion == "Strait of Georgia"] <- "Strait of\nGeorgia"
region.prop.long$region_label[region.prop.long$WinterRegion == "Western JDF"] <- "Western\nJDF"

region.prop.long$region_label <- factor(region.prop.long$region_label,
  levels = c("Northern\nSOG", "Strait of\nGeorgia", "Howe\nSound", 
             "Haro\nStrait", "Eastern\nJDF", "Central\nJDF", "Western\nJDF"))

region.n <- data.frame(table(winter.salish.ch$WinterRegion))
names(region.n) <- c("WinterRegion", "n")
region.n$WinterRegion <- as.character(region.n$WinterRegion)

region.n$region_label[region.n$WinterRegion == "Central JDF"] <- "Central\nJDF"
region.n$region_label[region.n$WinterRegion == "Eastern JDF"] <- "Eastern\nJDF"
region.n$region_label[region.n$WinterRegion == "Haro Strait"] <- "Haro\nStrait"
region.n$region_label[region.n$WinterRegion == "Howe Sound"] <- "Howe\nSound"
region.n$region_label[region.n$WinterRegion == "Northern SOG"] <- "Northern\nSOG"
region.n$region_label[region.n$WinterRegion == "Strait of Georgia"] <- "Strait of\nGeorgia"
region.n$region_label[region.n$WinterRegion == "Western JDF"] <- "Western\nJDF"
region.n$region_label <- factor(region.n$region_label,
  levels = c("Northern\nSOG", "Strait of\nGeorgia", "Howe\nSound", 
             "Haro\nStrait", "Eastern\nJDF", "Central\nJDF", "Western\nJDF"))

winter.region.diet <- ggplot() +
  geom_bar(data = region.prop.long, 
           aes(x = region_label, weight = Percentage, fill = PlotGroup),
           colour = "black", linewidth = 0.2) +
  scale_fill_manual(values = c("#0072B2", "#D55E00", "#009E73", "#C05D94",
                                "#5823AD", "#897d1c", "#4f4800",  "#deb887", "#9ec0f7", 
                                "#4CA914", "#105764",
                               "grey60")) +
  geom_text(data = region.n, aes(x = region_label, y = 95, label = n),
            colour = "black") +
  coord_cartesian(expand = FALSE) +
  xlab("Cluster") + ylab("Percentage (%)") + 
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.title = element_blank(),
        axis.title.y = element_text(size = 14),
        axis.text.x = element_text(size = 11.5, colour = "black",
                                   margin = margin(b = -10)),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 11.5, colour = "black"),
        legend.text = element_text(size = 10.5),
        legend.key.size = unit(0.48, "cm"),
        legend.position = "bottom",
        panel.border = element_rect(linewidth = 1)) +
  guides(fill = guide_legend(nrow = 3, byrow = TRUE))
winter.region.diet

# group together
ggarrange(winter.map, 
          ggarrange(winter.dendrogram, winter.region.diet, ncol = 1, labels = c("B", "C"),
                    vjust = 1.95), 
          nrow = 1, widths = c(0.58, 0.42), labels = c("A"), hjust = -1, vjust = 1.95)
ggsave("figures/winter/winter-regions.PNG", width = 12, height = 8, units = "in",
       bg = "white", dpi = 600)


## 9. Principal coordinates analysis ----

par(mfrow = c(1, 1))
lingoes.pcoa <- pcoa(grid.bray, correction = "lingoes", rn = NULL)
biplot.pcoa(lingoes.pcoa)
biplot.pcoa(lingoes.pcoa, winter.grid.prop.t)
corrected.vectors <- data.frame(lingoes.pcoa$vectors.cor[, 1:2])
corrected.eigen <- lingoes.pcoa$values[2]
plot(lingoes.pcoa$values[,3])
r.sq <- (corrected.eigen[1, ] + corrected.eigen[2, ]) / (sum(corrected.eigen))
r.sq # 48.1

corrected.eigen[1,] / sum(corrected.eigen)
corrected.eigen[2,] / sum(corrected.eigen)

corrected.vectors$GridID <- rownames(corrected.vectors)
corrected.vectors <- merge(x = corrected.vectors, y = grid.cluster5,
                           by = "GridID", all.x = FALSE, all.y = FALSE)

corrected.vectors$WinterRegion[corrected.vectors$Cluster == 1] <- "Western JDF"
corrected.vectors$WinterRegion[corrected.vectors$Grid == 16] <- "Central JDF"
corrected.vectors$WinterRegion[corrected.vectors$GridID == 17] <- "Eastern JDF"
corrected.vectors$WinterRegion[corrected.vectors$GridID == 27] <- "Haro Strait"
corrected.vectors$WinterRegion[corrected.vectors$Cluster == 2 &
                                !corrected.vectors$GridID == 16] <- "Strait of Georgia"
corrected.vectors$WinterRegion[corrected.vectors$GridID == 81] <- "Northern SOG"
corrected.vectors$WinterRegion[corrected.vectors$GridID == 56] <- "Howe Sound"

corrected.vectors %>% group_by(GridID, WinterRegion) %>% tally() %>% view()

ggplot() + 
  geom_point(data = corrected.vectors, 
             aes(x = Axis.1, y = Axis.2, colour = WinterRegion), 
             size = 2.5, alpha = 0.8) +
  xlab("Axis 1 (26.8%)") + ylab("Axis 1 (21.3%)") +
  scale_colour_manual(values = region.palette) +
  theme_bw() + labs(colour = "Region") +
  theme(legend.position = "inside", 
        legend.position.inside = c(0.86, 0.26),
        legend.text = element_text(size = 11),
        legend.title = element_text(hjust = 0.5, face = "bold"),
        legend.background = element_rect(fill = "transparent", colour = "black"),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 11))
ggsave("figures/winter/winter-PCoA.PNG", width = 17, height = 12, units = "cm")

## 10. Indicator species analysis ----
names(winter.salish.ch)
sort(unique(winter.salish.ch$WinterRegion))
set.seed(427)
winter.region.indval <- indval(winter.salish.ch[20:37], 
                               winter.salish.ch$WinterRegion, numitr = 1000)
summary(winter.region.indval, type = "short")
summary(winter.region.indval, type = "long")
winter.region.indval 
# Herring is indicator species for SOG, slightly higher than central JDF
# note that krill are not an indicator of western JDF (0.145, rounds to 0.15)

## 11. Single species maps ----
winter.grid.prop.t$GridID <- rownames(winter.grid.prop.t)
winter.grid.prop.melt <- melt(winter.grid.prop.t, id.vars = c("GridID"))
names(winter.grid.prop.melt) <- c("GridID", "Group", "Proportion")

winter.grid.prop.major <- subset(winter.grid.prop.melt, Group %in% c("Pacific herring",
  "Pacific sand lance", "Northern anchovy", "Surfperch", "Myctophid", "Gadiformes",
  "Euphausiid", "Mysid", "Pandalid"))
winter.grid.prop.major$Group <- as.character(winter.grid.prop.major$Group)
winter.grid.prop.major$Group[winter.grid.prop.major$Group == "Gadiformes"] <- "Cod and hake"
winter.grid.prop.major$Group[winter.grid.prop.major$Group == "Euphausiid"] <- "Krill"

winter.grid.prop.major$Percentage <- winter.grid.prop.major$Proportion * 100
winter.grid.prop.major$Group <- factor(winter.grid.prop.major$Group,
  levels = c("Pacific herring", "Pacific sand lance", "Northern anchovy", "Surfperch",
             "Myctophid", "Cod and hake", "Mysid", "Krill", "Pandalid"))
winter.grid.prop.major <- merge(x = grid, y = winter.grid.prop.major, by = "GridID",
                                all.x = FALSE, all.y = TRUE)

ggplot() +
  geom_sf(data = bc.coastUTM) +
  geom_sf(data = winter.grid.prop.major, aes(fill = Percentage), alpha = 0.8) +
  scale_fill_distiller(palette = "PuBu", direction = 1, trans = "sqrt",
                       limits = c(0, 100)) +
  coord_sf(datum = 4326,
           crs = st_crs(32610),
           xlim = c(300*10^3, 540*10^3), ylim = c(5325*10^3, 5600*10^3)) +
  theme_bw() + 
  theme(axis.title = element_blank(),
        panel.grid = element_blank(),
        axis.text = element_blank(), axis.ticks = element_blank(),
        strip.background = element_blank(), strip.text = element_text(size = 11),
        legend.title = element_text(hjust = 0.5, margin = margin(b = 7), size = 11), 
        legend.text = element_text(size = 10)) +
  facet_wrap(~ Group, nrow = 3)
ggsave("figures/winter/winter-single-species.PNG", 
       width = 17, height = 18, units = "cm", dpi = 600)

## 12. October - February cluster analysis ----
early.winter <- subset(winter.salish.ch, Month %in% c(10, 11, 12, 1, 2))
nrow(early.winter) # 480

early.winter.grid.tally <- early.winter %>% group_by(GridID) %>% tally()
range(early.winter.grid.tally) # n = 5 to 81

early.winter <- select(early.winter, -Cluster) # remove because cluster added later
names(early.winter)
colSums(early.winter[19:36])
early.winter.grid.prop <- data.frame(matrix(NA, nrow = 18, ncol = 0))
length(winter.n10.grids) # 18

for (i in 1:length(winter.n10.grids)) {
  grid.subset <- subset(early.winter, GridID == winter.n10.grids[i])
  mean.diet <- colMeans(grid.subset[19:36])
  early.winter.grid.prop <- cbind(early.winter.grid.prop, mean.diet)
  colnames(early.winter.grid.prop)[i] <- winter.n10.grids[i]
}

colSums(early.winter.grid.prop) # expected

# calculate Bray-Curtis dissimilarity matrix
early.winter.grid.prop.t <- as.data.frame(t(early.winter.grid.prop))
early.winter.grid.bray <- vegdist(early.winter.grid.prop.t, method = "bray")

# hierarchical cluster analysis
early.winter.grid.cluster <- hclust(early.winter.grid.bray, "average")
plot(early.winter.grid.cluster, main = "", sub = "", xlab = "",
     ylab = "Bray-Curtis dissimilarity", hang = -1)

# cophenetic correlation - does the dendrogram reflect the dissimilarity matrix
early.winter.coph <- cophenetic(early.winter.grid.cluster)
cor(early.winter.grid.bray, early.winter.coph) # 0.93

# determine optimal number of clusters
# calculate IndVal at stomach level
early.winter.indval.fill <- data.frame(k = 2:10, sum.sig.indicators = NA, 
                                      prop.sig.indicators = NA)

for (z in 2:10) {
  z.clusters <- data.frame(Cluster = cutree(early.winter.grid.cluster, k = z))
  z.clusters$GridID <- rownames(z.clusters)
  
  indiv.diet <- merge(x = early.winter, y = z.clusters, by = "GridID",
                      all.x = FALSE, all.y = FALSE)
  nrow(indiv.diet)
  
  set.seed(530)
  z.indval <- indval(indiv.diet[19:36], indiv.diet$Cluster, numitr = 1000)
  
  # calculate sum of significant indicators
  significant.indicators <- z.indval$indcls[z.indval$pval < 0.05]
  sum.sig.indicators <- sum(significant.indicators)
  
  # calculate prop with sig indicators
  clusters.with.sig.indval <- factor(z.indval$maxcls[z.indval$pval < 0.05])
  prop.sig.indicators <- length(levels(clusters.with.sig.indval)) / z
  
  early.winter.indval.fill$sum.sig.indicators[early.winter.indval.fill$k == z] <- sum.sig.indicators
  early.winter.indval.fill$prop.sig.indicators[early.winter.indval.fill$k == z] <- prop.sig.indicators
}

early.winter.indval.fill
early.winter.indval.fill[which.max(early.winter.indval.fill$sum.sig.indicators), ]
# k = 5 is best

# make map
early.winter.cluster5 <- data.frame(Cluster = cutree(early.winter.grid.cluster, 
                                                    k = 5))
early.winter.cluster5$GridID <- rownames(early.winter.cluster5)

early.winter.n10.grid <- subset(grid, GridID %in% winter.n10.grids)
early.winter.n10.grid <- merge(x = early.winter.n10.grid, y = early.winter.cluster5,
                               by = "GridID", all.x = FALSE, all.y = FALSE)
early.winter.n10.grid$Cluster <- factor(early.winter.n10.grid$Cluster)

early.winter.cluster5

early.winter.n10.grid$EarlyWinterRegion[early.winter.n10.grid$Cluster == 1] <- "Western JDF"
early.winter.n10.grid$EarlyWinterRegion[early.winter.n10.grid$Grid == 16] <- "Central JDF"
early.winter.n10.grid$EarlyWinterRegion[early.winter.n10.grid$GridID == 17] <- "Eastern JDF"
early.winter.n10.grid$EarlyWinterRegion[early.winter.n10.grid$Cluster == 4] <- "Haro Strait"
early.winter.n10.grid$EarlyWinterRegion[early.winter.n10.grid$Cluster == 2 & !early.winter.n10.grid$GridID == 16] <- "Strait of Georgia"
early.winter.n10.grid$EarlyWinterRegion[early.winter.n10.grid$GridID == 81] <- "Northern SOG"
early.winter.n10.grid$EarlyWinterRegion[early.winter.n10.grid$Cluster == 5] <- "Howe Sound"

early.winter.n10.grid$EarlyWinterRegion <- factor(early.winter.n10.grid$EarlyWinterRegion, 
  levels = c("Northern SOG", "Strait of Georgia", "Howe Sound", "Haro Strait",
             "Eastern JDF", "Central JDF", "Western JDF"))
early.winter.map <- ggplot() +
  geom_sf(data = bc.coastUTM) +
  geom_sf(data = early.winter.n10.grid, aes(fill = EarlyWinterRegion), alpha = 0.8) +
  
  geom_sf_text(data = n10.grid, aes(label = GridID), colour = "black", 
               size = 4.5) +
  scale_fill_manual(values = region.palette) +
  
  coord_sf(datum = 4326,
           crs = st_crs(32610),
           xlim = c(300*10^3, 540*10^3), ylim = c(5325*10^3, 5600*10^3)) +
  theme_bw() + 
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank(),
        panel.border = element_rect(colour = "black", linewidth = 0.8),
        legend.position = "inside",
        legend.position.inside = c(0.88, 0.88),
        legend.background = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(size = 11.5))
early.winter.map

# make dendrogram
ggdendrogram(early.winter.grid.cluster)
early.winter.df <- data.frame(x = 1:18, y = rep(-0.01, 18),
  region = c("Howe Sound", "Howe Sound", "Western JDF", rep("Strait of Georgia", 8),
             "Central JDF", "Strait of Georgia", "Strait of Georgia", "Strait of Georgia",
             "Haro Strait", "Eastern JDF", "Northern SOG"))
early.winter.df$region <- factor(early.winter.df$region, levels = c("Haro Strait", "Howe Sound", 
  "Western JDF", "Eastern JDF", "Northern SOG", "Strait of Georgia",
   "Central JDF"))

early.winter.order <- early.winter.grid.cluster$labels[early.winter.grid.cluster$order]
early.winter.order

early.winter.dendrogram <- ggdendrogram(early.winter.grid.cluster) +
  geom_point(data = early.winter.df, aes(x = x, y = y, colour = region),
             size = 7, shape = 15, show.legend = FALSE) +
  scale_colour_manual(values = region.palette) +

  annotate(geom = "rect", xmin = 0.55, xmax = 2.45, ymin = -0.045, ymax = -0.039,
           fill = "#009E73") +
  annotate(geom = "rect", xmin = 2.55, xmax = 3.45, ymin = -0.045, ymax = -0.039,
           fill = "#deb887") +
  annotate(geom = "rect", xmin = 3.55, xmax = 15.45, ymin = -0.045, ymax = -0.039,
           fill = "#0072B2") +
  annotate(geom = "rect", xmin = 15.55, xmax = 16.45, ymin = -0.045, ymax = -0.039,
           fill = "#D55E00") +
  annotate(geom = "rect", xmin = 16.55, xmax = 18.45, ymin = -0.045, ymax = -0.039,
           fill = "#5823AD") +

  annotate(geom = "text", x = 1.5, y = -0.06, label = "1", size = 4) +
  annotate(geom = "text", x = 3, y = -0.06, label = "2", size = 4) +
  annotate(geom = "text", x = 9.5, y = -0.06, label = "3", size = 4) +
  annotate(geom = "text", x = 16, y = -0.06, label = "4", size = 4) +
  annotate(geom = "text", x = 17.5, y = -0.06, label = "5", size = 4) +

  annotate(geom = "segment", x = 0.001, y = -0.02, yend = 0.52,
           colour = "black", linewidth = 1.2) +
  scale_x_continuous(limits = c(0, 19), expand = c(0, 0), breaks = 1:18,
                     labels = early.winter.order) +
  
  theme(axis.title.y = element_text(size = 14),
        axis.title.x = element_blank(),
        axis.ticks.y.left = element_line(),
        axis.text.y = element_text(size = 11.5, angle = 0, colour = "black"),
        axis.text.x = element_text(size = 11.5, angle = 0, hjust = 0.5,
                                   vjust = 15.3, colour = "black", margin = margin(b = -20)),
        legend.title = element_blank(),
        legend.text = element_text(size = 11.5),
        legend.position = "inside", 
        legend.position.inside = c(0.83, 0.8),
        plot.background = element_blank(),
        plot.margin = margin(t = -5, b = 2)) +
  ylab("Bray-Curtis dissimilarity")
early.winter.dendrogram

# mean diet composition in each region
early.winter <- merge(x = early.winter, y = early.winter.cluster5, by = "GridID",
                      all.x = FALSE, all.y = FALSE)

early.winter$EarlyWinterRegion[early.winter$Cluster == 1] <- "Western JDF"
early.winter$EarlyWinterRegion[early.winter$Grid == 16] <- "Central JDF"
early.winter$EarlyWinterRegion[early.winter$GridID == 17] <- "Eastern JDF"
early.winter$EarlyWinterRegion[early.winter$Cluster == 4] <- "Haro Strait"
early.winter$EarlyWinterRegion[early.winter$Cluster == 2 & !early.winter$GridID == 16] <- "Strait of Georgia"
early.winter$EarlyWinterRegion[early.winter$GridID == 81] <- "Northern SOG"
early.winter$EarlyWinterRegion[early.winter$Cluster == 5] <- "Howe Sound"

early.winter.region.levels <- sort(unique(early.winter$EarlyWinterRegion))
early.winter.region.prop <- data.frame(matrix(NA, nrow = 18, ncol = 0))

for (j in 1:length(early.winter.region.levels)) {
  region.subset <- subset(early.winter, EarlyWinterRegion == early.winter.region.levels[j])
  mean.diet <- colMeans(region.subset[19:36])
  early.winter.region.prop <- cbind(early.winter.region.prop, mean.diet)
  colnames(early.winter.region.prop)[j] <- early.winter.region.levels[j]
}

# account for the rare prey removed
colSums(early.winter.region.prop)

early.winter.rare.prey.adjust <- data.frame(t(data.frame(1 - colSums(early.winter.region.prop))))
rownames(early.winter.rare.prey.adjust) <- "Rare prey"
names(early.winter.rare.prey.adjust) <- names(early.winter.region.prop)
early.winter.rare.prey.adjust

early.winter.region.prop <- rbind(early.winter.region.prop, early.winter.rare.prey.adjust)
early.winter.region.prop <- early.winter.region.prop * 100 
early.winter.region.prop$Group <- rownames(early.winter.region.prop)

early.winter.region.prop.long <- melt(early.winter.region.prop)
names(early.winter.region.prop.long) <- c("Group", "EarlyWinterRegion", "Percentage")

early.winter.region.prop.long$PlotGroup <- early.winter.region.prop.long$Group

early.winter.major <- subset(early.winter.region.prop.long, Percentage >= 5)
sort(unique(early.winter.major$PlotGroup))

early.winter.region.prop.long$PlotGroup[!early.winter.region.prop.long$PlotGroup %in%
 early.winter.major$PlotGroup] <- "Other"
early.winter.region.prop.long$PlotGroup[early.winter.region.prop.long$PlotGroup == "Pacific herring"] <- "Herring"
early.winter.region.prop.long$PlotGroup[early.winter.region.prop.long$PlotGroup == "Pacific sand lance"] <- "Sand lance"
early.winter.region.prop.long$PlotGroup[early.winter.region.prop.long$PlotGroup == "Northern anchovy"] <- "Anchovy"
early.winter.region.prop.long$PlotGroup[early.winter.region.prop.long$PlotGroup == "Euphausiid"] <- "Krill"
early.winter.region.prop.long$PlotGroup[early.winter.region.prop.long$PlotGroup == "Threespine stickleback"] <- "Stickleback"
early.winter.region.prop.long$PlotGroup[early.winter.region.prop.long$PlotGroup == "Planktonic shrimp"] <- "Plank. shrimp"
early.winter.region.prop.long$PlotGroup[early.winter.region.prop.long$PlotGroup == "Gadiformes"] <- "Cod and hake"
early.winter.region.prop.long$PlotGroup[early.winter.region.prop.long$PlotGroup == "Deep-sea smelt"] <- "Bathylagid"
sort(unique(early.winter.region.prop.long$PlotGroup))

early.winter.region.prop.long$PlotGroup <- factor(early.winter.region.prop.long$PlotGroup,
  levels = c("Herring", "Sand lance", "Anchovy", "Surfperch", "Myctophid", 
             "Cod and hake", "Stickleback", "Bathylagid", "Krill", "Mysid", "Pandalid", 
             "Plank. shrimp", "Other"))
levels(early.winter.region.prop.long$PlotGroup)

early.winter.region.prop.long$region_label[early.winter.region.prop.long$EarlyWinterRegion == "Central JDF"] <- "Central\nJDF"
early.winter.region.prop.long$region_label[early.winter.region.prop.long$EarlyWinterRegion == "Eastern JDF"] <- "Eastern\nJDF"
early.winter.region.prop.long$region_label[early.winter.region.prop.long$EarlyWinterRegion == "Haro Strait"] <- "Haro\nStrait"
early.winter.region.prop.long$region_label[early.winter.region.prop.long$EarlyWinterRegion == "Howe Sound"] <- "Howe\nSound"
early.winter.region.prop.long$region_label[early.winter.region.prop.long$EarlyWinterRegion == "Northern SOG"] <- "Northern\nSOG"
early.winter.region.prop.long$region_label[early.winter.region.prop.long$EarlyWinterRegion == "Strait of Georgia"] <- "Strait of\nGeorgia"
early.winter.region.prop.long$region_label[early.winter.region.prop.long$EarlyWinterRegion == "Western JDF"] <- "Western\nJDF"

early.winter.region.prop.long$region_label <- factor(early.winter.region.prop.long$region_label,
                                        levels = c("Northern\nSOG", "Strait of\nGeorgia", "Howe\nSound", 
                                                   "Haro\nStrait", "Eastern\nJDF", "Central\nJDF", "Western\nJDF"))

early.winter.n <- data.frame(table(early.winter$EarlyWinterRegion))
names(early.winter.n) <- c("EarlyWinterRegion", "n")
early.winter.n$region_label[early.winter.n$EarlyWinterRegion == "Central JDF"] <- "Central\nJDF"
early.winter.n$region_label[early.winter.n$EarlyWinterRegion == "Eastern JDF"] <- "Eastern\nJDF"
early.winter.n$region_label[early.winter.n$EarlyWinterRegion == "Haro Strait"] <- "Haro\nStrait"
early.winter.n$region_label[early.winter.n$EarlyWinterRegion == "Howe Sound"] <- "Howe\nSound"
early.winter.n$region_label[early.winter.n$EarlyWinterRegion == "Northern SOG"] <- "Northern\nSOG"
early.winter.n$region_label[early.winter.n$EarlyWinterRegion == "Strait of Georgia"] <- "Strait of\nGeorgia"
early.winter.n$region_label[early.winter.n$EarlyWinterRegion == "Western JDF"] <- "Western\nJDF"

early.winter.n

early.winter.region.diet <- ggplot() +
  geom_bar(data = early.winter.region.prop.long, 
           aes(x = region_label, weight = Percentage, fill = PlotGroup),
           colour = "black", linewidth = 0.2) +
  scale_fill_manual(values = c("#0072B2", "#D55E00", "#009E73", "#C05D94",
                               "#5823AD", "#897d1c", "#4f4800", "black",
                               "#deb887", "#9ec0f7", "#4CA914", "#105764",
                               "grey60")) +
  geom_text(data = early.winter.n, 
            aes(x = region_label, y = 95, label = n)) +
  coord_cartesian(expand = FALSE) +
  ylab("Percentage (%)") + 
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.title = element_blank(),
        axis.title.y = element_text(size = 14),
        axis.text.x = element_text(size = 11.5, colour = "black",
                                   margin = margin(b = -10)),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 11.5, colour = "black"),
        legend.text = element_text(size = 9),
        legend.position = "bottom",
        panel.border = element_rect(linewidth = 1),
        legend.key.size = unit(0.4, "cm")) +
  guides(fill = guide_legend(nrow = 4, byrow = TRUE))
early.winter.region.diet

# group together
ggarrange(early.winter.map, 
          ggarrange(early.winter.dendrogram, early.winter.region.diet, ncol = 1, labels = c("B", "C"),
                    vjust = 1.95), 
          nrow = 1, widths = c(0.58, 0.42), labels = c("A"), hjust = -1, vjust = 1.95)
ggsave("figures/winter/early-winter-regions.PNG", width = 12, height = 8, units = "in",
       bg = "white", dpi = 600)
