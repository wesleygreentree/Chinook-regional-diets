# Summer analysis

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
# species (p <= 0.05) and where the sum of significant indicator species among
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

# 12. Seasonal sensitivity analyses: determine if results robust to fishery
# closures and the designation of summer as April to September.

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
table(salish.ch$Season) # 1589 summer, 1066 
nrow(subset(salish.ch, TotalWeight == 0)) # 575 # includes stomachs at PFMA spatial level of resolution

# calculate number of sampling days for discussion and the number of days where
# samples could be collected in the Salish Sea
exact.date <- subset(salish.ch, `Date Resolution` == "exact date")
exact.date$date <- as.Date(paste(exact.date$Year, exact.date$Month, exact.date$Day, sep = "-"))
length(unique(paste(exact.date$Year, exact.date$Month, exact.date$Day))) # 803
 
study.period.days <- data.frame(date = seq(as.Date("2017-04-01"), 
                                           as.Date("2022-03-31"), by = "days"))
length(study.period.days$date) # 1826
803/1826 # 44.0% of days in the study period

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

# plot seasonal sampling
summer.months <- summer.salish.ch %>% group_by(Year, Month) %>% tally()
# note April 2019 closure started on April 19
# subset(summer.salish.ch, Month == 5 & Year == 2020)) # Becher Bay, where retention allowed

ggplot() +
  geom_col(data = summer.months, aes(x = Month, y = n, fill = Month),
           show.legend = FALSE) +
  scale_fill_viridis_c() +
  scale_x_continuous(breaks = c(4, 5, 6, 7, 8, 9),
    labels = c("April", "May", "June", "July", "August", "September")) +
  ylab("Number of stomachs") +
  theme_bw() +
  theme(axis.title = element_text(size = 12),
        axis.text = element_text(size = 11),
        strip.text = element_text(size = 12),
        strip.background = element_blank()) +
  facet_wrap(~ Year, ncol = 1, scales = "free_y")
ggsave("figures/summer/summer-month-n.PNG", 
       width = 17, height = 14, units = "cm", dpi = 600)

# 2. Process diet composition ----

# save all coordinates for maps
all.summer <- summer.salish.ch

# remove empty stomachs
nrow(summer.salish.ch) # 1479
nrow(subset(summer.salish.ch, TotalProp == 0)) # 357
summer.salish.ch <- subset(summer.salish.ch, TotalProp != 0)
summer.salish.ch <- summer.salish.ch[, !names(summer.salish.ch) == "Empty"]
nrow(summer.salish.ch) # 1122
357/1479 # 24.1%

names(summer.salish.ch) # 16:54
sum(colMeans(summer.salish.ch[16:54])) # 1 # good

# remove stomachs with only unidentified material
nrow(subset(summer.salish.ch, `Unidentified material` == 1)) # 129
summer.salish.ch <- subset(summer.salish.ch, `Unidentified material` != 1)
nrow(summer.salish.ch) # 993

# calculate % identified weight
names(summer.salish.ch)
summer.salish.ch$IdentifiedTotal <- rowSums(summer.salish.ch[16:53])

summer.salish.ch[16:53] <- summer.salish.ch[16:53] / summer.salish.ch$IdentifiedTotal
unique(rowSums(summer.salish.ch[16:53])) # 1 # good

# Select grids with n >= 10
summer.tally <- summer.salish.ch %>% group_by(GridID) %>% tally()
summer.tally <- subset(summer.tally, n >= 10) 
summer.n10.grids <- summer.tally$GridID # 23 grid cells
summer.n10.grids

# save data for species accumulation curves before removal of rare prey
# and grid cells with fewer than 10 stomachs
saveRDS(summer.salish.ch, "data/out/summer-Salish-Chinook-for-specaccum.RDS")

summer.salish.ch <- subset(summer.salish.ch, GridID %in% summer.n10.grids)
nrow(summer.salish.ch) # 968

## 3. Remove rare prey ----
summer.major.prey <- readRDS("data/out/summer-major-prey-vector.RDS") # from scripts/data-processing.R
names(summer.salish.ch)

# exclude rare prey for analysis
summer.salish.ch <- select(summer.salish.ch, `Fish Code`, Species, Season, 
  `Date Resolution`, Year, Month, Day, Region, PFMA, `Site Resolution`,
  Longitude, Latitude, utm_x, utm_y, `Final Length (cm)`, Project, `Data flag`,
  GridID, all_of(summer.major.prey))

names(summer.salish.ch)
major.prey.row.sum <- rowSums(summer.salish.ch[19:34])
hist(major.prey.row.sum)
view(major.prey.row.sum) # only 4 < 90% major prey
# I don't normalize % identified to sum to 1 after removing major prey

# save for size sensitivity analysis
saveRDS(summer.salish.ch, "data/out/summer-diet-for-size-sensitivity.RDS")

# 4. Cluster analysis ----

ncol(summer.salish.ch[19:34]) # 16
summer.grid.prop <- data.frame(matrix(NA, nrow = 16, ncol = 0)) # need to change

for (i in 1:length(summer.n10.grids)) {
  grid.subset <- subset(summer.salish.ch, GridID == summer.n10.grids[i])
  mean.diet <- colMeans(grid.subset[19:34])
  summer.grid.prop <- cbind(summer.grid.prop, mean.diet)
  colnames(summer.grid.prop)[i] <- summer.n10.grids[i]
}

colSums(summer.grid.prop) # expected to not all sum to 1 since rare prey removed

# calculate Bray-Curtis dissimilarity matrix
summer.grid.prop.t <- as.data.frame(t(summer.grid.prop))
grid.bray <- vegdist(summer.grid.prop.t, method = "bray")

# hierarchical cluster analysis
grid.cluster <- hclust(grid.bray, "average")
plot(grid.cluster, main = "", sub = "", xlab = "",
     ylab = "Bray-Curtis dissimilarity", hang = -1)

# cophenetic correlation - does the dendrogram reflect the dissimilarity matrix
coph <- cophenetic(grid.cluster)
cor(grid.bray, coph) # 0.96, good fit

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
# the calculation of indicator species analysis for each region 

# Potentially, one might argue that the indicator species should be calculated
# using average grid diet composition. I disagree, as this method is likely to
# over-inflate indicator values. Still, I wanted to check that the number of 
# clusters suggested is the same (line 319).

# I follow the approach of Borcard et al. 2018. Numerical Ecology with R.

# calculate IndVal at stomach level
indiv.indval.fill <- data.frame(k = 2:10, sum.sig.indicators = NA, 
                          prop.sig.indicators = NA)

for (z in 2:10) {
  z.clusters <- data.frame(Cluster = cutree(grid.cluster, k = z))
  z.clusters$GridID <- rownames(z.clusters)

  indiv.diet <- merge(x = summer.salish.ch, y = z.clusters, by = "GridID",
                      all.x = FALSE, all.y = FALSE)
  nrow(indiv.diet) # should always be 968

  set.seed(530)
  z.indval <- indval(indiv.diet[19:34], indiv.diet$Cluster, numitr = 1000)

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
ggarrange(corrected.sil.plot, indiv.indval.plot, ncol = 1, labels = c("A", "B"))
ggsave("figures/summer/cluster-validation-summer.PNG", 
       width = 18.5, height = 20, units = "cm", bg = "white", dpi = 600)

## calculate IndVal at grid cell level, expect inflated values
grid.indval.fill <- data.frame(k = 2:10, sum.sig.indicators = NA, 
                                prop.sig.indicators = NA)

for (p in 2:10) {
  set.seed(256)
  p.indval <- indval(summer.grid.prop.t, cutree(grid.cluster, k = p))
  
  significant.indicators <- p.indval$indcls[p.indval$pval < 0.05]
  sum.sig.indicators <- sum(significant.indicators) # 1.395
  
  # calculate prop with significant indicators
  clusters.with.sig.indval <- factor(p.indval$maxcls[p.indval$pval < 0.05])
  prop.sig.indicators <- length(levels(clusters.with.sig.indval)) / p
  
  grid.indval.fill$sum.sig.indicators[grid.indval.fill$k == p] <- sum.sig.indicators
  grid.indval.fill$prop.sig.indicators[grid.indval.fill$k == p] <- prop.sig.indicators
}

grid.indval.fill

# 6. Map clusters ----
grid.cluster3 <- data.frame(Cluster = cutree(grid.cluster, k = 3))
grid.cluster3$GridID <- rownames(grid.cluster3)

n10.grid <- subset(grid, GridID %in% summer.n10.grids)
n10.grid <- merge(x = n10.grid, y = grid.cluster3, by = "GridID",
                  all.x = FALSE, all.y = FALSE)
n10.grid$Cluster <- factor(n10.grid$Cluster)

n10.grid$Region[n10.grid$Cluster == 1 & 
                !n10.grid$GridID %in% c(15, 16, 23)] <- "Strait of Georgia"
n10.grid$Region[n10.grid$GridID %in% c(15, 16, 23)] <- "Juan de Fuca Strait"
n10.grid$Region[n10.grid$Cluster == 2] <- "Haro Strait"
n10.grid$Region[n10.grid$Cluster == 3] <- "Howe Sound"

n10.grid$Region <- factor(n10.grid$Region, levels = c("Strait of Georgia",
  "Howe Sound", "Haro Strait", "Juan de Fuca Strait"))

n10.grid$GridID[n10.grid$GridID == 103] <- 99 # to improve readability of plots

summer.map <- ggplot() +
  geom_sf(data = bc.coastUTM) +
  geom_point(data = all.summer, aes(x = utm_x, y = utm_y), alpha = 0.3) +
  geom_sf(data = n10.grid, aes(fill = Region), alpha = 0.8,) +
  geom_sf_text(data = n10.grid, aes(label = GridID), colour = "black", size = 4.5) +
  scale_fill_manual(values = c("#0072B2", "#009E73", "#D55E00",  "#56b4e9")) +
  
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
        legend.position.inside = c(0.86, 0.93),
        legend.background = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(size = 11.5))
summer.map

## 7. Make dendrogram ----
df <- data.frame(x = 1:23, y = rep(-0.018, 23), 
                 region = c("Haro Strait", "Haro Strait", "Howe Sound",
                            rep("Strait of Georgia", 20)))
df$region[df$x %in% c(4, 5, 6)] <- "Juan de Fuca Strait"

ggdendrogram(grid.cluster)
grid.cell.order <- grid.cluster$labels[grid.cluster$order]
grid.cell.order[which(grid.cell.order == "103")] <- "99"
grid.cell.order

summer.dendrogram <- ggdendrogram(grid.cluster) +
  geom_point(data = df, aes(x = x, y = y, colour = region), 
             size = 5.5, shape = 15, show.legend = FALSE) +
  guides(colour = guide_legend(override.aes = list(size = 4))) +
  scale_colour_manual(values = c("#D55E00", "#009E73", "#56b4e9", "#0072B2")) +
  
  annotate(geom = "rect", xmin = 0.55, xmax = 2.45, ymin = -0.049, ymax = -0.043,
           fill = "#D55E00") +
  annotate(geom = "rect", xmin = 2.55, xmax = 3.45, ymin = -0.049, ymax = -0.043,
           fill = "#009E73") +
  annotate(geom = "rect", xmin = 3.55, xmax = 23.45, ymin = -0.049, ymax = -0.043,
           fill = "#0072B2") +
 
  annotate(geom = "text", x = 1.5, y = -0.063, label = "1", size = 4) +
  annotate(geom = "text", x = 3, y = -0.063, label = "2", size = 4) +
  annotate(geom = "text", x = 13.5, y = -0.063, label = "3", size = 4) +
  
  annotate(geom = "segment", x = 0.001, y = -0.02, yend = 0.64, colour = "black",
           linewidth = 1.2) +
  scale_x_continuous(limits = c(0, 24), expand = c(0, 0), breaks = c(1:23),
                     labels = grid.cell.order) +
  theme(axis.title.y = element_text(size = 14),
        axis.ticks.y.left = element_line(),
        axis.text.y = element_text(size = 11.5, angle = 0, colour = "black"),
        axis.text.x = element_text(size = 11, angle = 0, hjust = 0.5, 
                                   vjust = 15.3, colour = "black", margin = margin(b = -29)),
        plot.background = element_blank(),
        plot.margin = margin(t = -5, b = 2)) +
  ylab("Bray-Curtis dissimilarity")
summer.dendrogram

## 8. Stacked bar plot ----
summer.salish.ch <- merge(x = summer.salish.ch, y = grid.cluster3, by = "GridID",
                          all.x = FALSE, all.y = FALSE)
summer.salish.ch <- relocate(summer.salish.ch, GridID, .after = `Data flag`)
summer.salish.ch <- relocate(summer.salish.ch, Cluster, .after = GridID)

summer.salish.ch$SummerRegion[summer.salish.ch$Cluster == 1 & 
                  !summer.salish.ch$GridID %in% c(15, 16, 23)] <- "Strait of Georgia"
summer.salish.ch$SummerRegion[summer.salish.ch$GridID %in% c(15, 16, 23)] <- "Juan de Fuca Strait"
summer.salish.ch$SummerRegion[summer.salish.ch$Cluster == 2] <- "Haro Strait"
summer.salish.ch$SummerRegion[summer.salish.ch$Cluster == 3] <- "Howe Sound"
summer.salish.ch <- relocate(summer.salish.ch, SummerRegion, .after = Cluster)

summer.region.levels <- sort(unique(summer.salish.ch$SummerRegion))
names(summer.salish.ch)
region.prop <- data.frame(matrix(NA, nrow = 16, ncol = 0))

for (t in 1:length(summer.region.levels)) {
  region.subset <- subset(summer.salish.ch, SummerRegion == summer.region.levels[t])
  mean.diet <- colMeans(region.subset[21:36])
  region.prop <- cbind(region.prop, mean.diet)
  colnames(region.prop)[t] <- summer.region.levels[t]
}

# account for the rare prey removed:
colSums(region.prop)

rare.prey.adjust <- data.frame(t(data.frame(1 - colSums(region.prop))))
rownames(rare.prey.adjust) <- "Rare prey"
names(rare.prey.adjust) <- names(region.prop)
rare.prey.adjust

region.prop <- rbind(region.prop, rare.prey.adjust)
region.prop <- region.prop * 100
colSums(region.prop) # good
region.prop$Group <- rownames(region.prop)

region.prop.long <- melt(region.prop)
names(region.prop.long) <- c("Group", "SummerRegion", "Percentage")

region.prop.long$PlotGroup <- region.prop.long$Group
region.prop.long$PlotGroup[region.prop.long$PlotGroup == "Pacific herring"] <- "Herring"
region.prop.long$PlotGroup[region.prop.long$PlotGroup == "Pacific sand lance"] <- "Sand lance"
region.prop.long$PlotGroup[region.prop.long$PlotGroup == "Northern anchovy"] <- "Anchovy"
region.prop.long$PlotGroup[region.prop.long$PlotGroup == "Gadiformes"] <- "Cod and hake"
region.prop.long$PlotGroup[region.prop.long$PlotGroup == "Crustacean larvae"] <- "Crust. larvae"

region.prop.long$PlotGroup[!region.prop.long$PlotGroup %in% c("Herring",
  "Sand lance", "Anchovy", "Crust. larvae", "Surfperch",
  "Cod and hake")] <- "Other"
region.prop.long$PlotGroup <- factor(region.prop.long$PlotGroup,
  levels = c("Herring", "Sand lance", "Anchovy", "Cod and hake", "Surfperch",
             "Crust. larvae", "Other"))

region.prop.long$SummerRegion <- as.character(region.prop.long$SummerRegion)
region.prop.long$SummerRegion[region.prop.long$SummerRegion == "Strait of Georgia"] <- "Strait of\nGeorgia"
region.prop.long$SummerRegion[region.prop.long$SummerRegion == "Juan de Fuca Strait"] <- "Juan de Fuca\nStrait"
region.prop.long$SummerRegion <- factor(region.prop.long$SummerRegion,
      levels = c("Strait of\nGeorgia", "Howe Sound", "Haro Strait", 
                 "Juan de Fuca\nStrait"))

region.n <- data.frame(table(summer.salish.ch$SummerRegion))
names(region.n) <- c("SummerRegion", "n")
region.n$SummerRegion <- as.character(region.n$SummerRegion)
region.n$SummerRegion[region.n$SummerRegion == "Strait of Georgia"] <- "Strait of\nGeorgia"
region.n$SummerRegion[region.n$SummerRegion == "Juan de Fuca Strait"] <- "Juan de Fuca\nStrait"
region.n$SummerRegion <- factor(region.n$SummerRegion,
  levels = c("Strait of\nGeorgia", "Howe Sound", "Haro Strait", 
             "Juan de Fuca\nStrait"))

summer.region.diet <- ggplot() +
  geom_bar(data = region.prop.long, 
           aes(x = SummerRegion, weight = Percentage, fill = PlotGroup),
           colour = "black", linewidth = 0.2) +
  scale_fill_manual(values = c("#0072B2", "#D55E00", "#009E73", "#897d1c",
                               "#C05D94", "#daa21b", "grey60")) +
  geom_text(data = region.n,
            aes(x = SummerRegion, y = 95, label = n), colour = "black") +
  
  coord_cartesian(expand = FALSE) +
  ylab("Percentage (%)") + labs(fill = "Prey Category") +
  theme_bw() +
  theme(legend.title = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 14),
        axis.text.x = element_text(size = 11.5, colour = "black",
                                   margin = margin(b = -10)),
        axis.text.y = element_text(size = 11.5, colour = "black"),
        legend.text = element_text(size = 11),
        legend.position = "bottom",
        legend.key.size = unit(0.52, "cm"),
        panel.border = element_rect(linewidth = 1))
summer.region.diet

ggarrange(summer.map, 
          ggarrange(summer.dendrogram, summer.region.diet, ncol = 1, labels = c("B", "C"),
                    vjust = 1.95), 
          nrow = 1, widths = c(0.58, 0.42), labels = c("A"), hjust = -1, vjust = 1.95)
ggsave("figures/summer/summer-regions.PNG", width = 12, height = 8, units = "in",
       bg = "white", dpi = 600)


## 9. Principal coordinate analysis (PCoA) ----

# PCOA
par(mfrow = c(1, 1))
lingoes.pcoa <- pcoa(grid.bray, correction = "lingoes", rn = NULL)
biplot.pcoa(lingoes.pcoa)
biplot.pcoa(lingoes.pcoa, summer.grid.prop.t)
corrected.vectors <- data.frame(lingoes.pcoa$vectors.cor[, 1:2])
corrected.eigen <- lingoes.pcoa$values[2]
plot(lingoes.pcoa$values[,3])
r.sq <- (corrected.eigen[1, ] + corrected.eigen[2, ]) / (sum(corrected.eigen))
r.sq # 57.9

corrected.eigen[1,] / sum(corrected.eigen)
corrected.eigen[2,] / sum(corrected.eigen)

corrected.vectors$GridID <- rownames(corrected.vectors)

corrected.vectors <- merge(x = corrected.vectors, y = grid.cluster3,
                           by = "GridID", all.x = FALSE, all.y = FALSE)

corrected.vectors$SummerRegion[corrected.vectors$Cluster == 1 & 
                                !corrected.vectors$GridID %in% c(15, 16, 23)] <- "Strait of Georgia"
corrected.vectors$SummerRegion[corrected.vectors$GridID %in% c(15, 16, 23)] <- "Juan de Fuca Strait"
corrected.vectors$SummerRegion[corrected.vectors$Cluster == 2] <- "Haro Strait"
corrected.vectors$SummerRegion[corrected.vectors$Cluster == 3] <- "Howe Sound"
ggplot() + 
  geom_point(data = corrected.vectors, 
             aes(x = Axis.1, y = Axis.2, colour = SummerRegion), 
             size = 2.5, alpha = 0.8) +
  xlab("Axis 1 (38.4%)") + ylab("Axis 1 (19.6%)") +
  scale_colour_manual(values = c("#D55E00", "#009E73", "#56b4e9", "#0072B2")) +
  theme_bw() + labs(colour = "Region") +
  theme(legend.position = "inside", 
        legend.position.inside = c(0.16, 0.82),
        legend.text = element_text(size = 11),
        legend.title = element_text(hjust = 0.5, face = "bold"),
        legend.background = element_rect(fill = "transparent", colour = "black"),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 11))
ggsave("figures/summer/summer-PCoA.PNG", width = 17, height = 12, units = "cm")

# 10. Indicator species ----
set.seed(267)
summer.region.indval <- indval(summer.salish.ch[21:36], 
                               summer.salish.ch$SummerRegion, numitr = 1000)
summary(summer.region.indval, type = "short")
summary(summer.region.indval, type = "long")
summer.region.indval

## 11. Single species maps ----

summer.grid.prop.t$GridID <- rownames(summer.grid.prop.t)
summer.grid.prop.melt <- melt(summer.grid.prop.t, id.vars = c("GridID"))
names(summer.grid.prop.melt) <- c("GridID", "Group", "Proportion")

prop05 <- subset(region.prop.long, Percentage >= 5)
prop05.groups <- sort(unique(prop05$Group))
prop05.groups

summer.grid.prop.major <- subset(summer.grid.prop.melt, Group %in% prop05.groups)
summer.grid.prop.major$Percentage <- summer.grid.prop.major$Proportion * 100
summer.grid.prop.major$Group <- as.character(summer.grid.prop.major$Group)
summer.grid.prop.major$Group[summer.grid.prop.major$Group == "Gadiformes"] <- "Cod and hake"

summer.grid.prop.major$Group <- factor(summer.grid.prop.major$Group,
  levels = c("Pacific herring", "Pacific sand lance", "Northern anchovy",
             "Surfperch", "Cod and hake", "Crustacean larvae"))
summer.grid.prop.major <- merge(x = grid, y = summer.grid.prop.major, by = "GridID",
                                all.x = FALSE, all.y = TRUE)
class(summer.grid.prop.major) # sf

ggplot() +
  geom_sf(data = bc.coastUTM) +
  geom_sf(data = summer.grid.prop.major, aes(fill = Percentage), alpha = 0.8) +
  scale_fill_distiller(palette = "OrRd", direction = 1, trans = "sqrt",
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
  facet_wrap(~ Group, nrow = 2)
ggsave("figures/summer/summer-single-species.PNG", 
       width = 17, height = 12.5, units = "cm", dpi = 600)

## 12. Seasonal sensitivity analyses ----

## 12A. Summer 2017-2018 cluster analysis ----
summer2017.18 <- subset(summer.salish.ch, Year %in% c(2017, 2018))
nrow(summer2017.18) # 439
summer2017.18 <- select(summer2017.18, -Cluster)

summer2017.18.grid.tally <- summer2017.18 %>% group_by(GridID) %>% tally()
summer2017.18.grid.tally7 <- subset(summer2017.18.grid.tally, n >= 7)

summer2017.18 <- subset(summer2017.18, GridID %in% summer2017.18.grid.tally7$GridID)
nrow(summer2017.18) # 436

names(summer2017.18)
colSums(summer2017.18[20:35]) # remove rockfish and sculpin, never sampled in 2017, 2018
nrow(subset(summer2017.18, Rockfish > 0 | Sculpin > 0)) # 0
summer2017.18 <- select(summer2017.18, -Rockfish, -Sculpin)

summer2017.18.grid.prop <- data.frame(matrix(NA, nrow = 14, ncol = 0))

for (i in 1:length(summer2017.18.grid.tally7$GridID)) {
  grid.subset <- subset(summer2017.18, GridID == summer2017.18.grid.tally7$GridID[i])
  mean.diet <- colMeans(grid.subset[20:33])
  summer2017.18.grid.prop <- cbind(summer2017.18.grid.prop, mean.diet)
  colnames(summer2017.18.grid.prop)[i] <- summer2017.18.grid.tally7$GridID[i]
}
# grid cell 72 only has herring (Comox)
colSums(summer2017.18.grid.prop) # expected

# calculate Bray-Curtis dissimilarity matrix
summer2017.18.grid.prop.t <- as.data.frame(t(summer2017.18.grid.prop))
summer2017.18.grid.bray <- vegdist(summer2017.18.grid.prop.t, method = "bray")

# hierarchical cluster analysis
summer2017.18.grid.cluster <- hclust(summer2017.18.grid.bray, "average")
plot(summer2017.18.grid.cluster, main = "", sub = "", xlab = "",
     ylab = "Bray-Curtis dissimilarity", hang = -1)

summer2017.18.coph <- cophenetic(summer2017.18.grid.cluster)
cor(summer2017.18.grid.bray, summer2017.18.coph) # 0.98

# determine optimal number of clusters
# calculate IndVal at stomach level
summer2017.18.indval.fill <- data.frame(k = 2:10, sum.sig.indicators = NA, 
                                prop.sig.indicators = NA)

for (z in 2:10) {
  z.clusters <- data.frame(Cluster = cutree(summer2017.18.grid.cluster, k = z))
  z.clusters$GridID <- rownames(z.clusters)
  
  indiv.diet <- merge(x = summer2017.18, y = z.clusters, by = "GridID",
                      all.x = FALSE, all.y = FALSE)
  
  set.seed(530)
  z.indval <- indval(indiv.diet[20:33], indiv.diet$Cluster, numitr = 1000)
  
  # calculate sum of significant indicators
  significant.indicators <- z.indval$indcls[z.indval$pval < 0.05]
  sum.sig.indicators <- sum(significant.indicators)
  
  # calculate prop with sig indicators
  clusters.with.sig.indval <- factor(z.indval$maxcls[z.indval$pval < 0.05])
  prop.sig.indicators <- length(levels(clusters.with.sig.indval)) / z
  
  summer2017.18.indval.fill$sum.sig.indicators[summer2017.18.indval.fill$k == z] <- sum.sig.indicators
  summer2017.18.indval.fill$prop.sig.indicators[summer2017.18.indval.fill$k == z] <- prop.sig.indicators
}

summer2017.18.indval.fill

summer2017.18.cluster4 <- data.frame(Cluster = cutree(summer2017.18.grid.cluster, k = 4))
summer2017.18.cluster4$GridID <- rownames(summer2017.18.cluster4)

# map:
summer2017.18.n10.grid <- subset(grid, GridID %in% summer2017.18.grid.tally7$GridID)
summer2017.18.n10.grid <- merge(x = summer2017.18.n10.grid, y = summer2017.18.cluster4,
                                by = "GridID", all.x = FALSE, all.y = FALSE)
summer2017.18.n10.grid$Cluster <- factor(summer2017.18.n10.grid$Cluster)

summer2017.18.n10.grid$Region[summer2017.18.n10.grid$Cluster == 1] <- "Strait of Georgia"
summer2017.18.n10.grid$Region[summer2017.18.n10.grid$Cluster == 1 &
  summer2017.18.n10.grid$GridID %in% c(15, 16)] <- "Juan de Fuca Strait"
summer2017.18.n10.grid$Region[summer2017.18.n10.grid$Cluster == 2] <- "Haro Strait"
summer2017.18.n10.grid$Region[summer2017.18.n10.grid$Cluster == 3] <- "Haro Strait transition"
summer2017.18.n10.grid$Region[summer2017.18.n10.grid$Cluster == 4] <- "Howe Sound"

summer2017.18.n10.grid$Region <- factor(summer2017.18.n10.grid$Region, levels = c("Strait of Georgia",
  "Howe Sound", "Haro Strait transition", "Haro Strait", "Juan de Fuca Strait"))

summer2017.map <- ggplot() +
  geom_sf(data = bc.coastUTM) +
  geom_sf(data = summer2017.18.n10.grid, aes(fill = Region), alpha = 0.8,) +
  geom_sf_text(data = summer2017.18.n10.grid, aes(label = GridID), colour = "black", size = 4.5) +
  scale_fill_manual(values = c("#0072B2", "#009E73",  "#758000", "#D55E00",  "#56b4e9")) +
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
        legend.position.inside = c(0.82, 0.91),
        legend.background = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(size = 11.5))
summer2017.map

# dendrogram
df2017.18 <- data.frame(x = 1:18, y = rep(-0.025, 18),
  region = c("Haro Strait", "Haro Strait", "Howe Sound", "Haro Strait transition",
             rep("Strait of Georgia", 5), "Juan de Fuca Strait", "Juan de Fuca Strait",
             rep("Strait of Georgia", 7)))

ggdendrogram(summer2017.18.grid.cluster)
order2017.18 <- summer2017.18.grid.cluster$labels[summer2017.18.grid.cluster$order]
order2017.18

summer2017.dendrogram <- ggdendrogram(summer2017.18.grid.cluster) +
  geom_point(data = df2017.18, aes(x = x, y = y, colour = region),
             size = 5.5, shape = 15, show.legend = FALSE) +
  scale_colour_manual(values = c("#D55E00", "#758000", "#009E73", "#56b4e9", "#0072B2")) +
  annotate(geom = "segment", x = 0.001, y = -0.02, yend = 0.75, colour = "black",
           linewidth = 0.9) +
  annotate(geom = "rect", xmin = 0.55, xmax = 2.45, ymin = -0.055, ymax = -0.049,
           fill = "#D55E00") +
  annotate(geom = "rect", xmin = 2.55, xmax = 3.45, ymin = -0.055, ymax = -0.049,
           fill = "#009E73") +
  annotate(geom = "rect", xmin = 3.55, xmax = 4.45, ymin = -0.055, ymax = -0.049,
           fill = "#758000") +
  annotate(geom = "rect", xmin = 4.55, xmax = 18.45, ymin = -0.055, ymax = -0.049,
           fill = "#0072B2") +
  annotate("text", x = 1.5, y = -0.07, label = "1", size = 4) +
  annotate("text", x = 3, y = -0.07, label = "2", size = 4) +
  annotate("text", x = 4, y = -0.07, label = "3", size = 4) +
  annotate("text", x = 11.5, y = -0.07, label = "4", size = 4) +
  
  scale_x_continuous(limits = c(0, 19), expand = c(0, 0), breaks = 1:18,
                     labels = order2017.18) +
  scale_y_continuous(breaks = c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7)) +
  theme(axis.title.y = element_text(size = 14),
        axis.ticks.y.left = element_line(),
        axis.text.y = element_text(size = 11.5, angle = 0, colour = "black"),
        axis.text.x = element_text(size = 11, angle = 0, hjust = 0.5, 
                                   vjust = 14.2, colour = "black", margin = margin(b = -29)),
        plot.background = element_blank(),
        plot.margin = margin(t = -5, b = 2)) +
  ylab("Bray-Curtis dissimilarity")
summer2017.dendrogram

# calculate mean diet in the regions
summer2017.18 <- merge(x = summer2017.18, y = summer2017.18.cluster4, by = "GridID",
                          all.x = FALSE, all.y = FALSE)

summer2017.18$SummerRegion2017[summer2017.18$Cluster == 1] <- "Strait of Georgia"
summer2017.18$SummerRegion2017[summer2017.18$GridID %in% c(15, 16)] <- "Juan de Fuca Strait"
summer2017.18$SummerRegion2017[summer2017.18$Cluster == 2] <- "Haro Strait"
summer2017.18$SummerRegion2017[summer2017.18$Cluster == 3] <- "Haro Strait transition"
summer2017.18$SummerRegion2017[summer2017.18$Cluster == 4] <- "Howe Sound"

nrow(subset(summer2017.18, is.na(SummerRegion2017))) # 0

summer2017.region.levels <- sort(unique(summer2017.18$SummerRegion2017))
names(summer2017.18)
ncol(summer2017.18[20:33])
region.prop2017 <- data.frame(matrix(NA, nrow = 14, ncol = 0))

for (t in 1:length(summer2017.region.levels)) {
  region.subset <- subset(summer2017.18, SummerRegion2017 == summer2017.region.levels[t])
  mean.diet <- colMeans(region.subset[20:33])
  region.prop2017 <- cbind(region.prop2017, mean.diet)
  colnames(region.prop2017)[t] <- summer2017.region.levels[t]
}

# account for the rare prey removed (Juan de Fuca and Strait of Georgia only)
colSums(region.prop2017)

rare.prey.adjust2017 <- data.frame(t(data.frame(1 - colSums(region.prop2017))))
rownames(rare.prey.adjust2017) <- "Rare prey"
names(rare.prey.adjust2017) <- names(region.prop2017)
rare.prey.adjust2017

region.prop2017 <- rbind(region.prop2017, rare.prey.adjust2017)
region.prop2017 <- region.prop2017 * 100
colSums(region.prop2017) # all sum to 100
region.prop2017$Group <- rownames(region.prop2017)

region.prop.long2017 <- melt(region.prop2017)
names(region.prop.long2017) <- c("Group", "SummerRegion2017", "Percentage")

prop05.2017 <- subset(region.prop.long2017, Percentage >= 5)
unique(prop05.2017$Group)

unique(region.prop.long2017$PlotGroup)
region.prop.long2017$PlotGroup <- as.character(region.prop.long2017$Group)
region.prop.long2017$PlotGroup[region.prop.long2017$PlotGroup == "Threespine stickleback"] <- "Stickleback"
region.prop.long2017$PlotGroup[region.prop.long2017$PlotGroup == "Pacific herring"] <- "Herring"
region.prop.long2017$PlotGroup[region.prop.long2017$PlotGroup == "Pacific sand lance"] <- "Sand lance"
region.prop.long2017$PlotGroup[region.prop.long2017$PlotGroup == "Northern anchovy"] <- "Anchovy"
region.prop.long2017$PlotGroup[region.prop.long2017$PlotGroup == "Gadiformes"] <- "Cod and hake"
region.prop.long2017$PlotGroup[region.prop.long2017$PlotGroup == "Crustacean larvae"] <- "Crust. larvae"
sort(unique(region.prop.long2017$PlotGroup))

region.prop.long2017$PlotGroup[!region.prop.long2017$PlotGroup %in% c("Herring",
  "Sand lance", "Anchovy", "Crust. larvae", "Surfperch", "Cod and hake", "Pandalid",
  "Stickleback")] <- "Other"
region.prop.long2017$PlotGroup <- factor(region.prop.long2017$PlotGroup,
  levels = c("Herring", "Sand lance", "Anchovy", "Surfperch", "Cod and hake",
            "Stickleback",  "Crust. larvae", "Pandalid", "Other"))

region.prop.long2017$SummerRegion2017 <- as.character(region.prop.long2017$SummerRegion2017)
region.prop.long2017$SummerRegion2017[region.prop.long2017$SummerRegion2017 == "Strait of Georgia"] <- "Strait of\nGeorgia"
region.prop.long2017$SummerRegion2017[region.prop.long2017$SummerRegion2017 == "Juan de Fuca Strait"] <- "Juan de Fuca\nStrait"
region.prop.long2017$SummerRegion2017[region.prop.long2017$SummerRegion2017 == "Haro Strait transition"] <- "Haro Strait\ntransition"
region.prop.long2017$SummerRegion2017 <- factor(region.prop.long2017$SummerRegion2017,
  levels = c("Strait of\nGeorgia", "Howe Sound", "Haro Strait", "Haro Strait\ntransition",
             "Juan de Fuca\nStrait"))

region.n2017 <- data.frame(table(summer2017.18$SummerRegion2017))
names(region.n2017) <- c("SummerRegion2017", "n")
region.n2017$SummerRegion2017 <- as.character(region.n2017$SummerRegion2017)

region.n2017$SummerRegion2017[region.n2017$SummerRegion2017 == "Strait of Georgia"] <- "Strait of\nGeorgia"
region.n2017$SummerRegion2017[region.n2017$SummerRegion2017 == "Juan de Fuca Strait"] <- "Juan de Fuca\nStrait"
region.n2017$SummerRegion2017[region.n2017$SummerRegion2017 == "Haro Strait transition"] <- "Haro Strait\ntransition"
region.n2017$SummerRegion2017[region.n2017$SummerRegion2017 == "Haro Strait"] <- "Haro Strait"
region.n2017$SummerRegion2017[region.n2017$SummerRegion2017 == "Howe Sound"] <- "Howe Sound"
region.n2017$SummerRegion2017 <- factor(region.n2017$SummerRegion2017,
  levels = c("Strait of\nGeorgia", "Howe Sound", "Haro Strait", "Haro Strait\ntransition",
            "Juan de Fuca\nStrait"))
region.n2017

summer.region.diet2017 <- ggplot() +
  geom_bar(data = region.prop.long2017, 
           aes(x = SummerRegion2017, weight = Percentage, fill = PlotGroup),
           colour = "black", linewidth = 0.2) +
  scale_fill_manual(values = c("#0072B2", "#D55E00", "#009E73", "#C05D94", "#897d1c",
                               "#4f4800", "#daa21b", "#4CA914", "grey60")) +
  geom_text(data = region.n2017, aes(x = SummerRegion2017, y = 95, label = n)) +
  coord_cartesian(expand = FALSE) +
  ylab("Percentage (%)") + labs(fill = "Prey Category") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.title = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 14),
        axis.text.x = element_text(size = 11.5, colour = "black",
                                   margin = margin(b = -10)),
        axis.text.y = element_text(size = 11.5, colour = "black"),
        legend.text = element_text(size = 11),
        legend.key.size = unit(0.52, "cm"),
        legend.position = "bottom",
        panel.border = element_rect(linewidth = 1)) +
  guides(fill = guide_legend(nrow = 3, byrow = TRUE))
summer.region.diet2017

ggarrange(summer2017.map, 
          ggarrange(summer2017.dendrogram, summer.region.diet2017, ncol = 1, labels = c("B", "C"),
                    vjust = 1.95), 
          nrow = 1, widths = c(0.58, 0.42), labels = c("A"), hjust = -1, vjust = 1.95)
ggsave("figures/summer/summer-regions2017.PNG", width = 12, height = 8, units = "in",
       bg = "white", dpi = 600)

## 12B. July 15 - Sep 30 cluster analysis ----
late.summer <- subset(summer.salish.ch, Month %in% c(7, 8, 9))
nrow(late.summer) # 670

late.summer.grid.tally <- late.summer %>% group_by(GridID) %>% tally()
late.summer.grid.tally5 <- subset(late.summer.grid.tally, n >= 5)
late.summer <- select(late.summer, -Cluster)

late.summer <- subset(late.summer, GridID %in% late.summer.grid.tally5$GridID)
nrow(late.summer) # 664
colSums(late.summer[20:35]) # stickleback not sampled during late summer
late.summer <- select(late.summer, -`Threespine stickleback`)
ncol(late.summer[20:34])

late.summer.grid.prop <- data.frame(matrix(NA, nrow = 15, ncol = 0))

for (i in 1:length(late.summer.grid.tally5$GridID)) {
  grid.subset <- subset(late.summer, GridID == late.summer.grid.tally5$GridID[i])
  mean.diet <- colMeans(grid.subset[20:34])
  late.summer.grid.prop <- cbind(late.summer.grid.prop, mean.diet)
  colnames(late.summer.grid.prop)[i] <- late.summer.grid.tally5$GridID[i]
}

colSums(late.summer.grid.prop) # expected

# calculate Bray-Curtis dissimilarity matrix
late.summer.grid.prop.t <- as.data.frame(t(late.summer.grid.prop))
late.summer.grid.bray <- vegdist(late.summer.grid.prop.t, method = "bray")

# hierarchical cluster analysis
late.summer.grid.cluster <- hclust(late.summer.grid.bray, "average")
plot(late.summer.grid.cluster, main = "", sub = "", xlab = "",
     ylab = "Bray-Curtis dissimilarity", hang = -1)

late.summer.coph <- cophenetic(late.summer.grid.cluster)
cor(late.summer.grid.bray, late.summer.coph) # 0.98

# determine optimal number of clusters
# calculate IndVal at stomach level
late.summer.indval.fill <- data.frame(k = 2:10, sum.sig.indicators = NA, 
                                        prop.sig.indicators = NA)

for (z in 2:10) {
  z.clusters <- data.frame(Cluster = cutree(late.summer.grid.cluster, k = z))
  z.clusters$GridID <- rownames(z.clusters)
  
  indiv.diet <- merge(x = late.summer, y = z.clusters, by = "GridID",
                      all.x = FALSE, all.y = FALSE)
  nrow(indiv.diet)
  
  set.seed(530)
  z.indval <- indval(indiv.diet[20:34], indiv.diet$Cluster, numitr = 1000)
  
  # calculate sum of significant indicators
  significant.indicators <- z.indval$indcls[z.indval$pval <= 0.05]
  sum.sig.indicators <- sum(significant.indicators)
  
  # calculate prop with sig indicators
  clusters.with.sig.indval <- factor(z.indval$maxcls[z.indval$pval <= 0.05])
  prop.sig.indicators <- length(levels(clusters.with.sig.indval)) / z
  
  late.summer.indval.fill$sum.sig.indicators[late.summer.indval.fill$k == z] <- sum.sig.indicators
  late.summer.indval.fill$prop.sig.indicators[late.summer.indval.fill$k == z] <- prop.sig.indicators
}

late.summer.indval.fill

late.summer.cluster3 <- data.frame(Cluster = cutree(late.summer.grid.cluster, 
                                                    k = 3))
late.summer.cluster3$GridID <- rownames(late.summer.cluster3)
late.summer.cluster3

# map
late.summer.n10.grid <- subset(grid, GridID %in% late.summer.grid.tally5$GridID)
late.summer.n10.grid <- merge(x = late.summer.n10.grid, y = late.summer.cluster3,
                              by = "GridID", all.x = FALSE, all.y = FALSE)
late.summer.n10.grid$Cluster <- factor(late.summer.n10.grid$Cluster)

late.summer.n10.grid$LateSummerRegion[late.summer.n10.grid$Cluster == 1] <- "Strait of Georgia"
late.summer.n10.grid$LateSummerRegion[late.summer.n10.grid$GridID %in% c(15, 16, 17, 23)] <- "Juan de Fuca Strait"
late.summer.n10.grid$LateSummerRegion[late.summer.n10.grid$Cluster == 2] <- "Haro Strait"
late.summer.n10.grid$LateSummerRegion[late.summer.n10.grid$Cluster == 3] <- "Haro Strait transition"

late.summer.n10.grid$LateSummerRegion <- factor(late.summer.n10.grid$LateSummerRegion,
  levels = c("Strait of Georgia", "Haro Strait transition", "Haro Strait", "Juan de Fuca Strait"))

late.summer.map <- ggplot() +
  geom_sf(data = bc.coastUTM) +
  geom_sf(data = late.summer.n10.grid, aes(fill = LateSummerRegion), alpha = 0.8,) +
  geom_sf_text(data = late.summer.n10.grid, aes(label = GridID), colour = "black", size = 4.5) +
  scale_fill_manual(values = c("#0072B2", "#758000", "#D55E00",  "#56b4e9")) +
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
        legend.position.inside = c(0.82, 0.91),
        legend.background = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(size = 11.5))
late.summer.map

# dendrogram
df.late.summer <- data.frame(x = 1:19, y = rep(-0.015, 19),
  region = rev(c("Haro Strait transition", "Haro Strait", "Juan de Fuca Strait",
             "Juan de Fuca Strait", rep("Strait of Georgia", 3),
             "Juan de Fuca Strait", "Juan de Fuca Strait", rep("Strait of Georgia", 10))))
df.late.summer

ggdendrogram(late.summer.grid.cluster)
late.summer.order <-late.summer.grid.cluster$labels[late.summer.grid.cluster$order]
late.summer.order
  
late.summer.dendrogram <- ggdendrogram(late.summer.grid.cluster) +  
  geom_point(data = df.late.summer, aes(x = x, y = y, colour = region),
             size = 5.5, shape = 15, show.legend = FALSE) +
  scale_colour_manual(values = c("#D55E00", "#758000", "#56b4e9", "#0072B2")) +
  annotate(geom = "rect", xmin = 0.55, xmax = 17.45, ymin = -0.045, ymax = -0.039,
           fill = "#0072B2") +
  annotate(geom = "rect", xmin = 17.55, xmax = 18.45, ymin = -0.045, ymax = -0.039,
           fill = "#D55E00") +
  annotate(geom = "rect", xmin = 18.55, xmax = 19.45, ymin = -0.045, ymax = -0.039,
           fill = "#758000") +
  annotate("text", x = 8.5, y = -0.058, label = "1", size = 4) +
  annotate("text", x = 18, y = -0.058, label = "2", size = 4) +
  annotate("text", x = 19, y = -0.058, label = "3", size = 4) +
  
  annotate(geom = "segment", x = 0.001, y = -0.02, yend = 0.65, colour = "black",
           linewidth = 0.9) +
  scale_x_continuous(limits = c(0, 20), expand = c(0, 0), breaks = c(1:19),
                  labels = late.summer.order) +
  scale_y_continuous(breaks = c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6)) +
  theme(axis.title.y = element_text(size = 14),
        axis.ticks.y.left = element_line(),
        axis.text.y = element_text(size = 11.5, angle = 0, colour = "black"),
        axis.text.x = element_text(size = 11, angle = 0, hjust = 0.5, 
                                   vjust = 15, colour = "black", margin = margin(b = -29)),
        plot.background = element_blank(),
        plot.margin = margin(t = -5, b = 2)) +
  ylab("Bray-Curtis dissimilarity")
late.summer.dendrogram

# calculate mean diet in the regions
late.summer <- merge(x = late.summer, y = late.summer.cluster3, by = "GridID",
                       all.x = FALSE, all.y = FALSE)
late.summer$LateSummerRegion[late.summer$Cluster == 1] <- "Strait of Georgia"
late.summer$LateSummerRegion[late.summer$GridID %in% c(15, 16, 17, 23)] <- "Juan de Fuca Strait"
late.summer$LateSummerRegion[late.summer$Cluster == 2] <- "Haro Strait"
late.summer$LateSummerRegion[late.summer$Cluster == 3] <- "Haro Strait transition"

late.summer.region.levels <- sort(unique(late.summer$LateSummerRegion))
late.summer.region.prop <- data.frame(matrix(NA, nrow = 15, ncol = 0))

for (t in 1:length(late.summer.region.levels)) {
  region.subset <- subset(late.summer, LateSummerRegion == late.summer.region.levels[t])
  mean.diet <- colMeans(region.subset[20:34])
  late.summer.region.prop <- cbind(late.summer.region.prop, mean.diet)
  colnames(late.summer.region.prop)[t] <- late.summer.region.levels[t]
}

# account for the rare prey removed (Juan de Fuca and Strait of Georgia only)
colSums(late.summer.region.prop)

rare.prey.adjust.late <- data.frame(t(data.frame(1 - colSums(late.summer.region.prop))))
rownames(rare.prey.adjust.late) <- "Rare prey"
names(rare.prey.adjust.late) <- names(late.summer.region.prop)
rare.prey.adjust.late

late.summer.region.prop <- rbind(late.summer.region.prop, rare.prey.adjust.late)
late.summer.region.prop <- late.summer.region.prop * 100
colSums(late.summer.region.prop) # 100
late.summer.region.prop$Group <- rownames(late.summer.region.prop)

late.summer.region.prop <- melt(late.summer.region.prop)
names(late.summer.region.prop) <- c("Group", "LateSummerRegion", "Percentage")

# is there really only herring and sand lance in Haro Strait in late summer?
late.summer.major <- subset(late.summer.region.prop, Percentage >= 5)
unique(late.summer.major$Group)

late.summer.region.prop$PlotGroup <- late.summer.region.prop$Group
late.summer.region.prop$PlotGroup[!late.summer.region.prop$PlotGroup %in% late.summer.major$Group] <- "Other"
late.summer.region.prop$PlotGroup[late.summer.region.prop$PlotGroup == "Pacific herring"] <- "Herring"
late.summer.region.prop$PlotGroup[late.summer.region.prop$PlotGroup == "Pacific sand lance"] <- "Sand lance"
late.summer.region.prop$PlotGroup[late.summer.region.prop$PlotGroup == "Northern anchovy"] <- "Anchovy"
late.summer.region.prop$PlotGroup[late.summer.region.prop$PlotGroup == "Gadiformes"] <- "Cod and hake"

sort(unique(late.summer.region.prop$PlotGroup))
late.summer.region.prop$PlotGroup <- factor(late.summer.region.prop$PlotGroup,
  levels = c("Herring", "Sand lance", "Anchovy", "Cod and hake", "Pandalid", "Other"))

late.summer.region.prop$LateSummerRegion <- as.character(late.summer.region.prop$LateSummerRegion)
late.summer.region.prop$LateSummerRegion[late.summer.region.prop$LateSummerRegion == "Strait of Georgia"] <- "Strait of\nGeorgia"
late.summer.region.prop$LateSummerRegion[late.summer.region.prop$LateSummerRegion == "Juan de Fuca Strait"] <- "Juan de Fuca\nStrait"
late.summer.region.prop$LateSummerRegion[late.summer.region.prop$LateSummerRegion == "Haro Strait transition"] <- "Haro Strait\ntransition"
late.summer.region.prop$LateSummerRegion <- factor(late.summer.region.prop$LateSummerRegion,
  levels = c("Strait of\nGeorgia", "Haro Strait", "Haro Strait\ntransition", "Juan de Fuca\nStrait"))

late.summer.region.n <- data.frame(table(late.summer$LateSummerRegion))
names(late.summer.region.n) <- c("LateSummerRegion", "n")

late.summer.region.n$LateSummerRegion <- as.character(late.summer.region.n$LateSummerRegion)
late.summer.region.n$LateSummerRegion[late.summer.region.n$LateSummerRegion == "Strait of Georgia"] <- "Strait of\nGeorgia"
late.summer.region.n$LateSummerRegion[late.summer.region.n$LateSummerRegion == "Juan de Fuca Strait"] <- "Juan de Fuca\nStrait"
late.summer.region.n$LateSummerRegion[late.summer.region.n$LateSummerRegion == "Haro Strait transition"] <- "Haro Strait\ntransition"
late.summer.region.n$LateSummerRegion <- factor(late.summer.region.n$LateSummerRegion,
                                                   levels = c("Strait of\nGeorgia", "Haro Strait", "Haro Strait\ntransition", "Juan de Fuca\nStrait"))

late.summer.diet <- ggplot() +
  geom_bar(data = late.summer.region.prop, 
           aes(x = LateSummerRegion, weight = Percentage, fill = PlotGroup),
           colour = "black", linewidth = 0.2) +
  scale_fill_manual(values = c("#0072B2", "#D55E00", "#009E73", "#897d1c",
                               "#4CA914", "grey60")) +
  geom_text(data = late.summer.region.n,
            aes(x = LateSummerRegion, y = 95, label = n)) +
  coord_cartesian(expand = FALSE) +
  ylab("Percentage (%)") + labs(fill = "Prey Category") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.title = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 14),
        axis.text.x = element_text(size = 11.5, colour = "black",
                                   margin = margin(b = -10)),
        axis.text.y = element_text(size = 11.5, colour = "black"),
        legend.text = element_text(size = 11),
        legend.key.size = unit(0.52, "cm"),
        legend.position = "bottom",
        panel.border = element_rect(linewidth = 1))
late.summer.diet

ggarrange(late.summer.map, 
          ggarrange(late.summer.dendrogram, late.summer.diet, ncol = 1, labels = c("B", "C"),
                    vjust = 1.95), 
          nrow = 1, widths = c(0.58, 0.42), labels = c("A"), hjust = -1, vjust = 1.95)
ggsave("figures/summer/late-summer-regions.PNG", width = 12, height = 8, units = "in",
       bg = "white", dpi = 600)
