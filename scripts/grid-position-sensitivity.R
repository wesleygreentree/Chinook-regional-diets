## Sensitivity of results to grid position

## 0. Script description ----

# We used an iterative procedure to identify the most appropriate grid for the
# sampling coverage and used this grid for analyses. Two sensitivity analyses 
# were conducted to determine if the highest ranked grid was representative of
# the 625 candidate grids. First, we determined if there was consensus among the
# 625 candidate grids on the optimal number of clusters in summer and winter. We
# repeated the summer and winter cluster analyses with each grid and determined 
# the grid's optimal number of summer and winter clusters with the two goodness-
# of-clustering indices. Second, we tested that the different grids resulted in
# similar regional boundaries at the common number of clusters determined to be
# optimal across all grids in each season. To visualize clusters, we used
# indicator species with IndVal >= 0.15. Clusters were visualized by
# capture location and by grid centroid. The results of both sensitivity analyses
# were compared to the summer and winter cluster analyses conducted with the 
# highest-ranked grid.

# 1. Converts spatial data to UTM coordinates before joining with grid.

# 2. A loop determines the optimal number of clusters, based on the number of
# clusters most often determined to be best using indicator species analysis. 
# This is done separately in summer and winter.

# 3. Determines the general regionalization pattern for summer and winter to
# assess if results are robust to the position of the grid. 


# load packages
library(tidyverse)
library(sf)
library(vegan)
library(cluster)
library(labdsv)
library(reshape2)
library(ggpubr)

# load data
bc.coast <- read_sf("data/BC-coast-shapefile/bc-coast.shp")
summer.salish.ch <- readRDS("data/out/summer-Salish-Chinook.RDS") # from scripts/data-processing.R
winter.salish.ch <- readRDS("data/out/winter-Salish-Chinook.RDS") # from scripts/data-processing.R

## 1. Spatial processing ----

# convert to UTM
bc.coastUTM <- st_transform(bc.coast, crs = st_crs(32610))

summer.salish.ch.sf <- st_as_sf(summer.salish.ch, coords = c("Longitude", "Latitude"),
                         crs = st_crs(4326), remove = FALSE)
summer.salish.chUTM <- st_transform(summer.salish.ch.sf, crs = st_crs(32610))
summer.salish.chUTM.coords <- do.call(rbind, st_geometry(summer.salish.chUTM)) %>% 
  as_tibble() %>% setNames(c("utm_x", "utm_y"))
summer.salish.chUTM <- bind_cols(summer.salish.chUTM, summer.salish.chUTM.coords)
nrow(summer.salish.chUTM)

winter.salish.ch.sf <- st_as_sf(winter.salish.ch, coords = c("Longitude", "Latitude"),
                                crs = st_crs(4326), remove = FALSE)
winter.salish.chUTM <- st_transform(winter.salish.ch.sf, crs = st_crs(32610))
winter.salish.chUTM.coords <- do.call(rbind, st_geometry(winter.salish.chUTM)) %>% 
  as_tibble() %>% setNames(c("utm_x", "utm_y"))
winter.salish.chUTM <- bind_cols(winter.salish.chUTM, winter.salish.chUTM.coords)
nrow(winter.salish.chUTM)

ggplot() +
  geom_sf(data = bc.coastUTM) +
  geom_point(data = summer.salish.chUTM, aes(x = utm_x, y = utm_y), colour = "red", alpha = 0.5) +
  geom_point(data = winter.salish.chUTM, aes(x = utm_x, y = utm_y), colour = "blue", alpha = 0.5) +
  coord_sf(datum = 4326,
           crs = st_crs(32610),
           xlim = c(300*10^3, 540*10^3), ylim = c(5325*10^3, 5600*10^3)) +
  theme_bw()

# set up for grid
size <- 25*10^3 # 25 km x 25 km grid cell 
x <- 1*10^3 # x distance the grid is moved
y <- 1*10^3 # y distance the grid is moved

# bottom corner: 305, 5315 km
x.start <- 305*10^3
y.start <- 5315*10^3

x.seq <- seq(from = x.start, to = (x.start + size - x), by = x)
x.seq / 1000 # makes sense
y.seq <- seq(from = y.start, to = (y.start + size - y), by = y)
y.seq / 1000

## 2. Optimal number of clusters ----

# dataframe to fill:
grid.sensitivity <- data.frame(matrix(NA, nrow = 0, ncol = 13))

names(grid.sensitivity) <- c("grid_id", "summer.n10.grids", "summer.n10.total",
  "winter.n10.grids", "winter.n10.total", "summer.silhouette.best.k",
  "summer.indval.best.k", "summer.indval.prop", "winter.silhouette.best.k", 
  "winter.indval.best.k", "winter.indval.prop", 
  "summer_none", "winter_none")

for (i in 1:length(x.seq)) {
  for (j in 1:length(y.seq)) {
    
    grid_id <- paste("grid", i, j, sep = "_")
    print(grid_id)
    which <- ((i - 1) * 25) + j
    
    bounding.box <- st_bbox(c(xmin = x.seq[i],
                              ymin = y.seq[j],
                              xmax = x.seq[i] + (size * 10),
                              ymax = y.seq[j] + (size * 12)),
                            crs = st_crs(32610))
    x.range <- paste(bounding.box["xmin"], bounding.box["xmax"], sep = ", ")
    y.range <- paste(bounding.box["ymin"], bounding.box["ymax"], sep = ", ")
    
    grid <- st_make_grid(bounding.box, cellsize = size,
                         crs = st_crs(32610), what = "polygons")
    grid2 <- grid %>% 
      cbind(data.frame(GridID = seq(1, to = length(grid), by = 1))) %>% 
      st_as_sf()
    
    # summer
    summer.joined <- st_join(x = summer.salish.chUTM, y = grid2)
    summer.joined <- st_drop_geometry(summer.joined) # drop geometry column
    summer.n <- nrow(subset(summer.joined, !is.na(GridID))) # to make sure grid includes all stomachs
    summer.tally <- summer.joined %>% group_by(GridID) %>% tally()
    summer.n10.sub <- subset(summer.tally, n >= 10) # number of grid cells with >= 10 stomachs
    summer.n10 <- nrow(summer.n10.sub)
    
    summer.n10.stomachs <- subset(summer.joined, GridID %in% summer.n10.sub$GridID)
    summer.n10.total <- nrow(summer.n10.stomachs) # number of stomachs in a grid with >= 10 stomachs
    summer.n10.total
    
    # winter
    winter.joined <- st_join(x = winter.salish.chUTM, y = grid2)
    winter.joined <- st_drop_geometry(winter.joined) # drop geometry column
    winter.n <- nrow(subset(winter.joined, !is.na(GridID))) # to make sure grid includes all stomachs
    winter.tally <- winter.joined %>% group_by(GridID) %>% tally()
    winter.n10.sub <- subset(winter.tally, n >= 10) # number of grid cells with >= 10 stomachs
    winter.n10 <- nrow(winter.n10.sub)
    
    winter.n10.stomachs <- subset(winter.joined, GridID %in% winter.n10.sub$GridID)
    winter.n10.total <- nrow(winter.n10.stomachs) # number of stomachs in a grid with >= 10 stomachs
    winter.n10.total
    
    # what happens if a prey type isn't present ever? Check as this affects calculations
    summer.col.sum <- data.frame(sum = colSums(summer.n10.stomachs[16:31]))
    summer.none <- nrow(subset(summer.col.sum, sum == 0))
    
    winter.col.sum <- data.frame(sum = colSums(winter.n10.stomachs[16:33]))
    winter.none <- nrow(subset(winter.col.sum, sum == 0))
    
    ## Summer cluster analysis
    summer.grid.levels <- sort(unique(summer.n10.sub$GridID))
    names(summer.n10.stomachs)
    summer.grid.prop <- data.frame(matrix(NA, nrow = 16, ncol = 0))
    
    for (k in 1:length(summer.grid.levels)) {
      summer.grid.subset <- subset(summer.n10.stomachs, GridID == summer.grid.levels[k])
      summer.mean.diet <- colMeans(summer.grid.subset[16:31])
      summer.grid.prop <- cbind(summer.grid.prop, summer.mean.diet)
      colnames(summer.grid.prop)[k] <- summer.grid.levels[k]
    }
    
    # hierarchical cluster analysis on Bray-Curtis dissimilarity matrix
    # with average linkage
    summer.grid.prop.t <- as.data.frame(t(summer.grid.prop))
    summer.grid.bray <- vegdist(summer.grid.prop.t, method = "bray")
    summer.grid.cluster <- hclust(summer.grid.bray, "average")
    #plot(summer.grid.cluster, hang = -1) # compare to summer plot
    
    # silhouette width (using "corrected" method that is not biased against
    # single-object clusters)
    summer.sil.fill <- data.frame(k = 2:10, corrected_silhouette = NA)
    for (b in 2:10) {
      sil <- silhouette(cutree(summer.grid.cluster, k = b), summer.grid.bray)
      sil <- data.frame(sil[, 1:3])
      
      freq <- data.frame(table(sil$cluster))
      n1 <- subset(freq, Freq == 1)
      n1.clusters <- n1$Var1
      
      corrected.sil <- subset(sil, !cluster %in% n1.clusters)
      summer.sil.fill$corrected_silhouette[summer.sil.fill$k == b] <- mean(corrected.sil$sil_width)
    }
    summer.silhouette.best.k <- summer.sil.fill[which.max(summer.sil.fill$corrected_silhouette), ]
    summer.silhouette.best.k <- summer.silhouette.best.k$k
    
    # indicator species analysis
    summer.indiv.indval.fill <- data.frame(k = 2:10, sum.sig.indicators = NA, 
                                    prop.sig.indicators = NA)
    
    for (z in 2:10) {
      z.clusters <- data.frame(Cluster = cutree(summer.grid.cluster, k = z))
      z.clusters$GridID <- rownames(z.clusters)
      
      indiv.diet <- merge(x = summer.n10.stomachs, y = z.clusters, by = "GridID",
                          all.x = FALSE, all.y = FALSE)
      
      names(indiv.diet)
      set.seed(530)
      z.indval <- indval(indiv.diet[17:32], indiv.diet$Cluster, numitr = 1000)
      
      # calculate sum of significant indicators
      significant.indicators <- z.indval$indcls[z.indval$pval < 0.05]
      sum.sig.indicators <- sum(significant.indicators)
      
      # calculate prop with sig indicators
      clusters.with.sig.indval <- factor(z.indval$maxcls[z.indval$pval < 0.05])
      prop.sig.indicators <- length(levels(clusters.with.sig.indval)) / z
      
      summer.indiv.indval.fill$sum.sig.indicators[summer.indiv.indval.fill$k == z] <- sum.sig.indicators
      summer.indiv.indval.fill$prop.sig.indicators[summer.indiv.indval.fill$k == z] <- prop.sig.indicators
    }
    
    summer.indval.prop1 <- subset(summer.indiv.indval.fill, prop.sig.indicators == 1)

    if (nrow(summer.indval.prop1) >= 1) {
      summer.indval.best.k <- summer.indval.prop1[which.max(summer.indval.prop1$sum.sig.indicators), ]
      summer.indval.best.k <- summer.indval.best.k$k
      summer.indval.prop <- "one"
    }
    
    if (nrow(summer.indval.prop1) == 0) {
      summer.indval.best.k <- summer.indiv.indval.fill[which.max(summer.indiv.indval.fill$sum.sig.indicators), ]
      summer.indval.best.k <- summer.indval.best.k$k
      sumer.indval.prop <- "not one"
    }
    
    ## Winter cluster analysis
    winter.grid.levels <- sort(unique(winter.n10.sub$GridID))
    names(winter.n10.stomachs)
    winter.grid.prop <- data.frame(matrix(NA, nrow = 18, ncol = 0))
    
    for (k in 1:length(winter.grid.levels)) {
      winter.grid.subset <- subset(winter.n10.stomachs, GridID == winter.grid.levels[k])
      winter.mean.diet <- colMeans(winter.grid.subset[16:33])
      winter.grid.prop <- cbind(winter.grid.prop, winter.mean.diet)
      colnames(winter.grid.prop)[k] <- winter.grid.levels[k]
    }
    
    # hierarchical cluster analysis on Bray-Curtis dissimilarity matrix
    # with average linkage
    winter.grid.prop.t <- as.data.frame(t(winter.grid.prop))
    winter.grid.bray <- vegdist(winter.grid.prop.t, method = "bray")
    winter.grid.cluster <- hclust(winter.grid.bray, "average")
    #plot(winter.grid.cluster, hang = -1) # compare to winter plot
    
    # silhouette width
    winter.sil.fill <- data.frame(k = 2:10, corrected_silhouette = NA)
    for (b in 2:10) {
      sil <- silhouette(cutree(winter.grid.cluster, k = b), winter.grid.bray)
      sil <- data.frame(sil[, 1:3])
      
      freq <- data.frame(table(sil$cluster))
      n1 <- subset(freq, Freq == 1)
      n1.clusters <- n1$Var1
      
      corrected.sil <- subset(sil, !cluster %in% n1.clusters)
      winter.sil.fill$corrected_silhouette[winter.sil.fill$k == b] <- mean(corrected.sil$sil_width)
    }
    winter.silhouette.best.k <- winter.sil.fill[which.max(winter.sil.fill$corrected_silhouette), ]
    winter.silhouette.best.k <- winter.silhouette.best.k$k
    
    # indicator species analysis
    winter.indiv.indval.fill <- data.frame(k = 2:10, sum.sig.indicators = NA, 
                                           prop.sig.indicators = NA)
    
    for (z in 2:10) {
      z.clusters <- data.frame(Cluster = cutree(winter.grid.cluster, k = z))
      z.clusters$GridID <- rownames(z.clusters)
      
      indiv.diet <- merge(x = winter.n10.stomachs, y = z.clusters, by = "GridID",
                          all.x = FALSE, all.y = FALSE)
      names(indiv.diet)
      set.seed(530)
      z.indval <- indval(indiv.diet[17:34], indiv.diet$Cluster, numitr = 1000)
      
      # calculate sum of significant indicators
      significant.indicators <- z.indval$indcls[z.indval$pval < 0.05]
      sum.sig.indicators <- sum(significant.indicators)
      
      # calculate prop with sig indicators
      clusters.with.sig.indval <- factor(z.indval$maxcls[z.indval$pval < 0.05])
      prop.sig.indicators <- length(levels(clusters.with.sig.indval)) / z
      
      winter.indiv.indval.fill$sum.sig.indicators[winter.indiv.indval.fill$k == z] <- sum.sig.indicators
      winter.indiv.indval.fill$prop.sig.indicators[winter.indiv.indval.fill$k == z] <- prop.sig.indicators
    }
    
    winter.indval.prop1 <- subset(winter.indiv.indval.fill, prop.sig.indicators == 1)
    
    if (nrow(winter.indval.prop1) >= 1) {
      winter.indval.best.k <- winter.indval.prop1[which.max(winter.indval.prop1$sum.sig.indicators), ]
      winter.indval.best.k <- winter.indval.best.k$k
      winter.indval.prop <- "one"
    }
    
    if (nrow(winter.indval.prop1) == 0) {
      winter.indval.best.k <- winter.indiv.indval.fill[which.max(winter.indiv.indval.fill$sum.sig.indicators), ]
      winter.indval.best.k <- winter.indval.best.k$k
      winter.indval.prop <- "not one"
    }
    
    row <- data.frame(t(c(grid_id, summer.n10, summer.n10.total, winter.n10, winter.n10.total,
             summer.silhouette.best.k, summer.indval.best.k, summer.indval.prop,
             winter.silhouette.best.k, winter.indval.best.k, winter.indval.prop,
             summer.none, winter.none)))
    names(row) <- c("grid_id", "summer.n10.grids", "summer.n10.total",
       "winter.n10.grids", "winter.n10.total", "summer.silhouette.best.k",
       "summer.indval.best.k", "summer.indval.prop", "winter.silhouette.best.k", 
       "winter.indval.best.k", "winter.indval.prop",
       "summer_none", "winter_none")
    grid.sensitivity <- rbind(grid.sensitivity, row)
  }
}
saveRDS(grid.sensitivity, "data/out/grid-sensitivity-optimal-k.RDS")

unique(grid.sensitivity$summer_none) # always 0
unique(grid.sensitivity$winter_none) # always 0

## makes plots
#grid.sensitivity <- readRDS("data/out/grid-sensitivity-optimal-k.RDS")

table(grid.sensitivity$summer.silhouette.best.k) # 2 then 3
table(grid.sensitivity$winter.silhouette.best.k) # 2 then 5

table(grid.sensitivity$summer.indval.best.k) # 3
table(grid.sensitivity$winter.indval.best.k) # 5 then 6 and 4

summer.sil.table <- data.frame(table(grid.sensitivity$summer.silhouette.best.k))
names(summer.sil.table) <- c("Number of clusters", "Count")
summer.sil.table$Season <- "Summer"
summer.sil.table$Method <- "Silhouette width"

winter.sil.table <- data.frame(table(grid.sensitivity$winter.silhouette.best.k))
names(winter.sil.table) <- c("Number of clusters", "Count")
winter.sil.table$Season <- "Winter"
winter.sil.table$Method <- "Silhouette width"

summer.indval.table <- data.frame(table(grid.sensitivity$summer.indval.best.k))
names(summer.indval.table) <- c("Number of clusters", "Count")
summer.indval.table$Season <- "Summer"
summer.indval.table$Method <- "Sum of indicator species"

winter.indval.table <- data.frame(table(grid.sensitivity$winter.indval.best.k))
names(winter.indval.table) <- c("Number of clusters", "Count")
winter.indval.table$Season <- "Winter"
winter.indval.table$Method <- "Sum of indicator species"

best.k.table <- rbind(summer.sil.table, winter.sil.table, summer.indval.table,
                      winter.indval.table)
best.k.table$num_clusters <- as.numeric(as.character(best.k.table$`Number of clusters`))
range(best.k.table$num_clusters) 

ggplot() +
  geom_col(data = best.k.table,
           aes(x = num_clusters, y = Count, fill = Season),
           show.legend = FALSE) +
  scale_fill_manual(values = c("firebrick3", "dodgerblue3")) +
  xlab("Number of clusters") +
  scale_x_continuous(breaks = c(2:10)) +
  facet_grid(Season ~ Method) +
  theme_bw() +
  theme(strip.background = element_blank(), strip.text = element_text(size = 11))
ggsave("figures/supplement/best-k-sensitivity.PNG", 
       width = 17, height = 12, units = "cm", dpi = 600)

## 3. What is the general regionalization pattern? ----
summer.overall.best.k <- 3
winter.overall.best.k <- 5

summer.stomachs.fill <- data.frame(matrix(NA, ncol = 7, nrow = 0))
names(summer.stomachs.fill) <- c("Cluster", "grid", "Fish Code", "utm_x", "utm_y", 
                                  "GridID", "indicators")
winter.stomachs.fill <- data.frame(matrix(NA, ncol = 7, nrow = 0))
names(winter.stomachs.fill) <- c("Cluster", "grid", "Fish Code", "utm_x", "utm_y", 
                                 "GridID", "indicators")

summer.centroid.indicators <- data.frame(matrix(NA, ncol = 5, nrow = 0))
names(summer.centroid.indicators) <- c("Cluster", "GridID", "utm_x", "utm_y", 
                                       "indicators")
winter.centroid.indicators <- data.frame(matrix(NA, ncol = 5, nrow = 0))
names(winter.centroid.indicators) <- c("Cluster", "GridID", "utm_x", "utm_y", 
                                       "indicators")

for (i in 1:length(x.seq)) {
  for (j in 1:length(y.seq)) {
    
    grid_id <- paste("grid", i, j, sep = "_")
    print(grid_id)
    which <- ((i - 1) * 25) + j
    
    bounding.box <- st_bbox(c(xmin = x.seq[i],
                              ymin = y.seq[j],
                              xmax = x.seq[i] + (size * 10),
                              ymax = y.seq[j] + (size * 12)),
                            crs = st_crs(32610))
    x.range <- paste(bounding.box["xmin"], bounding.box["xmax"], sep = ", ")
    y.range <- paste(bounding.box["ymin"], bounding.box["ymax"], sep = ", ")
    
    grid <- st_make_grid(bounding.box, cellsize = size,
                         crs = st_crs(32610), what = "polygons")
    grid2 <- grid %>% 
      cbind(data.frame(GridID = seq(1, to = length(grid), by = 1))) %>% 
      st_as_sf()
    
    # calculate centroid of each grid cell
    suppressWarnings(grid.centroid <- st_centroid(grid2)) # not concerned about warning: st_centroid assumes attributes are constant over geometries
    grid.centroid.coords <- do.call(rbind, st_geometry(grid.centroid)) %>% 
      as_tibble() %>% setNames(c("utm_x", "utm_y"))
    grid.centroid <- bind_cols(grid.centroid, grid.centroid.coords)
    grid.centroid <- st_drop_geometry(grid.centroid)
    
    # summer
    summer.joined <- st_join(x = summer.salish.chUTM, y = grid2)
    summer.joined <- st_drop_geometry(summer.joined) # drop geometry column
    summer.n <- nrow(subset(summer.joined, !is.na(GridID))) # to make sure grid includes all stomachs
    summer.tally <- summer.joined %>% group_by(GridID) %>% tally()
    summer.n10.sub <- subset(summer.tally, n >= 10) # number of grid cells with >= 10 stomachs
    summer.n10 <- nrow(summer.n10.sub)
    
    summer.n10.stomachs <- subset(summer.joined, GridID %in% summer.n10.sub$GridID)
    summer.n10.total <- nrow(summer.n10.stomachs) # number of stomachs in a grid with >= 10 stomachs
    
    # winter
    winter.joined <- st_join(x = winter.salish.chUTM, y = grid2)
    winter.joined <- st_drop_geometry(winter.joined) # drop geometry column
    winter.n <- nrow(subset(winter.joined, !is.na(GridID))) # to make sure grid includes all stomachs
    winter.tally <- winter.joined %>% group_by(GridID) %>% tally()
    winter.n10.sub <- subset(winter.tally, n >= 10) # number of grid cells with >= 10 stomachs
    winter.n10 <- nrow(winter.n10.sub)
    
    winter.n10.stomachs <- subset(winter.joined, GridID %in% winter.n10.sub$GridID)
    winter.n10.total <- nrow(winter.n10.stomachs) # number of stomachs in a grid with >= 10 stomachs
    
    ## Summer cluster analysis
    summer.grid.levels <- sort(unique(summer.n10.sub$GridID))
    names(summer.n10.stomachs)
    summer.grid.prop <- data.frame(matrix(NA, nrow = 16, ncol = 0))
    
    for (k in 1:length(summer.grid.levels)) {
      summer.grid.subset <- subset(summer.n10.stomachs, GridID == summer.grid.levels[k])
      summer.mean.diet <- colMeans(summer.grid.subset[16:31])
      summer.grid.prop <- cbind(summer.grid.prop, summer.mean.diet)
      colnames(summer.grid.prop)[k] <- summer.grid.levels[k]
    }
    
    # hierarchical cluster analysis on Bray-Curtis dissimilarity matrix
    # with average linkage
    summer.grid.prop.t <- as.data.frame(t(summer.grid.prop))
    summer.grid.bray <- vegdist(summer.grid.prop.t, method = "bray")
    summer.grid.cluster <- hclust(summer.grid.bray, "average")
    #plot(summer.grid.cluster, hang = -1) # compare to summer plot
    
    # Indicator species analysis with summer best k
    summer.clusters <- data.frame(Cluster = cutree(summer.grid.cluster, k = summer.overall.best.k))
    summer.clusters$GridID <- rownames(summer.clusters)
    
    summer.n10.stomachs <- merge(x = summer.n10.stomachs, y = summer.clusters, by = "GridID",
                        all.x = FALSE, all.y = FALSE)
    
    names(summer.n10.stomachs)
    set.seed(530)
    summer.indval <- indval(summer.n10.stomachs[17:32], summer.n10.stomachs$Cluster, numitr = 1000)
    
    # visualize with indicators with IndVal >= 0.15
    summer.sig.indval <- data.frame(Cluster = factor(summer.indval$maxcls[summer.indval$indcls >= 0.15]))
    summer.sig.indval$Group <- rownames(summer.sig.indval)
  
    # Assign indicator species to individual stomachs
    summer.sig.indval.wide <- dcast(summer.sig.indval, Cluster ~ Group)
    summer.sig.indval.wide <- summer.sig.indval.wide %>% 
      unite(indicators, all_of(2:ncol(summer.sig.indval.wide)), 
            sep = ", ", na.rm = TRUE, remove = FALSE) # so it works if they are all NA?
    summer.sig.indval.wide <- select(summer.sig.indval.wide, Cluster, indicators)
    
    summer.n10.stomachs$grid <- grid_id
    
    summer.stomachs <- merge(x = select(summer.n10.stomachs, grid, `Fish Code`, utm_x, utm_y, GridID, Cluster),
                             y = summer.sig.indval.wide, by = "Cluster", 
                             all.x = TRUE, all.y = FALSE)
    # save to larger dataframe
    summer.stomachs.fill <-  rbind(summer.stomachs.fill, summer.stomachs)
    
    # Assign indicator species to grid centroids
    summer.grid.cluster.centroid <- merge(x = grid.centroid, y = summer.clusters,
                              by = "GridID", all.x = FALSE, all.y = FALSE)
    
    summer.grid.centroid.indicators <- merge(x = summer.grid.cluster.centroid,
      y = summer.sig.indval.wide, by = "Cluster", all.x = TRUE, all.y = FALSE)
    
    # save to larger dataframe
    summer.centroid.indicators <- rbind(summer.centroid.indicators, summer.grid.centroid.indicators)

    ## Winter cluster analysis
    winter.grid.levels <- sort(unique(winter.n10.sub$GridID))
    names(winter.n10.stomachs)
    winter.grid.prop <- data.frame(matrix(NA, nrow = 18, ncol = 0))
    
    for (k in 1:length(winter.grid.levels)) {
      winter.grid.subset <- subset(winter.n10.stomachs, GridID == winter.grid.levels[k])
      winter.mean.diet <- colMeans(winter.grid.subset[16:33])
      winter.grid.prop <- cbind(winter.grid.prop, winter.mean.diet)
      colnames(winter.grid.prop)[k] <- winter.grid.levels[k]
    }
    
    # hierarchical cluster analysis on Bray-Curtis dissimilarity matrix
    # with average linkage
    winter.grid.prop.t <- as.data.frame(t(winter.grid.prop))
    winter.grid.bray <- vegdist(winter.grid.prop.t, method = "bray")
    winter.grid.cluster <- hclust(winter.grid.bray, "average")
    #plot(winter.grid.cluster, hang = -1) # compare to winter plot
    
    # Indicator species analysis with winter best k
    winter.clusters <- data.frame(Cluster = cutree(winter.grid.cluster, k = winter.overall.best.k))
    winter.clusters$GridID <- rownames(winter.clusters)
    
    winter.n10.stomachs <- merge(x = winter.n10.stomachs, y = winter.clusters, by = "GridID",
                                 all.x = FALSE, all.y = FALSE)
    nrow(winter.n10.stomachs)
    
    names(winter.n10.stomachs)
    set.seed(530)
    winter.indval <- indval(winter.n10.stomachs[17:34], winter.n10.stomachs$Cluster, numitr = 1000)
    
    winter.sig.indval <- data.frame(Cluster = factor(winter.indval$maxcls[winter.indval$indcls >= 0.15])) # what min value
    winter.sig.indval$Group <- rownames(winter.sig.indval)
    
    # Assign indicator species to individual stomachs
    winter.sig.indval.wide <- dcast(winter.sig.indval, Cluster ~ Group)
    winter.sig.indval.wide <- winter.sig.indval.wide %>% 
      unite(indicators, all_of(2:ncol(winter.sig.indval.wide)), 
            sep = ", ", na.rm = TRUE, remove = FALSE) # so it works if they are all NA?
    winter.sig.indval.wide <- select(winter.sig.indval.wide, Cluster, indicators)
    
    winter.n10.stomachs$grid <- grid_id
    winter.stomachs <- merge(x = select(winter.n10.stomachs, grid, `Fish Code`, utm_x, utm_y, GridID, Cluster),
                             y = winter.sig.indval.wide, by = "Cluster", 
                             all.x = TRUE, all.y = FALSE)
    # save to larger dataframe
    winter.stomachs.fill <-  rbind(winter.stomachs.fill, winter.stomachs)
    
    # Assign indicator species to grid centroids
    winter.grid.cluster.centroid <- merge(x = grid.centroid, y = winter.clusters,
                                          by = "GridID", all.x = FALSE, all.y = FALSE)
    
    winter.grid.centroid.indicators <- merge(x = winter.grid.cluster.centroid,
                                             y = winter.sig.indval.wide, by = "Cluster", 
                                             all.x = TRUE, all.y = FALSE)
    
    # save to larger dataframe
    winter.centroid.indicators <- rbind(winter.centroid.indicators, winter.grid.centroid.indicators)
  }
}

nrow(summer.stomachs.fill) # 599,233
sum(as.numeric(grid.sensitivity$summer.n10.total)) # 599,233

nrow(summer.centroid.indicators) # 12,632
sum(as.numeric(grid.sensitivity$summer.n10.grids)) # 12,632

nrow(winter.stomachs.fill) # 444,526
sum(as.numeric(grid.sensitivity$winter.n10.total)) # 444,526

nrow(winter.centroid.indicators) # 9538
sum(as.numeric(grid.sensitivity$winter.n10.grids)) # 9538

## 3A Summer locations colour-coded by indicator species ----
# at the level of the location nested within grid # e.g., not affected by number of stomachs
nrow(summer.stomachs.fill) # 599,233
table(summer.stomachs.fill$indicators)

summer.location.indicators <- summer.stomachs.fill %>% 
  group_by(grid, utm_x, utm_y, indicators) %>% tally()
nrow(summer.location.indicators) # 75,703
table(summer.location.indicators$indicators)

summer.location.indicators.wide <- dcast(summer.location.indicators, utm_x + utm_y ~ indicators)
colSums(summer.location.indicators.wide[3:ncol(summer.location.indicators.wide)])
range(summer.location.indicators.wide$`Pacific herring`, na.rm = TRUE) # 0-625, good

summer.location.indicators.wide$location_total <- rowSums(summer.location.indicators.wide[3:ncol(summer.location.indicators.wide)])
range(summer.location.indicators.wide$location_total) # 32-65

summer.location.indicators.wide[is.na(summer.location.indicators.wide)] <- 0

names(summer.location.indicators.wide)
summer.location.indicators.wide$Myctophid <- summer.location.indicators.wide$`Myctophid, Pandalid, Squid` +
  summer.location.indicators.wide$`Myctophid, Squid`
summer.location.indicators.wide$Squid <- summer.location.indicators.wide$`Myctophid, Pandalid, Squid` +
  summer.location.indicators.wide$`Myctophid, Squid`
summer.location.indicators.wide$Pandalid <- summer.location.indicators.wide$`Myctophid, Pandalid, Squid`
summer.location.indicators.wide$`Northern anchovy2` <- summer.location.indicators.wide$`Northern anchovy` +
  summer.location.indicators.wide$`Northern anchovy, Surfperch`
summer.location.indicators.wide$Surfperch <- summer.location.indicators.wide$`Northern anchovy, Surfperch`

summer.location.indicators.wide <- select(summer.location.indicators.wide, -`Myctophid, Pandalid, Squid`,
  -`Myctophid, Squid`, -`Northern anchovy`, -`Northern anchovy, Surfperch`)
names(summer.location.indicators.wide)
names(summer.location.indicators.wide)[10] <- "Northern anchovy"

# calculate proportions
summer.location.indicators.prop <- melt(summer.location.indicators.wide, 
                                        id.vars = c("utm_x", "utm_y", "location_total"))
summer.location.indicators.prop$Proportion <- summer.location.indicators.prop$value / summer.location.indicators.prop$location_total

names(summer.location.indicators.prop)[4] <- "Group"
summer.location.indicators.prop <- subset(summer.location.indicators.prop, Group != "NA") # remove no indicator group

summer.location.indicators.prop$Proportion[summer.location.indicators.prop$Proportion == 0] <- NA

table(summer.location.indicators.prop$Group)
summer.location.indicators.prop$Group <- factor(summer.location.indicators.prop$Group,
  levels = c("Pacific herring", "Pacific sand lance", "Northern anchovy", "Surfperch",
             "Myctophid", "Squid", "Pandalid"))
view(subset(summer.location.indicators.prop, Group %in% c("Myctophid", "Squid", "Pandalid") &
              Proportion > 0))
# myctophid, squid, and pandalid look very similar because proportions are nearly identical
# (these are very rare indicator species in summer)

ggplot() +
  geom_sf(data = bc.coast) +
  geom_point(data = summer.location.indicators.prop, 
             aes(x = utm_x, utm_y, colour = Proportion)) + 
  scale_colour_viridis_c(option = "plasma", direction = 1, na.value = "grey40") +
  coord_sf(datum = 4326,
           crs = st_crs(32610),
           xlim = c(300*10^3, 540*10^3), ylim = c(5325*10^3, 5600*10^3)) +
  theme_bw() + 
  theme(strip.background = element_blank(), strip.text = element_text(size = 11),
        axis.title = element_blank(), axis.text = element_blank(),
        axis.ticks = element_blank(), panel.grid = element_blank(),
        legend.position = "inside",
        legend.position.inside = c(0.88, 0.22)) +
  facet_wrap(~ Group, nrow = 2)
ggsave("figures/summer/summer-sites-indicators.PNG", 
       width = 17, height = 12, units = "cm", dpi = 600)

## 3B Winter locations colour-coded by indicator species ----
# at the level of the location nested within grid # e.g., not affected by number of stomachs
nrow(winter.stomachs.fill) # 444,526
table(winter.stomachs.fill$indicators)
nrow(subset(winter.stomachs.fill, is.na(indicators))) # 5031

winter.location.indicators <- winter.stomachs.fill %>% 
  group_by(grid, utm_x, utm_y, indicators) %>% tally()
nrow(winter.location.indicators) # 57,919
table(winter.location.indicators$indicators)
nrow(subset(winter.location.indicators, is.na(indicators))) # 794

winter.location.indicators.wide <- dcast(winter.location.indicators, utm_x + utm_y ~ indicators)
colSums(winter.location.indicators.wide[3:ncol(winter.location.indicators.wide)]) # 794 NA
range(winter.location.indicators.wide$`Pacific herring`, na.rm = TRUE) # 0-625

winter.location.indicators.wide$location_total <- rowSums(winter.location.indicators.wide[3:ncol(winter.location.indicators.wide)])
range(winter.location.indicators.wide$location_total) # 20-625

names(winter.location.indicators.wide)
winter.location.indicators.wide$Krill <- winter.location.indicators.wide$Euphausiid +
  winter.location.indicators.wide$`Euphausiid, Mysid`
winter.location.indicators.wide$Mysid2 <- winter.location.indicators.wide$`Euphausiid, Mysid` +
  winter.location.indicators.wide$Mysid
winter.location.indicators.wide$`Northern anchovy2` <- winter.location.indicators.wide$`Northern anchovy` +
  winter.location.indicators.wide$`Northern anchovy, Surfperch`
winter.location.indicators.wide$Surfperch2 <- winter.location.indicators.wide$`Northern anchovy, Surfperch` +
  winter.location.indicators.wide$Surfperch

winter.location.indicators.wide <- select(winter.location.indicators.wide, -Euphausiid, -`Euphausiid, Mysid`,
                                          -Mysid, -`Northern anchovy`, -`Northern anchovy, Surfperch`,
                                          -Surfperch)
names(winter.location.indicators.wide)[9:11] <- c("Mysid", "Northern anchovy", "Surfperch")

# calculate proportions
winter.location.indicators.prop <- melt(winter.location.indicators.wide,
                                        id.vars = c("utm_x", "utm_y", "location_total"))
winter.location.indicators.prop$Proportion <- winter.location.indicators.prop$value / winter.location.indicators.prop$location_total

names(winter.location.indicators.prop)[4] <- "Group"
winter.location.indicators.prop <- subset(winter.location.indicators.prop, Group != "NA") # remove no indicator group

winter.location.indicators.prop$Proportion[winter.location.indicators.prop$Proportion == 0] <- NA

table(winter.location.indicators.prop$Group)
winter.location.indicators.prop$Group <- factor(winter.location.indicators.prop$Group,
  levels = c("Pacific herring", "Pacific sand lance", "Northern anchovy", "Surfperch",
             "Myctophid", "Mysid", "Krill"))

ggplot() +
  geom_sf(data = bc.coast) +
  geom_point(data = winter.location.indicators.prop, 
             aes(x = utm_x, utm_y, colour = Proportion)) + 
  scale_colour_viridis_c(option = "plasma", direction = 1, na.value = "grey40") +
  coord_sf(datum = 4326,
           crs = st_crs(32610),
           xlim = c(300*10^3, 540*10^3), ylim = c(5325*10^3, 5600*10^3)) +
  theme_bw() + 
  theme(strip.background = element_blank(), strip.text = element_text(size = 11),
        axis.title = element_blank(), axis.text = element_blank(),
        axis.ticks = element_blank(), panel.grid = element_blank(),
        legend.position = "inside",
        legend.position.inside = c(0.88, 0.22)) +
  facet_wrap(~ Group, nrow = 2)
ggsave("figures/winter/winter-sites-indicators.PNG", 
       width = 17, height = 12, units = "cm", dpi = 600)

## 3C Plot summer grid centroids by indicator species ----
head(summer.centroid.indicators)
table(summer.centroid.indicators$indicators)

summer.centroid.indicators$indicators <- factor(summer.centroid.indicators$indicators,
  levels = c("Pacific herring", "Pacific sand lance", "Northern anchovy",
             "Northern anchovy, Surfperch", "Myctophid, Pandalid, Squid",
             "Myctophid, Squid"))

ggplot() +
  geom_sf(data = bc.coast) +
  geom_point(data = summer.centroid.indicators,
             aes(x = utm_x, y = utm_y, colour = indicators),  
             size = 0.0001, alpha = 0.8) +
  scale_colour_manual(values = c("#0072B2", "#D55E00", "#009E73", "#009E73",
                                 "#5823AD", "#5823AD", "grey60"),
     guide = guide_legend(override.aes = list(size = 1.5))) +
  coord_sf(datum = 4326,
           crs = st_crs(32610),
           xlim = c(300*10^3, 540*10^3), ylim = c(5325*10^3, 5600*10^3)) +
  labs(colour = "Indicator species") +
  theme_bw() +
  theme(axis.title = element_blank(), axis.text = element_blank(),
        axis.ticks = element_blank(), panel.grid = element_blank(),
        legend.background = element_blank(), legend.key = element_blank(),
        legend.title = element_text(size = 11), legend.text = element_text(size = 10),
        legend.position = "right")
ggsave("figures/summer/summer-centroid-indicators.PNG",
       width = 17, height = 12, units = "cm", dpi = 600)

## 4D Plot winter grid centroids by indicator species ----
nrow(winter.centroid.indicators) # 9538
table(winter.centroid.indicators$indicators)

winter.centroid.indicators$indicators[winter.centroid.indicators$indicators == "Euphausiid"] <- "Krill"
winter.centroid.indicators$indicators[winter.centroid.indicators$indicators == "Euphausiid, Mysid"] <- "Krill, Mysid"

winter.centroid.indicators$indicators <- factor(winter.centroid.indicators$indicators,
  levels = c("Pacific herring", "Pacific sand lance", "Northern anchovy",
             "Northern anchovy, Surfperch", "Surfperch", "Myctophid",
             "Mysid", "Krill, Mysid", "Krill"))

ggplot() +
  geom_sf(data = bc.coast) +
  geom_point(data = winter.centroid.indicators,
             aes(x = utm_x, y = utm_y, colour = indicators),
             size = 0.0001, alpha = 0.8) +
  scale_colour_manual(values = c("#0072B2", "#D55E00", "#009E73", "#009E73", "#C05D94",
                                 "#5823AD", "#897d1c", "#897d1c", "#deb887", "grey60"),
                      guide = guide_legend(override.aes = list(size = 1.5))) +
  coord_sf(datum = 4326,
           crs = st_crs(32610),
           xlim = c(300*10^3, 540*10^3), ylim = c(5325*10^3, 5600*10^3)) +
  labs(colour = "Indicator species") +
  theme_bw() +
  theme(axis.title = element_blank(), axis.text = element_blank(),
        axis.ticks = element_blank(), panel.grid = element_blank(),
        legend.background = element_blank(), legend.key = element_blank(),
        legend.title = element_text(size = 11),
        legend.text = element_text(size = 10),
        legend.position = "right")
ggsave("figures/winter/winter-centroid-indicators.PNG", 
       width = 17, height = 12, units = "cm", dpi = 600)
