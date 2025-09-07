# Select the best grid for analysis

## 0. Script description ----
# We used a grid-based approach to analyze Chinook salmon diet composition, so
# that the cluster analysis of spatial units was not constrained by prior
# regionalizations of the Salish Sea. However, the placement of a grid over
# the sampling coverage may influence results. This script determines
# the grid position that best matches the sampling coverage.

# 1. Load data ----
library(tidyverse)
library(sf)
library(ggpubr)

bc.coast <- read_sf("data/BC-coast-shapefile/bc-coast.shp")

salish.ch <- readRDS("data/out/diet-data-for-grid-selection.RDS") # from scripts/data-processing.R
table(salish.ch$Season)

## 2. Spatial data wrangling ----

# transform to UTM zone 10N (EPSG 32610)
bc.coastUTM <- st_transform(bc.coast, crs = st_crs(32610))
st_crs(bc.coastUTM) # note that units are in metres

salish.ch.sf <- st_as_sf(salish.ch, coords = c("Longitude","Latitude"),
                         crs = st_crs(4326), remove = FALSE)
salish.chUTM <- st_transform(salish.ch.sf, crs = st_crs(32610))
salish.ch.summerUTM <- subset(salish.chUTM, Season == "Summer")
salish.ch.winterUTM <- subset(salish.chUTM, Season == "Winter")

size <- 25*10^3 # 25 km x 25 km grid cell 
x <- 1*10^3 # x distance the grid is moved
y <- 1*10^3 # y distance the grid is moved

# bottom corner of the grid: 305 km, 5315 km in UTM zone 10
x.start <- 305*10^3
y.start <- 5315*10^3

x.seq <- seq(from = x.start, to = (x.start + size - x), by = x)
y.seq <- seq(from = y.start, to = (y.start + size - y), by = y)

## 3. Iterative loop to assess each grid ----
grid.info <- data.frame(matrix(NA, nrow = 0, ncol = 10))
names(grid.info) <- c("grid_id", "which", "x.min", "y.min",
                      "summer.n", "summer.n10", "summer.n10.total",
                      "winter.n", "winter.n10", "winter.n10.total")

# make empty lists to store grid, summer_joined, winter_joined dataframes
grid.list <- vector(mode = "list", length = 625)
summer.joined.list <- vector(mode = "list", length = 625)
winter.joined.list <- vector(mode = "list", length = 625)

for (i in 1:length(x.seq)) {
  for (j in 1:length(y.seq)) {
    grid_id <- paste("grid", i, j, sep = "_")
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
    summer.joined <- st_join(x = salish.ch.summerUTM,
                             y = grid2)
    summer.n <- nrow(subset(summer.joined, !is.na(GridID))) # to make sure grid includes all stomachs
    summer.tally <- summer.joined %>% group_by(GridID) %>% tally()
    summer.n10.sub <- subset(summer.tally, n >= 10) # number of grid cells with >= 10 stomachs
    summer.n10 <- nrow(summer.n10.sub)
    
    summer.n10.stomachs <- subset(summer.joined, GridID %in% summer.n10.sub$GridID)
    summer.n10.total <- nrow(summer.n10.stomachs) # number of stomachs in a grid with >= 10 stomachs
    
    # winter
    winter.joined <- st_join(x = salish.ch.winterUTM,
                             y = grid2)
    winter.n <- nrow(subset(winter.joined, !is.na(GridID))) # to make sure grid includes all stomachs
    winter.tally <- winter.joined %>% group_by(GridID) %>% tally()
    winter.n10.sub <- subset(winter.tally, n >= 10) # number of grid cells with >= 10 stomachs
    winter.n10 <- nrow(winter.n10.sub)
    
    winter.n10.stomachs <- subset(winter.joined, GridID %in% winter.n10.sub$GridID)
    winter.n10.total <- nrow(winter.n10.stomachs) # number of stomachs in a grid with >= 10 stomachs
    
    # save an attribute that will be useful later
    attr(grid2, "grid") <- paste(i, j, sep = "_")
    attr(grid2, "which") <- which
    attr(summer.joined, "grid") <- paste(i, j, sep = "_")
    attr(summer.joined, "which") <- which
    attr(winter.joined, "grid") <- paste(i, j, sep = "_")
    attr(winter.joined, "which") <- which
    
    ## for plotting top grids later
    # save joined stomachs
    assign(paste("summer_joined", i, j, sep = "_"), summer.joined)
    assign(paste("winter_joined", i, j, sep = "_"), winter.joined)
    
    # save grid2
    assign(paste("grid", i, j, sep = "_"), grid2)
    
    row <- cbind(grid_id, which, x.range, y.range, 
                 summer.n, summer.n10, summer.n10.total, 
                 winter.n, winter.n10, winter.n10.total)
    grid.info <- rbind(grid.info, row)
    
    # add grid_x_y, summer_joined_x_y, winter_joined_x_y to their lists
    grid.list[[which]] <- grid2
    summer.joined.list[[which]] <- summer.joined
    winter.joined.list[[which]] <- winter.joined
  }
}

# metrics to rank grids
grid.info$summer.n <- as.numeric(grid.info$summer.n)
grid.info$winter.n <- as.numeric(grid.info$winter.n)
grid.info$summer.n10 <- as.numeric(grid.info$summer.n10)
grid.info$winter.n10 <- as.numeric(grid.info$winter.n10)
grid.info$summer.n10.total <- as.numeric(grid.info$summer.n10.total)
grid.info$winter.n10.total <- as.numeric(grid.info$winter.n10.total)

unique(grid.info$summer.n) # 993
unique(grid.info$winter.n) # 762

grid.info$mean.n10 <- (grid.info$summer.n10 + grid.info$winter.n10)/2
grid.info$mean.n10.total <- (grid.info$summer.n10.total + grid.info$winter.n10.total)/2

# alternatively grids could be ranked using the proportion of stomachs in
# each season. Since there are more stomachs in summer, this would avoid
# grid selection being biased to summer coverage. However, this does not
# affect which grid is identified as the best grid (line 314).
grid.info$summer.prop <- grid.info$summer.n10.total / grid.info$summer.n
grid.info$winter.prop <- grid.info$winter.n10.total / grid.info$winter.n
grid.info$mean.prop <- (grid.info$summer.prop + grid.info$winter.prop) / 2

# function to visualize the different grids
evaluate_grid <- function(g) {

  which <- as.numeric(grid.info$which[grid.info$grid_id == g])
  print(which)
  
  rank <- grid.info$rank[grid.info$grid_id == g]
  
  grid <- grid.list[[which]]
  summer.joined <- summer.joined.list[[which]]
  winter.joined <- winter.joined.list[[which]]
  
  i_j <- attr(grid, "grid")
  print(attr(grid, "grid"))
  print(attr(summer.joined, "grid"))
  print(attr(winter.joined, "grid"))

  summer.grid.tally <- summer.joined %>% group_by(GridID) %>% tally()
  summer.grid.tally10 <- subset(summer.grid.tally, n >= 10)
  summer.cells <- nrow(summer.grid.tally10)
  summer.stomachs <- nrow(subset(summer.joined, GridID %in% summer.grid.tally10$GridID)) # 964

  winter.grid.tally <- winter.joined %>% group_by(GridID) %>% tally()
  winter.grid.tally10 <- subset(winter.grid.tally, n >= 10)
  winter.cells <- nrow(winter.grid.tally10)
  winter.stomachs <- nrow(subset(winter.joined, GridID %in% winter.grid.tally10$GridID)) # 716
  
  grid$summer <- "no"
  grid$summer[grid$GridID %in% summer.grid.tally10$GridID] <- "yes"
  grid$winter <- "no"
  grid$winter[grid$GridID %in% winter.grid.tally10$GridID] <- "yes"
  
  grid$overall <- NA
  grid$overall[grid$summer == "no" & grid$winter == "no"] <- "no"
  grid$overall[grid$summer == "yes" & grid$winter == "no"] <- "summer only"
  grid$overall[grid$summer == "no" & grid$winter == "yes"] <- "winter only"
  grid$overall[grid$summer == "yes" & grid$winter == "yes"] <- "summer + winter"
  
  grid.title <- paste0("grid_", i_j)
  grid.subtitle <- paste0(summer.cells, " summer and ", 
                          winter.cells, " winter grid cells\n",
                          summer.stomachs, " summer and ", winter.stomachs,
                          " winter non-empty stomachs")
  
  main <- ggplot() +
    geom_sf(data = bc.coastUTM) + 
    geom_sf(data = salish.chUTM, alpha = 0.4, size = 0.8) +
    geom_sf(data = grid, aes(fill = overall)) +
    scale_fill_manual(
      values = c(`summer + winter` = alpha("#683968", 0.5), 
                `summer only` = alpha("#CA0020", 0.5),
                `winter only` = alpha("#0571B0", 0.5), 
                no = "transparent"),
      breaks = c("summer + winter", "summer only", "winter only")) +
    
    annotate("rect", xmin = 365*10^3, xmax = 495*10^3, ymin = 5325*10^3, 
             ymax = 5422.5*10^3, fill = "transparent", colour = "navy") +
    annotate("rect", xmin = 400*10^3, xmax = 500*10^3, ymin = 5425*10^3, 
             ymax = 5500*10^3, fill = "transparent", colour = "navy") +
    annotate("rect", xmin = 325*10^3, xmax = 390*10^3, ymin = 5510*10^3,
             ymax = 5560*10^3, fill = "transparent", colour = "navy") +
    coord_sf(datum = st_crs(32610),
             xlim = c(280*10^3, 600*10^3), ylim = c(5325*10^3, 5650*10^3)) +
    theme_bw() + labs(fill = "Sufficient sampling") +
    theme(legend.position = "inside", legend.position.inside = c(0.8, 0.84),
          legend.background = element_rect(fill = "grey90", colour = "black"),
          legend.title = element_text(size = 11), 
          legend.text = element_text(size = 10),
          axis.title = element_blank(), axis.text = element_blank(),
          axis.ticks = element_blank()) +
    ggtitle(grid.title, subtitle = grid.subtitle)
  
  jdf <- ggplot() +
    geom_sf(data = bc.coastUTM) + 
    geom_sf(data = salish.chUTM, alpha = 0.6, size = 0.8) +
    geom_sf(data = grid, aes(fill = overall), show.legend = FALSE) +
    scale_fill_manual(
      values = c(`summer + winter` = alpha("#683968", 0.5), 
                 `summer only` = alpha("#CA0020", 0.5),
                 `winter only` = alpha("#0571B0", 0.5), 
                 no = "transparent"),
      breaks = c("summer + winter", "summer only", "winter only")) +
    coord_sf(expand = FALSE, datum = st_crs(32610),
             xlim = c(365*10^3, 495*10^3), ylim = c(5325*10^3, 5422.5*10^3)) +
    theme_bw() + 
    theme(axis.text = element_blank(), axis.ticks = element_blank(),
          axis.title = element_blank(), plot.title = element_text(size = 8)) +
    ggtitle("Juan de Fuca Strait")

  csog <- ggplot() +
    geom_sf(data = bc.coastUTM) + 
    geom_sf(data = salish.chUTM, alpha = 0.6, size = 0.8) +
    geom_sf(data = grid, aes(fill = overall), show.legend = FALSE) +
    scale_fill_manual(
      values = c(`summer + winter` = alpha("#683968", 0.5), 
                 `summer only` = alpha("#CA0020", 0.5),
                 `winter only` = alpha("#0571B0", 0.5), 
                 no = "transparent"),
      breaks = c("summer + winter", "summer only", "winter only")) +
    coord_sf(expand = FALSE, datum = st_crs(32610),
             xlim = c(400*10^3, 500*10^3), ylim = c(5425*10^3, 5500*10^3)) +
    theme_bw() + 
    theme(axis.text = element_blank(), axis.ticks = element_blank(),
          axis.title = element_blank(), plot.title = element_text(size = 8)) +
    ggtitle("Central Strait of Georgia")
  
  nsog <- ggplot() +
    geom_sf(data = bc.coastUTM) + 
    geom_sf(data = salish.chUTM, alpha = 0.6, size = 0.8) +
    geom_sf(data = grid, aes(fill = overall), show.legend = FALSE) +
    scale_fill_manual(
      values = c(`summer + winter` = alpha("#683968", 0.5), 
                 `summer only` = alpha("#CA0020", 0.5),
                 `winter only` = alpha("#0571B0", 0.5), 
                 no = "transparent"),
      breaks = c("summer + winter", "summer only", "winter only")) +
    coord_sf(expand = FALSE, datum = st_crs(32610),
             xlim = c(325*10^3, 390*10^3), ylim = c(5510*10^3, 5560*10^3)) +
    theme_bw() + 
    theme(axis.text = element_blank(), axis.ticks = element_blank(),
          axis.title = element_blank(), plot.title = element_text(size = 8)) +
    ggtitle("Northern Strait of Georgia")
  nsog
  
  # arrange plots together
  zoom <- ggarrange(nsog, csog, jdf, nrow = 3, ncol = 1)
  grid.map <- ggarrange(main, zoom, widths = c(0.68, 0.36))

  file.name <- paste0("figures/select-grid/", rank, "-grid_", i_j, ".PNG")
  ggsave(filename = file.name, plot = grid.map, 
         width = 17, height = 13, units = "cm", bg = "white")
  return(grid.map)
}


# 4. Determine best grid ----
# Defined as highest mean number of cells with > 10 stomachs in summer and winter. 
# Ties ranked by mean total number of stomachs
grid.info <- grid.info[order(grid.info$mean.n10, grid.info$mean.n10.total, decreasing = TRUE), ]
grid.info$rank <- seq(1, 625, by = 1)
top.grids <- head(grid.info, 10) # 15_23, 12_23, 17_5, 16_5, 17_4
top.grids1 <- top.grids$grid_id
top.grids1

top.grids1 %>% purrr::map(evaluate_grid)

names(top.grids)
top.grids.table <- select(top.grids, rank, grid_id, summer.n10, summer.n10.total,
                          winter.n10, winter.n10.total, mean.n10, mean.n10.total)

# grid 15_23: not suitable, due to grid cell combining Juan de Fuca Strait and
# Saanich Inlet. Separation of Burrard Inlet stomachs into two grid cells is also
# not ideal.

# grid 12_23: not suitable, due to grid cell combining Juan de Fuca Strait and
# Saanich Inlet. This same grid cell disrupts spatially continuous sampling along
# Juan de Fuca Strait coastline. Separation of Burrard Inlet stomachs into two 
# grid cells is also not ideal.

# grid 17_5: suitable. Grid cell at Jervis Inlet is not ideal (separates 
# Malaspina Strait samples that are close together, with some being joined with
# Jervis Inlet), but acceptable. Grid cell at south Texada joining east Texada 
# and Lasqueti samples not ideal, but acceptable.

# grid 17_5 is the best grid.

# save grid_17_5 for use in analyses
saveRDS(grid_17_5, "data/out/grid_17_5.RDS") 

# Alternatively use mean proportion of stomachs in each season # no change
grid.info.prop <- grid.info[order(grid.info$mean.n10, grid.info$mean.prop, decreasing = TRUE), ]
head(grid.info.prop$grid_id, 5) 
# grid 17_5 remains the best grid, retains 96.2% of stomachs