## Map of the Salish Sea for main text

# Also makes maps showing previous regionalizations of the Salish Sea
# Most spatial data are not included in GitHub repository, but are 
# linked in this script or are only used to improve visualization
# (e.g., showing rivers and lakes is not necessary).

# load packages
library(tidyverse)
library(sf)
library(terra)
library(tidyterra)
library(rnaturalearth)
library(patchwork)
library(ggspatial)
library(ggnewscale)
library(ggrepel)

# load data
gebco <- readRDS("data/GEBCO/gebco.RDS") # topography file
gebcoUTM <- project(gebco, "epsg:32610") # UTM zone 10
plot(gebcoUTM)

# BC coast shapefile, derived from a shapefile from the US State Department
# (https://searchworks.stanford.edu/view/cq068zf3261)
bc.coast <- read_sf("data/BC-coast-shapefile/bc-coast.shp")
bc.coastUTM <- st_transform(bc.coast, crs = st_crs(32610))

# process hillshading for BC
gebco.clipped <- crop(gebcoUTM, bc.coastUTM, mask = TRUE)
gebco.slope <- terrain(gebco.clipped, v = "slope", unit = "radians") 
gebco.aspect <- terrain(gebco.clipped, v = "aspect", unit = "radians")
gebco.hillshade <- shade(slope = gebco.slope, aspect = gebco.aspect,
                         angle = 45, direction = 315)

# relevant rivers and lakes, provided by Pacific Salmon Foundation
psf.rivers <- read_sf("data/PSF-rivers/RL_Merge.shp")
psf.riversUTM <- st_transform(psf.rivers, crs = st_crs(32610))

# DFO Pacific Fishery Management Areas
pfma <- read_sf("data/DFO-PFMAs/DFO_BC_PFMA_SUBAREAS_CHS_V3_G.shp")
salish.sea.pfma <- subset(pfma, MGNT_AREA %in% c(13, 14, 15, 16, 17, 18, 19, 20,
                                                 28, 29))
salish.sea.pfmaUTM <- st_transform(salish.sea.pfma, crs = st_crs(32610))

## Figure 1: Study area map ----
salish.sea <- ggplot() +
  geom_sf(data = salish.sea.pfmaUTM, fill = "#CEF1F8", colour = "#CEF1F8") +
  geom_sf(data = bc.coastUTM) +
  geom_spatraster(data = gebco.hillshade, alpha = 0.35, show.legend = FALSE,
                  maxcell = 5e+15) +
  # set maxcell high, so that the raster is kept as original
  scale_fill_distiller(palette = "Greys", na.value = "transparent") +
  geom_sf(data = bc.coastUTM, fill = "transparent", linewidth = 0.05) +
  
  geom_sf(data = psf.riversUTM, fill = "#CEF1F8") +
  
  annotate("text",  x = 321*10^3, y = 5479*10^3, label = "Vancouver Island",
           size = 3.5) +
  annotate("text",  x = 548*10^3, y = 5455*10^3, label = "Fraser River",
           size = 3.5) +
  annotate("text",  x = 445*10^3, y = 5320*10^3, label = "Olympic Peninsula",
           size = 3.5) +
  annotate("text",  x = 595*10^3, y = 5435*10^3, label = "Canada",
           size = 3.5, fontface = "bold") +
  annotate("text",  x = 595*10^3, y = 5420*10^3, label = "USA",
           size = 3.5, fontface = "bold") +
  
  annotate("text", x = 483*10^3, y = 5378*10^3, label = "Haro\n      Strait", size = 3.5,
           colour = "dodgerblue3", fontface = "italic", lineheight = 0.75, vjust = 1) +
  annotate("text",  x = 448*10^3, y = 5460*10^3, label = "Strait of\nGeorgia",
           fontface = "italic", colour = "dodgerblue3", lineheight = 0.75, 
           size = 3.5) +
  annotate("text",  x = 484*10^3, y = 5496*10^3, label = "Howe Sound",
           fontface = "italic", colour = "dodgerblue3", size = 3.5) +
    annotate("text",  x = 405*10^3, y = 5545*10^3, label = "Mainland\ninlets",
             size = 3.5, colour = "dodgerblue3", lineheight = 0.75,
             fontface = "italic") +
  annotate("text",  x = 456*10^3, y = 5346*10^3, label = "Juan de Fuca Strait",
           fontface = "italic", colour = "dodgerblue3", size = 3.5) +
  annotate("text",  x = 305*10^3, y = 5365*10^3, label = "Pacific Ocean",
           fontface = "italic", colour = "dodgerblue3", size = 3.5) +
  
  annotation_scale(width_hint = 0.25, style = "bar") + 
  annotation_north_arrow(location = "br", 
                         height = unit(0.9, "cm"), width = unit(0.6, "cm"),
                         which_north = "true") +
  
  coord_sf(datum = 4326,
           crs = st_crs(32610),
           xlim = c(280*10^3, 600*10^3), ylim = c(5325*10^3, 5650*10^3)) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.title = element_blank(),
        axis.text = element_text(size = 10))
salish.sea

# inset
n.america <- ne_countries(scale = "large", returnclass = "sf")
n.americaUTM <- st_transform(n.america, crs = st_crs(32610))

inset <- ggplot() +
  geom_sf(data = n.americaUTM) +
  annotate(geom = "rect", xmin = 280*10^3, xmax = 600*10^3, 
           ymin = 5325*10^3, ymax = 5650*10^3, 
           colour = "black", fill = "transparent", linewidth = 0.8) +
  coord_sf(datum = 4326, crs = st_crs(32610),
           xlim = c(-2500*10^3, 1200*10^3), ylim = c(4000*10^3, 8200*10^3)) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.title = element_blank(), 
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "none",
        plot.background = element_blank(),
        plot.margin = unit(c(0, 0, 0, 0), "cm"),
        axis.ticks.length.x = unit(0,'pt'))
inset

study.area <- salish.sea + inset_element(inset, align_to = "panel",
                           left = 0.61, bottom = 0.57, right = 1, top = 1)
ggsave("figures/maps/study-area-map.PNG", 
       plot = study.area, width = 17, height = 17, units = "cm", dpi = 600)

## Supplemental figures ----

# DFO Pacific Fishery Management Areas:
salish.sea.pfmaUTM$PFMA <- as.character(salish.sea.pfmaUTM$MGNT_AREA)

pfma.centroid <- st_centroid(salish.sea.pfmaUTM)
pfma.centroid.coords <- do.call(rbind, st_geometry(pfma.centroid)) %>% 
  as_tibble() %>% setNames(c("utm_x", "utm_y"))
pfma.centroid <- bind_cols(pfma.centroid, pfma.centroid.coords)
pfma.centroid <- st_drop_geometry(pfma.centroid)

pfma.centroid.mean <- pfma.centroid %>% group_by(PFMA) %>% 
  summarize(mean_x = mean(utm_x),
            mean_y = mean(utm_y))

# define colour palette
pfma.pal <- c("#4e8397", "#a65a7b", "#0089ba", "#CBAF0D", "#7e6e9f", 
              "#a16807", "#00896f", "#78814a", "#4f50c3", "#0e4d5f")

pfma.map <- ggplot() +
  geom_sf(data = salish.sea.pfmaUTM, aes(fill = PFMA)) +
  geom_sf(data = bc.coastUTM) +

  scale_fill_manual(values = pfma.pal) +
  
  ggnewscale::new_scale_fill() +
  geom_spatraster(data = gebco.hillshade, alpha = 0.35, show.legend = FALSE,
                 maxcell = 5e+15) +
  # set maxcell high, so that the raster is kept as original
  scale_fill_distiller(palette = "Greys", na.value = "transparent") +
  
  geom_sf(data = bc.coastUTM, fill = "transparent", linewidth = 0.05,
          colour = "black") +
  geom_sf(data = psf.riversUTM, fill = "#CEF1F8") +
  geom_label(data = pfma.centroid.mean,
             aes(x = mean_x, y = mean_y, label = PFMA, colour = PFMA)) +
  scale_colour_manual(values = pfma.pal) +
  
  annotation_scale(width_hint = 0.25, style = "bar") + 
  annotation_north_arrow(location = "br", 
                         height = unit(0.9, "cm"), width = unit(0.6, "cm"),
                         which_north = "true") +
  
  coord_sf(datum = 4326,
           crs = st_crs(32610),
           xlim = c(280*10^3, 600*10^3), ylim = c(5325*10^3, 5650*10^3)) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.title = element_blank(),
        axis.text = element_text(size = 10),
        legend.position = "none")

pfma.w.inset <- pfma.map + inset_element(inset, align_to = "panel",
                         left = 0.61, bottom = 0.57, right = 1, top = 1)
ggsave("figures/maps/PFMA-map.PNG", 
       plot = pfma.w.inset, width = 17, height = 17, units = "cm", dpi = 600)

# DFO zooplankton regions: 
# downloaded from PSF Marine Data Centre
# https://soggy2.zoology.ubc.ca/geonetwork/srv/eng/catalog.search#/metadata/e1da07ca-353d-4b8a-a08b-643885e89e3b
dfo.zoop <- read_sf("data/DFO-zooplankton-regions/ssea_regions.shp")
dfo.zoop <- st_transform(dfo.zoop, crs = st_crs(32610)) # UTM zone 10

dfo.zoop <- st_make_valid(dfo.zoop)
zoop.centroid <- st_centroid(dfo.zoop)
zoop.centroid.coords <- do.call(rbind, st_geometry(zoop.centroid)) %>% 
  as_tibble() %>% setNames(c("utm_x", "utm_y"))
zoop.centroid <- bind_cols(zoop.centroid, zoop.centroid.coords)
zoop.centroid <- st_drop_geometry(zoop.centroid)

zoop.centroid$utm_x[zoop.centroid$region == "PugetSound"] <- 520*10^3
zoop.centroid$utm_y[zoop.centroid$region == "PugetSound"] <- 5346*10^3

zoop <- ggplot() +
  geom_sf(data = dfo.zoop, aes(fill = region), show.legend = FALSE) +
  geom_sf(data = bc.coastUTM) +
  new_scale_fill() +
  geom_spatraster(data = gebco.hillshade, alpha = 0.35, show.legend = FALSE,
                  maxcell = 5e+15) +
  # set maxcell high, so that the raster is kept as original
  scale_fill_distiller(palette = "Greys", na.value = "transparent") +
  
  geom_sf(data = bc.coastUTM, fill = "transparent", linewidth = 0.05) +
  geom_sf(data = psf.riversUTM, fill = "#CEF1F8") +
  
  geom_label_repel(data = zoop.centroid, 
                   aes(x = utm_x, y = utm_y, label = region, colour = region),
                   show.legend = FALSE, fill = alpha("white", 0.8)) +
  
  annotation_scale(width_hint = 0.25, style = "bar") + 
  annotation_north_arrow(location = "br", 
                         height = unit(0.9, "cm"), width = unit(0.6, "cm"),
                         which_north = "true") +
  
  coord_sf(datum = 4326,
           crs = st_crs(32610),
           xlim = c(280*10^3, 600*10^3), ylim = c(5325*10^3, 5650*10^3)) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.title = element_blank(),
        axis.text = element_text(size = 10))
zoop.map <- zoop + inset_element(inset, align_to = "panel",
                   left = 0.61, bottom = 0.57, right = 1, top = 1)
ggsave("figures/maps/DFO-zooplankton-map.PNG", 
       plot = zoop.map, width = 17, height = 17, units = "cm", dpi = 600)

# Oceanographic basins: drawn based on Pena et al. 2016
broad.regions <- read_sf("data/oceanographic-basins/Salish-Sea-basins.shp")
broad.regionsUTM <- st_transform(broad.regions, crs = st_crs(32610))

# should clip to the Canadian border
broad.regions <- st_intersection(broad.regionsUTM, salish.sea.pfmaUTM)

basins <- ggplot() +
  geom_sf(data = broad.regions, aes(fill = Region, colour = Region),
          show.legend = FALSE) +
  geom_sf(data = bc.coastUTM) +
  scale_fill_manual(values = c("#0072b2", "#e69f00", "#009e73", "#cc79a7", "#56b4e9"))+
  scale_colour_manual(values = c("#0072b2", "#e69f00", "#009e73", "#cc79a7", "#56b4e9")) +
  
  new_scale_fill() +
  geom_spatraster(data = gebco.hillshade, alpha = 0.35, show.legend = FALSE,
                  maxcell = 5e+15) +
  # set maxcell high, so that the raster is kept as original
  scale_fill_distiller(palette = "Greys", na.value = "transparent") +
  
  geom_sf(data = bc.coastUTM, fill = "transparent", linewidth = 0.05) +
  geom_sf(data = psf.riversUTM, fill = "#CEF1F8") +
  
  annotate("label", x = 315*10^3, y = 5520*10^3, label = "Northern Strait\nof Georgia",
           fontface = "bold", colour = "#cc79a7", lineheight = 0.75, size = 3.5,
           fill = alpha("white", 0.8)) +
  annotate("label", x = 390*10^3, y = 5452*10^3, label = "Central Strait\nof Georgia",
           fontface = "bold", colour = "#0072b2", lineheight = 0.75, size = 3.5,
           fill = alpha("white", 0.8)) +
  annotate("label", x = 505*10^3, y = 5420*10^3, label = "Southern Strait\nof Georgia",
           fontface = "bold", colour = "#56b4e9", lineheight = 0.75, size = 3.5,
           fill = alpha("white", 0.8)) +
  annotate("label", x = 445*10^3, y = 5376*10^3, label = "Haro Strait",
           fontface = "bold", colour = "#e69f00", lineheight = 0.75, size = 3.5,
           fill = alpha("white", 0.8)) +
  annotate("label", x = 400*10^3, y = 5345*10^3, label = "Juan de Fuca Strait",
           fontface = "bold", colour = "#009e73", lineheight = 0.75, size = 3.5,
           fill = alpha("white", 0.8)) +
  
  coord_sf(datum = 4326,
           crs = st_crs(32610),
           xlim = c(280*10^3, 600*10^3), ylim = c(5325*10^3, 5650*10^3)) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.title = element_blank(),
        axis.text = element_text(size = 10))
basins.map <- basins + inset_element(inset, align_to = "panel",
                     left = 0.61, bottom = 0.57, right = 1, top = 1)
ggsave("figures/maps/Salish-Sea-basins.PNG", 
       plot = basins.map, width = 17, height = 17, units = "cm", dpi = 600)
