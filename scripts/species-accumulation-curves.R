## Species accumulation curves for the supplemental material

## 0. Script description ----

# Species accumulation curves are used to assess if the minimum 10 stomachs
# per grid cell decision is appropriate.

# Uses diet data before the removal of rare prey. Empty stomachs and stomachs
# without identified mass were already removed. Data output from summer-grid-
# analysis.R and winter-grid-analysis.R

library(tidyverse)
library(vegan)
library(ggrepel)

## 1. Species accumulation curves ----

## 1A. Summer ----

summer.salish.ch <- readRDS("data/out/summer-Salish-Chinook-for-specaccum.RDS")
summer.major.prey <- readRDS("data/out/summer-major-prey-vector.RDS")

# in the main text, grid cell 103 is converted to 99 to fit into some figures
summer.salish.ch$GridID[summer.salish.ch$GridID == 102] <- 98
summer.salish.ch$GridID[summer.salish.ch$GridID == 103] <- 99
# 101 to 97 is not changed because only one stomach in grid cell 101 

summer.grid.table <- data.frame(table(summer.salish.ch$GridID))
summer.grid.table <- subset(summer.grid.table, Freq != 1)
summer.grids <- sort(unique(summer.grid.table$Var1))
summer.grids <- as.character(summer.grids)
summer.grids

summer.curve <- data.frame(matrix(NA, nrow = length(summer.grids), ncol = 3))
names(summer.curve) <- c("grid", "n", "slope")

summer.all.curves <- data.frame(matrix(NA, nrow = 0, ncol = 4))
names(summer.all.curves) <- c("Stomach", "Statistic", "Species", "grid")

for (i in 1:length(summer.grids)) {
  print(summer.grids[i])
  grid <- subset(summer.salish.ch, GridID == summer.grids[i])
  matrix <- select(grid, all_of(summer.major.prey))
  
  set.seed(214)
  curve <- specaccum(matrix, method = "random", permutations = 10000)
  
  curve.df <- as.data.frame(summary(curve))
  curve.df <- separate(curve.df, Freq, c("Statistic", "Species"), sep = ":",
                       remove = FALSE)
  curve.df <- separate(curve.df, Var2, c("Stomach", NA), sep = " ", remove = FALSE)
  curve.df$Species <- as.numeric(curve.df$Species)
  curve.df$Stomach <- as.numeric(curve.df$Stomach)
  
  curve.mean <- subset(curve.df, Statistic == "Mean   ")
  
  terminus <- tail(curve.mean, n = 4)
  
  ggplot() +
    geom_line(data = curve.mean, aes(x = Stomach, y = Species))
  
  regression <- lm(Species ~ Stomach, data = terminus)
  summary(regression)
  
  slope <- data.frame(coefficient = coefficients(regression))
  slope$var <- rownames(slope)
  slope <- subset(slope, var == "Stomach")
  grid.terminus.slope <- slope$coefficient
  
  summer.curve$grid[i] <- summer.grids[i]
  summer.curve$n[i] <- nrow(grid)
  summer.curve$slope[i] <- grid.terminus.slope
  
  curve.mean$grid <- summer.grids[i]
  curve.mean <- select(curve.mean, Stomach, Statistic, Species, grid)
  summer.all.curves <- rbind(summer.all.curves, curve.mean)
}

nrow(subset(summer.curve, slope < 0.05)) # 11
ggplot() +
  geom_hline(yintercept = 0.05, linetype = "dashed") +
  geom_vline(xintercept = 10, linetype = "dashed") +
  geom_point(data = summer.curve, aes(x = n, y = slope)) +
  geom_text_repel(data = summer.curve, aes(x = n, y = slope, label = grid),
                  min.segment.length = 0) +
  scale_y_continuous(breaks = c(0, 0.05, 0.25, 0.5, 0.75, 1)) +
  scale_x_continuous(breaks = c(0, 10, 30, 60, 90)) +
  xlab("Number of stomachs") + ylab("Slope") +
  theme_bw()

summer.all.curves.label <- summer.all.curves %>% group_by(grid) %>% 
  summarize(max_stomachs = max(Stomach, na.rm = FALSE),
            max_species = max(Species, na.rm = FALSE))

nrow(summer.all.curves)
summer.all.curves <- merge(x = summer.all.curves, y = summer.curve,
                           by = "grid", all.x = TRUE, all.y = TRUE)

ggplot() +
  geom_line(data = summer.all.curves,
            aes(x = Stomach, y = Species, group = grid, colour = slope),
            linewidth = 0.7) +
  geom_text_repel(data = summer.all.curves.label,
            aes(x = max_stomachs, y = max_species, label = grid),
            min.segment.length = 0) +
  scale_colour_viridis_c(direction = 1, trans = "sqrt", limits = c(0, 1)) +
  theme_bw()

## 1B. Winter ----

winter.salish.ch <- readRDS("data/out/winter-Salish-Chinook-for-specaccum.RDS")
winter.major.prey <- readRDS("data/out/winter-major-prey-vector.RDS")

winter.grid.table <- data.frame(table(winter.salish.ch$GridID))
winter.grid.table <- subset(winter.grid.table, Freq != 1)
winter.grids <- sort(unique(winter.grid.table$Var1))
winter.grids <- as.character(winter.grids)
winter.grids

winter.curve <- data.frame(matrix(NA, nrow = length(winter.grids), ncol = 3))
names(winter.curve) <- c("grid", "n", "slope")

winter.all.curves <- data.frame(matrix(NA, nrow = 0, ncol = 4))
names(winter.all.curves) <- c("Stomach", "Statistic", "Species", "grid")

for (i in 1:length(winter.grids)) {
  print(winter.grids[i])
  grid <- subset(winter.salish.ch, GridID == winter.grids[i])
  matrix <- select(grid, all_of(winter.major.prey))
  
  set.seed(214)
  curve <- specaccum(matrix, method = "random", permutations = 10000)
  
  curve.df <- as.data.frame(summary(curve))
  curve.df <- separate(curve.df, Freq, c("Statistic", "Species"), sep = ":",
                       remove = FALSE)
  curve.df <- separate(curve.df, Var2, c("Stomach", NA), sep = " ", remove = FALSE)
  curve.df$Species <- as.numeric(curve.df$Species)
  curve.df$Stomach <- as.numeric(curve.df$Stomach)
  curve.mean <- subset(curve.df, Statistic == "Mean   ")
  
  terminus <- tail(curve.mean, n = 4)
  
  ggplot() +
    geom_line(data = curve.mean, aes(x = Stomach, y = Species))
  
  regression <- lm(Species ~ Stomach, data = terminus)
  summary(regression)
  
  slope <- data.frame(coefficient = coefficients(regression))
  slope$var <- rownames(slope)
  slope <- subset(slope, var == "Stomach")
  grid.terminus.slope <- slope$coefficient
  
  winter.curve$grid[i] <- winter.grids[i]
  winter.curve$n[i] <- nrow(grid)
  winter.curve$slope[i] <- grid.terminus.slope
  
  curve.mean$grid <- winter.grids[i]
  curve.mean <- select(curve.mean, Stomach, Statistic, Species, grid)
  winter.all.curves <- rbind(winter.all.curves, curve.mean)
}

nrow(subset(winter.curve, slope < 0.05)) # 8 of 24
ggplot() +
  geom_hline(yintercept = 0.05, linetype = "dashed") +
  geom_vline(xintercept = 10, linetype = "dashed") +
  geom_point(data = winter.curve, aes(x = n, y = slope)) +
  geom_text_repel(data = winter.curve, aes(x = n, y = slope, label = grid),
                  min.segment.length = 0) +
  xlab("Number of stomachs") + ylab("Slope") +
  scale_y_continuous(breaks = c(0, 0.05, 0.25, 0.5, 0.75, 1), limits = c(0, 1)) +
  scale_x_continuous(breaks = c(0, 10, 25, 50, 75)) +
  theme_bw() 

winter.all.curves.label <- winter.all.curves %>% group_by(grid) %>% 
  summarize(max_stomachs = max(Stomach, na.rm = FALSE),
            max_species = max(Species, na.rm = FALSE))

winter.all.curves <- merge(x = winter.all.curves, y = winter.curve,
                           by = "grid", all.x = TRUE, all.y = TRUE)
ggplot() +
  geom_line(data = winter.all.curves,
            aes(x = Stomach, y = Species, group = grid, colour = slope),
            linewidth = 0.7) +
  geom_text_repel(data = winter.all.curves.label,
                  aes(x = max_stomachs, y = max_species, label = grid),
                  min.segment.length = 0) +
  scale_colour_viridis_c(direction = 1, trans = "sqrt", limits = c(0, 1)) +
  theme_bw()

# slope vs sample size
summer.curve$season <- "Summer"
winter.curve$season <- "Winter"
both.seasons.curves <- rbind(summer.curve, winter.curve)

ggplot() +
  geom_hline(yintercept = 0.05, linetype = "dashed") +
  geom_vline(xintercept = 10, linetype = "dashed") +
  geom_point(data = both.seasons.curves, aes(x = n, y = slope, colour = season),
             show.legend = FALSE, size = 2.5) +
  geom_text_repel(data = both.seasons.curves, aes(x = n, y = slope, label = grid),
                  min.segment.length = 0) +
  scale_colour_manual(values = c("firebrick3", "dodgerblue3")) +
  xlab("Number of stomachs") + ylab("Slope") + ylim(0, 1) +
  theme_bw() + 
  theme(axis.title = element_text(size = 12),
        axis.text = element_text(size = 11),
        strip.background = element_blank(),
        strip.text.x = element_text(size = 12)) +
  facet_wrap(~ season, nrow = 2)
ggsave("figures/supplement/cumulative-prey-curves-slope.PNG", 
       width = 17, height = 17, units = "cm", dpi = 600)

# show all curves
summer.all.curves$season <- "Summer"
winter.all.curves$season <- "Winter"

summer.all.curves.label$season <- "Summer"
winter.all.curves.label$season <- "Winter"

both.seasons.all.curves <- rbind(summer.all.curves, winter.all.curves)
both.seasons.all.curves.label <- rbind(summer.all.curves.label, 
                                       winter.all.curves.label)

ggplot() +
  geom_line(data = both.seasons.all.curves,
            aes(x = Stomach, y = Species, group = grid, colour = slope),
            linewidth = 0.7) +
  geom_text_repel(data = both.seasons.all.curves.label,
                  aes(x = max_stomachs, y = max_species, label = grid),
                  min.segment.length = 0) +
  labs(colour = "Slope") +
  scale_colour_viridis_b(direction = 1, trans = "sqrt", 
      breaks = c(0, 0.05, 0.25, 0.5, 0.75, 1), limits = c(0, 1)) +
  scale_y_continuous(breaks = c(3, 6, 9, 12)) +
  theme_bw() + 
  theme(axis.title = element_text(size = 12),
        axis.text = element_text(size = 11),
        strip.background = element_blank(),
        strip.text.x = element_text(size = 12),
        legend.title = element_text(hjust = 0.5)) +
  facet_wrap(~ season, nrow = 2)
ggsave("figures/supplement/cumulative-prey-curves-all.PNG", 
       width = 17, height = 17, units = "cm", dpi = 600)
