# ASDP database v. 1.1 (stomachs processed up to 2024-07-16)
# Produce data files for spatial analysis of adult Chinook salmon diet composition

## 0. Script description ----
# This script produces the data file used by analyses in Greentree et al. 
# to identify Salish Sea regions with distinct adult Chinook salmon diet 
# composition. The output is Salish-Chinook-diet.csv, which is used in analysis.
# Each row is a stomach, including capture data and proportional diet 
# composition (individual % weight, including unidentified material).

# The input files to the script include confidential information, such as 
# angler names and capture site names. As a result, we create a column called
# Site Resolution that indicates whether an angler's capture site name was
# successfully assigned to spatial coordinates or could only be assigned to 
# a fishery management area (typically also provided by the angler). 

# If interested in using Adult Salmon Diet Program data, please reach out
# (contact emails: uvicsalmondiet@gmail.com, wgreentree@outlook.com). Data
# collection is ongoing as this is a long-term monitoring program. Additional
# data may be available for your needs.


## 1. Load data/basic data wrangling ----
# load packages
library(tidyverse)
library(readxl)
library(reshape2)

# load data
nms <- names(read_excel("data/ASDP-database-v1.1/ASDP-salmondata-20240716.xlsx", n_max = 0))
nms
ct <- rep("guess", length(nms))
ct
ct[4] <- "text"
ct[5] <- "text"
ct[45] <- "text"
ct[46] <- "text"
ct[47] <- "text"
salmon.data <- read_excel("data/ASDP-database-v1.1/ASDP-salmondata-20240716.xlsx", col_types = ct)
warnings() # no warnings

stomach.data <- read_excel("data/ASDP-database-v1.1/ASDP-stomachdata-20240716.xlsx")

# A small subset of stomachs are month-only
month.only <- read_excel("data/ASDP-database-v1.1/ASDP-month-only-dates.xlsx")

# A small subset of stomachs are not appropriate for the proportional diet 
# composition analysis used in this paper.
stomach.rm <- read_excel("data/ASDP-database-v1.1/ASDP_stomachs_to_removeMarch2024.xlsx")

# change key column names/character types
salmon.data$`Fish Code` <- as.numeric(salmon.data$`Fish Code`)
stomach.data$`Fish Code` <- as.numeric(stomach.data$`Fish Code`)
stomach.rm$`Fish Code` <- as.numeric(stomach.rm$`Fish Code`)

salmon.data$Longitude <- as.numeric(salmon.data$Longitude)
salmon.data$Latitude <- as.numeric(salmon.data$Latitude)

names(salmon.data)[26] <- "PFMA"
names(stomach.data)[10] <- "Weight"
names(month.only)[3] <- "MultipleMonths"

# add column to indicate month-only date
month.only <- select(month.only, `Fish Code`, MultipleMonths)
names(month.only)
table(month.only$MultipleMonths)
month.only$MultipleMonths[month.only$MultipleMonths == 1] <- 2 
month.only$MultipleMonths[month.only$MultipleMonths == 0] <- 1
nrow(salmon.data)
salmon.data <- merge(x = salmon.data, y = month.only, by = "Fish Code", 
                     all.x = TRUE, all.y = TRUE)
nrow(salmon.data)
nrow(subset(salmon.data, is.na(MultipleMonths)))
salmon.data$MultipleMonths[is.na(salmon.data$MultipleMonths)] <- 0
table(salmon.data$MultipleMonths)

salmon.data$`Date Resolution`[salmon.data$MultipleMonths == 0] <- "exact date"
salmon.data$`Date Resolution`[salmon.data$MultipleMonths == 1] <- "date at month resolution"
salmon.data$`Date Resolution`[salmon.data$MultipleMonths == 2] <- "date at two-month resolution"
table(salmon.data$`Date Resolution`) 
# 46 stomachs were submitted by anglers only at the month resolution
# 8 winter stomachs were submitted by an angler in a mixed January and February batch

## 2. Salmon data file ----
# remove stomachs that can't be included for unique reasons
salmon.data <- subset(salmon.data, !`Fish Code` %in% c(17435, 20828))
# 17435 submitted without a data card, 20828 wasn't caught by recreational fishery

# remove stomachs of salmon sampled from trawls
salmon.data <- subset(salmon.data, Project != "Hake Bycatch" | is.na(Project))
salmon.data$Project[is.na(salmon.data$Project)] <- "Angler"

# samples from April 2017 - March 2022
salmon.data$YearMonth <- paste(salmon.data$Year, salmon.data$Month, sep = "-")
study.period <- c("2017-4", "2017-5", "2017-6", "2017-7", "2017-8", "2017-9",
                  "2017-10", "2017-11", "2017-12", "2018-1", "2018-2", "2018-3",
                  "2018-4", "2018-5", "2018-6", "2018-7", "2018-8", "2018-9", 
                  "2018-10", "2018-11", "2018-12", "2019-1", "2019-2", "2019-3",
                  "2019-4", "2019-5", "2019-6", "2019-7", "2019-8", "2019-9",
                  "2019-10", "2019-11", "2019-12", "2020-1", "2020-2", "2020-3",
                  "2020-4", "2020-5", "2020-6", "2020-7", "2020-8", "2020-9",
                  "2020-10", "2020-11", "2020-12", "2021-1", "2021-2", "2021-3",
                  "2021-4", "2021-5", "2021-6", "2021-7", "2021-8", "2021-9",
                  "2021-10", "2021-11", "2021-12", "2022-1", "2022-2", "2022-3")
salmon.data <- subset(salmon.data, YearMonth %in% study.period)
nrow(salmon.data) # 3753

# Chinook and coho only, removing pink salmon and where species unknown
table(salmon.data$Species) # 3192 Chinook, 483 coho
nrow(subset(salmon.data, is.na(Species))) # 75
salmon.data <- subset(salmon.data, Species %in% c("ch", "co"))
nrow(salmon.data) # 3675
salmon.data$Species[salmon.data$Species == "ch"] <- "Chinook"
salmon.data$Species[salmon.data$Species == "co"] <- "Coho"

# add season column
salmon.data$Season <- "Summer"
salmon.data$Season[salmon.data$Month %in% c(10, 11, 12, 1, 2, 3)] <- "Winter"

# add a region column
# remove samples without a known region
salmon.data <- subset(salmon.data, !is.na(PFMA))

salmon.data$Region <- rep(NA, nrow(salmon.data)) 
salmon.data$Region[salmon.data$PFMA %in% c(1, 101)] <- "Haida Gwaii"
salmon.data$Region[salmon.data$PFMA %in% c(13, 14, 15, 16, 17, 18, 19, 20, 28, 29)] <- "Salish Sea"
salmon.data$Region[salmon.data$PFMA %in% c(21, 23, 24, 25, 26, 27, 121, 123, 124, 125, 126)] <- "WCVI"
salmon.data$Region[salmon.data$PFMA %in% c(6, 9, 11, 12)] <- "Central Coast"
nrow(salmon.data) # 3651

# Create variable for capture site resolution. Sometimes anglers submit a stomach
# with the location recorded only to the PFMA or PFMA subarea level, or we are 
# not able to link the recorded capture site name to spatial coordinates, using
# ASDP standardized sites based on DFO creel monitoring (e.g., due to typos).
# Capture site names are considered private and not shared, but coordinates
# are shared.
creel.sites <- sort(unique(salmon.data$`Creel Site Name`))
creel.sites
area.resolution.sites <- creel.sites[1:28] # the first 28 are PFMA/PFMA subareas,
salmon.data$`Site Resolution` <- "creel site"
salmon.data$`Site Resolution`[salmon.data$`Creel Site Name` %in% area.resolution.sites] <- "PFMA or PFMA subarea"
table(salmon.data$`Site Resolution`) 
# 3356 at creel site resolution, 295 at PFMA or PFMA subarea resolution

# dates were entered as "YYYY-MM-DD", meaning that month-resolution dates
# were entered as "YYYY-MM-01" or "YYYY-MM-15". For these samples, Day should
# be NA
salmon.data$Day[salmon.data$`Date Resolution` != "exact date"] <- NA

# select columns needed
salmon.data <- select(salmon.data, `Fish Code`, Species, Season, `Date Resolution`, 
                     Year, Month, Day, Season, Region, PFMA, `Site Resolution`,
                      Longitude, Latitude, `Final Length (cm)`, Project)
nrow(salmon.data) # 3651
str(salmon.data)

## 3. Remove stomachs that don't meet analysis requirements ----
# stomachs to remove before analysis
stomach.rm <- read_excel("data/ASDP-database-v1.1/ASDP_stomachs_to_removeMarch2024.xlsx")
stomach.rm.prop <- subset(stomach.rm, `Keep for Proportion` == "no")
nrow(stomach.rm.prop) # 20

nrow(salmon.data) # 3651
length(unique(stomach.data$`Fish Code`))
salmon.data <- subset(salmon.data, !`Fish Code` %in% stomach.rm.prop$`Fish Code`)
stomach.data <- subset(stomach.data, !`Fish Code` %in% stomach.rm.prop$`Fish Code`)
nrow(salmon.data) # 3634
nrow(stomach.data) # 18,062

# weight was not measured for a small number of prey items. Replace 0 with 0.0001 g
# when those were the only prey item (i.e., % weight is 100%).
# If a prey item weight not being measured would likely change the proportion by
# more than a small amount, the stomach was removed above.
view(subset(stomach.rm, `Keep for Proportion` == "yes"))
keep.for.prop <- subset(stomach.rm, `Keep for Proportion` == "yes" & 
  `Keep for Fullness` == "no" & !`Fish Code` %in% c(19914, 20826, 20827))
# 19914: unidentified material not weighed, but still acceptable for % identified weight
# 20826, 20827 are dealt with separately - only add 0.0001 g to the Unidentified Material row
# Will be removed when stomachs that only contain Unidentified Material are removed
# prior to analysis, but they aren't empty stomachs of course

nrow(subset(keep.for.prop, `Fish Code` %in% stomach.data$`Fish Code`)) # 11
# 17015: don't add 0.0001 g to Pink salmonB since another pink salmon row with weight
# 17115: don't add 0.0001 g to SalmonChinook Salmon since another salmon row with weight
# 18705: don't add 0.0001 g to Unid FishHerring since another herring row with weight
keep.for.prop.no.weight <- subset(keep.for.prop, !`Fish Code` %in% c(17015, 17115, 18705))

# assign 0.0001 g so that proportions can be calculated
stomach.data$Weight[stomach.data$`Fish Code` %in% keep.for.prop.no.weight$`Fish Code`] <- 0.0001

# deal with 19914, 20826, 20827 separately, assign 0.0001 g to the Unidentified Material row
# so that unidentified material is dealt with appropriately. 
stomach.data$Weight[stomach.data$`Fish Code` %in% c(19914, 20826, 20827) &
                    stomach.data$FinalID == "Unidentified Material"] <- 0.0001

# add data flag to maximize usefulness for other researchers
stomach.rm$data_flag <- paste0("Keep for FO: ", stomach.rm$`Keep for FO`,
  "; Keep for Proportion: ", stomach.rm$`Keep for Proportion`, 
  "; Keep for Fullness: ", stomach.rm$`Keep for Fullness`)

stomach.rm <- select(stomach.rm, `Fish Code`, data_flag)

nrow(subset(stomach.rm, `Fish Code` %in% salmon.data$`Fish Code`)) # 11
salmon.data <- merge(x = salmon.data, y = stomach.rm, by = "Fish Code",
                     all.x = TRUE, all.y = FALSE)
nrow(salmon.data) # 3634

salmon.data$data_flag[is.na(salmon.data$data_flag)] <- "Keep for FO, Proportion, and Fullness"
table(salmon.data$data_flag)

# save salmon data
nrow(salmon.data) # 3634

## 4. Clean stomach data ----
stomach.data <- subset(stomach.data, `Fish Code` %in% salmon.data$`Fish Code`)
nrow(stomach.data) # 10,671
length(unique(stomach.data$`Fish Code`)) # 3634, good

# lowercase OriginalID, FinalID
stomach.data$OriginalID <- tolower(stomach.data$OriginalID)
stomach.data$FinalID <- tolower(stomach.data$FinalID)

# convert NA weights to 0
stomach.data$Weight[is.na(stomach.data$Weight)] <- 0

# double check there are no NAs in key variables
nrow(subset(stomach.data, is.na(`Fish Code`) | is.na(OriginalID) |
            is.na(FinalID) | is.na(`Confirmed?`) | is.na(Weight))) # 0
nrow(subset(stomach.data, is.na(`Master Code`))) # 0
subset(stomach.data, is.na(Weight)) # 0 # note that I just converted NA weights to 0

## 4A. Prey categories ----

# reduce to necessary columns
stomach.data <- select(stomach.data, `Fish Code`, FinalID, Weight)

stomach.data$Group <- rep(NA, nrow(stomach.data))

stomach.data$Group[stomach.data$FinalID %in% c("amphipod", "gammarid amphipod", "metacaprella",
                                               "hyperia medusarum", "caprellid",  
                                               "hyperiid amphipod", "themisto pacifica",
                                               "primno")] <- "Amphipod"
stomach.data$Group[stomach.data$FinalID %in% c("anchovy", "northern anchovy")] <- "Northern anchovy"
stomach.data$Group[stomach.data$FinalID %in% c("pandalus sp.")] <- "Pandalid"
stomach.data$Group[stomach.data$FinalID %in% c("octopus")] <- "Octopus"
stomach.data$Group[stomach.data$FinalID %in% c("squid", "stubby squid")] <- "Squid"

stomach.data$Group[stomach.data$FinalID %in% c("brachyuran megalopa", "brachyuran zoea",
                                               "brachyuran zoeae and megalopae", "cancrid megalopa",
                                               "cancrid zoea", "decapod", "dungeness crab megalopa",
                                               "unidentified brachyuran larvae",
                                               "pagurid megalopa")] <- "Crustacean larvae"

stomach.data$Group[stomach.data$FinalID %in% c("empty")] <- "Empty"
stomach.data$Group[stomach.data$FinalID %in% c("euphausiid")] <- "Euphausiid"
stomach.data$Group[stomach.data$FinalID %in% c("gadid", "pacific cod", "hake",
                                               "pacific hake", "walleye pollock")] <- "Gadiformes"
stomach.data$Group[stomach.data$FinalID %in% c("herring")] <- "Pacific herring"
stomach.data$Group[stomach.data$FinalID %in% c("myctophid")] <- "Myctophid"
stomach.data$Group[stomach.data$FinalID %in% c("mysid")] <- "Mysid"
stomach.data$Group[stomach.data$FinalID %in% c("algae", "gastropod", "insect", 
                                               "other", "zostera spp.")] <- "Other"
stomach.data$Group[stomach.data$FinalID %in% c("copepod")] <- "Copepod"
stomach.data$Group[stomach.data$FinalID %in% c("isopod")] <- "Isopod"
stomach.data$Group[stomach.data$FinalID %in% c("neotrypaea sp.")] <- "Neotrypaea sp."
stomach.data$Group[stomach.data$FinalID %in% c("pinnotherid")] <- "Pinnotherid"

stomach.data$Group[stomach.data$FinalID %in% c("barracudina")] <- "Barracudina"
stomach.data$Group[stomach.data$FinalID %in% c("deep-sea smelt", 
                                               "northern smoothtongue")] <- "Deep-sea smelt" 
stomach.data$Group[stomach.data$FinalID %in% c("larval flat fish", 
                                               "flatfish")] <- "Pleuronectiformes"
stomach.data$Group[stomach.data$FinalID %in% c("jack mackerel")] <- "Pacific jack mackerel" 
stomach.data$Group[stomach.data$FinalID %in% c("lingcod")] <- "Lingcod"
stomach.data$Group[stomach.data$FinalID %in% c("midshipman")] <- "Midshipman"
stomach.data$Group[stomach.data$FinalID %in% c("salmon", "pink salmon", 
                                               "chinook salmon", "chum salmon")] <- "Salmon"
stomach.data$Group[stomach.data$FinalID %in% c("poacher")] <- "Poacher"
stomach.data$Group[stomach.data$FinalID %in% c("prowfish")] <- "Prowfish"
stomach.data$Group[stomach.data$FinalID %in% c("rockfish")] <- "Rockfish"
stomach.data$Group[stomach.data$FinalID %in% c("sablefish")] <- "Sablefish"
stomach.data$Group[stomach.data$FinalID %in% c("threespine stickleback")] <- "Threespine stickleback"
stomach.data$Group[stomach.data$FinalID %in% c("tubesnout")] <- "Tubesnout"
stomach.data$Group[stomach.data$FinalID %in% c("eusergestes similis", "pasiphaea pacifica",
                                               "sergestid")] <- "Planktonic shrimp"
stomach.data$Group[stomach.data$FinalID %in% c("sand lance")] <- "Pacific sand lance"
stomach.data$Group[stomach.data$FinalID %in% c("surfperch")] <- "Surfperch"
stomach.data$Group[stomach.data$FinalID %in% c("polychaete")] <- "Polychaete"
stomach.data$Group[stomach.data$FinalID == "eulachon"] <- "Eulachon"
stomach.data$Group[stomach.data$FinalID == "surf smelt"] <- "Surf smelt"
stomach.data$Group[stomach.data$FinalID %in% c("pacific saury")] <- "Pacific saury"
stomach.data$Group[stomach.data$FinalID %in% c("staghorn sculpin")] <- "Sculpin"
stomach.data$Group[stomach.data$FinalID %in% c("crangonidae")] <- "Crangonidae"
# unid material
stomach.data$Group[stomach.data$FinalID %in% c("unidentified material", "gasterosteiformes",
                                               "brachyuran", "unidentified crustacean",
                                               "unidentified shrimp", "caridean shrimp",
                                               "unidentified fish", 
                                               "unidentified cephalopod")] <- "Unidentified material"

nrow(subset(stomach.data, is.na(Group))) # 0
view(stomach.data %>% group_by(Group, FinalID) %>% tally())
length(unique(stomach.data$Group)) # 40
# note that squid is at superorder Decapodiformes (cephalopods with 8 arms, 
# 2 tentacles), includes stubby squid

# FinalID = brachyura. It's not clear if this is a larval or settled crab, so
# it's unidentified material. If it was a settled crab, this would be in its own
# group and then not included in analyses (minor prey group, see below)
# FinalID = decapod: assumed decapod larvae since weighed less than the other larvae

## 4B. Calculate diet proportions ----
# aggregate weight by fish code and prey category
nrow(subset(stomach.data, is.na(Weight))) # 0 -> no NAs
nrow(stomach.data %>% group_by(`Fish Code`, Group) %>% tally()) # 5461
stomach.aggregate <- aggregate(stomach.data$Weight, 
                               by = list(stomach.data$`Fish Code`, stomach.data$Group),
                               FUN = sum)
nrow(stomach.aggregate) # 5461
names(stomach.aggregate) <- c("Fish Code", "Prey Category", "Weight")
length(unique(stomach.aggregate$`Fish Code`)) # 3634
stomach.aggregate <- arrange(stomach.aggregate, `Fish Code`)

# cast to wide, so there is one row per stomach
stomach.weight <- dcast(stomach.aggregate, `Fish Code` ~ `Prey Category`)
nrow(stomach.weight) # 3634
nrow(salmon.data) # 3634

stomach.weight[is.na(stomach.weight)] <- 0

last.col <- ncol(stomach.weight)
last.col # 41
stomach.weight$TotalWeight <- rowSums(stomach.weight[2:last.col])
head(stomach.weight)

# get proportions
stomach.prop <- stomach.weight
stomach.prop[2:last.col] <- stomach.prop[2:last.col] / stomach.prop$TotalWeight
head(stomach.prop)

# deals with NaN's (Not a Number) caused by dividing by zero for empty stomachs
stomach.prop <- stomach.prop %>% mutate_all(~replace(., is.nan(.), 0)) 
head(stomach.prop)

names(stomach.prop)
stomach.prop$TotalProp <- rowSums(stomach.prop[2:last.col])
unique(stomach.prop$TotalProp) # 0, 1, good

## 5. Merge salmon and stomach data ----
diet.data <- merge(x = salmon.data, y = stomach.prop, by = "Fish Code", 
                   all.x = FALSE, all.y = FALSE)
nrow(diet.data) # 3634, good

salish.ch <- subset(diet.data, Species == "Chinook" & Region == "Salish Sea")
nrow(salish.ch) # 2655

names(salish.ch)[15] <- "Data flag"

## 6. Save data ----
write_csv(salish.ch, "data/out/Salish-Sea-Chinook-salmon-diet.csv", na = "")
saveRDS(salish.ch, "data/out/Salish-Sea-Chinook-salmon-diet.RDS")

## 7. Further processing for use in grid analyses ----

## if you are running this code, read in salish.ch here, since raw data
# include confidential information
#salish.ch <- readRDS("data/out/Salish-Sea-Chinook-salmon-diet.RDS")
# if you read in the csv, convert PFMA from numeric to character to match

# remove stomachs without creel site
table(salish.ch$Season)
salish.ch <- subset(salish.ch, `Site Resolution` == "creel site")
table(salish.ch$Season) # 1479 summer, 1016 winter 
nrow(salish.ch) # 2495

# remove empty stomachs
salish.ch <- subset(salish.ch, TotalProp != 0)
salish.ch <- salish.ch[, !names(salish.ch) == "Empty"] 
nrow(salish.ch) # 1964
table(salish.ch$Season) # 1122 summer, 842 winter
# 1479-1122: 357 empty in summer
# 1016-842: 174 empty in winter

# remove stomachs with only unidentified material
salish.ch <- subset(salish.ch, `Unidentified material` != 1)
nrow(salish.ch) # 1755
table(salish.ch$Season) # 993 summer, 762 winter
# 1122-993: 129 only contained unid material in summer
# 842-762: 80 only contained unid material in winter

# % identified 
names(salish.ch)
salish.ch$IdentifiedTotal <- rowSums(salish.ch[16:53])

salish.ch[16:53] <- salish.ch[16:53] / salish.ch$IdentifiedTotal
unique(rowSums(salish.ch[16:53])) # 1

salish.ch <- salish.ch[, !names(salish.ch) %in% c("TotalWeight", "TotalProp",
                                                  "IdentifiedTotal")]

# save dataframe for select-grid.R. This dataframe is used because
# it includes summer and winter, and no stomachs contained only rare prey
# (i.e., no stomachs had to be removed).
saveRDS(salish.ch, "data/out/diet-data-for-grid-selection.RDS")

## 8. exclude rare prey, using Pacific Fishery Management Area spatial resolution ----
# I used 0.5% as the threshold, because at 1%, some stomachs only contained 
# rare prey and then would have to be removed
# Rare prey are determined separately for summer and winter

## 8A. Summer rare prey ----
summer.salish.ch <- subset(salish.ch, Season == "Summer")
summer.pfma.levels <- sort(unique(summer.salish.ch$PFMA))

names(summer.salish.ch) # don't include unidentified material
summer.pfma.prop <- data.frame(matrix(NA, nrow = 38, ncol = 0)) 
for (i in 1:length(summer.pfma.levels)) {
  summer.pfma.subset <- subset(summer.salish.ch, PFMA == summer.pfma.levels[i])
  summer.mean.diet <- colMeans(summer.pfma.subset[16:53])
  summer.pfma.prop <- cbind(summer.pfma.prop, summer.mean.diet)
  colnames(summer.pfma.prop)[i] <- summer.pfma.levels[i]
}
colSums(summer.pfma.prop) # all 1, good

# remove rare prey (<0.5% by individual mean weight)
summer.pfma.prop.major <- summer.pfma.prop[apply(summer.pfma.prop[1:ncol(summer.pfma.prop)],
                                                 1, function(x) any(x >= 0.005)), ]
summer.major.prey.vec <- rownames(summer.pfma.prop.major)
summer.major.prey.cols <- c("Fish Code", summer.major.prey.vec)
summer.major.prey.cols

summer.major.prey <- summer.salish.ch[, names(summer.salish.ch) %in% summer.major.prey.cols]
ncol(summer.major.prey) # 17
summer.major.prey$rowSum <- rowSums(summer.major.prey[2:ncol(summer.major.prey)])

summer.no.major.prey <- subset(summer.major.prey, rowSum == 0)
nrow(summer.no.major.prey) # 0 if 0.5%, 3 if 1% threshold to count as major prey

# output summer dataframe
names(summer.salish.ch)
summer.salish.ch <- select(summer.salish.ch, `Fish Code`, Species, Season, 
  `Date Resolution`, Year, Month, Day, Region, PFMA, `Site Resolution`,
  Longitude, Latitude, `Final Length (cm)`, Project, `Data flag`,
  all_of(summer.major.prey.vec))
saveRDS(summer.major.prey.vec, "data/out/summer-major-prey-vector.RDS")

# after removing rare prey, calculate % of each stomach's weight that was
# major prey
names(summer.salish.ch)
hist(rowSums(summer.salish.ch[16:31])) # good
range(rowSums(summer.salish.ch[16:31])) # fine
view(rowSums(summer.salish.ch[16:31])) # only 4 of 993 < 0.9, good

# save dataframe for grid sensitivity analysis
saveRDS(summer.salish.ch, "data/out/summer-Salish-Chinook.RDS") 

# 8B. Winter rare prey ----
winter.salish.ch <- subset(salish.ch, Season == "Winter") # 762
winter.pfma.levels <- sort(unique(winter.salish.ch$PFMA))

names(winter.salish.ch) # don't include unidentified material
winter.pfma.prop <- data.frame(matrix(NA, nrow = 38, ncol = 0)) # need to change
for (i in 1:length(winter.pfma.levels)) {
  winter.pfma.subset <- subset(winter.salish.ch, PFMA == winter.pfma.levels[i])
  winter.mean.diet <- colMeans(winter.pfma.subset[16:53])
  winter.pfma.prop <- cbind(winter.pfma.prop, winter.mean.diet)
  colnames(winter.pfma.prop)[i] <- winter.pfma.levels[i]
}
colSums(winter.pfma.prop) # all 1, good

# remove rare prey (<0.5% by individual mean weight)
winter.pfma.prop.major <- winter.pfma.prop[apply(winter.pfma.prop[1:ncol(winter.pfma.prop)],
                                                 1, function(x) any(x >= 0.005)), ]
winter.major.prey.vec <- rownames(winter.pfma.prop.major)
winter.major.prey.cols <- c("Fish Code", winter.major.prey.vec)
winter.major.prey.cols

winter.major.prey <- winter.salish.ch[, names(winter.salish.ch) %in% winter.major.prey.cols]
winter.major.prey$rowSum <- rowSums(winter.major.prey[2:ncol(winter.major.prey)])

winter.no.major.prey <- subset(winter.major.prey, rowSum == 0)
nrow(winter.no.major.prey) # 0 - good

names(winter.salish.ch)
winter.salish.ch <- select(winter.salish.ch, `Fish Code`, Species, Season, 
  `Date Resolution`, Year, Month, Day, Region, PFMA, `Site Resolution`,
  Longitude, Latitude, Project, `Final Length (cm)`, `Data flag`,
  all_of(winter.major.prey.vec))
saveRDS(winter.major.prey.vec, "data/out/winter-major-prey-vector.RDS")

# after removing rare prey, calculate % of each stomach's weight that was
# major prey
names(winter.salish.ch)
hist(rowSums(winter.salish.ch[16:33])) # good
min(rowSums(winter.salish.ch[16:33])) # 0.136
view(rowSums(winter.salish.ch[16:33])) # the second lowest is 0.9980063, so that is good

# save dataframe for grid sensitivity analysis
saveRDS(winter.salish.ch, "data/out/winter-Salish-Chinook.RDS")
