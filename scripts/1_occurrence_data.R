rm(list = ls()) # Reset R`s brain

library(dplyr)
library(tidyr)

# Uplad data from local file
occs <- read.csv("path-to-the-data-file.csv")

# Upload data from GBIF ####
# Method 1. {rgbif}
# See https://docs.ropensci.org/rgbif/articles/getting_occurrence_data.html
# There are two ways to get occurrence data from GBIF:
# 1. occ_download(): unlimited records. Useful for research and citation.
# 2. occ_search(): limited to 100K records. Useful primarily for testing.

# The function occ_search() (and related function occ_data()) should not be used for serious research. 
library(rgbif)
# Useful for exploration/training
occs <- occ_search(scientificName = 'Salamandra salamandra',
                   fields=c(
                    'scientificName',
                    'basisOfRecord',
                    'decimalLatitude',
                    'decimalLongitude',
                    'occurrenceStatus',
                    'coordinateUncertaintyInMeters',
                    'year'),
                   limit = 1000)$data

occs

#  A Very Simple Download in professional way
# remember to set up your GBIF credentials (https://docs.ropensci.org/rgbif/articles/gbif_credentials.html)
gbif_download <- occ_download(
                            pred("hasGeospatialIssue", FALSE),
                            pred("hasCoordinate", TRUE),
                            pred("occurrenceStatus","PRESENT"), 
                            pred_not(pred_in("basisOfRecord",c("FOSSIL_SPECIMEN","LIVING_SPECIMEN"))),
                            pred("taxonKey", 2431776),
                            format = "SIMPLE_CSV"
                            )

occ_download_wait(gbif_download) # checks if download is finished

occs <- occ_download_get(gbif_download) %>%
  occ_download_import()

# View data
View(occs)

# Save data to reuse it in the next session
# As Rdata file
save(occs, file = "./data/occs_salamandra.Rdata")
# As CSV file
write.csv(occs, file = "./data/occs_salamandra.csv")


# Method 2. {spocc}
# See https://docs.ropensci.org/spocc/
install.packages("spocc", dependencies = TRUE)

library(spocc)
# Query selected database for occurrence records
queryDb_Ss <- occs_queryDb(
  spNames = "Salamandra salamandra", 
  occDb = "gbif", 
  occNum = 1000,
  RmUncertain = TRUE)
occs_Ss <- queryDb_Ss$Salamandra_salamandra$cleaned

# Load occurrence data from the local file
# from Rdata
load(file = "./data/occs_salamandra.Rdata")
# alternatively, from CSV
occs <- read.csv("./data/occs_salamandra.csv")

# Occurrence data cleaning
occs <- occs %>% 
  # filter(conditions-you-specified) %>% 
  filter(coordinateUncertaintyInMeters < 1000) %>% 
  # select(columns-you-need) %>% 
  select(decimalLatitude, decimalLongitude) %>% 
  rename(x = decimalLongitude, y = decimalLatitude)

# Spatial thinning ####

# subset to just records with latitude and longitude
occs <- occs[!is.na(occs$x) & !is.na(occs$y),]

# round longitude and latitude with 5 digits
occs$x <- round(occs$x, 5)
occs$y <- round(occs$y, 5)

# remove spatial duplicates
occs <- unique(occs)

# Convert to `terra` spatVector
library(terra)
occv <- vect(occs, geom = c("x", "y"), crs = "epsg:4326") # "+proj=longlat +datum=WGS84"
plot(occv)

# Thinning the occurrences to given distance in km
thinDist <- 20 # the distance (in kilometers) that you want records to be separated by.
occ.df <- geom(occv, df = TRUE) %>% 
  select(x, y)
# Make a column with presence values (1)
occ.df$species <- 1

library(spThin)
output <- spThin::thin(loc.data = occ.df, lat.col = "y",
                       long.col = "x", spec.col = 'species',
                       thin.par = thinDist, reps = 10,
                       write.log.file = FALSE, 
                       verbose = FALSE,
                       locs.thinned.list.return = FALSE,
                       write.files = TRUE)

# pull thinned dataset with max records, not just the first in the list
maxThin <- which(sapply(output, nrow) == max(sapply(output, nrow)))
# if more than one max, pick first
maxThin <- output[[ifelse(length(maxThin) > 1, maxThin[1], maxThin)]]
occv.th <- occv[occ.df[as.numeric(rownames(maxThin)),]]
plot(occv.th, cex = 1)

