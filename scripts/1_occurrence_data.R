# Uplad data from local file
occs <- read.csv("path-to-the-data-file.csv")

# Upload data from GBIF.
# Method 1. {rgbif}
# See https://docs.ropensci.org/rgbif/articles/getting_occurrence_data.html
# There are two ways to get occurrence data from GBIF:
# 1. occ_download(): unlimited records. Useful for research and citation.
# 2. occ_search(): limited to 100K records. Useful primarily for testing.

# The function occ_search() (and related function occ_data()) should not be used for serious research. 
library(rgbif)
occs <- occ_search(scientificName = 'Ursus americanus',
                   fields=c(
                    'scientificName',
                    'basisOfRecord',
                    'decimalLatitude',
                    'decimalLongitude',
                    'occurrenceStatus',
                    'coordinateUncertaintyInMeters',
                    'year'),
                   limit = 10000)$data

# View data
View(occs)

# Method 2. {spocc}
# See https://docs.ropensci.org/spocc/
library(spocc)
# Query selected database for occurrence records
queryDb_Ss <- occs_queryDb(
  spNames = "Salamandra salamandra", 
  occDb = "gbif", 
  occNum = 1000,
  RmUncertain = TRUE)
occs_Ss <- queryDb_Ss$Salamandra_salamandra$cleaned

# Occurrence data cleaning
occs <- occs %>% 
  filter('conditions-you-specified') %>% 
  select('columns-you-need') %>% 
  rename('x' = "decimalLongitude", 'y' = "decimalLatitude")

# Spatial thinning
thinDist <- 2 # the distance (in kilometers) that you want records to be separated by.

# subset to just records with latitude and longitude
occs <- occs[!is.na(occs$x) & !is.na(occs$y),]

# round longitude and latitude with 5 digits
occs$x <- round(occs$x, 5)
occs$y <- round(occs$y, 5)

# remove spatial duplicates
occs <- unique(occs)

# Make a column with presence values (1)
occs$species <- 1


# Convert to `terra` spatVector
library(terra)
occv <- vect(occs, geom = c("x", "y"), crs = "epsg:4326") # "+proj=longlat +datum=WGS84"
plot(occv)


# Thinning the occurrences to given distance in km
occ.df <- geom(occv, df = TRUE)
occ.df$species <- 1

library(spThin)
output <- spThin::thin(loc.data = occ.df, lat.col = "y",
                       long.col = "x", spec.col = 'species',
                       thin.par = thinDist, reps = 100,
                       locs.thinned.list.return = TRUE, write.files = FALSE,
                       write.log.file = FALSE, verbose = FALSE)

# pull thinned dataset with max records, not just the first in the list
maxThin <- which(sapply(output, nrow) == max(sapply(output, nrow)))
# if more than one max, pick first
maxThin <- output[[ifelse(length(maxThin) > 1, maxThin[1], maxThin)]]
occv.th <- occv[occ.df[as.numeric(rownames(maxThin)),]]
plot(occv.th, cex = 1)