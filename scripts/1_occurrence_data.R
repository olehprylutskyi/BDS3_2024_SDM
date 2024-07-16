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