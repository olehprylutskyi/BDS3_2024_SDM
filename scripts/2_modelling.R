rm(list = ls()) # Reset R`s brain

# Set working directory (change it to yours)
setwd("~/git/BDS3_2024_SDM")

# Load libraries
library(tidyverse) # for data manipulations
library(sdm)       # for Species Distribution Modelling
library(terra)     # for raster operations (modern approach, use it when possible)


# We will use data on occurrences of Acer negundo L. - invasive tree species originated from N. America.
# Read presence data from files
presence <- read.csv("./data_raw/sampled_pixels.csv") %>%   # read data from CSV file
  select(bio02, bio05, bio06, bio12, bio15, dem, X, Y) %>%  # select variables with environmental data
  mutate(pb = 1) %>%                                        # make a new variable with "1" for all rows (indicates presence points)
  drop_na()                                                 # remove NAs

# Read background data from file (reduced to 300K points to meet GitHub file size limits)
background <- read.csv("./data_raw/background_pixels_reduced.csv") %>% 
  select(bio02, bio05, bio06, bio12, bio15, dem, X, Y) %>% 
  mutate(pb = 0) %>%                                       # make a new variable with "0" for all rows (indicates background points)
  drop_na()

# Randomly decrease the size of the background points if you face a 
# memory limits while working with large data
# background <- background[sample(nrow(background),size = 300000, replace = FALSE), ]

# Partitioning
# Then we need to randomly split the presence data on training (used for modelling) and testing 
# (used for testinmg model's predictive power) datasets
# make this example reproducible
set.seed(1) # the number might be any

# Use 70% of dataset as training set and 30% as testing set
sample <- sample(c(TRUE, FALSE), nrow(presence), replace = TRUE, prob = c(0.7,0.3))
train  <- presence[sample, ]
test   <- presence[!sample, ]

# Merge training data with background points, to use for training the model
sdmdata <- rbind(train, background)

# Modelling
m1 <- glm(
  # formula. . indicates all variables except "pb" are used as predictors
  pb ~ bio02 + bio05 + bio06 + bio12 + bio15 + dem,
  family = binomial, # since our response is binary (1 or 0)
  data = sdmdata     # data containing responce and predictor variables
)

# Check model class
class(m1)
# Check model output (evaluate explanatory side of the model as usual linear model)
summary(m1)
# We can see that all predictors are significant (p-value less than 0.05), except bio06.
# Effects themselves seem not so (very low values close to 0).
# Since we have a quite simple model, it would be to early to discuss this output without 
# at least careful checking of model assumptions.

# Prediction over the map
# Read raster predictors from the multi-band TIFF file
predictors <- terra::rast("./data_raw/combined_image.tif")
# Check predictors (must be of the same number and names as predictors from the model!)
names(predictors)

# Predict over the raster environmental data
p <- predict(predictors, # raster predictors
  m1,                    # model object
  type = "response"      # to make presence probability as the value between 0 and 1
)

# Plot the prediction
plot(p, main = "Acer negundo, probability of presence")

# Export the map
png("./figures/sdm_plot.png", width = 297, height = 210, units = "mm", res = 300)
plot(p, main = "Acer negundo, probability of presence")
dev.off()

# Evaluation
# To evaluate model's predictive performance, we will use validation (testing) data 
# set containing true presence points that were not used to train the model.
# We combine testing points with background points
test_w_env <- rbind(test, background)
# Then predict presence probabilities based on our model
test_w_env$predicted_values <- predict(
  m1,                                    # model object
  test_w_env,                            # testing points plus background
  type = "response"                      # to aquire probability of presence (from 0 to 1)
)


# Dedicated packages for SDM have in-built model evaluation toolkits.
# We will apply {sdm} (https://www.biogeoinformatics.org/) package for the example below.
library(sdm)

# Evaluation is resourse-demanding process, so below we first get a 
# smaller subset (30K points) to reduce computation time
test_w_env_subset <- test_w_env[sample(nrow(test_w_env), size = 30000), ]

ev <- evaluates(test_w_env_subset$pb, test_w_env_subset$predicted_values)

# Let's look at all the evaluation metrics
ev@statistics

# $Prevalence # the proportion of locations that are occupied
# [1] 0.007

# $AUC # Area under the curve, metric of sensitivity
# [1] 0.857

# $COR # Metric of specificity
#          cor      p.value 
# 1.210000e-01 2.075376e-98 

# $Deviance
# [1] 0.07274188

# We can access selected metric only. AUC for predictions
ev@statistics$AUC
AUC <- round(ev@statistics$AUC, 2) # usually rounded to 2nd decimal digit after the point
AUC

# To convert continuous presence probability map into a binary presence-ansence map,
# we first need to define a threshold of probability values.
# Pixels containing values above the threshold will be considered as "presence",
# below - as "absence".

# Get threshold for presence-absence (max(sensitivity + specificity), maximize TSS) binarization
th <- ev@threshold_based$threshold[2]

## Binary prediction maps
# Make empty raster with same properties as prediction raster
library(raster)
pa1 <- raster(p)

# Write binary response to the new raster
pa1[] <- ifelse(p[] >= th, 1, 0)

# Plot binary map
plot(pa1, main = "Acer negundo, modelled area of occupancy")



## Extrapolation ####
# To examine in which areas we can trust our model's prediction more,
# we have different technics. One called
# Multivariate environmental similarity surface (MESS)
# It reports the closeness of the point aty the map to the distribution of 
# reference points, gives negative values for dissimilar points and maps these 
# values across the whole prediction region 
# (Elith et al., 2010: https://doi.org/10.1111/j.2041-210X.2010.00036.x )

# Lower (negative) values indicate more dissimilar conditions

# First we need to prepare the data, and the data from known points must contains coordinates
# all_points <- data.frame(species = sdmdata$pb, coordinates(d))
# all_points <- all_points[,2:3]

# reference_points <- extract(covs, all_points)

library(dismo)
# Prepare data for MESS analysis
covariates <- raster::stack(predictors)
covariates <- subset(covariates, c("bio02", "bio05", "bio06", "bio12", "bio15", "dem"))
reference_points <- sdmdata[, c("bio02", "bio05", "bio06", "bio12", "bio15", "dem")]

mss <- dismo::mess(
  x = covariates,
  v = reference_points,
  # full = TRUE, # returns a RasterBrick w/ n layers corresponding to the layers of the input Raster
  full = FALSE
)

# Drop Inf values
mss[mss == Inf] <- NA
# Convert to terra spatRaster
mss <- rast(mss)

# Plot environmental similarity
plot(mss, main = "Multivariate environmental similarity surface (MESS)")





# It is a lot of dedicated R packages for SDM out there.
# I suggest giving a try to {sdm}, {ENMTools}, {ENMeval}, {biomod2}, {ENMTML}.

# Here we are trying {sdm} package approach, as an example.
# See videotutorial for the full tutorial:
# https://youtu.be/83dMS3bcjJM?si=SRbcKF8gT0k8D8Z4
library(sdm)

# To make model fitting faster, we selected only first 100 presence points.
# Background points will be selected randomly within predictor rasters
# Usually the first step is to prepare a specific data object  containig both 
# presence and background data, as well as data on covariates.

# Small subset of presence points
occurrences <- presence[1:100,]
# Convert to the {SpatialPointDataFrame} using {sp} package
library(sp)
coordinates(occurrences) <- c('X', 'Y') # define coordinates
# Assign coordinate reference system (WGS84 for decimal degrees)
proj4string(occurrences) = CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")
# Check it out
occurrences

# Plot training points upon one of the predictor raster
plot(predictors$dem)
plot(occurrences, add = T)

# {sdm} still applies {raster} package, so first we need to convert predictors in Raster Stack format
covs <- raster::stack(predictors)
# Select only those predictors used in the model building
covs <- subset(covs, c("bio02", "bio05", "bio06", "bio12", "bio15", "dem"))
# Technical detais
covs
# Preview predictors on maps
plot(covs)

# Construct sdmData object
sdmdata <- sdm::sdmData(
  formula = pb ~ bio02 + bio05 + bio06 + bio12 + bio15+ dem, # formula of response variable ~ predictors
  train = occurrences,      # training occurrence data
  predictors = covs,           # raster covariates (as RasterStack)
  bg = list(             # generating background
    method = "gRandom",  # selected randomly in geographic space
    remove = FALSE,      # whether points located in presence sites should be removed
    n = 10000)            # number of background records
)

# The next step is to cinstruct the model
m <- sdm::sdm(
  formula = pb ~ bio02 + bio05 + bio06 + bio12 + bio15+ dem,   # formula
  data = sdmdata,                         # sdm data object created with {sdm} function
  methods = c("glm", "maxent"),           # alternative modelling methods
  # replication = c("sub", "boot", "cv"), # method(s) for replications
  replications = "sub",
  test.p = 20,                            # test percentage for "sub", in %
  n = 2,                                  # replication for subsampling/bootstraping
  parallelSetting = list(ncore = 4,
    method = 'parallel')                  # parallel computing
)

# Since model fitting might take a while, it would be wise to save the model
save(m, file = paste0("./models/sdmModel.Rdata"))
# Load the model again
load(file = paste0("./models/sdmModel.Rdata"))
gc() # small useful command for memory cleaning

# We can explore model's preformance using a GUI tool
gui(m)
# May need manually stopping the process to proceed with the code

n_models <- 4 # specify the total number of models fitted

# Extract AUC and TSS vales for all models separately
modelEvaluations <- getEvaluation(
  m,                       # model object
  wtest = "test.dep",      # testing data
  stat = c("AUC", "TSS")   # statistics
)

# Look at AUC and TSS (True Skill Statistics) for all fitted models
modelEvaluations

# Variable importance ####
# You can get Mean variable importance (and confidence interval) for multiple models:
library(ggplot2)
vi <- getVarImp(m, # model object
                id = 1:n_models, # specify the modelIDs of the models
                wtest = 'test.dep')

vi

# Alternatively, with method specified
plot(getVarImp(m, method = "maxent"))

p_vi <- plot(vi, 'cor')
p_vi

# Prediction
# Predict the model over the raster stack
# The sdm::predict fiction write files in the working directory, be aware
p1 <- predict(
  m,                # model object
  newdata = covs,   # covariates to predict over
  overwrite = TRUE  # overwrite local files
)

# Prediction is a raster brick (multi-layer raster), made by {raster} package.
# One layer per model
p1
# Lets have a look
plot(p1)

# To get a single consensus map, we can apply ensembling.
# It is not just averaging individual predictions. We shell define weighting
# of each prediction, as well as the statisticts used for the binarization
en1 <- ensemble(
  m,                 # model object
  p1,                # prediction(s)
  setting = list(    # additional settings, see ?sdm::ensemble for details
    method = 'weighted',
    stat = 'tss', opt = 2
  )
)

# Plot ensemble prediction and save it
plot(en1)

# Apply the custom color palette and plot it again
cl <- colorRampPalette(c('#3E49BB','#3498DB','yellow','orange','red','darkred'))
plot(en1, col = cl(200), main = "Ensemble prediction (glm plus maxent)")

png("./figures/ensemble_prediction.png",
    width = 220, height = 140, units = "mm", res = 300)
plot(en1, col = cl(200), main = "Ensemble prediction (glm plus maxent)")
dev.off()

# Write ensemble prediction to a *.tiff file, if we are going to use it
# outside R
writeRaster(
  en1,                               # ensemble prediction as a raster
  filename = "./figures/ensemble_prediction.tif",           # filename
  overwrite = TRUE      # allow to overwrite the file if there was one
)

# Uncertainty estimation ####
# Binarization of the model
# Custom function to get threshold for each model and convert continual (from 0 
# to 1) prediction to the presence-absence binary raster (either 0 or 1).

# Use maximizing TSS to choose threshold (it is plenty of other methods of
# thresholding around there). Execute the code below to add the function 
# to the Environment.

binarize <- function(sdmMmodels, prediction) {
  # Find optimal threshold to particular model
  # Maximizing TSS
  th <- sdmMmodels[[1]]@models$pb[[1]]$`1`@evaluation$test.dep@threshold_based$threshold[2]
  # Make empty raster with same properties as prediction raster
  pa <- raster(prediction[[1]])
  # Write binary response to the new raster
  pa[] <- ifelse(prediction[[1]][] >= th, 1, 0)
  pa
}

# Create an empty list to write the binary map into
r_list <- list()

# Convert predictions to presence-absence rasters and write them in a list
for (i in 1:n_models) {
  pa <- binarize(m[[i]], p1[[i]])
  r_list[[i]] <- pa
}

# Convert a list to the raster stack
r_stack <- stack(r_list)
# Plot binarized prediction for each model (presence-absence)
plot(r_stack)


# Binary entropy
# The measure of variation in the expected values of the variable. Contrary to 
# the ususal entropy, vary from 0 (all models predict the same) to 1 (half of the 
# models predicts different outcome than another half).
# Read more: https://en.wikipedia.org/wiki/Binary_entropy_function
# Or watch: https://youtu.be/YtebGVx-Fxw?si=o8ylg___wCSvNGTB (orogoanl entropy conception)
# Code for calculation: https://stackoverflow.com/questions/27254550/calculating-entropy
# Formula: H(X) = −p log2 p − (1 − p) log2(1 − p)

# Custom function to calculate binary entropy
fun <- function(x) {
  # Frequency of outputs
  probabilities <- prop.table(table(x))
  # Binary entropy
  entropy <- -sum(probabilities * log2(probabilities))
  return(entropy)
}

# Calculate uncertainty (binary entropy) among binary predictions
un <- calc(r_stack, fun = fun)

# As the function return raster that covers all the extent, it might be useful 
# to mask pixels outside the area of interest (e.g., terrestrial areas or 
# specific region)

# Download World polygon from {spData} package
library(spData)
# Since {sdm} works with {sp}, convert the polygons to the Spatial Polygons
world_sp <- sp::as(world, "Spatial")

# Mask the area outside the spatial polygon
un_masked <- raster::mask(
  un,         # binary entropy raster
  world_sp    # Spatial Polygon (sp) of the mask
)

# Plot binary entropy
plot(un_masked, main = "Prediction uncertainty (binary entropy)")

# Save uncertainty plot as an image
png("./figures/uncertainty_binary_entropy.png",
    width = 220, height = 140, units = "mm", res = 300)
plot(un_masked, main = "Prediction uncertainty (binary entropy)")
dev.off()

# Prediction in environmental space ####
# Extract predictors' values for presence points from sdmData object, 
# generated with {sdm::sdmData}
occ.df <- as.data.frame(sdmdata)
occ.df <- occ.df[occ.df$pb == 1, ]
occ.df <- occ.df[, names(occ.df) != c("rID", "pb")]
str(occ.df)

# Visualize species probability of presence (favourable conditions) as a heatmap
# in the gradient of two selected environmental variables
p_bio02.bio12 <- niche(
  covs,                    # covariates
  en1,                     # ensemble prediction
  n = c("bio02", "bio12"), # choose any pair of covariates
  col = cl(200),           # colour palette
  plot = FALSE,
  out = TRUE
)

# Plot species' ecological niche based on the prediction, in 2-dimensional space of 
# 2 choosen environmental factors (covariates)
plot(p_bio02.bio12, col = cl(200))

# We shall change covariate combinations to explore species's niche
p_bio12.dem <- niche(
  covs,
  en1,
  n = c("bio12", "dem"),
  col = cl(200),
  plot = FALSE,
  out = TRUE
)

plot(p_bio12.dem, col = cl(200))

# This way we can explore how evenly (fundamental) species niche distributed
# across the environmental space (in contrast to the distribution of species 
# across geographical space that we can observe on the prediction map)

# Binary maps of ensemble prediction ####
## Define the threshold
# Extract coordinates from sdmData object
df <- as.data.frame(sdmdata)
coord <- coordinates(sdmdata)
df <- data.frame(species = df$pb, coordinates(sdmdata))
xy <- as.matrix(df[, c("X", "Y")])


# # Extract presence-background data as a vector
# species_presence <- c(
#   sdmdata@species$pb@presence,
#   sdmdata@species$pb@background
# )

# # Extract points coordinates as a data frame
# xy <- rbind(
#   presence[1:100, c("X", "Y")],
#   background[1:1000, c("X", "Y")]
# )


# Extract predicted values for data points
predicted_values <- raster::extract(en1, xy)

# Make prediction evaluation
ev <- evaluates(df$species, predicted_values)

# Alternativaly, we can use the sdmData object and the raster of ensemble prediction
ev <- evaluates(x = sdmdata, p = en1)

# Show statistics
ev@statistics

# Get threshold for presence-absence (max(se+sp), maximize TSS)
th <- ev@threshold_based$threshold[2]

## Binary ensemble maps ####
# Make empty raster with same properties as prediction raster
pa1 <- raster(en1)
# Write binary response to the new raster
pa1[] <- ifelse(en1[] >= th, 1, 0)

# Plot the binary ensemble prediction map
plot(pa1, main = "Expected area of occupancy for Acer negundo")
