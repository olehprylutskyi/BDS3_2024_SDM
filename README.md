# BDS^3 2024 - Species Distribution Modelling
=============================================

The repository contains materials for the projects related to the Species Distribution Modelling (SDM), carried out during the [Biological Data Science Summer School](https://www.bds3.org/), 7-20 July 2024, Uzhhorod, Ukraine.

__To be updated__

[Common mistakes in ecological niche models](https://doi.org/10.1080/13658816.2020.1798968)

## R packages for geoinformatics
[sf](https://r-spatial.github.io/sf/) - Simple Features for R

[terra](https://rspatial.github.io/terra/index.html) - modern raster and vector operations

[stars](https://r-spatial.github.io/stars/index.html) - Spatiotemporal Arrays: Raster and Vector Datacubes

[sp](http://edzer.github.io/sp/) - R Classes and Methods for Spatial Data

[raster](https://rspatial.org/raster/pkg/index.html) __outdated, not supported anymore!__

[A curated list of R packages for ecological niche modelling](https://www.sciencedirect.com/science/article/pii/S0304380022003404?via%3Dihub)

## Species occurrence data

[Global Biodiversity Inforamtion Facility](https://www.gbif.org/) - biodiversity data aggregator

[rgbif](https://docs.ropensci.org/rgbif/articles/getting_occurrence_data.html) - R package to access GBIF data

[spocc](https://docs.ropensci.org/spocc/) - Species occurrence data toolkit for R

[Selecting pseudo-absences for species distribution models: how, where and how many?](https://besjournals.onlinelibrary.wiley.com/doi/full/10.1111/j.2041-210X.2011.00172.x)

## Environmental data

[GEE Data Catalog](https://developers.google.com/earth-engine/datasets/) - the official data repository.

[GEE Community Catalog](https://developers.google.com/earth-engine/datasets/community/sat-io) - additional data.

[Sign Up](https://code.earthengine.google.com/register) to the GEE first, and [create a free account](https://youtu.be/3IwfRW8bjmo?si=_WbeUUSC7pQa_fnC).

[GEE basics tutorial](https://code.earthengine.google.com/2086c1681db9342b676b7a3243983508). See [official documentation](https://developers.google.com/earth-engine/guides) for more examples.

[Sample dataset workflow example](https://code.earthengine.google.com/ebe2c40e67fc1eb2e63af38aeca19960).

[Formules and explanations for spectral indices](https://awesome-ee-spectral-indices.readthedocs.io/en/latest/), and [GEE library for ready-to-use indices](https://ee-spectral.readthedocs.io/en/latest/).


- Bioclimatic
    * [CHELSA](https://chelsa-climate.org/)
    * [Worldclim](https://www.worldclim.org/)
- Edaphic
    * elevation
    * Slope
    * Exposition
    * Soil
- Landcover
- Primary productivity

## Data preparation in R/GEE
- Spatial thinning
- Resolution unification
- Dimension reduction
- Multicollinearity issue

## Modelling techniques
[A statistical explanation of MaxEnt for ecologist](https://onlinelibrary.wiley.com/doi/full/10.1111/j.1472-4642.2010.00725.x)

[ENMTools](https://nsojournals.onlinelibrary.wiley.com/doi/10.1111/ecog.05485)

[sdm](https://doi.org/10.1111/ecog.01881)

[ENMTML](https://doi.org/10.1016/j.envsoft.2019.104615)

[biomod2](https://cran.r-project.org/web/packages/biomod2/index.html)

## GUI tools
[QGIS Desktop](https://qgis.org/)

[Extra basemaps for QGIS](https://bnhr.xyz/2018/10/07/basemaps-in-qgis.html) - Python code to add convenient satellite/terrain/route webmaps to QGIS canva.

[MaxEnt](https://biodiversityinformatics.amnh.org/open_source/maxent/)

[Wallace Ecological Modelling App](https://wallaceecomod.github.io/)

## Prediction uncertainty
[The art of modelling range-shifting species](https://doi.org/10.1111/j.2041-210X.2010.00036.x) - Multivariate environmental similarity surface (MESS)


## Prediction binarization and thresholding


## Spatial autocorrelation (SAC)

## JointSDM
### Ensemble SDM

### HMSC
Basic paper: [Joint species distribution modelling with the r-package Hmsc](https://besjournals.onlinelibrary.wiley.com/doi/abs/10.1111/2041-210X.13345)

Important theory: [How to make more out of community data? A conceptual framework and its implementation as models and software](10.1111/ele.12757)

Hands-on: 
[EDS Seminar Series 2/22/22 - Joint Species Distribution Modeling in R with Hmsc](https://www.youtube.com/watch?v=u07eFE3Uqtg)

[Guide to using the Hmsc package for the production of Joint Species Distribtuion Models](https://www.r-bloggers.com/guide-to-using-the-hmsc-package-for-the-production-of-joint-species-distribtuion-models/)

[EDS Seminar Series 2/22/22 - Joint Species Distribution Modeling in R with Hmsc](https://www.youtube.com/watch?v=u07eFE3Uqtg)

## Reporting Species Distribution Models
[A standard protocol for reporting species distribution models](https://nsojournals.onlinelibrary.wiley.com/doi/full/10.1111/ecog.04960)

## Additional tutorials
[Species distribution modeling with terra package](https://rspatial.org/sdm/index.html#species-distribution-modeling)

[Introduction to species distribution modelling (SDM) in R](https://damariszurell.github.io/SDM-Intro/)

[ENM2020](https://www.youtube.com/watch?v=vj8qTo56rPA&list=PLCq9UxocboXPdulJteLT7MYj1WrW_tKcd) - the most comprehensive SDM/ENM video course ever. Best teachers from around the World, all possible SDM-related topics discussed. Recommend for those who already have a solid basis in SDM.