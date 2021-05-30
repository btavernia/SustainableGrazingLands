# SustainableGrazingLands

This repository stores codes related to assessing the effectiveness of grazing land managing practices in Colorado.  

MODISPhenologyMetrics.R takes, cleans (e.g., interpolates missing values), and smooths time series of MODIS-based NDVI values using a statistical process outlined in: 

Chen et al. 2004. A simple method for reconstructing a high-quality NDVI time-series data set based on Savitzky-Golay
filter. Remote Sensing of Environment. 91:332-344.

Once the data are smoothed, a series of phenology and productivity metrics are quantified using NDVI values from within the growing season. The metrics were the basis for judging the impact of sustainable grazing practices on ranchland productivity.
