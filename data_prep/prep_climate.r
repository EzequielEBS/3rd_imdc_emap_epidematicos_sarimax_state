# prep_climate.r — Step 1 of the data-prep pipeline.
#
# Gap-fills climate data for island municipalities that have no weather
# station / climate-grid coverage of their own, by copying the climate
# series of their nearest mainland municipality (with data) under the
# island's own geocode. Reads raw_data/climate.csv.gz and writes the
# completed table to processed_data/climate/climate.csv.gz, which is then
# consumed by data_prep/merge_dengue.r and data_prep/merge_chikungunya.r.
#
# Run before merge_dengue.r / merge_chikungunya.r.
library(tidyverse)
library(sf)
library(geobr)

climate <- read_csv("raw_data/climate.csv.gz")

# Download municipality boundaries (2020 IBGE mesh) to get coordinates for
# every geocode, since the raw climate data has no spatial reference of its
# own beyond the geocode itself.
munis <- read_municipality(year = 2020) %>%
  mutate(geocode = as.integer(code_muni))

# Get centroids
munis <- munis %>%
  mutate(centroid = st_centroid(geom),
         lon = st_coordinates(centroid)[,1],
         lat = st_coordinates(centroid)[,2])

# Known island municipalities without their own climate series: Fernando de
# Noronha (PE), Ilhabela (SP), and Florianópolis-area Ilha de
# Santa Catarina/Ilha do Mel-type codes. Hard-coded because there is no
# generic "is an island" flag in the IBGE mesh.
island_codes <- c(2916104, 2605