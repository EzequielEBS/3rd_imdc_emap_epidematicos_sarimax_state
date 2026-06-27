library(tidyverse)
library(sf)
library(geobr)

climate <- read_csv("raw_data/climate.csv.gz")

# Download all municipalities
munis <- read_municipality(year = 2020) %>%
  mutate(geocode = as.integer(code_muni))

# Get centroids
munis <- munis %>%
  mutate(centroid = st_centroid(geom),
         lon = st_coordinates(centroid)[,1],
         lat = st_coordinates(centroid)[,2])

island_codes <- c(2916104, 2605459, 2919926)

islands   <- munis %>% filter(geocode %in% island_codes)
mainland  <- munis %>% filter(!geocode %in% island_codes)

# Further restrict mainland to only geocodes present in your climate data
mainland  <- mainland %>% filter(geocode %in% unique(climate$geocode))

# st_nearest_feature finds the closest geometry
nearest_idx <- st_nearest_feature(islands, mainland)

island_to_mainland <- tibble(
  island_geocode   = islands$geocode,
  mainland_geocode = mainland$geocode[nearest_idx]
)

climate_islands <- climate %>%
  filter(geocode %in% island_to_mainland$mainland_geocode) %>%
  left_join(island_to_mainland, by = c("geocode" = "mainland_geocode")) %>%
  mutate(geocode = island_geocode) %>%
  select(-island_geocode)

climate_full <- bind_rows(climate, climate_islands)

if (!dir.exists("processed_data/climate")) {
  dir.create("processed_data/climate")
}

write_csv(climate_full, "processed_data/climate/climate.csv.gz")