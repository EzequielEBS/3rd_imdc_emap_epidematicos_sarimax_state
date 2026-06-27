# sel_cities.r — Optional data-prep step: extract a focal set of cities.
#
# Carves out a handful of geocode-level series (one CSV per city) from the
# already-merged dengue/chikungunya tables, used only for city-level model
# diagnostics in sarimax/src/model_sel.r (the IMDC state-level submission
# itself uses processed_data/<disease>/<disease>_<uf>_agg.csv.gz from
# agg_data_uf.r, not these files). Can be run any time after
# merge_dengue.r / merge_chikungunya.r.
library(tidyverse)

# IBGE geocodes of state-capital/large cities tracked for dengue diagnostics.
cities_dengue <- c(
  2931350,
  2933307,
  2302503,
  3119401,
  3549805,
  3541406,
  1200401,
  1200203,
  1716109,
  4113700,
  4103701,
  4104808,
  5201405,
  5102637,
  5215231
)

# IBGE geocodes of state-capital/large cities tracked for chikungunya diagnostics.
cities_chikungunya <- c(
  2211001,
  2931350,
  3143302,
  3119401,
  1721000,
  1716109,
  4104808,
  4219507,
  5103403,
  5102637
)

dengue_merged <- read_csv("processed_data/dengue/dengue_merged.csv.gz")
chikungunya_merged <- read_csv("processed_data/chikungunya/chikungunya_merged.csv.gz")

lapply(cities_dengue, function(city) {
  if (!dir.exists("processed_data/dengue/sel_cities")) {
    dir.create("processed_data/dengue/sel_cities")
  }
  dengue_merged %>%
    filter(geocode == city) %>%
    write_csv(paste0("processed_data/dengue/sel_cities/dengue_", city, ".csv.gz"))
})

lapply(cities_chikungunya, function(city) {
  if (!dir.exists("processed_data/chikungunya/sel_cities")) {
    dir.create("processed_data/chikungunya/sel_cities")
  }
  chikungunya_merged %>%
    filter(geocode == city) %>%
    write_csv(paste0("processed_data/chikungunya/sel_cities/chikungunya_", city, ".csv.gz"))
})
