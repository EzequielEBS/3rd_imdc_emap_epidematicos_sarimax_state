library(tidyverse)
library(lubridate)
library(aweek)

dengue <- read_csv("raw_data/dengue.csv.gz")
climate <- read_csv("processed_data/climate/climate.csv.gz")
env_vars <- read_csv("raw_data/environ_vars.csv.gz")
ocean <- read_csv("raw_data/ocean_climate_oscillations.csv.gz")
pop <- read_csv("raw_data/datasus_population_2001_2025.csv.gz")
map_regional_health <- read_csv("raw_data/map_regional_health.csv")

# Prepare data to merge
dengue <- dengue %>%
  mutate(year = year(date))
climate <- climate |> dplyr::select(-date)
ocean <- ocean %>%
  mutate(
    epiweek = as.integer(
      paste0(epiyear(date), sprintf("%02d", epiweek(date)))
    )
  ) |> 
  dplyr::select(-date)
pop <- rbind(
  pop,
  lapply(unique(pop$geocode), function(code) {
    data.frame(
      geocode = code,
      year = 2026,
      population = pop$population[pop$geocode == code & pop$year == 2025]
    )
  }) %>% bind_rows()
)

# Merge data
dengue_merged <- dengue %>%
  left_join(climate, by = c("epiweek", "geocode")) %>%
  left_join(env_vars, by = c("geocode", "uf_code")) %>%
  left_join(ocean, by = "epiweek") %>%
  left_join(pop, by = c("geocode", "year"))

# Count NAs by column
na_counts <- sapply(dengue_merged, function(x) sum(is.na(x)))

if (!dir.exists("processed_data/dengue")) {
  dir.create("processed_data/dengue")
}
write_csv(dengue_merged, "processed_data/dengue/dengue_merged.csv.gz")
