library(zoo)
library(tidyverse)
library(sf)
library(geobr)

dengue_merged <- read_csv("processed_data/dengue/dengue_merged.csv.gz")
chikungunya_merged <- read_csv("processed_data/chikungunya/chikungunya_merged.csv.gz")

dengue_merged <- dengue_merged %>%
  arrange(epiweek) %>%
  mutate(across(c(enso, iod, pdo), ~ na.approx(., na.rm = FALSE)))
chikungunya_merged <- chikungunya_merged %>%
  arrange(epiweek) %>%
  mutate(across(c(enso, iod, pdo), ~ na.approx(., na.rm = FALSE)))

# Count NAs by column
na_counts <- sapply(dengue_merged, function(x) sum(is.na(x)))
na_counts_chik <- sapply(chikungunya_merged, function(x) sum(is.na(x)))

if (!dir.exists("processed_data/dengue")) {
  dir.create("processed_data/dengue")
}
if (!dir.exists("processed_data/chikungunya")) {
  dir.create("processed_data/chikungunya")
}

write_csv(dengue_merged, "processed_data/dengue/dengue_merged.csv.gz")
write_csv(chikungunya_merged, "processed_data/chikungunya/chikungunya_merged.csv.gz")
