# handle_na.r — Step 3 of the data-prep pipeline.
#
# The only systematic missingness left after merge_dengue.r /
# merge_chikungunya.r is in the ocean-climate indices (ENSO/IOD/PDO), which
# are reported on the column's own irregular schedule and can have gaps once
# joined onto the weekly epiweek grid. This script linearly interpolates
# those three columns over time (`zoo::na.approx`, leaving any leading/
# trailing NA untouched via `na.rm = FALSE`) and overwrites the merged
# tables in place.
#
# Run after merge_dengue.r / merge_chikungunya.r; followed by
# data_prep/agg_data_uf.r.
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
