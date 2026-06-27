library(reticulate)
library(data.table)

py_config()

py_require(c("epiweeks","python-dotenv", "mosqlient"))

api_key <- "EzequielEBS:c74e2486-70d0-454a-983b-a12d55376324"
mosq <- import("mosqlient")

repository <- 'EzequielEBS/3rd_imdc_emap_epidematicos_sarimax_state' # fill with your repository name 
commit <- '6e1955aebd2c5874b2c489860186b458ecbbf217'
disease <- 'A90' # dengue prediction
adm_level <- 1 # state level prediction
adm_2 <- NULL
case_definition <- 'probable' # The IMDC uses probable cases 
published <- T

states <- c("AC", "AL", "AM", "AP", "BA", "CE", "DF", "GO", "MA",
            "MG", "MS", "MT", "PA", "PB", "PE", "PI", "PR", "RJ",
            "RN", "RO", "RR", "RS", "SC", "SE", "SP", "TO")

# submit dengue predictions for each state and target
for (st in states) {
  file_name <- paste0("processed_data/dengue/dengue_", st, "_agg.csv.gz")
  d <- read_csv(file_name, show_col_types = FALSE)
  adm_1 <- d$uf_code[1]
  target_ids <- paste0("target_", 1:4)

  for (target_id in target_ids) {
    pred_file <- paste0("sarimax/results/preds/pred_dengue_", st, "_", target_id, ".csv")
    pred <- read_csv(pred_file, show_col_types = FALSE)
    description <- paste0("Dengue prediction for state ", st, " and target ", target_id)
    mosq$upload_prediction(
                api_key = api_key,
                repository = repository,
                description = description,
                commit = commit,
                disease = disease,
                case_definition = case_definition,
                adm_level = adm_level,
                adm_1 = adm_1,
                published = published,
                prediction = pred
    )
  }
}