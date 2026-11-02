source("code/lib/load.R")
cfg <- read_run_config(Sys.getenv("CLIF_RUN_CONFIG",unset="config/config.json"))
prepared <- prepare_revision(cfg)
message("Prepared local cohort. Share only the share/ subfolder.")
