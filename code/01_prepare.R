source("code/lib/load.R")
cfg <- read_run_config()
prepared <- prepare_revision(cfg)
message("Prepared local cohort. Share only the share/ subfolder.")
