 OUTPUT_DIR <- "Models/RI-CLPM/html"

# RICLPM example ----
rmarkdown::render(
     "Models/RI-CLPM/RICLPM-lavaan-example.Rmd",
     output_dir = OUTPUT_DIR,
     encoding = "UTF-8"
)

# DGP tests ----
rmarkdown::render(
     "Models/RI-CLPM/CLPM-lavaan-simulation.Rmd",
     output_dir = OUTPUT_DIR,
     encoding = "UTF-8"
)
rmarkdown::render(
     "Models/RI-CLPM/RICLPM-lavaan-simulation.Rmd",
     output_dir = OUTPUT_DIR,
     encoding = "UTF-8"
)

# RICLPM vs OLS adjusted ----
rmarkdown::render(
     "Models/RI-CLPM/RICLPM-vs-DAG.Rmd",
     output_dir = OUTPUT_DIR,
     encoding = "UTF-8"
)

## WIP: brms, multilevel version ----
rmarkdown::render(
     "Models/RI-CLPM/CLPM-brms-simulation.Rmd",
     output_dir = OUTPUT_DIR,
     encoding = "UTF-8"
)
rmarkdown::render(
     "Models/RI-CLPM/RICLPM-brms-simulation.Rmd",
     output_dir = OUTPUT_DIR,
     encoding = "UTF-8"
)

## time-varying confounding
rmarkdown::render(
     "Models/RI-CLPM/RICLPM-time-varying-confounding.Rmd",
     output_dir = OUTPUT_DIR,
     encoding = "UTF-8"
)
