# Set global options here
library(ggplot2)

# Global knitr options
knitr::opts_chunk$set(
  echo = TRUE,
  warning = FALSE,
  cache = TRUE,
  error = FALSE,
  message = FALSE,
  fig.align = "center",
  fig.width = 8,
  dev = "ragg_png",
  fig.retina = 2
)

# If you've defined a custom font, we will use it in outputs
Font <- Sys.getenv("FONT")

# Plotting options
theme_set(
  theme_linedraw(
    base_family = Font,
    base_size = 9
  ) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_blank(),
      axis.text = element_text(size = rel(0.75))
    )
)
col1 <- "#2980b9"
W <- 17.4  # Figure width in cm
