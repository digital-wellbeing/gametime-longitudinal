# Set global options here
library(showtext)
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

# Plotting options
Font <- "Titillium Web"
font_add_google(Font, Font)
theme_set(
  theme_linedraw(
    base_family = Font,
    base_size = 9
  ) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank()
    )
)
col1 <- "#2980b9"
W <- 20  # Figure width in cm
