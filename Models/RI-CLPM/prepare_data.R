## Load and munge example data
library(tidyverse)
d <- read.table(
    "Models/RI-CLPM/data/RICLPM.dat",
    col.names = c("x1", "x2", "x3", "x4", "x5", "y1", "y2", "y3", "y4", "y5")
    )
d$id <- 1:nrow(d)

## Save subset (n = 100), for testing
set.seed(1337)
d_subset <- d[sample(1:nrow(d), 100), ]
saveRDS(d_subset, "Models/RI-CLPM/data/RICLPM-wide-subset.Rdata")


## Long format (subset, n = 100)
d_subset_long <- d_subset %>%
    pivot_longer(
        -id,
        names_to = c(".value", "time"),
        names_pattern = c("(.)(.)")
    )
saveRDS(d_subset_long, "Models/RI-CLPM/data/RICLPM-long-subset.Rdata")
