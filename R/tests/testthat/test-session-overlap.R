library(here)
library(lubridate)
library(tidyverse)
source(here("R/merge_intervals.R"))

## dummy data
tmp <- data.frame(
    interval = c(
        interval(
            as_datetime("2011-12-31 12:00:00"),
            as_datetime("2011-12-31 14:00:00")
            ),
        interval(
            as_datetime("2011-12-31 13:00:00"),
            as_datetime("2011-12-31 14:00:00")
        ),
        interval(
            as_datetime("2011-12-31 13:30:00"),
            as_datetime("2011-12-31 16:00:00")
        ),
        interval(
            as_datetime("2011-12-31 16:30:00"),
            as_datetime("2011-12-31 17:00:00")
        ),
        interval(
            as_datetime("2011-12-31 11:00:00"),
            as_datetime("2011-12-31 12:30:00")
        ),
        interval(
            as_datetime("2011-12-31 11:00:00"),
            as_datetime("2011-12-31 12:30:00")
        )
    )
) %>%
mutate(
  session_start = int_start(interval),
  session_end = int_end(interval)
) %>%
arrange(session_start) %>%
mutate(
  interval_merged = merge_intervals_all(interval)
)

test_that("merge_intervals method give expected results", {
  hours <- sum(as.duration(tmp$interval), na.rm = TRUE) / 3600
  hours_merged <- sum(as.duration(tmp$interval_merged), na.rm = TRUE) / 3600
  expect_equal(hours, 9)
  expect_equal(hours_merged, 5.5)
  expect_lt(hours_merged, hours)
})

test_that("merged and summed data agree with manual calculation", {
    # a participant with overlapping sessions and negative sessions
    # "raw" data
    tmp <- t_acnh %>%
    filter(pid == "05588dba956d2780") %>%
    arrange(session_start) %>%
    mutate(
        interval = interval(session_start, session_end)
    ) %>%
    filter(as.duration(interval) > 0) %>%
    mutate(
        interval_merged = merge_intervals_all(interval)
    )

    # same participant in cleaned data
    tmp_cleaned <- d %>%
        filter(pid == "05588dba956d2780") %>%
        select(StartDate, Hours)

    # only last session belong to wave 1
    # should be approx 1.12 hours
    hours <-  as.numeric(as.duration(tmp$interval_merged[9])) / 3600

    expect_equal(c(hours, 0, 0), tmp_cleaned$Hours)
})

test_that("duration of merged intervals is shorter or equal", {
    data_path <- here("Temp", "session-overlap-merged.rds")
    d_t <- read_rds(file = data_path)
    tmp <- d_t %>%
        group_by(pid) %>%
        summarise(
            duration = sum(as.duration(interval), na.rm = TRUE),
            duration_merged = sum(as.duration(interval_merged), na.rm = TRUE)
        )
    # all sum(duration_merged) should be less than or equal to sum(duration)
    # as all overlap is removed
    expect_true(all(tmp$duration_merged <= tmp$duration))
})
