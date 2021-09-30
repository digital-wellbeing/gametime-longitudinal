#' Check if int1 and int2 intersect
#' if they intersect:
#'      - set int1 to union of int1 & int2
#'      - set int2 to NA
#' if not intersect:
#'      - return original intervals
#' 
#' @return vector of int1, int2
merge_interval <- function(int2, int1) {
  if(!is.na(intersect(int1, int2))) {
    int1 <- union(int1, int2)
    int2 <- as.interval(NA)
    out <- c(int1, int2)
  } else {
    out <- c(int1, int2)
  }
  
  out
}

# Merge all overlapping intervals in df
merge_intervals_all <- function(interval, .pb = NULL) {
  if(!is.null(.pb)) .pb$tick()$print()
  interval2 <- interval
  tot_rows <- length(interval)
  for(i in seq_along(interval2)) {
    int1 <- interval2[i]
    if(is.na(int1)) next
    if(i < tot_rows) {
      for(j in (i + 1):tot_rows) {
        #cat("i: ", i, ", j: ", j, "\n")
        int2 <- interval2[j]
        merged <- merge_interval(int2, int1)
        int1 <- merged[1]
        interval2[i] <- int1
        interval2[j] <- merged[2]
        # we break the loop at the first non-overlapping interval
        # data must be sorted by `session_start` for this to work
        if(!is.na(merged[2])) break
      }
    }
  }
  interval2
}