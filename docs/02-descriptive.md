# Descriptives

Load the required packages.


```r
library(knitr)
library(scales)
library(gtsummary)
library(kableExtra)
library(here)
library(showtext)
library(tidyverse)
library(lubridate)
```

Figure options.


```r
# Plotting options
Font <- "Titillium Web"
font_add_google(Font, Font)
theme_set(
  theme_linedraw(
    base_family = Font,
    base_size = 12
  ) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank()
    )
)
```

We then load the previously cleaned data table.


```r
data_path <- here("Data", "cleaned_data.rds")
if (file.exists(data_path)) {
  d <- read_rds(file = data_path)
} else {
  stop(str_glue("{data_path} doesn't exist, run `01-clean.Rmd` to create it."))
}

# Make wave a nicely labelled factor
d <- d %>%
  mutate(Wave = factor(wid, levels = 1:3, labels = paste0("Wave ", 1:3)))

# Rename game to fit titles in plots
d <- d %>% 
  mutate(Game = if_else(Game == "Gran Turismo Sport", "GT Sport", Game))
```

## Demographics after exclusions


```r
d %>%
  filter(wid == 1) %>%
  distinct(pid, Game, Age, Experience, Gender) %>%
  select(-pid) %>%
  tbl_summary(by = Game, missing_text = "Missing") %>%
  add_overall() %>%
  bold_labels() %>%
  italicize_levels() %>%
  as_kable_extra(caption = "Sample demographics") %>% 
  kable_styling(full_width = FALSE, font_size = 12)
```

<table style="NAborder-bottom: 0; font-size: 12px; width: auto !important; margin-left: auto; margin-right: auto;" class="table">
<caption style="font-size: initial !important;">(\#tab:demographics-table)Sample demographics</caption>
 <thead>
  <tr>
   <th style="text-align:left;"> Characteristic </th>
   <th style="text-align:left;"> Overall, N = 38,935 </th>
   <th style="text-align:left;"> AC:NH, N = 13,646 </th>
   <th style="text-align:left;"> Apex Legends, N = 1,158 </th>
   <th style="text-align:left;"> EVE Online, N = 905 </th>
   <th style="text-align:left;"> Forza Horizon 4, N = 1,981 </th>
   <th style="text-align:left;"> GT Sport, N = 19,258 </th>
   <th style="text-align:left;"> Outriders, N = 1,530 </th>
   <th style="text-align:left;"> The Crew 2, N = 457 </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> Age </td>
   <td style="text-align:left;"> 34 (25, 42) </td>
   <td style="text-align:left;"> 32 (25, 41) </td>
   <td style="text-align:left;"> 25 (20, 32) </td>
   <td style="text-align:left;"> 39 (31, 50) </td>
   <td style="text-align:left;"> 33 (24, 42) </td>
   <td style="text-align:left;"> 35 (25, 43) </td>
   <td style="text-align:left;"> 38 (32, 45) </td>
   <td style="text-align:left;"> 25 (20, 35) </td>
  </tr>
  <tr>
   <td style="text-align:left;padding-left: 2em;" indentlevel="1"> Missing </td>
   <td style="text-align:left;"> 82 </td>
   <td style="text-align:left;"> 29 </td>
   <td style="text-align:left;"> 5 </td>
   <td style="text-align:left;"> 3 </td>
   <td style="text-align:left;"> 6 </td>
   <td style="text-align:left;"> 37 </td>
   <td style="text-align:left;"> 2 </td>
   <td style="text-align:left;"> 0 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Gender </td>
   <td style="text-align:left;">  </td>
   <td style="text-align:left;">  </td>
   <td style="text-align:left;">  </td>
   <td style="text-align:left;">  </td>
   <td style="text-align:left;">  </td>
   <td style="text-align:left;">  </td>
   <td style="text-align:left;">  </td>
   <td style="text-align:left;">  </td>
  </tr>
  <tr>
   <td style="text-align:left;padding-left: 2em;" indentlevel="1"> Man </td>
   <td style="text-align:left;"> 29,765 (77%) </td>
   <td style="text-align:left;"> 5,451 (40%) </td>
   <td style="text-align:left;"> 1,002 (87%) </td>
   <td style="text-align:left;"> 853 (95%) </td>
   <td style="text-align:left;"> 1,884 (95%) </td>
   <td style="text-align:left;"> 18,745 (98%) </td>
   <td style="text-align:left;"> 1,400 (92%) </td>
   <td style="text-align:left;"> 430 (94%) </td>
  </tr>
  <tr>
   <td style="text-align:left;padding-left: 2em;" indentlevel="1"> Non-binary / third gender </td>
   <td style="text-align:left;"> 705 (1.8%) </td>
   <td style="text-align:left;"> 557 (4.1%) </td>
   <td style="text-align:left;"> 25 (2.2%) </td>
   <td style="text-align:left;"> 7 (0.8%) </td>
   <td style="text-align:left;"> 11 (0.6%) </td>
   <td style="text-align:left;"> 87 (0.5%) </td>
   <td style="text-align:left;"> 15 (1.0%) </td>
   <td style="text-align:left;"> 3 (0.7%) </td>
  </tr>
  <tr>
   <td style="text-align:left;padding-left: 2em;" indentlevel="1"> Prefer not to say </td>
   <td style="text-align:left;"> 363 (0.9%) </td>
   <td style="text-align:left;"> 183 (1.3%) </td>
   <td style="text-align:left;"> 17 (1.5%) </td>
   <td style="text-align:left;"> 10 (1.1%) </td>
   <td style="text-align:left;"> 12 (0.6%) </td>
   <td style="text-align:left;"> 119 (0.6%) </td>
   <td style="text-align:left;"> 16 (1.0%) </td>
   <td style="text-align:left;"> 6 (1.3%) </td>
  </tr>
  <tr>
   <td style="text-align:left;padding-left: 2em;" indentlevel="1"> Woman </td>
   <td style="text-align:left;"> 8,017 (21%) </td>
   <td style="text-align:left;"> 7,426 (55%) </td>
   <td style="text-align:left;"> 109 (9.5%) </td>
   <td style="text-align:left;"> 32 (3.5%) </td>
   <td style="text-align:left;"> 68 (3.4%) </td>
   <td style="text-align:left;"> 267 (1.4%) </td>
   <td style="text-align:left;"> 97 (6.3%) </td>
   <td style="text-align:left;"> 18 (3.9%) </td>
  </tr>
  <tr>
   <td style="text-align:left;padding-left: 2em;" indentlevel="1"> Missing </td>
   <td style="text-align:left;"> 85 </td>
   <td style="text-align:left;"> 29 </td>
   <td style="text-align:left;"> 5 </td>
   <td style="text-align:left;"> 3 </td>
   <td style="text-align:left;"> 6 </td>
   <td style="text-align:left;"> 40 </td>
   <td style="text-align:left;"> 2 </td>
   <td style="text-align:left;"> 0 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Experience </td>
   <td style="text-align:left;"> 23 (16, 30) </td>
   <td style="text-align:left;"> 22 (15, 30) </td>
   <td style="text-align:left;"> 17 (10, 25) </td>
   <td style="text-align:left;"> 25 (20, 32) </td>
   <td style="text-align:left;"> 23 (15, 30) </td>
   <td style="text-align:left;"> 25 (16, 30) </td>
   <td style="text-align:left;"> 30 (22, 35) </td>
   <td style="text-align:left;"> 17 (11, 25) </td>
  </tr>
  <tr>
   <td style="text-align:left;padding-left: 2em;" indentlevel="1"> Missing </td>
   <td style="text-align:left;"> 181 </td>
   <td style="text-align:left;"> 62 </td>
   <td style="text-align:left;"> 8 </td>
   <td style="text-align:left;"> 5 </td>
   <td style="text-align:left;"> 9 </td>
   <td style="text-align:left;"> 91 </td>
   <td style="text-align:left;"> 5 </td>
   <td style="text-align:left;"> 1 </td>
  </tr>
</tbody>
<tfoot><tr><td style="padding: 0; " colspan="100%">
<sup>1</sup> Median (IQR); n (%)</td></tr></tfoot>
</table>

## Survey descriptives

These are for people after exclusions.

### Response rate & retention


```r
# Get data on invite dates and Ns
invites <- read_csv(here("Data", "invites.csv")) %>%
  rename(Game = game) %>%
  mutate(
    Game = str_replace(Game, "Animal Crossing: New Horizons", "AC:NH"),
    Game = str_replace(Game, "Gran Turismo Sport", "GT Sport")
    )
# Create a table where wave 0 are number of invites,
# then calculate response rate / retention at each wave.
# This assumes there are no new participants at wave 3
# (people who didn't participate in wave 2 showed up at wave 3).
bind_rows(
  select(invites, -date),
  d %>% filter(Responded) %>% count(Game, wid)
) %>%
  arrange(Game, wid) %>%
  group_by(Game) %>%
  mutate(
    R_rate = percent(n / lag(n), .01),
    n = comma(n)
  ) %>%
  pivot_wider(names_from = wid, values_from = c(n, R_rate)) %>%
  mutate(
    Invites = n_0,
    `Wave 1` = str_glue("{n_1} ({R_rate_1})"),
    `Wave 2` = str_glue("{n_2} ({R_rate_2})"),
    `Wave 3` = str_glue("{n_3} ({R_rate_3})")
  ) %>%
  select(Game, Invites:`Wave 3`) %>%
  mutate(across(everything(), ~ str_replace(., "NA", "0"))) %>%
  kbl(caption = "Number of people (response/retention rate) participating at each wave.") %>% 
  kable_styling(full_width = FALSE, font_size = 12)
```

<table class="table" style="font-size: 12px; width: auto !important; margin-left: auto; margin-right: auto;">
<caption style="font-size: initial !important;">(\#tab:response-rate-table)Number of people (response/retention rate) participating at each wave.</caption>
 <thead>
  <tr>
   <th style="text-align:left;"> Game </th>
   <th style="text-align:left;"> Invites </th>
   <th style="text-align:left;"> Wave 1 </th>
   <th style="text-align:left;"> Wave 2 </th>
   <th style="text-align:left;"> Wave 3 </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> AC:NH </td>
   <td style="text-align:left;"> 640,000 </td>
   <td style="text-align:left;"> 13,536 (2.11%) </td>
   <td style="text-align:left;"> 5,049 (37.30%) </td>
   <td style="text-align:left;"> 4,084 (80.89%) </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Apex Legends </td>
   <td style="text-align:left;"> 900,000 </td>
   <td style="text-align:left;"> 1,128 (0.13%) </td>
   <td style="text-align:left;"> 406 (35.99%) </td>
   <td style="text-align:left;"> 228 (56.16%) </td>
  </tr>
  <tr>
   <td style="text-align:left;"> EVE Online </td>
   <td style="text-align:left;"> 30,000 </td>
   <td style="text-align:left;"> 899 (3.00%) </td>
   <td style="text-align:left;"> 240 (26.70%) </td>
   <td style="text-align:left;"> 221 (92.08%) </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Forza Horizon 4 </td>
   <td style="text-align:left;"> 834,515 </td>
   <td style="text-align:left;"> 1,959 (0.23%) </td>
   <td style="text-align:left;"> 772 (39.41%) </td>
   <td style="text-align:left;"> 597 (77.33%) </td>
  </tr>
  <tr>
   <td style="text-align:left;"> GT Sport </td>
   <td style="text-align:left;"> 1,729,677 </td>
   <td style="text-align:left;"> 19,073 (1.10%) </td>
   <td style="text-align:left;"> 7,699 (40.37%) </td>
   <td style="text-align:left;"> 5,512 (71.59%) </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Outriders </td>
   <td style="text-align:left;"> 90,000.0 </td>
   <td style="text-align:left;"> 1,525.0 (1.69%) </td>
   <td style="text-align:left;"> 379.0 (24.85%) </td>
   <td style="text-align:left;"> 370.0 (97.63%) </td>
  </tr>
  <tr>
   <td style="text-align:left;"> The Crew 2 </td>
   <td style="text-align:left;"> 1,013,000 </td>
   <td style="text-align:left;"> 457 (0.05%) </td>
   <td style="text-align:left;"> 97 (21.23%) </td>
   <td style="text-align:left;"> 85 (87.63%) </td>
  </tr>
</tbody>
</table>

### Response dates

Only for actual responses (not rows where survey date was filled to be able to join telemetry)


```r
d %>%
  # Take only actually responded-to waves
  filter(Responded) %>%
  mutate(Date = as_date(StartDate)) %>%
  count(Game, Wave, Date) %>%
  ggplot(
    aes(Date, n, fill = Wave)
  ) +
  geom_col() +
  scale_y_continuous(
    "Responses",
    breaks = pretty_breaks(),
    expand = expansion(c(0, .1)),
  ) +
  scale_x_date(
    "Date",
    date_breaks = "7 day", date_labels = "%b\n%d", date_minor_breaks = "1 day"
  ) +
  facet_wrap("Game", scales = "free_y", ncol = 1)
```

<div class="figure" style="text-align: center">
<img src="02-descriptive_files/figure-html/response-date-histograms-1.png" alt="Histograms of response dates." width="672" />
<p class="caption">(\#fig:response-date-histograms)Histograms of response dates.</p>
</div>

Response times


```r
d %>%
  filter(Responded) %>%
  mutate(Hour = hour(StartDate)) %>%
  count(Game, Wave, Hour) %>%
  ggplot(aes(Hour, y = n, fill = Wave)) +
  scale_y_continuous(
    "Responses",
    breaks = pretty_breaks(),
    expand = expansion(c(0, .1)),
  ) +
  scale_x_continuous(
    breaks = seq(0, 21, by = 3),
    expand = expansion(c(0.01))
  ) +
  geom_col() +
  facet_wrap("Game", scales = "free", ncol = 2) +
  theme(legend.position = "bottom")
```

<div class="figure" style="text-align: center">
<img src="02-descriptive_files/figure-html/response-time-histograms-1.png" alt="Histograms of response times (in UTC)." width="672" />
<p class="caption">(\#fig:response-time-histograms)Histograms of response times (in UTC).</p>
</div>

#### Durations between waves

Participants could respond with variable delays due to variation in email schedules and late responding. So we also check the actual intervals between completing waves. Very small values are possible because a participant could have e.g. completed both waves 2 and 3 in succession after receiving wave 3 invitation. Note that negative values were also possible for this reason but they were excluded before. (This figure is restricted to 5-30 day intervals to display the bulk of the data.)


```r
# Table
d %>%
  select(Wave, Game, interval) %>%
  # group_by(Wave) %>%
  filter(Wave != "Wave 1") %>%
  summarise(
    Value = quantile(
      interval,
      probs = c(0, .10, .25, .5, .75, .90, 1),
      na.rm = T
    ) %>%
      round(3),
    Quantile = percent(c(0, .10, .25, .5, .75, .90, 1))
  ) %>%
  pivot_wider(names_from = Quantile, values_from = Value) %>%
  kbl(caption = "Interval duration percentiles preceding waves 2 and 3.") %>% 
  kable_styling(full_width = FALSE, font_size = 12)
```

<table class="table" style="font-size: 12px; width: auto !important; margin-left: auto; margin-right: auto;">
<caption style="font-size: initial !important;">(\#tab:intervals-between-waves-histogram)Interval duration percentiles preceding waves 2 and 3.</caption>
 <thead>
  <tr>
   <th style="text-align:right;"> 0% </th>
   <th style="text-align:right;"> 10% </th>
   <th style="text-align:right;"> 25% </th>
   <th style="text-align:right;"> 50% </th>
   <th style="text-align:right;"> 75% </th>
   <th style="text-align:right;"> 90% </th>
   <th style="text-align:right;"> 100% </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:right;"> 0.001 </td>
   <td style="text-align:right;"> 11.649 </td>
   <td style="text-align:right;"> 13.706 </td>
   <td style="text-align:right;"> 14.046 </td>
   <td style="text-align:right;"> 15.941 </td>
   <td style="text-align:right;"> 17.105 </td>
   <td style="text-align:right;"> 49.287 </td>
  </tr>
</tbody>
</table>

```r
# Figure
d %>%
  filter(Wave != "Wave 1") %>%
  filter(between(interval, 5, 30)) %>% 
  mutate(Wave = fct_drop(Wave)) %>%
  ggplot(aes(interval)) +
  geom_vline(xintercept = 14, size = .2) +
  geom_histogram(binwidth = 1, col = "white") +
  scale_y_continuous(
    "Count",
    expand = expansion(c(0, .1))
  ) +
  scale_x_continuous(
    "Days between responding",
    breaks = pretty_breaks()
  ) +
  facet_grid(Game ~ Wave, scales = "free_y")
```

<div class="figure" style="text-align: center">
<img src="02-descriptive_files/figure-html/intervals-between-waves-histogram-1.png" alt="Histograms of intervals between participants completing the survey waves (in days)." width="672" />
<p class="caption">(\#fig:intervals-between-waves-histogram)Histograms of intervals between participants completing the survey waves (in days).</p>
</div>

## System information


```r
sessionInfo()
```

```
## R version 4.1.1 (2021-08-10)
## Platform: x86_64-apple-darwin17.0 (64-bit)
## Running under: macOS Big Sur 10.16
## 
## Matrix products: default
## BLAS:   /Library/Frameworks/R.framework/Versions/4.1/Resources/lib/libRblas.0.dylib
## LAPACK: /Library/Frameworks/R.framework/Versions/4.1/Resources/lib/libRlapack.dylib
## 
## locale:
## [1] en_GB.UTF-8/en_GB.UTF-8/en_GB.UTF-8/C/en_GB.UTF-8/en_GB.UTF-8
## 
## attached base packages:
## [1] stats     graphics  grDevices datasets  utils     methods   base     
## 
## other attached packages:
##  [1] lubridate_1.8.0  forcats_0.5.1    stringr_1.4.0    dplyr_1.0.7     
##  [5] purrr_0.3.4      readr_2.0.2      tidyr_1.1.4      tibble_3.1.5    
##  [9] ggplot2_3.3.5    tidyverse_1.3.1  showtext_0.9-4   showtextdb_3.0  
## [13] sysfonts_0.8.5   here_1.0.1       kableExtra_1.3.4 gtsummary_1.4.2 
## [17] scales_1.1.1     knitr_1.36      
## 
## loaded via a namespace (and not attached):
##  [1] fs_1.5.0            bit64_4.0.5         webshot_0.5.2      
##  [4] httr_1.4.2          rprojroot_2.0.2     tools_4.1.1        
##  [7] backports_1.2.1     bslib_0.3.1         utf8_1.2.2         
## [10] R6_2.5.1            DBI_1.1.1           colorspace_2.0-2   
## [13] withr_2.4.2         tidyselect_1.1.1    downlit_0.2.1      
## [16] bit_4.0.4           curl_4.3.2          compiler_4.1.1     
## [19] textshaping_0.3.5   cli_3.0.1           rvest_1.0.1        
## [22] gt_0.3.1            xml2_1.3.2          labeling_0.4.2     
## [25] bookdown_0.24       sass_0.4.0          systemfonts_1.0.2  
## [28] digest_0.6.28       rmarkdown_2.11      svglite_2.0.0      
## [31] pkgconfig_2.0.3     htmltools_0.5.2     dbplyr_2.1.1       
## [34] fastmap_1.1.0       highr_0.9           rlang_0.4.11       
## [37] readxl_1.3.1        rstudioapi_0.13     farver_2.1.0       
## [40] jquerylib_0.1.4     generics_0.1.0      jsonlite_1.7.2     
## [43] vroom_1.5.5         magrittr_2.0.1      Matrix_1.3-4       
## [46] Rcpp_1.0.7          munsell_0.5.0       fansi_0.5.0        
## [49] lifecycle_1.0.1     stringi_1.7.5       yaml_2.2.1         
## [52] grid_4.1.1          parallel_4.1.1      crayon_1.4.1       
## [55] lattice_0.20-45     haven_2.4.3         splines_4.1.1      
## [58] hms_1.1.1           pillar_1.6.3        codetools_0.2-18   
## [61] reprex_2.0.1        glue_1.4.2          evaluate_0.14      
## [64] broom.helpers_1.4.0 renv_0.14.0         modelr_0.1.8       
## [67] vctrs_0.3.8         tzdb_0.1.2          cellranger_1.1.0   
## [70] gtable_0.3.0        assertthat_0.2.1    xfun_0.26          
## [73] broom_0.7.9         ragg_1.1.3          survival_3.2-13    
## [76] viridisLite_0.4.0   ellipsis_0.3.2
```
