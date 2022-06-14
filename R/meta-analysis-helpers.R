#' Get posterior samples from a meta-analysis
#'
#' @param x a brmsfit
#'
#' @return a tibble with posterior samples of the game-specific and mean parameters
get_ma_post <- function(x) {
  coef(x, summary = FALSE) %>%
    .[["Game"]] %>%
    .[, , 1] %>%
    as.data.frame() %>%
    cbind(fixef(x, summary = FALSE)) %>%
    as_tibble() %>%
    rename(Average = Intercept)
}

#' Function to draw forest plots
#' @param data
#' @param lavaan
#'
#' @return a ggplot2 object
forest_plot <- function(
    data,
    x_limits = c(-0.2, 0.2),
    lavaan = FALSE
) {
  out <- unnest(data, summary) %>%
    # We can use the xintercept mapping to control alpha below
    ggplot(aes(Mean, Game)) +
    # coord_cartesian(xlim = x_limits) +
    scale_x_continuous(
      "Estimated cross-lagged effect",
      breaks = pretty_breaks(7),
      expand = expansion(.0)
    ) +
    scale_y_discrete(
      expand = expansion(c(0.08, 0.1)),
      labels = function(x) ifelse(x == "Average", "**Average**", x)
    ) +
    scale_fill_brewer(
      palette = "Set1", direction = 1, aesthetics = c("fill", "color")
    ) +
    geom_vline(xintercept = 0, size = .1, col = "grey60") +
    scale_alpha_manual(values = c(.4, .8)) +
    # Posterior densities
    stat_halfeye(
      data = unnest(data, posterior) %>%
        select(-data, -summary) %>%
        pivot_longer(
          c(unique(d$Game), "Average"),
          names_to = "Game"
        ) %>%
        mutate(Game = fct_relevel(Game, "Average")),
      aes(
        value,
        fill = stat(x > 0),
        alpha = Game == "Average"
      ),
      height = .75, adjust = 1, point_interval = NULL,
      show.legend = FALSE, normalize = "panels"
    ) +
    # Summary geoms and texts
    geom_pointrangeh(
      aes(x = Mean, xmin = CI_low, xmax = CI_high),
      size = .175, fatten = 1, position = position_nudge(y = -.01)
    ) +
    # CIs in right margin
    geom_text(
      vjust = -0.5, size = 2.3, hjust = 1.05,
      aes(
        x = Inf,
        label = str_glue("{number(Mean, .01)} [{number(CI_low, .01)}, {number(CI_high, .01)}]")
        ),
      family = Font
    ) +
    # Posterior probability of direction for average
    geom_text(
      data = . %>% filter(Game == "Average"),
      vjust = 1.4, size = 2.1,
      aes(
        x = sign(Mean) * 0,
        hjust = ifelse(sign(Mean) == 1, 0, 1),
        label = str_glue("{percent(pd, .1)}"),
        col = as.logical(sign(Mean) == 1)
      ),
      family = Font,
      show.legend = FALSE
    ) +
    theme(
      axis.title.y = element_blank(),
      axis.text.y = element_markdown(),
      panel.spacing.x = unit(10, "pt")
    ) +
    facet_wrap("Panel_label", scales = "free_x", labeller = label_parsed)

  # Also display lavaan game-specific estimates?
  if (lavaan) {
    out +
      geom_pointrangeh(
        data = unnest(data, data),
        fatten = 1.25, size = .33, shape = 2,
        aes(x = est, xmin = ci.lower, xmax = ci.upper),
        position = position_nudge(y = -.075)
      )
  } else {
    out
  }
}
