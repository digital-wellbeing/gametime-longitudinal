#' Summarize meta-analyses
#' @param object a grouped df with all the brmsfits; e.g., fit_ma
#'
#' @return a list with
#'  - `summary`:
#'  - `draws`: the posterior draws
summarize_ma <- function(object) {

  # Data for plotting meta-analyses' posterior distributions
  p_data <- object %>%
    mutate(out = map(fit, get_ma_post)) %>%
    select(-fit) %>%
    unnest(out) %>%
    rename(Average = Intercept) %>%
    pivot_longer(`AC:NH`:Average, names_to = "Game") %>%
    mutate(Game = factor(Game, levels = c("Average", rev(unique(d$Game))))) %>%
    ungroup() %>%
    left_join(Order)

  # Summaries of posterior distributions
  fit_ma_sum <- p_data %>%
    group_by(y_var, x_var, Parameter, Type, Game, Panel_label) %>%
    summarise(
      describe_posterior(
        value,
        test = c("pd"),
        centrality = "mean"
      ) %>% select(-Parameter)
    ) %>%
    ungroup() %>%
    mutate(
      across(
        c(Mean, CI_low, CI_high),
        .fns = list(r = ~ format(round(., 2), nsmall = 2))
      )
    ) %>%
    mutate(Res = str_glue("{Mean_r} [{CI_low_r}, {CI_high_r}]")) %>%
    # Ordered labels
    left_join(Order) %>%
    mutate(Game = factor(Game, levels = c("Average", rev(unique(d$Game)))))

  list(
    summary = fit_ma_sum,
    draws = p_data
  )
}

#' Function to draw forest plots
#' @param object object created using `summarize_ma`
#' @param type
#' @param x_var
#' @param parameters
#' @param lavaan
#'
#' @return a ggplot2 object
forest_plot <- function(object,
                        type = "Unstandardized",
                        x_var = "Hours",
                        parameters = c("wy2 ~ wx1", "wx2 ~ wy1"),
                        x_limits = c(-0.2, 0.2),
                        lavaan = FALSE) {
  out <- object$summary %>%
    filter(
      Type == {{ type }},
      x_var == {{ x_var }},
      Parameter %in% {{ parameters }}
    ) %>%
    # We can use the xintercept mapping to control alpha below
    ggplot(aes(Mean, Game)) +
    coord_cartesian(xlim = x_limits) +
    scale_x_continuous(
      "Estimated cross-lagged effect",
      breaks = pretty_breaks(7),
      expand = expansion(.01)
    ) +
    scale_y_discrete(
      expand = expansion(c(0.15, 0.25)),
      labels = function(x) ifelse(x == "Average", "**Average**", x)
    ) +
    scale_fill_brewer(
      palette = "Set1", direction = 1, aesthetics = c("fill", "color")
    ) +
    scale_alpha_manual(values = c(.4, .8)) +
    # Vertical line for 0
    geom_vline(xintercept = 0, size = .1, col = "grey60") +
    # Posterior densities
    stat_halfeye(
      data = object$draws %>%
        filter(
          Type == {{ type }},
          x_var == {{ x_var }},
          Parameter %in% {{ parameters }}
        ),
      aes(
        value,
        fill = stat(x > 0),
        alpha = Game == "Average"
      ),
      height = .75, adjust = 1.5, point_interval = NULL,
      show.legend = FALSE, normalize = "panels"
    ) +
    # Summary geoms and texts
    geom_pointrangeh(
      aes(x = Mean, xmin = CI_low, xmax = CI_high),
      size = .2, fatten = 1.25, position = position_nudge(y = -.01)
    ) +
    # CIs in right margin
    geom_text(
      vjust = -0.5, size = 2.6, hjust = 1,
      aes(x = x_limits[2], label = Res),
      family = Font
    ) +
    # Posterior probability of direction for average
    geom_text(
      data = object$summary %>%
        filter(Game == "Average") %>%
        filter(
          Type == {{ type }},
          x_var == {{ x_var }},
          Parameter %in% {{ parameters }}
        ),
      vjust = 1.4, size = 3,
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
        data = d_ma %>%
          filter(
            Type == {{ type }},
            x_var == {{ x_var }},
            Parameter %in% {{ parameters }}
          ) %>%
          left_join(Order),
        col = "gray40", fatten = 1.25, size = .33,
        aes(x = est, xmin = ci.lower, xmax = ci.upper),
        position = position_nudge(y = -.075)
      )
  } else {
    out
  }
}