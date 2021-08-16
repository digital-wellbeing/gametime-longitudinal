

fixed_or_random <- function(b_mean, b_sd, n) {
    if(is.null(b_sd)) {
        return(b_mean)
    } else {
        return(rnorm(n, b_mean, b_sd))
    }
}

#' Simulate 3-wave RI-CLPM data
#' 
#' @param n number of participants
#'
#' @return a data.frame with the observed variables:
#' x1, x2, x3, y1, y2, y3
sim_RICLPM_data <- function(
    n,
    beta_2x, # beta = autoregressive params
    beta_2y,
    beta_3x,
    beta_3y,
    gamma_2x, # gamma = (averag) cross-lagged effects
    gamma_2y,
    gamma_3x,
    gamma_3y,
    gamma_2x_sd = NULL, # cross-lagged Random Effects
    gamma_2y_sd = NULL,
    gamma_3x_sd = NULL,
    gamma_3y_sd = NULL,
    x_mean,
    y_mean,
    sd_wx1, # covariances
    sd_wy1,
    cor_wx1_wy1,
    sd_wx2,
    sd_wy2,
    cor_wx2_wy2,
    sd_wx3,
    sd_wy3,
    cor_wx3_wy3,
    sd_x_RE, # random intercepts
    sd_y_RE,
    cor_xy_intercepts,
    keep_cross_effects = FALSE
) {
    # covariances
    # Time = 1
    sigma1 <- matrix(
        c(
            sd_wx1^2, sd_wx1 * sd_wy1 * cor_wx1_wy1,
            sd_wx1 * sd_wy1 * cor_wx1_wy1, sd_wy1^2
        ),
        nrow = 2, ncol = 2
    )
    # Time = 2
    sigma2 <- matrix(
        c(
            sd_wx2^2, sd_wx2 * sd_wy2 * cor_wx2_wy2,
            sd_wx2 * sd_wy2 * cor_wx2_wy2, sd_wy2^2
        ),
        nrow = 2, ncol = 2
    )
    # Time = 3
    sigma3 <- matrix(
        c(
            sd_wx3^2, sd_wx3 * sd_wy3 * cor_wx3_wy3,
            sd_wx3 * sd_wy3 * cor_wx3_wy3, sd_wy3^2
        ),
        nrow = 2, ncol = 2
    )
    ## Random intercepts
    sigmaRE <- matrix(
        c(
            sd_x_RE^2, sd_x_RE * sd_y_RE * cor_xy_intercepts,
            sd_x_RE * sd_y_RE * cor_xy_intercepts, sd_y_RE^2
        ),
        nrow = 2, ncol = 2
    )

    # Sample random variables
    e1 <- MASS::mvrnorm(n, c(0, 0), sigma1)
    e2 <- MASS::mvrnorm(n, c(0, 0), sigma2)
    e3 <- MASS::mvrnorm(n, c(0, 0), sigma3)
    RE <- MASS::mvrnorm(n, c(0, 0), sigmaRE)
    U_x <- RE[, 1]
    U_y <- RE[, 2]

    # within-subjects cross-lagged regressions
    wx1 <- e1[, 1]
    wy1 <- e1[, 2]

    gamma_2x <- fixed_or_random(gamma_2x, gamma_2x_sd, n)
    gamma_2y <- fixed_or_random(gamma_2y, gamma_2y_sd, n)
    gamma_3x <- fixed_or_random(gamma_3x, gamma_3x_sd, n)
    gamma_3y <- fixed_or_random(gamma_3y, gamma_3y_sd, n)
    ## Assume constant RE effect over time for now
    if(!is.null(gamma_3y_sd)) gamma_3y <- gamma_2y
    wx2 <- beta_2x * wx1 + gamma_2x * wy1 + e2[, 1]
    wy2 <- beta_2y * wy1 + gamma_2y * wx1 + e2[, 2]
    wx3 <- beta_3x * wx2 + gamma_3x * wy2 + e3[, 1]
    wy3 <- beta_3y * wy2 + gamma_3y * wx2 + e3[, 2]

    # observed variables
    x1 <- wx1 + x_mean + U_x
    x2 <- wx2 + x_mean + U_x
    x3 <- wx3 + x_mean + U_x
    y1 <- wy1 + y_mean + U_y
    y2 <- wy2 + y_mean + U_y
    y3 <- wy3 + y_mean + U_y

    out <- data.frame(x1, x2, x3, y1, y2, y3, U_x, U_y)
    if(keep_cross_effects == TRUE) {
        out <- cbind(
            out,
            data.frame(gamma_2x, gamma_2y, gamma_3x, gamma_3y)
            )
    }
    out
}