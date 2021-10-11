ARG R_VERSION
FROM rocker/r-ver:${R_VERSION}

## Install dependencies
RUN apt-get update \
  && apt-get install -y --no-install-recommends \
    pandoc \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    libfontconfig1-dev \
    libxt6 \
    libcairo2-dev \
    libv8-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    libfreetype6-dev \
    libpng-dev \
    libtiff5-dev \
    libjpeg-dev

## RENV
## Restore packages and install cmdstan
ARG RENV_VERSION
RUN R -e "install.packages('remotes', repos = c(CRAN = 'https://cloud.r-project.org'))"
RUN R -e "remotes::install_github('rstudio/renv@${RENV_VERSION}')"
WORKDIR /home
COPY renv.lock renv.lock
RUN R -e 'renv::restore()'
RUN R -e 'cmdstanr::install_cmdstan()'