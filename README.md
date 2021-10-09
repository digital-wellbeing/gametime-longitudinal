# Longitudinal study of game play and well-being

This repository ([GitHub](https://osf.io/fb38n/) / [OSF](https://osf.io/fb38n/)) contains the data and code required to reproduce all analyses reported in our manuscript, *Time spent playing video games is unlikely to impact well-being* (Vuorre, Johannes, Magnusson, & Przybylski, 2021). These analyses are presented in the [Online analysis supplement](https://digital-wellbeing.github.io/gametime-longitudinal).

## Materials

- [Preprint](todo)  
  - A publicly available version of our manuscript in advance of peer-review and formal publication
- [GitHub repository](https://github.com/digital-wellbeing/gametime-longitudinal)  
  - A version controlled repository containing all the raw data and code in this project
- [OSF repository](https://osf.io/fb38n/)  
  - An archived permanent copy of the GitHub repository
- [Online analysis supplement](https://digital-wellbeing.github.io/gametime-longitudinal)
  - The output document of our analyses

In addition to raw data, the repository contains many survey and telemetry variables that we did not analyse, but that may be of interest to further analyses of gameplay and well being. 

## Reproducibility

The raw data are in the `Data/` directory of this repository. The code that we used to clean and analyse the data are organised in R Markdown (`.Rmd`) files in this directory, which are meant to be run in the sequence indicated by their numeric prefixes. To run all the cleaning and analyses, and compile the resulting document ([the online analysis supplement](https://digital-wellbeing.github.io/gametime-longitudinal)), run `bookdown::render_book()` in R or click "Build Book" in the RStudio IDE.

### Docker

To ensure the reproducibility of our analyses, you can use Docker:

Build the Docker image
```
docker build \
    --build-arg R_VERSION=4.1.1 \
    --build-arg RENV_VERSION=0.14.0 \
    -t gametime-longitudinal .
```

Run the container and render output
```
docker run \
    --rm \
    -v "$(pwd):/home/" \
    -v "/home/renv/library" \
    -e MAX_CORES=1 \
    gametime-longitudinal \
    R -e 'renv::restore(prompt = FALSE); bookdown::render_book()'
```
