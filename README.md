# Longitudinal study of game play and well-being

The data and analytic code supporting the MS is here. The raw data from game publishers & qualtrics is stored elsewhere because of some potentially sensitive data (internal/external IDs, withdraw and no-consent responses.) Those were removed before data was placed in Data/.


# Docker
Build the Docker image
```
docker build \
    --no-cache \
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
    gametime-longitudinal \
    R -e 'renv::restore(prompt = FALSE); bookdown::render_book()'
```
