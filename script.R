#!/usr/bin/env Rscript

source("threemc_fit.R")

# load shell dataset
shell_dat <- readr::read_csv("data/shell_data.csv.gz")
# load shapefiles
areas <- sf::read_sf("data/areas.geojson")

# read in model
mod <- readLines("src/threemc.cpp")

fit <- threemc_fit(shell_dat, areas, mod)
