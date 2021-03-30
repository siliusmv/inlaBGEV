library(inlaBGEV)
library(ggplot2)
library(sf)

plot = plot_on_map(observations, use_tex = TRUE)

tikz_plot(file.path(here::here(), "results", "station-locations.pdf"),
          plot, width = 10, height = 7, view = TRUE)
