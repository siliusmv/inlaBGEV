library(inlaBGEV)
library(ggplot2)
library(magrittr)
library(sf)

# In this script we display all the available observation locations in a map

plot = plot_on_map(observations, response_name = "n_years") %>%
  style_map_plot(observations, use_tex = TRUE) +
  theme(text = element_text(size = 20)) +
  labs(col = "Years")

tikz_plot(file.path(here::here(), "results", "station-locations.pdf"),
          plot, width = 6, height = 9)
