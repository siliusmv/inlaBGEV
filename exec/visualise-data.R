library(inlaBGEV)
library(ggplot2)

plot = plot_on_map(observations) %>%
  style_map_plot(observations, use_tex = TRUE) +
  theme(text = element_text(size = 20))

tikz_plot(file.path(here::here(), "results", "station-locations.pdf"),
          plot, width = 6, height = 10, view = TRUE)
