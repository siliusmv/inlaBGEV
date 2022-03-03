library(inlaBGEV)
library(ggplot2)
library(magrittr)
library(sf)

# In this script we display all the available observation locations in a map

coords = observations %>%
  st_transform(get_proj_xy()) %>%
  dplyr::distinct(id, .keep_all = TRUE) %>%
  {cbind(., st_coordinates(.))}

plot = plot_on_map(coords, response_name = "n_years") %>%
  style_map_plot(prediction_grid, use_tex = TRUE) +
  theme(text = element_text(size = 17),
        legend.title = element_text(size = 17),
        legend.text = element_text(size = 17)) +
  labs(col = "Years")

tikz_plot(file.path(here::here(), "results", "station-locations.pdf"),
          plot, width = 6, height = 9)


plot1 = coords %>%
  dplyr::group_by(n_years) %>%
  dplyr::summarise(n = n()) %>%
  sf::st_drop_geometry() %>%
  ggplot() +
  ggpattern::geom_col_pattern(aes(x = n_years, y = n),
                              fill = NA,
                              col = "black",
                              pattern_fill = "black",
                              pattern_density = .03,
                              pattern_angle = 135,
                              pattern_spacing = .01) +
  theme_light() +
  theme(text = element_text(size = 17),
        legend.title = element_text(size = 17),
        legend.text = element_text(size = 17),
        panel.grid = element_blank()) +
  labs(y = "Number of stations", x = "Years", title = "a)") +
  scale_x_continuous(breaks = 1:16, expand = c(.02, .02)) +
  scale_y_continuous(expand = expansion(c(0, .02)))

plot2 = coords %>%
  {.[order(.$n_years), ]} %>%
  plot_on_map(response_name = "n_years") %>%
  style_map_plot(prediction_grid, use_tex = TRUE) +
  theme(text = element_text(size = 17),
        legend.title = element_text(size = 17),
        legend.text = element_text(size = 17)) +
  labs(col = "Years", title = "b)")

pp = (plot1 + theme(plot.margin = unit(c(0, 100, 0, 0), "pt"))) + plot2

tikz_plot(file.path(here::here(), "results", "station-locations-with-histogram.pdf"),
          pp, width = 15 * .75, height = 10 * .75)
