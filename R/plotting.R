#' @export
plot_stats = function(x, grid = TRUE,
                      upper = .975, lower = .025,
                      breaks = NULL, CI_breaks = NULL,
                      use_tex = FALSE,
                      ...) {
  x$CI = x[[paste0(as.character(upper * 100), "%")]] - x[[paste0(as.character(lower * 100), "%")]]
  if (grid) {
    gg1 = plot_grid(x, "mean", breaks, use_tex = use_tex, ...)
    gg2 = plot_grid(x, "CI", CI_breaks, use_tex = use_tex, ...)
  } else {
    gg1 = plot_on_map(x, "mean", breaks, use_tex = use_tex, ...)
    gg2 = plot_on_map(x, "CI", CI_breaks, use_tex = use_tex, ...)
  }
  gg1 = gg1 + labs(fill = "Posterior mean")
  if (use_tex) {
    gg2 = gg2 +
      labs(fill = paste0("Posterior $", as.character((upper - lower) * 100), "\\%$\nCI width"))
  } else {
    gg2 = gg2 +
      labs(fill = paste0("Posterior ", as.character((upper - lower) * 100), "%\nCI width"))
  }
  gg1 + gg2
}

#' @export
plot_grid = function(data, response_name, breaks = NULL, use_tex = FALSE, ...) {
  coords = as.data.frame(sf::st_coordinates(data))
  data$x = coords$X
  data$y = coords$Y
  if (!is.null(breaks)) {
    data[[response_name]] = get_binned_data(data[[response_name]], breaks, use_tex = use_tex)
  }
  gg = ggplot(data) +
    theme_light() +
    geom_raster(aes_string(x = "x", y = "y", fill = response_name)) +
    add_norway_map(sf::st_crs(data), sf::st_bbox(data), ...) +
    labs(x = "", y = "")
  gg = add_scico_fill(gg, binned = !is.null(breaks))
  if (use_tex) gg = latex_friendly_map_plot(gg)
  gg
}

#' @export
plot_on_map = function(data, response_name = NULL, breaks = NULL, use_tex = FALSE, ...) {
  #coords = as.data.frame(sf::st_coordinates(data))
  #data$x = coords$X
  #data$y = coords$Y
  if (!is.null(breaks)) {
    data[[response_name]] = get_binned_data(data[[response_name]], breaks, use_tex = use_tex)
  }
  if (is.null(response_name)) {
    gg = ggplot(data) + geom_sf()
  } else {
    gg = ggplot(data) + geom_sf(aes_string(col = response_name))
  }
  gg = gg +
    #geom_point(aes_string(x = "x", y = "y", col = response_name)) +
    add_norway_map(st_crs(data), st_bbox(data), ...) +
    scale_color_viridis_c() +
    theme_light() +
    labs(x = "", y = "")
  if (use_tex) gg = latex_friendly_map_plot(gg)
  gg
}

latex_friendly_map_plot = function(x) {
  info = ggplot_build(x)$layout$panel_params[[1]]$graticule
  east_ticks = info$degree[info$type == "E"]
  north_ticks = info$degree[info$type == "N"]
  x +
    scale_x_continuous(breaks = east_ticks, labels = paste0(east_ticks, "$^\\circ$E")) +
    scale_y_continuous(breaks = north_ticks, labels = paste0(north_ticks, "$^\\circ$N"))
}

add_scico_fill = function(x, binned = FALSE) {
  if (binned) {
    x = x + scico::scale_fill_scico_d(
      palette = "oslo", direction = -1, end = .9, begin = .1, drop = FALSE)
  } else {
    x = x + scico::scale_fill_scico(
      palette = "oslo", direction = -1, end = .9, begin = .1)
  }
  x
}

get_binned_data = function(x, breaks = NULL, digits = 1, use_tex = FALSE) {
  if (all(is.null(breaks))) breaks = seq(min(x), max(x), length.out = 5)
  if (!is.infinite(breaks[1])) breaks = c(-Inf, breaks)
  if (!is.infinite(breaks[length(breaks)])) breaks = c(breaks, Inf)
  breaks = round(breaks, digits)
  if (use_tex) {
    labels = c(
      paste("$<$", breaks[2]),
      paste(breaks[2:(length(breaks) - 2)], "-", breaks[3:(length(breaks) - 1)]),
      paste("$>$", breaks[length(breaks) - 1]))
  } else {
    labels = c(
      paste("<", breaks[2]),
      paste(breaks[2:(length(breaks) - 2)], "-", breaks[3:(length(breaks) - 1)]),
      paste(">", breaks[length(breaks) - 1]))
  }
  x = cut(x, breaks, labels = labels)
  x
}

add_norway_map = function(crs = NULL, bbox = NULL, ...) {
  map = rnaturalearth::ne_countries(scale = 50, country = "Norway", returnclass = "sf")
  c(layer_sf(
    geom = GeomSf, data = map, mapping = aes(),
    stat = "sf", position = "identity", show.legend = NA,
    inherit.aes = TRUE,
    params = list(na.rm = FALSE, fill = NA, ...)),
    coord_sf(default = TRUE, crs = crs,
             xlim = bbox[c(1, 3)], ylim = bbox[c(2, 4)]))
}
