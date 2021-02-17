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

fix_lengths = function(...) {
  call = match.call()
  varnames = sapply(call[-1], as.character)
  e = parent.frame()
  vars = lapply(varnames, get, envir = e)
  lengths = sapply(vars, length)
  max_length = max(lengths)
  if (any(max_length %% lengths != 0)) stop("Bad input lengths")
  for (i in seq_along(vars)) {
    if (lengths[i] < max_length) {
      assign(varnames[i], rep(vars[[i]], max_length / lengths[i]), envir = e)
    }
  }
  0
}

get_progress_bar = function(n) {
  pb = progress::progress_bar$new(
    format = ":percent [:bar] time elapsed: :elapsedfull, eta: :eta",
    total = n + 1, width = 70, clear = FALSE)
  pb$tick()
  pb
}

get_proj_xy = function() {
  sf::st_crs("+proj=utm +zone=33 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=km +no_defs")
}

#' @export
tikz_plot = function(file, expression, view = FALSE, ...) {
  operating_system = Sys.info()[["sysname"]]
  if (operating_system == "Windows") {
    proceed = readline(paste("This function was written on a Mac.",
                             "I have no idea if it will work on Windows.",
                             "Proceed? (y/n) "))
    if (proceed == "n") {
      return()
    } else if (proceed != "y") {
      stop("Invalid input")
    }
  }
  tmp = tempfile(tmpdir = getwd())
  tikzDevice::tikz(tmp, standAlone = TRUE, ...)
  eval(substitute(expression), envir = parent.frame())
  dev.off()
  system2("pdflatex", tmp)
  file.copy(paste0(tmp, ".pdf"), file, overwrite = TRUE)
  tmp_filename = tail(strsplit(tmp, "/")[[1]], 1)
  files_to_clean = grep(tmp_filename, list.files(full.names = TRUE), value = TRUE)
  file.remove(files_to_clean)
  if (view) {
    if (operating_system == "Darwin") {
      system2("open", file)
    } else {
      message("I don't know the command to open a pdf for your operating system")
    }
  }
}
