
# This script contains various utility functions

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

#' @export
get_proj_xy = function() {
  sf::st_crs("+proj=utm +zone=33 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=km +no_defs")
}

#' @export
tikz_plot = function(file, plot = NULL, expression = NULL, view = interactive(), ...) {
  # Ensure that you are on an operating system that you have tested
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

  # Create a temporary file for the tikz-output
  tmp = tempfile(tmpdir = getwd())
  # Clean up after yourself on early interrupt
  on.exit(suppressWarnings(file.remove(tmp)), add = TRUE)

  # Extract default tex usepackages and add the bm package for bold greek letters
  opt = options()
  on.exit(options(opt)) #Reset global options on exit
  tikzDevice::setTikzDefaults(overwrite = FALSE)
  tex_packages = options()$tikzLatexPackages
  if (!any(grepl("usepackage\\{bm\\}", tex_packages))) {
    tex_packages = c(tex_packages, "\\usepackage{bm}\n")
  }

  # Open a device for creating a tex-file
  tikzDevice::tikz(tmp, standAlone = TRUE, packages = tex_packages, ...)
  # Call dev.off() on exit in case of interruptions
  current_device = dev.cur()
  on.exit(dev.off(current_device))

  # Plot something into the tex-file
  if (!is.null(plot)) {
    if (any(class(plot) %in% c("gg", "ggplot", "patchwork"))) {
      print(plot)
    } else {
      for (p in plot) print(p)
    }
  } else {
    eval(substitute(expression), envir = parent.frame())
  }

  # Finish the creation of the tex-file
  dev.off()

  # Compile to pdf
  #system2("pdflatex", tmp)
  system2("lualatex", shQuote(tmp))

  # Copy pdf file to final destination
  file.copy(paste0(tmp, ".pdf"), file, overwrite = TRUE)

  # Clean up all temporary files
  tmp_filename = tail(strsplit(tmp, "/")[[1]], 1)
  files_to_clean = grep(tmp_filename, list.files(full.names = TRUE), value = TRUE)
  file.remove(files_to_clean)

  # Open the pdf with the final output
  if (view) {
    if (operating_system == "Darwin") {
      system2("open", shQuote(file))
    } else {
      message("I don't know the command to open a pdf for your operating system")
    }
  }
}
