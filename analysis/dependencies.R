
#' @return logical TRUE if package exists, FALSE otherwise
package_exists <- function(package) {
  tryCatch({find.package(package = package); TRUE},
           error = function(cnd) {
             rlang::cnd_muffle(cnd)
             FALSE
           })
}

#' @description
#' installs packages for this project. There may be an issue
#' with `plyxp`
#' 
require_packages <- function(packages, repo = c("CRAN", "Bioconductor")) {
  repo <- match.arg(repo, c("CRAN", "Bioconductor"))
  for (package in packages) {
    if (!package_exists(package)) {
      res <- readline(paste0("Package ", package, " not found. Install? [y/n]"))
      if (grepl("^y", res, ignore.case = TRUE)) {
        switch(repo,
               "Bioconductor" = {
                 if (!package_exists("BiocManager")) {
                   install.packages("BiocManager")
                 } 
                 BiocManager::install(package)
               },
               "CRAN" = {
                 install.packages(package)
               })
      } else {
        warning(paste0("Package ", package, " not installed."))
      }
    }
  }
}

require_packages(
  packages = c("dplyr", "readr", "purrr", "readxl", "stringr"),
  repo = "CRAN")
require_packages(
  packages = c("SummarizedExperiment", "plyxp",
               "plyranges", "rtracklayer", "GenomicRanges"),
  repo = "Bioconductor")


#' @name .__module__.
#' @description
#' Required packages that are exported with this module.
#' `plyxp` may not be available for the user, in which
#' case they should either upgrade their `BiocManager`
#' version to 3.20, or use `remotes::install_github("jtlandis/plyxp")`
#' 
#' @export
box::use(
  SummarizedExperiment[...],
  dplyr[...],
  plyranges[...],
  plyxp[...],
  GenomicRanges[...]
)

box::use(
  rlang[enquo, ensym, enquos, new_quosure,
        expr, quo_get_env, as_function,
        eval_tidy]
)

#' @export
plyxp_on <- function(.data, .f, ..., .on) {
  .data <- enquo(.data)
  .f <- enquo(.f)
  .on <- ensym(.on)
  dots <- enquos(...)
  quo <- new_quosure(
    expr({
      plyxp(!!.data, \(se, ...) {
        .f <- rlang::as_function(!!.f)
        obj <- (!!.on)(se)
        obj <- .f(obj, ...)
        (!!.on)(se) <- obj
        se
      }, !!!dots)
    }),
    rlang::quo_get_env(.data))
  eval_tidy(quo)
}
