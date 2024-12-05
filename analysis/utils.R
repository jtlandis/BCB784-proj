

box::use(
  ./dependencies[...]
)

#' @export
merge_cols <- function(x, y, ...) {
  # browser()
  quos <- plyxp:::plyxp_quos(..., .named = F, .ctx_default = "assays",
                             .ctx_opt = "cols")
  quos <- split(lapply(quos, rlang::quo_get_expr),
                vapply(quos, attr, "plyxp:::ctx", FUN.VALUE = ""))
  
  x <- select(x, !!!quos[["assays"]], rows(everything()),
              cols(!!!quos[["cols"]]))
  y <- select(y, !!!quos[["assays"]], cols(!!!quos[["cols"]])) |>
    plyxp(\(y, x) {
      rowRanges(y) <- rowRanges(x)
      y
    }, x = se(x))
  
  se <- cbind(
    se(x), se(y)
  )
  new_plyxp(se)
}