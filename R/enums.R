enum <- function(...) {
  nms <- eval(substitute(alist(...)))
  x <- as.list(setNames(seq_along(nms), nms))
  return(x)
}

#' @export
variance_types <- enum(homoskedastic, heteroskedastic)
