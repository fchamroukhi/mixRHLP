enum <- function(...) {
  nms <- eval(substitute(alist(...)))
  x <- as.list(setNames(seq_along(nms), nms))
  x
}

variance_types <- enum(free, common)
