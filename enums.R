enum <- function(...) {
  nms <- eval(substitute(alist(...)))
  x <- as.list(setNames(seq_along(nms), nms))
  x
}

# to enumerate the variance types
variance_types <- enum(free, common)

#algorithms <- enum(em, cem)
