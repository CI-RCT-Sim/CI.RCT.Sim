# Internal placeholders for survival::tmerge special terms.
# These are not meant to be called directly.

tdc <- function(...) {
  stop("tdc() is a special helper for survival::tmerge() and should not be called directly.")
}

event <- function(...) {
  stop("event() is a special helper for survival::tmerge() and should not be called directly.")
}
