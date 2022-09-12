qq <- function(...) {
  sapply(match.call()[-1], deparse)
}