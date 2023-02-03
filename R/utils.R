#' Title
#'
#' @param n
#' @param min
#' @param max
#' @param rate
#' @param power
#'
#' @return
#' @export
#'
#' @examples
make_radii <- function(n, min = 3, max = 50, rate = .1, power = 3) {

  exp <- rexp(n, rate = rate)^power |>
    sort() |>
    rev()

  radii <- rescale(exp, min = min,max = max) |>
    floor()

  radii
}




rescale <- function(x, min = 1, max = 100) {
  (x/max(x))*(max-min) +min
}
