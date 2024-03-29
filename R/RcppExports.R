# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#' Pack circles
#'
#' This function returns a dataframe of x and y coordinates and radii
#' of circles packed into an arbitrary polygon container.
#'
#' @param polygon A dataframe containing the x and y coordinates of vertices
#'   a polygon into which circles will be packed
#' @param radii A numeric vector giving the radii of circles which will be
#'   packed into the container specified in the polygon argument. This also
#'   determines the maximum number of circles which the algorithm will attempt
#'   to place.
#' @param existing_circles A dataframe of x and y coordinates and radii of any
#'   existing circles to be avoided when packing the new ones. This allows
#'   multiple iterations, for example to pack a shape within another shape.
#' @param max_attempts A scalar integer. How many times the algorithm will try
#'   to place a circle in a new spot before giving up.
#' @param seed A scalar integer. Used to seed the random number generator.
#' @param neat_edges Logical. If true, no part of a circle can extend beyond
#'   the container polygons boundary. If false, circles can extend beyond the
#'   boundary as long as any part of the circle remains within the container.
#' @export
pack_circles <- function(polygon, radii, existing_circles, max_attempts = 2000L, seed = 1L, neat_edges = TRUE) {
    .Call(`_pack_pack_circles`, polygon, radii, existing_circles, max_attempts, seed, neat_edges)
}

pack_circles_all <- function(polygon, radii, existing_circles, max_attempts = 200L, seed = 1L, neat_edges = TRUE) {
    .Call(`_pack_pack_circles_all`, polygon, radii, existing_circles, max_attempts, seed, neat_edges)
}

pack_polygons <- function(polygon, radii, sides, max_attempts = 2000L, seed = 1L, neat_edges = TRUE) {
    .Call(`_pack_pack_polygons`, polygon, radii, sides, max_attempts, seed, neat_edges)
}

