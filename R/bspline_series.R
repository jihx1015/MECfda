#' @title b-splines summation series.
#' @description
#' A s4 class that represents
#' the summation \eqn{\sum_{i=0}^{k}b_i B_{i,p}(x)} by a bspline_basis object
#' and coefficients \eqn{b_i} (\eqn{i = 0,\dots,k}).
#'
#' @slot coef coefficients of the b-splines, \eqn{b_i} (\eqn{i = 0,\dots,k}).
#' @slot bspline_basis a \code{\link{bspline_basis}} object,  represents the b-splines basis used, \eqn{\{B_{i,p}(x)\}_{i=-p}^{k}}.
#'
#' @export
#' @import methods
#' @author Heyang Ji
#' @examples bsb = bspline_basis(
#'             Boundary.knots = c(0,24),
#'             intercept      = TRUE,
#'             df             = NULL,
#'             degree         = 3
#' )
#' bss = bspline_series(
#'           coef = c(2,1,1.5,0.5),
#'           bspline_basis = bsb
#' )
bspline_series = setClass(
  "bspline_series",
  slots = c(
    coef = "numeric",
    bspline_basis = "bspline_basis")
)
