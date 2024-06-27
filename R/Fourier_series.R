#' @title s4 class of Fourier summation series
#' @description
#' A s4 class that represents the Fourier summation below:
#' \deqn{\frac{a_0}{2} +
#' \sum_{k=1}^{m_a} a_k \cos{(\frac{2\pi}{T}k(x-t_0))} +
#' \sum_{k=1}^{m_b} b_k \sin{(\frac{2\pi}{T}k(x-t_0))},
#' \qquad x\in[t_0,t_0+T]}
#'
#' @slot double_constant value of \eqn{a_0}.
#' @slot cos values of coefficients of \eqn{\cos} waves, \eqn{a_k}.
#' @slot sin values of coefficients of \eqn{\sin} waves, \eqn{b_k}.
#' @slot k_cos values of \eqn{k} corresponding to the coefficients of \eqn{\cos} waves
#' @slot k_sin values of \eqn{k} corresponding to the coefficients of \eqn{\sin} waves
#' @slot t_0 left end of the domain interval, \eqn{t_0}
#' @slot period length of the domain interval, \eqn{T}.
#' @details
#' If not assigned, \eqn{t_0 = 0}, \eqn{T = 2\pi}.
#' If not assigned, k_cos and k_sin equals 1, 2, 3, ...
#' @export
#' @import methods
#' @author Heyang Ji
#' @examples fsc = Fourier_series(
#'            double_constant = 0.5,
#'            cos = c(0,0.3),
#'            sin = c(1,0.7),
#'            k_cos = 1:2,
#'            )
Fourier_series = setClass(
  "Fourier_series",
  slots = c(
    double_constant = "numeric",
    cos             = "numeric",
    sin             = "numeric",
    k_cos           = "num.int.nul",
    k_sin           = "num.int.nul",
    # k_cos           = "ANY",
    # k_sin           = "ANY",
    t_0             = "numeric",
    period          = "numeric"
  )
)
