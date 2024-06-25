#' @title Compute the value of the Fourier summation series at certain points.
#'
#' @param object an object of \code{\link{Fourier_series}} class.
#' @param t Value of t.
#'
#' @return A numeric atomic vector
#' @export
#' @author Heyang Ji
#' @examples fsc = Fourier_series(
#'            double_constant = 0.5,
#'            cos = c(0,0.3),
#'            sin = c(1,0.7),
#'            k_cos = 1:2,
#'            )
#'           FourierSeries2fun(fsc,1:5)
#'
setGeneric("FourierSeries2fun",
           function(object,t) standardGeneric("FourierSeries2fun")
           # ,valueClass = "numeric"
)

#' @rdname FourierSeries2fun
#' @export
setMethod("FourierSeries2fun",
          signature(object="Fourier_series",
                    t = "numeric"),

          function(object, t){
            x1 = 2*pi*((object@k_cos) %*% t((t - object@t_0)/object@period))
            x2 = 2*pi*((object@k_sin) %*% t((t - object@t_0)/object@period))
            return(object@double_constant /2 + colSums(object@cos * cos(x1) + object@sin * sin(x2)))
          })
