#' @title num.int.nul
#' @import methods
#' @export
#' @exportClass num.int.nul
#' @author Heyang Ji
setClassUnion('num.int.nul',c("numeric","integer","NULL"))
