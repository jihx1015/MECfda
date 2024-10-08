% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fcRegression.R
\name{fcRegression}
\alias{fcRegression}
\title{Solve linear models with functional covariate(s)}
\usage{
fcRegression(
  Y,
  FC,
  Z,
  formula.Z,
  family = gaussian(link = "identity"),
  basis.type = c("Fourier", "Bspline"),
  basis.order = 6L,
  bs_degree = 3
)
}
\arguments{
\item{Y}{Response variable, can be an atomic vector, a one-column matrix or data frame,
recommended form is a one-column data frame with column name.}

\item{FC}{Functional covariate(s),
can be a "functional_variable" object or a matrix or a data frame or a list of these object(s).}

\item{Z}{Scalar covariate(s), can be NULL or not input (when there's no scalar covariate),
an atomic vector (when only one scalar covariate), a matrix or data frame,
recommended form is a data frame with column name(s).}

\item{formula.Z}{A formula without the response variable,
contains only scalar covariate(s) (or intercept),
use the format of lme4 package if random effects exist. e.g. ~ Z_1 + (1|Z_2).
(See \code{\link[lme4]{lmer}} and \code{\link[lme4]{glmer}})
If not assigned, include all scalar covariates and intercept term as fixed effects.}

\item{family}{A description of the error distribution and link function to be used in the model,
see \code{\link[stats]{family}}.}

\item{basis.type}{Type of funtion basis.
Can only be assigned as one type even if there is more than one functional covariates.
Available options: 'Fourier' or 'Bspline', represent Fourier basis and b-spline basis respectively.
For the detailed form for Fourier and b-splines basis,
see \code{\link{fourier_basis_expansion}} and \code{\link{bspline_basis_expansion}}.}

\item{basis.order}{Indicate number of the function basis.
When using Fourier basis \eqn{\frac{1}{2},\sin k t, \cos k t, k = 1,\dots,p_f},
basis.order is the number \eqn{p_f}.
When using b-splines basis \eqn{\{B_{i,p}(x)\}_{i=-p}^{k}},
basis.order is the number of splines, equal to \eqn{k+p+1}.
(same as arguement \code{df} in \code{\link[splines]{bs}}.)
May set a individual number for each functional covariate.
When the element of this argument is less than the number of functional covariates,
it will be used recursively.}

\item{bs_degree}{Degree of the piecewise polynomials if use b-splines basis,
default is 3. See \code{degree} in \code{\link[splines]{bs}}.}
}
\value{
fcRegression returns an object of class "fcRegression".
    It is a list that contains the following elements.
    \item{regression_result}{Result of the regression.}
    \item{FC.BasisCoefficient}{A list of Fourier_series or bspline_series object(s),
    represents the functional linear coefficient(s) of the functional covariates.}
    \item{function.basis.type}{Type of funtion basis used.}
    \item{basis.order}{Same as in the arguemnets.}
    \item{data}{Original data.}
    \item{bs_degree}{Degree of the splines, returned only if b-splines basis is used.}
}
\description{
Function to fit (generalized) linear model with functional covariate(s).
Model allows one or multiple functional covariate(s) as fixed effect(s), and zero,
one, or multiple scalar-valued fixed or random effect(s).
}
\details{
Solve linear models with functional covariates below
\deqn{g(E(Y_i|X_i,Z_i)) = \sum_{l=1}^{L} \int_{\Omega_l} \beta_l(t) X_{li}(t) dt + (1,Z_i^T)\gamma}
where the scalar-valued covariates can be fixed or random effect or doesn't exist
(may do not contain scalar-valued covariates).
}
\examples{
data(MECfda.data.sim.0.0)
res = fcRegression(FC = MECfda.data.sim.0.0$FC, Y=MECfda.data.sim.0.0$Y, Z=MECfda.data.sim.0.0$Z,
                   basis.order = 5, basis.type = c('Bspline'),
                   formula.Z = ~ Z_1 + (1|Z_2))
}
\author{
Heyang Ji
}
