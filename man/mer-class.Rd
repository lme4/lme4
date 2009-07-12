\name{mer-class}
\title{Mixed Model Representations and *mer Methods}
\docType{class}
% Classes
\alias{mer-class}
%\alias{coef.mer-class}
\alias{summary.mer-class}
%\alias{ranef.mer-class}
% Methods
\alias{anova,mer-method}
\alias{coef,mer-method}
\alias{coerce,mer,dtCMatrix-method}
%\alias{confint,mer-method}
\alias{deviance,mer-method}
\alias{expand,mer-method}
\alias{fitted,mer-method}
\alias{fitted,mer-method}
\alias{formula,mer-method}
\alias{logLik,mer-method}
\alias{model.frame,mer-method}
\alias{model.matrix,mer-method}
\alias{print,mer-method}
\alias{resid,mer-method}
\alias{residuals,mer-method}
\alias{show,mer-method}
\alias{simulate,mer-method}
\alias{summary,mer-method}
\alias{terms,mer-method}
\alias{update,mer-method}
\alias{vcov,mer-method}
\alias{with,mer-method}
%
%
\alias{deviance,summary.mer-method}
\alias{logLik,summary.mer-method}
\alias{vcov,summary.mer-method}
\alias{summary,summary.mer-method}

\description{The \code{mer} class represents linear or generalized
  linear or nonlinear mixed-effects models.  It incorporates
  sparse model matrices for the random effects and corresponding sparse
  Cholesky factors.  The \code{summary.mer} class represents the summary
  of these objects.
}

\section{Objects from the Class}{Objects can be created by calls of the
  form \code{new("mer", ...)}  or, more commonly, via the
  \code{\link{lmer}}, \code{\link{glmer}} or \code{\link{nlmer}}
  functions.
}

\section{Slots}{
  The class \code{"mer"} represents a linear or generalized linear or
  nonlinear or generalized nonlinear mixed model and contains the slots:

  \describe{

    \item{\code{env}:}{An environment (class \code{"environment"})
      created for the evaluation of the nonlinear model function. Not
      used except by \code{\link{nlmer}} models.}

    \item{\code{nlmodel}:}{The nonlinear model function as an object of
      class \code{"call"}.  Not used except by \code{\link{nlmer}}
      models.}

    \item{\code{frame}:}{The model frame (class \code{"data.frame"}).}

    \item{\code{call}:}{The matched call to the function that
      created the object. (class \code{"call"}).}

    \item{\code{flist}:}{The list of grouping factors for the random
      effects.}

    \item{\code{X}:}{Model matrix for the fixed effects. In an
      \code{\link{nlmer}} fitted model this matrix has \code{n * s} rows
      where \code{n} is the number of observations and \code{s} is the
      number of parameters in the nonlinear model.}

    \item{\code{Zt}:}{The transpose of model matrix for the random
      effects, stored as a compressed column-oriented sparse matrix (class
      \code{"\link[Matrix:dgCMatrix-class]{dgCMatrix}"}).}

    \item{\code{pWt}:}{Numeric prior weights vector. This may be of length
      zero (0), indicating unit prior weights.}

    \item{\code{offset}:}{Numeric offset vector. This may be of length
      zero (0), indicating no offset.}

    \item{\code{y}:}{The response vector (class \code{"numeric"}).}

    \item{\code{Gp}:}{Integer vector of group pointers within the random
      effects vector.  The elements of \code{Gp} are 0-based indices of
      the first element from each random-effects term.  Thus the first
      element is always 0.  The last element is the total length of the
      random effects vector.}

    \item{\code{dims}:}{A named integer vector of dimensions.  Some of
      the dimensions are \eqn{n}, the number of observations, \eqn{p}, the
      number of fixed effects, \eqn{q}, the total number of random
      effects, \eqn{s}, the number of parameters in the nonlinear model
      function and \eqn{nt}, the number of random-effects terms in the
      model.}

    \item{\code{ST}:}{A list of S and T factors in the TSST' Cholesky
      factorization of the relative variance matrices of the random
      effects associated with each random-effects term.  The unit lower
      triangular matrix, \eqn{T}, and the diagonal matrix, \eqn{S}, for
      each term are stored as a single matrix with diagonal elements
      from \eqn{S} and off-diagonal elements from \eqn{T}.}

    \item{\code{V}:}{Numeric gradient matrix (class \code{"matrix"}) of
      the nonlinear model function. Not used except by
      \code{\link{nlmer}} models.}

    \item{\code{A}:}{Scaled sparse model matrix (class
      \code{"\link[Matrix:dgCMatrix-class]{dgCMatrix}"}) for
      the the unit, orthogonal random effects, \eqn{U}.}

    \item{\code{Cm}:}{Reduced, weighted sparse model matrix (class
      \code{"\link[Matrix:dgCMatrix-class]{dgCMatrix}"}) for the
      unit, orthogonal random effects, U.  Not used except by
      \code{\link{nlmer}} models.}

    \item{\code{Cx}:}{The \code{"x"} slot in the weighted sparse model
      matrix (class \code{"\link[Matrix:dgCMatrix-class]{dgCMatrix}"})
      for the unit, orthogonal random effects, \eqn{U}, in generalized
      linear mixed models.  For these models the matrices \eqn{A} and
      \eqn{C} have the same sparsity pattern and only the \code{"x"}
      slot of \eqn{C} needs to be stored.}

    \item{\code{L}:}{The sparse lower Cholesky factor of \eqn{P(AA'+I)P'}
      (class \code{"\link[Matrix:CHMfactor-class]{dCHMfactor}"}) where \eqn{P}
      is the fill-reducing permutation calculated from the pattern of
      nonzeros in \eqn{A}.}

    \item{\code{deviance}:}{Named numeric vector containing the deviance
      corresponding to the maximum likelihood (the \code{"ML"} element)
      and \code{"REML"} criteria and various components. The
      \code{"ldL2"} element is twice the logarithm of the determinant of
      the Cholesky factor in the \code{L} slot. The \code{"usqr"}
      component is the value of the random-effects quadratic form.}

    \item{\code{fixef}:}{Numeric vector of fixed effects.}

    \item{\code{ranef}:}{Numeric vector of random effects on the
      original scale.}

    \item{\code{u}:}{Numeric vector of orthogonal, constant variance,
      random effects.}

    \item{\code{eta}:}{The linear predictor at the current values of the
      parameters and the random effects.}

    \item{\code{mu}:}{The means of the responses at the current parameter
      values.}

    \item{\code{muEta}:}{The diagonal of the Jacobian of \eqn{\mu}{mu}
      by \eqn{\eta}{eta}.  Has length zero (0) except for generalized
      mixed models.}

    \item{\code{var}:}{The diagonal of the conditional variance of
      \eqn{Y} given the random effects, up to prior weights.  In
      generalized mixed models this is the value of the variance
      function for the \code{\link{glm}} family.}

    \item{\code{resid}:}{The residuals, \eqn{y - \mu}{y-mu}, weighted by
      the \code{sqrtrWt} slot (when its length is \eqn{>0}).}

    \item{\code{sqrtXWt}:}{The square root of the weights applied to the
      model matrices \eqn{X} and \eqn{Z}. This may be of length zero
      (0), indicating unit weights.}

    \item{\code{sqrtrWt}:}{The square root of the weights applied to the
      residuals to obtain the weighted residual sum of squares. This may
      be of length zero (0), indicating unit weights.}

    \item{\code{RZX}:}{The dense solution (class \code{"matrix"}) to
      \eqn{L RZX = ST'Z'X = AX}.}

    \item{\code{RX}:}{The upper Cholesky factor (class \code{"matrix"})
      of the downdated \eqn{X'X}.}

  }

  The \code{"summary.mer"} class \emph{contains} the \code{"mer"},
  class and has additional slots,
  \describe{
    
    \item{\code{methTitle}:}{character string specifying a method title}
    \item{\code{logLik}:}{the same as \code{\link{logLik}(object)}.}
    \item{\code{ngrps}:}{the number of levels per grouping factor in the
      \code{flist} slot.}
    \item{\code{sigma}:}{the scale factor for the variance-covariance estimates}
    \item{\code{coefs}:}{the matrix of estimates, standard errors,
      etc. for the fixed-effects coefficients}
    \item{\code{vcov}:}{the same as \code{\link{vcov}(object)}.}
    \item{\code{REmat}:}{the formatted Random-Effects matrix}
    \item{\code{AICtab}:}{A named vector of values of AIC, BIC, log-likelihood
      and deviance}
  }% summary slots
}
\usage{
## Methods with "surprising" arguments
\S4method{deviance}{mer}(object, REML = NULL, \dots)
\S4method{expand}{mer}(x, sparse = TRUE, \dots)
\S4method{logLik}{mer}(object, REML = NULL, \dots)
\S4method{print}{mer}(x, digits, correlation, symbolic.cor, signif.stars, \dots)
}
\arguments{
  \item{object}{object of class \code{mer}.}
  \item{REML}{logical indicating if REML should be used. A value of
    \code{NULL}, the default, or \code{NA} indicates that the REML values
    should be returned if the model was fit by REML, otherwise the ML values.}
  \item{x}{object of class \code{mer} to expand.}
  \item{sparse}{logical scalar indicating if the sparse form of the
    expanded \code{T} and \code{S} matrices should be returned.}
  \item{digits}{number of digits to use when printing tables of
    parameter estimates.  Defaults to \code{max(3, getOption("digits") -
      3)}.}
  \item{correlation}{logical - should the correlation matrix of the
    fixed-effects parameter estimates be printed?  Defaults to \code{TRUE}.}
  \item{symbolic.cor}{logical - should a symbolic form of the
    correlation matrix be printed instead of the numeric form?  Defaults
    to \code{FALSE}.}
  \item{signif.stars}{logical - should the \sQuote{significance stars}
    be printed as part of the table of fixed-effects parameter
    estimates?  Defaults to \code{getOption("show.signif.stars")}.}
  \item{\dots}{potential further arguments passed to methods.}
}
\section{Methods}{
  \describe{
    \item{VarCorr}{\code{signature(x = "mer")}: Extract variance and
      correlation components.  See \code{\link{VarCorr}}}
    \item{anova}{\code{signature(object = "mer")}: returns the sequential
      decomposition of the contributions of fixed-effects terms or, for
      multiple arguments, model comparison statistics. See
      \code{\link{anova}}.}
    %% FIXME: say a bit more on lmer - model comparison ??
    \item{coef}{\code{signature(object = "mer")}: returns an object
      similar to the \code{\link{ranef}} method but incorporating the
      fixed-effects parameters, thereby forming a table of linear model
      coefficients (the columns) by level of the grouping factor (the rows).}
    \item{coerce}{\code{signature(from = "mer", to = "dtCMatrix")}:
      returns the \code{L} slot as a \code{"dtCMatrix"}
      %% FIXME - this produces a warning of a missing link in R-2.7.0
      %\code{"\linkS4class{dtCMatrix}"}
      (column-oriented, sparse, triangular matrix) object.}
    \item{deviance}{\code{signature(object = "mer")}: returns the
      \code{\link{deviance}} of the fitted model, or the \dQuote{REML
	deviance} (i.e. negative twice the REML criterion), according to
      the REML argument.  See the arguments section above for a description
      of the REML argument.}
    \item{expand}{\code{signature(object = "mer")}:

      returns a list of terms in the expansion of the \code{ST} slot.
      If \code{sparse} is \code{TRUE}, the default, the elements of the
      list are the numeric scalar \code{"sigma"}, the REML or ML
      estimate of the standard deviation in the model, and three sparse
      matrices: \code{"P"}, the permutation matrix, \code{"S"}, the
      diagonal scale matrix and \code{"T"}, the lower triangular matrix
      determining correlations.  When \code{sparse} is \code{FALSE} each
      element of the list is the expansions of the corresponding element
      of the \code{ST} slot into a list of \code{S}, the diagonal
      matrix, and \code{T}, the (dense) unit lower triangular matrix.

    }
    \item{fitted}{\code{signature(object = "mer")}:

      returns the fitted conditional means of the responses. See
      \code{\link{fitted}}.  The \code{\link{napredict}} function is
      called to align the result with the original data if the model was
      fit with \code{na.action = na.exclude}.

    }
    \item{fixef}{\code{signature(object = "mer")}:

      returns the estimates of the fixed-effects parameters.  See
      \code{\link{fixef}}.
    
    }
    \item{formula}{\code{signature(x = "mer")}:

      returns the model formula.  See \code{\link{formula}}.

    }
%%% some of the following are also in \usage{..} because they have extra args:
    \item{logLik}{\code{signature(object = "mer")}:

      returns the log-likelihood or the REML criterion, according to the
      optional \code{REML} argument (see the arguments section above),
      of the fitted model.  See also \code{\link{logLik}}.

    }
    \item{mcmcsamp}{\code{signature(object = "mer")}:

      Create a Markov chain Monte Carlo sample from a posterior
      distribution of the model's parameters.  See
      \code{\link{mcmcsamp}} for details.

    }
    \item{model.frame}{\code{signature(formula = "mer")}: returns the
      model frame (the \code{frame} slot).}
    \item{model.matrix}{\code{signature(object = "mer")}: returns the
      model matrix for the fixed-effects parameters (the \code{X}
      slot).} 
    \item{print}{\code{signature(x = "mer")}: print information about
      the fitted model.  See the arguments section above for a description
      of optional arguments.}
    \item{ranef}{\code{signature(object = "mer")}: returns the
      conditional modes of the random effects.  See \code{\link{ranef}}.}
    \item{refit}{\code{signature(object = "mer", newresp = "numeric")}:
      Update the response vector only and refit the model. See
      \code{\link{refit}}.}
    \item{resid}{\code{signature(object = "mer")}: returns the (raw)
      residuals.  This method calls \code{\link{napredict}}.  See the
      above description of the \code{fitted} method for details.  See
      also \code{\link{resid}}.}
    \item{residuals}{\code{signature(object = "mer")}: Another name
      for the \code{resid} method.}
    \item{show}{\code{signature(object = "mer")}: Same as the
      \code{print} method without the optional arguments.}
    \item{simulate}{\code{signature(object = "mer")}: simulate
      \code{nsim} (defaults to 1) responses from the theoretical
      distribution corresponding to the fitted model.  The \code{refit}
      method is particularly useful in combination with this method.
      See also \code{\link{simulate}}.}
    \item{terms}{\code{signature(x = "mer")}: Extract the 
      \code{\link{terms}} object for the fixed-effects terms in the
      model formula.}
    \item{update}{\code{signature(object = "mer")}: see
      \code{\link{update}} on how to update fitted models.}
    \item{vcov}{\code{signature(object = "mer")}: Calculate
      variance-covariance matrix of the \emph{fixed} effect terms,
      see also \code{\link{vcov}}.}
    \item{with}{\code{signature(data = "mer")}: Evaluate an R expression
      in an environment constructed from the \code{frame} slot.}
  }
}
\seealso{
  \code{\link{lmer}()}, \code{\link{glmer}()} and \code{\link{nlmer}()},
  which produce these objects.
  \cr
  \code{\link{VarCorr}} for extracting the variance and
  correlation components of the \emph{random}-effects terms.
}
\examples{
(fm2 <- lmer(Reaction ~ Days + (1|Subject) + (0+Days|Subject),
             data = sleepstudy))
print(fm2, digits = 10, corr = FALSE) # more precision; no corr.matrix

logLik(fm2)
(V2 <- vcov(fm2))
terms(fm2)
str(model.matrix(fm2))
str(model.frame(fm2))
str(resid(fm2))

VarCorr(fm2)
ee <- expand(fm2)
op <- options(digits = 3)
tcrossprod(ee$sigma * ee$P \%*\% ee$T \%*\% ee$S)
options(op)

\dontrun{
## Simulate 'Reaction' according to the fitted model:
dim(ss <- simulate(fm2, nsim = 200, seed = 101)) ## -> 180 x 200
}

}
\keyword{classes}
