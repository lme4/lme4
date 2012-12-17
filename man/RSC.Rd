\name{RSC}
\alias{RSC}
\title{Class "RSC" - regular sparse column representation}
\description{
  In the regular sparse column representation of a
  mixed-effects model predictor the Z and X matrices are
  represented by two dense matrices - a k by n integer
  matrix, where k is the number of non-zeros in each row of
  Z and n is the total number of observations, and a (k+p)
  by n numeric matrix where the first k rows are the
  non-zeros in Zt and the last p rows are Xt.  Also stored
  are variance-component parameter vector, theta, and the
  vector of lower bounds on these parameters, lower, a
  dsCMatrix of size (q+p) and the (q+p) numeric coefficient
  vector, ubeta.
}

