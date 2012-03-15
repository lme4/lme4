\docType{data}
\name{InstEval}
\alias{InstEval}
\title{University Lecture/Instructor Evaluations by Students at ETH}
\format{A data frame with 73421 observations on the following 7 variables.
  \describe{
    \item{\code{s}}{a factor with levels \code{1:2972} denoting
      individual students.}
    \item{\code{d}}{a factor with 1128 levels from \code{1:2160}, denoting
      individual professors or lecturers.}% ("d": \dQuote{Dozierende} in German)
    \item{\code{studage}}{an ordered factor with levels \code{2} <
      \code{4} < \code{6} < \code{8}, denoting student's \dQuote{age}
      measured in the \emph{semester} number the student has been enrolled.}
    \item{\code{lectage}}{an ordered factor with 6 levels, \code{1} <
      \code{2} < ... < \code{6}, measuring how many semesters back the
      lecture rated had taken place.}
    \item{\code{service}}{a binary factor with levels \code{0} and
      \code{1}; a lecture is a \dQuote{service}, if held for a
      different department than the lecturer's main one.}
    \item{\code{dept}}{a factor with 14 levels from \code{1:15}, using a
      random code for the department of the lecture.}

    \item{\code{y}}{a numeric vector of \emph{ratings} of lectures by
      the students, using the discrete scale \code{1:5}, with meanings
      of \sQuote{poor} to \sQuote{very good}.}
  }
  Each observation is one student's rating for a specific lecture
  (of one lecturer, during one semester in the past).}
\description{
  University lecture evaluations by students at ETH Zurich,
  anonymized for privacy protection.  This is an
  interesting \dQuote{medium} sized example of a
  \emph{partially} nested mixed effect model.
}
\details{
  The main goal of the survey is to find \dQuote{the best
  liked prof}, according to the lectures given.
  Statistical analysis of such data has been the basis for
  a (student) jury selecting the final winners.

  The present data set has been anonymized and slightly
  simplified on purpose.
}
\examples{
str(InstEval)

head(InstEval, 16)
xtabs(~ service + dept, InstEval)
}
\keyword{datasets}

