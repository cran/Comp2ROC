\name{roc.curves.boot}
\alias{roc.curves.boot}
\title{Compare curves}

\description{
This is the function which control the whole package.This uses all functions except the reading ones and \code{rocboot.summary} and \code{save.file.summary}.
}
\usage{
roc.curves.boot(data, nb = 1000, alfa = 0.05, name, mod1, mod2, related)
}

\arguments{
  \item{data}{
Data obtained throught \code{read.file} or \code{read.manually.introduced}
}

  \item{nb}{
Number of permutations
}
  \item{alfa}{
Confidance level for parametric methods
}
  \item{name}{
Name too show in graphs
}
  \item{mod1}{
Name of Modality 1
}
  \item{mod2}{
Name of Modality 2
}
  \item{related}{
Boolean parameter that represents if the two modalities are related or not
}
}

\value{
This function returns a list with:
\item{Area1 }{Area of Curve 1}
\item{SE1 }{Standard error of Curve 1}
\item{Area2 }{Area of Curve 2}
\item{SE2 }{Standard error of Curve 2}
\item{CorrCoef }{Correlation Coeficient}
\item{diff }{Difference Between Areas (TS)}
\item{zstats }{Z Statistic}
\item{pvalue1 }{p-value of Z Statistics}
\item{TrapArea1 }{Area of curve 1 using the Trapezoidal rule}
\item{TrapArea2 }{Area of curve 2 using the Trapezoidal rule}
\item{bootpvalue }{p-value of bootstrapping}
\item{nCross }{Number of Crossings}
\item{ICLB1 }{Confidance Interval: Lower Bound for Curve 1}
\item{ICUB1 }{Confidance Interval: Upper Bound for Curve 1}
\item{ICLB2 }{Confidance Interval: Lower Bound for Curve 2}
\item{ICUB2 }{Confidance Interval: Upper Bound for Curve 2}
\item{ICLBDiff }{Confidance Interval: Lower Bound for Difference between areas}
\item{ICUBDiff }{Confidance Interval: Upper Bound for Difference between areas}
}

\examples{
\donttest{
data(zhang)
nameE = "new_Zhang"
modality1DataColumn = "modality1"
modality2DataColumn = "modality2"
data = read.manually.introduced(zhang, moda1, TRUE, moda2, TRUE, "status", TRUE)
results = roc.curves.boot(zhang, 1000, 0.05, name=nameE,
                          mod1=modality1DataColumn, mod2=modality2DataColumn)
}
}

\keyword{ Permutation }
\keyword{ Comparation }
\keyword{ ROC }
\keyword{ Curves }
