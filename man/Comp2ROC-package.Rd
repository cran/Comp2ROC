\name{Comp2ROC-package}
\alias{Comp2ROC-package}
\alias{Comp2ROC}
\docType{package}
\title{Comparation of Two ROC Curves that Intersect}
\description{
Comaparation of ROC Curves using the methodology devoloped by Braga.
}
\details{
\tabular{ll}{
Package: \tab Comp2ROC\cr
Type: \tab Package\cr
Version: \tab 1.1.2\cr
Date: \tab 2016-05-18\cr
License: \tab GPL-2\cr
}
}
\author{
Ana C. Braga, with contributions from Hugo Frade, Sara Carvalho and Andre M Santiago.

Maintainer: Ana C. Braga <acb@dps.uminho.pt>; Andre M. Santiago <andreportugalsantiago@gmail.com>; 
}
\references{
BRAGA, A. AND COSTA, L. AND OLIVEIRA, P. 2011. An alternative method for global and partial comparasion of two diagnostic system based on ROC curves In Journal of Statistical Computation and Simulation.
}

\keyword{ package }
\keyword{Comp2ROC}
\keyword{compare}
\keyword{ROC}


\examples{
# This is a simple example how to use the package with the given dataset ZHANG:
nameE = "Zhang"
modality1DataColumn = "modality1"
modality2DataColumn = "modality2"
data(zhang)
results = roc.curves.boot(zhang, 10, 0.05, name=nameE,
                          mod1=modality1DataColumn, mod2=modality2DataColumn)
\donttest{rocboot.summary(results, "modality1", "modality2")}
}
