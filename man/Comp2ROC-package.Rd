\name{Comp2ROC-package}
\alias{Comp2ROC-package}
\alias{Comp2ROC}
\docType{package}
\title{
Comparation of two ROC curves that intersect.
}
\description{
Comparation of ROC curves using the methodology devoloped by Braga.
}
\details{
\tabular{ll}{
Package: \tab Comp2ROC\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2012-06-26\cr
License: \tab GPL-2\cr
}
}
\author{
Ana C. Braga, with contributions from Hugo Frade

Maintainer: Ana C. Braga <acb@dps.uminho.pt>
}
\references{
BRAGA, A. AND COSTA, L. AND OLIVEIRA, P. 2011. An alternative method for global and partial comparasion of two diagnostic system based on ROC curves In Journal of Statistical Computation and Simulation.
}

\keyword{ package }
\keyword{Comp2ROC}
\keyword{compare}
\keyword{ROC}

\examples{
	
#This is a simple example how to use the package:
#nameE="new_Zhang"
#moda1="mod1" 
#moda2="mod2"
#data(zhang)
#data=read.manually.introduced(zhang, moda1,TRUE, moda2, TRUE, "status", TRUE)
#results=roc.curves.boot(data,name=nameE,mod1= moda1,mod2=moda2) 
#rocboot.summary(results,moda1,moda2) 
#save.file.summary(results,nameE,app=TRUE,moda1,moda2) 

}
