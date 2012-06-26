\name{rocsampling}
\alias{rocsampling}
\title{
ROC Sampling
}
\description{
This function allows to calculate some statistical measures like extension and location.
}
\usage{
rocsampling(curve1.fpr, curve1.tpr, curve2.fpr, curve2.tpr, K = 100)
}

\arguments{
  \item{curve1.fpr}{
False positive rate vector with all points of the Curve 1
}
  \item{curve1.tpr}{
True positive rate vector with all points of the Curve 1
}
  \item{curve2.fpr}{
False positive rate vector with all points of the Curve 2
}
  \item{curve2.tpr}{
True positive rate vector with all points of the Curve 2
}
  \item{K}{
Number of sampling lines
}
}
\details{
This function uses functions like \code{areatriangles}, \code{curvesegslope}, \code{curvesegsloperef}, \code{diffareatriangles}, \code{linedistance} and \code{lineslope to calculate} that measures.
By default the number of sampling lines is 100, beacause it was proved by Braga that it was the optimal number.
}
\value{
This funcion returns a list with the following components:
\item{AUC1}{Total Area of Curve 1 (using triangles)}
\item{AUC2}{Total Area of Curve 2 (using triangles)}
\item{propc1}{Proportion of Curve1}
\item{propc2}{Proportion of Curve2}
\item{propties}{Proportion of ties}
\item{locc1}{Location of Curve 1}
\item{locc2}{Location of Curve 2}
\item{locties}{Location of Ties}
\item{K}{Number of sampling lines}
\item{lineslope}{Slopes of sampling lines}
\item{diffareas}{Difference of area of triangles}
\item{dist1}{Distance of the intersection points of Curve 1 to reference point}
\item{dist2}{Distance of the intersection points of Curve 2 to reference point}
}


\seealso{
\code{\link{areatriangles}}
\code{\link{curvesegslope}}
\code{\link{curvesegsloperef}}
\code{\link{diffareatriangles}}
\code{\link{linedistance}}
\code{\link{lineslope}}

}


\keyword{ ROC }
\keyword{ Sampling }
