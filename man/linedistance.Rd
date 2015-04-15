\name{linedistance}
\alias{linedistance}
\title{Intersection Points}

\description{
This function allows to calculate the intersection points between the ROC curve and the sampling lines. Also calculates the distance between this points and the reference point.
}
\usage{
linedistance(curve.fpr, curve.tpr, curve.segslope, curve.slope, line.slope, ref.point)
}

\arguments{
  \item{curve.fpr}{
False positive rate vector with all points of the given Curve
}
  \item{curve.tpr}{
True positive rate vector with all points of the given Curve
}
  \item{curve.segslope}{
Vector with all segments slope of the ROC curves
}
  \item{curve.slope}{
Vector with all the slope of all segments that connect the ROC curve with the reference point
}
  \item{line.slope}{
Vector with the slope of all sampling lines
}
  \item{ref.point}{
Reference point where we start drawing the sampling lines
}
}

\value{
This function returns a list with:
\item{dist}{Vector with distances between the intersection points and the reference points}
\item{x}{Vector with all x coordinates of intersection points}
\item{y}{Vector with all y coordinates of intersection points}
}

\seealso{
\code{\link{lineslope}}
\code{\link{curvesegslope}}
\code{\link{curvesegsloperef}}
}

\keyword{ Sampling }
\keyword{ Points }
\keyword{ Intersection }
\keyword{ Lines }
\keyword{ ROC }
