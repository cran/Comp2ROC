\name{curvesegsloperef}
\alias{curvesegsloperef}
\title{Segment Slopes to Reference Point}

\description{
This function allows to calculate the segments slope that connect the ROC curve segments with the reference point (1,0).
}
\usage{
curvesegsloperef(curve.fpr, curve.tpr, ref.point)
}

\arguments{
  \item{curve.fpr}{
False positive rate vector with all points of the given Curve
}
  \item{curve.tpr}{
True positive rate vector with all points of the given Curve
}
  \item{ref.point}{
Reference point where we start drawing the sample lines
}
}

\value{
This function returns a vector with all segments slopes that connect the ROC curve points to the reference point.
}

\keyword{ Curve }
\keyword{ Slope }
\keyword{ ROC }
\keyword{ Segment }
