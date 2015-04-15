\name{areatriangles}
\alias{areatriangles}
\title{Triangle Areas}
\description{
This function allows to calculate the triangles area formed with two points that was next to each other and the reference point. It also allows to calculate the total area based on the previous triangles.
}
\usage{
areatriangles(line.slope, line.dist1)
}

\arguments{
  \item{line.slope}{
Vector with all sampling lines slope
}
  \item{line.dist1}{
Vector with the ROC Curves and sampling lines intersection points, the distance between this points and the reference point
}
}

\value{
This function return a list with:
\item{auctri}{Total area}
\item{areatri}{Vector with all triangles areas}
}

\seealso{
\code{\link{lineslope}}
\code{\link{linedistance}}
\code{\link{curvesegslope}}
\code{\link{curvesegsloperef}}
}

\keyword{ Triangle }
\keyword{ Area }
