\name{diffareatriangles}
\alias{diffareatriangles}
\title{Difference Between Area Triangles}

\description{
This function allows to calculate the difference between triangles areas formed by the same sampling lines in two different ROC curves. It also allows to calculate the difference between total areas.
}
\usage{
diffareatriangles(area.triangle1, area.triangle2)
}

\arguments{
  \item{area.triangle1}{
Vector with all triangles areas of the Curve 1
}
  \item{area.triangle2}{
Vector with all triangles areas of the Curve 2
}
}

\value{
This function returns a list with:
\item{diffareas}{Difference between each triangle area}
\item{diffauc}{Difference between total areas}
}

\seealso{
\code{\link{areatriangles}}
}

\keyword{ Area }
\keyword{ Triangles }
\keyword{ Difference }

