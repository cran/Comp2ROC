\name{comp.roc.curves}
\alias{comp.roc.curves}
\title{Calculate distribution}

\description{
This funtion calculates by bootstrapping the real distribution for the entire length set.
}
\usage{
comp.roc.curves(result, ci.flag = FALSE, graph.flag = FALSE, nome)
}

\arguments{
  \item{result}{
List of statistical measures obtaind throught \code{rocsampling} 
}
  \item{ci.flag}{
Flag that indicates if the user wants to calculate the confidance intervals
}
  \item{graph.flag}{
Flag that indicates if the user wants to draw the graph
}
  \item{nome}{
Name to put on the graph
}
}
\details{
In this function \code{ci.flag} and \code{graph.flag} are set \code{FALSE} by defaut
}
\value{
\item{boot}{statistics test}
\item{p-value }{p-value for one-sided}
\item{p-value2 }{p-value for two-sided}
\item{ci}{confidance interval}
}

\seealso{
\code{\link{rocsampling}}
}

\keyword{ bootstrapping }
\keyword{ comparation }
\keyword{ distribution }
