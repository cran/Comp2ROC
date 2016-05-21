\name{comp.roc.delong}
\alias{comp.roc.delong}
\title{Calculate areas and stats}

\description{
This function allows to calculate the areas under the curve for each curve and some statistical measures.
}
\usage{
comp.roc.delong(sim1.ind, sim1.sta, sim2.ind, sim2.sta, related = TRUE)
}

\arguments{
  \item{sim1.ind}{
Vector with the data for Curve 1
}
  \item{sim1.sta}{
Vector with the status for Curve 1
}
  \item{sim2.ind}{
Vector with the data for Curve 2
}
  \item{sim2.sta}{
Vector with the status for Curve 2
}
  \item{related}{
Boolean parameter that represents if the two modalities are related or not
}
}
\details{
This function calculates the Wilcoxon Mann Whitney matrix for each modality, areas, standard deviations, variances and global correlations.
}
\value{
This function returns a list with:
\item{Z }{Hanley Z calculation}
\item{pvalue}{p-value for this Z}
\item{AUC}{Area under curve for each modality}
\item{SE}{Standard error}
\item{S}{Variance for each modality}
\item{R}{Correlation Coeficient}
}

\examples{
\donttest{
data(zhang)
modality1DataColumn = "modality1"
modality2DataColumn = "modality2"
data = read.manually.introduced(zhang, modality1DataColumn, TRUE,
                                modality2DataColumn, TRUE, "status", TRUE)
sim1.ind = unlist(data[1])
sim2.ind = unlist(data[2])  
sim1.sta = unlist(data[3])
sim2.sta = unlist(data[4])
comp.roc.delong(sim1.ind, sim1.sta, sim2.ind, sim2.sta)
}
}

\keyword{ DeLong }
\keyword{ Hanley }
\keyword{ Wilcoxon }
\keyword{ Whitney }
\keyword{ Mann }
\keyword{ AUC }
\keyword{ Areas }
\keyword{ Z-stats }
