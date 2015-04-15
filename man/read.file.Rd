\name{read.file}
\alias{read.file}
\title{Read data from file}

\description{
This function allows to read data from a file.
}
\usage{
read.file(name.file.csv, header.status = TRUE, separator = ";", decimal = ",", modality1,
testdirection1, modality2, testdirection2, status1, related = TRUE, status2 = NULL)
}

\arguments{
  \item{name.file.csv}{
Name of the file with data. The file must be in \code{csv} or \code{txt} format
}
  \item{header.status}{
Indicates if the file has a header row
}
  \item{separator}{
Indicates what is the column separator
}
  \item{decimal}{
Indicates what is the decimal separator
}
  \item{modality1}{
Name of the column of dataframe that represents the first modality
}
  \item{testdirection1}{
Indicates the direction of the test for modality 1. If \code{TRUE} means that larger test results represent more positive test
}
  \item{modality2}{
Name of the column of dataframe that represents the second modality
}
  \item{testdirection2}{
Indicates the direction of the test for modality 2. If \code{TRUE} means that larger test results represent more positive test
}
  \item{status1}{
Name of the column of dataframe that represents the Status 1
}
  \item{related}{
Boolean parameter that represents if the two modalities are related or not
}
  \item{status2}{
Name of the column of dataframe that represents the Status 2
}
}
\details{
The default column separator is \code{";"}. And the default decimal separator is \code{"."}. \code{header.status} has also a default value that is \code{TRUE}.
By default, the related parameter is set to \code{TRUE}. In this case the \code{status2} is not necessary (by default set to \code{(NULL)}, because in related modalities the status is the same. 
Otherwise, if related is set to \code{FALSE}, its necessary to indicate the name of \code{status2} column.
In the data must be listed first all values of the distribution of negative cases (0), followed by the positive ones (1).
}
\value{
This functions returns a list with the following data:
\item{sim1.ind}{Vector with the data for Curve 1}
\item{sim2.ind}{Vector with the data for Curve 2}
\item{sim1.sta}{Vector with the status for Curve 1}
\item{sim2.sta}{Vector with the status for Curve 2}
}

\seealso{
\code{\link{read.manually.introduced}}
}

\keyword{ Read }
\keyword{ Data }
\keyword{File }

