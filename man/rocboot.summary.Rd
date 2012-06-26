\name{rocboot.summary}
\alias{rocboot.summary}

\title{
Summary of Comparation}
\description{
This function allows to see the information obtained throught function \code{roc.curve.boot}.
}
\usage{
rocboot.summary(result,mod1,mod2)
}
\arguments{
  \item{result}{
List of statistical measures obtaind throught \code{roc.curves.boot}
}
  \item{mod1}{
Name of the column of dataframe that represents the first modality
 }
   \item{mod2}{
Name of the column of dataframe that represents the second modality
 }
}


\seealso{
\code{\link{save.file.summary}}
}

\examples{
#data(zhang)
#moda1="mod1" 
#moda2="mod2"
#nameE="new_zhang"
#data=read.manually.introduced(zhang, moda1, TRUE, moda2, TRUE, "status", TRUE)
#results=roc.curves.boot(data,name=nameE,mod1= moda1,mod2=moda2) 
#rocboot.summary(results,moda1,moda2)
}

\keyword{ Summary }
\keyword{ ROC }
