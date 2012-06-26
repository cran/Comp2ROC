\name{roc.curves.plot}
\alias{roc.curves.plot}
\title{
Plot ROC curves
}
\description{
This function allows to plot the two roc curves in comparasion.
}
\usage{
roc.curves.plot(sim1.curve, sim2.curve, mod1, mod2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{sim1.curve}{
Curve 1 created using the function \code{performance}.
}
  \item{sim2.curve}{
Curve 2 created using the function \code{performance}.
}
  \item{mod1}{
Name of Modality 1
}
  \item{mod2}{
Name of Modality 2
}
}

\seealso{
\code{\link{read.file}}
\code{\link{read.manually.introduced}}

}


\examples{
#data(zhang)
#moda1="mod1" 
#moda2="mod2"
#data=read.manually.introduced(zhang, moda1, TRUE, moda2, TRUE, "status", TRUE)
#sim1.ind=unlist(data[1])
#sim2.ind=unlist(data[2])  
#sim1.sta=unlist(data[3])
#sim2.sta=unlist(data[4])
#sim1.pred <- prediction(sim1.ind, sim1.sta)
#sim2.pred <- prediction(sim2.ind, sim2.sta)
#sim1.curve <- performance(sim1.pred,"tpr","fpr")
#sim2.curve <- performance(sim2.pred,"tpr","fpr")
#roc.curves.plot(sim1.curve, sim2.curve,mod1=moda1,mod2=moda2)
}

\keyword{ Plot }
\keyword{ ROC }
