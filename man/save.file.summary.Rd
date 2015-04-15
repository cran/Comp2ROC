\name{save.file.summary}
\alias{save.file.summary}
\title{Save File}

\description{
This functions allow to save the information on a file.
}
\usage{
save.file.summary(result, name, app = TRUE, mod1, mod2)
}

\arguments{
  \item{result}{
List of statistical measures obtaind throught roc.curves.boot
}
  \item{name}{
File name
}
  \item{app}{
Indicates if the user wants to append information on the same file
}
  \item{mod1}{
Name of the column of dataframe that represents the first modality
}
  \item{mod2}{
Name of the column of dataframe that represents the second modality
}
}
\details{
The user don't need to fill the \code{app} parameter, because by default it was set to \code{TRUE}. This parameter allow the user to choose if he wants the results of differents performances in the same file, or each time that he starts a new performance the file will be new.
}

\value{
This functions saves on the file with name \code{name} the performance parameters of the test.
}

\examples{
#data(zhang)
#moda1="mod1" 
#moda2="mod2"
#nameE="new_Zhang"
#data=read.manually.introduced(zhang, moda1, TRUE, moda2, TRUE, "status", TRUE)
#results=roc.curves.boot(data,name=nameE,mod1= moda1,mod2=moda2) 
#rocboot.summary(results,moda1,moda2)

#If the user wants to append the results
#save.file.summary(results,nameE,mod1=moda1,mod2=moda2)

#If the user doesn't want to append the results
#save.file.summary(results,nameE,app=FALSE,moda1,moda2)
}

\keyword{ Save }
\keyword{ File }
\keyword{ Rocboot }
