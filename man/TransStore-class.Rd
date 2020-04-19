\name{TransStore-class}
\Rdversion{1.1}
\docType{class}
\alias{TransStore-class}
\alias{describe,TransStore-method}
\alias{describe}

\title{Class \code{"TransStore"}}
\description{
Manage collection of related trans-gQTL results in BatchJobs registries, typically one per chromosome
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("TransStore", ...)}.
%%  ~~ describe objects here ~~ 
}
\section{Slots}{
  \describe{
    \item{\code{allRegistries}:}{Object of class \code{"list"} containing
\code{\link{Registry}} instances }
    \item{\code{numSubmitted}:}{Object of class \code{"numeric"} records
number of jobs submitted for each registry }
    \item{\code{numDone}:}{Object of class \code{"numeric"} records
number of jobs completed for each registry }
    \item{\code{nloci}:}{Object of class \code{"numeric"} records number of
loci with test results for each registry }
    \item{\code{jobinfos}:}{Object of class \code{"list"} records results
of \code{\link{getJobInfo}} for each registry }
  }
}
\section{Methods}{
  \describe{
    \item{describe}{\code{signature(object = "TransStore")}: summarize
information about a store }
	 }
}
%\references{
%%  ~~put references to the literature/web site here~~
%}
%\author{
%%  ~~who you are~~
%}
%\note{
%%%  ~~further notes~~
%}
%
%%% ~Make other sections like Warning with \section{Warning }{....} ~
%
%\seealso{
%%%  ~~objects to See Also as \code{\link{~~fun~~}}, ~~~
%%%  ~~or \code{\linkS4class{CLASSNAME}} for links to other classes ~~~
%}
\examples{
showClass("TransStore")
}
\keyword{classes}
