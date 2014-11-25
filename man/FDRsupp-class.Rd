\name{FDRsupp-class}
\Rdversion{1.1}
\docType{class}
\alias{FDRsupp-class}
\alias{getFDRfunc}
\alias{getTab}
\alias{getFDRfunc,FDRsupp-method}
\alias{getTab,FDRsupp-method}

\title{Class \code{"FDRsupp"}}
\description{
Support for FDR computations with ciseStore instances
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("FDRsupp", ...)}.
}
\section{Slots}{
  \describe{
    \item{\code{tab}:}{Object of class \code{"data.frame"} a table
with association scores and plug-in FDR estimates evaluated on 
selected score values }
    \item{\code{FDRfunc}:}{Object of class \code{"function"} a
function of one argument with input association score and output
the corresponding FDR estimate }
  }
}
\section{Methods}{
  \describe{
    \item{getFDRfunc}{\code{signature(x = "FDRsupp")}: ... }
    \item{getTab}{\code{signature(x = "FDRsupp")}: ... }
	 }
}
%\references{
%%  ~~put references to the literature/web site here~~
%}
%\author{
%%  ~~who you are~~
%}
\note{
Typically the \code{FDRfunc} function is constructed using
a smooth model relating the estimated FDR to association scores.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%%%  ~~objects to See Also as \code{\link{~~fun~~}}, ~~~
%%%  ~~or \code{\linkS4class{CLASSNAME}} for links to other classes ~~~
%}
\examples{
showClass("FDRsupp")
}
\keyword{classes}
