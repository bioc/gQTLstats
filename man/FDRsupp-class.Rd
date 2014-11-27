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
    \item{\code{FDRmodel}:}{Object of class \code{"gam"} 
that was fit to elements of \code{tab}}
    \item{\code{sessinfo}:}{sessionInfo() value at time
of construction}
    \item{\code{theCall}:}{instance of class \code{"call"}
showing call leading to construction}
  }
}
\section{Methods}{
  \describe{
    \item{getFDRfunc}{\code{signature(x = "FDRsupp")}: extract
the FDR approximating function, a function of one (vector) argument
assumed to represent association scores, evaluating to the
plug-in FDR estimates corresponding to these scores  }
    \item{getTab}{\code{signature(x = "FDRsupp")}: extract the
table of association scores and empirical FDR estimates }
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
