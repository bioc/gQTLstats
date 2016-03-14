\name{pifdr}
\alias{pifdr}
\title{
utility for computing plug-in FDR
}
\description{
utility for computing plug-in FDR
}
\usage{
pifdr( obs, perms, legacy=FALSE, trimToUnit = TRUE, ... ) 
}
\arguments{
  \item{obs}{
observed association scores
}
  \item{perms}{
vector of association scores under permutation; length should
be integer multiple of \code{length(obs)}
 }
\item{legacy}{logical, if TRUE, use the approximate version of
pifdr() available before 12/30/2013, with additional arguments if
desired}
\item{trimToUnit}{logical, if TRUE, values greater than 1 are
replaced by 1.  Such values can occur, for example, with relatively
small sample sizes.}
\item{\dots}{extra arguments passed if legacy is TRUE}
}
\details{
Revised 12/30/13 to employ hist() to rapidly bin the
permuted values.  

Use \code{legacy=TRUE} to obtain the
approximate implementation, for which the following remarks held:
``As currently implemented the algorithm is quadratic in \code{length(obs)}.
While it is possible to get a unique FDR value for every element of \code{obs},
an approximate approach yields practically identical precision and by default this
will be used for \code{obs} with length 2000 or more.  In this case,
\code{\link{approx}} is used with \code{rule=2} to interpolate from the
grid-based FDR estimates back to the data values.''
Additional parameters \code{npts} and \code{applier} may be supplied
if \code{legacy} is set to TRUE.

\code{npts} defined the number of points spanning the range
of \code{obs} to be used for a lossy grid-based computation
only used if \code{length(obs)>npts}.

\code{applier} is to be an \code{sapply}-like function.
}
\value{
vector of plug-in FDR estimates congruent to \code{obs}
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
Hastie Tibshirani and Friedman Elements of Statistical
Learning ch 18.7
%% ~put references to the literature/web site here ~
}
%\author{
%%  ~~who you are~~
%}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
\examples{
set.seed(1234)
op = par(no.readonly=TRUE)
par(mfrow=c(2,2))
X = c(rchisq(30000,1),rchisq(300,10))
Y = rchisq(30300*3,1)
qqplot(Y, X, xlab="null", ylab="observed")
hist(pp <- pifdr(X,Y), xlab="plug-in FDR", main=" ")
library(multtest)
rawp = 1-pchisq(X, 1)
MT <- mt.rawp2adjp(rawp) 
MT2 = MT[[1]][order(MT[[2]]),]
plot(MT2[,"BH"], pp, xlab="BH FDR", ylab="plug-in FDR")
par(op)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models }
