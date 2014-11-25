pifdr = function(obs,perms,trimToUnit=TRUE,...) {
#
# y binned into intervals using hist
# for computing plug-in FDR
#
 stopifnot((length(perms)%%length(obs)) == 0)
 y = perms
 x = obs
 fac = length(y)/length(x)
 h = hist(y, breaks=c(-Inf,x,Inf), plot=FALSE)
 oy = fac*length(x)-(cumsum(h$counts))  # cum. count up ordered binning of y from top
 rx = rank(x)
 brx = h$breaks[-c(1,length(x)+2)]
 ncalls = length(x):1
 fdr = oy[-length(oy)]/(fac*ncalls)
 # following can be used to verify current approach to handling bins
 #cbind(x=x, brx=brx, obrx=brx[rx], oy=oy[-length(oy)][rx], fdr=fdr[rx], ncalls=ncalls[rx])
 ans = fdr[rx]
 if (trimToUnit) ans = pmin(1, ans)
 ans
}

