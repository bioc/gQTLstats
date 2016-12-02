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

pifdr2 = function(obs,perms,trimToUnit=TRUE,expandPerms=TRUE,...) {
#
# y binned into intervals using hist
# for computing plug-in FDR
#
 if (!expandPerms) stopifnot((length(perms)%%length(obs)) == 0)
 lenperms = length(perms)
 lenobs = length(obs)
 approx_n_perms = lenperms/lenobs
 if (approx_n_perms > 1.5) message("possibly >1 permutations supplied, not yet handled, trimming")
 if (lenperms > lenobs) perms = sample(perms, size=lenobs, replace=FALSE)
 if (lenperms < lenobs) perms = sample(perms, size=lenobs, replace=TRUE)
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

pifdr3 = function(obs,perms,trimToUnit=TRUE,...) {
#
# y binned into intervals using hist
# for computing plug-in FDR
#
# previous implementation required length(perms) an integer multiple of length(obs)
# relaxing this to allow detailed filtering of permutation realizations
# in trans testing ... the numbers of observed and permuted 
# gene-snp pairs meeting a given
# trans condition (e.g., MAF and distance) for a given _rank_
# need not satisfy this condition
#
# stopifnot((length(perms)%%length(obs)) == 0)
 y = perms
 x = obs
 fac = length(y)/length(x)
 h = hist(y, breaks=c(-Inf,x,Inf), plot=FALSE)
 oy = ceiling(fac*length(x))-(cumsum(h$counts))  # cum. count up ordered binning of y from top
 rx = rank(x)
 brx = h$breaks[-c(1,length(x)+2)]
 ncalls = length(x):1
 fdr = oy[-length(oy)]/ceiling((fac*ncalls))
 # following can be used to verify current approach to handling bins
 #cbind(x=x, brx=brx, obrx=brx[rx], oy=oy[-length(oy)][rx], fdr=fdr[rx], ncalls=ncalls[rx])
 ans = fdr[rx]
 if (trimToUnit) ans = pmin(1, ans)
 ans
}
