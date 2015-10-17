
qqStore = function(st, ids=NULL, .probs = c(0,seq(.6,.8,.2),.9,.95,.99,.999,.9999,1),
   xlim.in=c(.2,75), lowfac = .5, xlab="Permutation distribution",
   ylab="Distribution of score statistic", countpos=50, plot.it=TRUE,
   doab=TRUE, scoreField="chisq", permField="permScore_1", ...) {
chisq = storeToFf(st, scoreField, ids=ids, ... )
ps1 = storeToFf(st, permField, ids = ids, ... )
qy = quantile(chisq, .probs)
qx = quantile(ps1, .probs)
myc = sapply(2:length(qy), function(x) 
       sum(chisq >= qy[x-1] & chisq < qy[x]))
myc[length(myc)] = myc[length(myc)] + length(chisq) - sum(myc) # highest vals
qx = qx[-length(qx)]
qy = qy[-length(qy)]
frac = paste0("(", round(100*myc/length(chisq),2), "%)")
qx[1] = lowfac*qx[2]
qy[1] = lowfac*qy[2]
if (plot.it) {
 plot(qx,qy,log="xy", pch="x", xlim=xlim.in,xlab="Permutation distribution",
     ylab="Distribution of score statistic")
 text(rep(countpos,length(qy)), qy, paste0(myc,frac))
 if (doab) abline(0,1, lty=2)
 }
ans = list(qx=qx, qy=qy, counts=myc, fracs=frac)
invisible(ans)
}

