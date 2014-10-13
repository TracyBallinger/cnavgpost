setwd("/Users/tracyballinger/Documents/cn-avg")
dat=read.table("allevents.dat", header=TRUE)

fpdat=dat[dat$true==0,]
tpdat=dat[dat$true==1,]

tpcost=dat[dat$true==1 & dat$avecost>0,]
tpzero=dat[dat$true==1 & dat$avecost==0,]
fpcost=dat[dat$true==0 & dat$avecost>0,]
fpzero=dat[dat$true==0 & dat$avecost==0,]

tpcosth=hist(tpcost$Lscore, plot=FALSE)
fpcosth=hist(fpcost$Lscore, plot=FALSE)
fpzeroh=hist(fpzero$Lscore, plot=FALSE)

mylwd=2
myylim=c(0,1000)
plot(tpcosth$mids, tpcosth$counts, type="S", col="blue", lwd=mylwd, ylim=myylim, xlab="Likelihood score", ylab="Number of events")
lines(fpcosth$mids, fpcosth$counts, type="S", col="grey", lwd=mylwd, ylim=myylim)
lines(fpzeroh$mids, fpzeroh$counts, type="S", col="black", lwd=mylwd, ylim=myylim)
legend("topleft", legend=c(
paste("TP, cost>0 (", nrow(tpcost), ")", sep=""), 
paste("FP, cost>0 (", nrow(fpcost),")", sep=""), 
paste("FP, cost=0 (", nrow(fpzero), ")", sep="")), 
lwd=mylwd, col=c("blue", "grey", "black"), bty="n")

par(mfrow=c(1,1))
par(mar=c(5,5,1,2)+0.1)
fph=hist(fpdat$Lscore, plot=FALSE)
tph=hist(tpdat$Lscore, plot=FALSE)
plot(tph$mids, tph$counts, type="S", col="black", lwd=mylwd, ylim=myylim, xlab="Likelihood score", ylab="Number of events", cex.lab=1.3)
lines(fph$mids, fph$counts, type="S", col="grey", lwd=mylwd, ylim=myylim)
legend("topleft", legend=c(
paste("TP (", prettyNum(nrow(tpdat), big.mark=","), ")", sep=""), 
paste("FP (", prettyNum(nrow(fpdat), big.mark=","),")", sep="")), 
lwd=mylwd, col=c("black", "grey"), bty="n")


h=hist(fpcost$avecost, breaks=seq(0, 60, 0.5), xlab="Average cost of event", ylab="frequency", col="black", main="Histogram of events with nonzero cost")
hist(tpcost$avecost, breaks=h$breaks, add=TRUE, col="red")
legend("topright", legend=c(
paste("TP, cost>0 (", nrow(tpcost), ")", sep=""), 
paste("FP, cost>0 (", nrow(fpcost),")", sep="")), 
fill=c("red", "black"), bty="n")

h=hist(fpcost$Lscore, xlab="Likelihood score of event", ylab="frequency", col="black", main="Histogram of events with nonzero cost", ylim=myylim)
hist(tpcost$Lscore, add=TRUE, col="red", border="red")
legend("topright", legend=c(
paste("TP, cost>0 (", nrow(tpcost), ")", sep=""), 
paste("FP, cost>0 (", nrow(fpcost),")", sep="")), 
fill=c("red", "black"), bty="n")

text(h$breaks[1], myylim[2], labels=c(h$counts[1]), pos=4, col="white", cex=0.7)


################################################
# Look at the average cost of the different types of events 
predevents=dat$true>-1
amps=(dat$true>-1) & dat$event_type=="amp"
dels=(dat$true>-1) & dat$event_type=="del"
adjs=(dat$true>-1) & dat$event_type=="adj"

hamps=hist(dat$avecost[amps])
hdels=hist(dat$avecost[dels])
hadjs=hist(dat$avecost[adjs])

plot(hamps$mids, hamps$counts, type="l", xlab="average cost of event", ylab="number of events", xlim=c(0,20))
lines(hdels$mids, hdels$counts, type="l", col="red")
lines(hadjs$mids, hadjs$counts, type="l", col="blue")
legend("topright", legend=c(
paste("Amplification (", sum(amps), ")", sep=""),
paste("Deletion (", sum(dels), ")", sep=""),
paste("Inversion (", sum(adjs), ")", sep="")), 
col=c("black", "red", "blue"), bty="n", lwd=2)

plot(dat$CNval, dat$avecost, xlab="CN change", ylab="Average cost of event", pch=20, col=gray(.2, alpha=0.2))