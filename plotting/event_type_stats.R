setwd("/Users/tracyballinger/Documents/cn-avg")

dat=read.table("allevents.dat", header=TRUE)
dat=read.table("alledges.dat", header=TRUE)

dat=dat[dat$avecost>0,]
eventype=dat$event_type
lscore=dat$Lscore
cnval=dat$CNval
truth=dat$true
length=dat$length
prevals=dat$prevals

###############################################
# Look at correlation between predicted and true orders of events.  Need to correct for the number of events per simulation.  Add this data in.  
x=matrix(as.numeric(unlist(strsplit(as.character(dat$prevals), ","))), ncol=length(prevals), nrow=3)
trueprevals=x[1,]
predprevals=x[2,]
prevalsd=x[3,]
x=matrix(as.numeric(unlist(strsplit(as.character(dat$orders), ","))), ncol=length(prevals), nrow=3)
trueorder=x[1,]
predorder=x[2,]
ordersd=x[3,]
myx=truth==1
x=cor.test(trueorder[myx], predorder[myx])
plot(trueorder[myx], predorder[myx], xlab="True Order", ylab="Predicted Order", cex=.5)
mtext(paste("correlation:", round(cor(trueorder[myx], predorder[myx]), digits=3), "p-value:", x$p.value))
# get the maximum order for each event 
simids=read.table("sim_ids_for_allevents.txt", skip=1)[,1]
unisimids=rle(as.character(simids))$values
maxtrueorders=matrix(nrow=length(simids), ncol=1)
maxpredorders=matrix(nrow=length(simids), ncol=1)
for (simid in unisimids){
	maxord=max(trueorder[simids==simid])
	maxtrueorders[simids==simid]=maxord
	maxord=max(predorder[simids==simid])
	maxpredorders[simids==simid]=maxord
}
truecorr=trueorder/maxtrueorders
predcorr=predorder/maxpredorders
x=cor.test(truecorr[myx], predcorr[myx])
plot(truecorr[myx], predcorr[myx], xlab="True Order", ylab="Predicted Order", cex=.5)
mtext(paste("correlation:", round(cor(truecorr[myx], predcorr[myx]), digits=3), "p-value:", signif(x$p.value, digits=3)))


##############################################################
# Look at accuracy of different types of events. 
types=c("amp","del", "adj")
recall=c(1:6)
precision=c(1:6)
total=c(1:6)
accuracy=c(1:6)
f1score=c(1:6)

for (i in c(1:3)){
	type=types[i]
	TP=sum(eventype==type & truth==1)
	FP=sum(eventype==type & truth==0)
	TN=sum(eventype==type & truth==2)
	FN=sum(eventype==type & truth==-1)
	recall[i]=TP/(TP+FN)
	precision[i]=TP/(TP+FP)
	total[i]=TP+FP+TN+FN
	accuracy[i]=(TP+TN)/total[i]
	f1score[i]=2*TP/(2*TP+FP+FN)
}
	
pvalue=0.5
for (i in c(4:6)){
	type=types[i-3]
	TP=sum(eventype==type & truth==1 & lscore> pvalue)
	FP=sum(eventype==type & truth==0 & lscore>pvalue)
	TN=sum(eventype==type & (truth==2 | (truth==0 & lscore<=pvalue)))
	FN=sum(eventype==type & (truth==-1 | (truth==1 & lscore <=pvalue)))
	recall[i]=TP/(TP+FN)
	precision[i]=TP/(TP+FP)
	total[i]=TP+FP+TN+FN
	accuracy[i]=(TP+TN)/total[i]
	f1score[i]=2*TP/(2*TP+FP+FN)
}


#################
# Barplot of accuracy for the various types of events. 
myy=rbind(recall[4:6], precision[4:6], f1score[4:6])
barplot(myy, beside=TRUE, names.arg=types, legend.text=c("recall", "precision", "F1 score"), args.legend=c(bty="n", x="topleft"))

mids=seq(0, max(length), by=100000)
blocklen=100000
x1=seq(0,10*blocklen, by=blocklen)
x2=seq(11*blocklen, 100*blocklen, by=10*blocklen)
x3=c(101*blocklen, max(length))
mids=c(x1,x2, x3)
recall=matrix(nrow=3, ncol=length(mids))
precision=matrix(nrow=3, ncol=length(mids))
accuracy=matrix(nrow=3, ncol=length(mids))
total=matrix(nrow=3, ncol=length(mids))
f1score=matrix(nrow=3, ncol=length(mids))

TP=matrix(nrow=3, ncol=length(mids))
FP=matrix(nrow=3, ncol=length(mids))
TN=matrix(nrow=3, ncol=length(mids))
FN=matrix(nrow=3, ncol=length(mids))

pvalue=0.5
for (i in c(1:length(mids))){
	l=mids[i]
	x=(length>(l-500) & length <= (l+500))
	for (t in c(1:3)){
		type=types[t]
		TP[t,i]=sum(eventype[x]==type & truth[x]==1 & lscore[x]> pvalue)
		FP[t,i]=sum(eventype[x]==type & truth[x]==0 & lscore[x]>pvalue)
		TN[t,i]=sum(eventype[x]==type & (truth[x]==2 | (truth[x]==0 & lscore[x]<=pvalue)))
		FN[t,i]=sum(eventype[x]==type & (truth[x]==-1 | (truth[x]==1 & lscore[x] <=pvalue)))
	}
}

recall=TP/(TP+FN)
precision=TP/(TP+FP)
total = TP + FP + TN + FN
accuracy=(TP+TN)/total
f1score=2*TP/(2*TP+FP+FN)
totaltrue=TP+FN

par(mfrow=c(5,1))
par(oma=c(4,4,2,2)+0.1)
nbin=20
nbin=length(mids)
margins=c(1,4,0,2)+0.1
par(mar=margins)
xvalues=seq(0, nbin)+0.5
mylwd=2
matplot(xvalues[1:nbin], t(totaltrue[,1:nbin]), type="s", xlab="size of event", ylab="Number true", xaxt="n", lwd=mylwd)
legend("topright", legend=types[1:3], col=c(1:3), lty=c(1:3), lwd=mylwd, bty="n")
par(mar=margins)
matplot(xvalues[1:nbin], t(total[,1:nbin]), type="s", lwd=mylwd, xlab="size of event", ylab="Number of Events", xaxt="n")
legend("topright", legend=types[1:3], col=c(1:3), lty=c(1:3), lwd=mylwd, bty="n")
par(mar=margins)
matplot(xvalues[1:nbin], t(recall[,1:nbin]), type="s", lwd=mylwd, xlab="size of event", ylab="Recall", xaxt="n")
legend("topright", legend=types[1:3], col=c(1:3), lty=c(1:3),lwd=mylwd, bty="n")
par(mar=margins)
matplot(xvalues[1:nbin], t(precision[,1:nbin]), type="s", xlab="size of event", ylab="Precision", lwd=mylwd, xaxt="n")
legend("topright", legend=types[1:3], col=c(1:3), lty=c(1:3),lwd=mylwd, bty="n")
par(mar=margins)
matplot(xvalues[1:nbin], t(f1score[,1:nbin]), type="s", lwd=mylwd, xlab="size of event", ylab="f1score", xaxt="n")
legend("topright", legend=types[1:3], col=c(1:3), lty=c(1:3), lwd=mylwd, bty="n")
axis(side=1, at=xvalues[1:nbin], labels=mids/blocklen, outer=TRUE)
mtext(side=1, text="Number of genomic blocks covered by event", line=4, cex=0.7)

# combine all the types
TP=apply(TP, 2, sum)
FP=apply(FP, 2, sum)
TN=apply(TN, 2, sum)
FN=apply(FN, 2, sum)
recall=TP/(TP+FN)
precision=TP/(TP+FP)
total = TP + FP + TN + FN
accuracy=(TP+TN)/total
f1score=2*TP/(2*TP+FP+FN)
totaltrue=TP+FN

par(mfrow=c(5,1))
par(oma=c(4,4,2,2)+0.1)
nbin=20
margins=c(1,4,0,2)+0.1
par(mar=margins)
plot(mids[1:nbin]-50000, totaltrue[1:nbin], type="s", xlab="size of event", ylab="Number true", xaxt="n")
par(mar=margins)
plot(mids[1:nbin]-50000,total[1:nbin], type="s", xlab="size of event", ylab="Number of Events", xaxt="n")
par(mar=margins)
plot(mids[1:nbin]-50000, recall[1:nbin], type="s", xlab="size of event", ylab="Recall", xaxt="n")
par(mar=margins)
plot(mids[1:nbin]-50000, precision[1:nbin], type="s", xlab="size of event", ylab="Precision", xaxt="n")
par(mar=margins)
plot(mids[1:nbin]-50000, f1score[1:nbin], type="s", xlab="size of event", ylab="f1score")

