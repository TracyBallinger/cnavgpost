setwd("/Users/tracyballinger/Documents/cn-avg")
#dat=read.table("simacc_filt50pc.dat", header=TRUE)
#dat=read.table("simint_events_pval50.dat", header=TRUE)
#edgedat=read.table("simint_edges_nofilt.dat", header=TRUE)
dat=read.table("intsims_events_pval50.txt", header=TRUE)
edgedat=read.table("intsims_edges_pval00.txt", header=TRUE)

totaledges=edgedat$TP+edgedat$FP+edgedat$TN+edgedat$FN
trueedges=edgedat$TP+edgedat$FN+edgedat$TN
connectivity=trueedges/edgedat$truenodes

TP=dat$TP
FP=dat$FP
TN=dat$TN
FN=dat$FN

total=TP+FP+TN+FN
totaltrue=dat$TP+dat$FN
precision=TP/(FP+TP)
recall=TP/totaltrue
f1score=2*(precision*recall)/(precision+recall)
f1score[is.nan(f1score)]=0
	
colors=c(rgb(27,158,119, maxColorValue=255), 
rgb(217,95,2,maxColorValue=255),
rgb(117,112,179, maxColorValue=255))

ptcolors=dat$blocks
ptcolors[dat$blocks==100]=colors[1]
ptcolors[dat$blocks==150]=colors[2]
ptcolors[dat$blocks==200]=colors[3]

shapes=c(15,16,17)
ptshapes=dat$blocks
ptshapes[dat$blocks==100]=shapes[1]
ptshapes[dat$blocks==150]=shapes[2]
ptshapes[dat$blocks==200]=shapes[3]

#################----------------
# Function for getting the pvalues out of a lm
lmp <- function (modelobject) {
    if (class(modelobject) != "lm") stop("Not an object of class 'lm' ")
    f <- summary(modelobject)$fstatistic
    p <- pf(f[1],f[2],f[3],lower.tail=F)
    attributes(p) <- NULL
    return(p)
}


par(mfrow=c(3,1))
margins=c(1,4,0,2)+0.1
par(oma=c(4,4,2,2)+0.1)
par(mar=margins)

legendtext=c(100,150,200)
x=100
legendtext[1]=paste(x, " (", min(totaltrue[dat$blocks==x]), "-", max(totaltrue[dat$blocks==x]), ")", sep="")
x=150
legendtext[2]=paste(x, " (", min(totaltrue[dat$blocks==x]), "-", max(totaltrue[dat$blocks==x]), ")", sep="")
x=200
legendtext[3]=paste(x, " (", min(totaltrue[dat$blocks==x]), "-", max(totaltrue[dat$blocks==x]), ")", sep="")

#betterpreds=dat$truescores>dat$minscores
mylabsize=1.5
y=recall
plot(connectivity, y, col=ptcolors, pch=ptshapes, xaxt="n", ylab="recall",cex.lab=mylabsize)
points(connectivity[betterpreds], y[betterpreds], pch=".")
mylm=lm(y~connectivity)
abline(mylm$coefficients)
mtext(paste("R-squared:",signif(summary(mylm)$r.squared, digits=3), ", p-value: ", signif(lmp(mylm), digits=3)), outer=FALSE, side=3, cex=0.7) 
legend("topright", legend=legendtext, title="number of blocks (#events)", pch=shapes, col=colors, bty="n")

y=precision
plot(connectivity, y, col=ptcolors, pch=ptshapes, xaxt="n", ylab="precision",cex.lab=mylabsize)
points(connectivity[betterpreds], y[betterpreds], pch=".")
mylm=lm(y~connectivity)
abline(mylm$coefficients)
mtext(paste("R-squared:",signif(summary(mylm)$r.squared, digits=3), ", p-value: ", signif(lmp(mylm), digits=3)), outer=FALSE, side=3, cex=0.7) 
legend("topright", legend=legendtext, title="number of blocks (#events)", pch=shapes, col=colors, bty="n")

y=f1score
plot(connectivity, y, col=ptcolors, pch=ptshapes, xlab="connectivity", ylab="F1 score", cex.lab=mylabsize)
points(connectivity[betterpreds], y[betterpreds], pch=".")
mylm=lm(y~connectivity)
abline(mylm$coefficients)
mtext(paste("R-squared:",signif(summary(mylm)$r.squared, digits=3), ", p-value: ", signif(lmp(mylm), digits=3)), outer=FALSE, side=3, cex=0.7) 
mtext("connectivity", outer=TRUE, side=1, line=2, cex.lab=mylabsize)
legend("topright", legend=legendtext, title="number of blocks (#events)", pch=shapes, col=colors, bty="n")




par(mfrow=c(1,1))
par(mar=c(5,4,4,2)+0.1)
par(xpd=FALSE)
plot(dat$truescores-dat$minscores, f1score, xlab="True Cost - Predicted Cost", pch=20)
z=line(dat$truescores-dat$minscores, f1score)
abline(coef(z))
z=cor.test(dat$truescores-dat$minscores, f1score)
mtext(paste("correlation: ", round(z$estimate, digits=4), " p-value: ", round(z$p.value, digits=4), sep=""), outer=TRUE, cex=0.8)