setwd("/Users/tracyballinger/Documents/cn-avg/lotadir")

plot_stats_over_time = function(TP, FP, Tots, total, legendtitle="", legendtext=c(1:4), linewidth=1, plot_title=""){
precision=TP/(FP+TP)
recall=TP/total
f1score=2*(precision*recall)/(precision+recall)
par(mfrow=c(3,1))
margins=c(1,4,0,8)+0.1
par(oma=c(4,4,2,2)+0.1)
par(xpd=TRUE)
par(mar=margins)
colors=c(rgb(217,240,163, maxColorValue=255),
rgb(217,240,163, maxColorValue=255),
rgb(217,240,163, maxColorValue=255),
rgb(173,221,142, maxColorValue=255),
rgb(120,198,121, maxColorValue=255),
rgb(65,171,93, maxColorValue=255),
rgb(35,132,67, maxColorValue=255),
rgb(0,104,55, maxColorValue=255),
rgb(0,69,41,maxColorValue=255),
rgb(0,69,41,maxColorValue=255))

nlines=ncol(TP)
mycolors=colors[(length(colors)-nlines):length(colors)]
legendcmd = function(){
	legend("topright", inset=c(-0.2,0),legend=legendtext, title=legendtitle, lty=c(1:length(legendtext)), lwd=linewidth, col=mycolors, bty="n")
	}

matplot(TP, type="l", lwd=linewidth, lty=c(1:nlines), col=mycolors, xlab="iteration", ylab="TP Count", xaxt="n")
title(outer=TRUE, main=plot_title)
legendcmd()
matplot(FP, type="l", lwd=linewidth, lty=c(1:nlines), col=mycolors, xlab="iteration", ylab="FP count", xaxt="n")
legendcmd()
#matplot(recall, type="l", lwd=linewidth, lty=c(1:nlines), col=mycolors, xlab="iteration", ylab="sensitivity", xaxt="n")
#legendcmd()
#matplot(precision, type="l", lwd=linewidth, lty=c(1:nlines), col=mycolors, xlab="iteration", ylab="specificity", xaxt="n")
#legendcmd()
matplot(f1score, type="l", lwd=linewidth, lty=c(1:nlines), col=mycolors, xlab="iteration", ylab="f1score")
legendcmd()
title(outer=TRUE, main=plot_title)

}

plot_counts_over_time = function(TP, FP, Tots, total, legendtitle="", legendtext=c(1:4), linewidth=1, plot_title=""){
par(mfrow=c(3,1))
margins=c(1,4,0,2)+0.1
par(oma=c(4,4,2,2)+0.1)
par(mar=margins)
legendcmd = function(){
legend("topleft", legend=legendtext, title=legendtitle, lty=c(1:length(legendtext)), lwd=linewidth, col=c(1:length(legendtext)), bty="n")}
matplot(TP, type="l", lwd=linewidth, xlab="iteration", ylab="TP Count", xaxt="n")
title(outer=TRUE, main=plot_title)
legendcmd()
matplot(FP, type="l", lwd=linewidth, xlab="iteration", ylab="FP count", xaxt="n")
legendcmd()
#matplot(Tots, type="l", lwd=linewidth, xlab="iteration", ylab="Total Positive")
#legendcmd()
}


###################################################################
##
pdf("accuracies_with_diff_k.pdf", width=7, height=8)

TP=read.table("counts_TP_ks.sum.dat")
FP=read.table("counts_FP_ks.sum.dat")
Tots=read.table("counts_Tot_ks.sum.dat")
total=1212
plot_stats_over_time(TP, FP, Tots, total, legendtitle="k value", legendtext=c(0,0.01, 0.1, 1), linewidth=1.5, plot_title="Accuracy at various values of k")
#plot_counts_over_time(TP, FP, Tots, total, legendtitle="k value", legendtext=c(0,0.01, 0.1, 1), linewidth=1.5, plot_title="Counts at various values of k")

TP=read.table("zero_TP_ks.sum.dat")
FP=read.table("zero_FP_ks.sum.dat")
Tots=read.table("zero_Tot_ks.sum.dat")
total=1212
plot_stats_over_time(TP, FP, Tots, total, legendtitle="k value", legendtext=c(0,0.01, 0.1, 1), linewidth=1.5, plot_title="Accuracy at various values of k with cost=0")
#plot_counts_over_time(TP, FP, Tots, total, legendtitle="k value", legendtext=c(0,0.01, 0.1, 1), linewidth=1.5, plot_title="Counts at various values of k with cost=0")

TP=read.table("nonzero_TP_ks.sum.dat")
FP=read.table("nonzero_FP_ks.sum.dat")
Tots=read.table("nonzero_Tot_ks.sum.dat")
total=1212
plot_stats_over_time(TP, FP, Tots, total, legendtitle="k value", legendtext=c(0,0.01, 0.1, 1), linewidth=1.5, plot_title="Accuracy at various values of k with cost>0")
#plot_counts_over_time(TP, FP, Tots, total, legendtitle="k value", legendtext=c(0,0.01, 0.1, 1), linewidth=1.5, plot_title="Counts at various values of k with cost>0")


linewidth=1
Tots=read.table("counts_Tot_ks.sum.dat")
zerotots=read.table("zero_Tot_ks.sum.dat")
nonzerotots=read.table("nonzero_Tot_ks.sum.dat")
matplot(Tots, type="l", lwd=linewidth, xlab="iteration", ylab="Total Events", col=1)
matplot(nonzerotots, type="l", lwd=linewidth, add=TRUE, col=2)
#matplot(zerotots, type="l", lwd=linewidth, add=TRUE, col=3)

legendtext=c("total events", "with cost >0")
legend("topleft", legend=legendtext, lwd=linewidth, col=c(1:length(legendtext)), bty="n", lty=1)

dev.off()

#################################_______________________________
pdf("accuracies_with_diff_n.pdf", width=7, height=8)

TP=read.table("counts_TP_ns.sum.dat")
FP=read.table("counts_FP_ns.sum.dat")
Tots=read.table("counts_Tot_ns.sum.dat")
total=1212
plot_stats_over_time(TP, FP, Tots, total, legendtitle="# runs", legendtext=c(1:ncol(TP)), linewidth=1.5, plot_title="Accuracy at various numbers of runs")
#plot_counts_over_time(TP, FP, Tots, total, legendtitle="k value", legendtext=c(0,0.01, 0.1, 1), linewidth=1.5, plot_title="Counts at various values of k")

TP=read.table("zero_TP_ns.sum.dat")
FP=read.table("zero_FP_ns.sum.dat")
Tots=read.table("zero_Tot_ns.sum.dat")
total=1212
plot_stats_over_time(TP, FP, Tots, total, legendtitle="# runs", legendtext=c(1:ncol(TP)), linewidth=1.5, plot_title="Accuracy at various numbers of runs with cost=0")
#plot_counts_over_time(TP, FP, Tots, total, legendtitle="k value", legendtext=c(0,0.01, 0.1, 1), linewidth=1.5, plot_title="Counts at various values of k with cost=0")

TP=read.table("nonzero_TP_ns.sum.dat")
FP=read.table("nonzero_FP_ns.sum.dat")
Tots=read.table("nonzero_Tot_ns.sum.dat")
total=1212
plot_stats_over_time(TP, FP, Tots, total, legendtitle="# runs", legendtext=c(1:ncol(TP)), linewidth=1.5, plot_title="Accuracy at various number of runs with cost>0")
#plot_counts_over_time(TP, FP, Tots, total, legendtitle="k value", legendtext=c(0,0.01, 0.1, 1), linewidth=1.5, plot_title="Counts at various values of k with cost>0")


linewidth=1
Tots=read.table("counts_Tot_ns.sum.dat")
zerotots=read.table("zero_Tot_ns.sum.dat")
nonzerotots=read.table("nonzero_Tot_ns.sum.dat")
matplot(Tots, type="l", lwd=linewidth, xlab="iteration", ylab="Total Events", col=1)
matplot(nonzerotots, type="l", lwd=linewidth, add=TRUE, col=2)
#matplot(zerotots, type="l", lwd=linewidth, add=TRUE, col=3)

legendtext=c("total events", "with cost >0")
legend("topleft", legend=legendtext, lwd=linewidth, col=c(1:length(legendtext)), bty="n", lty=1)

dev.off()

############################################
pdf("accuracies_with_diff_n_test.pdf", width=7, height=12)

TP=read.table("counts_TP_ns.sum.dat")
FP=read.table("counts_FP_ns.sum.dat")
Tots=read.table("counts_Tot_ns.sum.dat")
total=1212
plot_stats_over_time(TP, FP, Tots, total, legendtitle="# runs", legendtext=c(1:9), linewidth=1.5, plot_title="Accuracy at various numbers of runs")
#plot_counts_over_time(TP, FP, Tots, total, legendtitle="k value", legendtext=c(0,0.01, 0.1, 1), linewidth=1.5, plot_title="Counts at various values of k")

TP=read.table("zero_TP_ns.sum.dat")
FP=read.table("zero_FP_ns.sum.dat")
Tots=read.table("zero_Tot_ns.sum.dat")
total=1212
plot_stats_over_time(TP, FP, Tots, total, legendtitle="# runs", legendtext=c(1:9), linewidth=1.5, plot_title="Accuracy at various numbers of runs with cost=0")
#plot_counts_over_time(TP, FP, Tots, total, legendtitle="k value", legendtext=c(0,0.01, 0.1, 1), linewidth=1.5, plot_title="Counts at various values of k with cost=0")

TP=read.table("nonzero_TP_ns.sum.dat")
FP=read.table("nonzero_FP_ns.sum.dat")
Tots=read.table("nonzero_Tot_ns.sum.dat")
total=1212
plot_stats_over_time(TP, FP, Tots, total, legendtitle="# runs", legendtext=c(1:9), linewidth=1.5, plot_title="Accuracy at various number of runs with cost>0")
#plot_counts_over_time(TP, FP, Tots, total, legendtitle="k value", legendtext=c(0,0.01, 0.1, 1), linewidth=1.5, plot_title="Counts at various values of k with cost>0")


linewidth=1
Tots=read.table("counts_Tot_ns.dat")
zerotots=read.table("zero_Tot_ns.dat")
nonzerotots=read.table("nonzero_Tot_ns.dat")
matplot(Tots, type="l", lwd=linewidth, xlab="iteration", ylab="Total Events", col=1)
matplot(nonzerotots, type="l", lwd=linewidth, add=TRUE, col=2)
#matplot(zerotots, type="l", lwd=linewidth, add=TRUE, col=3)

legendtext=c("total events", "with cost >0")
legend("topleft", legend=legendtext, lwd=linewidth, col=c(1:length(legendtext)), bty="n", lty=1)

dev.off()

