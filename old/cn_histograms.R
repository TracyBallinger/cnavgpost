setwd("/Users/tracyballinger/Documents/cn-avg")

samplename="TCGA-06-0185"
directory="TCGA-06-0185_tracks"
require(lattice)

pdf(paste(samplename, "_histograms.pdf", sep=""), width=7, height=4)

#######################################################
## PLOT CGH Data from TCGA
# do this before reading in data: 
# cd HMS__HG-CGH-244A/Level_3
# cat hms.harvard.edu*.txt | grep -v barcode > all_copy_number_analysis.txt
# cd HMS__HG-CGH-244A/Level_2
### combine all the files into one
# cut -f2 file2 | paste file1 - > combined_sets.tsv 

cgh_level2_histogram_plots=function(data, dataname){
numsamples=ncol(data)-1
x = cor(data[,2:ncol(data)], use="na.or.complete")
correlationplot=levelplot(x, main=paste("Correlation matrix for", dataname), ylab="", xlab="", scales=list(x=list(rot=90)))
print(correlationplot)
nona=!(is.na(data[,2:ncol(data)]))
minval=min(data[,2:ncol(data)][nona])
maxval=max(data[,2:ncol(data)][nona])
numbins=20
mids=(seq(round(minval, digits=2) *100, round(maxval,digits=2)*100, length.out=numbins))/100
d=(mids[2]-mids[1])/2
mybreaks=c(mids-d, mids[numbins]+d)
histmat=matrix(nrow=numbins, ncol=numsamples)
for (i in 1:numsamples){
	h=hist(data[,(i+1)], breaks=mybreaks, plot=FALSE)
	histmat[,i]=h$counts
}
matplot(h$mids, histmat, type="l", xlab="Normalized Log2 ratio", ylab="Frequency", main=paste(dataname, "Level 2"), xlim=c(-5,5))
}

cgh_level3_histogram_plots= function(data, dataname){
minval=min(data[,5], na.rm=TRUE)
maxval=max(data[,5], na.rm=TRUE)
numbins=20
mids=(seq(round(minval, digits=2) *100, round(maxval,digits=2)*100, length.out=numbins))/100
d=(mids[2]-mids[1])/2
mybreaks=c(mids-d, mids[numbins]+d)
sampnames=unique(data[,1])
numsamples=length(sampnames)
histmat3=matrix(nrow=numbins, ncol=numsamples)
histbp=matrix(nrow=numbins, ncol=numsamples)
for (i in 1:numsamples){
	n = sampnames[i]
	subdata=data[data[,1]==n,]
	h=hist(subdata[,5], breaks=mybreaks, plot=FALSE)
	for (j in 1:numbins){
		keeps=subdata[,5]>=mybreaks[j] & subdata[,5]<mybreaks[j+1] & !is.na(subdata[,5])
		
		x = sum(as.numeric(subdata[keeps,4]-subdata[keeps,3]+sum(keeps)))
		histbp[j,i]=x
	}
	histmat3[,i]=h$counts
}
matplot(h$mids, histmat3, type="l", xlab="Normalized Log2 ratio", ylab="Frequency (segments)", main=paste(dataname, "Level 3"), xlim=c(-5,5))
legend("topright", legend=sampnames, lty=c(1:length(sampnames)), col=c(1:length(sampnames)), bty="n", cex=.7)
matplot(h$mids, histbp, type="l", xlab="Normalized Log2 ratio", ylab="Frequency (bp)", main=paste(dataname, "Level 3"), xlim=c(-5,5))
legend("topright", legend=sampnames, lty=c(1:length(sampnames)), col=c(1:length(sampnames)), bty="n", cex=.7)
}


hmslevel2=read.table(paste(directory, "/HMS__HG-CGH-244A/Level_2/combined_sets.tsv", sep=""), skip=1, header=TRUE, sep="\t")
cgh_level2_histogram_plots(data=hmslevel2, dataname="HMS__HG-CGH-244A")
hmslevel3=read.table(paste(directory, "/HMS__HG-CGH-244A/Level_3/all_copy_number_analysis.txt", sep=""))
cgh_level3_histogram_plots(data=hmslevel3, dataname="HMS__HG-CGH-244A")

msklevel2=read.table(paste(directory, "/MSKCC__HG-CGH-244A/Level_2/combined_sets.tsv", sep=""), skip=1, header=TRUE, sep="\t")
cgh_level2_histogram_plots(data=msklevel2, dataname="MSKCC__HG-CGH-244A")
msklevel3=read.table(paste(directory, "/MSKCC__HG-CGH-244A/Level_3/all_copy_number_analysis.txt", sep=""))
cgh_level3_histogram_plots(data=msklevel3[,c(1:4,6)], dataname="MSKCC__HG-CGH-244A")


######## Look at sequence-based CNV data  ##############################################
# Make histograms of cbs and cactus cnv
data=read.table(paste(directory, "/cactus.cnv.bg", sep=""))
cactuscnv=data[,3:4]
cactuscnv[,1]=data[,3]-data[,2]+1
data=read.table(paste(directory, "/cbs.cov.bg", sep=""))
cbscov=data[,3:4]
cbscov[,1]=data[,3]-data[,2]+1
data=read.table(paste(directory, "/majority.cnv.bg", sep=""))
majcnv=data[,3:4]
majcnv[,1]=data[,3]-data[,2]+1
data=read.table(paste(directory, "/minority.cnv.bg", sep=""))
mincnv=data[,3:4]
mincnv[,1]=data[,3]-data[,2]+1

minval=min(cactuscnv[,2], cbscov[,2], majcnv[,2], mincnv[,2])
maxval=max(cactuscnv[,2], cbscov[,2], majcnv[,2], mincnv[,2])
numbins=40
mids=(seq(round(minval, digits=2) *100, round(maxval,digits=2)*100, length.out=numbins))/100
d=(mids[2]-mids[1])/2
mybreaks=c(mids-d, mids[numbins]+d)
cnvhistmat=matrix(nrow=numbins, ncol=4)
cnvbpmat=matrix(nrow=numbins, ncol=4)
i=1
for (data in list(cactuscnv, cbscov, majcnv, mincnv)){
	h=hist(data[,2], breaks=mybreaks, plot=FALSE)
	cnvhistmat[,i]=h$counts
	for (j in 1:numbins){
		x = sum(data[data[,2]>=mybreaks[j] & data[,2]<mybreaks[j+1],1])
		cnvbpmat[j,i]=x
		}
	i=i+1
}

matplot(mids, cnvhistmat, type="l", xlab="Relative copy number", ylab="Frequency (by segment)", lwd=c(3,3,1,1))
legend("topright", c("cactus cnv", "CBS cnv", "majority cnv", "minority cnv"), lwd=c(3,3,1,1), col=c("black", "red", "green", "blue"), bty="n")

matplot(mids, cnvbpmat, type="l", xlab="Relative copy number", ylab="Frequency (by bp)", lwd=c(3,3,1,1))
legend("topright", c("cactus cnv", "CBS cnv", "majority cnv", "minority cnv"), lwd=c(3,3,1,1), col=c("black", "red", "green", "blue"), bty="n")


#########################################################
# LOOK at raw bambam input compared to CBS
# for the raw bambam input do this: 
# awk '{printf "%1.3f\n", $4}' bambam.cov.bg | sort | uniq -c > tmp
data=read.table(paste(directory, "/bambam.covhist.txt", sep=""))
bambamcnv=data
plot(bambamcnv[,2], bambamcnv[,1], type="l", main="bambam raw cnv", xlab="relative coverage", ylab="frequency", xlim=c(0,4))
 
data=read.table(paste(directory, "/cbs.cov.bg", sep=""))
cbscov=data[,3:4]
cbscov[,1]=data[,3]-data[,2]+1
minval=min(cbscov[,2])
maxval=max(cbscov[,2])

# Plot with lower resolution
numbins=40
mids=(seq(round(minval, digits=2) *100, round(maxval,digits=2)*100, length.out=numbins))/100
d=(mids[2]-mids[1])/2
mybreaks=c(mids-d, mids[numbins]+d)

h=hist(cbscov[,2], breaks=mybreaks, plot=FALSE)
## Bin the bambam raw data to compare to CBS segments
bambambinned=seq(1:numbins)
cbscovbp=bambambinned
for (j in 1:numbins){
	data=bambamcnv
	x=sum(data[data[,2]>=mybreaks[j] & data[,2]<mybreaks[j+1],1])
	bambambinned[j]=x
	data=cbscov
	cbscovbp[j]=sum(data[data[,2]>=mybreaks[j] & data[,2]<mybreaks[j+1],1])
}

bambamdens=bambamcnv[,1]/sum(bambamcnv[,1])

plot(bambamcnv[,2], bambamcnv[,1]/max(bambamcnv[,1]), type="l", main="bambam raw cnv", xlab="relative coverage (tumor/normal)", ylab="frequency (~bp)", xlim=c(0,4), col="grey")
lines(h$mids,  bambambinned/(max(bambambinned)), type="l", lwd=2)
#lines(h$mids, bambambinned/(max(bambambinned)), col="grey", lwd=2)
lines(h$mids, cbscovbp/max(cbscovbp), col="red", lwd=2)
#lines(h$mids, h$counts/(max(h$counts)), col="orange", lwd=2)
legend("topright", legend=c("raw bambam cnv", "CBS cnv"), lwd=c(1,2), col=c("black", "red"), bty="n")

################################################
#### NEW BAMBAM calls 
newbb = read.table(paste(directory,"/", samplename, "_whole_cnv.txt", sep=""))
newbbhist=seq(1:numbins)
for (j in 1:numbins){
	data=cbind(newbb[,3]-newbb[,2]+1, newbb[,4])
	x=sum(data[data[,2]>=mybreaks[j] & data[,2]<mybreaks[j+1],1])
	newbbhist[j]=x 
}
plot(bambamcnv[,2], bambamcnv[,1]/max(bambamcnv[,1]), type="l", main="bambam raw cnv", xlab="relative coverage (tumor/normal)", ylab="frequency (~bp)", xlim=c(0,4), col="grey")
lines(h$mids,  bambambinned/(max(bambambinned)), type="l", lwd=2)
lines(h$mids, cbscovbp/max(cbscovbp), col="red", lwd=2)
lines(h$mids, newbbhist/max(newbbhist), col="blue", lwd=2)
legend("topright", legend=c("raw bambam cnv", "CBS cnv", "new Bambam CNV"), lwd=c(2,2, 2), col=c("black", "red", "blue"), bty="n")

#plot with higher resolution
numbins=100
mids=(seq(round(minval, digits=2) *100, round(maxval,digits=2)*100, length.out=numbins))/100
d=(mids[2]-mids[1])/2
mybreaks=c(mids-d, mids[numbins]+d)
h=hist(cbscov[,2], breaks=mybreaks, plot=FALSE)
## Bin the bambam raw data to compare to CBS segments
bambambinned=seq(1:numbins)
cbscovbp=bambambinned
for (j in 1:numbins){
	data=bambamcnv
	x=sum(data[data[,2]>=mybreaks[j] & data[,2]<mybreaks[j+1],1])
	bambambinned[j]=x
	data=cbscov
	cbscovbp[j]=sum(data[data[,2]>=mybreaks[j] & data[,2]<mybreaks[j+1],1])
}
bambamdens=bambamcnv[,1]/sum(bambamcnv[,1])
newbbhist=seq(1:numbins)
for (j in 1:numbins){
	data=cbind(newbb[,3]-newbb[,2]+1, newbb[,4])
	x=sum(data[data[,2]>=mybreaks[j] & data[,2]<mybreaks[j+1],1])
	newbbhist[j]=x 
}
plot(bambamcnv[,2], bambamcnv[,1]/max(bambamcnv[,1]), type="l", main="bambam raw cnv", xlab="relative coverage (tumor/normal)", ylab="frequency (~bp)", xlim=c(0,4), col="grey")
lines(h$mids,  bambambinned/(max(bambambinned)), type="l", lwd=2)
lines(h$mids, cbscovbp/max(cbscovbp), col="red", lwd=2)
lines(h$mids, newbbhist/max(newbbhist), col="blue", lwd=2)
legend("topright", legend=c("raw bambam cnv", "CBS cnv", "new Bambam CNV"), lwd=c(2,2, 2), col=c("black", "red", "blue"), bty="n")


maxx=10
plot(newbb[,4], newbb[,7], main="new bambam cnv", xlab="Average relative coverage", ylab="Average allele fraction")

hist(newbb[newbb[,4]<maxx,4], main="new bambam cnv", xlab="Average relative coverage", ylab="Frequency")
 

dev.off()