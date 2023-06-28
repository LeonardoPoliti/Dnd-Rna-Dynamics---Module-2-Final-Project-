## DNA/RNA DYNAMICS - MODULE 2 FINAL PROJECT 
## Leonardo Politi - 28/6/2023 
## Master’s Degree in Bioinformatics, University of Bologna


suppressMessages(library(minfi))

data_path = "~/data/Input"
clean_Manifest_path = "~/Illumina450Manifest_clean.RData"


# 1. Load raw data with minfi and create an object called RGset storing the RGChannelSet object:  
#list.files(data_path)
targets <- read.metharray.sheet(data_path)  
RGset <- read.metharray.exp(targets = targets)
dim(RGset)
#save(RGset,file="RGset.RData")



# 2. Create the dataframes Red and Green to store the red and green fluorescence respectively: 
Green <- data.frame(getGreen(RGset))
#dim(Green)
Red <- data.frame(getRed(RGset))
#dim(Red)



# 3a. Red and Green fluorescence of address: 10715421 
Green[rownames(Green)=="10715421",] #samples intensities 
Red[rownames(Red)=="10715421",]

# 3b. Check in the manifest file if the address corresponds to a Type I or a Type II probe and, in case of Type I probe, report its color:
load(clean_Manifest_path)
Illumina450Manifest_clean[Illumina450Manifest_clean$AddressA_ID=="10715421",] #retrieve probe info from address



# 4. Create the object MSet.raw 
MSet.raw <- preprocessRaw(RGset)
#head(MSet.raw)
#save(MSet.raw,file="MSet_raw.RData")



# 5a. QCplot:
qc <- getQC(MSet.raw) #Consider the median of meth. and unmeth. channels for each sample
plotQC(qc)

# 5b. Check the intensity of negative controls using minfi:
controlStripPlot(RGset, controls="NEGATIVE") #plot the intensity values of each type of controls probes in our samples

# 5c. p-value:
detP <- detectionP(RGset) 
#head(detP)
#save(detP,file="detP.RData")

failed <- detP>0.01 # consider a detectionP threshold of 0.01
summary(failed) # returns the number of failed (TRUE) and not failed (FALSE) positions for each sample



# 6. Calculate raw beta and M values and plot the densities of mean methylation values, dividing the samples in WT and MUT

b <- getBeta(MSet.raw)
#summary(b)
M <- getM(MSet.raw)
#summary(M)

M_MUT <- M[,c(targets$Group == 'MUT')]
M_WT <- M[,c(targets$Group == 'WT')]

b_MUT <- b[,c(targets$Group == 'MUT')]
b_WT <- b[,c(targets$Group == 'WT')]

mean_of_M_MUT <- apply(M_MUT,1,mean)
mean_of_M_MUT <- na.omit(mean_of_M_MUT) # remove NaN
d_mean_of_M_MUT <- density(mean_of_M_MUT)

mean_of_M_WT <- apply(M_WT,1,mean)
mean_of_M_WT <- na.omit(mean_of_M_WT)
d_mean_of_M_WT <- density(mean_of_M_WT)

mean_of_b_MUT <- apply(b_MUT,1,mean)
mean_of_b_MUT <- na.omit(mean_of_b_MUT)
d_mean_of_b_MUT <- density(mean_of_b_MUT)

mean_of_b_WT <- apply(b_WT,1,mean)
mean_of_b_WT <- na.omit(mean_of_b_WT)
d_mean_of_b_WT <- density(mean_of_b_WT)

par(mfrow=c(2,2))
plot(d_mean_of_M_MUT,main="Density of M Values (MUT)",col="lightblue")
plot(d_mean_of_M_WT,main="Density of M Values (WT)",col="blue")
plot(d_mean_of_b_MUT,main="Density of b Values (MUT)",col="orange")
plot(d_mean_of_b_WT,main="Density of b Values (WT)",col="red")



# 7. NORMALIZATION - preprocessNoob

dfI <- Illumina450Manifest_clean[Illumina450Manifest_clean$Infinium_Design_Type=="I",] # Infinium 1
dfI <- droplevels(dfI) # removes any unused levels
#dim(dfI)
dfII <- Illumina450Manifest_clean[Illumina450Manifest_clean$Infinium_Design_Type=="II",] # Infinium 2
dfII <- droplevels(dfII)
#dim(dfII)

b_I <- b[rownames(b) %in% dfI$IlmnID,] #subsets the 'b' object 
b_II <- b[rownames(b) %in% dfII$IlmnID,]

# densities of the means
mean_of_b_I <- apply(b_I,1,mean)
mean_of_b_I <- na.omit(mean_of_b_I) # remove NaN
d_mean_of_b_I <- density(mean_of_b_I) 

mean_of_b_II <- apply(b_II,1,mean)
mean_of_b_II <- na.omit(mean_of_b_II)
d_mean_of_b_II <- density(mean_of_b_II)

# densities of the standard deviations
sd_of_b_I <- apply(b_I,1,sd) 
sd_of_b_I <- na.omit(sd_of_b_I)
d_sd_of_b_I <- density(sd_of_b_I,)

sd_of_b_II <- apply(b_II,1,sd)
sd_of_b_II <- na.omit(sd_of_b_II)
d_sd_of_b_II <- density(sd_of_b_II)

preprocessNoob_results <- preprocessNoob(RGset) # requires ‘IlluminaHumanMethylation450kanno.ilmn12.hg19'
#preprocessNoob_results

b_preprocessNoob <- getBeta(preprocessNoob_results)
#head(b_preprocessNoob)
#save(b_preprocessNoob,file="b_preprocessNoob.RData")

b_preprocessNoob_I <- b_preprocessNoob[rownames(b_preprocessNoob) %in% dfI$IlmnID,]
b_preprocessNoob_II <- b_preprocessNoob[rownames(b_preprocessNoob) %in% dfII$IlmnID,]

mean_of_b_preprocessNoob_I <- apply(b_preprocessNoob_I,1,mean)
d_mean_of_b_preprocessNoob_I <- density(mean_of_b_preprocessNoob_I,na.rm=T) #another way to remove NaN

mean_of_b_preprocessNoob_II <- apply(b_preprocessNoob_II,1,mean)
d_mean_of_b_preprocessNoob_II <- density(mean_of_b_preprocessNoob_II,na.rm=T)

sd_of_b_preprocessNoob_I <- apply(b_preprocessNoob_I,1,sd)
d_sd_of_b_preprocessNoob_I <- density(sd_of_b_preprocessNoob_I,na.rm=T)

sd_of_b_preprocessNoob_II <- apply(b_preprocessNoob_II,1,sd)
d_sd_of_b_preprocessNoob_II <- density(sd_of_b_preprocessNoob_II,na.rm=T)

# Plots:
Group = factor(targets$Group,levels = unique(targets$Group)) 

par(mfrow=c(2,3))
palette(c("orange","lightblue"))
plot(d_mean_of_b_I,col="blue",main="raw beta",xlim=c(-0.1,1.1),ylim=c(0,6))
lines(d_mean_of_b_II,col="red")
plot(d_sd_of_b_I,col="blue",main="raw sd",xlim=c(0,0.6),ylim=c(0,50))
lines(d_sd_of_b_II,col="red")
boxplot(b,col = Group,ylim=c(0,1))
plot(d_mean_of_b_preprocessNoob_I,col="blue",main="preprocessNoob beta",xlim=c(-0.1,1.1),ylim=c(0,6))
lines(d_mean_of_b_preprocessNoob_II,col="red")
plot(d_sd_of_b_preprocessNoob_I,col="blue",main="preprocessNoob sd",xlim=c(0,0.6),ylim=c(0,50))
lines(d_sd_of_b_preprocessNoob_II,col="red")
boxplot(b_preprocessNoob,col = Group,ylim=c(0,1)) 



# 8. PCA : 
pca_results <- prcomp(t(b_preprocessNoob),scale=T)
summary(pca_results)
samples = rownames(pca_results$x)
sample_labels = substr(samples, nchar(samples) - 5, nchar(samples)) #keep only the last part of the name for a cleaner visualization


# 8a.  Groups
#levels(Group)
palette(c("blue","red"))
plot(pca_results$x[,1], pca_results$x[,2],cex=2,pch=2,col=Group,xlab="PC1",ylab="PC2",xlim=c(-700,700),ylim=c(-1000,1000))
text(pca_results$x[,1], pca_results$x[,2],labels=sample_labels,cex=0.5,pos=1)
legend("bottomright",legend=levels(Group),col=c(1:nlevels(Group)),pch=2)

#8b. Sex
Gender = factor(targets$Sex,levels = unique(targets$Sex))
#levels(Gender)
palette(c("blue","red"))
plot(pca_results$x[,1], pca_results$x[,2],cex=2,pch=2,col=Gender,xlab="PC1",ylab="PC2",xlim=c(-700,700),ylim=c(-1000,1000))
text(pca_results$x[,1], pca_results$x[,2],labels=sample_labels,cex=0.5,pos=1)
legend("bottomright",legend=levels(Gender),col=c(1:nlevels(Gender)),pch=2)

#8c. Batch
sub_Sentrix_ID = substr(targets$Basename, start = nchar(targets$Basename) -2, stop = nchar(targets$Basename))
Batch = factor(sub_Sentrix_ID,levels = unique(sub_Sentrix_ID))

#levels(Batch)
palette(c("blue","red"))
plot(pca_results$x[,1], pca_results$x[,2],cex=2,pch=2,col=Batch,xlab="PC1",ylab="PC2",xlim=c(-700,700),ylim=c(-1000,1000))
text(pca_results$x[,1], pca_results$x[,2],labels=sample_labels,cex=0.5,pos=1)
legend("bottomright",legend=levels(Batch),col=c(1:nlevels(Batch)),pch=2)



#9. identify differentially methylated probes between group WT and group MUT using T-TEST
ttest_pvalue <- function(x) {
  t_test <- t.test(x~ Group)
  return(t_test$p.value)
} 

pValues_ttest <- apply(b_preprocessNoob,1, ttest_pvalue)
#length(pValues_ttest)

final_ttest <- data.frame(b_preprocessNoob, pValues_ttest)

final_ttest_0.05 <- final_ttest[final_ttest$pValues_ttest<=0.05,] # consider a threshold of 5%
dim(final_ttest_0.05)[1]



#10. Apply multiple test correction and set a significant threshold of 0.05
corrected_pValues_BH <- p.adjust(final_ttest$pValues_ttest,"BH")
corrected_pValues_Bonf <- p.adjust(final_ttest$pValues_ttest,"bonferroni")
final_ttest_corrected <- data.frame(final_ttest, corrected_pValues_BH, corrected_pValues_Bonf)
#head(final_ttest_corrected)

final_corrected_BH_0.05 <- final_ttest_corrected[final_ttest_corrected$corrected_pValues_BH<=0.05,]
dim(final_corrected_BH_0.05)[1]

final_corrected_Bonf_0.05 <- final_ttest_corrected[final_ttest_corrected$corrected_pValues_Bonf<=0.05,]
dim(final_corrected_Bonf_0.05)[1]



# 11a. Volcano Plot:
b_preprocessNoob_groupMUT <- b_preprocessNoob[,Group=="MUT"]
mean_b_preprocessNoob_groupMUT <- apply(b_preprocessNoob_groupMUT,1,mean)
b_preprocessNoob_groupWT <- b_preprocessNoob[,Group=="WT"]
mean_b_preprocessNoob_groupWT <- apply(b_preprocessNoob_groupWT,1,mean)

delta <- mean_b_preprocessNoob_groupMUT-mean_b_preprocessNoob_groupWT

toVolcPlot <- data.frame(delta, -log10(final_ttest_corrected$pValues_ttest))

plot(toVolcPlot[,1], toVolcPlot[,2],pch=16,cex=0.5, xlab = "mean methylation differences", ylab = "-log10(P)")
Highlight_left <- toVolcPlot[toVolcPlot[,1]< -0.1 & toVolcPlot[,2]>(-log10(0.01)),]
Highlight_right <- toVolcPlot[toVolcPlot[,1]> 0.1 & toVolcPlot[,2]>(-log10(0.01)),]
points(Highlight_left[,1], Highlight_left[,2],pch=16,cex=0.7,col=("red")) 
points(Highlight_right[,1], Highlight_right[,2],pch=16,cex=0.7,col=("green"))


# 11b. Manhattan Plot:

library(qqman)

final_ttest <- data.frame(rownames(final_ttest),final_ttest)
colnames(final_ttest)[1] <- "IlmnID"
#dim(final_ttest)      

final_ttest_annotated <- merge(final_ttest, Illumina450Manifest_clean,by="IlmnID")
#dim(final_ttest_annotated)

input_Manhattan <- final_ttest_annotated[colnames(final_ttest_annotated) %in% c("IlmnID","CHR","MAPINFO","pValues_ttest")]
#dim(input_Manhattan)

# reorder the levels of the CHR
order_chr <- c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X","Y")
input_Manhattan$CHR <- factor(input_Manhattan$CHR,levels=order_chr )

input_Manhattan$CHR <- as.numeric(input_Manhattan$CHR) # convert factors to numbers

manhattan(input_Manhattan, snp="IlmnID",chr="CHR", bp="MAPINFO", p="pValues_ttest",annotatePval = 0.00001,col= rainbow(24))



# 12. Produce an heatmap of the top 100 differentially methylated probes:
suppressMessages(library(gplots))

final_ttest <- final_ttest[order(final_ttest$pValues_ttest),] #order final_ttest
input_heatmap=as.matrix(final_ttest[1:100,2:9]) # considers only the top 100 most significant CpG probes
#head(final_ttest)

Group
targets
colorbar <- c("green","red","green","green","red","red","green","red")
col2=colorRampPalette(c("green","black","red"))(100)

# average linkage
heatmap.2(input_heatmap,col=col2,Rowv=T,Colv=T, hclustfun = function(x) hclust(x,method = 'average'), dendrogram="both",key=T,ColSideColors=colorbar,density.info="none",trace="none",scale="none",symm=F, main="Average linkage")


