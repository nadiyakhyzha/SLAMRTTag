# The following was used to generate Figure 2A. 

# Open up RNA-seq data
library(rtracklayer)

gtf<- import("~/Documents/AnalysisTools/hg19.gtf")
gtf_df<-as.data.frame(gtf)
gtf_df<-(gtf_df[gtf_df$type=="gene",])[,c(1:6,10,12,14)]

# Open up RNA-seq count tables
featureCounts_RNAseq_n1<- read.delim("~/Desktop/Lab Files/Sequencing/231128_VH00699_395_AAF5L7MM5/featureCounts/featureCounts_NK_hs_2100.exon.tabular", comment.char="#")[,c(1,7)]
featureCounts_RNAseq_n2<- read.delim("~/Desktop/Lab Files/Sequencing/231128_VH00699_395_AAF5L7MM5/featureCounts/featureCounts_NK_hs_2101.exon.tabular", comment.char="#")[,c(1,7)]
featureCounts_RNAseq_n3<- read.delim("~/Desktop/Lab Files/Sequencing/231128_VH00699_395_AAF5L7MM5/featureCounts/featureCounts_NK_hs_2102.exon.tabular", comment.char="#")[,c(1,7)]

#--------------------------------------------------------------------------------------------------------------------------------------
# Make a dataframe that contains all gene names and RNA-seq counts

RNAseq_merged<- cbind(gtf_df, featureCounts_RNAseq_n1[,2], featureCounts_RNAseq_n2[,2], featureCounts_RNAseq_n3[,2])

# Convert to CPM

RNAseq_merged$`featureCounts_RNAseq_n1[, 2]`<- (RNAseq_merged$`featureCounts_RNAseq_n1[, 2]`)/sum(RNAseq_merged$`featureCounts_RNAseq_n1[, 2]`)*1000000
RNAseq_merged$`featureCounts_RNAseq_n2[, 2]`<- (RNAseq_merged$`featureCounts_RNAseq_n2[, 2]`)/sum(RNAseq_merged$`featureCounts_RNAseq_n2[, 2]`)*1000000
RNAseq_merged$`featureCounts_RNAseq_n3[, 2]`<- (RNAseq_merged$`featureCounts_RNAseq_n3[, 2]`)/sum(RNAseq_merged$`featureCounts_RNAseq_n3[, 2]`)*1000000

# Find the average of the RNA-seq counts
average_CPM<- vector(length=nrow(RNAseq_merged))
RNAseq_merged <- cbind(RNAseq_merged, average_CPM)

for (i in 1:nrow(RNAseq_merged)){
  RNAseq_merged$average_CPM[i]<- mean(RNAseq_merged$`featureCounts_RNAseq_n1[, 2]`[i], RNAseq_merged$`featureCounts_RNAseq_n2[, 2]`[i])
}


#-----------------------------------------------------------------------------------------------------------------------------------
# Open up list of Polycomb-adjacent transcripts

K27me3_UP <- read.delim("~/Desktop/Lab Files/SLAM-RT&Tag/Figure 1/K562/K562_H3K27me3/K562_K27me3_IgG_UP_genebody.txt", header=FALSE)

# Assign RNA-seq counts for each Polycomb-adjacent transcript
for (i in 1:nrow(RNAseq_merged)){
K27me3_UP$RNAseq[i]<-RNAseq_merged$average_CPM[RNAseq_merged$gene_name==K27me3_UP$V5[i]]
}


#-----------------------------------------------------------------------------------------------------------------------------------
# Identify the top 50% expressed transcripts

# Reorder in descending order
top5<-RNAseq_merged[order(RNAseq_merged$average_CPM, decreasing=TRUE),]

# Remove everything that is not expressed (CPM = 0)
top5<-top5[top5$average_CPM!=0,]

# Isolate the top 50% expressed transcripts
top50<-top5[1:(0.5*nrow(top5)),]

#-----------------------------------------------------------------------------------------------------------------------------------
# Plotting RNA-seq counts as violin plots

library(ggplot2)
library(viridis)


RNAseq_UP<- K27me3_UP$RNAseq
RNAseq_top50<- top50$average_CPM

RNAseq_UP<- cbind(RNAseq_UP, "UP")
RNAseq_top50<- cbind(RNAseq_top50, "top50")

RNAseq<- as.data.frame(rbind(RNAseq_UP, RNAseq_top50))
colnames(RNAseq)<- c("CPM", "type")
RNAseq$type <- as.factor(RNAseq$type)
RNAseq$CPM <- as.numeric(RNAseq$CPM)

RNAseq$type <- factor(RNAseq$type, levels=c("UP", "top50"))

RNAseq_plot<- ggplot(RNAseq, aes(x=type, y=CPM, fill=type)) +
  geom_violin(width=1.4) +
  geom_boxplot(width=0.1, color="black") +
  scale_fill_viridis(discrete = TRUE) +
  theme_minimal()+ 
  ylim(0,200)+
  ylab("CPM") +
  xlab("")

RNAseq_plot<- RNAseq_plot +scale_fill_manual(values=c("lightblue", "grey"))

pdf("RNAseq_K27me3_UPvstop50.pdf", 4, 5)
RNAseq_plot
dev.off()

#-----------------------------------------------------------------------------------------------------------------------------------
# Performing statistical analysis

t.test(as.numeric(RNAseq_UP[,1]), as.numeric(RNAseq_top50[,1]), paired=FALSE, alternative="two.sided")
# Welch Two Sample t-test
# t = -19.055, df = 7710, p-value < 2.2e-16


median(as.numeric(RNAseq_UP[,1]))
# 3.257789
median(as.numeric(RNAseq_top50[,1]))
# 31.49196


#-----------------------------------------------------------------------------------------------------------------------------------
# Saving the coordinates of the top 50% expressed transcripts

# Function for getting gene body coordinates of top 50% expressed transcripts
library(GenomicRanges)

getFullLength <- function(x){
  x<-x
  x<- makeGRangesFromDataFrame(x, keep.extra.columns = TRUE)
  seqlevelsStyle(x) <- "UCSC"
  x<- as.data.frame(x)
  return(x)
}

# Apply to top 50% expressed transcripts

RNAseq_top50_genebody<- getFullLength(top50)

# The following function is to save the dataframe in bed format",
saveBed = function(x){
  #save file name
  file = deparse(substitute(x))
  file = paste(file,".txt", sep = "")
  #remove nonchromosomal postions
  x <- x[grep("chr", x$seqnames),]
  #rearrange columns in bed format, sort rows, and remove header
  x = x[ , c(1, 2, 3, 4, 9, 5)]
  x = x[order(x$seqnames, x$start),]
  names(x) = NULL
  
  #write to file
  write.table(x, file = file, row.names=FALSE, sep="\t", quote = FALSE)
}

# Save the bed file
saveBed(RNAseq_top50_genebody)
