# This code was used to generate Supplementary Figure 1B, C, and D

# Load the appropriate packages
library(rtracklayer)
library(DESeq2)
library(ggplot2)
library(ggrepel)
library(VennDiagram)

#-----------------------------------------------------------------------------------------------------------------------------------
# Open up gtf file and extract gene level information
gtf<- import("~/PATH/hg19.gtf")
gtf_df<-as.data.frame(gtf)
gtf_df<-(gtf_df[gtf_df$type=="gene",])[,c(1:6,10,12,14)]

#-----------------------------------------------------------------------------------------------------------------------------------
# Open up Subread count tables. Here perform analysis on read counts that fall within the whole length of the gene
featureCounts_IgG_n1<- read.delim("~/PATH/featureCounts_NK_hs_1876.gene.tabular", comment.char="#")[,c(1,7)]
featureCounts_H3K27me3_n1<- read.delim("~/PATH/featureCounts_NK_hs_1877.gene.tabular", comment.char="#")[,c(1,7)]
featureCounts_SC35_n1<- read.delim("~/PATH/featureCounts_NK_hs_1879.gene.tabular", comment.char="#")[,c(1,7)]


featureCounts_IgG_n2<- read.delim("~/PATH/featureCounts_NK_hs_1886.gene.tabular", comment.char="#")[,c(1,7)]
featureCounts_H3K27me3_n2<- read.delim("~/PATH/featureCounts_NK_hs_1887.gene.tabular", comment.char="#")[,c(1,7)]
featureCounts_SC35_n2<- read.delim("~/PATH/featureCounts_NK_hs_1888.gene.tabular", comment.char="#")[,c(1,7)]


featureCounts_IgG_n3<- read.delim("~/PATH/featureCounts_NK_hs_1942.gene.tabular", comment.char="#")[,c(1,7)]
featureCounts_H3K27me3_n3<- read.delim("~/PATH/featureCounts_NK_hs_1943.gene.tabular", comment.char="#")[,c(1,7)]
featureCounts_SC35_n3<- read.delim("~/PATH/featureCounts_NK_hs_1944.gene.tabular", comment.char="#")[,c(1,7)]

#-----------------------------------------------------------------------------------------------------------------------------------
# Set up the metadata for DESeq analysis

Counts<- cbind(featureCounts_IgG_n1[,2], featureCounts_H3K27me3_n1[,2], featureCounts_SC35_n1[,2],
               featureCounts_IgG_n2[,2], featureCounts_H3K27me3_n2[,2], featureCounts_SC35_n2[,2],
               featureCounts_IgG_n3[,2], featureCounts_H3K27me3_n3[,2], featureCounts_SC35_n3[,2])


colnames(Counts)<-c("IgG n1", "H3K27me3 n1", "SC35 n1", 
                    "IgG n2", "H3K27me3 n2", "SC35 n2", 
                    "IgG n3", "H3K27me3 n3", "SC35 n3")

rownames(Counts)<-featureCounts_IgG_n1[,1]
Protein<- c("IgG", "H3K27me3", "SC35", "IgG", "H3K27me3", "SC35", "IgG", "H3K27me3", "SC35")

Condition<- c("n1", "n1", "n1", "n2","n2", "n2", "n3", "n3", "n3")
Metadata<- cbind(Protein, Condition)
rownames(Metadata)<- c("IgG n1", "H3K27me3 n1", "SC35 n1", 
                       "IgG n2", "H3K27me3 n2", "SC35 n2", 
                       "IgG n3", "H3K27me3 n3", "SC35 n3")

#-----------------------------------------------------------------------------------------------------------------------------------
# Perform DESeq normalization

Count_Table <- DESeqDataSetFromMatrix(countData = Counts, colData = Metadata,
                                      design = ~ Protein + Condition)
Count_Table_DESeq<- DESeq(Count_Table)

Counts_normalized<- cbind(gtf_df,counts(Count_Table_DESeq, normalized=TRUE))

#-----------------------------------------------------------------------------------------------------------------------------------
# Perform PCA analysis and save PCA graph. Used in Figure 1C.
# Retrieve information regarding % variance for each principal component
PC_value<-plotPCA(rlog(Count_Table_DESeq), intgroup= "Protein")

PC_value$labels$y
# "PC2: 32% variance"

PC_value$labels$x
# "PC1: 48% variance"

# Generate and save the PCA plot
PCA<- plotPCA(rlog(Count_Table_DESeq), intgroup= "Protein", returnData=TRUE)
PCA$name<- c("n1", "n1", "n1", "n2", "n2", "n2", "n3", "n3","n3")
colnames(PCA)<- c("PC1", "PC2", "group","Target", "Replicate")

PCA_plot<-ggplot(PCA, aes(PC1, PC2, color=Target, shape=Replicate)) +
  geom_point(size=5) +
  scale_shape_manual(values=c(16, 17, 15, 18))+
  scale_color_manual(values=c('lightblue', 'black', "salmon"))+
  xlab(paste0("PC1: 48% variance")) +
  ylab(paste0("PC2: 32% variance")) + 
  coord_fixed() +
theme_classic()

pdf("K562_IgGvsH3K27me3vsSC35_PCAplot.pdf", 6,6)
PCA_plot
dev.off()


#-----------------------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------------------
# Set up the metadata for DESeq analysis for H3K27me3 RT&Tag. This is performed on read counts that fall within whole gene body
# due to the overwhelming amount of intronic signal within H3K27me3 RT&Tag libraries. Will be used for generating volcano plot
# in Supplementary Figure 1C. 

# Open up Subread count tables. Here perform analysis on read counts that fall within the whole genes

featureCounts_IgG_n1<- read.delim("~/PATH/featureCounts_NK_hs_1480.gene.tabular", comment.char="#")[,c(1,7)]
featureCounts_H3K27me3_n1<- read.delim("~/PATH/featureCounts_NK_hs_1481.gene.tabular", comment.char="#")[,c(1,7)]

featureCounts_IgG_n2<- read.delim("~/PATH/featureCounts_NK_hs_1509.gene.tabular", comment.char="#")[,c(1,7)]
featureCounts_H3K27me3_n2<- read.delim("~/PATH/featureCounts_NK_hs_1510.gene.tabular", comment.char="#")[,c(1,7)]

featureCounts_IgG_n3<- read.delim("~/PATH/featureCounts_NK_hs_1971.gene.tabular", comment.char="#")[,c(1,7)]
featureCounts_H3K27me3_n3<- read.delim("~/PATH/featureCounts_NK_hs_1972.gene.tabular", comment.char="#")[,c(1,7)]



Counts_H3K27me3<- cbind(featureCounts_IgG_n1[,2], featureCounts_H3K27me3_n1[,2],
                        featureCounts_IgG_n2[,2], featureCounts_H3K27me3_n2[,2],
                        featureCounts_IgG_n3[,2], featureCounts_H3K27me3_n3[,2])


colnames(Counts_H3K27me3)<-c("IgG n1", "H3K27me3 n1", 
                             "IgG n2", "H3K27me3 n2", 
                             "IgG n3", "H3K27me3 n3")

rownames(Counts)<-featureCounts_IgG_n1[,1]
Protein<- c("IgG", "H3K27me3", "IgG", "H3K27me3", "IgG", "H3K27me3")

Condition<- c("n1", "n1", "n2", "n2", "n3", "n3")
Metadata<- cbind(Protein, Condition)
rownames(Metadata)<- c("IgG n1", "H3K27me3 n1", 
                       "IgG n2", "H3K27me3 n2", 
                       "IgG n3", "H3K27me3 n3")

#-----------------------------------------------------------------------------------------------------------------------------------
# Perform DESeq normalization
Count_Table_H3K27me3 <- DESeqDataSetFromMatrix(countData = Counts_H3K27me3, colData = Metadata,
                                      design = ~ Protein + Condition)
Count_Table_DESeq_H3K27me3<- DESeq(Count_Table_H3K27me3)

Counts_normalized_H3K27me3<- cbind(gtf_df,counts(Count_Table_DESeq_H3K27me3, normalized=TRUE))

#-----------------------------------------------------------------------------------------------------------------------------------
# Function that performs DESeq differential expression and then assigns whether a transcript is
# enriched or depleted in H3K27me3 RT&Tag libraries relative to IgG RT&Tag libraries

RetrieveAll <- function(n1, n2, n3, n4) {
  n1 <- n1
  n2 <- n2
  n3 <- n3
  n4 <- n4
  
  Results<-results(n4, contrast= c("Protein", n1, n2), name = n3)
  Results_log2<-Results$log2FoldChange
  Results_pvalue<-Results$padj
  diffexpressed <- "NO"
  Results_combined<- cbind(gtf_df,Results_log2, Results_pvalue, diffexpressed)
  Results_combined<- Results_combined[!is.na(Results_log2),]
  Results_combined<-as.data.frame(Results_combined)
  Results_combined$diffexpressed[Results_combined$Results_log2 >1 &
                                   Results_combined$Results_pvalue <0.05] <- "UP"
  Results_combined$diffexpressed[Results_combined$Results_log2 < -1 &
                                   Results_combined$Results_pvalue <0.05] <- "DOWN"
  return(Results_combined)
}

#-----------------------------------------------------------------------------------------------------------------------------------
# Run the differential enrichment function
Results_H3K27me3_IgG_ALL<- RetrieveAll("H3K27me3", "IgG", "Protein_H3K27me3_vs_IgG", Count_Table_DESeq_H3K27me3)

# Save this file which will become Table 1
write.csv(Results_H3K27me3_IgG_ALL, "K562_H3K27me3_differentialenrichment_table.csv")

#-----------------------------------------------------------------------------------------------------------------------------------
# Plot the volcano plot
mycolors <- c("black", "lightblue", "grey70")
names(mycolors) <- c("DOWN", "UP", "NO")

# Label the top 4 H3K27me3 enriched transcripts
Results_H3K27me3_IgG_ALL$delabel <- NA
Results_H3K27me3_IgG_ALL$delabel[Results_H3K27me3_IgG_ALL$gene_name == "XIST"] <- "XIST"
Results_H3K27me3_IgG_ALL$delabel[Results_H3K27me3_IgG_ALL$gene_name == "WWOX"] <- "WWOX"
Results_H3K27me3_IgG_ALL$delabel[Results_H3K27me3_IgG_ALL$gene_name == "RP11-648K4.2"] <- "RP11-648K4.2"
Results_H3K27me3_IgG_ALL$delabel[Results_H3K27me3_IgG_ALL$gene_name == "IMMP2L"] <- "IMMP2L"


p<- ggplot(data=Results_H3K27me3_IgG_ALL, aes(x=Results_log2, y=-log10(Results_pvalue), col=diffexpressed,  label=delabel)) + geom_point(size=1) + theme_minimal() + geom_text (size=5, nudge_x= 1, nudge_y = 1) + 
  geom_vline(xintercept=c(-1, 1), col="black") +
  geom_hline(yintercept=-log10(0.05), col="black") +
  theme(text = element_text(size = 20)) +
  labs(x="Fold change (log2) H3K27me3 over IgG", y="FDR (-log10)") +
  theme_classic()
p2 <- p + scale_colour_manual(values = mycolors)


pdf("VolcanoPlot_K562_H3K27me3vsIgG.pdf", 7,5)
p2
dev.off()

#-----------------------------------------------------------------------------------------------------------------------------------
# Count the number of differentially enriched transcipts
H3K27me3_UP<- Results_H3K27me3_IgG_ALL[Results_H3K27me3_IgG_ALL$diffexpressed=="UP",]
nrow(H3K27me3_UP)
# 1390

#-----------------------------------------------------------------------------------------------------------------------------------
# Save the coordinates of H3K27me3 enriched transcripts 

# Function for saving the TES coordinates
get500TES <- function(x){
  x<- x
  colnames(x)<- c("chr", "start", "end", "width",  "strand", "source", "gene_id", "genetype", "name", "log2", "pvalue", "difexp")
  for (i in 1:nrow(x)){
    if (x$strand[i] == "+"){
      x$strand2[i] <- "-"
    } else {
      x$strand2[i] <- "+"
    }
  }
  x$strand<- x$strand2
  x<- makeGRangesFromDataFrame(x, keep.extra.columns = TRUE)
  x<- promoters(x, upstream=1, downstream=500)
  seqlevelsStyle(x) <- "UCSC"
  
  x<- as.data.frame(x)
  
  for (i in 1:nrow(x)){
    if (x$strand[i] == "+"){
      x$strand2[i] <- "-"
    } else {
      x$strand2[i] <- "+"
    }
  }
  x$strand<- x$strand2
  x<- x[,1:12]
  return(x)
}



# Function for saving the full gene length of transcripts
getFullLength <- function(x){
  x<-x
  x<- makeGRangesFromDataFrame(x, keep.extra.columns = TRUE)
  seqlevelsStyle(x) <- "UCSC"
  x<- as.data.frame(x)
  return(x)
}


# Function for saving the TSS of transcripts
getTSS <- function(x){
  x<- x
  x<- makeGRangesFromDataFrame(x, keep.extra.columns = TRUE)
  x<- promoters(x, upstream=1, downstream=1)
  seqlevelsStyle(x) <- "UCSC"
  x<- as.data.frame(x)
  return(x)
}


K562_H3K27me3_UP_500TES<- get500TES(H3K27me3_UP)

K562_H3K27me3_UP_genebody<- getFullLength(H3K27me3_UP)

K562_H3K27me3_UP_TSS<- getTSS(H3K27me3_UP)


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


saveBed(K562_H3K27me3_UP_500TES)

saveBed(K562_H3K27me3_UP_genebody)

saveBed(K562_H3K27me3_UP_TSS)



#-----------------------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------------------
# Set up the metadata for DESeq analysis for SC35 RT&Tag. This is performed on read counts that fall within exons as most signal is exonic
# Will be used for generating volcano plot in Supplementary Figure 1C. 

# Open up Subread count tables. Here perform analysis on read counts that fall within the exons of genes
featureCounts_IgG_n1<- read.delim("~/PATH/featureCounts_NK_hs_2088.exon.tabular", comment.char="#")[,c(1,7)]
featureCounts_SC35_n1<- read.delim("~/PATH/featureCounts_NK_hs_2089.exon.tabular", comment.char="#")[,c(1,7)]

featureCounts_IgG_n2<- read.delim("~/PATH/featureCounts_NK_hs_1918.exon.tabular", comment.char="#")[,c(1,7)]
featureCounts_SC35_n2<- read.delim("~/PATH/featureCounts_NK_hs_1920.exon.tabular", comment.char="#")[,c(1,7)]

featureCounts_IgG_n3<- read.delim("~/PATH/featureCounts_NK_hs_2084.exon.tabular", comment.char="#")[,c(1,7)]
featureCounts_SC35_n3<- read.delim("~/PATH/featureCounts_NK_hs_2085.exon.tabular", comment.char="#")[,c(1,7)]

#-----------------------------------------------------------------------------------------------------------------------------------
# Set up the metadata for DESeq analysis

Counts<- cbind(featureCounts_IgG_n1[,2], featureCounts_SC35_n1[,2],
               featureCounts_IgG_n2[,2], featureCounts_SC35_n2[,2],
               featureCounts_IgG_n3[,2], featureCounts_SC35_n3[,2])

colnames(Counts)<-c("IgG n1", "SC35 n1", 
                    "IgG n2", "SC35 n2", "IgG n3", "SC35 n3")

rownames(Counts)<-featureCounts_IgG_n1[,1]

Protein<- c("IgG", "SC35", "IgG",  "SC35", "IgG",  "SC35")
Condition<- c("n1", "n1", 
              "n2","n2", "n3", "n3")

Metadata<- cbind(Protein, Condition)
rownames(Metadata)<- c("IgG n1", "SC35 n1", 
                       "IgG n2", "SC35 n2", "IgG n3", "SC35 n3")


#-----------------------------------------------------------------------------------------------------------------------------------
# Perform DESeq normalization

Count_Table <- DESeqDataSetFromMatrix(countData = Counts, colData = Metadata,
                                      design = ~ Protein + Condition)

Count_Table_DESeq<- DESeq(Count_Table)

Counts_normalized<- cbind(gtf_df,counts(Count_Table_DESeq, normalized=TRUE))

#-----------------------------------------------------------------------------------------------------------------------------------
# Run the differential enrichment function
Results_SC35_IgG_ALL<- RetrieveAll("SC35", "IgG", "Protein_SC35_vs_IgG", Count_Table_DESeq)

# Save this file which will become Table 2
write.csv(Results_SC35_IgG_ALL, "K562_SC35_differentialenrichment_table.csv")

#-----------------------------------------------------------------------------------------------------------------------------------
# Plot the volcano plot
mycolors <- c("black", "#f49097", "grey70")
names(mycolors) <- c("DOWN", "UP", "NO")

# Label the top 4 H3K27me3 enriched transcripts
Results_SC35_IgG_ALL$delabel <- NA
Results_SC35_IgG_ALL$delabel[Results_SC35_IgG_ALL$gene_name == "MALAT1"] <- "MALAT1"
Results_SC35_IgG_ALL$delabel[Results_SC35_IgG_ALL$gene_name == "MAN2C1"] <- "MAN2C1"
Results_SC35_IgG_ALL$delabel[Results_SC35_IgG_ALL$gene_name == "PI4KAP2"] <- "PI4KAP2"
Results_SC35_IgG_ALL$delabel[Results_SC35_IgG_ALL$gene_name == "HDAC10"] <- "HDAC10"



p<- ggplot(data=Results_SC35_IgG_ALL, aes(x=Results_log2, y=-log10(Results_pvalue), col=diffexpressed,  label=delabel)) + geom_point(size=1) + theme_minimal() + geom_text (size=5, nudge_x= 1, nudge_y = 1) + 
  geom_vline(xintercept=c(-1, 1), col="black") +
  geom_hline(yintercept=-log10(0.05), col="black") +
  theme(text = element_text(size = 20)) +
  labs(x="Fold change (log2) SC35 over IgG", y="FDR (-log10)")+
  theme_classic()

p2 <- p + scale_colour_manual(values = mycolors)


pdf("VolcanoPlot_K562_SC35vsIgG.pdf", 7,5)
p2
dev.off()

#-----------------------------------------------------------------------------------------------------------------------------------
# Count the number of differentially enriched transcipts
SC35_UP<- Results_SC35_IgG_ALL[Results_SC35_IgG_ALL$diffexpressed=="UP",]
nrow(SC35_UP)
# 784

#-----------------------------------------------------------------------------------------------------------------------------------
# Save the coordinates of H3K27me3 enriched transcripts 

K562_SC35_UP_500TES<- get500TES(SC35_UP)
K562_SC35_UP_genebody<- getFullLength(SC35_UP)
K562_SC35_UP_TSS<- getTSS(SC35_UP)

saveBed(K562_SC35_UP_genebody)
saveBed(K562_SC35_UP_500TES)
saveBed(K562_SC35_UP_TSS)


#-----------------------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------------------
# Plot the Venn Diagram showing overlap between H3K27me3 and SC35 enriched transcripts. Used in Figure 1E. 

venn.diagram(
  x = list(H3K27me3_UP$gene_id, SC35_UP$gene_id),
  category.names = c("H3K27me3" , "SC35"),
  filename = 'K562_H3K27me3_vs_SC35_venn_diagramm.tiff',
  output=TRUE
)

# Count the number of all transcripts
all<-list(c(H3K27me3_UP$gene_id, SC35_UP$gene_id))
length(all[[1]])
#2174

common<-all[[1]][duplicated(all[[1]])]
length(common)
#40
