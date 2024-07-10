# This code was used to generate Figure 1C, D, E and Supplementary Figure 1A

# Load the appropriate packages
library(rtracklayer)
library(DESeq2)
library(ggplot2)
library(ggrepel)
library(Hmisc)
library(corrplot)
library(VennDiagram)

#-----------------------------------------------------------------------------------------------------------------------------------
# Open up gtf file and extract gene level information
gtf<- import("~/PATH/hg19.gtf")
gtf_df<-as.data.frame(gtf)
gtf_df<-(gtf_df[gtf_df$type=="gene",])[,c(1:6,10,12,14)]

#-----------------------------------------------------------------------------------------------------------------------------------
# Open up Subread count tables. Here perform analysis on read counts that fall within the whole length of the gene
featureCounts_IgG_n1<- read.delim("~/PATH/featureCounts_NK_hs_1983.gene.tabular", comment.char="#")[,c(1,7)]
featureCounts_H3K27me3_n1<- read.delim("~/PATH/featureCounts_NK_hs_1984.gene.tabular", comment.char="#")[,c(1,7)]
featureCounts_SC35_n1<- read.delim("~/PATH/featureCounts_NK_hs_1985.gene.tabular", comment.char="#")[,c(1,7)]

featureCounts_IgG_n2<- read.delim("~/PATH/featureCounts_NK_hs_2004.gene.tabular", comment.char="#")[,c(1,7)]
featureCounts_H3K27me3_n2<- read.delim("~/PATH/featureCounts_NK_hs_2005.gene.tabular", comment.char="#")[,c(1,7)]
featureCounts_SC35_n2<- read.delim("~/PATH/featureCounts_NK_hs_2006.gene.tabular", comment.char="#")[,c(1,7)]

featureCounts_IgG_n3<- read.delim("~/PATH/featureCounts_NK_hs_2030.gene.tabular", comment.char="#")[,c(1,7)]
featureCounts_H3K27me3_n3<- read.delim("~/PATH/featureCounts_NK_hs_2031.gene.tabular", comment.char="#")[,c(1,7)]
featureCounts_SC35_n3<- read.delim("~/PATH/featureCounts_NK_hs_2032.gene.tabular", comment.char="#")[,c(1,7)]

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
# "PC2: 35% variance"

PC_value$labels$x
# "PC1: 55% variance"

# Generate and save the PCA plot
PCA<- plotPCA(rlog(Count_Table_DESeq), intgroup= "Protein", returnData=TRUE)
PCA$name<- c("n1", "n1", "n1", "n2", "n2", "n2", "n3", "n3","n3")
colnames(PCA)<- c("PC1", "PC2", "group","Target", "Replicate")

PCA_plot<-ggplot(PCA, aes(PC1, PC2, color=Target, shape=Replicate)) +
  geom_point(size=5) +
  scale_shape_manual(values=c(16, 17, 15, 18))+
  scale_color_manual(values=c('black','lightblue', "salmon"))+
  xlab(paste0("PC1: 55% variance")) +
  ylab(paste0("PC2: 35% variance")) + 
  coord_fixed() +
theme_classic()

pdf("HEK293T_IgGvsH3K27me3vsSC35_PCAplot.pdf", 6,6)
PCA_plot
dev.off()


#-----------------------------------------------------------------------------------------------------------------------------------
# Determine the correlation values between IgG, H3K27me3, and SC35 RT&Tag libraries
rcor = rcorr(as.matrix(Counts_normalized[,c(10,13,16,11,14,17,12,15,18)]))$r

# Record these correlation values
rcor
#               IgG n1    IgG n2    IgG n3 H3K27me3 n1 H3K27me3 n2 H3K27me3 n3   SC35 n1   SC35 n2   SC35 n3
#IgG n1      1.0000000 0.9743646 0.9471214   0.4195537   0.4530245   0.2160021 0.8121837 0.6930691 0.6990008
#IgG n2      0.9743646 1.0000000 0.9365878   0.4031059   0.4446930   0.2059300 0.7844968 0.6866792 0.6859444
#IgG n3      0.9471214 0.9365878 1.0000000   0.4815796   0.5208806   0.3071915 0.8035978 0.7146212 0.7665111
#H3K27me3 n1 0.4195537 0.4031059 0.4815796   1.0000000   0.9869181   0.9581865 0.5889948 0.5992348 0.5680426
#H3K27me3 n2 0.4530245 0.4446930 0.5208806   0.9869181   1.0000000   0.9421749 0.6194535 0.6412085 0.6140189
#H3K27me3 n3 0.2160021 0.2059300 0.3071915   0.9581865   0.9421749   1.0000000 0.3831520 0.4201751 0.4229928
#SC35 n1     0.8121837 0.7844968 0.8035978   0.5889948   0.6194535   0.3831520 1.0000000 0.9423489 0.8907860
#SC35 n2     0.6930691 0.6866792 0.7146212   0.5992348   0.6412085   0.4201751 0.9423489 1.0000000 0.9465736
#SC35 n3     0.6990008 0.6859444 0.7665111   0.5680426   0.6140189   0.4229928 0.8907860 0.9465736 1.0000000

# Plot and save the correlation matrix. Used in Supplementary Figure 1A. 
pdf("HEK293T_IgGvsH3K27me3vsSC35_correlationplot.pdf", 9,9)
corrplot(rcor, method=c("shade"), tl.col="black", type="lower", col.lim=c(0,1))
dev.off()


#-----------------------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------------------
# Set up the metadata for DESeq analysis for H3K27me3 RT&Tag. This is performed on read counts that fall within whole gene body
# due to the overwhelming amount of intronic signal within H3K27me3 RT&Tag libraries. Will be used for generating volcano plot
# in Figure 1D. 

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
write.csv(Results_H3K27me3_IgG_ALL, "HEK293T_H3K27me3_differentialenrichment_table.csv")

#-----------------------------------------------------------------------------------------------------------------------------------
# Plot the volcano plot
mycolors <- c("black", "lightblue", "grey70")
names(mycolors) <- c("DOWN", "UP", "NO")

# Label the top 4 H3K27me3 enriched transcripts
Results_H3K27me3_IgG_ALL$delabel <- NA
Results_H3K27me3_IgG_ALL$delabel[Results_H3K27me3_IgG_ALL$gene_name == "XIST"] <- "XIST"
Results_H3K27me3_IgG_ALL$delabel[Results_H3K27me3_IgG_ALL$gene_name == "TSIX"] <- "TSIX"
Results_H3K27me3_IgG_ALL$delabel[Results_H3K27me3_IgG_ALL$gene_name == "ROBO2"] <- "ROBO2"
Results_H3K27me3_IgG_ALL$delabel[Results_H3K27me3_IgG_ALL$gene_name == "IMMP2L"] <- "IMMP2L"


p<- ggplot(data=Results_H3K27me3_IgG_ALL, aes(x=Results_log2, y=-log10(Results_pvalue), col=diffexpressed,  label=delabel)) + geom_point(size=1) + theme_minimal() + geom_text (size=5, nudge_x= 1, nudge_y = 1) + 
  geom_vline(xintercept=c(-1, 1), col="black") +
  geom_hline(yintercept=-log10(0.05), col="black") +
  theme(text = element_text(size = 20)) +
  labs(x="Fold change (log2) H3K27me3 over IgG", y="FDR (-log10)") +
  theme_classic()
p2 <- p + scale_colour_manual(values = mycolors)


pdf("VolcanoPlot_HEK293T_H3K27me3vsIgG.pdf", 7,5)
p2
dev.off()

#-----------------------------------------------------------------------------------------------------------------------------------
# Count the number of differentially enriched transcipts
H3K27me3_UP<- Results_H3K27me3_IgG_ALL[Results_H3K27me3_IgG_ALL$diffexpressed=="UP",]
nrow(H3K27me3_UP)
# 3893

#-----------------------------------------------------------------------------------------------------------------------------------
# Save the coordinates of H3K27me3 enriched transcripts 

# Function for saving the TES coordinates
getTES <- function(x){
  x<- x
  for (i in 1:nrow(x)){
    if (x$strand[i] == "+"){
      x$strand2[i] <- "-"
    } else {
      x$strand2[i] <- "+"
    }
  }
  x$strand<- x$strand2
  x<- makeGRangesFromDataFrame(x, keep.extra.columns = TRUE)
  x<- promoters(x, upstream=1, downstream=1)
  seqlevelsStyle(x) <- "UCSC"
  x<- as.data.frame(x)
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


HEK293T_H3K27me3_UP_TES<- getTES(H3K27me3_UP)

HEK293T_H3K27me3_UP_genebody<- getFullLength(H3K27me3_UP)


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


saveBed(HEK293T_H3K27me3_UP_TES)

saveBed(HEK293T_H3K27me3_UP_genebody)




#-----------------------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------------------
# Set up the metadata for DESeq analysis for SC35 RT&Tag. This is performed on read counts that fall within exons as most signal is exonic
# Will be used for generating volcano plot in Figure 1D. 

# Open up Subread count tables. Here perform analysis on read counts that fall within the exons of genes
featureCounts_IgG_n1<- read.delim("~/PATH/featureCounts_NK_hs_1983.exon.tabular", comment.char="#")[,c(1,7)]
featureCounts_SC35_n1<- read.delim("~/PATH/featureCounts_NK_hs_1985.exon.tabular", comment.char="#")[,c(1,7)]

featureCounts_IgG_n2<- read.delim("~/PATH/featureCounts_NK_hs_2004.exon.tabular", comment.char="#")[,c(1,7)]
featureCounts_SC35_n2<- read.delim("~/PATH/featureCounts_NK_hs_2006.exon.tabular", comment.char="#")[,c(1,7)]

featureCounts_IgG_n3<- read.delim("~/PATH/featureCounts_NK_hs_2030.exon.tabular", comment.char="#")[,c(1,7)]
featureCounts_SC35_n3<- read.delim("~/PATH/featureCounts_NK_hs_2032.exon.tabular", comment.char="#")[,c(1,7)]

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
write.csv(Results_SC35_IgG_ALL, "HEK293T_SC35_differentialenrichment_table.csv")

#-----------------------------------------------------------------------------------------------------------------------------------
# Plot the volcano plot
mycolors <- c("black", "#f49097", "grey70")
names(mycolors) <- c("DOWN", "UP", "NO")

# Label the top 4 H3K27me3 enriched transcripts
Results_SC35_IgG_ALL$delabel <- NA
Results_SC35_IgG_ALL$delabel[Results_SC35_IgG_ALL$gene_name == "MALAT1"] <- "MALAT1"
Results_SC35_IgG_ALL$delabel[Results_SC35_IgG_ALL$gene_name == "PNISR"] <- "PNISR"
Results_SC35_IgG_ALL$delabel[Results_SC35_IgG_ALL$gene_name == "LENG9"] <- "LENG9"
Results_SC35_IgG_ALL$delabel[Results_SC35_IgG_ALL$gene_name == "CHTF18"] <- "CHTF18"



p<- ggplot(data=Results_SC35_IgG_ALL, aes(x=Results_log2, y=-log10(Results_pvalue), col=diffexpressed,  label=delabel)) + geom_point(size=1) + theme_minimal() + geom_text (size=5, nudge_x= 1, nudge_y = 1) + 
  geom_vline(xintercept=c(-1, 1), col="black") +
  geom_hline(yintercept=-log10(0.05), col="black") +
  theme(text = element_text(size = 20)) +
  labs(x="Fold change (log2) SC35 over IgG", y="FDR (-log10)")+
  theme_classic()

p2 <- p + scale_colour_manual(values = mycolors)


pdf("VolcanoPlot_HEK293T_SC35vsIgG.pdf", 7,5)
p2
dev.off()

#-----------------------------------------------------------------------------------------------------------------------------------
# Count the number of differentially enriched transcipts
SC35_UP<- Results_SC35_IgG_ALL[Results_SC35_IgG_ALL$diffexpressed=="UP",]
nrow(SC35_UP)
# 1885

#-----------------------------------------------------------------------------------------------------------------------------------
# Save the coordinates of H3K27me3 enriched transcripts 

HEK293T_SC35_UP_TES<- getTES(SC35_UP)
HEK293T_SC35_UP_genebody<- getFullLength(SC35_UP)

saveBed(HEK293T_SC35_UP_genebody)
saveBed(HEK293T_SC35_UP_TES)


#-----------------------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------------------
# Plot the Venn Diagram showing overlap between H3K27me3 and SC35 enriched transcripts. Used in Figure 1E. 

venn.diagram(
  x = list(H3K27me3_UP$gene_id, SC35_UP$gene_id),
  category.names = c("H3K27me3" , "SC35"),
  filename = 'HEK293T_H3K27me3_vs_SC35_venn_diagramm.pdf',
  output=TRUE
)

# Count the number of all transcripts
all<-list(c(H3K27me3_UP$gene_id, SC35_UP$gene_id))
length(all[[1]])
#5778

common<-all[[1]][duplicated(all[[1]])]
length(common)
#351
