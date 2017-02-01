args <- commandArgs(trailingOnly=T)
cov_dir <- args[1]

library(gplots)

# Create boxplots of sequencing depth
cov <- read.delim(paste(cov_dir, "coverage.sample_summary", sep="/"))
cov <- subset(cov, sample_id != "Total")
cov <- subset(cov, mean < 500)

pdf(file=paste(cov_dir, "depth_boxplot.pdf", sep="/"), height=5, width=5)
boxplot(cov$mean,
    las=1,
    ylab="Depth")
dev.off()


# Create heatmap of sequencing depth across samples

sample.stats <- read.delim(paste(cov_dir, 
                                 "coverage.sample_statistics", 
                                 sep="/"))
rownames(sample.stats)<-sample.stats$Source_of_reads
sample.stats <- sample.stats[,-1]
sample.stats <- as.matrix(sample.stats)
sample.stats[sample.stats==0] <- 1
sample.stats.log10 <- log10(sample.stats)
pdf(file=paste(cov_dir, "depth_heatmap.pdf", sep="/"), height=17, width=20)
sample.stats.hist <- heatmap.2(
    sample.stats.log10,
    Rowv=TRUE,
    Colv=FALSE,
    dendrogram="row",
    symm=FALSE,
    trace="none",
    keysize = 0.75,
    cexRow=0.35,
    cexCol=0.2,
    xlab="Sequencing Depth Bins (1-500)",
    main="Depth of Sequencing Histograms (log10)")
dev.off()
