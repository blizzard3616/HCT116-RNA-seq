# |----------------------------------------------------------------------------------|
# | Project: HCT116 control vs PEITC treatment                                       |
# | Script: RNA-seq data analysis and visualization using DESeq2                     |
# | Author: Ran Yin                                                                  |
# | Created: 02/05/2019                                                              |
# |----------------------------------------------------------------------------------|
# sink(file = "feature")
date()

# Workflow: https://www.bioconductor.org/help/workflows/rnaseqGene/
# source("https://bioconductor.org/biocLite.R")
# biocLite("DESeq2")

require(data.table)
require(DESeq2)
require(readxl)
require(BiocParallel)
require(ggplot2)
require(knitr)




# NOTE on DESeq2 Output: 'baseMean' is the average of the normalized count values, 
# divided by the size factors, taken over all samples in the DESeqDataSet
# Load data----
dt1 <- fread("./docs/hct116.featurecounts.results.csv",
             skip = 1)
dt1
colnames(dt1)

# Remove unused columns----
# dt1 <- dt1[, c(1, 7:ncol(dt1)), with = FALSE]
# dt1
# colnames(dt1)

# Reorganize the data column to fit legend
dt1_1 <- subset(dt1, 
                select = c(1,7,8,10,11))
colnames(dt1_1)
cnames <- colnames(dt1_1)[-1]
cnames

# Renaming the column names
cnames <- gsub(x = cnames,
               pattern = ".sorted.dedup.bam",
               replacement = "")
cnames

#Select elements from column names
cnames <- gsub(x=cnames,
               pattern = "-to",
               replacement = "")
cnames

colnames(dt1_1)[-1] <- cnames
colnames(dt1_1)

# Read-in legend
library(readxl)
dt2 <- read_excel("./docs/HCT116 legend.xlsx")

# Specify treatment groups
mat <- data.frame(sample = cnames,
                  treatment = dt2$Description,
                  trt = rep(rep(c("np", "peitc"), 
                                each = 2),
                            1),
                  time = factor(rep(c(1), each = 2)),
                  repl = factor(rep(1, each = 2)))
mat

dtm <- as.matrix(dt1_1[, -1, with = FALSE])
rownames(dtm) <- dt1_1$Geneid
row.names(dtm)

# Part I: build a DESeq2 data----
dds <- DESeqDataSetFromMatrix(countData = dtm, 
                              colData = mat,
                              ~ trt + time)
dds

# Remove genes with low counts (<100 in all samples combined)
dds <- dds[rowSums(counts(dds)) >= 100, ]
nrow(counts(dds))
# 15,369 genes left, down from 24,421
# If all samples contain zeros, geometric means cannot be
# estimated. Change default 'type = "ratio"' to 'type = "poscounts"'.
# Type '?DESeq2::estimateSizeFactors' for more details.
dds <- estimateSizeFactors(dds,
                           type = "poscounts")
dds

# Clustering and PCA----
rld <- rlog(dds,
            blind = FALSE)
head(assay(rld))

#Clustering----
sampleDists <- dist(t(assay(rld)))
sampleDists
# 
tiff(filename = "tmp/skin_uvb_samples_cluster.tiff",
     height = 8,
     width = 8,
     units = 'in',
     res = 300,
     compression = "lzw+p")
heatmap(as.matrix(sampleDists), symm = TRUE)
graphics.off()

# PCA----
m1 <- prcomp(t(assay(rld)),
             center = TRUE,
             scale. = TRUE)
summary(m1)
plot(m1)
# 
# Biplot while keep only the most important variables (Javier)----
# Select PC-s to pliot (PC1 & PC2)
choices <- 1:2
# 
# Scores, i.e. points (df.u)
dt.scr <- data.table(m1$x[, choices])
# Add grouping variable
# dt.scr$grp <- colnames(assay(rld))
dt.scr$grp <- mat$treatment
dt.scr$sample <- mat$sample
dt.scr
# 
# Loadings, i.e. arrows (df.v)
dt.rot <- as.data.frame(m1$rotation[, choices])
dt.rot$feat <- rownames(dt.rot)
dt.rot <- data.table(dt.rot)
dt.rot
# 
# Axis labels
u.axis.labs <- paste(colnames(dt.rot)[1:2],
                     sprintf('(%0.1f%% explained var.)',
                             100*m1$sdev[choices]^2/sum(m1$sdev^2)))
u.axis.labs

# Based on Figure p0, keep only a few variables with high loadings in PC1 and PC2----
# var.keep.ndx <- which(dt.rot$feat %in% c("V1",
#                                          "V4",
#                                          "V6",
#                                          "V7",
#                                          "V8"))
gl1 <-dt.rot$feat[order(abs(dt.rot$PC1),
                        decreasing = TRUE)][1:5]
gl2 <-dt.rot$feat[order(abs(dt.rot$PC2),
                        decreasing = TRUE)][1:5]
var.keep.ndx <- which(dt.rot$feat %in% unique(c(gl1, gl2)))
# Or select all
# var.keep.ndx <- 1:ncol(assay(rld))

p1 <- ggplot(data = dt.rot[var.keep.ndx,],
             aes(x = PC1,
                 y = PC2)) +
  # coord_equal() +
  geom_point(data = dt.scr,
             aes(fill = grp),
             shape = 21,
             size = 3,
             alpha = 0.5) +
  # geom_segment(aes(x = 0,
  #                  y = 0,
  #                  xend = 10000*PC1,
  #                  yend = 10000*PC2),
  #              arrow = arrow(length = unit(1/2, 'picas')),
  #              color = "black",
  #              size = 0.5) +
  # geom_text(aes(x = 11000*PC1,
  #               y = 11000*PC2,
  #               label = dt.rot$feat[var.keep.ndx]),
  #           size = 2,
#           hjust = 0.2) +
geom_text(data = dt.scr,
          aes(x = PC1 + 20,
              y = PC2,
              label = dt.scr$sample),
          size = 2.0,
          hjust = 0.5) +
  scale_x_continuous(u.axis.labs[1]) +
  scale_y_continuous(u.axis.labs[2]) +
  scale_fill_manual(name = "Treatment",
                    values = c("light green", "red", "blue")) +
  ggtitle("") +
  theme(plot.title = element_text(hjust = 1.0,
                                  size = 2.75))
p1
tiff(filename = "./figures/SFN-no-tumor.tiff",
     height = 6,
     width = 7,
     units = 'in',
     res = 600,
     compression = "lzw+p")
print(p1)
graphics.off()
