# Perform a differential gene expression analysis of miRNAs with edgeR.

setwd("~/Dropbox/GitHub/RIS/ReedKM_Turkey_miRNA")

library("edgeR")
library("pheatmap")

# Read in the counts matrix
counts <- read.table("Results/miRNA_Expression/miRNA_Expression.txt", header=TRUE, sep="\t")
# Read in the experimental metadata
exp_meta <- read.csv("Resources/Experiment_Metadata.csv", header=TRUE)

# Save the miR names
mir_names <- as.character(counts$Geneid)

# Remove the two samples that are not meant for miRNA analysis
excl_samp <- c("NCT_48_38_38_R_smRNA_S37", "NCT_72_38_R_smRNA_S36")
counts <- counts[,!colnames(counts) %in% excl_samp]
# Remove the first six columns from the subread.txt file
counts <- counts[,-c(1:6)]
# Remove the two samples from the metadata csv too
exp_meta <- exp_meta[!exp_meta$Sample %in% excl_samp,]

# Build the group out of the treatment and time variables for now
exp_meta$Group <- as.factor(paste(exp_meta$Treatment, exp_meta$Time, sep="_"))

# Set up the DEGList object
dge_dat <- DGEList(
    counts=counts,
    genes=mir_names,
    group=exp_meta$Group)

dim(dge_dat)

# Calculate the normalization factors
dge_dat <- calcNormFactors(dge_dat)

# Let's make our first diagnostic plot - the total library size in fragments.
lib_sizes_M <- as.numeric(dge_dat$samples$lib.size) / 1000000
pdf(file="Results/miRNA_Expression/miRNA_Frag_Counts.pdf", height=6, width=6)
par(mar=c(9, 4, 2.5, 0.1))
at <- barplot(
    lib_sizes_M,
    ylim=c(0, max(lib_sizes_M)*1.1),
    col=as.numeric(dge_dat$samples$group),
    xlab="",
    ylab="Library Size (Millions of Fragments)",
    main="Library Sizes of Turkey miRNA Samples",
    axes=FALSE)
axis(side=2)
axis(side=1, at=at, labels=as.character(rownames(dge_dat$samples)), las=2, cex.axis=0.5)
dev.off()


# Do the counts-based filtering. We are being VERY permissive here.
min_cts <- 1
smallest_grp <- min(table(exp_meta$Group))
min_cpm <- log2((1 + min_cts) / min(dge_dat$samples$lib.size) * 1e6)
# Then if a gene has at least N samples with CPM greater than or equal to this
# CPM value, we keep it. N is the size of the smallest group.
cpms <- cpm(dge_dat, log=TRUE, prior.count=1)
keep <- apply(cpms, 1, function(x) {
    cvals <- as.numeric(x)
    pass_flt <- sum(cvals >= min_cpm)
    if(pass_flt >= smallest_grp) {
        return(TRUE)
    } else {
        return(FALSE)
    }
})
# Drop the genes that failed the filtering, and recalculate the library sizes
dge_dat <- dge_dat[keep, ,keep.lib.sizes=FALSE]
dim(dge_dat)
print(min_cpm)

# Make a model matrix based on the group assignments and estimate the
# dispersions
modmat <- model.matrix(~0 + Group, data=exp_meta)
dge_dat <- estimateDisp(dge_dat, design=modmat)

# Plot the biological coefficient of variation (BCV)
pdf(file="Results/miRNA_Expression/miRNA_BCV.pdf", height=6, width=6)
plotBCV(dge_dat)
dev.off()

# Plot a PCA of the samples based on the normalized CPMs for each miRNA
pdf(file="Results/miRNA_Expression/miRNA_PCA.pdf", height=6, width=6)
norm_cpm <- cpm(dge_dat, log=TRUE, prior.count=1, normalized=TRUE)
pc_dat <- prcomp(t(norm_cpm), scale=TRUE)
summary(pc_dat)
pc_1 <- pc_dat$x[,"PC1"]
pc_2 <- pc_dat$x[, "PC2"]
plot(
    pc_2 ~ pc_1,
    xlim=c(-20, 20),
    ylim=c(-10, 10),
    xlab="PC1 (29.4% Var. Exp.)",
    ylab="PC2 (13.9% Var. Exp.)",
    main="Turkey miRNA PCA",
    pch=19,
    col=as.numeric(dge_dat$samples$group))
grps <- as.character(unique(dge_dat$samples$group))
legend("top", grps, pch=19, col=seq_along(grps), ncol=4)
dev.off()
