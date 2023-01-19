# Perform a differential gene expression analysis of miRNAs with edgeR.

setwd("~/Dropbox/GitHub/RIS/ReedKM_Turkey_miRNA")

library("edgeR")
library("pheatmap")

# Read in the counts matrix
counts <- read.table("Results/miRNA_Expression/vs_NewGenome/miRNA_Expression.bt.NewGenome.txt", header=TRUE, sep="\t")
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

# Separate the data into 48 and 72 - these are separate experiments
counts_48 <- counts[,grepl("_48_", colnames(counts))]
counts_72 <- counts[,grepl("_72_", colnames(counts))]
exp_meta_48 <- exp_meta[grepl("_48_", exp_meta$Sample),]
exp_meta_72 <- exp_meta[grepl("_72_", exp_meta$Sample),]

# Make the group assignments based on genotype and incubation temp
exp_meta_48$Group <- as.factor(paste(exp_meta_48$Genotype, exp_meta_48$Incubation.Temp, sep="_"))
exp_meta_72$Group <- as.factor(paste(exp_meta_72$Genotype, exp_meta_72$Incubation.Temp, sep="_"))

# Summarize the 48h samples
dge_dat_48 <- DGEList(
    counts=counts_48,
    genes=mir_names,
    group=exp_meta_48$Group)
# # Calculate the normalization factors
dge_dat_48 <- calcNormFactors(dge_dat_48)
# Let's make our first diagnostic plot - the total library size in fragments.
lib_sizes_M <- as.numeric(dge_dat_48$samples$lib.size) / 1000000
pdf(file="Results/miRNA_Expression/vs_NewGenome/48_Fragment_Counts.pdf", height=6, width=6)
par(mar=c(9, 4, 2.5, 0.1))
at <- barplot(
    lib_sizes_M,
    ylim=c(0, max(lib_sizes_M)*1.25),
    col=as.numeric(dge_dat_48$samples$group),
    xlab="",
    ylab="Library Size (Millions of Fragments)",
    main="Library Sizes of 48h Samples",
    axes=FALSE)
axis(side=2)
axis(side=1, at=at, labels=as.character(rownames(dge_dat_48$samples)), las=2, cex.axis=0.5)
legend("top",
    legend=unique(dge_dat_48$samples$group),
    fill=seq_along(unique(dge_dat_48$samples$group)),
    ncol=length(unique(dge_dat_48$samples$group))/2)
dev.off()
# Do the counts-based filtering. We are being permissive here.
min_cts <- 3
smallest_grp <- min(table(exp_meta_48$Group))
min_cpm <- log2((1 + min_cts) / min(dge_dat_48$samples$lib.size) * 1e6)
# Then if a gene has at least N samples with CPM greater than or equal to this
# CPM value, we keep it. N is the size of the smallest group.
cpms <- cpm(dge_dat_48, log=TRUE, prior.count=1)
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
dge_dat_48 <- dge_dat_48[keep, ,keep.lib.sizes=FALSE]
dim(dge_dat_48)
print(min_cpm)

# Make a model matrix based on the group assignments and estimate the
# dispersions
modmat <- model.matrix(~0 + Group, data=exp_meta_48)
dge_dat_48 <- estimateDisp(dge_dat_48, design=modmat)

# Plot the biological coefficient of variation (BCV)
pdf(file="Results/miRNA_Expression/vs_NewGenome/48_BCV.pdf", height=6, width=6)
plotBCV(dge_dat_48)
dev.off()

# Plot a PCA of the samples based on the normalized CPMs for each miRNA
pdf(file="Results/miRNA_Expression/vs_NewGenome/48_PCA.pdf", height=6, width=6)
norm_cpm <- cpm(dge_dat_48, log=TRUE, prior.count=1, normalized=TRUE)
pc_dat <- prcomp(t(norm_cpm), scale=TRUE)
summary(pc_dat)
pc_1 <- pc_dat$x[,"PC1"]
pc_2 <- pc_dat$x[, "PC2"]
plot(
    pc_2 ~ pc_1,
    xlim=c(-10, 25),
    ylim=c(-15, 15),
    xlab="PC1 (35.3% Var. Exp.)",
    ylab="PC2 (20.2% Var. Exp.)",
    main="48h PCA",
    pch=19,
    col=as.numeric(dge_dat_48$samples$group))
grps <- as.character(unique(dge_dat_48$samples$group))
legend("top", grps, pch=19, col=seq_along(grps), ncol=3)
dev.off()

# Do the same for the 72h samples now
dge_dat_72 <- DGEList(
    counts=counts_72,
    genes=mir_names,
    group=exp_meta_72$Group)
# # Calculate the normalization factors
dge_dat_72 <- calcNormFactors(dge_dat_72)
# Let's make our first diagnostic plot - the total library size in fragments.
lib_sizes_M <- as.numeric(dge_dat_72$samples$lib.size) / 1000000
pdf(file="Results/miRNA_Expression/vs_NewGenome/72_Fragment_Counts.pdf", height=6, width=6)
par(mar=c(9, 4, 2.5, 0.1))
at <- barplot(
    lib_sizes_M,
    ylim=c(0, max(lib_sizes_M)*1.25),
    col=as.numeric(dge_dat_72$samples$group),
    xlab="",
    ylab="Library Size (Millions of Fragments)",
    main="Library Sizes of 72h Samples",
    axes=FALSE)
axis(side=2)
axis(side=1, at=at, labels=as.character(rownames(dge_dat_72$samples)), las=2, cex.axis=0.5)
legend("top",
    legend=unique(dge_dat_72$samples$group),
    fill=seq_along(unique(dge_dat_72$samples$group)),
    ncol=3)
dev.off()
# Do the counts-based filtering. We are being permissive here.
min_cts <- 3
smallest_grp <- min(table(exp_meta_72$Group))
min_cpm <- log2((1 + min_cts) / min(dge_dat_72$samples$lib.size) * 1e6)
# Then if a gene has at least N samples with CPM greater than or equal to this
# CPM value, we keep it. N is the size of the smallest group.
cpms <- cpm(dge_dat_72, log=TRUE, prior.count=1)
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
dge_dat_72 <- dge_dat_72[keep, ,keep.lib.sizes=FALSE]
dim(dge_dat_72)
print(min_cpm)

# Make a model matrix based on the group assignments and estimate the
# dispersions
modmat <- model.matrix(~0 + Group, data=exp_meta_72)
dge_dat_72 <- estimateDisp(dge_dat_72, design=modmat)

# Plot the biological coefficient of variation (BCV)
pdf(file="Results/miRNA_Expression/vs_NewGenome/72_BCV.pdf", height=6, width=6)
plotBCV(dge_dat_72)
dev.off()

# Plot a PCA of the samples based on the normalized CPMs for each miRNA
pdf(file="Results/miRNA_Expression/vs_NewGenome/72_PCA.pdf", height=6, width=6)
norm_cpm <- cpm(dge_dat_72, log=TRUE, prior.count=1, normalized=TRUE)
pc_dat <- prcomp(t(norm_cpm), scale=TRUE)
summary(pc_dat)
pc_1 <- pc_dat$x[,"PC1"]
pc_2 <- pc_dat$x[, "PC2"]
plot(
    pc_2 ~ pc_1,
    xlim=c(-20, 20),
    ylim=c(-20, 20),
    xlab="PC1 (36.6% Var. Exp.)",
    ylab="PC2 (19.2% Var. Exp.)",
    main="72h PCA",
    pch=19,
    col=as.numeric(dge_dat_72$samples$group))
grps <- as.character(unique(dge_dat_72$samples$group))
legend("top", grps, pch=19, col=seq_along(grps), ncol=3)
dev.off()
