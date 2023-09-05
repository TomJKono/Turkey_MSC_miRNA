# Perform a differential gene expression analysis of miRNAs with edgeR.

setwd("~/Dropbox/GitHub/RIS/ReedKM_Turkey_miRNA")

library("edgeR")

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

# OK, now let's start the miRNA DE analysis!
# Start with the 72h samples because these look pretty uniform and nice. We
# will do the simplistic approach of running pairwise comparisions for now
comp_1_meta <- exp_meta_72[exp_meta_72$Group %in% c("NCT_33", "NCT_38"),]
comp_1_meta$Group <- factor(comp_1_meta$Group, levels=c("NCT_33", "NCT_38"))
comp_1_dat <- dge_dat_72[,comp_1_meta$Sample]
comp_1_modmat <- model.matrix(~Group, data=comp_1_meta)
comp_1_fit <- glmQLFit(comp_1_dat, design=comp_1_modmat)
comp_1_dge <- glmQLFTest(comp_1_fit, coef=2)
comp_1_table <- topTags(comp_1_dge, n=nrow(comp_1_dge$genes))
write.csv(comp_1_table, file="Results/miRNA_Expression/vs_NewGenome/DEmiRs/72h_NCT38-vs-NCT33.csv", row.names=F, quote=F)

comp_2_meta <- exp_meta_72[exp_meta_72$Group %in% c("NCT_38", "NCT_43"),]
comp_2_meta$Group <- factor(comp_2_meta$Group, levels=c("NCT_38", "NCT_43"))
comp_2_dat <- dge_dat_72[,comp_2_meta$Sample]
comp_2_modmat <- model.matrix(~Group, data=comp_2_meta)
comp_2_fit <- glmQLFit(comp_2_dat, design=comp_2_modmat)
comp_2_dge <- glmQLFTest(comp_2_fit, coef=2)
comp_2_table <- topTags(comp_2_dge, n=nrow(comp_2_dge$genes))
write.csv(comp_2_table, file="Results/miRNA_Expression/vs_NewGenome/DEmiRs/72h_NCT43-vs-NCT38.csv", row.names=F, quote=F)

comp_3_meta <- exp_meta_72[exp_meta_72$Group %in% c("NCT_33", "NCT_43"),]
comp_3_meta$Group <- factor(comp_3_meta$Group, levels=c("NCT_33", "NCT_43"))
comp_3_dat <- dge_dat_72[,comp_3_meta$Sample]
comp_3_modmat <- model.matrix(~Group, data=comp_3_meta)
comp_3_fit <- glmQLFit(comp_3_dat, design=comp_3_modmat)
comp_3_dge <- glmQLFTest(comp_3_fit, coef=2)
comp_3_table <- topTags(comp_3_dge, n=nrow(comp_3_dge$genes))
write.csv(comp_3_table, file="Results/miRNA_Expression/vs_NewGenome/DEmiRs/72h_NCT43-vs-NCT33.csv", row.names=F, quote=F)

# OK, let's do the cross-genotype comparisons
comp_4_meta <- exp_meta_72[exp_meta_72$Group %in% c("NCT_38", "RBC2_38"),]
comp_4_meta$Group <- factor(comp_4_meta$Group, levels=c("NCT_38", "RBC2_38"))
comp_4_dat <- dge_dat_72[,comp_4_meta$Sample]
comp_4_modmat <- model.matrix(~Group, data=comp_4_meta)
comp_4_fit <- glmQLFit(comp_4_dat, design=comp_4_modmat)
comp_4_dge <- glmQLFTest(comp_4_fit, coef=2)
comp_4_table <- topTags(comp_4_dge, n=nrow(comp_4_dge$genes))
write.csv(comp_4_table, file="Results/miRNA_Expression/vs_NewGenome/DEmiRs/72h_RBC238-vs-NCT38.csv", row.names=F, quote=F)

comp_5_meta <- exp_meta_72[exp_meta_72$Group %in% c("NCT_43", "RBC2_43"),]
comp_5_meta$Group <- factor(comp_5_meta$Group, levels=c("NCT_43", "RBC2_43"))
comp_5_dat <- dge_dat_72[,comp_5_meta$Sample]
comp_5_modmat <- model.matrix(~Group, data=comp_5_meta)
comp_5_fit <- glmQLFit(comp_5_dat, design=comp_5_modmat)
comp_5_dge <- glmQLFTest(comp_5_fit, coef=2)
comp_5_table <- topTags(comp_5_dge, n=nrow(comp_5_dge$genes))
write.csv(comp_5_table, file="Results/miRNA_Expression/vs_NewGenome/DEmiRs/72h_RBC243-vs-NCT43.csv", row.names=F, quote=F)

# And now we will do a similar set of comparisons for the 48h samples.
comp_6_meta <- exp_meta_48[exp_meta_48$Group %in% c("NCT_33", "NCT_38"),]
comp_6_meta$Group <- factor(comp_6_meta$Group, levels=c("NCT_33", "NCT_38"))
comp_6_dat <- dge_dat_48[,comp_6_meta$Sample]
comp_6_modmat <- model.matrix(~Group, data=comp_6_meta)
comp_6_fit <- glmQLFit(comp_6_dat, design=comp_6_modmat)
comp_6_dge <- glmQLFTest(comp_6_fit, coef=2)
comp_6_table <- topTags(comp_6_dge, n=nrow(comp_6_dge$genes))
write.csv(comp_6_table, file="Results/miRNA_Expression/vs_NewGenome/DEmiRs/48h_NCT38-vs-NCT33.csv", row.names=F, quote=F)

comp_7_meta <- exp_meta_48[exp_meta_48$Group %in% c("NCT_38", "NCT_43"),]
comp_7_meta$Group <- factor(comp_7_meta$Group, levels=c("NCT_38", "NCT_43"))
comp_7_dat <- dge_dat_48[,comp_7_meta$Sample]
comp_7_modmat <- model.matrix(~Group, data=comp_7_meta)
comp_7_fit <- glmQLFit(comp_7_dat, design=comp_7_modmat)
comp_7_dge <- glmQLFTest(comp_7_fit, coef=2)
comp_7_table <- topTags(comp_7_dge, n=nrow(comp_7_dge$genes))
write.csv(comp_7_table, file="Results/miRNA_Expression/vs_NewGenome/DEmiRs/48h_NCT43-vs-NCT38.csv", row.names=F, quote=F)

comp_8_meta <- exp_meta_48[exp_meta_48$Group %in% c("NCT_33", "NCT_43"),]
comp_8_meta$Group <- factor(comp_8_meta$Group, levels=c("NCT_33", "NCT_43"))
comp_8_dat <- dge_dat_48[,comp_8_meta$Sample]
comp_8_modmat <- model.matrix(~Group, data=comp_8_meta)
comp_8_fit <- glmQLFit(comp_8_dat, design=comp_8_modmat)
comp_8_dge <- glmQLFTest(comp_8_fit, coef=2)
comp_8_table <- topTags(comp_8_dge, n=nrow(comp_8_dge$genes))
write.csv(comp_8_table, file="Results/miRNA_Expression/vs_NewGenome/DEmiRs/48h_NCT43-vs-NCT33.csv", row.names=F, quote=F)

# And the cross-genotype comparisons
comp_9_meta <- exp_meta_48[exp_meta_48$Group %in% c("NCT_33", "RBC2_33"),]
comp_9_meta$Group <- factor(comp_9_meta$Group, levels=c("NCT_33", "RBC2_33"))
comp_9_dat <- dge_dat_48[,comp_9_meta$Sample]
comp_9_modmat <- model.matrix(~Group, data=comp_9_meta)
comp_9_fit <- glmQLFit(comp_9_dat, design=comp_9_modmat)
comp_9_dge <- glmQLFTest(comp_9_fit, coef=2)
comp_9_table <- topTags(comp_9_dge, n=nrow(comp_9_dge$genes))
write.csv(comp_9_table, file="Results/miRNA_Expression/vs_NewGenome/DEmiRs/48h_RBC233-vs-NCT33.csv", row.names=F, quote=F)

comp_10_meta <- exp_meta_48[exp_meta_48$Group %in% c("NCT_38", "RBC2_38"),]
comp_10_meta$Group <- factor(comp_10_meta$Group, levels=c("NCT_38", "RBC2_38"))
comp_10_dat <- dge_dat_48[,comp_10_meta$Sample]
comp_10_modmat <- model.matrix(~Group, data=comp_10_meta)
comp_10_fit <- glmQLFit(comp_10_dat, design=comp_10_modmat)
comp_10_dge <- glmQLFTest(comp_10_fit, coef=2)
comp_10_table <- topTags(comp_10_dge, n=nrow(comp_10_dge$genes))
write.csv(comp_10_table, file="Results/miRNA_Expression/vs_NewGenome/DEmiRs/48h_RBC238-vs-NCT38.csv", row.names=F, quote=F)

comp_11_meta <- exp_meta_48[exp_meta_48$Group %in% c("NCT_43", "RBC2_43"),]
comp_11_meta$Group <- factor(comp_11_meta$Group, levels=c("NCT_43", "RBC2_43"))
comp_11_dat <- dge_dat_48[,comp_11_meta$Sample]
comp_11_modmat <- model.matrix(~Group, data=comp_11_meta)
comp_11_fit <- glmQLFit(comp_11_dat, design=comp_11_modmat)
comp_11_dge <- glmQLFTest(comp_11_fit, coef=2)
comp_11_table <- topTags(comp_11_dge, n=nrow(comp_11_dge$genes))
write.csv(comp_11_table, file="Results/miRNA_Expression/vs_NewGenome/DEmiRs/48h_RBC243-vs-NCT43.csv", row.names=F, quote=F)

# Oops - forgot to do the within-RBC2 comparisions! Do the 48h ones here:
comp_12_meta <- exp_meta_48[exp_meta_48$Group %in% c("RBC2_33", "RBC2_38"),]
comp_12_meta$Group <- factor(comp_12_meta$Group, levels=c("RBC2_33", "RBC2_38"))
comp_12_dat <- dge_dat_48[,comp_12_meta$Sample]
comp_12_modmat <- model.matrix(~Group, data=comp_12_meta)
comp_12_fit <- glmQLFit(comp_12_dat, design=comp_12_modmat)
comp_12_dge <- glmQLFTest(comp_12_fit, coef=2)
comp_12_table <- topTags(comp_12_dge, n=nrow(comp_12_dge$genes))
write.csv(comp_12_table, file="Results/miRNA_Expression/vs_NewGenome/DEmiRs/48h_RBC238-vs-RBC233.csv", row.names=F, quote=F)

comp_13_meta <- exp_meta_48[exp_meta_48$Group %in% c("RBC2_38", "RBC2_43"),]
comp_13_meta$Group <- factor(comp_13_meta$Group, levels=c("RBC2_38", "RBC2_43"))
comp_13_dat <- dge_dat_48[,comp_13_meta$Sample]
comp_13_modmat <- model.matrix(~Group, data=comp_13_meta)
comp_13_fit <- glmQLFit(comp_13_dat, design=comp_13_modmat)
comp_13_dge <- glmQLFTest(comp_13_fit, coef=2)
comp_13_table <- topTags(comp_13_dge, n=nrow(comp_13_dge$genes))
write.csv(comp_13_table, file="Results/miRNA_Expression/vs_NewGenome/DEmiRs/48h_RBC243-vs-RBC238.csv", row.names=F, quote=F)

comp_14_meta <- exp_meta_48[exp_meta_48$Group %in% c("RBC2_33", "RBC2_43"),]
comp_14_meta$Group <- factor(comp_14_meta$Group, levels=c("RBC2_33", "RBC2_43"))
comp_14_dat <- dge_dat_48[,comp_14_meta$Sample]
comp_14_modmat <- model.matrix(~Group, data=comp_14_meta)
comp_14_fit <- glmQLFit(comp_14_dat, design=comp_14_modmat)
comp_14_dge <- glmQLFTest(comp_14_fit, coef=2)
comp_14_table <- topTags(comp_14_dge, n=nrow(comp_14_dge$genes))
write.csv(comp_14_table, file="Results/miRNA_Expression/vs_NewGenome/DEmiRs/48h_RBC243-vs-RBC233.csv", row.names=F, quote=F)

# And the 72h ones here:
comp_15_meta <- exp_meta_72[exp_meta_72$Group %in% c("RBC2_38", "RBC2_43"),]
comp_15_meta$Group <- factor(comp_15_meta$Group, levels=c("RBC2_38", "RBC2_43"))
comp_15_dat <- dge_dat_72[,comp_15_meta$Sample]
comp_15_modmat <- model.matrix(~Group, data=comp_15_meta)
comp_15_fit <- glmQLFit(comp_15_dat, design=comp_15_modmat)
comp_15_dge <- glmQLFTest(comp_15_fit, coef=2)
comp_15_table <- topTags(comp_15_dge, n=nrow(comp_15_dge$genes))
write.csv(comp_15_table, file="Results/miRNA_Expression/vs_NewGenome/DEmiRs/72h_RBC243-vs-RBC238.csv", row.names=F, quote=F)

# Woof! OK, let's do a little more complicated analysis where we analyze the
# 72h and 48h samples as bigger experiments with incubation temp effects and
# genotype effects.
modmat_72 <- model.matrix(~0 + as.factor(Incubation.Temp) + Genotype, data=exp_meta_72)
colnames(modmat_72) <- c("Temp33", "Temp38", "Temp43", "GenotypeRBC2")
dge_dat_72 <- estimateDisp(dge_dat_72, design=modmat_72)
dge_72_fit <- glmQLFit(dge_dat_72, design=modmat_72)
contrasts_72 <- makeContrasts(
    Temp38_vs_Temp33=Temp38-Temp33,
    Temp43_vs_Temp38=Temp43-Temp38,
    Temp43_vs_Temp33=Temp43-Temp33,
    levels=modmat_72)
dge_72_38vs33 <- glmQLFTest(dge_72_fit, contrast=contrasts_72[,"Temp38_vs_Temp33"])
dge_72_38vs33_table <- topTags(dge_72_38vs33, n=nrow(dge_dat_72$genes))$table
dge_72_43vs38 <- glmQLFTest(dge_72_fit, contrast=contrasts_72[,"Temp43_vs_Temp38"])
dge_72_43vs38_table <- topTags(dge_72_43vs38, n=nrow(dge_dat_72$genes))$table
dge_72_43vs33 <- glmQLFTest(dge_72_fit, contrast=contrasts_72[,"Temp43_vs_Temp33"])
dge_72_43vs33_table <- topTags(dge_72_43vs33, n=nrow(dge_dat_72$genes))$table
# Let's merge these three together
colnames(dge_72_38vs33_table) <- c("miRNA.ID", "logFC.38vs33", "logCPM", "F.38vs33", "PValue.38vs33", "FDR.38vs33")
colnames(dge_72_43vs38_table) <- c("miRNA.ID", "logFC.43vs38", "logCPM", "F.43vs38", "PValue.43vs38", "FDR.43vs38")
colnames(dge_72_43vs33_table) <- c("miRNA.ID", "logFC.43vs33", "logCPM", "F.43vs33", "PValue.43vs33", "FDR.43vs33")
dge_72_tab_A <- merge(dge_72_38vs33_table, dge_72_43vs38_table, by=c("miRNA.ID", "logCPM"), all=TRUE)
dge_72_tab_B <- merge(dge_72_tab_A, dge_72_43vs33_table, by=c("miRNA.ID", "logCPM"), all=TRUE)
write.csv(dge_72_tab_B, file="Results/miRNA_Expression/vs_NewGenome/DEmiRs/72h_DE_miRs.csv", row.names=F, quote=F)


modmat_48 <- model.matrix(~0 + as.factor(Incubation.Temp) + Genotype, data=exp_meta_48)
colnames(modmat_48) <- c("Temp33", "Temp38", "Temp43", "GenotypeRBC2")
dge_dat_48 <- estimateDisp(dge_dat_48, design=modmat_48)
dge_48_fit <- glmQLFit(dge_dat_48, design=modmat_48)
contrasts_48 <- makeContrasts(
    Temp38_vs_Temp33=Temp38-Temp33,
    Temp43_vs_Temp38=Temp43-Temp38,
    Temp43_vs_Temp33=Temp43-Temp33,
    levels=modmat_48)
dge_48_38vs33 <- glmQLFTest(dge_48_fit, contrast=contrasts_48[,"Temp38_vs_Temp33"])
dge_48_38vs33_table <- topTags(dge_48_38vs33, n=nrow(dge_dat_48$genes))$table
dge_48_43vs38 <- glmQLFTest(dge_48_fit, contrast=contrasts_48[,"Temp43_vs_Temp38"])
dge_48_43vs38_table <- topTags(dge_48_43vs38, n=nrow(dge_dat_48$genes))$table
dge_48_43vs33 <- glmQLFTest(dge_48_fit, contrast=contrasts_48[,"Temp43_vs_Temp33"])
dge_48_43vs33_table <- topTags(dge_48_43vs33, n=nrow(dge_dat_48$genes))$table
# Let's merge these three together
colnames(dge_48_38vs33_table) <- c("miRNA.ID", "logFC.38vs33", "logCPM", "F.38vs33", "PValue.38vs33", "FDR.38vs33")
colnames(dge_48_43vs38_table) <- c("miRNA.ID", "logFC.43vs38", "logCPM", "F.43vs38", "PValue.43vs38", "FDR.43vs38")
colnames(dge_48_43vs33_table) <- c("miRNA.ID", "logFC.43vs33", "logCPM", "F.43vs33", "PValue.43vs33", "FDR.43vs33")
dge_48_tab_A <- merge(dge_48_38vs33_table, dge_48_43vs38_table, by=c("miRNA.ID", "logCPM"), all=TRUE)
dge_48_tab_B <- merge(dge_48_tab_A, dge_48_43vs33_table, by=c("miRNA.ID", "logCPM"), all=TRUE)
write.csv(dge_48_tab_B, file="Results/miRNA_Expression/vs_NewGenome/DEmiRs/48h_DE_miRs.csv", row.names=F, quote=F)

# Save the edgeR objects for the 72h and 48h experiments; we will use these for
# a variance partitioning analysis.
saveRDS(dge_dat_72, file="Results/miRNA_Expression/vs_NewGenome/72h_edgeR_Dat.rds")
saveRDS(dge_dat_48, file="Results/miRNA_Expression/vs_NewGenome/48h_edgeR_Dat.rds")

# Finally, let's just write the normalized CPMs for each of the experiments
# to do other analyses with them.
cpm_72 <- cpm(dge_dat_72, normalized=TRUE, log=TRUE, prior.count=1)
rownames(cpm_72) <- dge_dat_72$genes[,1]
write.csv(cpm_72, file="Results/miRNA_Expression/vs_NewGenome/CPMs/72h_NormCPM.csv", row.names=TRUE, quote=FALSE)
cpm_48 <- cpm(dge_dat_48, normalized=TRUE, log=TRUE, prior.count=1)
rownames(cpm_48) <- dge_dat_48$genes[,1]
write.csv(cpm_48, file="Results/miRNA_Expression/vs_NewGenome/CPMs/48h_NormCPM.csv", row.names=TRUE, quote=FALSE)

