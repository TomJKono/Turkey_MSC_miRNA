# Run a variancePartition analysis of the miRNA expression.

setwd("~/Dropbox/GitHub/RIS/ReedKM_Turkey_miRNA")
library(variancePartition)
library(edgeR)
library(doParallel)

# Run 8 worker threads
cl <- makeCluster(8)
registerDoParallel(cl)

# Read in the experimental metadata and the edgeR expression objects
exp_meta <- read.csv("Resources/Experiment_Metadata.csv", header=TRUE)
dge_dat_72 <- readRDS("Results/miRNA_Expression/vs_NewGenome/72h_edgeR_Dat.rds")
dge_dat_48 <- readRDS("Results/miRNA_Expression/vs_NewGenome/48h_edgeR_Dat.rds")

# Remove the two samples that are not meant for miRNA analysis
excl_samp <- c("NCT_48_38_38_R_smRNA_S37", "NCT_72_38_R_smRNA_S36")
exp_meta <- exp_meta[!exp_meta$Sample %in% excl_samp,]

# Subset the metadata sheet to separate the 72h and 48h experiments
exp_meta_48 <- exp_meta[grepl("_48_", exp_meta$Sample),]
exp_meta_72 <- exp_meta[grepl("_72_", exp_meta$Sample),]

# Build the model matrices
exp_meta_48$Incubation.Temp <- factor(exp_meta_48$Incubation.Temp, levels=c("33", "38", "43"))
modmat_48 <- model.matrix(~Incubation.Temp + Genotype + Incubation.Temp:Genotype, data=exp_meta_48)
exp_meta_72$Incubation.Temp <- factor(exp_meta_72$Incubation.Temp, levels=c("33", "38", "43"))
modmat_72 <- model.matrix(~Incubation.Temp + Genotype + Incubation.Temp:Genotype, data=exp_meta_72)

# use voom() to get the normalized expression values and the mean-variance
# relationships for fitting the varpart model
voom_72 <- voom(dge_dat_72, modmat_72)
voom_48 <- voom(dge_dat_48, modmat_48)

# Fit the varpart models
vp_mod_72 <- ~(1|Incubation.Temp) + (1|Genotype) + (1|Incubation.Temp:Genotype)
vp_fit_72 <- fitExtractVarPartModel(voom_72, vp_mod_72, exp_meta_72)
pdf(file="Results/miRNA_Expression/vs_NewGenome/72h_Variance_Partition.pdf", height=6, width=6)
plotVarPart(vp_fit_72)
dev.off()

vp_mod_48 <- ~(1|Incubation.Temp) + (1|Genotype) + (1|Incubation.Temp:Genotype)
vp_fit_48 <- fitExtractVarPartModel(voom_48, vp_mod_48, exp_meta_48)
pdf(file="Results/miRNA_Expression/vs_NewGenome/48h_Variance_Partition.pdf", height=6, width=6)
plotVarPart(vp_fit_48)
dev.off()

# Save the variance partition results to CSV
vp_fit_72_mat <- as.matrix(vp_fit_72)
rownames(vp_fit_72_mat) <- dge_dat_72$genes[,1]
write.csv(vp_fit_72_mat, file="Results/miRNA_Expression/vs_NewGenome/72h_Variance_Partition.csv", row.names=TRUE, quote=FALSE)
vp_fit_48_mat <- as.matrix(vp_fit_48)
rownames(vp_fit_48_mat) <- dge_dat_48$genes[,1]
write.csv(vp_fit_48_mat, file="Results/miRNA_Expression/vs_NewGenome/48h_Variance_Partition.csv", row.names=TRUE, quote=FALSE)
