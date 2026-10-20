### M. birostris analysis to produce publication figure ####
# M. Eppley, v4. current version uploaded to github 10/20/2026
library(vcfR)
library(pcadapt)
library(SNPRelate)
library(ggplot2)
library(pheatmap)
library(LEA)
library(adegenet)
library(ggrepel)
library(dplyr)
library(tidyr)
library(RColorBrewer)
library(dartR)
library(reshape2)
library(poppr)
library(ade4)
library(patchwork)
library(grid)
library(gridExtra)
library(maps)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(qvalue)

outdir <- "/Users/madelineeppley/Desktop/manta26pub"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
setwd(outdir)

pop_colors <- c("IND" = "#6DBAA4", "PER" = "#8C9FCB", "REV" = "#ED9A6C", "BYC" = "#DA8EC0")
structure_colors <- c("#264653", "#2a9d8f", "#8ab17d", "#e9c46a", "#f4a261", "#e76f51")

# construct the pop map
popmap <- data.frame(
  sample = c("BYC_RMB_01",
             "IN_1_MB","IN_2_MB","IN_4_MB","IN_5_MB","IN_6_MB","IN_7_MB",
             "PER_001_MB","PER_003_MB","PER_004_MB","PER_005_MB","PER_006_MB",
             "PER_007_MB","PER_DZW81_4_MB",
             "REV_10_MB","REV_13_MB","REV_14_MB_B","REV_15_MB",
             "REV_17_MB","REV_18_MB","REV_19_MB","REV_20_MB_B"),
  pop = c("BYC", rep("IND",6), rep("PER",7), rep("REV",8)))

# read original, subset, write out, then clean vcf
vcf_raw <- read.vcfR("/Users/madelineeppley/Desktop/manta/finalvcfs/birostris_final.vcf")
keep_cols <- c(TRUE, colnames(vcf_raw@gt)[-1] %in% popmap$sample)
vcf_raw@gt <- vcf_raw@gt[, keep_cols, drop = FALSE]
write.vcf(vcf_raw, "/Users/madelineeppley/Desktop/manta26pub/birostris_22samp.vcf.gz")
system("gunzip -f /Users/madelineeppley/Desktop/manta26pub/birostris_22samp.vcf.gz")
rm(vcf_raw)

# vcf back in
vcf_bir <- read.vcfR("/Users/madelineeppley/Desktop/manta26pub/birostris_22samp.vcf")

gl_bir <- vcfR2genlight(vcf_bir)
pop(gl_bir) <- popmap$pop[match(indNames(gl_bir), popmap$sample)]
n_snps_all <- nrow(vcf_bir@fix)
#print(nInd(gl_bir))
#print(n_snps_all)
print(table(pop(gl_bir)))

pop_labels_df <- data.frame(pop = c("IND","PER","REV","BYC"),
                            pop_label = c("India","Peru","Mexico","Bycatch"))

# sampling map, coordinates from bycatch df
sample_coords <- data.frame(
  sample = popmap$sample,
  lat = c(-4.467, 10.5,10.8,11.2,10.3,10.9,11.1,
          -6.0,-7.5,-8.0,-9.5,-10.0,-11.5,-13.5,
          18.7,18.9,18.6,19.0,18.5,18.8,19.1,18.7),
  lon = c(-99.017, 72.5,72.8,73.2,72.3,72.9,73.1,
          -81.0,-80.5,-81.2,-80.8,-81.5,-80.0,-80.5,
          -111.0,-110.5,-111.2,-110.3,-111.5,-110.9,-110.6,-111.1),
  pop = popmap$pop) %>%
  mutate(lon_pacific = ifelse(lon < 0, lon + 360, lon))

world_map <- map_data("world", wrap = c(0, 360))
pop_n <- sample_coords %>% group_by(pop) %>% summarise(n=n())

range_map <- ggplot() +
  geom_polygon(data=world_map, aes(x=long, y=lat, group=group), fill="gray90", color="gray70", linewidth=0.2) +
  geom_point(data=sample_coords, aes(x=lon_pacific, y=lat, fill=pop), shape=21, color="black", size=3.5, stroke=0.5, alpha=0.9) +
  scale_fill_manual(values=pop_colors, name="Sampling\nLocation", breaks=c("IND","PER","REV","BYC"),
                    labels=paste0(c("India","Peru","Mexico","Bycatch"), " (n=", pop_n$n[match(c("IND","PER","REV","BYC"), pop_n$pop)], ")")) +
  coord_fixed(ratio=1, xlim=c(30, 310), ylim=c(-50, 50), expand=FALSE) +
  scale_x_continuous(breaks=c(60, 120, 180, 240, 300), labels=c("60°E","120°E","180°","120°W","60°W")) +
  scale_y_continuous(breaks=c(-40, -20, 0, 20, 40), labels=c("40°S","20°S","0°","20°N","40°N")) +
  theme_minimal() + theme(legend.position="right", legend.title=element_text(size=13, face="bold"),
                          legend.text=element_text(size=12),
                          panel.grid.major=element_line(color="gray90", linewidth=0.3), panel.background=element_rect(fill="aliceblue", color=NA),
                          axis.title=element_blank(), axis.text=element_text(size=9)) + guides(fill=guide_legend(override.aes=list(size=5)))
ggsave("/Users/madelineeppley/Desktop/manta26pub/birostris_map.png", range_map, width=12, height=4, dpi=600)

# relatedness plot
gds.fn <- "/Users/madelineeppley/Desktop/manta26pub/bir.gds"
showfile.gds(closeall=TRUE)
snpgdsVCF2GDS("/Users/madelineeppley/Desktop/manta26pub/birostris_22samp.vcf", gds.fn, method="biallelic.only")
genofile <- snpgdsOpen(gds.fn)
rel <- snpgdsIBDKING(genofile, autosome.only=FALSE)
snpgdsClose(genofile)
kin_mat <- rel$kinship; sample_names_rel <- rel$sample.id
rownames(kin_mat) <- sample_names_rel; colnames(kin_mat) <- sample_names_rel
annotation_df <- data.frame(Population = popmap$pop[match(sample_names_rel, popmap$sample)])
rownames(annotation_df) <- sample_names_rel
relate <- pheatmap(kin_mat, labels_row=sample_names_rel, labels_col=sample_names_rel,
                   annotation_row=annotation_df, annotation_col=annotation_df,
                   annotation_colors=list(Population=pop_colors),
                   clustering_distance_rows="euclidean", clustering_distance_cols="euclidean",
                   main="Pairwise Genomic Relatedness, M. birostris", silent=TRUE)
write.csv(kin_mat, "/Users/madelineeppley/Desktop/manta26pub/birostris_kinship_matrix.csv")
diag(kin_mat) <- NA # set as NA so it doesnt return 0.5 as the highest relatedness
max_rel <- which(kin_mat == max(kin_mat, na.rm=TRUE), arr.ind=TRUE)[1,]
print(round(max(kin_mat, na.rm=TRUE), 4)) # highest relatedness 0.0936
print(round(min(kin_mat, na.rm=TRUE), 4)) # lowest relatedness 0.0197
print(round(mean(kin_mat, na.rm=TRUE), 4)) # mean relatedness 0.0469
print(rownames(kin_mat)[max_rel[1]]) # REV_15_MB
print(colnames(kin_mat)[max_rel[2]]) # REV_10_MB
ggsave("/Users/madelineeppley/Desktop/manta26pub/birostris_relatedness.png", relate, width=10, height=8, dpi=600)

# PCA with all SNPs
pca_all <- glPca(gl_bir, nf=10)
eig_all <- pca_all$eig / sum(pca_all$eig) * 100
pca_df_all <- data.frame(PC1=pca_all$scores[,1], PC2=pca_all$scores[,2],
                         sample=indNames(gl_bir), pop=pop(gl_bir)) %>% left_join(pop_labels_df, by="pop")
centroid_all <- pca_df_all %>% group_by(pop) %>% summarise(PC1_cen=mean(PC1), PC2_cen=mean(PC2), n=n())
pca_df_all <- pca_df_all %>% left_join(centroid_all, by="pop")
pops_ell_all <- centroid_all %>% filter(n >= 3) %>% pull(pop)

pca_plot_all <- ggplot(pca_df_all, aes(x=PC1, y=PC2)) +
  stat_ellipse(data=filter(pca_df_all, pop %in% pops_ell_all), aes(color=pop),
               type="norm", level=0.95, linewidth=1.2, show.legend=FALSE) +
  geom_segment(aes(xend=PC1_cen, yend=PC2_cen, color=pop), linewidth=0.5, alpha=0.6, show.legend=FALSE) +
  geom_point(aes(fill=pop, color=pop), size=4, alpha=0.8, shape=21, stroke=0.5) +
  scale_fill_manual(values=pop_colors, breaks=c("IND","PER","REV","BYC"), labels=c("India","Peru","Mexico","Bycatch")) +
  scale_color_manual(values=pop_colors, breaks=c("IND","PER","REV","BYC"), labels=c("India","Peru","Mexico","Bycatch")) +
  labs(x=paste0("PC1 (",round(eig_all[1],1),"%)"), y=paste0("PC2 (",round(eig_all[2],1),"%)"), subtitle="All SNPs") +
  theme_classic() +
  theme(legend.position="none", axis.text=element_text(size=13, color="black"),
        axis.title=element_text(size=14, face="bold"), plot.subtitle=element_text(size=13, hjust=0.5),
        panel.border=element_rect(color="black", fill=NA, linewidth=1),
        panel.grid.major=element_line(color="grey90", linewidth=0.3))
ggsave("/Users/madelineeppley/Desktop/manta26pub/birostris_PCA_allSNPs.png",
       pca_plot_all + labs(title="Mobula birostris") +
         theme(legend.position="bottom", legend.direction="horizontal", legend.title=element_blank(),
               plot.title=element_text(size=13, face="italic", hjust=0.5)) +
         guides(fill=guide_legend(override.aes=list(shape=21, size=5, alpha=0.8))),
       width=8, height=6, dpi=600)

# pcadapt for outlier loci
geno_bir <- read.pcadapt("/Users/madelineeppley/Desktop/manta26pub/birostris_22samp.vcf", type = "vcf")
obj <- pcadapt(geno_bir, K = 1)
pvals <- obj$pvalues
qvals <- qvalue(pvals)$qvalues
alpha <- 0.01
outliers <- which(qvals < alpha)
n_outliers <- length(outliers)
n_outliers # 4 outliers
outliers # 1081 1676 4995 8161

# updated 10/7/26
padj <- p.adjust(pvals, method = "bonferroni")
alphab <- 0.1
outliers_bonf <- which(padj < alphab)
length(outliers_bonf) #4 outliers
outliers_bonf # 1081 1676 4995 8161 same as above

# vcf with outliers
vcf_bir_outliers <- vcf_bir
vcf_bir_outliers@fix <- vcf_bir@fix[outliers, , drop = FALSE]
vcf_bir_outliers@gt  <- vcf_bir@gt[outliers, , drop = FALSE]
write.vcf(vcf_bir_outliers,
          "/Users/madelineeppley/Desktop/manta26pub/birostris_outliers.vcf.gz")
system("gunzip -f /Users/madelineeppley/Desktop/manta26pub/birostris_outliers.vcf.gz")

# PCA with outlier SNPs
vcf_outliers <- read.vcfR("/Users/madelineeppley/Desktop/manta26pub/birostris_outliers.vcf")
gl_outliers <- vcfR2genlight(vcf_outliers)
pop(gl_outliers) <- popmap$pop[match(indNames(gl_outliers), popmap$sample)]

pca_out <- glPca(gl_outliers, nf=10)
eig_out <- pca_out$eig / sum(pca_out$eig) * 100
pca_df_out <- data.frame(PC1=pca_out$scores[,1], PC2=pca_out$scores[,2],
                         sample=indNames(gl_outliers), pop=pop(gl_outliers)) %>% left_join(pop_labels_df, by="pop")
centroid_out <- pca_df_out %>% group_by(pop) %>% summarise(PC1_cen=mean(PC1), PC2_cen=mean(PC2), n=n())
pca_df_out <- pca_df_out %>% left_join(centroid_out, by="pop")
pops_ell_out <- centroid_out %>% filter(n >= 3) %>% pull(pop)

pca_plot_outliers <- ggplot(pca_df_out, aes(x=PC1, y=PC2)) +
  stat_ellipse(data=filter(pca_df_out, pop %in% pops_ell_out), aes(color=pop),
               type="norm", level=0.95, linewidth=1.2, show.legend=FALSE) +
  geom_segment(aes(xend=PC1_cen, yend=PC2_cen, color=pop), linewidth=0.5, alpha=0.6, show.legend=FALSE) +
  geom_point(aes(fill=pop, color=pop), size=4, alpha=0.8, shape=21, stroke=0.5) +
  scale_fill_manual(values=pop_colors, breaks=c("IND","PER","REV","BYC"), labels=c("India","Peru","Mexico","Bycatch")) +
  scale_color_manual(values=pop_colors, breaks=c("IND","PER","REV","BYC"), labels=c("India","Peru","Mexico","Bycatch")) +
  labs(x=paste0("PC1 (",round(eig_out[1],1),"%)"), y=paste0("PC2 (",round(eig_out[2],1),"%)"), subtitle="Outlier SNPs") +
  theme_classic() +
  theme(legend.position="none", axis.text=element_text(size=13, color="black"),
        axis.title=element_text(size=14, face="bold"), plot.subtitle=element_text(size=13, hjust=0.5),
        panel.border=element_rect(color="black", fill=NA, linewidth=1),
        panel.grid.major=element_line(color="grey90", linewidth=0.3))
ggsave("/Users/madelineeppley/Desktop/manta26pub/birostris_PCA_outlierSNPs.png",
       pca_plot_outliers + labs(title="Mobula birostris") +
         theme(legend.position="bottom", legend.direction="horizontal", legend.title=element_blank(),
               plot.title=element_text(size=13, face="italic", hjust=0.5)) +
         guides(fill=guide_legend(override.aes=list(shape=21, size=5, alpha=0.8))),
       width=8, height=6, dpi=600)

# DAPC analysis
dapc_result <- dapc(gl_bir, pop=pop(gl_bir), n.pca=4, n.da=3) # keeping 4 PCs and pcs - 1 = DA
dapc_df <- data.frame(LD1=dapc_result$ind.coord[,1], LD2=dapc_result$ind.coord[,2],
                      sample=indNames(gl_bir), pop=pop(gl_bir))
centroid_dapc <- dapc_df %>% group_by(pop) %>% summarise(LD1_cen=mean(LD1), LD2_cen=mean(LD2), n=n())
dapc_df <- dapc_df %>% left_join(centroid_dapc, by="pop")
pops_ell_dapc <- centroid_dapc %>% filter(n >= 3) %>% pull(pop)
pct_dapc <- dapc_result$eig / sum(dapc_result$eig) * 100

dapc_plot <- ggplot(dapc_df, aes(x=LD1, y=LD2)) +
  stat_ellipse(data=filter(dapc_df, pop %in% pops_ell_dapc), aes(color=pop),
               type="norm", level=0.95, linewidth=1.2, show.legend=FALSE) +
  geom_segment(aes(xend=LD1_cen, yend=LD2_cen, color=pop), linewidth=0.5, alpha=0.6, show.legend=FALSE) +
  geom_point(aes(fill=pop, color=pop), size=4, alpha=0.8, shape=21, stroke=0.5) +
  scale_fill_manual(values=pop_colors, breaks=c("IND","PER","REV","BYC"), labels=c("India","Peru","Mexico","Bycatch")) +
  scale_color_manual(values=pop_colors, breaks=c("IND","PER","REV","BYC"), labels=c("India","Peru","Mexico","Bycatch")) +
  labs(x=paste0("DAPC1 (",round(pct_dapc[1],1),"%)"), y=paste0("DAPC2 (",round(pct_dapc[2],1),"%)"), subtitle="All SNPs") +
  theme_classic() +
  theme(legend.position="none", axis.text=element_text(size=13, color="black"),
        axis.title=element_text(size=14, face="bold"),
        panel.border=element_rect(color="black", fill=NA, linewidth=1),
        panel.grid.major=element_line(color="grey90", linewidth=0.3))
ggsave("/Users/madelineeppley/Desktop/manta26pub/birostris_DAPC.png",
       dapc_plot + labs(title="Mobula birostris") +
         theme(legend.position="bottom", legend.direction="horizontal", legend.title=element_blank(),
               plot.title=element_text(size=13, face="italic", hjust=0.5)) +
         guides(fill=guide_legend(override.aes=list(shape=21, size=5, alpha=0.8))),
       width=8, height=6, dpi=600)

# ancestry proportions with SNMF on all snps
## JUST RELOAD
# needed to double-check k values, re-load
project_all <- load.snmfProject("/Users/madelineeppley/Desktop/manta26pub/birostris.snmfProject")
project_outliers <- load.snmfProject("/Users/madelineeppley/Desktop/manta26pub/birostris_outliers.snmfProject")
for(k in 1:5) cat("K =", k, ":", min(cross.entropy(project_all, K=k)), "\n") # lowest k = 1, 0.749
for(k in 1:5) cat("K =", k, ":", min(cross.entropy(project_outliers, K=k)), "\n") # lowest k = 3, 0.040

# to run again from the start:
#unlink("/Users/madelineeppley/Desktop/manta26pub/birostris.geno") # 'unlink' worked here when i had to re-run analyses otherwise snmf acts weird and overwrites
#unlink("/Users/madelineeppley/Desktop/manta26pub/birostris.snmfProject", recursive=TRUE)
#unlink("/Users/madelineeppley/Desktop/manta26pub/birostris.snmf", recursive=TRUE)
#vcf2geno("/Users/madelineeppley/Desktop/manta26pub/birostris_22samp.vcf",
#        "/Users/madelineeppley/Desktop/manta26pub/birostris.geno")
#project_all <- snmf("/Users/madelineeppley/Desktop/manta26pub/birostris.geno",
#                  K = 1:5, entropy = TRUE, repetitions = 10, project = "new")

# cross-entropy plot all SNPs
png("/Users/madelineeppley/Desktop/manta26pub/birostris_CE_allSNPs.png", width=800, height=600)
plot(project_all, col="blue", pch=19, cex=1.2, main="Cross-entropy: M. birostris - All SNPs")
dev.off()
print("cross entropy all snps")
for(k in 1:5) cat("K =", k, ":", min(cross.entropy(project_all, K=k)), "\n")

# outlier snps
#unlink("/Users/madelineeppley/Desktop/manta26pub/birostris_outliers.geno")
#unlink("/Users/madelineeppley/Desktop/manta26pub/birostris_outliers.snmfProject", recursive=TRUE)
#unlink("/Users/madelineeppley/Desktop/manta26pub/birostris_outliers.snmf", recursive=TRUE)
#vcf2geno("/Users/madelineeppley/Desktop/manta26pub/birostris_outliers.vcf",
#        "/Users/madelineeppley/Desktop/manta26pub/birostris_outliers.geno")
#project_outliers <- snmf("/Users/madelineeppley/Desktop/manta26pub/birostris_outliers.geno",
#                        K = 1:5, entropy = TRUE, repetitions = 10, project = "new")

# cross-entropy plot outlier snps
png("/Users/madelineeppley/Desktop/manta26pub/birostris_CE_outlierSNPs.png", width=800, height=600)
plot(project_outliers, col="red", pch=19, cex=1.2, main="Cross-entropy: M. birostris - Outlier SNPs")
dev.off()
print("cross entropy outlier snps")
for(k in 1:5) cat("K =", k, ":", min(cross.entropy(project_outliers, K=k)), "\n")

# structure plot with group names
sample_names <- indNames(gl_bir)
sample_pops <- pop(gl_bir)
pop_order <- c("IND","PER","REV","BYC")
pop_names_map <- c("IND"="India","PER"="Peru","REV"="Mexico","BYC"="Byc")

pop_label_df_str <- data.frame(Sample = sample_names,
                               Pop = factor(sample_pops, levels = pop_order)) %>%
  arrange(Pop) %>% mutate(x_pos = 1:n())
pop_boundaries <- pop_label_df_str %>% group_by(Pop) %>% summarise(x_end = max(x_pos)) %>%
  filter(x_end != max(pop_label_df_str$x_pos))
pop_spans <- pop_label_df_str %>% group_by(Pop) %>%
  summarise(x_start=min(x_pos), x_end=max(x_pos), x_mid=(min(x_pos)+max(x_pos))/2, n=n()) %>%
  mutate(display_name = pop_names_map[as.character(Pop)])
max_x <- max(pop_label_df_str$x_pos)

# K=2 all snps
ce_all <- cross.entropy(project_all, K = 2)
best_run_all <- which.min(ce_all)
qmat_all <- Q(project_all, K = 2, run = best_run_all)
qdf_all <- as.data.frame(qmat_all)
colnames(qdf_all) <- paste0("Cluster", 1:2)
qdf_all$Sample <- sample_names
qdf_all$Pop <- factor(sample_pops, levels = pop_order)
qdf_all <- qdf_all %>% arrange(Pop) %>% mutate(x_pos = 1:n())
qdf_all_long <- qdf_all %>% pivot_longer(cols = starts_with("Cluster"),
                                         names_to = "Cluster", values_to = "Proportion")

structure_all <- ggplot(qdf_all_long, aes(x = x_pos, y = Proportion, fill = Cluster)) +
  geom_bar(stat = "identity", width = 1, color = NA) +
  geom_vline(data = pop_boundaries, aes(xintercept = x_end + 0.5), color = "white", linewidth = 0.8) +
  scale_fill_manual(values = structure_colors[1:2]) +
  scale_y_continuous(expand = c(0,0), breaks = c(0,0.5,1), labels = c("0.0","0.5","1.0"), position = "left") +
  scale_x_continuous(expand = c(0,0), limits = c(0.5, max_x + 3)) +
  annotate("text", x = max_x + 2, y = 0.5, label = "K = 2", size = 6, fontface = "bold") +
  ggtitle("All SNPs") + coord_cartesian(clip = "off") +
  theme_minimal() +
  theme(plot.title=element_text(size=12, face="bold", hjust=0),
        axis.text.x=element_blank(), axis.ticks.x=element_blank(),
        axis.title=element_blank(), axis.text.y=element_text(size=12),
        legend.position="none", panel.grid=element_blank(),
        plot.margin=margin(5, 30, 2, 5))

# K=3 outlier snps
ce_out <- cross.entropy(project_outliers, K = 3)
best_run_out <- which.min(ce_out)
qmat_out <- Q(project_outliers, K = 3, run = best_run_out)
qdf_out_str <- as.data.frame(qmat_out)
colnames(qdf_out_str) <- paste0("Cluster", 1:3)
qdf_out_str$Sample <- sample_names
qdf_out_str$Pop <- factor(sample_pops, levels = pop_order)
qdf_out_str <- qdf_out_str %>% arrange(Pop) %>% mutate(x_pos = 1:n())
qdf_out_long <- qdf_out_str %>% pivot_longer(cols = starts_with("Cluster"),
                                             names_to = "Cluster", values_to = "Proportion")

structure_outlier <- ggplot(qdf_out_long, aes(x = x_pos, y = Proportion, fill = Cluster)) +
  geom_bar(stat = "identity", width = 1, color = NA) +
  geom_vline(data = pop_boundaries, aes(xintercept = x_end + 0.5), color = "white", linewidth = 0.8) +
  scale_fill_manual(values = structure_colors[c(3,4,5)]) +
  scale_y_continuous(expand = c(0,0), breaks = c(0,0.5,1), labels = c("0.0","0.5","1.0"), position = "left") +
  scale_x_continuous(expand = c(0,0), limits = c(0.5, max_x + 3)) +
  annotate("text", x = max_x + 2, y = 0.5, label = "K = 3", size = 6, fontface = "bold") +
  ggtitle("Outlier SNPs") + coord_cartesian(clip = "off") +
  theme_minimal() +
  theme(plot.title=element_text(size=12, face="bold", hjust=0),
        axis.text.x=element_blank(), axis.ticks.x=element_blank(),
        axis.title=element_blank(), axis.text.y=element_text(size=12),
        legend.position="none", panel.grid=element_blank(),
        plot.margin=margin(5, 30, 2, 5))

pop_label_plot <- ggplot() +
  geom_text(data = pop_spans, aes(x = x_mid, y = 0.5, label = display_name), size = 5, fontface = "bold") +
  geom_vline(data = pop_boundaries, aes(xintercept = x_end + 0.5), color = "black", linewidth = 1) +
  scale_x_continuous(limits = c(0.5, max_x + 3), expand = c(0,0)) +
  scale_y_continuous(limits = c(0,1), expand = c(0,0)) +
  coord_cartesian(clip = "off") + theme_void() + theme(plot.margin = margin(0, 30, 10, 5))

# amova FST
genI <- gl2gi(gl_bir)
pops_gi <- genI$pop; strata(genI) <- data.frame(pops = pops_gi)
p.amova <- poppr.amova(genI, ~pops)
amova.pvalues <- ade4::randtest(p.amova, nrepet = 9999)
overall_fst <- p.amova$statphi$Phi[length(p.amova$statphi$Phi)] # 0.003964456
# updated 10/7/26
print(amova.pvalues) # take row 3 from here for variation between pops 
overall_pval <- amova.pvalues$pvalue[length(amova.pvalues$pvalue)]
overall_pval # 9e-04

# fst but excluding BYC (n=1)
pop_sizes <- table(pop(gl_bir))
valid_pops <- names(pop_sizes[pop_sizes >= 2])
gl_valid <- gl_bir[pop(gl_bir) %in% valid_pops, ]
fst_result <- gl.fst.pop(gl_valid, nboots=1000, percent=95, nclusters=1)
fst_matrix <- fst_result$Fsts
fst_matrix[is.na(fst_matrix)] <- t(fst_matrix)[is.na(fst_matrix)]

# updated 10/7/26
pval_matrix <- fst_result$Pvalues
pval_matrix[is.na(pval_matrix)] <- t(pval_matrix)[is.na(pval_matrix)]
print("pairwise fst")
print(round(fst_matrix, 4))
print(pval_matrix)
## ok so noting the revision here, the results show significant comparisons between india-peru and india-mexico but not peru-mexico
# FST India vs Peru = 0.0067, p = 0 significant
# FST India vs Mexico = 0.0055, p = 0 significant 
# FST Peru vs Mexico = 0.0006, p= 0.278 not significant 

## outlier fst
table(ploidy(gl_outliers)) # 18 of the 22 inds are assigned ploidy 1 because there are no heterozygous calls across the 4 SNPs
ploidy(gl_outliers) <- 2 # set ploidy to 2 for all
nLoc(gl_outliers) # 4 outliers
pop(gl_outliers) <- popmap$pop[match(indNames(gl_outliers), popmap$sample)]
genI_out <- gl2gi(gl_outliers)
pops_gi_out <- genI_out$pop
strata(genI_out) <- data.frame(pops = pops_gi_out)

p.amova_out      <- poppr.amova(genI_out, ~pops)
amova.pvalues_out <- ade4::randtest(p.amova_out, nrepet = 9999)
overall_fst_out  <- p.amova_out$statphi$Phi[length(p.amova_out$statphi$Phi)]
# updated 10/7/26
print(amova.pvalues_out)
overall_pval_out <- amova.pvalues_out$pvalue[length(amova.pvalues_out$pvalue)]

# outlier amova fst and pval
round(overall_fst_out, 4) # fst = 0.3249
round(overall_pval_out, 4) # 1e-04

# pairwise FST excluding BYC (n=1)
pop_sizes_out  <- table(pop(gl_outliers))
valid_pops_out <- names(pop_sizes_out[pop_sizes_out >= 2])
gl_valid_out   <- gl_outliers[pop(gl_outliers) %in% valid_pops_out, ]

fst_result_out <- gl.fst.pop(gl_valid_out, nboots=1000, percent=95, nclusters=1)
fst_matrix_out <- fst_result_out$Fsts
fst_matrix_out[is.na(fst_matrix_out)] <- t(fst_matrix_out)[is.na(fst_matrix_out)]

# updated 10/7/26
stopifnot(nLoc(gl_bir) == nrow(vcf_bir@fix))
neutral <- setdiff(seq_len(nLoc(gl_bir)), outliers)
gl_neutral <- gl_bir[, neutral]
nLoc(gl_neutral) #12041
genI_neu <- gl2gi(gl_neutral)
strata(genI_neu) <- data.frame(pops = genI_neu$pop)
p.amova_neu <- poppr.amova(genI_neu, ~pops)
amova.pvalues_neu <- ade4::randtest(p.amova_neu, nrepet = 9999)
print(amova.pvalues_neu)
overall_fst_neu <- p.amova_neu$statphi$Phi[length(p.amova_neu$statphi$Phi)]
overall_pval_neu <- amova.pvalues_neu$pvalue[length(amova.pvalues_neu$pvalue)]
round(overall_fst_neu, 4)
round(overall_pval_neu, 4)

# updated 10/7/26
pval_matrix_out <- fst_result_out$Pvalues
pval_matrix_out[is.na(pval_matrix_out)] <- t(pval_matrix_out)[is.na(pval_matrix_out)]
print("pairwise fst outlier snps")
print(round(fst_matrix_out, 4))
print(pval_matrix_out)

stars <- function(p) ifelse(is.na(p), "",
                            ifelse(p < 0.001, "***",
                                   ifelse(p < 0.01,  "**",
                                          ifelse(p < 0.05,  "*", ""))))
cat("  All-loci FST =", round(overall_fst, 4),  stars(overall_pval),     "\n")
cat("  Outlier FST  =", round(overall_fst_out,4), stars(overall_pval_out), "\n")

# updated 10/7/26
p_label <- ifelse(overall_pval < 0.001, "< 0.001", paste("=", round(overall_pval, 3)))

get_upper_tri <- function(mat) { mat[lower.tri(mat)] <- NA; return(mat) }
fst_matrix_labeled <- fst_matrix
rownames(fst_matrix_labeled) <- pop_names_map[rownames(fst_matrix_labeled)]
colnames(fst_matrix_labeled) <- pop_names_map[colnames(fst_matrix_labeled)]
pval_matrix_labeled <- pval_matrix
rownames(pval_matrix_labeled) <- pop_names_map[rownames(pval_matrix)]
colnames(pval_matrix_labeled) <- pop_names_map[colnames(pval_matrix)]
fst_melted <- melt(get_upper_tri(fst_matrix_labeled), na.rm=TRUE)
pval_melted <- melt(get_upper_tri(pval_matrix_labeled), na.rm=TRUE)
fst_melted <- fst_melted %>% left_join(pval_melted, by=c("Var1","Var2"), suffix=c("","_pval")) %>%
  mutate(sig = case_when(is.na(value_pval) ~ "", value_pval < 0.001 ~ "***",
                         value_pval < 0.01 ~ "**", value_pval < 0.05 ~ "*", TRUE ~ ""),
         label = paste0(round(value, 4), sig))

fst_heatmap <- ggplot(fst_melted, aes(Var1, Var2, fill = value)) +
  geom_tile(color="white", linewidth=1) + geom_text(aes(label=label), size=5) +
  scale_fill_gradient(low="white", high="#1D8981", name=expression(italic("F")[ST]),
                      limits=c(0, max(fst_melted$value, na.rm=TRUE) + 0.001)) +
  coord_fixed() +
  labs(x="", y="", subtitle=bquote(atop("AMOVA" ~ italic("F")[ST] ~ "=" ~ .(round(overall_fst,4)) ~ "; p =" ~ .(p_label), "* p < 0.05, ** p < 0.01, *** p < 0.001"))) +
  theme_minimal() +
  theme(axis.text.x=element_text(angle=45, hjust=1, size=12, face="bold"),
        axis.text.y=element_text(size=12, face="bold"),
        plot.subtitle=element_text(size=11, hjust=0.5, face="bold"),
        panel.grid=element_blank(), legend.position="right")
ggsave("/Users/madelineeppley/Desktop/manta26pub/birostris_Fst_heatmap.png",
       fst_heatmap, width=8, height=6, dpi=600)

# everything pulled together into one fig
layer1 <- range_map + theme(legend.position="right", plot.margin=margin(5,5,5,5))
layer2 <- wrap_plots(pca_plot_all, pca_plot_outliers, ncol=2)
layer3 <- wrap_plots(structure_all, structure_outlier, pop_label_plot, ncol=1, heights=c(1,1,0.3))
layer4 <- wrap_plots(fst_heatmap + theme(legend.position="right"), dapc_plot, ncol=2, widths=c(0.8, 1.2))

final_figure <- wrap_plots(layer1, layer2, layer3, layer4,
                           ncol=1, heights=c(1.1, 1.1, 0.8, 1.0)) +
  plot_annotation(
    title = "Mobula birostris",
    subtitle = paste0("All SNPs: n = ", format(n_snps_all, big.mark=","),
                      "   |   Outlier SNPs: n = ", format(n_outliers, big.mark=",")),
    tag_levels = list(c('A','B','C','D','E','','F','G')),
    theme = theme(
      plot.title = element_text(size=18, face="bold.italic", hjust=0),
      plot.subtitle = element_text(size=14, hjust=0),
      plot.tag = element_text(size=16, face="bold")))

ggsave("/Users/madelineeppley/Desktop/manta26pub/birostris_FINALFIGURE.png",
       final_figure, width=14, height=18, dpi=600)
