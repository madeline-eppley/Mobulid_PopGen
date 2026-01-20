##### munkiana data ######
#### M.Eppley 01/19/26 ####
## this script covers VCF -> gl, PCA, LEA ancestry, Fst, outliers, pop structure figures ##

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

# read in filtered VCF (7188 SNPs, 13 individuals)
vcf_mun <- read.vcfR("/Users/madelineeppley/Desktop/manta/finalvcfs/munkiana_filtered_7188snps.vcf")

# population map updated for filtered individuals
popmap <- data.frame(
  sample = c("mu_3", "mu_4", "mu_7", "mu_11", "mu_12", "mu_21",
             "NIC_2728_01", "NIC_2740_01",
             "PER_DCLCW87_MM", "PER_DCLCW130_MM", "PER_DPPW019_MM", 
             "PER_DPPW28_MM", "PER_DPPW98_MM"),
  pop = c(rep("MEX", 6), rep("NIC", 2), rep("PER", 5)))

vcf_mun  #13 samples, 7188 SNPs
table(popmap$pop)  #MEX=6, NIC=2, PER=5

# convert to genlight
gl_mun <- vcfR2genlight(vcf_mun)

# initial data exploration - make a PCA
pca_mun <- glPca(gl_mun, nf=3)

pca_df_mun <- data.frame(
  PC1 = pca_mun$scores[,1],
  PC2 = pca_mun$scores[,2],
  Sample = indNames(gl_mun))

# plot 
ggplot(pca_df_mun, aes(x=PC1, y=PC2)) +
  geom_point(size=2, color="blue") +
  geom_text_repel(aes(label=Sample),
                  size=2.5,
                  max.overlaps=50,
                  segment.size=0.3) +
  theme_minimal() +
  labs(title="munkiana (filtered: 13 ind, 7188 SNPs)",
       x=paste0("PC1 (", round(pca_mun$eig[1]/sum(pca_mun$eig)*100,1), "%)"),
       y=paste0("PC2 (", round(pca_mun$eig[2]/sum(pca_mun$eig)*100,1), "%)"))

# check missingness per ind
gt_mun <- extract.gt(vcf_mun, element = "GT")
ind_missing_mun <- apply(gt_mun, 2, function(x) mean(is.na(x) | x == "./."))
snp_missing_mun <- apply(gt_mun, 1, function(x) mean(is.na(x) | x == "./."))

summary(ind_missing_mun)
summary(snp_missing_mun)

ggplot(data.frame(missing = ind_missing_mun), aes(x = missing)) +
  geom_histogram(binwidth = 0.05, fill = "steelblue", color = "black") +
  theme_minimal() +
  labs(
    title = "Individual Missingness (Filtered Dataset)",
    x = "Proportion of missing genotypes",
    y = "Number of individuals")

ggplot(data.frame(missing = snp_missing_mun), aes(x = missing)) +
  geom_histogram(binwidth = 0.05, fill = "darkorange", color = "black") +
  theme_minimal() +
  labs(
    title = "SNP Missingness (Filtered Dataset)",
    x = "Proportion of missing genotypes",
    y = "Number of SNPs")

#############
# RELATEDNESS
#############
vcf.fn_mun <- "/Users/madelineeppley/Desktop/manta/finalvcfs/munkiana_filtered_7188snps.vcf"
gds.fn_mun <- "/Users/madelineeppley/Desktop/manta/munkiana/mun_filtered.gds"  # NEW filename

# close any open gds
showfile.gds(closeall = TRUE)

snpgdsVCF2GDS(vcf.fn_mun, gds.fn_mun, method="biallelic.only")
snpgdsSummary(gds.fn_mun)

genofile_mun <- snpgdsOpen(gds.fn_mun)
rel_mun <- snpgdsIBDKING(genofile_mun, autosome.only=FALSE)
rel_mun$kinship[1:5, 1:5]

pheatmap(rel_mun$kinship,
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         main = "Pairwise Relatedness (KING) - munkiana filtered")

kin_mat_mun <- rel_mun$kinship
sample_names_mun <- rel_mun$sample.id

rel_df_mun <- as.data.frame(as.table(kin_mat_mun))
colnames(rel_df_mun) <- c("ind1", "ind2", "kinship")

# remove self comparisons
rel_df_mun <- rel_df_mun[rel_df_mun$ind1 != rel_df_mun$ind2, ]
rel_df_mun <- rel_df_mun[!duplicated(t(apply(rel_df_mun[,1:2], 1, sort))), ]

# sort by kinship
rel_df_mun <- rel_df_mun[order(-rel_df_mun$kinship), ]

# top related pairs
head(rel_df_mun, 10)

# add population info to relatedness
rel_df_mun <- rel_df_mun %>%
  mutate(ind1 = as.character(ind1), ind2 = as.character(ind2)) %>%
  left_join(popmap, by = c("ind1" = "sample")) %>%
  rename(pop1 = pop) %>%
  left_join(popmap, by = c("ind2" = "sample")) %>%
  rename(pop2 = pop)

head(rel_df_mun, 10)

# annotated heatmap
rownames(kin_mat_mun) <- sample_names_mun
colnames(kin_mat_mun) <- sample_names_mun

annotation_df <- data.frame(
  Population = popmap$pop[match(sample_names_mun, popmap$sample)])
rownames(annotation_df) <- sample_names_mun

pop_colors <- list(
  Population = c("MEX" = "#6DBAA4", "NIC" = "#ED9A6C", "PER" = "#8C9FCB"))

relate <- pheatmap(
  kin_mat_mun,
  labels_row = sample_names_mun,
  labels_col = sample_names_mun,
  annotation_row = annotation_df,
  annotation_col = annotation_df,
  annotation_colors = pop_colors,
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  main = "Pairwise Genomic Relatedness, M. munkiana")

ggsave("/Users/madelineeppley/Desktop/manta/munkiana/munkiana_relatedness_filtered.png", 
       relate, width = 10, height = 8, dpi = 600)

#############
# SNMF ancestry
#############
vcf2geno("/Users/madelineeppley/Desktop/manta/finalvcfs/munkiana_filtered_7188snps.vcf", 
         "/Users/madelineeppley/Desktop/manta/munkiana/munkiana_filtered.geno")  # NEW filename

# run
project <- snmf("/Users/madelineeppley/Desktop/manta/munkiana/munkiana_filtered.geno", 
                K = 1:5, entropy = TRUE, repetitions = 10, project = "new")

print(summary(project))
plot(project, col = "blue", pch = 19, cex = 1.2)

#############
# PCADAPT
#############
vcf_mun_path <- "/Users/madelineeppley/Desktop/manta/finalvcfs/munkiana_filtered_7188snps.vcf"
geno_mun <- read.pcadapt(vcf_mun_path, type = "vcf")

obj <- pcadapt(geno_mun, K = 1)
plot(obj, option = "manhattan")

pvals <- obj$pvalues
alpha <- 0.01  
outliers <- which(pvals < alpha)
length(outliers) #211 outliers
head(outliers)

out_meta <- vcf_mun@fix[outliers, c("CHROM", "POS", "ID")]
out_meta <- as.data.frame(out_meta, stringsAsFactors = FALSE)
out_meta$LocusID <- paste0(out_meta$CHROM, ":", out_meta$POS)

write.csv(out_meta, "/Users/madelineeppley/Desktop/manta/munkiana/munkiana_pcadapt_outlier_snps_meta_filtered.csv", 
          row.names = FALSE)
writeLines(out_meta$LocusID, "/Users/madelineeppley/Desktop/manta/munkiana/outlier_locus_ids_filtered.txt")

#############
# improving visualizations
#############

#############
## PCA
#############

pop(gl_mun) <- popmap$pop[match(indNames(gl_mun), popmap$sample)]

# pop labels
pop_labels_sep <- data.frame(
  pop = c("MEX", "NIC", "PER"),
  pop_label = c("Mexico", "Nicaragua", "Peru"))

pop_colors_sep <- c("MEX" = "#6DBAA4",
                    "NIC" = "#ED9A6C",                    
                    "PER" = "#8C9FCB")

# run pca, keep first 10 pcs
pca <- glPca(gl_mun, nf = 10)
eig <- pca$eig / sum(pca$eig) * 100

pca_df_sep <- data.frame(
  PC1 = pca$scores[,1],
  PC2 = pca$scores[,2],
  PC3 = pca$scores[,3],
  sample = indNames(gl_mun),
  pop = pop(gl_mun)) %>%
  left_join(pop_labels_sep, by = "pop")

# centroids for plotting
centroid_sep <- pca_df_sep %>%
  group_by(pop) %>%
  summarise(
    PC1_cen = mean(PC1),
    PC2_cen = mean(PC2),
    PC3_cen = mean(PC3),
    n = n())

pca_df_sep <- pca_df_sep %>%
  left_join(centroid_sep, by = "pop")

xlabel <- paste0("PC1 (", round(eig[1], 1), "%)")
ylabel <- paste0("PC2 (", round(eig[2], 1), "%)")

pca_plot_separate <- ggplot(data = pca_df_sep, aes(x = PC1, y = PC2)) +
  stat_ellipse(data = pca_df_sep %>% filter(pop != "NIC"),
               aes(color = pop),               
               type = "norm", 
               level = 0.95,
               linewidth = 1.2,
               show.legend = FALSE) +
  geom_segment(aes(xend = PC1_cen, yend = PC2_cen, color = pop),               
               linewidth = 0.5,
               alpha = 0.6,
               show.legend = FALSE) +
  geom_point(aes(fill = pop, color = pop),           
             size = 4,
             alpha = 0.8,
             shape = 21,
             stroke = 0.5) +
  scale_fill_manual(values = pop_colors_sep,                     
                    name = "Population",
                    breaks = c("MEX", "NIC", "PER"),
                    labels = c("Mexico", "Nicaragua", "Peru")) +
  scale_color_manual(values = pop_colors_sep,                       
                     name = "Population",
                     breaks = c("MEX", "NIC", "PER"),
                     labels = c("Mexico", "Nicaragua", "Peru")) +
  labs(x = xlabel, y = ylabel,
       title = "Mobula munkiana") +
  theme_classic() +
  theme(
    legend.position = "bottom",
    legend.direction = "horizontal",
    axis.text = element_text(size = 11, color = "black"),
    axis.title = element_text(size = 12, face = "bold"),
    plot.title = element_text(size = 13, face = "italic", hjust = 0.5),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    panel.grid.major = element_line(color = "grey90", linewidth = 0.3),
    legend.title = element_blank(),
    legend.text = element_text(size = 11)) +
  guides(
    fill = guide_legend(
      override.aes = list(
        shape = 21,
        size = 5,
        alpha = 0.8)))

pca_plot_separate

ggsave("/Users/madelineeppley/Desktop/manta/munkiana/munkiana_PCA_filtered.png", 
       pca_plot_separate, width = 8, height = 6, dpi = 600)

######### DAPC ########

pop(gl_mun) <- popmap$pop[match(indNames(gl_mun), popmap$sample)]

pop_colors <- c("MEX" = "#6DBAA4",
                "NIC" = "#ED9A6C",                
                "PER" = "#8C9FCB")

# run the dapc analysis
# 3 populations, so n.da = 2
dapc_result <- dapc(gl_mun,
                    pop = pop(gl_mun),
                    n.pca = 4,
                    n.da = 2)

# df for plotting DAPC
dapc_df <- data.frame(
  LD1 = dapc_result$ind.coord[,1],
  LD2 = dapc_result$ind.coord[,2],
  sample = indNames(gl_mun),
  pop = pop(gl_mun))

# centroids
centroid_dapc <- dapc_df %>%
  group_by(pop) %>%
  summarise(
    LD1_cen = mean(LD1),
    LD2_cen = mean(LD2),
    n = n())

dapc_df <- dapc_df %>%
  left_join(centroid_dapc, by = "pop")

# var exp
percent <- dapc_result$eig / sum(dapc_result$eig) * 100

dapc_plot <- ggplot(data = dapc_df, aes(x = LD1, y = LD2)) +
  stat_ellipse(data = dapc_df %>% filter(pop != "NIC"),
               aes(color = pop),               
               type = "norm", 
               level = 0.95,
               linewidth = 1.2,
               show.legend = FALSE) +
  geom_segment(aes(xend = LD1_cen, yend = LD2_cen, color = pop),               
               linewidth = 0.5,
               alpha = 0.6,
               show.legend = FALSE) +
  geom_point(aes(fill = pop, color = pop),           
             size = 4,
             alpha = 0.8,
             shape = 21,
             stroke = 0.5) +
  scale_fill_manual(values = pop_colors,                     
                    name = "Population",
                    breaks = c("MEX", "NIC", "PER"),
                    labels = c("Mexico", "Nicaragua", "Peru")) +
  scale_color_manual(values = pop_colors,                       
                     name = "Population",
                     breaks = c("MEX", "NIC", "PER"),
                     labels = c("Mexico", "Nicaragua", "Peru")) +
  labs(
    x = paste0("DAPC1 (", round(percent[1], 1), "%)"),
    y = paste0("DAPC2 (", round(percent[2], 1), "%)"),
    title = "Mobula munkiana") +
  theme_classic() +
  theme(
    legend.position = "bottom",
    legend.direction = "horizontal",
    axis.text = element_text(size = 11, color = "black"),
    axis.title = element_text(size = 12, face = "bold"),
    plot.title = element_text(size = 13, face = "italic", hjust = 0.5),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    panel.grid.major = element_line(color = "grey90", linewidth = 0.3),
    legend.title = element_blank(),
    legend.text = element_text(size = 11)) +
  guides(
    fill = guide_legend(
      override.aes = list(shape = 21, size = 5, alpha = 0.8)))

dapc_plot

ggsave("/Users/madelineeppley/Desktop/manta/munkiana/munkiana_DAPC_filtered.png", 
       dapc_plot, width = 8, height = 6, dpi = 600)

############### Fst heatmaps ###############

pop(gl_mun) <- popmap$pop[match(indNames(gl_mun), popmap$sample)]

fst_result <- gl.fst.pop(gl_mun, nboots = 1000, percent = 95, nclusters = 1)
fst_matrix <- fst_result$Fsts
fst_matrix[is.na(fst_matrix)] <- t(fst_matrix)[is.na(fst_matrix)]

get_upper_tri <- function(mat) {
  mat[lower.tri(mat)] <- NA
  return(mat)
}

upper_tri <- get_upper_tri(fst_matrix)
fst_melted <- melt(upper_tri, na.rm = TRUE)

fst_heatmap <- ggplot(fst_melted, aes(Var1, Var2, fill = value)) +
  geom_tile(color = "white", linewidth = 1) +
  scale_fill_gradient(
    low = "#F0F9F7",
    high = "#1D8981",
    name = expression(italic("F")[ST]),
    limits = c(0, max(fst_melted$value, na.rm = TRUE))) +
  coord_fixed() +
  labs(x = "", y = "", title = expression(italic("Mobula munkiana"))) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5, size = 11, face = "bold"),
    axis.text.y = element_text(size = 11, face = "bold"),
    plot.title = element_text(size = 13, face = "italic", hjust = 0.5),
    panel.grid = element_blank(),
    legend.position = "right")

fst_heatmap

ggsave("/Users/madelineeppley/Desktop/manta/munkiana/munkiana_Fstheatmap_filtered.png", 
       fst_heatmap, width = 8, height = 6, dpi = 600)

############## AMOVA ##############

pop(gl_mun) <- popmap$pop[match(indNames(gl_mun), popmap$sample)]

genI <- gl2gi(gl_mun)
pops <- genI$pop
strata(genI) <- data.frame(pops)

p.amova <- poppr.amova(genI, ~pops)
amova.pvalues <- ade4::randtest(p.amova, nrepet = 999)

overall_fst <- p.amova$statphi$Phi[3]
overall_pval <- amova.pvalues$pvalue

####################### plotting K 1-5 ########################

project <- load.snmfProject("/Users/madelineeppley/Desktop/manta/munkiana/munkiana_filtered.snmfProject")

sample_names <- indNames(gl_mun)
sample_pops <- pop(gl_mun)

pop_names <- c("MEX" = "Mexico", "NIC" = "Nicaragua", "PER" = "Peru")
pop_order <- c("MEX", "NIC", "PER")

# sort samples by pop
pop_label_df <- data.frame(
  Sample = sample_names,
  Pop = factor(sample_pops, levels = pop_order)) %>%
  arrange(Pop) %>%
  mutate(x_pos = 1:n())

pop_boundaries <- pop_label_df %>%
  group_by(Pop) %>%
  summarise(x_end = max(x_pos)) %>%
  ungroup() %>%
  filter(x_end != max(pop_label_df$x_pos))

pop_spans <- pop_label_df %>%
  group_by(Pop) %>%
  summarise(
    x_start = min(x_pos),
    x_end = max(x_pos),
    x_mid = (min(x_pos) + max(x_pos)) / 2,
    n = n()) %>%
  ungroup() %>%
  mutate(full_name = pop_names[as.character(Pop)])

max_x <- max(pop_label_df$x_pos)

structure_plots <- list()

for (k in 1:5) {
  ce <- cross.entropy(project, K = k)
  best_run <- which.min(ce)
  qmat <- Q(project, K = k, run = best_run)
  
  qdf <- as.data.frame(qmat)
  colnames(qdf) <- paste0("Cluster", 1:k)
  qdf$Sample <- sample_names
  qdf$Pop <- factor(sample_pops, levels = pop_order)
  
  qdf <- qdf %>%
    arrange(Pop) %>%
    mutate(
      Sample_ID = factor(Sample, levels = unique(Sample)),
      x_pos = 1:n())
  
  qdf_long <- qdf %>%
    pivot_longer(cols = starts_with("Cluster"),
                 names_to = "Cluster",
                 values_to = "Proportion")
  
  p <- ggplot(qdf_long, aes(x = x_pos, y = Proportion, fill = Cluster)) +
    geom_bar(stat = "identity", width = 1, color = NA) +
    geom_vline(data = pop_boundaries, aes(xintercept = x_end + 0.5),
               color = "white", linewidth = 0.8, linetype = "solid") +
    scale_fill_brewer(palette = "Set2") +
    scale_y_continuous(expand = c(0, 0),                       
                       breaks = c(0, 0.5, 1.0),
                       labels = c("0.0", "0.5", "1.0"),
                       position = "left") +
    scale_x_continuous(expand = c(0, 0),
                       limits = c(0.5, max_x + 3)) +
    annotate("text", x = max_x + 2, y = 0.5,
             label = paste0("K = ", k), size = 5, fontface = "bold") +
    coord_cartesian(clip = "off") +
    theme_minimal() +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      axis.text.y = element_text(size = 10),
      legend.position = "none",
      panel.grid = element_blank(),
      panel.border = element_blank(),
      plot.margin = margin(5, 30, 2, 5))
  
  structure_plots[[k]] <- p
}

pop_label_plot <- ggplot() +
  geom_text(data = pop_spans,             
            aes(x = x_mid, y = 0.5, label = full_name),
            size = 4, fontface = "bold") +
  geom_vline(data = pop_boundaries,          
             aes(xintercept = x_end + 0.5),
             color = "black", linewidth = 1) +
  scale_x_continuous(limits = c(0.5, max_x + 3),
                     expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
  coord_cartesian(clip = "off") +
  theme_void() +
  theme(plot.margin = margin(0, 30, 10, 5))

combined_structure <- wrap_plots(
  c(structure_plots, list(pop_label_plot)), 
  ncol = 1, 
  heights = c(rep(1, 5), 0.4))

ggsave("/Users/madelineeppley/Desktop/manta/munkiana/munkiana_STRUCTURE_filtered.png", combined_structure,       
       width = 12, height = 11, dpi = 600)

###################### multipanel plot! ######################

# FST heatmap with full names
fst_matrix_labeled <- fst_matrix
rownames(fst_matrix_labeled) <- c("Mexico", "Nicaragua", "Peru")
colnames(fst_matrix_labeled) <- c("Mexico", "Nicaragua", "Peru")

upper_tri_labeled <- get_upper_tri(fst_matrix_labeled)
fst_melted_labeled <- melt(upper_tri_labeled, na.rm = TRUE)

# with amova values 
fst_heatmap_final <- ggplot(fst_melted_labeled, aes(Var1, Var2, fill = value)) +
  geom_tile(color = "white", linewidth = 1) +
  scale_fill_gradient(
    low = "#F0F9F7",
    high = "#1D8981",
    name = expression(italic("F")[ST]),
    limits = c(0, max(fst_melted_labeled$value, na.rm = TRUE))) +
  coord_fixed() +
  labs(x = "", y = "",       
       subtitle = bquote("AMOVA" ~ italic("F")[ST] ~ "=" ~ .(round(overall_fst, 4)) ~ "; p =" ~ .(round(overall_pval, 4)))) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5, size = 11, face = "bold"),
    axis.text.y = element_text(size = 11, face = "bold"),
    plot.subtitle = element_text(size = 10, hjust = 0.5, face = "bold"),
    panel.grid = element_blank(),
    legend.position = "right")

# remove titles
pca_plot_final <- pca_plot_separate + labs(title = NULL)
dapc_plot_final <- dapc_plot + labs(title = NULL)

# top 2 panels will be structure runs K=2 and K=3
top_section <- wrap_plots(
  structure_plots[[2]],
  structure_plots[[3]], 
  pop_label_plot,
  ncol = 1,
  heights = c(1, 1, 0.3))

# middle row will have PCA and DAPC
middle_section <- wrap_plots(
  pca_plot_final,
  dapc_plot_final,
  ncol = 2)

# leave space for the manta picture
manta_space <- ggplot() + theme_void()

bottom_section <- wrap_plots(
  fst_heatmap_final,
  manta_space,
  ncol = 2,
  widths = c(0.8, 2))

final_figure <- wrap_plots(
  top_section,
  middle_section,
  bottom_section,
  ncol = 1,
  heights = c(1.3, 1.5, 0.8)) +
  plot_annotation(
    title = "Mobula munkiana",
    subtitle = "Neutral loci (n = 7,188 SNPs; 13 individuals)",
    tag_levels = list(c('A', 'B', '', 'C', 'D', 'E', '', '')),
    theme = theme(
      plot.title = element_text(size = 14, face = "bold.italic", hjust = 0),
      plot.subtitle = element_text(size = 12, face = "bold", hjust = 0),
      plot.tag = element_text(size = 14, face = "bold")))

final_figure

ggsave("/Users/madelineeppley/Desktop/manta/munkiana/munkiana_FULLFIGURE_filtered.png", final_figure,       
       width = 16, height = 15, dpi = 600)
