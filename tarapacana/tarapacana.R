##### tarapacana data #########
# M.Eppley, 1/19/26 ##
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

# read in filtered VCF and popmap
vcf_tar <- read.vcfR("/Users/madelineeppley/Desktop/manta/finalvcfs/tarapacana_final.vcf")
popmap <- read.table("/Users/madelineeppley/Desktop/manta/tarapacana/tarapacana_pop_map_3", 
                     header = FALSE, col.names = c("sample", "pop"))

# convert to genlight
gl_tar <- vcfR2genlight(vcf_tar)

# initial data exploration - make a PCA
pca_tar <- glPca(gl_tar, nf=3)
pca_df_tar <- data.frame(
  PC1 = pca_tar$scores[,1],
  PC2 = pca_tar$scores[,2],
  Sample = indNames(gl_tar))

# plot 
ggplot(pca_df_tar, aes(x=PC1, y=PC2)) +
  geom_point(size=2, color="blue") +
  geom_text_repel(aes(label=Sample),
                  size=2.5,
                  max.overlaps=50,
                  segment.size=0.3) +
  theme_minimal() +
  labs(title="tarapacana",
       x=paste0("PC1 (", round(pca_tar$eig[1]/sum(pca_tar$eig)*100,1), "%)"),
       y=paste0("PC2 (", round(pca_tar$eig[2]/sum(pca_tar$eig)*100,1), "%)"))

# check missingness per ind
gt_tar <- extract.gt(vcf_tar, element = "GT")
ind_missing_tar <- apply(gt_tar, 2, function(x) mean(is.na(x) | x == "./."))
snp_missing_tar <- apply(gt_tar, 1, function(x) mean(is.na(x) | x == "./."))

summary(ind_missing_tar)
summary(snp_missing_tar)
sort(ind_missing_tar, decreasing = TRUE)[1:22]

ggplot(data.frame(missing = ind_missing_tar), aes(x = missing)) +
  geom_histogram(binwidth = 0.05, fill = "steelblue", color = "black") +
  theme_minimal() +
  labs(
    title = "Individual Missingness",
    x = "Proportion of missing genotypes",
    y = "Number of individuals")

ggplot(data.frame(missing = snp_missing_tar), aes(x = missing)) +
  geom_histogram(binwidth = 0.05, fill = "darkorange", color = "black") +
  theme_minimal() +
  labs(
    title = "SNP Missingness",
    x = "Proportion of missing genotypes",
    y = "Number of SNPs")

# relatedness 
vcf.fn_tar <- "/Users/madelineeppley/Desktop/manta/finalvcfs/tarapacana_final.vcf"
gds.fn_tar <- "/Users/madelineeppley/Desktop/manta/tarapacana/tar.gds"
snpgdsVCF2GDS(vcf.fn_tar, gds.fn_tar, method="biallelic.only")
snpgdsSummary(gds.fn_tar)

genofile_tar <- snpgdsOpen(gds.fn_tar)
rel_tar <- snpgdsIBDKING(genofile_tar, autosome.only=FALSE)
rel_tar$kinship[1:5, 1:5]

pheatmap(rel_tar$kinship,
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         main = "Pairwise Relatedness (KING kinship coefficients) - tar")

kin_mat_tar <- rel_tar$kinship
sample_names_tar <- rel_tar$sample.id

rel_df_tar <- as.data.frame(as.table(kin_mat_tar))
colnames(rel_df_tar) <- c("ind1", "ind2", "kinship")

# remove self comparisons
rel_df_tar <- rel_df_tar[rel_df_tar$ind1 != rel_df_tar$ind2, ]
rel_df_tar <- rel_df_tar[!duplicated(t(apply(rel_df_tar[,1:2], 1, sort))), ]

# sort by kinship
rel_df_tar <- rel_df_tar[order(-rel_df_tar$kinship), ]

# top related pairs
head(rel_df_tar, 10)

relate <- pheatmap(
  kin_mat_tar,
  labels_row = sample_names_tar,
  labels_col = sample_names_tar,
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  main = "Pairwise Genomic Relatedness, M. tarapacana")

ggsave("/Users/madelineeppley/Desktop/manta/tarapacana/tarapacana_relatedness.png", 
       relate, width = 8, height = 6, dpi = 600)

# snmf - to LEA format
vcf2geno("/Users/madelineeppley/Desktop/manta/finalvcfs/tarapacana_final.vcf", 
         "/Users/madelineeppley/Desktop/manta/tarapacana/tarapacana.geno")

# run
project <- snmf("/Users/madelineeppley/Desktop/manta/tarapacana/tarapacana.geno", 
                K = 1:5, entropy = TRUE, repetitions = 10, project = "new")

print(summary(project))
plot(project, col = "blue", pch = 19, cex = 1.2)

# pcadapt
vcf_tar_path <- "/Users/madelineeppley/Desktop/manta/finalvcfs/tarapacana_final.vcf"
geno_tar <- read.pcadapt(vcf_tar_path, type = "vcf")

obj <- pcadapt(geno_tar, K = 1)
plot(obj, option = "manhattan")

pvals <- obj$pvalues
alpha <- 0.01  
outliers <- which(pvals < alpha)
length(outliers) #1157 outliers
head(outliers)

out_meta <- vcf_tar@fix[outliers, c("CHROM", "POS", "ID")]
out_meta <- as.data.frame(out_meta, stringsAsFactors = FALSE)
out_meta$LocusID <- paste0(out_meta$CHROM, ":", out_meta$POS)

write.csv(out_meta, "/Users/madelineeppley/Desktop/manta/tarapacana/tarapacana_pcadapt_outlier_snps_meta.csv", 
          row.names = FALSE)
writeLines(out_meta$LocusID, "/Users/madelineeppley/Desktop/manta/tarapacana/outlier_locus_ids.txt")

#############
# improving visualizations
#############

#############
## PCA
#############
pop(gl_tar) <- popmap$pop[match(indNames(gl_tar), popmap$sample)]

# pop labels
pop_labels_sep <- data.frame(
  pop = c("EAST", "OFF", "CTR", "WST"),
  pop_label = c("Eastern", "Offshore", "Central", "Western"))

pop_colors_sep <- c("EAST" = "#6DBAA4",
                    "OFF" = "#ED9A6C",                    
                    "CTR" = "#8C9FCB",
                    "WST" = "#DA8EC0")

# run pca, keep first 10 pcs
pca <- glPca(gl_tar, nf = 10)
eig <- pca$eig / sum(pca$eig) * 100

pca_df_sep <- data.frame(
  PC1 = pca$scores[,1],
  PC2 = pca$scores[,2],
  PC3 = pca$scores[,3],
  sample = indNames(gl_tar),
  pop = pop(gl_tar)) %>%
  left_join(pop_labels_sep, by = "pop")

# centroids for plotting
centroid_sep <- pca_df_sep %>%
  group_by(pop) %>%
  summarise(
    PC1_cen = mean(PC1),
    PC2_cen = mean(PC2),
    PC3_cen = mean(PC3),
    n = n())  # count samples per pop

pca_df_sep <- pca_df_sep %>%
  left_join(centroid_sep, by = "pop")

pops_with_ellipse <- centroid_sep %>% filter(n >= 3) %>% pull(pop)

xlabel <- paste0("PC1 (", round(eig[1], 1), "%)")
ylabel <- paste0("PC2 (", round(eig[2], 1), "%)")

pca_plot_separate <- ggplot(data = pca_df_sep, aes(x = PC1, y = PC2)) +
  stat_ellipse(data = filter(pca_df_sep, pop %in% pops_with_ellipse),
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
                    breaks = c("EAST", "OFF", "CTR", "WST"),
                    labels = c("Eastern", "Offshore", "Central", "Western")) +
  scale_color_manual(values = pop_colors_sep,                    
                     name = "Population",
                     breaks = c("EAST", "OFF", "CTR", "WST"),
                     labels = c("Eastern", "Offshore", "Central", "Western")) +
  labs(x = xlabel, y = ylabel,
       title = "Mobula tarapacana") +
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

ggsave("/Users/madelineeppley/Desktop/manta/tarapacana/tarapacana_PCA.png", 
       pca_plot_separate, width = 8, height = 6, dpi = 600)

########
# DAPC
########
pop(gl_tar) <- popmap$pop[match(indNames(gl_tar), popmap$sample)]

pop_colors <- c("EAST" = "#6DBAA4",
                "OFF" = "#ED9A6C",                    
                "CTR" = "#8C9FCB",
                "WST" = "#DA8EC0")


# run the dapc analysis
dapc_result <- dapc(gl_tar,
                    pop = pop(gl_tar),
                    n.pca = 4,
                    n.da = 3)

# df for plotting DAPC
dapc_df <- data.frame(
  LD1 = dapc_result$ind.coord[,1],
  LD2 = dapc_result$ind.coord[,2],
  sample = indNames(gl_tar),
  pop = pop(gl_tar))

# centroids with sample counts
centroid_dapc <- dapc_df %>%
  group_by(pop) %>%
  summarise(
    LD1_cen = mean(LD1),
    LD2_cen = mean(LD2),
    n = n())

dapc_df <- dapc_df %>%
  left_join(centroid_dapc, by = "pop")

# identify pops with enough samples for ellipses
pops_with_ellipse_dapc <- centroid_dapc %>% filter(n >= 3) %>% pull(pop)

# var exp
percent <- dapc_result$eig / sum(dapc_result$eig) * 100

dapc_plot <- ggplot(data = dapc_df, aes(x = LD1, y = LD2)) +
  stat_ellipse(data = filter(dapc_df, pop %in% pops_with_ellipse_dapc),
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
                    breaks = c("EAST", "OFF", "CTR", "WST"),
                    labels = c("Eastern", "Offshore", "Central", "Western")) +
  scale_color_manual(values = pop_colors,                    
                     name = "Population",
                     breaks = c("EAST", "OFF", "CTR", "WST"),
                     labels = c("Eastern", "Offshore", "Central", "Western")) +
  labs(
    x = paste0("DAPC1 (", round(percent[1], 1), "%)"),
    y = paste0("DAPC2 (", round(percent[2], 1), "%)"),
    title = "Mobula tarapacana") +
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

ggsave("/Users/madelineeppley/Desktop/manta/tarapacana/tarapacana_DAPC.png", 
       dapc_plot, width = 8, height = 6, dpi = 600)

###############
# Fst heatmaps
###############
pop(gl_tar) <- popmap$pop[match(indNames(gl_tar), popmap$sample)]

fst_result <- gl.fst.pop(gl_tar, nboots = 1000, percent = 95, nclusters = 1)
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
  labs(x = "", y = "", title = expression(italic("Mobula tarapacana"))) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 11, face = "bold"),
    axis.text.y = element_text(size = 11, face = "bold"),
    plot.title = element_text(size = 13, face = "italic", hjust = 0.5),
    panel.grid = element_blank(),
    legend.position = "right")

fst_heatmap

ggsave("/Users/madelineeppley/Desktop/manta/tarapacana/tarapacana_Fstheatmap.png", 
       fst_heatmap, width = 8, height = 6, dpi = 600)

##############
# AMOVA
##############
pop(gl_tar) <- popmap$pop[match(indNames(gl_tar), popmap$sample)]

genI <- gl2gi(gl_tar)
pops <- genI$pop
strata(genI) <- data.frame(pops)

p.amova <- poppr.amova(genI, ~pops)
amova.pvalues <- ade4::randtest(p.amova, nrepet = 999)

overall_fst <- p.amova$statphi$Phi[3]
overall_pval <- amova.pvalues$pvalue
overall_fst
overall_pval

#####################
## plotting K 1-5 ###
#####################
project <- load.snmfProject("/Users/madelineeppley/Desktop/manta/tarapacana/tarapacana.snmfProject")

sample_names <- indNames(gl_tar)
sample_pops <- pop(gl_tar)

pop_names <- c("EAST" = "Eastern", "OFF" = "Offshore", "CTR" = "Central", "WST" = "Western")
pop_order <- c("EAST", "OFF", "CTR", "WST")

# sort by pop
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
  mutate(
    full_name = pop_names[as.character(Pop)],
    # kind of need abbreviations here to avoid overlap
    display_name = case_when(
      n <= 2 ~ as.character(Pop),
      TRUE ~ full_name))

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
            aes(x = x_mid, y = 0.5, label = display_name),
            size = 3, fontface = "bold") +
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

ggsave("/Users/madelineeppley/Desktop/manta/tarapacana/tarapacana_STRUCTURE.png", combined_structure,      
       width = 15, height = 11, dpi = 600)

#####################
# multipanel plot!
#####################

# using abbreviations here too to avoid overlap
fst_matrix_labeled <- fst_matrix
rownames(fst_matrix_labeled) <- c("E", "Off", "C", "W")
colnames(fst_matrix_labeled) <- c("E", "Off", "C", "W")

upper_tri_labeled <- get_upper_tri(fst_matrix_labeled)
fst_melted_labeled <- melt(upper_tri_labeled, na.rm = TRUE)

fst_heatmap_final <- ggplot(fst_melted_labeled, aes(Var1, Var2, fill = value)) +
  geom_tile(color = "white", linewidth = 1) +
  scale_fill_gradient(
    low = "#F0F9F7",
    high = "#1D8981",
    name = expression(italic("F")[ST]),
    limits = c(0, max(fst_melted_labeled$value, na.rm = TRUE))) +
  coord_fixed() +
  labs(x = "", y = "",      
       subtitle = bquote("AMOVA" ~ italic("F")[ST] ~ "= 0.0032; p = 0.008")) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5, size = 11, face = "bold"),
    axis.text.y = element_text(size = 11, face = "bold"),
    plot.subtitle = element_text(size = 10, hjust = 0.5, face = "bold"),
    panel.grid = element_blank(),
    legend.position = "right")

# remove the titles from the plots
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
    title = "Mobula tarapacana",
    subtitle = "Neutral loci (n=43,165)",
    tag_levels = list(c('A', 'B', '', 'C', 'D', 'E', '', '')),
    theme = theme(
      plot.title = element_text(size = 14, face = "bold.italic", hjust = 0),
      plot.subtitle = element_text(size = 12, face = "bold", hjust = 0),
      plot.tag = element_text(size = 14, face = "bold")))

final_figure

ggsave("/Users/madelineeppley/Desktop/manta/tarapacana/tarapacana_FULLFIGURE.png", final_figure,      
       width = 16, height = 15, dpi = 600)
