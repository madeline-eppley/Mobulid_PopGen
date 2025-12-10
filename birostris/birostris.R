##### birostris data #####
#### author M.Eppley, date: 11/24/25 ####
## this script covers VCF -> gl, PCA, LEA ancestry, and Fst outliers ##

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
vcf_bir <- read.vcfR("/Users/madelineeppley/Desktop/manta/birostris/minDP10_maxmiss0.8_filtInd.recode.vcf")
popmap <- read.table("/Users/madelineeppley/Desktop/manta/birostris/pop_map_birostris", header = FALSE, col.names = c("sample", "pop"))

# convert to genlight
gl_bir <- vcfR2genlight(vcf_bir)

# initial data exploration - make a PCA
pca_bir <- glPca(gl_bir, nf=3)
#scatter(pca_bir)

pca_df_bir <- data.frame(
  PC1 = pca_bir$scores[,1],
  PC2 = pca_bir$scores[,2],
  Sample = indNames(gl_bir))

# plot 
ggplot(pca_df_bir, aes(x=PC1, y=PC2)) +
  geom_point(size=2, color="blue") +
  geom_text_repel(aes(label=Sample),
                  size=2.5,
                  max.overlaps=50,
                  segment.size=0.3) +
  theme_minimal() +
  labs(title="birostris",
       x=paste0("PC1 (", round(pca_bir$eig[1]/sum(pca_bir$eig)*100,1), "%)"),
       y=paste0("PC2 (", round(pca_bir$eig[2]/sum(pca_bir$eig)*100,1), "%)"))

# check missingness per ind
gt_bir <- extract.gt(vcf_bir, element = "GT")
ind_missing_bir <- apply(gt_bir, 2, function(x) mean(is.na(x) | x == "./."))
snp_missing_bir <- apply(gt_bir, 1, function(x) mean(is.na(x) | x == "./."))
summary(ind_missing_bir)
summary(snp_missing_bir)

sort(ind_missing_bir, decreasing = TRUE)[1:22]

ggplot(data.frame(missing = ind_missing_bir), aes(x = missing)) +
  geom_histogram(binwidth = 0.05, fill = "steelblue", color = "black") +
  theme_minimal() +
  labs(
    title = "Individual Missingness",
    x = "Proportion of missing genotypes",
    y = "Number of individuals")

ggplot(data.frame(missing = snp_missing_bir), aes(x = missing)) +
  geom_histogram(binwidth = 0.05, fill = "darkorange", color = "black") +
  theme_minimal() +
  labs(
    title = "SNP Missingness",
    x = "Proportion of missing genotypes",
    y = "Number of SNPs")


# relatedness 
vcf.fn_bir <- "/Users/madelineeppley/Desktop/manta/birostris/minDP10_maxmiss0.8_filtInd.recode.vcf"
gds.fn_bir <- "/Users/madelineeppley/Desktop/manta/bir.gds"

snpgdsVCF2GDS(vcf.fn_bir, gds.fn_bir, method="biallelic.only")  # keep only biallelic SNPs
snpgdsSummary(gds.fn_bir)

genofile_bir <- snpgdsOpen(gds.fn_bir)
rel_bir <- snpgdsIBDKING(genofile_bir, autosome.only=FALSE)
rel_bir$kinship[1:5, 1:5]

pheatmap(rel_bir$kinship,
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         main = "Pairwise Relatedness (KING kinship coefficients) - bir")

kin_mat_bir <- rel_bir$kinship
sample_names_bir <- rel_bir$sample.id

rel_df_bir <- as.data.frame(as.table(kin_mat_bir))
colnames(rel_df_bir) <- c("ind1", "ind2", "kinship")

# remove self comparisions
rel_df_bir <- rel_df_bir[rel_df_bir$ind1 != rel_df_bir$ind2, ]
rel_df_bir <- rel_df_bir[!duplicated(t(apply(rel_df_bir[,1:2], 1, sort))), ]

# sort by kinship
rel_df_bir <- rel_df_bir[order(-rel_df_bir$kinship), ]

# top related pairs
head(rel_df_bir, 10) # highest value is 0.1, great

relate <- pheatmap(
  kin_mat_bir,
  labels_row = sample_names_bir,
  labels_col = sample_names_bir,
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  main = "Pairwise Genomic Relatedness, M. birostris")

ggsave("/Users/madelineeppley/Desktop/manta/birostris/birostris_relatedness.png", relate, width = 8, height = 6, dpi = 600)

## looks great here - now no relatedness between samples 

# snmf - to LEA format
vcf2geno("/Users/madelineeppley/Desktop/manta/birostris/minDP10_maxmiss0.8_filtInd.recode.vcf", "/Users/madelineeppley/Desktop/manta/birostris/birostris.geno")

# run
project <- snmf("/Users/madelineeppley/Desktop/manta/birostris/birostris.geno", K = 1:5, entropy = TRUE, repetitions = 10, project = "new")

# find the best K value and best run
print(summary(project)) # our highest support is K=1
plot(project, col = "blue", pch = 19, cex = 1.2) # highest support K=1

# just check K=2 to see
best_k <- 2
k2_runs <- which(project@K == 2)
best_run <- which.min(cross.entropy(project, K = best_k))
qmatrix <- Q(project, K = best_k, run = best_run)

sample_names <- colnames(vcf_bir@gt)[-1] # need to remove first col here
sample_names

# nope looks kind of artificial to me
barplot(t(qmatrix), col = rainbow(best_k), border = NA, 
        names.arg = sample_names, las = 2, cex.names = 0.8,
        main = paste("birostris", best_k))

# pcadapt
vcf_bir_path <- "/Users/madelineeppley/Desktop/manta/birostris/minDP10_maxmiss0.8_filtInd.recode.vcf"
geno_bir <- read.pcadapt(vcf_bir_path, type = "vcf")
obj <- pcadapt(geno_bir, K = 1)
plot(obj, option = "manhattan")

pvals <- obj$pvalues
alpha <- 0.01  
outliers <- which(pvals < alpha)
length(outliers) #179 outliers
head(outliers)

out_meta <- vcf_bir@fix[outliers, c("CHROM", "POS", "ID")]
out_meta <- as.data.frame(out_meta, stringsAsFactors = FALSE)

# create a LocusID col
out_meta$LocusID <- paste0(out_meta$CHROM, ":", out_meta$POS)

write.csv(out_meta, "/Users/madelineeppley/Desktop/manta/birostris/birostris_pcadapt_outlier_snps_meta.csv", row.names = FALSE)
writeLines(out_meta$LocusID, "/Users/madelineeppley/Desktop/manta/birostris/outlier_locus_ids.txt")


################################
##### improving visualizations
#################################
################################

### PCA
pop_separate <- popmap$pop
pop_separate[popmap$sample == "BYC_RMB_01"] <- "BYC"
pop(gl_bir) <- pop_separate[match(indNames(gl_bir), popmap$sample)]

# pop labels
pop_labels_sep <- data.frame(
  pop = c("PER", "REV", "IND", "BYC"),
  pop_label = c("Peru", "Mexico", "India", "Bycatch"))

pop_colors_sep <- c("IND" = "#6DBAA4",
                    "REV" = "#ED9A6C",  
                    "PER" = "#8C9FCB",
                    "BYC" = "#DA8EC0")  

# run pca, keep first 10 pcs
pca <- glPca(gl_bir, nf = 10)
eig <- pca$eig / sum(pca$eig) * 100

pca_df_sep <- data.frame(
  PC1 = pca$scores[,1],
  PC2 = pca$scores[,2],
  PC3 = pca$scores[,3],
  sample = indNames(gl_bir),
  pop = pop(gl_bir)) %>%
  left_join(pop_labels_sep, by = "pop")

# centriods for plotting
centroid_sep <- pca_df_sep %>%
  group_by(pop) %>%
  summarise(
    PC1_cen = mean(PC1),
    PC2_cen = mean(PC2),
    PC3_cen = mean(PC3))

pca_df_sep <- pca_df_sep %>%
  left_join(centroid_sep, by = "pop")

xlabel <- paste0("PC1 (", round(eig[1], 1), "%)")
ylabel <- paste0("PC2 (", round(eig[2], 1), "%)")

pca_plot_separate <- ggplot(data = pca_df_sep, aes(x = PC1, y = PC2)) +
  stat_ellipse(aes(color = pop), 
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
                    breaks = c("IND", "PER", "REV", "BYC"),
                    labels = c("India", "Peru", "Mexico", "Bycatch")) +
  scale_color_manual(values = pop_colors_sep, 
                     name = "Population",
                     breaks = c("IND", "PER", "REV", "BYC"),
                     labels = c("India", "Peru", "Mexico", "Bycatch")) +
  labs(x = xlabel, y = ylabel,
       title = "Mobula birostris") +
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

ggsave("/Users/madelineeppley/Desktop/manta/birostris/birostris_PCA.png", pca_plot_separate, width = 8, height = 6, dpi = 600)

#######
# DAPC
#######
pop(gl_bir) <- popmap$pop[match(indNames(gl_bir), popmap$sample)]

# set a nice color palette :o
pop_colors <- c("IND" = "#6DBAA4",
                "REV" = "#ED9A6C",
                "PER" = "#8C9FCB",
                "BYC" = "#DA8EC0")

set.seed(999)
grp <- find.clusters(gl_bir,
                     scale = FALSE,
                     n.pca = 10,      
                     choose.n.clust = TRUE,
                     max.n.clust = 10) # optimal 1 cluster

grp

#plot(grp$Kstat, type="b", xlab="K", ylab="BIC")

# run the dapc analysis, keep 4 PCs --- revisited this and decreased from original 8
dapc_result <- dapc(gl_bir,
                    pop = pop(gl_bir),
                    n.pca = 4, # based on sampling locations, K is set to 4 (IND, REV, PER, BYC)
                    n.da = 3) # one minus n PCA retained

# df for plotting DAPC
dapc_df <- data.frame(
  LD1 = dapc_result$ind.coord[,1],
  LD2 = dapc_result$ind.coord[,2],
  sample = indNames(gl_bir),
  pop = pop(gl_bir),
  is_bycatch = indNames(gl_bir) == "BYC_RMB_01",
  display_pop = ifelse(indNames(gl_bir) == "BYC_RMB_01", "BYC", as.character(pop(gl_bir))))

# centriods - not for bycatch since there's only 1 sample 
centroid_dapc <- dapc_df %>%
  filter(!is_bycatch) %>%
  group_by(pop) %>%
  summarise(
    LD1_cen = mean(LD1),
    LD2_cen = mean(LD2))

dapc_df <- dapc_df %>%
  left_join(centroid_dapc, by = "pop")

# var exp
percent <- dapc_result$eig / sum(dapc_result$eig) * 100

# DAPC plot!~!
dapc_plot <- ggplot(data = dapc_df, aes(x = LD1, y = LD2)) +
  stat_ellipse(data = filter(dapc_df, !is_bycatch), # bycatch just 1 sample, no ellipses
               aes(color = pop), 
               type = "norm",
               level = 0.95,
               linewidth = 1.2,
               show.legend = FALSE) +
  geom_segment(data = filter(dapc_df, !is_bycatch), # bycatch just has 1 sample here, so we will not get a centriod
               aes(xend = LD1_cen, yend = LD2_cen, color = pop), 
               linewidth = 0.5,
               alpha = 0.6,
               show.legend = FALSE) +
  geom_point(aes(fill = display_pop, color = display_pop), 
             size = 4, 
             alpha = 0.8,
             shape = 21,
             stroke = 0.5) +
  scale_fill_manual(values = pop_colors, 
                    name = "Population",
                    breaks = c("IND", "PER", "REV", "BYC"),
                    labels = c("India", "Peru", "Mexico", "Bycatch")) +
  scale_color_manual(values = pop_colors, 
                     name = "Population",
                     breaks = c("IND", "PER", "REV", "BYC"),
                     labels = c("India", "Peru", "Mexico", "Bycatch")) +
  labs(
    x = paste0("DAPC1 (", round(percent[1], 1), "%)"),
    y = paste0("DAPC2 (", round(percent[2], 1), "%)"),
    title = "Mobula birostris") +
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


ggsave("/Users/madelineeppley/Desktop/manta/birostris/birostris_DAPC.png", dapc_plot, width = 8, height = 6, dpi = 600)

###############
# Fst heatmaps
###############
pop(gl_bir) <- popmap$pop[match(indNames(gl_bir), popmap$sample)]

# calculate pairwise Fst
# weir & cockerham Fst, 1000 bootstrap replicates, 95% confidence intervals
fst_result <- gl.fst.pop(gl_bir, nboots = 1000, percent = 95, nclusters = 1)
fst_matrix <- fst_result$Fsts

# lower plot
fst_matrix[is.na(fst_matrix)] <- t(fst_matrix)[is.na(fst_matrix)]

# upper plot
get_upper_tri <- function(mat) {
  mat[lower.tri(mat)] <- NA
  return(mat)}

upper_tri <- get_upper_tri(fst_matrix)

# melt for ggplot
fst_melted <- melt(upper_tri, na.rm = TRUE)

# heatmap!!
fst_heatmap <- ggplot(fst_melted, aes(Var1, Var2, fill = value)) +
  geom_tile(color = "white", linewidth = 1) +
  scale_fill_gradient(
    low = "#F0F9F7",
    high = "#1D8981",
    name = expression(italic("F")[ST]),
    limits = c(0, max(fst_melted$value, na.rm = TRUE))) +
  coord_fixed() +
  labs(x = "", y = "", title = expression(italic("Mobula birostris"))) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5, size = 11, face = "bold"),
    axis.text.y = element_text(size = 11, face = "bold"),
    plot.title = element_text(size = 13, face = "italic", hjust = 0.5),
    panel.grid = element_blank(),
    legend.position = "right")

fst_heatmap
ggsave("/Users/madelineeppley/Desktop/manta/birostris/birostris_Fstheatmap.png", fst_heatmap, width = 8, height = 6, dpi = 600)

############
### AMOVA
############
# pop assignment check again
pop(gl_bir) <- popmap$pop[match(indNames(gl_bir), popmap$sample)]

# genid
genI <- gl2gi(gl_bir)
pops <- genI$pop
strata(genI) <- data.frame(pops)

# AMOVA
p.amova <- poppr.amova(genI, ~pops)

# p vals
amova.pvalues <- ade4::randtest(p.amova, nrepet = 999)

# overall fst  = Phi statistic and return both that and p-value
overall_fst <- p.amova$statphi$Phi[3]
overall_pval <- amova.pvalues$pvalue
overall_fst # 0.004515259
overall_pval
# AMOVA FST = 0.0045; p = 0.001


#####################
## plotting K 1-5 ###
####################
### load in SNMF
project <- load.snmfProject("/Users/madelineeppley/Desktop/manta/birostris/birostris.snmfProject")

# sample names again to organize
sample_names <- indNames(gl_bir)
sample_pops <- pop(gl_bir)

# pop names replace with readable
pop_names <- c("IND" = "India", "PER" = "Peru", "REV" = "Mexico")

pop_label_df <- data.frame(
  Sample = sample_names,
  Pop = sample_pops) %>%
  arrange(Pop, Sample) %>%
  mutate(x_pos = 1:n())

# we need to figure out pop boundaries and max x values
pop_boundaries <- pop_label_df %>%
  group_by(Pop) %>%
  summarise(x_end = max(x_pos)) %>%
  filter(x_end != max(pop_label_df$x_pos))

pop_spans <- pop_label_df %>%
  group_by(Pop) %>%
  summarise(
    x_start = min(x_pos),
    x_end = max(x_pos),
    x_mid = mean(x_pos),
    full_name = pop_names[Pop])

max_x <- max(pop_label_df$x_pos)

# let's make stacked barplots of K 1-5
structure_plots <- list()

for (k in 1:5) {
  ce <- cross.entropy(project, K = k)
  best_run <- which.min(ce)
  qmat <- Q(project, K = k, run = best_run)
  qdf <- as.data.frame(qmat)
  colnames(qdf) <- paste0("Cluster", 1:k)
  qdf$Sample <- sample_names
  qdf$Pop <- sample_pops
  
  qdf <- qdf %>% # change the format of the q dataframe
    arrange(Pop, Sample) %>%
    mutate(
      Sample_ID = factor(Sample, levels = Sample),
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

# we need some white lines to separate the populations, add lines
pop_label_plot <- ggplot() +
  geom_text(data = pop_spans, 
            aes(x = x_mid, y = 0.5, label = full_name),
            size = 5, fontface = "bold") +
  geom_vline(data = pop_boundaries,
             aes(xintercept = x_end + 0.5),
             color = "black", linewidth = 1) +
  scale_x_continuous(limits = c(0.5, max_x + 3),
                     expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
  coord_cartesian(clip = "off") +
  theme_void() +
  theme(plot.margin = margin(0, 30, 10, 5))

# plot!
combined_structure <- wrap_plots(
  c(structure_plots, list(pop_label_plot)), 
  ncol = 1, 
  heights = c(rep(1, 5), 0.4))

ggsave("/Users/madelineeppley/Desktop/manta/birostris/birostris_STRUCTURE.png", combined_structure, 
       width = 12, height = 11, dpi = 600)

###################
## multipanel plot!
####################
# we probably want to do K=2 and K=3 here since they're the most interesting
structure_k2 <- structure_plots[[2]]
structure_k3 <- structure_plots[[3]]

# fix our fst matrix labels to be more readable for composite image
fst_matrix_labeled <- fst_matrix
rownames(fst_matrix_labeled) <- c("India", "Peru", "Mexico")
colnames(fst_matrix_labeled) <- c("India", "Peru", "Mexico")

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
       subtitle = bquote("AMOVA" ~ italic("F")[ST] ~ "= 0.005; p = 0.001")) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5, size = 11, face = "bold"),
    axis.text.y = element_text(size = 11, face = "bold"),
    plot.subtitle = element_text(size = 10, hjust = 0.5, face = "bold"),
    panel.grid = element_blank(),
    legend.position = "right")

# remove the titles from the plots - we don't need species names repeated
pca_plot_final <- pca_plot_separate + labs(title = NULL)

dapc_plot_final <- dapc_plot + labs(title = NULL)

#birostris_img <- readPNG("/Users/madelineeppley/Desktop/manta/birostris/birostris_image.png")
# this didn't work bc of scaling issues, but it's located here to try again if it makes sense later

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

# leave some space for the manta picture
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
    title = "Mobula birostris",
    subtitle = "Neutral loci (n=11,843)",
    tag_levels = list(c('A', 'B', '', 'C', 'D', 'E', '', '')),
    theme = theme(
      plot.title = element_text(size = 14, face = "bold.italic", hjust = 0),
      plot.subtitle = element_text(size = 12, face = "bold", hjust = 0),
      plot.tag = element_text(size = 14, face = "bold")))

final_figure

# I instead saved it as a PDF in US LETTER size because the font text was scaled better, this size was too small
#ggsave("/Users/madelineeppley/Desktop/manta/birostris/birostris_FULLFIGURE.png", final_figure, 
       #width = 16, height = 15, dpi = 600)
