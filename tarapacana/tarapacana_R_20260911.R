### M. tarapacana analysis to produce publication figure ####
# M. Eppley, v4. updated 10/20/26

library(vcfR); library(pcadapt); library(SNPRelate); library(ggplot2)
library(pheatmap); library(LEA); library(adegenet); library(ggrepel)
library(dplyr); library(tidyr); library(RColorBrewer); library(dartR)
library(reshape2); library(poppr); library(ade4); library(patchwork)
library(grid); library(gridExtra); library(maps); library(sf)
library(qvalue) # updated 10/20/26 - qvalue() is used below but was never loaded in this script

outdir <- "/Users/madelineeppley/Desktop/manta26pub"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
setwd(outdir)  # so pcadapt can write positions.txt here

pop_colors <- c("EAST" = "#6DBAA4", "GAL" = "#ED9A6C", "CTR" = "#8C9FCB", "WST" = "#DA8EC0")
structure_colors <- c("#264653", "#2a9d8f", "#8ab17d", "#e9c46a", "#f4a261", "#e76f51")

vcf_tar <- read.vcfR("/Users/madelineeppley/Desktop/manta/finalvcfs/tarapacana_final.vcf")
popmap <- data.frame(
  sample = c("BYC_RMB_57","BYC_RMO_45","BYC_RMT_27","BYC_RMT_29","BYC_RMT_46",
             "BYC_RMM_30","BYC_RMT_07","BYC_RMT_28","BYC_RMT_49",
             "BYC_RMT_04","BYC_RMT_06","BYCI_RMT_69","BYCI_RMT_71"),
  pop = c(rep("EAST",5), rep("GAL",4), rep("CTR",3), "WST"))

gl_tar <- vcfR2genlight(vcf_tar)
pop(gl_tar) <- popmap$pop[match(indNames(gl_tar), popmap$sample)]
n_snps_all <- nrow(vcf_tar@fix)
print(nInd(gl_tar))
print(n_snps_all)
print(table(pop(gl_tar)))

pop_labels_df <- data.frame(pop = c("EAST","GAL","CTR","WST"),
                            pop_label = c("Eastern","Galapagos","Central","Western"))

# sampling map global
world_map <- map_data("world", wrap = c(0, 360))
sample_coords <- data.frame(
  sample = popmap$sample,
  lat = c(1.433, 2.217, 2.65, 0.35, 1.767, 2.717, -0.9, 3.35, 2.4, 3.817, 4.283, 4.317, -4.6),
  lon = c(-82.433,-85.317,-84.617,-82.333,-84.8,-91.983,-92.233,-92.1,-90.083,-110.083,-110.083,-144.4,179.883),
  pop = popmap$pop) %>% mutate(lon_pacific = ifelse(lon < 0, lon + 360, lon))

range_map <- ggplot() +
  geom_polygon(data=world_map, aes(x=long,y=lat,group=group), fill="gray90", color="gray70", linewidth=0.2) +
  geom_point(data=sample_coords, aes(x=lon_pacific,y=lat,fill=pop), shape=21, color="black", size=5.5, stroke=0.5) +
  scale_fill_manual(values=pop_colors, name="Sampling\nLocation", breaks=c("EAST","GAL","CTR","WST"),
                    labels=c("Eastern (n=5)","Galapagos (n=4)","Central (n=3)","Western (n=1)")) +
  coord_fixed(ratio = 1, xlim=c(120,300), ylim=c(-35,35), expand=FALSE) +
  scale_x_continuous(breaks=c(140,180,220,260), labels=c("140°E","180°","140°W","100°W")) +
  scale_y_continuous(breaks=c(-30,0,30), labels=c("30°S","0°","30°N")) +
  theme_minimal() +
  theme(legend.position="right", legend.title=element_text(size=13,face="bold"), legend.text=element_text(size=12),
        panel.grid.major=element_line(color="gray80",linewidth=0.3), panel.background=element_rect(fill="aliceblue",color=NA),
        axis.title=element_blank(), axis.text=element_text(size=9)) +
  guides(fill=guide_legend(override.aes=list(size=5)))
ggsave("/Users/madelineeppley/Desktop/manta26pub/tarapacana_map.png", range_map, width=12, height=4, dpi=600)

# relatedness
gds.fn <- "/Users/madelineeppley/Desktop/manta26pub/tar.gds"
snpgdsVCF2GDS("/Users/madelineeppley/Desktop/manta/finalvcfs/tarapacana_final.vcf", gds.fn, method="biallelic.only")
genofile <- snpgdsOpen(gds.fn); rel <- snpgdsIBDKING(genofile, autosome.only=FALSE); snpgdsClose(genofile)
kin_mat <- rel$kinship; sample_names_rel <- rel$sample.id
rownames(kin_mat) <- sample_names_rel; colnames(kin_mat) <- sample_names_rel
annotation_df <- data.frame(Population=popmap$pop[match(sample_names_rel,popmap$sample)]); rownames(annotation_df) <- sample_names_rel
relatedness_plot <- pheatmap(kin_mat, labels_row=sample_names_rel, labels_col=sample_names_rel,
                             clustering_distance_rows="euclidean", clustering_distance_cols="euclidean",
                             main="Pairwise Genomic Relatedness, M. tarapacana", silent=TRUE)
write.csv(kin_mat, "/Users/madelineeppley/Desktop/manta26pub/tarapacana_kinship_matrix.csv")
diag(kin_mat) <- NA # set as NA so it doesnt return 0.5 as the highest relatedness
max_rel <- which(kin_mat == max(kin_mat, na.rm=TRUE), arr.ind=TRUE)[1,]
print(round(max(kin_mat, na.rm=TRUE), 4)) # highest relatedness 0.0108
print(rownames(kin_mat)[max_rel[1]]) # BYC_RMT_49
print(colnames(kin_mat)[max_rel[2]]) # BYC_RMT_06
ggsave("/Users/madelineeppley/Desktop/manta26pub/tarapacana_relatedness.png", relatedness_plot, width=8, height=6, dpi=600)

# pca with all snps
pca_all <- glPca(gl_tar, nf=10); eig_all <- pca_all$eig/sum(pca_all$eig)*100
pca_df_all <- data.frame(PC1=pca_all$scores[,1],PC2=pca_all$scores[,2],sample=indNames(gl_tar),pop=pop(gl_tar)) %>% left_join(pop_labels_df,by="pop")
centroid_all <- pca_df_all %>% group_by(pop) %>% summarise(PC1_cen=mean(PC1),PC2_cen=mean(PC2),n=n())
pca_df_all <- pca_df_all %>% left_join(centroid_all,by="pop")
pops_ell_all <- centroid_all %>% filter(n>=3) %>% pull(pop)
pca_plot_all <- ggplot(pca_df_all, aes(x=PC1,y=PC2)) +
  stat_ellipse(data=filter(pca_df_all,pop%in%pops_ell_all),aes(color=pop),type="norm",level=0.95,linewidth=1.2,show.legend=FALSE) +
  geom_segment(aes(xend=PC1_cen,yend=PC2_cen,color=pop),linewidth=0.5,alpha=0.6,show.legend=FALSE) +
  geom_point(aes(fill=pop,color=pop),size=4,alpha=0.8,shape=21,stroke=0.5) +
  scale_fill_manual(values=pop_colors,breaks=c("EAST","GAL","CTR","WST"),labels=c("Eastern","Galapagos","Central","Western")) +
  scale_color_manual(values=pop_colors,breaks=c("EAST","GAL","CTR","WST"),labels=c("Eastern","Galapagos","Central","Western")) +
  labs(x=paste0("PC1 (",round(eig_all[1],1),"%)"),y=paste0("PC2 (",round(eig_all[2],1),"%)"),subtitle="All SNPs") +
  theme_classic() + theme(legend.position="none",axis.text=element_text(size=13,color="black"),axis.title=element_text(size=14,face="bold"),
                          plot.subtitle=element_text(size=13,hjust=0.5),panel.border=element_rect(color="black",fill=NA,linewidth=1),
                          panel.grid.major=element_line(color="grey90",linewidth=0.3))
ggsave("/Users/madelineeppley/Desktop/manta26pub/tarapacana_PCA_allSNPs.png",
       pca_plot_all+labs(title="Mobula tarapacana",subtitle=paste0("All SNPs (n = ",format(n_snps_all,big.mark=","),")"))+
         theme(legend.position="bottom",legend.direction="horizontal",legend.title=element_blank(),plot.title=element_text(size=13,face="italic",hjust=0.5))+
         guides(fill=guide_legend(override.aes=list(shape=21,size=5,alpha=0.8))), width=8, height=6, dpi=600)

# pcadapt
geno_tar <- read.pcadapt("/Users/madelineeppley/Desktop/manta/finalvcfs/tarapacana_final.vcf", type="vcf")
obj <- pcadapt(geno_tar, K=1)
pvals <- obj$pvalues
qvals <- qvalue(pvals)$qvalues
alpha <- 0.01
outliers <- which(qvals < alpha)
n_outliers <- length(outliers)
n_outliers # 116

# updated 10/20/26 - bonferroni check, same as birostris
padj <- p.adjust(pvals, method = "bonferroni")
alphab <- 0.1
outliers_bonf <- which(padj < alphab)
length(outliers_bonf) # re-run and record
outliers_bonf #97 bonferroni outliers

vcf_out <- vcf_tar; vcf_out@fix <- vcf_tar@fix[outliers,,drop=FALSE]; vcf_out@gt <- vcf_tar@gt[outliers,,drop=FALSE]
write.vcf(vcf_out, "/Users/madelineeppley/Desktop/manta26pub/tarapacana_outliers.vcf.gz")
system(paste0("gunzip -f ","/Users/madelineeppley/Desktop/manta26pub/tarapacana_outliers.vcf.gz"))

# pca with outlier snps
vcf_outliers <- read.vcfR("/Users/madelineeppley/Desktop/manta26pub/tarapacana_outliers.vcf")
gl_outliers <- vcfR2genlight(vcf_outliers); pop(gl_outliers) <- popmap$pop[match(indNames(gl_outliers),popmap$sample)]
pca_out <- glPca(gl_outliers,nf=10); eig_out <- pca_out$eig/sum(pca_out$eig)*100
pca_df_out <- data.frame(PC1=pca_out$scores[,1],PC2=pca_out$scores[,2],sample=indNames(gl_outliers),pop=pop(gl_outliers)) %>% left_join(pop_labels_df,by="pop")
centroid_out <- pca_df_out %>% group_by(pop) %>% summarise(PC1_cen=mean(PC1),PC2_cen=mean(PC2),n=n())
pca_df_out <- pca_df_out %>% left_join(centroid_out,by="pop")
pops_ell_out <- centroid_out %>% filter(n>=3) %>% pull(pop)
pca_plot_outliers <- ggplot(pca_df_out,aes(x=PC1,y=PC2)) +
  stat_ellipse(data=filter(pca_df_out,pop%in%pops_ell_out),aes(color=pop),type="norm",level=0.95,linewidth=1.2,show.legend=FALSE) +
  geom_segment(aes(xend=PC1_cen,yend=PC2_cen,color=pop),linewidth=0.5,alpha=0.6,show.legend=FALSE) +
  geom_point(aes(fill=pop,color=pop),size=4,alpha=0.8,shape=21,stroke=0.5) +
  scale_fill_manual(values=pop_colors,breaks=c("EAST","GAL","CTR","WST"),labels=c("Eastern","Galapagos","Central","Western")) +
  scale_color_manual(values=pop_colors,breaks=c("EAST","GAL","CTR","WST"),labels=c("Eastern","Galapagos","Central","Western")) +
  labs(x=paste0("PC1 (",round(eig_out[1],1),"%)"),y=paste0("PC2 (",round(eig_out[2],1),"%)"),subtitle="Outlier SNPs") +
  theme_classic() + theme(legend.position="none",axis.text=element_text(size=13,color="black"),axis.title=element_text(size=14,face="bold"),
                          plot.subtitle=element_text(size=13,hjust=0.5),panel.border=element_rect(color="black",fill=NA,linewidth=1),
                          panel.grid.major=element_line(color="grey90",linewidth=0.3))
ggsave("/Users/madelineeppley/Desktop/manta26pub/tarapacana_PCA_outlierSNPs.png",
       pca_plot_outliers+labs(title="Mobula tarapacana",subtitle=paste0("Outlier SNPs (n = ",format(n_outliers,big.mark=","),")"))+
         theme(legend.position="bottom",legend.direction="horizontal",legend.title=element_blank(),plot.title=element_text(size=13,face="italic",hjust=0.5))+
         guides(fill=guide_legend(override.aes=list(shape=21,size=5,alpha=0.8))), width=8, height=6, dpi=600)

# dapc
dapc_result <- dapc(gl_tar,pop=pop(gl_tar),n.pca=4,n.da=3)
dapc_df <- data.frame(LD1=dapc_result$ind.coord[,1],LD2=dapc_result$ind.coord[,2],sample=indNames(gl_tar),pop=pop(gl_tar))
centroid_dapc <- dapc_df %>% group_by(pop) %>% summarise(LD1_cen=mean(LD1),LD2_cen=mean(LD2),n=n())
dapc_df <- dapc_df %>% left_join(centroid_dapc,by="pop")
pops_ell_dapc <- centroid_dapc %>% filter(n>=3) %>% pull(pop)
pct_dapc <- dapc_result$eig/sum(dapc_result$eig)*100
dapc_plot <- ggplot(dapc_df,aes(x=LD1,y=LD2)) +
  stat_ellipse(data=filter(dapc_df,pop%in%pops_ell_dapc),aes(color=pop),type="norm",level=0.95,linewidth=1.2,show.legend=FALSE) +
  geom_segment(aes(xend=LD1_cen,yend=LD2_cen,color=pop),linewidth=0.5,alpha=0.6,show.legend=FALSE) +
  geom_point(aes(fill=pop,color=pop),size=4,alpha=0.8,shape=21,stroke=0.5) +
  scale_fill_manual(values=pop_colors,breaks=c("EAST","GAL","CTR","WST"),labels=c("Eastern","Galapagos","Central","Western")) +
  scale_color_manual(values=pop_colors,breaks=c("EAST","GAL","CTR","WST"),labels=c("Eastern","Galapagos","Central","Western")) +
  labs(x=paste0("DAPC1 (",round(pct_dapc[1],1),"%)"),y=paste0("DAPC2 (",round(pct_dapc[2],1),"%)"),subtitle="All SNPs") +
  theme_classic() + theme(legend.position="none",axis.text=element_text(size=13,color="black"),axis.title=element_text(size=14,face="bold"),
                          panel.border=element_rect(color="black",fill=NA,linewidth=1),panel.grid.major=element_line(color="grey90",linewidth=0.3))
ggsave("/Users/madelineeppley/Desktop/manta26pub/tarapacana_DAPC.png",
       dapc_plot+labs(title="Mobula tarapacana")+theme(legend.position="bottom",legend.direction="horizontal",legend.title=element_blank(),
                                                       plot.title=element_text(size=13,face="italic",hjust=0.5))+guides(fill=guide_legend(override.aes=list(shape=21,size=5,alpha=0.8))),
       width=8, height=6, dpi=600)

# snmf analysis
#unlink("/Users/madelineeppley/Desktop/manta26pub/tarapacana.geno")
#unlink("/Users/madelineeppley/Desktop/manta26pub/tarapacana.snmfProject", recursive=TRUE)
#unlink("/Users/madelineeppley/Desktop/manta26pub/tarapacana.snmf", recursive=TRUE)
#vcf2geno("/Users/madelineeppley/Desktop/manta/finalvcfs/tarapacana_final.vcf", "/Users/madelineeppley/Desktop/manta26pub/tarapacana.geno")
#project_all <- snmf("/Users/madelineeppley/Desktop/manta26pub/tarapacana.geno", K=1:5, entropy=TRUE, repetitions=10, project="new")

#unlink("/Users/madelineeppley/Desktop/manta26pub/tarapacana_outliers.geno")
#unlink("/Users/madelineeppley/Desktop/manta26pub/tarapacana_outliers.snmfProject", recursive=TRUE)
#unlink("/Users/madelineeppley/Desktop/manta26pub/tarapacana_outliers.snmf", recursive=TRUE)
#vcf2geno("/Users/madelineeppley/Desktop/manta26pub/tarapacana_outliers.vcf",
#        "/Users/madelineeppley/Desktop/manta26pub/tarapacana_outliers.geno")
#project_outliers <- snmf("/Users/madelineeppley/Desktop/manta26pub/tarapacana_outliers.geno", K=1:10, entropy=TRUE, repetitions=10, project="new")

# ancestry proportions with SNMF on all snps
## JUST RELOAD
# needed to double-check k values, re-load
project_all <- load.snmfProject("/Users/madelineeppley/Desktop/manta26pub/tarapacana.snmfProject")
project_outliers <- load.snmfProject("/Users/madelineeppley/Desktop/manta26pub/tarapacana_outliers.snmfProject")
for(k in 1:5) cat("K =", k, ":", min(cross.entropy(project_all, K=k)), "\n") # lowest k = 1, 0.749
for(k in 1:5) cat("K =", k, ":", min(cross.entropy(project_outliers, K=k)), "\n") # lowest k = 3, 0.040

# ce plot for all snps
png("/Users/madelineeppley/Desktop/manta26pub/tarapacana_CE_allSNPs.png", width=800, height=600)
plot(project_all, col="blue", pch=19, cex=1.2, main="Cross-entropy: M. tarapacana - All SNPs")
dev.off()
print("cross entropy all snps")
for(k in 1:5) cat("K =", k, ":", min(cross.entropy(project_all, K=k)), "\n")

# ce plot with outlier snps
png("/Users/madelineeppley/Desktop/manta26pub/tarapacana_CE_outlierSNPs.png", width=800, height=600)
plot(project_outliers, col="red", pch=19, cex=1.2, main="Cross-entropy: M. tarapacana - Outlier SNPs")
dev.off()
print("cross entropy outlier snps")
for(k in 1:10) cat("K =", k, ":", min(cross.entropy(project_outliers, K=k)), "\n")

sample_names <- indNames(gl_tar); sample_pops <- pop(gl_tar)
pop_order <- c("WST","CTR", "GAL", "EAST")
pop_names_map <- c("WST"="West","CTR"="Central","GAL"="Galapagos","EAST"="Eastern")
pop_label_df_str <- data.frame(Sample=sample_names,Pop=factor(sample_pops,levels=pop_order)) %>% arrange(Pop) %>% mutate(x_pos=1:n())
pop_boundaries <- pop_label_df_str %>% group_by(Pop) %>% summarise(x_end=max(x_pos)) %>% filter(x_end!=max(pop_label_df_str$x_pos))
pop_spans <- pop_label_df_str %>% group_by(Pop) %>% summarise(x_start=min(x_pos),x_end=max(x_pos),x_mid=(min(x_pos)+max(x_pos))/2,n=n()) %>%
  mutate(display_name=pop_names_map[as.character(Pop)])
max_x <- max(pop_label_df_str$x_pos)

ce_all <- cross.entropy(project_all,K=2); qmat_all <- Q(project_all,K=2,run=which.min(ce_all))
qdf_all <- as.data.frame(qmat_all); colnames(qdf_all) <- paste0("Cluster",1:2)
qdf_all$Sample <- sample_names; qdf_all$Pop <- factor(sample_pops,levels=pop_order)
qdf_all <- qdf_all %>% arrange(Pop) %>% mutate(x_pos=1:n())
qdf_all_long <- qdf_all %>% pivot_longer(cols=starts_with("Cluster"),names_to="Cluster",values_to="Proportion")
structure_all_k2 <- ggplot(qdf_all_long,aes(x=x_pos,y=Proportion,fill=Cluster)) +
  geom_bar(stat="identity",width=1,color=NA) + geom_vline(data=pop_boundaries,aes(xintercept=x_end+0.5),color="white",linewidth=0.8) +
  scale_fill_manual(values=structure_colors[1:2]) + scale_y_continuous(expand=c(0,0),breaks=c(0,0.5,1),labels=c("0.0","0.5","1.0"),position="left") +
  scale_x_continuous(expand=c(0,0),limits=c(0.5,max_x+3)) + annotate("text",x=max_x+2,y=0.5,label="K = 2",size=6,fontface="bold") +
  ggtitle("All SNPs") + coord_cartesian(clip="off") + theme_minimal() +
  theme(plot.title=element_text(size=12,face="bold",hjust=0),axis.text.x=element_blank(),axis.ticks.x=element_blank(),
        axis.title=element_blank(),axis.text.y=element_text(size=12),legend.position="none",panel.grid=element_blank(),plot.margin=margin(5,30,2,5))

ce_out <- cross.entropy(project_outliers,K=3); qmat_out_str <- Q(project_outliers,K=3,run=which.min(ce_out))
qdf_out_str <- as.data.frame(qmat_out_str); colnames(qdf_out_str) <- paste0("Cluster",1:3)
qdf_out_str$Sample <- sample_names; qdf_out_str$Pop <- factor(sample_pops,levels=pop_order)
qdf_out_str <- qdf_out_str %>% arrange(Pop) %>% mutate(x_pos=1:n())
qdf_out_long <- qdf_out_str %>% pivot_longer(cols=starts_with("Cluster"),names_to="Cluster",values_to="Proportion")
structure_outlier_k3 <- ggplot(qdf_out_long,aes(x=x_pos,y=Proportion,fill=Cluster)) +
  geom_bar(stat="identity",width=1,color=NA) + geom_vline(data=pop_boundaries,aes(xintercept=x_end+0.5),color="white",linewidth=0.8) +
  scale_fill_manual(values=structure_colors[c(1,2,3)]) + scale_y_continuous(expand=c(0,0),breaks=c(0,0.5,1),labels=c("0.0","0.5","1.0"),position="left") +
  scale_x_continuous(expand=c(0,0),limits=c(0.5,max_x+3)) + annotate("text",x=max_x+2,y=0.5,label="K = 3",size=6,fontface="bold") +
  ggtitle("Outlier SNPs") + coord_cartesian(clip="off") + theme_minimal() +
  theme(plot.title=element_text(size=12,face="bold",hjust=0),axis.text.x=element_blank(),axis.ticks.x=element_blank(),
        axis.title=element_blank(),axis.text.y=element_text(size=12),legend.position="none",panel.grid=element_blank(),plot.margin=margin(5,30,2,5))

pop_label_plot <- ggplot() + geom_text(data=pop_spans,aes(x=x_mid,y=0.5,label=display_name),size=5,fontface="bold") +
  geom_vline(data=pop_boundaries,aes(xintercept=x_end+0.5),color="black",linewidth=1) +
  scale_x_continuous(limits=c(0.5,max_x+3),expand=c(0,0)) + scale_y_continuous(limits=c(0,1),expand=c(0,0)) +
  coord_cartesian(clip="off") + theme_void() + theme(plot.margin=margin(0,30,10,5))

# amova FST
# updated 10/20/26: pairwise_amova() removed, pairwise p-values now come from gl.fst.pop
genI <- gl2gi(gl_tar); pops_gi <- genI$pop; strata(genI) <- data.frame(pops=pops_gi)
p.amova <- poppr.amova(genI,~pops); amova.pvalues <- ade4::randtest(p.amova,nrepet=9999)
overall_fst <- p.amova$statphi$Phi[length(p.amova$statphi$Phi)]
# updated 10/20/26
print(amova.pvalues) # take row 3 from here for variation between pops
overall_pval <- amova.pvalues$pvalue[length(amova.pvalues$pvalue)]
overall_fst # -0.0002 
overall_pval # was 0.5395 with pvalue[1], is now 0.5932

# pairwise fst excluding WST (n=1)
pop_sizes <- table(pop(gl_tar))
valid_pops <- names(pop_sizes[pop_sizes >= 2])
gl_valid <- gl_tar[pop(gl_tar) %in% valid_pops, ]
fst_result <- gl.fst.pop(gl_valid, nboots=1000, percent=95, nclusters=1)
fst_matrix <- fst_result$Fsts
fst_matrix[is.na(fst_matrix)] <- t(fst_matrix)[is.na(fst_matrix)]
fst_matrix[fst_matrix < 0] <- 0

# updated 10/20/26
pval_matrix <- fst_result$Pvalues
pval_matrix[is.na(pval_matrix)] <- t(pval_matrix)[is.na(pval_matrix)]
print("pairwise fst")
print(round(fst_matrix, 4))
print(pval_matrix)

## outlier fst
table(ploidy(gl_outliers))
ploidy(gl_outliers) <- 2 # set ploidy to 2 for all
nLoc(gl_outliers) # 116 outliers
pop(gl_outliers) <- popmap$pop[match(indNames(gl_outliers), popmap$sample)]
genI_out <- gl2gi(gl_outliers)
pops_gi_out <- genI_out$pop
strata(genI_out) <- data.frame(pops = pops_gi_out)

p.amova_out       <- poppr.amova(genI_out, ~pops)
amova.pvalues_out <- ade4::randtest(p.amova_out, nrepet = 9999)
overall_fst_out   <- p.amova_out$statphi$Phi[length(p.amova_out$statphi$Phi)]
# updated 10/20/26
print(amova.pvalues_out)
overall_pval_out <- amova.pvalues_out$pvalue[length(amova.pvalues_out$pvalue)]

# outlier amova fst and pval
round(overall_fst_out, 4) # fst = 0.1932
round(overall_pval_out, 4) # updated 0.1071

# pairwise outlier FST excluding WST (n=1)
pop_sizes_out  <- table(pop(gl_outliers))
valid_pops_out <- names(pop_sizes_out[pop_sizes_out >= 2])
gl_valid_out   <- gl_outliers[pop(gl_outliers) %in% valid_pops_out, ]

fst_result_out <- gl.fst.pop(gl_valid_out, nboots=1000, percent=95, nclusters=1)
fst_matrix_out <- fst_result_out$Fsts
fst_matrix_out[is.na(fst_matrix_out)] <- t(fst_matrix_out)[is.na(fst_matrix_out)]

# updated 10/20/26
pval_matrix_out <- fst_result_out$Pvalues
pval_matrix_out[is.na(pval_matrix_out)] <- t(pval_matrix_out)[is.na(pval_matrix_out)]
print("pairwise fst outlier snps")
print(round(fst_matrix_out, 4))
print(pval_matrix_out)

# neutral loci fst
# updated 10/20/26
stopifnot(nLoc(gl_tar) == nrow(vcf_tar@fix))
neutral <- setdiff(seq_len(nLoc(gl_tar)), outliers)
gl_neutral <- gl_tar[, neutral]
nLoc(gl_neutral) # 43622 = n_snps_all - 116
genI_neu <- gl2gi(gl_neutral)
strata(genI_neu) <- data.frame(pops = genI_neu$pop)
p.amova_neu <- poppr.amova(genI_neu, ~pops)
amova.pvalues_neu <- ade4::randtest(p.amova_neu, nrepet = 9999)
print(amova.pvalues_neu)
overall_fst_neu  <- p.amova_neu$statphi$Phi[length(p.amova_neu$statphi$Phi)]
overall_pval_neu <- amova.pvalues_neu$pvalue[length(amova.pvalues_neu$pvalue)]
round(overall_fst_neu, 4) # -9e-04
round(overall_pval_neu, 4) # 0.8337

stars <- function(p) ifelse(is.na(p), "",
                            ifelse(p < 0.001, "***",
                                   ifelse(p < 0.01,  "**",
                                          ifelse(p < 0.05,  "*", ""))))
cat("  All-loci FST =", round(overall_fst, 4),  stars(overall_pval),     "\n")
cat("  Outlier FST  =", round(overall_fst_out,4), stars(overall_pval_out), "\n")

# updated 10/20/26
p_label <- ifelse(overall_pval < 0.001, "< 0.001", paste("=", round(overall_pval, 3)))

pop_name_map <- c("EAST"="Eastern","GAL"="Galapagos","CTR"="Central")
fst_matrix_labeled <- fst_matrix; rownames(fst_matrix_labeled) <- pop_name_map[rownames(fst_matrix_labeled)]
colnames(fst_matrix_labeled) <- pop_name_map[colnames(fst_matrix_labeled)]
pval_matrix_labeled <- pval_matrix; rownames(pval_matrix_labeled) <- pop_name_map[rownames(pval_matrix)]
colnames(pval_matrix_labeled) <- pop_name_map[colnames(pval_matrix)]
get_upper_tri <- function(mat){mat[lower.tri(mat)]<-NA;return(mat)}
fst_melted <- melt(get_upper_tri(fst_matrix_labeled),na.rm=TRUE)
pval_melted <- melt(get_upper_tri(pval_matrix_labeled),na.rm=TRUE)
fst_melted <- fst_melted %>% left_join(pval_melted,by=c("Var1","Var2"),suffix=c("","_pval")) %>%
  mutate(sig=case_when(is.na(value_pval)~"",value_pval<0.001~"***",value_pval<0.01~"**",value_pval<0.05~"*",TRUE~""),
         label=paste0(round(value,4),sig))

fst_heatmap <- ggplot(fst_melted,aes(Var1,Var2,fill=value)) +
  geom_tile(color="white",linewidth=1) + geom_text(aes(label=label),size=5) +
  scale_fill_gradient(low="white",high="#1D8981",name=expression(italic("F")[ST]),limits=c(0,max(fst_melted$value,na.rm=TRUE)+0.001)) +
  coord_fixed() +
  labs(x="",y="",subtitle=bquote(atop("AMOVA"~italic("F")[ST]~"="~.(round(overall_fst,4))~"; p"~.(p_label),"* p < 0.05, ** p < 0.01, *** p < 0.001"))) +
  theme_minimal() + theme(axis.text.x=element_text(angle=45,hjust=1,size=12,face="bold"),axis.text.y=element_text(size=12,face="bold"),
                          plot.subtitle=element_text(size=11,hjust=0.5,face="bold"),panel.grid=element_blank(),legend.position="right")
ggsave("/Users/madelineeppley/Desktop/manta26pub/tarapacana_Fst_heatmap.png", fst_heatmap, width=8, height=6, dpi=600)

# create the final figure!
layer1 <- range_map+theme(legend.position="right",plot.margin=margin(5,5,5,5))
layer2 <- wrap_plots(pca_plot_all,pca_plot_outliers,ncol=2)
layer3 <- wrap_plots(structure_all_k2,structure_outlier_k3,pop_label_plot,ncol=1,heights=c(1,1,0.3))
layer4 <- wrap_plots(fst_heatmap+theme(legend.position="right"),dapc_plot,ncol=2,widths=c(0.8,1.2))
final_figure <- wrap_plots(layer1,layer2,layer3,layer4,ncol=1,heights=c(1.3,1.1,0.8,1.0)) +
  plot_annotation(title="Mobula tarapacana",
                  subtitle=paste0("All SNPs: n = ",format(n_snps_all,big.mark=","),"   |   Outlier SNPs: n = ",format(n_outliers,big.mark=",")),
                  tag_levels=list(c('A','B','C','D','E','','F','G')),
                  theme=theme(plot.title=element_text(size=18,face="bold.italic",hjust=0),plot.subtitle=element_text(size=14,hjust=0),plot.tag=element_text(size=16,face="bold")))
ggsave("/Users/madelineeppley/Desktop/manta26pub/tarapacana_FINALFIGURE.png", final_figure, width=14, height=18, dpi=600)
