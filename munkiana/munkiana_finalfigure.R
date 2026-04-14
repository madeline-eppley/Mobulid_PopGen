### M. munkiana analysis to produce publication figure ####
# M. Eppley, v2. current version uploaded to github 4/14/2026

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

outdir <- "/Users/madelineeppley/Desktop/manta26pub"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
setwd(outdir)  # so pcadapt can write positions.txt here

pop_colors <- c("MEX" = "#6DBAA4", "NIC" = "#ED9A6C", "PER" = "#8C9FCB")
structure_colors <- c("#264653", "#2a9d8f", "#8ab17d", "#e9c46a", "#f4a261", "#e76f51")

vcf_mun <- read.vcfR("/Users/madelineeppley/Desktop/manta/finalvcfs/munkiana_filtered_7188snps.vcf")
popmap <- data.frame(
  sample = c("mu_11","mu_12","mu_21","mu_3","mu_4","mu_7",
             "NIC_2728_01","NIC_2740_01",
             "PER_DCLCW130_MM","PER_DCLCW87_MM","PER_DPPW019_MM","PER_DPPW28_MM","PER_DPPW98_MM"),
  pop = c(rep("MEX",6), rep("NIC",2), rep("PER",5)))

gl_mun <- vcfR2genlight(vcf_mun)
pop(gl_mun) <- popmap$pop[match(indNames(gl_mun), popmap$sample)]
n_snps_all <- nrow(vcf_mun@fix)
cat("Munkiana:", nInd(gl_mun), "individuals,", n_snps_all, "SNPs\n")
print(table(pop(gl_mun)))

pop_labels_df <- data.frame(pop=c("MEX","NIC","PER"), pop_label=c("Mexico","Nicaragua","Peru"))

# map with sampling coords
sample_coords <- data.frame(
  sample = popmap$sample,
  lat = c(24.1,24.3,23.6,24.5,23.8,24.0, 11.2,11.5, -6.0,-7.5,-8.0,-9.5,-10.0),
  lon = c(-109.8,-110.0,-110.3,-109.5,-110.2,-109.7, -87.0,-87.3, -81.0,-80.5,-81.2,-80.8,-81.5),
  pop = popmap$pop)
world_map <- map_data("world")

# IUCN range since this is not a globally distributed species
tryCatch({
  iucn_data <- st_read("/Users/madelineeppley/Desktop/manta/rangemap/SHARKS_RAYS_CHIMAERAS/SHARKS_RAYS_CHIMAERAS.shp", quiet=TRUE)
  munkiana_range <- iucn_data %>% filter(sci_name == "Mobula munkiana") %>% st_make_valid()
  munkiana_range_df <- fortify(as(munkiana_range, "Spatial"))
  has_range <- TRUE
}, error = function(e) { has_range <<- FALSE })

range_map <- ggplot() +
  geom_polygon(data=world_map, aes(x=long,y=lat,group=group), fill="gray90", color="gray70", linewidth=0.2)
if(exists("has_range") && has_range) range_map <- range_map + geom_polygon(data=munkiana_range_df, aes(x=long,y=lat,group=group), fill="#6DBAA4", alpha=0.3, color=NA)
range_map <- range_map +
  geom_point(data=sample_coords, aes(x=lon,y=lat,fill=pop), shape=21, color="black", size=5.5, stroke=0.5) +
  scale_fill_manual(values=pop_colors, name="Sampling\nLocation", breaks=c("MEX","NIC","PER"),
                    labels=c("Mexico (n=6)","Nicaragua (n=2)","Peru (n=5)")) +
  coord_cartesian(xlim=c(-170,-20), ylim=c(-25,35), expand=FALSE) +
  theme_minimal() + theme(legend.position="right",legend.title=element_text(size=13,face="bold"),legend.text=element_text(size=12),
                          panel.grid.major=element_line(color="gray80",linewidth=0.3),panel.background=element_rect(fill="aliceblue",color=NA),
                          axis.title=element_blank(),axis.text=element_text(size=9)) + guides(fill=guide_legend(override.aes=list(size=5)))
ggsave("/Users/madelineeppley/Desktop/manta26pub/munkiana_map.png", range_map, width=12, height=4, dpi=600)

# relatedness
gds.fn <- "/Users/madelineeppley/Desktop/manta26pub/mun.gds"; showfile.gds(closeall=TRUE)
snpgdsVCF2GDS("/Users/madelineeppley/Desktop/manta/finalvcfs/munkiana_filtered_7188snps.vcf",gds.fn,method="biallelic.only")
genofile <- snpgdsOpen(gds.fn); rel <- snpgdsIBDKING(genofile,autosome.only=FALSE); snpgdsClose(genofile)
kin_mat <- rel$kinship; sample_names_rel <- rel$sample.id
rownames(kin_mat) <- sample_names_rel; colnames(kin_mat) <- sample_names_rel
annotation_df <- data.frame(Population=popmap$pop[match(sample_names_rel,popmap$sample)]); rownames(annotation_df) <- sample_names_rel
relate <- pheatmap(kin_mat,labels_row=sample_names_rel,labels_col=sample_names_rel,annotation_row=annotation_df,annotation_col=annotation_df,
                   annotation_colors=list(Population=pop_colors),clustering_distance_rows="euclidean",clustering_distance_cols="euclidean",
                   main="Pairwise Genomic Relatedness, M. munkiana",silent=TRUE)
ggsave("/Users/madelineeppley/Desktop/manta26pub/munkiana_relatedness.png", relate, width=10, height=8, dpi=600)

# pca using all snps
pca_all <- glPca(gl_mun,nf=10); eig_all <- pca_all$eig/sum(pca_all$eig)*100
pca_df_all <- data.frame(PC1=pca_all$scores[,1],PC2=pca_all$scores[,2],sample=indNames(gl_mun),pop=pop(gl_mun)) %>% left_join(pop_labels_df,by="pop")
centroid_all <- pca_df_all %>% group_by(pop) %>% summarise(PC1_cen=mean(PC1),PC2_cen=mean(PC2),n=n())
pca_df_all <- pca_df_all %>% left_join(centroid_all,by="pop")
pops_ell_all <- centroid_all %>% filter(n>=3) %>% pull(pop)
pca_plot_all <- ggplot(pca_df_all,aes(x=PC1,y=PC2)) +
  stat_ellipse(data=filter(pca_df_all,pop%in%pops_ell_all),aes(color=pop),type="norm",level=0.95,linewidth=1.2,show.legend=FALSE) +
  geom_segment(aes(xend=PC1_cen,yend=PC2_cen,color=pop),linewidth=0.5,alpha=0.6,show.legend=FALSE) +
  geom_point(aes(fill=pop,color=pop),size=4,alpha=0.8,shape=21,stroke=0.5) +
  scale_fill_manual(values=pop_colors,breaks=c("MEX","NIC","PER"),labels=c("Mexico","Nicaragua","Peru")) +
  scale_color_manual(values=pop_colors,breaks=c("MEX","NIC","PER"),labels=c("Mexico","Nicaragua","Peru")) +
  labs(x=paste0("PC1 (",round(eig_all[1],1),"%)"),y=paste0("PC2 (",round(eig_all[2],1),"%)"),subtitle="All SNPs") +
  theme_classic() + theme(legend.position="none",axis.text=element_text(size=13,color="black"),axis.title=element_text(size=14,face="bold"),
                          plot.subtitle=element_text(size=13,hjust=0.5),panel.border=element_rect(color="black",fill=NA,linewidth=1),
                          panel.grid.major=element_line(color="grey90",linewidth=0.3))
ggsave("/Users/madelineeppley/Desktop/manta26pub/munkiana_PCA_allSNPs.png",
       pca_plot_all+labs(title="Mobula munkiana")+theme(legend.position="bottom",legend.direction="horizontal",legend.title=element_blank(),
                                                        plot.title=element_text(size=13,face="italic",hjust=0.5))+guides(fill=guide_legend(override.aes=list(shape=21,size=5,alpha=0.8))),width=8,height=6,dpi=600)

# pcadapt
geno_mun <- read.pcadapt("/Users/madelineeppley/Desktop/manta/finalvcfs/munkiana_filtered_7188snps.vcf",type="vcf")
obj <- pcadapt(geno_mun,K=1); pvals <- obj$pvalues; outliers <- which(pvals<0.01); n_outliers <- length(outliers)
n_outliers
vcf_out <- vcf_mun; vcf_out@fix <- vcf_mun@fix[outliers,,drop=FALSE]; vcf_out@gt <- vcf_mun@gt[outliers,,drop=FALSE]
write.vcf(vcf_out,"/Users/madelineeppley/Desktop/manta26pub/munkiana_outliers.vcf.gz"); system(paste0("gunzip -f ","/Users/madelineeppley/Desktop/manta26pub/munkiana_outliers.vcf.gz"))

# pca with outliers
vcf_outliers <- read.vcfR("/Users/madelineeppley/Desktop/manta26pub/munkiana_outliers.vcf")
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
  scale_fill_manual(values=pop_colors,breaks=c("MEX","NIC","PER"),labels=c("Mexico","Nicaragua","Peru")) +
  scale_color_manual(values=pop_colors,breaks=c("MEX","NIC","PER"),labels=c("Mexico","Nicaragua","Peru")) +
  labs(x=paste0("PC1 (",round(eig_out[1],1),"%)"),y=paste0("PC2 (",round(eig_out[2],1),"%)"),subtitle="Outlier SNPs") +
  theme_classic() + theme(legend.position="none",axis.text=element_text(size=13,color="black"),axis.title=element_text(size=14,face="bold"),
                          plot.subtitle=element_text(size=13,hjust=0.5),panel.border=element_rect(color="black",fill=NA,linewidth=1),
                          panel.grid.major=element_line(color="grey90",linewidth=0.3))
ggsave("/Users/madelineeppley/Desktop/manta26pub/munkiana_PCA_outlierSNPs.png",
       pca_plot_outliers+labs(title="Mobula munkiana")+theme(legend.position="bottom",legend.direction="horizontal",legend.title=element_blank(),
                                                             plot.title=element_text(size=13,face="italic",hjust=0.5))+guides(fill=guide_legend(override.aes=list(shape=21,size=5,alpha=0.8))),width=8,height=6,dpi=600)

# dapc
dapc_result <- dapc(gl_mun,pop=pop(gl_mun),n.pca=4,n.da=2)
dapc_df <- data.frame(LD1=dapc_result$ind.coord[,1],LD2=dapc_result$ind.coord[,2],sample=indNames(gl_mun),pop=pop(gl_mun))
centroid_dapc <- dapc_df %>% group_by(pop) %>% summarise(LD1_cen=mean(LD1),LD2_cen=mean(LD2),n=n())
dapc_df <- dapc_df %>% left_join(centroid_dapc,by="pop")
pops_ell_dapc <- centroid_dapc %>% filter(n>=3) %>% pull(pop)
pct_dapc <- dapc_result$eig/sum(dapc_result$eig)*100
dapc_plot <- ggplot(dapc_df,aes(x=LD1,y=LD2)) +
  stat_ellipse(data=filter(dapc_df,pop%in%pops_ell_dapc),aes(color=pop),type="norm",level=0.95,linewidth=1.2,show.legend=FALSE) +
  geom_segment(aes(xend=LD1_cen,yend=LD2_cen,color=pop),linewidth=0.5,alpha=0.6,show.legend=FALSE) +
  geom_point(aes(fill=pop,color=pop),size=4,alpha=0.8,shape=21,stroke=0.5) +
  scale_fill_manual(values=pop_colors,breaks=c("MEX","NIC","PER"),labels=c("Mexico","Nicaragua","Peru")) +
  scale_color_manual(values=pop_colors,breaks=c("MEX","NIC","PER"),labels=c("Mexico","Nicaragua","Peru")) +
  labs(x=paste0("DAPC1 (",round(pct_dapc[1],1),"%)"),y=paste0("DAPC2 (",round(pct_dapc[2],1),"%)"),subtitle="All SNPs") +
  theme_classic() + theme(legend.position="none",axis.text=element_text(size=13,color="black"),axis.title=element_text(size=14,face="bold"),
                          panel.border=element_rect(color="black",fill=NA,linewidth=1),panel.grid.major=element_line(color="grey90",linewidth=0.3))
ggsave("/Users/madelineeppley/Desktop/manta26pub/munkiana_DAPC.png",
       dapc_plot+labs(title="Mobula munkiana")+theme(legend.position="bottom",legend.direction="horizontal",legend.title=element_blank(),
                                                     plot.title=element_text(size=13,face="italic",hjust=0.5))+guides(fill=guide_legend(override.aes=list(shape=21,size=5,alpha=0.8))),width=8,height=6,dpi=600)

# snmf ancestry
unlink("/Users/madelineeppley/Desktop/manta26pub/munkiana.geno")
unlink("/Users/madelineeppley/Desktop/manta26pub/munkiana.snmfProject", recursive=TRUE)
unlink("/Users/madelineeppley/Desktop/manta26pub/munkiana.snmf", recursive=TRUE)
vcf2geno("/Users/madelineeppley/Desktop/manta/finalvcfs/munkiana_filtered_7188snps.vcf","/Users/madelineeppley/Desktop/manta26pub/munkiana.geno")
project_all <- snmf("/Users/madelineeppley/Desktop/manta26pub/munkiana.geno", K=1:5, entropy=TRUE, repetitions=10, project="new")

# ce plot for all snps
png("/Users/madelineeppley/Desktop/manta26pub/munkiana_CE_allSNPs.png", width=800, height=600)
plot(project_all, col="blue", pch=19, cex=1.2, main="Cross-entropy: M. munkiana - All SNPs")
dev.off()
for(k in 1:5) cat("K =", k, ":", min(cross.entropy(project_all, K=k)), "\n")
unlink("/Users/madelineeppley/Desktop/manta26pub/munkiana_outliers.geno")
unlink("/Users/madelineeppley/Desktop/manta26pub/munkiana_outliers.snmfProject", recursive=TRUE)
unlink("/Users/madelineeppley/Desktop/manta26pub/munkiana_outliers.snmf", recursive=TRUE)
vcf2geno("/Users/madelineeppley/Desktop/manta26pub/munkiana_outliers.vcf",
         "/Users/madelineeppley/Desktop/manta26pub/munkiana_outliers.geno")
project_outliers <- snmf("/Users/madelineeppley/Desktop/manta26pub/munkiana_outliers.geno",K=1:5,entropy=TRUE,repetitions=10,project="new")

# ce plot for outlier snps
png("/Users/madelineeppley/Desktop/manta26pub/munkiana_CE_outlierSNPs.png", width=800, height=600)
plot(project_outliers, col="red", pch=19, cex=1.2, main="Cross-entropy: M. munkiana - Outlier SNPs")
dev.off()
for(k in 1:5) cat("K =", k, ":", min(cross.entropy(project_outliers, K=k)), "\n")

sample_names <- indNames(gl_mun); sample_pops <- pop(gl_mun)
pop_order <- c("MEX","NIC","PER"); pop_names_map <- c("MEX"="Mexico","NIC"="Nicaragua","PER"="Peru")
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
structure_all <- ggplot(qdf_all_long,aes(x=x_pos,y=Proportion,fill=Cluster)) +
  geom_bar(stat="identity",width=1,color=NA) + geom_vline(data=pop_boundaries,aes(xintercept=x_end+0.5),color="white",linewidth=0.8) +
  scale_fill_manual(values=structure_colors[1:2]) + scale_y_continuous(expand=c(0,0),breaks=c(0,0.5,1),labels=c("0.0","0.5","1.0"),position="left") +
  scale_x_continuous(expand=c(0,0),limits=c(0.5,max_x+3)) + annotate("text",x=max_x+2,y=0.5,label="K = 2",size=6,fontface="bold") +
  ggtitle("All SNPs") + coord_cartesian(clip="off") + theme_minimal() +
  theme(plot.title=element_text(size=12,face="bold",hjust=0),axis.text.x=element_blank(),axis.ticks.x=element_blank(),
        axis.title=element_blank(),axis.text.y=element_text(size=12),legend.position="none",panel.grid=element_blank(),plot.margin=margin(5,30,2,5))

ce_out <- cross.entropy(project_outliers,K=2); qmat_out_str <- Q(project_outliers,K=2,run=which.min(ce_out))
qdf_out_str <- as.data.frame(qmat_out_str); colnames(qdf_out_str) <- paste0("Cluster",1:2)
qdf_out_str$Sample <- sample_names; qdf_out_str$Pop <- factor(sample_pops,levels=pop_order)
qdf_out_str <- qdf_out_str %>% arrange(Pop) %>% mutate(x_pos=1:n())
qdf_out_long <- qdf_out_str %>% pivot_longer(cols=starts_with("Cluster"),names_to="Cluster",values_to="Proportion")
structure_outlier <- ggplot(qdf_out_long,aes(x=x_pos,y=Proportion,fill=Cluster)) +
  geom_bar(stat="identity",width=1,color=NA) + geom_vline(data=pop_boundaries,aes(xintercept=x_end+0.5),color="white",linewidth=0.8) +
  scale_fill_manual(values=structure_colors[c(3,4)]) + scale_y_continuous(expand=c(0,0),breaks=c(0,0.5,1),labels=c("0.0","0.5","1.0"),position="left") +
  scale_x_continuous(expand=c(0,0),limits=c(0.5,max_x+3)) + annotate("text",x=max_x+2,y=0.5,label="K = 2",size=6,fontface="bold") +
  ggtitle("Outlier SNPs") + coord_cartesian(clip="off") + theme_minimal() +
  theme(plot.title=element_text(size=12,face="bold",hjust=0),axis.text.x=element_blank(),axis.ticks.x=element_blank(),
        axis.title=element_blank(),axis.text.y=element_text(size=12),legend.position="none",panel.grid=element_blank(),plot.margin=margin(5,30,2,5))

pop_label_plot <- ggplot() + geom_text(data=pop_spans,aes(x=x_mid,y=0.5,label=display_name),size=5,fontface="bold") +
  geom_vline(data=pop_boundaries,aes(xintercept=x_end+0.5),color="black",linewidth=1) +
  scale_x_continuous(limits=c(0.5,max_x+3),expand=c(0,0)) + scale_y_continuous(limits=c(0,1),expand=c(0,0)) +
  coord_cartesian(clip="off") + theme_void() + theme(plot.margin=margin(0,30,10,5))

# AMOVA, fst and pairwise sig
pairwise_amova <- function(genind_obj, pop1, pop2, nrepet=9999) {
  keep <- pop(genind_obj) %in% c(pop1,pop2); gi_sub <- genind_obj[keep,]; pop(gi_sub) <- droplevels(pop(gi_sub))
  strata(gi_sub) <- data.frame(pops=pop(gi_sub))
  amova_r <- poppr.amova(gi_sub,~pops); amova_t <- ade4::randtest(amova_r,nrepet=nrepet)
  list(phi_st=amova_r$statphi$Phi[length(amova_r$statphi$Phi)], p_value=amova_t$pvalue[1])}
genI <- gl2gi(gl_mun); pops_gi <- genI$pop; strata(genI) <- data.frame(pops=pops_gi)
p.amova <- poppr.amova(genI,~pops); amova.pvalues <- ade4::randtest(p.amova,nrepet=9999)
overall_fst <- p.amova$statphi$Phi[length(p.amova$statphi$Phi)]; overall_pval <- amova.pvalues$pvalue[1]

fst_result <- gl.fst.pop(gl_mun,nboots=1000,percent=95,nclusters=1)
fst_matrix <- fst_result$Fsts; fst_matrix[is.na(fst_matrix)] <- t(fst_matrix)[is.na(fst_matrix)]
pop_sizes <- table(pop(gl_mun))

pval_matrix <- matrix(NA,nrow=nrow(fst_matrix),ncol=ncol(fst_matrix),dimnames=dimnames(fst_matrix))
for(i in 1:(nrow(fst_matrix)-1)){for(j in (i+1):ncol(fst_matrix)){
  p1 <- rownames(fst_matrix)[i]; p2 <- colnames(fst_matrix)[j]
  if(pop_sizes[p1]>=2 & pop_sizes[p2]>=2){
    res <- pairwise_amova(genI,p1,p2); pval_matrix[i,j] <- res$p_value; pval_matrix[j,i] <- res$p_value
    cat(p1,"vs",p2,": FST =",round(fst_matrix[i,j],4),", p =",round(res$p_value,4),"\n")}}}

get_upper_tri <- function(mat){mat[lower.tri(mat)]<-NA;return(mat)}
pop_name_map_mun <- c("MEX"="Mexico","NIC"="Nicaragua","PER"="Peru")
fst_matrix_labeled <- fst_matrix; rownames(fst_matrix_labeled) <- pop_name_map_mun[rownames(fst_matrix_labeled)]
colnames(fst_matrix_labeled) <- pop_name_map_mun[colnames(fst_matrix_labeled)]
pval_matrix_labeled <- pval_matrix; rownames(pval_matrix_labeled) <- pop_name_map_mun[rownames(pval_matrix)]
colnames(pval_matrix_labeled) <- pop_name_map_mun[colnames(pval_matrix)]
fst_melted <- melt(get_upper_tri(fst_matrix_labeled),na.rm=TRUE)
pval_melted <- melt(get_upper_tri(pval_matrix_labeled),na.rm=TRUE)
fst_melted <- fst_melted %>% left_join(pval_melted,by=c("Var1","Var2"),suffix=c("","_pval")) %>%
  mutate(sig=case_when(is.na(value_pval)~"",value_pval<0.001~"***",value_pval<0.01~"**",value_pval<0.05~"*",TRUE~""),
         label=paste0(round(value,4),sig))

fst_heatmap <- ggplot(fst_melted,aes(Var1,Var2,fill=value)) +
  geom_tile(color="white",linewidth=1) + geom_text(aes(label=label),size=5) +
  scale_fill_gradient(low="white",high="#1D8981",name=expression(italic("F")[ST]),limits=c(0,max(fst_melted$value,na.rm=TRUE)+0.001)) +
  coord_fixed() + labs(x="",y="",subtitle=bquote(atop("AMOVA"~italic("F")[ST]~"="~.(round(overall_fst,4))~"; p ="~.(round(overall_pval,3)),"* p < 0.05, ** p < 0.01, *** p < 0.001"))) +
  theme_minimal() + theme(axis.text.x=element_text(angle=45,hjust=1,size=12,face="bold"),axis.text.y=element_text(size=12,face="bold"),
                          plot.subtitle=element_text(size=11,hjust=0.5,face="bold"),panel.grid=element_blank(),legend.position="right")
ggsave("/Users/madelineeppley/Desktop/manta26pub/munkiana_Fst_heatmap.png", fst_heatmap, width=8, height=6, dpi=600)

# all layers together in final fig
layer1 <- range_map+theme(legend.position="right",plot.margin=margin(5,5,5,5))
layer2 <- wrap_plots(pca_plot_all,pca_plot_outliers,ncol=2)
layer3 <- wrap_plots(structure_all,structure_outlier,pop_label_plot,ncol=1,heights=c(1,1,0.3))
layer4 <- wrap_plots(fst_heatmap+theme(legend.position="right"),dapc_plot,ncol=2,widths=c(0.8,1.2))
final_figure <- wrap_plots(layer1,layer2,layer3,layer4,ncol=1,heights=c(0.8,1.1,0.8,1.0)) +
  plot_annotation(title="Mobula munkiana",
                  subtitle=paste0("All SNPs: n = ",format(n_snps_all,big.mark=","),"   |   Outlier SNPs: n = ",format(n_outliers,big.mark=",")),
                  tag_levels=list(c('A','B','C','D','E','','F','G')),
                  theme=theme(plot.title=element_text(size=18,face="bold.italic",hjust=0),plot.subtitle=element_text(size=14,hjust=0),plot.tag=element_text(size=16,face="bold")))
ggsave("/Users/madelineeppley/Desktop/manta26pub/munkiana_FINALFIGURE.png", final_figure, width=14, height=18, dpi=600)
