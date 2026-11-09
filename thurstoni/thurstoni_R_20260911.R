### M. thurstoni analysis to produce publication figure ####
# M. Eppley, v4. updated 10/20/26

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
dir.create(outdir, showWarnings=FALSE, recursive=TRUE)
setwd(outdir)  # so pcadapt can write positions.txt here

pop_colors <- c("IND"="#6DBAA4","MEX"="#ED9A6C","NIC"="#8C9FCB","ECU"="#DA8EC0","BYC"="#B8B8B8")
structure_colors <- c("#264653","#2a9d8f","#8ab17d","#e9c46a","#f4a261","#e76f51")

vcf_thur <- read.vcfR("/Users/madelineeppley/Desktop/manta/finalvcfs/thurstoni_final_v2.vcf")
popmap <- data.frame(
  sample = c("mthur_446","IN_10_MB_B","IN_12_MB_B",
             "MT_07_ELP","MT_04_ELP","MT_02_ELP","MT_05_ELP",
             "NIC_6003_E","NIC_6048","NIC_6003_F","NIC_6016_1",
             "MT_01_EC","MT_05_EC","MT_04_EC_conc","MT_07_EC_conc",
             "MT_10_EC","MT_09_EC_conc","MT_08_EC_conc","MT_03_EC",
             "BYC_RMO_40","BYC_RMO_47","BYC_RMO_48","BYC_RMU_16",
             "BYC_RMU_17","BYC_RMU_18","BYC_RMU_26","BYC_DESC_60","BYC_DESC_62"),
  pop = c(rep("IND",3),rep("MEX",4),rep("NIC",4),rep("ECU",8),rep("BYC",9)))

gl_thur <- vcfR2genlight(vcf_thur)
pop(gl_thur) <- popmap$pop[match(indNames(gl_thur),popmap$sample)]
n_snps_all <- nrow(vcf_thur@fix)
nInd(gl_thur) #28
print(table(pop(gl_thur))) # BYC 9, ECU 8, IND 3, MEX 4, NIC 4

pop_labels_df <- data.frame(pop=c("IND","MEX","NIC","ECU","BYC"),pop_label=c("India","Mexico","Nicaragua","Ecuador","Bycatch"))

# range wide sampling map
# IND=India, MEX=Pacific Mexico, NIC=Nicaragua, ECU=Ecuador, BYC=Bycatch
sample_coords <- data.frame(
  sample = popmap$sample,
  lat = c(10.5,11.0, 11.2,18.5,18.7,18.3,18.6, 11.2,11.5,11.3,11.0,
          -1.0,-0.8,-1.2,-0.5,-1.5,-1.3,-0.7,-0.9,
          -1.333,1.067,2.45,-9.217,-9.217,-9.217,-9.183,NA,NA),
  lon = c(72.5,73.0, 73.2,-105.0,-105.2,-104.8,-105.1, -87.0,-87.3,-87.1,-86.8,
          -80.5,-80.3,-80.7,-80.1,-81.0,-80.8,-80.2,-80.6,
          -92.5,-86.15,-90.333,-83.667,-83.667,-83.667,-85.85,NA,NA),
  pop = popmap$pop)
sample_coords_map <- sample_coords %>% filter(!is.na(lat)) %>%
  mutate(lon_pacific = ifelse(lon < 0, lon + 360, lon))
world_map <- map_data("world", wrap = c(0, 360))
pop_n <- popmap %>% group_by(pop) %>% summarise(n=n())

range_map <- ggplot() +
  geom_polygon(data=world_map, aes(x=long,y=lat,group=group), fill="gray90", color="gray70", linewidth=0.2)
if(exists("has_range") && has_range) range_map <- range_map + geom_polygon(data=thurstoni_range_df, aes(x=long,y=lat,group=group), fill="#6DBAA4", alpha=0.3, color=NA)
range_map <- range_map +
  geom_point(data=sample_coords_map, aes(x=lon_pacific,y=lat,fill=pop), shape=21, color="black", size=4, stroke=0.5) +
  scale_fill_manual(values=pop_colors, name="Sampling\nLocation", breaks=c("IND","MEX","NIC","ECU","BYC"),
                    labels=paste0(c("India","Mexico","Nicaragua","Ecuador","Bycatch")," (n=",pop_n$n[match(c("IND","MEX","NIC","ECU","BYC"),pop_n$pop)],")")) +
  coord_fixed(ratio = 1, xlim=c(30, 310), ylim=c(-50, 50), expand=FALSE) +
  scale_x_continuous(breaks=c(60,120,180,240,300), labels=c("60°E","120°E","180°","120°W","60°W")) +
  scale_y_continuous(breaks=c(-40,-20,0,20,40), labels=c("40°S","20°S","0°","20°N","40°N")) +
  theme_minimal() + theme(legend.position="right",legend.title=element_text(size=13,face="bold"),legend.text=element_text(size=12),
                          panel.grid.major=element_line(color="gray80",linewidth=0.3),panel.background=element_rect(fill="aliceblue",color=NA),
                          axis.title=element_blank(),axis.text=element_text(size=9)) + guides(fill=guide_legend(override.aes=list(size=5)))
ggsave("/Users/madelineeppley/Desktop/manta26pub/thurstoni_map.png",range_map,width=12,height=4,dpi=600)

# relatedness analysis
gds.fn <- "/Users/madelineeppley/Desktop/manta26pub/thur.gds"; showfile.gds(closeall=TRUE)
snpgdsVCF2GDS("/Users/madelineeppley/Desktop/manta/finalvcfs/thurstoni_final_v2.vcf",gds.fn,method="biallelic.only")
genofile <- snpgdsOpen(gds.fn); rel <- snpgdsIBDKING(genofile,autosome.only=FALSE); snpgdsClose(genofile)
kin_mat <- rel$kinship; sample_names_rel <- rel$sample.id
rownames(kin_mat) <- sample_names_rel; colnames(kin_mat) <- sample_names_rel
annotation_df <- data.frame(Population=popmap$pop[match(sample_names_rel,popmap$sample)]); rownames(annotation_df) <- sample_names_rel
relate <- pheatmap(kin_mat,labels_row=sample_names_rel,labels_col=sample_names_rel,annotation_row=annotation_df,annotation_col=annotation_df,
                   annotation_colors=list(Population=pop_colors),clustering_distance_rows="euclidean",clustering_distance_cols="euclidean",
                   main="Pairwise Genomic Relatedness, M. thurstoni",fontsize_row=6,fontsize_col=6,silent=TRUE)
write.csv(kin_mat, "/Users/madelineeppley/Desktop/manta26pub/thurstoni_kinship_matrix.csv")
diag(kin_mat) <- NA # set as NA so it doesnt return 0.5 as the highest relatedness
max_rel <- which(kin_mat == max(kin_mat, na.rm=TRUE), arr.ind=TRUE)[1,]
print(round(max(kin_mat, na.rm=TRUE), 4)) # highest relatedness 0.387
print(round(min(kin_mat, na.rm=TRUE), 4)) # lowest relatedness 0.0822
print(round(mean(kin_mat, na.rm=TRUE), 4)) # mean relatedness 0.1884
print(rownames(kin_mat)[max_rel[1]]) # BYC_RMU_17
print(colnames(kin_mat)[max_rel[2]]) # BYC_RMU_16
ggsave("/Users/madelineeppley/Desktop/manta26pub/thurstoni_relatedness.png",relate,width=12,height=10,dpi=600)

# pca using all snps
pca_all <- glPca(gl_thur,nf=10); eig_all <- pca_all$eig/sum(pca_all$eig)*100
pca_df_all <- data.frame(PC1=pca_all$scores[,1],PC2=pca_all$scores[,2],sample=indNames(gl_thur),pop=pop(gl_thur)) %>% left_join(pop_labels_df,by="pop")
centroid_all <- pca_df_all %>% group_by(pop) %>% summarise(PC1_cen=mean(PC1),PC2_cen=mean(PC2),n=n())
pca_df_all <- pca_df_all %>% left_join(centroid_all,by="pop")
pops_ell_all <- centroid_all %>% filter(n>=3) %>% pull(pop)
pca_plot_all <- ggplot(pca_df_all,aes(x=PC1,y=PC2)) +
  stat_ellipse(data=filter(pca_df_all,pop%in%pops_ell_all),aes(color=pop),type="norm",level=0.95,linewidth=1.2,show.legend=FALSE) +
  geom_segment(aes(xend=PC1_cen,yend=PC2_cen,color=pop),linewidth=0.5,alpha=0.6,show.legend=FALSE) +
  geom_point(aes(fill=pop,color=pop),size=4,alpha=0.8,shape=21,stroke=0.5) +
  scale_fill_manual(values=pop_colors,breaks=c("IND","MEX","NIC","ECU","BYC"),labels=c("India","Mexico","Nicaragua","Ecuador","Bycatch")) +
  scale_color_manual(values=pop_colors,breaks=c("IND","MEX","NIC","ECU","BYC"),labels=c("India","Mexico","Nicaragua","Ecuador","Bycatch")) +
  labs(x=paste0("PC1 (",round(eig_all[1],1),"%)"),y=paste0("PC2 (",round(eig_all[2],1),"%)"),subtitle="All SNPs") +
  theme_classic() + theme(legend.position="none",axis.text=element_text(size=13,color="black"),axis.title=element_text(size=14,face="bold"),
                          plot.subtitle=element_text(size=13,hjust=0.5),panel.border=element_rect(color="black",fill=NA,linewidth=1),
                          panel.grid.major=element_line(color="grey90",linewidth=0.3))
ggsave("/Users/madelineeppley/Desktop/manta26pub/thurstoni_PCA_allSNPs.png",
       pca_plot_all+labs(title="Mobula thurstoni")+theme(legend.position="bottom",legend.direction="horizontal",legend.title=element_blank(),
                                                         plot.title=element_text(size=13,face="italic",hjust=0.5))+guides(fill=guide_legend(override.aes=list(shape=21,size=5,alpha=0.8))),width=8,height=6,dpi=600)

# pcadapt
geno_thur <- read.pcadapt("/Users/madelineeppley/Desktop/manta/finalvcfs/thurstoni_final_v2.vcf",type="vcf")
obj <- pcadapt(geno_thur,K=1)
pvals <- obj$pvalues
qvals <- qvalue(pvals)$qvalues
alpha <- 0.01
outliers <- which(qvals < alpha)
n_outliers <- length(outliers)
n_outliers #13

# updated 10/20/26
padj <- p.adjust(pvals, method = "bonferroni")
alphab <- 0.1
outliers_bonf <- which(padj < alphab)
length(outliers_bonf) # 11 bonf outliers
outliers_bonf

vcf_out <- vcf_thur; vcf_out@fix <- vcf_thur@fix[outliers,,drop=FALSE]; vcf_out@gt <- vcf_thur@gt[outliers,,drop=FALSE]
write.vcf(vcf_out,"/Users/madelineeppley/Desktop/manta26pub/thurstoni_outliers.vcf.gz"); system(paste0("gunzip -f ","/Users/madelineeppley/Desktop/manta26pub/thurstoni_outliers.vcf.gz"))

# pca on outlier snps
vcf_outliers <- read.vcfR("/Users/madelineeppley/Desktop/manta26pub/thurstoni_outliers.vcf")
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
  scale_fill_manual(values=pop_colors,breaks=c("IND","MEX","NIC","ECU","BYC"),labels=c("India","Mexico","Nicaragua","Ecuador","Bycatch")) +
  scale_color_manual(values=pop_colors,breaks=c("IND","MEX","NIC","ECU","BYC"),labels=c("India","Mexico","Nicaragua","Ecuador","Bycatch")) +
  labs(x=paste0("PC1 (",round(eig_out[1],1),"%)"),y=paste0("PC2 (",round(eig_out[2],1),"%)"),subtitle="Outlier SNPs") +
  theme_classic() + theme(legend.position="none",axis.text=element_text(size=13,color="black"),axis.title=element_text(size=14,face="bold"),
                          plot.subtitle=element_text(size=13,hjust=0.5),panel.border=element_rect(color="black",fill=NA,linewidth=1),
                          panel.grid.major=element_line(color="grey90",linewidth=0.3))
ggsave("/Users/madelineeppley/Desktop/manta26pub/thurstoni_PCA_outlierSNPs.png",
       pca_plot_outliers+labs(title="Mobula thurstoni")+theme(legend.position="bottom",legend.direction="horizontal",legend.title=element_blank(),
                                                              plot.title=element_text(size=13,face="italic",hjust=0.5))+guides(fill=guide_legend(override.aes=list(shape=21,size=5,alpha=0.8))),width=8,height=6,dpi=600)

# dapc
dapc_result <- dapc(gl_thur,pop=pop(gl_thur),n.pca=4,n.da=4)
dapc_df <- data.frame(LD1=dapc_result$ind.coord[,1],LD2=dapc_result$ind.coord[,2],sample=indNames(gl_thur),pop=pop(gl_thur))
centroid_dapc <- dapc_df %>% group_by(pop) %>% summarise(LD1_cen=mean(LD1),LD2_cen=mean(LD2),n=n())
dapc_df <- dapc_df %>% left_join(centroid_dapc,by="pop")
pops_ell_dapc <- centroid_dapc %>% filter(n>=3) %>% pull(pop)
pct_dapc <- dapc_result$eig/sum(dapc_result$eig)*100
dapc_plot <- ggplot(dapc_df,aes(x=LD1,y=LD2)) +
  stat_ellipse(data=filter(dapc_df,pop%in%pops_ell_dapc),aes(color=pop),type="norm",level=0.95,linewidth=1.2,show.legend=FALSE) +
  geom_segment(aes(xend=LD1_cen,yend=LD2_cen,color=pop),linewidth=0.5,alpha=0.6,show.legend=FALSE) +
  geom_point(aes(fill=pop,color=pop),size=4,alpha=0.8,shape=21,stroke=0.5) +
  scale_fill_manual(values=pop_colors,breaks=c("IND","MEX","NIC","ECU","BYC"),labels=c("India","Mexico","Nicaragua","Ecuador","Bycatch")) +
  scale_color_manual(values=pop_colors,breaks=c("IND","MEX","NIC","ECU","BYC"),labels=c("India","Mexico","Nicaragua","Ecuador","Bycatch")) +
  labs(x=paste0("DAPC1 (",round(pct_dapc[1],1),"%)"),y=paste0("DAPC2 (",round(pct_dapc[2],1),"%)"),subtitle="All SNPs") +
  theme_classic() + theme(legend.position="none",axis.text=element_text(size=13,color="black"),axis.title=element_text(size=14,face="bold"),
                          panel.border=element_rect(color="black",fill=NA,linewidth=1),panel.grid.major=element_line(color="grey90",linewidth=0.3))
ggsave("/Users/madelineeppley/Desktop/manta26pub/thurstoni_DAPC.png",
       dapc_plot+labs(title="Mobula thurstoni")+theme(legend.position="bottom",legend.direction="horizontal",legend.title=element_blank(),
                                                      plot.title=element_text(size=13,face="italic",hjust=0.5))+guides(fill=guide_legend(override.aes=list(shape=21,size=5,alpha=0.8))),width=8,height=6,dpi=600)

# snmf ancestry analysis
# ancestry proportions with SNMF on all snps
## JUST RELOAD
# needed to double-check k values, re-load
project_all <- load.snmfProject("/Users/madelineeppley/Desktop/manta26pub/thurstoni.snmfProject")
project_outliers <- load.snmfProject("/Users/madelineeppley/Desktop/manta26pub/thurstoni_outliers.snmfProject")
for(k in 1:5) cat("K =", k, ":", min(cross.entropy(project_all, K=k)), "\n")
#K = 1 : 0.6466615
#K = 2 : 0.6623363
#K = 3 : 0.7291858
#K = 4 : 0.7488793
#K = 5 : 0.7856974
for(k in 1:5) cat("K =", k, ":", min(cross.entropy(project_outliers, K=k)), "\n")
#K = 1 : 0.1643111
#K = 2 : 0.1135736
#K = 3 : 0.07721085
#K = 4 : 0.4698943
#K = 5 : 0.2952678

#unlink("/Users/madelineeppley/Desktop/manta26pub/thurstoni.geno")
#unlink("/Users/madelineeppley/Desktop/manta26pub/thurstoni.snmfProject", recursive=TRUE)
#unlink("/Users/madelineeppley/Desktop/manta26pub/thurstoni.snmf", recursive=TRUE)
#vcf2geno("/Users/madelineeppley/Desktop/manta/finalvcfs/thurstoni_final_v2.vcf","/Users/madelineeppley/Desktop/manta26pub/thurstoni.geno")
#project_all <- snmf("/Users/madelineeppley/Desktop/manta26pub/thurstoni.geno",K=1:5,entropy=TRUE,repetitions=10,project="new")

# cross entropy plot with all snps
png("/Users/madelineeppley/Desktop/manta26pub/thurstoni_CE_allSNPs.png", width=800, height=600)
plot(project_all, col="blue", pch=19, cex=1.2, main="Cross-entropy: M. thurstoni - All SNPs")
dev.off()
print("all snps ce values")
for(k in 1:5) cat("K =", k, ":", min(cross.entropy(project_all, K=k)), "\n")

#unlink("/Users/madelineeppley/Desktop/manta26pub/thurstoni_outliers.geno")
#unlink("/Users/madelineeppley/Desktop/manta26pub/thurstoni_outliers.snmfProject", recursive=TRUE)
#unlink("/Users/madelineeppley/Desktop/manta26pub/thurstoni_outliers.snmf", recursive=TRUE)
#vcf2geno("/Users/madelineeppley/Desktop/manta26pub/thurstoni_outliers.vcf",
#        "/Users/madelineeppley/Desktop/manta26pub/thurstoni_outliers.geno")
#project_outliers <- snmf("/Users/madelineeppley/Desktop/manta26pub/thurstoni_outliers.geno",K=1:5,entropy=TRUE,repetitions=10,project="new")

# Cross-entropy plot - Outlier SNPs
png("/Users/madelineeppley/Desktop/manta26pub/thurstoni_CE_outlierSNPs.png", width=800, height=600)
plot(project_outliers, col="red", pch=19, cex=1.2, main="Cross-entropy: M. thurstoni - Outlier SNPs")
dev.off()
print("ce values from outlier snps")
for(k in 1:5) cat("K =", k, ":", min(cross.entropy(project_outliers, K=k)), "\n")

sample_names <- indNames(gl_thur); sample_pops <- pop(gl_thur)
pop_order <- c("IND","MEX","NIC","ECU","BYC"); pop_names_map <- c("IND"="India","MEX"="Mexico","NIC"="Nicaragua","ECU"="Ecuador","BYC"="Bycatch")
pop_label_df_str <- data.frame(Sample=sample_names,Pop=factor(sample_pops,levels=pop_order)) %>% arrange(Pop) %>% mutate(x_pos=1:n())
pop_boundaries <- pop_label_df_str %>% group_by(Pop) %>% summarise(x_end=max(x_pos)) %>% filter(x_end!=max(pop_label_df_str$x_pos))
pop_spans <- pop_label_df_str %>% group_by(Pop) %>% summarise(x_start=min(x_pos),x_end=max(x_pos),x_mid=(min(x_pos)+max(x_pos))/2,n=n()) %>%
  mutate(display_name=pop_names_map[as.character(Pop)])
max_x <- max(pop_label_df_str$x_pos)

ce_all <- cross.entropy(project_all,K=3); qmat_all <- Q(project_all,K=3,run=which.min(ce_all))
qdf_all <- as.data.frame(qmat_all); colnames(qdf_all) <- paste0("Cluster",1:3)
qdf_all$Sample <- sample_names; qdf_all$Pop <- factor(sample_pops,levels=pop_order)
qdf_all <- qdf_all %>% arrange(Pop) %>% mutate(x_pos=1:n())
qdf_all_long <- qdf_all %>% pivot_longer(cols=starts_with("Cluster"),names_to="Cluster",values_to="Proportion")
structure_all <- ggplot(qdf_all_long,aes(x=x_pos,y=Proportion,fill=Cluster)) +
  geom_bar(stat="identity",width=1,color=NA) + geom_vline(data=pop_boundaries,aes(xintercept=x_end+0.5),color="white",linewidth=0.8) +
  scale_fill_manual(values=structure_colors[1:3]) + scale_y_continuous(expand=c(0,0),breaks=c(0,0.5,1),labels=c("0.0","0.5","1.0"),position="left") +
  scale_x_continuous(expand=c(0,0),limits=c(0.5,max_x+3)) + annotate("text",x=max_x+2,y=0.5,label="K = 3",size=6,fontface="bold") +
  ggtitle("All SNPs") + coord_cartesian(clip="off") + theme_minimal() +
  theme(plot.title=element_text(size=12,face="bold",hjust=0),axis.text.x=element_blank(),axis.ticks.x=element_blank(),
        axis.title=element_blank(),axis.text.y=element_text(size=12),legend.position="none",panel.grid=element_blank(),plot.margin=margin(5,30,2,5))

ce_out <- cross.entropy(project_outliers,K=3); qmat_out_str <- Q(project_outliers,K=3,run=which.min(ce_out))
qdf_out_str <- as.data.frame(qmat_out_str); colnames(qdf_out_str) <- paste0("Cluster",1:3)
qdf_out_str$Sample <- sample_names; qdf_out_str$Pop <- factor(sample_pops,levels=pop_order)
qdf_out_str <- qdf_out_str %>% arrange(Pop) %>% mutate(x_pos=1:n())
qdf_out_long <- qdf_out_str %>% pivot_longer(cols=starts_with("Cluster"),names_to="Cluster",values_to="Proportion")
structure_outlier <- ggplot(qdf_out_long,aes(x=x_pos,y=Proportion,fill=Cluster)) +
  geom_bar(stat="identity",width=1,color=NA) + geom_vline(data=pop_boundaries,aes(xintercept=x_end+0.5),color="white",linewidth=0.8) +
  scale_fill_manual(values=structure_colors[c(3,4,5)]) + scale_y_continuous(expand=c(0,0),breaks=c(0,0.5,1),labels=c("0.0","0.5","1.0"),position="left") +
  scale_x_continuous(expand=c(0,0),limits=c(0.5,max_x+3)) + annotate("text",x=max_x+2,y=0.5,label="K = 3",size=6,fontface="bold") +
  ggtitle("Outlier SNPs") + coord_cartesian(clip="off") + theme_minimal() +
  theme(plot.title=element_text(size=12,face="bold",hjust=0),axis.text.x=element_blank(),axis.ticks.x=element_blank(),
        axis.title=element_blank(),axis.text.y=element_text(size=12),legend.position="none",panel.grid=element_blank(),plot.margin=margin(5,30,2,5))

pop_label_plot <- ggplot() + geom_text(data=pop_spans,aes(x=x_mid,y=0.5,label=display_name),size=5,fontface="bold") +
  geom_vline(data=pop_boundaries,aes(xintercept=x_end+0.5),color="black",linewidth=1) +
  scale_x_continuous(limits=c(0.5,max_x+3),expand=c(0,0)) + scale_y_continuous(limits=c(0,1),expand=c(0,0)) +
  coord_cartesian(clip="off") + theme_void() + theme(plot.margin=margin(0,30,10,5))

# amova FST
genI <- gl2gi(gl_thur); pops_gi <- genI$pop; strata(genI) <- data.frame(pops=pops_gi)
p.amova <- poppr.amova(genI,~pops); amova.pvalues <- ade4::randtest(p.amova,nrepet=9999)
overall_fst <- p.amova$statphi$Phi[length(p.amova$statphi$Phi)] # 0.002854827
# updated 10/20/26
print(amova.pvalues) # take row 3 from here for variation between pops
overall_pval <- amova.pvalues$pvalue[length(amova.pvalues$pvalue)]
overall_pval # was 1 with pvalue[1], 0.1369

# pairwise fst (all pops have n >= 2 here)
pop_sizes <- table(pop(gl_thur))
valid_pops <- names(pop_sizes[pop_sizes >= 2])
gl_valid <- gl_thur[pop(gl_thur) %in% valid_pops, ]
fst_result <- gl.fst.pop(gl_valid, nboots=1000, percent=95, nclusters=1)
fst_matrix <- fst_result$Fsts
fst_matrix[is.na(fst_matrix)] <- t(fst_matrix)[is.na(fst_matrix)]

# updated 10/20/26
pval_matrix <- fst_result$Pvalues
pval_matrix[is.na(pval_matrix)] <- t(pval_matrix)[is.na(pval_matrix)]
print("pairwise fst")
print(round(fst_matrix, 4))
print(pval_matrix)

## outlier fst
table(ploidy(gl_outliers))
ploidy(gl_outliers) <- 2 # set ploidy to 2 for all
nLoc(gl_outliers) # 13 outliers
pop(gl_outliers) <- popmap$pop[match(indNames(gl_outliers), popmap$sample)]
genI_out <- gl2gi(gl_outliers)
pops_gi_out <- genI_out$pop
strata(genI_out) <- data.frame(pops = pops_gi_out)

p.amova_out      <- poppr.amova(genI_out, ~pops)
amova.pvalues_out <- ade4::randtest(p.amova_out, nrepet = 9999)
overall_fst_out  <- p.amova_out$statphi$Phi[length(p.amova_out$statphi$Phi)]
# updated 10/20/26
print(amova.pvalues_out)
overall_pval_out <- amova.pvalues_out$pvalue[length(amova.pvalues_out$pvalue)]

# outlier amova fst and pval
round(overall_fst_out, 4)
round(overall_pval_out, 4)

# pairwise outlier FST
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
stopifnot(nLoc(gl_thur) == nrow(vcf_thur@fix))
neutral <- setdiff(seq_len(nLoc(gl_thur)), outliers)
gl_neutral <- gl_thur[, neutral]
nLoc(gl_neutral) #1178
genI_neu <- gl2gi(gl_neutral)
strata(genI_neu) <- data.frame(pops = genI_neu$pop)
p.amova_neu <- poppr.amova(genI_neu, ~pops)
amova.pvalues_neu <- ade4::randtest(p.amova_neu, nrepet = 9999)
print(amova.pvalues_neu)
overall_fst_neu  <- p.amova_neu$statphi$Phi[length(p.amova_neu$statphi$Phi)]
overall_pval_neu <- amova.pvalues_neu$pvalue[length(amova.pvalues_neu$pvalue)]
round(overall_fst_neu, 4) #0.0029
round(overall_pval_neu, 4) #0.1376

stars <- function(p) ifelse(is.na(p), "",
                            ifelse(p < 0.001, "***",
                                   ifelse(p < 0.01,  "**",
                                          ifelse(p < 0.05,  "*", ""))))
cat("  All-loci FST =", round(overall_fst, 4),  stars(overall_pval),     "\n")
cat("  Outlier FST  =", round(overall_fst_out,4), stars(overall_pval_out), "\n")

# updated 10/20/26 - poppr's default cutoff (0.05) removed 12/13 outliers and 1038/1178 neutral loci
amova_cutoff <- 0.5 # max proportion missing per locus, use the same value for every AMOVA and report it in methods

amova_filter <- function(gl, cutoff) {
  miss <- colMeans(is.na(as.matrix(gl)))
  gl <- gl[, miss <= cutoff]
  gl <- gl[rowSums(!is.na(as.matrix(gl))) > 0, ]  # drop inds with no data so poppr doesn't break
  cat("AMOVA using", nLoc(gl), "loci and", nInd(gl), "individuals\n")
  gl
}

gl_out_amova <- amova_filter(gl_outliers, amova_cutoff)
genI_out <- gl2gi(gl_out_amova)
strata(genI_out) <- data.frame(pops = genI_out$pop)
p.amova_out <- poppr.amova(genI_out, ~pops, cutoff = amova_cutoff)

# updated 10/20/26
p_label <- ifelse(overall_pval < 0.001, "< 0.001", paste("=", round(overall_pval, 3)))

get_upper_tri <- function(mat){mat[lower.tri(mat)]<-NA;return(mat)}
fst_matrix_labeled <- fst_matrix; rownames(fst_matrix_labeled) <- pop_names_map[rownames(fst_matrix_labeled)]
colnames(fst_matrix_labeled) <- pop_names_map[colnames(fst_matrix_labeled)]
pval_matrix_labeled <- pval_matrix; rownames(pval_matrix_labeled) <- pop_names_map[rownames(pval_matrix)]
colnames(pval_matrix_labeled) <- pop_names_map[colnames(pval_matrix)]
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
ggsave("/Users/madelineeppley/Desktop/manta26pub/thurstoni_Fst_heatmap.png",fst_heatmap,width=8,height=6,dpi=600)

# create the final figure :)
layer1 <- range_map+theme(legend.position="right",plot.margin=margin(5,5,5,5))
layer2 <- wrap_plots(pca_plot_all,pca_plot_outliers,ncol=2)
layer3 <- wrap_plots(structure_all,structure_outlier,pop_label_plot,ncol=1,heights=c(1,1,0.3))
layer4 <- wrap_plots(fst_heatmap+theme(legend.position="right"),dapc_plot,ncol=2,widths=c(0.8,1.2))
final_figure <- wrap_plots(layer1,layer2,layer3,layer4,ncol=1,heights=c(1.1,1.1,0.8,1.0)) +
  plot_annotation(title="Mobula thurstoni",
                  subtitle=paste0("All SNPs: n = ",format(n_snps_all,big.mark=","),"   |   Outlier SNPs: n = ",format(n_outliers,big.mark=",")),
                  tag_levels=list(c('A','B','C','D','E','','F','G')),
                  theme=theme(plot.title=element_text(size=18,face="bold.italic",hjust=0),plot.subtitle=element_text(size=14,hjust=0),plot.tag=element_text(size=16,face="bold")))
ggsave("/Users/madelineeppley/Desktop/manta26pub/thurstoni_FINALFIGURE.png",final_figure,width=14,height=18,dpi=600)
